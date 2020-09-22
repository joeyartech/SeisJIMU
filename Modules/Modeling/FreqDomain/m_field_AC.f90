module m_field
use m_global
use m_string
use m_mpienv
use m_arrayop
use m_sysio
use m_model, only:m
use m_geometry, only:geom
use m_computebox, only:cb, computebox_dealloc_model

! use, intrinsic :: ieee_arithmetic

    ! private
    ! public init_field_mumps, init_field_extmodel
    ! public field_matrix_factorize, field_RHS_substitute, field_RHS_substitute_adjoint
    ! public field_write

    complex :: oga

    real :: h2

    include 'cmumps_struc.h'
    type(cmumps_struc) :: mumps

    logical :: if_hicks=.true.

    !extended model
    real,dimension(:,:),allocatable :: buo !buoyancy = 1/rho
    complex,dimension(:,:),allocatable :: kpa !complex bulk modulus
    
    contains

    subroutine init_field_mumps
    
        if_hicks=get_setup_logical('IF_HICKS',default=.true.)

        !initialize MUMPS
        mumps%COMM = MPI_COMM_WORLD ! Parallel session 
        mumps%SYM =  0 ! unsymmetric matrix
        mumps%PAR =  1 ! master processor also involved in solving Ax=B
        
        mumps%JOB = -1 ! initialize mumps
        call cmumps(mumps)
                
        if(mpiworld%is_master) then
            mumps%ICNTL(20)  =  0     !dense RHS
            mumps%ICNTL(21)  =  0     !0->centralized solution 
            mumps%ICNTL(7)   =  get_setup_int('MUMPS_ORDERING','ICNTL(7)',default=7)  !choice of ordering for analysis (default: automatic choice)
            mumps%ICNTL(14)  =  get_setup_int('MUMPS_EXTRA_RAM','ICNTL(14)',default=60)  !percentage of increasing of estimated workspace for facto (default=20). No longer really used since MUMPS v4.9.X
            mumps%ICNTL(23)  =  get_setup_int('MUMPS_MAM_RAM','ICNTL(23)',default=760) !memory (MB) per processor, introduced since MUMPS v4.9.X
            mumps%KEEP(84)   =  get_setup_int('MUMPS_BLOCKING_FACTOR','KEEP(84)',default=16)   !blocking factor for multiple RHS
        end if

    end subroutine

    subroutine init_field_extmodel(freq)
        complex :: freq

        oga=freq*2.*r_pi
        h2=m%h*m%h

        if(mpiworld%is_master) then
            call alloc(kpa,[0,cb%nz-2+1],[0,cb%nx-2+1])
            call alloc(buo,[0,cb%nz-2+1],[0,cb%nx-2+1])

            ! !Kolsky-Futterman
            ! kpa=rho* vp*vp*(1.-0.5*c_i/qp)*(1.-0.5*c_i/qp) !dispersion-free

            ! vc=1./cmplx( (1./vp)+(1./(r_pi*vp*qp))*log(oga_ref/oga) , 0.5/(vp*qp) )
            ! kpa=rho* vc*vc !with dispersion


            kpa=cb%rho*cb%vp*cb%vp !ACoustic
            buo=1./cb%rho

    !         !left & right edges
    !         kpa(1:cb%nz,0)      =kpa(1:cb%nz,1)
    !         buo(1:cb%nz,0)      =buo(1:cb%nz,1)
    !         kpa(1:cb%nz,cb%nx+1)=kpa(1:cb%nz,cb%nx)       
    !         buo(1:cb%nz,cb%nx+1)=buo(1:cb%nz,cb%nx)
    ! 
    !         
    !         !top & bottom edges
    !         kpa(0,1:cb%nx)      =kpa(1,1:cb%nx)
    !         buo(0,1:cb%nx)      =buo(1,1:cb%nx)
    !         kpa(cb%nz+1,1:cb%nx)=kpa(cb%nz,1:cb%nx)       
    !         buo(cb%nz+1,1:cb%nx)=buo(cb%nz,1:cb%nx)
    ! 
    !         !4 corners
    !         kpa(0,0)            =kpa(1,1)
    !         buo(0,0)            =buo(1,1)
    !         
    !         kpa(0,cb%nx+1)      =kpa(1,cb%nx)
    !         buo(0,cb%nx+1)      =buo(1,cb%nx)
    !         
    !         kpa(cb%nz+1,0)      =kpa(cb%nz,1)
    !         buo(cb%nz+1,0)      =buo(cb%nz,1)
    ! 
    !         kpa(cb%nz+1,cb%nx+1)=kpa(cb%nz,cb%nx)
    !         buo(cb%nz+1,cb%nx+1)=buo(cb%nz,cb%nx)

            !save some memory
            call computebox_dealloc_model
            
        endif

    end subroutine


    subroutine field_matrix_factorize
        logical,save :: if_analyze=.true.

        if (mpiworld%is_master) then

            !allocate mumps table
            !estimate maximum size of matrix
            mumps%N=(cb%nz-2)*(cb%nx-2)
            n=9*mumps%N

            if(.not.associated(mumps%A))   allocate(mumps%A(n))
            if(.not.associated(mumps%IRN)) allocate(mumps%IRN(n))
            if(.not.associated(mumps%JCN)) allocate(mumps%JCN(n))
            
            !build matrix A
            call build_matrix(mumps%A,mumps%IRN,mumps%JCN,n,mumps%Nz)

        end if

        !matrix analysis
        !only needs to be done once
        if (if_analyze) then
            mumps%JOB = 1
            call cmumps(mumps)

            if_analyze=.false.
        end if

        !matrix factorization
        mumps%JOB = 2
        call cmumps(mumps)
  
    end subroutine

    subroutine field_RHS_substitute
        !fill RHS
        if (mpiworld%is_master) then
            if(.not.associated(mumps%RHS)) allocate(mumps%RHS((cb%nz-2)*(cb%nx-2)*geom%nsrc))
            
            call fill_RHS_source(mumps%RHS)
            mumps%ICNTL(9) = 1
            mumps%NRHS = geom%nsrc
            mumps%LRHS = mumps%N

        endif
        
        !solve Ax=b
        mumps%JOB = 3
        !if(mpiworld%is_master) call cpu_time(t1)
        call cmumps(mumps)
        !if(mpiworld%is_master) call cpu_time(t2)

        ! if(mpiworld%is_master) then
        !     write(*,*) 'time for mumps solution :',t2-t1, ' s' 
        ! endif

    end subroutine

    subroutine field_RHS_substitute_adjoint(dres)
        complex,dimension(geom%nsrc) :: dres

        !fill RHS
        if (mpiworld%is_master) then
            if(.not.associated(mumps%RHS)) allocate(mumps%RHS(cb%n*geom%nsrc))

            call fill_RHS_receiver(mumps%RHS,conjg(dres))

            mumps%icntl(9) = 0
            mumps%NRHS = geom%nsrc
            mumps%LRHS = mumps%N

        endif
        
        !solve A^Tx=b
        mumps%JOB = 3
        !if(mpiworld%is_master) call cpu_time(t1)
        call cmumps(mumps)
        !if(mpiworld%is_master) call cpu_time(t2)

        ! if(mpiworld%is_master) then
        !     write(*,*) 'time for mumps solution :',t2-t1, ' s' 
        ! endif

    end subroutine

    subroutine field_write(ifreq)
        complex,dimension(:),allocatable :: seismo

        if (mpiworld%is_master) then
            call alloc(seismo,geom%ntr)

            call get_RHS(mumps%RHS,seismo)

            open(20,file='synth_data_'//num2str(ifreq,'(04i)'),access='direct',recl=4*geom%ntr)
            write(20,rec=1) real(seismo(:))
            write(20,rec=2) aimag(seismo(:))
            close(20)
            
            open(20,file='wavefield',access='direct',recl=4*cb%n)
            write(20,rec=1) real(mumps%RHS(:))
            write(20,rec=2) aimag(mumps%RHS(:))
            close(20)

        end if

    end subroutine

    !========= set up Ax=b =================

    subroutine build_matrix(A,IRN,ICN,n,k)
        complex,dimension(n) :: A
        integer,dimension(n) :: IRN, ICN

        real,parameter :: damp=0.

        real,parameter :: w1=0.4382634,  w2=1.-w1
        real,parameter :: wm1=0.6287326, wm2=0.3712667, wm3=1.-wm1-wm2

        real,parameter :: quarter_wm2=0.25*wm2, quarter_wm3=0.25*wm3
        
        real,parameter :: az=0., bz=1. !; cz=0.; dz=0.; ez=0.; fz=0.; gz=0.; hz=0.
        real,parameter :: ax=1., bx=0. !; cx=0.; dx=0.; ex=0.; fx=0.; gx=0.; hx=0.

        complex :: dz0, dx0, dzp, dzm, dxp, dxm
        complex :: oga2

        w1_p2 = w1/(h2*h2)
        w2_p2 = w2/(h2*h2)
        quarter_w1_p2 = 0.25*w1_p2
        quarter_w2_p2 = 0.25*w2_p2

        oga2=oga*oga
        
        ! initialize sparse impedance matrix
        IRN=0
        ICN=0
        A=cmplx(0.,0.)

        k=0
        jr=0


        do ix=1,cb%nx-2
        do iz=1,cb%nz-2

            !interpolate buoyancy at half grid points
            !arithmetic averaging is used
            bpm=0.25*(buo(iz,ix)+buo(iz+1,ix)+buo(iz,ix-1)+buo(iz+1,ix-1))
            bmp=0.25*(buo(iz,ix)+buo(iz-1,ix)+buo(iz,ix+1)+buo(iz-1,ix+1))
            bpp=0.25*(buo(iz,ix)+buo(iz+1,ix)+buo(iz,ix+1)+buo(iz+1,ix+1))
            bmm=0.25*(buo(iz,ix)+buo(iz-1,ix)+buo(iz,ix-1)+buo(iz-1,ix-1))
            b00=buo(iz,ix)
            bp0=0.5*(buo(iz,ix)+buo(iz+1,ix))
            bm0=0.5*(buo(iz,ix)+buo(iz-1,ix))
            b0p=0.5*(buo(iz,ix)+buo(iz,ix+1))
            b0m=0.5*(buo(iz,ix)+buo(iz,ix-1))

            !damping functions for PML
            dz0=cb%damp_z(iz)
            dx0=cb%damp_x(ix)
            dxp=cb%damp_x_half(ix)
            dxm=cb%damp_x_half(ix-1)
            dzp=cb%damp_z_half(iz)
            dzm=cb%damp_z_half(iz-1)
            

            !Node  0 0
            jr=(ix-1)*cb%nz+iz
            jc=(ix-1)*cb%nz+iz

            k=k+1
            irn(k) = jr
            icn(k) = jc
            
            A(k) = wm1*oga2/kpa(iz,ix) &
                    +quarter_w1_p2*( &
                        -dx0*ax*( bmp*dxp +bpm*dxm +bpp*dxp +bmm*dxm)  &
                        +dx0*bx*( bmp*dzm +bpm*dzp -bpp*dzp -bmm*dzm)  &
                        +dz0*az*( bmp*dxp +bpm*dxm -bpp*dxp -bmm*dxm)  &
                        -dz0*bz*( bmp*dzm +bpm*dzp +bpp*dzp +bmm*dzm)) &
                    +w2_p2*( &
                        -dx0*ax*( b0p*dxp +b0m*dxm) &
                        -dz0*bz*( bp0*dzp +bm0*dzm))

            A(k)=A(k)+damp*4.*A(k)

            !Node  0+1
            jc=0
            if (ix<cb%nx) then
                jc=(ix-1+1)*cb%nz+iz
                
                k=k+1
                irn(k)=jr
                icn(k)=jc
                
                A(k)= quarter_wm2*oga2/kpa(iz,ix+1) &
                    +quarter_w1_p2*( &
                         dx0*ax*( bmp*dxp +bpp*dxp)   &
                        +dx0*bx*( bmp*dzm -bpp*dzp)   &
                        +dz0*az*(-bmp*dxp +bpp*dxp)   &
                        +dz0*bz*(-bmp*dzm -bpp*dzp) ) &
                    +w2_p2*( &
                              dx0*ax* b0p*dxp &
                        +0.25*dz0*az*(bp0*dx0 -bm0*dx0))

                A(k)=A(k)-damp*A(k)
            end if
           
            !Node  0-1
            if (ix>1) then
                jc=(ix-1-1)*cb%nz+iz

                k=k+1
                irn(k)=jr
                icn(k)=jc

                A(k)= quarter_wm2*oga2/kpa(iz,ix-1) &
                    +quarter_w1_p2*( &
                         dx0*ax*( bpm*dxm +bmm*dxm)  &
                        +dx0*bx*( bpm*dzp -bmm*dzm)  &
                        +dz0*az*(-bpm*dxm +bmm*dxm)  &
                        +dz0*bz*(-bpm*dzp -bmm*dzm)) &
                    +w2_p2*( &
                              dx0*ax*  b0m*dxm &
                        +0.25*dz0*az*(-bp0*dx0 +bm0*dx0) )

                A(k)=A(k)-damp*A(k)
            end if

            !Node -1 0
            if (iz > 1) then
                jc=(ix-1)*cb%nz+iz-1

                k=k+1
                irn(k)=jr
                icn(k)=jc

                A(k)= quarter_wm2*oga2/kpa(iz-1,ix) &
                    +quarter_w1_p2*( &
                         dx0*ax*(-bmp*dxp -bmm*dxm)  &
                        +dx0*bx*(-bmp*dzm +bmm*dzm)  &
                        +dz0*az*( bmp*dxp -bmm*dxm)  &
                        +dz0*bz*( bmp*dzm +bmm*dzm)) &
                    +w2_p2*( &
                              dz0*bz*  bm0*dzm &
                        +0.25*dx0*bx*(-b0p*dz0 +b0m*dz0) )

                A(k)=A(k)-damp*A(k)
            end if

            !Node +1 0
            if (iz < cb%nz) then
                jc=(ix-1)*cb%nz+iz+1

                k=k+1
                irn(k)=jr
                icn(k)=jc

                A(k)= quarter_wm2*oga2/kpa(iz+1,ix) &
                    +quarter_w1_p2*( &
                         dx0*ax*(-bpm*dxm -bpp*dxp)  &
                        +dx0*bx*(-bpm*dzp +bpp*dzp)  &
                        +dz0*az*( bpm*dxm -bpp*dxp)  &
                        +dz0*bz*( bpm*dzp +bpp*dzp)) &
                    +w2_p2*( &
                              dz0*bz*  bp0*dzp & 
                        +0.25*dx0*bx*(-b0m*dz0 +b0p*dz0) )

                A(k)=A(k)-damp*A(k)
            end if

            !Node -1+1
            if (iz > 1 .and. ix < cb%nx) then
                jc=(ix-1+1)*cb%nz+iz-1

                k=k+1
                irn(k)=jr
                icn(k)=jc
                
                A(k)= quarter_wm3*oga2/kpa(iz-1,ix+1) &
                    +quarter_w1_p2*( &
                         dx0*ax*bmp*dxp  &
                        -dx0*bx*bmp*dzm  &
                        -dz0*az*bmp*dxp  &
                        +dz0*bz*bmp*dzm) &
                    +quarter_w2_p2*( &
                        -dx0*bx*b0p*dz0 &
                        -dz0*az*bm0*dx0 )

            end if

            !Node +1-1
            if (iz < cb%nz .and. ix >1) then
                jc=(ix-1-1)*cb%nz+iz+1

                k=k+1
                irn(k)=jr
                icn(k)=jc
                
                A(k)= quarter_wm3*oga2/kpa(iz+1,ix-1) &
                    +quarter_w1_p2*( &
                         dx0*ax*bpm*dxm  &
                        -dx0*bx*bpm*dzp  &
                        -dz0*az*bpm*dxm  &
                        +dz0*bz*bpm*dzp) &
                    +quarter_w2_p2*( &
                        -dx0*bx*b0m*dz0 &
                        -dz0*az*bp0*dx0)

            end if

            !Node +1+1
            if (iz < cb%nz .and. ix < cb%nx) then
                jc=(ix-1+1)*cb%nz+iz+1

                k=k+1
                irn(k)=jr
                icn(k)=jc
                
                A(k)= quarter_wm3*oga2/kpa(iz+1,ix+1) &
                   +quarter_w1_p2*( &
                         dx0*ax*bpp*dxp  &
                        +dx0*bx*bpp*dzp  &
                        +dz0*az*bpp*dxp  &
                        +dz0*bz*bpp*dzp) &
                   +quarter_w2_p2*( &
                         dx0*bx*b0p*dz0 &
                        +dz0*az*bp0*dx0 )

            end if
           
            !Node -1-1
            if (iz > 1 .and. ix > 1) then
                jc=(ix-1-1)*cb%nz+iz-1

                k=k+1
                irn(k)=jr
                icn(k)=jc
                
                A(k)= quarter_wm3*oga2/kpa(iz-1,ix-1) &
                    +quarter_w1_p2*( &
                         dx0*ax*bmm*dxm  &
                        +dx0*bx*bmm*dzm  &
                        +dz0*az*bmm*dxm  &
                        +dz0*bz*bmm*dzm) &
                    +quarter_w2_p2*( &
                         dx0*bx*b0m*dz0 &
                        +dz0*az*bm0*dx0 )

            end if

        end do
        end do


        if(m%if_freesurface) then
            do ix=1,cb%nx
            do iz=1,cb%npml
                k=k+1  !PROBLEMATIC
                irn(k)= (ix-1)*cb%nz+iz
                icn(k)= (ix-1)*cb%nz+iz
                A(k)= oga2/kpa(iz,ix)
            enddo
            enddo
        endif

    end subroutine


    subroutine fill_RHS_source(RHS)
        complex,dimension(cb%nz,cb%nx,geom%nsrc) :: RHS

        complex :: source_term
        
        RHS=cmplx(0.,0.)

        source_term=1./(m%h*m%h)

        do i=1,geom%nsrc

            ifz=geom%src(i)%ifz+cb%npml; iz=geom%src(i)%iz+cb%npml; ilz=geom%src(i)%ilz+cb%npml
            ifx=geom%src(i)%ifx+cb%npml; ix=geom%src(i)%ix+cb%npml; ilx=geom%src(i)%ilx+cb%npml


            if(if_hicks) then
                RHS(ifz:ilz,ifx:ilx,i) = RHS(ifz:ilz,ifx:ilx,i) + source_term *geom%src(i)%interp_coeff
            else
                RHS(iz,ix,i) = RHS(iz,ix,i) + source_term
            endif

        enddo

    end subroutine

    subroutine fill_RHS_receiver(RHS,adjsource)
        complex,dimension(cb%nz,cb%nx,geom%nsrc) :: RHS
        complex,dimension(*) :: adjsource
        
        complex :: source_term

        RHS=cmplx(0.,0.)

        k=1

        do i=1,geom%nsrc
        do j=1,geom%src(i)%nrcv

            ifz=geom%src(i)%rcv(j)%ifz+cb%npml; iz=geom%src(i)%rcv(j)%iz+cb%npml; ilz=geom%src(i)%rcv(j)%ilz+cb%npml
            ifx=geom%src(i)%rcv(j)%ifx+cb%npml; ix=geom%src(i)%rcv(j)%ix+cb%npml; ilx=geom%src(i)%rcv(j)%ilx+cb%npml

            source_term=adjsource(k)/(m%h*m%h)
            k=k+1

            if(if_hicks) then
                ! select case (shot%rcv(ircv)%icomp)
                !     case (1) !pressure adjsource
                    RHS(ifz:ilz,ifx:ilx,i) = RHS(ifz:ilz,ifx:ilx,i) + source_term *geom%src(i)%rcv(j)%interp_coeff
                ! end select

            else
                ! select case (shot%rcv(ircv)%icomp)
                !     case (1) !pressure adjsource
                    !p[iz,ix,iy]
                    RHS(iz,ix,i) = RHS(iz,ix,i) + source_term
                ! end select

            endif

        enddo
        enddo

    end subroutine

    subroutine get_RHS(RHS,seismo)
        complex,dimension(cb%nz-2,cb%nx-2,geom%nsrc) :: RHS
        complex,dimension(*) :: seismo

        k=1

        do i=1,geom%nsrc
        do j=1,geom%src(i)%nrcv

            ifz=geom%src(i)%rcv(j)%ifz+cb%npml; iz=geom%src(i)%rcv(j)%iz+cb%npml; ilz=geom%src(i)%rcv(j)%ilz+cb%npml
            ifx=geom%src(i)%rcv(j)%ifx+cb%npml; ix=geom%src(i)%rcv(j)%ix+cb%npml; ilx=geom%src(i)%rcv(j)%ilx+cb%npml

            if(if_hicks) then
                seismo(k)=sum(RHS(ifz:ilz,ifx:ilx,i) *geom%src(i)%rcv(j)%interp_coeff )
            else
                seismo(k)= RHS(iz,ix,i)
            endif
            k=k+1

        enddo
        enddo

    end subroutine

end
