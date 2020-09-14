module m_field
use m_mpienv
use m_arrayop
use m_sysio
use m_model, only:m

! use, intrinsic :: ieee_arithmetic

    ! private fdcoeff_o4,fdcoeff_o8,c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z,c4x,c4y,c4z
    ! private k_x,k_y,k_z,npower,rcoef
    ! public
    
    ! !FD coeff
    ! real,dimension(2),parameter :: fdcoeff_o4 = [1.125,      -1./24.]
    ! real,dimension(4),parameter :: fdcoeff_o8 = [1225./1024, -245./3072., 49./5120., -5./7168]
    
    ! real :: c1x, c1y, c1z
    ! real :: c2x, c2y, c2z
    ! real :: c3x, c3y, c3z
    ! real :: c4x, c4y, c4z
    
    ! !CPML
    ! real,parameter :: k_x = 1.
    ! real,parameter :: k_y = 1.
    ! real,parameter :: k_z = 1.
    ! real,parameter :: npower=2.
    ! real,parameter :: rcoef=0.001
    
    ! !local models in computebox
    ! type t_localmodel
    !     sequence
    !     real,dimension(:,:,:),allocatable ::  buox, buoy, buoz, kpa, invkpa
    ! end type
    
    ! type(t_localmodel) :: lm
    
    ! !fields
    ! type t_field
    !     sequence
    !     !flat arrays
    !     real,dimension(:,:,:),allocatable :: vx,vy,vz,p
    !     real,dimension(:,:,:),allocatable :: cpml_dvx_dx,cpml_dvy_dy,cpml_dvz_dz !for cpml
    !     real,dimension(:,:,:),allocatable :: cpml_dp_dz, cpml_dp_dx, cpml_dp_dy
    ! end type
    
    ! !hicks interpolation
    ! logical :: if_hicks
    
    ! !info
    ! character(*),parameter :: waveeq_info='time-domain isotropic 2D/3D acoustic'
    ! character(*),parameter :: gradient_info='kpa-rho'
    ! integer,parameter :: ncorr=2

    real :: omega

    include 'cmumps_struc.h'
    type(cmumps_struc) :: mumps

    logical :: if_firstcall_mumps=.true.

    !extended model
    real,dimension(:,:),allocatable :: buo !buoyancy = 1/rho
    complex,dimension(:,:),allocatable :: kpa !complex bulk modulus
    
    contains

    subroutine init_field_mumps

        if_freesurf=get_setup_logical('IF_FREESURFACE','IF_FREESURF',default=.true.)

        !initialize MUMPS
        mumps%COMM = MPI_COMM_WORLD ! Parallel session 
        mumps%SYM =  0 ! Unsymmetric matrix
        mumps%PAR =  1 ! Host working
        mumps%JOB = -1 ! Initilization of MUMPS

        call cmumps(mumps)
        
        if(mpiworld%is_master) then
            mumps%icntl(20)  =  0     !dense rhs
            mumps%icntl(21)  =  0     !0->centralized solution 
            mumps%icntl(7)   =  get_setup_int('icntl_7',default=7)  !choice of ordering for analysis (default: automatic choice)
            mumps%icntl(14)  =  get_setup_int('icntl_14',default=60)  !percentage of increasing of estimated workspace for facto (default=20). No longer really used since MUMPS v4.9.X
            mumps%icntl(23)  =  get_setup_int('icntl_23',default=760) !memory (MB) per processor, introduced since MUMPS v4.9.X
            mumps%keep(84)   =  get_setup_int('keep_84',default=16)   !blocking factor for multiple RHS
        end if

        if_firstcall_mumps = .false.

    end subroutine

    subroutine init_field_extmodel(r_omega)

        omega=r_omega

        call alloc(kpa,[0,cb%n1+1],[0,cb%n2+1])
        call alloc(buo,[0,cb%n1+1],[0,cb%n2+1])

        ! !Kolsky-Futterman
        ! kpa=rho* vp*vp*(1.-0.5*c_i/qp)*(1.-0.5*c_i/qp) !dispersion-free

        ! vc=1./cmplx( (1./vp)+(1./(r_pi*vp*qp))*log(oga_ref/oga) , 0.5/(vp*qp) )
        ! kpa=rho* vc*vc !with dispersion

        kpa=m%rho*m%vp*m%vp !ACoustic

        do iz=1,cb%nz
            kpa(iz,0)=kpa(iz,1)
            kpa(iz,cb%nx+1)=kpa(iz,cb%nx)
        end do

        do ix=1,cb%nx
            kpa(0,ix)=kpa(1,ix)
            kpa(cb%nz+1,ix)=kpa(cb%nz,ix)
        end do

        kpa(0,0)            =kpa(1,1)
        kpa(0,cb%nx+1)      =kpa(1,cb%nx)
        kpa(cb%nz+1,0)      =kpa(cb%nz,1)
        kpa(cb%nz+1,cb%nx+1)=kpa(cb%nz,cb%nx)


        !buoyancy
        buo=1./m%rho

        buo(0,1:m%nx)     =buo(1,1:m%nx)
        buo(m%nz+1,1:m%nx)=buo(m%nz,1:m%nx)
        
        buo(1:m%nz,0)     =buo(1:m%nz,1)
        buo(1:m%nz,m%nx+1)=buo(1:m%nz,m%nx)
    
        buo(0,0)          =buo(1,1)
        buo(0,m%nx+1)     =buo(1,m%nx)
        buo(m%nz+1,0)     =buo(m%nz,1)
        buo(m%nz+1,m%nx+1)=buo(m%nz,m%nx)

    end subroutine

    subroutine field_matrix

        if (mpiworld%is_master) then

            !allocate mumps table
            !estimate maximum size of matrix
            n=9*cb%nz*cb%nx

            call alloc(mumps%A,n)
            call alloc(mumps%IRN,n)
            call alloc(mumps%JCN,n)

            !build matrix A
            call build_matrix(mumps%A,mumps%IRN,mumps%JCN,n)

        endif

        !LU decomposition
        if (mpiworld%is_master) then 
            mumps%n = cb%nz*cb%nx
            mumps%nz = nnz     
        end if

        if (if_firstcall_mumps==.false.) then
            call hud('MATRIX ANALYSIS')
            mumps%job = 1
            call cmumps(mumps)

            if_firstcall_mumps = .true.

        end if

        call hud('MATRIX FACTORIZATION')
        mumps%job = 2
        call cmumps(mumps)
  
    end subroutine



    subroutine field_solve(rhs)
        !fill RHS
        call cpu_time(t1)
        
        if (mpiworld%is_master) then
            allocate(mumps%rhs(cb%n1e*cb%n2e*acqui%nsrc))

            if(rhs=='source') then
                if(if_hicks) then
                    call fill_rhs_source_hicks(pbdir,acqui%nsrc,pbdir%id%rhs)
                else
                    call fill_rhs_source(mumps%rhs,cb%n1e,cb%n2e,cb%npml,m%h,acqui%nsrc,acqui%coord_src)
                endif

            elseif(rhs=='receiver') then
                if(if_hicks) then
                    call fill_rhs_receiver_hicks(pbdir,acqui%nsrc,acqui%ntrace,acqui%nrecmax,&
                        acqui%nbr_rec,acqui%associated_rec,inv%residual,mumps%rhs)
                else
                    call fill_rhs_receiver(mumps%rhs,cb%n1e,cb%n2e,cb%npml,&
                        m%h,acqui%nsrc,acqui%nrec_tot,acqui%ntrace,acqui%nrecmax,&
                        acqui%coord_rec,acqui%nbr_rec,acqui%associated_rec,inv%residual)
                endif

            endif

            mumps%nrhs = acqui%nsrc
            mumps%lrhs = mumps%n

        endif
        
        call cpu_time(t2)
        if(mpiworld%is_master) write(*,*) '1st part solve acqui :',t2-t1, ' s'

        if(rhs=='source') then ! src sources -> solve Ax=b
            mumps%icntl(9) = 1
        elseif(rhs=='receiver') then ! rec sources -> adjoint system -> solve A^tx=b
            mumps%icntl(9) = 0
        end if

        !solve 
        mumps%job = 3
        
        if(mpiworld%is_master) call cpu_time(t1)
        call cmumps(mumps)
        if(mpiworld%is_master) call cpu_time(t2)

        if(mpiworld%is_master) then
            write(*,*) 'time for mumps solution :',t2-t1, ' s' 
        endif

    end subroutine


    subroutine field_write
        complex,dimension(:),allocatable ::data_cal
        complex :: data_cal_k

        if (mpiworld%is_master) then
            allocate(dcal(acqui%ntrace)); dcal(:)=0.
            k=1

            do isrc = 1, acqui%nsrc
            do irec = 1, acqui%nbr_rec(isrc)

                irec_true = acqui%associated_rec(isrc,irec)
                
                if(if_hicks==.false.) then
                    !position of the source in the grid
                    xr=acqui%coord_rec(irec_true,1)
                    zr=acqui%coord_rec(irec_true,2)
                    ixr=nint(xr/pbdir%h)+1+pbdir%npml
                    izr=nint(zr/pbdir%h)+1+pbdir%npml

                    !ISO
                    dcal(k) = mumps%rhs((isrc-1)*(cb%n1e*cb%n2e) +(ixr-1)*cb%n1e+izr)
                    
                    !VTI: dcal(k) = mumps%rhs((isrc-1)*(cb%n1e*cb%n2e*2) +(ixr-1)*cb%n1e+izr)

                else
                    call hicks_extraction(pbdir,acqui,isrc,irec_true,dcal_k)
                    dcal(k)=data_cal_k
                endif

                k=k+1

            end do
            end do

            open(20,file='data_modeling',access='direct',recl=8*acqui%ntrace)
            write(20,rec=ifreq) dcal(:)
            close(20)
            write(*,*) 'data modeling writing ok'
              
            open(20,file='wavefield',access='direct',recl=4*pbdir%n1e*pbdir%n2e)
            write(20,rec=ifreq) real(mumps%rhs(1:cb%n1e*cb%n2e))
            close(20)

            write(*,*) 'wavefield modeling writing ok'

        end if

    end subroutine


    subroutine build_matrix(A,IRN,ICN,n)
        complex,dimension(n) :: A
        integer,dimension(n) :: IRN, ICN

        real,parameter :: damp=0.

        real,parameter :: w1=0.4382634,  w2=1.-w1
        real,parameter :: wm1=0.6287326, wm2=0.3712667, wm3=1.-wm1-wm2

        real,parameter :: quarter_wm2=0.25*wm2, quarter_wm3=0.25*wm3
        
        real,parameter :: az=0., bz=1. !; cz=0.; dz=0.; ez=0.; fz=0.; gz=0.; hz=0.
        real,parameter :: ax=1., bx=0. !; cx=0.; dx=0.; ex=0.; fx=0.; gx=0.; hx=0.

        complex :: dz0, dx0, dzp, dzm, dxp, dxm


        w1_p2 = w1/h2/h2
        w2_p2 = w2/h2/h2
        quarter_w1_p2 = 0.25*w1_p2
        quarter_w2_p2 = 0.25*w2_p2

        omega2=omega*omega
        
        ! initialize sparse impedance matrix
        IRN=0.
        ICN=0.
        A=cmplx(0.,0.)

        nnz=0
        jr=0


        do ix=1,cb%nx
        do iz=1,cb%nz

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
            jc=jr

            nnz=nnz+1
            irn(nnz) = jr
            icn(nnz) = jc
            
            A(nnz) = wm1*omega2/kpa(iz,ix) &
                    +quarter_w1_p2*( &
                        -dx0*ax*( bmp*dxp +bpm*dxm +bpp*dxp +bmm*dxm)  &
                        +dx0*bx*( bmp*dzm +bpm*dzp -bpp*dzp -bmm*dzm)  &
                        +dz0*az*( bmp*dxp +bpm*dxm -bpp*dxp -bmm*dxm)  &
                        -dz0*bz*( bmp*dzm +bpm*dzp +bpp*dzp +bmm*dzm)) &
                    +w2_p2*( &
                        -dx0*ax*( b0p*dxp +b0m*dxm) &
                        -dz0*bz*( bp0*dzp +bm0*dzm))

            A(nnz)=A(nnz)+damp*4.*A(nnz)

            !Node  0+1
            jc=0
            if (ix<cb%nx) then
                jc=(ix-1+1)*cb%nz+iz
                
                nnz=nnz+1
                irn(nnz)=jr
                icn(nnz)=jc
                
                A(nnz)= quarter_wm2*omega2/kpa(iz,ix+1) &
                    +quarter_w1_p2*( &
                         dx0*ax*( bmp*dxp +bpp*dxp)   &
                        +dx0*bx*( bmp*dzm -bpp*dzp)   &
                        +dz0*az*(-bmp*dxp +bpp*dxp)   &
                        +dz0*bz*(-bmp*dzm -bpp*dzp) ) &
                    +w2_p2*( &
                              dx0*ax* b0p*dxp &
                        +0.25*dz0*az*(bp0*dx  -bm0*dx) )

                A(nnz)=A(nnz)-damp*A(nnz)
            end if
           
            !Node  0-1
            if (ix>1) then
                jc=(ix-1-1)*cb%nz+iz

                nnz=nnz+1
                irn(nnz)=jr
                icn(nnz)=jc

                A(nnz)= quarter_wm2*omega2/kpa(iz,ix-1) &
                    +quarter_w1_p2*( &
                         dx0*ax*( bpm*dxm +bmm*dxm)  &
                        +dx0*bx*( bpm*dzp -bmm*dzm)  &
                        +dz0*az*(-bpm*dxm +bmm*dxm)  &
                        +dz0*bz*(-bpm*dzp -bmm*dzm)) &
                    +w2_p2*( &
                              dx0*ax*  b0m*dxm &
                        +0.25*dz0*az*(-bp0*dx0 +bm0*dx0) )

                A(nnz)=A(nnz)-damp*A(nnz)
            end if

            !Node -1 0
            if (iz > 1) then
                jc=(ix-1)*cb%nz+iz-1

                nnz=nnz+1
                irn(nnz)=jr
                icn(nnz)=jc

                A(nnz)= quarter_wm2*omega2/kpa(iz-1,ix) &
                    +quarter_w1_p2*( &
                         dx0*ax*(-bmp*dxp -bmm*dxm)  &
                        +dx0*bx*(-bmp*dzm +bmm*dzm)  &
                        +dz0*az*( bmp*dxp -bmm*dxm)  &
                        +dz0*bz*( bmp*dzm +bmm*dzm)) &
                    +w2_p2*( &
                              dz0*bz*  bm0*dzm &
                        +0.25*dx0*bx*(-b0p*dz0 +b0m*dz0) )

                A(nnz)=A(nnz)-damp*A(nnz)
            end if

            !Node +1 0
            if (iz < cb%nz) then
                jc=(ix-1)*cb%nz+iz+1

                nnz=nnz+1
                irn(nnz)=jr
                icn(nnz)=jc

                A(nnz)= quarter_wm2*omega2/kpa(iz+1,ix) &
                    +quarter_w1_p2*( &
                         dx0*ax*(-bpm*dxm -bpp*dxp)  &
                        +dx0*bx*(-bpm*dzp +bpp*dzp)  &
                        +dz0*az*( bpm*dxm -bpp*dxp)  &
                        +dz0*bz*( bpm*dzp +bpp*dzp)) &
                    +w2_p2*( &
                              dz0*bz*  bp0*dzp & 
                        +0.25*dx0*bx*(-b0m*dz0 +b0p*dz0) )

                A(nnz)=A(nnz)-damp*A(nnz)
            end if

            !Node -1+1
            if (iz > 1 .and. ix < cb%nx) then
                jc=(ix-1+1)*cb%nz+iz-1

                nnz=nnz+1
                irn(nnz)=jr
                icn(nnz)=jc
                
                A(nnz)= quarter_wm3*omega2/kpa(iz-1,ix+1) &
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

                nnz=nnz+1
                irn(nnz)=jr
                icn(nnz)=jc
                
                A(nnz)= quarter_wm3*omega2/kpa(iz+1,ix-1) &
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

                nnz=nnz+1
                irn(nnz)=jr
                icn(nnz)=jc
                
                A(nnz)= quarter_wm3*omega2/kpa(iz+1,ix+1) &
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

                nnz=nnz+1
                irn(nnz)=jr
                icn(nnz)=jc
                
                A(nnz)= quarter_wm3*omega2/kpa(iz-1,ix-1) &
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


        if(if_freesurf) then
            do ix=1,cb%nx
            do iz=1,cb%npml then
                nnz=nnz+1
                irn(nnz)= (ix-1)*cb%nz+iz
                icn(nnz)= (ix-1)*cb%nz+iz
                A(nnz)= omega2/kpa(iz,ix)
            enddo
            enddo
        endif

    end subroutine



!     !========= use before propagation =================
!     subroutine field_print_info
!         !modeling method
!         call hud('WaveEq : Time-domain isotropic 2D/3D ACoustic system')
!         call hud('1st-order Velocity-Stress formulation')
!         call hud('Staggered-grid Finite-difference method')
!         call hud('Cartesian O(x4,t2) stencil')
        
!         !stencil constant
!         if(mpiworld%is_master) then
!             write(*,*) 'Coeff:',fdcoeff_o4
!         endif
!         c1x=fdcoeff_o4(1)/m%dx; c1y=fdcoeff_o4(1)/m%dy; c1z=fdcoeff_o4(1)/m%dz
!         c2x=fdcoeff_o4(2)/m%dx; c2y=fdcoeff_o4(2)/m%dy; c2z=fdcoeff_o4(2)/m%dz
                
!     end subroutine
    
!     subroutine check_model  !not required
        
!     end
    
!     subroutine check_discretization
!         !grid dispersion condition
!         if (5.*m%cell_diagonal > cb%velmin/shot%src%fpeak/2.) then  !FDTDo4 rule
!             write(*,*) 'WARNING: Shot# '//shot%cindex//' can have grid dispersion!'
!             write(*,*) 'Shot# '//shot%cindex//' 5*dx, velmin, fpeak:',5.*m%cell_diagonal, cb%velmin,shot%src%fpeak
!         endif
        
!         !CFL condition
!         cfl = cb%velmax*shot%src%dt*m%cell_inv_diagonal*sum(abs(fdcoeff_o4))! ~0.494 (3D); ~0.606 (2D)
!         if(mpiworld%is_master) write(*,*) 'CFL value:',CFL
        
!         if(cfl>1.) then
!             write(*,*) 'ERROR: CFL > 1 on shot# '//shot%cindex//'!'
!             write(*,*) 'Shot# '//shot%cindex//' velmax, dt, 1/dx:',cb%velmax,shot%src%dt,m%cell_inv_diagonal
!             stop
!         endif
        
!     end subroutine
    
!     !========= use inside propagation =================
    
!     subroutine init_field_localmodel
    
!         call alloc(lm%buox,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],initialize=.false.)
!         call alloc(lm%buoy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],initialize=.false.)
!         call alloc(lm%buoz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],initialize=.false.)
!         call alloc(lm%kpa, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],initialize=.false.)
!         call alloc(lm%invkpa,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],initialize=.false.)
        
!         lm%kpa(:,:,:)=cb%vp(:,:,:)*cb%vp(:,:,:)*cb%rho(:,:,:)
!         lm%invkpa(:,:,:)=1./lm%kpa(:,:,:)

!         lm%buoz(cb%ifz,:,:)=1./cb%rho(cb%ifz,:,:)
!         lm%buox(:,cb%ifx,:)=1./cb%rho(:,cb%ifx,:)
!         lm%buoy(:,:,cb%ify)=1./cb%rho(:,:,cb%ify)

!         do iz=cb%ifz+1,cb%ilz
!             lm%buoz(iz,:,:)=0.5/cb%rho(iz,:,:)+0.5/cb%rho(iz-1,:,:)
!         enddo
        
!         do ix=cb%ifx+1,cb%ilx
!             lm%buox(:,ix,:)=0.5/cb%rho(:,ix,:)+0.5/cb%rho(:,ix-1,:)
!         enddo

!         do iy=cb%ify+1,cb%ily
!             lm%buoy(:,:,iy)=0.5/cb%rho(:,:,iy)+0.5/cb%rho(:,:,iy-1)
!         enddo
        
!     end subroutine
    
!     subroutine init_field(f)
!         type(t_field) :: f
        
!         call alloc(f%vx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
!         call alloc(f%vy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
!         call alloc(f%vz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
!         call alloc(f%p, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        
!         call alloc(f%cpml_dvx_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
!         call alloc(f%cpml_dvy_dy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
!         call alloc(f%cpml_dvz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
!         call alloc(f%cpml_dp_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
!         call alloc(f%cpml_dp_dy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
!         call alloc(f%cpml_dp_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        
!         !hicks point interpolation
!         if_hicks=get_setup_logical('IF_HICKS',default=.true.)
        
!     end subroutine
    
!     subroutine check_field(f,name)
!         type(t_field) :: f
!         character(*) :: name
        
!         if(mpiworld%is_master) write(*,*) name//' sample values:',minval(f%p),maxval(f%p)
        
!         if(any(.not. ieee_is_finite(f%vz))) then
!             write(*,*) 'ERROR: '//name//' values become Infinity on Shot# '//shot%cindex//' !!'
!             stop
!         endif
!         if(any(ieee_is_nan(f%vz))) then
!            write(*,*) 'ERROR: '//name//' values become NaN on Shot# '//shot%cindex//' !!'
!            stop
!         endif
        
!     end subroutine
    
!     subroutine field_cpml_reinitialize(f)
!         type(t_field) :: f
!         f%cpml_dvx_dx=0.
!         f%cpml_dvy_dy=0.
!         f%cpml_dvz_dz=0.
!         f%cpml_dp_dx=0.
!         f%cpml_dp_dy=0.
!         f%cpml_dp_dz=0.
!     end subroutine
    
!     subroutine write_field(iunit,f)
!         type(t_field) :: f
!         write(iunit) f%vz
!     end subroutine
    
!     !========= forward propagation =================
!     !WE: du_dt = MDu
!     !u=[vx vy vz p]^T, p=(sxx+syy+szz)/3=sxx+syy+szz
!     !  [diag3(b)  0  ]    [0   0   0   dx]
!     !M=[   0     kpa ], D=|0   0   0   dy|
!     !  (b=1/rho)          |0   0   0   dz|
!     !                     [dx  dy  dz  0 ]
!     !
!     !Discretization (staggered grid in space and time):
!     !  (grid index)     (real index)
!     !  s(iz,ix,iy):=   s[iz,ix,iy]^it+0.5,it+1.5,...
!     ! vx(iz,ix,iy):=  vx[iz,ix-0.5,iy]^it,it+1,...
!     ! vy(iz,ix,iy):=  vy[iz,ix,iy-0.5]^it,it+1,...
!     ! vz(iz,ix,iy):=  vy[iz-0.5,ix,iy]^it,it+1,...
    
!     !add RHS to v^it
!     subroutine put_velocities(time_dir,it,w,f)
!         integer :: time_dir,it
!         real :: w
!         type(t_field) :: f
        
!         ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
!         ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
!         ify=shot%src%ify; iy=shot%src%iy; ily=shot%src%ily
        
!         source_term=time_dir*w
        
!         if(if_hicks) then
!             select case (shot%src%icomp)
!                 case (2)
!                 f%vx(ifz:ilz,ifx:ilx,ify:ily) = f%vx(ifz:ilz,ifx:ilx,ify:ily) + source_term *lm%buox(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff
                
!                 case (3)
!                 f%vy(ifz:ilz,ifx:ilx,ify:ily) = f%vy(ifz:ilz,ifx:ilx,ify:ily) + source_term *lm%buoy(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff
                
!                 case (4)
!                 f%vz(ifz:ilz,ifx:ilx,ify:ily) = f%vz(ifz:ilz,ifx:ilx,ify:ily) + source_term *lm%buoz(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff
                
!             end select
            
!         else
!             select case (shot%src%icomp)
!                 case (2) !horizontal x force on vx[iz,ix-0.5,iy]
!                 f%vx(iz,ix,iy) = f%vx(iz,ix,iy) + source_term*lm%buox(iz,ix,iy)
                
!                 case (3) !horizontal y force on vy[iz,ix,iy-0.5]
!                 f%vy(iz,ix,iy) = f%vy(iz,ix,iy) + source_term*lm%buoy(iz,ix,iy)
                
!                 case (4) !vertical force     on vz[iz-0.5,ix,iy]
!                 f%vz(iz,ix,iy) = f%vz(iz,ix,iy) + source_term*lm%buoz(iz,ix,iy)
                
!             end select
            
!         endif
        
!     end subroutine
    
!     !v^it -> v^it+1 by FD of s^it+0.5
!     subroutine update_velocities(time_dir,it,f,bl)
!         integer :: time_dir, it
!         type(t_field) :: f
!         integer,dimension(6) :: bl
        
!         ifz=bl(1)+2
!         ilz=bl(2)-1
!         ifx=bl(3)+2
!         ilx=bl(4)-1
!         ify=bl(5)+2
!         ily=bl(6)-1
        
!         if(m%if_freesurface) ifz=max(ifz,1)
        
!         if(m%is_cubic) then
!             call fd3d_flat_velocities(f%vx,f%vy,f%vz,f%p,                         &
!                                       f%cpml_dp_dx,f%cpml_dp_dy,f%cpml_dp_dz,     &
!                                       lm%buox,lm%buoy,lm%buoz,                    &
!                                       ifz,ilz,ifx,ilx,ify,ily,time_dir*shot%src%dt)
!         else
!             call fd2d_flat_velocities(f%vx,f%vz,f%p,                      &
!                                       f%cpml_dp_dx,f%cpml_dp_dz,          &
!                                       lm%buox,lm%buoz,                    &
!                                       ifz,ilz,ifx,ilx,time_dir*shot%src%dt)
!         endif
        
!         !apply free surface boundary condition if needed
!         !free surface is located at [1,ix,iy] level
!         !so symmetric mirroring: vz[0.5]=vz[1.5], ie. vz(1,ix,iy)=vz(2,ix,iy) -> dp(1,ix,iy)=0.
!         if(m%if_freesurface) then
!             f%vz(1,:,:)=f%vz(2,:,:)
! !             !$omp parallel default (shared)&
! !             !$omp private(ix,iy,i)
! !             !$omp do schedule(dynamic)
! !             do iy=ify,ily
! !             do ix=ifx,ilx
! !                 i=(1-cb%ifz) + (ix-cb%ifx)*nz + (iy-cb%ify)*nz*nx +1 !iz=1,ix,iy
! !                 
! !                 f%vz(i)=f%vz(i+1)
! !             enddo
! !             enddo
! !             !$omp enddo
! !             !$omp end parallel
!         endif
        
!     end subroutine
    
!     !add RHS to s^it+0.5
!     subroutine put_stresses(time_dir,it,w,f)
!         integer :: time_dir
!         real :: w
!         type(t_field) :: f
        
!         ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
!         ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
!         ify=shot%src%ify; iy=shot%src%iy; ily=shot%src%ily
        
!         source_term=time_dir*w
        
!         if(if_hicks) then
!             select case (shot%src%icomp)
!                 case (1)
!                 f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) + source_term *shot%src%interp_coeff
!             end select
            
!         else
!             select case (shot%src%icomp)
!                 case (1) !explosion on s[iz,ix,iy]
!                 f%p(iz,ix,iy) = f%p(iz,ix,iy) + source_term
!             end select
            
!         endif
        
!     end subroutine
    
!     !s^it+0.5 -> s^it+1.5 by FD of v^it+1
!     subroutine update_stresses(time_dir,it,f,bl)
!         integer :: time_dir, it
!         type(t_field) :: f
!         integer,dimension(6) :: bl
        
!         ifz=bl(1)+1
!         ilz=bl(2)-2
!         ifx=bl(3)+1
!         ilx=bl(4)-2
!         ify=bl(5)+1
!         ily=bl(6)-2
        
!         if(m%if_freesurface) ifz=max(ifz,1)
        
!         if(m%is_cubic) then
!             call fd3d_flat_stresses(f%vx,f%vy,f%vz,f%p,                         &
!                                     f%cpml_dvx_dx,f%cpml_dvy_dy,f%cpml_dvz_dz,  &
!                                     lm%kpa,                                     &
!                                     ifz,ilz,ifx,ilx,ify,ily,time_dir*shot%src%dt)
!         else
!             call fd2d_flat_stresses(f%vx,f%vz,f%p,                      &
!                                     f%cpml_dvx_dx,f%cpml_dvz_dz,        &
!                                     lm%kpa,                             &
!                                     ifz,ilz,ifx,ilx,time_dir*shot%src%dt)
!         endif
        
!         !apply free surface boundary condition if needed
!         !free surface is located at [1,ix,iy] level
!         !so explicit boundary condition: p(1,ix,iy)=0
!         !and antisymmetric mirroring: p(0,ix,iy)=-p(2,ix,iy) -> vz(2,ix,iy)=vz(1,ix,iy)
!         if(m%if_freesurface) then
!             f%p(1,:,:)=0.
!             f%p(0,:,:)=-f%p(2,:,:)
! !             !$omp parallel default (shared)&
! !             !$omp private(ix,iy,i)
! !             !$omp do schedule(dynamic)
! !             do iy=ify,ily
! !             do ix=ifx,ilx
! !                 i=(1-cb%ifz) + (ix-cb%ifx)*nz + (iy-cb%ify)*nz*nx +1 !iz=1,ix,iy 
! !                 
! !                 f%p(i)=0.
! !                 
! !                 f%p(i-1)=-f%p(i+1)
! !             enddo
! !             enddo
! !             !$omp enddo
! !             !$omp end parallel
!         endif
        
!     end subroutine
    
!     !get v^it+1 or s^it+1.5
!     subroutine get_field(f,seismo)
!         type(t_field) :: f
!         real,dimension(*) :: seismo
        
!         do ircv=1,shot%nrcv
!             ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
!             ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
!             ify=shot%rcv(ircv)%ify; iy=shot%rcv(ircv)%iy; ily=shot%rcv(ircv)%ily
            
!             if(if_hicks) then
!                 select case (shot%rcv(ircv)%icomp)
!                     case (1)
!                     seismo(ircv)=sum(  f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
!                     case (2)
!                     seismo(ircv)=sum( f%vx(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
!                     case (3)
!                     seismo(ircv)=sum( f%vy(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
!                     case (4)
!                     seismo(ircv)=sum( f%vz(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
!                 end select
                
!             else
!                 select case (shot%rcv(ircv)%icomp)
!                     case (1) !p[iz,ix,iy]
!                     seismo(ircv)= f%p(iz,ix,iy)
!                     case (2) !vx[iz,ix-0.5,iy]
!                     seismo(ircv)=f%vx(iz,ix,iy)
!                     case (3) !vy[iz,ix,iy-0.5]
!                     seismo(ircv)=f%vy(iz,ix,iy)
!                     case (4) !vz[iz-0.5,ix,iy]
!                     seismo(ircv)=f%vz(iz,ix,iy)
!                 end select
                
!             endif
            
!         enddo
        
!     end subroutine
    
    
!     !========= adjoint propagation =================
!     !WE:        du_dt   = MDu   + src
!     !->        Ndu_dt   = Du    + Nsrc, N=M^-1
!     !Adjoint:  Ndv_dt^T = D^Tv  + adjsrc
!     !->        -dv_dt   =-MDv   + Mr
!     !     (reverse time) (negative derivatives)
!     !
!     !Discrete form (staggered grid in space and time, 2D as example):
!     ! N  [vx^it+1  ] [vx^it    ]   [ 0     0    dx(-)][vx^it+1  ]
!     !----|vz^it+1  |-|vz^it    | = | 0     0    dz(-)||vz^it+1  |  +Nsrc
!     ! dt [ s^it+1.5] [ s^it+0.5]   [dx(+) dz(+)  0   ][ s^it+0.5]
!     !dx(-):=a1(s(ix)-s(ixm1))+a2(s(ixp1)-s(ixm2))  (O(x4))
!     !dx(+):=a1(v(ixp1)-v(ix))+a2(v(ixp2)-v(ixm1))  (O(x4))
!     !Adjoint:
!     ! N  [vx^it    ] [vx^it+1  ]   [ 0       0      dx(+)^T][vx^it+1  ]
!     !----|vz^it    |-|vz^it+1  | = | 0       0      dz(+)^T||vz^it+1  |  +Nsrc
!     ! dt [ s^it+0.5] [ s^it+1.5]   [dx(-)^T dz(-)^T  0     ][ s^it+0.5]
!     !dx(+)^T=a1(s(ixm1)-s(ix))+a2(s(ixm2)-s(ixp1))=-dx(-)
!     !dx(-)^T=a1(v(ix)-v(ixp1))+a2(v(ixm1)-v(ixp2))=-dx(+)
!     !-dx,-dz can be regarded as -dt, thus allowing to use same code with flip of dt sign.
!     !transposes of time derivatives centered on v^it+0.5 and s^it+1
!     !so v^it+1   - v^it     => v^it - v^it+1
!     !   s^it+1.5 - s^it+0.5 => s^it+0.5 - s^it+1.5
!     !
!     !In each time step:
!     !d=RGAsrc, src:source term, A:put source, G:propagator, R:sampling at receivers
!     !d=[vx vy vz p]^T, R=[I  0   0 ], A=[diag3(b) 0], src=[fx fy fz p]^T
!     !                    |0 2/3 1/3|    [   0     1]
!     !                    [   0    1]
!     !Adjoint: src=A^T N G^T M R^T d
!     !M R^T:put adjoint sources, 
!     !G^T:  adjoint propagator,
!     !A^T N:get adjoint fields
    
    
!     !add RHS to s^it+1.5
!     subroutine put_stresses_adjoint(time_dir,it,adjsource,f)
!         integer :: time_dir
!         real,dimension(*) :: adjsource
!         type(t_field) :: f
        
!         do ircv=1,shot%nrcv
!             ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
!             ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
!             ify=shot%rcv(ircv)%ify; iy=shot%rcv(ircv)%iy; ily=shot%rcv(ircv)%ily
            
!             tmp=adjsource(ircv)
            
!             if(if_hicks) then
!                 select case (shot%rcv(ircv)%icomp)
!                     case (1) !pressure adjsource
!                     f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) &
!                         +tmp* lm%kpa(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
!                 end select
                
!             else
!                 select case (shot%rcv(ircv)%icomp)
!                     case (1) !pressure adjsource
!                     !p[iz,ix,iy]
!                     f%p(iz,ix,iy) = f%p(iz,ix,iy)  &
!                         +tmp* lm%kpa(iz,ix,iy)
!                 end select
                
!             endif
            
!         enddo
        
!     end subroutine
    
!     !s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
!     subroutine update_stresses_adjoint(time_dir,it,f,bl)
!         integer :: time_dir, it
!         type(t_field) :: f
!         integer,dimension(6) :: bl
        
!         !same code with -dt
!         call update_stresses(time_dir,it,f,bl)
        
!     end subroutine
    
!     !add RHS to v^it+1
!     subroutine put_velocities_adjoint(time_dir,it,adjsource,f)
!         integer :: time_dir
!         real,dimension(*) :: adjsource
!         type(t_field) :: f
        
!         do ircv=1,shot%nrcv
!             ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
!             ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
!             ify=shot%rcv(ircv)%ify; iy=shot%rcv(ircv)%iy; ily=shot%rcv(ircv)%ily
            
!             tmp=adjsource(ircv)
            
!             if(if_hicks) then
!                 select case (shot%rcv(ircv)%icomp)
!                     case (2) !horizontal x adjsource
!                     f%vx(ifz:ilz,ifx:ilx,ify:ily) = f%vx(ifz:ilz,ifx:ilx,ify:ily) + tmp*lm%buox(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
                    
!                     case (3) !horizontal y adjsource
!                     f%vy(ifz:ilz,ifx:ilx,ify:ily) = f%vy(ifz:ilz,ifx:ilx,ify:ily) + tmp*lm%buoy(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
                    
!                     case (4) !horizontal z adjsource
!                     f%vz(ifz:ilz,ifx:ilx,ify:ily) = f%vz(ifz:ilz,ifx:ilx,ify:ily) + tmp*lm%buoz(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
!                 end select
                
!             else
!                 select case (shot%rcv(ircv)%icomp)
!                     case (2) !horizontal x adjsource
!                     !vx[ix-0.5,iy,iz]
!                     f%vx(iz,ix,iy) = f%vx(iz,ix,iy) + tmp*lm%buox(iz,ix,iy)
                    
!                     case (3) !horizontal y adjsource
!                     !vy[ix,iy-0.5,iz]
!                     f%vy(iz,ix,iy) = f%vy(iz,ix,iy) + tmp*lm%buoy(iz,ix,iy)
                    
!                     case (4) !horizontal z adjsource
!                     !vz[ix,iy,iz-0.5]
!                     f%vz(iz,ix,iy) = f%vz(iz,ix,iy) + tmp*lm%buoz(iz,ix,iy)
!                 end select
                
!             endif
            
!         enddo
        
!     end subroutine
    
!     !v^it+1 -> v^it by FD^T of s^it+0.5
!     subroutine update_velocities_adjoint(time_dir,it,f,bl)
!         integer :: time_dir, it
!         type(t_field) :: f
!         integer,dimension(6) :: bl
        
!         !same code with -dt
!         call update_velocities(time_dir,it,f,bl)
        
!     end subroutine
    
!     !get v^it or s^it+0.5
!     subroutine get_field_adjoint(f,w)
!         type(t_field) :: f
!         real w
        
!         ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
!         ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
!         ify=shot%src%ify; iy=shot%src%iy; ily=shot%src%ily
        
!         if(if_hicks) then
!             select case (shot%src%icomp)
!                 case (1)
!                 w = sum(lm%invkpa(ifz:ilz,ifx:ilx,ify:ily)*f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
!                 case (2)
!                 w = sum(     f%vx(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
!                 case (3)
!                 w = sum(     f%vy(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
!                 case (4)
!                 w = sum(     f%vz(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
!             end select
            
!         else
!             select case (shot%src%icomp)
!                 case (1) !p[iz,ix,iy]
!                 w=lm%invkpa(iz,ix,iy)*f%p(iz,ix,iy)
                
!                 case (2) !vx[iz,ix-0.5,iy]
!                 w=f%vx(iz,ix,iy)
                
!                 case (3) !vy[iz,ix,iy-0.5]
!                 w=f%vy(iz,ix,iy)
                
!                 case (4) !vz[iz-0.5,ix,iy]
!                 w=f%vz(iz,ix,iy)
!             end select
            
!         endif
        
!     end subroutine
    
!     !========= for wavefield correlation ===================   
!     !gkpa = -1/kpa2    dp_dt \dot adjp
!     !     = -1/kpa \nabla.v  \dot adjp
!     !v^it+1, p^it+0.5, adjp^it+0.5
!     !
!     !grho = dv_dt      \dot adjv
!     !     = b \nabla p \dot adjv
!     !p^it+0.5, v^it, adjv^it
!     !use (v[i+1]+v[i])/2 to approximate v[i+0.5], so is adjv
    
!     subroutine field_correlation_moduli(it,sf,rf,sb,rb,corr)
!         type(t_field) :: sf,rf
!         integer,dimension(6) :: sb,rb
!         real,dimension(cb%mz,cb%mx,cb%my,ncorr) :: corr
        
!         !nonzero only when sf touches rf
!         ifz=max(sb(1),rb(1),2)
!         ilz=min(sb(2),rb(2),cb%mz-2)
!         ifx=max(sb(3),rb(3),2)
!         ilx=min(sb(4),rb(4),cb%mx-2)
!         ify=max(sb(5),rb(5),2)
!         ily=min(sb(6),rb(6),cb%my-2)
        
!         if(m%is_cubic) then
!             call corr3d_flat_moduli(sf%vx,sf%vy,sf%vz,rf%p,&
!                                     corr(:,:,:,1),         &
!                                     ifz,ilz,ifx,ilx,ify,ily)
!         else
!             call corr2d_flat_moduli(sf%vx,sf%vz,rf%p,&
!                                     corr(:,:,:,1),   &
!                                     ifz,ilz,ifx,ilx)
!         endif
        
!     end subroutine
    
!     subroutine field_correlation_density(it,sf,rf,sb,rb,corr)
!         type(t_field) :: sf,rf
!         integer,dimension(6) :: sb,rb
!         real,dimension(cb%mz,cb%mx,cb%my,ncorr) :: corr
        
!         !nonzero only when sf touches rf
!         ifz=max(sb(1),rb(1),2)
!         ilz=min(sb(2),rb(2),cb%mz-2)
!         ifx=max(sb(3),rb(3),2)
!         ilx=min(sb(4),rb(4),cb%mx-2)
!         ify=max(sb(5),rb(5),2)
!         ily=min(sb(6),rb(6),cb%my-2)
        
!         if(m%is_cubic) then
!             call corr3d_flat_density(sf%p,rf%vx,rf%vy,rf%vz,&
!                                      corr(:,:,:,2),         &
!                                      ifz,ilz,ifx,ilx,ify,ily)
!         else
!             call corr2d_flat_density(sf%p,rf%vx,rf%vz,&
!                                      corr(:,:,:,2),   &
!                                      ifz,ilz,ifx,ilx)
!         endif
        
!     end subroutine
    
!     subroutine field_correlation_scaling(corr)
!         real,dimension(cb%mz,cb%mx,cb%my,ncorr) :: corr
        
!         corr(:,:,:,1)=corr(:,:,:,1) * (-lm%invkpa(1:cb%mz,1:cb%mx,1:cb%my))

!         corr(:,:,:,2)=corr(:,:,:,2) / cb%rho(1:cb%mz,1:cb%mx,1:cb%my)

!         corr(1,:,:,:) = corr(2,:,:,:)
!         corr(cb%mz-1,:,:,:) = corr(cb%mz-2,:,:,:)
!         corr(cb%mz,  :,:,:) = corr(cb%mz-2,:,:,:)

!         corr(:,1,:,:) = corr(:,2,:,:)
!         corr(:,cb%mx-1,:,:) = corr(:,cb%mx-2,:,:)
!         corr(:,cb%mx  ,:,:) = corr(:,cb%mx-2,:,:)

!         if(m%is_cubic) then
!             corr(:,:,1,:) = corr(:,2,:,:)
!             corr(:,:,cb%my-1,:) = corr(:,:,cb%my-2,:)
!             corr(:,:,cb%my  ,:) = corr(:,:,cb%my-2,:)
!         endif
        
!         !set unit of gkpa to be [m3], grho to be [m5/s2]
!         !such that after multiplied by (kpa_max-kpa_min) or (rho_max-rho_min) (will be done in m_parameterization.f90)
!         !the unit of parameter update is [Nm], same as Lagrangian
!         !and the unit of gradient scaling factor is [1/N/m] (in m_scaling.f90)
!         !therefore parameters become unitless
!         corr=corr*m%cell_volume*shot%src%dt
        
!     end subroutine
    
!     !========= Finite-Difference on flattened arrays ==================
    
!     subroutine fd3d_flat_velocities(vx,vy,vz,p,                       &
!                                     cpml_dp_dx,cpml_dp_dy,cpml_dp_dz, &
!                                     buox,buoy,buoz,                   &
!                                     ifz,ilz,ifx,ilx,ify,ily,dt)
!         real,dimension(*) :: vx,vy,vz,p
!         real,dimension(*) :: cpml_dp_dx,cpml_dp_dy,cpml_dp_dz
!         real,dimension(*) :: buox,buoy,buoz
        
!         nz=cb%nz
!         nx=cb%nx
!         ny=cb%ny
        
!         dp_dx=0.;dp_dy=0.;dp_dz=0.
        
!         !$omp parallel default (shared)&
!         !$omp private(iz,ix,iy,i,&
!         !$omp         izm2_ix_iy,izm1_ix_iy,iz_ix_iy,izp1_ix_iy,&
!         !$omp         iz_ixm2_iy,iz_ixm1_iy,iz_ixp1_iy,&
!         !$omp         iz_ix_iym2,iz_ix_iym1,iz_ix_iyp1,&
!         !$omp         dp_dx,dp_dy,dp_dz)
!         !$omp do schedule(dynamic) collapse(2)
!         do iy=ify,ily
!         do ix=ifx,ilx
        
!             !dir$ simd
!             do iz=ifz,ilz
                
!                 i=(iz-cb%ifz)+(ix-cb%ifx)*nz+(iy-cb%ify)*nz*nx+1
                
!                 izm2_ix_iy=i-2  !iz-2,ix,iy
!                 izm1_ix_iy=i-1  !iz-1,ix,iy
!                 iz_ix_iy  =i    !iz,ix,iy
!                 izp1_ix_iy=i+1  !iz+1,ix,iy
                
!                 iz_ixm2_iy=i  -2*nz  !iz,ix-2,iy
!                 iz_ixm1_iy=i    -nz  !iz,ix-1,iy
!                 iz_ixp1_iy=i    +nz  !iz,ix+1,iy
                
!                 iz_ix_iym2=i  -2*nz*nx  !iz,ix,iy-2
!                 iz_ix_iym1=i    -nz*nx  !iz,ix,iy-1
!                 iz_ix_iyp1=i    +nz*nx  !iz,ix,iy+1
                
!                 dp_dx= c1x*(p(iz_ix_iy)-p(iz_ixm1_iy)) +c2x*(p(iz_ixp1_iy)-p(iz_ixm2_iy))
!                 dp_dy= c1y*(p(iz_ix_iy)-p(iz_ix_iym1)) +c2y*(p(iz_ix_iyp1)-p(iz_ix_iym2))
!                 dp_dz= c1z*(p(iz_ix_iy)-p(izm1_ix_iy)) +c2z*(p(izp1_ix_iy)-p(izm2_ix_iy))
                
!                 !cpml
!                 cpml_dp_dx(iz_ix_iy)= cb%b_x_half(ix)*cpml_dp_dx(iz_ix_iy) + cb%a_x_half(ix)*dp_dx
!                 cpml_dp_dy(iz_ix_iy)= cb%b_y_half(iy)*cpml_dp_dy(iz_ix_iy) + cb%a_y_half(iy)*dp_dy
!                 cpml_dp_dz(iz_ix_iy)= cb%b_z_half(iz)*cpml_dp_dz(iz_ix_iy) + cb%a_z_half(iz)*dp_dz

!                 dp_dx=dp_dx*cb%kappa_x_half(ix) + cpml_dp_dx(iz_ix_iy)
!                 dp_dy=dp_dy*cb%kappa_y_half(iy) + cpml_dp_dy(iz_ix_iy)
!                 dp_dz=dp_dz*cb%kappa_z_half(iz) + cpml_dp_dz(iz_ix_iy)
                
!                 !velocity
!                 vx(iz_ix_iy)=vx(iz_ix_iy) + dt*buox(iz_ix_iy)*dp_dx
!                 vy(iz_ix_iy)=vy(iz_ix_iy) + dt*buoy(iz_ix_iy)*dp_dy
!                 vz(iz_ix_iy)=vz(iz_ix_iy) + dt*buoz(iz_ix_iy)*dp_dz
                
!             enddo
            
!         enddo
!         enddo
!         !$omp end do
!         !$omp end parallel
        
!     end subroutine
    
!     subroutine fd2d_flat_velocities(vx,vz,p,              &
!                                     cpml_dp_dx,cpml_dp_dz,&
!                                     buox,buoz,            &
!                                     ifz,ilz,ifx,ilx,dt)
!         real,dimension(*) :: vx,vz,p
!         real,dimension(*) :: cpml_dp_dx,cpml_dp_dz
!         real,dimension(*) :: buox,buoz
        
!         nz=cb%nz
!         nx=cb%nx
        
!         dp_dx=0.; dp_dz=0.
        
!         !$omp parallel default (shared)&
!         !$omp private(iz,ix,i,&
!         !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
!         !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
!         !$omp         dp_dx,dp_dz)
!         !$omp do schedule(dynamic)
!         do ix=ifx,ilx
        
!             !dir$ simd
!             do iz=ifz,ilz
                
!                 i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1
                
!                 izm2_ix=i-2  !iz-2,ix
!                 izm1_ix=i-1  !iz-1,ix
!                 iz_ix  =i    !iz,ix
!                 izp1_ix=i+1  !iz+1,ix
                
!                 iz_ixm2=i  -2*nz  !iz,ix-2
!                 iz_ixm1=i    -nz  !iz,ix-1
!                 iz_ixp1=i    +nz  !iz,ix+1
                
!                 dp_dx= c1x*(p(iz_ix)-p(iz_ixm1)) +c2x*(p(iz_ixp1)-p(iz_ixm2))
!                 dp_dz= c1z*(p(iz_ix)-p(izm1_ix)) +c2z*(p(izp1_ix)-p(izm2_ix))
                
!                 !cpml
!                 cpml_dp_dx(iz_ix)= cb%b_x_half(ix)*cpml_dp_dx(iz_ix) + cb%a_x_half(ix)*dp_dx
!                 cpml_dp_dz(iz_ix)= cb%b_z_half(iz)*cpml_dp_dz(iz_ix) + cb%a_z_half(iz)*dp_dz

!                 dp_dx=dp_dx*cb%kappa_x_half(ix) + cpml_dp_dx(iz_ix)
!                 dp_dz=dp_dz*cb%kappa_z_half(iz) + cpml_dp_dz(iz_ix)
                
!                 !velocity
!                 vx(iz_ix)=vx(iz_ix) + dt*buox(iz_ix)*dp_dx
!                 vz(iz_ix)=vz(iz_ix) + dt*buoz(iz_ix)*dp_dz
                
!             enddo
            
!         enddo
!         !$omp end do
!         !$omp end parallel
        
!     end subroutine
    
!     subroutine fd3d_flat_stresses(vx,vy,vz,p,                         &
!                                   cpml_dvx_dx,cpml_dvy_dy,cpml_dvz_dz,&
!                                   kpa,                                &
!                                   ifz,ilz,ifx,ilx,ify,ily,dt)
!         real,dimension(*) :: vx,vy,vz,p
!         real,dimension(*) :: cpml_dvx_dx,cpml_dvy_dy,cpml_dvz_dz
!         real,dimension(*) :: kpa
        
!         nz=cb%nz
!         nx=cb%nx
!         ny=cb%ny
        
!         dvx_dx=0.;dvy_dy=0.;dvz_dz=0.
        
!         !$omp parallel default (shared)&
!         !$omp private(iz,ix,iy,i,&
!         !$omp         izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
!         !$omp         iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
!         !$omp         iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
!         !$omp         dvx_dx,dvy_dy,dvz_dz)
!         !$omp do schedule(dynamic) collapse(2)
!         do iy=ify,ily
!         do ix=ifx,ilx
        
!             !dir$ simd
!             do iz=ifz,ilz
            
!                 i=(iz-cb%ifz)+(ix-cb%ifx)*nz+(iy-cb%ify)*nz*nx+1
                
!                 izm1_ix_iy=i-1  !iz-1,ix,iy
!                 iz_ix_iy  =i    !iz,ix,iy
!                 izp1_ix_iy=i+1  !iz+1,ix,iy
!                 izp2_ix_iy=i+2  !iz+2,ix,iy
                
!                 iz_ixm1_iy=i       -nz  !iz,ix-1,iy
!                 iz_ixp1_iy=i       +nz  !iz,ix+1,iy
!                 iz_ixp2_iy=i     +2*nz  !iz,ix+2,iy
                
!                 iz_ix_iym1=i    -nz*nx  !iz,ix,iy-1
!                 iz_ix_iyp1=i    +nz*nx  !iz,ix,iy+1
!                 iz_ix_iyp2=i  +2*nz*nx  !iz,ix,iy+2
                
!                 dvx_dx= c1x*(vx(iz_ixp1_iy)-vx(iz_ix_iy))  +c2x*(vx(iz_ixp2_iy)-vx(iz_ixm1_iy))
!                 dvy_dy= c1y*(vy(iz_ix_iyp1)-vy(iz_ix_iy))  +c2y*(vy(iz_ix_iyp2)-vy(iz_ix_iym1))
!                 dvz_dz= c1z*(vz(izp1_ix_iy)-vz(iz_ix_iy))  +c2z*(vz(izp2_ix_iy)-vz(izm1_ix_iy))
                
!                 !cpml
!                 cpml_dvx_dx(iz_ix_iy)=cb%b_x(ix)*cpml_dvx_dx(iz_ix_iy)+cb%a_x(ix)*dvx_dx
!                 cpml_dvy_dy(iz_ix_iy)=cb%b_y(iy)*cpml_dvy_dy(iz_ix_iy)+cb%a_y(iy)*dvy_dy
!                 cpml_dvz_dz(iz_ix_iy)=cb%b_z(iz)*cpml_dvz_dz(iz_ix_iy)+cb%a_z(iz)*dvz_dz

!                 dvx_dx=dvx_dx*cb%kappa_x(ix) + cpml_dvx_dx(iz_ix_iy)
!                 dvy_dy=dvy_dy*cb%kappa_y(iy) + cpml_dvy_dy(iz_ix_iy)
!                 dvz_dz=dvz_dz*cb%kappa_z(iz) + cpml_dvz_dz(iz_ix_iy)
                
!                 !pressure
!                 p(iz_ix_iy) = p(iz_ix_iy) + dt * kpa(iz_ix_iy)*(dvx_dx+dvy_dy+dvz_dz)
                
!             enddo
            
!         enddo
!         enddo
!         !$omp enddo 
!         !$omp end parallel
        
!     end subroutine
    
!     subroutine fd2d_flat_stresses(vx,vz,p,                &
!                                   cpml_dvx_dx,cpml_dvz_dz,&
!                                   kpa,                    &
!                                   ifz,ilz,ifx,ilx,dt)
!         real,dimension(*) :: vx,vz,p
!         real,dimension(*) :: cpml_dvx_dx,cpml_dvz_dz
!         real,dimension(*) :: kpa
        
!         nz=cb%nz
!         nx=cb%nx
        
!         dvx_dx=0.;dvz_dz=0.
        
!         !$omp parallel default (shared)&
!         !$omp private(iz,ix,i,&
!         !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
!         !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
!         !$omp         dvx_dx,dvz_dz)
!         !$omp do schedule(dynamic)
!         do ix=ifx,ilx
        
!             !dir$ simd
!             do iz=ifz,ilz
            
!                 i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1
                
!                 izm1_ix=i-1  !iz-1,ix
!                 iz_ix  =i    !iz,ix
!                 izp1_ix=i+1  !iz+1,ix
!                 izp2_ix=i+2  !iz+2,ix
                
!                 iz_ixm1=i  -nz  !iz,ix-1
!                 iz_ixp1=i  +nz  !iz,ix+1
!                 iz_ixp2=i  +2*nz !iz,ix+2
                
!                 dvx_dx= c1x*(vx(iz_ixp1)-vx(iz_ix))  +c2x*(vx(iz_ixp2)-vx(iz_ixm1))
!                 dvz_dz= c1z*(vz(izp1_ix)-vz(iz_ix))  +c2z*(vz(izp2_ix)-vz(izm1_ix))
                
!                 !cpml
!                 cpml_dvx_dx(iz_ix)=cb%b_x(ix)*cpml_dvx_dx(iz_ix)+cb%a_x(ix)*dvx_dx
!                 cpml_dvz_dz(iz_ix)=cb%b_z(iz)*cpml_dvz_dz(iz_ix)+cb%a_z(iz)*dvz_dz

!                 dvx_dx=dvx_dx*cb%kappa_x(ix) + cpml_dvx_dx(iz_ix)
!                 dvz_dz=dvz_dz*cb%kappa_z(iz) + cpml_dvz_dz(iz_ix)
                
!                 !pressure
!                 p(iz_ix) = p(iz_ix) + dt * kpa(iz_ix)*(dvx_dx+dvz_dz)
                
!             enddo
            
!         enddo
!         !$omp enddo 
!         !$omp end parallel
        
!     end subroutine
    
!     subroutine corr3d_flat_moduli(sf_vx,sf_vy,sf_vz,rf_p,&
!                                   corr,                  &
!                                   ifz,ilz,ifx,ilx,ify,ily)
!         real,dimension(*) :: sf_vx,sf_vy,sf_vz,rf_p
!         real,dimension(*) :: corr
        
!         nz=cb%nz
!         nx=cb%nx
        
!         dvx_dx=0.
!         dvy_dx=0.
!         dvz_dz=0.
        
!         dsp=0.
!          rp=0.
        
!         !$omp parallel default (shared)&
!         !$omp private(iz,ix,iy,i,j,&
!         !$omp         izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
!         !$omp         iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
!         !$omp         iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
!         !$omp         dvx_dx,dvy_dy,dvz_dz,&
!         !$omp         dsp,rp)
!         !$omp do schedule(dynamic) collapse(2)
!         do iy=ify,ily
!         do ix=ifx,ilx
        
!             !dir$ simd
!             do iz=ifz,ilz
                
!                 i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(iy-cb%ify)*cb%nz*cb%nx+1 !field has boundary layers
!                 j=(iz-1)     +(ix-1)     *cb%mz+(iy-1)     *cb%mz*cb%mx+1 !corr has no boundary layers
                
!                 izm1_ix_iy=i-1  !iz-1,ix,iy
!                 iz_ix_iy  =i    !iz,ix,iy
!                 izp1_ix_iy=i+1  !iz+1,ix,iy
!                 izp2_ix_iy=i+2  !iz+2,ix,iy
                
!                 iz_ixm1_iy=i    -nz  !iz,ix-1,iy
!                 iz_ixp1_iy=i    +nz  !iz,ix+1,iy
!                 iz_ixp2_iy=i  +2*nz  !iz,ix+2,iy
                
!                 iz_ix_iym1=i    -nz*nx  !iz,ix,iy-1
!                 iz_ix_iyp1=i    +nz*nx  !iz,ix,iy+1
!                 iz_ix_iyp2=i  +2*nz*nx  !iz,ix,iy+2
                
!                 dvx_dx = c1x*(sf_vx(iz_ixp1_iy)-sf_vx(iz_ix_iy)) +c2x*(sf_vx(iz_ixp2_iy)-sf_vx(iz_ixm1_iy))
!                 dvy_dy = c1y*(sf_vy(iz_ix_iyp1)-sf_vy(iz_ix_iy)) +c2y*(sf_vy(iz_ix_iyp2)-sf_vy(iz_ix_iym1))
!                 dvz_dz = c1z*(sf_vz(izp1_ix_iy)-sf_vz(iz_ix_iy)) +c2z*(sf_vz(izp2_ix_iy)-sf_vz(izm1_ix_iy))
                
!                 dsp = dvx_dx +dvy_dy +dvz_dz
!                  rp = rf_p(i)
                
!                 corr(j)=corr(j) + dsp*rp
                
!             end do
            
!         end do
!         end do
!         !$omp end do
!         !$omp end parallel
        
!     end subroutine

!     subroutine corr2d_flat_moduli(sf_vx,sf_vz,rf_p,&
!                                   corr,            &
!                                   ifz,ilz,ifx,ilx)
!         real,dimension(*) :: sf_vx,sf_vz,rf_p
!         real,dimension(*) :: corr
        
!         nz=cb%nz
        
!         dvx_dx=0.
!         dvz_dz=0.
        
!         dsp=0.
!          rp=0.
        
!         !$omp parallel default (shared)&
!         !$omp private(iz,ix,i,j,&
!         !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
!         !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
!         !$omp         dvx_dx,dvz_dz,&
!         !$omp         dsp,rp)
!         !$omp do schedule(dynamic)
!         do ix=ifx,ilx
        
!             !dir$ simd
!             do iz=ifz,ilz
                
!                 i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
!                 j=(iz-1)     +(ix-1)     *cb%mz !corr has no boundary layers
                
!                 izm1_ix=i-1  !iz-1,ix
!                 iz_ix  =i    !iz,ix
!                 izp1_ix=i+1  !iz+1,ix
!                 izp2_ix=i+2  !iz+2,ix
                
!                 iz_ixm1=i    -nz  !iz,ix-1
!                 iz_ixp1=i    +nz  !iz,ix+1
!                 iz_ixp2=i  +2*nz  !iz,ix+2
                
!                 dvx_dx = c1x*(sf_vx(iz_ixp1)-sf_vx(iz_ix)) +c2x*(sf_vx(iz_ixp2)-sf_vx(iz_ixm1))
!                 dvz_dz = c1z*(sf_vz(izp1_ix)-sf_vz(iz_ix)) +c2z*(sf_vz(izp2_ix)-sf_vz(izm1_ix))
                
!                 dsp = dvx_dx +dvz_dz
!                  rp = rf_p(i)
                
!                 corr(j)=corr(j) + dsp*rp
                
!             end do
            
!         end do
!         !$omp end do
!         !$omp end parallel

!     end subroutine
    
!     subroutine corr3d_flat_density(sf_p,rf_vx,rf_vy,rf_vz,&
!                                    corr,                  &
!                                    ifz,ilz,ifx,ilx,ify,ily)
!         real,dimension(*) :: sf_p,rf_vx,rf_vy,rf_vz
!         real,dimension(*) :: corr
        
!         nz=cb%nz
!         nx=cb%nx
        
!         dsvx=0.; dsvy=0.; dsvz=0.
!          rvx=0.;  rvy=0.;  rvz=0.
        
!         !$omp parallel default (shared)&
!         !$omp private(iz,ix,iy,i,j,&
!         !$omp         izm2_ix_iy,izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
!         !$omp         iz_ixm2_iy,iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
!         !$omp         iz_ix_iym2,iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
!         !$omp         dsvx,dsvy,dsvz,&
!         !$omp          rvx, rvy, rvz)
!         !$omp do schedule(dynamic) collapse(2)
!         do iy=ify,ily
!         do ix=ifx,ilx
        
!             !dir$ simd
!             do iz=ifz,ilz
            
!                 i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(iy-cb%ify)*cb%nz*cb%nx+1 !field has boundary layers
!                 j=(iz-1)     +(ix-1)     *cb%mz+(iy-1)     *cb%mz*cb%mx+1 !corr has no boundary layers

!                 izm2_ix_iy=i-2  !iz-2,ix,iy
!                 izm1_ix_iy=i-1  !iz-1,ix,iy
!                 iz_ix_iy  =i    !iz,ix,iy
!                 izp1_ix_iy=i+1  !iz+1,ix,iy
!                 izp2_ix_iy=i+2  !iz+2,ix,iy
                
!                 iz_ixm2_iy=i  -2*nz  !iz,ix-2,iy
!                 iz_ixm1_iy=i    -nz  !iz,ix-1,iy
!                 iz_ixp1_iy=i    +nz  !iz,ix+1,iy
!                 iz_ixp2_iy=i  +2*nz  !iz,ix+2,iy
                
!                 iz_ix_iym2=i  -2*nz*nx  !iz,ix,iy-2
!                 iz_ix_iym1=i    -nz*nx  !iz,ix,iy-1
!                 iz_ix_iyp1=i    +nz*nx  !iz,ix,iy+1
!                 iz_ix_iyp2=i  +2*nz*nx  !iz,ix,iy+2               

!                 dsvx = (c1x*(sf_p(iz_ix_iy  )-sf_p(iz_ixm1_iy)) +c2x*(sf_p(iz_ixp1_iy)-sf_p(iz_ixm2_iy))) &
!                       +(c1x*(sf_p(iz_ixp1_iy)-sf_p(iz_ix_iy  )) +c2x*(sf_p(iz_ixp2_iy)-sf_p(iz_ixm1_iy)))
!                 dsvy = (c1y*(sf_p(iz_ix_iy  )-sf_p(iz_ix_iym1)) +c2y*(sf_p(iz_ix_iyp1)-sf_p(iz_ix_iym2))) &
!                       +(c1y*(sf_p(iz_ix_iyp1)-sf_p(iz_ix_iy  )) +c2y*(sf_p(iz_ix_iyp2)-sf_p(iz_ix_iym1)))
!                 dsvz = (c1z*(sf_p(iz_ix_iy  )-sf_p(izm1_ix_iy)) +c2z*(sf_p(izp1_ix_iy)-sf_p(izm2_ix_iy))) &
!                       +(c1z*(sf_p(izp1_ix_iy)-sf_p(iz_ix_iy  )) +c2z*(sf_p(izp2_ix_iy)-sf_p(izm1_ix_iy)))
!                 !complete equation with unnecessary terms e.g. sf_p(iz_ix) for better understanding
!                 !with flag -O, the compiler should automatically detect such possibilities of simplification
                
!                 rvx = rf_vx(iz_ixp1_iy) +rf_vx(iz_ix_iy)
!                 rvy = rf_vy(iz_ix_iyp1) +rf_vy(iz_ix_iy)
!                 rvz = rf_vz(izp1_ix_iy) +rf_vz(iz_ix_iy)
                
!                 corr(j)=corr(j) + 0.25*( dsvx*rvx + dsvy*rvy + dsvz*rvz )
                
!             enddo
            
!         enddo
!         enddo
!         !$omp end do
!         !$omp end parallel
        
!     end subroutine
    
!     subroutine corr2d_flat_density(sf_p,rf_vx,rf_vz,&
!                                    corr,            &
!                                    ifz,ilz,ifx,ilx)
!         real,dimension(*) :: sf_p,rf_vx,rf_vz
!         real,dimension(*) :: corr
        
!         nz=cb%nz
        
!         dsvx=0.; dsvz=0.
!          rvx=0.; rvz=0.
        
!         !$omp parallel default (shared)&
!         !$omp private(iz,ix,i,j,&
!         !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
!         !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
!         !$omp         dsvx,dsvz,&
!         !$omp          rvx, rvz)
!         !$omp do schedule(dynamic)
!         do ix=ifx,ilx
        
!             !dir$ simd
!             do iz=ifz,ilz
            
!                 i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
!                 j=(iz-1)     +(ix-1)     *cb%mz+1 !corr has no boundary layers
                
!                 izm2_ix=i-2  !iz-2,ix
!                 izm1_ix=i-1  !iz-1,ix
!                 iz_ix  =i    !iz,ix
!                 izp1_ix=i+1  !iz+1,ix
!                 izp2_ix=i+2  !iz+2,ix
                
!                 iz_ixm2=i  -2*nz  !iz,ix-2
!                 iz_ixm1=i    -nz  !iz,ix-1
!                 iz_ixp1=i    +nz  !iz,ix+1
!                 iz_ixp2=i  +2*nz  !iz,ix+2
                
!                 dsvx = (c1x*(sf_p(iz_ix  )-sf_p(iz_ixm1)) +c2x*(sf_p(iz_ixp1)-sf_p(iz_ixm2))) &
!                       +(c1x*(sf_p(iz_ixp1)-sf_p(iz_ix  )) +c2x*(sf_p(iz_ixp2)-sf_p(iz_ixm1)))
!                 dsvz = (c1z*(sf_p(iz_ix  )-sf_p(izm1_ix)) +c2z*(sf_p(izp1_ix)-sf_p(izm2_ix))) &
!                       +(c1z*(sf_p(izp1_ix)-sf_p(iz_ix  )) +c2z*(sf_p(izp2_ix)-sf_p(izm1_ix)))
!                 !complete equation with unnecessary terms e.g. sf_p(iz_ix) for better understanding
!                 !with flag -O, the compiler should automatically detect such possibilities of simplification
                
!                 rvx = rf_vx(iz_ix) +rf_vx(iz_ixp1)
!                 rvz = rf_vz(iz_ix) +rf_vz(izp1_ix)
                
!                 corr(j)=corr(j) + 0.25*( dsvx*rvx + dsvz*rvz )
                
!             enddo
            
!         enddo
!         !$omp end do
!         !$omp end parallel
        
!     end subroutine
    

end
