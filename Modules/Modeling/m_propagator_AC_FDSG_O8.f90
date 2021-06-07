module m_propagator_fdsg_o8
use m_sysio
use m_arrayop
use m_shot
use m_field
use m_boundarystore

    !info
    character(*),parameter :: info = &
    'Time-domain Isotropic 2D/3D ACoustic system'//s_return// &
    '1st-order Velocity-Stress formulation'//s_return// &
    'Staggered-Grid Finite-Difference (FDSG) method'//s_return// &
    'Cartesian O(x4,t2) stencil'//s_return// &
    'gradients for: kpa rho'
    integer,parameter :: ngrad=2

    !FD coeff
    real,dimension(4),parameter :: coef = [1225./1024, -245./3072., 49./5120., -5./7168]
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z
    real :: c3x, c3y, c3z
    real :: c4x, c4y, c4z

    !blooming
    integer,dimension(:,:),allocatable :: bloom
    integer,parameter :: initial_half_bloomwidth=5 !half bloom width at ift, should involve hicks points
    
    contains
    
    subroutine init(orig,if_will_restore)
        logical,optional :: if_will_restore

        !propagator - local model
        call alloc(self%buox,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],initialize=.false.)
        call alloc(self%buoy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],initialize=.false.)
        call alloc(self%buoz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],initialize=.false.)
        call alloc(self%kpa, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],initialize=.false.)
        call alloc(self%_kpa,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],initialize=.false.)
        
        self%kpa=cb%rho*cb%vp**2
        self%_kpa=1./self%kpa

        self%buoz(cb%ifz,:,:)=1./cb%rho(cb%ifz,:,:)
        self%buox(:,cb%ifx,:)=1./cb%rho(:,cb%ifx,:)
        self%buoy(:,:,cb%ify)=1./cb%rho(:,:,cb%ify)

        do iz=cb%ifz+1,cb%ilz
            self%buoz(iz,:,:)=0.5/cb%rho(iz,:,:)+0.5/cb%rho(iz-1,:,:)
        enddo
        
        do ix=cb%ifx+1,cb%ilx
            self%buox(:,ix,:)=0.5/cb%rho(:,ix,:)+0.5/cb%rho(:,ix-1,:)
        enddo

        do iy=cb%ify+1,cb%ily
            self%buoy(:,:,iy)=0.5/cb%rho(:,:,iy)+0.5/cb%rho(:,:,iy-1)
        enddo

        call alloc(f%vx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%vy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%vz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%p, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        
        call alloc(f%cpml_dvx_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%cpml_dvy_dy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%cpml_dvz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%cpml_dp_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%cpml_dp_dy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%cpml_dp_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        

        !blooming
        call alloc(bloom,6,nt) !sfield blooming
        
        !directional maximum propagation distance per time step
        distz = cb%velmax * dt / m%dz
        distx = cb%velmax * dt / m%dx
        disty = cb%velmax * dt / m%dy
        
        bloom(1,nt)=minval(orig(:)%iz) -initial_half_bloomwidth !shot%rcv(:)%iz
        bloom(2,nt)=maxval(orig(:)%iz) +initial_half_bloomwidth
        bloom(3,nt)=minval(orig(:)%ix) -initial_half_bloomwidth
        bloom(4,nt)=maxval(orig(:)%ix) +initial_half_bloomwidth
        bloom(5,nt)=minval(orig(:)%iy) -initial_half_bloomwidth
        bloom(6,nt)=maxval(orig(:)%iy) +initial_half_bloomwidth
        do it=nt-1,1,-1
            it_fwd=nt-it+1
            bloom(1,it)=max(nint(bloom(1,nt)-it_fwd*distz),cb%ifz) !bloombox ifz
            bloom(2,it)=min(nint(bloom(2,nt)+it_fwd*distz),cb%ilz) !bloombox ilz
            bloom(3,it)=max(nint(bloom(3,nt)-it_fwd*distx),cb%ifx) !bloombox ifx
            bloom(4,it)=min(nint(bloom(4,nt)+it_fwd*distx),cb%ilx) !bloombox ilx
            bloom(5,it)=max(nint(bloom(5,nt)-it_fwd*disty),cb%ify) !bloombox ify
            bloom(6,it)=min(nint(bloom(6,nt)+it_fwd*disty),cb%ily) !bloombox ily
        enddo
        
        if(.not.m%is_cubic) bloom(5:6,:)=1


        !boundarystore
        if(if_will_restore) then
            call init_boundarystore
        endif

    end subroutine
    
    !========= for FDSG O(dx8,dt2) ===================  

    subroutine print_info

        call hud('Invoked field & propagator modules info : '//s_return//self%info)
        if(mpiworld%is_master) then
            write(*,*) 'Coeff:',coef
        endif
        c1x=coef(1)/m%dx; c1y=coef(1)/m%dy; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz
        c3x=coef(3)/m%dx; c3y=coef(3)/m%dy; c3z=coef(3)/m%dz
        c4x=coef(4)/m%dx; c4y=coef(4)/m%dy; c4z=coef(4)/m%dz
        
    end subroutine

    !========= forward propagation =================
    !WE: du_dt = MDu
    !u=[vx vy vz p]^T, p=(sxx+syy+szz)/3=sxx+syy+szz
    !  [diag3(b)  0  ]    [0   0   0   dx]
    !M=[   0     kpa ], D=|0   0   0   dy|
    !  (b=1/rho)          |0   0   0   dz|
    !                     [dx  dy  dz  0 ]
    !
    !Discretization (staggered grid in space and time):
    !  (grid index)     (real index)
    !  s(iz,ix,iy):=   s[iz,ix,iy]^it+0.5,it+1.5,...
    ! vx(iz,ix,iy):=  vx[iz,ix-0.5,iy]^it,it+1,...
    ! vy(iz,ix,iy):=  vy[iz,ix,iy-0.5]^it,it+1,...
    ! vz(iz,ix,iy):=  vy[iz-0.5,ix,iy]^it,it+1,...

    subroutine propagate(self)
        type(t_field) :: self

        real,dimension(:,:),allocatable :: seismo
        
        integer,parameter :: time_dir=1 !time direction
        
        call alloc(self%seismo,shot%nrcv,nt,initialize=.false.)
        
        tt1=0.; tt2=0.; tt3=0.; tt4=0.; tt5=0.; tt6=0.

        if(if_snapshot) open(16,file='snapshot_field%vz',access='stream')
        
        do it=ift,ilt
            if(mod(it,500)==0 .and. mpiworld%is_master) write(*,*) 'it----',it
            if(mod(it,500)==0) call check_field(field,'field')
            
            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add forces to v^it
            call cpu_time(tic)
            call field%inject_velocities(time_dir,it,wavelet(it))
            call cpu_time(toc)
            tt1=tt1+tic-toc
                        
            !step 2: from v^it to v^it+1 by differences of s^it+0.5
            call cpu_time(tic)
            call field%update_velocities(time_dir,it,sbloom(:,it))
            call cpu_time(toc)
            tt2=tt2+tic-toc
            
            !step 3: add pressure to s^it+0.5
            call cpu_time(tic)
            call field%inject_stresses(time_dir,it,wavelet(it))
            call cpu_time(toc)
            tt3=tt3+tic-toc
            
            !step 4: from s^it+0.5 to s^it+1.5 by differences of v^it+1
            call cpu_time(tic)
            call field%update_stresses(time_dir,it,bloom(:,it))
            call cpu_time(toc)
            tt4=tt4+tic-toc
            
            !step 5: sample v^it+1 or s^it+1.5 at receivers
            call cpu_time(tic)
            call field%extract_field(self%seismo(:,it))
            call cpu_time(toc)
            tt5=tt5+tic-toc
            
            !snapshot
            if(if_snapshot) then
            if(mod(it-1,it_delta_snapshot)==0 .or. it==ilt) then
                call field%write(16)
            endif
            endif
            
            !step 6: save v^it+1 in boundary layers
            if(field%if_will_restore) then
                call cpu_time(tic)
                call boundarystore_transport('save',it,sfield)
                call cpu_time(toc)
                tt6=tt6+t6-t5
            endif
            
        enddo
        
        if(if_snapshot) then
            close(16)
            write(*,'(a,i0.5,a)') 'ximage < snapshot_sfield%vz  n1=',cb%nz,' perc=99'
            write(*,'(a,i0.5,a,i0.5,a)') 'xmovie < snapshot_sfield%vz  n1=',cb%nz,' n2=',cb%nx,' clip=?e-?? loop=2 title=%g'
        endif
        
        if(mpiworld%is_master) then
            write(*,*) 'time add source velocities',tt1/mpiworld%max_threads
            write(*,*) 'time update velocities    ',tt2/mpiworld%max_threads
            write(*,*) 'time add source stresses  ',tt3/mpiworld%max_threads
            write(*,*) 'time update stresses      ',tt4/mpiworld%max_threads
            write(*,*) 'time extract field        ',tt5/mpiworld%max_threads
            write(*,*) 'time save boundary        ',tt6/mpiworld%max_threads
        endif
        
        !synthetic data
        call alloc(dsyn,shot%rcv(1)%nt,shot%nrcv,shot%ncomp)
        do i=1,shot%ncomp
            dsyn()=transpose(seismo())
        enddo
        
        deallocate(seismo)

    end subroutine

    !add RHS to v^it
    subroutine inject_velocities(f,time_dir,it,w)
        class(t_field), intent(in) :: f
        integer :: time_dir,it
        real :: w
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        ify=shot%src%ify; iy=shot%src%iy; ily=shot%src%ily
        
        source_term=time_dir*w
        
        if(if_hicks) then
            select case (shot%src%comp)
                case ('vx')
                f%vx(ifz:ilz,ifx:ilx,ify:ily) = f%vx(ifz:ilz,ifx:ilx,ify:ily) + source_term*shot%src%interp_coeff!source_term *lm%buox(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff
                
                case ('vy')
                f%vy(ifz:ilz,ifx:ilx,ify:ily) = f%vy(ifz:ilz,ifx:ilx,ify:ily) + source_term*shot%src%interp_coeff!source_term *lm%buoy(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff
                
                case ('vz')
                f%vz(ifz:ilz,ifx:ilx,ify:ily) = f%vz(ifz:ilz,ifx:ilx,ify:ily) + source_term*shot%src%interp_coeff!source_term *lm%buoz(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff
                
            end select
            
        else
            select case (shot%src%comp)
                case ('vx') !horizontal x force on vx[iz,ix-0.5,iy]
                f%vx(iz,ix,iy) = f%vx(iz,ix,iy) + source_term*f%buox(iz,ix,iy)
                
                case ('vy') !horizontal y force on vy[iz,ix,iy-0.5]
                f%vy(iz,ix,iy) = f%vy(iz,ix,iy) + source_term*f%buoy(iz,ix,iy)
                
                case ('vz') !vertical force     on vz[iz-0.5,ix,iy]
                f%vz(iz,ix,iy) = f%vz(iz,ix,iy) + source_term*f%buoz(iz,ix,iy)
                
            end select
            
        endif
        
    end subroutine
    
    !v^it -> v^it+1 by FD of s^it+0.5
    subroutine update_velocities(f,time_dir,it,bl)
        class(t_field), intent(in) :: f
        integer :: time_dir, it
        
        ifz=self%bloom(1)+2
        ilz=self%bloom(2)-1
        ifx=self%bloom(3)+2
        ilx=self%bloom(4)-1
        ify=self%bloom(5)+2
        ily=self%bloom(6)-1
        
        if(m%if_freesurface) ifz=max(ifz,1)
        
        if(m%is_cubic) then
            call fd3d_velocities(f%vx,f%vy,f%vz,f%p,                         &
                                 f%cpml_dp_dx,f%cpml_dp_dy,f%cpml_dp_dz,     &
                                 f%buox,f%buoy,f%buoz,                       &
                                 ifz,ilz,ifx,ilx,ify,ily,time_dir*shot%src%dt)
        else
            call fd2d_velocities(f%vx,f%vz,f%p,                      &
                                 f%cpml_dp_dx,f%cpml_dp_dz,          &
                                 f%buox,f%buoz,                      &
                                 ifz,ilz,ifx,ilx,time_dir*shot%src%dt)
        endif
        
        !apply free surface boundary condition if needed
        if(m%if_freesurface) call fd_freesurface_velocities(f%vz)
        
    end subroutine
    
    !add RHS to s^it+0.5
    subroutine inject_stresses(f,time_dir,it,w)
        class(t_field), intent(in) :: f
        integer :: time_dir
        real :: w
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        ify=shot%src%ify; iy=shot%src%iy; ily=shot%src%ily
        
        source_term=time_dir*w
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (1)
                f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) + source_term *shot%src%interp_coeff
            end select
            
        else
            select case (shot%src%icomp)
                case (1) !explosion on s[iz,ix,iy]
                f%p(iz,ix,iy) = f%p(iz,ix,iy) + source_term
            end select
            
        endif
        
    end subroutine
    
    !s^it+0.5 -> s^it+1.5 by FD of v^it+1
    subroutine update_stresses(f,time_dir,it,bl)
        class(t_field), intent(in) :: f
        integer :: time_dir, it
        integer,dimension(6) :: bl
        
        ifz=bl(1)+1
        ilz=bl(2)-2
        ifx=bl(3)+1
        ilx=bl(4)-2
        ify=bl(5)+1
        ily=bl(6)-2
        
        if(m%if_freesurface) ifz=max(ifz,1)
        
        if(m%is_cubic) then
            call fd3d_stresses(f%vx,f%vy,f%vz,f%p,                         &
                                    f%cpml_dvx_dx,f%cpml_dvy_dy,f%cpml_dvz_dz,  &
                                    lm%kpa,                                     &
                                    ifz,ilz,ifx,ilx,ify,ily,time_dir*shot%src%dt)
        else
            call fd2d_stresses(f%vx,f%vz,f%p,                      &
                                    f%cpml_dvx_dx,f%cpml_dvz_dz,        &
                                    lm%kpa,                             &
                                    ifz,ilz,ifx,ilx,time_dir*shot%src%dt)
        endif
        
        !apply free surface boundary condition if needed
        if(m%if_freesurface) call fd_freesurface_stresses(f%p)

    end subroutine
    
    !========= adjoint propagation =================
    !WE:        du_dt   = MDu   + src
    !->        Ndu_dt   = Du    + Nsrc, N=M^-1
    !Adjoint:  Ndv_dt^T = D^Tv  + adjsrc
    !->        -dv_dt   =-MDv   + Mr
    !     (reverse time) (negative derivatives)
    !
    !Discrete form (staggered grid in space and time, 2D as example):
    ! N  [vx^it+1  ] [vx^it    ]   [ 0     0    dx(-)][vx^it+1  ]
    !----|vz^it+1  |-|vz^it    | = | 0     0    dz(-)||vz^it+1  |  +Nsrc
    ! dt [ s^it+1.5] [ s^it+0.5]   [dx(+) dz(+)  0   ][ s^it+0.5]
    !dx(-):=a1(s(ix)-s(ixm1))+a2(s(ixp1)-s(ixm2))  (O(x4))
    !dx(+):=a1(v(ixp1)-v(ix))+a2(v(ixp2)-v(ixm1))  (O(x4))
    !Adjoint:
    ! N  [vx^it    ] [vx^it+1  ]   [ 0       0      dx(+)^T][vx^it+1  ]
    !----|vz^it    |-|vz^it+1  | = | 0       0      dz(+)^T||vz^it+1  |  +Nsrc
    ! dt [ s^it+0.5] [ s^it+1.5]   [dx(-)^T dz(-)^T  0     ][ s^it+0.5]
    !dx(+)^T=a1(s(ixm1)-s(ix))+a2(s(ixm2)-s(ixp1))=-dx(-)
    !dx(-)^T=a1(v(ix)-v(ixp1))+a2(v(ixm1)-v(ixp2))=-dx(+)
    !-dx,-dz can be regarded as -dt, thus allowing to use same code with flip of dt sign.
    !transposes of time derivatives centered on v^it+0.5 and s^it+1
    !so v^it+1   - v^it     => v^it - v^it+1
    !   s^it+1.5 - s^it+0.5 => s^it+0.5 - s^it+1.5
    !
    !In each time step:
    !d=RGAsrc, src:source term, A:put source, G:propagator, R:sampling at receivers
    !d=[vx vy vz p]^T, R=[I  0   0 ], A=[diag3(b) 0], src=[fx fy fz p]^T
    !                    |0 2/3 1/3|    [   0     1]
    !                    [   0    1]
    !Adjoint: src=A^T N G^T M R^T d
    !M R^T:put adjoint sources, 
    !G^T:  adjoint propagator,
    !A^T N:get adjoint fields
        
    subroutine propagate_reverse(o_sfield,o_imag,o_grad,o_rdt,dout)
        real,dimension(nt),optional :: dout
        real,dimension(cb%mz,cb%mx,cb%my,3),optional :: gradient
        
        real,dimension(:,:),allocatable :: seismo
        
        integer,parameter :: time_dir=-1 !time direction
        
        !reinitialize memory for incident wavefield reconstruction
        sfield%cpml_dvx_dx=0.
        sfield%cpml_dvy_dy=0.
        sfieldf%cpml_dvz_dz=0.
        sfieldf%cpml_dp_dx=0.
        sfieldf%cpml_dp_dy=0.
        sfieldf%cpml_dp_dz=0.

        call alloc(seismo,shot%nrcv,nt,initialize=.false.)
        seismo=transpose(dres) !to save mem space, seismo could be just time slices..
        
        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.

        if(if_snapshot) then
            open(20,file='snapshot_sfield%vz_back',access='stream')
            !open(22,file='snapshot_sfield%vz_deri',access='stream')
            open(24,file='snapshot_rfield%vz',access='stream')
            open(26,file='snapshot_scorr',access='stream')
            open(28,file='snapshot_vcorr',access='stream')
        endif
        
        
        do it=ilt,ift,time_dir
            
            if(mod(it,500)==0 .and. mpiworld%is_master) write(*,*) 'it----',it
            if(mod(it,500)==0) call check_field(sfield,'sfield')
            if(mod(it,500)==0) call check_field(rfield,'rfield')
            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            if(present(o_sfield)) then
                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call boundarystore_transport('load',it,sfield)
                call cpu_time(toc)
                tt1=tt1+tic-toc
                
                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call sfield%update_stresses(time_dir,it,sbloom(:,it))
                call cpu_time(toc)
                tt2=tt2+tic-toc

                !backward step 3: rm pressure from s^it+0.5
                call cpu_time(tic)
                call sfield%inject_stresses(time_dir,it,wavelet(it))
                call cpu_time(toc)
                tt3=tt3+tic-toc
            endif

            !--------------------------------------------------------!

            !adjoint step 5: fill s^it+1.5 at receivers
            call cpu_time(tic)
            call rfield%inject_stresses_adjoint(time_dir,it,seismo(:,it))
            call cpu_time(toc)
            tt4=tt4+tic-toc

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call rfield%update_stresses_adjoint(time_dir,it,rbloom(:,it))
            call cpu_time(toc)
            tt5=tt5+tic-toc

            !gkpa: sfield%s_dt^it+0.5 \dot rfield%s^it+0.5
            !use sfield%v^it+1 to compute sfield%s_dt^it+0.5, as backward step 4
            if(present(o_grad)) then
                call cpu_time(tic)
                call rfield%grad_moduli(sfield,it,sbloom(:,it),rbloom(:,it),gradient)
                call cpu_time(toc)
            endif
            tt6=tt6+tic-toc

            if(present(o_imag)) then
                call cpu_time(tic)
                call rfield%imag(sfield,it,sbloom(:,it),rbloom(:,it),image)
                call cpu_time(toc)
            endif
            tt6=tt6+tic-toc

            !========================================================!

            if(present(o_sfield)) then
                !backward step 2: v^it+1 -> v^it by FD of s^it+0.5
                call cpu_time(tic)
                call sfield%update_velocities(time_dir,it,sbloom(:,it))
                call cpu_time(toc)
                tt7=tt7+tic-toc

                !backward step 1: rm forces from v^it
                call cpu_time(tic)
                call sfield%inject_velocities(time_dir,it,wavelet(it))
                call cpu_time(toc)
                tt8=tt8+tic-toc
            endif

            !--------------------------------------------------------!

            !adjoint step 5: fill v^it+1 at receivers
            call cpu_time(tic)
            call rfield%inject_velocities_adjoint(time_dir,it,seismo(:,it))
            call cpu_time(toc)
            tt9=tt9+tic-toc

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call rfield%update_velocities_adjoint(time_dir,it,rbloom(:,it))
            call cpu_time(toc)
            tt10=tt10+tic-toc
            
            !adjoint step 1: sample v^it or s^it+0.5 at source position
            if(present(dout)) then
                call cpu_time(tic)
                call rfield%extract_adjoint(dout(it))
                call cpu_time(toc)
                tt11=tt11+tic-toc
            endif
            
            !grho: sfield%v_dt^it \dot rfield%v^it
            !use sfield%s^it+0.5 to compute sfield%v_dt^it, as backward step 2
            if(present(o_grad)) then
                call cpu_time(tic)
                call rfield%grad_density(sfield,it,sbloom(:,it),rbloom(:,it),gradient)
                call cpu_time(toc)
                tt12=tt12+t12-t11
            endif          
            
            
            !snapshot
            if(if_snapshot) then
            if(mod(it-1,it_delta_snapshot)==0 .or. it==ift) then
                call write_field(20,sfield)
                call write_field(24,rfield)
                write(26)gradient(:,:,:,1)
                write(28)gradient(:,:,:,3)
            endif
            endif
                        
        enddo
        
        !scale gradient
        if(present(o_grad)) then
            call rfield%gradient_scaling(o_grad)
        endif
        
        if(if_snapshot) then
            close(20)
            !close(22)
            close(24)
            close(26)
            close(28)
            write(*,'(a,i0.5,a)') 'ximage < snapshot_rfield%vz  n1=',cb%nz,' perc=99'
            write(*,'(a,i0.5,a,i0.5,a)') 'xmovie < snapshot_rfield%vz  n1=',cb%nz,' n2=',cb%nx,' clip=?e-?? loop=2 title=%g'
            
            write(*,'(a,i0.5,a,i0.5,a)') 'xmovie < snapshot_grad  n1=',cb%mz,' n2=',cb%mx,' clip=?e-?? loop=2 title=%g'
            write(*,'(a,i0.5,a,i0.5,a)') 'xmovie < snapshot_grad  n1=',cb%mz,' n2=',cb%mx,' clip=?e-?? loop=2 title=%g'
        endif
        
        
        if(mpiworld%is_master) then
            write(*,*) 'time load boundary            ',tt1/mpiworld%max_threads
            write(*,*) 'time update stresses          ',tt2/mpiworld%max_threads
            write(*,*) 'time rm source stresses       ',tt3/mpiworld%max_threads
            write(*,*) 'time update velocities        ',tt7/mpiworld%max_threads
            write(*,*) 'time rm source velocities     ',tt8/mpiworld%max_threads
            write(*,*) 'time -------------------------'
            write(*,*) 'time add adjsource stresses   ',tt4/mpiworld%max_threads
            write(*,*) 'time update adj stresses      ',tt5/mpiworld%max_threads
            write(*,*) 'time add adjsource velocities ',tt9/mpiworld%max_threads
            write(*,*) 'time update adj velocities    ',tt10/mpiworld%max_threads
            write(*,*) 'time extract&write fields     ',tt11/mpiworld%max_threads
            write(*,*) 'time correlation              ',(tt6+tt12)/mpiworld%max_threads

        endif
        
        deallocate(seismo)

    end subroutine

    !add RHS to s^it+1.5
    subroutine inject_stresses_adjoint(f,time_dir,it,adjsource)
        class(t_field), intent(in) :: f
        integer :: time_dir
        real,dimension(*) :: adjsource
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            ify=shot%rcv(ircv)%ify; iy=shot%rcv(ircv)%iy; ily=shot%rcv(ircv)%ily
            
            tmp=adjsource(ircv)
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (1) !pressure adjsource
                    f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) &
                        +tmp* lm%kpa(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (1) !pressure adjsource
                    !p[iz,ix,iy]
                    f%p(iz,ix,iy) = f%p(iz,ix,iy)  &
                        +tmp* lm%kpa(iz,ix,iy)
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    !s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
    subroutine update_stresses_adjoint(f,time_dir,it,bl)
        class(t_field), intent(in) :: f
        integer :: time_dir, it
        integer,dimension(6) :: bl
        
        !same code with -dt
        call f%update_stresses(time_dir,it,bl)
        
    end subroutine
    
    !add RHS to v^it+1
    subroutine inject_velocities_adjoint(f,time_dir,it,adjsource)
        class(t_field), intent(in) :: f
        integer :: time_dir
        real,dimension(*) :: adjsource
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            ify=shot%rcv(ircv)%ify; iy=shot%rcv(ircv)%iy; ily=shot%rcv(ircv)%ily
            
            tmp=adjsource(ircv)
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (2) !horizontal x adjsource
                    f%vx(ifz:ilz,ifx:ilx,ify:ily) = f%vx(ifz:ilz,ifx:ilx,ify:ily) + tmp*lm%buox(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
                    
                    case (3) !horizontal y adjsource
                    f%vy(ifz:ilz,ifx:ilx,ify:ily) = f%vy(ifz:ilz,ifx:ilx,ify:ily) + tmp*lm%buoy(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
                    
                    case (4) !horizontal z adjsource
                    f%vz(ifz:ilz,ifx:ilx,ify:ily) = f%vz(ifz:ilz,ifx:ilx,ify:ily) + tmp*lm%buoz(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (2) !horizontal x adjsource
                    !vx[ix-0.5,iy,iz]
                    f%vx(iz,ix,iy) = f%vx(iz,ix,iy) + tmp*lm%buox(iz,ix,iy)
                    
                    case (3) !horizontal y adjsource
                    !vy[ix,iy-0.5,iz]
                    f%vy(iz,ix,iy) = f%vy(iz,ix,iy) + tmp*lm%buoy(iz,ix,iy)
                    
                    case (4) !horizontal z adjsource
                    !vz[ix,iy,iz-0.5]
                    f%vz(iz,ix,iy) = f%vz(iz,ix,iy) + tmp*lm%buoz(iz,ix,iy)
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    !v^it+1 -> v^it by FD^T of s^it+0.5
    subroutine update_velocities_adjoint(f,time_dir,it,bl)
        class(t_field), intent(in) :: f
        integer :: time_dir, it
        integer,dimension(6) :: bl
        
        !same code with -dt
        call f%update_velocities(time_dir,it,bl)
        
    end subroutine
    
    !========= for wavefield correlation ===================   
    !gkpa = -1/kpa2    dp_dt \dot adjp
    !     = -1/kpa \nabla.v  \dot adjp
    !v^it+1, p^it+0.5, adjp^it+0.5
    !
    !grho = dv_dt      \dot adjv
    !     = b \nabla p \dot adjv
    !p^it+0.5, v^it, adjv^it
    !use (v[i+1]+v[i])/2 to approximate v[i+0.5], so is adjv
    
    subroutine gkpa(sf,rf,it,sb,rb)
        class(t_field), intent(in) :: sf, rf
        integer,dimension(6) :: sb,rb
        real,dimension(cb%mz,cb%mx,cb%my,ngrad) :: grad
        
        !nonzero only when sf touches rf
        ifz=max(sb(1),rb(1),2)
        ilz=min(sb(2),rb(2),cb%mz-2)
        ifx=max(sb(3),rb(3),2)
        ilx=min(sb(4),rb(4),cb%mx-2)
        ify=max(sb(5),rb(5),2)
        ily=min(sb(6),rb(6),cb%my-2)
        
        if(m%is_cubic) then
            call grad3d_moduli(sf%vx,sf%vy,sf%vz,rf%p,&
                               grad(:,:,:,1),         &
                               ifz,ilz,ifx,ilx,ify,ily)
        else
            call grad2d_moduli(sf%vx,sf%vz,rf%p,&
                               grad(:,:,:,1),   &
                               ifz,ilz,ifx,ilx)
        endif
        
    end subroutine
    
    subroutine grho(sf,rf,it,sb,rb)
        class(t_field), intent(in) :: sf,rf
        integer,dimension(6) :: sb,rb
        real,dimension(cb%mz,cb%mx,cb%my,ngrad) :: grad
        
        !nonzero only when sf touches rf
        ifz=max(sb(1),rb(1),2)
        ilz=min(sb(2),rb(2),cb%mz-2)
        ifx=max(sb(3),rb(3),2)
        ilx=min(sb(4),rb(4),cb%mx-2)
        ify=max(sb(5),rb(5),2)
        ily=min(sb(6),rb(6),cb%my-2)
        
        if(m%is_cubic) then
            call grad3d_density(sf%p,rf%vx,rf%vy,rf%vz,&
                                grad(:,:,:,2),         &
                                ifz,ilz,ifx,ilx,ify,ily)
        else
            call grad2d_density(sf%p,rf%vx,rf%vz,&
                                grad(:,:,:,2),   &
                                ifz,ilz,ifx,ilx)
        endif
        
    end subroutine
    
    subroutine grad_scaling(grad)
        real,dimension(cb%mz,cb%mx,cb%my,ngrad) :: grad
        
        grad(:,:,:,1)=grad(:,:,:,1) * (-lm%invkpa(1:cb%mz,1:cb%mx,1:cb%my))

        grad(:,:,:,2)=grad(:,:,:,2) / cb%rho(1:cb%mz,1:cb%mx,1:cb%my)

        grad(1,:,:,:) = grad(2,:,:,:)
        grad(cb%mz-1,:,:,:) = grad(cb%mz-2,:,:,:)
        grad(cb%mz,  :,:,:) = grad(cb%mz-2,:,:,:)

        grad(:,1,:,:) = grad(:,2,:,:)
        grad(:,cb%mx-1,:,:) = grad(:,cb%mx-2,:,:)
        grad(:,cb%mx  ,:,:) = grad(:,cb%mx-2,:,:)

        if(m%is_cubic) then
            grad(:,:,1,:) = grad(:,2,:,:)
            grad(:,:,cb%my-1,:) = grad(:,:,cb%my-2,:)
            grad(:,:,cb%my  ,:) = grad(:,:,cb%my-2,:)
        endif
        
        !set unit of gkpa to be [m3], grho to be [m5/s2]
        !such that after multiplied by (kpa_max-kpa_min) or (rho_max-rho_min) (will be done in m_parameterization.f90)
        !the unit of parameter update is [Nm], same as Lagrangian
        !and the unit of gradient scaling factor is [1/N/m] (in m_scaling.f90)
        !therefore parameters become unitless
        grad=grad*m%cell_volume*shot%src%dt
        
    end subroutine

    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd3d_flat_velocities(vx,vy,vz,p,                       &
                                    cpml_dp_dx,cpml_dp_dy,cpml_dp_dz, &
                                    buox,buoy,buoz,                   &
                                    ifz,ilz,ifx,ilx,ify,ily,dt)
        real,dimension(*) :: vx,vy,vz,p
        real,dimension(*) :: cpml_dp_dx,cpml_dp_dy,cpml_dp_dz
        real,dimension(*) :: buox,buoy,buoz
        
        nz=cb%nz
        nx=cb%nx
        ny=cb%ny
        
        dp_dx=0.;dp_dy=0.;dp_dz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,&
        !$omp         izm2_ix_iy,izm1_ix_iy,iz_ix_iy,izp1_ix_iy,&
        !$omp         iz_ixm2_iy,iz_ixm1_iy,iz_ixp1_iy,&
        !$omp         iz_ix_iym2,iz_ix_iym1,iz_ix_iyp1,&
        !$omp         dp_dx,dp_dy,dp_dz)
        !$omp do schedule(dynamic) collapse(2)
        do iy=ify,ily
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+(iy-cb%ify)*nz*nx+1
                
                izm2_ix_iy=i-2  !iz-2,ix,iy
                izm1_ix_iy=i-1  !iz-1,ix,iy
                iz_ix_iy  =i    !iz,ix,iy
                izp1_ix_iy=i+1  !iz+1,ix,iy
                
                iz_ixm2_iy=i  -2*nz  !iz,ix-2,iy
                iz_ixm1_iy=i    -nz  !iz,ix-1,iy
                iz_ixp1_iy=i    +nz  !iz,ix+1,iy
                
                iz_ix_iym2=i  -2*nz*nx  !iz,ix,iy-2
                iz_ix_iym1=i    -nz*nx  !iz,ix,iy-1
                iz_ix_iyp1=i    +nz*nx  !iz,ix,iy+1
                
                dp_dx= c1x*(p(iz_ix_iy)-p(iz_ixm1_iy)) +c2x*(p(iz_ixp1_iy)-p(iz_ixm2_iy))
                dp_dy= c1y*(p(iz_ix_iy)-p(iz_ix_iym1)) +c2y*(p(iz_ix_iyp1)-p(iz_ix_iym2))
                dp_dz= c1z*(p(iz_ix_iy)-p(izm1_ix_iy)) +c2z*(p(izp1_ix_iy)-p(izm2_ix_iy))
                
                !cpml
                cpml_dp_dx(iz_ix_iy)= cb%b_x_half(ix)*cpml_dp_dx(iz_ix_iy) + cb%a_x_half(ix)*dp_dx
                cpml_dp_dy(iz_ix_iy)= cb%b_y_half(iy)*cpml_dp_dy(iz_ix_iy) + cb%a_y_half(iy)*dp_dy
                cpml_dp_dz(iz_ix_iy)= cb%b_z_half(iz)*cpml_dp_dz(iz_ix_iy) + cb%a_z_half(iz)*dp_dz

                dp_dx=dp_dx*cb%kappa_x_half(ix) + cpml_dp_dx(iz_ix_iy)
                dp_dy=dp_dy*cb%kappa_y_half(iy) + cpml_dp_dy(iz_ix_iy)
                dp_dz=dp_dz*cb%kappa_z_half(iz) + cpml_dp_dz(iz_ix_iy)
                
                !velocity
                vx(iz_ix_iy)=vx(iz_ix_iy) + dt*buox(iz_ix_iy)*dp_dx
                vy(iz_ix_iy)=vy(iz_ix_iy) + dt*buoy(iz_ix_iy)*dp_dy
                vz(iz_ix_iy)=vz(iz_ix_iy) + dt*buoz(iz_ix_iy)*dp_dz
                
            enddo
            
        enddo
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd2d_flat_velocities(vx,vz,p,              &
                                    cpml_dp_dx,cpml_dp_dz,&
                                    buox,buoz,            &
                                    ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vx,vz,p
        real,dimension(*) :: cpml_dp_dx,cpml_dp_dz
        real,dimension(*) :: buox,buoz
        
        nz=cb%nz
        nx=cb%nx
        
        dp_dx=0.; dp_dz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dp_dx,dp_dz)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1
                
                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                
                iz_ixm2=i  -2*nz  !iz,ix-2
                iz_ixm1=i    -nz  !iz,ix-1
                iz_ixp1=i    +nz  !iz,ix+1
                
                dp_dx= c1x*(p(iz_ix)-p(iz_ixm1)) +c2x*(p(iz_ixp1)-p(iz_ixm2))
                dp_dz= c1z*(p(iz_ix)-p(izm1_ix)) +c2z*(p(izp1_ix)-p(izm2_ix))
                
                !cpml
                cpml_dp_dx(iz_ix)= cb%b_x_half(ix)*cpml_dp_dx(iz_ix) + cb%a_x_half(ix)*dp_dx
                cpml_dp_dz(iz_ix)= cb%b_z_half(iz)*cpml_dp_dz(iz_ix) + cb%a_z_half(iz)*dp_dz

                dp_dx=dp_dx*cb%kappa_x_half(ix) + cpml_dp_dx(iz_ix)
                dp_dz=dp_dz*cb%kappa_z_half(iz) + cpml_dp_dz(iz_ix)
                
                !velocity
                vx(iz_ix)=vx(iz_ix) + dt*buox(iz_ix)*dp_dx
                vz(iz_ix)=vz(iz_ix) + dt*buoz(iz_ix)*dp_dz
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd3d_flat_stresses(vx,vy,vz,p,                         &
                                  cpml_dvx_dx,cpml_dvy_dy,cpml_dvz_dz,&
                                  kpa,                                &
                                  ifz,ilz,ifx,ilx,ify,ily,dt)
        real,dimension(*) :: vx,vy,vz,p
        real,dimension(*) :: cpml_dvx_dx,cpml_dvy_dy,cpml_dvz_dz
        real,dimension(*) :: kpa
        
        nz=cb%nz
        nx=cb%nx
        ny=cb%ny
        
        dvx_dx=0.;dvy_dy=0.;dvz_dz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,&
        !$omp         izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
        !$omp         iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
        !$omp         iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
        !$omp         dvx_dx,dvy_dy,dvz_dz)
        !$omp do schedule(dynamic) collapse(2)
        do iy=ify,ily
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+(iy-cb%ify)*nz*nx+1
                
                izm1_ix_iy=i-1  !iz-1,ix,iy
                iz_ix_iy  =i    !iz,ix,iy
                izp1_ix_iy=i+1  !iz+1,ix,iy
                izp2_ix_iy=i+2  !iz+2,ix,iy
                
                iz_ixm1_iy=i       -nz  !iz,ix-1,iy
                iz_ixp1_iy=i       +nz  !iz,ix+1,iy
                iz_ixp2_iy=i     +2*nz  !iz,ix+2,iy
                
                iz_ix_iym1=i    -nz*nx  !iz,ix,iy-1
                iz_ix_iyp1=i    +nz*nx  !iz,ix,iy+1
                iz_ix_iyp2=i  +2*nz*nx  !iz,ix,iy+2
                
                dvx_dx= c1x*(vx(iz_ixp1_iy)-vx(iz_ix_iy))  +c2x*(vx(iz_ixp2_iy)-vx(iz_ixm1_iy))
                dvy_dy= c1y*(vy(iz_ix_iyp1)-vy(iz_ix_iy))  +c2y*(vy(iz_ix_iyp2)-vy(iz_ix_iym1))
                dvz_dz= c1z*(vz(izp1_ix_iy)-vz(iz_ix_iy))  +c2z*(vz(izp2_ix_iy)-vz(izm1_ix_iy))
                
                !cpml
                cpml_dvx_dx(iz_ix_iy)=cb%b_x(ix)*cpml_dvx_dx(iz_ix_iy)+cb%a_x(ix)*dvx_dx
                cpml_dvy_dy(iz_ix_iy)=cb%b_y(iy)*cpml_dvy_dy(iz_ix_iy)+cb%a_y(iy)*dvy_dy
                cpml_dvz_dz(iz_ix_iy)=cb%b_z(iz)*cpml_dvz_dz(iz_ix_iy)+cb%a_z(iz)*dvz_dz

                dvx_dx=dvx_dx*cb%kappa_x(ix) + cpml_dvx_dx(iz_ix_iy)
                dvy_dy=dvy_dy*cb%kappa_y(iy) + cpml_dvy_dy(iz_ix_iy)
                dvz_dz=dvz_dz*cb%kappa_z(iz) + cpml_dvz_dz(iz_ix_iy)
                
                !pressure
                p(iz_ix_iy) = p(iz_ix_iy) + dt * kpa(iz_ix_iy)*(dvx_dx+dvy_dy+dvz_dz)
                
            enddo
            
        enddo
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine
    
    subroutine fd2d_flat_stresses(vx,vz,p,                &
                                  cpml_dvx_dx,cpml_dvz_dz,&
                                  kpa,                    &
                                  ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vx,vz,p
        real,dimension(*) :: cpml_dvx_dx,cpml_dvz_dz
        real,dimension(*) :: kpa
        
        nz=cb%nz
        nx=cb%nx
        
        dvx_dx=0.;dvz_dz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dvx_dx,dvz_dz)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1
                
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1
                iz_ixp2=i  +2*nz !iz,ix+2
                
                dvx_dx= c1x*(vx(iz_ixp1)-vx(iz_ix))  +c2x*(vx(iz_ixp2)-vx(iz_ixm1))
                dvz_dz= c1z*(vz(izp1_ix)-vz(iz_ix))  +c2z*(vz(izp2_ix)-vz(izm1_ix))
                
                !cpml
                cpml_dvx_dx(iz_ix)=cb%b_x(ix)*cpml_dvx_dx(iz_ix)+cb%a_x(ix)*dvx_dx
                cpml_dvz_dz(iz_ix)=cb%b_z(iz)*cpml_dvz_dz(iz_ix)+cb%a_z(iz)*dvz_dz

                dvx_dx=dvx_dx*cb%kappa_x(ix) + cpml_dvx_dx(iz_ix)
                dvz_dz=dvz_dz*cb%kappa_z(iz) + cpml_dvz_dz(iz_ix)
                
                !pressure
                p(iz_ix) = p(iz_ix) + dt * kpa(iz_ix)*(dvx_dx+dvz_dz)
                
            enddo
            
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine


    subroutine fd_freesurface_velocities(vz)
        !free surface is located at [1,ix,iy] level
        !so symmetric mirroring: vz[0.5]=vz[1.5], ie. vz(1,ix,iy)=vz(2,ix,iy) -> dp(1,ix,iy)=0.
            vz(1,:,:)=vz(2,:,:)
!             !$omp parallel default (shared)&
!             !$omp private(ix,iy,i)
!             !$omp do schedule(dynamic)
!             do iy=ify,ily
!             do ix=ifx,ilx
!                 i=(1-cb%ifz) + (ix-cb%ifx)*nz + (iy-cb%ify)*nz*nx +1 !iz=1,ix,iy
!                 
!                 f%vz(i)=f%vz(i+1)
!             enddo
!             enddo
!             !$omp enddo
!             !$omp end parallel
    endif

    subroutine fd_freesurface_stresses(p)
        !free surface is located at [1,ix,iy] level
        !so explicit boundary condition: p(1,ix,iy)=0
        !and antisymmetric mirroring: p(0,ix,iy)=-p(2,ix,iy) -> vz(2,ix,iy)=vz(1,ix,iy)
            f%p(1,:,:)=0.
            f%p(0,:,:)=-f%p(2,:,:)
!             !$omp parallel default (shared)&
!             !$omp private(ix,iy,i)
!             !$omp do schedule(dynamic)
!             do iy=ify,ily
!             do ix=ifx,ilx
!                 i=(1-cb%ifz) + (ix-cb%ifx)*nz + (iy-cb%ify)*nz*nx +1 !iz=1,ix,iy 
!                 
!                 f%p(i)=0.
!                 
!                 f%p(i-1)=-f%p(i+1)
!             enddo
!             enddo
!             !$omp enddo
!             !$omp end parallel
    endif

    
    subroutine grad3d_flat_moduli(sf_vx,sf_vy,sf_vz,rf_p,&
                                  grad,                  &
                                  ifz,ilz,ifx,ilx,ify,ily)
        real,dimension(*) :: sf_vx,sf_vy,sf_vz,rf_p
        real,dimension(*) :: grad
        
        nz=cb%nz
        nx=cb%nx
        
        dvx_dx=0.
        dvy_dx=0.
        dvz_dz=0.
        
        dsp=0.
         rp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,j,&
        !$omp         izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
        !$omp         iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
        !$omp         iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
        !$omp         dvx_dx,dvy_dy,dvz_dz,&
        !$omp         dsp,rp)
        !$omp do schedule(dynamic) collapse(2)
        do iy=ify,ily
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(iy-cb%ify)*cb%nz*cb%nx+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+(iy-1)     *cb%mz*cb%mx+1 !grad has no boundary layers
                
                izm1_ix_iy=i-1  !iz-1,ix,iy
                iz_ix_iy  =i    !iz,ix,iy
                izp1_ix_iy=i+1  !iz+1,ix,iy
                izp2_ix_iy=i+2  !iz+2,ix,iy
                
                iz_ixm1_iy=i    -nz  !iz,ix-1,iy
                iz_ixp1_iy=i    +nz  !iz,ix+1,iy
                iz_ixp2_iy=i  +2*nz  !iz,ix+2,iy
                
                iz_ix_iym1=i    -nz*nx  !iz,ix,iy-1
                iz_ix_iyp1=i    +nz*nx  !iz,ix,iy+1
                iz_ix_iyp2=i  +2*nz*nx  !iz,ix,iy+2
                
                dvx_dx = c1x*(sf_vx(iz_ixp1_iy)-sf_vx(iz_ix_iy)) +c2x*(sf_vx(iz_ixp2_iy)-sf_vx(iz_ixm1_iy))
                dvy_dy = c1y*(sf_vy(iz_ix_iyp1)-sf_vy(iz_ix_iy)) +c2y*(sf_vy(iz_ix_iyp2)-sf_vy(iz_ix_iym1))
                dvz_dz = c1z*(sf_vz(izp1_ix_iy)-sf_vz(iz_ix_iy)) +c2z*(sf_vz(izp2_ix_iy)-sf_vz(izm1_ix_iy))
                
                dsp = dvx_dx +dvy_dy +dvz_dz
                 rp = rf_p(i)
                
                grad(j)=grad(j) + dsp*rp
                
            end do
            
        end do
        end do
        !$omp end do
        !$omp end parallel
        
    end subroutine

    subroutine grad2d_flat_moduli(sf_vx,sf_vz,rf_p,&
                                  grad,            &
                                  ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_vx,sf_vz,rf_p
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        dvx_dx=0.
        dvz_dz=0.
        
        dsp=0.
         rp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dvx_dx,dvz_dz,&
        !$omp         dsp,rp)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm1=i    -nz  !iz,ix-1
                iz_ixp1=i    +nz  !iz,ix+1
                iz_ixp2=i  +2*nz  !iz,ix+2
                
                dvx_dx = c1x*(sf_vx(iz_ixp1)-sf_vx(iz_ix)) +c2x*(sf_vx(iz_ixp2)-sf_vx(iz_ixm1))
                dvz_dz = c1z*(sf_vz(izp1_ix)-sf_vz(iz_ix)) +c2z*(sf_vz(izp2_ix)-sf_vz(izm1_ix))
                
                dsp = dvx_dx +dvz_dz
                 rp = rf_p(i)
                
                grad(j)=grad(j) + dsp*rp
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine
    
    subroutine grad3d_flat_density(sf_p,rf_vx,rf_vy,rf_vz,&
                                   grad,                  &
                                   ifz,ilz,ifx,ilx,ify,ily)
        real,dimension(*) :: sf_p,rf_vx,rf_vy,rf_vz
        real,dimension(*) :: grad
        
        nz=cb%nz
        nx=cb%nx
        
        dsvx=0.; dsvy=0.; dsvz=0.
         rvx=0.;  rvy=0.;  rvz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,j,&
        !$omp         izm2_ix_iy,izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
        !$omp         iz_ixm2_iy,iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
        !$omp         iz_ix_iym2,iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
        !$omp         dsvx,dsvy,dsvz,&
        !$omp          rvx, rvy, rvz)
        !$omp do schedule(dynamic) collapse(2)
        do iy=ify,ily
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(iy-cb%ify)*cb%nz*cb%nx+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+(iy-1)     *cb%mz*cb%mx+1 !grad has no boundary layers

                izm2_ix_iy=i-2  !iz-2,ix,iy
                izm1_ix_iy=i-1  !iz-1,ix,iy
                iz_ix_iy  =i    !iz,ix,iy
                izp1_ix_iy=i+1  !iz+1,ix,iy
                izp2_ix_iy=i+2  !iz+2,ix,iy
                
                iz_ixm2_iy=i  -2*nz  !iz,ix-2,iy
                iz_ixm1_iy=i    -nz  !iz,ix-1,iy
                iz_ixp1_iy=i    +nz  !iz,ix+1,iy
                iz_ixp2_iy=i  +2*nz  !iz,ix+2,iy
                
                iz_ix_iym2=i  -2*nz*nx  !iz,ix,iy-2
                iz_ix_iym1=i    -nz*nx  !iz,ix,iy-1
                iz_ix_iyp1=i    +nz*nx  !iz,ix,iy+1
                iz_ix_iyp2=i  +2*nz*nx  !iz,ix,iy+2               

                dsvx = (c1x*(sf_p(iz_ix_iy  )-sf_p(iz_ixm1_iy)) +c2x*(sf_p(iz_ixp1_iy)-sf_p(iz_ixm2_iy))) &
                      +(c1x*(sf_p(iz_ixp1_iy)-sf_p(iz_ix_iy  )) +c2x*(sf_p(iz_ixp2_iy)-sf_p(iz_ixm1_iy)))
                dsvy = (c1y*(sf_p(iz_ix_iy  )-sf_p(iz_ix_iym1)) +c2y*(sf_p(iz_ix_iyp1)-sf_p(iz_ix_iym2))) &
                      +(c1y*(sf_p(iz_ix_iyp1)-sf_p(iz_ix_iy  )) +c2y*(sf_p(iz_ix_iyp2)-sf_p(iz_ix_iym1)))
                dsvz = (c1z*(sf_p(iz_ix_iy  )-sf_p(izm1_ix_iy)) +c2z*(sf_p(izp1_ix_iy)-sf_p(izm2_ix_iy))) &
                      +(c1z*(sf_p(izp1_ix_iy)-sf_p(iz_ix_iy  )) +c2z*(sf_p(izp2_ix_iy)-sf_p(izm1_ix_iy)))
                !complete equation with unnecessary terms e.g. sf_p(iz_ix) for better understanding
                !with flag -O, the compiler should automatically detect such possibilities of simplification
                
                rvx = rf_vx(iz_ixp1_iy) +rf_vx(iz_ix_iy)
                rvy = rf_vy(iz_ix_iyp1) +rf_vy(iz_ix_iy)
                rvz = rf_vz(izp1_ix_iy) +rf_vz(iz_ix_iy)
                
                grad(j)=grad(j) + 0.25*( dsvx*rvx + dsvy*rvy + dsvz*rvz )
                
            enddo
            
        enddo
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine grad2d_flat_density(sf_p,rf_vx,rf_vz,&
                                   grad,            &
                                   ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_p,rf_vx,rf_vz
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        dsvx=0.; dsvz=0.
         rvx=0.; rvz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dsvx,dsvz,&
        !$omp          rvx, rvz)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers
                
                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2=i  -2*nz  !iz,ix-2
                iz_ixm1=i    -nz  !iz,ix-1
                iz_ixp1=i    +nz  !iz,ix+1
                iz_ixp2=i  +2*nz  !iz,ix+2
                
                dsvx = (c1x*(sf_p(iz_ix  )-sf_p(iz_ixm1)) +c2x*(sf_p(iz_ixp1)-sf_p(iz_ixm2))) &
                      +(c1x*(sf_p(iz_ixp1)-sf_p(iz_ix  )) +c2x*(sf_p(iz_ixp2)-sf_p(iz_ixm1)))
                dsvz = (c1z*(sf_p(iz_ix  )-sf_p(izm1_ix)) +c2z*(sf_p(izp1_ix)-sf_p(izm2_ix))) &
                      +(c1z*(sf_p(izp1_ix)-sf_p(iz_ix  )) +c2z*(sf_p(izp2_ix)-sf_p(izm1_ix)))
                !complete equation with unnecessary terms e.g. sf_p(iz_ix) for better understanding
                !with flag -O, the compiler should automatically detect such possibilities of simplification
                
                rvx = rf_vx(iz_ix) +rf_vx(iz_ixp1)
                rvz = rf_vz(iz_ix) +rf_vz(izp1_ix)
                
                grad(j)=grad(j) + 0.25*( dsvx*rvx + dsvz*rvz )
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine


    subroutine imag3d_flat(sf_p,rf_p,&
                           imag,                  &
                           ifz,ilz,ifx,ilx,ify,ily)
        real,dimension(*) :: sf_p,rf_p
        real,dimension(*) :: imag
        
        nz=cb%nz
        nx=cb%nx
        
        dvx_dx=0.
        dvy_dx=0.
        dvz_dz=0.
        
        sp=0.
        rp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,j,&
        !$omp         sp,rp)
        !$omp do schedule(dynamic) collapse(2)
        do iy=ify,ily
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(iy-cb%ify)*cb%nz*cb%nx+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+(iy-1)     *cb%mz*cb%mx+1 !grad has no boundary layers
                
                sp = sf_p(i)
                rp = rf_p(i)
                
                imag(j)=imag(j) + sp*rp
                
            end do
            
        end do
        end do
        !$omp end do
        !$omp end parallel
        
    end subroutine

    subroutine imag2d_flat(sf_p,rf_p,&
                           imag,            &
                           ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_p,rf_p
        real,dimension(*) :: imag
        
        nz=cb%nz
        
        sp=0.
        rp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         sp,rp)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
                sp = sf_p(i)
                rp = rf_p(i)
                
                imag(j)=imag(j) + sp*rp
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

end
