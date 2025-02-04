module m_propagator
use m_System
use m_hicks, only : hicks_r
use m_resampler
use m_model
use m_shotlist, only : shls
use m_shot
use m_computebox
use m_field
use m_correlate
use m_cpml

    private
    public :: GBRW, GFLOP

    !FD coef
    real,dimension(1),parameter :: coef = 1.
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z

    !local const
    real :: dt2,inv_2dz,inv_2dx

    !scaling source wavelet
    real :: wavelet_scaler

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain ISOtropic 2D constant-density ACoustic propagation'//s_NL// &
            '2nd-order Pressure formulation'//s_NL// &
            'Regular-Grid Finite-Difference (FDRG) method'//s_NL// &
            'Cartesian O(x²,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.5*Vmax/dx'//s_NL// &
            'Required model attributes: vp'//s_NL// &
            'Required field components: p, p_prev, p_next'

        integer :: nbndlayer=0

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:,:),allocatable :: scaled_vpz, scaled_vpx

        !time frames
        integer :: nt
        real :: dt

        !!reference value
        !real :: invsqEref

        contains
        procedure :: print_info
        procedure :: estim_RAM
        procedure :: check_model
        procedure :: check_discretization
        procedure :: init
        procedure :: init_field
        procedure :: init_correlate
        procedure :: init_abslayer

        procedure :: forward
        
        procedure :: inject_pressure
        ! procedure :: set_pressure
        procedure :: update_pressure
        procedure :: evolve_pressure
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    contains

   !========= for FDSG O(dx2,dt2) ===================  

    subroutine print_info(self)
        class(t_propagator) :: self

        call hud('Invoked field & propagator modules info : '//s_NL//self%info)
        call hud('FDSG Coef : 1') !//num2str(coef(1))//', '//num2str(coef(2)))
        
    end subroutine
    
    subroutine estim_RAM(self)
        class(t_propagator) :: self
    end subroutine
    
    subroutine check_model(self)
        class(t_propagator) :: self
        
        if(index(self%info,'vp')>0  .and. .not. allocated(m%vp)) then
            !call error('vp model is NOT given.')
            call alloc(m%vp,m%nz,m%nx,m%ny,o_init=1500.)
            call warn('Constant vp model (1500 m/s) is allocated by propagator.')
        endif

        if(.not. allocated(m%rho)) then
            call alloc(m%rho,m%nz,m%nx,m%ny,o_init=1000.)
            call warn('Constant rho model (1000 kg/m³) is allocated by propgator, but will NOT be used.')
        endif
                
    end subroutine

    subroutine check_discretization(self)
        class(t_propagator) :: self

        !grid dispersion condition
        if (5.*m%dmin > cb%velmin/shot%fmax) then  !O(x4) rule: 5 points per wavelength
            call warn(shot%sindex//' can have grid dispersion!'//s_NL// &
                ' 5*dz, velmin, fmax = '//num2str(5.*m%dmin)//', '//num2str(cb%velmin)//', '//num2str(shot%fmax))
        endif
        
        !time frames
        self%nt=shot%nt
        self%dt=shot%dt
        time_window=(shot%nt-1)*shot%dt

        sumcoef=1. !sum(abs(coef))

        CFL = sumcoef*cb%velmax*self%dt*m%rev_cell_diagonal !R. Courant, K. O. Friedrichs & H. Lewy (1928)

        call hud('CFL value: '//num2str(CFL))
        
        if(CFL>1.) then
            self%dt = setup%get_real('CFL',o_default='0.9')/(sumcoef*cb%velmax*m%rev_cell_diagonal)
            self%nt=nint(time_window/self%dt)+1

            call warn('CFL > 1 on '//shot%sindex//'!'//s_NL//&
                'vmax, dt, 1/dx = '//num2str(cb%velmax)//', '//num2str(self%dt)//', '//num2str(m%rev_cell_diagonal) //s_NL//&
                'Adjusted dt, nt = '//num2str(self%dt)//', '//num2str(self%nt))

        endif

    end subroutine

    subroutine init(self)
        class(t_propagator) :: self

        character(:),allocatable :: file

        ! c1x=coef(1)/(m%dx**2); c1z=coef(1)/(m%dz**2)

        dt2=self%dt**2
        ! inv_2dz=1/(2*m%dz)
        ! inv_2dx=1/(2*m%dx)

        wavelet_scaler=dt2/m%cell_volume

        if_hicks=shot%if_hicks

        call alloc(self%scaled_vpz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%scaled_vpx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        self%scaled_vpz = cb%vp*self%dt/m%dz
        self%scaled_vpx = cb%vp*self%dt/m%dx

        !initialize m_field
        call field_init(.false.,self%nt,self%dt)

        ! !initialize m_correlate
        ! call correlate_init(ppg%nt,ppg%dt)

        !rectified interval for time integration
        !default to Nyquist, and must be a multiple of dt
        rdt=setup%get_real('REF_RECT_TIME_INTEVAL','RDT',o_default=num2str(0.5/shot%fmax))
        irdt=floor(rdt/self%dt)
        if(irdt==0) irdt=1
        rdt=irdt*self%dt
        call hud('rdt, irdt = '//num2str(rdt)//', '//num2str(irdt))

    end subroutine

    subroutine init_field(self,f,name,ois_adjoint,oif_will_reconstruct)
        class(t_propagator) :: self
        type(t_field) :: f
        character(*) :: name
        logical,optional :: oif_will_reconstruct
        logical,optional :: ois_adjoint

        !field
        ! call f%init(name)
        f%name=name

        f%is_adjoint=either(ois_adjoint,.false.,present(ois_adjoint))

        call f%init_bloom

        !f%if_will_reconstruct=either(oif_will_reconstruct,.not.f%is_adjoint,present(oif_will_reconstruct))
        !if(f%if_will_reconstruct) call f%init_boundary
        call f%init_boundary_pressure

        call alloc(f%p     , [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%p_prev, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%p_next, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

    end subroutine

    subroutine init_correlate(self,corr,name,purpose)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name,purpose

    end subroutine

    subroutine init_abslayer(self)
        class(t_propagator) :: self

    end subroutine

    !========= Derivations =================
    !PDE:      A u = ϰ∂ₜ²u - ∇·b∇u = f
    !Adjoint:  Aᵀa = ϰ∂ₜ²a - ∇·b∇a = d
    !where
    !u=p=tr(s)=szz=sxx=syy is (hydrostatic) pressure
    !f=fp*δ(x-xs), d is recorded data
    !b=ρ⁻¹ is buoyancy, ϰ=κ⁻¹ is bulk compliance (inverse of modulus)
    !and a=pᵃ is the adjoint field
    !
    !Discrete case:
    !Meshing with staggered grids in time and space (2D example):
    !                     |       |   -½ ∂zp      |       |
    !                     |       |       bz      |       |
    !                     |       |       |       |       |
    !                     κ   bx  κ   bx  κ   bx  κ   bx  κ
    !  --p--p--p--→ t    -p--∂ₓp--p--∂ₓp--p--∂ₓp--p--∂ₓp--p-→ x
    !   -1  0  1         -2  -1½ -1  -½   0   ½   1   1½  2    
    !                     |       |       |       |       | 
    !                     |       |    ½ ∂zp      |       | 
    !                     |       |       bz      |       | 
    !                     |       |       |       |       | 
    !                    -|-------|-----1-p-------|-------|-
    !                     |       |       κ       |       | 
    !                     |       |       |       |       | 
    !                     |       |   1½ ∂zp      |       | 
    !                     |       |       bz      |       | 
    !                     |       |       |       |       | 
    !                                   z ↓
    !
    !Convention for half-integer index:
    !(array index)  =>    (real index)     
    !∂zp(iz,ix,iy)  => ∂zp[iz-½,ix,  iy  ]^n   :=vz((iz-½)*dz,ix*dx,iy*dy,n*dt)
    !∂ₓp(iz,ix,iy)  => ∂ₓp[iz,  ix-½,iy  ]^n  
    !∂yp(iz,ix,iy)  => ∂yp[iz,  ix,  iy-½]^n  
    !  p(iz,ix,iy)  =>   p[iz,  ix,  iy  ]^n+½ :=p(iz*dz,ix*dx,iy*dy,(n+½)*dt)
    !
    !Forward:
    !FD eqn:
    !ϰ*∂ₜ²p = ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp) +f
    !ϰ*∂ₜ²p = [∂zᶠ ∂ₓᶠ][bz  ][∂zᵇ]p
    !                  [  bx][∂ₓᵇ]
    !where
    !∂ₜ²*dt² := p^n+1 -2p^n +p^n-1  ~O(t²)
    !∂zᵇ*dz  := p(iz  )-p(iz-1)     ~O(x¹)
    !∂zᶠ*dz  := p(iz+1)-p(iz  )     ~O(x¹)
    !Step #1: p^n  += src
    !Step #2: sample p^n at receivers
    !Step #3: save p^n to boundary values
    !Step #4: p^n+1 = 2p^n -p^n-1 +laplacian of p^n
    !Step #5: (p^n-1,p^n) = (p^n,p^n+1)
    !in reverse time:
    !Step #5: (p^n,p^n+1) = (p^n-1,p^n)
    !Step #4: p^n-1 = 2p^n -p^n+1 +laplacian of p^n
    !Step #3: load boundary values for p^n+1
    !Step #1: p^n -= src
    !
    !Adjoint:
    !since
    !∂ₜ²ᵀ = ∂ₜ²
    !∂zᵇᵀ = p(iz)-p(iz+1) = -∂zᶠ, ∂zᶠᵀ = -∂zᵇ
    !FD eqn:
    !ϰ*∂ₜ²ᵀpᵃ = [∂zᵇᵀ ∂ₓᵇᵀ][bz  ][∂zᶠᵀ]pᵃ = [∂zᶠ ∂ₓᶠ][bz  ][∂zᵇ]pᵃ
    !                      [  bx][∂ₓᶠᵀ]              [  bx][∂ₓᵇ]  
    !ie. ϰ*∂ₜ²pᵃ = ∂zᶠbz*(∂zᵇpᵃ) + ∂ₓᶠbx*(∂ₓᵇpᵃ) +d
    !SAME as the discretized FD eqn!
    !
    !Time marching (in reverse time):
    !Step #1: pᵃ^n += adjsrc
    !Step #2: sample pᵃ^n at source
    !Step #4: pᵃ^n-1 = 2pᵃ^n -pᵃ^n+1 +laplacian of pᵃ^n
    !Step #5: (pᵃ^n,pᵃ^n+1) = (pᵃ^n-1,pᵃ^n)    

    subroutine forward(self,fld_u, o_u_star_u)
        class(t_propagator) :: self
        type(t_field) :: fld_u
        type(t_correlate),optional :: o_u_star_u

        real,parameter :: time_dir=1. !time direction

        !seismo
        call alloc(fld_u%seismo, shot%nrcv,self%nt)
#define RWs_init    cb%n

!         tt1=0.; tt2=0.; tt3=0.; tt4=0.; tt5=0.; tt6=0.

        ift=1; ilt=self%nt

!         call cpu_time(tic)

        do it=ift,ilt
!             if(mod(it,500)==0 .and. mpiworld%is_master) then
!                 write(*,*) 'it----',it
!                 call fld_u%check_value
!             endif

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add pressure
!             call cpu_time(tic)
            call self%inject_pressure(fld_u,time_dir,it)
!             call cpu_time(toc)
!             tt1=tt1+toc-tic

            !step 2: save p^it in boundary layers
            ! if(fld_u%if_will_reconstruct) then
                ! call cpu_time(tic)
                ! call fld_u%boundary_transport_pressure('save',it)
                ! call cpu_time(toc)
                ! tt2=tt2+toc-tic
            ! endif

            ! !step 3: set hardBC
            ! call cpu_time(tic)
            ! call self%set_pressure(fld_u,time_dir,it)
            ! call cpu_time(toc)
            ! tt3=tt3+toc-tic

            !step 4: update pressure
!             call cpu_time(tic)
            call self%update_pressure(fld_u,time_dir,it)
!             call cpu_time(toc)
!             tt4=tt4+toc-tic

            !step 5: evolve pressure, it -> it+1
!             call cpu_time(tic)
            call self%evolve_pressure(fld_u,time_dir,it)
!             call cpu_time(toc)
!             tt5=tt5+toc-tic

            !step 6: sample p^it+1 at receivers
!             call cpu_time(tic)
            call self%extract(fld_u,it)
!             call cpu_time(toc)
!             tt6=tt6+toc-tic

!             !snapshot
!             call fld_u%write(it)

        enddo

!         call cpu_time(toc)
!         tt1=tt1+toc-tic

!         if(mpiworld%is_master) then
!             write(*,*) 'Elapsed time to add source   ',tt1/mpiworld%max_threads
!             write(*,*) 'Elapsed time to save boundary',tt2/mpiworld%max_threads
!             ! write(*,*) 'Elapsed time to set field    ',tt3/mpiworld%max_threads
!             write(*,*) 'Elapsed time to update field ',tt4/mpiworld%max_threads
!             write(*,*) 'Elapsed time to evolve field ',tt5/mpiworld%max_threads
!             write(*,*) 'Elapsed time to extract field',tt6/mpiworld%max_threads
!         endif


        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')


        !postprocess
        if(present(o_u_star_u)) then
            !scale by m%cell_volume*rdt tobe an energy distribution in the discretized world
            call o_u_star_u%scale(m%cell_volume*rdt)
        endif

    end subroutine

    subroutine inject_pressure(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        if(.not. f%is_adjoint) then

!             if(if_hicks) then
!                 ifz=shot%src%ifz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
!                 ifx=shot%src%ifx-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
!                 ify=shot%src%ify-cb%ioy+1; ily=shot%src%ily-cb%ioy+1
!             else
                iz=shot%src%iz-cb%ioz+1
                ix=shot%src%ix-cb%iox+1
                iy=shot%src%iy-cb%ioy+1
!             endif
!
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler
!
!             !explosion
!             if(if_hicks) then
!                 f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) + wl*cb%vp(ifz:ilz,ifx:ilx,ify:ily)**2 *shot%src%interp_coef
!             else
                f%p(iz,ix,iy)                = f%p(iz,ix,iy)                + wl*cb%vp(iz,ix,iy)**2
!             endif

            return

        endif
#define RWs_inject    (9+3+4)
#define FLOPs_inject  (2+5)
        ! do i=1,shot%nrcv

        !     if(if_hicks) then
        !         ifz=shot%rcv(i)%ifz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
        !         ifx=shot%rcv(i)%ifx-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
        !         ify=shot%rcv(i)%ify-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1
        !     else
        !         iz=shot%rcv(i)%iz-cb%ioz+1
        !         ix=shot%rcv(i)%ix-cb%iox+1
        !         iy=shot%rcv(i)%iy-cb%ioy+1
        !     endif

        !     !adjsource for pressure
        !     wl = f%wavelet(i,it)*wavelet_scaler    !no time_dir needed!

        !     if(if_hicks) then 
        !         f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) +wl*cb%vp(ifz:ilz,ifx:ilx,ify:ily)**2*shot%rcv(i)%interp_coef
        !     else
        !         f%p(iz,ix,iy)                = f%p(iz,ix,iy)                +wl*cb%vp(iz,ix,iy)**2
        !     endif

        ! enddo
        
    end subroutine

    ! subroutine set_pressure(self,f,time_dir,it)
    !     class(t_propagator) :: self
    !     type(t_field) :: f

    !     if(.not. f%is_adjoint) then

    !         if(if_hicks) then
    !             ifz=shot%src%ifz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
    !             ifx=shot%src%ifx-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
    !             ify=shot%src%ify-cb%ioy+1; ily=shot%src%ily-cb%ioy+1
    !         else
    !             ifz=shot%src%iz-cb%ioz+1; ilz=ifz
    !             ifx=shot%src%ix-cb%iox+1; ilx=ifx
    !             ify=shot%src%iy-cb%ioy+1; ily=ify
    !         endif
            
    !         wl=f%wavelet(1,it)
            
    !         if(shot%src%comp=='pbnd') then !hard BC
    !             f%p(ifz:ilz,ifx:ilx,ify:ily) =                                wl!                                  *shot%src%interp_coef
    !         endif
                
    !         return

    !     endif

    !     do i=1,shot%nrcv

    !         if(if_hicks) then
    !             ifz=shot%rcv(i)%ifz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
    !             ifx=shot%rcv(i)%ifx-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
    !             ify=shot%rcv(i)%ify-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1
    !         else
    !             ifz=shot%rcv(i)%iz-cb%ioz+1; ilz=ifz
    !             ifx=shot%rcv(i)%ix-cb%iox+1; ilx=ifx
    !             ify=shot%rcv(i)%iy-cb%ioy+1; ily=ify
    !         endif

    !         ! !adjsource for pressure
    !         ! wl=f%wavelet(i,it)
                
    !         if(shot%rcv(i)%comp=='pbnd') then
    !             f%p(ifz:ilz,ifx:ilx,ify:ily) =                               0. !wl                                  !*shot%rcv(i)%interp_coef !no time_dir needed!
    !         endif

    !     enddo
        
    ! end subroutine


    subroutine update_pressure(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        real,dimension(:,:,:),allocatable :: tmp_p

        ifz=f%bloom(1,it)
        if(m%is_freesurface) ifz=max(ifz,1)

        if(m%is_cubic) then
            ! call fd3d_pressure(f%p,                                      &
            !                    f%dp_dz,f%dp_dx,f%dp_dy,                  &
            !                    self%buoz,self%buox,self%buoy,self%kpa,   &
            !                    ifz,f%bloom(2,it),f%bloom(3,it),f%bloom(4,it))
        else
            call fd2d(f%p_next,f%p,f%p_prev, &
                self%scaled_vpz, self%scaled_vpx, &
                ifz,f%bloom(2,it),f%bloom(3,it),f%bloom(4,it))

        endif


    end subroutine

    subroutine evolve_pressure(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        real,dimension(:,:,:),pointer :: tmp

        if(time_dir>0.) then !in forward time
            tmp=>f%p_prev
            f%p_prev => f%p
            f%p      => f%p_next
            f%p_next => tmp

        else !in reverse time
            tmp=>f%p_next
            f%p_next => f%p
            f%p      => f%p_prev
            f%p_prev => tmp

        endif
        
    end subroutine

    subroutine extract(self,f,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        if(.not.f%is_adjoint) then

            do i=1,shot%nrcv
!                 ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
!                 ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
!                 ify=shot%rcv(i)%ify-cb%ioy+1; iy=shot%rcv(i)%iy-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1
                iz=shot%rcv(i)%iz-cb%ioz+1;
                ix=shot%rcv(i)%ix-cb%iox+1;
                iy=shot%rcv(i)%iy-cb%ioy+1;

!                 if(if_hicks) then
!                     select case (shot%rcv(i)%comp)
!                         case default
!                         !case ('p')
!                         f%seismo(i,it)=sum(f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(i)%interp_coef)
!
!                         ! case ('vz')
!                         ! f%seismo(i,it)=sum(f%vz(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
!                         ! case ('vx')
!                         ! f%seismo(i,it)=sum(f%vx(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
!                         ! case ('vy')
!                         ! f%seismo(i,it)=sum(f%vy(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
!                     end select
!
!                 else
!                     select case (shot%rcv(i)%comp)
!                         case default
                        !case ('p') !p[iz,ix,iy]
                        f%seismo(i,it)=f%p(iz,ix,iy)

                        ! case ('vz') !vz[iz-0.5,ix,iy]
                        ! f%seismo(i,it)=f%vz(iz,ix,iy)
                        ! case ('vx') !vx[iz,ix-0.5,iy]
                        ! f%seismo(i,it)=f%vx(iz,ix,iy)
                        ! case ('vy') !vy[iz,ix,iy-0.5]
                        ! f%seismo(i,it)=f%vy(iz,ix,iy)
!                     end select
                    
!                 endif

            enddo

            return

        endif

#define RWs_extract   (9+2)*shot%nrcv
#define FLOPs_extract  0

!             ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
!             ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
!             ify=shot%src%ify-cb%ioy+1; iy=shot%src%iy-cb%ioy+1; ily=shot%src%ily-cb%ioy+1
!
!             if(if_hicks) then
!                 select case (shot%src%comp)
!                     case default
!                     !case ('p')
!                     f%seismo(1,it)=sum(f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coef)
!
!                     ! case ('vz')
!                     ! f%seismo(1,it)=sum(f%vz(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
!
!                     ! case ('vx')
!                     ! f%seismo(1,it)=sum(f%vx(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
!
!                     ! case ('vy')
!                     ! f%seismo(1,it)=sum(f%vy(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
!
!                 end select
!
!             else
!                 select case (shot%src%comp)
!                     case default
!                     !case ('p') !p[iz,ix,iy]
!                     f%seismo(1,it)=f%p(iz,ix,iy)
!
!                     ! case ('vz') !vz[iz-0.5,ix,iy]
!                     ! f%seismo(1,it)=f%vz(iz,ix,iy)
!
!                     ! case ('vx') !vx[iz,ix-0.5,iy]
!                     ! f%seismo(1,it)=f%vx(iz,ix,iy)
!
!                     ! case ('vy') !vy[iz,ix,iy-0.5]
!                     ! f%seismo(1,it)=f%vy(iz,ix,iy)
!
!                 end select
!
!             endif
        
    end subroutine
    
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%scaled_vpz, self%scaled_vpx)
    end subroutine



    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd2d(pn, pc, pp, scaled_vpz, scaled_vpx, ifz,ilz,ifx,ilx)
        real,dimension(*) :: pn, pc, pp
        real,dimension(*) :: scaled_vpz, scaled_vpx

        real lapz, lapx
        
        nz=cb%nz
        nx=cb%nx
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm1,      iz_ixp1,&
        !$omp         pc_izm1,pp_izm1,pc_izp1,pp_izp1,pc_ixm1,pp_ixm1,pc_ixp1,pp_ixp1,&
        !$omp         lapz,lapx)
        !$omp do schedule(dynamic)
        do ix = ifx,ilx
            !dir$ simd
            do iz = ifz,ilz

                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1

                iz_ixm1=i    -nz  !iz,ix-1
                izm1_ix=i-1       !iz+1,ix
                iz_ix  =i         !iz,ix
                izp1_ix=i+1       !iz+1,ix
                iz_ixp1=i    +nz  !iz,ix+1

                !branching
                if(iz>ifz) then
                    pc_izm1 = pc(izm1_ix)
                    pp_izm1 = pp(izm1_ix)
                else
                    pc_izm1 = 0.
                    pp_izm1 = 0.
                endif

                if(iz<cb%ilz) then
                    pc_izp1 = pc(izp1_ix)
                    pp_izp1 = pp(izp1_ix)
                else
                    pc_izp1 = 0.
                    pp_izp1 = 0.
                endif

                if(ix>cb%ifx) then
                    pc_ixm1 = pc(iz_ixm1)
                    pp_ixm1 = pp(iz_ixm1)
                else
                    pc_ixm1 = 0.
                    pp_ixm1 = 0.
                endif

                if(ix<cb%ilx) then
                    pc_ixp1 = pc(iz_ixp1)
                    pp_ixp1 = pp(iz_ixp1)
                else
                    pc_ixp1 = 0.
                    pp_ixp1 = 0.
                endif

#define RWs_fd2d_1    ( 2*cb%nz + (4*cb%n-(cb%nx+cb%nz-2)*2)*2 )
#define FLOPs_fd2d_1  0

                lapz=scaled_vpz(iz_ix)*scaled_vpz(iz_ix)*( pc_izp1 +pc_izm1 -2.*pc(iz_ix) )
                lapx=scaled_vpx(iz_ix)*scaled_vpx(iz_ix)*( pc_ixp1 +pc_ixm1 -2.*pc(iz_ix) )

#define RWs_fd2d_2    4*cb%n
#define FLOPs_fd2d_2  5*cb%n*2

                !! top boundary */
                if(m%is_freesurface) then
                if(iz==ifz) then
                    lapz=scaled_vpz(iz_ix)*(-pc(iz_ix)+pc_izp1+pp(iz_ix)-pp_izp1);
                    if(ix>cb%ifx .and. ix<cb%ilx) lapx=0.5*lapx;
                endif
                endif

                !! bottom boundary */
                if(iz==cb%ilz) then
                    lapz=scaled_vpz(iz_ix)*(pc_izm1-pc(iz_ix)-pp_izm1+pp(iz_ix));
                    if(ix>cb%ifx .and. ix<cb%ilx) lapx=0.5*lapx;
                endif

                !! left boundary */
                if(ix==cb%ifx) then
                    if(iz>ifz .and. iz<cb%ilz) lapz=0.5*lapz;
                    lapx=scaled_vpx(iz_ix)*(-pc(iz_ix)+pc_ixp1+pp(iz_ix)-pp_ixp1);
                endif

                !! right boundary */
                if(ix==cb%ilx) then
                    if(iz>ifz .and. iz<cb%ilz) lapz=0.5*lapz;
                    lapx=scaled_vpx(iz_ix)*(pc_iz_ixm1-pc(iz_ix)-pp_iz_ixm1+pp(iz_ix));
                endif

#define RWs_fd2d_3    3*(cb%nx+cb%nz-2)*2
#define FLOPs_fd2d_3  (4+1)*(cb%nx+cb%nz-2)*2

                pn(i) = 2.*pc(i) -pp(i) +lapz +lapx; !forward in time

#define RWs_fd2d_4    3*cb%n
#define FLOPs_fd2d_4  4*cb%n


            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine


    !giga bytes of data transferred between global & local memory on device
    pure real function GBRW()

        integer(8) :: total

        total = ( &
             int(RWs_fd2d_1+RWs_fd2d_2+RWs_fd2d_3+RWs_fd2d_4, 8) *ppg%nt*shls%nlists &
            +int(RWs_inject ,8) *ppg%nt*shls%nlists &
            +int(RWs_extract,8) *ppg%nt*shls%nlists &
            +int(RWs_init   ,8)        *shls%nlists &
            )

        GBRW = total*4/1e9

    end function

    !giga FP32 arithmics(+-*/) operations
    pure real function GFLOP()

        integer(8) :: total

        total = ( &
             int(FLOPs_fd2d_1+FLOPs_fd2d_2+FLOPs_fd2d_3+FLOPs_fd2d_4, 8) *ppg%nt*shls%nlists &
            +int(FLOPs_inject ,8) *ppg%nt*shls%nlists &
            +int(FLOPs_extract,8) *ppg%nt*shls%nlists &
            )

        GFLOP = total/1e9

    end function

end
