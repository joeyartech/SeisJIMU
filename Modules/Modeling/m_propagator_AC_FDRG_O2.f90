module m_propagator

use m_System
use m_hicks, only : hicks_r
use m_resampler
use m_model
use m_shot
use m_computebox
use m_field
use m_cpml

    private

    !FD coef
    real,dimension(1),parameter :: coef = 1.
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain ISOtropic 2D/3D ACoustic propagation'//s_NL// &
            '1st-order Velocity-Stress formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.5*Vmax/dx'//s_NL// &
            'Required model attributes: vp, rho'//s_NL// &
            'Required field components: vz, vx, vy(3D), p'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Imaging conditions: P-Pxcorr'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: grho gkpa'

        integer :: nbndlayer=max(1,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=2 !number of basic gradients
        integer :: nimag=1 !number of basic images
        integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:,:),allocatable :: buoz, buox, buoy, kpa, inv_kpa

        !time frames
        integer :: nt
        real :: dt

        contains
        procedure :: print_info
        procedure :: estim_RAM
        procedure :: check_model
        procedure :: check_discretization
        procedure :: init
        procedure :: init_field
        procedure :: init_abslayer

        procedure :: forward
        procedure :: adjoint
        
        procedure :: inject_pressure
        procedure :: update_pressure
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    !these procedures will be contained in an m_correlate module in future release
    public :: gradient_density, gradient_moduli, imaging, energy, gradient_postprocess, imaging_postprocess

    ! real,dimension(:,:,:),allocatable :: sf_p_save

    contains


   !========= for FDSG O(dx2,dt2) ===================  

    subroutine print_info(self)
        class(t_propagator) :: self

        call hud('Invoked field & propagator modules info : '//s_NL//self%info)
        call hud('FDRG Coef : 1') !//num2str(coef(1))//', '//num2str(coef(2)))
        
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

        if(index(self%info,'rho')>0 .and. .not. allocated(m%rho)) then
            call alloc(m%rho,m%nz,m%nx,m%ny,o_init=1000.)
            call warn('Constant rho model (1000 kg/m³) is allocated by propagator.')
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

        c1x=coef(1)/m%dx; c1y=coef(1)/m%dy; c1z=coef(1)/m%dz
        !c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz
        
        if_hicks=shot%if_hicks

        !call alloc(self%buoz,   [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%buox,   [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%buoy,   [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%kpa,    [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%inv_kpa,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        
        self%kpa=cb%rho*cb%vp**2
        self%inv_kpa=1./self%kpa

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

        !initialize m_field
        call field_init(.false.,self%nt,self%dt)

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
        call f%init_boundary

        call alloc(f%p_prev, [cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%p_curr, [cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%p_next, [cb%ifx,cb%ilx],[cb%ify,cb%ily])

        !call alloc(f%dp_dz, [cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dp_dx, [cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dp_dy, [cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%dpxx_dx, [cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dpyy_dy, [cb%ifx,cb%ilx],[cb%ify,cb%ily])

    end subroutine

    subroutine init_abslayer(self)
        class(t_propagator) :: self

        call cpml%init

    end subroutine

    !========= Derivations =================
    !PDE:      A u = ∂ₜ²u - κ∇·b∇u = f
    !Adjoint:  Aᵀa = Aa = d
    !where
    !u=p=tr(s)=szz=sxx=syy is (hydrostatic) pressure
    !f=fp*δ(x-xs), d is recorded data
    !b=ρ⁻¹ is buoyancy, and a=pᵃ is the adjoint field
    !
    !Discrete case:
    !Meshing with staggered grids in time and space (2D example):
    !                            |       |   -½ ∂zp      |       |
    !                            |       |       bz      |       |
    !                            |       |       |       |       |
    !                            κ   bx  κ   bx  κ   bx  κ   bx  κ
    !  -∂ₜp--p-∂ₜp-p-∂ₜp-→ t    -p--∂ₓp--p--∂ₓp--p--∂ₓp--p--∂ₓp--p-→ x
    !    -1 -½  0  ½  1         -2  -1½ -1  -½   0   ½   1   1½  2    
    !                            |       |       |       |       | 
    !                            |       |    ½ ∂zp      |       | 
    !                            |       |       bz      |       | 
    !                            |       |       |       |       | 
    !                           -|-------|-----1-p-------|-------|-
    !                            |       |       κ       |       | 
    !                            |       |       |       |       | 
    !                            |       |   1½ ∂zp      |       | 
    !                            |       |       bz      |       | 
    !                            |       |       |       |       | 
    !                                          z ↓
    !
    !Forward:
    !FD eqn:
    !∂ₜ²p = ∂zᵇ(bz*∂zᶠp) + ∂ₓᵇ(bx*∂ₓᶠp) +f
    !where
    !∂ₜ²*dt² := p^n+1 -2p^n +p^n-1  ~O(t²)
    !∂zᶠ*dz  := p(iz+1)-p(iz  )     ~O(x¹)
    !∂zᵇ*dz  := p(iz  )-p(iz-1)     ~O(x¹)
    !Step #1: p^n+½ += src
    !Step #2: p^n+1½ = p^n+½ + spatial FD of p^n
    !Step #3: sample p^n+1½ at receivers
    !Step #4: save p^n+1½ to boundary values
    !in reverse time:
    !Step #4: load boundary values for p^n+1½
    !Step #2: p^n+½ = p^n+1½ - spatial FD of p^n
    !Step #1: p^n+½ -= src
    !
    !Adjoint:
    !FD eqn:
    !∂ₜ²ᵀpᵃ = ∂zᵇᵀ(bz*∂zᶠᵀpᵃ) + ∂ₓᵇᵀ(bx*∂ₓᶠᵀpᵃ) +d
    !∂ₜ²ᵀ = ∂ₜ²
    !∂zᶠᵀ = p(iz-1)-p(iz) = -∂zᵇ, ∂zᵇᵀ = -∂zᶠ
    !so
    !∂ₜ²pᵃ = ∂zᶠ(bz*∂zᵇpᵃ) + ∂ₓᶠ(bx*∂ₓᵇpᵃ) +d
    !NOT exactly same as the discretized FD eqn!
    !
    !Time marching (in reverse time):
    !Step #5: pᵃ^n+1½ += adjsrc
    !Step #4: pᵃ^n+½ = pᵃ^n+1½ + spatial FD of pᵃ^n+1
    !
    !For adjoint test:
    !In each step of forward time marching: dsyn=RGA
    !!f:source wavelet, N=M⁻¹dt: diagonal 
    !A:inject source into field, G=UK∂ᵇB∂ᶠ:propagator, with
    !∂ᶠ: forward FD, ∂ᵇ: backward FD, U:2nd-ord time integration B:buoyancy, K:bulk modulus, 
    !R:extract field at receivers
    !while in each step of reverse-time adjoint marching: dadj=AᵀGᵀRᵀdsyn
    !Rᵀ:inject adjoint sources, Gᵀ:adjoint propagator, Aᵀ:extract adjoint fields
    !

    subroutine forward(self,fld_u)
        class(t_propagator) :: self
        type(t_field) :: fld_u

        real,parameter :: time_dir=1. !time direction

        !seismo
        call alloc(fld_u%seismo,shot%nrcv,self%nt)
            
        tt1=0.; tt2=0.; tt3=0.; tt4=0.; tt5=0.; tt6=0.

        ift=1; ilt=self%nt

        do it=ift,ilt
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_u%check_value
            endif

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add pressure
            call cpu_time(tic)
            call self%inject_stresses(fld_u,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2: update pressure
            call cpu_time(tic)
            call self%update_stresses(fld_u,time_dir,it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            !step 3: sample pressure at receivers
            call cpu_time(tic)
            call self%extract(fld_u,it)
            tt3=tt3+toc-tic

            !snapshot
            call fld_u%write(it)

            !step 6: save v^it+1 in boundary layers
            ! if(fld_u%if_will_reconstruct) then
                call cpu_time(tic)
                call fld_u%boundary_transport('save',it)
                call cpu_time(toc)
                tt4=tt4+toc-tic
            ! endif

            ! move new values to old values (the present becomes the past, the future becomes the present)
            f%p_prev = f%p_curr
            f%p_curr = f%p_next


        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add source',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary',tt4/mpiworld%max_threads
        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')

    end subroutine

    

    !forward: add RHS to s^it+0.5
    !adjoint: add RHS to s^it+1.5
    subroutine inject_pressure(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        if(.not. f%is_adjoint) then

            if(if_hicks) then
                ifz=shot%src%ifz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
                ifx=shot%src%ifx-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
                ify=shot%src%ify-cb%ioy+1; ily=shot%src%ily-cb%ioy+1
            else
                ifz=shot%src%iz-cb%ioz+1; ilz=ifz
                ifx=shot%src%ix-cb%iox+1; ilx=ifx
                ify=shot%src%iy-cb%ioy+1; ily=ify
            endif
            
            wl=time_dir*f%wavelet(1,it)
            
            if(shot%src%comp=='p') then !explosion
                f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) + wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef
            endif
            if(shot%src%comp=='pbnd') then !hard BC
                f%p(ifz:ilz,ifx:ilx,ify:ily) =                                wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef
            endif
                
            return

        endif

        do i=1,shot%nrcv

            if(if_hicks) then
                ifz=shot%rcv(i)%ifz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                ify=shot%rcv(i)%ify-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1
            else
                ifz=shot%rcv(i)%iz-cb%ioz+1; ilz=ifz
                ifx=shot%rcv(i)%ix-cb%iox+1; ilx=ifx
                ify=shot%rcv(i)%iy-cb%ioy+1; ily=ify
            endif

            !adjsource for pressure
            wl=f%wavelet(i,it)
                
            if(shot%rcv(i)%comp=='p') then
                f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) +wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef !no time_dir needed!
            endif
            if(shot%rcv(i)%comp=='pbnd') then
                f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) +wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef !no time_dir needed!
            endif

        enddo
        
    end subroutine

    !forward: s^it+0.5 -> s^it+1.5 by FD of v^it+1
    !adjoint: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
    subroutine update_stresses(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ifz=f%bloom(1,it)+1
        ilz=f%bloom(2,it)-2
        ifx=f%bloom(3,it)+1
        ilx=f%bloom(4,it)-2
        ify=f%bloom(5,it)+1
        ily=f%bloom(6,it)-2
        
        if(m%is_freesurface) ifz=max(ifz,1)

        if(m%is_cubic) then
            ! call fd3d_stresses(f%vz,f%vx,f%vy,f%p,                     &
            !                    f%dvz_dz,f%dvx_dx,f%dvy_dy,             &
            !                    self%kpa,                               &
            !                    ifz,ilz,ifx,ilx,ify,ily,time_dir*self%dt)
        else
            call fd2d_pressure(f%p,                            &
                               f%dp_dx,f%dp_dy,                &
                               self%buox,self%buoy,self%kpa,   &
                               ifx,ilx,time_dir*self%dt        )
        endif
        

        ! apply the time evolution scheme
        ! we apply it everywhere, including at some points on the edges of the domain that have not be calculated above,
        ! which is of course wrong (or more precisely undefined), but this does not matter because these values
        ! will be erased by the Dirichlet conditions set on these edges below
        f%p_next = - f%p_prev + 2*f%p_curr + &
              DELTAT*DELTAT * ((f%dpxx_dx + f%dpyy_dy) * kpa ) !+ &
              !4 * r_pi * vp**2 * source_term * Kronecker_source)

        ! !apply free surface boundary condition if needed
        ! if(m%is_freesurface) call fd_freesurface_stresses(f%p)

        ! apply Dirichlet conditions at the bottom of the C-PML layers,
        ! which is the right condition to implement in order for C-PML to remain stable at long times
        
        ! Dirichlet condition for pressure on the left boundary
        f%p_next(1,:) = 0.

        ! Dirichlet condition for pressure on the right boundary
        f%p_next(NX,:) = 0.

        ! Dirichlet condition for pressure on the bottom boundary
        f%p_next(:,1) = 0.

        ! Dirichlet condition for pressure on the top boundary
        f%p_next(:,NY) = 0.

    end subroutine

    subroutine extract(self,f,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        if(.not.f%is_adjoint) then

            do i=1,shot%nrcv
                ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                ify=shot%rcv(i)%ify-cb%ioy+1; iy=shot%rcv(i)%iy-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1

                if(if_hicks) then
                    select case (shot%rcv(i)%comp)
                        
                        case ('p')
                        f%seismo(i,it)=sum(f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(i)%interp_coef)
                        ! case ('vz')
                        ! f%seismo(i,it)=sum(f%vz(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
                        ! case ('vx')
                        ! f%seismo(i,it)=sum(f%vx(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
                        ! case ('vy')
                        ! f%seismo(i,it)=sum(f%vy(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                        case ('p') !p[iz,ix,iy]
                        f%seismo(i,it)=f%p(iz,ix,iy)
                        ! case ('vz') !vz[iz-0.5,ix,iy]
                        ! f%seismo(i,it)=f%vz(iz,ix,iy)
                        ! case ('vx') !vx[iz,ix-0.5,iy]
                        ! f%seismo(i,it)=f%vx(iz,ix,iy)
                        ! case ('vy') !vy[iz,ix,iy-0.5]
                        ! f%seismo(i,it)=f%vy(iz,ix,iy)
                    end select
                    
                endif

            enddo

            return

        endif

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            ify=shot%src%ify-cb%ioy+1; iy=shot%src%iy-cb%ioy+1; ily=shot%src%ily-cb%ioy+1
            
            if(if_hicks) then
                select case (shot%src%comp)
                    case ('p')
                    f%seismo(1,it)=sum(f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coef)
                    
                    case ('vz')
                    f%seismo(1,it)=sum(f%vz(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
                    
                    case ('vx')
                    f%seismo(1,it)=sum(f%vx(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
                    
                    case ('vy')
                    f%seismo(1,it)=sum(f%vy(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
                    
                end select
                
            else
                select case (shot%src%comp)
                    case ('p') !p[iz,ix,iy]
                    f%seismo(1,it)=f%p(iz,ix,iy)
                    
                    case ('vz') !vz[iz-0.5,ix,iy]
                    f%seismo(1,it)=f%vz(iz,ix,iy)
                    
                    case ('vx') !vx[iz,ix-0.5,iy]
                    f%seismo(1,it)=f%vx(iz,ix,iy)
                    
                    case ('vy') !vy[iz,ix,iy-0.5]
                    f%seismo(1,it)=f%vy(iz,ix,iy)
                    
                end select
                
            endif
        
    end subroutine
    
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%buoz, self%buox, self%buoy, self%kpa, self%inv_kpa)
    end subroutine



    subroutine fd2d_pressure(p,                &
                             dp_dx,dp_dy,      &
                             buox,buoy,kpa,    &
                             ifx,ilx,ify,ily,dt)
        real,dimension(*) :: p
        real,dimension(*) :: dp_dx,dp_dy
        real,dimension(*) :: buox,buoy
        
        nz=cb%nz
        nx=cb%nx
        
        dp_dx_=0.; dp_dy_=0.
        dpxx_dx_=0.; dpyy_dy_=0.

        !dir$ simd
        do j = 1,NY
            do i = 1,NX-1
                dp_dx_ = c1x*(f%p_curr(i+1,j) - f%p_curr(i,j))

                f%dp_dx(i,j) = cmpl%b_x_half(i)*f%dp_dx(i,j) + cmpl%a_x_half(i)*dp_dx_

                dp_dx_ = dp_dx_ / cmpl%kpa_x_half(i) + f%dp_dx(i,j)

                rho_half_x = 0.5 * (rho(i+1,j) + rho(i,j))
                f%pxx(i,j) = dp_dx_ / rho_half_x
            enddo
        enddo

        do j = 1,NY-1
            do i = 1,NX
                dp_dy_ = c1y*(f%p_curr(i,j+1) - f%p_curr(i,j))

                f%dp_dy(i,j) = cmpl%b_y_half(j)*f%dp_dy(i,j) + cmpl%a_y_half(j)*dp_dy_

                dp_dy_ = dp_dy_ / cmpl%kpa_y_half(j) + f%dp_dy(i,j)

                rho_half_y = 0.5 * (rho(i,j+1) + rho(i,j))
                f%pyy(i,j) = dp_dy_ / rho_half_y
            enddo
        enddo

        ! compute the second spatial derivatives

        do j = 1,NY
            do i = 2,NX
                dpxx_dx_ = c1x*(f%pxx(i,j) - f%pxx(i-1,j))

                f%dpxx_dx(i,j) = cmpl%b_x(i)*f%dpxx_dx(i,j) + cmpl%a_x(i)*dpxx_dx_

                dpxx_dx_ = dpxx_dx_ / cmpl%kpa_x(i) + f%dpxx_dx(i,j)

                f%dpxx_dx(i,j) = dpxx_dx_
            enddo
        enddo

        do j = 2,NY
            do i = 1,NX
                dpyy_dy_ = c1y*(f%pyy(i,j) - f%pyy(i,j-1))

                f%dpyy_dy(i,j) = cmpl%b_y(j) * f%dpyy_dy(i,j) + cmpl%a_y(j) * dpyy_dy_

                dpyy_dy_ = dpyy_dy_ / cmpl%kpa_y(j) + f%dpyy_dy(i,j)

                f%dpyy_dy(i,j) = dpyy_dy_
            enddo
        enddo

    end subroutine

end