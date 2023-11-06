module m_propagator
use m_System
use m_hicks, only : hicks_r
use m_resampler
use m_model
use m_shot
use m_computebox
use m_field
use m_correlate
use m_cpml

    private

    !FD coef
    real,dimension(2),parameter :: coef = [9./8.,-1./24.] !Fornberg, 1988, Generation of Finite-Difference Formulas on Arbitrary Spaced Grids.
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z

    !local const
    real :: dt2, inv_2dt, inv_2dz, inv_2dx


    !scaling source wavelet
    real :: wavelet_scaler

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain ISOtropic 2D/3D Acoustic propagation'//s_NL// &
            '2nd-order Pressure formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.5*Vmax/dx'//s_NL// &
            'Required model attributes: vp, rho'//s_NL// &
            'Required field components: p, p_prev, p_next'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Poynting definitions: Esq_gradphi'//s_NL// &
            'Imaging conditions: ipp ibksc ifwsc (P-Pxcorr of backward & forward scattering)'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: gikpa, gbuo'

        integer :: nbndlayer=max(1,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=2 !number of basic gradients
        integer :: nimag=3 !number of basic images
        !integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:,:),allocatable :: buoz, buox, buoy, kpa !r2, tilD_vp2

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
        procedure :: assemble

        procedure :: forward
        procedure :: adjoint_poynting
        
        procedure :: inject_pressure
        procedure :: update_pressure
        procedure :: evolve_pressure
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    character(:),allocatable :: s_poynting_def

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    ! logical :: is_absolute_virtual

    contains
    
    !========= for FDSG O(dx4,dt2) ===================  

    subroutine print_info(self)
        class(t_propagator) :: self

        call hud('Invoked field & propagator modules info : '//s_NL//self%info)
        call hud('FDSG Coef : '//num2str(coef(1))//', '//num2str(coef(2)))
        
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

        sumcoef=sum(abs(coef))

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

        c1x=coef(1)/m%dx; c1y=coef(1)/m%dy; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz

        dt2=self%dt**2
        inv_2dt =1./2/self%dt
        inv_2dz =1./2/m%dz
        inv_2dx =1./2/m%dx
        
        wavelet_scaler=dt2/m%cell_volume

        if_hicks=shot%if_hicks

        call alloc(self%buoz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%buox,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%buoy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%kpa, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        self%kpa=cb%rho*cb%vp**2

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

        ! call alloc(self%r2,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        ! self%r2 = (cb%vp*self%dt/m%dx)**2

        !initialize m_field
        call field_init(.false.,self%nt,self%dt)

        !initialize m_correlate
        call correlate_init(ppg%nt,ppg%dt)

        !rectified interval for time integration
        !default to Nyquist, and must be a multiple of dt
        rdt=setup%get_real('REF_RECT_TIME_INTEVAL','RDT',o_default=num2str(0.5/shot%fmax))
        irdt=floor(rdt/self%dt)
        if(irdt==0) irdt=1
        rdt=irdt*self%dt
        call hud('rdt, irdt = '//num2str(rdt)//', '//num2str(irdt))

        s_poynting_def=setup%get_str('POYNTING_DEF',o_default='Esq_gradphi')
        if(s_poynting_def/='Esq_gradphi') call error('Sorry, other Poynting definitions have not yet implemented.')

    end subroutine

    subroutine init_field(self,f,name,ois_simple,ois_adjoint,oif_will_reconstruct)
        class(t_propagator) :: self
        type(t_field) :: f
        character(*) :: name
        logical,optional :: oif_will_reconstruct
        logical,optional :: ois_adjoint, ois_simple

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

        ! if(.not.either(ois_simple,.false.,present(ois_simple))) then
            call alloc(f%dp_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
            call alloc(f%dp_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
            call alloc(f%dp_dy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

            call alloc(f%dpzz_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
            call alloc(f%dpxx_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
            call alloc(f%dpyy_dy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

            call alloc(f%lap, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        ! endif

        call alloc(f%poynz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%poynx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

    end subroutine

    subroutine init_correlate(self,corr,name)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name
        
        corr%name=name

        if(name(1:1)=='g') then !gradient components
            call alloc(corr%gikpa,m%nz,m%nx,m%ny)
            call alloc(corr%gbuo, m%nz,m%nx,m%ny)
        else !image components
            call alloc(corr%ipp,m%nz,m%nx,m%ny)
            call alloc(corr%ibksc,m%nz,m%nx,m%ny)
            call alloc(corr%ifwsc,m%nz,m%nx,m%ny)
        endif

    end subroutine

    subroutine init_abslayer(self)
        class(t_propagator) :: self

        call cpml%init

    end subroutine

    subroutine assemble(self,corr)
        class(t_propagator) :: self
        type(t_correlate) :: corr

        if(allocated(correlate_image)) then
            call correlate_assemble(corr%ipp,   correlate_image(:,:,:,1))
            call correlate_assemble(corr%ibksc, correlate_image(:,:,:,2))
            call correlate_assemble(corr%ifwsc, correlate_image(:,:,:,3))
        endif

        if(allocated(correlate_gradient)) then
            call correlate_assemble(corr%gikpa, correlate_gradient(:,:,:,1))
            call correlate_assemble(corr%gbuo,  correlate_gradient(:,:,:,2))
        endif        
        
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
        

    subroutine forward(self,fld_u)
        class(t_propagator) :: self
        type(t_field) :: fld_u

        real,parameter :: time_dir=1. !time direction

        !seismo
        call alloc(fld_u%seismo,shot%nrcv,self%nt)
            
        tt1=0.; tt2=0.; tt3=0.; tt4=0.; tt5=0.; tt6=0.; tt7=0.

        ift=1; ilt=self%nt

        do it=ift,ilt
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_u%check_value
            endif

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add pressure
            call cpu_time(tic)
            call self%inject_pressure(fld_u,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2: save p^it+1 in boundary layers
            ! if(fld_E0%if_will_reconstruct) then
                call cpu_time(tic)
                call fld_u%boundary_transport_pressure('save',it)
                call cpu_time(toc)
                tt2=tt2+toc-tic
            ! endif

            ! !step 3: set hardBC
            ! call cpu_time(tic)
            ! call self%set_pressure(fld_E0,time_dir,it)
            ! call cpu_time(toc)
            ! tt3=tt3+toc-tic

            !step 4: update pressure
            call cpu_time(tic)
            call self%update_pressure(fld_u,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !step 5: evolve pressure, it -> it+1
            call cpu_time(tic)
            call self%evolve_pressure(fld_u,time_dir,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !step 6: sample p^it+1 at receivers
            call cpu_time(tic)
            call self%extract(fld_u,it)
            call cpu_time(toc)
            tt7=tt7+toc-tic

            !snapshot
            call fld_u%write(it)

        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add source              ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary           ',tt2/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set field             ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field            ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve field            ',tt6/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field           ',tt7/mpiworld%max_threads
        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')

    end subroutine

    subroutine adjoint_poynting(self,fld_q,fld_p, fld_v,fld_u, a_star_u)
        class(t_propagator) :: self
        type(t_field) :: fld_q,fld_p,fld_v,fld_u
        type(t_correlate) :: a_star_u

        real,parameter :: time_dir=-1. !time direction

        !reinitialize absorbing boundary for incident wavefield reconstruction
        call fld_v%reinit
        call fld_u%reinit
        
        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.
        
        ift=1; ilt=self%nt

        do it=ilt,ift,int(time_dir)
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_q%check_value
                call fld_p%check_value
                call fld_v%check_value
                call fld_u%check_value
            endif

            ! if(present(o_sf)) then
                !backward step 5: it+1 -> it
                call cpu_time(tic)
                call self%evolve_pressure(fld_v,time_dir,it)
                call self%evolve_pressure(fld_u,time_dir,it)
                call cpu_time(toc)
                tt1=tt1+toc-tic

                ! !backward step 2: retrieve p^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_v%boundary_transport_pressure('load',it)
                call fld_u%boundary_transport_pressure('load',it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 4:
                call cpu_time(tic)
                call self%update_pressure(fld_v,time_dir,it)
                call self%update_pressure(fld_u,time_dir,it)
                call cpu_time(toc)
                tt4=tt4+toc-tic

                !backward step 1: rm p^it at source
                call cpu_time(tic)
                call self%inject_pressure(fld_v,time_dir,it)
                call self%inject_pressure(fld_u,time_dir,it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            ! endif

            !adjoint step 6: inject to p^it+1 at receivers
            call cpu_time(tic)
            call self%inject_pressure(fld_q,time_dir,it)
            call self%inject_pressure(fld_p,time_dir,it)
            call cpu_time(toc)
            tt8=tt8+toc-tic

            !adjoint step 4:
            call cpu_time(tic)
            call self%update_pressure(fld_q,time_dir,it)
            call self%update_pressure(fld_p,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !image: rf%p^it star sf%p^it
            if(mod(it,irdt)==0) then
                call cpu_time(tic)
                call compute_poynting(fld_q,fld_p)
                call compute_poynting(fld_v,fld_u)
                call cpu_time(toc)
                tt3=tt3+toc-tic

                call cpu_time(tic)
                call cross_correlate_image(fld_p,fld_u,a_star_u,it)
                call cpu_time(toc)
                tt10=tt10+toc-tic
            endif

            !adjoint step 5
            ! this step is moved to update_pressure for easier management
            call cpu_time(tic)
            call self%evolve_pressure(fld_q,time_dir,it)
            call self%evolve_pressure(fld_p,time_dir,it)
            call cpu_time(toc)
            tt11=tt11+toc-tic

            !adjoint step 1: sample p^it at source position
            ! if(if_record_adjseismo) then
                ! call cpu_time(tic)
                ! call self%extract(fld_q,it)
                ! call self%extract(fld_p,it)
                ! call cpu_time(toc)
                ! tt12=tt12+toc-tic
            ! endif


            !--------------------------------------------------------!
            
            !snapshot
            call fld_p%write(it,o_suffix='_rev')
            call fld_u%write(it,o_suffix='_rev')

            call a_star_u%write(it,o_suffix='_rev')

        enddo

        call a_star_u%scale(m%cell_volume*rdt)


        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to evolve field        ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to load boundary       ',tt2/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set field           ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field        ',tt4/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to separate field      ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source           ',tt6/mpiworld%max_threads        
            write(*,*) 'Elapsed time -----------------------'
            write(*,*) 'Elapsed time to add adj & virtual src    ',tt8/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj field         ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve adj field         ',tt11/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set adj field          ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract fields           ',tt12/mpiworld%max_threads
            write(*,*) 'Elapsed time to compute Poynting vectors ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate                ',tt10/mpiworld%max_threads

        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_rfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_rfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')
        call hud('ximage < snap_*  n1='//num2str(cb%mz)//' perc=99')
        call hud('xmovie < snap_*  n1='//num2str(cb%mz)//' n2='//num2str(cb%mx)//' clip=?e-?? loop=2 title=%g')

    end subroutine


    subroutine inject_pressure(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        if(.not. f%is_adjoint) then
!this loop takes more time when nthreads>1.
!e.g. 0.1s vs 0.046s from Elapased time to rm source (tt5)
!tobe fixed..

            if(if_hicks) then
                ifz=shot%src%ifz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
                ifx=shot%src%ifx-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
                ify=shot%src%ify-cb%ioy+1; ily=shot%src%ily-cb%ioy+1
            else
                iz=shot%src%iz-cb%ioz+1
                ix=shot%src%ix-cb%iox+1
                iy=shot%src%iy-cb%ioy+1
            endif
            
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler

            !explosion
            if(if_hicks) then
                f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) + wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef
            else
                f%p(iz,ix,iy)                = f%p(iz,ix,iy)                + wl*self%kpa(iz,ix,iy)
            endif

            return

        endif

!this loop takes more time when nthreads>1.
!e.g. 38s vs 8.7s from Elapased time to add adj source  (tt6)
!tobe fixed..
        do i=1,shot%nrcv

            if(if_hicks) then
                ifz=shot%rcv(i)%ifz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                ify=shot%rcv(i)%ify-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1
            else
                iz=shot%rcv(i)%iz-cb%ioz+1
                ix=shot%rcv(i)%ix-cb%iox+1
                iy=shot%rcv(i)%iy-cb%ioy+1
            endif

            !adjsource for pressure
            wl = f%wavelet(i,it)*wavelet_scaler    !no time_dir needed!

            if(if_hicks) then 
                f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) +wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef
            else
                f%p(iz,ix,iy)                = f%p(iz,ix,iy)                +wl*self%kpa(iz,ix,iy)
            endif

        enddo
        
    end subroutine

    subroutine update_pressure(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ! !necessary after computing the secondary source
        ! f%lap=0.

        ifz=f%bloom(1,it)
        if(m%is_freesurface) ifz=max(ifz,1)
        ilz=f%bloom(2,it)
        ifx=f%bloom(3,it)
        ilx=f%bloom(4,it)

        if(m%is_cubic) then
            ! call fd3d_pressure(f%p,                                      &
            !                    f%dp_dz,f%dp_dx,f%dp_dy,                  &
            !                    self%buoz,self%buox,self%buoy,self%kpa,   &
            !                    ifz,f%bloom(2,it),f%bloom(3,it),f%bloom(4,it))
        else
            call fd2d_laplacian(f%p,                            &
                                f%dp_dz,f%dp_dx,f%dpzz_dz,f%dpxx_dx,&
                                f%lap,&
                                self%buoz,self%buox,            &
                                ifz,ilz,ifx,ilx)
        endif


        if(time_dir>0.) then !in forward time
            f%p_next(ifz:ilz,ifx:ilx,:) = 2*f%p(ifz:ilz,ifx:ilx,:) -f%p_prev(ifz:ilz,ifx:ilx,:) & 
                +dt2*self%kpa(ifz:ilz,ifx:ilx,:)*f%lap(ifz:ilz,ifx:ilx,:)
        else !in reverse time
            f%p_prev(ifz:ilz,ifx:ilx,:) = 2*f%p(ifz:ilz,ifx:ilx,:) -f%p_next(ifz:ilz,ifx:ilx,:) &
                +dt2*self%kpa(ifz:ilz,ifx:ilx,:)*f%lap(ifz:ilz,ifx:ilx,:)
        endif

        ! !apply free surface boundary condition if needed
        ! if(m%is_freesurface) call fd_freesurface_stresses(f%p)

        ! ! apply Dirichlet conditions at the bottom of the C-PML layers,
        ! ! the right condition to keep C-PML stable at long time
        ! f%p_next(cb%ifz,:,:)=0.
        ! f%p_next(cb%ilz,:,:)=0.
        ! f%p_next(:,cb%ifx,:)=0.
        ! f%p_next(:,cb%ilx,:)=0.
        ! f%p_next(:,:,cb%ify)=0.
        ! f%p_next(:,:,cb%ily)=0.

    end subroutine

    subroutine evolve_pressure(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        real,dimension(:,:,:),pointer :: tmp

        if(time_dir>0.) then !in forward time
            tmp=>f%p_prev
            f%p_prev=>f%p
            f%p     =>f%p_next
            f%p_next=>tmp

            ! f%p_prev = f%p
            ! f%p      = f%p_next
            ! !f%p_next = f%p_prev

        else !in reverse time
            tmp=>f%p_next
            f%p_next=>f%p
            f%p     =>f%p_prev
            f%p_prev=>tmp

            ! f%p_next = f%p
            ! f%p      = f%p_prev
            ! !f%p_prev = f%p_next

        endif

        nullify(tmp)
        
    end subroutine

    subroutine compute_poynting(v,u)
        type(t_field) :: v, u

        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: E2, ph !envelope squared & inst phase
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: dph_dz, dph_dx

    	!E=sqrt(u%p*u%p+v%p*v%p)
        E2=u%p*u%p+v%p*v%p
        ph=atan2(v%p,u%p)

    	do ix=cb%ifx+1,cb%ilx-1
    	do iz=cb%ifz+1,cb%ilz-1
    	    dph_dz(iz,ix,1) = asin(sin(ph(iz+1,ix,1) - ph(iz-1,ix,1)))*inv_2dz
    	    dph_dx(iz,ix,1) = asin(sin(ph(iz,ix+1,1) - ph(iz,ix-1,1)))*inv_2dx
    	enddo
    	enddo

        u%poynz=E2*dph_dz
        u%poynx=E2*dph_dx

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
                        case default
                        !case ('p')
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
                        case default
                        !case ('p') !p[iz,ix,iy]
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
                    case default
                    !case ('p')
                    f%seismo(1,it)=sum(f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coef)
                    
                    ! case ('vz')
                    ! f%seismo(1,it)=sum(f%vz(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
                    
                    ! case ('vx')
                    ! f%seismo(1,it)=sum(f%vx(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
                    
                    ! case ('vy')
                    ! f%seismo(1,it)=sum(f%vy(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
                    
                end select
                
            else
                select case (shot%src%comp)
                    case default
                    !case ('p') !p[iz,ix,iy]
                    f%seismo(1,it)=f%p(iz,ix,iy)
                    
                    ! case ('vz') !vz[iz-0.5,ix,iy]
                    ! f%seismo(1,it)=f%vz(iz,ix,iy)
                    
                    ! case ('vx') !vx[iz,ix-0.5,iy]
                    ! f%seismo(1,it)=f%vx(iz,ix,iy)
                    
                    ! case ('vy') !vy[iz,ix,iy-0.5]
                    ! f%seismo(1,it)=f%vy(iz,ix,iy)
                    
                end select
                
            endif
        
    end subroutine
        
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%buoz, self%buox, self%buoy, self%kpa)
    end subroutine


    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !Kₘ<a|Au> = Kₘ<a|ϰ∂ₜ²u - ∇·b∇u>
    !for ϰ: Kₘ<a|Au> = ∫ a ∂ₜ²u dt =-∫ ∂ₜa ∂ₜu dt, or = ∫ a κ∇·b∇u dt 
    !for b: Kₘ<a|Au> = -Kₘ<a|∇·b∇u> = ∫ ∇a·∇u dt

    subroutine auto_correlate(f,corr,it)
        type(t_field), intent(in) :: f
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=f%bloom(1,it)
        ilz=f%bloom(2,it)
        ifx=f%bloom(3,it)
        ilx=f%bloom(4,it)
        ify=f%bloom(5,it)
        ily=f%bloom(6,it)
        
        ! if(m%is_cubic) then
        ! else
            
            ! call imag2d_inverse_scattering(rf%p_next,rf%p,rf%p_prev,sf%p_next,sf%p,sf%p_prev,&
            !                       imag%drp_dt_dsp_dt,imag%nab_rp_nab_sp,                     &
            !                       ifz,ilz,ifx,ilx)
        ! endif


    end subroutine

    subroutine cross_correlate(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ! ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ! ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)

        !for gikpa
        corr%gikpa = corr%gikpa + &
            !rf%p(1:m%nz,1:m%nx,1:m%ny) * ( sf%p_next(1:m%nz,1:m%nx,1:m%ny)-2*sf%p(1:m%nz,1:m%nx,1:m%ny)+sf%p_prev(1:m%nz,1:m%nx,1:m%ny) )/dt2
            rf%p(1:m%nz,1:m%nx,1:m%ny) * ppg%kpa(1:m%nz,1:m%nx,1:m%ny)*sf%lap(1:m%nz,1:m%nx,1:m%ny)

        !for gbuo
        call fd2d_grho(rf%p,sf%p,corr%gbuo,   ifz,ilz,ifx,ilx)

    end subroutine

    subroutine cross_correlate_image(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ! ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ! ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)

        call imag2d(rf%p,sf%p,&
                    rf%poynz,rf%poynx,sf%poynz,sf%poynx, &
                    corr%ipp,corr%ibksc,corr%ifwsc, &
                    ifz,ilz,ifx,ilx)

    end subroutine

    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd2d_laplacian(p,dp_dz,dp_dx,&
                                dpzz_dz,dpxx_dx,&
                                lap,&
                                buoz,buox,&
                                ifz,ilz,ifx,ilx)
        real,dimension(*) :: p,dp_dz,dp_dx
        real,dimension(*) :: dpzz_dz,dpxx_dx
        real,dimension(*) :: lap
        real,dimension(*) :: buoz,buox

        real,dimension(:),allocatable :: pzz,pxx
        call alloc(pzz,cb%n)
        call alloc(pxx,cb%n)

        nz=cb%nz
        nx=cb%nx
        
        !flux: b∇u ~= ( bz*∂zᵇp , bx*∂ₓᵇp )
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dp_dz_,dp_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx+2,ilx-1
            !dir$ simd
            do iz = ifz+2,ilz-1

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                
                iz_ixm2=i  -2*nz !iz,ix-2
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1

                dp_dz_ = c1z*(p(iz_ix) - p(izm1_ix)) +c2z*(p(izp1_ix)-p(izm2_ix))
                dp_dx_ = c1x*(p(iz_ix) - p(iz_ixm1)) +c2x*(p(iz_ixp1)-p(iz_ixm2))

                dp_dz(iz_ix) = cpml%b_z_half(iz)*dp_dz(iz_ix) + cpml%a_z_half(iz)*dp_dz_
                dp_dx(iz_ix) = cpml%b_x_half(ix)*dp_dx(iz_ix) + cpml%a_x_half(ix)*dp_dx_

                dp_dz_ = dp_dz_/cpml%kpa_z_half(iz) + dp_dz(iz_ix)
                dp_dx_ = dp_dx_/cpml%kpa_x_half(ix) + dp_dx(iz_ix)

                pzz(iz_ix) = buoz(iz_ix)*dp_dz_
                pxx(iz_ix) = buox(iz_ix)*dp_dx_

            enddo
        enddo
        !$omp end do
        !$omp end parallel
        
        !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dpzz_dz_,dpxx_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx+1,ilx-2
            !dir$ simd
            do iz = ifz+1,ilz-2

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1
                iz_ixp2=i  +2*nz !iz,ix+2
                
                dpzz_dz_ = c1z*(pzz(izp1_ix) - pzz(iz_ix))  +c2z*(pzz(izp2_ix) - pzz(izm1_ix))
                dpxx_dx_ = c1x*(pxx(iz_ixp1) - pxx(iz_ix))  +c2x*(pxx(iz_ixp2) - pxx(iz_ixm1))

                dpzz_dz(iz_ix) = cpml%b_z(iz)*dpzz_dz(iz_ix) + cpml%a_z(iz)*dpzz_dz_
                dpxx_dx(iz_ix) = cpml%b_x(ix)*dpxx_dx(iz_ix) + cpml%a_x(ix)*dpxx_dx_

                dpzz_dz_ = dpzz_dz_/cpml%kpa_z(iz) + dpzz_dz(iz_ix)
                dpxx_dx_ = dpxx_dx_/cpml%kpa_x(ix) + dpxx_dx(iz_ix)

                lap(iz_ix) = dpzz_dz_ + dpxx_dx_

            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine fd2d_laplacian_nocpml(p,&
                                    lap,&
                                    buoz,buox,&
                                    ifz,ilz,ifx,ilx)
        real,dimension(*) :: p
        real,dimension(*) :: lap
        real,dimension(*) :: buoz,buox

        real,dimension(:),allocatable :: pzz,pxx
        call alloc(pzz,cb%n)
        call alloc(pxx,cb%n)

        nz=cb%nz
        nx=cb%nx
        
        !flux: b∇u ~= ( bz*∂zᵇp , bx*∂ₓᵇp )
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dp_dz_,dp_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx+2,ilx-1
            !dir$ simd
            do iz = ifz+2,ilz-1

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                
                iz_ixm2=i  -2*nz !iz,ix-2
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1

                dp_dz_ = c1z*(p(iz_ix) - p(izm1_ix)) +c2z*(p(izp1_ix)-p(izm2_ix))
                dp_dx_ = c1x*(p(iz_ix) - p(iz_ixm1)) +c2x*(p(iz_ixp1)-p(iz_ixm2))

                pzz(iz_ix) = buoz(iz_ix)*dp_dz_
                pxx(iz_ix) = buox(iz_ix)*dp_dx_

            enddo
        enddo
        !$omp end do
        !$omp end parallel
        
        !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dpzz_dz_,dpxx_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx+1,ilx-2
            !dir$ simd
            do iz = ifz+1,ilz-2

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1
                iz_ixp2=i  +2*nz !iz,ix+2
                
                dpzz_dz_ = c1z*(pzz(izp1_ix) - pzz(iz_ix))  +c2z*(pzz(izp2_ix) - pzz(izm1_ix))
                dpxx_dx_ = c1x*(pxx(iz_ixp1) - pxx(iz_ix))  +c2x*(pxx(iz_ixp2) - pxx(iz_ixm1))

                lap(iz_ix) = dpzz_dz_ + dpxx_dx_

            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

    ! !2nd order
    ! subroutine fd2d_grho(rp,sp,corr,&
    !                     ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: rp, sp
    !     real,dimension(*) :: corr

    !     real,dimension(:),allocatable :: pzz, pxx
    !     call alloc(pzz,cb%n)
    !     call alloc(pxx,cb%n)

    !     nz=cb%nz
    !     nx=cb%nx
        
        
    !     !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
    !     !$omp parallel default (shared)&
    !     !$omp private(iz,ix,i,&
    !     !$omp         iz_ixm1,izm1_ix,iz_ix,izp1_ix,iz_ixp1)
    !     !$omp do schedule(dynamic)
    !     do ix = ifx,ilx-1
    !         !dir$ simd
    !         do iz = ifz,ilz-1

    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
    !             j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers

    !             iz_ixm1=i    -nz  !iz,ix-1
    !             izm1_ix=i-1       !iz-1,ix
    !             iz_ix  =i
    !             izp1_ix=i+1       !iz+1,ix
    !             iz_ixp1=i    +nz  !iz,ix+1

                
    !             corr(j) = corr(j) &
    !                 + (rp(izp1_ix)-rp(izm1_ix))*(sp(izp1_ix)-sp(izm1_ix))*inv_2dz*inv_2dz &
    !                 + (rp(iz_ixp1)-rp(iz_ixm1))*(sp(iz_ixp1)-sp(iz_ixm1))*inv_2dx*inv_2dx

    !         enddo
    !     enddo
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

    !4th order
    subroutine fd2d_grho(rp,sp,corr,&
                        ifz,ilz,ifx,ilx)
        real,dimension(*) :: rp, sp
        real,dimension(*) :: corr

        real c1z, c2z, c1x, c2x !no confusion with c1z etc defined in the modules

        real,dimension(:),allocatable :: pzz, pxx
        call alloc(pzz,cb%n)
        call alloc(pxx,cb%n)

        nz=cb%nz
        nx=cb%nx
        
        c1z=2./3./m%dz; c2z=-1./12./m%dz
        c1x=2./3./m%dx; c2x=-1./12./m%dx
        
        !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         iz_ixm1,izm1_ix,iz_ix,izp1_ix,iz_ixp1)
        !$omp do schedule(dynamic)
        do ix = ifx+2,ilx-2
            !dir$ simd
            do iz = ifz+2,ilz-2

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers

                iz_ixm2=i    -2*nz
                iz_ixm1=i    -nz  !iz,ix-1
                izm2_ix=i-2
                izm1_ix=i-1       !iz-1,ix
                iz_ix  =i
                izp1_ix=i+1       !iz+1,ix
                izp2_ix=i+2
                iz_ixp1=i    +nz  !iz,ix+1
                iz_ixp2=i    +2*nz

                drp_dz = c1z*(rp(izp1_ix)-rp(izm1_ix)) +c2z*(rp(izp2_ix)-rp(izm2_ix))
                drp_dx = c1x*(rp(iz_ixp1)-rp(iz_ixm1)) +c2x*(rp(iz_ixp2)-rp(iz_ixm2))
                
                dsp_dz = c1z*(sp(izp1_ix)-sp(izm1_ix)) +c2z*(sp(izp2_ix)-sp(izm2_ix))
                dsp_dx = c1x*(sp(iz_ixp1)-sp(iz_ixm1)) +c2x*(sp(iz_ixp2)-sp(iz_ixm2))
                
                corr(j) = corr(j) + (drp_dz*dsp_dz + drp_dx*dsp_dx)
            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine imag2d(rf_p,sf_p,&
                        rf_poynz,rf_poynx,sf_poynz,sf_poynx,&
                        ipp, ibksc, ifwsc,&
                        ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_p,sf_p
        real,dimension(*) :: rf_poynz,rf_poynx,sf_poynz,sf_poynx
        real,dimension(*) :: ipp, ibksc, ifwsc
        
        nz=cb%nz
        
        rp=0.
        sp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         rp,sp)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers
                
                ipp(j)=ipp(j) + rf_p(i)*sf_p(i)

                if(rf_poynz(i)*sf_poynz(i)+rf_poynx(i)*sf_poynx(i) < 0.) then !backward scattering
                    ibksc(j)=ibksc(j) + rf_p(i)*sf_p(i)
                else
                    ifwsc(j)=ifwsc(j) + rf_p(i)*sf_p(i)
                endif
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

end