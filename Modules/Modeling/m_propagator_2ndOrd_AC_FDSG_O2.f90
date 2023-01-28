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
            '2nd-order Pressure formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x²,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.5*Vmax/dx'//s_NL// &
            'Required model attributes: vp, rho'//s_NL// &
            'Required field components: vz, vx, vy(3D), p'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Imaging conditions: P-Pxcorr'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: grho gkpa gD'

        integer :: nbndlayer=max(1,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=3 !number of basic gradients
        integer :: nimag=1 !number of basic images
        integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:,:),allocatable :: buoz, buox, buoy, kpa, inv_kpa

        !time frames
        integer :: nt
        real :: dt, dt2

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
        procedure :: set_pressure
        procedure :: update_pressure
        procedure :: evolve_pressure
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    real,dimension(:,:,:),allocatable :: D

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

        self%dt2=self%dt*self%dt
        
    end subroutine

    subroutine init(self)
        class(t_propagator) :: self

        character(:),allocatable :: file

        c1x=coef(1)/m%dx; c1y=coef(1)/m%dy; c1z=coef(1)/m%dz
        !c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz
        
        if_hicks=shot%if_hicks

        call alloc(self%buoz,   [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
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


        !if the "D" model exists
        file=setup%get_file('FILE_D')
        if(file/='') then
            call alloc(D,m%nz,m%nx,m%ny)
            call sysio_read(file,D,m%n)
            call hud('D model is read.')
        endif


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
        call f%init_boundary_pressure

        call alloc(f%p     , [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%p_prev, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%p_next, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%dp_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dp_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dp_dy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%dpzz_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dpxx_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dpyy_dy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%lapz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%lapx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%lapy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

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
    !∂ₜ²p = κ*∂zᶠ(bz*∂zᵇp) + κ*∂ₓᶠ(bx*∂ₓᵇp) +f
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
    !FD eqn:
    !∂ₜ²ᵀpᵃ = ∂zᶠᵀbz*(∂zᵇᵀκ*pᵃ) + ∂ₓᶠᵀ(bx*∂ₓᵇᵀκ*pᵃ) +d
    !∂ₜ²ᵀ = ∂ₜ²
    !∂zᵇᵀ = p(iz)-p(iz+1) = -∂zᶠ, ∂zᶠᵀ = -∂zᵇ
    !so
    !∂ₜ²pᵃ = ∂zᵇbz*(∂zᶠκ*pᵃ) + ∂ₓᵇbx*(∂ₓᶠκ*pᵃ) +d
    !NOT same as the discretized FD eqn!
    !
    !Time marching (in reverse time):
    !Step #1: pᵃ^n += adjsrc
    !Step #2: sample pᵃ^n at source
    !Step #4: pᵃ^n-1 = 2pᵃ^n -pᵃ^n+1 +laplacian of pᵃ^n
    !Step #5: (pᵃ^n,pᵃ^n+1) = (pᵃ^n-1,pᵃ^n)
    !
    !For adjoint test:
    !In each step of forward time marching: dsyn=REGBndA
    !A:inject source into field, Bnd:apply boundary condition, G=IK∂ᵇB∂ᶠ:propagator, with
    !∂ᶠ: forward FD, ∂ᵇ: backward FD, I:2nd-ord time integration B:buoyancy, K:bulk modulus, 
    !E:evolution, R:extract field at receivers
    !while in each step of reverse-time adjoint marching: dadj=AᵀBndᵀGᵀEᵀRᵀdsyn
    !Rᵀ:inject adjoint sources, Eᵀ:reverse-time evolution, Gᵀ:adjoint propagator, with
    !Gᵀ=∂ᶠᵀBᵀ∂ᵇᵀKᵀIᵀ=∂ᵇBᵀ∂ᶠKIᵀ
    !Bndᵀ:adjoint BC: pᵃ=0 where ∇p≠0 (Neumann BC); ∇pᵃ=0 where p≠0 (Dirichlet BC))
    !Aᵀ:extract adjoint fields

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
            call self%inject_pressure(fld_u,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2: save v^it+1 in boundary layers
            ! if(fld_u%if_will_reconstruct) then
                call cpu_time(tic)
                call fld_u%boundary_transport_pressure('save',it)
                call cpu_time(toc)
                tt2=tt2+toc-tic
            ! endif

            !step 3: set hardBC
            call cpu_time(tic)
            call self%set_pressure(fld_u,time_dir,it)
            call cpu_time(toc)
            tt3=tt3+toc-tic

            !step 4: update pressure
            call cpu_time(tic)
            call self%update_pressure(fld_u,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !step 5: evolve pressure
            call cpu_time(tic)
            call self%evolve_pressure(fld_u,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !step 6: sample pressure at receivers
            call cpu_time(tic)
            call self%extract(fld_u,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !snapshot
            call fld_u%write(it)

        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add source   ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to set field    ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve field ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field',tt6/mpiworld%max_threads
        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')

    end subroutine

    subroutine adjoint(self,fld_a,fld_u,oif_record_adjseismo,oif_compute_imag,oif_compute_grad)
        class(t_propagator) :: self
        type(t_field) :: fld_a,fld_u
        logical,optional :: oif_record_adjseismo, oif_compute_imag, oif_compute_grad
        
        real,parameter :: time_dir=-1. !time direction
        logical :: if_record_adjseismo, if_compute_imag, if_compute_grad

        if_record_adjseismo =either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))
        if_compute_imag     =either(oif_compute_imag,    .false.,present(oif_compute_imag))
        if_compute_grad     =either(oif_compute_grad,    .false.,present(oif_compute_grad))
        self%if_compute_engy=self%if_compute_engy.and.(if_compute_imag.or.if_compute_grad)
        
        if(if_record_adjseismo)  call alloc(fld_a%seismo,1,self%nt)
        if(if_compute_imag)      call alloc(cb%imag,cb%mz,cb%mx,cb%my,self%nimag)
        if(if_compute_grad)      call alloc(cb%grad,cb%mz,cb%mx,cb%my,self%ngrad)
        if(self%if_compute_engy) call alloc(cb%engy,cb%mz,cb%mx,cb%my,self%nengy)

        !reinitialize absorbing boundary for incident wavefield reconstruction
        call fld_u%reinit
                    
        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.
        
        ift=1; ilt=self%nt

        ! call alloc(sf_p_save,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        
        do it=ilt,ift,int(time_dir)
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_a%check_value
                call fld_u%check_value
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then
                !backward step 5:
                call cpu_time(tic)
                call self%evolve_pressure(fld_u,time_dir,it)
                call cpu_time(toc)
                tt1=tt1+toc-tic

                !backward step 2: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_u%boundary_transport_pressure('load',it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%set_pressure(fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic

                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%update_pressure(fld_u,time_dir,it)
                call cpu_time(toc)
                tt4=tt4+toc-tic

                !backward step 1: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_pressure(fld_u,time_dir,it)
                call cpu_time(toc)
                tt5=tt5+toc-tic
            ! endif


            !--------------------------------------------------------!

            !adjoint step 6: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_pressure(fld_a,time_dir,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !adjoint step 5
            call cpu_time(tic)
            call self%evolve_pressure(fld_a,time_dir,it)
            call cpu_time(toc)
            tt7=tt7+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_pressure(fld_a,time_dir,it)
            call cpu_time(toc)
            tt8=tt8+toc-tic

            !adjoint step 5: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%set_pressure(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 1: sample v^it or s^it+0.5 at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_a,it)
                call cpu_time(toc)
                tt10=tt10+toc-tic
            endif

            !gkpa: rf%s^it+0.5 star D sf%s_dt^it+0.5
            !use sf%v^it+1 to compute sf%s_dt^it+0.5, as backward step 4
            if(if_compute_grad.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call gradient_moduli(fld_a,fld_u,it,cb%grad(:,:,:,2))
                call gradient_D     (fld_a,fld_u,it,cb%grad(:,:,:,3))
                call cpu_time(toc)
                tt11=tt11+toc-tic
            endif

            ! if(if_compute_imag.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call imaging(fld_a,fld_u,it,cb%imag)
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif

            ! !energy term of sfield
            ! if(self%if_compute_engy.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call energy(fld_u,it,cb%engy)
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
                
            !========================================================!
            
            
            ! !grho: sfield%v_dt^it \dot rfield%v^it
            ! !use sfield%s^it+0.5 to compute sfield%v_dt^it, as backward step 2
            ! if(if_compute_grad.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call gradient_density(fld_a,fld_u,it,cb%grad(:,:,:,1))
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
            
            !snapshot
            call fld_a%write(it,o_suffix='_rev')
            call fld_u%write(it,o_suffix='_rev')
            if(if_compute_imag) then
                call fld_a%write_ext(it,'imag' ,cb%imag,size(cb%imag))
            endif
            if(if_compute_grad) then
                ! call fld_a%write_ext(it,'grad_density',cb%grad(:,:,:,1),size(cb%grad(:,:,:,1)))
                call fld_a%write_ext(it,'grad_moduli' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
                ! call fld_a%write_ext(it,'grad_a_star_Du' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
            endif
            if(self%if_compute_engy) then
                call fld_a%write_ext(it,'engy',cb%engy(:,:,:,1),size(cb%engy(:,:,:,1)))
            endif

        enddo
        
        !postprocess
        ! if(if_compute_imag) call imaging_postprocess
        if(if_compute_grad) call gradient_postprocess
        
        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to evolve field        ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to load boundary       ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to set field           ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field        ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source           ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time -----------------------'
            write(*,*) 'Elapsed time to add adj source      ',tt6/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve adj field    ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj field    ',tt8/mpiworld%max_threads
            write(*,*) 'Elapsed time to set adj field       ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract&write fields',tt10/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate           ',tt11/mpiworld%max_threads

        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_rfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_rfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')
        call hud('ximage < snap_*  n1='//num2str(cb%mz)//' perc=99')
        call hud('xmovie < snap_*  n1='//num2str(cb%mz)//' n2='//num2str(cb%mx)//' clip=?e-?? loop=2 title=%g')
        
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
                f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) + wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)!*shot%src%interp_coef
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
                f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) +wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)!*shot%rcv(i)%interp_coef !no time_dir needed!
            endif

        enddo
        
    end subroutine

    subroutine set_pressure(self,f,time_dir,it)
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
            
            wl=f%wavelet(1,it)
            
            if(shot%src%comp=='pbnd') then !hard BC
                f%p(ifz:ilz,ifx:ilx,ify:ily) =                                wl!                                  *shot%src%interp_coef
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

            ! !adjsource for pressure
            ! wl=f%wavelet(i,it)
                
            if(shot%rcv(i)%comp=='pbnd') then
                f%p(ifz:ilz,ifx:ilx,ify:ily) =                               0. !wl                                  !*shot%rcv(i)%interp_coef !no time_dir needed!
            endif

        enddo
        
    end subroutine

    !forward: s^it+0.5 -> s^it+1.5 by FD of v^it+1
    !adjoint: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
    subroutine update_pressure(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ifz=f%bloom(1,it)
        if(m%is_freesurface) ifz=max(ifz,1)

        if(m%is_cubic) then
            ! call fd3d_pressure(f%p,                                      &
            !                    f%dp_dz,f%dp_dx,f%dp_dy,                  &
            !                    self%buoz,self%buox,self%buoy,self%kpa,   &
            !                    ifz,f%bloom(2,it),f%bloom(3,it),f%bloom(4,it))
        else
            call fd2d_laplacian(f%p,                            &
                                f%dp_dz,f%dp_dx,f%dpzz_dz,f%dpxx_dx,&
                                f%lapz,f%lapx,&
                                self%buoz,self%buox,self%kpa,           &
                                ifz,f%bloom(2,it),f%bloom(3,it),f%bloom(4,it),f%is_adjoint)
        endif

        if(time_dir>0.) then !in forward time
                f%p_next = 2*f%p -f%p_prev +self%dt2*self%kpa*(f%lapz + f%lapx)
        else !in reverse time
            if(.not.f%is_adjoint) then
                f%p_prev = 2*f%p -f%p_next +self%dt2*self%kpa*(f%lapz + f%lapx)
            else
                f%p_prev = 2*f%p -f%p_next +self%dt2*         (f%lapz + f%lapx) !kpa has been multiplied in fd2d_laplacian
            endif
        endif

        ! if(allocated(D)) then
        !     if(time_dir>0.) then !in forward time
        !         f%p_next(1:cb%mz,1:cb%mx,1:cb%my) = f%p_next(1:cb%mz,1:cb%mx,1:cb%my) + D*f%p(1:cb%mz,1:cb%mx,1:cb%my)
        !     else !in reverse time
        !         f%p_prev(1:cb%mz,1:cb%mx,1:cb%my) = f%p_prev(1:cb%mz,1:cb%mx,1:cb%my) + D*f%p(1:cb%mz,1:cb%mx,1:cb%my)
        !     endif
        ! endif


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

        if(time_dir>0.) then !in forward time
            f%p_prev = f%p
            f%p      = f%p_next
        else !in reverse time
            f%p_next = f%p
            f%p      = f%p_prev
        endif

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
        call dealloc(self%buoz, self%buox, self%buoy, self%kpa, self%inv_kpa)
    end subroutine


    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !Kₘ<a|Au> = Kₘ<a|∂ₜ²u-κ∇·b∇u> ≐ Kₘ -∫ aᵀκ∇·b∇u dt
    !for κ: Kₘ<a|Au> = -aᵀ∇·b∇u
    !for b: Kₘ<a|Au> = ∇(κa)·∇u
    !
    !For imaging:
    !I = ∫ a u dt =: a★u

    subroutine gradient_moduli(rf,sf,it,grad)
        type(t_field), intent(in) :: rf, sf
        real,dimension(cb%mz,cb%mx,cb%my) :: grad

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        if(m%is_cubic) then
            ! call grad3d_moduli(rf%p,sf%vz,sf%vx,sf%vy,&
            !                    grad,                  &
            !                    ifz,ilz,ifx,ilx,ify,ily)
        else
            ! call grad2d_moduli(rf%p,sf%p,sf_p_save,&
            !                    grad,            &
            !                    ifz,ilz,ifx,ilx)
            ! sf_p_save = sf%p
            
            !inexact greadient
            call grad2d_moduli(rf%p,sf%lapz,sf%lapx,&
                               grad,            &
                               ifz,ilz,ifx,ilx)
        endif
        
    end subroutine

    subroutine gradient_D(rf,sf,it,grad)
        type(t_field), intent(in) :: rf, sf
        real,dimension(cb%mz,cb%mx,cb%my) :: grad

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        if(m%is_cubic) then
            ! call grad3d_moduli(rf%p,sf%vz,sf%vx,sf%vy,&
            !                    grad,                  &
            !                    ifz,ilz,ifx,ilx,ify,ily)
        else
            call grad2d_D(rf%p,sf%p,&
                        grad,            &
                       ifz,ilz,ifx,ilx)
            
        endif
        
    end subroutine
    
    subroutine gradient_postprocess

        !scale the kernel tobe a gradient in the discretized world
        cb%grad = cb%grad*m%cell_volume*rdt
        
        ! !grho
        ! cb%grad(:,:,:,1) = cb%grad(:,:,:,1) / cb%rho(1:cb%mz,1:cb%mx,1:cb%my)
        !gkpa
        cb%grad(:,:,:,2) = cb%grad(:,:,:,2) * (-ppg%inv_kpa(1:cb%mz,1:cb%mx,1:cb%my))
                
        !preparing for cb%project_back
        cb%grad(1,:,:,:) = cb%grad(2,:,:,:)

    end subroutine

    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd2d_laplacian(p,                &
                              dp_dz,dp_dx,dpzz_dz,dpxx_dx,&
                              lapz,lapx,&
                              buoz,buox,kpa,&
                              ifz,ilz,ifx,ilx,is_adjoint)
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx) :: p
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx) :: dp_dz,dp_dx,dpzz_dz,dpxx_dx
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx) :: lapz,lapx,kpa
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx) :: buoz,buox
        logical :: is_adjoint

        real,dimension(:,:),allocatable :: pzz,pxx
        call alloc(pzz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(pxx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        nz=cb%nz
        nx=cb%nx
        
        dp_dz_=0.; dp_dx_=0.
        dpzz_dz_=0.; dpxx_dx_=0.

        if(.not. is_adjoint) then

            !flux: b∇u ~= ( bz*∂zᵇp , bx*∂ₓᵇp )
            !dir$ simd
            do ix = ifx,ilx
                do iz = ifz+1,ilz
                    dp_dz_ = c1z*(p(iz,ix) - p(iz-1,ix))

                    dp_dz(iz,ix) = cpml%b_z_half(iz)*dp_dz(iz,ix) + cpml%a_z_half(iz)*dp_dz_

                    dp_dz_ = dp_dz_/cpml%kpa_z_half(iz) + dp_dz(iz,ix)

                    pzz(iz,ix) =  buoz(iz,ix)*dp_dz_
                enddo
            enddo

            do ix = ifx+1,ilx
                do iz = ifz,ilz
                    dp_dx_ = c1x*(p(iz,ix) - p(iz,ix-1))

                    dp_dx(iz,ix) = cpml%b_x_half(ix)*dp_dx(iz,ix) + cpml%a_x_half(ix)*dp_dx_

                    dp_dx_ = dp_dx_/cpml%kpa_x_half(ix) + dp_dx(iz,ix)

                    pxx(iz,ix) = buox(iz,ix)*dp_dx_
                enddo
            enddo

            !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
            do ix = ifx,ilx
                do iz = ifz,ilz-1
                    dpzz_dz_ = c1z*(pzz(iz+1,ix) - pzz(iz,ix))

                    dpzz_dz(iz,ix) = cpml%b_z(iz)*dpzz_dz(iz,ix) + cpml%a_z(iz)*dpzz_dz_

                    dpzz_dz_ = dpzz_dz_/cpml%kpa_z(iz) + dpzz_dz(iz,ix)

                    lapz(iz,ix) = dpzz_dz_ !*kpa(iz,ix) !keep kpa later because lapz is needed for gradient computatin
                enddo
            enddo

            do ix = ifx,ilx-1
                do iz = ifz,ilz
                    dpxx_dx_ = c1x*(pxx(iz,ix+1) - pxx(iz,ix))

                    dpxx_dx(iz,ix) = cpml%b_x(ix)*dpxx_dx(iz,ix) + cpml%a_x(ix)*dpxx_dx_

                    dpxx_dx_ = dpxx_dx_/cpml%kpa_x(ix) + dpxx_dx(iz,ix)

                    lapx(iz,ix) = dpxx_dx_ !*kpa(iz,ix) !keep kpa later because lapz is needed for gradient computatin
                enddo
            enddo

        else

            !flux: b∇κa ~= ( bz*∂zᶠκ*pᵃ , bx*∂ₓᶠκ*pᵃ )
            !dir$ simd
            do ix = ifx,ilx
                do iz = ifz,ilz-1
                    dp_dz_ = c1z*(kpa(iz+1,ix)*p(iz+1,ix) - kpa(iz,ix)*p(iz,ix))

                    dp_dz(iz,ix) = cpml%b_z_half(iz)*dp_dz(iz,ix) + cpml%a_z_half(iz)*dp_dz_

                    dp_dz_ = dp_dz_/cpml%kpa_z_half(iz) + dp_dz(iz,ix)

                    pzz(iz,ix) =  buoz(iz,ix)*dp_dz_
                enddo
            enddo

            do ix = ifx,ilx-1
                do iz = ifz,ilz
                    dp_dx_ = c1x*(kpa(iz,ix+1)*p(iz,ix+1) - kpa(iz,ix)*p(iz,ix))

                    dp_dx(iz,ix) = cpml%b_x_half(ix)*dp_dx(iz,ix) + cpml%a_x_half(ix)*dp_dx_

                    dp_dx_ = dp_dx_/cpml%kpa_x_half(ix) + dp_dx(iz,ix)

                    pxx(iz,ix) = buox(iz,ix)*dp_dx_
                enddo
            enddo

            !laplacian: ∇·b∇κa ~= ∂zᵇ(bz*∂zᶠpᵃ) + ∂ₓᵇ(bx*∂ₓᶠpᵃ)
            do ix = ifx,ilx
                do iz = ifz+1,ilz
                    dpzz_dz_ = c1z*(pzz(iz,ix) - pzz(iz-1,ix))

                    dpzz_dz(iz,ix) = cpml%b_z(iz)*dpzz_dz(iz,ix) + cpml%a_z(iz)*dpzz_dz_

                    dpzz_dz_ = dpzz_dz_/cpml%kpa_z(iz) + dpzz_dz(iz,ix)

                    lapz(iz,ix) = dpzz_dz_
                enddo
            enddo

            do ix = ifx+1,ilx
                do iz = ifz,ilz
                    dpxx_dx_ = c1x*(pxx(iz,ix) - pxx(iz,ix-1))

                    dpxx_dx(iz,ix) = cpml%b_x(ix)*dpxx_dx(iz,ix) + cpml%a_x(ix)*dpxx_dx_

                    dpxx_dx_ = dpxx_dx_/cpml%kpa_x(ix) + dpxx_dx(iz,ix)

                    lapx(iz,ix) = dpxx_dx_
                enddo
            enddo

        endif

    end subroutine

    subroutine grad2d_moduli(rf_p,sf_lapz,sf_lapx,&
                             grad,          &
                             ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_p,sf_lapz,sf_lapx
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
                grad(j)=grad(j) + rf_p(i)*(sf_lapz(i)+sf_lapx(i))
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine grad2d_D(rf_p,sf_p,&
                             grad,          &
                             ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_p,sf_p
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
                grad(j)=grad(j) + rf_p(i)*sf_p(i)
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

end