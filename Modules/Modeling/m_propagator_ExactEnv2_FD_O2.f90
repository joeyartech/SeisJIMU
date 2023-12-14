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
    real,dimension(1),parameter :: coef = 1.
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z

    !local const
    real :: dt2, inv_2dt, inv_2dz, inv_2dx


    !scaling source wavelet
    real :: wavelet_scaler, virtualsrc_scaler

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain ISOtropic 2D/3D Exact Envelope propagation (the adjoint is NOT)'//s_NL// &
            '2nd-order Pressure mono tilD formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.5*Vmax/dx'//s_NL// &
            'Required model attributes: vp, rho, tilD'//s_NL// &
            'Required field components: p, p_prev, p_next'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Imaging conditions: P-Pxcorr'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: gikpa, gbuo'

        integer :: nbndlayer=max(1,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=2  !number of basic gradients
        !integer :: nimag=1 !number of basic images
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

        procedure :: forward
        procedure :: adjoint
        
        procedure :: inject_pressure
        procedure :: inject_pressure_scattering
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

    ! logical :: is_absolute_virtual

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

        if(index(self%info,'tilD')>0 .and. .not. allocated(m%tilD)) then
            call alloc(m%tilD,m%nz,m%nx,m%ny,o_init=0.)
            call warn('Constant tilD model (0) is allocated by propagator.')
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

        c1x=coef(1)/m%dx; c1y=coef(1)/m%dy; c1z=coef(1)/m%dz
        !c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz

        dt2=self%dt**2
        inv_2dt =1./2/self%dt
        inv_2dz =1./2/m%dz
        inv_2dx =1./2/m%dx
        
        wavelet_scaler=dt2/m%cell_volume

        if(setup%get_bool('IF_VIRTUAL_SRC_SCALED_BY_VOLUME_SIZE',o_default='T')) then
            virtualsrc_scaler=dt2/m%cell_volume
        else
            virtualsrc_scaler=dt2
        endif

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

        !self%invsqEref=setup%get_real('REF_ENVELOPE','RE0',o_default=num2str(1))**(-2)

        ! call alloc(self%r2,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        ! self%r2 = (cb%vp*self%dt/m%dx)**2

        ! call alloc(self%tilD_vp2,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        ! self%tilD_vp2=cb%tilD*cb%vp**2


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

        ! is_absolute_virtual=setup%get_bool('IS_ABSOLUTE_VIRTUAL','IS_ABS_VIRT',o_default='F')

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

        call alloc(f%lap, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

    end subroutine

    subroutine init_correlate(self,corr,name,purpose)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name,purpose
        
        corr%name=name
        corr%purpose=purpose

        select case (purpose)
        case ('gradient')
            call alloc(corr%rp_ddsp,        m%nz,m%nx,m%ny)
            call alloc(corr%grad_rp_grad_sp,m%nz,m%nx,m%ny)
        end select

    end subroutine

    subroutine init_abslayer(self)
        class(t_propagator) :: self

        call cpml%init

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
    !
    !!!!!!!!!!!!!!
    !! ARCHIVED !!
    !!!!!!!!!!!!!!
    !Possibility I:
    !For adjoint test:
    !forward time marching:
    ![pprev]        [pprev]   [1 ][0 1][-1 2+κ∂ᵇb∂ᶠ][1 ][1 ][pprev]
    ![pcurr] = REGBA[pcurr] = [ R][1 0][ 0 1       ][ B][ A][pcurr]
    !A:inject source into field, B:boundary load, G:propagator, with
    !∂ᶠ: forward FD, ∂ᵇ: backward FD, b:buoyancy, κ:bulk modulus, 
    !E:evolution, R:extract field at receivers.
    !While for adjoint marching (∂ᶠᵀBᵀ∂ᵇᵀKᵀ=∂ᵇBᵀ∂ᶠK):
    ![pᵃprev] [1  ][1  ][ -1        0][0 1][1  ][pᵃprev]
    ![pᵃcurr]=[ Aᵀ][ Bᵀ][2+∂ᵇbᵀ∂ᶠκ  1][1 0][ Rᵀ][pᵃcurr]
    !in reverse time:
    ![pᵃcurr]=[Aᵀ ][Bᵀ ][1 2+∂ᵇbᵀ∂ᶠκ][0 1][Rᵀ ][pᵃcurr]
    ![pᵃprev] [  1][  1][0 -1       ][1 0][  1][pᵃprev]
    !since
    ![1 2+∂ᵇbᵀ∂ᶠκ][0 1]=[0 1][ -1        0]
    ![0 -1       ][1 0] [1 0][2+∂ᵇbᵀ∂ᶠκ  1]
    !so
    ![pᵃcurr]=[Aᵀ ][Bᵀ ][0 1][ -1        0][Rᵀ ][pᵃcurr]
    ![pᵃprev] [  1][  1][1 0][2+∂ᵇbᵀ∂ᶠκ  1][  1][pᵃprev]
    !Aᵀ:extract adjoint fields, Bᵀ:boundary condition

    ![ pprev]        [ pprev]   [1   ][0 1   ][-1 2+κ∂ᵇb∂ᶠ          ][1   ][1   ][ 1       ][ pprev]
    ![ pcurr] = REGBA[ pcurr] = [ R  ][1 0   ][ 0 1                 ][ B  ][ A  ][    1    ][ pcurr]
    ![δpprev]        [δpprev]   [  1 ][   0 1][          -1 2+κ∂ᵇb∂ᶠ][  1 ][  1 ][      1  ][δpprev]
    ![δpcurr]        [δpcurr]   [   R][   1 0][           0 1       ][   B][   1][Sp Sc 0 1][δpcurr]

    ![δpᵃprev]   [1  Spᵀ][1   ][1    ][1 2+∂ᵇbᵀ∂ᶠκ           ][0 1   ][1    ][δpᵃprev]
    ![δpᵃcurr] = [ 1 Scᵀ][ Aᵀ ][ Bᵀ  ][0 -1                  ][1 0   ][ Rᵀ  ][δpᵃcurr]
    ![ pᵃprev]   [  1 0 ][  1 ][  1  ][           1 2+∂ᵇbᵀ∂ᶠκ][   0 1][  1  ][ pᵃprev]
    ![ pᵃcurr]   [    1 ][   1][   Bᵀ][           0 -1       ][   1 0][   Rᵀ][ pᵃcurr]
    !in reverse time:
    ![ pᵃcurr]   [ 1    ][1   ][Bᵀ  ][0 1   ][ -1                  ][Rᵀ  ][ pᵃcurr]
    ![ pᵃprev]   [ 0 1  ][ 1  ][ 1  ][1 0   ][2+∂ᵇbᵀ∂ᶠκ            ][ 1  ][ pᵃprev]
    ![δpᵃcurr] = [Scᵀ 1 ][  Aᵀ][  Bᵀ][   0 1][          1 2+∂ᵇbᵀ∂ᶠκ][  Rᵀ][δpᵃcurr]
    ![δpᵃprev]   [Spᵀ  1][   1][   1][   1 0][          0 -1       ][   1][δpᵃprev]
    !
    !
    !Possibility II:
    ![pprev]      [pprev] [1  ] [0 1 0] [ 1           ] [1  ] [1  ] [pprev]
    !|pcurr|=REGBA|pcurr|=| R | |0 0 1| |      1      | | B | | A | |pcurr|
    ![pnext]      [pnext] [  1] [1 0 0] [-1 2+κ∂ᵇb∂ᶠ 0] [  1] [  1] [pnext]
    ![pᵃprev]           [pᵃprev] [1   ] [1   ] [1    -1        ] [0 0 1] [1   ] [pᵃprev]!
    !|pᵃcurr|=BᵀAᵀGᵀEᵀRᵀ|pᵃcurr|=| Bᵀ | | Aᵀ | |   1  2+∂ᵇbᵀ∂ᶠκ| |1 0 0| | Rᵀ | |pᵃcurr|
    ![pᵃnext]           [pᵃnext] [   1] [   1] [      0        ] [0 1 0] [   1] [pᵃnext]
    !in reverse time:
    ![pᵃnext] [1   ] [1   ] [             ] [0 1 0] [1   ] [pᵃnext]
    !|pᵃcurr|=| Bᵀ | | Aᵀ | |2+∂ᵇbᵀ∂ᶠκ 1 0| |0 0 1| | Rᵀ | |pᵃcurr|
    ![pᵃprev] [   1] [   1] [ -1       1  ] [1 0 0] [   1] [pᵃprev]
    !                       [             ] [pᵃcurr]
    !                      ~|2+∂ᵇbᵀ∂ᶠκ 1 0| |pᵃprev|
    !                       [ -1       0 1] [pᵃnext]
    !
    !Possibility III:
    !For adjoint test:
    !forward time marching:
    ![pprev]      [pprev] [1  ] [0 1 0] [ 1           ] [1  ] [1  ] [pprev]
    !|pcurr|=REGBA|pcurr|=| R | |0 0 1| |      1      | | B | | A | |pcurr|
    ![pnext]      [pnext] [  1] [1 0 0] [-1 2+κ∂ᵇb∂ᶠ 0] [  1] [  1] [pnext]
    !A:inject source into field, B:boundary load, G:propagator, with
    !∂ᶠ: forward FD, ∂ᵇ: backward FD, b:buoyancy, κ:bulk modulus, 
    !E:evolution, R:extract field at receivers.
    !While for adjoint marching (∂ᶠᵀBᵀ∂ᵇᵀKᵀ=∂ᵇBᵀ∂ᶠK):
    ![pᵃprev]           [pᵃprev] [1   ] [1   ] [1    -1        ] [0 0 1] [1   ] [pᵃprev]
    !|pᵃcurr|=BᵀAᵀGᵀEᵀRᵀ|pᵃcurr|=| Bᵀ | | Aᵀ | |   1  2+∂ᵇbᵀ∂ᶠκ| |1 0 0| | Rᵀ | |pᵃcurr|
    ![pᵃnext]           [pᵃnext] [   1] [   1] [      0        ] [0 1 0] [   1] [pᵃnext]
    !In reverse time:
    ![pᵃnext] [1   ] [1   ] [             ] [0 1 0] [1   ] [pᵃnext]
    !|pᵃcurr|=| Bᵀ | | Aᵀ | |2+∂ᵇbᵀ∂ᶠκ 1 0| |0 0 1| | Rᵀ | |pᵃcurr|
    ![pᵃprev] [   1] [   1] [-1        0 1] [1 0 0] [   1] [pᵃprev]
    !since
    ![  0       0 0][0 1 0]   [0   0       0]   [0 1 0][1  -1        0]
    !|2+∂ᵇbᵀ∂ᶠκ 1 0||0 0 1| = |0 2+∂ᵇbᵀ∂ᶠκ 1| = |0 0 1||0   0        0|
    ![ -1       0 1][1 0 0]   [1  -1       0]   [1 0 0][0  2+∂ᵇbᵀ∂ᶠκ 1]
    ![pᵃnext] [1   ] [1   ] [0 1 0] [1  -1        0] [1   ] [pᵃnext]
    !|pᵃcurr|=| Bᵀ | | Aᵀ | |0 0 1| |0   0        0| | Rᵀ | |pᵃcurr|
    ![pᵃprev] [   1] [   1] [1 0 0] [0  2+∂ᵇbᵀ∂ᶠκ 1] [   1] [pᵃprev]
    !Aᵀ:extract adjoint fields, Bᵀ:boundary condition
    
    ! [pprev] [0 1 0][ 1           ][pprev]
    ! |pcurr|=|0 0 1||      1      ||pcurr|
    ! [pnext] [1 0 0][-1 2+κ∂ᵇb∂ᶠ 0][pnext]
    ! [pᵃprev] [1    -1        ][0 0 1][pᵃnext]
    ! |pᵃcurr|=|   1  2+∂ᵇbᵀ∂ᶠκ||1 0 0||pᵃprev|
    ! [pᵃnext] [      0        ][0 1 0][pᵃcurr]

    ! [ pprev] [0 1 0      ][           ][ 1                         ][ pprev]
    ! | pcurr|=|0 0 1      |[           ]|      1                    || pcurr|
    ! [ pnext] [1 0 0      ][           ][-1 2+κ∂ᵇb∂ᶠ 0              ][ pnext]
    ! [δpprev] [      0 1 0][0  0   0   ][               1           ][δpprev]
    ! |δpcurr|=|      0 0 1||1 -2+D 1   ||                    1      ||δpcurr|
    ! [δpnext] [      1 0 0][0  0   0   ][              -1 2+κ∂ᵇb∂ᶠ 0][δpnext]
    
    ! [ pᵃprev] [1    -1                        ][   0  1    0][0 0 1      ][ pᵃprev]
    ! | pᵃcurr|=|   1  2+∂ᵇbᵀ∂ᶠκ                |[   0 -2+Dᵀ 0]|1 0 0      || pᵃcurr|
    ! [ pᵃnext] [      0                        ][   0  1    0][0 1 0      ][ pᵃnext]
    ! [δpᵃprev] [                1    -1        ][            ][      0 0 1][δpᵃprev]
    ! |δpᵃcurr|=|                   1  2+∂ᵇbᵀ∂ᶠκ||            ||      1 0 0||δpᵃcurr|
    ! [δpᵃnext] [                      0        ][            ][      0 1 0][δpᵃnext]

    

    subroutine forward(self,fld_u)
        class(t_propagator) :: self
        type(t_field) :: fld_u

        real,parameter :: time_dir=1. !time direction

        !seismo
        call alloc(fld_u%seismo, shot%nrcv,self%nt)
            
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
            tt5=tt5+toc-tic

            !step 6: sample p^it+1 at receivers
            call cpu_time(tic)
            call self%extract(fld_u,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            ! if(present(o_E0_star_E0)) then
            !     call cpu_time(tic)
            !     call auto_correlate(fld_E0,o_E0_star_E0,it)
            !     call cpu_time(toc)
            !     tt7=tt7+toc-tic
            ! endif

            !snapshot
            call fld_u%write(it)

        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add source   ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary',tt2/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set field    ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve field ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field',tt6/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to correlate    ',tt7/mpiworld%max_threads
        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')


        ! !postprocess
        ! if(present(o_E0_star_E0)) then
        !     !scale by m%cell_volume*rdt tobe an energy distribution in the discretized world
        !     call o_E0_star_E0%scale(m%cell_volume*rdt)
        ! endif

    end subroutine


    subroutine adjoint(self,fld_q,fld_p, fld_v,fld_u, q_star_v,p_star_u)
        class(t_propagator) :: self
        type(t_field) :: fld_q,fld_p,fld_v,fld_u
        type(t_correlate) :: q_star_v,  p_star_u

        real,parameter :: time_dir=-1. !time direction

        !seismo
        call alloc(fld_q%seismo, shot%nrcv,self%nt)
        call alloc(fld_p%seismo, shot%nrcv,self%nt)

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

                !set hardBC
                !call self%set_pressure(fld_dE,         it)
                ! fld_dE%p=surD(:,:,:,it)*fld_E0%p

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

            !gradient: rf%p^it star sf%p^it
            ! if(mod(it,irdt)==0) then
                call cpu_time(tic)
                call cross_correlate(fld_q,fld_v,q_star_v,it)
                call cross_correlate(fld_p,fld_u,p_star_u,it)
                call cpu_time(toc)
                tt10=tt10+toc-tic
            ! endif

            !adjoint step 5
            ! this step is moved to update_pressure for easier management
            call cpu_time(tic)
            call self%evolve_pressure(fld_q,time_dir,it)
            call self%evolve_pressure(fld_p,time_dir,it)
            call cpu_time(toc)
            tt11=tt11+toc-tic

            ! !adjoint step 5: inject to s^it+1.5 at receivers
            ! call cpu_time(tic)
            ! call self%set_pressure(fld_a,time_dir,it)
            ! fld_dF%p = surD(:,:,:,it)*fld_F0%p
            ! call cpu_time(toc)
            ! tt9=tt9+toc-tic

            !adjoint step 1: sample p^it at source position
            ! if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_q,it)
                call self%extract(fld_p,it)
                call cpu_time(toc)
                tt12=tt12+toc-tic
            ! endif

            ! !gkpa: rf%s^it+0.5 star tilD sf%s_dt^it+0.5
            ! !if(if_compute_grad.and.mod(it,irdt)==0) then
            ! if(mod(it,irdt)==0 .and. present(o_F2_star_E0)) then
            !     call cpu_time(tic)
            !     call gradient_tilD(fld_F2,fld_E0,o_F2_star_E0,self%dt,it)
            !     call cpu_time(toc)
            !     tt11=tt11+toc-tic
            ! endif


            !--------------------------------------------------------!
            
            !snapshot
            call fld_v%write(it,o_suffix='_adj')
            call fld_u%write(it,o_suffix='_adj')
            call fld_q%write(it,o_suffix='_adj')
            call fld_p%write(it,o_suffix='_adj')
            
            call q_star_v%write(it,o_suffix='_adj')
            call p_star_u%write(it,o_suffix='_adj')
            
        enddo


        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to evolve field        ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to load boundary       ',tt2/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set field           ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field        ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source           ',tt6/mpiworld%max_threads        
            write(*,*) 'Elapsed time -----------------------'
            write(*,*) 'Elapsed time to add adj & virtual src',tt8/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj field    ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve adj field    ',tt11/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set adj field       ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract fields      ',tt12/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate           ',tt10/mpiworld%max_threads

        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_rfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_rfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')
        call hud('ximage < snap_*  n1='//num2str(cb%mz)//' perc=99')
        call hud('xmovie < snap_*  n1='//num2str(cb%mz)//' n2='//num2str(cb%mx)//' clip=?e-?? loop=2 title=%g')


        ! if(if_corr) then
            call q_star_v%scale(m%cell_volume*rdt)
            call p_star_u%scale(m%cell_volume*rdt)
            
            select case (setup%get_str('GRADIENT_TERMS',o_default='pu+qv'))
            case ('pu+qv')
                call correlate_stack(q_star_v%rp_ddsp        +p_star_u%rp_ddsp        , correlate_gradient(:,:,:,1)) !gikpa
                call correlate_stack(q_star_v%grad_rp_grad_sp+p_star_u%grad_rp_grad_sp, correlate_gradient(:,:,:,2)) !gbuo

            case ('pu')
                call correlate_stack(p_star_u%rp_ddsp        , correlate_gradient(:,:,:,1)) !gikpa
                call correlate_stack(p_star_u%grad_rp_grad_sp, correlate_gradient(:,:,:,2)) !gbuo

            case ('qv')
                call correlate_stack(q_star_v%rp_ddsp        , correlate_gradient(:,:,:,1)) !gikpa
                call correlate_stack(q_star_v%grad_rp_grad_sp, correlate_gradient(:,:,:,2)) !gbuo
            end select
        ! endif

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

    subroutine set_pressure(self,f,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ! if(.not. f%is_adjoint) then
!this loop takes more time when nthreads>1.
!e.g. 0.1s vs 0.046s from Elapased time to rm source (tt5)
!tobe fixed..

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
            
            wl=f%wavelet(i,it)!*wavelet_scaler

            if(if_hicks) then 
                ! f%p(ifz:ilz,ifx:ilx,ify:ily) = wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef
            else
                f%p(ifz:ilz,ifx:ilx,ify:ily)                = wl!*self%kpa(iz,ix,iy)
            endif

        enddo

    end subroutine

    subroutine inject_pressure_scattering(self,fld_d,fld,D,it)
        class(t_propagator) :: self
        type(t_field) :: fld_d, fld
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: D

        ifz=fld%bloom(1,it)
        ilz=fld%bloom(2,it)
        ifx=fld%bloom(3,it)
        ilx=fld%bloom(4,it)

        !fld_d%p = fld_d%p +1.0*dt2*self%kpa*D*fld%p /m%cell_volume !I don't understand why should divide by the cell_volume..
        fld_d%p(2:m%nz-1,2:m%nz-1,1) = fld_d%p(2:m%nz-1,2:m%nz-1,1) &
            ! +1.0*dt2*self%kpa(2:m%nz-1,2:m%nz-1,1)*D(2:m%nz-1,2:m%nz-1,1)*fld%p(2:m%nz-1,2:m%nz-1,1)
            !+1.0*dt2*self%kpa(2:m%nz-1,2:m%nz-1,1)*D(2:m%nz-1,2:m%nz-1,1)*fld%p(2:m%nz-1,2:m%nz-1,1) /m%cell_volume !I don't understand why should divide by the cell_volume..
            +1.0*self%kpa(2:m%nz-1,2:m%nz-1,1)*D(2:m%nz-1,2:m%nz-1,1)*fld%p(2:m%nz-1,2:m%nz-1,1) *virtualsrc_scaler !w cell_volume by default
            
        ! call flattend_inject_scattering(fld_d%p,fld%p,self%kpa,D,ifz,ilz,ifx,ilx)

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

    subroutine cross_correlate(rf,sE,corr,it)
        type(t_field), intent(in) :: rf, sE
        type(t_correlate) :: corr

        ! real,dimension(:,:,:),allocatable :: tmp_sf_p

        !nonzero only when sf touches rf
        ifz=max(sE%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sE%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sE%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sE%bloom(4,it),rf%bloom(4,it),cb%mx)
        ! ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ! ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        ! if(m%is_cubic) then
        !     ! call grad3d_moduli(rf%p,sf%vz,sf%vx,sf%vy,&
        !     !                    grad,                  &
        !     !                    ifz,ilz,ifx,ilx,ify,ily)
        ! else
        !     call grad2d_tilD((rf%p_next-rf%p_prev)/2/dt,(sf%p_next-sf%p_prev)/2/dt,&
        !                       grad%sp_rp,dt,     &
        !                       ifz,ilz,ifx,ilx)
            
        ! endif

        if(m%is_cubic) then
        
        else if(corr%purpose=='gradient') then

            !for gikpa
            corr%rp_ddsp = corr%rp_ddsp + &
                rf%p(1:m%nz,1:m%nx,1:m%ny) * ( sE%p_next(1:m%nz,1:m%nx,1:m%ny)-2*sE%p(1:m%nz,1:m%nx,1:m%ny)+sE%p_prev(1:m%nz,1:m%nx,1:m%ny) )/dt2

            !for gbuo
            call grad2d(rf%p,sE%p,corr%grad_rp_grad_sp,   ifz,ilz,ifx,ilx)

        endif

    end subroutine

    subroutine cross_correlate_Eph(rf,sE,sph,corr,it)
        type(t_field), intent(in) :: rf, sE, sph
        type(t_correlate) :: corr

        ! real,dimension(:,:,:),allocatable :: tmp_sf_p

        !nonzero only when sf touches rf
        ifz=max(sE%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sE%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sE%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sE%bloom(4,it),rf%bloom(4,it),cb%mx)
        ! ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ! ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        ! if(m%is_cubic) then
        !     ! call grad3d_moduli(rf%p,sf%vz,sf%vx,sf%vy,&
        !     !                    grad,                  &
        !     !                    ifz,ilz,ifx,ilx,ify,ily)
        ! else
        !     call grad2d_tilD((rf%p_next-rf%p_prev)/2/dt,(sf%p_next-sf%p_prev)/2/dt,&
        !                       grad%sp_rp,dt,     &
        !                       ifz,ilz,ifx,ilx)
            
        ! endif

        if(m%is_cubic) then
        
        else if(corr%purpose=='gradient') then

            corr%rp_ddsp = corr%rp_ddsp + &
                rf%p(1:m%nz,1:m%nx,1:m%ny) * sE%p(1:m%nz,1:m%nx,1:m%ny) *((asin(sin(sph%p_next(1:m%nz,1:m%nx,1:m%ny)-sph%p_prev(1:m%nz,1:m%nx,1:m%ny))))*inv_2dt)**2

            !for gbuo [wait]
            !call grad2d_ph(rf%p,sph%p,corr%grad_rp_grad_sp,   ifz,ilz,ifx,ilx)

        endif

    end subroutine


    !========= Finite-Difference on flattened arrays ==================
    
    ! subroutine flattend_inject_scattering(dp,p,kpa,D,ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: dp,p,kpa,D
        
    !     !$omp parallel default (shared)&
    !     !$omp private(i)
    !     !$omp do schedule(dynamic)
    !     do ix = ifx,ilx
    !         !dir$ simd
    !         do iz = ifz,ilz

    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

    !             dp(i) = dp(i) +1.*dt2*kpa(i)*D(i)*p(i) /m%cell_volume

    !         enddo
    !     enddo
    !     !$omp end do
    !     !$omp end parallel
        
    ! end subroutine

    subroutine fd2d_grad_ph_sq(ph, grad_ph_sq, ifz,ilz,ifx,ilx)
        real,dimension(*) :: ph, grad_ph_sq

        nz=cb%nz
        nx=cb%nx
        
        !$omp parallel default (shared)&
        !$omp private(i)
        !$omp do schedule(dynamic)
        do ix = ifx+1,ilx-1
            !dir$ simd
            do iz = ifz+1,ilz-1

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm1_ix=i-1       !iz-1,ix
                iz_ix  =i         !iz,ix
                izp1_ix=i+1       !iz+1,ix

                iz_ixm1=i    -nz  !iz,ix-1
                iz_ixp1=i    +nz  !iz,ix+1

                grad_ph_sq(iz_ix) = ( asin(sin(ph(izp1_ix)-ph(izm1_ix))) *inv_2dx)**2 &
                                   +( asin(sin(ph(iz_ixp1)-ph(iz_ixm1))) *inv_2dz)**2

            enddo
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine

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
        !$omp         izm1_ix,iz_ix,iz_ixm1,&
        !$omp         dp_dz_,dp_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx+1,ilx
            !dir$ simd
            do iz = ifz+1,ilz

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm1_ix=i-1       !iz-1,ix
                iz_ix  =i         !iz,ix
                iz_ixm1=i    -nz  !iz,ix-1

                dp_dz_ = c1z*(p(iz_ix) - p(izm1_ix))
                dp_dx_ = c1x*(p(iz_ix) - p(iz_ixm1))

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
        !$omp         iz_ix,izp1_ix,iz_ixp1,&
        !$omp         dpzz_dz_,dpxx_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx,ilx-1
            !dir$ simd
            do iz = ifz,ilz-1

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                iz_ix  =i         !iz,ix
                izp1_ix=i+1       !iz+1,ix
                iz_ixp1=i    +nz  !iz,ix+1

                dpzz_dz_ = c1z*(pzz(izp1_ix) - pzz(iz_ix))
                dpxx_dx_ = c1x*(pxx(iz_ixp1) - pxx(iz_ix))

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


    subroutine grad2d(rp,sE,grad_rp_grad_sp,&
                        ifz,ilz,ifx,ilx)
        real,dimension(*) :: rp, sE
        real,dimension(*) :: grad_rp_grad_sp

        real,dimension(:),allocatable :: pzz, pxx
        call alloc(pzz,cb%n)
        call alloc(pxx,cb%n)

        nz=cb%nz
        nx=cb%nx
        
        
        !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         iz_ixm1,izm1_ix,iz_ix,izp1_ix,iz_ixp1)
        !$omp do schedule(dynamic)
        do ix = ifx,ilx-1
            !dir$ simd
            do iz = ifz,ilz-1

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers

                iz_ixm1=i    -nz  !iz,ix-1
                izm1_ix=i-1       !iz-1,ix
                iz_ix  =i
                izp1_ix=i+1       !iz+1,ix
                iz_ixp1=i    +nz  !iz,ix+1

                
                grad_rp_grad_sp(j) = grad_rp_grad_sp(j) &
                    + (rp(izp1_ix)-rp(izm1_ix))*(sE(izp1_ix)-sE(izm1_ix))*inv_2dz*inv_2dz &
                    + (rp(iz_ixp1)-rp(iz_ixm1))*(sE(iz_ixp1)-sE(iz_ixm1))*inv_2dx*inv_2dx

            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

end
