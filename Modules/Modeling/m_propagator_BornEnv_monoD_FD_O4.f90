module m_propagator
use m_System
use m_hicks, only : hicks_r
use m_resampler
use m_model
use m_shot
use m_computebox
use m_field
use m_correlation
use m_cpml

    private

    !FD coef
    real,dimension(2),parameter :: coef = [9./8., -1./24.]
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z

    !local const
    real :: inv_dt
    real :: dt2, inv_2dz, inv_2dx, inv_2dt, inv_dz2,inv_dx2

    !scaling source wavelet
    real :: wavelet_scaler

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain ISOtropic 2D/3D Born-based Envelope propagation'//s_NL// &
            '2nd-order Pressure mono tilD formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x²,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.5*Vmax/dx'//s_NL// &
            'Required model attributes: vp, rho, tilD'//s_NL// &
            'Required field components: p, p_prev, p_next'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Imaging conditions: P-Pxcorr'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: gvp2 gtilD'
            !'Basic gradients: grho gkpa gtilD'

        integer :: nbndlayer=max(1,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=1
        !integer :: ngrad=3 !number of basic gradients
        !integer :: nimag=1 !number of basic images
        !integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:,:),allocatable :: buoz, buox, buoy, kpa, inv_kpa, r2, tilD_vp2

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
        procedure :: init_abslayer

        procedure :: forward
        procedure :: forward_scattering
        procedure :: adjoint
        procedure :: adjoint_tilD
        procedure :: adjoint_vp2
        
        procedure :: inject_pressure
        procedure :: inject_pressure_scattering
        procedure :: inject_pressure_adjoint_scattering
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
        c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz

        rgc1z = 2./3. /m%dz; rgc1x = 2./3. /m%dx
        rgc2z =-1./12./m%dz; rgc2x =-1./12./m%dx

        inv_dz2 =1./m%dz/m%dz
        inv_dx2 =1./m%dx/m%dx

        inv_dt=1./self%dt
        inv_2dt=1./2/self%dt
        dt2=self%dt**2
        
        wavelet_scaler=dt2/m%cell_volume

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

        !self%invsqEref=setup%get_real('REF_ENVELOPE','RE0',o_default=num2str(1))**(-2)

        call alloc(self%r2,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        self%r2 = (cb%vp*self%dt/m%dx)**2

        call alloc(self%tilD_vp2,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        self%tilD_vp2=cb%tilD*cb%vp**2


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
    !Adjoint:  Aᵀa = d
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

    

    subroutine forward(self,fld_E0)
        class(t_propagator) :: self
        type(t_field) :: fld_E0

        real,parameter :: time_dir=1. !time direction

        !seismo
        call alloc(fld_E0%seismo, shot%nrcv,self%nt)
            
        tt1=0.; tt2=0.; tt3=0.; tt4=0.; tt5=0.; tt6=0.

        ift=1; ilt=self%nt

        do it=ift,ilt
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_E0%check_value
            endif

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add pressure
            call cpu_time(tic)
            call self%inject_pressure(fld_E0,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2: save p^it+1 in boundary layers
            ! if(fld_E0%if_will_reconstruct) then
                call cpu_time(tic)
                call fld_E0%boundary_transport_pressure('save',it)
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
            call self%update_pressure(fld_E0,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !step 5: evolve pressure
            call cpu_time(tic)
            call self%evolve_pressure(fld_E0,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !step 6: sample pressure at receivers
            call cpu_time(tic)
            call self%extract(fld_E0,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !snapshot
            call fld_E0%write(it)

        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add source   ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary',tt2/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set field    ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve field ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field',tt6/mpiworld%max_threads
        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')

    end subroutine

    subroutine forward_scattering(self,fld_dE,fld_E0)
        class(t_propagator) :: self
        type(t_field) :: fld_dE, fld_E0

        real,parameter :: time_dir=1. !time direction

        !seismo
        call alloc(fld_dE%seismo,shot%nrcv,self%nt)
        call alloc(fld_E0%seismo, shot%nrcv,self%nt)
            
        tt1=0.; tt2=0.; tt3=0.; tt4=0.; tt5=0.; tt6=0.

        ift=1; ilt=self%nt

call warn('no dE modeling to save time while testing..')

        do it=ift,ilt
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_E0%check_value
                call fld_dE%check_value
            endif

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add pressure
            call cpu_time(tic)
            call self%inject_pressure(fld_E0,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2: save p^it+1 in boundary layers
            ! if(fld_E0%if_will_reconstruct) then
                call cpu_time(tic)
                call fld_E0%boundary_transport_pressure('save',it)
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
            call self%update_pressure(fld_E0,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic


!            !do forward time stepping (step# conforms with backward & adjoint time stepping)
!            !step 1: add pressure
!            call cpu_time(tic)
!            call self%inject_pressure_scattering(fld_dE,fld_E0,time_dir,it)
!            call cpu_time(toc)
!            tt1=tt1+toc-tic
!
!             !step 2: save p^it+1 in boundary layers
!             ! if(fld_E0%if_will_reconstruct) then
!                 call cpu_time(tic)
!                 call fld_dE%boundary_transport_pressure('save',it)
!                 call cpu_time(toc)
!                 tt2=tt2+toc-tic
!             ! endif
!
!             !step 4: update pressure
!             call cpu_time(tic)
!             call self%update_pressure(fld_dE,time_dir,it)
!             call cpu_time(toc)
!             tt4=tt4+toc-tic



            !step 5: evolve pressure
            call cpu_time(tic)
            call self%evolve_pressure(fld_E0,time_dir,it)
!             call self%evolve_pressure(fld_dE,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !step 6: sample pressure at receivers
            call cpu_time(tic)
            call self%extract(fld_E0,it)
!             call self%extract(fld_dE,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !snapshot
            call fld_E0%write(it)
!             call fld_dE%write(it)

        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add source   ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary',tt2/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set field    ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve field ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field',tt6/mpiworld%max_threads
        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')

    end subroutine

    subroutine adjoint(self,fld_a,fld_u,oif_record_adjseismo,o_a_star_u)
        class(t_propagator) :: self
        type(t_field) :: fld_a,fld_u
        logical,optional :: oif_record_adjseismo
        type(t_correlation),optional :: o_a_star_u
        
        real,parameter :: time_dir=-1. !time direction
        logical :: if_record_adjseismo

        if_record_adjseismo =either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))
        if(if_record_adjseismo)  call alloc(fld_a%seismo,1,self%nt)

        ! call alloc(cb%grad,cb%mz,cb%mx,cb%my,3)
        
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

                !backward step 2: retrieve p^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_u%boundary_transport_pressure('load',it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                ! !backward step 3: s^it+1.5 -> s^it+0.5
                ! call cpu_time(tic)
                ! call self%set_pressure(fld_u,time_dir,it)
                ! call cpu_time(toc)
                ! tt3=tt3+toc-tic

                !backward step 4: s^it+1.5 -> s^it+0.5
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

            !gkpa: rf%s^it+0.5 star tilD sf%s_dt^it+0.5
            !if(if_compute_grad.and.mod(it,irdt)==0) then
            if(mod(it,irdt)==0 .and. present(o_a_star_u)) then
                call cpu_time(tic)
                call gradient_vp2_nab_rp_nab_sp(fld_a,fld_u,o_a_star_u,it)
                call cpu_time(toc)
                tt11=tt11+toc-tic
            endif

            !--------------------------------------------------------!

            !adjoint step 6: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_pressure(fld_a,time_dir,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5
            call cpu_time(tic)
            call self%update_pressure(fld_a,time_dir,it)
            call cpu_time(toc)
            tt8=tt8+toc-tic

            !adjoint step 5
            ! this step is moved to update_pressure for easier management
            call cpu_time(tic)
            call self%evolve_pressure(fld_a,time_dir,it)
            call cpu_time(toc)
            tt7=tt7+toc-tic

            ! !adjoint step 5: inject to s^it+1.5 at receivers
            ! call cpu_time(tic)
            ! call self%set_pressure(fld_a,time_dir,it)
            ! call cpu_time(toc)
            ! tt9=tt9+toc-tic

            !adjoint step 1: sample v^it or s^it+0.5 at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_a,it)
                call cpu_time(toc)
                tt10=tt10+toc-tic
            endif
            
            !snapshot
            call fld_a%write(it,o_suffix='_rev')
            call fld_u%write(it,o_suffix='_rev')
            
            if(present(o_a_star_u)) then
                call fld_u%write_ext(it,'corr_a_star_u_nab2',o_a_star_u%nab_rp_nab_sp,m%n)
            endif
            
        enddo
        
        if(present(o_a_star_u)) then
            !call gradient_postprocess
            ! !scale by m%cell_volume*rdt tobe a gradient in the discretized world
            ! cb%corr = cb%corr*m%cell_volume*rdt
            ! !preparing for cb%project_back
            ! cb%corr(1,:,:,:) = cb%corr(2,:,:,:)
            o_a_star_u%nab_rp_nab_sp = o_a_star_u%nab_rp_nab_sp*m%cell_volume*rdt
            o_a_star_u%nab_rp_nab_sp(1,:,:,:) = o_a_star_u%nab_rp_nab_sp(2,:,:,:)
        endif

        
        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to evolve field        ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to load boundary       ',tt2/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set field           ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field        ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source           ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time -----------------------'
            write(*,*) 'Elapsed time to add adj source      ',tt6/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve adj field    ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj field    ',tt8/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set adj field       ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract&write fields',tt10/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate           ',tt11/mpiworld%max_threads

        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_rfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_rfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')
        call hud('ximage < snap_*  n1='//num2str(cb%mz)//' perc=99')
        call hud('xmovie < snap_*  n1='//num2str(cb%mz)//' n2='//num2str(cb%mx)//' clip=?e-?? loop=2 title=%g')
        
    end subroutine


    subroutine adjoint_tilD(self,fld_F2,fld_E0,oif_record_adjseismo,o_F2_star_E0)
        class(t_propagator) :: self
        type(t_field) :: fld_F2,fld_E0
        logical,optional :: oif_record_adjseismo
        type(t_correlation),optional :: o_F2_star_E0
        
        real,parameter :: time_dir=-1. !time direction
        logical :: if_record_adjseismo

        if_record_adjseismo =either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))
        if(if_record_adjseismo)  call alloc(fld_F2%seismo,1,self%nt)

        ! call alloc(cb%grad,cb%mz,cb%mx,cb%my,3)
        
        !reinitialize absorbing boundary for incident wavefield reconstruction
        call fld_E0%reinit
                    
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
                call fld_F2%check_value
                call fld_E0%check_value
            endif            


! !gkpa: rf%s^it+0.5 star tilD sf%s_dt^it+0.5
! !if(if_compute_grad.and.mod(it,irdt)==0) then
! if(mod(it,irdt)==0 .and. present(o_F2_star_E0)) then
!     call cpu_time(tic)
!     call gradient_tilD(fld_F2,fld_E0,o_F2_star_E0,self%dt,it)! 1st cond   1.00000005E-03   1.21089315      0.921757936       289.135193      -289.541931      -1.00140679     F
!     call cpu_time(toc)
!     tt11=tt11+toc-tic
! endif

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then
                !backward step 5:
                call cpu_time(tic)
                call self%evolve_pressure(fld_E0,time_dir,it)
                call cpu_time(toc)
                tt1=tt1+toc-tic

                !backward step 2: retrieve p^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_E0%boundary_transport_pressure('load',it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                ! !backward step 3: s^it+1.5 -> s^it+0.5
                ! call cpu_time(tic)
                ! call self%set_pressure(fld_E0,time_dir,it)
                ! call cpu_time(toc)
                ! tt3=tt3+toc-tic

                !backward step 4: s^it+1.5 -> s^it+0.5
                call cpu_time(tic)
                call self%update_pressure(fld_E0,time_dir,it)
                call cpu_time(toc)
                tt4=tt4+toc-tic

                !backward step 1: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_pressure(fld_E0,time_dir,it)
                call cpu_time(toc)
                tt5=tt5+toc-tic
            ! endif

            !gkpa: rf%s^it+0.5 star tilD sf%s_dt^it+0.5
            !if(if_compute_grad.and.mod(it,irdt)==0) then
            if(mod(it,irdt)==0 .and. present(o_F2_star_E0)) then
                call cpu_time(tic)
                call gradient_tilD(fld_F2,fld_E0,o_F2_star_E0,self%dt,it) ! 1st cond   1.00000005E-03   1.28531706      0.913123488       372.193542      -288.561340     -0.775299132     F
                call cpu_time(toc)
                tt11=tt11+toc-tic
            endif

            !--------------------------------------------------------!

            !adjoint step 6: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_pressure(fld_F2,time_dir,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5
            call cpu_time(tic)
            call self%update_pressure(fld_F2,time_dir,it)
            call cpu_time(toc)
            tt8=tt8+toc-tic

! !gkpa: rf%s^it+0.5 star tilD sf%s_dt^it+0.5
! !if(if_compute_grad.and.mod(it,irdt)==0) then
! if(mod(it,irdt)==0 .and. present(o_F2_star_E0)) then
!     call cpu_time(tic)
!     call gradient_tilD(fld_F2,fld_E0,o_F2_star_E0,self%dt,it) !1st cond   1.00000005E-03  0.624549747      0.459028095       165.521637      -288.936798      -1.74561346     F
!     call cpu_time(toc)
!     tt11=tt11+toc-tic
! endif

            !adjoint step 5
            ! this step is moved to update_pressure for easier management
            call cpu_time(tic)
            call self%evolve_pressure(fld_F2,time_dir,it)
            call cpu_time(toc)
            tt7=tt7+toc-tic

            ! !adjoint step 5: inject to s^it+1.5 at receivers
            ! call cpu_time(tic)
            ! call self%set_pressure(fld_a,time_dir,it)
            ! call cpu_time(toc)
            ! tt9=tt9+toc-tic

            !adjoint step 1: sample v^it or s^it+0.5 at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_F2,it)
                call cpu_time(toc)
                tt10=tt10+toc-tic
            endif

            ! !gkpa: rf%s^it+0.5 star tilD sf%s_dt^it+0.5
            ! !if(if_compute_grad.and.mod(it,irdt)==0) then
            ! if(mod(it,irdt)==0 .and. present(o_F2_star_E0)) then
            !     call cpu_time(tic)
            !     call gradient_tilD(fld_F2,fld_E0,o_F2_star_E0,self%dt,it)
            !     call cpu_time(toc)
            !     tt11=tt11+toc-tic
            ! endif
            
            !snapshot
            call fld_F2%write(it,o_suffix='_rev')
            call fld_E0%write(it,o_suffix='_rev')
            
            if(present(o_F2_star_E0)) then
                call fld_E0%write_ext(it,'corr_F2_star_E0_dt2', o_F2_star_E0%drp_dt_dsp_dt,m%n)
                call fld_E0%write_ext(it,'corr_F2_star_E0_nab2',o_F2_star_E0%nab_rp_nab_sp,m%n)
            endif
            
        enddo
        
        if(present(o_F2_star_E0)) then
            !call gradient_postprocess
            ! !scale by m%cell_volume*rdt tobe a gradient in the discretized world
            ! cb%corr = cb%corr*m%cell_volume*rdt
            ! !preparing for cb%project_back
            ! cb%corr(1,:,:,:) = cb%corr(2,:,:,:)
            o_F2_star_E0%drp_dt_dsp_dt = o_F2_star_E0%drp_dt_dsp_dt*m%cell_volume*rdt !*self%invsqEref
            o_F2_star_E0%drp_dt_dsp_dt(1,:,:,:) = o_F2_star_E0%drp_dt_dsp_dt(2,:,:,:)

            o_F2_star_E0%nab_rp_nab_sp = o_F2_star_E0%nab_rp_nab_sp*m%cell_volume*rdt
            o_F2_star_E0%nab_rp_nab_sp(1,:,:,:) = o_F2_star_E0%nab_rp_nab_sp(2,:,:,:)
        endif

        
        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to evolve field        ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to load boundary       ',tt2/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set field           ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field        ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source           ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time -----------------------'
            write(*,*) 'Elapsed time to add adj source      ',tt6/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve adj field    ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj field    ',tt8/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set adj field       ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract&write fields',tt10/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate           ',tt11/mpiworld%max_threads

        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_rfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_rfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')
        call hud('ximage < snap_*  n1='//num2str(cb%mz)//' perc=99')
        call hud('xmovie < snap_*  n1='//num2str(cb%mz)//' n2='//num2str(cb%mx)//' clip=?e-?? loop=2 title=%g')
        
    end subroutine

    subroutine adjoint_vp2(self,fld_F1,fld_F2,fld_dE,fld_E0,oif_record_adjseismo,o_F1_star_E0,o_F2_star_dE,o_F2_star_E0)
        class(t_propagator) :: self
        type(t_field) :: fld_F1,fld_F2,fld_dE,fld_E0
        logical,optional :: oif_record_adjseismo
        type(t_correlation),optional :: o_F1_star_E0,o_F2_star_dE,o_F2_star_E0

        real,parameter :: time_dir=-1. !time direction
        logical :: if_record_adjseismo
        logical :: if_corr

        if_record_adjseismo =either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))
        if(if_record_adjseismo) then
            call alloc(fld_F1%seismo,1,self%nt)
            call alloc(fld_F2%seismo,1,self%nt)
        endif
        if_corr = present(o_F1_star_E0).and.present(o_F2_star_dE).and.present(o_F2_star_E0)

        ! call alloc(cb%grad,cb%mz,cb%mx,cb%my,3)

call warn('no dE and F2 modeling to save time while testing..')

        !reinitialize absorbing boundary for incident wavefield reconstruction
        call fld_E0%reinit
!         call fld_dE%reinit

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
!                 call fld_F1%check_value
                call fld_F2%check_value
                call fld_E0%check_value
!                 call fld_dE%check_value
            endif

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then
                !backward step 5:
                call cpu_time(tic)
!                 call self%evolve_pressure(fld_dE,time_dir,it)
                call self%evolve_pressure(fld_E0,time_dir,it)
                call cpu_time(toc)
                tt1=tt1+toc-tic

                !backward step 2: retrieve p^it+1 at boundary layers (BC)
                call cpu_time(tic)
!                 call fld_dE%boundary_transport_pressure('load',it)
                call fld_E0%boundary_transport_pressure('load',it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 4: s^it+1.5 -> s^it+0.5
                call cpu_time(tic)
!                 call self%update_pressure(fld_dE,time_dir,it)
                call self%update_pressure(fld_E0,time_dir,it)
                call cpu_time(toc)
                tt4=tt4+toc-tic

                !backward step 1: rm pressure from s^it+0.5
                call cpu_time(tic)
!                 call self%inject_pressure_scattering(fld_dE,fld_E0,time_dir,it)
                call self%inject_pressure(fld_E0,time_dir,it)
                call cpu_time(toc)
                tt5=tt5+toc-tic

            ! endif

            !gkpa: rf%s^it+0.5 star tilD sf%s_dt^it+0.5
            !if(if_compute_grad.and.mod(it,irdt)==0) then
            if(mod(it,irdt)==0 .and. if_corr) then
                call cpu_time(tic)
                ! call gradient_vp2_oneterm(fld_F1,fld_E0,o_F1_star_E0,it)
                ! call gradient_vp2_oneterm(fld_F2,fld_dE,o_F2_star_dE,it)
                ! call gradient_vp2_3rdterm(fld_F2,fld_dE,o_F2_star_E0,it)

                call gradient_vp2_nab_rp_nab_sp(fld_F1,fld_E0,o_F1_star_E0,it)
!                 call gradient_vp2_nab_rp_nab_sp(fld_F2,fld_dE,o_F2_star_dE,it)
!                 call gradient_vp2_nab_rp_nab_sp(fld_F2,fld_E0,o_F2_star_E0,it)
                call cpu_time(toc)
                tt11=tt11+toc-tic
            endif

            !--------------------------------------------------------!

            !adjoint step 6: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_pressure(fld_F1,time_dir,it)
!             call self%inject_pressure(fld_F2,time_dir,it)
!             call self%inject_pressure_adjoint_scattering(fld_F1,fld_F2,fld_E0,time_dir,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic


            !adjoint step 4: s^it+1.5 -> s^it+0.5
            call cpu_time(tic)
            call self%update_pressure(fld_F1,time_dir,it)
!             call self%update_pressure(fld_F2,time_dir,it)
            call cpu_time(toc)
            tt8=tt8+toc-tic
            
            !adjoint step 5
            call cpu_time(tic)
            call self%evolve_pressure(fld_F1,time_dir,it)
!             call self%evolve_pressure(fld_F2,time_dir,it)
            call cpu_time(toc)
            tt7=tt7+toc-tic

            !adjoint step 1: sample v^it or s^it+0.5 at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_F1,it)
!                 call self%extract(fld_F2,it)
                call cpu_time(toc)
                tt10=tt10+toc-tic
            endif

            ! !adjoint step 6: inject to s^it+1.5 at receivers
            ! call cpu_time(tic)
            ! call self%inject_pressure_adjoint_scattering(fld_F1,fld_F2,fld_E0,time_dir,it)
            ! call cpu_time(toc)
            ! tt6=tt6+toc-tic
            
            !snapshot
            call fld_F1%write(it,o_suffix='_rev')
!             call fld_F2%write(it,o_suffix='_rev')
            call fld_E0%write(it,o_suffix='_rev')
!             call fld_dE%write(it,o_suffix='_rev')
            
            if(if_corr) then
                call fld_E0%write_ext(it,'corr_F1_star_E0',o_F1_star_E0%nab_rp_nab_sp,m%n)
!                 call fld_E0%write_ext(it,'corr_F2_star_dE',o_F2_star_dE%nab_rp_nab_sp,m%n)
!                 call fld_E0%write_ext(it,'corr_F2_star_E0',o_F2_star_E0%nab_rp_nab_sp,m%n)
            endif
            
        enddo
        
        if(if_corr) then
            ! call gradient_postprocess
            ! !scale by m%cell_volume*rdt tobe a gradient in the discretized world
            ! cb%corr = cb%corr*m%cell_volume*rdt
            ! !preparing for cb%project_back
            ! cb%corr(1,:,:,:) = cb%corr(2,:,:,:)
            o_F1_star_E0%nab_rp_nab_sp = o_F1_star_E0%nab_rp_nab_sp*m%cell_volume*rdt
            o_F1_star_E0%nab_rp_nab_sp(1,:,:,:) = o_F1_star_E0%nab_rp_nab_sp(2,:,:,:)

!             o_F2_star_dE%nab_rp_nab_sp = o_F2_star_dE%nab_rp_nab_sp*m%cell_volume*rdt
!             o_F2_star_dE%nab_rp_nab_sp(1,:,:,:) = o_F2_star_dE%nab_rp_nab_sp(2,:,:,:)
!
!             o_F2_star_E0%nab_rp_nab_sp = o_F2_star_E0%nab_rp_nab_sp*m%cell_volume*rdt
!             o_F2_star_E0%nab_rp_nab_sp(1,:,:,:) = o_F2_star_E0%nab_rp_nab_sp(2,:,:,:)
        endif

        
        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to evolve field        ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to load boundary       ',tt2/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set field           ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update field        ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source           ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time -----------------------'
            write(*,*) 'Elapsed time to add adj source      ',tt6/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve adj field    ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj field    ',tt8/mpiworld%max_threads
            ! write(*,*) 'Elapsed time to set adj field       ',tt9/mpiworld%max_threads
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
            
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler
            
            ! if(shot%src%comp=='p') then !explosion
                f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) + wl!*shot%src%interp_coef
            ! endif
                
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
            wl = f%wavelet(i,it)*wavelet_scaler 
                
            ! if(shot%rcv(i)%comp=='p') then
                f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) +wl!*shot%rcv(i)%interp_coef !no time_dir needed!
            ! endif

        enddo
        
    end subroutine

    subroutine inject_pressure_scattering(self,fld_d,fld,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: fld_d, fld
        real :: time_dir

        ! if(.not.fld%is_adjoint) then
        !     fld_d%p(1:cb%mz,1:cb%mx,1:cb%my) = fld_d%p(1:cb%mz,1:cb%mx,1:cb%my) + time_dir*W2Idt*fld%p(1:cb%mz,1:cb%mx,1:cb%my)*self%kpa(1:cb%mz,1:cb%mx,1:cb%my)
        !     return
        ! endif

        if(m%is_cubic) then !1/3D
            ! fld_d%p(1:cb%mz,1:cb%mx,1:cb%my) = fld_d%p(1:cb%mz,1:cb%mx,1:cb%my) +         self%tilD_dt2*fld%p(1:cb%mz,1:cb%mx,1:cb%my)
        else !2D
            !fld_d%p = fld_d%p + self%tilD_dt*dt*fld%p

            !方案1 w/ backward FD
            !fld_d%p(1:cb%mz,1:cb%mx,1:1) = fld_d%p(1:cb%mz,1:cb%mx,1:1) +   self%invsqEref*self%tilD_dt(1:cb%mz,1:cb%mx,1:1)*(fld%p(1:cb%mz,1:cb%mx,1:1)**2-fld%p_prev(1:cb%mz,1:cb%mx,1:1)**2)
            !fld_d%p = fld_d%p + self%invsqEref*self%tilD_dt*(fld%p**2-fld%p_prev**2)

            !方案1 w/ averaged FD, the forward_scattering subroutine should be updated
            !fld_d%p = fld_d%p + self%invsqEref*self%tilD_dt*((fld%p**2-fld%p_prev**2)+(fld%p_next**2-fld%p**2))/2

            !方案2, the forward_scattering subroutine should be updated
            !fld_d%p = fld_d%p + self%invsqEref*self%tilD_dt*dt*(fld%p_next**2-2*fld%p**2+fld%p_prev**2)

            !方案1.2
            ! fld_d%p = fld_d%p + time_dir*cb%tilD*(abs(fld%p_next)-2*abs(fld%p)+abs(fld%p_prev))

            ! fld_d%p = fld_d%p + time_dir*cb%tilD*fld%p
            !fld_d%p = fld_d%p + time_dir*cb%tilD*( fld%p -fld%p_prev )
            fld_d%p = fld_d%p + time_dir*cb%tilD*(fld%p_next-2*fld%p+fld%p_prev)


            ! ! do ix=cb%ifx,cb%ilx
            ! ! do iz=cb%ifz,cb%ilz
            do ix=cb%ifx+1,cb%ilx-1
            do iz=cb%ifz+1,cb%ilz-1
                ixm1=either(ix-1,ix,ix>cb%ifx)
                ixp1=either(ix+1,ix,ix<cb%ilx)
                izm1=either(iz-1,iz,iz>cb%ifz)
                izp1=either(iz+1,iz,iz<cb%ilz)

                dp1=(self%tilD_vp2(izp1,ix  ,1)+self%tilD_vp2(iz  ,ix  ,1)) * ((fld%p(izp1,ix  ,1))-(fld%p(iz  ,ix  ,1))) !(abs(fld%p(izp1,ix  ,1))-abs(fld%p(iz  ,ix  ,1)))
                dm1=(self%tilD_vp2(iz  ,ix  ,1)+self%tilD_vp2(izm1,ix  ,1)) * ((fld%p(iz  ,ix  ,1))-(fld%p(izm1,ix  ,1))) !(abs(fld%p(iz  ,ix  ,1))-abs(fld%p(izm1,ix  ,1)))
                ddz=(dp1-dm1)*0.5*inv_dz2

                dp1=(self%tilD_vp2(iz  ,ixp1,1)+self%tilD_vp2(iz  ,ix  ,1)) * ((fld%p(iz  ,ixp1,1))-(fld%p(iz  ,ix  ,1))) !(abs(fld%p(iz  ,ixp1,1))-abs(fld%p(iz  ,ix  ,1)))
                dm1=(self%tilD_vp2(iz  ,ix  ,1)+self%tilD_vp2(iz  ,ixm1,1)) * ((fld%p(iz  ,ix  ,1))-(fld%p(iz  ,ixm1,1))) !(abs(fld%p(iz  ,ix  ,1))-abs(fld%p(iz  ,ixm1,1)))
                ddx=(dp1-dm1)*0.5*inv_dx2

                fld_d%p(iz,ix,1) = fld_d%p(iz,ix,1) -time_dir*(ddz+ddx)*dt2

            enddo; enddo

        endif


    end subroutine

    subroutine inject_pressure_adjoint_scattering(self,fld_F1,fld_F2,fld_E0,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: fld_F1,fld_F2,fld_E0
        real :: time_dir

        ! if(m%is_cubic) then !1/3D
        !     ! fld_F1%p(1:cb%mz,1:cb%mx,1:cb%my) = fld_F1%p(1:cb%mz,1:cb%mx,1:cb%my) +         self%tilD_dt2*fld_F1%p(1:cb%mz,1:cb%mx,1:cb%my)
        ! else !2D
        !     fld_F1%p(1:cb%mz,1:cb%mx,1:1) = fld_F1%p(1:cb%mz,1:cb%mx,1:1) -2*self%invsqEref*fld_E0%p(1:cb%mz,1:cb%mx,1:1)*self%tilD_dt(1:cb%mz,1:cb%mx,1:1)*(fld_F2%p(1:cb%mz,1:cb%mx,1:1)-fld_F2%p_prev(1:cb%mz,1:cb%mx,1:1))
        ! endif


            !方案1.2
            ! fld_F1%p = fld_F1%p + cb%tilD*(fld_F2%p_prev-2*fld_F2%p+fld_F2%p_next)*sgns(fld_E0%p)

            !fld_F1%p_prev = fld_F1%p_prev + cb%tilD*fld_F2%p
            ! fld_F1%p      = fld_F1%p      + cb%tilD*(fld_F2%p-fld_F2%p_next)
            fld_F1%p = fld_F1%p + cb%tilD*(fld_F2%p_prev-2*fld_F2%p+fld_F2%p_next)

            ! do ix=cb%ifx,cb%ilx
            ! do iz=cb%ifz,cb%ilz
            do ix=cb%ifx+1,cb%ilx-1
            do iz=cb%ifz+1,cb%ilz-1
                ixm1=either(ix-1,ix,ix>cb%ifx)
                ixp1=either(ix+1,ix,ix<cb%ilx)
                izm1=either(iz-1,iz,iz>cb%ifz)
                izp1=either(iz+1,iz,iz<cb%ilz)

                dp1=(self%tilD_vp2(izp1,ix  ,1)+self%tilD_vp2(iz  ,ix  ,1)) * ((fld_F2%p(izp1,ix  ,1))-(fld_F2%p(iz  ,ix  ,1)))!(abs(fld_F2%p(izp1,ix  ,1))-abs(fld_F2%p(iz  ,ix  ,1)))
                dm1=(self%tilD_vp2(iz  ,ix  ,1)+self%tilD_vp2(izm1,ix  ,1)) * ((fld_F2%p(iz  ,ix  ,1))-(fld_F2%p(izm1,ix  ,1)))!(abs(fld_F2%p(iz  ,ix  ,1))-abs(fld_F2%p(izm1,ix  ,1)))
                ddz=(dp1-dm1)*0.5*inv_dz2

                dp1=(self%tilD_vp2(iz  ,ixp1,1)+self%tilD_vp2(iz  ,ix  ,1)) * ((fld_F2%p(iz  ,ixp1,1))-(fld_F2%p(iz  ,ix  ,1)))
                dm1=(self%tilD_vp2(iz  ,ix  ,1)+self%tilD_vp2(iz  ,ixm1,1)) * ((fld_F2%p(iz  ,ix  ,1))-(fld_F2%p(iz  ,ixm1,1)))
                ddx=(dp1-dm1)*0.5*inv_dx2

                fld_F1%p(iz,ix,1) = fld_F1%p(iz,ix,1) -(ddz+ddx)*dt2!*sgn(fld_E0%p(iz,ix,1))

            enddo; enddo

    end subroutine

    pure function sgns(a) result(s)
    use m_math, only: r_eps
        real,dimension(:,:,:),intent(in) :: a
        real,dimension(:,:,:),allocatable :: s
        s=a
        where (a>r_eps)
            s=1.
        elsewhere (a<-r_eps)
            s=-1.
        elsewhere
            s=0.
        endwhere
    end function

    pure function sgn(a) result(s)
    use m_math, only: r_eps
        real,intent(in) :: a
        if (a>r_eps) then
            s=1.
        elseif (a<-r_eps) then
            s=-1.
        else
            s=0.
        endif
    end function

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

    !forward: s^it+0.5 -> s^it+1.5 by FD of v^it+1
    !adjoint: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
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
            call fd2d_laplacian(f%p,                            &
                                f%dp_dz,f%dp_dx,f%dpzz_dz,f%dpxx_dx,&
                                f%lapz,f%lapx,&
                                self%buoz,self%buox,self%kpa,           &
                                ifz,f%bloom(2,it),f%bloom(3,it),f%bloom(4,it),f%is_adjoint)
        endif


        if(.not.f%is_adjoint) then
            if(time_dir>0.) then !in forward time
                f%p_next = 2*f%p -f%p_prev +dt2*self%kpa*(f%lapz + f%lapx)
            else !in reverse time
                f%p_prev = 2*f%p -f%p_next +dt2*self%kpa*(f%lapz + f%lapx)
            endif

            return

        endif


        !is adjoint        
        !f%p_prev=f%p_prev+2*f%p+dt2*       (f%lapz + f%lapx)  !kpa has been multiplied in fd2d_laplacian
        !replace f%p_prev by f%p_next
        !since f%p_prev=f%p_next in the previous time step
        !f%p_prev=f%p_next+2*f%p+dt2*       (f%lapz + f%lapx)  !kpa has been multiplied in fd2d_laplacian
        !f%p=-f%p

        f%p_prev=-f%p_next+2*f%p+dt2*       (f%lapz + f%lapx)  !kpa has been multiplied in fd2d_laplacian

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

        if(.not.f%is_adjoint) then
            if(time_dir>0.) then !in forward time
                f%p_prev = f%p
                f%p      = f%p_next
                !f%p_next = f%p_prev

            else !in reverse time
                f%p_next = f%p
                f%p      = f%p_prev
                !f%p_prev = f%p_next

            endif
            
            return

        endif

        !is adjoint
        f%p_next=f%p
        f%p      = f%p_prev
        !f%p_prev = f%p_next
        
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
    !if Au=∂ₜ²u-vp²∇²u
    !for vp²: Kₘ<a|Au> = -aᵀ∇²u

    subroutine gradient_tilD(rf,sf,grad,dt,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlation) :: grad
        real :: dt

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        ! if(m%is_cubic) then
        !     ! call grad3d_moduli(rf%p,sf%vz,sf%vx,sf%vy,&
        !     !                    grad,                  &
        !     !                    ifz,ilz,ifx,ilx,ify,ily)
        ! else
        !     call grad2d_tilD((rf%p_next-rf%p_prev)/2/dt,(sf%p_next-sf%p_prev)/2/dt,&
        !                       grad%sp_rp,dt,     &
        !                       ifz,ilz,ifx,ilx)
            
        ! endif

        call grad2d_inverse_scattering(rf%p_next,rf%p,rf%p_prev,sf%p_next,sf%p,sf%p_prev,&
                              grad%drp_dt_dsp_dt,grad%nab_rp_nab_sp,                     &
                              ifz,ilz,ifx,ilx)

    end subroutine

    subroutine gradient_vp2_nab_rp_nab_sp(rf,sf,grad,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlation) :: grad

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        if(m%is_cubic) then
        else
            call grad2d_vp2_nab_rp_nab_sp(rf%p,sf%p,&
                            grad%nab_rp_nab_sp,&
                            ifz,ilz,ifx,ilx)
        endif

    end subroutine


    ! subroutine gradient_vp2_oneterm(rf,sf,grad,it)
    !     type(t_field), intent(in) :: rf, sf
    !     type(t_correlation) :: grad

    !     !nonzero only when sf touches rf
    !     ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
    !     ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
    !     ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
    !     ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
    !     ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
    !     ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
    !     if(m%is_cubic) then
    !         ! call grad3d_moduli(rf%p,sf%vz,sf%vx,sf%vy,&
    !         !                    grad,                  &
    !         !                    ifz,ilz,ifx,ilx,ify,ily)
    !     else
    !         ! call grad2d_moduli(rf%p,sf%p,sf_p_save,&
    !         !                    grad,            &
    !         !                    ifz,ilz,ifx,ilx)
    !         ! sf_p_save = sf%p
            
    !         !inexact greadient
    !         call grad2d_vp2_oneterm(rf%p,sf%lapz,sf%lapx,&
    !                         grad%sp_rp,          &
    !                         ifz,ilz,ifx,ilx)
    !     endif
        
    ! end subroutine

    ! subroutine gradient_vp2_3rdterm(rf,sf,grad,it)
    !     type(t_field), intent(in) :: rf, sf
    !     type(t_correlation) :: grad

    !     !nonzero only when sf touches rf
    !     ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
    !     ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
    !     ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
    !     ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
    !     ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
    !     ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
    !     if(m%is_cubic) then
    !         ! call grad3d_moduli(rf%p,sf%vz,sf%vx,sf%vy,&
    !         !                    grad,                  &
    !         !                    ifz,ilz,ifx,ilx,ify,ily)
    !     else
    !         ! call grad2d_moduli(rf%p,sf%p,sf_p_save,&
    !         !                    grad,            &
    !         !                    ifz,ilz,ifx,ilx)
    !         ! sf_p_save = sf%p
            
    !         !inexact greadient
    !         call grad2d_vp2_3rdterm(rf%p,sf%p,&
    !                         grad%nab_rp_nab_sp,          &
    !                         ifz,ilz,ifx,ilx)
    !     endif

    ! end subroutine
    
    ! subroutine gradient_postprocess

    !     !scale the kernel tobe a gradient in the discretized world
    !     cb%grad = cb%grad*m%cell_volume*rdt
        
    !     ! !grho
    !     ! cb%grad(:,:,:,1) = cb%grad(:,:,:,1) / cb%rho(1:cb%mz,1:cb%mx,1:cb%my)
    !     !gkpa
    !     !cb%grad(:,:,:,2) = cb%grad(:,:,:,2) * (-ppg%inv_kpa(1:cb%mz,1:cb%mx,1:cb%my))
    !     !gvp2
                
    !     !preparing for cb%project_back
    !     cb%grad(1,:,:,:) = cb%grad(2,:,:,:)

    ! end subroutine

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
                do iz = ifz+2,ilz-1
                    dp_dz_ = c1z*(p(iz,ix) - p(iz-1,ix)) + c2z*(p(iz+1,ix) - p(iz-2,ix))

                    dp_dz(iz,ix) = cpml%b_z_half(iz)*dp_dz(iz,ix) + cpml%a_z_half(iz)*dp_dz_

                    dp_dz_ = dp_dz_/cpml%kpa_z_half(iz) + dp_dz(iz,ix)

                    pzz(iz,ix) =  buoz(iz,ix)*dp_dz_
                enddo
            enddo

            do ix = ifx+2,ilx-1
                do iz = ifz,ilz
                    dp_dx_ = c1x*(p(iz,ix) - p(iz,ix-1)) + c2x*(p(iz,ix+1) - p(iz,ix-2))

                    dp_dx(iz,ix) = cpml%b_x_half(ix)*dp_dx(iz,ix) + cpml%a_x_half(ix)*dp_dx_

                    dp_dx_ = dp_dx_/cpml%kpa_x_half(ix) + dp_dx(iz,ix)

                    pxx(iz,ix) = buox(iz,ix)*dp_dx_
                enddo
            enddo

            !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
            do ix = ifx,ilx
                do iz = ifz+1,ilz-2
                    dpzz_dz_ = c1z*(pzz(iz+1,ix) - pzz(iz,ix)) + c2z*(pzz(iz+2,ix) - pzz(iz-1,ix))

                    dpzz_dz(iz,ix) = cpml%b_z(iz)*dpzz_dz(iz,ix) + cpml%a_z(iz)*dpzz_dz_

                    dpzz_dz_ = dpzz_dz_/cpml%kpa_z(iz) + dpzz_dz(iz,ix)

                    lapz(iz,ix) = dpzz_dz_ !*kpa(iz,ix) !keep kpa later because lapz is needed for gradient computatin
                enddo
            enddo

            do ix = ifx+1,ilx-2
                do iz = ifz,ilz
                    dpxx_dx_ = c1x*(pxx(iz,ix+1) - pxx(iz,ix)) + c2x*(pxx(iz,ix+2) - pxx(iz,ix-1))

                    dpxx_dx(iz,ix) = cpml%b_x(ix)*dpxx_dx(iz,ix) + cpml%a_x(ix)*dpxx_dx_

                    dpxx_dx_ = dpxx_dx_/cpml%kpa_x(ix) + dpxx_dx(iz,ix)

                    lapx(iz,ix) = dpxx_dx_ !*kpa(iz,ix) !keep kpa later because lapz is needed for gradient computatin
                enddo
            enddo

        else

            !flux: b∇κa ~= ( bz*∂zᶠκ*pᵃ , bx*∂ₓᶠκ*pᵃ )
            !dir$ simd
            do ix = ifx,ilx
                do iz = ifz+1,ilz-2
                    dp_dz_ = c1z*(kpa(iz+1,ix)*p(iz+1,ix) - kpa(iz,ix)*p(iz,ix)) + c2z*(kpa(iz+2,ix)*p(iz+2,ix) - kpa(iz-1,ix)*p(iz-1,ix))

                    dp_dz(iz,ix) = cpml%b_z_half(iz)*dp_dz(iz,ix) + cpml%a_z_half(iz)*dp_dz_

                    dp_dz_ = dp_dz_/cpml%kpa_z_half(iz) + dp_dz(iz,ix)

                    pzz(iz,ix) =  buoz(iz,ix)*dp_dz_
                enddo
            enddo

            do ix = ifx+1,ilx-2
                do iz = ifz,ilz
                    dp_dx_ = c1x*(kpa(iz,ix+1)*p(iz,ix+1) - kpa(iz,ix)*p(iz,ix)) + c2x*(kpa(iz,ix+2)*p(iz,ix+2) - kpa(iz,ix-1)*p(iz,ix-1))

                    dp_dx(iz,ix) = cpml%b_x_half(ix)*dp_dx(iz,ix) + cpml%a_x_half(ix)*dp_dx_

                    dp_dx_ = dp_dx_/cpml%kpa_x_half(ix) + dp_dx(iz,ix)

                    pxx(iz,ix) = buox(iz,ix)*dp_dx_
                enddo
            enddo

            !laplacian: ∇·b∇κa ~= ∂zᵇ(bz*∂zᶠpᵃ) + ∂ₓᵇ(bx*∂ₓᶠpᵃ)
            do ix = ifx,ilx
                do iz = ifz+2,ilz-1
                    dpzz_dz_ = c1z*(pzz(iz,ix) - pzz(iz-1,ix)) + c2z*(pzz(iz+1,ix) - pzz(iz-2,ix))

                    dpzz_dz(iz,ix) = cpml%b_z(iz)*dpzz_dz(iz,ix) + cpml%a_z(iz)*dpzz_dz_

                    dpzz_dz_ = dpzz_dz_/cpml%kpa_z(iz) + dpzz_dz(iz,ix)

                    lapz(iz,ix) = dpzz_dz_
                enddo
            enddo

            do ix = ifx+2,ilx-1
                do iz = ifz,ilz
                    dpxx_dx_ = c1x*(pxx(iz,ix) - pxx(iz,ix-1)) + c2x*(pxx(iz,ix+1) - pxx(iz,ix-2))

                    dpxx_dx(iz,ix) = cpml%b_x(ix)*dpxx_dx(iz,ix) + cpml%a_x(ix)*dpxx_dx_

                    dpxx_dx_ = dpxx_dx_/cpml%kpa_x(ix) + dpxx_dx(iz,ix)

                    lapx(iz,ix) = dpxx_dx_
                enddo
            enddo

        endif

    end subroutine

    ! subroutine grad2d_tilD(rf_p_next,rf_p,sf_p,&
    !                        grad,            &
    !                        ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: rf_p_next,rf_p,sf_p
    !     real,dimension(*) :: grad
        
    !     nz=cb%nz
        
    !     !$omp parallel default (shared)&
    !     !$omp private(iz,ix,i,j)
    !     !$omp do schedule(dynamic)
    !     do ix=ifx,ilx
        
    !         !dir$ simd
    !         do iz=ifz,ilz
                
    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
    !             j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
    !             grad(j)=grad(j) + (rf_p_next(i)-rf_p(i))*inv_dt *sf_p(i)*sf_p(i)
                                
    !         end do
            
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

    subroutine grad2d_inverse_scattering(rf_p_next,rf_p,rf_p_prev,sf_p_next,sf_p,sf_p_prev,&
                      grad_drp_dt_dsp_dt,grad_nab_rp_nab_sp,&
                      ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_p_next,rf_p,rf_p_prev
        real,dimension(*) :: sf_p_next,sf_p,sf_p_prev
        real,dimension(*) :: grad_drp_dt_dsp_dt, grad_nab_rp_nab_sp

        !double precision :: drp_dt, dsp_dt, drp_dz, dsp_dz, drp_dx, dsp_dx

        nz=cb%nz

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers

                
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                
                iz_ixm1=i    -nz  !iz,ix-1
                iz_ixp1=i    +nz  !iz,ix+1
                
                !for time diff, we have 3 choices
                !which choice is correct depends on where we form the crosscorrelation
                !Choice I:
                drp_dt = (rf_p_next(i)-rf_p_prev(i)) *inv_2dt
                dsp_dt = (sf_p_next(i)-sf_p_prev(i)) *inv_2dt
                grad_drp_dt_dsp_dt(j)=grad_drp_dt_dsp_dt(j) + drp_dt*dsp_dt

                !!Choice II:
                ! drp_dt = -rf_p(i)
                ! d2sp_dt2 = (sf_p_next(i)-2*sf_p(i)+sf_p_prev(i)) /dt2
                ! grad_drp_dt_dsp_dt(j)=grad_drp_dt_dsp_dt(j) - rf_p(i)*d2sp_dt2

                !!Choice III:
                ! d2rp_dt2 = (rf_p_next(i)-2*rf_p(i)+rf_p_prev(i)) /dt2
                ! grad_drp_dt_dsp_dt(j)=grad_drp_dt_dsp_dt(j) - d2rp_dt2*sf_p(i)


                !drp_dz = (rf_p(izp1_ix)-rf_p(izm1_ix)) *inv_2dz
                !dsp_dz = (sf_p(izp1_ix)-sf_p(izm1_ix)) *inv_2dz
                !drp_dx = (rf_p(iz_ixp1)-rf_p(iz_ixm1)) *inv_2dx
                !dsp_dx = (sf_p(iz_ixp1)-sf_p(iz_ixm1)) *inv_2dx

                drp_dz = rgc1z*(rf_p(izp1_ix)-rf_p(izm1_ix)) + rgc2z*(rf_p(izp2_ix)-rf_p(izm2_ix))
                dsp_dz = rgc1z*(sf_p(izp1_ix)-sf_p(izm1_ix)) + rgc2z*(sf_p(izp2_ix)-sf_p(izm2_ix))

                drp_dx = rgc1x*(rf_p(iz_ixp1)-rf_p(iz_ixm1)) + rgc2x*(rf_p(iz_ixp2)-rf_p(iz_ixm2))
                dsp_dx = rgc1x*(sf_p(iz_ixp1)-sf_p(iz_ixm1)) + rgc2x*(sf_p(iz_ixp2)-sf_p(iz_ixm2))

                grad_nab_rp_nab_sp(j)=grad_nab_rp_nab_sp(j) + drp_dz*dsp_dz + drp_dx*dsp_dx

            end do
            
        end do
        !$omp end do
        !$omp end parallel


    end subroutine

    subroutine grad2d_vp2_nab_rp_nab_sp(rf_p,sf_p,     &
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
                
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                
                iz_ixm1=i    -nz  !iz,ix-1
                iz_ixp1=i    +nz  !iz,ix+1
                
                drp_dz = rgc1z*(rf_p(izp1_ix)-rf_p(izm1_ix)) + rgc2z*(rf_p(izp2_ix)-rf_p(izm2_ix))
                dsp_dz = rgc1z*(sf_p(izp1_ix)-sf_p(izm1_ix)) + rgc2z*(sf_p(izp2_ix)-sf_p(izm2_ix))

                drp_dx = rgc1x*(rf_p(iz_ixp1)-rf_p(iz_ixm1)) + rgc2x*(rf_p(iz_ixp2)-rf_p(iz_ixm2))
                dsp_dx = rgc1x*(sf_p(iz_ixp1)-sf_p(iz_ixm1)) + rgc2x*(sf_p(iz_ixp2)-sf_p(iz_ixm2))

                grad(j)=grad(j) + drp_dz*dsp_dz + drp_dx*dsp_dx

            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

    ! subroutine grad2d_vp2_oneterm(rf_p,sf_lapz,sf_lapx,&
    !                               grad,          &
    !                               ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: rf_p,sf_lapz,sf_lapx
    !     real,dimension(*) :: grad

    !     nz=cb%nz
        
    !     !$omp parallel default (shared)&
    !     !$omp private(iz,ix,i,j)
    !     !$omp do schedule(dynamic)
    !     do ix=ifx,ilx
        
    !         !dir$ simd
    !         do iz=ifz,ilz
                
    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
    !             j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
    !             grad(j)=grad(j) + rf_p(i)*(sf_lapz(i)+sf_lapx(i))
                
    !         end do
            
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

    ! subroutine grad2d_vp2_3rdterm(rf_p,sf_p,     &
    !                               grad,          &
    !                               ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: rf_p,sf_p
    !     real,dimension(*) :: grad
        
    !     ! double precision :: drp_dz, dsp_dz, drp_dx, dsp_dx

    !     nz=cb%nz
        
    !     !$omp parallel default (shared)&
    !     !$omp private(iz,ix,i,j)
    !     !$omp do schedule(dynamic)
    !     do ix=ifx,ilx
        
    !         !dir$ simd
    !         do iz=ifz,ilz
                
    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
    !             j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
    !             izm1_ix=i-1  !iz-1,ix
    !             iz_ix  =i    !iz,ix
    !             izp1_ix=i+1  !iz+1,ix
                
    !             iz_ixm1=i    -nz  !iz,ix-1
    !             iz_ixp1=i    +nz  !iz,ix+1
                
    !             drp_dz = (    rf_p(izp1_ix) -    rf_p(izm1_ix)) !*inv_2dz
    !             dsp_dz = (abs(sf_p(izp1_ix))-abs(sf_p(izm1_ix)))!*inv_2dz

    !             drp_dx = (    rf_p(iz_ixp1) -    rf_p(iz_ixm1)) !*inv_2dx
    !             dsp_dx = (abs(sf_p(iz_ixp1))-abs(sf_p(iz_ixm1)))!*inv_2dx

    !             grad(j)=grad(j) + drp_dz*dsp_dz + drp_dx*dsp_dx
                
    !         end do
            
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

    ! subroutine grad2d_moduli(rf_p,sf_lapz,sf_lapx,&
    !                          grad,          &
    !                          ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: rf_p,sf_lapz,sf_lapx
    !     real,dimension(*) :: grad
        
    !     nz=cb%nz
        
    !     !$omp parallel default (shared)&
    !     !$omp private(iz,ix,i,j)
    !     !$omp do schedule(dynamic)
    !     do ix=ifx,ilx
        
    !         !dir$ simd
    !         do iz=ifz,ilz
                
    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
    !             j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
    !             grad(j)=grad(j) + rf_p(i)*(sf_lapz(i)+sf_lapx(i))
                
    !         end do
            
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

    ! subroutine grad2d_D(rf_p,sf_p,&
    !                          grad,          &
    !                          ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: rf_p,sf_p
    !     real,dimension(*) :: grad
        
    !     nz=cb%nz
        
    !     !$omp parallel default (shared)&
    !     !$omp private(iz,ix,i,j)
    !     !$omp do schedule(dynamic)
    !     do ix=ifx,ilx
        
    !         !dir$ simd
    !         do iz=ifz,ilz
                
    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
    !             j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
    !             grad(j)=grad(j) + rf_p(i)*sf_p(i)
                
    !         end do
            
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

end
