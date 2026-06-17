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
    real :: inv_dt, inv_4dt

    !old sfield
    real,dimension(:,:,:),allocatable :: old_Hx, old_Hz, old_Ey

    !scaling source wavelet
    real :: wavelet_scaler

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain isotropic Transverse Electric (TE) wave propagation'//s_NL// &
            '1st-order H-E formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.606(for 2D) or 0.494(3D) *Vmax/dx'//s_NL// &
            'Required model attributes: eps, sgma, mu'//s_NL// &
            'Required field components: Hx, Hy, Ey'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Imaging conditions: iEyEy'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%Ey² dt'//s_NL// &
            'Basic gradients: gmu geps gsgma'

        integer :: nbndlayer=max(2,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=3 !number of basic gradients
        integer :: nimag=1 !number of basic images
        integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:),allocatable :: imuz, imux, eps, sgma

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
        procedure :: init_correlate
        procedure :: init_abslayer
        procedure :: assemble

        procedure :: forward
        procedure :: adjoint
        
        !procedure :: inject_H
        procedure :: inject_E
        procedure :: update_H
        procedure :: update_E
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    logical,public :: propagator_if_record_adjseismo=.false.

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
        
        ! if(index(self%info,'vs')>0  .and. .not. allocated(m%vs)) then
        !     call alloc(m%vs,m%nz,m%nx,1,o_init=866.)
        !     call warn('Constant vs model (866 m/s) is allocated by propagator.')
        ! endif

        ! if(index(self%info,'rho')>0 .and. .not. allocated(m%rho)) then
        !     call alloc(m%rho,m%nz,m%nx,m%ny,o_init=1000.)
        !     call warn('Constant rho model (1000 kg/m³) is allocated by propagator.')
        ! endif
                
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

        inv_dt  = 1./dt
        inv_4dt = 1./4./dt
        
    end subroutine

    subroutine init(self)
        class(t_propagator) :: self

        real,dimension(:,:),allocatable :: temp_imu

        c1x=coef(1)/m%dx; c1y=coef(1)/m%dy; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz

        ! inv_2dz =1./2/m%dz
        ! inv_2dx =1./2/m%dx
        
        wavelet_scaler=self%dt/m%cell_volume

        if_hicks=shot%if_hicks

        call alloc(self%eps, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%sgma,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%imuz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%imux,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        
        call alloc(temp_imu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        temp_imu(:,:) = 1./cb%mu(:,:,1)

        do iz=cb%ifz+1,cb%ilz
            self%imuz(iz,:)=(temp_imu(iz,:)+temp_imu(iz-1,:))/2.
        enddo
        
        do ix=cb%ifx+1,cb%ilx
            self%imux(:,ix)=(temp_imu(:,ix)+temp_imu(:,ix-1))/2.
        enddo

        deallocate(temp_imu)

        self%eps =cb%eps (:,:,1)
        self%sgma=cb%sgma(:,:,1)

        !initialize m_field
        call field_init(.true.,self%nt,self%dt)

        !initialize m_correlate
        call correlate_init(self%nt,self%dt)

        !rectified interval for time integration
        !default to Nyquist, and must be a multiple of dt
        rdt=setup%get_real('REF_RECT_TIME_INTEVAL','RDT',o_default=num2str(0.5/shot%fmax))
        irdt=floor(rdt/self%dt)
        if(irdt==0) irdt=1
        rdt=irdt*self%dt
        call hud('rdt, irdt = '//num2str(rdt)//', '//num2str(irdt))

    end subroutine

    subroutine init_field(self,f,name,ois_adjoint)
        class(t_propagator) :: self
        type(t_field) :: f
        character(*) :: name
        logical,optional :: ois_adjoint

        !field
        ! call f%init(name)
        f%name=name

        f%is_adjoint=either(ois_adjoint,.false.,present(ois_adjoint))

        call f%init_bloom

        call f%init_boundary_magnetic

        call alloc(f%Hx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%Hz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%Ey,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])

        call alloc(f%dHx_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dHz_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dEy_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dEy_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
                
    end subroutine

    subroutine init_correlate(self,corr,name)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name
        
        corr%name=name

        ! if(name(1:1)=='g') then !gradient components
            call alloc(corr%gmu,m%nz,m%nx,m%ny)
            call alloc(corr%geps,m%nz,m%nx,m%ny)
            call alloc(corr%gsgma,m%nz,m%nx,m%ny)
        !else !image components
        !    call alloc(corr%ipp,m%nz,m%nx,m%ny)
        !    call alloc(corr%ibksc,m%nz,m%nx,m%ny)
        !    call alloc(corr%ifwsc,m%nz,m%nx,m%ny)
        ! endif

    end subroutine

    subroutine init_abslayer(self)
        class(t_propagator) :: self

        call cpml%init

    end subroutine

    subroutine assemble(self,corr)
        class(t_propagator) :: self
        type(t_correlate) :: corr

        !if(allocated(correlate_image)) then
        !    call correlate_assemble(corr%ipp, correlate_image(:,:,:,1))
        !    call correlate_assemble(corr%ibksc, correlate_image(:,:,:,2))
        !    call correlate_assemble(corr%ifwsc, correlate_image(:,:,:,3))
        !endif

        if(allocated(correlate_gradient)) then
            call correlate_assemble(corr%gmu,  correlate_gradient(:,:,:,1))
            call correlate_assemble(corr%geps, correlate_gradient(:,:,:,2))
            call correlate_assemble(corr%gsgma,correlate_gradient(:,:,:,3))
        endif        
        
    end subroutine

    
    !========= Derivations =================
    !PDE:      A u = Mₚ∂ₜ u + Md u - D u - f = 0
    !Adjoint:  Aᵀa = Mₚ∂ₜᵀa - Md a - Dᵀa - d = 0
    !where
    !u=[Hx,Hz,Ey]ᵀ, where H=[Hx,0,-Hz]ᵀ is the magnetic field, E=[0 Ey 0]ᵀ is the electric field
    !f=[0,0,Jy]ᵀδ(x-xs) is the source term, where J=[0,-Jy,0] is the injection current density vector,
    !and d is recorded data
    !Mₚ=diag[μ μ ε], μ is magnetic permeability, ε is electric permittivity,
    !Md=diag[0 0 σ], σ is electric conductivity,
    !  [0  0  ∂z]
    !D=|0  0  ∂ₓ|
    !  [∂z ∂ₓ 0 ]
    !and a=[Hxᵃ Hzᵃ Eyᵃ]ᵀ is the adjoint field
    !
    !Continuous case:
    !<a|Au> = ∫ a(x,t) (Mₚ∂ₜ+Md-D)u(x,t) dx³dt
    !Integration by parts, eg.:
    !∫aᵀMₚ∂ₜu dt = aMₚu|ₜ₌₀ᵀ - ∫(∂ₜa)ᵀMₚu dt, and freely choosing a(t=T)=0 (final condition),
    !∫aᵀMₚ∂ₜu dt = -∫(∂ₜa)ᵀMₚu dt
    !Similar procedure on spatial derivatives, we have
    !∫aᵀDu dx³ = -∫(Dᵀa)ᵀ u dx³, with same boundary conditions on a
    !Therefore, Aᵀa = Mₚ∂ₜa - Md a -Dᵀa
    !However, this method (finding the adjoint FD eqn by integration by parts)
    !is NOT accurate enough in the discrete world to pass the adjoint test.
    !
    !Discrete case:
    !Meshing with staggered grids in time and space:
    !                       |        |    -½ Hx        |        |
    !                       |        |       μz⁻¹      |        |
    !                       |        |        |        |        |
    !                      εσ   μx  εσ   μx  εσ   μx  εσ   μx  εσ  ⁻¹
    !  -H-Ey-H-Ey-H-→ t   -Ey---Hz--Ey---Hz--Ey---Hz--Ey---Hz--Ey-→ x
    !  -1 -½ 0  ½ 1        -2  -1½  -1   -½   0    ½   1   1½   2    
    !                       |        |        |        |        | 
    !                       |        |     ½ Hx        |        | 
    !                       |        |       μz⁻¹      |        |
    !                       |        |        |        |        | 
    !                      -|--------|-----1-Ey--------|--------|-
    !                       |        |       εσ        |        | 
    !                       |        |        |        |        | 
    !                       |        |    1½ Hx        |        | 
    !                       |        |       μz⁻¹      |        | 
    !                       |        |        |        |        | 
    !                                       z ↓
    !
    !Convention for half-integer index:
    !(array index)  =>    (real index)     
    !Ey(iz,ix)  => Ey[iz,  ix, ]^n+½ :=Ey(iz*dz,ix*dx,(n+½)*dt)
    !Hx(iz,ix)  => Hx[iz-½,ix, ]^n   :=Hx((iz-½)*dz,ix*dx,n*dt)
    !Hz(iz,ix)  => Hz[iz,  ix-½]^n   
    !
    ! Forward:
    ! FD eqn:
    !           [Hx^n  ]   [ 0    0   ∂zᵇ][Hx^n+1]      [∂zᵇ Ey^n+½             ]
    ! (Mₚ∂ₜᶠ+Md)|Hz^n  | = | 0    0   ∂ₓᵇ||Hz^n+1| +f = |∂ₓᵇ Ey^n+½             | +f
    !           [Ey^n+1]   [∂zᶠ  ∂ₓᶠ   0 ][Ey^n+½]      [∂zᶠ Hx^n+1 + ∂ₓᶠ Hx^n+1]
    ! where
    ! ∂ₜᶠ*dt := v^n+1 - v^n                                ~O(t²)
    ! ∂zᵇ*dz := c₁(v(iz  )-v(iz-1)) +c₂(v(iz+1)-v(iz-2))  ~O(x⁴)
    ! ∂zᶠ*dz := c₁(s(iz+1)-s(iz  )) +c₂(s(iz+2)-s(iz-1))  ~O(x⁴)
    
    ! Time marching:    
    !           [Hx^n  ]               [Hx^n+1 ]               [Hx^n  ])
    ! (Mₚ∂ₜᶠ+Md)|Hz^n  | = (Mₚ/dt+Md/2)|Hz^n+1 | - (Mₚ/dt-Md/2)|Hz^n  |)
    !           [Ey^n+1]               [Ey^n+1½]               [Ey^n+½])
    !                     [ 0    0   ∂zᵇ][Hx^n+1]      [∂zᵇ Ey^n+½             ]
    !                    =| 0    0   ∂ₓᵇ||Hz^n+1| +f = |∂ₓᵇ Ey^n+½             | +f
    !                     [∂zᶠ  ∂ₓᶠ   0 ][Ey^n+½]      [∂zᶠ Hx^n+1 + ∂ₓᶠ Hz^n+1]
    !that is,
    !            [Hx^n+1 ]               [Hx^n  ]   [∂zᵇ Ey^n+½             ]
    !(Mₚ/dt+Md/2)|Hz^n+1 | - (Mₚ/dt-Md/2)|Hz^n  | = |∂ₓᵇ Ey^n+½             | +f
    !            [Ey^n+1½]               [Ey^n+½]   [∂zᶠ Hx^n+1 + ∂ₓᶠ Hz^n+1]
    !
    ! [Hx^n+1 ]                 [            [Hx^n  ]   [∂zᵇ Ey^n+½             ]   ]
    ! |Hz^n+1 | = (Mₚ/dt+Md/2)⁻¹[(Mₚ/dt-Md/2)|Hz^n  | + |∂ₓᵇ Ey^n+½             | +f]
    ! [Ey^n+1½]                 [            [Ey^n+½]   [∂zᶠ Hy^n+1 + ∂ₓᶠ Hy^n+1]   ]
    !compared w/ SH propagator, swap Steps 1 & 2 and 3 & 4 here
    ! Step #1: H^n+1 = H^n + spatial FD(E^n+½)
    !! Step #2: H^n += src
    ! Step #3: E^n+1½ = E^n+½ + spatial FD(H^n+1)
    ! Step #4: E^n+½ += src
    ! Step #5: sample H^n & E^n+½ at receivers
    ! Step #6: save H^n+1 to boundary values

    ! Reverse time marching (for wavefield reconstruction)
    !            [Hx^n  ]               [Hx^n+1 ]   [∂zᵇ Ey^n+½             ]
    !(Mₚ/dt-Md/2)|Hz^n  | = (Mₚ/dt+Md/2)|Hz^n+1 | - |∂ₓᵇ Ey^n+½             | -f
    !            [Ey^n+½]               [Ey^n+1½]   [∂zᶠ Hx^n+1 + ∂ₓᶠ Hz^n+1]
    !
    !            [Ey^n+½]               [Ey^n+1½]   [∂zᶠ Hx^n+1 + ∂ₓᶠ Hz^n+1]
    !(Mₚ/dt-Md/2)|Hz^n  | = (Mₚ/dt+Md/2)|Hz^n+1 | - |∂ₓᵇ Ey^n+½             | -f
    !            [Hx^n  ]               [Hx^n+1 ]   [∂zᵇ Ey^n+½             ]
    !
    ![Ey^n+½]                 [            [Ey^n+1½]   [∂zᶠ Hx^n+1 + ∂ₓᶠ Hz^n+1]    ]
    !|Hz^n  | = (Mₚ/dt-Md/2)⁻¹[(Mₚ/dt+Md/2)|Hz^n+1 | - |∂ₓᵇ Ey^n+½             | -f ]
    ![Hx^n  ]                 [            [Hx^n+1 ]   [∂zᵇ Ey^n+½             ]    ]
    !
    ! Step #6: load boundary values for v^n+1
    ! Step #4: E^n+½ -= src
    ! Step #3: E^n+½ = E^n+1½ - spatial FD(H^n+1)
    ! Step #2: H^n -= src
    ! Step #1: H^n+1 = E^n - spatial FD(E^n+½)
    ! N.B. Same codes for spatial FDs as in forward time marching, with a negated dt.
    
    ! Adjoint:
    ! FD eqn:
    !           [Hxᵃ^n  ]   [∂zᵇ Eyᵃ^n+½              ]
    ! (Mₚ∂ₜᶠ-Md)|Hzᵃ^n  | = |∂ₓᵇ Eyᵃ^n+½              | -d
    !           [Eyᵃ^n+1]   [∂zᶠ Hxᵃ^n+1 + ∂ₓᶠ Hxᵃ^n+1]
    !because ∂ₜᶠᵀ = -∂ₜᵇ, ∂zᵇᵀ = -∂zᶠ, ∂zᶠᵀ = -∂zᵇ, Dᵀ=-D (antisymmetric)
    
    ! Time marching:
    !             [Hxᵃ^n+1 ]               [Hxᵃ^n  ]   [∂zᵇ Eyᵃ^n+½              ]
    ! (Mₚ/dt-Md/2)|Hzᵃ^n+1 | - (Mₚ/dt+Md/2)|Hzᵃ^n  | = |∂ₓᵇ Eyᵃ^n+½              | -d
    !             [Eyᵃ^n+1½]               [Eyᵃ^n+½]   [∂zᶠ Hxᵃ^n+1 + ∂ₓᶠ Hzᵃ^n+1]
    !but we evolve it in reverse time:
    !             [Eyᵃ^n+½]               [Eyᵃ^n+1½]   [∂zᶠ Hxᵃ^n+1 + ∂ₓᶠ Hzᵃ^n+1]
    ! (Mₚ/dt+Md/2)|Hzᵃ^n  | - (Mₚ/dt-Md/2)|Hzᵃ^n+1 | =-|∂ₓᵇ Eyᵃ^n+½              | +d
    !             [Hxᵃ^n  ]               [Hxᵃ^n+1 ]   [∂zᵇ Eyᵃ^n+½              ]
    !
    ! [Eyᵃ^n+½]                 [            [Eyᵃ^n+1½]  [∂zᶠ Hxᵃ^n+1 + ∂ₓᶠ Hzᵃ^n+1]   ]
    ! |Hzᵃ^n  | = (Mₚ/dt+Md/2)⁻¹[(Mₚ/dt-Md/2)|Hzᵃ^n+1 | -|∂ₓᵇ Eyᵃ^n+½              | +d]
    ! [Hxᵃ^n  ]                 [            [Hxᵃ^n+1 ]  [∂zᵇ Eyᵃ^n+½              ]   ]
    !
    ! Step #5: Eᵃ^n+½ = Eᵃ^n+1½ - spatial FD(Hᵃ^n+1)
    ! Step #4: Eᵃ^n+1½ += adjsrc
    ! Step #3: Hᵃ^n = Hᵃ^n+1 - spatial FD(Hᵃ^n+½)
    !! Step #2: Hᵃ^n+1 += adjsrc
    ! N.B. Same codes for spatial FDs as in forward time marching, just negate the sign of the Laplacian
    ! and do NOT negate the sign of the RHS
    
!!    ! For adjoint test:
!!    ! In each step of forward time marching: dsyn=RGANf
!!    ! f:source wavelet, N=M⁻¹dt: diagonal
!!    ! A:inject source into field, G:propagator, R:extract field at receivers
!!    ! while in each step of reverse-time adjoint marching: dadj=NAᵀGᵀRᵀdsyn=AᵀGᵀRᵀNdsyn
!!    ! Rᵀ:inject adjoint sources, Gᵀ:adjoint propagator, Aᵀ:extract adjoint fields
    

    subroutine forward(self,fld_u)
        class(t_propagator) :: self
        type(t_field) :: fld_u

        real,parameter :: time_dir=1. !time direction

        !seismo
        call alloc(fld_u%seismo,shot%nrcv,self%nt)
            
        tt1=0.; tt2=0.; tt3=0.; tt4=0.; tt5=0.; tt6=0.; tt7=0.

        ift=1; ilt=self%nt

        do it=ift,ilt
            if(mod(it,100)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_u%check_value(fld_u%Ey)
            endif

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !unlike m_propagator_SH_FDSG_O4
            !here, let's update the field before injection..

            !step 2: from H^it to H^it+1 by differences of E^it+0.5
            call cpu_time(tic)
            call self%update_H(fld_u,time_dir,it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            !!step 1: add forces to v^it, no need
            !call cpu_time(tic)
            !call self%inject_H(fld_u,time_dir,it)
            !call cpu_time(toc)
            !tt1=tt1+toc-tic

            !step 4: from E^it+0.5 to E^it+1.5 by differences of H^it+1
            call cpu_time(tic)
            call self%update_E(fld_u,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !step 3: add current density vector to E^it+0.5
            call cpu_time(tic)
            call self%inject_E(fld_u,time_dir,it)
            call cpu_time(toc)
            tt3=tt3+toc-tic

            !step 5: sample E^it+1.5 at receivers
            call cpu_time(tic)
            call self%extract(fld_u,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !snapshot
            call fld_u%write(it)

            !step 6: save H^it+1 in boundary layers
            ! if(fld_u%if_will_reconstruct) then
                call cpu_time(tic)
                call fld_u%boundary_transport_magnetic('save',it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            ! endif

        enddo

        if(mpiworld%is_master) then
            !write(*,*) 'Elapsed time to add stress source  ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update H    ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to add current density vector',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update E  ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field      ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary      ',tt6/mpiworld%max_threads
            write(*,*) 'Total elapsed time (min):',(tt1+tt2+tt3+tt4+tt5+tt6)/60./mpiworld%max_threads
        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')

    end subroutine

    subroutine adjoint(self,fld_a,fld_u, a_star_u)
    !adjoint_a_star_Du
        class(t_propagator) :: self
        type(t_field) :: fld_a,fld_u
        type(t_correlate) :: a_star_u
        
        real,parameter :: time_dir=-1. !time direction

        !reinitialize absorbing boundary for incident wavefield reconstruction
        call fld_u%reinit

        call alloc(old_Ey, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(old_Hx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(old_Hz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        
        !for adjoint test
        if(propagator_if_record_adjseismo)  call alloc(fld_a%seismo,1,self%nt)

        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.
        
        ift=1; ilt=self%nt

        do it=ilt,ift,int(time_dir)
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_a%check_value(fld_a%Ey)
                call fld_u%check_value(fld_u%Ey)
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            old_Ey = fld_u%Ey

            !backward step 6: retrieve H^it+1 at boundary layers (BC)
            call cpu_time(tic)
            call fld_u%boundary_transport_magnetic('load',it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !backward step 3: rm force from E^it+0.5
            call cpu_time(tic)
            call self%inject_E(fld_u,time_dir,it)
            call cpu_time(toc)
            tt3=tt3+toc-tic
        
            !backward step 4: E^it+1.5 -> E^it+0.5 by FD of H^it+1
            call cpu_time(tic)
            call self%update_E(fld_u,time_dir,it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            !--------------------------------------------------------!

            !adjoint step 4: E^it+1.5 -> E^it+0.5 by FD^T of H^it+1
            call cpu_time(tic)
            call self%update_E(fld_a,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !adjoint step 5: inject to E^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_E(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !gkpa: rf%s^it+0.5 star D sf%s_dt^it+0.5
            !use sf%v^it+1 to compute sf%s_dt^it+0.5, as backward step

            if(mod(it,irdt)==0) then
                call cpu_time(tic)
                call cross_correlate_geps_gsgma(fld_a,fld_u,a_star_u,it)
                ! call cross_correlate_image(fld_a,fld_u,a_star_u,it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif
                            
            !========================================================!
            old_Hx = fld_u%Hx
            old_Hz = fld_u%Hz

            !!backward step 1: rm  from H^it
            !call cpu_time(tic)
            !call self%inject_H(fld_u,time_dir,it)
            !call cpu_time(toc)
            !tt8=tt8+toc-tic

            !backward step 2: H^it+1 -> H^it by FD of E^it+0.5
            call cpu_time(tic)
            call self%update_H(fld_u,time_dir,it)
            call cpu_time(toc)
            tt7=tt7+toc-tic


            !--------------------------------------------------------!

            !adjoint step 2: H^it+1 -> H^it by FD^T of E^it+0.5
            call cpu_time(tic)
            call self%update_H(fld_a,time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            !!adjoint step 3: inject to v^it+1 at receivers
            !call cpu_time(tic)
            !call self%inject_H(fld_a,time_dir,it)
            !call cpu_time(toc)
            !tt9=tt9+toc-tic

            !adjoint step 1: sample E^it+0.5 at source position
            if(propagator_if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_a,it)
                call cpu_time(toc)
                tt11=tt11+toc-tic
            endif
            
!            !grho: sfield%v_dt^it \dot rfield%v^it
!            !use sfield%s^it+0.5 to compute sfield%v_dt^it, as backward step 2
           if(mod(it,irdt)==0) then
               call cpu_time(tic)
               call cross_correlate_gmu(fld_a,fld_u,a_star_u,it)
               call cpu_time(toc)
               tt6=tt6+toc-tic
           endif
            
            !snapshot
            call fld_a%write(it,o_suffix='_rev')
            call fld_u%write(it,o_suffix='_rev')

            call a_star_u%write(it,o_suffix='_rev')

        enddo
        
        !postprocess
        call cross_correlate_postprocess(a_star_u)
        call a_star_u%scale(m%cell_volume*rdt)
        
        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to load boundary         ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update E     ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source E  ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update H       ',tt7/mpiworld%max_threads
!            write(*,*) 'Elapsed time to rm source stresses    ',tt8/mpiworld%max_threads
            write(*,*) 'Total elapsed time for forward (min)',(tt1+tt2+tt3+tt7+tt8)/60./mpiworld%max_threads
            write(*,*) ' ---------------------------- '
            write(*,*) 'Elapsed time to add adjsource E ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj E    ',tt5/mpiworld%max_threads
!            write(*,*) 'Elapsed time to add adjsource stresses   ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj H      ',tt10/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract&write fields     ',tt11/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate                ',tt6/mpiworld%max_threads
            write(*,*) 'Total elapsed time for adjoint&correlate (min)',(tt4+tt5+tt9+tt10+tt11+tt6)/60./mpiworld%max_threads
            write(*,*) 'Total elapsed time (min):',(tt1+tt2+tt3+tt7+tt8+tt4+tt5+tt9+tt10+tt11+tt6)/60./mpiworld%max_threads

        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_rfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_rfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')
        call hud('ximage < snap_*  n1='//num2str(cb%mz)//' perc=99')
        call hud('xmovie < snap_*  n1='//num2str(cb%mz)//' n2='//num2str(cb%mx)//' clip=?e-?? loop=2 title=%g')
        
    end subroutine


    !forward: add RHS to v^it
    !adjoint: add RHS to v^it+1
    subroutine inject_E(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        if(.not. f%is_adjoint) then

            if(shot%src%comp=='Ey') then

                ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
                ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
                
                wl=time_dir*f%wavelet(1,it)*wavelet_scaler
                ! if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl !required to pass adjointtest.
                
            
                if(if_hicks) then
                    f%Ey(ifz:ilz,ifx:ilx,1) = f%Ey(ifz:ilz,ifx:ilx,1) + wl/self%eps(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)
                    
                else
                    f%Ey(iz,ix,1) = f%Ey(iz,ix,1) + wl/self%eps(iz,ix)
                
                endif
            
            endif

            return

        endif

            do i=1,shot%nrcv

                if(shot%rcv(i)%comp=='Ey') then  !horizontal y adjsource !Ey[ix,iy-0.5,iz]
                    
                    ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                    ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                    
                    wl=f%wavelet(i,it)*wavelet_scaler
                    ! if(m%is_freesurface.and.shot%rcv(i)%iz==1) wl=2*wl !required to pass adjointtest.
                    
                    if(if_hicks) then
                        f%Ey(ifz:ilz,ifx:ilx,1) = f%Ey(ifz:ilz,ifx:ilx,1) + wl/self%eps(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                        
                    else
                        f%Ey(iz,ix,1) = f%Ey(iz,ix,1) + wl/self%eps(iz,ix) !no time_dir needed!
                    
                    endif
               
                endif

            enddo
        
    end subroutine
    
    !forward: v^it -> v^it+1 by FD  of s^it+0.5
    !adjoint: v^it+1 -> v^it by FDᵀ of s^it+0.5
    subroutine update_H(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ifz=f%bloom(1,it)+2
        ilz=f%bloom(2,it)-1
        ifx=f%bloom(3,it)+2
        ilx=f%bloom(4,it)-1
        ify=f%bloom(5,it)+2
        ily=f%bloom(6,it)-1

        if(m%is_freesurface) ifz=max(ifz,1)

        call fd2d_H(f%Hx,f%Hz,f%Ey,                 &
                    f%dEy_dz,f%dEy_dx,              &
                    self%imuz,self%imux,            &
                    ifz,ilz,ifx,ilx,time_dir*self%dt)

!        if(m%is_freesurface) then
!            !apply free surface boundary condition if needed
!            !Levandar & Roberttson's stress image method
!            f%sxy(1,:,1)=0.
!            f%sxy(0:cb%ifz:-1, :,1)=-f%sxy(2:2+0-cb%ifz, :,1)
!
!            !image szy
!            f%szy(1:cb%ifz:-1, :,1)=-f%szy(2:2+1-cb%ifz, :,1)
            
!        endif

    end subroutine

    ! !forward: add RHS to s^it+0.5
    ! !adjoint: add RHS to s^it+1.5
    ! subroutine inject_H(self,f,time_dir,it)
    !     class(t_propagator) :: self
    !     type(t_field) :: f

    !     if(.not. f%is_adjoint) then

    !         ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
    !         ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
    !         wl=time_dir*f%wavelet(1,it)*wavelet_scaler
            
    !         if(if_hicks) then
    !             if(shot%src%comp=='szy') then
    !                 f%szy(ifz:ilz,ifx:ilx,1) = f%szy(ifz:ilz,ifx:ilx,1) + wl*self%muz(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)

    !             else if(shot%src%comp=='sxy') then
    !                 f%sxy(ifz:ilz,ifx:ilx,1) = f%sxy(ifz:ilz,ifx:ilx,1) + wl*self%mux(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)
                
    !             endif
                
    !         else
    !             if(shot%src%comp=='szy') then
    !                 f%szy(iz,ix,1) = f%szy(iz,ix,1) + wl*self%muz(iz,ix)
                
    !             else if(shot%src%comp=='sxy') then
    !                 f%sxy(iz,ix,1) = f%sxy(iz,ix,1) + wl*self%mux(iz,ix)
                
    !             endif
                
    !         endif

    !         return

    !     endif

    !         do i=1,shot%nrcv

    !             ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
    !             ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                
    !             !adjsource for pressure
    !             wl=f%wavelet(i,it)*wavelet_scaler
                
    !             if(if_hicks) then 

    !                 if(shot%rcv(i)%comp=='szy') then
    !                     f%szy(ifz:ilz,ifx:ilx,1) = f%szy(ifz:ilz,ifx:ilx,1) +wl*self%muz(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
    !                 elseif(shot%rcv(i)%comp=='sxy') then
    !                     f%sxy(ifz:ilz,ifx:ilx,1) = f%sxy(ifz:ilz,ifx:ilx,1) +wl*self%mux(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
    !                 endif

    !             else

    !                 if(shot%rcv(i)%comp=='szy') then
    !                     !szy[iz,ix,1]
    !                     f%szy(iz,ix,1) = f%szy(iz,ix,1) +wl*self%muz(iz,ix) !no time_dir needed!
    !                 elseif(shot%rcv(i)%comp=='sxy') then
    !                     !sxy[iz,ix,1]
    !                     f%sxy(iz,ix,1) = f%sxy(iz,ix,1) +wl*self%mux(iz,ix) !no time_dir needed!
    !                 endif

    !             endif

    !         enddo
        
    ! end subroutine

    !forward: s^it+0.5 -> s^it+1.5 by FD of v^it+1
    !adjoint: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
    subroutine update_E(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ifz=f%bloom(1,it)+1
        ilz=f%bloom(2,it)-2
        ifx=f%bloom(3,it)+1
        ilx=f%bloom(4,it)-2
        
        if(m%is_freesurface) ifz=max(ifz,1)

        call fd2d_E(f%Hx,f%Hz,f%Ey,                 &
                    f%dHx_dz,f%dHz_dx,              &
                    self%eps,self%sgma,             &
                    ifz,ilz,ifx,ilx,time_dir*self%dt)
        
        !if(m%is_freesurface) then
        !    !apply free surface boundary condition if needed
        !    !Levandar & Roberttson's stress image method
        !    f%vy(cb%ifz:0,:,1)=0.
	!
        !endif

    end subroutine

    subroutine extract(self,f,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        if(.not.f%is_adjoint) then

            do i=1,shot%nrcv
                ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1

                if(if_hicks) then
                    select case (shot%rcv(i)%comp)
                        
                        case ('Ey')
                        f%seismo(i,it)=sum(f%Ey(ifz:ilz,ifx:ilx,1) *shot%rcv(i)%interp_coef(:,:,1))
                        case ('Hx')
                        f%seismo(i,it)=sum(f%Hx(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1))
                        case ('Hz')
                        f%seismo(i,it)=sum(f%Hz(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1))
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                        case ('Ey') !p[iz,ix,1]
                        f%seismo(i,it)=f%Ey(iz,ix,1)
                        case ('Hx') !vz[iz-0.5,ix,1]
                        f%seismo(i,it)=f%Hx(iz,ix,1)
                        case ('Hz') !vx[iz,ix-0.5,1]
                        f%seismo(i,it)=f%Hz(iz,ix,1)
                    end select
                    
                endif

            enddo

            return

        endif

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            if(if_hicks) then
                select case (shot%src%comp)
                    case ('Ey')
                    f%seismo(1,it)=sum(f%Ey(ifz:ilz,ifx:ilx,1) *shot%src%interp_coef(:,:,1))
                    
                    case ('Hx')
                    f%seismo(1,it)=sum(f%Hx(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                    
                    case ('Hz')
                    f%seismo(1,it)=sum(f%Hz(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                    
                end select
                
            else
                select case (shot%src%comp)
                    case ('Ey') !p[iz,ix,1]
                    f%seismo(1,it)=f%Ey(iz,ix,1)
                    
                    case ('Hx') !vz[iz-0.5,ix,1]
                    f%seismo(1,it)=f%Hx(iz,ix,1)
                    
                    case ('Hz') !vx[iz,ix-0.5,1]
                    f%seismo(1,it)=f%Hz(iz,ix,1)
                    
                end select
                
            endif
        
    end subroutine
    
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%imuz, self%imux, self%eps, self%sgma)
        call dealloc(old_Hx, old_Hz, old_Ey)
    end subroutine


    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !Kₘ<a|Au> = Kₘ<a|Mₚ∂ₜu +Mdu -Du>
    !Since it's cumbersome to get ∂ₜu by time marching,
    !replace ∂ₜu by -Mₚ⁻¹(Mdu-Du) and neglect f

    subroutine cross_correlate_geps_gsgma(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)    
        
        call grad2d_geps_gsgma(rf%Ey,sf%Ey,old_Ey,&
                         corr%geps, corr%gsgma,   &
                         ifz,ilz,ifx,ilx)
        
    end subroutine

    subroutine cross_correlate_gmu(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        
        !inexact greadient
        call grad2d_gmu(rf%Hx,rf%Hz,sf%Hx,sf%Hz,old_Hx,old_Hz,&
                        corr%gmu,                             &
                        ifz,ilz,ifx,ilx                       )

    end subroutine

    subroutine cross_correlate_image(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr
        
        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        
        ! if(m%is_cubic) then
        !     call imag3d_xcorr(rf%p,sf%p,&
        !                       imag,                  &
        !                       ifz,ilz,ifx,ilx,ify,ily)
        ! else
        !    call imag2d(rf%p,sf%p,&
        !                rf%poynz,rf%poynx,sf%poynz,sf%poynx, &
        !                corr%ipp,corr%ibksc,corr%ifwsc, &
        !                ifz,ilz,ifx,ilx)
        ! endif

        ! call imag2d_xcorr(rf%p,rf%vz,rf%vx,&
        !                   sf%p,sf%vz,sf%vx,&
        !                   imag,            &
        !                   ifz,ilz,ifx,ilx)

    end subroutine

    subroutine cross_correlate_postprocess(corr)
        type(t_correlate) :: corr
        
        if(allocated(correlate_gradient)) then
            
            !remove singular point at the src position,
            !because we didn't consider src when deriving the gradient formula
            iz=shot%src%iz-cb%ioz+1
            ix=shot%src%ix-cb%iox+1
            !for point source
            !safeguards
            ifz=either(iz-2,iz,iz>=3); ilz=either(iz+2,iz,iz<=m%nz-2)
            ifx=either(ix-2,ix,ix>=3); ilx=either(ix+2,ix,ix<=m%nx-2)
            
            !first remove otherwise will appear in the sum below
            corr%geps (iz,ix,1)=0
            corr%gsgma(iz,ix,1)=0
            corr%gmu  (iz,ix,1)=0

            !then replace singular point by sum
            ncells=size(corr%geps(ifz:ilz,ifx:ilx,1))-1
            corr%geps (iz,ix,1) = sum(corr%geps (ifz:ilz,ifx:ilx,1))/ncells
            corr%gsgma(iz,ix,1) = sum(corr%gsgma(ifz:ilz,ifx:ilx,1))/ncells
            corr%gmu  (iz,ix,1) = sum(corr%gmu  (ifz:ilz,ifx:ilx,1))/ncells
            
            !remove singular top boundary
            corr%geps (1,:,:) = corr%geps (2,:,:)
            corr%gsgma(1,:,:) = corr%gsgma(2,:,:)
            corr%gmu  (1,:,:) = corr%gmu  (2,:,:)

        endif

        ! if(allocated(correlate_image)) then
        !     corr%ipp (1,:,:) = corr%ipp (2,:,:)
        !     corr%ibksc(1,:,:) = corr%ibksc(2,:,:)
        !     corr%ifwsc(1,:,:) = corr%ifwsc(2,:,:)
        ! endif

    end subroutine

    !========= Finite-Difference on flattened arrays ==================

    subroutine fd2d_H(Hx,Hz,Ey,       &
                      dEy_dz,dEy_dx,    &
                      imuz,imux,        &
                      ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: Hx,Hz,Ey
        real,dimension(*) :: dEy_dz,dEy_dx
        real,dimension(*) :: imuz,imux
        
        nz=cb%nz
        nx=cb%nx
        
        dEy_dz_=0.; dEy_dx_=0.

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dEy_dz_,dEy_dx_)
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

                dEy_dz_= c1z*(Ey(iz_ix)-Ey(izm1_ix)) +c2z*(Ey(izp1_ix)-Ey(izm2_ix))
                dEy_dx_= c1x*(Ey(iz_ix)-Ey(iz_ixm1)) +c2x*(Ey(iz_ixp1)-Ey(iz_ixm2))

                !cpml
                dEy_dz(iz_ix)= cpml%b_z_half(iz)*dEy_dz(iz_ix) + cpml%a_z_half(iz)*dEy_dz_
                dEy_dx(iz_ix)= cpml%b_x_half(ix)*dEy_dx(iz_ix) + cpml%a_x_half(ix)*dEy_dx_

                dEy_dz_=dEy_dz_*cpml%kpa_z_half(iz) + dEy_dz(iz_ix)
                dEy_dx_=dEy_dx_*cpml%kpa_x_half(ix) + dEy_dx(iz_ix)

                !H
                Hx(iz_ix)=Hx(iz_ix) + dt*imuz(iz_ix)*dEy_dz_
                Hz(iz_ix)=Hz(iz_ix) + dt*imux(iz_ix)*dEy_dx_

            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd2d_E(Hx,Hz,Ey,       &
                      dHx_dz,dHz_dx,  &
                      eps,sgma,              &
                      ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: Hx,Hz,Ey
        real,dimension(*) :: dHx_dz,dHz_dx
        real,dimension(*) :: eps,sgma
        
        nz=cb%nz
        nx=cb%nx
        
        dHx_dz_=0.; dHz_dx_=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dHx_dz_,dHz_dx_)
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
                
                dHx_dz_= c1z*(Hx(izp1_ix)-Hx(iz_ix)) +c2z*(Hx(izp2_ix)-Hx(izm1_ix))
                dHz_dx_= c1x*(Hz(iz_ixp1)-Hz(iz_ix)) +c2x*(Hz(iz_ixp2)-Hz(iz_ixm1))
                
                !cpml
                dHx_dz(iz_ix)=cpml%b_z(iz)*dHx_dz(iz_ix)+cpml%a_z(iz)*dHx_dz_
                dHz_dx(iz_ix)=cpml%b_x(ix)*dHz_dx(iz_ix)+cpml%a_x(ix)*dHz_dx_

                dHx_dz_=dHx_dz_*cpml%kpa_z(iz) + dHx_dz(iz_ix)
                dHz_dx_=dHz_dx_*cpml%kpa_x(ix) + dHz_dx(iz_ix)
                
                !εdE + σE = ∂zHx + ∂ₓHz
                !ε(E^n+1-E^n)/dt + σ(E^n+1+E^n)/2 = ∂zHx + ∂ₓHz
                !(ε/dt+σ/2) (E^n+1) -(ε/dt-σ/2)(E^n)   = ∂zHx + ∂ₓHz
                Ey(iz_ix) = ((eps(iz_ix)/dt-sgma(iz_ix)/2)*Ey(iz_ix) + dHx_dz_+dHz_dx_) &
                           / (eps(iz_ix)/dt+sgma(iz_ix)/2)
                
            enddo
            
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine
    

    subroutine grad2d_gmu(rf_Hx,rf_Hz,sf_Hx,sf_Hz,old_Hx,old_Hz,&
                          grad,                                 &
                          ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_Hx,rf_Hz,sf_Hx,sf_Hz,old_Hx,old_Hz
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        dvy_dz=0.
        dvy_dx=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         iz_ix,izp1_ix,iz_ixp1)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers

                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix            
                iz_ixp1=i  +nz  !iz,ix+1
                
                grad(j)=grad(j) +( (rf_Hx(izp1_ix)+rf_Hx(iz_ix))*(sf_Hx(izp1_ix)+sf_Hx(iz_ix)-old_Hx(izp1_ix)-old_Hx(iz_ix)) &
                                  +(rf_Hz(iz_ixp1)+rf_Hx(iz_ix))*(sf_Hx(iz_ixp1)+sf_Hx(iz_ix)-old_Hx(iz_ixp1)-old_Hx(iz_ix)) &
                                 )*inv_4dt

            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine
    
    subroutine grad2d_geps_gsgma(rf_Ey,sf_Ey,old_Ey,&
                                 geps,gsgma,        &
                                 ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_Ey,sf_Ey,old_Ey
        real,dimension(*) :: geps,gsgma
        
        nz=cb%nz
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers
                
                geps (j)=geps (j) + rf_Ey(i)*(sf_Ey(i)-old_Ey(i)) *invdt
                gsgma(j)=gsgma(j) + rf_Ey(i)* sf_Ey(i)
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine

end
