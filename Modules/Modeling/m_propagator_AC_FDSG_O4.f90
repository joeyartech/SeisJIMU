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
            'Time-domain ISOtropic 2D/3D ACoustic propagation'//s_NL// &
            '1st-order Velocity-Stress formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.606(for 2D) or 0.494(3D) *Vmax/dx'//s_NL// &
            'Required model attributes: vp, rho'//s_NL// &
            'Required field components: vz, vx, vy(3D), p'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Poynting definitions: p_dotp_gradp, dotp_gradp, p_v, Esq_gradphi'//s_NL// &
            'Imaging conditions: ipp ibksc ifwsc (P-Pxcorr of backward & forward scattering)'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: grho gkpa'

        integer :: nbndlayer=max(2,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=2 !number of basic gradients
        integer :: nimag=3 !number of basic images
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
        procedure :: init_correlate
        procedure :: init_abslayer
        procedure :: assemble

        procedure :: forward
        procedure :: adjoint_poynting
        ! procedure :: adjoint_3terms
        
        procedure :: inject_velocities
        procedure :: inject_stresses
        procedure :: update_velocities
        procedure :: update_stresses
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    character(:),allocatable :: s_poynting_def

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    ! real,dimension(:,:,:),allocatable :: sf_p_save

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

        c1x=coef(1)/m%dx; c1y=coef(1)/m%dy; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz

        inv_2dz =1./2/m%dz
        inv_2dx =1./2/m%dx
        
        wavelet_scaler=self%dt/m%cell_volume

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

        s_poynting_def=setup%get_str('POYNTING_DEF',o_default='p_v')
        if(s_poynting_def/='p_dotp_gradp'.and.s_poynting_def/='dotp_gradp'.and.s_poynting_def/='p_v'.and.s_poynting_def/='Esq_gradphi') then
            call error('Sorry, other Poynting definitions have not yet implemented.')
        endif

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

        call alloc(f%vz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%vx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%vy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%p, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%dvz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dvx_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dvy_dy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dp_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dp_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dp_dy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%poynz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%poynx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
                
    end subroutine

    subroutine init_correlate(self,corr,name)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name
        
        corr%name=name

        if(name(1:1)=='g') then !gradient components
            call alloc(corr%gkpa,m%nz,m%nx,m%ny)
            call alloc(corr%grho,m%nz,m%nx,m%ny)
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
            call correlate_assemble(corr%ipp, correlate_image(:,:,:,1))
            call correlate_assemble(corr%ibksc, correlate_image(:,:,:,2))
            call correlate_assemble(corr%ifwsc, correlate_image(:,:,:,3))
        endif

        if(allocated(correlate_gradient)) then
            call correlate_assemble(corr%grho, correlate_gradient(:,:,:,1))
            call correlate_assemble(corr%gkpa, correlate_gradient(:,:,:,2))
        endif        
        
    end subroutine

    
    !========= Derivations =================
    !PDE:      A u = M ∂ₜ u - D u = f
    !Adjoint:  Aᵀa = M ∂ₜᵀa - Dᵀa = d
    !where
    !u=[vz vx vy p]ᵀ, [vz vx vy] are velocities, p=tr(s)=szz=sxx=syy is (hydrostatic) pressure
    !f=[fz fx fy fp]ᵀδ(x-xs), d is recorded data
    !M=[diag(ρ) κ⁻¹], N=M⁻¹=[diag(b) κ], b=ρ⁻¹ is buoyancy,
    !  [0   0   0   ∂z]
    !D=|0   0   0   ∂ₓ|
    !  |0   0   0   ∂y|
    !  [∂z  ∂ₓ  ∂y  0 ]
    !and a=[vzᵃ vxᵃ vyᵃ pᵃ]ᵀ is the adjoint field
    !
    !Continuous case:
    !<a|Au> = ∫ a(x,t) (M∂ₜ-D)u(x,t) dx³dt
    !Integration by parts, eg.:
    !∫aᵀM∂ₜu dt = aMu|ₜ₌₀ᵀ - ∫(∂ₜa)ᵀMu dt, and freely choosing a(t=T)=0 (final condition),
    !∫aᵀM∂ₜu dt = -∫(∂ₜa)ᵀMu dt
    !Similar procedure on spatial derivatives, we have
    !∫aᵀDu dx³ = -∫(Dᵀa)ᵀ u dx³, with same boundary conditions on a
    !Therefore, Aᵀa = ∂ₜa - Dᵀa
    !However, this method (finding the adjoint FD eqn by integration by parts)
    !is NOT accurate enough in the discrete world to pass the adjoint test.
    !
    !Discrete case:
    !Meshing with staggered grids in time and space (2D example):
    !                      |      |   -½ vz     |      |
    !                      |      |      bz     |      |
    !                      |      |      |      |      |
    !                      κ  bx  κ  bx  κ  bx  κ  bx  κ
    !  -v--p-v-p-v-→ t    -p--vx--p--vx--p--vx--p--vx--p-→ x
    !  -1 -½ 0 ½ 1        -2 -1½ -1  -½  0   ½  1  1½  2    
    !                      |      |      |      |      | 
    !                      |      |    ½ vz     |      | 
    !                      |      |      bz     |      | 
    !                      |      |      |      |      | 
    !                     -|------|----1-p------|------|-
    !                      |      |      κ      |      | 
    !                      |      |      |      |      | 
    !                      |      |   1½ vz     |      | 
    !                      |      |      bz     |      | 
    !                      |      |      |      |      | 
    !                                  z ↓
    !
    !Convention for half-integer index:
    !(array index)  =>    (real index)     
    ! vz(iz,ix,iy)  =>  vz[iz-½,ix,  iy  ]^n   :=vz((iz-½)*dz,ix*dx,iy*dy,n*dt)
    ! vx(iz,ix,iy)  =>  vx[iz,  ix-½,iy  ]^n  
    ! vy(iz,ix,iy)  =>  vy[iz,  ix,  iy-½]^n  
    !  p(iz,ix,iy)  =>   p[iz,  ix,  iy  ]^n+½ :=p(iz*dz,ix*dx,iy*dy,(n+½)*dt)
    !
    !Forward:
    !FD eqn:
    !      [vz^n  ]   [ 0     0    ∂zᵇ][vz^n+1]      [∂zᵇ p^n+½              ] 
    !M ∂ₜᶠ |vx^n  | = | 0     0    ∂ₓᵇ||vx^n+1| +f = |∂ₓᵇ p^n+½              | +f
    !      [ p^n+1]   [∂zᶠ   ∂ₓᶠ    0 ][ p^n+½]      [∂zᶠ vz^n+1 + ∂ₓᶠ vx^n+1]
    !where
    !∂ₜᶠ*dt := v^n+1 - v^n                               ~O(t²)
    !∂zᵇ*dz := c₁(p(iz  )-p(iz-1)) +c₂(p(iz+1)-p(iz-2))  ~O(x⁴)
    !∂zᶠ*dz := c₁(v(iz+1)-v(iz  )) +c₂(v(iz+2)-v(iz-1))  ~O(x⁴)
    !
    !Time marching:
    ![vz^n+1 ]   [vz^n  ]      [∂zᵇ p^n+½              ]
    !|vx^n+1 | = |vx^n  | + M⁻¹|∂ₓᵇ p^n+½              |dt  +M⁻¹f*dt
    ![ p^n+1½]   [ p^n+½]      [∂zᶠ vz^n+1 + ∂ₓᶠ vx^n+1]
    !Step #1: v^n += src
    !Step #2: v^n+1 = v^n + spatial FD(p^n+½)
    !Step #3: p^n+½ += src
    !Step #4: p^n+1½ = p^n+½ + spatial FD(v^n+1)
    !Step #5: sample v^n & p^n++½ at receivers
    !Step #6: save v^n+1 to boundary values
    !
    !Reverse time marching (for wavefield reconstruction)
    ![ p^n+½]   [ p^n+1½]      [∂zᶠ vz^n+1 + ∂ₓᶠ vx^n+1]
    !|vx^n  | = |vx^n+1 | - M⁻¹|∂ₓᵇ p^n+½              |dt  -M⁻¹f*dt
    ![vz^n  ]   [vz^n+1 ]      [∂zᵇ p^n+½              ]
    !Step #6: load boundary values for v^n+1
    !Step #4: p^n+½ = p^n+1½ - spatial FD(v^n+1)
    !Step #3: p^n+½ -= src
    !Step #2: v^n+1 = v^n - spatial FD(p^n+½)
    !Step #1: v^n -= src
    !N.B. Same codes for spatial FDs as in forward time marching, with a negated dt.
    !
    !Adjoint:
    !FD eqn:
    !      [vzᵃ^n  ]   [ 0      0     ∂zᶠᵀ][vzᵃ^n+1]
    !M ∂ₜᶠᵀ|vxᵃ^n  | = | 0      0     ∂ₓᶠᵀ||vxᵃ^n+1|  +d
    !      [ pᵃ^n+1]   [∂zᵇᵀ   ∂ₓᵇᵀ    0  ][ pᵃ^n+½]
    !∂ₜᶠᵀ = v^n-1 -v^n   = -∂ₜᵇ
    !∂zᵇᵀ = c₁(v[i  ]-v[i+½]) +c₂(v[i- ½]-v[i+1½]) = -∂zᶠ
    !∂zᶠᵀ = c₁(p[i-½]-p[i  ]) +c₂(p[i-1½]-p[i+ ½]) = -∂zᵇ
    !      [vzᵃ^n  ]   [ 0      0     ∂zᵇ ][vzᵃ^n+1]
    !M ∂ₜᵇ |vxᵃ^n  | = | 0      0     ∂ₓᵇ ||vxᵃ^n+1|  -d
    !      [ pᵃ^n+1]   [∂zᶠ    ∂ₓᶠ    0   ][ pᵃ^n+½]
    !ie. Dᵀ=-D, antisymmetric
    !
    !Time marching:
    ![vzᵃ^n+1 ]   [vzᵃ^n  ]      [∂zᵇ pᵃ^n+½               ]
    !|vxᵃ^n+1 | = |vxᵃ^n  | + M⁻¹|∂ₓᵇ pᵃ^n+½               |dt  -M⁻¹d*dt
    ![ pᵃ^n+1½]   [ pᵃ^n+½]      [∂zᶠ vzᵃ^n+1 + ∂ₓᶠ vxᵃ^n+1]
    !but we have to do it in reverse time:
    ![ pᵃ^n+½]   [ pᵃ^n+1½]      [∂zᶠ vzᵃ^n+1 + ∂ₓᶠ vxᵃ^n+1]
    !|vxᵃ^n  | = |vxᵃ^n+1 | - M⁻¹|∂ₓᵇ pᵃ^n+½               |dt  +M⁻¹d*dt
    ![vzᵃ^n  ]   [vzᵃ^n+1 ]      [∂zᵇ pᵃ^n+½               ]
    !Step #5: pᵃ^n+1½ += adjsrc
    !Step #4: pᵃ^n+½ = pᵃ^n+1½ - spatial FD(vᵃ^n+1)
    !Step #3: vᵃ^n+1 += adjsrc
    !Step #2: vᵃ^n = vᵃ^n+1 - spatial FD(pᵃ^n+½)
    !N.B. Same codes for spatial FDs as in forward time marching, with a negated dt, but the RHS should use a "+" sign (regardless of the reverse time direction).
    !
    !For adjoint test:
    !In each step of forward time marching: dsyn=RGANf
    !f:source wavelet, N=M⁻¹dt: diagonal
    !A:inject source into field, G:propagator, R:extract field at receivers
    !while in each step of reverse-time adjoint marching: dadj=NAᵀGᵀRᵀdsyn=AᵀGᵀRᵀNdsyn
    !Rᵀ:inject adjoint sources, Gᵀ:adjoint propagator, Aᵀ:extract adjoint fields
    !

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
            !step 1: add forces to v^it
            call cpu_time(tic)
            call self%inject_velocities(fld_u,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2: from v^it to v^it+1 by differences of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities(fld_u,time_dir,it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            !step 3: add pressure to s^it+0.5
            call cpu_time(tic)
            call self%inject_stresses(fld_u,time_dir,it)
            call cpu_time(toc)
            tt3=tt3+toc-tic

            !step 4: from s^it+0.5 to s^it+1.5 by differences of v^it+1
            call cpu_time(tic)
            call self%update_stresses(fld_u,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !step 5: sample v^it+1 or s^it+1.5 at receivers
            call cpu_time(tic)
            call self%extract(fld_u,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !snapshot
            call fld_u%write(it)

            !step 6: save v^it+1 in boundary layers
            ! if(fld_u%if_will_reconstruct) then
                call cpu_time(tic)
                call fld_u%boundary_transport('save',it)
                call cpu_time(toc)
                tt7=tt7+toc-tic
            ! endif

        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add source velocities',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update velocities    ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to add source stresses  ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update stresses      ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field        ',tt6/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary        ',tt7/mpiworld%max_threads
        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')

    end subroutine

    subroutine adjoint_poynting(self,fld_q,fld_p,fld_v,fld_u, a_star_u)
    !adjoint_a_star_Du
        class(t_propagator) :: self
        type(t_field) :: fld_q,fld_p
        type(t_field) :: fld_v,fld_u
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

        ! call alloc(sf_p_save,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        do it=ilt,ift,int(time_dir)
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_q%check_value
                call fld_p%check_value
                call fld_v%check_value
                call fld_u%check_value
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then

                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_v%boundary_transport('load',it)
                call fld_u%boundary_transport('load',it)
                call cpu_time(toc)
                tt1=tt1+toc-tic
                
                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%update_stresses(fld_v,time_dir,it)
                call self%update_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_stresses(fld_v,time_dir,it)
                call self%inject_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_stresses(fld_q,time_dir,it)
            call self%inject_stresses(fld_p,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_stresses(fld_q,time_dir,it)
            call self%update_stresses(fld_p,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !compute Poynting vectors before imaging
            call cpu_time(tic)
            call compute_poynting(fld_q,fld_p,it)
            call compute_poynting(fld_v,fld_u,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic
            
            !gkpa: rf%s^it+0.5 star D sf%s_dt^it+0.5
            !use sf%v^it+1 to compute sf%s_dt^it+0.5, as backward step

            if(mod(it,irdt)==0) then
                call cpu_time(tic)
                call cross_correlate_image(fld_p,fld_u,a_star_u,it)
                call cpu_time(toc)
                tt7=tt7+toc-tic
            endif

            ! if(mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call cross_correlate_gkpa(fld_a,fld_u,a_star_u,it)
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

            ! if(present(o_sf)) then
                !backward step 2: v^it+1 -> v^it by FD of s^it+0.5
                call cpu_time(tic)
                call self%update_velocities(fld_v,time_dir,it)
                call self%update_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt8=tt8+toc-tic

                !backward step 1: rm forces from v^it
                call cpu_time(tic)
                call self%inject_velocities(fld_v,time_dir,it)
                call self%inject_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt9=tt9+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 3: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_velocities(fld_q,time_dir,it)
            call self%inject_velocities(fld_p,time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities(fld_q,time_dir,it)
            call self%update_velocities(fld_p,time_dir,it)
            call cpu_time(toc)
            tt11=tt11+toc-tic
            
            ! !adjoint step 1: sample v^it or s^it+0.5 at source position
            ! if(if_record_adjseismo) then
            !     call cpu_time(tic)
            !     call self%extract(fld_a,it)
            !     call cpu_time(toc)
            !     tt12=tt12+toc-tic
            ! endif
            
            ! !grho: sfield%v_dt^it \dot rfield%v^it
            ! !use sfield%s^it+0.5 to compute sfield%v_dt^it, as backward step 2
            ! if(mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call cross_correlate_grho(fld_a,fld_u,a_star_u,it)
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
            
            !snapshot
            call fld_p%write(it,o_suffix='_rev')
            call fld_u%write(it,o_suffix='_rev')

            call a_star_u%write(it,o_suffix='_rev')

        enddo

        !postprocess
        call cross_correlate_postprocess(a_star_u)
        call a_star_u%scale(m%cell_volume*rdt)

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to load boundary            ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update stresses          ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source stresses       ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update velocities        ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source velocities     ',tt8/mpiworld%max_threads
            write(*,*) 'Elapsed time ----------------------------'
            write(*,*) 'Elapsed time to add adjsource stresses   ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj stresses      ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to add adjsource velocities ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj velocities    ',tt10/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract&write fields     ',tt11/mpiworld%max_threads
            write(*,*) 'Elapsed time to compute Poynting vectors ',tt6/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate                ',tt7/mpiworld%max_threads

        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_rfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_rfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')
        call hud('ximage < snap_*  n1='//num2str(cb%mz)//' perc=99')
        call hud('xmovie < snap_*  n1='//num2str(cb%mz)//' n2='//num2str(cb%mx)//' clip=?e-?? loop=2 title=%g')

    end subroutine


    !forward: add RHS to v^it
    !adjoint: add RHS to v^it+1
    subroutine inject_velocities(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        if(.not. f%is_adjoint) then

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            ify=shot%src%ify-cb%ioy+1; iy=shot%src%iy-cb%ioy+1; ily=shot%src%ily-cb%ioy+1
            
            wl=time_dir*f%wavelet(1,it)
            
            if(if_hicks) then
                select case (shot%src%comp)
                    case ('vz')
                    f%vz(ifz:ilz,ifx:ilx,ify:ily) = f%vz(ifz:ilz,ifx:ilx,ify:ily) + wl*self%buoz(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coef
                    
                    case ('vx')
                    f%vx(ifz:ilz,ifx:ilx,ify:ily) = f%vx(ifz:ilz,ifx:ilx,ify:ily) + wl*self%buox(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coef
                    
                    case ('vy')
                    f%vy(ifz:ilz,ifx:ilx,ify:ily) = f%vy(ifz:ilz,ifx:ilx,ify:ily) + wl*self%buoy(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coef
                    
                end select
                
            else
                select case (shot%src%comp)
                    case ('vz') !vertical force     on vz[iz-0.5,ix,iy]
                    f%vz(iz,ix,iy) = f%vz(iz,ix,iy) + wl*self%buoz(iz,ix,iy)
                    
                    case ('vx') !horizontal x force on vx[iz,ix-0.5,iy]
                    f%vx(iz,ix,iy) = f%vx(iz,ix,iy) + wl*self%buox(iz,ix,iy)
                    
                    case ('vy') !horizontal y force on vy[iz,ix,iy-0.5]
                    f%vy(iz,ix,iy) = f%vy(iz,ix,iy) + wl*self%buoy(iz,ix,iy)
                    
                end select
                
            endif

            return

        endif


            do i=1,shot%nrcv
                ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                ify=shot%rcv(i)%ify-cb%ioy+1; iy=shot%rcv(i)%iy-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1
                
                wl=f%wavelet(i,it)
                
                if(if_hicks) then
                    select case (shot%rcv(i)%comp)
                        case ('vz') !vertical z adjsource
                        f%vz(ifz:ilz,ifx:ilx,ify:ily) = f%vz(ifz:ilz,ifx:ilx,ify:ily) + wl*self%buoz(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef !no time_dir needed!

                        case ('vx') !horizontal x adjsource
                        f%vx(ifz:ilz,ifx:ilx,ify:ily) = f%vx(ifz:ilz,ifx:ilx,ify:ily) + wl*self%buox(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef !no time_dir needed!
                        
                        case ('vy') !horizontal y adjsource
                        f%vy(ifz:ilz,ifx:ilx,ify:ily) = f%vy(ifz:ilz,ifx:ilx,ify:ily) + wl*self%buoy(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef !no time_dir needed!
                        
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                        case ('vz') !vertical z adjsource
                        !vz[ix,iy,iz-0.5]
                        f%vz(iz,ix,iy) = f%vz(iz,ix,iy) + wl*self%buoz(iz,ix,iy) !no time_dir needed!

                        case ('vx') !horizontal x adjsource
                        !vx[ix-0.5,iy,iz]
                        f%vx(iz,ix,iy) = f%vx(iz,ix,iy) + wl*self%buox(iz,ix,iy) !no time_dir needed!
                        
                        case ('vy') !horizontal y adjsource
                        !vy[ix,iy-0.5,iz]
                        f%vy(iz,ix,iy) = f%vy(iz,ix,iy) + wl*self%buoy(iz,ix,iy) !no time_dir needed!
                        
                    end select
                    
                endif
                
            enddo
        
    end subroutine
    
    !forward: v^it -> v^it+1 by FD of s^it+0.5
    !adjoint: v^it+1 -> v^it by FD^T of s^it+0.5
    subroutine update_velocities(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ifz=f%bloom(1,it)+2
        ilz=f%bloom(2,it)-1
        ifx=f%bloom(3,it)+2
        ilx=f%bloom(4,it)-1
        ify=f%bloom(5,it)+2
        ily=f%bloom(6,it)-1

        if(m%is_freesurface) ifz=max(ifz,1)

        if(m%is_cubic) then
            call fd3d_velocities(f%vz,f%vx,f%vy,f%p,                     &
                                 f%dp_dz,f%dp_dx,f%dp_dy,                &
                                 self%buoz,self%buox,self%buoy,          &
                                 ifz,ilz,ifx,ilx,ify,ily,time_dir*self%dt)
        else

            call fd2d_velocities(f%vz,f%vx,f%p,                  &
                                 f%dp_dz,f%dp_dx,                &
                                 self%buoz,self%buox,            &
                                 ifz,ilz,ifx,ilx,time_dir*self%dt)

        endif

        !apply free surface boundary condition if needed
        if(m%is_freesurface) call fd_freesurface_velocities(f%vz)

    end subroutine
    
    !forward: add RHS to s^it+0.5
    !adjoint: add RHS to s^it+1.5
    subroutine inject_stresses(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        if(.not. f%is_adjoint) then

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            ify=shot%src%ify-cb%ioy+1; iy=shot%src%iy-cb%ioy+1; ily=shot%src%ily-cb%ioy+1
            
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler
            
            if(if_hicks) then
                if(shot%src%comp=='p') then
                    f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) + wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef
                endif
                
            else
                if(shot%src%comp=='p') then
                    !explosion on s[iz,ix,iy]
                    f%p(iz,ix,iy) = f%p(iz,ix,iy) + wl*self%kpa(iz,ix,iy)
                endif
                
            endif

            return

        endif

            do i=1,shot%nrcv

                if(shot%rcv(i)%comp=='p') then

                    ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                    ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                    ify=shot%rcv(i)%ify-cb%ioy+1; iy=shot%rcv(i)%iy-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1
                    
                    !adjsource for pressure
                    wl=f%wavelet(i,it)*wavelet_scaler
                    
                    if(if_hicks) then 

                        f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) +wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef !no time_dir needed!

                    else           
                        !p[iz,ix,iy]
                        f%p(iz,ix,iy) = f%p(iz,ix,iy) +wl*self%kpa(iz,ix,iy) !no time_dir needed!

                    endif

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
            call fd3d_stresses(f%vz,f%vx,f%vy,f%p,                     &
                               f%dvz_dz,f%dvx_dx,f%dvy_dy,             &
                               self%kpa,                               &
                               ifz,ilz,ifx,ilx,ify,ily,time_dir*self%dt)
        else
            call fd2d_stresses(f%vz,f%vx,f%p,                  &
                               f%dvz_dz,f%dvx_dx,              &
                               self%kpa,                       &
                               ifz,ilz,ifx,ilx,time_dir*self%dt)
        endif
        
        !apply free surface boundary condition if needed
        if(m%is_freesurface) call fd_freesurface_stresses(f%p)

    end subroutine

    subroutine compute_poynting(v,u,it)
        type(t_field) :: v, u

        real,dimension(:,:,:),allocatable :: E2, ph !envelope squared & inst phase

        !nonzero only when sf touches rf
        ifz=u%bloom(1,it)+2
        ilz=u%bloom(2,it)-2
        ifx=u%bloom(3,it)+2
        ilx=u%bloom(4,it)-2

        ! ify=max(f%bloom(5,it),1)
        ! ily=min(f%bloom(6,it),cb%my)
        
        ! if(m%is_cubic) then
        !     ! call grad3d_moduli(rf%p,sf%vz,sf%vx,sf%vy,&
        !     !                    grad,                  &
        !     !                    ifz,ilz,ifx,ilx,ify,ily)
        ! else
            select case(s_poynting_def)

            case('p_dotp_gradp')! s:=p*dotp*∇p
                call poyn_dotp_gradp(u%p,u%vz,u%vx,ppg%kpa,u%poynz,u%poynx,ifz,ilz,ifx,ilx)
                u%poynz=u%p*u%poynz
                u%poynx=u%p*u%poynx

            case('dotp_gradp')! s:=dotp*∇p
                call poyn_dotp_gradp(u%p,u%vz,u%vx,ppg%kpa,u%poynz,u%poynx,ifz,ilz,ifx,ilx)

            case('p_v')! s:=p*v
                !call poyn_p_v(u%p,u%v,u%poynz,u%poynx,ifz,ilz,ifx,ilx)
                u%poynz=u%p*u%vz
                u%poynx=u%p*u%vx

            case('Esq_gradphi')! s:=E²*∇ϕ
                !E=sqrt(u%p*u%p+v%p*v%p)
                E2=u%p*u%p+v%p*v%p
                ph=atan2(v%p,u%p)

                !$omp parallel default (shared)&
                !$omp private(iz,ix,&
                !$omp         dph_dz,dph_dx)
                !$omp do schedule(dynamic)
                do ix=ifx,ilx
                do iz=ifz,ilz
                    dph_dz = asin(sin(ph(iz+1,ix,1) - ph(iz-1,ix,1)))*inv_2dz
                    dph_dx = asin(sin(ph(iz,ix+1,1) - ph(iz,ix-1,1)))*inv_2dx

                    u%poynz(iz,ix,1)=E2(iz,ix,1)*dph_dz
                    u%poynx(iz,ix,1)=E2(iz,ix,1)*dph_dx

                enddo
                enddo
                !$omp end do
                !$omp end parallel

            end select
            
        ! endif

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
                        case ('vz')
                        f%seismo(i,it)=sum(f%vz(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
                        case ('vx')
                        f%seismo(i,it)=sum(f%vx(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
                        case ('vy')
                        f%seismo(i,it)=sum(f%vy(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                        case ('p') !p[iz,ix,iy]
                        f%seismo(i,it)=f%p(iz,ix,iy)
                        case ('vz') !vz[iz-0.5,ix,iy]
                        f%seismo(i,it)=f%vz(iz,ix,iy)
                        case ('vx') !vx[iz,ix-0.5,iy]
                        f%seismo(i,it)=f%vx(iz,ix,iy)
                        case ('vy') !vy[iz,ix,iy-0.5]
                        f%seismo(i,it)=f%vy(iz,ix,iy)
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

    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !Kₘ<a|Au> = Kₘ<a|M∂ₜu-Du> = ∫ aᵀ KₘM ∂ₜu dt
    !Since it's cumbersome to get ∂ₜu by time marching,
    !replace ∂ₜu by M⁻¹Du and neglect f
    !ie. M∂ₜu=Du+f -> ∂ₜu=M⁻¹Du+M⁻¹f ≐ M⁻¹Du
    !This simplification introduces singularities in the gradient only at source positions,
    !which are probably removed by gradient masking.
    !
    !Therefore, ∫ aᵀ KₘM ∂ₜu dt ≐ ∫ aᵀ KₘM M⁻¹Du dt =: a★Du
    !where
    !                            [ 0   0   ∂ₓᵇ] [vz]
    !a★Du = [vzᵃ vxᵃ pᵃ] Kₘln(M) | 0   0   ∂zᵇ| |vx|
    !                            [∂ₓᶠ  ∂zᶠ  0 ] [p ]
    !In particular, we compute
    !  grho = vᵃ ∂ₜv = vᵃ b∇p
    !  gkpa = pᵃ (-κ⁻²) ∂ₜp = = pᵃ (-κ⁻¹) ∇·v
    !
    !For imaging:
    !I = ∫ a u dt =: a★u

    subroutine cross_correlate_grho(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        if(m%is_cubic) then
            call grad3d_grho(rf%vz,rf%vx,rf%vy,sf%p,&
                             corr%grho,             &
                             ifz,ilz,ifx,ilx,ify,ily)
        else
            call grad2d_grho(rf%vz,rf%vx,sf%p,&
                             corr%grho,       &
                             ifz,ilz,ifx,ilx)
        endif
        
    end subroutine

    subroutine cross_correlate_gkpa(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        if(m%is_cubic) then
            call grad3d_gkpa(rf%p,sf%vz,sf%vx,sf%vy,&
                             corr%gkpa,             &
                             ifz,ilz,ifx,ilx,ify,ily)
        else
            ! call grad2d_moduli(rf%p,sf%p,sf_p_save,&
            !                    grad,            &
            !                    ifz,ilz,ifx,ilx)
            ! sf_p_save = sf%p
            
            !inexact greadient
            call grad2d_gkpa(rf%p,sf%vz,sf%vx,&
                             corr%gkpa,       &
                             ifz,ilz,ifx,ilx)
        endif
        
    end subroutine

    subroutine cross_correlate_image(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr
        
        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        ! if(m%is_cubic) then
        !     call imag3d_xcorr(rf%p,sf%p,&
        !                       imag,                  &
        !                       ifz,ilz,ifx,ilx,ify,ily)
        ! else
            call imag2d(rf%p,sf%p,&
                        rf%poynz,rf%poynx,sf%poynz,sf%poynx, &
                        corr%ipp,corr%ibksc,corr%ifwsc, &
                        ifz,ilz,ifx,ilx)
        ! endif

        ! call imag2d_xcorr(rf%p,rf%vz,rf%vx,&
        !                   sf%p,sf%vz,sf%vx,&
        !                   imag,            &
        !                   ifz,ilz,ifx,ilx)

    end subroutine

    subroutine cross_correlate_postprocess(corr)
        type(t_correlate) :: corr
        
        if(allocated(correlate_gradient)) then
            !scaling gradients by model parameters
            corr%grho = corr%grho / cb%rho(1:cb%mz,1:cb%mx,1:cb%my)
            corr%gkpa = corr%gkpa * (-ppg%inv_kpa(1:cb%mz,1:cb%mx,1:cb%my))
                    
            !preparing for projection back
            corr%gkpa(1,:,:) = corr%gkpa(2,:,:)
            corr%grho(1,:,:) = corr%grho(2,:,:)
        endif

        if(allocated(correlate_image)) then
            corr%ipp (1,:,:) = corr%ipp (2,:,:)
            corr%ibksc(1,:,:) = corr%ibksc(2,:,:)
            corr%ifwsc(1,:,:) = corr%ifwsc(2,:,:)
        endif

    end subroutine

    subroutine energy(sf,it,engy)
        type(t_field),intent(in) :: sf
        real,dimension(cb%mz,cb%mx,cb%my) :: engy

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),cb%my)
        
        if(m%is_cubic) then
            call engy3d_xcorr(sf%p,&
                              engy,                  &
                              ifz,ilz,ifx,ilx,ify,ily)
        else
            call engy2d_xcorr(sf%p,&
                              engy,            &
                              ifz,ilz,ifx,ilx)
        endif

    end subroutine

    !========= Finite-Difference on flattened arrays ==================
    
    subroutine poyn_dotp_gradp(p,vz,vx,kpa,&
                            poynz,poynx,   &
                            ifz,ilz,ifx,ilx)
        real,dimension(*) :: p,vz,vx,kpa,poynz,poynx
        
        nz=cb%nz
        nx=cb%nx
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dp_dz_,dp_dx_)
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
                iz_ixp2=i  +2*nz  !iz,ix+2

                dvz_dz= c1z*(vz(izp1_ix)-vz(iz_ix))  +c2z*(vz(izp2_ix)-vz(izm1_ix))
                dvx_dx= c1x*(vx(iz_ixp1)-vx(iz_ix))  +c2x*(vx(iz_ixp2)-vx(iz_ixm1))

                dp_dz= c1z*(p(iz_ix)-p(izm1_ix)) +c2z*(p(izp1_ix)-p(izm2_ix))
                dp_dx= c1x*(p(iz_ix)-p(iz_ixm1)) +c2x*(p(iz_ixp1)-p(iz_ixm2))

                !s := dotP*nabP = kpa*divv *nabP
                poynz(iz_ix)= kpa(iz_ix)*(dvz_dz+dvx_dx) *dp_dz
                poynx(iz_ix)= kpa(iz_ix)*(dvz_dz+dvx_dx) *dp_dx

            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine

    subroutine fd3d_velocities(vz,vx,vy,p,               &
                               dp_dz,dp_dx,dp_dy,        &
                               buoz,buox,buoy,           &
                               ifz,ilz,ifx,ilx,ify,ily,dt)
        real,dimension(*) :: vz,vx,vy,p
        real,dimension(*) :: dp_dz,dp_dx,dp_dy
        real,dimension(*) :: buoz,buox,buoy
        
        nz=cb%nz
        nx=cb%nx
        ny=cb%ny
        
        dp_dz_=0.;dp_dx_=0.;dp_dy_=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,&
        !$omp         izm2_ix_iy,izm1_ix_iy,iz_ix_iy,izp1_ix_iy,&
        !$omp         iz_ixm2_iy,iz_ixm1_iy,iz_ixp1_iy,&
        !$omp         iz_ix_iym2,iz_ix_iym1,iz_ix_iyp1,&
        !$omp         dp_dz_,dp_dx_,dp_dy_)
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
                
                dp_dz_= c1z*(p(iz_ix_iy)-p(izm1_ix_iy)) +c2z*(p(izp1_ix_iy)-p(izm2_ix_iy))
                dp_dx_= c1x*(p(iz_ix_iy)-p(iz_ixm1_iy)) +c2x*(p(iz_ixp1_iy)-p(iz_ixm2_iy))
                dp_dy_= c1y*(p(iz_ix_iy)-p(iz_ix_iym1)) +c2y*(p(iz_ix_iyp1)-p(iz_ix_iym2))
                
                !cpml
                dp_dz(iz_ix_iy)= cpml%b_z_half(iz)*dp_dz(iz_ix_iy) + cpml%a_z_half(iz)*dp_dz_
                dp_dx(iz_ix_iy)= cpml%b_x_half(ix)*dp_dx(iz_ix_iy) + cpml%a_x_half(ix)*dp_dx_
                dp_dy(iz_ix_iy)= cpml%b_y_half(iy)*dp_dy(iz_ix_iy) + cpml%a_y_half(iy)*dp_dy_

                dp_dz_ = dp_dz_*cpml%kpa_z_half(iz) + dp_dz(iz_ix_iy)
                dp_dx_ = dp_dx_*cpml%kpa_x_half(ix) + dp_dx(iz_ix_iy)
                dp_dy_ = dp_dy_*cpml%kpa_y_half(iy) + dp_dy(iz_ix_iy)
                
                !velocity
                vz(iz_ix_iy)=vz(iz_ix_iy) + dt*buoz(iz_ix_iy)*dp_dz_
                vx(iz_ix_iy)=vx(iz_ix_iy) + dt*buox(iz_ix_iy)*dp_dx_
                vy(iz_ix_iy)=vy(iz_ix_iy) + dt*buoy(iz_ix_iy)*dp_dy_
                
            enddo
            
        enddo
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd2d_velocities(vz,vx,p,          &
                               dp_dz,dp_dx,      &
                               buoz,buox,        &
                               ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vz,vx,p
        real,dimension(*) :: dp_dz,dp_dx
        real,dimension(*) :: buoz,buox
        
        nz=cb%nz
        nx=cb%nx
        
        dp_dz_=0.; dp_dx_=0.

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dp_dz_,dp_dx_)
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

                dp_dz_= c1z*(p(iz_ix)-p(izm1_ix)) +c2z*(p(izp1_ix)-p(izm2_ix))
                dp_dx_= c1x*(p(iz_ix)-p(iz_ixm1)) +c2x*(p(iz_ixp1)-p(iz_ixm2))

                !cpml
                dp_dz(iz_ix)= cpml%b_z_half(iz)*dp_dz(iz_ix) + cpml%a_z_half(iz)*dp_dz_
                dp_dx(iz_ix)= cpml%b_x_half(ix)*dp_dx(iz_ix) + cpml%a_x_half(ix)*dp_dx_

                dp_dz_=dp_dz_*cpml%kpa_z_half(iz) + dp_dz(iz_ix)
                dp_dx_=dp_dx_*cpml%kpa_x_half(ix) + dp_dx(iz_ix)

                !velocity
                vz(iz_ix)=vz(iz_ix) + dt*buoz(iz_ix)*dp_dz_
                vx(iz_ix)=vx(iz_ix) + dt*buox(iz_ix)*dp_dx_

            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd3d_stresses(vz,vx,vy,p,               &
                             dvz_dz,dvx_dx,dvy_dy,     &
                             kpa,                      &
                             ifz,ilz,ifx,ilx,ify,ily,dt)
        real,dimension(*) :: vz,vx,vy,p
        real,dimension(*) :: dvz_dz,dvx_dx,dvy_dy
        real,dimension(*) :: kpa
        
        nz=cb%nz
        nx=cb%nx
        ny=cb%ny
        
        dvz_dz_=0.;dvx_dx_=0.;dvy_dy_=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,&
        !$omp         izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
        !$omp         iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
        !$omp         iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
        !$omp         dvz_dz_,dvx_dx_,dvy_dy_)
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
                
                dvx_dx_= c1x*(vx(iz_ixp1_iy)-vx(iz_ix_iy))  +c2x*(vx(iz_ixp2_iy)-vx(iz_ixm1_iy))
                dvy_dy_= c1y*(vy(iz_ix_iyp1)-vy(iz_ix_iy))  +c2y*(vy(iz_ix_iyp2)-vy(iz_ix_iym1))
                dvz_dz_= c1z*(vz(izp1_ix_iy)-vz(iz_ix_iy))  +c2z*(vz(izp2_ix_iy)-vz(izm1_ix_iy))
                
                !cpml
                dvz_dz(iz_ix_iy)=cpml%b_z(iz)*dvz_dz(iz_ix_iy)+cpml%a_z(iz)*dvz_dz_
                dvx_dx(iz_ix_iy)=cpml%b_x(ix)*dvx_dx(iz_ix_iy)+cpml%a_x(ix)*dvx_dx_
                dvy_dy(iz_ix_iy)=cpml%b_y(iy)*dvy_dy(iz_ix_iy)+cpml%a_y(iy)*dvy_dy_

                dvz_dz_=dvz_dz_*cpml%kpa_z(iz) + dvz_dz(iz_ix_iy)
                dvx_dx_=dvx_dx_*cpml%kpa_x(ix) + dvx_dx(iz_ix_iy)
                dvy_dy_=dvy_dy_*cpml%kpa_y(iy) + dvy_dy(iz_ix_iy)
                
                !pressure
                p(iz_ix_iy) = p(iz_ix_iy) + dt * kpa(iz_ix_iy)*(dvz_dz_+dvx_dx_+dvy_dy_)
                
            enddo
            
        enddo
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine
    
    subroutine fd2d_stresses(vz,vx,p,          &
                             dvz_dz,dvx_dx,    &
                             kpa,              &
                             ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vz,vx,p
        real,dimension(*) :: dvz_dz,dvx_dx
        real,dimension(*) :: kpa
        
        nz=cb%nz
        nx=cb%nx
        
        dvz_dz_=0.;dvx_dx_=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dvz_dz_,dvx_dx_)
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
                
                dvz_dz_= c1z*(vz(izp1_ix)-vz(iz_ix))  +c2z*(vz(izp2_ix)-vz(izm1_ix))
                dvx_dx_= c1x*(vx(iz_ixp1)-vx(iz_ix))  +c2x*(vx(iz_ixp2)-vx(iz_ixm1))
                
                !cpml
                dvz_dz(iz_ix)=cpml%b_z(iz)*dvz_dz(iz_ix)+cpml%a_z(iz)*dvz_dz_
                dvx_dx(iz_ix)=cpml%b_x(ix)*dvx_dx(iz_ix)+cpml%a_x(ix)*dvx_dx_

                dvz_dz_=dvz_dz_*cpml%kpa_z(iz) + dvz_dz(iz_ix)
                dvx_dx_=dvx_dx_*cpml%kpa_x(ix) + dvx_dx(iz_ix)
                
                !pressure
                p(iz_ix) = p(iz_ix) + dt * kpa(iz_ix)*(dvz_dz_+dvx_dx_)
                
            enddo
            
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine


    subroutine fd_freesurface_velocities(vz)
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: vz

        !free surface is located at [1,ix,iy] level
        !so symmetric mirroring: vz[0.5]=vz[1.5], ie. vz(1,ix,iy)=vz(2,ix,iy) -> dp(1,ix,iy)=0.
            vz(1,:,:)=vz(2,:,:)
            ! !$omp parallel default (shared)&
            ! !$omp private(ix,iy,i)
            ! !$omp do schedule(dynamic)
            ! do iy=ify,ily
            ! do ix=ifx,ilx
            !     i=(1-cb%ifz) + (ix-cb%ifx)*nz + (iy-cb%ify)*nz*nx +1 !iz=1,ix,iy
                
            !     f%vz(i)=f%vz(i+1)
            ! enddo
            ! enddo
            ! !$omp enddo
            ! !$omp end parallel

    end subroutine

    subroutine fd_freesurface_stresses(p)
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: p

        !free surface is located at [1,ix,iy] level
        !so explicit boundary condition: p(1,ix,iy)=0
        !and antisymmetric mirroring: p(0,ix,iy)=-p(2,ix,iy) -> vz(2,ix,iy)=vz(1,ix,iy)
            p(1,:,:)=0.
            p(0,:,:)=-p(2,:,:)
            ! !$omp parallel default (shared)&
            ! !$omp private(ix,iy,i)
            ! !$omp do schedule(dynamic)
            ! do iy=ify,ily
            ! do ix=ifx,ilx
            !     i=(1-cb%ifz) + (ix-cb%ifx)*nz + (iy-cb%ify)*nz*nx +1 !iz=1,ix,iy 
                
            !     f%p(i)=0.
                
            !     f%p(i-1)=-f%p(i+1)
            ! enddo
            ! enddo
            ! !$omp enddo
            ! !$omp end parallel

    end subroutine

    
    subroutine grad3d_gkpa(rf_p,sf_vz,sf_vx,sf_vy,&
                           grad,                  &
                           ifz,ilz,ifx,ilx,ify,ily)
        real,dimension(*) :: rf_p,sf_vz,sf_vx,sf_vy
        real,dimension(*) :: grad
        
        nz=cb%nz
        nx=cb%nx
        
        dvz_dz=0.
        dvx_dx=0.
        dvy_dx=0.
        
         rp=0.
        dsp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,j,&
        !$omp         izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
        !$omp         iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
        !$omp         iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
        !$omp         dvz_dz,dvx_dx,dvy_dy,&
        !$omp         rp,dsp)
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
                
                dvz_dz = c1z*(sf_vz(izp1_ix_iy)-sf_vz(iz_ix_iy)) +c2z*(sf_vz(izp2_ix_iy)-sf_vz(izm1_ix_iy))
                dvx_dx = c1x*(sf_vx(iz_ixp1_iy)-sf_vx(iz_ix_iy)) +c2x*(sf_vx(iz_ixp2_iy)-sf_vx(iz_ixm1_iy))
                dvy_dy = c1y*(sf_vy(iz_ix_iyp1)-sf_vy(iz_ix_iy)) +c2y*(sf_vy(iz_ix_iyp2)-sf_vy(iz_ix_iym1))
                
                 rp = rf_p(i)
                dsp = dvz_dz +dvx_dx +dvy_dy
                
                grad(j)=grad(j) + rp*dsp
                
            end do
            
        end do
        end do
        !$omp end do
        !$omp end parallel
        
    end subroutine

    ! subroutine grad2d_moduli(rf_p,sf_p,sf_p_save,&
    !                          grad,            &
    !                          ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: rf_p,sf_p,sf_p_save
    !     real,dimension(*) :: grad
        
    !     nz=cb%nz
        
    !      rp=0.
    !     dsp=0.
        
    !     !$omp parallel default (shared)&
    !     !$omp private(iz,ix,i,j,&
    !     !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
    !     !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
    !     !$omp         dvz_dz,dvx_dx,&
    !     !$omp         dsp,rp)
    !     !$omp do schedule(dynamic)
    !     do ix=ifx,ilx
        
    !         !dir$ simd
    !         do iz=ifz,ilz
                
    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
    !             j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
    !              rp = rf_p(i)
    !             dsp = sf_p_save(i) - sf_p(i)
                
    !             grad(j)=grad(j) + rp*dsp
                
    !         end do
            
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

    subroutine grad2d_gkpa(rf_p,sf_vz,sf_vx,&
                           grad,            &
                           ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_p,sf_vz,sf_vx
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        dvz_dz=0.
        dvx_dx=0.
        
         rp=0.
        dsp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dvz_dz,dvx_dx,&
        !$omp         rp,dsp)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers
                
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm1=i    -nz  !iz,ix-1
                iz_ixp1=i    +nz  !iz,ix+1
                iz_ixp2=i  +2*nz  !iz,ix+2
                
                dvz_dz = c1z*(sf_vz(izp1_ix)-sf_vz(iz_ix)) +c2z*(sf_vz(izp2_ix)-sf_vz(izm1_ix))
                dvx_dx = c1x*(sf_vx(iz_ixp1)-sf_vx(iz_ix)) +c2x*(sf_vx(iz_ixp2)-sf_vx(iz_ixm1))
                
                 rp = rf_p(i)
                dsp = dvz_dz +dvx_dx
                
                grad(j)=grad(j) + rp*dsp
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine
    
    subroutine grad3d_grho(rf_vz,rf_vx,rf_vy,sf_p,&
                           grad,                  &
                           ifz,ilz,ifx,ilx,ify,ily)
        real,dimension(*) :: rf_vz,rf_vx,rf_vy,sf_p
        real,dimension(*) :: grad
        
        nz=cb%nz
        nx=cb%nx
        
        dsvz=0.; dsvx=0.; dsvy=0.
         rvz=0.;  rvx=0.;  rvy=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,j,&
        !$omp         izm2_ix_iy,izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
        !$omp         iz_ixm2_iy,iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
        !$omp         iz_ix_iym2,iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
        !$omp          rvz, rvx, rvy,&
        !$omp         dsvz,dsvx,dsvy)
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
                
                rvz = rf_vz(izp1_ix_iy) +rf_vz(iz_ix_iy)
                rvx = rf_vx(iz_ixp1_iy) +rf_vx(iz_ix_iy)
                rvy = rf_vy(iz_ix_iyp1) +rf_vy(iz_ix_iy)

                dsvz = (c1z*(sf_p(iz_ix_iy  )-sf_p(izm1_ix_iy)) +c2z*(sf_p(izp1_ix_iy)-sf_p(izm2_ix_iy))) &
                      +(c1z*(sf_p(izp1_ix_iy)-sf_p(iz_ix_iy  )) +c2z*(sf_p(izp2_ix_iy)-sf_p(izm1_ix_iy)))
                dsvx = (c1x*(sf_p(iz_ix_iy  )-sf_p(iz_ixm1_iy)) +c2x*(sf_p(iz_ixp1_iy)-sf_p(iz_ixm2_iy))) &
                      +(c1x*(sf_p(iz_ixp1_iy)-sf_p(iz_ix_iy  )) +c2x*(sf_p(iz_ixp2_iy)-sf_p(iz_ixm1_iy)))
                dsvy = (c1y*(sf_p(iz_ix_iy  )-sf_p(iz_ix_iym1)) +c2y*(sf_p(iz_ix_iyp1)-sf_p(iz_ix_iym2))) &
                      +(c1y*(sf_p(iz_ix_iyp1)-sf_p(iz_ix_iy  )) +c2y*(sf_p(iz_ix_iyp2)-sf_p(iz_ix_iym1)))
                
                !complete equation with unnecessary terms e.g. sf_p(iz_ix) for better understanding
                !with flag -O, the compiler should automatically detect such possibilities of simplification
                
                grad(j)=grad(j) + 0.25*( rvz*dsvz + rvx*dsvx + rvy*dsvy )
                
            enddo
            
        enddo
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine grad2d_grho(rf_vz,rf_vx,sf_p,&
                           grad,            &
                           ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_vz,rf_vx,sf_p
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        dsvz=0.; dsvx=0.
         rvz=0.; rvx=0.

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         rvz, rvx,&
        !$omp         dsvz,dsvx)
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
                
                rvz = rf_vz(iz_ix) +rf_vz(izp1_ix)
                rvx = rf_vx(iz_ix) +rf_vx(iz_ixp1)
                
                dsvz = (c1z*(sf_p(iz_ix  )-sf_p(izm1_ix)) +c2z*(sf_p(izp1_ix)-sf_p(izm2_ix))) &
                      +(c1z*(sf_p(izp1_ix)-sf_p(iz_ix  )) +c2z*(sf_p(izp2_ix)-sf_p(izm1_ix)))
                dsvx = (c1x*(sf_p(iz_ix  )-sf_p(iz_ixm1)) +c2x*(sf_p(iz_ixp1)-sf_p(iz_ixm2))) &
                      +(c1x*(sf_p(iz_ixp1)-sf_p(iz_ix  )) +c2x*(sf_p(iz_ixp2)-sf_p(iz_ixm1)))
                !complete equation with unnecessary terms e.g. sf_p(iz_ix) for better understanding
                !with flag -Ox, the compiler should automatically detect such possible simplification
                
                grad(j)=grad(j) + 0.25*( rvz*dsvz + rvx*dsvx )
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine


    ! subroutine imag3d_ipp(rf_p,sf_p,             &
    !                       imag,                  &
    !                       ifz,ilz,ifx,ilx,ify,ily)
    !     real,dimension(*) :: rf_p,sf_p
    !     real,dimension(*) :: imag
        
    !     nz=cb%nz
    !     nx=cb%nx
        
    !     rp=0.
    !     sp=0.
        
    !     !$omp parallel default (shared)&
    !     !$omp private(iz,ix,iy,i,j,&
    !     !$omp         rp,sp)
    !     !$omp do schedule(dynamic) collapse(2)
    !     do iy=ify,ily
    !     do ix=ifx,ilx
        
    !         !dir$ simd
    !         do iz=ifz,ilz
                
    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(iy-cb%ify)*cb%nz*cb%nx+1 !field has boundary layers
    !             j=(iz-1)     +(ix-1)     *cb%mz+(iy-1)     *cb%mz*cb%mx+1 !grad has no boundary layers
                
    !             rp = rf_p(i)
    !             sp = sf_p(i)
                
    !             imag(j)=imag(j) + rp*sp
                
    !         end do
            
    !     end do
    !     end do
    !     !$omp end do
    !     !$omp end parallel
        
    ! end subroutine

    ! subroutine imag2d_xcorr(rf_p,rf_vz,rf_vx,&
    !                         sf_p,sf_vz,sf_vx,&
    !                         imag,          &
    !                         ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: rf_p,rf_vz,rf_vx
    !     real,dimension(*) :: sf_p,sf_vz,sf_vx
    !     real,dimension(*) :: imag
        
    !     nz=cb%nz
        
    !     !$omp parallel default (shared)&
    !     !$omp private(iz,ix,i,j)
    !     !$omp do schedule(dynamic)
    !     do ix=ifx,ilx
        
    !         !dir$ simd
    !         do iz=ifz,ilz
                
    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
    !             j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
    !             imag(j)=imag(j) + &
    !                 rf_p(i)*sf_p(i) + rf_vz(i)*sf_vz(i) + rf_vx(i)*sf_vx(i)
                
    !         end do
            
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

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

    subroutine engy3d_xcorr(sf_p,          &
                            engy,          &
                            ifz,ilz,ifx,ilx,ify,ily)
        real,dimension(*) :: sf_p
        real,dimension(*) :: engy
        
        nz=cb%nz
        nx=cb%nx
        
        sp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,j,&
        !$omp         sp)
        !$omp do schedule(dynamic) collapse(2)
        do iy=ify,ily
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(iy-cb%ify)*cb%nz*cb%nx+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+(iy-1)     *cb%mz*cb%mx+1 !grad has no boundary layers
                
                sp = sf_p(i)
                
                engy(j)=engy(j) + sp*sp
                
            end do
            
        end do
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine engy2d_xcorr(sf_p,          &
                            engy,          &
                            ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_p
        real,dimension(*) :: engy
        
        nz=cb%nz
        
        sp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         sp)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers
                
                sp = sf_p(i)
                
                engy(j)=engy(j) + sp*sp
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

end
