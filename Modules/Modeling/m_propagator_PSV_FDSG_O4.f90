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

use, intrinsic :: ieee_arithmetic

    private

    !FD coef
    real,dimension(2),parameter :: coef = [9./8.,-1./24.] !Fornberg, 1988, Generation of Finite-Difference Formulas on Arbitrary Spaced Grids.
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z

    !local const
    real :: dt2, inv_2dt, inv_2dz, inv_2dx

    !scaling source wavelet
    real :: wavelet_scaler

    character(:),allocatable :: FS_method

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain ISOtropic 2D P-SV (ELastic) propagation'//s_NL// &
            '1st-order Velocity-Stress formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.606 *Vmax/dx'//s_NL// &
            'Required model attributes: vp, vs, rho'//s_NL// &
            'Required field components: vz, vx'//s_NL// &
            'Required field components: szz, sxx, szx'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Imaging conditions: P-Pxcorr'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: grho glda gmu'

        integer :: nbndlayer=max(2,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=3 !number of basic gradients
        integer :: nimag=1 !number of basic images
        integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:),allocatable :: buoz, buox, ldap2mu, lda, mu
        real,dimension(:,:),allocatable :: ldapmu, two_ldapmu, inv_ldapmu_4mu

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
        
        procedure :: inject_velocities
        procedure :: inject_stresses
        procedure :: update_velocities
        procedure :: update_stresses
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    !conversion from normal stresses to pressure
    real,parameter :: sn_p=0.5  !2D

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    logical :: if_record_adjseismo=.false.

    contains
    
    !========= for FDSG O(dx4,dt2) ===================  

    subroutine print_info(self)
        class(t_propagator) :: self

        call hud('Invoked field & propagator modules info : '//s_NL//self%info)
        call hud('FDGS Coef : '//num2str(coef(1))//', '//num2str(coef(2)))
        
    end subroutine
    
    subroutine estim_RAM(self)
        class(t_propagator) :: self
    end subroutine
    
    subroutine check_model(self)
        class(t_propagator) :: self
        
        if(index(self%info,'vp')>0  .and. .not. allocated(m%vp)) then
            !call error('vp model is NOT given.')
            call alloc(m%vp,m%nz,m%nx,1,o_init=1500.)
            call warn('Constant vp model (1500 m/s) is allocated by propagator.')
        endif

        if(index(self%info,'vs')>0  .and. .not. allocated(m%vs)) then
            call alloc(m%vs,m%nz,m%nx,1)
            m%vs=m%vp/sqrt(3.)
            call warn('Poisson solid (vs=vp/√3) is assumed by propagator.')
        endif

        if(index(self%info,'rho')>0 .and. .not. allocated(m%rho)) then
            call alloc(m%rho,m%nz,m%nx,1,o_init=1000.)
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

        CFL = sumcoef*cb%velmax*self%dt*m%rev_cell_diagonal

        call hud('CFL value: '//num2str(CFL))
        
        if(CFL>1.) then
            self%dt = setup%get_real('CFL',o_default='0.9')/(sumcoef*cb%velmax*m%rev_cell_diagonal)
            self%nt=nint(time_window/self%dt)+1

            call warn('CFL > 1 on '//shot%sindex//'!'//s_NL//&
                'vmax, dt, 1/dx = '//num2str(cb%velmax)//', '//num2str(self%dt)//', '//num2str(m%rev_cell_diagonal) //s_NL//&
                'Adjusted dt, nt = '//num2str(self%dt)//', '//num2str(self%nt))

        endif       
        
    end subroutine

    subroutine init(self,oif_record_adjseismo)
        class(t_propagator) :: self

        logical,optional :: oif_record_adjseismo

        real,dimension(:,:),allocatable :: temp_mu

        c1x=coef(1)/m%dx; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2z=coef(2)/m%dz
        
        wavelet_scaler=self%dt/m%cell_volume

        if_hicks=shot%if_hicks

        if_record_adjseismo=either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))

        FS_method=setup%get_str('FS_METHOD',o_default='stress_image')

        call alloc(self%buoz,   [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%buox,   [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%ldap2mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%lda,    [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%mu,     [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%inv_ldapmu_4mu, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%ldapmu, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%two_ldapmu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])



        call alloc(temp_mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        self%ldap2mu(:,:)=cb%rho(:,:,1)*cb%vp(:,:,1)**2
             temp_mu(:,:)=cb%rho(:,:,1)*cb%vs(:,:,1)**2

        self%lda=self%ldap2mu-2.*temp_mu
        ! if(mpiworld%is_master) then
        ! write(*,*) 'self%ldap2mu sanity:', minval(self%ldap2mu),maxval(self%ldap2mu)
        ! write(*,*) 'self%lda     sanity:', minval(self%lda),maxval(self%lda)
        ! endif

        self%ldapmu=self%lda+temp_mu
        self%two_ldapmu=2.*self%ldapmu

! print*,'!!!!!!!!!!!!!!!!!!!!!!!'
        self%inv_ldapmu_4mu=0.25/(self%lda+temp_mu)/temp_mu
! print*,'!!!!!!!!!!!!!!!!!!!!!!!'
        !interpolat mu by harmonic average
        temp_mu=1./temp_mu

        do ix=cb%ifx+1,cb%ilx
        do iz=cb%ifz+1,cb%ilz
                self%mu(iz,ix)=4./(temp_mu(iz-1,ix-1) &
                                  +temp_mu(iz-1,ix  ) &
                                  +temp_mu(iz  ,ix-1) &
                                  +temp_mu(iz  ,ix  ))
        end do
        end do

        where( ieee_is_nan(self%mu) .or. .not.ieee_is_finite(self%mu) )
            self%mu=0.
        endwhere

        ! open(8,file='self%mu',access='stream')
        ! write(8) self%mu
        ! close(8)

        ! !interpolat mu by arithmic average
        ! do ix=cb%ifx+1,cb%ilx
        ! do iz=cb%ifz+1,cb%ilz
        !         self%mu(iz,ix)=( temp_mu(iz-1,ix-1) &
        !                       +temp_mu(iz-1,ix  ) &
        !                       +temp_mu(iz  ,ix-1) &
        !                       +temp_mu(iz  ,ix  ))
        ! end do
        ! end do
        ! self%mu=self%mu*0.25

        deallocate(temp_mu)


        self%mu(cb%ifz,:)=self%mu(cb%ifz+1,:)
        self%mu(:,cb%ifx)=self%mu(:,cb%ifx+1)

        !check mu values
        if(mpiworld%is_master) then
            write(*,*) 'ppg%mu sanity:', minval(self%mu),maxval(self%mu), any(ieee_is_nan(self%mu)), any(.not. ieee_is_finite(self%mu))
        endif


        self%buoz(cb%ifz,:)=1./cb%rho(cb%ifz,:,1)
        self%buox(:,cb%ifx)=1./cb%rho(:,cb%ifx,1)

        do iz=cb%ifz+1,cb%ilz
            self%buoz(iz,:)=0.5/cb%rho(iz,:,1)+0.5/cb%rho(iz-1,:,1)
        enddo
        
        do ix=cb%ifx+1,cb%ilx
            self%buox(:,ix)=0.5/cb%rho(:,ix,1)+0.5/cb%rho(:,ix-1,1)
        enddo



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
        call f%init_boundary_velocities

        call alloc(f%vz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%vx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%szz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%sxx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%szx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])

        call alloc(f%dvz_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dvz_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dvx_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dvx_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dszz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dsxx_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dszx_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dszx_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
                
    end subroutine

    subroutine init_correlate(self,corr,name)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name
        
        corr%name=name

        ! if(name(1:1)=='g') then !gradient components
            call alloc(corr%grho,m%nz,m%nx,m%ny)
            call alloc(corr%glda,m%nz,m%nx,m%ny)
            call alloc(corr%gmu, m%nz,m%nx,m%ny)
        ! else !image components
        !     call alloc(corr%ipp,m%nz,m%nx,m%ny)
        !     call alloc(corr%ibksc,m%nz,m%nx,m%ny)
        !     call alloc(corr%ifwsc,m%nz,m%nx,m%ny)
        ! endif

    end subroutine

    subroutine init_abslayer(self)
        class(t_propagator) :: self

        call cpml%init

    end subroutine

    subroutine assemble(self,corr)
        class(t_propagator) :: self
        type(t_correlate) :: corr

        if(allocated(correlate_gradient)) then
            call correlate_assemble(corr%grho, correlate_gradient(:,:,:,1))
            call correlate_assemble(corr%glda, correlate_gradient(:,:,:,2))
            call correlate_assemble(corr%gmu,  correlate_gradient(:,:,:,3))
        endif        
        
    end subroutine
    
    !========= Derivations =================
    !PDE:      A u = M ∂ₜ u - D u = f
    !Adjoint:  Aᵀa = M ∂ₜᵀa - Dᵀa = d
    !where
    !u=[vz vx szz sxx szx]ᵀ, [vz vx] are velocities, [szz sxx] are normal stresses, szx is shear stress
    !f=[fz fx fzz fxx 0]ᵀδ(x-xs) with xs source position, d is recorded data
    !p=tr(s)=½(szz+sxx) is (hydrostatic) pressure, fzz=fxx=½fp is pressure source (explosion)
    !  [ρ                  ]    [0  0  ∂z 0  ∂ₓ]
    !  |  ρ                |    |0  0  0  ∂ₓ ∂z|
    !M=|    [λ+2μ  λ    ]⁻¹|, D=|∂z 0  0  0  0 |
    !  |    | λ   λ+2μ  |  |    |0  ∂ₓ 0  0  0 |
    !  [    [          μ]  ]    [∂ₓ ∂z 0  0  0 ]
    !a=[vzᵃ vxᵃ szzᵃ, sxxᵃ, szxᵃ]ᵀ is the adjoint field
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
    !                      |    ss    |   ss -½ vz  ss    |         |
    !                      |    μ     |   μ     bz  μ     |         |
    !                      |          |         |         |         |
    !                     λ,μ   bx   λ,μ  bx   λ,μ  bx   λ,μ  bx   λ,μ
    !  -v--s-v-s-v-→ t    -sn---vx---sn---vx---sn---vx---sn---vx---sn-→ x
    !  -1 -½ 0 ½ 1        -2   -1½   -1   -½    0    ½    1   1½    2    
    !                      |          |         |         |         | 
    !                      |    ss    |   ss  ½ vz  ss    |         | 
    !                      |    μ     |   μ     bz  μ     |         | 
    !                      |          |         |         |         | 
    !                     -|----------|-------1-sn--------|---------|-
    !                      |          |         κ         |         | 
    !                      |          |         |         |         | 
    !                      |          |      1½ vz        |         | 
    !                      |          |         bz        |         | 
    !                      |          |         |         |         | 
    !                                         z ↓
    !
    !  s: szz, sxx or szx
    ! sn: szz or sxx (normal stresses)
    ! ss: szx (shear stresses)
    ! λ,μ: λ and λ+2μ
    !
    !Convention for half-integer index:
    !(array index)  =>     (real index)     
    !  vz(iz,ix)    =>   vz[iz-½,ix  ]^n   :=vz((iz-½)*dz,ix*dx,n*dt)
    !  vx(iz,ix)    =>   vx[iz,  ix-½]^n  
    !  sn(iz,ix)    =>   sn[iz,  ix  ]^n+½ :=sn(iz*dz,    ix*dx,(n+½)*dt)
    !  ss(iz,ix)    =>   ss[iz-½,ix-½]^n  
    !
    !Forward:
    !FD eqn:
    !      [ vz^n  ]   [ 0     0    ∂zᵇ  0  ∂ₓᶠ][ vz^n+1]      [∂zᵇ szz^n+½ + ∂ₓᶠ szx^n+½]  
    !      | vx^n  |   | 0     0     0  ∂ₓᵇ ∂zᶠ|| vx^n+1|      |∂ₓᵇ sxx^n+½ + ∂zᶠ szx^n+½|  
    !M ∂ₜᶠ |szz^n+1| = |∂zᶠ    0     0   0   0 ||szz^n+½| +f = |∂zᶠ  vz^n+1              | +f
    !      |sxx^n+1|   | 0    ∂ₓᶠ    0   0   0 ||sxx^n+½|      |∂ₓᶠ  vx^n+1              |  
    !      [szx^n+1]   [∂ₓᵇ   ∂zᵇ    0   0   0 ][szx^n+½]      [∂ₓᵇ  vz^n+1 + ∂zᵇ  vx^n+1]  
    !where
    !∂ₜᶠ*dt := v^n+1 - v^n                             ~O(t²)
    !∂zᵇ*dz := c₁(s(iz  )-s(iz-1) +c₂(s(iz+1)-s(iz-2)  ~O(x⁴)
    !∂zᶠ*dz := c₁(v(iz+1)-v(iz  ) +c₂(v(iz+2)-v(iz-1)  ~O(x⁴)
    !
    !Time marching:
    ![ vz^n+1 ]   [ vz^n  ]      [∂zᵇ szz^n+½ + ∂ₓᶠ szx^n+½]
    !| vx^n+1 |   | vx^n  |      |∂ₓᵇ sxx^n+½ + ∂zᶠ szx^n+½|
    !|szz^n+1½| = |szz^n+½| + M⁻¹|∂zᶠ  vz^n+1              |dt  +M⁻¹f*dt
    !|sxx^n+1½|   |sxx^n+½|      |∂ₓᶠ  vx^n+1              |
    ![szx^n+1½]   [szx^n+½]      [∂ₓᵇ  vz^n+1 + ∂zᵇ  vx^n+1]
    !where M⁻¹f=[b*fz b*fx 2(λ+μ)fp 2(λ+μ)fp 0]
    !Step #1: v^n += src
    !Step #2: v^n+1 = v^n + spatial FD(p^n+½)
    !Step #3: s^n+½ += src
    !Step #4: s^n+1½ = s^n+½ + spatial FD(v^n+1)
    !Step #5: sample v^n & s^n++½ at receivers
    !Step #6: save v^n+1 to boundary values
    !
    !Reverse time marching (for wavefield reconstruction)
    ![szx^n+½]   [szx^n+1½]      [∂ₓᵇ  vz^n+1 + ∂zᵇ  vx^n+1]
    !|sxx^n+½|   |sxx^n+1½|      |∂ₓᶠ  vx^n+1              |
    !|szz^n+½| = |szz^n+1½| - M⁻¹|∂zᶠ  vz^n+1              |dt  -M⁻¹f*dt
    !| vx^n  |   | vx^n+1 |      |∂ₓᵇ sxx^n+½ + ∂zᶠ szx^n+½|
    ![ vz^n  ]   [ vz^n+1 ]      [∂zᵇ szz^n+½ + ∂ₓᶠ szx^n+½]
    !Step #6: load boundary values for v^n+1
    !Step #4: s^n+½ = s^n+1½ - spatial FD(v^n+1)
    !Step #3: s^n+½ -= src
    !Step #2: v^n+1 = v^n - spatial FD(s^n+½)
    !Step #1: v^n -= src
    !N.B. Same codes for spatial FDs as in forward time marching, with a negated dt.
    !
    !Adjoint:
    !FD eqn:
    !      [ vzᵃ^n  ]   [ 0     0    ∂zᶠᵀ  0   ∂ₓᵇᵀ][ vzᵃ^n+1]
    !      | vxᵃ^n  |   | 0     0     0   ∂ₓᶠᵀ ∂zᵇᵀ|| vxᵃ^n+1|
    !M ∂ₜᶠᵀ|szzᵃ^n+1| = |∂zᵇᵀ   0     0    0    0  ||szzᵃ^n+½| +d
    !      |sxxᵃ^n+1|   | 0    ∂ₓᵇᵀ   0    0    0  ||sxxᵃ^n+½|
    !      [szxᵃ^n+1]   [∂ₓᶠᵀ  ∂zᶠᵀ   0    0    0  ][szxᵃ^n+½]
    !∂ₜᶠᵀ = v^n-1 -v^n   = -∂ₜᵇ
    !∂zᵇᵀ = c₁(v[i  ]-v[i+½]) +c₂(v[i- ½]-v[i+1½]) = -∂zᶠ
    !∂zᶠᵀ = c₁(s[i-½]-s[i  ]) +c₂(s[i-1½]-s[i+ ½]) = -∂zᵇ
    !      [ vzᵃ^n  ]   [ 0    0   ∂zᵇ  0  ∂ₓᶠ][ vzᵃ^n+1]
    !      | vxᵃ^n  |   | 0    0    0  ∂ₓᵇ ∂zᶠ|| vxᵃ^n+1|
    !M ∂ₜᵇ |szzᵃ^n+1| = |∂zᶠ   0    0   0   0 ||szzᵃ^n+½| -d
    !      |sxxᵃ^n+1|   | 0   ∂ₓᶠ   0   0   0 ||sxxᵃ^n+½|
    !      [szxᵃ^n+1]   [∂ₓᵇ  ∂zᵇ   0   0   0 ][szxᵃ^n+½]
    !ie. Dᵀ=-D, antisymmetric
    !
    !Time marching:
    ![ vzᵃ^n+1 ]   [ vzᵃ^n  ]      [∂zᵇ szzᵃ^n+½ + ∂ₓᵇ szxᵃ^n+½]
    !| vxᵃ^n+1 |   | vxᵃ^n  |      |∂ₓᵇ sxxᵃ^n+½ + ∂zᵇ szxᵃ^n+½|
    !|szzᵃ^n+1½| = |szzᵃ^n+½| + M⁻¹|∂zᶠ  vzᵃ^n+1               |dt  -M⁻¹d*dt
    !|sxxᵃ^n+1½|   |sxxᵃ^n+½|      |∂ₓᶠ  vxᵃ^n+1               |
    ![szxᵃ^n+1½]   [szxᵃ^n+½]      [∂ₓᵇ  vzᵃ^n+1 + ∂zᵇ  vxᵃ^n+1]
    !but we have to do it in reverse time:
    ![szxᵃ^n+½] = [szxᵃ^n+1½]      [∂ₓᵇ  vzᵃ^n+1 + ∂zᵇ  vxᵃ^n+1]
    !|sxxᵃ^n+½| = |sxxᵃ^n+1½|      |∂ₓᶠ  vxᵃ^n+1               |
    !|szzᵃ^n+½| = |szzᵃ^n+1½| - M⁻¹|∂zᶠ  vzᵃ^n+1               |dt  +M⁻¹d*dt
    !| vxᵃ^n  | = | vxᵃ^n+1 |      |∂ₓᵇ sxxᵃ^n+½ + ∂zᵇ szxᵃ^n+½|
    ![ vzᵃ^n  ] = [ vzᵃ^n+1 ]      [∂zᵇ szzᵃ^n+½ + ∂ₓᵇ szxᵃ^n+½]
    !Step #5: sᵃ^n+1½ += adjsrc
    !Step #4: sᵃ^n+½ = sᵃ^n+1½ - spatial FD(vᵃ^n+1)
    !Step #3: vᵃ^n+1 += adjsrc
    !Step #2: vᵃ^n = vᵃ^n+1 - spatial FD(pᵃ^n+½)
    !N.B. Same codes for spatial FDs as in forward time marching, with a negated dt, but the RHS should use a "+" sign (regardless of the reverse time direction).
    !
    !For adjoint test:
    !Forward time marching: dsyn=RGNf
    !f:source wavelet, N=M⁻¹dt: diagonal, G:propagator, R:extract field at receivers
    !Reverse-time adjoint marching: Ndadj=GᵀRᵀdsyn
    !Rᵀ:inject adjoint sources, Gᵀ:adjoint propagator

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
                call fld_u%check_value(fld_u%vz)
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
            tt5=tt5+toc-tic

            !snapshot
            call fld_u%write(it)

            !step 6: save v^it+1 in boundary layers
            ! if(fld_u%if_will_reconstruct) then
                call cpu_time(tic)
                call fld_u%boundary_transport_velocities('save',it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            ! endif

        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add source velocities',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update velocities    ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to add source stresses  ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update stresses      ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field        ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary        ',tt6/mpiworld%max_threads
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
        
        !for adjoint test
        if(if_record_adjseismo)  call alloc(fld_a%seismo,1,self%nt)

        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.
        
        ift=1; ilt=self%nt

        ! call alloc(sf_p_save,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%1,cb%1])
        
        do it=ilt,ift,int(time_dir)
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_a%check_value(fld_a%vz)
                call fld_u%check_value(fld_u%vz)
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then

                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_u%boundary_transport_velocities('load',it)
                call cpu_time(toc)
                tt1=tt1+toc-tic
                
                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%update_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_stresses(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_stresses(fld_a,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: rf%s^it+0.5 star D sf%s_dt^it+0.5
            !use sf%v^it+1 to compute sf%s_dt^it+0.5, as backward step 4
            if(mod(it,irdt)==0) then
                call cpu_time(tic)
                call cross_correlate_glda_gmu(fld_a,fld_u,a_star_u,it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif
                            
            !========================================================!

            ! if(present(o_sf)) then
                !backward step 2: v^it+1 -> v^it by FD of s^it+0.5
                call cpu_time(tic)
                call self%update_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt7=tt7+toc-tic

                !backward step 1: rm forces from v^it
                call cpu_time(tic)
                call self%inject_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt8=tt8+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 3: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_velocities(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities(fld_a,time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            !adjoint step 1: sample v^it or s^it+0.5 at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_a,it)
                call cpu_time(toc)
                tt11=tt11+toc-tic
            endif
            
!            !grho: sfield%v_dt^it \dot rfield%v^it
!            !use sfield%s^it+0.5 to compute sfield%v_dt^it, as backward step 2
!            if(if_compute_grad.and.mod(it,irdt)==0) then
!                call cpu_time(tic)
!                call gradient_density(fld_a,fld_u,it,cb%grad(:,:,1,1))
!                call cpu_time(toc)
!                tt6=tt6+toc-tic
!            endif
            
            !snapshot
            call fld_a%write(it,o_suffix='_rev')
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
            write(*,*) 'Total elapsed time for forward (min)',(tt1+tt2+tt3+tt7+tt8)/60./mpiworld%max_threads
            write(*,*) ' ---------------------------- '
            write(*,*) 'Elapsed time to add adjsource stresses   ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj stresses      ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to add adjsource velocities ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj velocities    ',tt10/mpiworld%max_threads
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
    subroutine inject_velocities(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        if(.not. f%is_adjoint) then

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler
            
            if(if_hicks) then
                select case (shot%src%comp)
                case ('vz')
                    f%vz(ifz:ilz,ifx:ilx,1) = f%vz(ifz:ilz,ifx:ilx,1) + wl*self%buoz(ifz:ilz,ifx:ilx) *shot%src%interp_coef(:,:,1)
                
                case ('vx')
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl !required to pass adjointtest. Why weaker when vx as src?
                    f%vx(ifz:ilz,ifx:ilx,1) = f%vx(ifz:ilz,ifx:ilx,1) + wl*self%buox(ifz:ilz,ifx:ilx) *shot%src%interp_coef(:,:,1)
                    
                end select
                
            else
                select case (shot%src%comp)
                case ('vz') !vertical force     on vz[iz-0.5,ix]
                    f%vz(iz,ix,1) = f%vz(iz,ix,1) + wl*self%buoz(iz,ix)
                    
                case ('vx') !horizontal x force on vx[iz,ix-0.5]
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl !required to pass adjointtest. Why weaker when vx as src?
                    f%vx(iz,ix,1) = f%vx(iz,ix,1) + wl*self%buox(iz,ix)
                    
                end select
                
            endif

            return

        endif


            do i=1,shot%nrcv
                ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                
                wl=f%wavelet(i,it)*wavelet_scaler
                
                if(if_hicks) then
                    select case (shot%rcv(i)%comp)
                    case ('vz') !vertical z adjsource
                        f%vz(ifz:ilz,ifx:ilx,1) = f%vz(ifz:ilz,ifx:ilx,1) + wl*self%buoz(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!

                    case ('vx') !horizontal x adjsource
                        if(m%is_freesurface.and.shot%rcv(i)%iz==1) wl=2*wl !required to pass adjointtest. Why weaker when vx as src?
                        f%vx(ifz:ilz,ifx:ilx,1) = f%vx(ifz:ilz,ifx:ilx,1) + wl*self%buox(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                        
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                    case ('vz') !vertical z adjsource
                        !vz[ix,1,iz-0.5]
                        f%vz(iz,ix,1) = f%vz(iz,ix,1) + wl*self%buoz(iz,ix) !no time_dir needed!

                    case ('vx') !horizontal x adjsource
                        !vx[ix-0.5,1,iz]
                        if(m%is_freesurface.and.shot%rcv(i)%iz==1) wl=2*wl !required to pass adjointtest. Why weaker when vx as src?
                        f%vx(iz,ix,1) = f%vx(iz,ix,1) + wl*self%buox(iz,ix) !no time_dir needed!
                        
                    end select
                    
                endif
                
            enddo
        
    end subroutine
    
    !forward: v^it -> v^it+1 by FD  of s^it+0.5
    !adjoint: v^it+1 -> v^it by FDᵀ of s^it+0.5
    subroutine update_velocities(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ifz=f%bloom(1,it)+2
        ilz=f%bloom(2,it)-2  !-1
        ifx=f%bloom(3,it)+2
        ilx=f%bloom(4,it)-2  !-1

        if(m%is_freesurface) ifz=max(ifz,1)

        call fd2d_velocities(f%vz,f%vx,f%szz,f%sxx,f%szx,            &
                             f%dszz_dz,f%dsxx_dx,f%dszx_dz,f%dszx_dx,&
                             self%buoz,self%buox,                    &
                             ifz,ilz,ifx,ilx,time_dir*self%dt)

        
        !apply free surface boundary condition if needed
        !free surface is located at [1,ix,1] level
        if(m%is_freesurface) then
            if (FS_method=='zero_stress') then
                !Δszz=0 : -(λ+2μ)∂_z vz = λ ∂ₓvx
                !         -(λ+2μ)[1,ix]*(vz[1.5,ix]-vz[0.5,ix])/dz = λ[1,ix]*(vx[1,ix-0.5]-vx[1,ix+0.5])/dx
                !         -(λ+2μ)(1,ix)*(vz(2  ,ix)-vz(1  ,ix))/dz = λ(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
                !          (λ+2μ)(1,ix)*(vz(1  ,ix)-vz(2  ,ix))/dz = λ(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
                !Δszx=0 : -∂_z vx = ∂ₓvz
                !         -(vx[2,ix+0.5]-vx[0,ix+0.5])/2dz = ( (vz[0.5,ix]+vz[1.5,ix])/2 - (vz[0.5,ix+1]+vz[1.5,ix+1])/2 )/dx
                !         -(vx(2,ix+1  )-vx(0,ix+1  ))/2dz = ( (vz(1  ,ix)+vz(2  ,ix))/2 - (vz(1  ,ix+1)+vz(2  ,ix+1))/2 )/dx
                !          (vx(0,ix+1  )-vx(2,ix+1  ))/ dz = ( (vz(1  ,ix)+vz(2  ,ix))   -  vz(1  ,ix+1)-vz(2  ,ix+1)    )/dx
                dz_dx = m%dz/m%dx

                do ix=ifx,ilx
                    f%vz(1,ix,1)= f%vz(2,ix,1) + self%lda(1,ix)*(f%vx(1,ix,1)-f%vx(1,ix+1,1))*dz_dx/self%ldap2mu(1,ix)
                    !f%vx(0,ix,1)= f%vx(2,ix,1) + (f%vz(1,ix,1)+f%vz(2,ix,1)-f%vz(1,ix+1,1)-f%vz(2,ix+1,1))*dz_dx/self%ldap2mu(1,ix) !toy2del says this condition is not needed
                enddo

            elseif (FS_method=='stress_image') then !Levandar & Roberttson
                !Roberttson's 3rd method
                f%vz(cb%ifz:1,:,1)=0.
                f%vx(cb%ifz:0,:,1)=0.

! do ix=cb%ifx+2,cb%ilx-2

!     ! dszz_dz_= c1z*(f%szz(2,ix  ,1)-f%szz(1,ix,1)) +c2z*(f%szz(3,ix  ,1)-f%szz(0,ix  ,1))
!     ! dszx_dx_= c1x*(f%szx(2,ix+1,1)-f%szx(2,ix,1)) +c2x*(f%szx(2,ix+2,1)-f%szx(2,ix-1,1))
    
!     ! f%dszz_dz(2,ix,1)= cpml%b_z_half(2)*f%dszz_dz(2,ix,1) + cpml%a_z_half(2)*dszz_dz_
!     ! f%dszx_dx(2,ix,1)= cpml%b_x(ix)    *f%dszx_dx(2,ix,1) + cpml%a_x(ix)    *dszx_dx_

!     ! dszz_dz_=dszz_dz_*cpml%kpa_z_half(2) + f%dszz_dz(2,ix,1)
!     ! dszx_dx_=dszx_dx_*cpml%kpa_x(2)      + f%dszx_dx(2,ix,1)

!     ! f%vz(2,ix,1)=f%vz(2,ix,1) + time_dir*self%dt*self%buoz(2,ix)*(dszz_dz_+dszx_dx_)



!     dsxx_dx_= c1x*(f%sxx(1,ix,1)-f%sxx(1,ix-1,1)) +c2x*(f%sxx(1,ix+1,1)-f%sxx(1,ix-2,1))
!     dszx_dz_= c1z*(f%szx(2,ix,1)-f%szx(1,ix  ,1)) +c2z*(f%szx(3,ix  ,1)-f%szx(0,ix  ,1))
!     !dszx_dz_=f%szx(2,ix,1)
                    
!     ! !cpml
!     ! f%dsxx_dx(1,ix,1)= cpml%b_x_half(ix)*f%dsxx_dx(1,ix,1) + cpml%a_x_half(ix)*dsxx_dx_
!     ! f%dszx_dz(1,ix,1)= cpml%b_z(1)      *f%dszx_dz(1,ix,1) + cpml%a_z(1)      *dszx_dz_

!     ! dsxx_dx_=dsxx_dx_*cpml%kpa_x_half(ix) + f%dsxx_dx(1,ix,1)  !kappa's should have been inversed in m_computebox.f90
!     ! dszx_dz_=dszx_dz_*cpml%kpa_z(1)       + f%dszx_dz(1,ix,1)
    
!     !velocity
!     f%vx(1,ix,1)=f%vx(1,ix,1) + time_dir*self%dt*self%buox(1,ix)*(dszx_dz_+dsxx_dx_)

! enddo




            else  !FS_method='effective_medium' (Mittet, Cao & Chen)            

                !required for high-ord FD
                f%vz(cb%ifz:1,:,1)=0.
                f%vx(cb%ifz:0,:,1)=0.

                do ix=ifx,ilx
                    dszz_dz_= (f%szz(2,ix,1)                )/m%dz !c1z*(f%szz(2,ix,1)-f%szz(1,ix  ,1)) +c2z*(f%szz(3,ix  ,1)-f%szz(0,ix,1))
                    dsxx_dx_= (f%sxx(1,ix,1)-f%sxx(1,ix-1,1))/m%dx !c1x*(f%sxx(1,ix,1)-f%sxx(1,ix-1,1)) +c2x*(f%sxx(1,ix+1,1)-f%sxx(1,ix-2,1))

                    ! dszx_dz_= c1z*(f%szx(2,ix,1)-f%szx(1,ix,1)) +c2z*(f%szx(3,ix,1)-f%szx(0,ix,1))
                    dszx_dx_= (f%szx(2,ix+1,1)-f%szx(2,ix,1))/m%dx !c1x*(f%szx(2,ix+1,1)-f%szx(2,ix,1)) +c2x*(f%szx(2,ix+2,1)-f%szx(2,ix-1,1))
                                    
                    !velocity
                    f%vz(2,ix,1)=f%vz(2,ix,1) + self%dt*   self%buoz(2,ix)*(dszz_dz_           +dszx_dx_)

                    !f%vx(i)=f%vx(i) + self%dt*2.*self%buox(i)*(dszx_dz_+dsxx_dx_)
                    f%vx(1,ix,1)=f%vx(1,ix,1) + self%dt*2.*self%buox(1,ix)*(f%szx(2,ix,1)/m%dz +dsxx_dx_)

                enddo

                ! do ix=ifx,ilx

                                    
                !     do iz=cb%ifz,2
                !         dszz_dz_= (f%szz(iz,ix,1)                 )/m%dz !c1z*(f%szz(2,ix,1)-f%szz(1,ix  ,1)) +c2z*(f%szz(3,ix  ,1)-f%szz(0,ix,1))
                !         dszx_dx_= (f%szx(iz,ix+1,1)-f%szx(iz,ix,1))/m%dx !c1x*(f%szx(2,ix+1,1)-f%szx(2,ix,1)) +c2x*(f%szx(2,ix+2,1)-f%szx(2,ix-1,1))
                !         !velocity
                !         f%vz(iz,ix,1)=f%vz(iz,ix,1) + self%dt*   self%buoz(iz,ix)*(dszz_dz_           +dszx_dx_)
                !     enddo

                !     do iz=cb%ifz,1
                !         ! dszx_dz_= c1z*(f%szx(2,ix,1)-f%szx(1,ix,1)) +c2z*(f%szx(3,ix,1)-f%szx(0,ix,1))
                !         dsxx_dx_= (f%sxx(iz,ix,1)-f%sxx(iz,ix-1,1))/m%dx !c1x*(f%sxx(1,ix,1)-f%sxx(1,ix-1,1)) +c2x*(f%sxx(1,ix+1,1)-f%sxx(1,ix-2,1))
                !         !f%vx(i)=f%vx(i) + self%dt*2.*self%buox(i)*(dszx_dz_+dsxx_dx_)
                !         f%vx(iz,ix,1)=f%vx(iz,ix,1) + self%dt*2.*self%buox(iz,ix)*(f%szx(iz+1,ix,1)/m%dz +dsxx_dx_)
                !     enddo

                ! enddo

            endif

        endif

    end subroutine

    !forward: add RHS to s^it+0.5
    !adjoint: add RHS to s^it+1.5
    subroutine inject_stresses(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        if(.not. f%is_adjoint) then

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler
            
            if(if_hicks) then
                select case (shot%src%comp)
                case ('p')
                    f%szz(ifz:ilz,ifx:ilx,1) = f%szz(ifz:ilz,ifx:ilx,1) +wl*sn_p*( self%ldap2mu(ifz:ilz,ifx:ilx)*shot%src%interp_coef_anti(:,:,1) &
                                                                                  +self%lda    (ifz:ilz,ifx:ilx)*shot%src%interp_coef_symm(:,:,1) )
                    f%sxx(ifz:ilz,ifx:ilx,1) = f%sxx(ifz:ilz,ifx:ilx,1) +wl*sn_p*( self%lda    (ifz:ilz,ifx:ilx)*shot%src%interp_coef_anti(:,:,1) &
                                                                                  +self%ldap2mu(ifz:ilz,ifx:ilx)*shot%src%interp_coef_symm(:,:,1) )
                case ('szz')
                    f%szz(ifz:ilz,ifx:ilx,1) = f%szz(ifz:ilz,ifx:ilx,1) + wl*self%ldap2mu(ifz:ilz,ifx:ilx)*shot%src%interp_coef_anti(:,:,1)
                    f%sxx(ifz:ilz,ifx:ilx,1) = f%sxx(ifz:ilz,ifx:ilx,1) + wl*self%lda    (ifz:ilz,ifx:ilx)*shot%src%interp_coef_anti(:,:,1)
                case ('sxx')
                    f%szz(ifz:ilz,ifx:ilx,1) = f%szz(ifz:ilz,ifx:ilx,1) + wl*self%lda    (ifz:ilz,ifx:ilx)*shot%src%interp_coef_symm(:,:,1)
                    f%sxx(ifz:ilz,ifx:ilx,1) = f%sxx(ifz:ilz,ifx:ilx,1) + wl*self%ldap2mu(ifz:ilz,ifx:ilx)*shot%src%interp_coef_symm(:,:,1)

                case ('szx')
                    f%szx(ifz:ilz,ifx:ilx,1) = f%szx(ifz:ilz,ifx:ilx,1) + wl*self%mu(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)
                
                end select
                
            else
                select case (shot%src%comp)
                case ('p')
                    !explosion on s[iz,ix,1]
                    f%szz(iz,ix,1) = f%szz(iz,ix,1) + wl*sn_p*self%two_ldapmu(iz,ix)
                    f%sxx(iz,ix,1) = f%sxx(iz,ix,1) + wl*sn_p*self%two_ldapmu(iz,ix)

                case ('szz')
                    f%szz(iz,ix,1) = f%szz(iz,ix,1) + wl*self%ldap2mu(iz,ix)
                    f%sxx(iz,ix,1) = f%sxx(iz,ix,1) + wl*    self%lda(iz,ix)
                case ('sxx')
                    f%szz(iz,ix,1) = f%szz(iz,ix,1) + wl*    self%lda(iz,ix)
                    f%sxx(iz,ix,1) = f%sxx(iz,ix,1) + wl*self%ldap2mu(iz,ix)

                case ('szx')
                    f%szx(iz,ix,1) = f%szx(iz,ix,1) + wl*self%mu(iz,ix)
                
                end select
            
            endif

            return

        endif

            do i=1,shot%nrcv

                ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                
                !adjsource for pressure
                wl=f%wavelet(i,it)*wavelet_scaler
                
                if(if_hicks) then 

                    select case (shot%rcv(i)%comp)
                    case ('p')
                        f%szz(ifz:ilz,ifx:ilx,1) = f%szz(ifz:ilz,ifx:ilx,1) +wl*sn_p*( self%ldap2mu(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef_anti(:,:,1) &
                                                                                      +self%lda    (ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef_symm(:,:,1) )
                        f%sxx(ifz:ilz,ifx:ilx,1) = f%sxx(ifz:ilz,ifx:ilx,1) +wl*sn_p*( self%lda    (ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef_anti(:,:,1) &
                                                                                      +self%ldap2mu(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef_symm(:,:,1) )

                    case ('szz')
                        f%szz(ifz:ilz,ifx:ilx,1) = f%szz(ifz:ilz,ifx:ilx,1) + wl*self%ldap2mu(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef_anti(:,:,1)
                        f%sxx(ifz:ilz,ifx:ilx,1) = f%sxx(ifz:ilz,ifx:ilx,1) + wl*    self%lda(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef_anti(:,:,1)
                    case ('sxx')
                        f%szz(ifz:ilz,ifx:ilx,1) = f%szz(ifz:ilz,ifx:ilx,1) + wl*    self%lda(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef_symm(:,:,1)
                        f%sxx(ifz:ilz,ifx:ilx,1) = f%sxx(ifz:ilz,ifx:ilx,1) + wl*self%ldap2mu(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef_symm(:,:,1)

                    case ('szx')
                        f%szx(ifz:ilz,ifx:ilx,1) = f%szx(ifz:ilz,ifx:ilx,1) + wl*self%mu(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1)

                    end select

                else           
                    select case (shot%rcv(i)%comp)
                    case ('p')
                        !s[iz,ix,1]
                        f%szz(iz,ix,1) = f%szz(iz,ix,1) +wl*sn_p*self%two_ldapmu(iz,ix) !no time_dir needed!
                        f%sxx(iz,ix,1) = f%sxx(iz,ix,1) +wl*sn_p*self%two_ldapmu(iz,ix) !no time_dir needed!

                    case ('szz')
                        f%szz(iz,ix,1) = f%szz(iz,ix,1) + wl*self%ldap2mu(iz,ix)
                        f%sxx(iz,ix,1) = f%sxx(iz,ix,1) + wl*    self%lda(iz,ix)
                    case ('sxx')
                        f%szz(iz,ix,1) = f%szz(iz,ix,1) + wl*    self%lda(iz,ix)
                        f%sxx(iz,ix,1) = f%sxx(iz,ix,1) + wl*self%ldap2mu(iz,ix)

                    case ('szx')
                        f%szx(iz,ix,1) = f%szx(iz,ix,1) + wl*self%mu(iz,ix)

                    end select

                endif

            enddo
        
    end subroutine

    !forward: s^it+0.5 -> s^it+1.5 by FD of v^it+1
    !adjoint: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
    subroutine update_stresses(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        real factor

        ifz=f%bloom(1,it)+2  !1
        ilz=f%bloom(2,it)-2
        ifx=f%bloom(3,it)+2  !1
        ilx=f%bloom(4,it)-2
        
        if(m%is_freesurface) ifz=max(ifz,2)

        call fd2d_stresses(f%vz,f%vx,f%szz,f%sxx,f%szx,    &
                           f%dvz_dz,f%dvx_dx,f%dvz_dx,f%dvx_dz,&
                           self%ldap2mu,self%lda,self%mu,  &
                           ifz,ilz,ifx,ilx,time_dir*self%dt)
        

        !apply free surface boundary condition if needed
        !free surface is located at [1,ix] level
        if(m%is_freesurface) then
            if (FS_method=='zero_stress') then
                !so explicit boundary condition: szz(1,ix)=0
                !and antisymmetric mirroring: szx[0.5,ix-0.5]=-szx[1.5,ix-0.5] -> szx(1,ix)=-szx(2,ix)
                f%szz(1,:,1)=0.
                f%szx(1,:,1)=-f%szx(2,:,1)

            elseif (FS_method=='stress_image') then !Levandar & Roberttson

                !image szz
                ! f%szz(-2,:,1)=-f%szz(4,:,1)
                ! f%szz(-1,:,1)=-f%szz(3,:,1)
                ! f%szz( 0,:,1)=-f%szz(2,:,1)
                f%szz( 1,:,1)=0.
                f%szz(0:cb%ifz:-1, :,1)=-f%szz(2:2+0-cb%ifz, :,1)
                ! if(f%is_adjoint) f%szz(cb%ifz:0, :,1)=0.

                !not image on sxx
                ! f%sxx(0:cb%ifz:-1,:,1)=0. !no needed
                
                ! do ix=cb%ifx+1,cb%ilx-2
                ! do iz=cb%ifz+1,0
                !     dvz_dz_= c1z*(f%vz(iz+1,ix,1)-f%vz(iz,ix,1))  +c2z*(f%vz(iz+2,ix,1)-f%vz(iz-1,ix,1))
                !     dvx_dx_= c1x*(f%vx(iz,ix+1,1)-f%vx(iz,ix,1))  +c2x*(f%vx(iz,ix+2,1)-f%vx(iz,ix-1,1))                
                !     !normal stresses
                !     f%sxx(iz,ix,1)  = f%sxx(iz,ix,1) + time_dir*self%dt * (self%lda(iz,ix)*dvz_dz_ + self%ldap2mu(iz,ix)*dvx_dx_)
                ! enddo
                ! enddo
                do ix=cb%ifx+1,cb%ilx-2
                    dvx_dx_= c1x*(f%vx(1,ix+1,1)-f%vx(1,ix,1))  +c2x*(f%vx(1,ix+2,1)-f%vx(1,ix-1,1))
                    
                    factor=-self%lda(1,ix)**2/self%ldap2mu(1,ix) + self%ldap2mu(1,ix)
                    f%sxx(1,ix,1)  = f%sxx(1,ix,1) + time_dir*self%dt * factor*dvx_dx_
                enddo
                
                ! if(f%is_adjoint) then
                !     f%sxx( 1,:,1)=0.
                !     nnz=0-cb%ifz
                !     f%sxx(cb%ifz:0, :,1)=-f%sxx(2+nnz:2:-1, :,1)
                ! endif

! do ix=cb%ifx+1,cb%ilx-2

!     dvx_dx_= c1x*(f%vx(1,ix+1,1)-f%vx(1,ix,1))  +c2x*(f%vx(1,ix+2,1)-f%vx(1,ix-1,1))
    
!     ! !cpml
!     ! f%dvx_dx(1,ix,1)=cpml%b_x(ix)*f%dvx_dx(1,ix,1)+cpml%a_x(ix)*dvx_dx_

!     ! dvx_dx_=dvx_dx_*cpml%kpa_x(ix) + f%dvx_dx(iz,ix,1)
    
!     factor=-self%lda(1,ix)**2/self%ldap2mu(1,ix) + self%ldap2mu(1,ix)

!     !normal stresses
!     f%sxx(1,ix,1) = f%sxx(1,ix,1) + time_dir*self%dt * factor*dvx_dx_

! enddo


                !image szx
                ! f%szx(-1,:,1)=-f%szx(4,:,1)
                ! f%szx( 0,:,1)=-f%szx(3,:,1)
                ! f%szx( 1,:,1)=-f%szx(2,:,1)            
                !nnz=1-cb%ifz
                f%szx(1:cb%ifz:-1, :,1)=-f%szx(2:2+1-cb%ifz, :,1)
                ! if(f%is_adjoint) f%szx(cb%ifz:1, :,1)=0.

            else  !FS_method='effective_medium' (Mittet, Cao & Chen)            

                !required for high-ord FD
                ! f%szx(cb%ifz:1,:,1)=0.
                ! f%szz(cb%ifz:0,:,1)=0.
                f%szz( 1,:,1)=0.
                nnz=0-cb%ifz
                f%szz(cb%ifz:0, :,1)=-f%szz(2+nnz:2:-1, :,1)
                nnz=1-cb%ifz
                f%szx(cb%ifz:1, :,1)=-f%szx(2+nnz:2:-1, :,1)
                f%sxx(cb%ifz:0,:,1)=0.

                do ix=ifx,ilx
                
                    dvx_dx_= (f%vx(1,ix+1,1)-f%vx(1,ix,1))/m%dx !c1x*(f%vx(1,ix+1,1)-f%vx(1,ix,1))  +c2x*(f%vx(1,ix+2,1)-f%vx(1,ix-1,1))
                    !dvx_dx_= c1x*(f%vx(1,ix+1,1)-f%vx(1,ix,1))  +c2x*(f%vx(1,ix+2,1)-f%vx(1,ix-1,1))
                    
                    !normal stresses
                    f%szz(1,ix,1) = 0.

                    factor=-self%lda(1,ix)**2/self%ldap2mu(1,ix) + self%ldap2mu(1,ix)
                    factor=factor/2.
                    f%sxx(1,ix,1) = f%sxx(1,ix,1) + time_dir*self%dt * factor*dvx_dx_


                    dvz_dx_= (f%vz(2,ix,1)-f%vz(2,ix-1,1))/m%dx
                    !dvz_dx_= c1x*(f%vz(2,ix,1)-f%vz(2,ix-1,1))  +c2x*(f%vz(2,ix+1,1)-f%vz(2,ix-2,1))
                    dvx_dz_= (f%vx(2,ix,1)-f%vx(1,ix  ,1))/m%dz
                    !dvx_dz_= c1z*(f%vx(2,ix,1)-f%vx(1,ix  ,1))  +c2z*(f%vx(3,ix  ,1)-f%vx(0,ix  ,1))
                    
                    !shear stress
                    f%szx(2,ix,1) = f%szx(2,ix,1) + time_dir*self%dt * self%mu(2,ix)*(dvz_dx_+dvx_dz_)
                        
                enddo

                
                ! do ix=ifx,ilx
                
                !     do iz=cb%ifz,1
                !         !normal stresses
                !         f%szz(iz,ix,1) = 0.
                        
                !         dvx_dx_= (f%vx(iz,ix+1,1)-f%vx(iz,ix,1))/m%dx

                !         factor=-self%lda(iz,ix)**2/self%ldap2mu(iz,ix) + self%ldap2mu(iz,ix)
                !         factor=factor/2.
                !         f%sxx(iz,ix,1) = f%sxx(iz,ix,1) + self%dt * factor*dvx_dx_

                !     enddo

                !     do iz=cb%ifz+1,2
                !         dvz_dx_= (f%vz(iz,ix,1)-f%vz(iz  ,ix-1,1))/m%dx
                !         dvx_dz_= (f%vx(iz,ix,1)-f%vx(iz-1,ix  ,1))/m%dz
                        
                !         !shear stress
                !         f%szx(iz,ix,1) = f%szx(iz,ix,1) + self%dt * self%mu(iz,ix)*(dvz_dx_+dvx_dz_)

                !     enddo
                        
                ! enddo

            endif

        endif

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
                    case ('vz')
                        f%seismo(i,it)=sum(f%vz(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1))
                    case ('vx')
                        f%seismo(i,it)=sum(f%vx(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1))

                    case ('p')
                        f%seismo(i,it)=sum(sn_p*(f%szz(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef_anti (:,:,1)&
                                                +f%sxx(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef_trunc(:,:,1)))
                    case ('szz')
                        f%seismo(i,it)=sum(f%szz(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef_anti (:,:,1))
                    case ('sxx')
                        f%seismo(i,it)=sum(f%sxx(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef_trunc(:,:,1))

                    case ('szx')
                        f%seismo(i,it)=sum(f%szx(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )

                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                    case ('vz') !vz[iz-0.5,ix]
                        f%seismo(i,it)=f%vz(iz,ix,1)
                    case ('vx') !vx[iz,ix-0.5]
                        f%seismo(i,it)=f%vx(iz,ix,1)

                    case ('p') !p[iz,ix]
                        f%seismo(i,it)=(f%szz(iz,ix,1)+f%sxx(iz,ix,1))*sn_p
                    case ('szz')
                        f%seismo(i,it)=f%szz(iz,ix,1)
                    case ('sxx')
                        f%seismo(i,it)=f%sxx(iz,ix,1)

                    case ('szx')
                        f%seismo(i,it)=f%szx(iz,ix,1)

                    end select
                    
                endif

            enddo

            return

        endif

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            if(if_hicks) then
                select case (shot%src%comp)
                case ('vz')
                    f%seismo(1,it)=sum(f%vz(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                case ('vx')
                    f%seismo(1,it)=sum(f%vx(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))

                case ('p')
                    f%seismo(1,it)=sum(sn_p*(f%szz(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef_anti (:,:,1)&
                                            +f%sxx(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef_trunc(:,:,1)))
                case ('szz')
                    f%seismo(1,it)=sum(f%szz(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef_anti (:,:,1))
                case ('sxx')
                    f%seismo(1,it)=sum(f%sxx(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef_trunc(:,:,1))

                case ('szx')
                    f%seismo(1,it)=sum(f%szx(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1) )
                    
                end select
                
            else
                select case (shot%src%comp)                    
                case ('vz') !vz[iz-0.5,ix,1]
                    f%seismo(1,it)=f%vz(iz,ix,1)
                    
                case ('vx') !vx[iz,ix-0.5,1]
                    f%seismo(1,it)=f%vx(iz,ix,1)

                case ('p') !p[iz,ix,1]
                    f%seismo(1,it)=(f%szz(iz,ix,1)+f%sxx(iz,ix,1))*sn_p
                case ('szz')
                    f%seismo(1,it)=f%szz(iz,ix,1)
                case ('sxx')
                    f%seismo(1,it)=f%sxx(iz,ix,1)

                case ('szx')
                    f%seismo(1,it)=f%szx(iz,ix,1)
                                        
                end select
                
            endif
        
    end subroutine
    
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%buoz, self%buox)
        call dealloc(self%ldap2mu, self%lda, self%mu)
        call dealloc(self%two_ldapmu, self%inv_ldapmu_4mu)
    end subroutine

    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !Kₘ<a|Au> = Kₘ<a|M∂ₜu-Du> = ∫ aᵀ KₘM ∂ₜu dt
    !Since it's cumbersome to get ∂ₜu by time marching,
    !replace ∂ₜu by M⁻¹Du and neglect f
    !ie. M∂ₜu=Du+f -> ∂ₜu=M⁻¹Du+M⁻¹f ≐ M⁻¹Du
    !This simplification introduces singularities in the gradient only at source positions, which are probably removed by gradient masking.
    !
    !Therefore, ∫ aᵀ KₘM ∂ₜu dt ≐ ∫ aᵀ KₘM M⁻¹Du dt
    !
    !M₂ = diag₂(ρ)
    !K_ρM₂ M₂⁻¹ = diag₂(1) diag₂(b) = diag₂(b)
    !
    !     [λ+2μ  λ    ]⁻¹   [___1____[(λ+2μ)     -λ]    ]
    !M₃ = | λ   λ+2μ  |   = |4(λμ+μ²)[    -λ (λ+2μ)]    |
    !     [          μ]     [                        1/μ]
    !
    !        __-1___[1 1 0]
    !K_λM₃ = 4(λ+μ)²|1 1 0|
    !               [0 0 0]
    !
    !             __-1___[1 1 0][λ+2μ  λ    ]   __-1__[1 1 0]
    !K_λM₃ M₃⁻¹ = 4(λ+μ)²|1 1 0|| λ   λ+2μ  | = 2(λ+μ)|1 1 0|
    !                    [0 0 0][          μ]         [0 0 0]
    !
    !        [___-1____[ λ²+2λμ+2μ² -λ²-2λμ    ]  0  ]
    !K_µM₃ = |4(λμ+μ²)²[-λ²-2λμ      λ²+2λμ+2μ²]  0  |
    !        [              0           0        -μ⁻²]
    !
    !             [___-1____[ λ²+2λμ+2μ² -λ²-2λμ    ][λ+2μ  λ    ]     ]   [__-1___[λ+2μ  -λ ]     ]
    !K_μM₃ M₃⁻¹ = |4(λμ+μ²)²[-λ²-2λμ      λ²+2λμ+2μ²][ λ   λ+2μ  ]     | = |2(λ+μ)μ[ -λ  λ+2μ]     |
    !             [                                                -μ⁻¹]   [                   -μ⁻¹]
    !
    !       [∂zᵇszz + ∂ₓᶠszx]
    !       |∂ₓᵇsxx + ∂zᶠszx|
    !and Du=|∂zᶠvz          |
    !       |∂ₓᶠvx          |
    !       [∂ₓᵇvz  + ∂zᵇvx ]
    !
    !Therefore,
    !grho = [vzᵃ vxᵃ] diag₂(b) [∂zᵇszz + ∂ₓᶠszx]
    !                          |∂ₓᵇsxx + ∂zᶠszx|
    !     = b[ vzᵃ(∂zᵇszz +∂ₓᶠszx) + vxᵃ(∂ₓᵇsxx+ ∂zᶠszx) ]
    !
    !glda = [szzᵃ sxxᵃ szxᵃ] __-1__[1 1 0] [∂zᶠvz          ]
    !                        2(λ+μ)|1 1 0| |∂ₓᶠvx          |
    !                              [0 0 0] [∂ₓᵇvz + ∂zᵇvx  ]
    !     = __-1__[ szzᵃ(∂zᶠvz +∂ₓᶠvx) + sxxᵃ(∂zᶠvz +∂ₓᶠvx) ]
    !       2(λ+μ)  
    !
    !gmu = [szzᵃ sxxᵃ szxᵃ] __-1___[λ+2μ  -λ        ] [∂zᶠvz        ]
    !                       2(λ+μ)μ[ -λ  λ+2μ       | |∂ₓᶠvx        |
    !                              [          2(λ+μ)] [∂ₓᵇvz + ∂zᵇvx]
    !    = __-1___[ szzᵃ((λ+2μ)∂zᶠvz -λ*∂ₓᶠvx) + sxxᵃ(-λ*∂zᶠvz +(λ+2μ)∂ₓᶠvx) + szxᵃ*2(λ+μ)(∂ₓᵇvz +∂zᵇvx) ]
    !      2(λ+μ)μ 


    !On the free surface
    !PDE:      A u = M ∂ₜ u - D u = f
    !where
    !u=[vx sxx]ᵀ,
    !f=[fz fx fzz fxx 0]ᵀδ(x-xs) with xs source position, d is recorded data
    !  [ρ        ]    [0  ∂ₓ]
    !M=|  _λ+2μ__|, D=[∂ₓ 0 ]
    !  [  4μ(λ+μ)]
    !since 0=sz=(λ+2μ)∂zvz+λ∂ₓvx => sx=λ∂zvz+(λ+2μ)∂ₓvx=4μ(λ+μ)/(λ+2μ)*∂ₓvx
    !
    !K_λM M⁻¹ =  [0        ][ρ          ] = [0            ]
    !            |  __-1__ ||  _4μ(λ+μ)_|   |  ____-μ_____|
    !            [  4(λ+μ)²][    λ+2μ   ]   [  (λ+μ)(λ+2μ)]
    !
    !compare w/ the 'full-elastic' case:
    !        __-1___[1 1 0]
    !K_λM₃ = 4(λ+μ)²|1 1 0|
    !               [0 0 0]
    !
    !             __-1___[1 1 0][λ+2μ  λ    ]   __-1__[1 1 0]
    !K_λM₃ M₃⁻¹ = 4(λ+μ)²|1 1 0|| λ   λ+2μ  | = 2(λ+μ)|1 1 0|
    !                    [0 0 0][          μ]         [0 0 0]
    !
    !
    !K_µM M⁻¹ =  [0             ][ρ          ] = [0               ]
    !            |  _-(λ+μ)²-μ²_||  _4μ(λ+μ)_|   |  _-4(λ+μ)²-4μ²_|
    !            [    μ²(λ+μ)²  ][    λ+2μ   ]   [   μ(λ+μ)(λ+2μ) ]
    !
    !compare w/ the 'full-elastic' case
    !        [___-1____[ λ²+2λμ+2μ² -λ²-2λμ    ]  0  ]
    !K_µM₃ = |4(λμ+μ²)²[-λ²-2λμ      λ²+2λμ+2μ²]  0  |
    !        [              0           0        -μ⁻²]
    !
    !             [___-1____[ λ²+2λμ+2μ² -λ²-2λμ    ][λ+2μ  λ    ]     ]   [__-1___[λ+2μ  -λ ]     ]
    !K_μM₃ M₃⁻¹ = |4(λμ+μ²)²[-λ²-2λμ      λ²+2λμ+2μ²][ λ   λ+2μ  ]     | = |2(λ+μ)μ[ -λ  λ+2μ]     |
    !             [                                                -μ⁻¹]   [                   -μ⁻¹]
    !
    !Therefore, 
    !glda = ____-μ_____[ sxxᵃ * ∂ₓᶠvx ]
    !       (λ+μ)(λ+2μ)  
    !
    !gmu = _-4(λ+μ)²-4μ²_[ sxxᵃ * ∂ₓᶠvx ]
    !       μ(λ+μ)(λ+2μ)  


    
   ! subroutine gradient_density(rf,sf,it,grad)
   !     type(t_field), intent(in) :: rf, sf
   !     real,dimension(cb%mz,cb%mx,1) :: grad
       
   !     !nonzero only when sf touches rf
   !     ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
   !     ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
   !     ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
   !     ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
       
   !     call grad2d_density(rf%vz,rf%vx,sf%szz,sf%sxx,sf%szx,&
   !                         grad,                            &
   !                         ifz,ilz,ifx,ilx)
       
   ! end subroutine

    subroutine cross_correlate_glda_gmu(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2) !
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
                    
        !inexact greadient
        call grad2d_glda_gmu(rf%szz,rf%sxx,rf%szx,sf%vz,sf%vx,  &
                           ppg%ldap2mu,ppg%lda,ppg%two_ldapmu,&
                           corr%glda,corr%gmu,            &
                           ifz,ilz,ifx,ilx)

        ! if(m%is_freesurface) call grad2d_freesurf_glda_gmu(rf%sxx,sf%vx,&
        !                     corr%glda,corr%gmu,            &
        !                     ifx,ilx)
        
    end subroutine
    
    subroutine cross_correlate_postprocess(corr)
        type(t_correlate) :: corr

        real,dimension(:,:),allocatable :: mhalf_inv_ldapmu
        
        !grho
        corr%grho(:,:,1) = corr%grho(:,:,1) / cb%rho(1:cb%mz,1:cb%mx,1)

        mhalf_inv_ldapmu = -0.5/(ppg%ldap2mu(1:cb%mz,1:cb%mx)-ppg%mu(1:cb%mz,1:cb%mx))

        !glda
        corr%glda(:,:,1) = corr%glda(:,:,1) *mhalf_inv_ldapmu
        !gmu
        corr%gmu(:,:,1) = corr%gmu(:,:,1) *mhalf_inv_ldapmu/ppg%mu(1:cb%mz,1:cb%mx)
        

        ! if(m%is_freesurface) then
        !     corr%glda(1,:,1) = corr%glda(1,:,1)/mhalf_inv_ldapmu(1,:) *& !canceling
        !         (-ppg%mu(1,1:cb%mx))                                /ppg%ldapmu(1,1:cb%mx)/ppg%ldap2mu(1,1:cb%mx)
        !     corr%gmu(1,:,1)  = corr%gmu(1,:,1) /mhalf_inv_ldapmu(1,:)*ppg%mu(1,1:cb%mx) *& !canceling
        !         (-4*ppg%ldapmu(1,1:cb%mx)**2-4*ppg%mu(1,1:cb%mx)**2)/ppg%mu(1,1:cb%mx)/ppg%ldapmu(1,1:cb%mx)/ppg%ldap2mu(1,1:cb%mx)

        ! else
            !preparing for projection back
            corr%grho(1,:,:) = corr%grho(2,:,:)
            corr%glda(1,:,:) = corr%glda(2,:,:)
            corr%gmu (1,:,:) = corr%gmu (2,:,:)

        ! endif

        where( ieee_is_nan(corr%gmu(:,:,1)) .or. .not.ieee_is_finite(corr%gmu(:,:,1)) )
            corr%gmu(:,:,1)=0.
        endwhere


        deallocate(mhalf_inv_ldapmu)

    end subroutine

    ! subroutine imaging(rf,sf,it,imag)
    !     type(t_field), intent(in) :: rf, sf
    !     real,dimension(cb%mz,cb%mx,1) :: imag
        
    !     !nonzero only when sf touches rf
    !     ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
    !     ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
    !     ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
    !     ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        
    !     call imag2d_xcorr(rf%szz,rf%sxx, &
    !                       sf%szz,sf%sxx, &
    !                       imag,          &
    !                       ifz,ilz,ifx,ilx)
        
    ! end subroutine

    ! subroutine imaging_postprocess

    !     !scale by rdt tobe an image in the discretized world
    !     cb%imag = cb%imag*rdt

    !     !for cb%project_back
    !     cb%imag(1,:,:,:) = cb%imag(2,:,:,:)

    ! end subroutine

    ! subroutine energy(sf,it,engy)
    !     type(t_field),intent(in) :: sf
    !     real,dimension(cb%mz,cb%mx,1) :: engy

    !     !nonzero only when sf touches rf
    !     ifz=max(sf%bloom(1,it),2)
    !     ilz=min(sf%bloom(2,it),cb%mz)
    !     ifx=max(sf%bloom(3,it),1)
    !     ilx=min(sf%bloom(4,it),cb%mx)
        
    !     call engy2d_xcorr(sf%szz,sf%sxx, &
    !                       engy,          &
    !                       ifz,ilz,ifx,ilx)
        
    ! end subroutine

    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd2d_velocities(vz,vx,szz,sxx,szx,              &
                               dszz_dz,dsxx_dx,dszx_dz,dszx_dx,&
                               buoz,buox,                      &
                               ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vz,vx,szz,sxx,szx
        real,dimension(*) :: dszz_dz,dsxx_dx,dszx_dz,dszx_dx
        real,dimension(*) :: buoz,buox
        
        nz=cb%nz
        nx=cb%nx
        
        dszz_dz_=0.; dsxx_dx_=0.; dszx_dz_=0.; dszx_dx_=0.

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dszz_dz_,dsxx_dx_,dszx_dz_,dszx_dx_)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx

            !dir$ simd
            do iz=ifz,ilz

                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1
                
                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2=i  -2*nz  !iz,ix-2
                iz_ixm1=i    -nz  !iz,ix-1
                iz_ixp1=i    +nz  !iz,ix+1
                iz_ixp2=i  +2*nz  !iz,ix+2
                

                dszz_dz_= c1z*(szz(iz_ix)-szz(izm1_ix)) +c2z*(szz(izp1_ix)-szz(izm2_ix))
                dsxx_dx_= c1x*(sxx(iz_ix)-sxx(iz_ixm1)) +c2x*(sxx(iz_ixp1)-sxx(iz_ixm2))

                dszx_dz_= c1z*(szx(izp1_ix)-szx(iz_ix)) +c2z*(szx(izp2_ix)-szx(izm1_ix))
                dszx_dx_= c1x*(szx(iz_ixp1)-szx(iz_ix)) +c2x*(szx(iz_ixp2)-szx(iz_ixm1))
                                
                !cpml
                dszz_dz(i)= cpml%b_z_half(iz)*dszz_dz(i) + cpml%a_z_half(iz)*dszz_dz_
                dsxx_dx(i)= cpml%b_x_half(ix)*dsxx_dx(i) + cpml%a_x_half(ix)*dsxx_dx_
                dszx_dz(i)= cpml%b_z(iz)     *dszx_dz(i) + cpml%a_z(iz)     *dszx_dz_
                dszx_dx(i)= cpml%b_x(ix)     *dszx_dx(i) + cpml%a_x(ix)     *dszx_dx_

                dszz_dz_=dszz_dz_*cpml%kpa_z_half(iz) + dszz_dz(i)
                dsxx_dx_=dsxx_dx_*cpml%kpa_x_half(ix) + dsxx_dx(i)  !kappa's should have been inversed in m_computebox.f90
                dszx_dz_=dszx_dz_*cpml%kpa_z(iz)      + dszx_dz(i)
                dszx_dx_=dszx_dx_*cpml%kpa_x(ix)      + dszx_dx(i)
                
                !velocity
                vz(i)=vz(i) + dt*buoz(i)*(dszz_dz_+dszx_dx_)
                vx(i)=vx(i) + dt*buox(i)*(dszx_dz_+dsxx_dx_)

            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        

! do ix=ifx,ilx

!     iz=1

!         i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1
        
!         izm2_ix=i-2  !iz-2,ix
!         izm1_ix=i-1  !iz-1,ix
!         iz_ix  =i    !iz,ix
!         izp1_ix=i+1  !iz+1,ix
!         izp2_ix=i+2  !iz+2,ix
        
!         iz_ixm2=i  -2*nz  !iz,ix-2
!         iz_ixm1=i    -nz  !iz,ix-1
!         iz_ixp1=i    +nz  !iz,ix+1
!         iz_ixp2=i  +2*nz  !iz,ix+2
        
!         dszz_dz_= (szz(izp1_ix)-szz(iz_ix    ))/m%dz
!         dsxx_dx_= (sxx(izp1_ix)-sxx(izp1_ixm1))/m%dx

!         dszx_dz_= (szx(izp1_ix)           )/m%dz
!         dszx_dx_= (szx(iz_ixp1)-szx(iz_ix))/m%dx
        
!         !velocity
!         vz(iz_ixp1)=vz(iz_ixp1) + dt*   buoz(iz_ixp1)*(dszz_dz_+dszx_dx_)
!         vx(iz_ix  )=vx(iz_ix  ) + dt*2.*buox(iz_ix  )*(dszx_dz_+dsxx_dx_)

! enddo

    end subroutine
    
    subroutine fd2d_stresses(vz,vx,szz,sxx,szx,          &
                             dvz_dz,dvx_dx,dvz_dx,dvx_dz,&
                             ldap2mu,lda,mu,             &
                             ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vz,vx,szz,sxx,szx
        real,dimension(*) :: dvz_dz,dvx_dx,dvz_dx,dvx_dz
        real,dimension(*) :: ldap2mu,lda,mu
        
        nz=cb%nz
        nx=cb%nx
        
        dvz_dz_=0.
        dvx_dx_=0.
        dvz_dx_=0.
        dvx_dz_=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dvz_dz_,dvx_dx_,dvz_dx_,dvx_dz_)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1
                
                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2=i  -2*nz !iz,ix-2
                iz_ixm1=i    -nz !iz,ix-1
                iz_ixp1=i    +nz !iz,ix+1
                iz_ixp2=i  +2*nz !iz,ix+2
                

                dvz_dz_= c1z*(vz(izp1_ix)-vz(iz_ix))  +c2z*(vz(izp2_ix)-vz(izm1_ix))
                dvx_dx_= c1x*(vx(iz_ixp1)-vx(iz_ix))  +c2x*(vx(iz_ixp2)-vx(iz_ixm1))
                
                !cpml
                dvz_dz(i)=cpml%b_z(iz)*dvz_dz(i)+cpml%a_z(iz)*dvz_dz_
                dvx_dx(i)=cpml%b_x(ix)*dvx_dx(i)+cpml%a_x(ix)*dvx_dx_

                dvz_dz_=dvz_dz_*cpml%kpa_z(iz) + dvz_dz(iz_ix)
                dvx_dx_=dvx_dx_*cpml%kpa_x(ix) + dvx_dx(iz_ix)
                
                !normal stresses
                szz(i) = szz(i) + dt * (ldap2mu(i)*dvz_dz_ + lda(i)    *dvx_dx_)
                sxx(i) = sxx(i) + dt * (lda(i)    *dvz_dz_ + ldap2mu(i)*dvx_dx_)


                dvz_dx_= c1x*(vz(iz_ix)-vz(iz_ixm1))  +c2x*(vz(iz_ixp1)-vz(iz_ixm2))
                dvx_dz_= c1z*(vx(iz_ix)-vx(izm1_ix))  +c2z*(vx(izp1_ix)-vx(izm2_ix))

                !cpml
                dvz_dx(i)=cpml%b_x_half(ix)*dvz_dx(i)+cpml%a_x_half(ix)*dvz_dx_
                dvx_dz(i)=cpml%b_z_half(iz)*dvx_dz(i)+cpml%a_z_half(iz)*dvx_dz_

                dvz_dx_=dvz_dx_*cpml%kpa_x_half(ix) + dvz_dx(i)
                dvx_dz_=dvx_dz_*cpml%kpa_z_half(iz) + dvx_dz(i)

                !shear stress
                szx(i) = szx(i) + dt * mu(i)*(dvz_dx_+dvx_dz_)
                
            enddo
            
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine

    subroutine grad2d_glda_gmu(rf_szz,rf_sxx,rf_szx,&
                             sf_vz,sf_vx,         &
                             ldap2mu,lda,two_ldapmu,&
                             grad_lda,grad_mu,    &
                             ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_szz,rf_sxx,rf_szx,sf_vz,sf_vx
        real,dimension(*) :: ldap2mu,lda,two_ldapmu
        real,dimension(*) :: grad_lda,grad_mu
        
        nz=cb%nz
        
        sf_dvx_dx=0.
        sf_dvz_dz=0.
        sf_4_dvzdx_p_dvxdz=0.
        rf_4_szx=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,izp1_ixm2,iz_ixm1,izp1_ixm1,&
        !$omp         izm2_ixp1,izm1_ixp1,iz_ixp1,izp1_ixp1,izp2_ixp1,&
        !$omp         iz_ixp2,izp1_ixp2,&
        !$omp         sf_dvx_dx,sf_dvz_dz,&
        !$omp         sf_4_dvzdx_p_dvxdz,rf_4_szx)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers

                izm2_ix= i-2  !iz-2,ix
                izm1_ix= i-1  !iz-1,ix
                iz_ix  = i    !iz  ,ix
                izp1_ix= i+1  !iz+1,ix
                izp2_ix= i+2  !iz+2,ix

                iz_ixm2  = i    -2*nz  !iz  ,ix-2
                izp1_ixm2= i+1  -2*nz  !iz+1,ix-2
                iz_ixm1  = i      -nz  !iz  ,ix-1
                izp1_ixm1= i+1    -nz  !iz+1,ix-1

                izm2_ixp1= i-2    +nz  !iz-2,ix+1
                izm1_ixp1= i-1    +nz  !iz-1,ix+1
                iz_ixp1  = i      +nz  !iz  ,ix+1
                izp1_ixp1= i+1    +nz  !iz+1,ix+1
                izp2_ixp1= i+2    +nz  !iz+2,ix+1

                iz_ixp2  = i    +2*nz  !iz  ,ix+2
                izp1_ixp2= i+1  +2*nz  !iz+1,ix+2
                

                sf_dvx_dx = c1x*(sf_vx(iz_ixp1)-sf_vx(iz_ix)) +c2x*(sf_vx(iz_ixp2)-sf_vx(iz_ixm1))
                sf_dvz_dz = c1z*(sf_vz(izp1_ix)-sf_vz(iz_ix)) +c2z*(sf_vz(izp2_ix)-sf_vz(izm1_ix))

                grad_lda(j)=grad_lda(j) + rf_szz(i)*(sf_dvz_dz + sf_dvx_dx) &
                                        + rf_sxx(i)*(sf_dvz_dz + sf_dvx_dx)


                sf_4_dvzdx_p_dvxdz = &   !(dvz_dx+dvx_dz)(iz,ix)     [iz-0.5,ix-0.5]
                                         c1x*(sf_vz(iz_ix    )-sf_vz(iz_ixm1  )) +c2x*(sf_vz(iz_ixp1  )-sf_vz(iz_ixm2  )) &
                                       + c1z*(sf_vx(iz_ix    )-sf_vx(izm1_ix  )) +c2z*(sf_vx(izp1_ix  )-sf_vx(izm2_ix  )) &
                                     & &
                                     & & !(dvz_dx+dvx_dz)(iz,ix+1)   [iz-0.5,ix+0.5]
                                       + c1x*(sf_vz(iz_ixp1  )-sf_vz(iz_ix    )) +c2x*(sf_vz(iz_ixp2  )-sf_vz(iz_ixm1  )) &
                                       + c1z*(sf_vx(iz_ixp1  )-sf_vx(izm1_ixp1)) +c2z*(sf_vx(izp1_ixp1)-sf_vx(izm2_ixp1)) &
                                     & &
                                     & & !(dvz_dx+dvx_dz)(iz+1,ix)   [iz+0.5,ix-0.5]
                                       + c1x*(sf_vz(izp1_ix  )-sf_vz(izp1_ixm1)) +c2x*(sf_vz(izp1_ixp1)-sf_vz(izp1_ixm2)) &
                                       + c1z*(sf_vx(izp1_ix  )-sf_vx(iz_ix    )) +c2z*(sf_vx(izp2_ix  )-sf_vx(izm1_ix  )) &
                                     & &
                                     & & !(dvz_dx+dvx_dz)(iz+1,ix+1) [iz+0.5,ix+0.5]
                                       + c1x*(sf_vz(izp1_ixp1)-sf_vz(izp1_ix  )) +c2x*(sf_vz(izp1_ixp2)-sf_vz(izp1_ixm1)) &
                                       + c1z*(sf_vx(izp1_ixp1)-sf_vx(iz_ixp1  )) +c2z*(sf_vx(izp2_ixp1)-sf_vx(izm1_ixp1))

                        ![iz-0.5,ix-0.5]   [iz+0.5,ix-0.5]   [iz-0.5,ix+0.5]     [iz+0.5,ix+0.5]
                rf_4_szx = rf_szx(iz_ix) + rf_szx(izp1_ix) + rf_szx(iz_ixp1) + rf_szx(izp1_ixp1)

                grad_mu(j)=grad_mu(j) + rf_szz(i)*( ldap2mu(i)*sf_dvz_dz -    lda(i)*sf_dvx_dx) &
                                      + rf_sxx(i)*(-lda(i)    *sf_dvz_dz +ldap2mu(i)*sf_dvx_dx) &
                                      + two_ldapmu(i)*0.0625*rf_4_szx*sf_4_dvzdx_p_dvxdz   !0.0625=1/16

            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine grad2d_freesurf_glda_gmu(rf_sxx,sf_vx,         &
                                         grad_lda,grad_mu,    &
                                         ifx,ilx)
        real,dimension(*) :: rf_sxx,sf_vx
        real,dimension(*) :: grad_lda,grad_mu
        
        nz=cb%nz
        
        sf_dvx_dx=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         iz_ix,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         sf_dvx_dx)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !!!dir$ simd
            iz=1
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !corr has no boundary layers       
                
                iz_ix  = i    !iz  ,ix
                iz_ixm1  = i      -nz  !iz  ,ix-1
                iz_ixp1  = i      +nz  !iz  ,ix+1
                iz_ixp2  = i    +2*nz  !iz  ,ix+2

                sf_dvx_dx = c1x*(sf_vx(iz_ixp1)-sf_vx(iz_ix)) +c2x*(sf_vx(iz_ixp2)-sf_vx(iz_ixm1))

                grad_lda(j)=grad_lda(j) + rf_sxx(i)*sf_dvx_dx

                grad_mu(j) =grad_mu(j)  + rf_sxx(i)*sf_dvx_dx
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine
    
    
    subroutine grad2d_density(rf_vz,rf_vx,         &
                              sf_szz,sf_sxx,sf_szx,&
                              grad,                &
                              ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_vz,rf_vx
        real,dimension(*) :: sf_szz,sf_sxx,sf_szx
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        sf_2_dszzdz_p_dszxdx=0.
        sf_2_dsxxdx_p_dszxdz=0.
        rf_2vz=0.
        rf_2vx=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,izp1_ixm1,&
        !$omp         izm1_ixp1,iz_ixp1,izp1_ixp1,izp2_ixp1,&
        !$omp         iz_ixp2,izp1_ixp2,&
        !$omp         sf_2_dszzdz_p_dszxdx,sf_2_dsxxdx_p_dszxdz,&
        !$omp         rf_2vz,rf_2vx)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !corr has no boundary layers
                
                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2  = i    -2*nz  !iz  ,ix-2
                iz_ixm1  = i      -nz  !iz  ,ix-1
                izp1_ixm1= i+1    -nz  !iz-1,ix-1

                izm1_ixp1= i-1    +nz  !iz-1,ix+1
                iz_ixp1  = i      +nz  !iz  ,ix+1
                izp1_ixp1= i+1    +nz  !iz+1,ix+1
                izp2_ixp1= i+2    +nz  !iz+2,ix+1

                iz_ixp2  = i    +2*nz  !iz ,ix+2
                izp1_ixp2= i+1  +2*nz  !iz+1,ix+2

                sf_2_dszzdz_p_dszxdx = &   !(dszx_dx+dszz_dz)(iz,ix)   [iz-0.5,ix]
                                           c1z*(sf_szz(iz_ix    )-sf_szz(izm1_ix)) +c2z*(sf_szz(izp1_ix  )-sf_szz(izm2_ix  )) &
                                         + c1x*(sf_szx(iz_ixp1  )-sf_szx(iz_ix  )) +c2x*(sf_szx(iz_ixp2  )-sf_szx(iz_ixm1  )) &
                                       & &
                                       & & !(dszx_dx+dszz_dz)(iz+1,ix) [iz+0.5,ix]
                                         + c1z*(sf_szz(izp1_ix  )-sf_szz(iz_ix  )) +c2z*(sf_szz(izp2_ix  )-sf_szz(izm1_ix  )) &
                                         + c1x*(sf_szx(izp1_ixp1)-sf_szx(izp1_ix)) +c2x*(sf_szx(izp1_ixp2)-sf_szx(izp1_ixm1))
                
                sf_2_dsxxdx_p_dszxdz = &   !(dsxx_dx+dszx_dz)(iz,ix)   [iz,ix-0.5]
                                           c1x*(sf_sxx(iz_ix    )-sf_sxx(iz_ixm1)) +c2x*(sf_sxx(iz_ixp1  )-sf_sxx(iz_ixm2  )) &
                                         + c1z*(sf_szx(izp1_ix  )-sf_szx(iz_ix  )) +c2z*(sf_szx(izp2_ix  )-sf_szx(izm1_ix  )) &
                                       & &
                                       & & !(dsxx_dx+dszx_dz)(iz,ix+1) [iz,ix+0.5]
                                         + c1x*(sf_sxx(iz_ixp1  )-sf_sxx(iz_ix  )) +c2x*(sf_sxx(iz_ixp2  )-sf_sxx(iz_ixm1  )) &
                                         + c1z*(sf_szx(izp1_ixp1)-sf_szx(iz_ixp1)) +c2z*(sf_szx(izp2_ixp1)-sf_szx(izm1_ixp1))
                                         
                rf_2vz = rf_vz(iz_ix) + rf_vz(izp1_ix)
                         ![iz-0.5,ix]      [iz+0.5,ix]
                rf_2vx = rf_vx(iz_ix) + rf_vx(iz_ixp1)
                         ![iz,ix-0.5]      [iz,ix+0.5]
                
                grad(j)=grad(j) + 0.25*( rf_2vz*sf_2_dszzdz_p_dszxdx + rf_2vx*sf_2_dsxxdx_p_dszxdz )
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine

    subroutine imag2d_xcorr(rf_szz,rf_sxx, &
                            sf_szz,sf_sxx, &
                            imag,          &
                            ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_szz,rf_sxx
        real,dimension(*) :: sf_szz,sf_sxx
        real,dimension(*) :: imag
        
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
                
                rp = sn_p*( rf_szz(i) + rf_sxx(i) )
                sp = sn_p*( sf_szz(i) + sf_sxx(i) )
                
                imag(j)=imag(j) + rp*sp !+rf_vz(i)*sf_vz(i) + rf_vx(i)*sf_vx(i)
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine engy2d_xcorr(sf_szz,sf_sxx, &
                            engy,          &
                            ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_szz,sf_sxx
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
                
                sp = sn_p* ( sf_szz(i) + sf_sxx(i) )
                
                engy(j)=engy(j) + sp*sp
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

end
