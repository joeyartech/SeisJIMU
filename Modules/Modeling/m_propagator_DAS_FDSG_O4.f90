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

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain ISOtropic 2D P-SV (ELastic) propagation'//s_NL// &
            '1st-order Momemtum-Strain formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.606 *Vmax/dx'//s_NL// &
            'Required model attributes: vp, vs, rho'//s_NL// &
            'Required field components: pz,px,'//s_NL// &
            'Required field components: ez,ex,es'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Basic gradients: glda gmu'

        integer :: nbndlayer=max(2,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=2 !number of basic gradients

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:),allocatable :: buoz, buox, ldap2mu, lda, mu
        real,dimension(:,:),allocatable :: two_ldapmu, ldapmu!, inv_ladpmu_4mu

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
        
        procedure :: inject_momenta
        procedure :: inject_strains
        procedure :: update_momenta
        procedure :: update_strains
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    !conversion from shear strain to es
    real,parameter :: ezx_es=2

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

        call alloc(self%buoz,   [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%buox,   [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%ldap2mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%lda,    [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%mu,     [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%ldapmu, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        call alloc(temp_mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        self%ldap2mu(:,:)=cb%rho(:,:,1)*cb%vp(:,:,1)**2
             temp_mu(:,:)=cb%rho(:,:,1)*cb%vs(:,:,1)**2

        self%lda=self%ldap2mu-2.*temp_mu
        ! if(mpiworld%is_master) then
        ! write(*,*) 'self%ldap2mu sanity:', minval(self%ldap2mu),maxval(self%ldap2mu)
        ! write(*,*) 'self%lda     sanity:', minval(self%lda),maxval(self%lda)
        ! endif

        self%ldapmu=(self%lda+temp_mu)!why?

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
        call f%init_boundary_momenta

        call alloc(f%pz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%px, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%ez, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%ex, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%es, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])

        call alloc(f%dpz_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dpz_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dpx_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dpx_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dez_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dez_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dex_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dex_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%des_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%des_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        
    end subroutine

    subroutine init_correlate(self,corr,name)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name
        
        corr%name=name

        if(name(1:1)=='g') then !gradient components
            call alloc(corr%grho,m%nz,m%nx,m%ny)
            call alloc(corr%glda,m%nz,m%nx,m%ny)
            call alloc(corr%gmu, m%nz,m%nx,m%ny)
        ! else !image components
        !     call alloc(corr%ipp,m%nz,m%nx,m%ny)
        !     call alloc(corr%ibksc,m%nz,m%nx,m%ny)
        !     call alloc(corr%ifwsc,m%nz,m%nx,m%ny)
        endif

    end subroutine

    subroutine init_abslayer(self)
        class(t_propagator) :: self

        call cpml%init

    end subroutine

    subroutine assemble(self,corr)
        class(t_propagator) :: self
        type(t_correlate) :: corr

        if(allocated(correlate_gradient)) then
!            call correlate_assemble(corr%grho, correlate_gradient(:,:,:,1))
            call correlate_assemble(corr%glda, correlate_gradient(:,:,:,2))
            call correlate_assemble(corr%gmu,  correlate_gradient(:,:,:,3))
        endif        
        
    end subroutine
    
    !========= Derivations =================
    !PDE:      A u = M ∂ₜ u - MD Mu = Mf
    !Adjoint:  Aᵀa = M ∂ₜᵀa - MDᵀMa = Md
    !where
    !u=[pz px ez ex es]ᵀ, [pz px] are momenta, [ez ex] are normal strains, es is shear strain
    !f=[fz fx sz sx ss]ᵀδ(x-xs) with xs source position, d is recorded data
    !  [b              ]    [0  0  ∂z 0  ∂ₓ]
    !  |  b            |    |0  0  0  ∂ₓ ∂z|
    !M=|    λ+2μ  λ    |, D=|∂z 0  0  0  0 |
    !  |     λ   λ+2μ  |    |0  ∂ₓ 0  0  0 |
    !  [              μ]    [∂ₓ ∂z 0  0  0 ]
    !a=[pzᵃ pxᵃ ezᵃ, exᵃ, esᵃ]ᵀ is the adjoint field
    !
    !Continuous case:
    !<a|Au-Mf> = ∫ a(x,t) (M∂ₜ-MD)u(x,t) -Mf dx³dt
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
    !                         |    es    |   es -½ pz  es    |   es    |
    !                         |    μ     |   μ     bz  μ     |   μ     |
    !                         |          |         |         |         |
    !                        λ,μ   bx   λ,μ  bx   λ,μ  bx   λ,μ  bx   λ,μ
    !  -pz--e-pz-e-pz-→ t    -sn---px---sn---px---en---px---en---px---sn-→ x
    !   -1 -½  0 ½ 1          -2  -1½   -1   -½    0    ½    1   1½    2    
    !                         |          |         |         |         | 
    !                         |    es    |   es  ½ pz  es    pz  es    | 
    !                         |    μ     |   μ     bz  μ     bz  μ     | 
    !                         |          |         |         |         | 
    !                        -|----------|-------1-en--------en--------|-
    !                         |          |        λ,μ       λ,μ        | 
    !                         |          |         |         |         | 
    !                         |          |      1½ pz        |         | 
    !                         |          |         bz        |         | 
    !                         |          |         |         |         | 
    !                                            z ↓
    !
    ! en: ez or ex (normal strains)
    ! es: shear strains
    ! λ,μ: λ and λ+2μ
    !
    !Convention for half-integer index:
    !(array index)  =>     (real index)     
    !  pz(iz,ix)    =>   pz[iz-½,ix  ]^n   :=pz((iz-½)*dz,ix*dx,n*dt)
    !  px(iz,ix)    =>   px[iz,  ix-½]^n  
    !  en(iz,ix)    =>   en[iz,  ix  ]^n+½ :=en(iz*dz,    ix*dx,(n+½)*dt)
    !  es(iz,ix)    =>   es[iz-½,ix-½]^n  
    !
    !Forward:
    !FD eqn:
    !    [pz^n  ]   [ 0     0    ∂zᵇ  0  ∂ₓᶠ] [pz^n+1]      [∂zᵇ(λ+2μ)ez^n+½ + ∂zᵇ λez^n+½ + ∂ₓᶠ μes^n+½]   
    !    |px^n  |   | 0     0     0  ∂ₓᵇ ∂zᶠ| |px^n+1|      |∂ₓᵇ(λ+2μ)ex^n+½ + ∂ₓᵇ λex^n+½ + ∂zᶠ μes^n+½|   
    !∂ₜᶠ |ez^n+1| = |∂zᶠ    0     0   0   0 |M|ez^n+½| +f = |∂zᶠ bpz^n+1                                | +f
    !    |ex^n+1|   | 0    ∂ₓᶠ    0   0   0 | |ex^n+½|      |∂ₓᶠ bpx^n+1                                |   
    !    [es^n+1]   [∂ₓᵇ   ∂zᵇ    0   0   0 ] [es^n+½]      [∂ₓᵇ bpz^n+1 + ∂zᵇ bpx^n+1                  ]   
    !where
    !∂ₜᶠ*dt := pz^n+1 - pz^n                             ~O(t²)
    !∂zᵇ*dz := c₁(ez(iz  )-ez(iz-1) +c₂(ez(iz+1)-ez(iz-2)  ~O(x⁴)
    !∂zᶠ*dz := c₁(pz(iz+1)-pz(iz  ) +c₂(pz(iz+2)-pz(iz-1)  ~O(x⁴)
    !
    !Time marching:
    ![pz^n+1 ]   [pz^n  ]   [∂zᵇ(λ+2μ)ez^n+½ + ∂zᵇ λez^n+½ + ∂ₓᶠ μes^n+½]
    !|px^n+1 |   |px^n  |   |∂ₓᵇ(λ+2μ)ex^n+½ + ∂ₓᵇ λex^n+½ + ∂zᶠ μes^n+½|
    !|ez^n+1½| = |ez^n+½| + |∂zᶠ bpz^n+1                                |dt +f*dt
    !|ex^n+1½|   |ex^n+½|   |∂ₓᶠ bpx^n+1                                |
    ![es^n+1½]   [es^n+½]   [∂ₓᵇ bpz^n+1 + ∂zᵇ bpx^n+1                  ]
    !Step #1: pz^n += src
    !Step #2: pz^n+1 = pz^n + spatial FD(e^n+½)
    !Step #3: e^n+½ += src
    !Step #4: e^n+1½ = e^n+½ + spatial FD(pz,px^n+1)
    !Step #5: sample pz,px^n & e^n+1½ at receivers
    !Step #6: save pz,px^n+1 to boundary values
    !
    !Reverse time marching (for wavefield reconstruction)
    ![es^n+½]   [es^n+1½]   [∂zᵇ(λ+2μ)ez^n+½ + ∂zᵇ λez^n+½ + ∂ₓᶠ μes^n+½]
    !|ex^n+½|   |ex^n+1½|   |∂ₓᵇ(λ+2μ)ex^n+½ + ∂ₓᵇ λex^n+½ + ∂zᶠ μes^n+½|
    !|ez^n+½| = |ez^n+1½| - |∂zᶠ bpz^n+1                                |dt -f*dt
    !|px^n  |   |px^n+1 |   |∂ₓᶠ bpx^n+1                                |
    ![pz^n  ]   [pz^n+1 ]   [∂ₓᵇ bpz^n+1 + ∂zᵇ bpx^n+1                  ]
    !Step #6: load boundary values for pz,px^n+1
    !Step #4: e^n+½ = e^n+1½ - spatial FD(pz,px^n+1)
    !Step #3: e^n+½ -= src
    !Step #2: pz^n+1 = pz^n - spatial FD(e^n+½)
    !Step #1: pz^n -= src
    !N.B. Same codes for spatial FDs as in forward time marching, with a negated dt.
    !
    !Adjoint:
    !FD eqn:
    !    [pzᵃ^n  ]   [ 0     0    ∂zᶠᵀ  0   ∂ₓᵇᵀ] [pz^n+1]
    !    |pxᵃ^n  |   | 0     0     0   ∂ₓᶠᵀ ∂zᵇᵀ| |px^n+1|
    !∂ₜᶠᵀ|ezᵃ^n+1| = |∂zᵇᵀ   0     0    0    0  |M|ez^n+½| +d
    !    |exᵃ^n+1|   | 0    ∂ₓᵇᵀ   0    0    0  | |ex^n+½|
    !    [esᵃ^n+1]   [∂ₓᶠᵀ  ∂zᶠᵀ   0    0    0  ] [es^n+½]
    !where
    !∂ₜᶠᵀ = pz^n-1 -pz^n   = -∂ₜᵇ
    !∂zᵇᵀ = c₁(ez[i  ]-ez[i+½]) +c₂(ez[i- ½]-ez[i+1½]) = -∂zᶠ
    !∂zᶠᵀ = c₁(pz[i-½]-pz[i  ]) +c₂(pz[i-1½]-pz[i+ ½]) = -∂zᵇ
    !    [pzᵃ^n  ]   [ 0    0   ∂zᵇ  0  ∂ₓᶠ] [pzᵃ^n+1]      [∂zᵇ(λ+2μ)ez^n+½ + ∂zᵇ λez^n+½ + ∂ₓᶠ μes^n+½]   
    !    |pxᵃ^n  |   | 0    0    0  ∂ₓᵇ ∂zᶠ| |pxᵃ^n+1|      |∂ₓᵇ(λ+2μ)ex^n+½ + ∂ₓᵇ λex^n+½ + ∂zᶠ μes^n+½|   
    !∂ₜᵇ |ezᵃ^n+1| = |∂zᶠ   0    0   0   0 |M|ezᵃ^n+½| -d = |∂zᶠ bpz^n+1                                | -d
    !    |exᵃ^n+1|   | 0   ∂ₓᶠ   0   0   0 | |exᵃ^n+½|      |∂ₓᶠ bpx^n+1                                |   
    !    [esᵃ^n+1]   [∂ₓᵇ  ∂zᵇ   0   0   0 ] [esᵃ^n+½]      [∂ₓᵇ bpz^n+1 + ∂zᵇ bpx^n+1                  ]   
    !ie. Dᵀ=-D, antisymmetric
    !
    !Time marching:
    ![pzᵃ^n+1 ]   [pzᵃ^n  ]   [∂zᵇ(λ+2μ)ez^n+½ + ∂zᵇ λez^n+½ + ∂ₓᶠ μes^n+½]        
    !|pxᵃ^n+1 |   |pxᵃ^n  |   |∂ₓᵇ(λ+2μ)ex^n+½ + ∂ₓᵇ λex^n+½ + ∂zᶠ μes^n+½|        
    !|ezᵃ^n+1½| = |ezᵃ^n+½| + |∂zᶠ bpz^n+1                                |dt -d*dt
    !|exᵃ^n+1½|   |exᵃ^n+½|   |∂ₓᶠ bpx^n+1                                |        
    ![esᵃ^n+1½]   [esᵃ^n+½]   [∂ₓᵇ bpz^n+1 + ∂zᵇ bpx^n+1                  ]        
    !but we have to do it in reverse time:
    ![esᵃ^n+½] = [esᵃ^n+1½]   [∂zᵇ(λ+2μ)ez^n+½ + ∂zᵇ λez^n+½ + ∂ₓᶠ μes^n+½]        
    !|exᵃ^n+½| = |exᵃ^n+1½|   |∂ₓᵇ(λ+2μ)ex^n+½ + ∂ₓᵇ λex^n+½ + ∂zᶠ μes^n+½|        
    !|ezᵃ^n+½| = |ezᵃ^n+1½| - |∂zᶠ bpz^n+1                                |dt +d*dt
    !|pxᵃ^n  | = |pxᵃ^n+1 |   |∂ₓᶠ bpx^n+1                                |        
    ![pzᵃ^n  ] = [pzᵃ^n+1 ]   [∂ₓᵇ bpz^n+1 + ∂zᵇ bpx^n+1                  ]        
    !Step #5: eᵃ^n+1½ += adjsrc
    !Step #4: eᵃ^n+½ = eᵃ^n+1½ - spatial FD(pzᵃ,pxᵃ^n+1)
    !Step #3: pzᵃ^n+1 += adjsrc
    !Step #2: pzᵃ^n = pzᵃ^n+1 - spatial FD(eᵃ^n+½)
    !N.B. Same codes for spatial FDs as in forward time marching, with a negated dt, but the RHS should use a "+" sign (regardless of the reverse time direction).
    !
    !For adjoint test:
    !In each step of forward time marching: dsyn=RGAf
    !f:source wavelet, A:inject source into field, G:propagator, R:extract field at receivers
    !while in each step of reverse-time adjoint marching: dadj=AᵀGᵀRᵀdsyn=AᵀGᵀRᵀdsyn
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
                call fld_u%check_value(fld_u%pz)
            endif

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add forces to pz^it
            call cpu_time(tic)
            call self%inject_momenta(fld_u,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2: from pz^it to pz^it+1 by differences of e^it+0.5
            call cpu_time(tic)
            call self%update_momenta(fld_u,time_dir,it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            !step 3: add pressure to e^it+0.5
            call cpu_time(tic)
            call self%inject_strains(fld_u,time_dir,it)
            call cpu_time(toc)
            tt3=tt3+toc-tic

            !step 4: from e^it+0.5 to e^it+1.5 by differences of pz^it+1
            call cpu_time(tic)
            call self%update_strains(fld_u,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !step 5: sample pz^it+1 or e^it+1.5 at receivers
            call cpu_time(tic)
            call self%extract(fld_u,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !snapshot
            call fld_u%write(it)

            !step 6: save pz^it+1 in boundary layers
            ! if(fld_u%if_will_reconstruct) then
                call cpu_time(tic)
                call fld_u%boundary_transport_momenta('save',it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            ! endif

        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add source momenta',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update momenta    ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to add source strains  ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update strains      ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field        ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary        ',tt6/mpiworld%max_threads
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
                call fld_a%check_value(fld_a%pz)
                call fld_u%check_value(fld_u%pz)
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then

                !backward step 6: retrieve pz^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_u%boundary_transport_momenta('load',it)
                call cpu_time(toc)
                tt1=tt1+toc-tic
                
                !backward step 4: e^it+1.5 -> e^it+0.5 by FD of pz^it+1
                call cpu_time(tic)
                call self%update_strains(fld_u,time_dir,it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: rm pressure from e^it+0.5
                call cpu_time(tic)
                call self%inject_strains(fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to e^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_strains(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: e^it+1.5 -> e^it+0.5 by FD^T of pz^it+1
            call cpu_time(tic)
            call self%update_strains(fld_a,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: rf%e^it+0.5 star D sf%s_dt^it+0.5
            !use sf%pz^it+1 to compute sf%s_dt^it+0.5, as backward step 4
            if(mod(it,irdt)==0) then
                call cpu_time(tic)
                call cross_correlate_glda_gmu(fld_a,fld_u,a_star_u,it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif
                            
            !========================================================!

            ! if(present(o_sf)) then
                !backward step 2: pz^it+1 -> pz^it by FD of e^it+0.5
                call cpu_time(tic)
                call self%update_momenta(fld_u,time_dir,it)
                call cpu_time(toc)
                tt7=tt7+toc-tic

                !backward step 1: rm forces from pz^it
                call cpu_time(tic)
                call self%inject_momenta(fld_u,time_dir,it)
                call cpu_time(toc)
                tt8=tt8+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 3: inject to pz^it+1 at receivers
            call cpu_time(tic)
            call self%inject_momenta(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: pz^it+1 -> pz^it by FD^T of e^it+0.5
            call cpu_time(tic)
            call self%update_momenta(fld_a,time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            !adjoint step 1: sample pz^it or e^it+0.5 at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_a,it)
                call cpu_time(toc)
                tt11=tt11+toc-tic
            endif
            
            ! !grho: sfield%pz_dt^it \dot rfield%pz^it
            ! !use sfield%e^it+0.5 to compute sfield%pz_dt^it, as backward step 2
            ! if(if_compute_grad.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call gradient_density(fld_a,fld_u,it,cb%grad(:,:,1,1))
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
            
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
            write(*,*) 'Elapsed time to update strains           ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source strains        ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update momenta           ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source momenta        ',tt8/mpiworld%max_threads
            write(*,*) 'Elapsed time ----------------------------'
            write(*,*) 'Elapsed time to add adjsource strains    ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj strains       ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to add adjsource momenta    ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj momenta       ',tt10/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract&write fields     ',tt11/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate                ',tt6/mpiworld%max_threads

        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_rfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_rfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')
        call hud('ximage < snap_*  n1='//num2str(cb%mz)//' perc=99')
        call hud('xmovie < snap_*  n1='//num2str(cb%mz)//' n2='//num2str(cb%mx)//' clip=?e-?? loop=2 title=%g')
        
    end subroutine


    !forward: add RHS to pz^it
    !adjoint: add RHS to pz^it+1
    subroutine inject_momenta(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        if(.not. f%is_adjoint) then

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler
            
            if(if_hicks) then
                select case (shot%src%comp)
                    case ('pz')
                    !f%pz(ifz:ilz,ifx:ilx,1) = f%pz(ifz:ilz,ifx:ilx,1) + wl*shot%src%interp_coef(:,:,1)
                    f%pz(ifz:ilz,ifx:ilx,1) = f%pz(ifz:ilz,ifx:ilx,1) + wl/self%buoz(ifz:ilz,ifx:ilx) *shot%src%interp_coef(:,:,1)
                    
                    case ('px')
                    !f%px(ifz:ilz,ifx:ilx,1) = f%px(ifz:ilz,ifx:ilx,1) + wl*shot%src%interp_coef(:,:,1)
                    f%px(ifz:ilz,ifx:ilx,1) = f%px(ifz:ilz,ifx:ilx,1) + wl/self%buox(ifz:ilz,ifx:ilx) *shot%src%interp_coef(:,:,1)
                    
                end select
                
            else
                select case (shot%src%comp)
                    case ('pz') !vertical force     on pz[iz-0.5,ix]
                    !f%pz(iz,ix,1) = f%pz(iz,ix,1) + wl
                    f%pz(iz,ix,1) = f%pz(iz,ix,1) + wl/self%buoz(iz,ix)
                    
                    case ('px') !horizontal x force on px[iz,ix-0.5]
                    !f%px(iz,ix,1) = f%px(iz,ix,1) + wl
                    f%px(iz,ix,1) = f%px(iz,ix,1) + wl/self%buox(iz,ix)
                    
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
                        case ('pz') !vertical z adjsource
                        !f%pz(ifz:ilz,ifx:ilx,1) = f%pz(ifz:ilz,ifx:ilx,1) + wl*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                        f%pz(ifz:ilz,ifx:ilx,1) = f%pz(ifz:ilz,ifx:ilx,1) + wl/self%buoz(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!

                        case ('px') !horizontal x adjsource
                        !f%px(ifz:ilz,ifx:ilx,1) = f%px(ifz:ilz,ifx:ilx,1) + wl*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                        f%px(ifz:ilz,ifx:ilx,1) = f%px(ifz:ilz,ifx:ilx,1) + wl/self%buox(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                        
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                        case ('pz') !vertical z adjsource
                        !pz[ix,1,iz-0.5]
                        !f%pz(iz,ix,1) = f%pz(iz,ix,1) + wl !no time_dir needed!
                        f%pz(iz,ix,1) = f%pz(iz,ix,1) + wl/self%buoz(iz,ix)  !self%buoz(iz,ix) !no time_dir needed!

                        case ('px') !horizontal x adjsource
                        !px[ix-0.5,1,iz]
                        !f%px(iz,ix,1) = f%px(iz,ix,1) + wl !no time_dir needed!
                        f%px(iz,ix,1) = f%px(iz,ix,1) + wl/self%buox(iz,ix) !no time_dir needed!
                        
                    end select
                    
                endif
                
            enddo
        
    end subroutine
    
    !forward: pz^it -> pz^it+1 by FD  of e^it+0.5
    !adjoint: pz^it+1 -> pz^it by FDᵀ of e^it+0.5
    subroutine update_momenta(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ifz=f%bloom(1,it)+2
        ilz=f%bloom(2,it)-2  !-1
        ifx=f%bloom(3,it)+2
        ilx=f%bloom(4,it)-2  !-1

        if(m%is_freesurface) ifz=max(ifz,1)

        call fd2d_momenta(f%pz,f%px,f%ez,f%ex,f%es,            &
                        f%dez_dz,f%dex_dx,f%dex_dz,f%dez_dx,f%des_dz,f%des_dx,&
                        self%ldap2mu,self%lda,self%mu,             &
                        ifz,ilz,ifx,ilx,time_dir*self%dt)

        !apply free surface boundary condition if needed
        !free surface is located at [1,ix,1] level
        !Δszz=0 : -(λ+2μ)∂_z vz = λ ∂ₓvx
        !         -(λ+2μ)[1,ix]*(vz[1.5,ix]-vz[0.5,ix])/dz = λ[1,ix]*(vx[1,ix-0.5]-vx[1,ix+0.5])/dx
        !         -(λ+2μ)(1,ix)*(vz(2  ,ix)-vz(1  ,ix))/dz = λ(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
        !          (λ+2μ)(1,ix)*(vz(1  ,ix)-vz(2  ,ix))/dz = λ(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
        !Δszx=0 : -∂_z vx = ∂ₓvz
        !         -(vx[2,ix+0.5]-vx[0,ix+0.5])/2dz = ( (vz[0.5,ix]+vz[1.5,ix])/2 - (vz[0.5,ix+1]+vz[1.5,ix+1])/2 )/dx
        !         -(vx(2,ix+1  )-vx(0,ix+1  ))/2dz = ( (vz(1  ,ix)+vz(2  ,ix))/2 - (vz(1  ,ix+1)+vz(2  ,ix+1))/2 )/dx
        !          (vx(0,ix+1  )-vx(2,ix+1  ))/ dz = ( (vz(1  ,ix)+vz(2  ,ix))   -  vz(1  ,ix+1)-vz(2  ,ix+1)    )/dx
        if(m%is_freesurface) then
            dz_dx = m%dz/m%dx

            do ix=ifx,ilx
                f%pz(1,ix,1)= f%pz(2,ix,1) + 1/cb%rho(1,ix,1)* &
                self%lda(1,ix)*(self%buox(1,ix)*f%px(1,ix,1)-self%buox(1,ix+1)*f%px(1,ix+1,1))*dz_dx/self%ldap2mu(1,ix)
                !f%px(0,ix)= f%px(2,ix) + (f%pz(1,ix)+f%pz(2,ix)-f%pz(1,ix+1)-f%pz(2,ix+1))*dz_dx/lm%ldap2mu(1,ix) !toy2del says this condition is not needed
            enddo
        endif

    end subroutine
    
    !forward: add RHS to e^it+0.5
    !adjoint: add RHS to e^it+1.5
    subroutine inject_strains(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        if(.not. f%is_adjoint) then

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler
            
            if(if_hicks) then
                select case (shot%src%comp)
                    case ('ez')
                    f%ez(ifz:ilz,ifx:ilx,1) = f%ez(ifz:ilz,ifx:ilx,1) + wl/self%ldapmu(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)
                    case ('ex')
                    f%ex(ifz:ilz,ifx:ilx,1) = f%ex(ifz:ilz,ifx:ilx,1) + wl/self%ldapmu(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)
                    case ('es')
                    f%es(ifz:ilz,ifx:ilx,1) = f%es(ifz:ilz,ifx:ilx,1) + wl/self%mu(ifz:ilz,ifx:ilx)        *shot%src%interp_coef(:,:,1)
                endselect
                
            else
                select case (shot%src%comp)
                    case ('ez')
                    !f%ez(iz,ix,1) = f%ez(iz,ix,1) + wl
                    f%ez(iz,ix,1) = f%ez(iz,ix,1) + wl/(self%lda(iz,ix)+self%mu(iz,ix))*2!/self%two_ldapmu(iz,ix)
                    case ('ex')
                    f%ex(iz,ix,1) = f%ex(iz,ix,1) + wl/self%ldapmu(iz,ix)
                    case ('es')
                    f%es(iz,ix,1) = f%es(iz,ix,1) + wl/self%mu(iz,ix)
                endselect
                
            endif

            return

        endif

            do i=1,shot%nrcv
                ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                
                wl=f%wavelet(i,it)*wavelet_scaler
                    
                if(if_hicks) then
                    select case (shot%rcv(i)%comp)
                        case ('ez')
                        f%ez(ifz:ilz,ifx:ilx,1) = f%ez(ifz:ilz,ifx:ilx,1) +wl/self%ldapmu(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                        case ('ex')
                        f%ex(ifz:ilz,ifx:ilx,1) = f%ex(ifz:ilz,ifx:ilx,1) +wl/self%ldapmu(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                        case ('es')
                        f%es(ifz:ilz,ifx:ilx,1) = f%es(ifz:ilz,ifx:ilx,1) +wl/self%mu(ifz:ilz,ifx:ilx)        *shot%rcv(i)%interp_coef(:,:,1)
                    endselect
                else           
                    select case (shot%rcv(i)%comp)
                        case ('ez')
                        !f%ez(iz,ix,1) = f%ez(iz,ix,1) +wl
                        f%ez(iz,ix,1) = f%ez(iz,ix,1) +wl/self%ldapmu(iz,ix)
                        case ('ex')
                        f%ex(iz,ix,1) = f%ex(iz,ix,1) +wl/self%ldapmu(iz,ix)
                        case ('es')
                        f%es(iz,ix,1) = f%es(iz,ix,1) +wl/self%mu(iz,ix)
                    endselect
                endif

            enddo
        
    end subroutine

    !forward: e^it+0.5 -> e^it+1.5 by FD of pz^it+1
    !adjoint: e^it+1.5 -> e^it+0.5 by FD^T of pz^it+1
    subroutine update_strains(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ifz=f%bloom(1,it)+2  !1
        ilz=f%bloom(2,it)-2
        ifx=f%bloom(3,it)+2  !1
        ilx=f%bloom(4,it)-2
        
        if(m%is_freesurface) ifz=max(ifz,1)

        call fd2d_strains(f%pz,f%px,f%ez,f%ex,f%es,    &
                           f%dpz_dz,f%dpx_dx,f%dpz_dx,f%dpx_dz,&
                           self%buoz,self%buox,&
                           ifz,ilz,ifx,ilx,time_dir*self%dt)
        
        !apply free surface boundary condition if needed
        !free surface is located at [1,ix] level
        !so explicit boundary condition: 0=sz(1,ix)=ldap2mu*ez(1,ix)+lda*ex(1,ix) -> ez(1,ix)=-lda(1,ix)*ex(1,ix)/ldap2mu(1,ix)
        !and antisymmetric mirroring: mu*es[0.5,ix-0.5]=-mu*es[1.5,ix-0.5] -> mu(1,ix)*es(1,ix)=-mu(2,ix)*es(2,ix)
        if(m%is_freesurface) then
            f%ez(1,:,1)=-self%lda(1,:)*f%ex(1,:,1)/self%ldap2mu(1,:)
            f%es(1,:,1)=-self%mu(2,:) *f%es(2,:,1)/self%mu(1,:)
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
                        
                        case ('ez')
                        f%seismo(i,it)=sum(f%ez(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                        case ('ex')
                        f%seismo(i,it)=sum(f%ex(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                        case ('es')
                        f%seismo(i,it)=sum(f%es(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                        case ('pz')
                        f%seismo(i,it)=sum(f%pz(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                        case ('px')
                        f%seismo(i,it)=sum(f%px(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                        case ('ez') !ez[iz,ix]
                        f%seismo(i,it)=f%ez(iz,ix,1)
                        case ('ex') !ez[iz,ix]
                        f%seismo(i,it)=f%ex(iz,ix,1)
                        case ('es') !ez[iz,ix]
                        f%seismo(i,it)=f%es(iz,ix,1)
                        case ('pz') !pz[iz-0.5,ix]
                        f%seismo(i,it)=f%pz(iz,ix,1)
                        case ('px') !px[iz,ix-0.5]
                        f%seismo(i,it)=f%px(iz,ix,1)
                    end select
                    
                endif

            enddo

            return

        endif

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            if(if_hicks) then
                select case (shot%src%comp)
                    case ('ez')
                    f%seismo(1,it)=sum(f%ez(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))

                    case ('ex')
                    f%seismo(1,it)=sum(f%ex(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))

                    case ('es')
                    f%seismo(1,it)=sum(f%es(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                    
                    case ('pz')
                    f%seismo(1,it)=sum(f%pz(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                    
                    case ('px')
                    f%seismo(1,it)=sum(f%px(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                    
                end select
                
            else
                select case (shot%src%comp)
                    case ('ez') !ez[iz-0.5,ix,1]
                    f%seismo(1,it)=f%ez(iz,ix,1)

                    case ('ex') !ex[iz-0.5,ix,1]
                    f%seismo(1,it)=f%ex(iz,ix,1)

                    case ('es') !ex[iz-0.5,ix,1]
                    f%seismo(1,it)=f%es(iz,ix,1)
                    
                    case ('pz') !pz[iz-0.5,ix,1]
                    f%seismo(1,it)=f%pz(iz,ix,1)
                    
                    case ('px') !px[iz,ix-0.5,1]
                    f%seismo(1,it)=f%px(iz,ix,1)
                    
                end select
                
            endif
        
    end subroutine
    
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%buoz, self%buox)
        call dealloc(self%ldap2mu, self%lda, self%mu)
        ! call dealloc(self%two_ldapmu, self%inv_ladpmu_4mu)
    end subroutine

    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !<a|Au-Mf> = <a|M∂ₜu-MDMu-Mf> = <a|M∂ₜu> -<a|MDMu> -<a|Mf>
    !Kₘ<a|Au-Mf> = <a|KₘM∂ₜu> -<a|(KₘM)DMu> -<a|MD(KₘM)u> -<a|KₘMf>
    !            = <a|KₘM(DMu+Mf))> -<a|(KₘM)DMu> -<a|MD(KₘM)u> -<a|KₘMf>
    !            =                                -<a|MD(KₘM)u>
    !            =                                -<a|MD(KₘM)u>
    !            =-∫ aᵀ M D KₘM u dt
    !No approximations!
    !
    !  [ 0     0    ∂zᵇ  0  ∂ₓᶠ]    [b              ]
    !  | 0     0     0  ∂ₓᵇ ∂zᶠ|    |  b            |
    !D=|∂zᶠ    0     0   0   0 |, M=|    λ+2μ  λ    |
    !  | 0    ∂ₓᶠ    0   0   0 |    |     λ   λ+2μ  |
    !  [∂ₓᵇ   ∂zᵇ    0   0   0 ]    [              μ]
    !                            [0     ]
    !                            | 0    |
    !K_bM=diag{1,1,0,0,0}, K_λ M=|  1 1 |, K_μM=diag{0,0,2,2,1}
    !                            |  1 1 |
    !                            [     0]
    !
    !Therefore,
    !
    !       [pzᵃ]ᵀ    [  0  ]   [pzᵃ]ᵀ  [∂zᵇ(ez+ex)]
    !       |pxᵃ|     |  0  |   |pxᵃ|   |∂ₓᵇ(ez+ex)|
    !glda =-|ezᵃ| M D |ez+ex| =-|ezᵃ| M |          | =-b[ pzᵃ*∂zᵇ(ez+ex) + pxᵃ*∂ₓᵇ(ez+ex) ]  ie. -b\vec{pᵃ}\cdot\nabla(ez+ex)
    !       |exᵃ|     |ez+ex|   |exᵃ|   |          |
    !       [esᵃ]     [  0  ]   [esᵃ]   [          ]
    !
    !       [pzᵃ]ᵀ    [ 0 ]   [pzᵃ]ᵀ  [2∂zᵇez+∂ₓᶠes]
    !       |pxᵃ|     | 0 |   |pxᵃ|   |2∂ₓᵇex+∂zᶠes|
    ! gmu =-|ezᵃ| M D |2ez| =-|ezᵃ| M |            | =-b[ pzᵃ*(2∂zᵇez+∂ₓᶠes) + pxᵃ*(2∂ₓᵇex+∂zᶠes) ]
    !       |exᵃ|     |2ex|   |exᵃ|   |            |
    !       [esᵃ]     [ es]   [esᵃ]   [            ]

    
    subroutine cross_correlate_glda_gmu(rf,sf,corr,it)
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
            ! call grad3d_gkpa(rf%p,sf%pz,sf%px,sf%vy,&
            !                  corr%gkpa,             &
            !                  ifz,ilz,ifx,ilx,ify,ily)
        else
            ! call grad2d_moduli(rf%p,sf%p,sf_p_save,&
            !                    grad,            &
            !                    ifz,ilz,ifx,ilx)
            ! sf_p_save = sf%p
            
            !inexact greadient
            call grad2d_glda_gmu(rf%pz,rf%px,sf%ez,sf%ex,sf%es,&
                                ppg%buoz,ppg%buox,             &
                                corr%glda,corr%gmu,            &
                                ifz,ilz,ifx,ilx)
        endif
        
    end subroutine
    
    subroutine cross_correlate_postprocess(corr)
        type(t_correlate) :: corr

        if(allocated(correlate_gradient)) then
            !scaling gradients by model parameters
            corr%glda = corr%glda / (-m%rho)
            corr%gmu  = corr%gmu  / (-m%rho)
                    
            !preparing for projection back
            corr%glda(1,:,:) = corr%glda(2,:,:)
            corr%gmu (1,:,:) = corr%gmu (2,:,:)
        endif

    end subroutine

    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd2d_momenta(pz,px,ez,ex,es,              &
                               dez_dz,dex_dx,dex_dz,dez_dx,des_dz,des_dx,&
                               ldap2mu,lda,mu,             &
                               ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: pz,px,ez,ex,es
        real,dimension(*) :: dez_dz,dex_dx,dex_dz,dez_dx,des_dz,des_dx
        real,dimension(*) :: ldap2mu,lda,mu
        
        nz=cb%nz
        nx=cb%nx
        
        dez_dz_=0.; dex_dx_=0.; dex_dz_=0.; dez_dx_=0.; des_dz_=0.; des_dx_=0.

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dez_dz_,dex_dx_,dex_dz_,dez_dx_,des_dz_,des_dx_)
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
                

                dez_dz_= c1z*(ldap2mu(iz_ix)*ez(iz_ix)-ldap2mu(izm1_ix)*ez(izm1_ix)) +c2z*(ldap2mu(izp1_ix)*ez(izp1_ix)-ldap2mu(izm2_ix)*ez(izm2_ix))
                dex_dx_= c1x*(ldap2mu(iz_ix)*ex(iz_ix)-ldap2mu(iz_ixm1)*ex(iz_ixm1)) +c2x*(ldap2mu(iz_ixp1)*ex(iz_ixp1)-ldap2mu(iz_ixm2)*ex(iz_ixm2))

                dex_dz_= c1z*(lda(iz_ix)*ex(iz_ix)-lda(izm1_ix)*ex(izm1_ix)) +c2z*(lda(izp1_ix)*ex(izp1_ix)-lda(izm2_ix)*ex(izm2_ix))
                dez_dx_= c1x*(lda(iz_ix)*ez(iz_ix)-lda(iz_ixm1)*ez(iz_ixm1)) +c2x*(lda(iz_ixp1)*ez(iz_ixp1)-lda(iz_ixm2)*ez(iz_ixm2))

                des_dz_= c1z*(mu(izp1_ix)*es(izp1_ix)-mu(iz_ix)*es(iz_ix)) +c2z*(mu(izp2_ix)*es(izp2_ix)-mu(izm1_ix)*es(izm1_ix))
                des_dx_= c1x*(mu(iz_ixp1)*es(iz_ixp1)-mu(iz_ix)*es(iz_ix)) +c2x*(mu(iz_ixp2)*es(iz_ixp2)-mu(iz_ixm1)*es(iz_ixm1))
                                
                !cpml
                dez_dz(i)= cpml%b_z_half(iz)*dez_dz(i) + cpml%a_z_half(iz)*dez_dz_
                dex_dx(i)= cpml%b_x_half(ix)*dex_dx(i) + cpml%a_x_half(ix)*dex_dx_
                dex_dz(i)= cpml%b_z_half(iz)*dex_dz(i) + cpml%a_z_half(iz)*dex_dz_
                dez_dx(i)= cpml%b_x_half(ix)*dez_dx(i) + cpml%a_x_half(ix)*dez_dx_
                des_dz(i)= cpml%b_z(iz)     *des_dz(i) + cpml%a_z(iz)     *des_dz_
                des_dx(i)= cpml%b_x(ix)     *des_dx(i) + cpml%a_x(ix)     *des_dx_

                dez_dz_=dez_dz_*cpml%kpa_z_half(iz) + dez_dz(i)
                dex_dx_=dex_dx_*cpml%kpa_x_half(ix) + dex_dx(i)  !kappa's should have been inversed in m_computebox.f90
                dex_dz_=dex_dz_*cpml%kpa_z_half(iz) + dex_dz(i)
                dez_dx_=dez_dx_*cpml%kpa_x_half(ix) + dez_dx(i)
                des_dz_=des_dz_*cpml%kpa_z(iz)      + des_dz(i)
                des_dx_=des_dx_*cpml%kpa_x(ix)      + des_dx(i)
                
                !momenta
                pz(i)=pz(i) + dt*(dez_dz_+dex_dz_+des_dx_)
                px(i)=px(i) + dt*(dez_dx_+dex_dx_+des_dz_)

            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd2d_strains(pz,px,ez,ex,es,              &
                             dpz_dz,dpx_dx,dpz_dx,dpx_dz,&
                             buoz,buox,                  &
                             ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: pz,px,ez,ex,es
        real,dimension(*) :: dpz_dz,dpx_dx,dpz_dx,dpx_dz
        real,dimension(*) :: buoz,buox
        
        nz=cb%nz
        nx=cb%nx
        
        dpz_dz_=0.
        dpx_dx_=0.
        dpz_dx_=0.
        dpx_dz_=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dpz_dz_,dpx_dx_,dpz_dx_,dpx_dz_)
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
                

                dpz_dz_= c1z*(buoz(izp1_ix)*pz(izp1_ix)-buoz(iz_ix)*pz(iz_ix))  +c2z*(buoz(izp2_ix)*pz(izp2_ix)-buoz(izm1_ix)*pz(izm1_ix))
                dpx_dx_= c1x*(buox(iz_ixp1)*px(iz_ixp1)-buox(iz_ix)*px(iz_ix))  +c2x*(buox(iz_ixp2)*px(iz_ixp2)-buox(iz_ixm1)*px(iz_ixm1))
                
                !cpml
                dpz_dz(i)=cpml%b_z(iz)*dpz_dz(i)+cpml%a_z(iz)*dpz_dz_
                dpx_dx(i)=cpml%b_x(ix)*dpx_dx(i)+cpml%a_x(ix)*dpx_dx_

                dpz_dz_=dpz_dz_*cpml%kpa_z(iz) + dpz_dz(iz_ix)
                dpx_dx_=dpx_dx_*cpml%kpa_x(ix) + dpx_dx(iz_ix)
                
                !normal strains
                ez(i) = ez(i) + dt * dpz_dz_
                ex(i) = ex(i) + dt * dpx_dx_

                dpz_dx_= c1x*(buoz(iz_ix)*pz(iz_ix)-buoz(iz_ixm1)*pz(iz_ixm1))  +c2x*(buoz(iz_ixp1)*pz(iz_ixp1)-buoz(iz_ixm2)*pz(iz_ixm2))
                dpx_dz_= c1z*(buox(iz_ix)*px(iz_ix)-buox(izm1_ix)*px(izm1_ix))  +c2z*(buox(izp1_ix)*px(izp1_ix)-buox(izm2_ix)*px(izm2_ix))

                !cpml
                dpz_dx(i)=cpml%b_x_half(ix)*dpz_dx(i)+cpml%a_x_half(ix)*dpz_dx_
                dpx_dz(i)=cpml%b_z_half(iz)*dpx_dz(i)+cpml%a_z_half(iz)*dpx_dz_

                dpz_dx_=dpz_dx_*cpml%kpa_x_half(ix) + dpz_dx(i)
                dpx_dz_=dpx_dz_*cpml%kpa_z_half(iz) + dpx_dz(i)

                !shear stress
                es(i) = es(i) + dt * (dpz_dx_+dpx_dz_)
                
            enddo
            
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine

    subroutine grad2d_glda_gmu(rf_pz,rf_px,&
                             sf_ez,sf_ex,sf_es,&
                             buoz,buox,&
                             grad_lda,grad_mu,    &
                             ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_pz,rf_px,sf_ez,sf_ex,sf_es
        real,dimension(*) :: buoz,buox
        real,dimension(*) :: grad_lda,grad_mu
        
        nz=cb%nz
        
        sf_dpx_dx=0.
        sf_dpz_dz=0.
        sf_4_dpzdx_p_dpxdz=0.
        rf_4_es=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,izp1_ixm2,iz_ixm1,izp1_ixm1,&
        !$omp         izm2_ixp1,izm1_ixp1,iz_ixp1,izp1_ixp1,izp2_ixp1,&
        !$omp         iz_ixp2,izp1_ixp2,&
        !$omp         sf_dez_dz,sf_dex_dz,sf_des_dz,sf_dez_dx,sf_dex_dx,sf_des_dx)
        !$omp do schedule(dynamic)
!!        !$omp         sf_dpx_dx,sf_dpz_dz,&
!!        !$omp         sf_4_dpzdx_p_dpxdz,rf_4_es)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers

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
                

                sf_dez_dz = c1z*(sf_ez(iz_ix)-sf_ez(izm1_ix)) +c2z*(sf_ez(izp1_ix)-sf_ez(izm2_ix))
                sf_dex_dz = c1z*(sf_ex(iz_ix)-sf_ex(izm1_ix)) +c2z*(sf_ex(izp1_ix)-sf_ex(izm2_ix))

                sf_dez_dx = c1x*(sf_ez(iz_ix)-sf_ez(izm1_ix)) +c2x*(sf_ez(iz_ixp1)-sf_ez(iz_ixm2))
                sf_dex_dx = c1x*(sf_ex(iz_ix)-sf_ex(izm1_ix)) +c2x*(sf_ex(iz_ixp1)-sf_ex(iz_ixm2))

                !Let's try this simplified version
                !(neglecting the fact rf_pz & rf_px are living on two diff grids -- staggered grids --
                !so do sf_dez_dz & sf_dez_dx.)
                !If tests really show sth diff, then we can modify here further..
                grad_lda(j)=grad_lda(j) + rf_pz(i)*(sf_dez_dz + sf_dex_dz) &
                                        + rf_px(i)*(sf_dez_dx + sf_dex_dx)

                sf_des_dz = c1z*(sf_es(izp1_ix)-sf_es(iz_ix)) +c2z*(sf_es(izp2_ix)-sf_es(izm1_ix))
                sf_des_dx = c1x*(sf_es(iz_ixp1)-sf_es(iz_ix)) +c2x*(sf_es(iz_ixp2)-sf_es(iz_ixm1))

                grad_mu(j)=grad_mu(j) + rf_pz(i)*(2*sf_dez_dz + sf_des_dx) &
                                      + rf_px(i)*(2*sf_dex_dx + sf_des_dz)


                ! sf_4_dpzdx_p_dpxdz = &   !(dpz_dx+dpx_dx)(iz,ix)     [iz-0.5,ix-0.5]
                !                          c1x*(sf_pz(iz_ix    )-sf_pz(iz_ixm1  )) +c2x*(sf_pz(iz_ixp1  )-sf_pz(iz_ixm2  )) &
                !                        + c1z*(sf_px(iz_ix    )-sf_px(izm1_ix  )) +c2z*(sf_px(izp1_ix  )-sf_px(izm2_ix  )) &
                !                      & &
                !                      & & !(dpz_dx+dpx_dx)(iz,ix+1)   [iz-0.5,ix+0.5]
                !                        + c1x*(sf_pz(iz_ixp1  )-sf_pz(iz_ix    )) +c2x*(sf_pz(iz_ixp2  )-sf_pz(iz_ixm1  )) &
                !                        + c1z*(sf_px(iz_ixp1  )-sf_px(izm1_ixp1)) +c2z*(sf_px(izp1_ixp1)-sf_px(izm2_ixp1)) &
                !                      & &
                !                      & & !(dpz_dx+dpx_dx)(iz+1,ix)   [iz+0.5,ix-0.5]
                !                        + c1x*(sf_pz(izp1_ix  )-sf_pz(izp1_ixm1)) +c2x*(sf_pz(izp1_ixp1)-sf_pz(izp1_ixm2)) &
                !                        + c1z*(sf_px(izp1_ix  )-sf_px(iz_ix    )) +c2z*(sf_px(izp2_ix  )-sf_px(izm1_ix  )) &
                !                      & &
                !                      & & !(dpz_dx+dpx_dx)(iz+1,ix+1) [iz+0.5,ix+0.5]
                !                        + c1x*(sf_pz(izp1_ixp1)-sf_pz(izp1_ix  )) +c2x*(sf_pz(izp1_ixp2)-sf_pz(izp1_ixm1)) &
                !                        + c1z*(sf_px(izp1_ixp1)-sf_px(iz_ixp1  )) +c2z*(sf_px(izp2_ixp1)-sf_px(izm1_ixp1))

                !         ![iz-0.5,ix-0.5]   [iz+0.5,ix-0.5]   [iz-0.5,ix+0.5]     [iz+0.5,ix+0.5]
                ! rf_4_es = rf_es(iz_ix) + rf_es(izp1_ix) + rf_es(iz_ixp1) + rf_es(izp1_ixp1)

                grad_mu(j)=grad_mu(j) + rf_pz(i)*(2*sf_dez_dz +sf_des_dx) &
                                      + rf_px(i)*(2*sf_dex_dx +sf_des_dz)
                                      ! *0.0625*rf_4_es*sf_4_dpzdx_p_dpxdz   !0.0625=1/16

            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine
    
    subroutine grad2d_density(rf_pz,rf_px,         &
                              sf_ez,sf_ex,sf_es,&
                              grad,                &
                              ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_pz,rf_px
        real,dimension(*) :: sf_ez,sf_ex,sf_es
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        sf_2_dezdz_p_desdx=0.
        sf_2_dexdx_p_desdz=0.
        rf_2pz=0.
        rf_2px=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,izp1_ixm1,&
        !$omp         izm1_ixp1,iz_ixp1,izp1_ixp1,izp2_ixp1,&
        !$omp         iz_ixp2,izp1_ixp2,&
        !$omp         sf_2_dezdz_p_desdx,sf_2_dexdx_p_desdz,&
        !$omp         rf_2pz,rf_2px)
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

                sf_2_dezdz_p_desdx = &   !(des_dx+dez_dz)(iz,ix)   [iz-0.5,ix]
                                           c1z*(sf_ez(iz_ix    )-sf_ez(izm1_ix)) +c2z*(sf_ez(izp1_ix  )-sf_ez(izm2_ix  )) &
                                         + c1x*(sf_es(iz_ixp1  )-sf_es(iz_ix  )) +c2x*(sf_es(iz_ixp2  )-sf_es(iz_ixm1  )) &
                                       & &
                                       & & !(des_dx+dez_dz)(iz+1,ix) [iz+0.5,ix]
                                         + c1z*(sf_ez(izp1_ix  )-sf_ez(iz_ix  )) +c2z*(sf_ez(izp2_ix  )-sf_ez(izm1_ix  )) &
                                         + c1x*(sf_es(izp1_ixp1)-sf_es(izp1_ix)) +c2x*(sf_es(izp1_ixp2)-sf_es(izp1_ixm1))
                
                sf_2_dexdx_p_desdz = &   !(dex_dx+des_dz)(iz,ix)   [iz,ix-0.5]
                                           c1x*(sf_ex(iz_ix    )-sf_ex(iz_ixm1)) +c2x*(sf_ex(iz_ixp1  )-sf_ex(iz_ixm2  )) &
                                         + c1z*(sf_es(izp1_ix  )-sf_es(iz_ix  )) +c2z*(sf_es(izp2_ix  )-sf_es(izm1_ix  )) &
                                       & &
                                       & & !(dex_dx+des_dz)(iz,ix+1) [iz,ix+0.5]
                                         + c1x*(sf_ex(iz_ixp1  )-sf_ex(iz_ix  )) +c2x*(sf_ex(iz_ixp2  )-sf_ex(iz_ixm1  )) &
                                         + c1z*(sf_es(izp1_ixp1)-sf_es(iz_ixp1)) +c2z*(sf_es(izp2_ixp1)-sf_es(izm1_ixp1))
                                         
                rf_2pz = rf_pz(iz_ix) + rf_pz(izp1_ix)
                         ![iz-0.5,ix]      [iz+0.5,ix]
                rf_2px = rf_px(iz_ix) + rf_px(iz_ixp1)
                         ![iz,ix-0.5]      [iz,ix+0.5]
                
                grad(j)=grad(j) + 0.25*( rf_2pz*sf_2_dezdz_p_desdx + rf_2px*sf_2_dexdx_p_desdz )
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine

    subroutine imag2d_xcorr(rf_ez,rf_ex, &
                            sf_ez,sf_ex, &
                            imag,          &
                            ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_ez,rf_ex
        real,dimension(*) :: sf_ez,sf_ex
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
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
                rp = sn_p*( rf_ez(i) + rf_ex(i) )
                sp = sn_p*( sf_ez(i) + sf_ex(i) )
                
                imag(j)=imag(j) + rp*sp !+rf_pz(i)*sf_pz(i) + rf_px(i)*sf_px(i)
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine engy2d_xcorr(sf_ez,sf_ex, &
                            engy,          &
                            ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_ez,sf_ex
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
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
                sp = sn_p* ( sf_ez(i) + sf_ex(i) )
                
                engy(j)=engy(j) + sp*sp
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

end
