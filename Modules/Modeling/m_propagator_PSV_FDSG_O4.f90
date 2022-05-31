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
    real,dimension(2),parameter :: coef = [9./8.,-1./24.] !Fornberg, 1988, Generation of Finite-Difference Formulas on Arbitrary Spaced Grids.
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain ISOtropic 2D P-SV (ELastic) propagation'//s_NL// &
            '1st-order Velocity-Stress formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x4,t2) stencil'//s_NL// &
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
        real,dimension(:,:),allocatable :: two_ldapmu, inv_ladpmu_4mu

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
        procedure :: forward_scattering
        
        procedure :: adjoint
        
        procedure :: inject_velocities
        procedure :: inject_stresses
        procedure :: inject_stresses_scattering
        procedure :: update_velocities
        procedure :: update_stresses
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    !factors from normal stresses to pressure
    !real,parameter :: factor_xx=0.6666667 , factor_zz=0.3333333
    real,parameter :: factor_xx=0.5 , factor_zz=0.5

    logical :: if_hicks
    integer :: irdt
    real :: rdt

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

    subroutine init(self)
        class(t_propagator) :: self

        real,dimension(:,:),allocatable :: temp_mu

        c1x=coef(1)/m%dx; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2z=coef(2)/m%dz
        
        if_hicks=shot%if_hicks

        call alloc(self%buoz,   [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%buox,   [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%ldap2mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%lda,    [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%mu,     [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%two_ldapmu,    [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%inv_ldapmu_4mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        call alloc(temp_mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        !Staggered grid:
        ! <---*---<  
        ! |   |   |  * integer point   - lda+2mu,lda 
        ! x---v---x  < half in x dir   - buox  
        ! |   |   |  v half in z dir   - buoz
        ! <---*---<  x half in x&z dir - mu
        !     (grid index)      (real index)
        ! lda+2mu,lda(iz,ix) := lda+2mu,lda[iz,ix]
        !        buox(iz,ix) := buox[iz,ix-0.5]
        !        buoz(iz,ix) := buoz[iz-0.5,ix]
        !          mu(iz,ix) :=   mu[iz-0.5,ix-0.5]

        self%ldap2mu(:,:)=cb%rho(:,:,1)*cb%vp(:,:,1)**2
           temp_mu(:,:)=cb%rho(:,:,1)*cb%vs(:,:,1)**2

        self%lda=self%ldap2mu-2.*temp_mu
        ! if(mpiworld%is_master) then
        ! write(*,*) 'self%ldap2mu sanity:', minval(self%ldap2mu),maxval(self%ldap2mu)
        ! write(*,*) 'self%lda     sanity:', minval(self%lda),maxval(self%lda)
        ! endif

        self%two_ldapmu=2.*(self%lda+temp_mu)

        where(temp_mu>100.)
            self%inv_ldapmu_4mu=0.25/(self%lda+temp_mu)/temp_mu
        endwhere


        !interpolat mu by harmonic average
        temp_mu=1./temp_mu

        do ix=cb%ifx+1,cb%ilx
        do iz=cb%ifz+1,cb%ilz
                self%mu(iz,ix)=4./( temp_mu(iz-1,ix-1) &
                                 +temp_mu(iz-1,ix  ) &
                                 +temp_mu(iz  ,ix-1) &
                                 +temp_mu(iz  ,ix  ))
        end do
        end do

        where( ieee_is_nan(self%mu) .or. .not. ieee_is_finite(self%mu) )
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

        deallocate(temp_mu)


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
        call f%init_boundary

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

    subroutine init_abslayer(self)
        class(t_propagator) :: self

        call cpml%init

    end subroutine
    
    !========= Derivations =================
    !PDE:      A u = M ∂ₜ u - D u = f
    !Adjoint:  Aᵀa = M ∂ₜᵀa - Dᵀa = d
    !where
    !u=[vz vx szz sxx szx]ᵀ, [vz vx] are velocities, p=(szz+sxx)/2 is pressure
    !f=[fz fx fzz fxx]ᵀδ(x-xs), d is recorded data
    !  [ρ                ]    [0  0  ∂ₓ 0  ∂z]
    !  |  ρ              |    |0  0  0  ∂z ∂ₓ|
    !M=|    [λ+2μ  λ  ]⁻¹|, D=|∂ₓ 0  0  0  0 |
    !  |    | λ   λ+2μ|  |    |0  ∂z 0  0  0 |
    !  [    [      μ  ]  ]    [∂z ∂ₓ 0  0  0 ]
    !N=M⁻¹=[diag(b) κ], b=ρ⁻¹ is buoyancy,
    !a=[vzᵃ vxᵃ szzᵃ, sxxᵃ, szxᵃ]ᵀ is the adjoint field
    !
    !Continuous case:
    !<a|Au> = ∫ a(x,t) (M∂ₜ-D)u(x,t) dx³dt
    !Integration by parts, eg.:
    !∫aᵀM∂ₜu dt = aMu|t=0,T - ∫(∂ₜa)ᵀMu dt, and freely choosing a(t=T)=0 (final condition),
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
    !                      |      |     buoz    |      |
    !                      |      |      |      |      |
    !                      κ buox κ buox λ+2μ buox κ buox κ
    !  -v--p-v-p-v-→ t    -p--vx--p--vx--szz--vx--p--vx--p-→ x
    !  -1 -½ 0 ½ 1        -2 -1½ -1  -½  0   ½  1  1½  2    
    !                      |      |      |      |      | 
    !                      |      |    ½ vz     |      | 
    !                      |      |     buoz    |      | 
    !                      |      |      |      |      | 
    !                     -|------|----1-p------|------|-
    !                      |      |      κ      |      | 
    !                      |      |      |      |      | 
    !                      |      |   1½ vz     |      | 
    !                      |      |     buoz    |      | 
    !                      |      |      |      |      | 
    !                                  z ↓
    ! in space:
    ! <---*---<  
    ! |   |   |  * integer point   - sxx,szz - lda+2mu,lda 
    ! x---v---x  < half in x dir   - vx      - buox  
    ! |   |   |  v half in z dir   - vz      - buoz
    ! <---*---<  x half in x&z dir - sxz     - mu
    !
    !Forward:
    !FD eqn:
    !      [vz^n  ]   [ 0     0    ∂zᵇ][vz^n+1]      [∂zᵇ p^n+½              ] 
    !M ∂ₜᶠ |vx^n  | = | 0     0    ∂ₓᵇ||vx^n+1| +f = |∂ₓᵇ p^n+½              | +f
    !      [ p^n+1]   [∂zᶠ   ∂ₓᶠ    0 ][ p^n+½]      [∂zᶠ vz^n+1 + ∂ₓᶠ vx^n+1]
    !where
    !∂ₜᶠ*dt := v^n+1 - v^n                             ~O(t²)
    !∂zᵇ*dz := c₁(p[i  ]-p[i-½]) +c₂(p[i+ ½]-p[i-1½])  ~O(x⁴)
    !∂zᶠ*dz := c₁(v[i+½]-v[i  ]) +c₂(v[i+1½]-v[i- ½])  ~O(x⁴)
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
    !∂zᵇᵀ = c₁(v[i  ]-v[i+½]) +c₂(v[i- ½]-v[i+1½]) = -∂ₓᶠ
    !∂zᶠᵀ = c₁(p[i-½]-p[i  ]) +c₂(p[i-1½]-p[i+ ½]) = -∂ₓᵇ
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
    !Convention for half-integer index:
    !(array index)  =>    (real index)     
    !  vz(iz,ix)    =>   vz[i-½,j  ]^n  
    !  vx(iz,ix)    =>   vx[i,  j-½]^n  
    !  vy(iz,ix)    =>   vy[i,  j  ]^n  
    ! szz(iz,ix)    =>  szz[i,  j  ]^n+½
    ! szx(iz,ix)    =>  szx[i-½,j-½]^n
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
                call fld_u%boundary_transport('save',it)
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
        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')

    end subroutine

    subroutine adjoint(self,fld_a,fld_u,oif_record_adjseismo,oif_compute_imag,oif_compute_grad)
    !adjoint_a_star_Du
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
        ! if(present(o_sf)) then
                                         fld_u%dvz_dz=0.
                                         fld_u%dvx_dx=0.
            if(allocated(fld_u%dvy_dy))  fld_u%dvy_dy=0.
                                         fld_u%dp_dz=0.
                                         fld_u%dp_dx=0.
            if(allocated(fld_u%dp_dy))   fld_u%dp_dy=0.
        ! endif
                    
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
                call fld_a%check_value
                call fld_u%check_value
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then

                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_u%boundary_transport('load',it)
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
            if(if_compute_grad.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call gradient_moduli(fld_a,fld_u,it,cb%grad(:,:,:,2))
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif

            if(if_compute_imag.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call imaging(fld_a,fld_u,it,cb%imag)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif

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
            
            ! !grho: sfield%v_dt^it \dot rfield%v^it
            ! !use sfield%s^it+0.5 to compute sfield%v_dt^it, as backward step 2
            ! if(if_compute_grad.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call gradient_density(o_sf,rf,it,cb%grad(:,:,:,1))
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
                ! call fld_a%write_ext(it,'grad_moduli' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
                call fld_a%write_ext(it,'grad_a_star_Du' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
            endif
            if(self%if_compute_engy) then
                call fld_a%write_ext(it,'engy',cb%engy(:,:,:,1),size(cb%engy(:,:,:,1)))
            endif

        enddo
        
        !postprocess
        if(if_compute_imag) call imaging_postprocess
        if(if_compute_grad) call gradient_postprocess
        
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
            write(*,*) 'Elapsed time to correlate                ',tt6/mpiworld%max_threads

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
            
            wl=time_dir*f%wavelet(1,it)
            
            if(if_hicks) then
                select case (shot%src%comp)
                    case ('vz')
                    f%vz(ifz:ilz,ifx:ilx,1) = f%vz(ifz:ilz,ifx:ilx,1) + wl*self%buoz(ifz:ilz,ifx:ilx,1) *shot%src%interp_coef
                    
                    case ('vx')
                    f%vx(ifz:ilz,ifx:ilx,1) = f%vx(ifz:ilz,ifx:ilx,1) + wl*self%buox(ifz:ilz,ifx:ilx,1) *shot%src%interp_coef
                    
                end select
                
            else
                select case (shot%src%comp)
                    case ('vz') !vertical force     on vz[iz-0.5,ix]
                    f%vz(iz,ix,1) = f%vz(iz,ix,1) + wl*self%buoz(iz,ix,1)
                    
                    case ('vx') !horizontal x force on vx[iz,ix-0.5]
                    f%vx(iz,ix,1) = f%vx(iz,ix,1) + wl*self%buox(iz,ix,1)
                    
                end select
                
            endif

            return

        endif


            do i=1,shot%nrcv
                ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                
                wl=f%wavelet(i,it)
                
                if(if_hicks) then
                    select case (shot%rcv(i)%comp)
                        case ('vz') !vertical z adjsource
                        f%vz(ifz:ilz,ifx:ilx,1) = f%vz(ifz:ilz,ifx:ilx,1) + wl*self%buoz(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef !no time_dir needed!

                        case ('vx') !horizontal x adjsource
                        f%vx(ifz:ilz,ifx:ilx,1) = f%vx(ifz:ilz,ifx:ilx,1) + wl*self%buox(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef !no time_dir needed!
                        
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                        case ('vz') !vertical z adjsource
                        !vz[ix,1,iz-0.5]
                        f%vz(iz,ix,1) = f%vz(iz,ix,1) + wl*self%buoz(iz,ix,1) !no time_dir needed!

                        case ('vx') !horizontal x adjsource
                        !vx[ix-0.5,1,iz]
                        f%vx(iz,ix,1) = f%vx(iz,ix,1) + wl*self%buox(iz,ix,1) !no time_dir needed!
                        
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

        if(m%is_freesurface) ifz=max(ifz,1)

        call fd2d_velocities(f%vz,f%vx,f%szz,f%sxx,f%szx,            &
                             f%dszz_dz,f%dsxx_dx,f%dszx_dz,f%dszx_dx,&
                             self%buoz,self%buox,                    &
                             ifz,ilz,ifx,ilx,time_dir*self%dt)

        !apply free surface boundary condition if needed
        !free surface is located at [1,ix,1] level
        !dszz=0 -> -(lda+2mu)dvz_dz = lda dvx_dx
        !       -> -ldap2mu[1,ix]*(vz[1.5,ix]-vz[0.5,ix])/dz = lda[1,ix]*(vx[1,ix-0.5]-vx[1,ix+0.5])/dx
        !       -> -ldap2mu(1,ix)*(vz(2  ,ix)-vz(1  ,ix))/dz = lda(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
        !       ->  ldap2mu(1,ix)*(vz(1  ,ix)-vz(2  ,ix))/dz = lda(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
        !dsxz=0 -> -dvx_dz = dvz_dx
        !dsxz=0 -> -(vx[2,ix+0.5]-vx[0,ix+0.5])/2dz = ( (vz[0.5,ix]+vz[1.5,ix])/2 - (vz[0.5,ix+1]+vz[1.5,ix+1])/2 )/dx
        !       -> -(vx(2,ix+1  )-vx(0,ix+1  ))/2dz = ( (vz(1  ,ix)+vz(2  ,ix))/2 - (vz(1  ,ix+1)+vz(2  ,ix+1))/2 )/dx
        !       ->  (vx(0,ix+1  )-vx(2,ix+1  ))/ dz = ( (vz(1  ,ix)+vz(2  ,ix))   -  vz(1  ,ix+1)-vz(2  ,ix+1)    )/dx
        if(m%is_freesurface) then
            dz_dx = m%dz/m%dx

            do ix=ifx,ilx
                f%vz(1,ix,1)= f%vz(2,ix,1) + self%lda(1,ix,1)*(f%vx(1,ix,1)-f%vx(1,ix+1,1))*dz_dx/self%ldap2mu(1,ix,1)
                !f%vx(0,ix)= f%vx(2,ix) + (f%vz(1,ix)+f%vz(2,ix)-f%vz(1,ix+1)-f%vz(2,ix+1))*dz_dx/lm%ldap2mu(1,ix) !toy2del says this condition is not needed
            enddo
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
            
            wl=time_dir*f%wavelet(1,it)
            
            if(if_hicks) then
                if(shot%src%comp=='p') then
                    f%szz(ifz:ilz,ifx:ilx,1) = f%szz(ifz:ilz,ifx:ilx,1) + wl*self%ldap2mu(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)
                    f%sxx(ifz:ilz,ifx:ilx,1) = f%sxx(ifz:ilz,ifx:ilx,1) + wl*self%ldap2mu(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)
                endif
                
            else
                if(shot%src%comp=='p') then
                    !explosion on s[iz,ix,1]
                    f%szz(iz,ix,1) = f%szz(iz,ix,1) + wl*self%ldap2mu(iz,ix,1)
                    f%sxx(iz,ix,1) = f%sxx(iz,ix,1) + wl*self%ldap2mu(iz,ix,1)
                endif
                
            endif

            return

        endif

            do i=1,shot%nrcv

                if(shot%rcv(i)%comp=='p') then

                    ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                    ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                    
                    !adjsource for pressure
                    wl=f%wavelet(i,it)
                    
                    if(if_hicks) then 

                        f%szz(ifz:ilz,ifx:ilx,1) = f%szz(ifz:ilz,ifx:ilx,1) +wl*self%ldap2mu(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                        f%sxx(ifz:ilz,ifx:ilx,1) = f%sxx(ifz:ilz,ifx:ilx,1) +wl*self%ldap2mu(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!

                    else           
                        !s[iz,ix,1]
                        f%szz(iz,ix,1) = f%szz(iz,ix,1) +wl*self%ldap2mu(iz,ix,1) !no time_dir needed!
                        f%sxx(iz,ix,1) = f%sxx(iz,ix,1) +wl*self%ldap2mu(iz,ix,1) !no time_dir needed!

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
        1=f%bloom(5,it)+1
        1=f%bloom(6,it)-2
        
        if(m%is_freesurface) ifz=max(ifz,1)

        call fd2d_stresses(f%vz,f%vx,f%szz,f%sxx,f%szx,    &
                           f%dvz_dz,f%dvx_dx,              &
                           self%ldap2mu,self%lda,self%mu,  &
                           ifz,ilz,ifx,ilx,time_dir*self%dt)
        
        !apply free surface boundary condition if needed
        !free surface is located at [1,ix] level
        !so explicit boundary condition: szz(1,ix)=0
        !and antisymmetric mirroring: sxz[0.5,ix-0.5]=-sxz[1.5,ix-0.5] -> sxz(1,ix)=-sxz(2,ix)
        if(m%is_freesurface) then
            f%szz(1,:,1)=0.
            f%sxz(1,:,1)=-f%sxz(2,:,1)
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
                        
                        case ('p')
                        f%seismo(i,it)=sum( (factor_zz*f%szz(ifz:ilz,ifx:ilx,1) + factor_xx*f%sxx(ifz:ilz,ifx:ilx,1)) *shot%rcv(i)%interp_coef(:,:,1) )
                        case ('vz')
                        f%seismo(i,it)=sum(f%vz(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                        case ('vx')
                        f%seismo(i,it)=sum(f%vx(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                        case ('p') !p[iz,ix]
                        f%seismo(i,it)=factor_zz*f%szz(iz,ix,1) + factor_xx*f%sxx(iz,ix,1)
                        case ('vz') !vz[iz-0.5,ix]
                        f%seismo(i,it)=f%vz(iz,ix,1)
                        case ('vx') !vx[iz,ix-0.5]
                        f%seismo(i,it)=f%vx(iz,ix,1)
                    end select
                    
                endif

            enddo

            return

        endif

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            if(if_hicks) then
                select case (shot%src%comp)
                    case ('p')
                    f%seismo(1,it)=sum(f%p(ifz:ilz,ifx:ilx,1) *shot%src%interp_coef)
                    
                    case ('vz')
                    f%seismo(1,it)=sum(f%vz(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef)
                    
                    case ('vx')
                    f%seismo(1,it)=sum(f%vx(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef)
                    
                end select
                
            else
                select case (shot%src%comp)
                    case ('p') !p[iz,ix,1]
                    f%seismo(1,it)=f%p(iz,ix,1)
                    
                    case ('vz') !vz[iz-0.5,ix,1]
                    f%seismo(1,it)=f%vz(iz,ix,1)
                    
                    case ('vx') !vx[iz,ix-0.5,1]
                    f%seismo(1,it)=f%vx(iz,ix,1)
                    
                end select
                
            endif
        
    end subroutine
    
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%buox, self%buoy, self%buoz, self%kpa, self%inv_kpa)
    end subroutine

    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !∇ₘ<a|Au> = ∇ₘ<a|M∂ₜu-Du> = ∫ aᵀ ∇ₘM ∂ₜu dt
    !Since it's cumbersome to get ∂ₜu by time marching,
    !replace ∂ₜu by M⁻¹Du and neglect f
    !ie. M∂ₜu=Du+f -> ∂ₜu=M⁻¹Du+M⁻¹f ≐ M⁻¹Du
    !This simplification introduces singularities in the gradient only at source positions, which are probably removed by gradient masking.
    !
    !Therefore, ∫ aᵀ ∇ₘM ∂ₜu dt ≐ ∫ aᵀ ∇ₘln(M) Du dt =: a★Du
    !where
    !                            [ 0   0   ∂ₓᵇ] [vx]
    !a★Du = [vxᵃ vzᵃ pᵃ] ∇ₘln(M) | 0   0   ∂zᵇ| |vz|
    !                            [∂ₓᶠ  ∂zᶠ  0 ] [p ]
    !In particular, we compute
    !  grho = vᵃ ∂ₜv = vᵃ b∇p
    !  gkpa = pᵃ (-κ⁻²) ∂ₜp = = pᵃ (-κ⁻¹) ∇·v
    !
    !For imaging:
    !I = ∫ a u dt =: a★u

    subroutine gradient_density(rf,sf,it,grad)
        type(t_field), intent(in) :: rf, sf
        real,dimension(cb%mz,cb%mx,cb%my) :: grad
        
        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        
        call grad2d_density(rf%vz,rf%vx,sf%szz,sf%sxx,sf%szx,&
                            grad,                            &
                            ifz,ilz,ifx,ilx)
        
    end subroutine

    subroutine gradient_moduli(rf,sf,it,grad)
        type(t_field), intent(in) :: rf, sf
        real,dimension(cb%mz,cb%mx,cb%my) :: grad

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
                    
        !inexact greadient
        call grad2d_moduli(rf%szz,rf%sxx,rf%szx,sf%vz,sf%vx,&
                           grad,                            &
                           ifz,ilz,ifx,ilx)
        
    end subroutine
    
    subroutine gradient_postprocess

        !scale the kernel tobe a gradient in the discretized world
        cb%grad = cb%grad*m%cell_volume*rdt
        
        !grho
        cb%grad(:,:,:,1) = cb%grad(:,:,:,1) / cb%rho(1:cb%mz,1:cb%mx,1:cb%my)
        !gkpa
        cb%grad(:,:,:,2) = cb%grad(:,:,:,2) * (-ppg%inv_kpa(1:cb%mz,1:cb%mx,1:cb%my))
                
        !preparing for cb%project_back
        cb%grad(1,:,:,:) = cb%grad(2,:,:,:)

    end subroutine

    ! subroutine imaging(rf,sf,it,imag)
    !     type(t_field), intent(in) :: rf, sf
    !     real,dimension(cb%mz,cb%mx,cb%my) :: imag
        
    !     !nonzero only when sf touches rf
    !     ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
    !     ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
    !     ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
    !     ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
    !     1=max(sf%bloom(5,it),rf%bloom(5,it),1)
    !     1=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
    !     ! if(m%is_cubic) then
    !     !     call imag3d_xcorr(rf%p,sf%p,&
    !     !                       imag,                  &
    !     !                       ifz,ilz,ifx,ilx,1,1)
    !     ! else
    !         call imag2d_xcorr(rf%p,sf%p,&
    !                           imag,            &
    !                           ifz,ilz,ifx,ilx)
    !     ! endif

    !     ! call imag2d_xcorr(rf%p,rf%vz,rf%vx,&
    !     !                   sf%p,sf%vz,sf%vx,&
    !     !                   imag,            &
    !     !                   ifz,ilz,ifx,ilx)

    ! end subroutine

    ! subroutine imaging_postprocess

    !     !scale by rdt tobe an image in the discretized world
    !     cb%imag = cb%imag*rdt

    !     !for cb%project_back
    !     cb%imag(1,:,:,:) = cb%imag(2,:,:,:)

    ! end subroutine

    ! subroutine energy(sf,it,engy)
    !     type(t_field),intent(in) :: sf
    !     real,dimension(cb%mz,cb%mx,cb%my) :: engy

    !     !nonzero only when sf touches rf
    !     ifz=max(sf%bloom(1,it),2)
    !     ilz=min(sf%bloom(2,it),cb%mz)
    !     ifx=max(sf%bloom(3,it),1)
    !     ilx=min(sf%bloom(4,it),cb%mx)
    !     1=max(sf%bloom(5,it),1)
    !     1=min(sf%bloom(6,it),cb%my)
        
    !     if(m%is_cubic) then
    !         call engy3d_xcorr(sf%p,&
    !                           engy,                  &
    !                           ifz,ilz,ifx,ilx,1,1)
    !     else
    !         call engy2d_xcorr(sf%p,&
    !                           engy,            &
    !                           ifz,ilz,ifx,ilx)
    !     endif

    ! end subroutine

    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd2d_velocities(vz,vx,szz,sxx,sxz,&
                               dszz_dz,dsxx_dx,dszx_dz,dszx_dx,&
                               buoz,buox,        &
                               ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vz,vx,p
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
                

                dszz_dz= c1z*(szz(iz_ix)-szz(izm1_ix)) +c2z*(szz(izp1_ix)-szz(izm2_ix))
                dsxx_dx= c1x*(sxx(iz_ix)-sxx(iz_ixm1)) +c2x*(sxx(iz_ixp1)-sxx(iz_ixm2))

                dszx_dz= c1z*(szx(izp1_ix)-szx(iz_ix)) +c2z*(szx(izp2_ix)-szx(izm1_ix))
                dszx_dx= c1x*(szx(iz_ixp1)-szx(iz_ix)) +c2x*(szx(iz_ixp2)-szx(iz_ixm1))
                                
                !cpml
                cpml%dsxx_dx(i)= cpml%b_x_half(ix)*cpml%dsxx_dx(i) + cpml%a_x_half(ix)*dsxx_dx
                cpml%dsxz_dz(i)= cpml%b_z(iz)     *cpml%dsxz_dz(i) + cpml%a_z(iz)     *dsxz_dz 
                cpml%dsxz_dx(i)= cpml%b_x(ix)     *cpml%dsxz_dx(i) + cpml%a_x(ix)     *dsxz_dx
                cpml%dszz_dz(i)= cpml%b_z_half(iz)*cpml%dszz_dz(i) + cpml%a_z_half(iz)*dszz_dz

                dszz_dz=dszz_dz*cb%kappa_z_half(iz) + cpml%dszz_dz(i)
                dsxx_dx=dsxx_dx*cb%kappa_x_half(ix) + cpml%dsxx_dx(i)  !kappa's should have been inversed in m_computebox.f90
                dszx_dz=dszx_dz*cb%kappa_z(iz)      + cpml%dszx_dz(i)
                dszx_dx=dszx_dx*cb%kappa_x(ix)      + cpml%dszx_dx(i)
                
                !velocity
                vz(i)=vz(i) + dt*buoz(i)*(dszz_dz+dszx_dx)
                vx(i)=vx(i) + dt*buox(i)*(dszx_dz+dsxx_dx)

            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd3d_stresses(vz,vx,vy,p,               &
                             dvz_dz,dvx_dx,dvy_dy,     &
                             kpa,                      &
                             ifz,ilz,ifx,ilx,1,1,dt)
        real,dimension(*) :: vz,vx,vy,p
        real,dimension(*) :: dvz_dz,dvx_dx,dvy_dy
        real,dimension(*) :: kpa
        
        nz=cb%nz
        nx=cb%nx
        ny=cb%ny
        
        dvz_dz_=0.;dvx_dx_=0.;dvy_dy_=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,1,i,&
        !$omp         izm1_ix_1,iz_ix_1,izp1_ix_1,izp2_ix_1,&
        !$omp         iz_ixm1_1,iz_ixp1_1,iz_ixp2_1,&
        !$omp         iz_ix_1m1,iz_ix_1p1,iz_ix_1p2,&
        !$omp         dvz_dz_,dvx_dx_,dvy_dy_)
        !$omp do schedule(dynamic) collapse(2)
        do 1=1,1
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+(1-cb%1)*nz*nx+1
                
                izm1_ix_1=i-1  !iz-1,ix,1
                iz_ix_1  =i    !iz,ix,1
                izp1_ix_1=i+1  !iz+1,ix,1
                izp2_ix_1=i+2  !iz+2,ix,1
                
                iz_ixm1_1=i       -nz  !iz,ix-1,1
                iz_ixp1_1=i       +nz  !iz,ix+1,1
                iz_ixp2_1=i     +2*nz  !iz,ix+2,1
                
                iz_ix_1m1=i    -nz*nx  !iz,ix,1-1
                iz_ix_1p1=i    +nz*nx  !iz,ix,1+1
                iz_ix_1p2=i  +2*nz*nx  !iz,ix,1+2
                
                dvx_dx_= c1x*(vx(iz_ixp1_1)-vx(iz_ix_1))  +c2x*(vx(iz_ixp2_1)-vx(iz_ixm1_1))
                dvy_dy_= c1y*(vy(iz_ix_1p1)-vy(iz_ix_1))  +c2y*(vy(iz_ix_1p2)-vy(iz_ix_1m1))
                dvz_dz_= c1z*(vz(izp1_ix_1)-vz(iz_ix_1))  +c2z*(vz(izp2_ix_1)-vz(izm1_ix_1))
                
                !cpml
                dvz_dz(iz_ix_1)=cpml%b_z(iz)*dvz_dz(iz_ix_1)+cpml%a_z(iz)*dvz_dz_
                dvx_dx(iz_ix_1)=cpml%b_x(ix)*dvx_dx(iz_ix_1)+cpml%a_x(ix)*dvx_dx_
                dvy_dy(iz_ix_1)=cpml%b_y(1)*dvy_dy(iz_ix_1)+cpml%a_y(1)*dvy_dy_

                dvz_dz_=dvz_dz_*cpml%kpa_z(iz) + dvz_dz(iz_ix_1)
                dvx_dx_=dvx_dx_*cpml%kpa_x(ix) + dvx_dx(iz_ix_1)
                dvy_dy_=dvy_dy_*cpml%kpa_y(1) + dvy_dy(iz_ix_1)
                
                !pressure
                p(iz_ix_1) = p(iz_ix_1) + dt * kpa(iz_ix_1)*(dvz_dz_+dvx_dx_+dvy_dy_)
                
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
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%1:cb%1) :: vz

        !free surface is located at [1,ix,1] level
        !so symmetric mirroring: vz[0.5]=vz[1.5], ie. vz(1,ix,1)=vz(2,ix,1) -> dp(1,ix,1)=0.
            vz(1,:,:)=vz(2,:,:)
            ! !$omp parallel default (shared)&
            ! !$omp private(ix,1,i)
            ! !$omp do schedule(dynamic)
            ! do 1=1,1
            ! do ix=ifx,ilx
            !     i=(1-cb%ifz) + (ix-cb%ifx)*nz + (1-cb%1)*nz*nx +1 !iz=1,ix,1
                
            !     f%vz(i)=f%vz(i+1)
            ! enddo
            ! enddo
            ! !$omp enddo
            ! !$omp end parallel

    end subroutine

    subroutine fd_freesurface_stresses(p)
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%1:cb%1) :: p

        !free surface is located at [1,ix,1] level
        !so explicit boundary condition: p(1,ix,1)=0
        !and antisymmetric mirroring: p(0,ix,1)=-p(2,ix,1) -> vz(2,ix,1)=vz(1,ix,1)
            p(1,:,:)=0.
            p(0,:,:)=-p(2,:,:)
            ! !$omp parallel default (shared)&
            ! !$omp private(ix,1,i)
            ! !$omp do schedule(dynamic)
            ! do 1=1,1
            ! do ix=ifx,ilx
            !     i=(1-cb%ifz) + (ix-cb%ifx)*nz + (1-cb%1)*nz*nx +1 !iz=1,ix,1 
                
            !     f%p(i)=0.
                
            !     f%p(i-1)=-f%p(i+1)
            ! enddo
            ! enddo
            ! !$omp enddo
            ! !$omp end parallel

    end subroutine

    
    subroutine grad3d_moduli(rf_p,sf_vz,sf_vx,sf_vy,&
                             grad,                  &
                             ifz,ilz,ifx,ilx,1,1)
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
        !$omp private(iz,ix,1,i,j,&
        !$omp         izm1_ix_1,iz_ix_1,izp1_ix_1,izp2_ix_1,&
        !$omp         iz_ixm1_1,iz_ixp1_1,iz_ixp2_1,&
        !$omp         iz_ix_1m1,iz_ix_1p1,iz_ix_1p2,&
        !$omp         dvz_dz,dvx_dx,dvy_dy,&
        !$omp         rp,dsp)
        !$omp do schedule(dynamic) collapse(2)
        do 1=1,1
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(1-cb%1)*cb%nz*cb%nx+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+(1-1)     *cb%mz*cb%mx+1 !grad has no boundary layers
                
                izm1_ix_1=i-1  !iz-1,ix,1
                iz_ix_1  =i    !iz,ix,1
                izp1_ix_1=i+1  !iz+1,ix,1
                izp2_ix_1=i+2  !iz+2,ix,1
                
                iz_ixm1_1=i    -nz  !iz,ix-1,1
                iz_ixp1_1=i    +nz  !iz,ix+1,1
                iz_ixp2_1=i  +2*nz  !iz,ix+2,1
                
                iz_ix_1m1=i    -nz*nx  !iz,ix,1-1
                iz_ix_1p1=i    +nz*nx  !iz,ix,1+1
                iz_ix_1p2=i  +2*nz*nx  !iz,ix,1+2
                
                dvz_dz = c1z*(sf_vz(izp1_ix_1)-sf_vz(iz_ix_1)) +c2z*(sf_vz(izp2_ix_1)-sf_vz(izm1_ix_1))
                dvx_dx = c1x*(sf_vx(iz_ixp1_1)-sf_vx(iz_ix_1)) +c2x*(sf_vx(iz_ixp2_1)-sf_vx(iz_ixm1_1))
                dvy_dy = c1y*(sf_vy(iz_ix_1p1)-sf_vy(iz_ix_1)) +c2y*(sf_vy(iz_ix_1p2)-sf_vy(iz_ix_1m1))
                
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

    subroutine grad2d_moduli(rf_p,sf_vz,sf_vx,&
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
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
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
    
    subroutine grad3d_density(rf_vz,rf_vx,rf_vy,sf_p,&
                              grad,                  &
                              ifz,ilz,ifx,ilx,1,1)
        real,dimension(*) :: rf_vz,rf_vx,rf_vy,sf_p
        real,dimension(*) :: grad
        
        nz=cb%nz
        nx=cb%nx
        
        dsvz=0.; dsvx=0.; dsvy=0.
         rvz=0.;  rvx=0.;  rvy=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,1,i,j,&
        !$omp         izm2_ix_1,izm1_ix_1,iz_ix_1,izp1_ix_1,izp2_ix_1,&
        !$omp         iz_ixm2_1,iz_ixm1_1,iz_ixp1_1,iz_ixp2_1,&
        !$omp         iz_ix_1m2,iz_ix_1m1,iz_ix_1p1,iz_ix_1p2,&
        !$omp          rvz, rvx, rvy,&
        !$omp         dsvz,dsvx,dsvy)
        !$omp do schedule(dynamic) collapse(2)
        do 1=1,1
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(1-cb%1)*cb%nz*cb%nx+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+(1-1)     *cb%mz*cb%mx+1 !grad has no boundary layers

                izm2_ix_1=i-2  !iz-2,ix,1
                izm1_ix_1=i-1  !iz-1,ix,1
                iz_ix_1  =i    !iz,ix,1
                izp1_ix_1=i+1  !iz+1,ix,1
                izp2_ix_1=i+2  !iz+2,ix,1
                
                iz_ixm2_1=i  -2*nz  !iz,ix-2,1
                iz_ixm1_1=i    -nz  !iz,ix-1,1
                iz_ixp1_1=i    +nz  !iz,ix+1,1
                iz_ixp2_1=i  +2*nz  !iz,ix+2,1
                
                iz_ix_1m2=i  -2*nz*nx  !iz,ix,1-2
                iz_ix_1m1=i    -nz*nx  !iz,ix,1-1
                iz_ix_1p1=i    +nz*nx  !iz,ix,1+1
                iz_ix_1p2=i  +2*nz*nx  !iz,ix,1+2               
                
                rvz = rf_vz(izp1_ix_1) +rf_vz(iz_ix_1)
                rvx = rf_vx(iz_ixp1_1) +rf_vx(iz_ix_1)
                rvy = rf_vy(iz_ix_1p1) +rf_vy(iz_ix_1)

                dsvz = (c1z*(sf_p(iz_ix_1  )-sf_p(izm1_ix_1)) +c2z*(sf_p(izp1_ix_1)-sf_p(izm2_ix_1))) &
                      +(c1z*(sf_p(izp1_ix_1)-sf_p(iz_ix_1  )) +c2z*(sf_p(izp2_ix_1)-sf_p(izm1_ix_1)))
                dsvx = (c1x*(sf_p(iz_ix_1  )-sf_p(iz_ixm1_1)) +c2x*(sf_p(iz_ixp1_1)-sf_p(iz_ixm2_1))) &
                      +(c1x*(sf_p(iz_ixp1_1)-sf_p(iz_ix_1  )) +c2x*(sf_p(iz_ixp2_1)-sf_p(iz_ixm1_1)))
                dsvy = (c1y*(sf_p(iz_ix_1  )-sf_p(iz_ix_1m1)) +c2y*(sf_p(iz_ix_1p1)-sf_p(iz_ix_1m2))) &
                      +(c1y*(sf_p(iz_ix_1p1)-sf_p(iz_ix_1  )) +c2y*(sf_p(iz_ix_1p2)-sf_p(iz_ix_1m1)))
                
                !complete equation with unnecessary terms e.g. sf_p(iz_ix) for better understanding
                !with flag -O, the compiler should automatically detect such possibilities of simplification
                
                grad(j)=grad(j) + 0.25*( rvz*dsvz + rvx*dsvx + rvy*dsvy )
                
            enddo
            
        enddo
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine grad2d_density(rf_vz,rf_vx,sf_p,&
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
        !$omp          rvz, rvx,&
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


    subroutine imag3d_xcorr(rf_p,sf_p,             &
                            imag,                  &
                            ifz,ilz,ifx,ilx,1,1)
        real,dimension(*) :: rf_p,sf_p
        real,dimension(*) :: imag
        
        nz=cb%nz
        nx=cb%nx
        
        rp=0.
        sp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,1,i,j,&
        !$omp         rp,sp)
        !$omp do schedule(dynamic) collapse(2)
        do 1=1,1
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(1-cb%1)*cb%nz*cb%nx+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+(1-1)     *cb%mz*cb%mx+1 !grad has no boundary layers
                
                rp = rf_p(i)
                sp = sf_p(i)
                
                imag(j)=imag(j) + rp*sp
                
            end do
            
        end do
        end do
        !$omp end do
        !$omp end parallel
        
    end subroutine

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

    subroutine imag2d_xcorr(rf_p,sf_p,     &
                            imag,          &
                            ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_p,sf_p
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
                
                rp = rf_p(i)
                sp = sf_p(i)
                
                imag(j)=imag(j) + rp*sp
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine engy3d_xcorr(sf_p,          &
                            engy,          &
                            ifz,ilz,ifx,ilx,1,1)
        real,dimension(*) :: sf_p
        real,dimension(*) :: engy
        
        nz=cb%nz
        nx=cb%nx
        
        sp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,1,i,j,&
        !$omp         sp)
        !$omp do schedule(dynamic) collapse(2)
        do 1=1,1
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(1-cb%1)*cb%nz*cb%nx+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+(1-1)     *cb%mz*cb%mx+1 !grad has no boundary layers
                
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
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
                sp = sf_p(i)
                
                engy(j)=engy(j) + sp*sp
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

end








module m_field
use m_mpienv
use m_arrayop
use m_sysio
use m_model, only:m
use m_shot, only:shot
use m_computebox, only: cb

use, intrinsic :: ieee_arithmetic

    private fdcoeff_o4,fdcoeff_o8,c1x,c1z,c2x,c2z,c3x,c3z,c4x,c4z
    !private kappa,npower,rcoef
    public
    
    !FD coeff
    real,dimension(2),parameter :: fdcoeff_o4 = [1.125,      -1./24.]
    real,dimension(4),parameter :: fdcoeff_o8 = [1225./1024, -245./3072., 49./5120., -5./7168]
    
    real :: c1x, c1z
    real :: c2x, c2z
    real :: c3x, c3z
    real :: c4x, c4z
    
    !!CPML
    !real,parameter :: kappa=1.
    !real,parameter :: npower=2.
    !real,parameter :: rcoef=0.001
    
    !local models in computebox
    type t_localmodel
        sequence
        !shorthand for greek letters
        !alp bta gma  del(dta) eps zta 
        !eta tht iota kpa lda mu
        !nu xi omi pi rho sgm
        !tau ups phi chi psi oga
        !etc:
        !buo=buoyancy
        real,dimension(:,:),allocatable :: buox, buoz, ldap2mu, lda, mu
        real,dimension(:,:),allocatable :: two_ldapmu,inv_ldapmu_4mu
    end type
    
    type(t_localmodel) :: lm
    
    !fields
    type t_field
        sequence
        real,dimension(:,:),allocatable :: vx,vy,vz,sxx,szz,sxz  !vy is just for fun
        real,dimension(:,:),allocatable :: cpml_dvx_dx,cpml_dvz_dz,cpml_dvx_dz,cpml_dvz_dx !for cpml
        real,dimension(:,:),allocatable :: cpml_dsxx_dx,cpml_dszz_dz,cpml_dsxz_dx,cpml_dsxz_dz
    end type

    !sampling factor
    !real,parameter :: factor_xx=0.6666667 , factor_zz=0.3333333
    real,parameter :: factor_xx=0.5 , factor_zz=0.5
    
    !hicks interpolation
    logical :: if_hicks
    
    !info
    character(*),parameter :: waveeq_info='time-domain isotropic P-SV (2D elastic)'
    character(*),parameter :: gradient_info='lda-mu-rho'
    integer,parameter :: ncorr=3
    
    contains
    
    !========= use before propagation =================
    subroutine field_print_info
        !modeling method
        call hud('WaveEq : Time-domain isotropic P-SV (2D elastic) system')
        call hud('1st-order Velocity-Stress formulation')
        call hud('Staggered-grid Finite-difference method')
        call hud('Cartesian O(x4,t2) stencil')
        
        !stencil constant
        if(mpiworld%is_master) then
            write(*,*) 'Coeff:',fdcoeff_o4
        endif
        c1x=fdcoeff_o4(1)/m%dx; c1z=fdcoeff_o4(1)/m%dz
        c2x=fdcoeff_o4(2)/m%dx; c2z=fdcoeff_o4(2)/m%dz
        
    end subroutine
    
    subroutine check_model  !not required
        
    end
    
    subroutine check_discretization
        !grid dispersion condition
        if (5.*m%cell_diagonal > cb%velmin/shot%src%fpeak/2.) then  !FDTDo4 rule
            write(*,*) 'WARNING: Shot# '//shot%cindex//' can have grid dispersion!'
            write(*,*) 'Shot# '//shot%cindex//' 5*dx, velmin, fpeak:',5.*m%cell_diagonal, cb%velmin,shot%src%fpeak
        endif
        
        !CFL condition
        cfl = cb%velmax*shot%src%dt*m%cell_inv_diagonal*sum(abs(fdcoeff_o4))! ~0.606
        if(mpiworld%is_master) write(*,*) 'CFL value:',CFL
        
        if(cfl>1.) then
            write(*,*) 'ERROR: CFL > 1 on shot# '//shot%cindex//'!'
            write(*,*) 'Shot# '//shot%cindex//' velmax, dt, 1/dx:',cb%velmax,shot%src%dt,m%cell_inv_diagonal
            stop
        endif
        
    end subroutine
    
    !========= use inside propagation =================
    
    subroutine init_field_localmodel
        real,dimension(:,:),allocatable :: temp_mu
        
        ! (Lamé par)  (Rock phy)  (Seismic wave)     (Voigt)
        ! lda+2mu   = K+4/3G    = rho*Vp^2         = c11=c33
        ! lda       = K-2/3G    = rho*(Vp^2-2Vs^2) = c13
        ! mu        = G         = rho*Vs^2         = c44
        !
        ! Poisson ratio nv = (3K-2G)/(6k+2G) = 0.5-3/(6sqrt(Vp/Vs)+2)

        call alloc(lm%buox,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(lm%buoz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(lm%ldap2mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(lm%lda, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(lm%mu,  [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(lm%two_ldapmu,    [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(lm%inv_ldapmu_4mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        call alloc(temp_mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        
        !Staggered grid:
        ! <---*---<  
        ! |   |   |  * integer point   - lda+2mu,lda 
        ! x---v---x  < half in x dir   - buox  
        ! |   |   |  v half in z dir   - buoz
        ! <---*---<  x half in x&z dir - mu
        !     (grid index)      (real index)
        ! lda+2mu,lda(iz,ix) := lda+2mu,lda[iz,ix]
        !        buox(iz,ix) := buox[iz,ix-0.5]
        !        buoz(iz,ix) := buoz[iz-0.5,ix]
        !          mu(iz,ix) :=   mu[iz-0.5,ix-0.5]

        lm%ldap2mu(:,:)=cb%rho(:,:,1)*cb%vp(:,:,1)**2
           temp_mu(:,:)=cb%rho(:,:,1)*cb%vs(:,:,1)**2

        lm%lda=lm%ldap2mu-2.*temp_mu
        ! if(mpiworld%is_master) then
        ! write(*,*) 'lm%ldap2mu sanity:', minval(lm%ldap2mu),maxval(lm%ldap2mu)
        ! write(*,*) 'lm%lda     sanity:', minval(lm%lda),maxval(lm%lda)
        ! endif

        lm%two_ldapmu=2.*(lm%lda+temp_mu)

        where(temp_mu>100.)
            lm%inv_ldapmu_4mu=0.25/(lm%lda+temp_mu)/temp_mu
        endwhere


        !interpolat mu by harmonic average
        temp_mu=1./temp_mu

        do ix=cb%ifx+1,cb%ilx
        do iz=cb%ifz+1,cb%ilz
                lm%mu(iz,ix)=4./( temp_mu(iz-1,ix-1) &
                                 +temp_mu(iz-1,ix  ) &
                                 +temp_mu(iz  ,ix-1) &
                                 +temp_mu(iz  ,ix  ))
        end do
        end do

        where( ieee_is_nan(lm%mu) .or. .not. ieee_is_finite(lm%mu) )
            lm%mu=0.
        endwhere

        ! open(8,file='lm%mu',access='stream')
        ! write(8) lm%mu
        ! close(8)

        ! !interpolat mu by arithmic average
        ! do ix=cb%ifx+1,cb%ilx
        ! do iz=cb%ifz+1,cb%ilz
        !         lm%mu(iz,ix)=( temp_mu(iz-1,ix-1) &
        !                       +temp_mu(iz-1,ix  ) &
        !                       +temp_mu(iz  ,ix-1) &
        !                       +temp_mu(iz  ,ix  ))
        ! end do
        ! end do
        ! lm%mu=lm%mu*0.25

        lm%mu(cb%ifz,:)=lm%mu(cb%ifz+1,:)
        lm%mu(:,cb%ifx)=lm%mu(:,cb%ifx+1)

        !check mu values
        if(mpiworld%is_master) then
            write(*,*) 'lm%mu sanity:', minval(lm%mu),maxval(lm%mu), any(ieee_is_nan(lm%mu)), any(.not. ieee_is_finite(lm%mu))
        endif
        

        do iz=cb%ifz+1,cb%ilz
            lm%buoz(iz,:)=0.5/cb%rho(iz,:,1)+0.5/cb%rho(iz-1,:,1)
        enddo
        
        do ix=cb%ifx+1,cb%ilx
            lm%buox(:,ix)=0.5/cb%rho(:,ix,1)+0.5/cb%rho(:,ix-1,1)
        enddo

        lm%buoz(cb%ifz,:)=1./cb%rho(cb%ifz,:,1)
        lm%buox(:,cb%ifx)=1./cb%rho(:,cb%ifx,1)

        deallocate(temp_mu)

    end subroutine
    
    subroutine init_field(f)
        type(t_field) :: f
        
        call alloc(f%vx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%vz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%sxx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%szz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%sxz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        
        call alloc(f%cpml_dvx_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dvz_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dvx_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dvz_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dsxx_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dszz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dsxz_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dsxz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        
        !hicks point interpolation
        if_hicks=get_setup_logical('IF_HICKS',default=.true.)

    end subroutine
    
    subroutine check_field(f,name)
        type(t_field) :: f
        character(*) :: name
        
        if(mpiworld%is_master) write(*,*) name//' sample values:',minval(f%vz),maxval(f%vz)
        
        if(any(.not. ieee_is_finite(f%vz))) then
            write(*,*) 'ERROR: '//name//' values become Infinity on Shot# '//shot%cindex//' !!'
            stop
        endif
        if(any(ieee_is_nan(f%vz))) then
           write(*,*) 'ERROR: '//name//' values become NaN on Shot# '//shot%cindex//' !!'
           stop
        endif
        
    end subroutine
    
    subroutine field_cpml_reinitialize(f)
        type(t_field) :: f
        f%cpml_dvx_dx=0.
        f%cpml_dvz_dz=0.
        f%cpml_dvx_dz=0.
        f%cpml_dvz_dx=0.
        f%cpml_dsxx_dx=0.
        f%cpml_dszz_dz=0.
        f%cpml_dsxz_dx=0.
        f%cpml_dsxz_dz=0.
    end subroutine
    
    subroutine write_field(iunit,f)
        type(t_field) :: f
        write(iunit) f%vz
    end subroutine
    
    !========= forward propagation =================
    !WE: du_dt = MDu
    !u=[vx vz sxx szz sxz]^T, p=(2sxx+szz)/3
    !  [b                     ]    [      dx 0  dz]
    !  |  b                   |    |      0  dz dx|
    !M=|   lda+2mu   lda      |, D=|dx 0          |
    !  |      lda   lda+2mu   |    |0  dz         |
    !  [                    mu]    [dz dx         ]
    !  (b=1/rho)
    !
    !Discretization (staggered grid in space and time):
    !  (grid index)     (real index)
    ! sxx,szz(iz,ix) := sxx,szz[iz,ix]^it+0.5,it+1.5,...
    !      vx(iz,ix) := vx[iz,ix-0.5]^it,it+1,...
    !      vz(iz,ix) := vy[iz-0.5,ix]^it,it+1,...
    !     sxz(iz,ix) := sxz[iz-0.5,ix-0.5]^it,it+1,...
    ! in space:
    ! <---*---<  
    ! |   |   |  * integer point   - sxx,szz - lda+2mu,lda 
    ! x---v---x  < half in x dir   - vx      - buox  
    ! |   |   |  v half in z dir   - vz      - buoz
    ! <---*---<  x half in x&z dir - sxz     - mu
    
    !add RHS to v^it
    subroutine put_velocities(time_dir,it,w,f)
        integer :: time_dir,it
        real :: w
        type(t_field) :: f
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        
        source_term=time_dir*w
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (2)
                f%vx(ifz:ilz,ifx:ilx) = f%vx(ifz:ilz,ifx:ilx) + source_term *lm%buox(ifz:ilz,ifx:ilx) *shot%src%interp_coeff(:,:,1)
                
                case (4)
                f%vz(ifz:ilz,ifx:ilx) = f%vz(ifz:ilz,ifx:ilx) + source_term *lm%buoz(ifz:ilz,ifx:ilx) *shot%src%interp_coeff(:,:,1)
                
            end select
            
        else
            select case (shot%src%icomp)
                case (2) !horizontal x force on vx[iz,ix-0.5,1]
                f%vx(iz,ix) = f%vx(iz,ix) + source_term*lm%buox(iz,ix)
                
                case (4) !vertical force     on vz[iz-0.5,ix,1]
                f%vz(iz,ix) = f%vz(iz,ix) + source_term*lm%buoz(iz,ix)
                
            end select
            
        endif
        
    end subroutine
    
    !v^it -> v^it+1 by FD of s^it+0.5
    subroutine update_velocities(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl

        ifz=bl(1)+2
        ilz=bl(2)-2
        ifx=bl(3)+2
        ilx=bl(4)-2
        !1=bl(5)+2
        !1=bl(6)-2

        if(m%if_freesurface) ifz=max(ifz,1)
        
        call fd2d_flat_velocities(f%vx,f%vz,f%sxx,f%szz,f%sxz,                                &
                                  f%cpml_dsxx_dx,f%cpml_dszz_dz,f%cpml_dsxz_dx,f%cpml_dsxz_dz,&
                                  lm%buox,lm%buoz,                                            &
                                  ifz,ilz,ifx,ilx,time_dir*shot%src%dt)
        

        dz_dx = m%dz/m%dx

        !apply free surface boundary condition if needed
        !free surface is located at [1,ix,1] level
        !dszz=0 -> -(lda+2mu)dvz_dz = lda dvx_dx
        !       -> -ldap2mu[1,ix]*(vz[1.5,ix]-vz[0.5,ix])/dz = lda[1,ix]*(vx[1,ix-0.5]-vx[1,ix+0.5])/dx
        !       -> -ldap2mu(1,ix)*(vz(2  ,ix)-vz(1  ,ix))/dz = lda(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
        !       ->  ldap2mu(1,ix)*(vz(1  ,ix)-vz(2  ,ix))/dz = lda(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
        !dsxz=0 -> -dvx_dz = dvz_dx
        !dsxz=0 -> -(vx[2,ix+0.5]-vx[0,ix+0.5])/2dz = ( (vz[0.5,ix]+vz[1.5,ix])/2 - (vz[0.5,ix+1]+vz[1.5,ix+1])/2 )/dx
        !       -> -(vx(2,ix+1  )-vx(0,ix+1  ))/2dz = ( (vz(1  ,ix)+vz(2  ,ix))/2 - (vz(1  ,ix+1)+vz(2  ,ix+1))/2 )/dx
        !       ->  (vx(0,ix+1  )-vx(2,ix+1  ))/ dz = ( (vz(1  ,ix)+vz(2  ,ix))   -  vz(1  ,ix+1)-vz(2  ,ix+1)    )/dx
        if(m%if_freesurface) then
            do ix=ifx,ilx
                f%vz(1,ix)= f%vz(2,ix) + lm%lda(1,ix)*(f%vx(1,ix)-f%vx(1,ix+1))*dz_dx/lm%ldap2mu(1,ix)
                !f%vx(0,ix)= f%vx(2,ix) + (f%vz(1,ix)+f%vz(2,ix)-f%vz(1,ix+1)-f%vz(2,ix+1))*dz_dx/lm%ldap2mu(1,ix) !toy2del says this condition is not needed
            enddo
        endif
        
    end subroutine
    
    !add RHS to s^it+0.5
    subroutine put_stresses(time_dir,it,w,f)
        integer :: time_dir
        real :: w
        type(t_field) :: f
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        
        source_term=time_dir*w
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (1)
                f%sxx(ifz:ilz,ifx:ilx) = f%sxx(ifz:ilz,ifx:ilx) + source_term *shot%src%interp_coeff(:,:,1)
                f%szz(ifz:ilz,ifx:ilx) = f%szz(ifz:ilz,ifx:ilx) + source_term *shot%src%interp_coeff(:,:,1)
            end select
            
        else
            select case (shot%src%icomp)
                case (1) !explosion on s[iz,ix,1]
                f%sxx(iz,ix) = f%sxx(iz,ix) + source_term
                f%szz(iz,ix) = f%szz(iz,ix) + source_term
            end select
            
        endif
        
    end subroutine
    
    !s^it+0.5 -> s^it+1.5 by FD of v^it+1
    subroutine update_stresses(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl
        
        ifz=bl(1)+2
        ilz=bl(2)-2
        ifx=bl(3)+2
        ilx=bl(4)-2
        !1=bl(5)+2
        !1=bl(6)-2

        if(m%if_freesurface) ifz=max(ifz,1)
        
        call fd2d_flat_stresses(f%vx,f%vz,f%sxx,f%szz,f%sxz,                            &
                                f%cpml_dvx_dx,f%cpml_dvz_dz,f%cpml_dvx_dz,f%cpml_dvz_dx,&
                                lm%ldap2mu,lm%lda,lm%mu,                                &
                                ifz,ilz,ifx,ilx,time_dir*shot%src%dt)

        !apply free surface boundary condition if needed
        !free surface is located at [1,ix] level
        !so explicit boundary condition: szz(1,ix)=0
        !and antisymmetric mirroring: sxz[0.5,ix-0.5]=-sxz[1.5,ix-0.5] -> sxz(1,ix)=-sxz(2,ix)
        if(m%if_freesurface) then
            f%szz(1,:)=0.
            f%sxz(1,:)=-f%sxz(2,:)
        endif
        
    end subroutine
    
    !get v^it+1 or s^it+1.5
    subroutine get_field(f,seismo)
        type(t_field) :: f
        real,dimension(*) :: seismo
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (1)
                    seismo(ircv)=sum( (factor_xx*f%sxx(ifz:ilz,ifx:ilx) + factor_zz*f%szz(ifz:ilz,ifx:ilx)) *shot%rcv(ircv)%interp_coeff(:,:,1) )
                    case (2)
                    seismo(ircv)=sum( f%vx(ifz:ilz,ifx:ilx) *shot%rcv(ircv)%interp_coeff(:,:,1) )
                    case (4)
                    seismo(ircv)=sum( f%vz(ifz:ilz,ifx:ilx) *shot%rcv(ircv)%interp_coeff(:,:,1) )
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (1) !sxx,szz[iz,ix,1]
                    seismo(ircv) = factor_xx * f%sxx(iz,ix) + factor_zz * f%szz(iz,ix)
                    case (2) !vx[iz,ix-0.5,1]
                    seismo(ircv)=f%vx(iz,ix)
                    case (4) !vz[iz-0.5,ix,1]
                    seismo(ircv)=f%vz(iz,ix)
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    
    !========= adjoint propagation =================
    !WE:        du_dt   = MDu   + src
    !->        Ndu_dt   = Du    + Nsrc, N=M^-1
    !Adjoint:  Ndv_dt^T = D^Tv  + adjsrc
    !->        -dv_dt   =-MDv   + Mr
    !     (reverse time) (negative derivatives)
    !
    !Discrete form (staggered grid in space and time, 2D as example):
    !    [ vx^it+1  ] [ vx^it    ]   [            dx(-)      dz(-)][ vx^it+1  ]
    ! N  | vz^it+1  | | vz^it    |   |                 dz(-) dx(-)|| vz^it+1  |
    !----|sxx^it+1.5|-|sxx^it+0.5| = |dx(+)                       ||sxx^it+0.5|  +Nsrc
    ! dt |szz^it+1.5| |szz^it+0.5|   |      dz(+)                 ||szz^it+0.5|
    !    [sxz^it+1.5] [sxz^it+0.5]   [dz(+) dx(+)                 ][sxz^it+0.5]
    !dx(-):=a1(s(ix)-s(ixm1))+a2(s(ixp1)-s(ixm2))  (O(x4))
    !dx(+):=a1(v(ixp1)-v(ix))+a2(v(ixp2)-v(ixm1))  (O(x4))
    !Adjoint:
    !    [ vx^it    ] [ vx^it+1  ]   [              dx(+)^T        dz(+)^T][ vx^it+1  ]
    ! N  | vz^it    |-| vz^it+1  |   |                     dz(+)^T dx(+)^T|| vz^it+1  |
    !----|sxx^it+0.5| [sxx^it+1.5] = |dx(-)^T                             ||sxx^it+0.5]  +Nsrc
    ! dt |szz^it+0.5| [szz^it+1.5]   |        dz(-)^T                     ||szz^it+0.5]
    !    [sxz^it+0.5] [sxz^it+1.5]   [dz(-)^T dx(-)^T                     ][sxz^it+0.5]
    !dx(+)^T=a1(s(ixm1)-s(ix))+a2(s(ixm2)-s(ixp1))=-dx(-)
    !dx(-)^T=a1(v(ix)-v(ixp1))+a2(v(ixm1)-v(ixp2))=-dx(+)
    !-dx,-dz can be regarded as -dt, thus allowing to use same code with flip of dt sign.
    !transposes of time derivatives centered on v^it+0.5 and s^it+1
    !so v^it+1   - v^it     => v^it - v^it+1
    !   s^it+1.5 - s^it+0.5 => s^it+0.5 - s^it+1.5
    !
    !In each time step:
    !d=RGAsrc, src:source term, A:put source, G:propagator, R:sampling at receivers
    !d=[vx vz p]^T, R=[I  0   0 ], A=[diag2(b) 0], src=[fx fz p]^T
    !                 |0 2/3 1/3|    [   0     1]
    !                 [   0    1]
    !Adjoint: src=A^T N G^T M R^T d
    !M R^T:put adjoint sources, 
    !G^T:  adjoint propagator,
    !A^T N:get adjoint fields
    
    
    !add RHS to s^it+1.5
    subroutine put_stresses_adjoint(time_dir,it,adjsource,f)
        integer :: time_dir
        real,dimension(*) :: adjsource
        type(t_field) :: f
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            
            tmp=adjsource(ircv)
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (1) !pressure adjsource
                    f%sxx(ifz:ilz,ifx:ilx) = f%sxx(ifz:ilz,ifx:ilx)  &
                        +tmp*( factor_xx*lm%ldap2mu(ifz:ilz,ifx:ilx) &
                              +factor_zz*lm%lda(ifz:ilz,ifx:ilx) )   &
                        *shot%rcv(ircv)%interp_coeff(:,:,1)
        
                    f%szz(ifz:ilz,ifx:ilx) = f%szz(ifz:ilz,ifx:ilx)    &
                        +tmp*( factor_xx*lm%lda(ifz:ilz,ifx:ilx)       &
                              +factor_zz*lm%ldap2mu(ifz:ilz,ifx:ilx) ) &
                        *shot%rcv(ircv)%interp_coeff(:,:,1)
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (1) !pressure adjsource
                    !sxx[iz,ix]
                    f%sxx(iz,ix) = f%sxx(iz,ix)            &
                        +tmp*( factor_xx*lm%ldap2mu(iz,ix) &
                              +factor_zz*lm%lda(iz,ix) )
                    
                    !szz[iz,ix]
                    f%szz(iz,ix) = f%szz(iz,ix)            &
                        +tmp*( factor_xx*lm%lda(iz,ix)     &
                              +factor_zz*lm%ldap2mu(iz,ix) )
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    !s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
    subroutine update_stresses_adjoint(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl
        
        !same code with -dt
        call update_stresses(time_dir,it,f,bl)
        
    end subroutine
    
    !add RHS to v^it+1
    subroutine put_velocities_adjoint(time_dir,it,adjsource,f)
        integer :: time_dir
        real,dimension(*) :: adjsource
        type(t_field) :: f
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            
            tmp=adjsource(ircv)
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (2) !horizontal x adjsource
                    f%vx(ifz:ilz,ifx:ilx) = f%vx(ifz:ilz,ifx:ilx) + tmp*lm%buox(ifz:ilz,ifx:ilx) *shot%rcv(ircv)%interp_coeff(:,:,1)
                    
                    case (4) !horizontal z adjsource
                    f%vz(ifz:ilz,ifx:ilx) = f%vz(ifz:ilz,ifx:ilx) + tmp*lm%buoz(ifz:ilz,ifx:ilx) *shot%rcv(ircv)%interp_coeff(:,:,1)
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (2) !horizontal x adjsource
                    !vx[ix-0.5,1,iz]
                    f%vx(iz,ix) = f%vx(iz,ix) + tmp*lm%buox(iz,ix)
                    
                    case (4) !horizontal z adjsource
                    !vz[ix,1,iz-0.5]
                    f%vz(iz,ix) = f%vz(iz,ix) + tmp*lm%buoz(iz,ix)
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    !v^it+1 -> v^it by FD^T of s^it+0.5
    subroutine update_velocities_adjoint(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl
        
        !same code with -dt
        call update_velocities(time_dir,it,f,bl)
        
    end subroutine
    
    !get v^it or s^it+0.5
    subroutine get_field_adjoint(f,w)
        type(t_field) :: f
        real w
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (1)
                w = sum( &
        ( lm%ldap2mu(ifz:ilz,ifx:ilx)*f%sxx(ifz:ilz,ifx:ilx) -lm%lda    (ifz:ilz,ifx:ilx)*f%szz(ifz:ilz,ifx:ilx)   &
         -lm%lda    (ifz:ilz,ifx:ilx)*f%sxx(ifz:ilz,ifx:ilx) +lm%ldap2mu(ifz:ilz,ifx:ilx)*f%szz(ifz:ilz,ifx:ilx) ) &
        *lm%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)  *shot%src%interp_coeff(:,:,1) )
                
                case (2)
                w = sum( f%vx(ifz:ilz,ifx:ilx) *shot%src%interp_coeff(:,:,1) )
                
                case (4)
                w = sum( f%vz(ifz:ilz,ifx:ilx) *shot%src%interp_coeff(:,:,1) )
                
            end select
            
        else
            select case (shot%src%icomp)
                case (1) !s[iz,ix,1]
                w = ( lm%ldap2mu(iz,ix)*f%sxx(iz,ix) -lm%lda    (iz,ix)*f%szz(iz,ix)   &
                     -lm%lda    (iz,ix)*f%sxx(iz,ix) +lm%ldap2mu(iz,ix)*f%szz(iz,ix) ) &
                    *lm%inv_ldapmu_4mu(iz,ix)
                
                case (2) !vx[iz,ix-0.5,1]
                w=f%vx(iz,ix)
                
                case (4) !vz[iz-0.5,ix,1]
                w=f%vz(iz,ix)
            end select
            
        endif
        
    end subroutine
    
    !========= for wavefield correlation ===================
    !dC/dm = \int a  dN/dm du/dt
    !      = \int a  dN/dm M Du
    !      = \int a -N dM/dm Du
    !                            [0                    ]
    !                            |  0                  |
    !-NdM/dlda = -1/(lda+mu)/4mu |    lda+2mu  -lda    |
    !                            |     -lda   lda+2mu  |
    !                            [                    0]
    ! => corr_lda := [adjsxx adjszz][lda+2mu  -lda  ][dvx_dx]
    !                               [ -lda   lda+2mu][dvz_dz]
    !                           [0                          ]
    !                           |  0                        |
    !-NdM/dmu = -2/(lda+mu)/4mu |    lda+2mu                |
    !                           |           lda+2mu         |
    !                           [                  2(lda+mu)]
    !                                     [lda+2mu                ][   dvx_dx    ]
    ! => corr_mu := [adjsxx adjszz adjsxz][       lda+2mu         ][   dvz_dz    ]
    !                                     [              2(lda+mu)][dvx_dz+dvz_dx]
    !            [b    ]
    !dN/drho M = |  b  |
    !            [    0]
    ! => corr_rho := [adjvx adjvz][sxx_dx+sxz_dz]
    !                             [szz_dz+sxz_dx]
    !
    !For gmu, we need to spatially interpolate adjsxz and dvx_dz+dvz_dx.
    !For grho, we need to spatially interpolate adjvx, adjvz and sxx_dx+sxz_dz, szz_dz+sxz_dx
    !
    !For time, we just use v^it+1, s^it+0.5, adjs^it+0.5
    
    subroutine field_correlation_moduli(it,sf,rf,sb,rb,corr)
        type(t_field) :: sf,rf
        integer,dimension(6) :: sb,rb
        real,dimension(cb%mz,cb%mx,ncorr) :: corr
        
        !nonzero only when sf touches rf
        ifz=max(sb(1),rb(1),2)
        ilz=min(sb(2),rb(2),cb%mz-2)
        ifx=max(sb(3),rb(3),2)
        ilx=min(sb(4),rb(4),cb%mx-2)
        1=max(sb(5),rb(5),2)
        1=min(sb(6),rb(6),cb%my-2)
        
        call corr2d_flat_moduli(sf%vx,sf%vz,             &
                                rf%sxx,rf%szz,rf%sxz,    &
                                lm%ldap2mu,lm%two_ldapmu,&
                                corr(:,:,1),corr(:,:,2), &
                                ifz,ilz,ifx,ilx)
        
    end subroutine
    
    subroutine field_correlation_density(it,sf,rf,sb,rb,corr)
        type(t_field) :: sf,rf
        integer,dimension(6) :: sb,rb
        real,dimension(cb%mz,cb%mx,ncorr) :: corr
        
        !nonzero only when sf touches rf
        ifz=max(sb(1),rb(1),2)
        ilz=min(sb(2),rb(2),cb%mz-2)
        ifx=max(sb(3),rb(3),2)
        ilx=min(sb(4),rb(4),cb%mx-2)
        1=max(sb(5),rb(5),2)
        1=min(sb(6),rb(6),cb%my-2)
        
        call corr2d_flat_density(sf%sxx,sf%szz,sf%sxz,&
                                 rf%vx,rf%vz,         &
                                 corr(:,:,3),         &
                                 ifz,ilz,ifx,ilx)

    end subroutine
    
    subroutine field_correlation_scaling(corr)
        real,dimension(cb%mz,cb%mx,ncorr) :: corr
        
        !corr_lda
        corr(:,:,1)=corr(:,:,1) / (-lm%two_ldapmu(1:cb%mz,1:cb%mx))

        !corr_mu
        corr(:,:,2)=corr(:,:,2) * (-2.*lm%inv_ldapmu_4mu(1:cb%mz,1:cb%mx))
        
        !corr_rho
        corr(:,:,3)=corr(:,:,3) / cb%rho(1:cb%mz,1:cb%mx,1)

        corr(1,:,:) = corr(3,:,:)
        corr(2,:,:) = corr(3,:,:)
        corr(cb%mz-1,:,:) = corr(cb%mz-2,:,:)
        corr(cb%mz,  :,:) = corr(cb%mz-2,:,:)

        corr(:,1,:) = corr(:,3,:)
        corr(:,2,:) = corr(:,3,:)
        corr(:,cb%mx-1,:) = corr(:,cb%mx-2,:)
        corr(:,cb%mx  ,:) = corr(:,cb%mx-2,:)
        
        !the unit of glda,gmu should be [m3], grho should be [m5/s2]
        !after multiplied by (lda_max-lda_min) or (rho_max-rho_min) (to be done in m_parameterization.f90)
        !the unit of parameter update is [Nm], same as Lagrangian
        !and since the unit of gradient scaling factor is [1/N/m] (in m_scaling.f90)
        !then the parameter update becomes unitless
        corr=corr*m%cell_volume*shot%src%dt
        
    end subroutine
    
    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd2d_flat_velocities(vx,vz,sxx,szz,sxz,                                  &
                                    cpml_dsxx_dx,cpml_dszz_dz,cpml_dsxz_dx,cpml_dsxz_dz,&
                                    buox,buoz,                                          &
                                    ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vx,vz,sxx,szz,sxz
        real,dimension(*) :: cpml_dsxx_dx,cpml_dszz_dz,cpml_dsxz_dx,cpml_dsxz_dz
        real,dimension(*) :: buox,buoz
        
        nz=cb%nz
        nx=cb%nx
        
        dsxx_dx=0.
        dsxz_dz=0.
        dsxz_dx=0.
        dszz_dz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dsxx_dx,dsxz_dz,dsxz_dx,dszz_dz)
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
                

                dsxx_dx= c1x*(sxx(iz_ix)-sxx(iz_ixm1)) +c2x*(sxx(iz_ixp1)-sxx(iz_ixm2))
                dsxz_dz= c1z*(sxz(izp1_ix)-sxz(iz_ix)) +c2z*(sxz(izp2_ix)-sxz(izm1_ix))

                dsxz_dx= c1x*(sxz(iz_ixp1)-sxz(iz_ix)) +c2x*(sxz(iz_ixp2)-sxz(iz_ixm1))
                dszz_dz= c1z*(szz(iz_ix)-szz(izm1_ix)) +c2z*(szz(izp1_ix)-szz(izm2_ix))
                                
                !cpml
                cpml_dsxx_dx(i)= cb%b_x_half(ix)*cpml_dsxx_dx(i) + cb%a_x_half(ix)*dsxx_dx
                cpml_dsxz_dz(i)= cb%b_z(iz)     *cpml_dsxz_dz(i) + cb%a_z(iz)     *dsxz_dz 
                cpml_dsxz_dx(i)= cb%b_x(ix)     *cpml_dsxz_dx(i) + cb%a_x(ix)     *dsxz_dx
                cpml_dszz_dz(i)= cb%b_z_half(iz)*cpml_dszz_dz(i) + cb%a_z_half(iz)*dszz_dz

                dsxx_dx=dsxx_dx*cb%kappa_x_half(ix) + cpml_dsxx_dx(i)  !kappa's should have been inversed in m_computebox.f90
                dsxz_dz=dsxz_dz*cb%kappa_z(iz)      + cpml_dsxz_dz(i)
                dsxz_dx=dsxz_dx*cb%kappa_x(ix)      + cpml_dsxz_dx(i)
                dszz_dz=dszz_dz*cb%kappa_z_half(iz) + cpml_dszz_dz(i)
                
                !velocity
                vx(i)=vx(i) + dt*buox(i)*(dsxx_dx+dsxz_dz)
                vz(i)=vz(i) + dt*buoz(i)*(dsxz_dx+dszz_dz)
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd2d_flat_stresses(vx,vz,sxx,szz,sxz,                              &
                                  cpml_dvx_dx,cpml_dvz_dz,cpml_dvx_dz,cpml_dvz_dx,&
                                  ldap2mu,lda,mu,                                 &
                                  ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vx,vz,sxx,szz,sxz
        real,dimension(*) :: cpml_dvx_dx,cpml_dvz_dz,cpml_dvx_dz,cpml_dvz_dx
        real,dimension(*) :: ldap2mu,lda,mu
        
        nz=cb%nz
        nx=cb%nx
        
        dvx_dx=0.
        dvz_dz=0.
        dvx_dz=0.
        dvz_dx=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dvx_dx,dvz_dz,&
        !$omp         dvx_dz,dvz_dx)
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
                

                dvx_dx= c1x*(vx(iz_ixp1)-vx(iz_ix))  +c2x*(vx(iz_ixp2)-vx(iz_ixm1))
                dvz_dz= c1z*(vz(izp1_ix)-vz(iz_ix))  +c2z*(vz(izp2_ix)-vz(izm1_ix))
                
                !cpml
                cpml_dvx_dx(i)=cb%b_x(ix)*cpml_dvx_dx(i)+cb%a_x(ix)*dvx_dx
                cpml_dvz_dz(i)=cb%b_z(iz)*cpml_dvz_dz(i)+cb%a_z(iz)*dvz_dz

                dvx_dx=dvx_dx*cb%kappa_x(ix) + cpml_dvx_dx(iz_ix)
                dvz_dz=dvz_dz*cb%kappa_z(iz) + cpml_dvz_dz(iz_ix)
                
                !normal stresses
                sxx(i) = sxx(i) + dt * (ldap2mu(i)*dvx_dx+lda(i)    *dvz_dz)
                szz(i) = szz(i) + dt * (lda(i)    *dvx_dx+ldap2mu(i)*dvz_dz)


                dvx_dz= c1z*(vx(iz_ix)-vx(izm1_ix))  +c2z*(vx(izp1_ix)-vx(izm2_ix))
                dvz_dx= c1x*(vz(iz_ix)-vz(iz_ixm1))  +c2x*(vz(iz_ixp1)-vz(iz_ixm2))

                !cpml
                cpml_dvx_dz(i)=cb%b_z_half(iz)*cpml_dvx_dz(i)+cb%a_z_half(iz)*dvx_dz
                cpml_dvz_dx(i)=cb%b_x_half(ix)*cpml_dvz_dx(i)+cb%a_x_half(ix)*dvz_dx

                dvx_dz=dvx_dz*cb%kappa_z_half(iz) + cpml_dvx_dz(i)
                dvz_dx=dvz_dx*cb%kappa_x_half(ix) + cpml_dvz_dx(i)

                !shear stress
                sxz(i) = sxz(i) + dt * mu(i)*(dvx_dz+dvz_dx)
                
            enddo
            
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine
    
    subroutine corr2d_flat_moduli(sf_vx,sf_vz,         &
                                  rf_sxx,rf_szz,rf_sxz,&
                                  ldap2mu,two_ldapmu,  &
                                  corr_lda,corr_mu,    &
                                  ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_vx,sf_vz
        real,dimension(*) :: rf_sxx,rf_szz,rf_sxz
        real,dimension(*) :: ldap2mu,two_ldapmu
        real,dimension(*) :: corr_lda,corr_mu
        
        nz=cb%nz
        
        sf_dvx_dx=0.
        sf_dvz_dz=0.
        sf_4_dvzdx_p_dvxdz=0.
        rf_4_sxz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,izp1_ixm2,iz_ixm1,izp1_ixm1,&
        !$omp         izm2_ixp1,izm1_ixp1,iz_ixp1,izp1_ixp1,izp2_ixp1,&
        !$omp         iz_ixp2,izp1_ixp2,&
        !$omp         sf_dvx_dx,sf_dvz_dz,&
        !$omp         sf_4_dvzdx_p_dvxdz,rf_4_sxz)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !corr has no boundary layers

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

                corr_lda(j)=corr_lda(j) + rf_sxx(i)*sf_dvx_dx + rf_sxx(i)*sf_dvz_dz &
                                        + rf_szz(i)*sf_dvx_dx + rf_szz(i)*sf_dvz_dz


                sf_4_dvzdx_p_dvxdz = &  !(dvz_dx+dvx_dz)(iz,ix)     [iz-0.5,ix-0.5]
                                         c1x*(sf_vz(iz_ix)-sf_vz(iz_ixm1)) +c2x*(sf_vz(iz_ixp1)-sf_vz(iz_ixm2)) &
                                       + c1z*(sf_vx(iz_ix)-sf_vx(izm1_ix)) +c2z*(sf_vx(izp1_ix)-sf_vx(izm2_ix)) &
                                     & &
                                     & & !(dvz_dx+dvx_dz)(iz,ix+1)   [iz-0.5,ix+0.5]
                                       + c1x*(sf_vz(iz_ixp1)-sf_vz(iz_ix    )) +c2x*(sf_vz(iz_ixp2  )-sf_vz(iz_ixm1  )) &
                                       + c1z*(sf_vx(iz_ixp1)-sf_vx(izm1_ixp1)) +c2z*(sf_vx(izp1_ixp1)-sf_vx(izm2_ixp1)) &
                                     & &
                                     & & !(dvz_dx+dvx_dz)(iz+1,ix)   [iz+0.5,ix-0.5]
                                       + c1x*(sf_vz(izp1_ix)-sf_vz(izp1_ixm1)) +c2x*(sf_vz(izp1_ixp1)-sf_vz(izp1_ixm2)) &
                                       + c1z*(sf_vx(izp1_ix)-sf_vx(iz_ix    )) +c2z*(sf_vx(izp2_ix  )-sf_vx(izm1_ix  )) &
                                     & &
                                     & & !(dvz_dx+dvx_dz)(iz+1,ix+1) [iz+0.5,ix+0.5]
                                       + c1x*(sf_vz(izp1_ixp1)-sf_vz(izp1_ix)) +c2x*(sf_vz(izp1_ixp2)-sf_vz(izp1_ixm1)) &
                                       + c1z*(sf_vx(izp1_ixp1)-sf_vx(iz_ixp1)) +c2z*(sf_vx(izp2_ixp1)-sf_vx(izm1_ixp1))

                        ![iz-0.5,ix-0.5]   [iz+0.5,ix-0.5]   [iz-0.5,ix+0.5]     [iz+0.5,ix+0.5]
                rf_4_sxz = rf_sxz(iz_ix) + rf_sxz(izp1_ix) + rf_sxz(iz_ixp1) + rf_sxz(izp1_ixp1)

                corr_mu(j)=corr_mu(j) + ldap2mu(i)*rf_sxx(i)*sf_dvx_dx &
                                      + ldap2mu(i)*rf_szz(i)*sf_dvz_dz &
                                      + two_ldapmu(i)*0.0625*rf_4_sxz*sf_4_dvzdx_p_dvxdz   !0.0625=1/16
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine corr2d_flat_density(sf_sxx,sf_szz,sf_sxz,&
                                   rf_vx,rf_vz,         &
                                   corr,                &
                                   ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_sxx,sf_szz,sf_sxz
        real,dimension(*) :: rf_vx,rf_vz
        real,dimension(*) :: corr
        
        nz=cb%nz
        
        sf_2_dsxxdx_p_dsxzdz=0.
        sf_2_dsxzdx_p_dszzdz=0.
        rf_2vx=0.
        rf_2vz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,izp1_ixm1,&
        !$omp         izm1_ixp1,iz_ixp1,izp1_ixp1,izp2_ixp1,&
        !$omp         iz_ixp2,izp1_ixp2,&
        !$omp         sf_2_dsxxdx_p_dsxzdz,sf_2_dsxzdx_p_dszzdz,&
        !$omp         rf_2vx,rf_2vz)
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
                
                sf_2_dsxxdx_p_dsxzdz = &   !(dsxx_dx+dsxz_dz)(iz,ix)   [iz,ix-0.5]
                                           c1x*(sf_sxx(iz_ix  )-sf_sxx(iz_ixm1)) +c2x*(sf_sxx(iz_ixp1)-sf_sxx(iz_ixm2)) &
                                         + c1z*(sf_sxz(izp1_ix)-sf_sxz(iz_ix  )) +c2z*(sf_sxz(izp2_ix)-sf_sxz(izm1_ix)) &
                                       & &
                                       & & !(dsxx_dx+dsxz_dz)(iz,ix+1) [iz,ix+0.5]
                                         + c1x*(sf_sxx(iz_ixp1  )-sf_sxx(iz_ix  )) +c2x*(sf_sxx(iz_ixp2  )-sf_sxx(iz_ixm1  )) &
                                         + c1z*(sf_sxz(izp1_ixp1)-sf_sxz(iz_ixp1)) +c2z*(sf_sxz(izp2_ixp1)-sf_sxz(izm1_ixp1))

                sf_2_dsxzdx_p_dszzdz = &  !(dsxz_dx+dszz_dz)(iz,ix)   [iz-0.5,ix]
                                           c1x*(sf_sxz(iz_ixp1)-sf_sxz(iz_ix  )) +c2x*(sf_sxz(iz_ixp2)-sf_sxz(iz_ixm1)) &
                                         + c1z*(sf_szz(iz_ix  )-sf_szz(izm1_ix)) +c2z*(sf_szz(izp1_ix)-sf_szz(izm1_ix)) &
                                       & &
                                       & & !(dsxz_dx+dszz_dz)(iz+1,ix) [iz+0.5,ix]
                                         + c1x*(sf_sxz(izp1_ixp1)-sf_sxz(izp1_ix)) +c2x*(sf_sxz(izp1_ixp2)-sf_sxz(izp1_ixm1)) &
                                         + c1z*(sf_szz(izp1_ix  )-sf_szz(iz_ix  )) +c2z*(sf_szz(izp2_ix  )-sf_szz(iz_ix    ))
                                         
                         ![iz,ix-0.5]      [iz,ix+0.5]
                rf_2vx = rf_vx(iz_ix) + rf_vx(iz_ixp1)
                         ![iz-0.5,ix]      [iz+0.5,ix]
                rf_2vz = rf_vz(iz_ix) + rf_vz(izp1_ix)
                
                corr(j)=corr(j) + 0.25*( rf_2vx*sf_2_dsxxdx_p_dsxzdz + rf_2vz*sf_2_dsxzdx_p_dszzdz )
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    

end
