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
            'Time-domain ISOtropic 2D PSV (ELastic) propagation'//s_NL// &
            '2nd-order Displacement-Strain formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.5*Vmax/dx'//s_NL// &
            'Required model attributes: vp, vs, rho'//s_NL// &
            'Required field components: uz, ux, uz_prev, ux_prev, uz_next, ux_next'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: grho(wait) glda gmu'

        integer :: nbndlayer=max(2,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=3 !number of basic gradients
        ! integer :: nimag=3 !number of basic images
        !integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:),allocatable :: buoz, buox, ldap2mu, lda, mu
        real,dimension(:,:),allocatable :: inv_ldapmu_4mu!, ldapmu

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
        
        procedure :: inject_displacement
        procedure :: inject_strains
        procedure :: update_displacement
        procedure :: update_strains
        procedure :: evolve
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    logical,public :: if_propagator_record_adjseismo=.false.

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
            call alloc(m%vp,m%nz,m%nx,1,o_init=1500.)
            call warn('Constant vp model (1500 m/s) is allocated by propagator.')
        endif

        if(index(self%info,'vs')>0  .and. .not. allocated(m%vs)) then
            call alloc(m%vs,m%nz,m%nx,1)
            m%vs=m%vp/sqrt(3.)
            call warn('Poisson solid (vs=vp/√3) is assumed by propagator.')
        endif

        if(index(self%info,'rho')>0 .and. .not. allocated(m%rho)) then
            call alloc(m%rho,m%nz,m%nx,1,o_init=1.)
            call warn('Constant rho model (1 kg/m³) is allocated by propagator.')
        endif

        if(m%is_freesurface) then
            call warn('Sorry, free surface has NOT yet considered in this propagator. Switch off free-surface condition.')
            m%is_freesurface=.false.
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

        real,dimension(:,:),allocatable :: temp_mu

        c1x=coef(1)/m%dx; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2z=coef(2)/m%dz

        dt2=self%dt**2
        inv_2dt =1./2/self%dt
        inv_2dz =1./2/m%dz
        inv_2dx =1./2/m%dx
        
        wavelet_scaler=dt2/m%cell_volume

        if_hicks=shot%if_hicks

        call alloc(self%buoz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%buox,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%ldap2mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%lda,    [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%mu,             [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%inv_ldapmu_4mu, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        call alloc(temp_mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        self%ldap2mu(:,:)=cb%rho(:,:,1)*cb%vp(:,:,1)**2
             temp_mu(:,:)=cb%rho(:,:,1)*cb%vs(:,:,1)**2

        self%lda=self%ldap2mu-2.*temp_mu
        ! if(mpiworld%is_master) then
        ! write(*,*) 'self%ldap2mu sanity:', minval(self%ldap2mu),maxval(self%ldap2mu)
        ! write(*,*) 'self%lda     sanity:', minval(self%lda),maxval(self%lda)
        ! endif

        self%inv_ldapmu_4mu=0.25/(self%lda+temp_mu)/temp_mu

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
        call f%init_boundary_displacement

        call alloc(f%uz     , [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%uz_prev, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%uz_next, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%ux     , [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%ux_prev, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%ux_next, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%duz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dux_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dux_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%duz_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
    
        call alloc(f%duz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dux_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%duz_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dux_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
    
        call alloc(f%lapz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%lapx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%ez,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%ex,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%es,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%dez_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dex_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dex_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dez_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%des_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%des_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

    end subroutine

    subroutine init_correlate(self,corr,name)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name
        
        corr%name=name

        ! if(name(1:1)=='g') then !gradient components
            ! call alloc(corr%gikpa,m%nz,m%nx,m%ny)
            ! call alloc(corr%gbuo, m%nz,m%nx,m%ny)
            call alloc(corr%glda, m%nz,m%nx,m%ny)
            call alloc(corr%gmu,  m%nz,m%nx,m%ny)
        ! else !image components
            ! call alloc(corr%ipp,m%nz,m%nx,m%ny)
            ! call alloc(corr%ibksc,m%nz,m%nx,m%ny)
            ! call alloc(corr%ifwsc,m%nz,m%nx,m%ny)
        ! endif

    end subroutine

    subroutine init_abslayer(self)
        class(t_propagator) :: self

        call cpml%init

    end subroutine

    subroutine assemble(self,corr)
        class(t_propagator) :: self
        type(t_correlate) :: corr

        ! if(allocated(correlate_image)) then
        !     call correlate_assemble(corr%ipp,   correlate_image(:,:,:,1))
        !     call correlate_assemble(corr%ibksc, correlate_image(:,:,:,2))
        !     call correlate_assemble(corr%ifwsc, correlate_image(:,:,:,3))
        ! endif

        if(allocated(correlate_gradient)) then
            ! call correlate_assemble(corr%gikpa, correlate_gradient(:,:,:,1))
            ! call correlate_assemble(corr%gbuo,  correlate_gradient(:,:,:,2))
            call correlate_assemble(corr%glda,  correlate_gradient(:,:,:,2))
            call correlate_assemble(corr%gmu,   correlate_gradient(:,:,:,3))

        endif        
        
    end subroutine

    
    !========= Derivations =================
    !PDE:      A u = CP u - CDC u - f
    !Adjoint:  Aᵀa = CP a + CDC a - d
    !where
    !u=[uz ux ez ex es]ᵀ, f=[fz fx]ᵀδ(x-xs) with xs source position, d is recorded data
    !P=diag[ρ∂ₜₜ ρ∂ₜₜ 1 1 1]
    !  [1             ]     [ 0  0 ∂z  0 ∂ₓ]
    !  |  1           |     | 0  0  0 ∂ₓ ∂z|
    !C=|    λ+2μ  λ   |, D =|∂z  0  0  0  0|
    !  |    λ   λ+2μ  |     | 0 ∂ₓ  0  0  0|
    !  [             μ]     [∂ₓ ∂z  0  0  0]
    !a=[uzᵃ uxᵃ ezᵃ exᵃ esᵃ]ᵀ is the adjoint field
    !
    !Discrete case:
    !Meshing with staggered grids in time and space (2D example):
    !                   |    es    |   es -½ uz  es    |         |
    !                   |    μ     |   μ     bz  μ     |         |
    !                   |          |         |         |         |
    !   e  e  e        λ,μ   bx   λ,μ  bx   λ,μ  bx   λ,μ  bx   λ,μ
    !  -u--u--u-→ t    -en---ux---en---ux---en---ux---en---ux---en-→ x
    !  -1  0  1        -2   -1½   -1   -½    0    ½    1   1½    2    
    !                   |          |         |         |         | 
    !                   |    es    |   es  ½ uz  es    |         | 
    !                   |    μ     |   μ     bz  μ     |         | 
    !                   |          |         |         |         | 
    !                  -|----------|-------1-en--------|---------|-
    !                   |          |         κ         |         | 
    !                   |          |         |         |         | 
    !                   |          |      1½ uz        |         | 
    !                   |          |         bz        |         | 
    !                   |          |         |         |         | 
    !                                      z ↓
    !
    !Convention for half-integer index:
    !(array index)  =>     (real index)     
    !  uz(iz,ix)    =>   uz[iz-½,ix  ]^n   :=uz((iz-½)*dz,ix*dx,n*dt)
    !  ux(iz,ix)    =>   ux[iz,  ix-½]^n  
    !  en(iz,ix)    =>   en[iz,  ix  ]^n+½ :=en(iz*dz,    ix*dx,(n+½)*dt)
    !  es(iz,ix)    =>   es[iz-½,ix-½]^n  
    !
    !Forward:
    !FD eqn:
    !  [uz]      [ 0   0 ∂zᵇ  0 ∂ₓᶠ][      uz          ]
    !  |ux|      | 0   0  0  ∂ₓ ∂zᶠ||      ux          |
    !CP|ez|^n = C|∂zᶠ  0  0   0  0 ||(λ+2μ)ez +     λex|^n
    !  |ex|      | 0  ∂ₓᶠ 0   0  0 ||     λez +(λ+2μ)ex|
    !  [es]      [∂ₓᵇ ∂zᵇ 0   0  0 ][     μes          ]
    !where
    !∂ₜᶠ*dt := v^n+1 - v^n                             ~O(t²)
    !∂zᵇ*dz := c₁(s(iz  )-s(iz-1) +c₂(s(iz+1)-s(iz-2)  ~O(x⁴)
    !∂zᶠ*dz := c₁(v(iz+1)-v(iz  ) +c₂(v(iz+2)-v(iz-1)  ~O(x⁴)
    
    !Step #1: u^n  += src
    !Step #2: save u^n to boundary values
    !Step #3: e^n = spatial FD(u^n)
    !Step #4: e^n += src
    !Step #5: lap^n = spatial FD(e^n)
    !Step #6: u^n+1 = 2u^n -u^n-1 + lap^n
    !Step #7: (u^n-1,u^n) = (u^n,u^n+1)
    !Step #8: sample u^n at receivers
    !
    !in reverse time:
    !Step #7: (u^n,u^n+1) = (u^n-1,u^n)
    !Step #2: load boundary values for u^n
    !Step #3: e^n = spatial FD(u^n)
    !Step #4: e^n -= src
    !Step #5: lap^n = spatial FD(e^n)
    !Step #6: u^n-1 = 2u^n -u^n+1 + lap^n
    !Step #1: u^n -= src
    !
    !Adjoint:
    !since
    !∂ₜ²ᵀ = ∂ₜ²
    !∂zᵇᵀ = u(iz)-u(iz+1) = -∂zᶠ, ∂zᶠᵀ = -∂zᵇ
    !FD eqn:
    !  [uz]      [ 0   0 ∂zᵇ  0 ∂ₓᶠ][      uz          ]
    !  |ux|      | 0   0  0  ∂ₓ ∂zᶠ||      ux          |
    !CP|ez|^n =-C|∂zᶠ  0  0   0  0 ||(λ+2μ)ez +     λex|^n
    !  |ex|      | 0  ∂ₓᶠ 0   0  0 ||     λez +(λ+2μ)ex|
    !  [es]      [∂ₓᵇ ∂zᵇ 0   0  0 ][     μes          ]
    !SAME as the discretized forward FD eqn!
    !
    !Time marching (in reverse time):
    !Step #1: uᵃ^n += adjsrc
    !Step #3: eᵃ^n = spatial FD(uᵃ^n)
    !Step #4: eᵃ^n += src
    !Step #5: lapᵃ^n = spatial FD(eᵃ^n)
    !Step #6: uᵃ^n+1 = 2uᵃ^n -uᵃ^n-1 + lapᵃ^n
    !Step #7: cross correlate
    !Step #8: (uᵃ^n-1,uᵃ^n) = (uᵃ^n,uᵃ^n+1)
    !Step #9: sample uᵃ^n at source      


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
                call fld_u%check_value(fld_u%uz)
            endif

            !step 1
            call cpu_time(tic)
            call self%inject_displacement(fld_u,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2
            call cpu_time(tic)
            call fld_u%boundary_transport_displacement('save',it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            !step 3
            call cpu_time(tic)
            call self%update_strains(fld_u,time_dir,it)
            call cpu_time(toc)
            tt3=tt3+toc-tic

            !step 4
            call cpu_time(tic)
            call self%inject_strains(fld_u,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 5 & 6
            call cpu_time(tic)
            call self%update_displacement(fld_u,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !step 7: evolve, it -> it+1
            call cpu_time(tic)
            call self%evolve(fld_u,time_dir,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !step 8
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
            write(*,*) 'Elapsed time to update strain           ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update displacement     ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to evolve field            ',tt6/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field           ',tt7/mpiworld%max_threads
        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_sfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_sfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')

    end subroutine

    subroutine adjoint(self,fld_a,fld_u,a_star_u)
    !adjoint_a_star_Du
        class(t_propagator) :: self
        type(t_field) :: fld_a,fld_u
        type(t_correlate) :: a_star_u

        real,parameter :: time_dir=-1. !time direction

        !reinitialize absorbing boundary for incident wavefield reconstruction
        call fld_u%reinit
        
        !for adjoint test
        if(if_propagator_record_adjseismo)  call alloc(fld_a%seismo,1,self%nt)

        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.
        
        ift=1; ilt=self%nt

        do it=ilt,ift,int(time_dir)
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_a%check_value(fld_a%uz)
                call fld_u%check_value(fld_u%uz)
            endif   

            !backward step 4
            call cpu_time(tic)
            call self%evolve(fld_u,time_dir,it)
            call cpu_time(toc)
            tt11=tt11+toc-tic

            !backward step 2
            call cpu_time(tic)
            call fld_u%boundary_transport_displacement('load',it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !backward step 3
            call cpu_time(tic)
            call self%update_strains(fld_u,time_dir,it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            !backward step 4
            call cpu_time(tic)
            call self%inject_strains(fld_u,time_dir,it)
            call cpu_time(toc)
            tt3=tt3+toc-tic

            !backward step 5 & 6
            call cpu_time(tic)
            call self%update_displacement(fld_u,time_dir,it)
            call cpu_time(toc)
            tt7=tt7+toc-tic

            !backward step 1
            call cpu_time(tic)
            call self%inject_displacement(fld_u,time_dir,it)
            call cpu_time(toc)
            tt8=tt8+toc-tic

            !--------------------------------------------------------!

            !adjoint step 1: inject to pz^it+1 at receivers
            call cpu_time(tic)
            call self%inject_displacement(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 3: e^it+1.5 -> e^it+0.5 by FD^T of pz^it+1
            call cpu_time(tic)
            call self%update_strains(fld_a,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !adjoint step 4: inject to e^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_strains(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 5 & 6: pz^it+1 -> pz^it by FD^T of e^it+0.5
            call cpu_time(tic)
            call self%update_displacement(fld_a,time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic

            !adjoint step 7: gkpa: rf%e^it+0.5 star D sf%s_dt^it+0.5
            !use sf%pz^it+1 to compute sf%s_dt^it+0.5, as backward step 4
            if(mod(it,irdt)==0) then
                call cpu_time(tic)
                call cross_correlate_glda_gmu(fld_a,fld_u,a_star_u,it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif

            !adjoint step 8
            ! this step is moved to update_pressure for easier management
            call cpu_time(tic)
            call self%evolve(fld_a,time_dir,it)
            call cpu_time(toc)
            tt11=tt11+toc-tic
            
            !adjoint step 9: sample pz^it or e^it+0.5 at source position
            if(if_propagator_record_adjseismo) then
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


    !forward: add RHS to pz^it
    !adjoint: add RHS to pz^it+1
    subroutine inject_displacement(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        if(.not. f%is_adjoint) then

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler
            
            if(if_hicks) then
                select case (shot%src%comp)
                case ('uz')
                    f%uz(ifz:ilz,ifx:ilx,1) = f%uz(ifz:ilz,ifx:ilx,1) + wl*self%buoz(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)

                case ('ux')
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                    f%ux(ifz:ilz,ifx:ilx,1) = f%ux(ifz:ilz,ifx:ilx,1) + wl*self%buox(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)
                
                end select
                
            else
                select case (shot%src%comp)
                case ('uz')
                    f%uz(iz,ix,1) = f%uz(iz,ix,1) + wl*self%buoz(iz,ix)

                case ('ux')
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                    f%ux(iz,ix,1) = f%ux(iz,ix,1) + wl*self%buox(iz,ix)
                    
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
                    case ('uz') !vertical z adjsource
                        f%uz(ifz:ilz,ifx:ilx,1) = f%uz(ifz:ilz,ifx:ilx,1) + wl*self%buoz(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1)

                    case ('ux') !horizontal x adjsource
                        f%ux(ifz:ilz,ifx:ilx,1) = f%ux(ifz:ilz,ifx:ilx,1) + wl*self%buox(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1)
                        
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                    case ('uz')
                        f%uz(iz,ix,1) = f%uz(iz,ix,1) + wl*self%buoz(iz,ix)

                    case ('ux')
                        if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                        f%ux(iz,ix,1) = f%ux(iz,ix,1) + wl*self%buox(iz,ix)

                    end select
                    
                endif
                
            enddo
        
    end subroutine

    subroutine update_displacement(self,f,time_dir,it)
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

            call fd2d_laplacian(f%ez,f%ex,f%es,&
                                f%dez_dz,f%dex_dx,f%dex_dz,f%dez_dx,f%des_dz,f%des_dx,&
                                f%lapz,f%lapx,&
                                self%ldap2mu,self%lda,self%mu,&
                                ifz,ilz,ifx,ilx)

        endif

        if(f%is_adjoint) then !flip
            f%lapz=-f%lapz
            f%lapx=-f%lapx
        endif

        if(time_dir>0.) then !in forward time
            f%uz_next(ifz:ilz,ifx:ilx,1) = 2*f%uz(ifz:ilz,ifx:ilx,1) -f%uz_prev(ifz:ilz,ifx:ilx,1) +dt2*self%buoz(ifz:ilz,ifx:ilx)*f%lapz(ifz:ilz,ifx:ilx,1)
            f%ux_next(ifz:ilz,ifx:ilx,1) = 2*f%ux(ifz:ilz,ifx:ilx,1) -f%ux_prev(ifz:ilz,ifx:ilx,1) +dt2*self%buox(ifz:ilz,ifx:ilx)*f%lapx(ifz:ilz,ifx:ilx,1)
        else !in reverse time
            f%uz_prev(ifz:ilz,ifx:ilx,1) = 2*f%uz(ifz:ilz,ifx:ilx,1) -f%uz_next(ifz:ilz,ifx:ilx,1) +dt2*self%buoz(ifz:ilz,ifx:ilx)*f%lapz(ifz:ilz,ifx:ilx,1)
            f%ux_prev(ifz:ilz,ifx:ilx,1) = 2*f%ux(ifz:ilz,ifx:ilx,1) -f%ux_next(ifz:ilz,ifx:ilx,1) +dt2*self%buox(ifz:ilz,ifx:ilx)*f%lapx(ifz:ilz,ifx:ilx,1)
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


    !forward: add RHS to e^it+0.5
    !adjoint: add RHS to e^it+1.5
    subroutine inject_strains(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        if(.not. f%is_adjoint) then

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1

            wl=time_dir*f%wavelet(1,it)/m%cell_volume !as no time derivative on strains
            
            if(if_hicks) then
                select case (shot%src%comp)
                case ('ez')
                    if(m%is_freesurface.and.shot%src%iz==1) then !on FS
                        f%ez(ifz:ilz,ifx:ilx,1) = f%ez(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*( self%lda(ifz:ilz,ifx:ilx)**2/self%ldap2mu(ifz:ilz,ifx:ilx))*shot%src%interp_coef_symm(:,:,1)
                        f%ex(ifz:ilz,ifx:ilx,1) = f%ex(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*(-self%lda(ifz:ilz,ifx:ilx)) *shot%src%interp_coef_symm(:,:,1)

                    else !no interfaction with FS
                        f%ez(ifz:ilz,ifx:ilx,1) = f%ez(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*  self%ldap2mu(ifz:ilz,ifx:ilx) *shot%src%interp_coef_full(:,:,1)
                        f%ex(ifz:ilz,ifx:ilx,1) = f%ex(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*(-self%lda    (ifz:ilz,ifx:ilx))*shot%src%interp_coef_full(:,:,1)

                    endif

                case ('ex')
                    f%ez(ifz:ilz,ifx:ilx,1) = f%ez(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*(-self%lda    (ifz:ilz,ifx:ilx))*shot%src%interp_coef_symm(:,:,1)
                    f%ex(ifz:ilz,ifx:ilx,1) = f%ex(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*  self%ldap2mu(ifz:ilz,ifx:ilx) *shot%src%interp_coef_symm(:,:,1)

                case ('es')
                    f%es(ifz:ilz,ifx:ilx,1) = f%es(ifz:ilz,ifx:ilx,1) +   wl/self%mu(ifz:ilz,ifx:ilx) *shot%src%interp_coef(:,:,1)

                endselect
                
            else
                select case (shot%src%comp)
                case ('ez')
                    !if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                    f%ez(iz,ix,1) = f%ez(iz,ix,1) + wl*self%inv_ldapmu_4mu(iz,ix)*  self%ldap2mu(iz,ix)
                    f%ex(iz,ix,1) = f%ex(iz,ix,1) + wl*self%inv_ldapmu_4mu(iz,ix)*(-self%lda    (iz,ix))

                case ('ex')
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                    f%ez(iz,ix,1) = f%ez(iz,ix,1) + wl*self%inv_ldapmu_4mu(iz,ix)*(-self%lda    (iz,ix))
                    f%ex(iz,ix,1) = f%ex(iz,ix,1) + wl*self%inv_ldapmu_4mu(iz,ix)*  self%ldap2mu(iz,ix)

                case ('es')
                    f%es(iz,ix,1) = f%es(iz,ix,1) +   wl/self%mu(iz,ix)

                endselect

            endif

            return

        endif

            do i=1,shot%nrcv
                ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                
                wl=f%wavelet(i,it)/m%cell_volume !as no time derivative on strains
                    
                if(if_hicks) then
                    select case (shot%rcv(i)%comp)
                    case ('ez')
                        if(m%is_freesurface.and.shot%rcv(i)%iz==1) then !on FS
                            f%ez(ifz:ilz,ifx:ilx,1) = f%ez(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*( self%lda(ifz:ilz,ifx:ilx)**2/self%ldap2mu(ifz:ilz,ifx:ilx))*shot%rcv(i)%interp_coef_symm(:,:,1)
                            f%ex(ifz:ilz,ifx:ilx,1) = f%ex(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*(-self%lda(ifz:ilz,ifx:ilx)) *shot%rcv(i)%interp_coef_symm(:,:,1)

                        else !no interfaction with FS
                            f%ez(ifz:ilz,ifx:ilx,1) = f%ez(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*  self%ldap2mu(ifz:ilz,ifx:ilx) *shot%rcv(i)%interp_coef_full(:,:,1)
                            f%ex(ifz:ilz,ifx:ilx,1) = f%ex(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*(-self%lda    (ifz:ilz,ifx:ilx))*shot%rcv(i)%interp_coef_full(:,:,1)

                        endif


                    case ('ex')
                        f%ez(ifz:ilz,ifx:ilx,1) = f%ez(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*(-self%lda    (ifz:ilz,ifx:ilx))*shot%rcv(i)%interp_coef_symm(:,:,1)
                        f%ex(ifz:ilz,ifx:ilx,1) = f%ex(ifz:ilz,ifx:ilx,1) + wl*self%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)*  self%ldap2mu(ifz:ilz,ifx:ilx) *shot%rcv(i)%interp_coef_symm(:,:,1)

                    case ('es')
                        f%es(ifz:ilz,ifx:ilx,1) = f%es(ifz:ilz,ifx:ilx,1) + wl/self%mu(ifz:ilz,ifx:ilx) *shot%rcv(i)%interp_coef(:,:,1)

                    endselect

                else           
                    select case (shot%rcv(i)%comp)
                    case ('ez')
                        !if(m%is_freesurface.and.shot%rcv(i)%iz==1) wl=2*wl
                        f%ez(iz,ix,1) = f%ez(iz,ix,1) + wl*self%inv_ldapmu_4mu(iz,ix)*  self%ldap2mu(iz,ix)
                        f%ex(iz,ix,1) = f%ex(iz,ix,1) + wl*self%inv_ldapmu_4mu(iz,ix)*(-self%lda    (iz,ix))

                    case ('ex')
                        if(m%is_freesurface.and.shot%rcv(i)%iz==1) wl=2*wl
                        f%ez(iz,ix,1) = f%ez(iz,ix,1) + wl*self%inv_ldapmu_4mu(iz,ix)*(-self%lda    (iz,ix))
                        f%ex(iz,ix,1) = f%ex(iz,ix,1) + wl*self%inv_ldapmu_4mu(iz,ix)*  self%ldap2mu(iz,ix)

                    case ('es')
                        f%es(iz,ix,1) = f%es(iz,ix,1) +wl/self%mu(iz,ix)

                    endselect

                endif

            enddo
        
    end subroutine


    subroutine update_strains(self,f,time_dir,it)
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

            call fd2d_strains(f%uz,f%ux,f%ez,f%ex,f%es,&
                              f%duz_dz,f%dux_dx,f%duz_dx,f%dux_dz,&
                              ifz,ilz,ifx,ilx)

        endif

        if(f%is_adjoint) then !flip
            f%ez=-f%ez
            f%ex=-f%ex
            f%es=-f%es
        endif

    end subroutine

    subroutine evolve(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        real,dimension(:,:,:),pointer :: tmp

        if(time_dir>0.) then !in forward time
                  tmp=>f%uz_prev
            f%uz_prev=>f%uz
            f%uz     =>f%uz_next
            f%uz_next=>tmp

                  tmp=>f%ux_prev
            f%ux_prev=>f%ux
            f%ux     =>f%ux_next
            f%ux_next=>tmp

        else !in reverse time
                  tmp=>f%uz_next
            f%uz_next=>f%uz
            f%uz     =>f%uz_prev
            f%uz_prev=>tmp

                  tmp=>f%ux_next
            f%ux_next=>f%ux
            f%ux     =>f%ux_prev
            f%ux_prev=>tmp
            
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
                
                if(if_hicks) then
                    select case (shot%rcv(i)%comp)
                    case ('uz')
                        f%seismo(i,it)=sum(f%uz(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                    case ('ux')
                        f%seismo(i,it)=sum(f%ux(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )


                    case ('ez')
                        f%seismo(i,it)=sum(f%ez(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef_full(:,:,1) )
                    case ('ex')
                        f%seismo(i,it)=sum(f%ex(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef_full(:,:,1) )

                    case ('es')
                        f%seismo(i,it)=sum(f%es(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )

                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                    case ('uz') !pz[iz-0.5,ix]
                        f%seismo(i,it)=f%uz(iz,ix,1)
                    case ('ux') !px[iz,ix-0.5]
                        f%seismo(i,it)=f%ux(iz,ix,1)

                    case ('ez') !ez[iz,ix]
                        f%seismo(i,it)=f%ez(iz,ix,1)
                    case ('ex') !ez[iz,ix]
                        f%seismo(i,it)=f%ex(iz,ix,1)

                    case ('es') !ez[iz,ix]
                        f%seismo(i,it)=f%es(iz,ix,1)

                    end select
                    
                endif

            enddo

            return

        endif

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            if(if_hicks) then
                select case (shot%src%comp)
                case ('uz')
                    f%seismo(1,it)=sum(f%uz(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                case ('ux')
                    f%seismo(1,it)=sum(f%ux(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))

                case ('ez')
                    f%seismo(1,it)=sum(f%ez(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef_full(:,:,1))
                case ('ex')
                    f%seismo(1,it)=sum(f%ex(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef_full(:,:,1))

                case ('es')
                    f%seismo(1,it)=sum(f%es(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                
                end select
                
            else
                select case (shot%src%comp)
                case ('uz') !pz[iz-0.5,ix,1]
                    f%seismo(1,it)=f%uz(iz,ix,1)
                case ('ux') !px[iz,ix-0.5,1]
                    f%seismo(1,it)=f%ux(iz,ix,1)

                case ('ez') !ez[iz-0.5,ix,1]
                    f%seismo(1,it)=f%ez(iz,ix,1)
                case ('ex') !ex[iz-0.5,ix,1]
                    f%seismo(1,it)=f%ex(iz,ix,1)

                case ('es') !ex[iz-0.5,ix,1]
                    f%seismo(1,it)=f%es(iz,ix,1)
                
                end select
                
            endif
        
    end subroutine
        
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%buoz, self%buox, self%ldap2mu, self%lda, self%mu)
    end subroutine


    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !
    !PDE:      A u = CP u - CDC u - f
    !Adjoint:  Aᵀa = CP a + CDC a - d
    !
    !<a|Au> = <a|CPu-CDCu-f> = <a|CPu> - <a|CDCu> - <a|f>
    !K_ρ<a|CPu> = <a|C(K_ρP)u>
    !Kₘ<a|Au> = <a|(KₘC)Pu> - <a|(KₘC)DCu> - <a|C(KₘDC)u>
    !         = <a|(KₘC)C⁻¹f>              - <a|C(KₘDC)u>
    !         ≐-<a|C(KₘDC)u> = <DCa|(KₘC)u>
    !         ≐-<Pa|(KₘC)u> =-∫ (Pa)ᵀ KₘC u dt
    !
    !Pa = [ρ∂ₜₜuzᵃ ρ∂ₜₜuxᵃ ezᵃ exᵃ esᵃ]ᵀ
    ! u = [    uz      ux  ez  ex  es ]ᵀ
    !
    !     [0     ]
    !     | 0    |
    !K_λC=|  1 1 |, K_μC=diag{0,0,2,2,1}
    !     |  1 1 |
    !     [     0]
    !
    !Therefore,
    !glda =  (ezᵃ+exᵃ)*(ez+ex)
    !gmu  = 2(ezᵃ*ez+exᵃ*ex) + esᵃ*es 

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

    subroutine cross_correlate_glda_gmu(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ! ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ! ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)

        if(m%is_cubic) then
        else
            call grad2d_glda_gmu(rf%ez(:,:,1),rf%ex(:,:,1),rf%es(:,:,1),&
                                 sf%ez(:,:,1),sf%ex(:,:,1),sf%es(:,:,1),&
                                 corr%glda,corr%gmu, &
                                 ifz,ilz,ifx,ilx)
        endif

    end subroutine

    subroutine cross_correlate_postprocess(corr)
        type(t_correlate) :: corr

        ! if(allocated(correlate_gradient)) then

            !preparing for projection back
            corr%glda(1,:,:) = corr%glda(2,:,:)
            corr%gmu (1,:,:) = corr%gmu (2,:,:)

            call interp2D(corr%gmu(:,:,1),[1,1])
            ! corr%gmu(1,:,:) = corr%gmu(2,:,:)

        ! endif

    end subroutine

    subroutine interp2D(array,ishifts)
        real,dimension(:,:) :: array
        integer :: ishifts(2)

        real,dimension(:,:),allocatable :: tmp

        n1 = size(array,dim=1); ish1 = ishifts(1)
        n2 = size(array,dim=2); ish2 = ishifts(2)

        allocate(tmp(n1,n2))

        do i2=1,n2-ish2 !may cause out-of-bounds if ish2<0
        do i1=1,n1-ish1 !may cause out-of-bounds if ish1<0
            tmp(i1,i2)=( array(i1     ,i2     ) &
                        +array(i1+ish1,i2     ) &
                        +array(i1     ,i2+ish2) &
                        +array(i1+ish1,i2+ish2) )/4
        enddo
        enddo

        tmp(:,n2)=tmp(:,n2-ish2)
        tmp(n1,:)=tmp(n1-ish1,:)

        array=tmp

        deallocate(tmp)

    end subroutine

    ! subroutine roll2D(array,ishifts)
    !     real,dimension(:,:) :: array
    !     integer :: ishifts(2)

    !     real,dimension(:,:),allocatable :: tmp

    !     n1 = size(array,dim=1)
    !     n2 = size(array,dim=2)

    !     allocate(tmp(n1,n2))

    !     tmp(1:n1,1:n2)=array(1:n1,1:n2)

    !     array=tmp

    !     deallocate(tmp)

    ! end subroutine

    ! subroutine cross_correlate_image(rf,sf,corr,it)
    !     type(t_field), intent(in) :: rf, sf
    !     type(t_correlate) :: corr

    !     !nonzero only when sf touches rf
    !     ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
    !     ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
    !     ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
    !     ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
    !     ! ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
    !     ! ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)

    !     call imag2d(rf%p,sf%p,&
    !                 rf%poynz,rf%poynx,sf%poynz,sf%poynx, &
    !                 corr%ipp,corr%ibksc,corr%ifwsc, &
    !                 ifz,ilz,ifx,ilx)

    ! end subroutine

    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd2d_laplacian(ez,ex,es,&
                              dez_dz,dex_dx,dex_dz,dez_dx,des_dz,des_dx,&
                              lapz,lapx,&
                              ldap2mu,lda,mu,&
                              ifz,ilz,ifx,ilx)
        real,dimension(*) :: ez,ex,es
        real,dimension(*) :: dez_dz,dex_dx,dex_dz,dez_dx,des_dz,des_dx
        real,dimension(*) :: lapz,lapx
        real,dimension(*) :: ldap2mu,lda,mu

        nz=cb%nz
        nx=cb%nx

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dez_dz_,dex_dx_,dex_dz_,dez_dx_,des_dz_,des_dx_)
        !$omp do schedule(dynamic)
        do ix=ifx+2,ilx-2

            !dir$ simd
            do iz=ifz+2,ilz-2

                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1

                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2=i  -2*nz !iz,ix-2
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1
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
                dex_dx_=dex_dx_*cpml%kpa_x_half(ix) + dex_dx(i)
                dex_dz_=dex_dz_*cpml%kpa_z_half(iz) + dex_dz(i)
                dez_dx_=dez_dx_*cpml%kpa_x_half(ix) + dez_dx(i)
                des_dz_=des_dz_*cpml%kpa_z(iz)      + des_dz(i)
                des_dx_=des_dx_*cpml%kpa_x(ix)      + des_dx(i)
                
                !displacement
                lapz(i)= dez_dz_+dex_dz_+des_dx_
                lapx(i)= dez_dx_+dex_dx_+des_dz_

            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine fd2d_strains(uz,ux,ez,ex,es,&
                            duz_dz,dux_dx,duz_dx,dux_dz,&
                            ifz,ilz,ifx,ilx)
        real,dimension(*) :: uz,ux,ez,ex,es
        real,dimension(*) :: duz_dz,dux_dx,duz_dx,dux_dz
        
        nz=cb%nz
        nx=cb%nx
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         duz_dz_,dux_dx_,duz_dx_,dux_dz_)
        !$omp do schedule(dynamic)
        do ix=ifx+2,ilx-2
        
            !dir$ simd
            do iz=ifz+2,ilz-2
            
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
                

                duz_dz_= c1z*(uz(izp1_ix)-uz(iz_ix))  +c2z*(uz(izp2_ix)-uz(izm1_ix))
                dux_dx_= c1x*(ux(iz_ixp1)-ux(iz_ix))  +c2x*(ux(iz_ixp2)-ux(iz_ixm1))
                
                !cpml
                duz_dz(i)=cpml%b_z(iz)*duz_dz(i)+cpml%a_z(iz)*duz_dz_
                dux_dx(i)=cpml%b_x(ix)*dux_dx(i)+cpml%a_x(ix)*dux_dx_

                duz_dz_=duz_dz_*cpml%kpa_z(iz) + duz_dz(iz_ix)
                dux_dx_=dux_dx_*cpml%kpa_x(ix) + dux_dx(iz_ix)
                
                !normal strains
                ez(i) = duz_dz_
                ex(i) = dux_dx_

                duz_dx_= c1x*(uz(iz_ix)-uz(iz_ixm1))  +c2x*(uz(iz_ixp1)-uz(iz_ixm2))
                dux_dz_= c1z*(ux(iz_ix)-ux(izm1_ix))  +c2z*(ux(izp1_ix)-ux(izm2_ix))

                !cpml
                duz_dx(i)=cpml%b_x_half(ix)*duz_dx(i)+cpml%a_x_half(ix)*duz_dx_
                dux_dz(i)=cpml%b_z_half(iz)*dux_dz(i)+cpml%a_z_half(iz)*dux_dz_

                duz_dx_=duz_dx_*cpml%kpa_x_half(ix) + duz_dx(i)
                dux_dz_=dux_dz_*cpml%kpa_z_half(iz) + dux_dz(i)

                !shear stress
                es(i) = duz_dx_+dux_dz_
                
            enddo
            
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine

    subroutine grad2d_glda_gmu(rf_ez,rf_ex,rf_es,&
                           sf_ez,sf_ex,sf_es,&
                           glda,gmu,&
                           ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_ez,rf_ex,rf_es
        real,dimension(*) :: sf_ez,sf_ex,sf_es
        real,dimension(*) :: glda,gmu
        
        nz=cb%nz
        nx=cb%nx
        
        !glda = -(ezᵃ+exᵃ)*(ez+ex)
        !gmu  =-2(ezᵃ*ez+exᵃ*ex) + esᵃ*es 
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j)
        !$omp do schedule(dynamic)
        do ix = ifx,ilx
            !dir$ simd
            do iz = ifz,ilz

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers
                
                glda(j) = glda(j) -  (rf_ez(i)+rf_ex(i))*(sf_ez(i)+sf_ex(i))

                gmu (j) = gmu (j) -2*(rf_ez(i)*sf_ez(i) +rf_ex(i)*sf_ex(i)) &
                                    - rf_es(i)*sf_es(i)

            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

end
