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
            '2nd-order Displacement formulation'//s_NL// &
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
        integer :: ngrad=66 !number of basic gradients
        ! integer :: nimag=3 !number of basic images
        !integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:),allocatable :: buoz, buox!, buoy
        real,dimension(:,:),allocatable :: ldap2mu, lda, mu

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
        procedure :: update
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
            call alloc(m%rho,m%nz,m%nx,1,o_init=1000.)
            call warn('Constant rho model (1000 kg/m³) is allocated by propagator.')
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
        call alloc(self%mu,     [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        call alloc(temp_mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])

        self%ldap2mu(:,:)=cb%rho(:,:,1)*cb%vp(:,:,1)**2
             temp_mu(:,:)=cb%rho(:,:,1)*cb%vs(:,:,1)**2

        self%lda=self%ldap2mu-2.*temp_mu
        ! if(mpiworld%is_master) then
        ! write(*,*) 'self%ldap2mu sanity:', minval(self%ldap2mu),maxval(self%ldap2mu)
        ! write(*,*) 'self%lda     sanity:', minval(self%lda),maxval(self%lda)
        ! endif

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

        ! call f%init_bloom
        call warn("Disable field%bloom, otherwise back propagation of incident field is unstable (don't know why..)")
        call alloc(f%bloom,6,self%nt)
        f%bloom(1,:)=cb%ifz
        f%bloom(2,:)=cb%ilz
        f%bloom(3,:)=cb%ifx
        f%bloom(4,:)=cb%ilx
        f%bloom(5,:)=cb%ify
        f%bloom(6,:)=cb%ily


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
        
        call alloc(f%dz_ldap2mu_duzdz_p_lda_duxdx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dx_lda_duzdz_p_ldap2mu_duxdx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dz_mu_duxdz_p_duzdx         ,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%dx_mu_duxdz_p_duzdx         ,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        call alloc(f%lapz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%lapx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

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

call alloc(corr%g11,m%nz,m%nx,m%ny)
call alloc(corr%g12,m%nz,m%nx,m%ny)
call alloc(corr%g16,m%nz,m%nx,m%ny)
call alloc(corr%g21,m%nz,m%nx,m%ny)
call alloc(corr%g22,m%nz,m%nx,m%ny)
call alloc(corr%g26,m%nz,m%nx,m%ny)
call alloc(corr%g61,m%nz,m%nx,m%ny)
call alloc(corr%g62,m%nz,m%nx,m%ny)
call alloc(corr%g66,m%nz,m%nx,m%ny)

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

call correlate_assemble(corr%g11, correlate_gradient(:,:,:,11))
call correlate_assemble(corr%g12, correlate_gradient(:,:,:,12))
call correlate_assemble(corr%g16, correlate_gradient(:,:,:,16))
call correlate_assemble(corr%g21, correlate_gradient(:,:,:,21))
call correlate_assemble(corr%g22, correlate_gradient(:,:,:,22))
call correlate_assemble(corr%g26, correlate_gradient(:,:,:,26))
call correlate_assemble(corr%g61, correlate_gradient(:,:,:,61))
call correlate_assemble(corr%g62, correlate_gradient(:,:,:,62))
call correlate_assemble(corr%g66, correlate_gradient(:,:,:,66))

        endif        
        
    end subroutine

    
    !========= Derivations =================
    !PDE:      A u = ρ∂ₜₜ u - DᵀCD u = f
    !Adjoint:  Aᵀa = ρ∂ₜₜᵀa - DᵀCD a = d
    !where
    !u=[uz ux]ᵀ is displacement, f=[fz fx]ᵀδ(x-xs) with xs source position, d is recorded data
    !  [λ+2μ  λ      ]
    !C=| λ   λ+2μ    |, Dᵀ=[∂z 0  0  ∂ₓ]
    !  |          μ μ|     [0  ∂ₓ ∂z 0 ]
    !  [          μ μ]
    !a=[uzᵃ uxᵃ]ᵀ is the adjoint field
    !
    !Discrete case:
    !Meshing with staggered grids in time and space (2D example):
    !                   |    ss    |   ss -½ uz  ss    |         |
    !                   |    μ     |   μ     bz  μ     |         |
    !                   |          |         |         |         |
    !                  λ,μ   bx   λ,μ  bx   λ,μ  bx   λ,μ  bx   λ,μ
    !  -u--u--u-→ t    -sn---ux---sn---ux---sn---ux---sn---ux---sn-→ x
    !  -1  0  1        -2   -1½   -1   -½    0    ½    1   1½    2    
    !                   |          |         |         |         | 
    !                   |    ss    |   ss  ½ uz  ss    |         | 
    !                   |    μ     |   μ     bz  μ     |         | 
    !                   |          |         |         |         | 
    !                  -|----------|-------1-sn--------|---------|-
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
    !  sn(iz,ix)    =>   sn[iz,  ix  ]^n+½ :=sn(iz*dz,    ix*dx,(n+½)*dt)
    !  ss(iz,ix)    =>   ss[iz-½,ix-½]^n  
    !
    !Forward:
    !FD eqn:
    !                                [∂zᶠ  0 ]
    !ρ∂ₜ² [uz^n] = [∂zᵇ 0   0 ∂ₓᶠ] C | 0  ∂ₓᶠ|[uz^n]
    !     [ux^n]   [0  ∂ₓᵇ ∂zᶠ 0 ]   | 0  ∂zᵇ|[ux^n]
    !                                [∂ₓᵇ  0 ]
    !where
    !∂ₜᶠ*dt := v^n+1 - v^n                             ~O(t²)
    !∂zᵇ*dz := c₁(s(iz  )-s(iz-1) +c₂(s(iz+1)-s(iz-2)  ~O(x⁴)
    !∂zᶠ*dz := c₁(v(iz+1)-v(iz  ) +c₂(v(iz+2)-v(iz-1)  ~O(x⁴)
    !Step #1: u^n  += src
    !Step #2: save u^n to boundary values
    !Step #3: u^n+1 = 2u^n -u^n-1 +laplacian of u^n
    !Step #4: (u^n-1,u^n) = (u^n,u^n+1)
    !Step #5: sample u^n at receivers
    !
    !in reverse time:
    !Step #4: (u^n,u^n+1) = (u^n-1,u^n)
    !Step #2: load boundary values for u^n
    !Step #3: u^n-1 = 2u^n -u^n+1 +laplacian of u^n
    !Step #1: u^n -= src
    !
    !Adjoint:
    !since
    !∂ₜ²ᵀ = ∂ₜ²
    !∂zᵇᵀ = u(iz)-u(iz+1) = -∂zᶠ, ∂zᶠᵀ = -∂zᵇ
    !FD eqn:
    !                                [∂zᶠ  0 ]
    !ρ∂ₜ² [uz^n] = [∂zᵇ 0   0 ∂ₓᶠ] C | 0  ∂ₓᶠ|[uz^n]
    !     [ux^n]   [0  ∂ₓᵇ ∂zᶠ 0 ]   | 0  ∂zᵇ|[ux^n]
    !                                [∂ₓᵇ  0 ]
    !SAME as the discretized FD eqn!
    !
    !Time marching (in reverse time):
    !Step #1: uᵃ^n += adjsrc
    !Step #2: uᵃ^n-1 = 2uᵃ^n -uᵃ^n+1 +laplacian of uᵃ^n
    !Step #3: cross correlate
    !Step #4: (uᵃ^n,uᵃ^n+1) = (uᵃ^n-1,uᵃ^n)
    !Step #5: sample uᵃ^n at source        

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

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add force
            call cpu_time(tic)
            call self%inject_displacement(fld_u,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2: save p^it+1 in boundary layers
            call cpu_time(tic)
            call fld_u%boundary_transport_displacement('save',it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            ! !step 3: set hardBC
            ! call cpu_time(tic)
            ! call self%set_pressure(fld_E0,time_dir,it)
            ! call cpu_time(toc)
            ! tt3=tt3+toc-tic

            !step 3: update
            call cpu_time(tic)
            call self%update(fld_u,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !step 4: evolve, it -> it+1
            call cpu_time(tic)
            call self%evolve(fld_u,time_dir,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !step 5: sample p^it+1 at receivers
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

            ! if(present(o_sf)) then
                !backward step 4: it+1 -> it
                call cpu_time(tic)
                call self%evolve(fld_u,time_dir,it)
                call cpu_time(toc)
                tt1=tt1+toc-tic

                ! !backward step 2: retrieve p^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_u%boundary_transport_displacement('load',it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3:
                call cpu_time(tic)
                call self%update(fld_u,time_dir,it)
                call cpu_time(toc)
                tt4=tt4+toc-tic

                !backward step 1: rm p^it at source
                call cpu_time(tic)
                call self%inject_displacement(fld_u,time_dir,it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            ! endif

            !adjoint step 1: inject to p^it+1 at receivers
            call cpu_time(tic)
            call self%inject_displacement(fld_a,time_dir,it)
            call cpu_time(toc)
            tt8=tt8+toc-tic

            !adjoint step 2:
            call cpu_time(tic)
            call self%update(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            ! !image: rf%p^it star sf%p^it
            ! if(mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call cross_correlate_image(fld_p,fld_u,a_star_u,it)
            !     call cpu_time(toc)
            !     tt10=tt10+toc-tic
            ! endif

            !adjoint step 3 gradient: rf%p^it star sf%p^it
            if(mod(it,irdt)==0) then
                call cpu_time(tic)
                call cross_correlate(fld_a,fld_u,a_star_u,it)
                call cross_correlate_gij(fld_a,fld_u,a_star_u,it)
                call cpu_time(toc)
                tt10=tt10+toc-tic
            endif

            !adjoint step 4
            ! this step is moved to update_pressure for easier management
            call cpu_time(tic)
            call self%evolve(fld_a,time_dir,it)
            call cpu_time(toc)
            tt11=tt11+toc-tic

            !adjoint step 5: sample p^it at source position
            if(if_propagator_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_a,it)
                call cpu_time(toc)
                tt12=tt12+toc-tic
            endif


            !--------------------------------------------------------!
            
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
            ! write(*,*) 'Elapsed time to compute Poynting vectors ',tt3/mpiworld%max_threads
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

    subroutine update(self,f,time_dir,it)
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

            call fd2d_laplacian(f%uz,f%ux,&
                                f%duz_dz,f%dux_dx,f%dux_dz,f%duz_dx,&
                                f%dz_ldap2mu_duzdz_p_lda_duxdx,&
                                f%dx_lda_duzdz_p_ldap2mu_duxdx,&
                                f%dz_mu_duxdz_p_duzdx,&
                                f%dx_mu_duxdz_p_duzdx,&
                                f%lapz,f%lapx,&
                                self%ldap2mu,self%lda,self%mu,&
                                ifz,ilz,ifx,ilx)

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

                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                    case ('uz') !pz[iz-0.5,ix]
                        f%seismo(i,it)=f%uz(iz,ix,1)
                    case ('ux') !px[iz,ix-0.5]
                        f%seismo(i,it)=f%ux(iz,ix,1)

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

                end select
                
            else
                select case (shot%src%comp)
                case ('uz') !pz[iz-0.5,ix,1]
                    f%seismo(1,it)=f%uz(iz,ix,1)
                case ('ux') !px[iz,ix-0.5,1]
                    f%seismo(1,it)=f%ux(iz,ix,1)

                end select
                
            endif
        
    end subroutine
        
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%buoz, self%buox, self%ldap2mu, self%lda, self%mu)
    end subroutine


    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !<a|Au> = <a|ρ∂ₜₜu-DᵀCDu> = <a|ρ∂ₜₜu> - <a|DᵀCDu>
    !   K_ρ <a|ρ∂ₜₜu> = <a|∂ₜₜu> ≐ <a|bDᵀCDu>
    !   Kₘ  <a|DᵀCDu> =-<Da|(KₘC)Du>
    !
    !  [λ+2μ  λ      ]         [1 1    ]         [2 0    ]
    !C=| λ   λ+2μ    |, K_λC = |1 1    |, K_μC = |0 2    |
    !  |          μ μ|         |    0 0|         |    1 1|
    !  [          μ μ]         [    0 0]         [    1 1]
    !
    !Therefore,
    !glda = ...
    !gmu  = ...

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

        if(m%is_cubic) then
        else
            call grad2d_glda_gmu(rf%uz(:,:,1),rf%ux(:,:,1),&
                                 sf%uz(:,:,1),sf%ux(:,:,1),&
                                 corr%glda,corr%gmu, &
                                 ifz,ilz,ifx,ilx)
        endif

    end subroutine

    subroutine cross_correlate_gij(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        real,dimension(:,:,:),allocatable :: sf_duz_dz, sf_dux_dx, sf_duz_dx, sf_dux_dz
        real,dimension(:,:,:),allocatable :: rf_duz_dz, rf_dux_dx, rf_duz_dx, rf_dux_dz

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2) !
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        
        call alloc(rf_duz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(rf_dux_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(rf_duz_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(rf_dux_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        
        call alloc(sf_duz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(sf_dux_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(sf_duz_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(sf_dux_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])

        rf_duz_dz(ifz:ilz,ifx:ilx,1) = (rf%uz(ifz+1:ilz+1,ifx:ilx,1)-rf%uz(ifz:ilz,ifx:ilx,1))/m%dz
        rf_dux_dx(ifz:ilz,ifx:ilx,1) = (rf%ux(ifz:ilz,ifx+1:ilx+1,1)-rf%ux(ifz:ilz,ifx:ilx,1))/m%dx
        rf_duz_dx(ifz:ilz,ifx:ilx,1) = (rf%uz(ifz:ilz,ifx:ilx,1)-rf%uz(ifz:ilz,ifx-1:ilx-1,1))/m%dx
        rf_dux_dz(ifz:ilz,ifx:ilx,1) = (rf%uz(ifz:ilz,ifx:ilx,1)-rf%uz(ifz-1:ilz-1,ifx:ilx,1))/m%dz

        sf_duz_dz(ifz:ilz,ifx:ilx,1) = (sf%uz(ifz+1:ilz+1,ifx:ilx,1)-sf%uz(ifz:ilz,ifx:ilx,1))/m%dz
        sf_dux_dx(ifz:ilz,ifx:ilx,1) = (sf%ux(ifz:ilz,ifx+1:ilx+1,1)-sf%ux(ifz:ilz,ifx:ilx,1))/m%dx
        sf_duz_dx(ifz:ilz,ifx:ilx,1) = (sf%uz(ifz:ilz,ifx:ilx,1)-sf%uz(ifz:ilz,ifx-1:ilx-1,1))/m%dx
        sf_dux_dz(ifz:ilz,ifx:ilx,1) = (sf%uz(ifz:ilz,ifx:ilx,1)-sf%uz(ifz-1:ilz-1,ifx:ilx,1))/m%dz

        corr%g11(ifz:ilz,ifx:ilx,1) = corr%g11(ifz:ilz,ifx:ilx,1) + rf_duz_dz(ifz:ilz,ifx:ilx,1) *  sf_duz_dz(ifz:ilz,ifx:ilx,1)
        corr%g12(ifz:ilz,ifx:ilx,1) = corr%g12(ifz:ilz,ifx:ilx,1) + rf_duz_dz(ifz:ilz,ifx:ilx,1) *  sf_dux_dx(ifz:ilz,ifx:ilx,1)
        corr%g16(ifz:ilz,ifx:ilx,1) = corr%g16(ifz:ilz,ifx:ilx,1) + rf_duz_dz(ifz:ilz,ifx:ilx,1) * (sf_duz_dx(ifz:ilz,ifx:ilx,1)+sf_dux_dz(ifz:ilz,ifx:ilx,1))
        
        corr%g21(ifz:ilz,ifx:ilx,1) = corr%g21(ifz:ilz,ifx:ilx,1) + rf_dux_dx(ifz:ilz,ifx:ilx,1) *  sf_duz_dz(ifz:ilz,ifx:ilx,1)
        corr%g22(ifz:ilz,ifx:ilx,1) = corr%g22(ifz:ilz,ifx:ilx,1) + rf_dux_dx(ifz:ilz,ifx:ilx,1) *  sf_dux_dx(ifz:ilz,ifx:ilx,1)
        corr%g26(ifz:ilz,ifx:ilx,1) = corr%g26(ifz:ilz,ifx:ilx,1) + rf_dux_dx(ifz:ilz,ifx:ilx,1) * (sf_duz_dx(ifz:ilz,ifx:ilx,1)+sf_dux_dz(ifz:ilz,ifx:ilx,1))
        
        corr%g61(ifz:ilz,ifx:ilx,1) = corr%g61(ifz:ilz,ifx:ilx,1) + (rf_duz_dx(ifz:ilz,ifx:ilx,1)+rf_dux_dz(ifz:ilz,ifx:ilx,1)) &
                                                                  *  sf_duz_dz(ifz:ilz,ifx:ilx,1) 
        corr%g62(ifz:ilz,ifx:ilx,1) = corr%g62(ifz:ilz,ifx:ilx,1) + (rf_duz_dx(ifz:ilz,ifx:ilx,1)+rf_dux_dz(ifz:ilz,ifx:ilx,1)) &
                                                                  *  sf_dux_dx(ifz:ilz,ifx:ilx,1) 
        corr%g66(ifz:ilz,ifx:ilx,1) = corr%g66(ifz:ilz,ifx:ilx,1) + (rf_duz_dx(ifz:ilz,ifx:ilx,1)+rf_dux_dz(ifz:ilz,ifx:ilx,1)) &
                                                                  * (sf_duz_dx(ifz:ilz,ifx:ilx,1)+sf_dux_dz(ifz:ilz,ifx:ilx,1)) 

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
    
    subroutine fd2d_laplacian(uz,ux,&
                              duz_dz,dux_dx,dux_dz,duz_dx,&
                              dz_ldap2mu_duzdz_p_lda_duxdx,&
                              dx_lda_duzdz_p_ldap2mu_duxdx,&
                              dz_mu_duxdz_p_duzdx,&
                              dx_mu_duxdz_p_duzdx,&
                              lapz,lapx,&
                              ldap2mu,lda,mu,&
                              ifz,ilz,ifx,ilx)
        real,dimension(*) :: uz,ux
        real,dimension(*) :: duz_dz,dux_dx,dux_dz,duz_dx
        real,dimension(*) :: dz_ldap2mu_duzdz_p_lda_duxdx
        real,dimension(*) :: dx_lda_duzdz_p_ldap2mu_duxdx
        real,dimension(*) :: dz_mu_duxdz_p_duzdx
        real,dimension(*) :: dx_mu_duxdz_p_duzdx
        real,dimension(*) :: ldap2mu,lda,mu,lapz,lapx

        real,dimension(:),allocatable :: ldap2mu_duzdz_p_lda_duxdx
        real,dimension(:),allocatable :: lda_duzdz_p_ldap2mu_duxdx
        real,dimension(:),allocatable :: mu_duxdz_p_duz_dx

        call alloc(ldap2mu_duzdz_p_lda_duxdx, cb%n)
        call alloc(lda_duzdz_p_ldap2mu_duxdx, cb%n)
        call alloc(mu_duxdz_p_duz_dx,         cb%n)

        nz=cb%nz
        nx=cb%nx

        
        !       [λ+2μ  λ      ][∂zᶠ  0 ]       [λ+2μ  λ      ][∂zᶠuz]   [(λ+2μ)∂zᶠuz +  λ    ∂ₓᶠux]
        !flux = | λ   λ+2μ    || 0  ∂ₓᶠ|[uz] = | λ   λ+2μ    ||∂ₓᶠux| = | λ    ∂zᶠuz + (λ+2μ)∂ₓᶠux|
        !       |          μ μ|| 0  ∂zᵇ|[ux]   |          μ μ||∂zᵇux|   |    μ ∂zᵇux + μ     ∂ₓᵇuz|
        !       [          μ μ][∂ₓᵇ  0 ]       [          μ μ][∂ₓᵇuz]   [    μ ∂zᵇux + μ     ∂ₓᵇuz]
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         duz_dz_,dux_dx_,dux_dz_,duz_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx+2,ilx-2
            !dir$ simd
            do iz = ifz+2,ilz-2

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2=i  -2*nz !iz,ix-2
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1
                iz_ixp2=i  +2*nz  !iz,ix+2

                duz_dz_ = c1z*(uz(izp1_ix)-uz(iz_ix)) +c2z*(uz(izp2_ix)-uz(izm1_ix)) !∂zᶠ
                dux_dx_ = c1x*(ux(iz_ixp1)-ux(iz_ix)) +c2x*(ux(iz_ixp2)-ux(iz_ixm1)) !∂ₓᶠ

                duz_dz(i)=cpml%b_z(iz)*duz_dz(i)+cpml%a_z(iz)*duz_dz_
                dux_dx(i)=cpml%b_x(ix)*dux_dx(i)+cpml%a_x(ix)*dux_dx_

                duz_dz_ = duz_dz_*cpml%kpa_z(iz) + duz_dz(iz_ix)
                dux_dx_ = dux_dx_*cpml%kpa_x(ix) + dux_dx(iz_ix)


                dux_dz_ = c1z*(ux(iz_ix)-ux(izm1_ix)) +c2z*(ux(izp1_ix)-ux(izm2_ix)) !∂zᵇ
                duz_dx_ = c1x*(uz(iz_ix)-uz(iz_ixm1)) +c2x*(uz(iz_ixp1)-uz(iz_ixm2)) !∂ₓᵇ

                dux_dz(i)=cpml%b_z_half(iz)*dux_dz(i)+cpml%a_z_half(iz)*dux_dz_
                duz_dx(i)=cpml%b_x_half(ix)*duz_dx(i)+cpml%a_x_half(ix)*duz_dx_

                dux_dz_ = dux_dz_*cpml%kpa_z_half(iz) + dux_dz(iz_ix)
                duz_dx_ = duz_dx_*cpml%kpa_x_half(ix) + duz_dx(iz_ix)


                ldap2mu_duzdz_p_lda_duxdx(iz_ix) = &
                    ldap2mu(iz_ix)*duz_dz_ +lda    (iz_ix)*dux_dx_
                lda_duzdz_p_ldap2mu_duxdx(iz_ix) = &
                    lda    (iz_ix)*duz_dz_ +ldap2mu(iz_ix)*dux_dx_

                mu_duxdz_p_duz_dx(iz_ix) = mu(iz_ix)*(duz_dx_+dux_dz_)
                
            enddo
        enddo
        !$omp end do
        !$omp end parallel


        !                          [ldap2mu_duzdz_p_lda_duxdx    ]
        !Laplacian= [∂zᵇ 0   0 ∂ₓᶠ]|    lda_duzdz_p_ldap2mu_duxdx|
        !           [0  ∂ₓᵇ ∂zᶠ 0 ]|     mu_duxdz_p_duz_dx       |
        !                          [     mu_duxdz_p_duz_dx       ]
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dz_ldap2mu_duzdz_p_lda_duxdx_,&
        !$omp         dx_lda_duzdz_p_ldap2mu_duxdx_,&
        !$omp         dz_mu_duxdz_p_duzdx_,&
        !$omp         dx_mu_duxdz_p_duzdx_)
        !$omp do schedule(dynamic)
        do ix = ifx+2,ilx-2
            !dir$ simd
            do iz = ifz+2,ilz-2

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2=i  -2*nz !iz,ix-2
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1
                iz_ixp2=i  +2*nz  !iz,ix+2
                
                dz_ldap2mu_duzdz_p_lda_duxdx_ = & !∂zᵇ
                     c1z*(ldap2mu_duzdz_p_lda_duxdx(iz_ix)  -ldap2mu_duzdz_p_lda_duxdx(izm1_ix)) &
                    +c2z*(ldap2mu_duzdz_p_lda_duxdx(izp1_ix)-ldap2mu_duzdz_p_lda_duxdx(izm2_ix))

                dx_lda_duzdz_p_ldap2mu_duxdx_ = & !∂ₓᵇ
                     c1x*(lda_duzdz_p_ldap2mu_duxdx(iz_ix)  -lda_duzdz_p_ldap2mu_duxdx(iz_ixm1)) &
                    +c2x*(lda_duzdz_p_ldap2mu_duxdx(iz_ixp1)-lda_duzdz_p_ldap2mu_duxdx(iz_ixm2))

                dz_ldap2mu_duzdz_p_lda_duxdx(i)=cpml%b_z_half(iz)*dz_ldap2mu_duzdz_p_lda_duxdx(i)+cpml%a_z_half(iz)*dz_ldap2mu_duzdz_p_lda_duxdx_
                dx_lda_duzdz_p_ldap2mu_duxdx(i)=cpml%b_x_half(ix)*dx_lda_duzdz_p_ldap2mu_duxdx(i)+cpml%a_x_half(ix)*dx_lda_duzdz_p_ldap2mu_duxdx_

                dz_ldap2mu_duzdz_p_lda_duxdx_ = dz_ldap2mu_duzdz_p_lda_duxdx_*cpml%kpa_z_half(iz) + dz_ldap2mu_duzdz_p_lda_duxdx(iz_ix)
                dx_lda_duzdz_p_ldap2mu_duxdx_ = dx_lda_duzdz_p_ldap2mu_duxdx_*cpml%kpa_x_half(ix) + dx_lda_duzdz_p_ldap2mu_duxdx(iz_ix)


                dz_mu_duxdz_p_duzdx_ = & !∂zᶠ
                     c1z*(mu_duxdz_p_duz_dx(izp1_ix)-mu_duxdz_p_duz_dx(iz_ix)  ) &
                    +c2z*(mu_duxdz_p_duz_dx(izp2_ix)-mu_duxdz_p_duz_dx(izm1_ix))

                dx_mu_duxdz_p_duzdx_ = & !∂ₓᶠ
                     c1x*(mu_duxdz_p_duz_dx(iz_ixp1)-mu_duxdz_p_duz_dx(iz_ix)  ) &
                    +c2x*(mu_duxdz_p_duz_dx(iz_ixp2)-mu_duxdz_p_duz_dx(iz_ixm1))

                dz_mu_duxdz_p_duzdx(i)=cpml%b_z(iz)*dz_mu_duxdz_p_duzdx(i)+cpml%a_z(iz)*dz_mu_duxdz_p_duzdx_
                dx_mu_duxdz_p_duzdx(i)=cpml%b_x(ix)*dx_mu_duxdz_p_duzdx(i)+cpml%a_x(ix)*dx_mu_duxdz_p_duzdx_
                
                dz_mu_duxdz_p_duzdx_ = dz_mu_duxdz_p_duzdx_*cpml%kpa_z(iz) + dz_mu_duxdz_p_duzdx(iz_ix)
                dx_mu_duxdz_p_duzdx_ = dx_mu_duxdz_p_duzdx_*cpml%kpa_x(ix) + dx_mu_duxdz_p_duzdx(iz_ix)
                

                lapz(iz_ix) = dz_ldap2mu_duzdz_p_lda_duxdx_ + dx_mu_duxdz_p_duzdx_
                lapx(iz_ix) = dx_lda_duzdz_p_ldap2mu_duxdx_ + dz_mu_duxdz_p_duzdx_

            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine


    subroutine grad2d_glda_gmu(rf_uz,rf_ux,&
                               sf_uz,sf_ux,&
                               glda,gmu,&
                               ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_uz,rf_ux, sf_uz, sf_ux
        real,dimension(*) :: glda,gmu
        
        nz=cb%nz
        nx=cb%nx
        
        !     [∂zᶠ  0 ]       [∂zᶠuz]         [1 1    ]         [2 0    ]
        !Du = | 0  ∂ₓᶠ|[uz] = |∂ₓᶠux|, K_λC = |1 1    |, K_μC = |0 2    |
        !     | 0  ∂zᵇ|[ux]   |∂zᵇux|         |    0 0|         |    1 1|
        !     [∂ₓᵇ  0 ]       [∂ₓᵇuz]         [    0 0]         [    1 1]
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         rf_duz_dz,rf_dux_dx,rf_dux_dz,rf_duz_dx,&
        !$omp         sf_duz_dz,sf_dux_dx,sf_dux_dz,sf_duz_dx)
        !$omp do schedule(dynamic)
        do ix = ifx,ilx
            !dir$ simd
            do iz = ifz,ilz

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers
                
                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2=i  -2*nz !iz,ix-2
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1
                iz_ixp2=i  +2*nz !iz,ix+2

                rf_duz_dz = c1z*(rf_uz(izp1_ix)-rf_uz(iz_ix)) +c2z*(rf_uz(izp2_ix)-rf_uz(izm1_ix))
                rf_dux_dx = c1x*(rf_ux(iz_ixp1)-rf_ux(iz_ix)) +c2x*(rf_ux(iz_ixp2)-rf_ux(iz_ixm1))

                sf_duz_dz = c1z*(sf_uz(izp1_ix)-sf_uz(iz_ix)) +c2z*(sf_uz(izp2_ix)-sf_uz(izm1_ix))
                sf_dux_dx = c1x*(sf_ux(iz_ixp1)-sf_ux(iz_ix)) +c2x*(sf_ux(iz_ixp2)-sf_ux(iz_ixm1))

                rf_dux_dz = c1z*(rf_ux(iz_ix)-rf_ux(izm1_ix)) +c2z*(rf_ux(izp1_ix)-rf_ux(izm2_ix))
                rf_duz_dx = c1x*(rf_uz(iz_ix)-rf_uz(iz_ixm1)) +c2x*(rf_uz(iz_ixp1)-rf_uz(iz_ixm2))

                sf_dux_dz = c1z*(sf_ux(iz_ix)-sf_ux(izm1_ix)) +c2z*(sf_ux(izp1_ix)-sf_ux(izm2_ix))
                sf_duz_dx = c1x*(sf_uz(iz_ix)-sf_uz(iz_ixm1)) +c2x*(sf_uz(iz_ixp1)-sf_uz(iz_ixm2))
                
                glda(j) = glda(j) + (rf_duz_dz+rf_dux_dx)*(sf_duz_dz+sf_dux_dx)
                                    ! rf_duz_dz*sf_duz_dz +rf_duz_dz*sf_dux_dx &
                                    !+rf_dux_dx*sf_duz_dz +rf_dux_dx*sf_dux_dx

                gmu (j) = gmu (j) +2*(rf_duz_dz*sf_duz_dz+rf_dux_dx*sf_dux_dx) &
                                    +(rf_dux_dz+rf_duz_dx)*(sf_dux_dz+sf_duz_dx)
                                  ! +2*rf_duz_dz*sf_duz_dz +2*rf_dux_dx*sf_dux_dx &
                                  !   +rf_dux_dz*sf_dux_dz   +rf_dux_dz*sf_duz_dx &
                                  !   +rf_duz_dx*sf_dux_dz   +rf_duz_dx*sf_duz_dx

            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

end
