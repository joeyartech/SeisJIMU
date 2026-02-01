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

    ! !local const
    ! real :: dt2, inv_2dt, inv_2dz, inv_2dx


    !scaling source wavelet
    real :: wavelet_scaler

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain ISOtropic 2D SH wave propagation'//s_NL// &
            '1st-order Velocity-Stress formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.606(for 2D) or 0.494(3D) *Vmax/dx'//s_NL// &
            'Required model attributes: vs, rho'//s_NL// &
            'Required field components: szy, sxy, vy'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Imaging conditions: ipp'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: grho gmu'

        integer :: nbndlayer=max(2,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=2 !number of basic gradients
        integer :: nimag=1 !number of basic images
        integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:),allocatable :: muz, mux, buo

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
        
        procedure :: inject_stresses
        procedure :: inject_velocities
        procedure :: update_stresses
        procedure :: update_velocities
        procedure :: extract


        final :: final

    end type

    type(t_propagator),public :: ppg

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    logical :: if_record_adjseismo=.false.

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
        
        if(index(self%info,'vs')>0  .and. .not. allocated(m%vs)) then
            call alloc(m%vs,m%nz,m%nx,1,o_init=866.)
            call warn('Constant vs model (866 m/s) is allocated by propagator.')
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

        real,dimension(:,:),allocatable :: temp_mu

        c1x=coef(1)/m%dx; c1y=coef(1)/m%dy; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz

        ! inv_2dz =1./2/m%dz
        ! inv_2dx =1./2/m%dx
        
        wavelet_scaler=self%dt/m%cell_volume

        if_hicks=shot%if_hicks

        call alloc(self%buo,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%muz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%mux,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        
        call alloc(temp_mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        temp_mu = cb%rho(:,:,1)*cb%vs(:,:,1)**2

        do iz=cb%ifz+1,cb%ilz
            self%muz(iz,:)=(temp_mu(iz,:)+temp_mu(iz-1,:))/2.
        enddo
        
        do ix=cb%ifx+1,cb%ilx
            self%mux(:,ix)=(temp_mu(:,ix)+temp_mu(:,ix-1))/2.
        enddo

        deallocate(temp_mu)

        self%buo=1./cb%rho(:,:,1)

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
        call f%init_boundary_stresses

        call alloc(f%szy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%szx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%vy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])

        call alloc(f%dszy_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dsxy_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dvy_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dvy_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
                
    end subroutine

    subroutine init_correlate(self,corr,name)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name
        
        corr%name=name

        ! if(name(1:1)=='g') then !gradient components
            call alloc(corr%grho,m%nz,m%nx,m%ny)
            call alloc(corr%gimu,m%nz,m%nx,m%ny)
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
            call correlate_assemble(corr%grho, correlate_gradient(:,:,:,1))
            call correlate_assemble(corr%gimu, correlate_gradient(:,:,:,2))
        endif        
        
    end subroutine

    
    !========= Derivations =================
    !PDE:      A u = M ∂ₜ u - D u = f
    !Adjoint:  Aᵀa = M ∂ₜᵀa - Dᵀa = d
    !where
    !u=[szy sxy vy]ᵀ, [szy sxy] are shear stresses, vy is velocity
    !f=[mzy mxy fy]ᵀδ(x-xs) is the source term, d is recorded data
    !M=[diag(μ⁻¹) ρ], N=M⁻¹=[diag(μ) b], b=ρ⁻¹ is buoancy,
    !  [0  0  ∂z]
    !D=|0  0  ∂ₓ|
    !  [∂z ∂ₓ 0 ]
    !and a=[szyᵃ sxyᵃ vyᵃ]ᵀ is the adjoint field
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
    !Meshing with staggered grids in time and space:
    !                       |        |    -½ szy       |        |
    !                       |        |       μz        |        |
    !                       |        |        |        |        |
    !                      buo  μx  buo  μx  buo  μx  buo  μx  buo
    !  -s-vy-s-vy-s-→ t   -vy--sxy--vy--sxy--vy--sxy--vy--sxy--vy-→ x
    !  -1 -½ 0  ½ 1        -2  -1½  -1   -½   0   ½    1   1½   2    
    !                       |        |        |        |        | 
    !                       |        |     ½ szy       |        | 
    !                       |        |       μz        |        |
    !                       |        |        |        |        | 
    !                      -|--------|-----1-vy--------|--------|-
    !                       |        |       buo       |        | 
    !                       |        |        |        |        | 
    !                       |        |    1½ szy       |        | 
    !                       |        |       μz        |        | 
    !                       |        |        |        |        | 
    !                                       z ↓
    !
    !Convention for half-integer index:
    !(array index)  =>    (real index)     
    ! vy(iz,ix)  =>  vy[iz,  ix, ]^n+½ := vy(iz*dz,ix*dx,(n+½)*dt)
    !szy(iz,ix)  => szy[iz-½,ix, ]^n   :=szy((iz-½)*dz,ix*dx,n*dt)
    !sxy(iz,ix)  => sxy[iz,  ix-½]^n   
    !
    ! Forward:
    ! FD eqn:
    !      [vy^n+1]   [∂zᶠ   ∂ₓᶠ    0 ][ p^n+½]      [∂zᶠ vz^n+1 + ∂ₓᶠ vx^n+1]
    ! M ∂ₜᶠ |sxy^n  | = | 0     0    ∂ₓᵇ||vx^n+1| +f = |∂ₓᵇ p^n+½              | +f
    !      [szy^n  ]   [ 0     0    ∂zᵇ][vz^n+1]      [∂zᵇ p^n+½              ]
    ! where
    ! ∂ₜᶠ*dt := v^n+1 - v^n                               ~O(t²)
    ! ∂zᵇ*dz := c₁(p(iz  )-p(iz-1)) +c₂(p(iz+1)-p(iz-2))  ~O(x⁴)
    ! ∂zᶠ*dz := c₁(v(iz+1)-v(iz  )) +c₂(v(iz+2)-v(iz-1))  ~O(x⁴)
    
!!    ! Time marching:
!!    ! [vz^n+1 ]   [vz^n  ]      [∂zᵇ p^n+½              ]
!!    ! |vx^n+1 | = |vx^n  | + M⁻¹|∂ₓᵇ p^n+½              |dt  +M⁻¹f*dt
!!    ! [ p^n+1½]   [ p^n+½]      [∂zᶠ vz^n+1 + ∂ₓᶠ vx^n+1]
!!    ! Step #1: v^n += src
!!    ! Step #2: v^n+1 = v^n + spatial FD(p^n+½)
!!    ! Step #3: p^n+½ += src
!!    ! Step #4: p^n+1½ = p^n+½ + spatial FD(v^n+1)
!!    ! Step #5: sample v^n & p^n++½ at receivers
!!    ! Step #6: save v^n+1 to boundary values

!!    ! Reverse time marching (for wavefield reconstruction)
!!    ! [ p^n+½]   [ p^n+1½]      [∂zᶠ vz^n+1 + ∂ₓᶠ vx^n+1]
!!    ! |vx^n  | = |vx^n+1 | - M⁻¹|∂ₓᵇ p^n+½              |dt  -M⁻¹f*dt
!!    ! [vz^n  ]   [vz^n+1 ]      [∂zᵇ p^n+½              ]
!!    ! Step #6: load boundary values for v^n+1
!!    ! Step #4: p^n+½ = p^n+1½ - spatial FD(v^n+1)
!!    ! Step #3: p^n+½ -= src
!!    ! Step #2: v^n+1 = v^n - spatial FD(p^n+½)
!!    ! Step #1: v^n -= src
!!    ! N.B. Same codes for spatial FDs as in forward time marching, with a negated dt.
    
    ! Adjoint:
    ! FD eqn:
!!    !      [vzᵃ^n  ]   [ 0      0     ∂zᶠᵀ][vzᵃ^n+1]
!!    ! M ∂ₜᶠᵀ|vxᵃ^n  | = | 0      0     ∂ₓᶠᵀ||vxᵃ^n+1|  +d
!!    !      [ pᵃ^n+1]   [∂zᵇᵀ   ∂ₓᵇᵀ    0  ][ pᵃ^n+½]
!!    ! ∂ₜᶠᵀ = v^n-1 -v^n   = -∂ₜᵇ
!!    ! ∂zᵇᵀ = c₁(v[i  ]-v[i+½]) +c₂(v[i- ½]-v[i+1½]) = -∂zᶠ
!!    ! ∂zᶠᵀ = c₁(p[i-½]-p[i  ]) +c₂(p[i-1½]-p[i+ ½]) = -∂zᵇ
!!    !      [vzᵃ^n  ]   [ 0      0     ∂zᵇ ][vzᵃ^n+1]
!!    ! M ∂ₜᵇ |vxᵃ^n  | = | 0      0     ∂ₓᵇ ||vxᵃ^n+1|  -d
!!    !      [ pᵃ^n+1]   [∂zᶠ    ∂ₓᶠ    0   ][ pᵃ^n+½]
!!    ! ie. Dᵀ=-D, antisymmetric
    
!!    ! Time marching:
!!    ! [vzᵃ^n+1 ]   [vzᵃ^n  ]      [∂zᵇ pᵃ^n+½               ]
!!    ! |vxᵃ^n+1 | = |vxᵃ^n  | + M⁻¹|∂ₓᵇ pᵃ^n+½               |dt  -M⁻¹d*dt
!!    ! [ pᵃ^n+1½]   [ pᵃ^n+½]      [∂zᶠ vzᵃ^n+1 + ∂ₓᶠ vxᵃ^n+1]
!!    ! but we have to do it in reverse time:
!!    ! [ pᵃ^n+½]   [ pᵃ^n+1½]      [∂zᶠ vzᵃ^n+1 + ∂ₓᶠ vxᵃ^n+1]
!!    ! |vxᵃ^n  | = |vxᵃ^n+1 | - M⁻¹|∂ₓᵇ pᵃ^n+½               |dt  +M⁻¹d*dt
!!    ! [vzᵃ^n  ]   [vzᵃ^n+1 ]      [∂zᵇ pᵃ^n+½               ]
!!    ! Step #5: pᵃ^n+1½ += adjsrc
!!    ! Step #4: pᵃ^n+½ = pᵃ^n+1½ - spatial FD(vᵃ^n+1)
!!    ! Step #3: vᵃ^n+1 += adjsrc
!!    ! Step #2: vᵃ^n = vᵃ^n+1 - spatial FD(pᵃ^n+½)
!!    ! N.B. Same codes for spatial FDs as in forward time marching, with a negated dt, but the RHS should use a "+" sign (regardless of the reverse time direction).
    
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
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_u%check_value(fld_u%vz)
            endif

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add forces to v^it
            call cpu_time(tic)
            call self%inject_stresses(fld_u,time_dir,it)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2: from v^it to v^it+1 by differences of s^it+0.5
            call cpu_time(tic)
            call self%update_stresses(fld_u,time_dir,it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            !step 3: add pressure to s^it+0.5
            call cpu_time(tic)
            call self%inject_velocities(fld_u,time_dir,it)
            call cpu_time(toc)
            tt3=tt3+toc-tic

            !step 4: from s^it+0.5 to s^it+1.5 by differences of v^it+1
            call cpu_time(tic)
            call self%update_velocities(fld_u,time_dir,it)
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
                call fld_u%boundary_transport_stresses('save',it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            ! endif

        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add stress source  ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update stresses    ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to add source stresses',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update velocities  ',tt4/mpiworld%max_threads
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
        
        !for adjoint test
        if(if_record_adjseismo)  call alloc(fld_a%seismo,1,self%nt)

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
                call fld_a%check_value(fld_a%vz)
                call fld_u%check_value(fld_u%vz)
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then

                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_u%boundary_transport_stresses('load',it)
                call cpu_time(toc)
                tt1=tt1+toc-tic
                
                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%update_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_velocities(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_velocities(fld_a,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: rf%s^it+0.5 star D sf%s_dt^it+0.5
            !use sf%v^it+1 to compute sf%s_dt^it+0.5, as backward step

            if(mod(it,irdt)==0) then
                call cpu_time(tic)
                call cross_correlate_grho(fld_a,fld_u,a_star_u,it)
                ! call cross_correlate_image(fld_a,fld_u,a_star_u,it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif
                            
            !========================================================!

            ! if(present(o_sf)) then
                !backward step 2: v^it+1 -> v^it by FD of s^it+0.5
                call cpu_time(tic)
                call self%update_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt7=tt7+toc-tic

                !backward step 1: rm forces from v^it
                call cpu_time(tic)
                call self%inject_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt8=tt8+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 3: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_stresses(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_stresses(fld_a,time_dir,it)
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
           if(mod(it,irdt)==0) then
               call cpu_time(tic)
               call cross_correlate_gimu(fld_a,fld_u,a_star_u,it)
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
            write(*,*) 'Elapsed time to update velocities     ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source velocities  ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update stresses       ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source stresses    ',tt8/mpiworld%max_threads
            write(*,*) 'Total elapsed time for forward (min)',(tt1+tt2+tt3+tt7+tt8)/60./mpiworld%max_threads
            write(*,*) ' ---------------------------- '
            write(*,*) 'Elapsed time to add adjsource velocities ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj velocities    ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to add adjsource stresses   ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj stresses      ',tt10/mpiworld%max_threads
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

            if(shot%src%comp=='vy') then

                ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
                ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
                
                wl=time_dir*f%wavelet(1,it)*wavelet_scaler
            
                if(if_hicks) then
                    f%vy(ifz:ilz,ifx:ilx,1) = f%vy(ifz:ilz,ifx:ilx,1) + wl*self%buo(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)
                    
                else
                    f%vy(iz,ix,1) = f%vy(iz,ix,1) + wl*self%buo(iz,ix)
                
                endif
            
            endif

            return

        endif

            do i=1,shot%nrcv

                if(shot%rcv(i)%comp=='vy') then  !horizontal y adjsource !vy[ix,iy-0.5,iz]
                    
                    ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                    ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                    
                    wl=f%wavelet(i,it)*wavelet_scaler
                    
                    if(if_hicks) then
                        f%vy(ifz:ilz,ifx:ilx,1) = f%vy(ifz:ilz,ifx:ilx,1) + wl*self%buo(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                        
                    else
                        f%vy(iz,ix,1) = f%vy(iz,ix,1) + wl*self%buo(iz,ix) !no time_dir needed!
                    
                    endif
               
                endif

            enddo
        
    end subroutine
    
    !forward: v^it -> v^it+1 by FD  of s^it+0.5
    !adjoint: v^it+1 -> v^it by FDᵀ of s^it+0.5
    subroutine update_stresses(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ifz=f%bloom(1,it)+2
        ilz=f%bloom(2,it)-1
        ifx=f%bloom(3,it)+2
        ilx=f%bloom(4,it)-1
        ify=f%bloom(5,it)+2
        ily=f%bloom(6,it)-1

        if(m%is_freesurface) ifz=max(ifz,1)

        call fd2d_stresses(f%szy,f%sxy,f%vy,               &
                           f%dszy_dz,f%dsxy_dx,            &
                           self%muz,self%mux,              &
                           ifz,ilz,ifx,ilx,time_dir*self%dt)

        !apply free surface boundary condition if needed
        !Levandar & Roberttson's stress image method
        f%sxy(1,:,1)=0.
        f%sxy(0:cb%ifz:-1, :,1)=-f%sxy(2:2+0-cb%ifz, :,1)

        !image szx
        f%szy(1:cb%ifz:-1, :,1)=-f%szy(2:2+1-cb%ifz, :,1)
        
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
                if(shot%src%comp=='szy') then
                    f%szy(ifz:ilz,ifx:ilx,1) = f%szy(ifz:ilz,ifx:ilx,1) + wl*self%muz(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)

                else if(shot%src%comp=='sxy') then
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl !required to pass adjointtest.
                    f%sxy(ifz:ilz,ifx:ilx,1) = f%sxy(ifz:ilz,ifx:ilx,1) + wl*self%mux(ifz:ilz,ifx:ilx)*shot%src%interp_coef(:,:,1)
                
                endif
                
            else
                if(shot%src%comp=='szy') then
                    f%szy(iz,ix,1) = f%szy(iz,ix,1) + wl*self%muz(iz,ix)
                
                else if(shot%src%comp=='sxy') then
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl !required to pass adjointtest.
                    f%sxy(iz,ix,1) = f%sxy(iz,ix,1) + wl*self%mux(iz,ix)
                
                endif
                
            endif

            return

        endif

            do i=1,shot%nrcv

                ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                
                !adjsource for pressure
                wl=f%wavelet(i,it)*wavelet_scaler
                
                if(if_hicks) then 

                    if(shot%rcv(i)%comp=='szy') then
                        f%szy(ifz:ilz,ifx:ilx,1) = f%szy(ifz:ilz,ifx:ilx,1) +wl*self%muz(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                    elseif(shot%rcv(i)%comp=='sxy') then
                        f%sxy(ifz:ilz,ifx:ilx,1) = f%sxy(ifz:ilz,ifx:ilx,1) +wl*self%mux(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!
                    endif

                else

                    if(shot%rcv(i)%comp=='szy') then
                        !szy[iz,ix,1]
                        f%szy(iz,ix,1) = f%szy(iz,ix,1) +wl*self%muz(iz,ix) !no time_dir needed!
                    elseif(shot%rcv(i)%comp=='sxy') then
                        !sxy[iz,ix,1]
                        f%sxy(iz,ix,1) = f%sxy(iz,ix,1) +wl*self%mux(iz,ix) !no time_dir needed!
                    endif

                endif

            enddo
        
    end subroutine

    !forward: s^it+0.5 -> s^it+1.5 by FD of v^it+1
    !adjoint: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
    subroutine update_velocities(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        ifz=f%bloom(1,it)+1
        ilz=f%bloom(2,it)-2
        ifx=f%bloom(3,it)+1
        ilx=f%bloom(4,it)-2
        
        if(m%is_freesurface) ifz=max(ifz,1)

        call fd2d_velocities(f%szy,f%sxy,f%vy,              &
                            f%dvy_dz,f%dvy_dx,              &
                            self%buo,                       &
                            ifz,ilz,ifx,ilx,time_dir*self%dt)
        
        !apply free surface boundary condition if needed
        !Levandar & Roberttson's stress image method
        f%vy(cb%ifz:0,:,1)=0.

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
                        
                        case ('vy')
                        f%seismo(i,it)=sum(f%vy(ifz:ilz,ifx:ilx,1) *shot%rcv(i)%interp_coef(:,:,1))
                        case ('szy')
                        f%seismo(i,it)=sum(f%szy(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1))
                        case ('sxy')
                        f%seismo(i,it)=sum(f%sxy(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1))
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                        case ('vy') !p[iz,ix,1]
                        f%seismo(i,it)=f%vy(iz,ix,1)
                        case ('szy') !vz[iz-0.5,ix,1]
                        f%seismo(i,it)=f%szy(iz,ix,1)
                        case ('sxy') !vx[iz,ix-0.5,1]
                        f%seismo(i,it)=f%sxy(iz,ix,1)
                    end select
                    
                endif

            enddo

            return

        endif

            ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
            ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
            
            if(if_hicks) then
                select case (shot%src%comp)
                    case ('vy')
                    f%seismo(1,it)=sum(f%vy(ifz:ilz,ifx:ilx,1) *shot%src%interp_coef(:,:,1))
                    
                    case ('szy')
                    f%seismo(1,it)=sum(f%szy(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                    
                    case ('sxy')
                    f%seismo(1,it)=sum(f%sxy(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                    
                end select
                
            else
                select case (shot%src%comp)
                    case ('vy') !p[iz,ix,1]
                    f%seismo(1,it)=f%vy(iz,ix,1)
                    
                    case ('szy') !vz[iz-0.5,ix,1]
                    f%seismo(1,it)=f%szy(iz,ix,1)
                    
                    case ('sxy') !vx[iz,ix-0.5,1]
                    f%seismo(1,it)=f%sxy(iz,ix,1)
                    
                end select
                
            endif
        
    end subroutine
    
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%muz, self%mux, self%buo)
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
        
        call grad2d_grho(rf%vy,sf%szy,sf%sxy,&
                         corr%grho,          &
                         ifz,ilz,ifx,ilx)
        
    end subroutine

    subroutine cross_correlate_gimu(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        
        !inexact greadient
        call grad2d_gimu(rf%szy,rf%sxy,sf%vy,&
                         corr%gimu,          &
                         ifz,ilz,ifx,ilx     )

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
            !scale gradients by model parameters
            corr%grho(:,:,1) = corr%grho(:,:,1) / cb%rho(1:cb%mz,1:cb%mx,1)
            corr%gimu(:,:,1) = corr%gimu(:,:,1) * cb%rho(1:cb%mz,1:cb%mx,1)*cb%vp(1:cb%mz,1:cb%mx,1)**2 !mu
            
            !remove singular point at the src position,
            !because we didn't consider src when deriving the gradient formula    
            iz=shot%src%iz-cb%ioz+1
            ix=shot%src%ix-cb%iox+1
            !for point source
            !safeguards
            ifz=either(iz-2,iz,iz>=3); ilz=either(iz+2,iz,iz<=m%nz-2)
            ifx=either(ix-2,ix,ix>=3); ilx=either(ix+2,ix,ix<=m%nx-2)
            
            corr%grho(iz,ix,1)=0 !first remove otherwise will appear in the sum below
            corr%gimu(iz,ix,1)=0 !first remove otherwise will appear in the sum below
            ncells=size(corr%grho(ifz:ilz,ifx:ilx,1))-1
            corr%grho(iz,ix,1) = sum(corr%grho(ifz:ilz,ifx:ilx,1))/ncells
            corr%gimu(iz,ix,1) = sum(corr%gimu(ifz:ilz,ifx:ilx,1))/ncells
            
            !remove singular top boundary..
            corr%grho(1,:,:) = corr%grho(2,:,:)
            corr%gimu(1,:,:) = corr%gimu(2,:,:)

        endif

        ! if(allocated(correlate_image)) then
        !     corr%ipp (1,:,:) = corr%ipp (2,:,:)
        !     corr%ibksc(1,:,:) = corr%ibksc(2,:,:)
        !     corr%ifwsc(1,:,:) = corr%ifwsc(2,:,:)
        ! endif

    end subroutine

    !========= Finite-Difference on flattened arrays ==================

    subroutine fd2d_stresses(szy,sxy,vy,       &
                             dszy_dz,dsxy_dx,  &
                             muz,mux,          &
                             ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: szy,sxy,vy
        real,dimension(*) :: dszy_dz,dsxy_dx
        real,dimension(*) :: muz,mux
        
        nz=cb%nz
        nx=cb%nx
        
        dszy_dz_=0.; dsxy_dx_=0.

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dszy_dz_,dsxy_dx_)
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

                dszy_dz_= c1z*(szy(iz_ix)-szy(izm1_ix)) +c2z*(szy(izp1_ix)-szy(izm2_ix))
                dsxy_dx_= c1x*(sxy(iz_ix)-sxy(iz_ixm1)) +c2x*(sxy(iz_ixp1)-sxy(iz_ixm2))

                !cpml
                dszy_dz(iz_ix)= cpml%b_z_half(iz)*dszy_dz(iz_ix) + cpml%a_z_half(iz)*dszy_dz_
                dsxy_dx(iz_ix)= cpml%b_x_half(ix)*dsxy_dx(iz_ix) + cpml%a_x_half(ix)*dsxy_dx_

                dszy_dz_=dszy_dz_*cpml%kpa_z_half(iz) + dszy_dz(iz_ix)
                dsxy_dx_=dsxy_dx_*cpml%kpa_x_half(ix) + dsxy_dx(iz_ix)

                !velocity
                szy(iz_ix)=szy(iz_ix) + dt*muz(iz_ix)*dszy_dz_
                sxy(iz_ix)=sxy(iz_ix) + dt*mux(iz_ix)*dsxy_dx_

            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd2d_velocities(szy,sxy,vy,       &
                               dvy_dz,dvy_dx,    &
                               buo,              &
                               ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: szy,sxy,vy
        real,dimension(*) :: dvy_dz,dvy_dx
        real,dimension(*) :: buo
        
        nz=cb%nz
        nx=cb%nx
        
        dvy_dz_=0.;dvy_dx_=0.
        
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
                
                dvy_dz_= c1z*(vy(izp1_ix)-vy(iz_ix)) +c2z*(vy(izp2_ix)-vy(izm1_ix))
                dvy_dx_= c1x*(vy(iz_ixp1)-vy(iz_ix)) +c2x*(vy(iz_ixp2)-vy(iz_ixm1))
                
                !cpml
                dvy_dz(iz_ix)=cpml%b_z(iz)*dvy_dz(iz_ix)+cpml%a_z(iz)*dvy_dz_
                dvy_dx(iz_ix)=cpml%b_x(ix)*dvy_dx(iz_ix)+cpml%a_x(ix)*dvy_dx_

                dvy_dz_=dvy_dz_*cpml%kpa_z(iz) + dvy_dz(iz_ix)
                dvy_dx_=dvy_dx_*cpml%kpa_x(ix) + dvy_dx(iz_ix)
                
                !pressure
                vy(iz_ix) = vy(iz_ix) + dt*buo(iz_ix)*(dvy_dz_+dvy_dx_)
                
            enddo
            
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine
    

    subroutine grad2d_gimu(rf_szy,rf_sxy,sf_vy,&
                           grad,               &
                           ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_szy,rf_sxy,sf_vy
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        dvy_dz=0.
        dvy_dx=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dvy_dz,dvy_dx)
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
                
                dvy_dz = c1z*(sf_vy(izp1_ix)-sf_vy(iz_ix)) +c2z*(sf_vy(izp2_ix)-sf_vy(izm1_ix))
                dvy_dx = c1x*(sf_vy(iz_ixp1)-sf_vy(iz_ix)) +c2x*(sf_vy(iz_ixp2)-sf_vy(iz_ixm1))

                grad(j)=grad(j) + rf_szy(i)*dvy_dz + rf_sxy(i)*dvy_dx
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine
    
    subroutine grad2d_grho(rf_vy,sf_szy,sf_sxy,&
                           grad,               &
                           ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_vy,sf_szy,sf_sxy
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        dz_dszy=0.; dx_dsxy=0.

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dz_dszy, dx_dsxy)
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
                
                ! rvz = rf_vz(iz_ix) +rf_vz(izp1_ix)
                ! rvx = rf_vx(iz_ix) +rf_vx(iz_ixp1)
                
                ! dsvz = (c1z*(sf_p(iz_ix  )-sf_p(izm1_ix)) +c2z*(sf_p(izp1_ix)-sf_p(izm2_ix))) &
                !       +(c1z*(sf_p(izp1_ix)-sf_p(iz_ix  )) +c2z*(sf_p(izp2_ix)-sf_p(izm1_ix)))
                ! dsvx = (c1x*(sf_p(iz_ix  )-sf_p(iz_ixm1)) +c2x*(sf_p(iz_ixp1)-sf_p(iz_ixm2))) &
                !       +(c1x*(sf_p(iz_ixp1)-sf_p(iz_ix  )) +c2x*(sf_p(iz_ixp2)-sf_p(iz_ixm1)))
                !complete equation with unnecessary terms e.g. sf_p(iz_ix) for better understanding
                !with flag -Ox, the compiler should automatically detect such possible simplification
                
                dz_dszy = (c1z*(sf_szy(iz_ix)-sf_szy(izm1_ix)) +c2z*(sf_szy(izp1_ix)-sf_szy(izm2_ix)))
                dx_dsxy = (c1x*(sf_sxy(iz_ix)-sf_sxy(iz_ixm1)) +c2x*(sf_sxy(iz_ixp1)-sf_sxy(iz_ixm2)))

                grad(j)=grad(j) + rf_vy(i)*(dz_dszy + dx_dsxy) !0.25*( rvz*dsvz + rvx*dsvx )
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine

end
