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
            'Time-domain ISOtropic 2D ACoustic propagation'//s_NL// &
            '1st-order Momemtum-Strain formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.606 *Vmax/dx'//s_NL// &
            'Required model attributes: vp, rho'//s_NL// &
            'Required field components: pz, px'//s_NL// &
            'Required field components: thta'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Basic gradients: grho(wait) gkpa'

        integer :: nbndlayer=max(2,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=3 !number of basic gradients

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:),allocatable :: buoz, buox, kpa
        real,dimension(:,:),allocatable :: ikpa

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

        c1x=coef(1)/m%dx; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2z=coef(2)/m%dz
        
        wavelet_scaler=self%dt/m%cell_volume

        if_hicks=shot%if_hicks

        if(if_hicks.and.m%is_freesurface) then
            if(shot%src%comp=='ez') then
            if(shot%src%iz>1.and.shot%src%ifz<1) then 
                write(*,*) shot%sindex//"'s iz, ifz = "//num2str(shot%src%iz)//', '//num2str(shot%src%ifz)
                call error("ez src below FS but some Hick interp coef are above FS. Code stop as I don't know how to do Hicks interpolation.")
            endif
            endif

            do i=1,shot%nrcv
                if(shot%rcv(i)%comp=='ez') then
                if(shot%rcv(i)%iz>1.and.shot%rcv(i)%ifz<1) then 
                    write(*,*) shot%sindex//"has receiver iz, ifz = "//num2str(shot%src%iz)//', '//num2str(shot%src%ifz)
                    call error("ez rcv below FS but some Hick interp coef are above FS. Code stop as I don't know how to do Hicks interpolation.")
                endif
                endif

            enddo

        endif

        if_record_adjseismo=either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))

        FS_method=setup%get_str('FS_METHOD',o_default='stress_image')

        call alloc(self%buoz,           [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%buox,           [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(self%kpa,            [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        
        self%kpa(:,:)=cb%rho(:,:,1)*cb%vp(:,:,1)**2

        self%kpa(cb%ifz,:)=self%kpa(cb%ifz+1,:)
        self%kpa(:,cb%ifx)=self%kpa(:,cb%ifx+1)

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
        call f%init_boundary_momenta

        call alloc(f%pz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%px, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%thta, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])

        call alloc(f%dpz_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dpx_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dthta_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        call alloc(f%dthta_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[1,1])
        
        ! !needed for free surface BC
        ! nnz=0-cb%ifz+1
        ! call alloc(f%sz,[cb%ifz,1+nnz],[cb%ifx,cb%ilx],[1,1])
        ! call alloc(f%sx,[cb%ifz,1+nnz],[cb%ifx,cb%ilx],[1,1])

    end subroutine

    subroutine init_correlate(self,corr,name)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name
        
        corr%name=name

        ! if(name(1:1)=='g') then !gradient components
            call alloc(corr%grho,m%nz,m%nx,m%ny)
            call alloc(corr%gkpa,m%nz,m%nx,m%ny)
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
!            call correlate_assemble(corr%grho, correlate_gradient(:,:,:,1))
            call correlate_assemble(corr%gkpa, correlate_gradient(:,:,:,2))
        endif        
        
    end subroutine
    
    !========= Derivations =================
    !PDE:      A u = M ∂ₜ u - MD Mu = f
    !Adjoint:  Aᵀa = M ∂ₜᵀa - MDᵀMa = d
    !where
    !u=[pz px θ]ᵀ, [pz px] are momenta, θ are volumetric strain
    !f=[fz fx p]ᵀ  δ(x-xs) with xs source position, d is recorded data
    !  [b    ]    [0  0  ∂z]
    !M=|  b  |, D=|0  0  ∂ₓ|
    !  [    κ]    [∂z ∂ₓ 0 ]
    !a=[pzᵃ pxᵃ θᵃ]ᵀ is the adjoint field
    !
    !Continuous case:
    !<a|Au> = ∫ a(x,t) (M∂ₜ-MD)u(x,t) dx³dt
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
    !                         |        |     -½ pz       |        |
    !                         |        |        bz       |        |
    !                         |        |        |        |        |
    !                         κ   bx   κ   bx   κ   bx   κ        κ
    !  -pz--θ-pz-θ-pz-→ t    -θ---px---θ---px---θ---px---θ---px---θ-→ x
    !   -1 -½  0 ½ 1          -2  -1½ -1   -½   0    ½   1   1½   2    
    !                         |        |        |        |        | 
    !                         |        |      ½ pz       pz       | 
    !                         |        |        bz       bz       | 
    !                         |        |        |        |        | 
    !                        -|--------|-------1-θ-------θ--------|-
    !                         |        |        κ        κ        | 
    !                         |        |        |        |        | 
    !                         |        |     1½ pz       |        | 
    !                         |        |        bz       |        | 
    !                         |        |        |        |        | 
    !                                         z ↓
    !
    !Convention for half-integer index:
    !(array index)  =>     (real index)     
    !  pz(iz,ix)    =>   pz[iz-½,ix  ]^n   :=pz((iz-½)*dz,ix*dx,n*dt)
    !  px(iz,ix)    =>   px[iz,  ix-½]^n  
    !   θ(iz,ix)    =>    θ[iz,  ix  ]^n+½ := θ(iz*dz,    ix*dx,(n+½)*dt)
    !
    !Forward:
    !FD eqn:
    !    [pz^n  ]   [ 0    0   ∂zᵇ] [pz^n+1]         [∂zᵇκ θ^n+½]   [  ρfz]
    !∂ₜᶠ |px^n  |   | 0    0   ∂ₓᵇ|M|px^n+1| +M⁻¹f = |∂ₓᵇκ θ^n+½| + |  ρfx|
    !    [θ^n+1 ] = [∂zᶠ  ∂ₓᶠ  0  ] [ θ^n+½]         [∂zᶠbpz^n+1]   [κ⁻¹p ]
    !∂ₜᶠ*dt := pz^n+1 - pz^n                               ~O(t²)
    !∂zᵇ*dz := c₁( θ(iz  )- θ(iz-1) +c₂( θ(iz+1)- θ(iz-2)  ~O(x⁴)
    !∂zᶠ*dz := c₁(pz(iz+1)-pz(iz  ) +c₂(pz(iz+2)-pz(iz-1)  ~O(x⁴)
    !
    !Time marching:
    ![pz^n+1 ]   [pz^n  ]   [∂zᵇκθ^n+½              ]
    !|px^n+1 | = |px^n  | + |∂ₓᵇκθ^n+½              |dt +M⁻¹f*dt    
    ![ θ^n+1½]   [ θ^n+½]   [∂zᶠbpz^n+1 + ∂ₓᶠbpx^n+1]
    !Step #1: pz^n += src
    !Step #2: pz^n+1 = pz^n + spatial FD(θ^n+½)
    !Step #3: θ^n+½ += src
    !Step #4: θ^n+1½ = θ^n+½ + spatial FD(pz,px^n+1)
    !Step #5: sample pz,px^n & θ^n+1½ at receivers
    !Step #6: save pz,px^n+1 to boundary values
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
                call cross_correlate_gkpa(fld_a,fld_u,a_star_u,it)
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
            write(*,*) 'Total elapsed time for forward (min)',(tt1+tt2+tt3+tt7+tt8)/60./mpiworld%max_threads
            write(*,*) ' ---------------------------- '
            write(*,*) 'Elapsed time to add adjsource strains    ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj strains       ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to add adjsource momenta    ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj momenta       ',tt10/mpiworld%max_threads
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
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl !required to pass adjointtest. Why weaker when vx as src?
                    f%px(ifz:ilz,ifx:ilx,1) = f%px(ifz:ilz,ifx:ilx,1) + wl/self%buox(ifz:ilz,ifx:ilx) *shot%src%interp_coef(:,:,1)

                case ('vz')
                    f%pz(ifz:ilz,ifx:ilx,1) = f%pz(ifz:ilz,ifx:ilx,1) + wl*shot%src%interp_coef(:,:,1)

                case ('vx')
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                    f%px(ifz:ilz,ifx:ilx,1) = f%px(ifz:ilz,ifx:ilx,1) + wl*shot%src%interp_coef(:,:,1)
                    
                end select
                
            else
                select case (shot%src%comp)
                case ('pz') !vertical force     on pz[iz-0.5,ix]
                    !f%pz(iz,ix,1) = f%pz(iz,ix,1) + wl
                    f%pz(iz,ix,1) = f%pz(iz,ix,1) + wl/self%buoz(iz,ix)
                    
                case ('px') !horizontal x force on px[iz,ix-0.5]
                    !f%px(iz,ix,1) = f%px(iz,ix,1) + wl
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                    f%px(iz,ix,1) = f%px(iz,ix,1) + wl/self%buox(iz,ix)
                    
                case ('vz')
                    f%pz(iz,ix,1) = f%pz(iz,ix,1) + wl

                case ('vx')
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                    f%px(iz,ix,1) = f%px(iz,ix,1) + wl
                    
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
                        f%pz(ifz:ilz,ifx:ilx,1) = f%pz(ifz:ilz,ifx:ilx,1) + wl/self%buoz(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!

                    case ('px') !horizontal x adjsource
                        if(m%is_freesurface.and.shot%rcv(i)%iz==1) wl=2*wl !required to pass adjointtest. Why weaker when vx as src?
                        f%px(ifz:ilz,ifx:ilx,1) = f%px(ifz:ilz,ifx:ilx,1) + wl/self%buox(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef(:,:,1) !no time_dir needed!

                    case ('vz') !vertical z adjsource
                        f%pz(ifz:ilz,ifx:ilx,1) = f%pz(ifz:ilz,ifx:ilx,1) + wl*shot%rcv(i)%interp_coef(:,:,1)
                    case ('vx') !horizontal x adjsource
                        f%px(ifz:ilz,ifx:ilx,1) = f%px(ifz:ilz,ifx:ilx,1) + wl*shot%rcv(i)%interp_coef(:,:,1)
                        
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                    case ('pz') !vertical z adjsource
                        !pz[ix,1,iz-0.5]
                        f%pz(iz,ix,1) = f%pz(iz,ix,1) + wl/self%buoz(iz,ix)  !self%buoz(iz,ix) !no time_dir needed!

                    case ('px') !horizontal x adjsource
                        !px[ix-0.5,1,iz]
                        if(m%is_freesurface.and.shot%rcv(i)%iz==1) wl=2*wl
                        f%px(iz,ix,1) = f%px(iz,ix,1) + wl/self%buox(iz,ix) !no time_dir needed!
                        
                    case ('vz')
                        f%pz(iz,ix,1) = f%pz(iz,ix,1) + wl

                    case ('vx')
                        if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                        f%px(iz,ix,1) = f%px(iz,ix,1) + wl

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

        call fd2d_momenta(f%pz,f%px,f%thta,             &
                        f%dthta_dz,f%dthta_dx,          &
                        self%kpa,                       &
                        ifz,ilz,ifx,ilx,time_dir*self%dt)

        if(m%is_freesurface) then
            if (FS_method=='stress_image') then !Levandar & Roberttson
                !Roberttson's 3rd method
                f%pz(cb%ifz:1,:,1)=0.
                f%px(cb%ifz:0,:,1)=0.

            endif

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
                    f%thta(ifz:ilz,ifx:ilx,1) = f%thta(ifz:ilz,ifx:ilx,1) + wl*self%ikpa(ifz:ilz,ifx:ilx)*shot%src%interp_coef_symm(:,:,1)

                case ('ex')
                    f%thta(ifz:ilz,ifx:ilx,1) = f%thta(ifz:ilz,ifx:ilx,1) + wl*self%ikpa(ifz:ilz,ifx:ilx)*shot%src%interp_coef_symm(:,:,1)

                endselect
                
            else
                select case (shot%src%comp)
                case ('ez')
                    !if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                    f%thta(iz,ix,1) = f%thta(iz,ix,1) + wl*self%ikpa(iz,ix)

                case ('ex')
                    if(m%is_freesurface.and.shot%src%iz==1) wl=2*wl
                    f%thta(iz,ix,1) = f%thta(iz,ix,1) + wl*self%ikpa(iz,ix)
                    
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
                        f%thta(ifz:ilz,ifx:ilx,1) = f%thta(ifz:ilz,ifx:ilx,1) + wl*self%ikpa(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef_symm(:,:,1)
                            

                    case ('ex')
                        f%thta(ifz:ilz,ifx:ilx,1) = f%thta(ifz:ilz,ifx:ilx,1) + wl*self%ikpa(ifz:ilz,ifx:ilx)*shot%rcv(i)%interp_coef_symm(:,:,1)
                    
                    endselect

                else           
                    select case (shot%rcv(i)%comp)
                    case ('ez')
                        !if(m%is_freesurface.and.shot%rcv(i)%iz==1) wl=2*wl
                        f%thta(iz,ix,1) = f%thta(iz,ix,1) + wl*self%ikpa(iz,ix)

                    case ('ex')
                        f%thta(iz,ix,1) = f%thta(iz,ix,1) + wl*self%ikpa(iz,ix)
                        
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
        
        if(m%is_freesurface) ifz=max(ifz,2)

        call fd2d_strains(self%buoz,f%pz,self%buox,f%px,f%thta,&
                           f%dpz_dz,f%dpx_dx,&
                           ifz,ilz,ifx,ilx,time_dir*self%dt)

        if(m%is_freesurface) then !FS condition: sz=ss=0.

            ! if (FS_method=='stress_image') then !Levandar & Roberttson
                
            !     nnz=0-(-4)+1 !0-cb%ifz+1
            !     f%sz(2:1+nnz,:,1)=self%ldap2mu(2:1+nnz,:)*f%ez(2:1+nnz,:,1)+self%kpa    (2:1+nnz,:)*f%ex(2:1+nnz,:,1)
            !     f%sx(2:1+nnz,:,1)=self%kpa    (2:1+nnz,:)*f%ez(2:1+nnz,:,1)+self%ldap2mu(2:1+nnz,:)*f%ex(2:1+nnz,:,1)

            !     !image szz
            !     f%sz(1,:,1)=0.
            !     f%sz(0:-4:-1, :,1)=-f%sz(2:2+0-(-4), :,1)
            !     !f%sz(0:cb%ifz:-1, :,1)=-f%sz(2:2+0-cb%ifz, :,1)

            !     !not image on sxx
            !     f%sx(-4:0,:,1)=0.
            !     !f%sx(cb%ifz:0,:,1)=0.

            !     do ix=cb%ifx+1,cb%ilx-2
            !         dvx_dx_= c1x*(self%buox(1,ix+1)*f%px(1,ix+1,1)-self%buox(1,ix  )*f%px(1,ix  ,1)) &
            !                 +c2x*(self%buox(1,ix+2)*f%px(1,ix+2,1)-self%buox(1,ix-1)*f%px(1,ix-1,1))
            !         !factor = 1./(self%ldap2mu(1,ix)/self%inv_ldapmu_4mu(1,ix)) !less precise
            !         factor = -self%kpa(1,ix)**2/self%ldap2mu(1,ix) + self%ldap2mu(1,ix)
            !         f%sx(1,ix,1) = f%sx(1,ix,1) + time_dir*self%dt * factor*dvx_dx_
            !     enddo

            !     !convert to ez & ex
            !     ![sz]=[λ+2μ λ   ][ez] => [ez]=   1   [λ+2μ -λ  ][sz]
            !     ![sx] [λ    λ+2μ][ex]    [ex] (λ+μ)4μ[-λ   λ+2μ][sx]
            !     f%ez(-4:1,:,1) = self%inv_ldapmu_4mu(-4:1,:)*( self%ldap2mu(-4:1,:)*f%sz(-4:1,:,1)-self%kpa    (-4:1,:)*f%sx(-4:1,:,1) )
            !     f%ex(-4:1,:,1) = self%inv_ldapmu_4mu(-4:1,:)*(-self%kpa    (-4:1,:)*f%sz(-4:1,:,1)+self%ldap2mu(-4:1,:)*f%sx(-4:1,:,1) )
            !     !f%ez(cb%ifz:1,:,1) = self%inv_ldapmu_4mu(cb%ifz:1,:)*( self%ldap2mu(cb%ifz:1,:)*f%sz(cb%ifz:1,:,1)-self%kpa    (cb%ifz:1,:)*f%sx(cb%ifz:1,:,1) )
            !     !f%ex(cb%ifz:1,:,1) = self%inv_ldapmu_4mu(cb%ifz:1,:)*(-self%kpa    (cb%ifz:1,:)*f%sz(cb%ifz:1,:,1)+self%ldap2mu(cb%ifz:1,:)*f%sx(cb%ifz:1,:,1) )
            !
            !
            !     !image on szx = μ*es
            !     f%es(1:-4:-1, :,1)=-f%es(2:2+1-(-4), :,1)
            !     !f%es(1:cb%ifz:-1, :,1)=-f%es(2:2+1-cb%ifz, :,1)


            if (FS_method=='stress_image') then !replicates Levandar & Roberttson method under strain system

                !image szz
                f%thta( 1,:,1)=0.
                f%thta(0:cb%ifz:-1, :,1)=-f%thta(2:2+0-cb%ifz, :,1)

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
                        f%seismo(i,it)=sum(self%buoz(ifz:ilz,ifx:ilx)*f%pz(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                    case ('vx')
                        f%seismo(i,it)=sum(self%buox(ifz:ilz,ifx:ilx)*f%px(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )

                    case ('pz')
                        f%seismo(i,it)=sum(f%pz(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                    case ('px')
                        f%seismo(i,it)=sum(f%px(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef(:,:,1) )
                        
                    case ('ez')
                        f%seismo(i,it)=sum(f%thta(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef_full(:,:,1) )
                    case ('ex')
                        f%seismo(i,it)=sum(f%thta(ifz:ilz,ifx:ilx,1)*shot%rcv(i)%interp_coef_full(:,:,1) )

                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                    case ('vz') !pz[iz-0.5,ix]
                        f%seismo(i,it)=self%buoz(iz,ix)*f%pz(iz,ix,1)
                    case ('vx') !px[iz,ix-0.5]
                        f%seismo(i,it)=self%buox(iz,ix)*f%px(iz,ix,1)

                    case ('pz') !pz[iz-0.5,ix]
                        f%seismo(i,it)=f%pz(iz,ix,1)
                    case ('px') !px[iz,ix-0.5]
                        f%seismo(i,it)=f%px(iz,ix,1)

                    case ('ez') !ez[iz,ix]
                        f%seismo(i,it)=f%thta(iz,ix,1)
                    case ('ex') !ez[iz,ix]
                        f%seismo(i,it)=f%thta(iz,ix,1)

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
                    f%seismo(1,it)=sum(f%pz(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                case ('vx')
                    f%seismo(1,it)=sum(f%px(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))

                case ('pz')
                    f%seismo(1,it)=sum(f%pz(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                case ('px')
                    f%seismo(1,it)=sum(f%px(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef(:,:,1))
                
                case ('ez')
                    f%seismo(1,it)=sum(f%thta(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef_full(:,:,1))
                case ('ex')
                    f%seismo(1,it)=sum(f%thta(ifz:ilz,ifx:ilx,1)*shot%src%interp_coef_full(:,:,1))

                end select
                
            else
                select case (shot%src%comp)
                case ('pz') !pz[iz-0.5,ix,1]
                    f%seismo(1,it)=f%pz(iz,ix,1)
                case ('px') !px[iz,ix-0.5,1]
                    f%seismo(1,it)=f%px(iz,ix,1)
                
                case ('ez') !ez[iz-0.5,ix,1]
                    f%seismo(1,it)=f%thta(iz,ix,1)
                case ('ex') !ex[iz-0.5,ix,1]
                    f%seismo(1,it)=f%thta(iz,ix,1)

                end select
                
            endif
        
    end subroutine
    
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%buoz, self%buox)
        call dealloc(self%ikpa, self%kpa)
        ! call dealloc(self%two_ldapmu, self%inv_ldapmu_4mu)
    end subroutine

    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !<a|Au> = <a|M∂ₜu-MDMu> = <a|M∂ₜu> - <a|MDMu>
    !   Kₘ<a|M∂ₜu> = <a|(KₘM)∂ₜu> = <a|(KₘM)(DMu+f)> ≐ <a|(KₘM)DMu>
    !   Kₘ<a|MDMu> = <a|(KₘM)DMu> + <a|M(KₘDM)u> = <a|(KₘM)DMu> - <DMa|(KₘM)u>
    !So Kₘ<a|Au> ≐ <DMa|(KₘM)u> =∫ (DMa)ᵀ KₘM u dt
    !
    !  [b     ]      [∂zᵇκθ          ]  
    !M=|  b   |, DMa=|∂ₓᵇκθ          |,  
    !  |    κ |      [∂zᶠbpz + ∂ₓᶠbpx]
    !
    !K_bM=diag{1,1,0}, K_κM=diag{0,0,1}
    !
    !Therefore,
    !gkpa =  ∂zᶠbpzᵃ*θ + ∂ₓᶠbpxᵃ*θ  ie. ∇b\vec{pᵃ}⋅θ

    subroutine cross_correlate_gkpa(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2) !
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ! ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ! ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
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
            call grad2d_gkpa(ppg%buoz*rf%pz(:,:,1),ppg%buox*rf%px(:,:,1),&
                                sf%thta,&
                                corr%gkpa,&
                                ifz,ilz,ifx,ilx)
        endif
        
    end subroutine
    
    subroutine cross_correlate_postprocess(corr)
        type(t_correlate) :: corr

        if(allocated(correlate_gradient)) then        
            !preparing for projection back
            corr%gkpa(1,:,:) = corr%gkpa(2,:,:)
        endif

    end subroutine

    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd2d_momenta(pz,px,thta,      &
                           dthta_dz,dthta_dx,&
                           kpa,              &
                           ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: pz,px,thta
        real,dimension(*) :: dthta_dz,dthta_dx
        real,dimension(*) :: kpa
        
        nz=cb%nz
        nx=cb%nx
        
        dthta_dz_=0.
        dthta_dx_=0.

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dthta_dz_,dthta_dx_)
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
                
                

                dthta_dz_= c1z*(kpa(iz_ix)*thta(iz_ix)-kpa(izm1_ix)*thta(izm1_ix)) +c2z*(kpa(izp1_ix)*thta(izp1_ix)-kpa(izm2_ix)*thta(izm2_ix))
                dthta_dx_= c1x*(kpa(iz_ix)*thta(iz_ix)-kpa(iz_ixm1)*thta(iz_ixm1)) +c2x*(kpa(iz_ixp1)*thta(iz_ixp1)-kpa(iz_ixm2)*thta(iz_ixm2))
                                
                !cpml
                dthta_dz(i)= cpml%b_z_half(iz)*dthta_dz(i) + cpml%a_z_half(iz)*dthta_dz_
                dthta_dx(i)= cpml%b_x_half(ix)*dthta_dx(i) + cpml%a_x_half(ix)*dthta_dx_

                dthta_dz_=dthta_dz_*cpml%kpa_z_half(iz) + dthta_dz(i)
                dthta_dx_=dthta_dx_*cpml%kpa_x_half(ix) + dthta_dx(i)  !kappa's should have been inversed in m_computebox.f90
                
                !momenta
                pz(i)=pz(i) + dt*dthta_dz_
                px(i)=px(i) + dt*dthta_dx_

            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd2d_strains(buoz,pz,buox,px,thta,  &
                             dpz_dz,dpx_dx,     &
                             ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: buoz,pz,buox,px,thta
        real,dimension(*) :: dpz_dz,dpx_dx
        
        nz=cb%nz
        nx=cb%nx
        
        dpz_dz_=0.
        dpx_dx_=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dpz_dz_,dpx_dx_)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1
                
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                
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
                thta(i) = thta(i) + dt * (dpz_dz_+dpx_dx_)

            enddo
            
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine

    subroutine grad2d_gkpa(rf_bpz,rf_bpx,&
                            sf_thta,&
                            grad_kpa,&
                            ifz,ilz,ifx,ilx)
        real,dimension(*) :: rf_bpz,rf_bpx
        real,dimension(*) :: sf_thta
        real,dimension(*) :: grad_kpa
        
        nz=cb%nz
        
        rf_dz_bpz=0.; rf_dx_bpx=0.

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         rf_dz_bpz,rf_dx_bpx)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers

                izm1_ix= i-1  !iz-1,ix
                iz_ix  = i    !iz  ,ix
                izp1_ix= i+1  !iz+1,ix
                izp2_ix= i+2  !iz+2,ix

                iz_ixm1  = i      -nz  !iz  ,ix-1        
                iz_ixp1  = i      +nz  !iz  ,ix+1
                iz_ixp2  = i    +2*nz  !iz  ,ix+2
                
                rf_dz_bpz = c1z*(rf_bpz(izp1_ix)-rf_bpz(iz_ix)) + c2z*(rf_bpz(izp2_ix)-rf_bpz(izm1_ix))
                rf_dx_bpx = c1x*(rf_bpx(iz_ixp1)-rf_bpx(iz_ix)) + c2x*(rf_bpx(iz_ixp2)-rf_bpx(iz_ixm1))

                grad_kpa(j) = grad_kpa(j) + (rf_dz_bpz+rf_dx_bpx) * sf_thta(i)

            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

end
