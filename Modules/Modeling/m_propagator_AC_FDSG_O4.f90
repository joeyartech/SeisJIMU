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
            'Time-domain ISOtropic 2D/3D ACoustic propagation'//s_NL// &
            '1st-order Velocity-Stress formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x4,t2) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt <= 0.606(2D) or 0.494(3D) *Vmax/dx'//s_NL// &
            'Required model attributes: vp, rho'//s_NL// &
            'Required field components: vz, vx, vy(for 3D), p'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Basic gradients: grho gkpa'//s_NL// &
            'Imaging conditions: P-Pxcorr'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'

        integer :: nbndlayer=max(2,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=2 !number of basic gradients
        integer :: nimag=1 !number of basic images
        integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !shorthand for greek letters
        !alfa bta gma del(dta) eps
        !zta eta thta iota 
        !kpa lda mu nu
        !xi omi pi rho
        !sgma tau ups phi chi psi oga
        !
        !buo : buoyancy
        !vx,vy : horizontal velocities
        !leading _ : inverse

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
        procedure :: init_abslayer

        procedure :: forward
        procedure :: forward_scattering
        procedure :: inject_velocities
        procedure :: update_velocities
        procedure :: inject_stresses
        procedure :: inject_stresses_scattering
        procedure :: update_stresses
        procedure :: extract

        procedure :: adjoint
        procedure :: adjoint_du_star_a
        procedure :: adjoint_scattering_u_star_da
        procedure :: adjoint_scattering_du_star_da
        procedure :: inject_stresses_adjoint
        procedure :: update_stresses_adjoint
        procedure :: inject_velocities_adjoint
        procedure :: update_velocities_adjoint
        procedure :: extract_adjoint

        final :: final

    end type

    type(t_propagator),public :: ppg

    logical :: if_hicks
    integer :: irdt
    real :: rdt

! real,dimension(:,:,:),allocatable :: sf_p_save

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
            call warn('Shot# '//shot%sindex//' can have grid dispersion!'//s_NL// &
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
            call warn('CFL > 1 on Shot# '//shot%sindex//'!'//s_NL//&
                'vmax, dt, 1/dx = '//num2str(cb%velmax)//', '//num2str(self%dt)//', '//num2str(m%rev_cell_diagonal) //s_NL//&
                'Will adjust dt s.t. CFL < 1')

            self%dt = setup%get_real('CFL',o_default='0.9')/(sumcoef*cb%velmax*m%rev_cell_diagonal)
            self%nt=nint(time_window/self%dt)+1

            write(*,*) 'Shot# '//shot%sindex//': Adjusted dt, nt =',self%dt,self%nt

        endif       
        
    end subroutine

    subroutine init(self)
        class(t_propagator) :: self

        c1x=coef(1)/m%dx; c1y=coef(1)/m%dy; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz
        
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
                
    end subroutine

    subroutine init_abslayer(self)
        class(t_propagator) :: self

        call cpml%init(o_kpa_max=1.)

    end subroutine
    
    !========= forward modeling =================
    !WE: M du_dt = Du + s
    !u=[vx vy vz p]^T, p=(sxx+syy+szz)/3=sxx+syy+szz, s=[fvx fvy fvz fp]^T*delta(x-xs)
    !
    !  [rho   0  ]    [0   0   0   dx]
    !M=[ 0  1/kpa], D=|0   0   0   dy|
    !                 |0   0   0   dz|
    !                 [dx  dy  dz  0 ]
    !
    !Discretization (staggered grid in space and time):
    !  (grid index)     (real index)
    !  s(iz,ix,iy):=   s[iz,ix,iy]^it+0.5,it+1.5,...
    ! vz(iz,ix,iy):=  vy[iz-0.5,ix,iy]^it,it+1,...
    ! vx(iz,ix,iy):=  vx[iz,ix-0.5,iy]^it,it+1,...
    ! vy(iz,ix,iy):=  vy[iz,ix,iy-0.5]^it,it+1,...

    !========= adjoint propagation =================
    !WEq:      Mdu_dt   = Du   + src
    !Adjoint:  Mdv_dt^T = D^Tv + adjsrc
    !
    !Discrete form (staggered grid in space and time, 2D as example):
    ! M  [ [vx^it+1  ] [vx^it    ] ]   [ 0     0    dx(-)][vx^it+1  ]
    !----| |vz^it+1  |-|vz^it    | | = | 0     0    dz(-)||vz^it+1  |  +src
    ! dt [ [ s^it+1.5] [ s^it+0.5] ]   [dx(+) dz(+)  0   ][ s^it+0.5]
    !dx(-):=c1(s(ix)-s(ixm1))+c2(s(ixp1)-s(ixm2))  (O(x4))
    !dx(+):=c1(v(ixp1)-v(ix))+c2(v(ixp2)-v(ixm1))  (O(x4))
    !
    !Adjoint:
    ! M  [ [vx^it    ] [vx^it+1  ] ]   [ 0       0      dx(+)^T][vx^it+1  ]
    !----| |vz^it    |-|vz^it+1  | | = | 0       0      dz(+)^T||vz^it+1  |  +src
    ! dt [ [ s^it+0.5] [ s^it+1.5] ]   [dx(-)^T dz(-)^T  0     ][ s^it+0.5]
    !dx(+)^T=c1(s(ixm1)-s(ix))+c2(s(ixm2)-s(ixp1))=-dx(-)
    !dx(-)^T=c1(v(ix)-v(ixp1))+c2(v(ixm1)-v(ixp2))=-dx(+)
    !i.e. D^T=-D, allowing to use same code with a negated sign for dt.
    !
    !transposes of time derivatives centered on v^it+0.5 and s^it+1
    !so v^it+1   - v^it     => v^it - v^it+1
    !   s^it+1.5 - s^it+0.5 => s^it+0.5 - s^it+1.5
    !
    !In each time step:
    !d=RGAsrc, src:source term, A:put source, G:propagator, R:sampling at receivers
    !d=[vx vy vz p]^T, R=[I  0   0 ], A=[diag3(b) 0], src=[fx fy fz p]^T
    !                    |0 2/3 1/3|    [   0     1]
    !                    [   0    1]
    !Adjoint: src=A^T N G^T M R^T d
    !M R^T:put adjoint sources, 
    !G^T:  adjoint propagator,
    !A^T N:get adjoint fields

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

    subroutine forward_scattering(self,fld_du,fld_u,W2imag)
        class(t_propagator) :: self
        type(t_field) :: fld_du,fld_u
        real,dimension(cb%mz,cb%mx,cb%my) :: W2imag

        real,parameter :: time_dir=1. !time direction

        !seismo
        call alloc(fld_du%seismo,shot%nrcv,self%nt)
        
        tt1=0.; tt2=0.; tt3=0.; tt4=0.; tt5=0.; tt6=0.
        
        ift=1; ilt=self%nt

        do it=ift,ilt
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_du%check_value
                 call fld_u%check_value
            endif

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            ! !step 1: add forces to v^it
            ! call cpu_time(tic)
            ! call self%inject_velocities_ext(fld_mu,-fld_a)
            ! call cpu_time(toc)
            ! tt1=tt1+toc-tic

            !step 2: from v^it to v^it+1 by differences of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities(fld_u, time_dir,it)
            call self%update_velocities(fld_du,time_dir,it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            !step 3: add pressure to s^it+0.5
            call cpu_time(tic)
            call self%inject_stresses(fld_u,time_dir,it)
            call self%inject_stresses_scattering(fld_du,fld_u%p(1:cb%mz,1:cb%mx,1:cb%my)*W2imag,time_dir)
            call cpu_time(toc)
            tt3=tt3+toc-tic

            !step 4: from s^it+0.5 to s^it+1.5 by differences of v^it+1
            call cpu_time(tic)
            call self%update_stresses(fld_u, time_dir,it)
            call self%update_stresses(fld_du,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !step 5: sample v^it+1 or s^it+1.5 at receivers
            call cpu_time(tic)
            call self%extract(fld_du,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !snapshot
            ! call fld_u%write(it,o_suffix='_for2')
            call fld_du%write(it)

            !step 6: save v^it+1 in boundary layers
            ! if(fld_u%if_will_reconstruct) then
                call cpu_time(tic)
                call  fld_u%boundary_transport('save',it)
                call fld_du%boundary_transport('save',it)
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

    end subroutine

    subroutine adjoint(self,fld_a,fld_u,oif_record_adjseismo,oif_compute_imag,oif_compute_grad)
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
            call self%inject_stresses_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_stresses_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: sfield%s_dt^it+0.5 \dot rfield%s^it+0.5
            !use sfield%v^it+1 to compute sfield%s_dt^it+0.5, as backward step 4
            if(if_compute_grad.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call gradient_moduli(fld_u,fld_a,it,cb%grad(:,:,:,2))
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif

            if(if_compute_imag.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call imaging(fld_u,fld_a,it,cb%imag)
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

            !adjoint step 5: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_velocities_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            ! !adjoint step 1: sample v^it or s^it+0.5 at source position
            ! if(if_record_adjseismo) then
            !     call cpu_time(tic)
            !     call self%extract_adjoint(fld_a,it)
            !     call cpu_time(toc)
            !     tt11=tt11+toc-tic
            ! endif
            
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
                call fld_a%write_ext(it,'grad_u_dot_a' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
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

    subroutine adjoint_du_star_a(self,fld_a,fld_du,fld_u,W2imag,oif_compute_imag,oif_compute_grad)
        class(t_propagator) :: self
        type(t_field) :: fld_a,fld_du,fld_u
        real,dimension(cb%mz,cb%mx,cb%my) :: W2imag
        logical,optional :: oif_compute_imag,oif_compute_grad
        
        real,parameter :: time_dir=-1. !time direction
        logical :: if_compute_imag,if_compute_grad

        ! if_record_adjseismo =either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))
        if_compute_imag     =either(oif_compute_imag,    .false.,present(oif_compute_imag))
        if_compute_grad     =either(oif_compute_grad,    .false.,present(oif_compute_grad))
        ! self%if_compute_engy=self%if_compute_engy.and.(if_compute_imag.or.if_compute_grad)
        
        ! if(if_record_adjseismo)  call alloc(fld_a%seismo,1,self%nt)
        if(if_compute_imag)      call alloc(cb%imag,cb%mz,cb%mx,cb%my,self%nimag)
        if(if_compute_grad)      call alloc(cb%grad,cb%mz,cb%mx,cb%my,self%ngrad)
        ! if(self%if_compute_engy) call alloc(cb%engy,cb%mz,cb%mx,cb%my,self%nengy)

        !reinitialize absorbing boundary for incident wavefield reconstruction
        ! if(present(o_sf)) then
                                          fld_du%dvz_dz=0.
                                           fld_u%dvz_dz=0.
                                          fld_du%dvx_dx=0.
                                           fld_u%dvx_dx=0.
            ! if(allocated(fld_du%dvy_dy))  fld_du%dvy_dy=0.
            ! if(allocated( fld_u%dvy_dy))   fld_u%dvy_dy=0.
                                          fld_du%dp_dz=0.
                                           fld_u%dp_dz=0.
                                          fld_du%dp_dx=0.
                                           fld_u%dp_dx=0.
            ! if(allocated(fld_du%dp_dy))   fld_du%dp_dy=0.
            ! if(allocated( fld_u%dp_dy))    fld_u%dp_dy=0.
        ! endif
                    
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
                call fld_du%check_value
                call fld_u%check_value
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then

                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_du%boundary_transport('load',it)
                call  fld_u%boundary_transport('load',it)
                call cpu_time(toc)
                tt1=tt1+toc-tic
                
                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%update_stresses(fld_du,time_dir,it)
                call self%update_stresses( fld_u,time_dir,it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_stresses_scattering(fld_du,fld_u%p(1:cb%mz,1:cb%mx,1:cb%my)*W2imag,time_dir)
                call self%inject_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_stresses_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_stresses_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: sfield%s_dt^it+0.5 \dot rfield%s^it+0.5
            !use sfield%v^it+1 to compute sfield%s_dt^it+0.5, as backward step 4
            if(if_compute_grad.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call gradient_moduli(fld_du,fld_a,it,cb%grad(:,:,:,2))
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif

            if(if_compute_imag.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call imaging(fld_du,fld_a,it,cb%imag)
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
                call self%update_velocities(fld_du,time_dir,it)
                call self%update_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt7=tt7+toc-tic

                !backward step 1: rm forces from v^it
                call cpu_time(tic)
                ! call self%inject_velocities_ext(fld_du,fld_u,W2imag)
                call self%inject_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt8=tt8+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_velocities_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            ! !adjoint step 1: sample v^it or s^it+0.5 at source position
            ! if(if_record_adjseismo) then
            !     call cpu_time(tic)
            !     call self%extract_adjoint(fld_a,it)
            !     call cpu_time(toc)
            !     tt11=tt11+toc-tic
            ! endif
            
            ! !grho: sfield%v_dt^it \dot rfield%v^it
            ! !use sfield%s^it+0.5 to compute sfield%v_dt^it, as backward step 2
            ! if(if_compute_grad.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call gradient_density(o_sf,rf,it,cb%grad(:,:,:,1))
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
            
            !snapshot
            ! call fld_a%write(it,o_suffix='_back2')
            call fld_du%write(it,o_suffix='_rev')
            ! call fld_u%write(it,o_suffix='_back2')
            if(if_compute_imag) then
                call fld_a%write_ext(it,'imag_du_dot_a' ,cb%imag,size(cb%imag))
            endif
            if(if_compute_grad) then
                ! call fld_a%write_ext(it,'grad_density',cb%grad(:,:,:,1),size(cb%grad(:,:,:,1)))
                ! call fld_a%write_ext(it,'grad_moduli' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
                call fld_a%write_ext(it,'grad_du_dot_a' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
            endif
            ! if(self%if_compute_engy) then
            !     call fld_a%write_ext(it,'engy',cb%engy(:,:,:,1),size(cb%engy(:,:,:,1)))
            ! endif

        enddo
        
        !postprocess gradient
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

    subroutine adjoint_scattering_u_star_da(self,fld_da,fld_a,W2imag,fld_u,oif_record_adjseismo,oif_compute_grad)
        class(t_propagator) :: self
        type(t_field) :: fld_da,fld_a,fld_u
        real,dimension(cb%mz,cb%mx,cb%my) :: W2imag
        logical,optional :: oif_record_adjseismo, oif_compute_grad
        
        real,parameter :: time_dir=-1. !time direction
        logical :: if_record_adjseismo, if_compute_imag, if_compute_grad

        if_record_adjseismo =either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))
        ! if_compute_imag     =either(oif_compute_imag,    .false.,present(oif_compute_imag))
        if_compute_grad     =either(oif_compute_grad,    .false.,present(oif_compute_grad))
        ! self%if_compute_engy=self%if_compute_engy.and.(if_compute_imag.or.if_compute_grad)
        
        if(if_record_adjseismo)  call alloc(fld_da%seismo,1,self%nt)
        ! if(if_compute_imag)      call alloc(cb%imag,cb%mz,cb%mx,cb%my,self%nimag)
        if(if_compute_grad)      call alloc(cb%grad,cb%mz,cb%mx,cb%my,self%ngrad)
        ! if(self%if_compute_engy) call alloc(cb%engy,cb%mz,cb%mx,cb%my,self%nengy)

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
        
! call alloc(sf_p_save,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        do it=ilt,ift,int(time_dir)
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_da%check_value
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
            call self%inject_stresses_scattering(fld_da,fld_a%p(1:cb%mz,1:cb%mx,1:cb%my)*W2imag,time_dir)
            call self%inject_stresses_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_stresses_adjoint(fld_da,time_dir,it)
            call self%update_stresses_adjoint(fld_a, time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: sfield%s_dt^it+0.5 \dot rfield%s^it+0.5
            !use sfield%v^it+1 to compute sfield%s_dt^it+0.5, as backward step 4
            if(if_compute_grad.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call gradient_moduli(fld_da,fld_u,it,cb%grad(:,:,:,2))
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif

            ! if(if_compute_imag.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call imaging(fld_u,fld_a,it,cb%imag)
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

            !adjoint step 5: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_velocities_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities_adjoint(fld_da,time_dir,it)
            call self%update_velocities_adjoint(fld_a, time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            !adjoint step 1: sample v^it or s^it+0.5 at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract_adjoint(fld_da,it)
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
            call fld_da%write(it,o_suffix='_rev')
            ! call fld_a%write(it,o_suffix='_back3')
            ! call fld_u%write(it,o_suffix='_back3')
            ! if(if_compute_imag) then
                ! call fld_lda%write_ext(it,'imag' ,cb%imag,size(cb%imag))
            ! endif
            if(if_compute_grad) then
                ! call fld_da%write_ext(it,'grad_density',cb%grad(:,:,:,1),size(cb%grad(:,:,:,1)))
                ! call fld_da%write_ext(it,'grad_moduli' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
                call fld_da%write_ext(it,'grad_u_dot_da' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
            endif
            ! if(self%if_compute_engy) then
                ! call fld_a%write_ext(it,'engy',cb%engy(:,:,:,1),size(cb%engy(:,:,:,1)))
            ! endif

        enddo
        
        !postprocess gradient
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

    end subroutine

    subroutine adjoint_scattering_du_star_da(self,fld_da,fld_a,fld_du,fld_u,W2imag,oif_record_adjseismo,oif_compute_grad)
        class(t_propagator) :: self
        type(t_field) :: fld_da,fld_a
        type(t_field) :: fld_du,fld_u
        real,dimension(cb%mz,cb%mx,cb%my) :: W2imag
        logical,optional :: oif_record_adjseismo, oif_compute_grad
        
        real,parameter :: time_dir=-1. !time direction
        logical :: if_record_adjseismo, if_compute_imag, if_compute_grad

        ! if_record_adjseismo =either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))
        ! if_compute_imag     =either(oif_compute_imag,    .false.,present(oif_compute_imag))
        if_compute_grad     =either(oif_compute_grad,    .false.,present(oif_compute_grad))
        ! self%if_compute_engy=self%if_compute_engy.and.(if_compute_imag.or.if_compute_grad)
        
        ! if(if_record_adjseismo)  call alloc(fld_da%seismo,1,self%nt)
        ! if(if_compute_imag)      call alloc(cb%imag,cb%mz,cb%mx,cb%my,self%nimag)
        if(if_compute_grad)      call alloc(cb%grad,cb%mz,cb%mx,cb%my,self%ngrad)
        ! if(self%if_compute_engy) call alloc(cb%engy,cb%mz,cb%mx,cb%my,self%nengy)

        !reinitialize absorbing boundary for incident wavefield reconstruction
        ! if(present(o_sf)) then
                                          fld_du%dvz_dz=0.
                                           fld_u%dvz_dz=0.
                                          fld_du%dvx_dx=0.
                                           fld_u%dvx_dx=0.
            ! if(allocated(fld_du%dvy_dy))  fld_du%dvy_dy=0.
            ! if(allocated( fld_u%dvy_dy))   fld_u%dvy_dy=0.
                                          fld_du%dp_dz=0.
                                           fld_u%dp_dz=0.
                                          fld_du%dp_dx=0.
                                           fld_u%dp_dx=0.
        ! endif

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
                call fld_da%check_value
                call fld_a%check_value
                call fld_du%check_value
                call fld_u%check_value
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then
                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_du%boundary_transport('load',it)
                call  fld_u%boundary_transport('load',it)
                call cpu_time(toc)
                tt1=tt1+toc-tic
                
                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%update_stresses(fld_du,time_dir,it)
                call self%update_stresses( fld_u,time_dir,it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_stresses_scattering(fld_du,fld_u%p(1:cb%mz,1:cb%mx,1:cb%my)*W2imag,time_dir)
                call self%inject_stresses( fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_stresses_scattering(fld_da,fld_a%p(1:cb%mz,1:cb%mx,1:cb%my)*W2imag,time_dir)
            call self%inject_stresses_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_stresses_adjoint(fld_da,time_dir,it)
            call self%update_stresses_adjoint(fld_a, time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: sfield%s_dt^it+0.5 \dot rfield%s^it+0.5
            !use sfield%v^it+1 to compute sfield%s_dt^it+0.5, as backward step 4
            if(if_compute_grad.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call gradient_moduli(fld_da,fld_du,it,cb%grad(:,:,:,2))
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif

            ! if(if_compute_imag.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call imaging(fld_u,fld_a,it,cb%imag)
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
                call self%update_velocities(fld_du,time_dir,it)
                call self%update_velocities( fld_u,time_dir,it)
                call cpu_time(toc)
                tt7=tt7+toc-tic

                !backward step 1: rm forces from v^it
                call cpu_time(tic)
                ! call self%inject_velocities(fld_du,fld_u,it)
                call self%inject_velocities( fld_u,time_dir,it)
                call cpu_time(toc)
                tt8=tt8+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_velocities_adjoint(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities_adjoint(fld_da,time_dir,it)
            call self%update_velocities_adjoint(fld_a, time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            !adjoint step 1: sample v^it or s^it+0.5 at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract_adjoint(fld_da,it)
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
            call fld_da%write(it,o_suffix='_rev')
            ! call fld_a%write(it,o_suffix='_back3')
            ! call fld_u%write(it,o_suffix='_back3')
            ! if(if_compute_imag) then
                ! call fld_lda%write_ext(it,'imag' ,cb%imag,size(cb%imag))
            ! endif
            if(if_compute_grad) then
                ! call fld_da%write_ext(it,'grad_density',cb%grad(:,:,:,1),size(cb%grad(:,:,:,1)))
                ! call fld_da%write_ext(it,'grad_moduli' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
                call fld_da%write_ext(it,'grad_du_dot_da' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
            endif
            ! if(self%if_compute_engy) then
                ! call fld_a%write_ext(it,'engy',cb%engy(:,:,:,1),size(cb%engy(:,:,:,1)))
            ! endif

        enddo
        
        !postprocess gradient
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

    end subroutine

    !add RHS to v^it
    subroutine inject_velocities(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
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
        
    end subroutine

    ! subroutine inject_velocities_ext(self,f,RHS)
    !     class(t_propagator) :: self
    !     type(t_field) :: f, RHS

    !     f%vz = f%vz + RHS%vz
    !     f%vx = f%vx + RHS%vx
    !     f%vy = f%vy + RHS%vy
    ! end subroutine
    
    !v^it -> v^it+1 by FD of s^it+0.5
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
    
    !add RHS to s^it+0.5
    subroutine inject_stresses(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        ifz=shot%src%ifz-cb%ioz+1; iz=shot%src%iz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
        ifx=shot%src%ifx-cb%iox+1; ix=shot%src%ix-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
        ify=shot%src%ify-cb%ioy+1; iy=shot%src%iy-cb%ioy+1; ily=shot%src%ily-cb%ioy+1
        
        wl=time_dir*f%wavelet(1,it)
        
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
        
    end subroutine
    
    subroutine inject_stresses_scattering(self,f,RHS,time_dir)
        class(t_propagator) :: self
        type(t_field) :: f
        real,dimension(cb%mz,cb%mx,cb%my) :: RHS
        real :: time_dir

        f%p(1:cb%mz,1:cb%mx,1:cb%my) = f%p(1:cb%mz,1:cb%mx,1:cb%my) + time_dir*RHS*self%kpa(1:cb%mz,1:cb%mx,1:cb%my)
    end subroutine

    !s^it+0.5 -> s^it+1.5 by FD of v^it+1
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

    subroutine extract(self,f,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
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
        
    end subroutine

    !add RHS to s^it+1.5
    subroutine inject_stresses_adjoint(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        do i=1,shot%nrcv

            if(shot%rcv(i)%comp=='p') then

                ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                ify=shot%rcv(i)%ify-cb%ioy+1; iy=shot%rcv(i)%iy-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1
                
                !adjsource for pressure
                wl=f%wavelet(i,it)
                
                if(if_hicks) then 

                    f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) +wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef

                else           
                    !p[iz,ix,iy]
                    f%p(iz,ix,iy) = f%p(iz,ix,iy) +wl*self%kpa(iz,ix,iy)

                endif

            endif

        enddo
        
    end subroutine
    
    !s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
    subroutine update_stresses_adjoint(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        !same code with -dt
        call self%update_stresses(f,time_dir,it)
        
    end subroutine
    
    !add RHS to v^it+1
    subroutine inject_velocities_adjoint(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        do i=1,shot%nrcv
            ifz=shot%rcv(i)%ifz-cb%ioz+1; iz=shot%rcv(i)%iz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
            ifx=shot%rcv(i)%ifx-cb%iox+1; ix=shot%rcv(i)%ix-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
            ify=shot%rcv(i)%ify-cb%ioy+1; iy=shot%rcv(i)%iy-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1
            
            wl=f%wavelet(i,it)
            
            if(if_hicks) then
                select case (shot%rcv(i)%comp)
                    case ('vz') !vertical z adjsource
                    f%vz(ifz:ilz,ifx:ilx,ify:ily) = f%vz(ifz:ilz,ifx:ilx,ify:ily) + wl*self%buoz(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef

                    case ('vx') !horizontal x adjsource
                    f%vx(ifz:ilz,ifx:ilx,ify:ily) = f%vx(ifz:ilz,ifx:ilx,ify:ily) + wl*self%buox(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef
                    
                    case ('vy') !horizontal y adjsource
                    f%vy(ifz:ilz,ifx:ilx,ify:ily) = f%vy(ifz:ilz,ifx:ilx,ify:ily) + wl*self%buoy(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef
                    
                end select
                
            else
                select case (shot%rcv(i)%comp)
                    case ('vz') !vertical z adjsource
                    !vz[ix,iy,iz-0.5]
                    f%vz(iz,ix,iy) = f%vz(iz,ix,iy) + wl*self%buoz(iz,ix,iy)

                    case ('vx') !horizontal x adjsource
                    !vx[ix-0.5,iy,iz]
                    f%vx(iz,ix,iy) = f%vx(iz,ix,iy) + wl*self%buox(iz,ix,iy)
                    
                    case ('vy') !horizontal y adjsource
                    !vy[ix,iy-0.5,iz]
                    f%vy(iz,ix,iy) = f%vy(iz,ix,iy) + wl*self%buoy(iz,ix,iy)
                    
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    !v^it+1 -> v^it by FD^T of s^it+0.5
    subroutine update_velocities_adjoint(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
        !same code with -dt
        call self%update_velocities(f,time_dir,it)
        
    end subroutine
    
    subroutine extract_adjoint(self,f,it)
        class(t_propagator) :: self
        type(t_field) :: f
        
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
        call dealloc(self%buox, self%buoy, self%buoz, self%kpa, self%inv_kpa)
    end subroutine

    !========= for wavefield correlation ===================   
    !grho = -adjv · dv_dt
    !     = -adjv · b ∇p
    !
    !gkpa = -adjp · -1/kpa2 dp_dt
    !     = -adjp · -1/kpa  ∇·v
    !
    !v^it+1, p^it+0.5, adjp^it+0.5
    !p^it+0.5, v^it, adjv^it
    !use (v[i+1]+v[i])/2 to approximate v[i+0.5], so is adjv
    
    subroutine gradient_density(sf,rf,it,grad)
        type(t_field), intent(in) :: sf, rf
        real,dimension(cb%mz,cb%mx,cb%my) :: grad
        
        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        if(m%is_cubic) then
            call grad3d_density(sf%p,rf%vz,rf%vx,rf%vy,&
                                grad,                  &
                                ifz,ilz,ifx,ilx,ify,ily)
        else
            call grad2d_density(sf%p,rf%vz,rf%vx,&
                                grad,            &
                                ifz,ilz,ifx,ilx)
        endif
        
    end subroutine

    subroutine gradient_moduli(sf,rf,it,grad)
        type(t_field), intent(in) :: sf, rf
        real,dimension(cb%mz,cb%mx,cb%my) :: grad
        
        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        if(m%is_cubic) then
            call grad3d_moduli(sf%vz,sf%vx,sf%vy,rf%p,&
                               grad,                  &
                               ifz,ilz,ifx,ilx,ify,ily)
        else
            ! call grad2d_moduli(sf%p,sf_p_save,rf%p,&
            !                    grad,            &
            !                    ifz,ilz,ifx,ilx)
            ! sf_p_save = sf%p
            
            !inexact greadient
            call grad2d_moduli(sf%vz,sf%vx,rf%p,&
                               grad,            &
                               ifz,ilz,ifx,ilx)
        endif
        
    end subroutine
    
    subroutine gradient_postprocess
        !grho
        cb%grad(:,:,:,1)=-cb%grad(:,:,:,1) / cb%rho(1:cb%mz,1:cb%mx,1:cb%my)         *m%cell_volume*rdt !scaled by m%cell_volume*rdt tobe a gradient in the discretized world
        !gkpa
        cb%grad(:,:,:,2)=-cb%grad(:,:,:,2) * (-ppg%inv_kpa(1:cb%mz,1:cb%mx,1:cb%my)) *m%cell_volume*rdt
                
        !preparing for cb%project_back
        cb%grad(1,:,:,:) = cb%grad(2,:,:,:)

    end subroutine

    subroutine imaging(sf,rf,it,imag)
        type(t_field), intent(in) :: sf, rf
        real,dimension(cb%mz,cb%mx,cb%my) :: imag
        
        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)
        
        ! if(m%is_cubic) then
        !     call imag3d_xcorr(sf%p,rf%p,&
        !                       imag,                  &
        !                       ifz,ilz,ifx,ilx,ify,ily)
        ! else
            call imag2d_xcorr(sf%p,rf%p,&
                              imag,            &
                              ifz,ilz,ifx,ilx)
        ! endif

        ! call imag2d_xcorr(sf%p,sf%vz,sf%vx,&
        !                   rf%p,rf%vz,rf%vx,&
        !                   imag,            &
        !                   ifz,ilz,ifx,ilx)

    end subroutine

    subroutine imaging_postprocess
        cb%imag = cb%imag * rdt

        !for cb%project_back
        cb%imag(1,:,:,:) = cb%imag(2,:,:,:)
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

    
    subroutine grad3d_moduli(sf_vz,sf_vx,sf_vy,rf_p,&
                             grad,                  &
                             ifz,ilz,ifx,ilx,ify,ily)
        real,dimension(*) :: sf_vz,sf_vx,sf_vy,rf_p
        real,dimension(*) :: grad
        
        nz=cb%nz
        nx=cb%nx
        
        dvz_dz=0.
        dvx_dx=0.
        dvy_dx=0.
        
        dsp=0.
         rp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,j,&
        !$omp         izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
        !$omp         iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
        !$omp         iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
        !$omp         dvz_dz,dvx_dx,dvy_dy,&
        !$omp         dsp,rp)
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
                
                dsp = dvz_dz +dvx_dx +dvy_dy
                 rp = rf_p(i)
                
                grad(j)=grad(j) - dsp*rp
                
            end do
            
        end do
        end do
        !$omp end do
        !$omp end parallel
        
    end subroutine

    ! subroutine grad2d_moduli(sf_p,sf_p_save,rf_p,&
    !                          grad,            &
    !                          ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: sf_p,sf_p_save,rf_p
    !     real,dimension(*) :: grad
        
    !     nz=cb%nz
        
    !     dsp=0.
    !      rp=0.
        
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
                
    !             dsp = sf_p_save(i) - sf_p(i)
    !              rp = rf_p(i)
                
    !             grad(j)=grad(j) - dsp*rp
                
    !         end do
            
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

    subroutine grad2d_moduli(sf_vz,sf_vx,rf_p,&
                             grad,            &
                             ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_vz,sf_vx,rf_p
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        dvz_dz=0.
        dvx_dx=0.
        
        dsp=0.
         rp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dvz_dz,dvx_dx,&
        !$omp         dsp,rp)
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
                
                dsp = dvz_dz +dvx_dx
                 rp = rf_p(i)
                
                grad(j)=grad(j) - dsp*rp
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine
    
    subroutine grad3d_density(sf_p,rf_vz,rf_vx,rf_vy,&
                              grad,                  &
                              ifz,ilz,ifx,ilx,ify,ily)
        real,dimension(*) :: sf_p,rf_vz,rf_vx,rf_vy
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
        !$omp         dsvz,dsvx,dsvy,&
        !$omp          rvz, rvx, rvy)
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

                dsvz = (c1z*(sf_p(iz_ix_iy  )-sf_p(izm1_ix_iy)) +c2z*(sf_p(izp1_ix_iy)-sf_p(izm2_ix_iy))) &
                      +(c1z*(sf_p(izp1_ix_iy)-sf_p(iz_ix_iy  )) +c2z*(sf_p(izp2_ix_iy)-sf_p(izm1_ix_iy)))
                dsvx = (c1x*(sf_p(iz_ix_iy  )-sf_p(iz_ixm1_iy)) +c2x*(sf_p(iz_ixp1_iy)-sf_p(iz_ixm2_iy))) &
                      +(c1x*(sf_p(iz_ixp1_iy)-sf_p(iz_ix_iy  )) +c2x*(sf_p(iz_ixp2_iy)-sf_p(iz_ixm1_iy)))
                dsvy = (c1y*(sf_p(iz_ix_iy  )-sf_p(iz_ix_iym1)) +c2y*(sf_p(iz_ix_iyp1)-sf_p(iz_ix_iym2))) &
                      +(c1y*(sf_p(iz_ix_iyp1)-sf_p(iz_ix_iy  )) +c2y*(sf_p(iz_ix_iyp2)-sf_p(iz_ix_iym1)))
                
                !complete equation with unnecessary terms e.g. sf_p(iz_ix) for better understanding
                !with flag -O, the compiler should automatically detect such possibilities of simplification
                
                rvz = rf_vz(izp1_ix_iy) +rf_vz(iz_ix_iy)
                rvx = rf_vx(iz_ixp1_iy) +rf_vx(iz_ix_iy)
                rvy = rf_vy(iz_ix_iyp1) +rf_vy(iz_ix_iy)
                
                grad(j)=grad(j) - 0.25*( dsvz*rvz + dsvx*rvx + dsvy*rvy )
                
            enddo
            
        enddo
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine grad2d_density(sf_p,rf_vz,rf_vx,&
                              grad,            &
                              ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_p,rf_vz,rf_vx
        real,dimension(*) :: grad
        
        nz=cb%nz
        
        dsvz=0.; dsvx=0.
         rvz=0.; rvx=0.

        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dsvz,dsvx,&
        !$omp          rvz, rvx)
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
                
                dsvz = (c1z*(sf_p(iz_ix  )-sf_p(izm1_ix)) +c2z*(sf_p(izp1_ix)-sf_p(izm2_ix))) &
                      +(c1z*(sf_p(izp1_ix)-sf_p(iz_ix  )) +c2z*(sf_p(izp2_ix)-sf_p(izm1_ix)))
                dsvx = (c1x*(sf_p(iz_ix  )-sf_p(iz_ixm1)) +c2x*(sf_p(iz_ixp1)-sf_p(iz_ixm2))) &
                      +(c1x*(sf_p(iz_ixp1)-sf_p(iz_ix  )) +c2x*(sf_p(iz_ixp2)-sf_p(iz_ixm1)))
                !complete equation with unnecessary terms e.g. sf_p(iz_ix) for better understanding
                !with flag -Ox, the compiler should automatically detect such possible simplification
                
                rvz = rf_vz(iz_ix) +rf_vz(izp1_ix)
                rvx = rf_vx(iz_ix) +rf_vx(iz_ixp1)
                
                grad(j)=grad(j) - 0.25*( dsvz*rvz + dsvx*rvx )
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine


    subroutine imag3d_xcorr(sf_p,rf_p,             &
                            imag,                  &
                            ifz,ilz,ifx,ilx,ify,ily)
        real,dimension(*) :: sf_p,rf_p
        real,dimension(*) :: imag
        
        nz=cb%nz
        nx=cb%nx
        
        sp=0.
        rp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,j,&
        !$omp         sp,rp)
        !$omp do schedule(dynamic) collapse(2)
        do iy=ify,ily
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(iy-cb%ify)*cb%nz*cb%nx+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+(iy-1)     *cb%mz*cb%mx+1 !grad has no boundary layers
                
                sp = sf_p(i)
                rp = rf_p(i)
                
                imag(j)=imag(j) + sp*rp
                
            end do
            
        end do
        end do
        !$omp end do
        !$omp end parallel
        
    end subroutine

    ! subroutine imag2d_xcorr(sf_p,sf_vz,sf_vx,&
    !                         rf_p,rf_vz,rf_vx,&
    !                         imag,          &
    !                         ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: sf_p,sf_vz,sf_vx
    !     real,dimension(*) :: rf_p,rf_vz,rf_vx
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
    !                 sf_p(i)*rf_p(i) + sf_vz(i)*rf_vz(i) + sf_vx(i)*rf_vx(i)
                
    !         end do
            
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

    subroutine imag2d_xcorr(sf_p,rf_p,     &
                            imag,          &
                            ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_p,rf_p
        real,dimension(*) :: imag
        
        nz=cb%nz
        
        sp=0.
        rp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         sp,rp)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !grad has no boundary layers
                
                sp = sf_p(i)
                rp = rf_p(i)
                
                imag(j)=imag(j) + sp*rp
                
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
