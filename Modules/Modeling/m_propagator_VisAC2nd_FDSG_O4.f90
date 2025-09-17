module m_propagator
use m_System
use m_hicks, only : hicks_r
use m_hilbert
use m_resampler
use m_model
use m_shot
use m_computebox
use m_field
use m_correlate
use m_cpml
use singleton

    private

    !FD coef
    real,dimension(2),parameter :: coef = [9./8.,-1./24.] !Fornberg, 1988, Generation of Finite-Difference Formulas on Arbitrary Spaced Grids.
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z

    !ViscoAcoustic coef
    ! real :: visac_a = 5.7356
    ! real :: visac_b = -762.1606
    ! real :: visac_d = 46054.0
    ! real :: visac_e = 1

    real :: visac_a = 3.7801
    real :: visac_b = -110.687
    real :: visac_d = 1111.5
    real :: visac_e = 1

    !local const
    real :: dt2, inv_2dt, inv_2dz, inv_2dx

    !scaling source wavelet
    real :: wavelet_scaler

    type,public :: t_propagator
        !info
        character(i_str_xxlen) :: info = &
            'Time-domain ISOtropic 2D/3D VIScoACoustic propagation'//s_NL// &
            '2nd-order Pressure formulation'//s_NL// &
            'Vireux-Levandar Staggered-Grid Finite-Difference (FDSG) method'//s_NL// &
            'Cartesian O(x⁴,t²) stencil'//s_NL// &
            'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
            '   -> dt ≤ 0.5*Vmax/dx'//s_NL// &
            'Required model attributes: vp, rho, qp'//s_NL// &
            'Required field components: p, p_prev, p_next'//s_NL// &
            'Required boundary layer thickness: 2'//s_NL// &
            'Poynting definitions: Esq_gradphi'//s_NL// &
            'Imaging conditions: ipp ibksc ifwsc (P-Pxcorr of backward & forward scattering)'//s_NL// &
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: gbuo(wait), gikpa, gqp'

        integer :: nbndlayer=max(1,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=2 !number of basic gradients
        integer :: nimag=3 !number of basic images
        !integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:,:),allocatable :: buoz, buox, buoy, kpa!, qp
        complex,dimension(:,:,:),allocatable :: invC2,  C1n,  C0n
        complex,dimension(:,:,:),allocatable :: invC2H, C1nH, C0nH

        !time frames
        integer :: nt
        real :: dt

        !!reference value
        !real :: invsqEref

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
        
        procedure :: inject_pressure
        procedure :: update_pressure
        procedure :: evolve_pressure
        procedure :: extract
        ! procedure :: gaussian_smooth

        ! procedure :: dft2d
        ! procedure :: idft2d
        ! procedure :: fftshift2d
        ! procedure :: ifftshift2d
        ! procedure :: fft_filter

        final :: final

    end type

    type(t_propagator),public :: ppg

    character(:),allocatable :: s_poynting_def

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    logical :: if_record_adjseismo=.true.

    ! logical :: is_absolute_virtual

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

        if(index(self%info,'qp')>0 .and. .not. allocated(m%qp)) then
            call alloc(m%qp,m%nz,m%nx,m%ny,o_init=1000.)
            call warn('Constant Qp model (1000) is allocated by propagator.')
        endif

        ! if(.not.setup%get_bool('IS_Q_DISPERSION',o_default='T')) then
        !     visac_a = 0.
        !     visac_b = 0.
        !     visac_d = 0.
        ! endif
        ! if(setup%get_bool('IS_Q_DISSIPATION',o_default='T')) then
        !     visac_e = 0.
        ! endif
        if(.not.setup%get_bool('IS_Q_ATTENUATION',o_default='T')) then
            visac_a = 0.
            visac_b = 0.
            visac_d = 0.
            visac_e = 0.
        endif

        call hud('visac coefs:'//num2str(visac_a)//', '//num2str(visac_b)//', '//num2str(visac_d)//', '//num2str(visac_e))
                
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
        complex,dimension(:,:,:),allocatable :: C2,C1,C0

        c1x=coef(1)/m%dx; c1y=coef(1)/m%dy; c1z=coef(1)/m%dz
        c2x=coef(2)/m%dx; c2y=coef(2)/m%dy; c2z=coef(2)/m%dz

        dt2=self%dt**2
        inv_2dt =1./2/self%dt
        inv_2dz =1./2/m%dz
        inv_2dx =1./2/m%dx
        
        wavelet_scaler=dt2/m%cell_volume

        if_hicks=shot%if_hicks

        call alloc(self%buoz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%buox,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%buoy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(self%kpa, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        !call alloc(self%qp,  [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        self%kpa=cb%rho*(cb%vp)**2
        !self%qp =cb%qp

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

        ! call alloc(self%r2,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        ! self%r2 = (cb%vp*self%dt/m%dx)**2

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

        ! s_poynting_def=setup%get_str('POYNTING_DEF',o_default='Esq_gradphi')
        ! if(s_poynting_def/='Esq_gradphi') call error('Sorry, other Poynting definitions have not yet implemented.')

        !coef
        C2 = 1 -2*visac_a/r_pi/cb%qp -c_i*visac_e/cb%qp; self%invC2 = 1/C2
        C1 =    2*visac_b/r_pi/cb%qp;                    self%C1n = C1*self%invC2
        C0 =    2*visac_d/r_pi/cb%qp;                    self%C0n = C0*self%invC2

        C2 = 1 -2*visac_a/r_pi/cb%qp +c_i*visac_e/cb%qp; self%invC2H = 1/C2
        C1 =    2*visac_b/r_pi/cb%qp;                    self%C1nH= C1*self%invC2H
        C0 =    2*visac_d/r_pi/cb%qp;                    self%C0nH= C0*self%invC2H

        deallocate(C2,C1,C0) !save some RAM

    end subroutine

    subroutine init_field(self,f,name,ois_simple,ois_adjoint,oif_will_reconstruct)
        class(t_propagator) :: self
        type(t_field) :: f
        character(*) :: name
        logical,optional :: oif_will_reconstruct
        logical,optional :: ois_adjoint, ois_simple

        !field
        ! call f%init(name)
        f%name=name

        f%is_adjoint=either(ois_adjoint,.false.,present(ois_adjoint))

        call f%init_bloom

        !f%if_will_reconstruct=either(oif_will_reconstruct,.not.f%is_adjoint,present(oif_will_reconstruct))
        !if(f%if_will_reconstruct) call f%init_boundary
        call f%init_boundary_pressure

        call alloc(f%p     , [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%p_prev, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%p_next, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        ! if(.not.either(ois_simple,.false.,present(ois_simple))) then
            call alloc(f%dp_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
            call alloc(f%dp_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
            call alloc(f%dp_dy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

            call alloc(f%dpzz_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
            call alloc(f%dpxx_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
            call alloc(f%dpyy_dy, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

            call alloc(f%lap, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        ! endif

        call alloc(f%poynz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(f%poynx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

    end subroutine

    subroutine init_correlate(self,corr,name)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name
        
        corr%name=name

        ! if(name(1:1)=='g') then !gradient components
            call alloc(corr%gbuo, m%nz,m%nx,m%ny)
            call alloc(corr%gikpa,m%nz,m%nx,m%ny)
            ! call alloc(corr%giqp,m%nz,m%nx,m%ny)
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

        if(allocated(correlate_image)) then
            call correlate_assemble(corr%ipp,   correlate_image(:,:,:,1))
            call correlate_assemble(corr%ibksc, correlate_image(:,:,:,2))
            call correlate_assemble(corr%ifwsc, correlate_image(:,:,:,3))
        endif

        if(allocated(correlate_gradient)) then
            call correlate_assemble(corr%gbuo,  correlate_gradient(:,:,:,1))
            call correlate_assemble(corr%gikpa, correlate_gradient(:,:,:,2))
            ! call correlate_assemble(corr%giqp, correlate_gradient(:,:,:,3))
        endif        
        
    end subroutine

    
    !========= Derivations =================
    !PDE:      A U  = ϰ₀(C₂ ∂ₜ² -iC₁∂ₜ +C₀)U - ∇·b∇U = f
    !Adjoint:  AᴴUᵃ = ϰ₀(C₂ᴴ∂ₜ² -iC₁∂ₜ +C₀)Uᵃ- ∇·b∇Uᵃ= d
    !where
    !U=u+i*Hilbert(u) is the complex-valued wavefield (pressure component)
    !Uᵃ is the associated complex-valued adjoint field
    !f=fp*δ(x-xs), d is recorded data
    !b=ρ⁻¹ is buoyancy, 
    !ϰ₀=κ₀⁻¹ is bulk compliance (inverse of modulus) at the reference frequency
          
    !PDE:    (C₂∂ₜ² -iC₁∂ₜ +C₀)U = κ₀∇·b∇U + κ₀f
    !Discretized:
    !    C₂(Uⁿ⁺¹-2Uⁿ+Uⁿ⁻¹)/dt² -iC₁(Uⁿ-Uⁿ⁻¹)/dt +C₀Uⁿ =κ₀Lap
    !Forward:
    !    Uⁿ⁺¹ -2Uⁿ +Uⁿ⁻¹ -iC₁/C₂(Uⁿ-Uⁿ⁻¹)dt +C₀/C₂Uⁿdt² =1/C₂*dt²*κ₀Lap
    !    Uⁿ⁺¹ =2Uⁿ -Uⁿ⁻¹ +iC₁/C₂(Uⁿ-Uⁿ⁻¹)dt -C₀/C₂Uⁿdt² +1/C₂*dt²*κ₀Lap
    !Backward:
    !    Uⁿ⁻¹ -2Uⁿ +Uⁿ⁺¹ -iC₁/C₂(Uⁿ⁺¹-Uⁿ)dt +C₀/C₂Uⁿdt² =1/C₂*dt²*κ₀Lap
    !    Uⁿ⁻¹ =2Uⁿ -Uⁿ⁺¹ +iC₁/C₂(Uⁿ⁺¹-Uⁿ)dt -C₀/C₂Uⁿdt² +1/C₂*dt²*κ₀Lap
     
    !Adjoint:(C₂ᴴ∂ₜ² -iC₁∂ₜ +C₀)Uᵃ = κ₀∇·b∇Uᵃ + κ₀d
    !Discretized:
    !    C₂ᴴ(Uᵃⁿ⁺¹-2Uᵃⁿ+Uᵃⁿ⁻¹)/dt² +iC₁(Uᵃⁿ-Uᵃⁿ⁻¹)/dt +C₀Uᵃⁿ =κ₀Lap
    !Forward:
    !    Uᵃⁿ⁺¹ -2Uᵃⁿ +Uᵃⁿ⁻¹ -iC₁/C₂ᴴ(Uᵃⁿ-Uᵃⁿ⁻¹)dt +C₀/C₂ᴴUᵃⁿdt² =1/C₂ᴴ*dt²*κ₀Lap
    !    Uᵃⁿ⁺¹ =2Uᵃⁿ -Uᵃⁿ⁻¹ +iC₁/C₂ᴴ(Uᵃⁿ-Uᵃⁿ⁻¹)dt -C₀/C₂ᴴUᵃⁿdt² +1/C₂ᴴ*dt²*κ₀Lap
    !Backward:
    !    Uᵃⁿ⁻¹ -2Uᵃⁿ +Uᵃⁿ⁺¹ -iC₁/C₂ᴴ(Uᵃⁿ⁺¹-Uᵃⁿ)dt +C₀/C₂ᴴUᵃⁿdt² =1/C₂ᴴ*dt²*κ₀Lap
    !    Uᵃⁿ⁻¹ =2Uᵃⁿ -Uᵃⁿ⁺¹ +iC₁/C₂ᴴ(Uᵃⁿ⁺¹-Uᵃⁿ)dt -C₀/C₂ᴴUᵃⁿdt² +1/C₂ᴴ*dt²*κ₀Lap

    subroutine forward(self,fld_reU,fld_imU)
        class(t_propagator) :: self
        type(t_field) :: fld_reU,fld_imU
        real, allocatable :: mask(:,:,:)
        real,parameter :: time_dir=1. !time direction
        
        !seismo
        call alloc(fld_reU%seismo,shot%nrcv,self%nt)
        call alloc(fld_imU%seismo,shot%nrcv,self%nt)

        ! call hilbert_transform(fld_reU%wavelet,fld_imU%wavelet,1,self%nt,o_axis=2)

        tt1=0.; tt2=0.; tt3=0.; tt4=0.; tt5=0.; tt6=0.; tt7=0.

        ift=1; ilt=self%nt

        do it=ift,ilt
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_reU%check_value
                call fld_imU%check_value
            endif
            
            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add pressure
            call cpu_time(tic)
            call self%inject_pressure(fld_reU,time_dir,it)
            call self%inject_pressure(fld_imU,time_dir,it)
            ! print *, 'time=',it,'wavelet for u:',size(fld_reU%wavelet,1),size(fld_reU%wavelet,2),'wavelet for v:',size(fld_imU%wavelet,1),size(fld_imU%wavelet,2)
            call cpu_time(toc)
            tt1=tt1+toc-tic

            !step 2: save p^it+1 in boundary layers
            ! if(fld_E0%if_will_reconstruct) then
                call cpu_time(tic)
                call fld_reU%boundary_transport_pressure('save',it)
                call fld_imU%boundary_transport_pressure('save',it)
                call cpu_time(toc)
                tt2=tt2+toc-tic
            ! endif

            ! !step 3: set hardBC
            ! call cpu_time(tic)
            ! call self%set_pressure(fld_E0,time_dir,it)
            ! call cpu_time(toc)
            ! tt3=tt3+toc-tic

            !step 4: update pressure
            call cpu_time(tic)
            call self%update_pressure(fld_reU,fld_imU,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            ! if(mod(it,1000)==0) then
            !     call self%gaussian_smooth(fld_reU)
            !     call self%gaussian_smooth(fld_imU)
            ! endif
            

            ! if(mod(it,20)==0) then
            !     call build_mask(mask, fld_reU, it)
            !     ! print *,'mask shape', size(mask,1),size(mask,2),size(mask,3), 'fld_reU%p shape',size(fld_reU%p,1),size(fld_reU%p,2),size(fld_reU%p,3)
            !     ! open(unit=10, file='../../Demo/11_Marmousi/mask.bin', form='unformatted', access='stream', status='replace')
            !     ! write(10) mask
            !     ! close(10)
            !     ! stop
            !     call fft_gassian_filt(fld_reU, mask, it)
            !     call fft_gassian_filt(fld_imU, mask, it)
            !     deallocate(mask)
            ! endif
            

            !step 5: evolve pressure, it -> it+1
            call cpu_time(tic)
            call self%evolve_pressure(fld_reU,time_dir,it)
            call self%evolve_pressure(fld_imU,time_dir,it)
            call cpu_time(toc)
            tt6=tt6+toc-tic

            !step 6: sample p^it+1 at receivers
            call cpu_time(tic)
            call self%extract(fld_reU,it)
            call self%extract(fld_imU,it)
            call cpu_time(toc)
            tt7=tt7+toc-tic

            !snapshot
            call fld_reU%write(it)
            call fld_imU%write(it)

        enddo
        
        ! deallocate(wavelet_hilb)
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

    subroutine adjoint(self, fld_reA,fld_imA, fld_reU,fld_imU, A_star_U)
        class(t_propagator) :: self
        type(t_field) :: fld_reA,fld_imA, fld_reU,fld_imU
        type(t_correlate) :: A_star_U

        real, allocatable :: mask(:,:,:)

        real,parameter :: time_dir=-1. !time direction

        !reinitialize absorbing boundary for incident wavefield reconstruction
        call fld_reU%reinit
        call fld_imU%reinit

        if(if_record_adjseismo)  call alloc(fld_reA%seismo,1,self%nt)
        if(if_record_adjseismo)  call alloc(fld_imA%seismo,1,self%nt)
        
        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.
        
        ift=1; ilt=self%nt
        do it=ilt,ift,int(time_dir)
            if(mod(it,500)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_reA%check_value
                call fld_imA%check_value
                call fld_reU%check_value
                call fld_imU%check_value
                
            endif

            ! if(present(o_sf)) then
                !backward step 5: it+1 -> it
                call cpu_time(tic)
                call self%evolve_pressure(fld_reU,time_dir,it)
                call self%evolve_pressure(fld_imU,time_dir,it)
                call cpu_time(toc)
                tt1=tt1+toc-tic

                ! !backward step 2: retrieve p^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_reU%boundary_transport_pressure('load',it)
                call fld_imU%boundary_transport_pressure('load',it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 4:
                call cpu_time(tic)
                call self%update_pressure(fld_reU,fld_imU,time_dir,it)
                call cpu_time(toc)
                tt4=tt4+toc-tic

                ! if(mod(it,100)==0) then
                !     call build_mask(mask, fld_reU, it)
                !     ! print *,'mask size=',size(mask,1),size(mask,2),'it=',it
                !     ! if(it==4400) then
                !     ! open(unit=10, file='../../Demo/11_Marmousi/mask.bin', form='unformatted', access='stream', status='replace')
                !     ! write(10) mask
                !     ! close(10)
                !     ! stop
                !     ! endif 
                    
                !     call fft_gassian_filt(fld_reU, mask, it)
                !     call fft_gassian_filt(fld_imU, mask, it)
                !     deallocate(mask)
                !     ! print *,'fld_reU shape=',size(fld_reU%p,1),size(fld_reU%p,2)
                !  endif
                 

                !backward step 1: rm p^it at source
                call cpu_time(tic)
                call self%inject_pressure(fld_reU,time_dir,it)
                call self%inject_pressure(fld_imU,time_dir,it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            ! endif

            !adjoint step 6: inject to p^it+1 at receivers
            call cpu_time(tic)
            call self%inject_pressure(fld_reA,time_dir,it)
            call self%inject_pressure(fld_imA,time_dir,it)
            call cpu_time(toc)
            tt8=tt8+toc-tic

            !adjoint step 4:
            call cpu_time(tic)
            call self%update_pressure(fld_reA,fld_imA,time_dir,it)
            ! call self%update_pressure(fld_reA,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            ! if(mod(it,100)==0) then
            !     call build_mask(mask, fld_reA, it)
            !     call fft_gassian_filt(fld_reA, mask, it)
            !     call fft_gassian_filt(fld_imA, mask, it)
            !     deallocate(mask)
            ! endif

            !image: rf%p^it star sf%p^it
            if(mod(it,irdt)==0) then
                ! call cpu_time(tic)
                ! call compute_poynting(fld_imA,fld_reA)
                ! call compute_poynting(fld_imU,fld_reU)
                ! call cpu_time(toc)
                ! tt3=tt3+toc-tic

                call cpu_time(tic)
                call cross_correlate_gradient(fld_reA,fld_imA,fld_reU,fld_imU,A_star_U,it)
                call cpu_time(toc)
                tt10=tt10+toc-tic
            endif

            !adjoint step 5
            ! this step is moved to update_pressure for easier management
            call cpu_time(tic)
            call self%evolve_pressure(fld_reA,time_dir,it)
            call self%evolve_pressure(fld_imA,time_dir,it)
            call cpu_time(toc)
            tt11=tt11+toc-tic

            !adjoint step 1: sample p^it at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_reA,it)
                call self%extract(fld_imA,it)
                call cpu_time(toc)
                tt12=tt12+toc-tic
            endif


            !--------------------------------------------------------!
         
            !snapshot
            call fld_reA%write(it,o_suffix='_rev')
            call fld_imA%write(it,o_suffix='_rev')
            call fld_reU%write(it,o_suffix='_rev')
            call fld_imU%write(it,o_suffix='_rev')

            call A_star_U%write(it,o_suffix='_rev')

        enddo

        call A_star_U%scale(m%cell_volume*rdt)


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


    subroutine inject_pressure(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        if(.not. f%is_adjoint) then
        
            if(if_hicks) then
                ifz=shot%src%ifz-cb%ioz+1; ilz=shot%src%ilz-cb%ioz+1
                ifx=shot%src%ifx-cb%iox+1; ilx=shot%src%ilx-cb%iox+1
                ify=shot%src%ify-cb%ioy+1; ily=shot%src%ily-cb%ioy+1
            else
                iz=shot%src%iz-cb%ioz+1
                ix=shot%src%ix-cb%iox+1
                iy=shot%src%iy-cb%ioy+1
            endif

            ! call hilbert_transform2(f%wavelet,wavelet_hilb,1,self%nt)
            ! wl=time_dir*f%wavelet(1,it)*wavelet_scaler
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler

            if(if_hicks) then
                select case (shot%src%comp)
                case ('p') !explosion
                    f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) + wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef
                
                case ('dpdz') !vertical force
                    f%p(ifz+1:ilz+1,ifx:ilx,ify:ily) = f%p(ifz+1:ilz+1,ifx:ilx,ify:ily) + wl*self%kpa(ifz+1:ilz+1,ifx:ilx,ify:ily)*inv_2dz*shot%src%interp_coef
                    f%p(ifz-1:ilz-1,ifx:ilx,ify:ily) = f%p(ifz-1:ilz-1,ifx:ilx,ify:ily) - wl*self%kpa(ifz-1:ilz-1,ifx:ilx,ify:ily)*inv_2dz*shot%src%interp_coef
                
                endselect

            else
                select case (shot%src%comp)
                case ('p') !explosion
                    f%p(iz,ix,iy)                = f%p(iz,ix,iy)                + wl*self%kpa(iz,ix,iy)
                
                case ('dpdz') !vertical force
                    f%p(iz+1,ix,iy)                = f%p(iz+1,ix,iy)            + wl*self%kpa(iz+1,ix,iy)*inv_2dz
                    f%p(iz-1,ix,iy)                = f%p(iz-1,ix,iy)            - wl*self%kpa(iz-1,ix,iy)*inv_2dz
                
                endselect

            endif

            return

        endif

        do i=1,shot%nrcv

            if(if_hicks) then
                ifz=shot%rcv(i)%ifz-cb%ioz+1; ilz=shot%rcv(i)%ilz-cb%ioz+1
                ifx=shot%rcv(i)%ifx-cb%iox+1; ilx=shot%rcv(i)%ilx-cb%iox+1
                ify=shot%rcv(i)%ify-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1
            else
                iz=shot%rcv(i)%iz-cb%ioz+1
                ix=shot%rcv(i)%ix-cb%iox+1
                iy=shot%rcv(i)%iy-cb%ioy+1
            endif

            wl = f%wavelet(i,it)*wavelet_scaler    !no time_dir needed!

            if(if_hicks) then 

                select case (shot%rcv(i)%comp)
                case ('p') !adjsource for pressure
                    f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) +wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef

                case ('dpdz') !adjsource for vertical force 
                    f%p(ifz+1:ilz+1,ifx:ilx,ify:ily) = f%p(ifz+1:ilz+1,ifx:ilx,ify:ily) + wl*self%kpa(ifz+1:ilz+1,ifx:ilx,ify:ily)*inv_2dz*shot%rcv(i)%interp_coef
                    f%p(ifz-1:ilz-1,ifx:ilx,ify:ily) = f%p(ifz-1:ilz-1,ifx:ilx,ify:ily) - wl*self%kpa(ifz-1:ilz-1,ifx:ilx,ify:ily)*inv_2dz*shot%rcv(i)%interp_coef

                endselect

            else
                select case (shot%rcv(i)%comp)
                case ('p') !adjsource for pressure
                    f%p(iz,ix,iy)                = f%p(iz,ix,iy)                +wl*self%kpa(iz,ix,iy)
            
                case ('dpdz') !adjsource for vertical force 
                    f%p(iz+1,ix,iy)              = f%p(iz+1,ix,iy)              + wl*self%kpa(iz+1,ix,iy)*inv_2dz
                    f%p(iz-1,ix,iy)              = f%p(iz-1,ix,iy)              - wl*self%kpa(iz-1,ix,iy)*inv_2dz

                endselect

            endif

        enddo
        
    end subroutine

    subroutine update_pressure(self,f_re,f_im,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f_re, f_im

        complex,dimension(:,:,:),allocatable :: Uprev, U, Unext, Lap

                
        ! !necessary after computing the secondary source
        ! f%lap=0.

        ifz=f_re%bloom(1,it)
        if(m%is_freesurface) ifz=max(ifz,1)
        ilz=f_re%bloom(2,it)
        ifx=f_re%bloom(3,it)
        ilx=f_re%bloom(4,it)
        
        if(m%is_cubic) then
            ! call fd3d_pressure(f%p,                                      &
            !                    f%dp_dz,f%dp_dx,f%dp_dy,                  &
            !                    self%buoz,self%buox,self%buoy,self%kpa,   &
            !                    ifz,f%bloom(2,it),f%bloom(3,it),f%bloom(4,it))
        else
            call fd2d_laplacian(f_re%p,                                         &
                                f_re%dp_dz,f_re%dp_dx,f_re%dpzz_dz,f_re%dpxx_dx,&
                                f_re%lap,                                       &
                                self%buoz,self%buox,                            &
                                ifz,ilz,ifx,ilx)
            call fd2d_laplacian(f_im%p,                                         &
                                f_im%dp_dz,f_im%dp_dx,f_im%dpzz_dz,f_im%dpxx_dx,&
                                f_im%lap,                                       &
                                self%buoz,self%buox,                            &
                                ifz,ilz,ifx,ilx)
        endif

        if(.not. f_re%is_adjoint) then !PDE
            if(time_dir>0. ) then !forward in time
                Uprev = cmplx(f_re%p_prev,f_im%p_prev)
                U     = cmplx(f_re%p     ,f_im%p     )
                Lap   = cmplx(f_re%lap   ,f_im%lap   )

               !Uⁿ⁺¹  =2Uⁿ -Uⁿ⁻¹  +  i C₁/C₂    (Uⁿ-Uⁿ⁻¹) dt      -   C₀/C₂ Uⁿdt² +      1/C₂*dt²*Lap
                Unext =2*U -Uprev +c_i*self%C1n*(U-Uprev)*self%dt -self%C0n*U*dt2 +self%invC2*dt2*self%kpa*Lap
               ! Unext =2*U -Uprev +dt2*self%kpa*Lap
                f_re%p_next =  real(Unext)
                f_im%p_next = aimag(Unext)

            else !backward in time
            !for 1st-order time derivative, use forward FD instead of backward FD like forward PDE
            !this is fine but cannot make the field doesn't disappear after backpropagation
                Unext = cmplx(f_re%p_next,f_im%p_next)
                U     = cmplx(f_re%p     ,f_im%p     )
                Lap   = cmplx(f_re%lap   ,f_im%lap   )

               !Uⁿ⁻¹  =2Uⁿ -Uⁿ⁺¹  +  i C₁/C₂    (Uⁿ⁺¹-Uⁿ) dt      -   C₀/C₂ Uⁿdt² +      1/C₂*dt²*Lap
                Uprev =2*U -Unext +c_i*self%C1n*(Unext-U)*self%dt -self%C0n*U*dt2 +self%invC2*dt2*self%kpa*Lap
                f_re%p_prev =  real(Uprev)
                f_im%p_prev = aimag(Uprev)

            endif

        else !Adjoint, backward in time
                Unext = cmplx(f_re%p_next,f_im%p_next)
                U     = cmplx(f_re%p     ,f_im%p     )
                Lap   = cmplx(f_re%lap   ,f_im%lap   )

               !Uᵃⁿ⁻¹ =2Uᵃⁿ-Uᵃⁿ⁺¹ +  i    C₁/C₂ᴴ (Uᵃⁿ⁺¹-Uᵃⁿ)dt     -   C₀/C₂ᴴ Uᵃⁿdt² +     1/C₂ᴴ*dt²*Lap
                Uprev =2*U -Unext +c_i*self%C1nH*(Unext-U)*self%dt -self%C0nH*U*dt2 +self%invC2H*dt2*self%kpa*Lap
               ! Uprev =2*U -Unext +dt2*self%kpa*Lap
                f_re%p_prev =  real(Uprev)
                f_im%p_prev = aimag(Uprev)
        endif

        deallocate(Uprev, U, Unext, Lap) !save some RAM

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

    subroutine evolve_pressure(self,f,time_dir,it)
        class(t_propagator) :: self
        type(t_field) :: f

        real,dimension(:,:,:),pointer :: tmp

        if(time_dir>0.) then !in forward time
            tmp=>f%p_prev
            f%p_prev=>f%p
            f%p     =>f%p_next
            f%p_next=>tmp

            ! f%p_prev = f%p
            ! f%p      = f%p_next
            ! !f%p_next = f%p_prev

        else !in reverse time
            tmp=>f%p_next
            f%p_next=>f%p
            f%p     =>f%p_prev
            f%p_prev=>tmp

            ! f%p_next = f%p
            ! f%p      = f%p_prev
            ! !f%p_prev = f%p_next

        endif

        nullify(tmp)
        
    end subroutine

    subroutine compute_poynting(v,u)
        type(t_field) :: v, u

        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: E2, ph !envelope squared & inst phase
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: dph_dz, dph_dx

    	!E=sqrt(u%p*u%p+v%p*v%p)
        E2=u%p*u%p+v%p*v%p
        ph=atan2(v%p,u%p)

    	do ix=cb%ifx+1,cb%ilx-1
    	do iz=cb%ifz+1,cb%ilz-1
    	    dph_dz(iz,ix,1) = asin(sin(ph(iz+1,ix,1) - ph(iz-1,ix,1)))*inv_2dz
    	    dph_dx(iz,ix,1) = asin(sin(ph(iz,ix+1,1) - ph(iz,ix-1,1)))*inv_2dx
    	enddo
    	enddo

        u%poynz=E2*dph_dz
        u%poynx=E2*dph_dx

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
                        f%seismo(i,it)=sum(f%p(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)

                    case ('dpdz')
                        f%seismo(i,it)=( sum(f%p(ifz+1:ilz+1,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef) &
                                        -sum(f%p(ifz-1:ilz-1,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef) )*inv_2dz 

                        ! case ('vx')
                        ! f%seismo(i,it)=sum(f%vx(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
                        ! case ('vy')
                        ! f%seismo(i,it)=sum(f%vy(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef)
                    end select
                    
                else
                    select case (shot%rcv(i)%comp)
                    case ('p') !p[iz,ix,iy]
                        f%seismo(i,it)=f%p(iz,ix,iy)

                    case ('dpdz')
                        f%seismo(i,it)=( f%p(iz+1,ix,iy) &
                                        -f%p(iz-1,ix,iy) )*inv_2dz 

                        ! case ('vz') !vz[iz-0.5,ix,iy]
                        ! f%seismo(i,it)=f%vz(iz,ix,iy)
                        ! case ('vx') !vx[iz,ix-0.5,iy]
                        ! f%seismo(i,it)=f%vx(iz,ix,iy)
                        ! case ('vy') !vy[iz,ix,iy-0.5]
                        ! f%seismo(i,it)=f%vy(iz,ix,iy)
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
                    
                case ('dpdz')
                    f%seismo(1,it)=( sum(f%p(ifz+1:ilz+1,ifx:ilx,ify:ily)*shot%src%interp_coef) &
                                    -sum(f%p(ifz-1:ilz-1,ifx:ilx,ify:ily)*shot%src%interp_coef) )*inv_2dz 

                    ! case ('vz')
                    ! f%seismo(1,it)=sum(f%vz(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
                    
                    ! case ('vx')
                    ! f%seismo(1,it)=sum(f%vx(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
                    
                    ! case ('vy')
                    ! f%seismo(1,it)=sum(f%vy(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef)
                    
                end select
                
            else
                select case (shot%src%comp)
                case ('p') !p[iz,ix,iy]
                    f%seismo(1,it)=f%p(iz,ix,iy)

                case ('dpdz')
                    f%seismo(1,it)=( f%p(iz+1,ix,iy) &
                                    -f%p(iz-1,ix,iy) )*inv_2dz 
                    
                    ! case ('vz') !vz[iz-0.5,ix,iy]
                    ! f%seismo(1,it)=f%vz(iz,ix,iy)
                    
                    ! case ('vx') !vx[iz,ix-0.5,iy]
                    ! f%seismo(1,it)=f%vx(iz,ix,iy)
                    
                    ! case ('vy') !vy[iz,ix,iy-0.5]
                    ! f%seismo(1,it)=f%vy(iz,ix,iy)
                    
                end select
                
            endif
        
    end subroutine
        
    subroutine final(self)
        type(t_propagator) :: self
        call dealloc(self%buoz, self%buox, self%buoy, self%kpa)
    end subroutine


    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !Kₘ<a|Au> = Kₘ<a|ϰ∂ₜ²u - ∇·b∇u>
    !for ϰ: Kₘ<a|Au> = ∫ a ∂ₜ²u dt =-∫ ∂ₜa ∂ₜu dt, or = ∫ a κ∇·b∇u dt 
    !for b: Kₘ<a|Au> = -Kₘ<a|∇·b∇u> = ∫ ∇a·∇u dt

    subroutine cross_correlate_gradient(reA,imA, reU,imU, corr,it)
        type(t_field), intent(in) :: reA, imA, reU, imU
        type(t_correlate) :: corr

        real,dimension(:,:,:),allocatable,save :: re_gikpa_rere, re_gikpa_imim, im_gikpa

        complex,dimension(:,:,:),allocatable :: Ulap, Aconj

        !nonzero only when sf touches rf
        ifz=max(reU%bloom(1,it),reA%bloom(1,it),2)
        ilz=min(reU%bloom(2,it),reA%bloom(2,it),cb%mz)
        ifx=max(reU%bloom(3,it),reA%bloom(3,it),1)
        ilx=min(reU%bloom(4,it),reA%bloom(4,it),cb%mx)
        ! ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ! ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)


        ! Ulap  = cmplx( reU%lap(1:m%nz,1:m%nx,1:m%ny), imU%lap(1:m%nz,1:m%nx,1:m%ny) )
        ! Aconj = cmplx( reA%p  (1:m%nz,1:m%nx,1:m%ny),-imA%p  (1:m%nz,1:m%nx,1:m%ny) )

        ! !for gikpa
        ! corr%gikpa = corr%gikpa + &
        !     real ( Aconj*ppg%kpa(1:m%nz,1:m%nx,1:m%ny)*Ulap )


        !other parts
        call alloc(re_gikpa_rere,cb%mz,cb%mx,cb%my, oif_protect=.true.)
        call alloc(re_gikpa_imim,cb%mz,cb%mx,cb%my, oif_protect=.true.)
        call alloc(im_gikpa,     cb%mz,cb%mx,cb%my, oif_protect=.true.)
        
        re_gikpa_rere = re_gikpa_rere + &
                   reA%p(1:m%nz,1:m%nx,1:m%ny)*ppg%kpa(1:m%nz,1:m%nx,1:m%ny)*reU%lap(1:m%nz,1:m%nx,1:m%ny)
        re_gikpa_imim = re_gikpa_imim + &
                   imA%p(1:m%nz,1:m%nx,1:m%ny)*ppg%kpa(1:m%nz,1:m%nx,1:m%ny)*imU%lap(1:m%nz,1:m%nx,1:m%ny)

        corr%gikpa  = re_gikpa_rere + re_gikpa_imim

        ! ! im_gikpa = im_gikpa + &
        !     ! aimag( Aconj(1:m%nz,1:m%nx,1:m%ny)*ppg%kpa(1:m%nz,1:m%nx,1:m%ny)*Ulap(1:m%nz,1:m%nx,1:m%ny) )
        ! im_gikpa = im_gikpa + ppg%kpa(1:m%nz,1:m%nx,1:m%ny)* (&
        !      reA%p(1:m%nz,1:m%nx,1:m%ny)*imU%lap(1:m%nz,1:m%nx,1:m%ny) &
        !     -imA%p(1:m%nz,1:m%nx,1:m%ny)*reU%lap(1:m%nz,1:m%nx,1:m%ny) &
        !     )

        ! call sysio_write('re_gikpa_rere',re_gikpa_rere,m%n)
        ! call sysio_write('re_gikpa_imim',re_gikpa_imim,m%n)
        ! call sysio_write('re_gikpa',     corr%gikpa,   m%n)
        ! call sysio_write('im_gikpa',     im_gikpa,     m%n)

        ! !for gbuo
        ! call fd2d_grho(reA%p,reU%p,corr%gbuo,   ifz,ilz,ifx,ilx)

        ! !for giqp
        ! corr%giqp = corr%giqp + &
        !     A(1:m%nz,1:m%nx,1:m%ny) * ppg%kpa(1:m%nz,1:m%nx,1:m%ny)*Ulap(1:m%nz,1:m%nx,1:m%ny)        

    end subroutine

    subroutine cross_correlate_image(rf,sf,corr,it)
        type(t_field), intent(in) :: rf, sf
        type(t_correlate) :: corr

        !nonzero only when sf touches rf
        ifz=max(sf%bloom(1,it),rf%bloom(1,it),2)
        ilz=min(sf%bloom(2,it),rf%bloom(2,it),cb%mz)
        ifx=max(sf%bloom(3,it),rf%bloom(3,it),1)
        ilx=min(sf%bloom(4,it),rf%bloom(4,it),cb%mx)
        ! ify=max(sf%bloom(5,it),rf%bloom(5,it),1)
        ! ily=min(sf%bloom(6,it),rf%bloom(6,it),cb%my)

        call imag2d(rf%p,sf%p,&
                    rf%poynz,rf%poynx,sf%poynz,sf%poynx, &
                    corr%ipp,corr%ibksc,corr%ifwsc, &
                    ifz,ilz,ifx,ilx)

    end subroutine

    !========= Finite-Difference on flattened arrays ==================

    subroutine fd2d_laplacian(p,dp_dz,dp_dx,&
                                dpzz_dz,dpxx_dx,&
                                lap,&
                                buoz,buox,&
                                ifz,ilz,ifx,ilx)
        real,dimension(*) :: p,dp_dz,dp_dx
        real,dimension(*) :: dpzz_dz,dpxx_dx
        real,dimension(*) :: lap
        real,dimension(*) :: buoz,buox

        real,dimension(:),allocatable :: pzz,pxx
        call alloc(pzz,cb%n)
        call alloc(pxx,cb%n)

        nz=cb%nz
        nx=cb%nx
        
        !flux: b∇u ~= ( bz*∂zᵇp , bx*∂ₓᵇp )
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dp_dz_,dp_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx+2,ilx-1
            !dir$ simd
            do iz = ifz+2,ilz-1

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                
                iz_ixm2=i  -2*nz !iz,ix-2
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1

                dp_dz_ = c1z*(p(iz_ix) - p(izm1_ix)) +c2z*(p(izp1_ix)-p(izm2_ix))
                dp_dx_ = c1x*(p(iz_ix) - p(iz_ixm1)) +c2x*(p(iz_ixp1)-p(iz_ixm2))

                dp_dz(iz_ix) = cpml%b_z_half(iz)*dp_dz(iz_ix) + cpml%a_z_half(iz)*dp_dz_
                dp_dx(iz_ix) = cpml%b_x_half(ix)*dp_dx(iz_ix) + cpml%a_x_half(ix)*dp_dx_

                dp_dz_ = dp_dz_/cpml%kpa_z_half(iz) + dp_dz(iz_ix)
                dp_dx_ = dp_dx_/cpml%kpa_x_half(ix) + dp_dx(iz_ix)

                pzz(iz_ix) = buoz(iz_ix)*dp_dz_
                pxx(iz_ix) = buox(iz_ix)*dp_dx_

            enddo
        enddo
        !$omp end do
        !$omp end parallel
        
        !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dpzz_dz_,dpxx_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx+1,ilx-2
            !dir$ simd
            do iz = ifz+1,ilz-2

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1
                iz_ixp2=i  +2*nz !iz,ix+2
                
                dpzz_dz_ = c1z*(pzz(izp1_ix) - pzz(iz_ix))  +c2z*(pzz(izp2_ix) - pzz(izm1_ix))
                dpxx_dx_ = c1x*(pxx(iz_ixp1) - pxx(iz_ix))  +c2x*(pxx(iz_ixp2) - pxx(iz_ixm1))

                dpzz_dz(iz_ix) = cpml%b_z(iz)*dpzz_dz(iz_ix) + cpml%a_z(iz)*dpzz_dz_
                dpxx_dx(iz_ix) = cpml%b_x(ix)*dpxx_dx(iz_ix) + cpml%a_x(ix)*dpxx_dx_

                dpzz_dz_ = dpzz_dz_/cpml%kpa_z(iz) + dpzz_dz(iz_ix)
                dpxx_dx_ = dpxx_dx_/cpml%kpa_x(ix) + dpxx_dx(iz_ix)

                lap(iz_ix) = dpzz_dz_ + dpxx_dx_

            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine fd2d_laplacian_nocpml(p,&
                                    lap,&
                                    buoz,buox,&
                                    ifz,ilz,ifx,ilx)
        real,dimension(*) :: p
        real,dimension(*) :: lap
        real,dimension(*) :: buoz,buox

        real,dimension(:),allocatable :: pzz,pxx
        call alloc(pzz,cb%n)
        call alloc(pxx,cb%n)

        nz=cb%nz
        nx=cb%nx
        
        !flux: b∇u ~= ( bz*∂zᵇp , bx*∂ₓᵇp )
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,&
        !$omp         dp_dz_,dp_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx+2,ilx-1
            !dir$ simd
            do iz = ifz+2,ilz-1

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                
                iz_ixm2=i  -2*nz !iz,ix-2
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1

                dp_dz_ = c1z*(p(iz_ix) - p(izm1_ix)) +c2z*(p(izp1_ix)-p(izm2_ix))
                dp_dx_ = c1x*(p(iz_ix) - p(iz_ixm1)) +c2x*(p(iz_ixp1)-p(iz_ixm2))

                pzz(iz_ix) = buoz(iz_ix)*dp_dz_
                pxx(iz_ix) = buox(iz_ix)*dp_dx_

            enddo
        enddo
        !$omp end do
        !$omp end parallel
        
        !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dpzz_dz_,dpxx_dx_)
        !$omp do schedule(dynamic)
        do ix = ifx+1,ilx-2
            !dir$ simd
            do iz = ifz+1,ilz-2

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1

                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm1=i  -nz  !iz,ix-1
                iz_ixp1=i  +nz  !iz,ix+1
                iz_ixp2=i  +2*nz !iz,ix+2
                
                dpzz_dz_ = c1z*(pzz(izp1_ix) - pzz(iz_ix))  +c2z*(pzz(izp2_ix) - pzz(izm1_ix))
                dpxx_dx_ = c1x*(pxx(iz_ixp1) - pxx(iz_ix))  +c2x*(pxx(iz_ixp2) - pxx(iz_ixm1))

                lap(iz_ix) = dpzz_dz_ + dpxx_dx_

            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

    ! !2nd order
    ! subroutine fd2d_grho(rp,sp,corr,&
    !                     ifz,ilz,ifx,ilx)
    !     real,dimension(*) :: rp, sp
    !     real,dimension(*) :: corr

    !     real,dimension(:),allocatable :: pzz, pxx
    !     call alloc(pzz,cb%n)
    !     call alloc(pxx,cb%n)

    !     nz=cb%nz
    !     nx=cb%nx
        
        
    !     !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
    !     !$omp parallel default (shared)&
    !     !$omp private(iz,ix,i,&
    !     !$omp         iz_ixm1,izm1_ix,iz_ix,izp1_ix,iz_ixp1)
    !     !$omp do schedule(dynamic)
    !     do ix = ifx,ilx-1
    !         !dir$ simd
    !         do iz = ifz,ilz-1

    !             i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
    !             j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers

    !             iz_ixm1=i    -nz  !iz,ix-1
    !             izm1_ix=i-1       !iz-1,ix
    !             iz_ix  =i
    !             izp1_ix=i+1       !iz+1,ix
    !             iz_ixp1=i    +nz  !iz,ix+1

                
    !             corr(j) = corr(j) &
    !                 + (rp(izp1_ix)-rp(izm1_ix))*(sp(izp1_ix)-sp(izm1_ix))*inv_2dz*inv_2dz &
    !                 + (rp(iz_ixp1)-rp(iz_ixm1))*(sp(iz_ixp1)-sp(iz_ixm1))*inv_2dx*inv_2dx

    !         enddo
    !     enddo
    !     !$omp end do
    !     !$omp end parallel

    ! end subroutine

    !4th order
    subroutine fd2d_grho(rp,sp,corr,&
                        ifz,ilz,ifx,ilx)
        real,dimension(*) :: rp, sp
        real,dimension(*) :: corr

        real c1z, c2z, c1x, c2x !no confusion with c1z etc defined in the modules

        real,dimension(:),allocatable :: pzz, pxx
        call alloc(pzz,cb%n)
        call alloc(pxx,cb%n)

        nz=cb%nz
        nx=cb%nx
        
        c1z=2./3./m%dz; c2z=-1./12./m%dz
        c1x=2./3./m%dx; c2x=-1./12./m%dx
        
        !laplacian: ∇·b∇u ~= ∂zᶠ(bz*∂zᵇp) + ∂ₓᶠ(bx*∂ₓᵇp)
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         iz_ixm1,izm1_ix,iz_ix,izp1_ix,iz_ixp1)
        !$omp do schedule(dynamic)
        do ix = ifx+2,ilx-2
            !dir$ simd
            do iz = ifz+2,ilz-2

                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !grad has no boundary layers

                iz_ixm2=i    -2*nz
                iz_ixm1=i    -nz  !iz,ix-1
                izm2_ix=i-2
                izm1_ix=i-1       !iz-1,ix
                iz_ix  =i
                izp1_ix=i+1       !iz+1,ix
                izp2_ix=i+2
                iz_ixp1=i    +nz  !iz,ix+1
                iz_ixp2=i    +2*nz

                drp_dz = c1z*(rp(izp1_ix)-rp(izm1_ix)) +c2z*(rp(izp2_ix)-rp(izm2_ix))
                drp_dx = c1x*(rp(iz_ixp1)-rp(iz_ixm1)) +c2x*(rp(iz_ixp2)-rp(iz_ixm2))
                
                dsp_dz = c1z*(sp(izp1_ix)-sp(izm1_ix)) +c2z*(sp(izp2_ix)-sp(izm2_ix))
                dsp_dx = c1x*(sp(iz_ixp1)-sp(iz_ixm1)) +c2x*(sp(iz_ixp2)-sp(iz_ixm2))
                
                corr(j) = corr(j) + (drp_dz*dsp_dz + drp_dx*dsp_dx)
            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

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


    ! subroutine gaussian_smooth(self, f, dx_in)
    !     class(t_propagator) :: self
    !     type(t_field), intent(inout) :: f
    !     real(kind=4) :: sigma
    !     real(kind=4), intent(in), optional :: dx_in
        
    !     integer :: n, i, j, il, idx
    !     real(kind=4), allocatable :: temp(:), temp_prev(:), temp_next(:), kernel(:)
    !     real(kind=4) :: sig2, weight_sum, dist

    !     n = size(f%p)
    !     allocate(temp(n), temp_next(n),temp_prev(n))
    !     temp = 0.0

    !     if (present(dx_in)) then
    !         dx = dx_in
    !     else
    !         dx = 12.5 
    !     end if

    !     sigma = 200

    !     ! 高斯核范围（±4σ）
    !     il = int(4.0 * sigma / dx)
    !     allocate(kernel(-il:il))
    !     sig2 = 2.0 * sigma * sigma

    !     ! 构建高斯核
    !     do j = -il, il
    !         dist = real(j) * dx
    !         kernel(j) = exp(- (dist * dist) / sig2)
    !     end do

    !     ! 应用高斯平滑
    !     do i = 1, n
    !         weight_sum = 0.0
    !         do j = -il, il
    !             idx = i + j
    !             if (idx >= 1 .and. idx <= n) then
    !             temp(i) = temp(i) + kernel(j) * f%p(idx)
    !             weight_sum = weight_sum + kernel(j)
    !             end if
    !         end do
    !         if (weight_sum > 0.0) temp(i) = temp(i) / weight_sum
    !     end do
    !     f%p = temp
            
    !     do i = 1, n
    !         weight_sum = 0.0
    !         do j = -il, il
    !             idx = i + j
    !             if (idx >= 1 .and. idx <= n) then
    !             temp_next(i) = temp_next(i) + kernel(j) * f%p_next(idx)
    !             weight_sum = weight_sum + kernel(j)
    !             end if
    !         end do
    !         if (weight_sum > 0.0) temp_next(i) = temp_next(i) / weight_sum
    !     end do
    !     f%p_next = temp_next

    !     do i = 1, n
    !         weight_sum = 0.0
    !         do j = -il, il
    !             idx = i + j
    !             if (idx >= 1 .and. idx <= n) then
    !             temp_prev(i) = temp_prev(i) + kernel(j) * f%p_prev(idx)
    !             weight_sum = weight_sum + kernel(j)
    !             end if
    !         end do
    !         if (weight_sum > 0.0) temp_prev(i) = temp_prev(i) / weight_sum
    !     end do
    !     f%p_prev = temp_prev

    !     deallocate(temp, temp_next, temp_prev, kernel)
    ! end subroutine gaussian_smooth

    ! subroutine gaussian_smooth(self, f, dx_in, dy_in)
    !     class(t_propagator) :: self
    !     type(t_field), intent(inout) :: f
    !     real(kind=4), intent(in), optional :: dx_in, dy_in
    !     real(kind=4) :: sigma_x, sigma_y, dx, dy

    !     integer :: nx, ny, i, j, ii, jj, idx, idy
    !     integer :: ilx, ily
    !     real(kind=4), allocatable :: kernel_x(:), kernel_y(:)
    !     real(kind=4), allocatable :: temp(:,:), temp_next(:,:), temp_prev(:,:)
    !     real(kind=4) :: dist, sig2x, sig2y, weight_sum

    !     nz = size(f%p, 1)
    !     nx = size(f%p, 2)
    !     allocate(temp(nz, nx), temp_next(nz, nx), temp_prev(nz, nx))
    !     temp = 0.0

    !     ! 默认网格间距
    !     if (present(dx_in)) then
    !         dx = dx_in
    !     else
    !         dx = 12.5
    !     end if
    !     if (present(dy_in)) then
    !         dy = dy_in
    !     else
    !         dy = 12.5
    !     end if

    !     ! 默认标准差
    !     sigma_x = 1
    !     sigma_y = 1

    !     ! 核半径（±4σ）
    !     ilx = int(4.0 * sigma_x / dx)
    !     ily = int(4.0 * sigma_y / dy)
    !     allocate(kernel_x(-ilx:ilx))
    !     allocate(kernel_y(-ily:ily))

    !     sig2x = 2.0 * sigma_x * sigma_x
    !     sig2y = 2.0 * sigma_y * sigma_y

    !     ! 构造高斯核
    !     do i = -ilx, ilx
    !         dist = real(i) * dx
    !         kernel_x(i) = exp(- (dist * dist) / sig2x)
    !     end do
    !     do j = -ily, ily
    !         dist = real(j) * dy
    !         kernel_y(j) = exp(- (dist * dist) / sig2y)
    !     end do

    !     ! 平滑 f%p：二维卷积（可分离核，x后y）
    !     call apply_gaussian_2d(f%p, temp, kernel_x, ilx, dx, kernel_y, ily, dy)
    !     f%p(:,:,1) = temp

    !     call apply_gaussian_2d(f%p_next, temp_next, kernel_x, ilx, dx, kernel_y, ily, dy)
    !     f%p_next(:,:,1) = temp_next

    !     call apply_gaussian_2d(f%p_prev, temp_prev, kernel_x, ilx, dx, kernel_y, ily, dy)
    !     f%p_prev(:,:,1) = temp_prev

    !     deallocate(kernel_x, kernel_y, temp, temp_next, temp_prev)
    ! end subroutine gaussian_smooth


    ! subroutine apply_gaussian_2d(input, output, kernel_x, ilx, dx, kernel_y, ily, dy)
    !     implicit none
    !     real(kind=4), intent(in) :: input(:,:,:)
    !     real(kind=4), intent(in) :: kernel_x(-ilx:ilx), dx
    !     real(kind=4), intent(in) :: kernel_y(-ily:ily), dy
    !     real(kind=4), intent(out) :: output(:,:)
    !     integer, intent(in) :: ilx, ily

    !     integer :: nz, nx, i, j, ii, jj
    !     real(kind=4) :: sum, weight_sum

    !     ! 临时缓冲
    !     real(kind=4), allocatable :: tmp(:,:)

    !     nz = size(input,1)
    !     nx = size(input,2)

    !     allocate(tmp(nz, nx))
    !     tmp = 0.0

    !     ! 一维 z 向平滑
    !     do j = 1, nx
    !         do i = 1, nz
    !             sum = 0.0
    !             weight_sum = 0.0
    !             do ii = -ilx, ilx
    !                 if (i+ii >= 1 .and. i+ii <= nz) then
    !                     sum = sum + kernel_x(ii) * input(i+ii, j,1)
    !                     weight_sum = weight_sum + kernel_x(ii)
    !                 end if
    !             end do
    !             tmp(i,j) = sum / weight_sum
    !         end do
    !     end do

    !     ! 一维 x 向平滑
    !     do i = 1, nz
    !         do j = 1, nx
    !             sum = 0.0
    !             weight_sum = 0.0
    !             do jj = -ily, ily
    !                 if (j+jj >= 1 .and. j+jj <= nx) then
    !                     sum = sum + kernel_y(jj) * tmp(i, j+jj)
    !                     weight_sum = weight_sum + kernel_y(jj)
    !                 end if
    !             end do
    !             output(i,j) = sum / weight_sum
    !         end do
    !     end do

    !     deallocate(tmp)
    ! end subroutine apply_gaussian_2d

    subroutine fft_gassian_filt(f, mask, it)
        type(t_field) :: f
        complex(fftkind), allocatable :: spectrum_ppre(:,:,:), spectrum_p(:,:,:), spectrum_pnext(:,:,:), spectrum_filted_ppre(:,:,:), spectrum_filted_p(:,:,:), spectrum_filted_pnext(:,:,:)
        real, intent(in) :: mask(:,:,:)  
        character(len=200) :: filename
        ! if(time_dir>0.) then
        ifz=f%bloom(1,it)
        if(m%is_freesurface) ifz=max(ifz,1)
        ilz=f%bloom(2,it)
        ifx=f%bloom(3,it)
        ilx=f%bloom(4,it)
        ! else
        !     ifz=1
        !     ilz=size(f%p,1)
        !     ifx=1
        !     ilx=size(f%p,2)
        ! end if

        ! print *,'ifz,ilz=',ifz,ilz,'ifx,ilx=',ifx,ilx
        

        allocate(spectrum_ppre(ifz:ilz,ifx:ilx,1), spectrum_p(ifz:ilz,ifx:ilx,1), spectrum_pnext(ifz:ilz,ifx:ilx,1))
        allocate(spectrum_filted_ppre(ifz:ilz,ifx:ilx,1), spectrum_filted_p(ifz:ilz,ifx:ilx,1), spectrum_filted_pnext(ifz:ilz,ifx:ilx,1))
        ! print *,'spectrum_p shape:',size(spectrum_p,1),size(spectrum_p,2)
        spectrum_ppre=fft(cmplx(f%p_prev(ifz:ilz,ifx:ilx,:),0.0,kind=fftkind),inv=.TRUE.)
        spectrum_p=fft(cmplx(f%p(ifz:ilz,ifx:ilx,:),0.0,kind=fftkind),inv=.TRUE.)
        spectrum_pnext=fft(cmplx(f%p_next(ifz:ilz,ifx:ilx,:),0.0,kind=fftkind),inv=.TRUE.)
        ! write(filename,'("../../Demo/11_Marmousi/spectrum_p_",I0,".bin")') it
        ! open(unit=12, file=trim(filename), form='unformatted', access='stream', status='replace')
        ! write(12) spectrum_p
        ! close(12)
        ! stop
        spectrum_filted_ppre=spectrum_ppre*cmplx(mask,mask, kind=fftkind)
        spectrum_filted_p=spectrum_p*cmplx(mask,mask, kind=fftkind)
        spectrum_filted_pnext=spectrum_pnext*cmplx(mask,mask, kind=fftkind)
        ! write(filename,'("../../Demo/11_Marmousi/spectrum_p_filted_",I0,".bin")') it
        ! open(unit=13, file=trim(filename), form='unformatted', access='stream', status='replace')
        ! write(13) spectrum_filted_p
        ! close(13)
        ! print *,'ifz ilz ifx ilx', ifz,ilz,ifx,ilx, 'spectrum_p shape:',size(spectrum_p,1),size(spectrum_p,2),'mask shape',size(mask,1),size(mask,2)
        ! print *,'spectrum shape',size(spectrum_p,1),size(spectrum_p,2),'spectrum_filted',size(spectrum_filted_p,1),size(spectrum_filted_p,2)
        f%p_prev(ifz:ilz,ifx:ilx,:)=real(fft(spectrum_filted_ppre,inv=.False.))
        f%p(ifz:ilz,ifx:ilx,:)=real(fft(spectrum_filted_p,inv=.False.))
        f%p_next(ifz:ilz,ifx:ilx,:)=real(fft(spectrum_filted_pnext,inv=.False.))
        
        deallocate(spectrum_ppre, spectrum_p, spectrum_pnext)
        deallocate(spectrum_filted_ppre, spectrum_filted_p, spectrum_filted_pnext)
    end subroutine fft_gassian_filt

    subroutine fft_gassian_filt2(f, mask, it)
        type(t_field) :: f
        complex(fftkind), allocatable :: spectrum_ppre(:,:,:), spectrum_p(:,:,:), spectrum_pnext(:,:,:), spectrum_filted_ppre(:,:,:), spectrum_filted_p(:,:,:), spectrum_filted_pnext(:,:,:)
        real, intent(in) :: mask(:,:,:)  
        character(len=200) :: filename
        ! if(time_dir>0.) then
        ifz=f%bloom(1,it)
        if(m%is_freesurface) ifz=max(ifz,1)
        ilz=f%bloom(2,it)
        ifx=f%bloom(3,it)
        ilx=f%bloom(4,it)
        ! else
        !     ifz=1
        !     ilz=size(f%p,1)
        !     ifx=1
        !     ilx=size(f%p,2)
        ! end if

        ! print *,'ifz,ilz=',ifz,ilz,'ifx,ilx=',ifx,ilx
        

        allocate(spectrum_ppre(ifz:ilz,ifx:ilx,1), spectrum_p(ifz:ilz,ifx:ilx,1), spectrum_pnext(ifz:ilz,ifx:ilx,1))
        allocate(spectrum_filted_ppre(ifz:ilz,ifx:ilx,1), spectrum_filted_p(ifz:ilz,ifx:ilx,1), spectrum_filted_pnext(ifz:ilz,ifx:ilx,1))
        ! print *,'spectrum_p shape:',size(spectrum_p,1),size(spectrum_p,2)
        spectrum_ppre=fft(cmplx(f%p_prev(ifz:ilz,ifx:ilx,:),0.0,kind=fftkind),inv=.TRUE.)
        spectrum_p=fft(cmplx(f%p(ifz:ilz,ifx:ilx,:),0.0,kind=fftkind),inv=.TRUE.)
        spectrum_pnext=fft(cmplx(f%p_next(ifz:ilz,ifx:ilx,:),0.0,kind=fftkind),inv=.TRUE.)
        write(filename,'("../../Demo/11_Marmousi/spectrum_p_",I0,".bin")') it
        open(unit=12, file=trim(filename), form='unformatted', access='stream', status='replace')
        write(12) spectrum_p
        close(12)
        ! stop
        spectrum_filted_ppre=spectrum_ppre*cmplx(mask,0.0, kind=fftkind)
        spectrum_filted_p=spectrum_p*cmplx(mask,0.0, kind=fftkind)
        spectrum_filted_pnext=spectrum_pnext*cmplx(mask,0.0, kind=fftkind)
        write(filename,'("../../Demo/11_Marmousi/spectrum_p_filted_",I0,".bin")') it
        open(unit=13, file=trim(filename), form='unformatted', access='stream', status='replace')
        write(13) spectrum_filted_p
        close(13)
        ! print *,'ifz ilz ifx ilx', ifz,ilz,ifx,ilx, 'spectrum_p shape:',size(spectrum_p,1),size(spectrum_p,2),'mask shape',size(mask,1),size(mask,2)
        ! print *,'spectrum shape',size(spectrum_p,1),size(spectrum_p,2),'spectrum_filted',size(spectrum_filted_p,1),size(spectrum_filted_p,2)
        f%p_prev(ifz:ilz,ifx:ilx,:)=real(fft(spectrum_filted_ppre,inv=.False.))
        f%p(ifz:ilz,ifx:ilx,:)=real(fft(spectrum_filted_p,inv=.False.))
        f%p_next(ifz:ilz,ifx:ilx,:)=real(fft(spectrum_filted_pnext,inv=.False.))
        deallocate(spectrum_ppre, spectrum_p, spectrum_pnext)
        deallocate(spectrum_filted_ppre, spectrum_filted_p, spectrum_filted_pnext)
    end subroutine fft_gassian_filt2

    subroutine build_mask(mask, f, it, freq_cut, h_in, sigma_filter_in)
        type(t_field)        :: f
        real, intent(in),  optional      :: h_in, sigma_filter_in, freq_cut
        real, allocatable, intent(out)   :: mask(:,:,:) 

        real                :: h, sigma_filter, vp_min, k_cut
        integer             :: nz, nx, i, j, nz2, nx2
        integer             :: idx_x, idx_z
        real, allocatable   :: kx(:), kz(:), kx_grid(:,:), kz_grid(:,:), K(:,:)
        
        ! fp     = shot%fpeak
        vp_min = cb%velmin
        ! if(time_dir>0) then
        ifz=f%bloom(1,it)
        if(m%is_freesurface) ifz=max(ifz,1)
        ilz=f%bloom(2,it)
        ifx=f%bloom(3,it)
        ilx=f%bloom(4,it)
        nz = ilz-ifz+1
        nx = ilx-ifx+1
        ! else
        !     nz = size(f%p,1)
        !     nx = size(f%p,2)
        ! end if
        
        ! print *,'nz=',nz,'nx=',nx

        if (present(h_in)) then
            h = h_in
        else
            h = 12.5
        end if

        if (present(sigma_filter_in)) then
            sigma_filter = sigma_filter_in
        else
            sigma_filter = 2
        end if

        allocate(kx(nx), kz(nz), kx_grid(nz,nx), kz_grid(nz,nx), K(nz,nx))
        allocate(mask(nz,nx,1))
        ! allocate(kx(ilx-ifx+1), kz(ilz-ifz+1), kx_grid(ilz-ifz+1,ilx-ifx+1), kz_grid(ilz-ifz+1,ilx-ifx+1), K(nilz-ifz+1,ilx-ifx+1))
        ! allocate(mask(ilz-ifz+1,ilx-ifx+1,1))
        ! print *,'For mask:ifz ilz ifx ilx', ifz,ilz,ifx,ilx, 'K shape:',size(K,1),size(K,2), 'nz nx=',nz,nx

        do i = 1, nx/2
            kx(i) = 2.0*r_pi*real(i-1)/(h*real(nx))    ! [0 ... nx/2-1]
        end do
        do i = nx/2+1, nx
            kx(i) = 2.0d0*r_pi*real(i-nx-1)/(h*real(nx)) ! [-nx/2 ... -1]
        end do

        ! ---- kz 向量 ----
        do i = 1, nz/2
            kz(i) = 2.0d0*r_pi*real(i-1)/(h*real(nz))
        end do
        do i = nz/2+1, nz
            kz(i) = 2.0d0*r_pi*real(i-nz-1)/(h*real(nz))
        end do

        do j = 1, nz
        do i = 1, nx
            kx_grid(j,i) = kx(i)
            kz_grid(j,i) = kz(j)
            K(j,i) = sqrt(kx_grid(j,i)**2 + kz_grid(j,i)**2)
        end do
        end do


        ! 截止波数（单位 rad/m）
        if(present(freq_cut)) then
            k_cut = 2.0*r_pi*freq_cut / vp_min
        else
            k_cut = 2.0*r_pi*shot%fmax / vp_min
        endif

        where (K <= k_cut)
            mask(:,:,1) = 1.0
        elsewhere
            mask(:,:,1) = exp( - ((K - k_cut)/(sigma_filter*k_cut))**2 )
        end where

        deallocate(kx, kz, kx_grid, kz_grid, K)
    end subroutine build_mask

    subroutine build_mask2(mask, f, it, freq_cut, h_in, sigma_filter_in)
        type(t_field)        :: f
        real, intent(in),  optional      :: h_in, sigma_filter_in, freq_cut
        real, allocatable, intent(out)   :: mask(:,:,:) 

        real                :: h, sigma_filter, vp_min, k_cut
        integer             :: nz, nx, i, j, nz2, nx2
        integer             :: idx_x, idx_z
        real, allocatable   :: kx(:), kz(:), kx_grid(:,:), kz_grid(:,:), K(:,:)
        
        ! fp     = shot%fpeak
        vp_min = cb%velmin
        ! if(time_dir>0) then
        ifz=f%bloom(1,it)
        if(m%is_freesurface) ifz=max(ifz,1)
        ilz=f%bloom(2,it)
        ifx=f%bloom(3,it)
        ilx=f%bloom(4,it)
        nz = ilz-ifz+1
        nx = ilx-ifx+1
        ! else
        !     nz = size(f%p,1)
        !     nx = size(f%p,2)
        ! end if
        
        ! print *,'nz=',nz,'nx=',nx

        if (present(h_in)) then
            h = h_in
        else
            h = 12.5
        end if

        if (present(sigma_filter_in)) then
            sigma_filter = sigma_filter_in
        else
            sigma_filter = 2
        end if

        allocate(kx(nx), kz(nz), kx_grid(nz,nx), kz_grid(nz,nx), K(nz,nx))
        allocate(mask(nz,nx,1))
        ! allocate(kx(ilx-ifx+1), kz(ilz-ifz+1), kx_grid(ilz-ifz+1,ilx-ifx+1), kz_grid(ilz-ifz+1,ilx-ifx+1), K(nilz-ifz+1,ilx-ifx+1))
        ! allocate(mask(ilz-ifz+1,ilx-ifx+1,1))
        ! print *,'For mask:ifz ilz ifx ilx', ifz,ilz,ifx,ilx, 'K shape:',size(K,1),size(K,2), 'nz nx=',nz,nx

        do i = 1, nx/2
            kx(i) = 2.0*r_pi*real(i-1)/(h*real(nx))    ! [0 ... nx/2-1]
        end do
        do i = nx/2+1, nx
            kx(i) = 2.0d0*r_pi*real(i-nx-1)/(h*real(nx)) ! [-nx/2 ... -1]
        end do

        ! ---- kz 向量 ----
        do i = 1, nz/2
            kz(i) = 2.0d0*r_pi*real(i-1)/(h*real(nz))
        end do
        do i = nz/2+1, nz
            kz(i) = 2.0d0*r_pi*real(i-nz-1)/(h*real(nz))
        end do

        do j = 1, nz
        do i = 1, nx
            kx_grid(j,i) = kx(i)
            kz_grid(j,i) = kz(j)
            K(j,i) = sqrt(kx_grid(j,i)**2 + kz_grid(j,i)**2)
        end do
        end do


        ! 截止波数（单位 rad/m）
        ! if(present(freq_cut)) then
        !     k_cut = 2.0*r_pi*freq_cut / vp_min
        ! else
        k_cut = 2.0*r_pi*shot%fmax / vp_min
        ! endif

        where (K <= k_cut)
            mask(:,:,1) = 1.0
        elsewhere
            mask(:,:,1) = exp( - ((K - k_cut)/(sigma_filter*k_cut))**2 )
        end where

        deallocate(kx, kz, kx_grid, kz_grid, K)
    end subroutine build_mask2

    ! subroutine build_mask(mask, f, h_in, sigma_filter_in)
    !     type(t_field), intent(in)        :: f
    !     real, intent(in),  optional      :: h_in, sigma_filter_in
    !     real, allocatable, intent(out)   :: mask(:,:,:) 

    !     ! 局部
    !     real                :: h, sigma_filter, fp, vp_min, k_cut!, r_pi
    !     integer             :: nz, nx, i, j
    !     integer             :: idx_x, idx_z
    !     real, allocatable   :: kx(:), kz(:), kx_grid(:,:), kz_grid(:,:), K(:,:)

    !     r_pi = acos(-1.0)
    !     nz = size(f%p,1)
    !     nx = size(f%p,2)

    !     ! 这两个来自你的外部环境；确保在本模块 use 到它们
    !     fp     = shot%fpeak
    !     vp_min = cb%velmin

    !     ifz=f%bloom(1,it)
    !     if(m%is_freesurface) ifz=max(ifz,1)
    !     ilz=f%bloom(2,it)
    !     ifx=f%bloom(3,it)
    !     ilx=f%bloom(4,it)

    !     if (present(h_in)) then
    !         h = h_in
    !     else
    !         h = 12.5
    !     end if

    !     if (present(sigma_filter_in)) then
    !         sigma_filter = sigma_filter_in
    !     else
    !         sigma_filter = 0.5
    !     end if
    !     print *,'22222222222222'
    !     allocate(kx(ilx-ifx+1), kz(ilz-ifz+1), kx_grid(ilz-ifz+1,ilx-ifx+1), kz_grid(ilz-ifz+1,ilx-ifx+1), K(ilz-ifz+1,ilx-ifx+1))
    !     allocate(mask(ilz-ifz+1,ilx-ifx+1,1))
    !     print *,'333333333333'
    !     ! ---- 关键修正：环形折返的索引（奇偶长度都正确）----
    !     ! ---- kx 向量 ----
    !     do i = ifx, (ilx-ifx+1)/2
    !         kx(i) = 2.0*r_pi*real(i-1)/(h*real((ilx-ifx+1)))    ! [0 ... nx/2-1]
    !     end do
    !     do i = (ilx-ifx+1)/2+1, (ilx-ifx+1)
    !         kx(i) = 2.0d0*r_pi*real(i-(ilx-ifx+1)-1)/(h*real(ilx-ifx+1)) ! [-nx/2 ... -1]
    !     end do

    !     ! ---- kz 向量 ----
    !     do i = ifz, (ilz-ifz+1)/2
    !         kz(i) = 2.0d0*r_pi*real(i-1)/(h*real(ilz-ifz+1))
    !     end do
    !     do i = (ilz-ifz+1)/2+1, (ilz-ifz+1)
    !         kz(i) = 2.0d0*r_pi*real(i-(ilz-ifz+1)-1)/(h*real(ilz-ifz+1))
    !     end do

    !     do j = ifz, (ilz-ifz+1)
    !     do i = ifx, (ilx-ifx+1)
    !         kx_grid(j,i) = kx(i)
    !         kz_grid(j,i) = kz(j)
    !         K(j,i) = sqrt(kx_grid(j,i)**2 + kz_grid(j,i)**2)
    !     end do
    !     end do


    !     ! 截止波数（单位 rad/m）
    !     k_cut = 2.0*r_pi*fp / vp_min * 2

    !     where (K <= k_cut)
    !         mask(:,:,1) = 1.0
    !     elsewhere
    !         mask(:,:,1) = exp( - ((K - k_cut)/(sigma_filter*k_cut))**2 )
    !     end where

    !     deallocate(kx, kz, kx_grid, kz_grid, K)
    ! end subroutine build_mask

end
