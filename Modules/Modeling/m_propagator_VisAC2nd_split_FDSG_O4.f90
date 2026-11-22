module m_propagator
use m_System
use m_hicks, only : hicks_r
use m_hilbert
use m_resampler
use m_model
use m_shot
use m_computebox
use m_freesurface
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

    logical :: if_freesurface_ppg

    ! logical :: is_Q_attenuation
    character(:),allocatable :: stablize_method
    real split_fpeak
    
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
            'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
            'Basic gradients: gbuo(wait), gikpa, giqp'

        integer :: nbndlayer=max(1,hicks_r) !minimum absorbing layer thickness
        integer :: ngrad=3 !number of basic gradients
        integer :: nimag=3 !number of basic images
        !integer :: nengy=1 !number of energy terms

        logical :: if_compute_engy=.false.

        !local models shared between fields
        real,dimension(:,:,:),allocatable :: buoz, buox, buoy, kpa!, qp
        complex,dimension(:,:,:),allocatable :: invC2,  C1n,  C0n
        complex,dimension(:,:,:),allocatable :: invC2H, C1nH, C0nH
        
        complex,dimension(:,:,:), allocatable :: Ap, Am, Bcoef
        real,dimension(:,:,:)   , allocatable :: damp,dampH,scale
        complex,dimension(:,:,:), allocatable :: w_disp

        complex :: dC2dQ
        ! real :: dC2dQ
        complex :: dC1dQ, dC0dQ

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
        
        procedure :: inject_pressure
        procedure :: update_pressure
        procedure :: evolve_pressure
        procedure :: extract

        ! procedure :: dft2d
        ! procedure :: idft2d
        ! procedure :: fftshift2d
        ! procedure :: ifftshift2d
        ! procedure :: fft_filter

        final :: final

    end type

    type(t_propagator),public :: ppg

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    logical :: if_record_adjseismo=.false.

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

        if(.not.setup%get_bool('IS_Q_DISPERSION',o_default='T')) then
            visac_a = 0.
            visac_b = 0.
            visac_d = 0.
        endif
        if(.not.setup%get_bool('IS_Q_DISSIPATION',o_default='T')) then
            visac_e = 0.
        endif
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
        real,   dimension(:,:,:),allocatable :: C2r   !【split】Re(C2)，Step1 实刚度（前向/伴随共用）
        complex :: Cq, disc, w1, w2
        complex :: disch, w1h, w2h                  !【split】伴随求根临时量

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

        !coef
        C2 = 1 -2*visac_a/r_pi/cb%qp -c_i*visac_e/cb%qp; self%invC2 = 1/C2
        C1 =    2*visac_b/r_pi/cb%qp;                    self%C1n = C1*self%invC2
        C0 =    2*visac_d/r_pi/cb%qp;                    self%C0n = C0*self%invC2

        C2 = 1 -2*visac_a/r_pi/cb%qp +c_i*visac_e/cb%qp; self%invC2H = 1/C2
        C1 =    2*visac_b/r_pi/cb%qp;                    self%C1nH= C1*self%invC2H
        C0 =    2*visac_d/r_pi/cb%qp;                    self%C0nH= C0*self%invC2H

        stablize_method=setup%get_str('PPG_STABLIZE_METHOD',o_default='split')
        split_fpeak=setup%get_real('SPLIT_FPEAK',o_default=num2str(shot%fpeak))

        ! # ===========================================================================
        ! # 2) 逐点衰减 d(x,z)：完整色散在主频 w_ref 处的物理支
        ! #    k_ref = w_ref / v  =>  v^2 k_ref^2 = w_ref^2（与空间无关）
        ! # ===========================================================================        
        C2 = 1 -2*visac_a/r_pi/cb%qp -c_i*visac_e/cb%qp
        C1 =    2*visac_b/r_pi/cb%qp
        C0 =    2*visac_d/r_pi/cb%qp

        ! ---- CFL：Q_min 决定 C2r_min ----
        dt_cfl_out = m%dz * sqrt(minval(real(C2))) / (sqrt(2.) * maxval(m%vp))
        write(*,*) 'Q range = [', minval(m%qp), ', ', maxval(m%qp), '];'
        write(*,*) 'C2r_min = ', minval(real(C2)), '; vp_max = ', maxval(m%vp), ' m/s'
        write(*,*) 'dtCFL = ', dt_cfl_out, ' s;  dt = ', self%dt, ' s;  dt/dtCFL = ', self%dt/dt_cfl_out
        if (self%dt > dt_cfl_out) call error('dt 超过有效 CFL，本方案不保证稳定，请减小 dt')

        ! ---- 逐点衰减 ----
        call alloc(self%damp,  cb%nz,cb%nx,1)
        call alloc(self%dampH, cb%nz,cb%nx,1)
        call alloc(self%w_disp,cb%nz,cb%nx,1)
        call alloc(C2r,        cb%nz,cb%nx,1)

        !【split】Step1 用的实刚度 C2r=Re(C2)。C2 的虚部(-i/Q)全部交给 Step2 的 d 补回。
        C2r = real(C2)

        ! 参考频率 w_ref=2π·split_fpeak：v²k_ref²=w_ref²，与空间位置无关，故 Cq=C0+w_ref²。
        ! 本模块 PDE 为 C2·P̈ - i·C1·Ṗ + C0·P = κ∇·b∇P，一阶项已带 -i，
        ! 故 C1=2b/(πQ) 是实数；色散 C2·w² + C1·w - Cq = 0
        !   =>  disc=sqrt(C1²+4C2Cq)，物理支 w=(C1+disc)/(2C2) (Re>0, Im>0)。
        ! 关键：不是 (-C1±disc)/(2C2)，也不是 notebook 里把 C1 当初纯虚数的
        ! (iC1±sqrt(-C1²+4C2Cq))/(2C2)。后两者代入本模块实数 C1 会让 Im(w)<0 被
        ! max(.,0) 截成 d=0 —— "没跑出耗散"的根因。数值核对(Q=36,f=20Hz): d=1.9e-3。
        wref2 = (2.0*r_pi*split_fpeak)**2

        do j = 1, cb%nx
        do i = 1, cb%nz
            Cq   = C0(i,j,1) + wref2

            disc = sqrt(C1(i,j,1)**2 + 4.*C2(i,j,1)*Cq)
            w1   = (C1(i,j,1) + disc) / (2.*C2(i,j,1))
            w2   = (C1(i,j,1) - disc) / (2.*C2(i,j,1))
            if (real(w1) >= real(w2)) then
                self%w_disp(i,j,1) = w1
            else
                self%w_disp(i,j,1) = w2
            end if
            self%damp(i,j,1) = max(aimag(self%w_disp(i,j,1)),0.) *self%dt

            ! 伴随：C2H=conjg(C2)=C2r+i/Q(虚部 +i/Q)，色散 C2H·w² + C1·w - Cq=0。
            ! 因 C1,Cq 实数，wH=conjg(w)，Im(wH)=-Im(w)<0；倒传衰减率取
            ! d_H=-Im(wH)dt=Im(w)dt(与正演相等)。若误取 +Im(wH)dt 会得 0、耗散丢失。
            disch = sqrt(C1(i,j,1)**2 + 4.*conjg(C2(i,j,1))*Cq)
            w1h   = (C1(i,j,1) + disch) / (2.*conjg(C2(i,j,1)))
            w2h   = (C1(i,j,1) - disch) / (2.*conjg(C2(i,j,1)))
            if (real(w1h) >= real(w2h)) then
                self%dampH(i,j,1) = max(-aimag(w1h),0.) *self%dt
            else
                self%dampH(i,j,1) = max(-aimag(w2h),0.) *self%dt
            end if

        end do
        end do

        ! ---- Step1 差分系数（逐点数组，实刚度 C2r + 一阶项照模块格式带 c_i）----
        ! Ap=1-i·C1·dt/(2·C2r) 作分母；Am=1+i·C1·dt/(2·C2r) 乘 Uprev；
        ! Bcoef=2-C0·dt²/C2r；scale=dt²/C2r。对应推进 U*=(Bcoef·U-Am·Uprev+scale·κ·Lap)/Ap。
        call alloc(self%Ap,   cb%nz,cb%nx,1)
        call alloc(self%Am,   cb%nz,cb%nx,1)
        call alloc(self%Bcoef,cb%nz,cb%nx,1)
        call alloc(self%scale,cb%nz,cb%nx,1)
        self%Ap    = 1. - c_i*C1*self%dt/(2.*C2r)
        self%Am    = 1. + c_i*C1*self%dt/(2.*C2r)
        self%Bcoef = 2. - C0*self%dt**2/C2r
        self%scale = self%dt**2 / C2r


        ! self%dC2dQ = 2*visac_a/r_pi - c_i*visac_e
        ! self%dC1dQ = 2*visac_b/r_pi
        ! self%dC0dQ = 2*visac_d/r_pi

        if_freesurface_ppg=setup%get_bool('FREE_SURFACE_PPG',o_default='T')

        deallocate(C2,C1,C0,C2r) !save some RAM

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

    end subroutine

    subroutine init_correlate(self,corr,name)
        class(t_propagator) :: self
        type(t_correlate) :: corr
        character(*) :: name
        
        corr%name=name

        ! if(name(1:1)=='g') then !gradient components
            call alloc(corr%gbuo, m%nz,m%nx,m%ny)
            call alloc(corr%gikpa,m%nz,m%nx,m%ny)
            call alloc(corr%giqp,  m%nz,m%nx,m%ny)
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
            call correlate_assemble(corr%gbuo,  correlate_gradient(:,:,:,1))
            call correlate_assemble(corr%gikpa, correlate_gradient(:,:,:,2))
            call correlate_assemble(corr%giqp,  correlate_gradient(:,:,:,3))

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
        real,parameter :: time_dir=1. !time direction

        !seismo
        call alloc(fld_reU%seismo,shot%nrcv,self%nt)
        call alloc(fld_imU%seismo,shot%nrcv,self%nt)

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

            !step 4: update pressure
            call cpu_time(tic)
            call self%update_pressure(fld_reU,fld_imU,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

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

            !image: rf%p^it star sf%p^it
            if(mod(it,irdt)==0) then
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

        !postprocess
        call cross_correlate_postprocess(A_star_U)
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
            ! write(*,*) 'Elapsed time to compute Poynting vectors ',tt3/mpiworld%max_threads
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
            
            wl=time_dir*f%wavelet(1,it)*wavelet_scaler

            if(if_hicks) then
                select case (shot%src%comp)
                case ('p') !explosion
                    f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) + wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%src%interp_coef_full
                    
                case ('dpdz') !vertical force
                    f%p(ifz+1:ilz+1,ifx:ilx,ify:ily) = f%p(ifz+1:ilz+1,ifx:ilx,ify:ily) + wl*self%kpa(ifz+1:ilz+1,ifx:ilx,ify:ily)*inv_2dz*shot%src%interp_coef
                    f%p(ifz-1:ilz-1,ifx:ilx,ify:ily) = f%p(ifz-1:ilz-1,ifx:ilx,ify:ily) - wl*self%kpa(ifz-1:ilz-1,ifx:ilx,ify:ily)*inv_2dz*shot%src%interp_coef
                
                endselect
            
            else
                select case (shot%src%comp)
                case ('p') !explosion
                    f%p(iz,ix,iy)                  = f%p(iz,ix,iy)                + wl*self%kpa(iz,ix,iy)
                    
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
                    f%p(ifz:ilz,ifx:ilx,ify:ily) = f%p(ifz:ilz,ifx:ilx,ify:ily) +wl*self%kpa(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef_full

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
        complex,dimension(:,:,:),allocatable :: Ustar
                
        ! !necessary after computing the secondary source
        ! f%lap=0.

        ifz=f_re%bloom(1,it)
        if(m%is_freesurface) ifz=max(ifz,1)
        ilz=f_re%bloom(2,it)
        ifx=f_re%bloom(3,it)
        ilx=f_re%bloom(4,it)
        
        if(m%is_cubic) then

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
            
        if(if_freesurface_ppg) then
            ! call freesurface_stress(f_re%lap)
            ! call freesurface_stress(f_im%lap)
            call freesurface_stress(f_re%p)
            call freesurface_stress(f_im%p)
        endif

        if(.not. f_re%is_adjoint) then !PDE
            if(time_dir>0. ) then !forward in time
                Uprev = cmplx(f_re%p_prev,f_im%p_prev)
                U     = cmplx(f_re%p     ,f_im%p     )
                Lap   = cmplx(f_re%lap   ,f_im%lap   )

                !【split 正演】Step1: 实刚度 C2r 推进得 U*（临界稳定，无耗散）
                !   Ap·U* = Bcoef·U - Am·Uprev + scale·κ·Lap
                !   （lap=∇·(1/ρ∇P)，方程 (1/κ)(...)=∇·(1/ρ∇P) 故 Lap 必乘 κ=kpa）
                Ustar = ( self%Bcoef*U - self%Am*Uprev + self%scale*self%kpa*Lap ) / self%Ap
                ! Step2: 逐点二阶分裂衰减（纯空间域，无需 FFT/DST）
                Unext = ( Ustar + self%damp*Uprev ) / (1 + self%damp)

                f_re%p_next =  real(Unext)
                f_im%p_next = aimag(Unext)

            else !backward in time —— 入射场反向重建（必须用正向 split 的精确逆）
                Unext = cmplx(f_re%p_next,f_im%p_next)
                U     = cmplx(f_re%p     ,f_im%p     )
                Lap   = cmplx(f_re%lap   ,f_im%lap   )

                !【split 反向精确逆】正演两步联立：
                !   Ap·U* = Bcoef·U - Am·Uprev + scale·κ·Lap ;  Unext·(1+d)=U*+d·Uprev
                ! 消去 U* 反解 Uprev：
                !   Uprev = [Bcoef·U + scale·κ·Lap - (1+d)·Ap·Unext] / (Am - d·Ap)
                ! 注意不能用"把 d 对称写一遍"的伪逆，否则重建不精确（往返误差大）。
                Uprev = ( self%Bcoef*U + self%scale*self%kpa*Lap - (1.+self%damp)*self%Ap*Unext ) &
                      / ( self%Am - self%damp*self%Ap )

                f_re%p_prev =  real(Uprev)
                f_im%p_prev = aimag(Uprev)

            endif

        else !Adjoint, backward in time
                Unext = cmplx(f_re%p_next,f_im%p_next)
                U     = cmplx(f_re%p     ,f_im%p     )
                Lap   = cmplx(f_re%lap   ,f_im%lap   )

                !【split 伴随倒传】C2H=C2r+i/Q（虚部与正演相反）。倒传同样有 |g|>1 的
                ! 不稳定根，必须分裂。Step1 用实刚度 C2r（与正演共用 Ap/Am/Bcoef/scale），
                ! 但一阶导数项按本模块伴随格式取相反的 i 号 —— 分母用 Am、Unext 系数用 Ap：
                !   Am·U* = Bcoef·U + scale·κ·Lap - Ap·Unext
                ! Step2 用伴随专属 dampH=-Im(wH)dt 做对称分裂（倒传不能用精确逆，会发散）。
                Ustar = ( self%Bcoef*U + self%scale*self%kpa*Lap - self%Ap*Unext ) / self%Am
                Uprev = ( Ustar + self%dampH*Unext ) / (1 + self%dampH)

                f_re%p_prev =  real(Uprev)
                f_im%p_prev = aimag(Uprev)

        endif

        deallocate(Uprev, U, Unext, Lap) !save some RAM

    end subroutine


    subroutine fd_freesurface_stresses(p)
        ! real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: p
        real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: p

        p(1,:,:)=0.
        p(0,:,:)=-p(2,:,:)

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

        else !in reverse time
            tmp=>f%p_next
            f%p_next=>f%p
            f%p     =>f%p_prev
            f%p_prev=>tmp

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
                ify=shot%rcv(i)%ify-cb%ioy+1; iy=shot%rcv(i)%iy-cb%ioy+1; ily=shot%rcv(i)%ily-cb%ioy+1

                if(if_hicks) then
                    select case (shot%rcv(i)%comp)
                    case ('p')
                        f%seismo(i,it)=sum(f%p(ifz:ilz,ifx:ilx,ify:ily)*shot%rcv(i)%interp_coef_full)

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
                    f%seismo(1,it)=sum(f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coef_full)
                    
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
        !【split】释放分裂法新增的逐点系数/衰减数组
        if(allocated(self%Ap))    deallocate(self%Ap, self%Am, self%Bcoef, self%scale)
        if(allocated(self%damp))  deallocate(self%damp, self%dampH, self%w_disp)
    end subroutine


    !========= gradient, imaging or other correlations ===================
    !For gradient:
    !Kₘ<a|Au> = Kₘ<a|ϰ∂ₜ²u - ∇·b∇u>
    !for ϰ: Kₘ<a|Au> = ∫ a ∂ₜ²u dt =-∫ ∂ₜa ∂ₜu dt, or = ∫ a κ∇·b∇u dt 
    !for b: Kₘ<a|Au> = -Kₘ<a|∇·b∇u> = ∫ ∇a·∇u dt
    !for qp: Kₘ<a|Au> = ∫ a 1/κ(dC2dQ ∂ₜₜu +dC1dQ ∂ₜu +dC0dQ u) dt 

    subroutine cross_correlate_gradient(reA,imA, reU,imU, corr,it)
        type(t_field), intent(in) :: reA, imA, reU, imU
        type(t_correlate) :: corr

        logical,dimension(:,:,:),allocatable :: mask
        complex,dimension(:,:,:),allocatable :: Uprev, U, Unext, A

        !nonzero only when sf touches rf
        ifz=max(reU%bloom(1,it),reA%bloom(1,it),2)
        ilz=min(reU%bloom(2,it),reA%bloom(2,it),cb%mz)
        ifx=max(reU%bloom(3,it),reA%bloom(3,it),1)
        ilx=min(reU%bloom(4,it),reA%bloom(4,it),cb%mx)
        ify=max(reU%bloom(5,it),reA%bloom(5,it),1)
        ily=min(reU%bloom(6,it),reA%bloom(6,it),cb%my)

        call alloc(mask,m%nz,m%nx,m%ny)
        do ix = ifx,ilx
        do iz = ifz,ilz
            mask(iz,ix,1) = .true.
        enddo
        enddo

        where(mask)

            corr%gikpa = corr%gikpa + reA%p(1:m%nz,1:m%nx,1:m%ny)*reU%lap(1:m%nz,1:m%nx,1:m%ny) !should use this one because we just use the real part in the misfit function, only the re-re term is needed here.

            corr%giqp = corr%giqp + reA%p(1:m%nz,1:m%nx,1:m%ny)/ppg%kpa(1:m%nz,1:m%nx,1:m%ny)* &
                    real(ppg%dC2dQ*(reU%p_next(1:m%nz,1:m%nx,1:m%ny) + reU%p_prev(1:m%nz,1:m%nx,1:m%ny)- 2*reU%p(1:m%nz,1:m%nx,1:m%ny))/dt2 &
                        -ppg%dC1dQ*(reU%p_next(1:m%nz,1:m%nx,1:m%ny) - reU%p_prev(1:m%nz,1:m%nx,1:m%ny))/(2*ppg%dt) &
                        +ppg%dC0dQ* reU%p(1:m%nz,1:m%nx,1:m%ny) &
                        )

        end where

    end subroutine

    subroutine cross_correlate_postprocess(corr)
        type(t_correlate) :: corr

        if(allocated(correlate_gradient)) then
            !scaling gradients by model parameters
            ! corr%grho = corr%grho / cb%rho(1:cb%mz,1:cb%mx,1:cb%my)
            corr%gikpa=corr%gikpa*ppg%kpa(1:m%nz,1:m%nx,1:m%ny)

            !remove singular point at the src position,
            !because we didn't consider src when deriving the gradient formula    
            iz=shot%src%iz-cb%ioz+1
            ix=shot%src%ix-cb%iox+1

            if(shot%src%comp=='dpdz') then !dipole source   
                corr%gikpa(iz-1,ix,1) = corr%gikpa(iz,ix,1)
                corr%gikpa(iz+1,ix,1) = corr%gikpa(iz,ix,1)

                corr%giqp(iz-1,ix,1) = corr%giqp(iz,ix,1)
                corr%giqp(iz+1,ix,1) = corr%giqp(iz,ix,1)
            
            elseif(shot%src%comp=='dpdx') then !dipole source
                corr%gikpa(iz,ix-1,1) = corr%gikpa(iz,ix,1)
                corr%gikpa(iz,ix+1,1) = corr%gikpa(iz,ix,1)

                corr%giqp(iz,ix-1,1) = corr%giqp(iz,ix,1)
                corr%giqp(iz,ix+1,1) = corr%giqp(iz,ix,1)
            
            else !point source
                !safeguards
                ifz=either(iz-2,iz,iz>=3); ilz=either(iz+2,iz,iz<=m%nz-2)
                ifx=either(ix-2,ix,ix>=3); ilx=either(ix+2,ix,ix<=m%nx-2)
                
                corr%gikpa(iz,ix,1)=0 !first remove otherwise will appear in the sum below
                corr%giqp (iz,ix,1)=0 !first remove otherwise will appear in the sum below
                ncells=size(corr%gikpa(ifz:ilz,ifx:ilx,1))-1
                corr%gikpa(iz,ix,1) = sum(corr%gikpa(ifz:ilz,ifx:ilx,1))/ncells
                corr%giqp (iz,ix,1) = sum(corr%giqp (ifz:ilz,ifx:ilx,1))/ncells
            
            endif

            !removing singular top boundary..
            ! corr%grho(1,:,:) = corr%grho(2,:,:)
            corr%gikpa(1,:,:) = corr%gikpa(2,:,:)
            corr%giqp (1,:,:) = corr%giqp (2,:,:)

        endif

        ! if(allocated(correlate_image)) then
        !     corr%ipp (1,:,:) = corr%ipp (2,:,:)
        !     corr%ibksc(1,:,:) = corr%ibksc(2,:,:)
        !     corr%ifwsc(1,:,:) = corr%ifwsc(2,:,:)
        ! endif

    end subroutine

    !========= Finite-Difference on flattened arrays ==================

    subroutine fd2d_flux(p,dp_dz,dp_dx,&
                         az,ax,&
                         buoz,buox,&
                         ifz,ilz,ifx,ilx)
        real,dimension(*) :: p,dp_dz,dp_dx
        real,dimension(*) :: az,ax
        real,dimension(*) :: buoz,buox

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

                az(iz_ix) = buoz(iz_ix)*dp_dz_
                ax(iz_ix) = buox(iz_ix)*dp_dx_

            enddo
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine

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

end
