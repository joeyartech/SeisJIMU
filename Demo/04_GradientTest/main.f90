program main
use m_System
use m_Modeling
use m_Kernel
use m_linesearcher

!     type(t_checkpoint) :: chp_shls, chp_qp
    type(t_querypoint),target :: qp0

    !mpiworld lives in t_mpienv
    call mpiworld%init(name='MPIWorld')

    call hud('======================================'//s_NL// &
             '       WELCOME TO SeisJIMU PFEI       '//s_NL// &
             '======================================')

    call setup%init
    call sysio_init
    
    if(.not. setup%exist) then
        if(mpiworld%is_master) then
            call hud('No input setup file given. Print manual..')
!             call print_manual
            call ppg%print_info
        endif
        call mpiworld%final
        stop
    endif

    !print propagator info
    call ppg%print_info

    ! !estimate required memory
    ! call m%estim_RAM
    ! call cb%estim_RAM
    ! call sfield%estim_RAM
    ! call rfield%estim_RAM
    
    !checkpoint
!     call checkpoint_init

    !model
    call m%init
    call m%read
    call ppg%check_model

    !shotlist
    call shls%read_from_data
    call shls%build
!     call chp_shls%init('FWI_shotlist_gradient',oif_fuse=.true.)
!     if(.not.shls%is_registered(chp_shls,'sampled_shots')) then
        call shls%sample
!         call shls%register(chp_shls,'sampled_shots')
!     endif
    call shls%assign

    !if preconditioner needs energy terms
    if(index(setup%get_str('PRECONDITIONING','PRECO'),'energy')>0) then
        ppg%if_compute_engy=.true.
    endif

    !parametrizer
    call param%init

    !current querypoint
    call qp0%init('qp0')

    !objective function and gradient
    call fobj%init
!     call chp_qp%init('FWI_querypoint_gradient')
!     if(.not.qp0%is_registered(chp_qp)) then
        call fobj%eval(qp0,oif_update_m=.false.)
!         call qp0%register(chp_qp)
!     endif

    call sysio_write('correlation_gradient',correlation_gradient,m%n*ppg%ngrad)
    
    call sysio_write('qp0%g',qp0%g,size(qp0%g))
    call sysio_write('qp0%pg',qp0%pg,size(qp0%pg))
    
    !linesearch init and scale
    call ls%init
    call ls%scale(qp0)
    
    call hud('qp0%f, ║g║₁ = '//num2str(qp0%f)//', '//num2str(sum(abs(qp0%g))))
    
    !if just estimate the wavelet or compute the gradient then this is it.
    if(setup%get_str('JOB')=='gradient') then
        call mpiworld%final
        stop
    endif
    
    call optimizer_init_loop(qp0)
        
    call mpiworld%final

end
