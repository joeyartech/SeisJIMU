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
             '       WELCOME TO SeisJIMU FWI        '//s_NL// &
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
    
    call sysio_write('qp0%g',qp0%g,size(qp0%g))
    call sysio_write('qp0%pg',qp0%pg,size(qp0%pg))
    
    !linesearch init and scale
    call ls%init
    call ls%scale(qp0)
    
    call hud('qp0%f, ║g║₂² = '//num2str(qp0%f)//', '//num2str(norm2(qp0%g)))
    
    !if just estimate the wavelet or compute the gradient then this is it.
    if(setup%get_str('JOB')=='gradient') then
        call mpiworld%final
        stop
    endif
    
    call optimizer_init_loop(qp0)
        
    call mpiworld%final

end

!copy-paste from FWI/main.f90
subroutine modeling_gradient(is_fitting_data)!(oif_gradient)
use mpi
use m_System
use m_Modeling
use m_weighter
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse

    !logical,optional :: oif_gradient
    logical :: is_fitting_data

    logical,save :: is_first_in=.true.
    type(t_field) :: fld_u, fld_a
    character(:),allocatable :: update_wavelet

    is_fitting_data=.true.

    fobj%dnorms=0.
    fobj%xnorms=0.

    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor

        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project

        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer

        call ppg%init_field(fld_u,name='fld_u');  call fld_u%ignite
        
        !forward modeling
        call ppg%forward(fld_u);   call fld_u%acquire

        if(mpiworld%is_master) call shot%write('draw_',shot%dsyn)

        update_wavelet=setup%get_str('UPDATE_WAVELET',o_default='per shot')

        if(update_wavelet/='no') call shot%update_wavelet !call gradient_matchfilter_data
        
        !write synthetic data
        call shot%write('dsyn_',shot%dsyn)

        !data weighting
        call wei%update
        
        !Adjoint state method with Lagrangian formulation
        !to compute the gradient (L2 norm example)
        !C = ½║u║² = ½∫ (u-d)² δ(x-xr) dtdx³
        !∇ᵤC = (u-d)δ(x-xr)
        !L = C + <a|Au-s> ≐ C + <Aᴴa|u>
        !0 = ∇ᵤL = ∇ᵤC + Aᴴa => Aᴴa = -∇ᵤC
        !∇ₘL = a ∇ₘA u =: a★Du

        !objective function and adjoint source
        call fobj%stack_dnorms
        shot%dadj=-shot%dadj
        
        if(mpiworld%is_master) call fobj%print_dnorms('Shotloop-stacked','upto Shot#'//shot%sindex)
        ! if(either(oif_gradient,.true.,present(oif_gradient))) then
        
            !adjoint source
            if(update_wavelet/='no') call shot%update_adjsource
            
            call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)
            
            call fld_a%ignite
            
            !adjoint modeling
            call ppg%adjoint(fld_a,fld_u,oif_compute_grad=.true.)
            
            call cb%project_back
            
        ! endif

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    !collect global objective function value
    call mpi_allreduce(mpi_in_place, fobj%dnorms, fobj%n_dnorms, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    !scale by shotlist
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)

    call fobj%print_dnorms('Shotloop-stacked, shotlist-scaled, but yet linesearch-scaled','')

    ! if(either(oif_gradient,.true.,present(oif_gradient))) then
        !collect global gradient
        call mpi_allreduce(mpi_in_place, m%gradient,  m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        !scale
        call shls%scale(m%n*ppg%ngrad,o_from_sampled=m%gradient)

        if(ppg%if_compute_engy) then
            call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        endif
    ! endif

    call mpiworld%barrier

end subroutine
