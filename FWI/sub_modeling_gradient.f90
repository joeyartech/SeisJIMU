subroutine modeling_gradient_ip
end subroutine


subroutine modeling_gradient_vp
use mpi
use m_System
use m_Modeling
use m_weighter
use m_Lpnorm
use m_Envnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    logical,save :: is_first_in=.true.

    character(:),allocatable :: s_update_wavelet
    type(t_weighter) :: wei_wl

    type(t_field) :: fld_u,fld_a
    type(t_correlate) :: a_star_u

    character(:),allocatable :: s_dnorm
    real,dimension(:,:),allocatable :: tmp,Eobs
    
    ! type :: t_S
    !     real,dimension(:),allocatable :: scale
    ! end type
    ! type(t_S),dimension(:),allocatable,save :: S
    ! real,dimension(:,:),allocatable :: tmp_dsyn, gmwindow

    ! if(is_first_in) allocate(S(shls%nshots_per_processor)) !then can NOT randomly sample shots..


    !misfit
    fobj%misfit=0.

    if(.not.allocated(s_dnorm)) s_dnorm=setup%get_str('DATA_NORM','DNORM',o_default='L2sq')

    call alloc(correlate_gradient,m%nz,m%nx,m%ny,ppg%ngrad)
    
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        if(s_dnorm=='Envsq'.or.s_dnorm=='Env2sq') then
            call alloc(Eobs,shot%nt,shot%nrcv)
            call hilbert_envelope(shot%dobs,Eobs,shot%nt,shot%nrcv)
            call shot%write('Eobs_',Eobs)
        endif

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project
        
        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer

        call hud('----  Solving Au=s  ----')
        call ppg%init_field(fld_u, name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire; call shot%write('Ru_',shot%dsyn)


        if(setup%get_str('JOB')=='forward') cycle

        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)

        ! s_update_wavelet=setup%get_str('UPDATE_WAVELET')
        ! if(s_update_wavelet/='') then
        !     call hud('----  Update Wavelet  ----')    
        !     call wei_wl%update(o_suffix='_4WAVELET')
        !     call shot%update_wavelet(wei_wl%weight) !call gradient_matchfilter_data    
        !     call shot%write('updated_Ru_',shot%dsyn)
        !     call suformat_write('updated_wavelet_'//shot%sindex,shot%wavelet,shot%nt,1,shot%dt)
        ! endif

        call hud('----  Computing obj func & dadj  ----')

            call hud('Using DNORM '//s_dnorm)
            select case (s_dnorm)
                case ('L2sq')
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
                    
                call kernel_L2sq(shot%dadj)


                !check if S is changing..
    !            if(shot%index==1)   print*, 'on '//shot%sindex,i,S(i)%scale
    !            if(shot%index==112) print*, 'on '//shot%sindex,i,S(i)%scale

                case('Envsq')
                fobj%misfit = fobj%misfit &
                    + Envsq(0.5, shot%nt, shot%nrcv, wei%weight, shot%dsyn, Eobs, shot%dt)
                call kernel_Envsq(shot%dadj,shot%nt,shot%nrcv)


                case('Qsq')
                fobj%misfit = fobj%misfit &
                    + Qsq(0.5, shot%nt, shot%nrcv, wei%weight, shot%dsyn, shot%dobs, shot%dt)
                call kernel_Qsq(shot%dadj,shot%nt,shot%nrcv,shot%dt)
                

                case default
                call error('No DNORM specified!')

            end select

            call shot%write('dadj_',shot%dadj)

        call hud('----  Solving adjoint eqn & xcorrelate  ----')
        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.); call fld_a%ignite
        call ppg%init_correlate(a_star_u,'a_star_u')
        call ppg%adjoint(fld_a,fld_u,a_star_u)

        call hud('----  Assemble  ----')
        call ppg%assemble(a_star_u)

        call hud('---------------------------------')

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    
    !allreduce misfit values
    call mpi_allreduce(mpi_in_place, [fobj%misfit], 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked misfit '//num2str(fobj%misfit))

    fobj%dnorms=fobj%misfit

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
    
    !scale by shotlist
    call shls%scale(1,o_from_sampled=[fobj%misfit])
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)


    !write correlate
    if(mpiworld%is_master) then
        call a_star_u%write
    endif

    !allreduce energy, gradient
    ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale by shotlist
    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)

    if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)

    call mpiworld%barrier
    
    is_first_in=.false.

end subroutine
