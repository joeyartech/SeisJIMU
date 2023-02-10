subroutine modeling_gradient
use mpi
use m_System
use m_Modeling
use m_weighter
use m_Lpnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_E0, fld_dE, fld_F1, fld_F2
    type(t_correlation) :: F2_star_E0, F1_star_lapE0, F2_star_lapdE, gtDt, gtDs, gvp2
    ! real,dimension(:,:),allocatable :: Wdres
    !real,dimension(:,:,:),allocatable :: Ddt2
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff
    
    call gtDt%init('gtDt',shape123_from='model')
    ! call gtDs%init('gtDs',shape123_from='model')
    call gvp2%init('gvp2',shape123_from='model')
    
    !PFEI misfit
    fobj%misfit=0.    
    
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

        call ppg%init_field(fld_E0,name='fld_E0');    call fld_E0%ignite
        call ppg%init_field(fld_dE,name='fld_dE')

        !forward modeling
        !A E₀= s
        !AδE = DE₀
        call ppg%forward_scattering(fld_dE,fld_E0)
        call fld_dE%acquire; call shot%write('RdE_',shot%dsyn); shot%dsyn_aux=shot%dsyn
        call fld_E0%acquire; call shot%write('RE0_',shot%dsyn)

        if(setup%get_str('JOB')=='forward modeling') cycle

        if(index(setup%get_str('JOB'),'build tilD')>0) then
            call hud('--------------------------------')
            call hud('        Build tilD model         ')
            call hud('--------------------------------')

            call F2_star_E0%init('F2_star_E0',shape123_from='model') !F₂★E₀ for gtDt

            call wei%update!('_4IMAGING')
            fobj%misfit = fobj%misfit &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-(shot%dsyn+shot%dsyn_aux), shot%dt)
            
            call alloc(shot%dadj,shot%nt,shot%nrcv)
            call kernel_L2sq(shot%dadj)
            call shot%write('dadj_',shot%dadj)
            
            !adjoint modeling
            !AᴴF₂ = -RᴴN(E+δE)
            call ppg%init_field(fld_F2,name='fld_F2',ois_adjoint=.true.); call fld_F2%ignite
            call ppg%adjoint_tilD(fld_F2,fld_E0,F2_star_E0)
            call hud('---------------------------------')

            !call cb%project_back
            gtDt%sp_rp(cb%ioz:cb%ioz+cb%mz-1,&
                       cb%iox:cb%iox+cb%mx-1,&
                       cb%ioy:cb%ioy+cb%my-1,:) = &
            gtDt%sp_rp(cb%ioz:cb%ioz+cb%mz-1,&
                       cb%iox:cb%iox+cb%mx-1,&
                       cb%ioy:cb%ioy+cb%my-1,:) + F2_star_E0%sp_rp(:,:,:,:)
            ! gtDs%sp_rp(cb%ioz:cb%ioz+cb%mz-1,&
            !            cb%iox:cb%iox+cb%mx-1,&
            !            cb%ioy:cb%ioy+cb%my-1,shot%index) = F2_star_E0%sp_rp(:,:,:,1)

        endif
        
        if(index(setup%get_str('JOB'),'update velocity')>0) then
            call hud('---------------------------------')
            call hud('     Update velocity model       ')
            call hud('---------------------------------')

            call F1_star_lapE0%init('F1_star_lapE0',shape123_from='model') !F₁★∇²E₀ for gvp2 (one RE)
            call F2_star_lapdE%init('F2_star_lapdE',shape123_from='model') !F₂★∇²δE for gvp2 (the other RE)

            call wei%update
            fobj%misfit = fobj%misfit &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-(shot%dsyn+shot%dsyn_aux), shot%dt)

            call alloc(shot%dadj,shot%nt,shot%nrcv)
            call kernel_L2sq(shot%dadj)
            call shot%write('dadj_',shot%dadj)
            
            !adjoint modeling
            !AᴴF₂ =      +Rᴴ(d-(E₀+δE))
            !AᴴF₁ = -DF₂ +Rᴴ(d-(E₀+δE))
            call ppg%init_field(fld_F1,name='fld_F1',ois_adjoint=.true.)
            call ppg%init_field(fld_F2,name='fld_F2',ois_adjoint=.true.); call fld_F2%ignite
            call ppg%adjoint_vp2(fld_F1,fld_F2,fld_dE,fld_E0,F1_star_lapE0,F2_star_lapdE)
            call hud('---------------------------------')
    
            !project back    
            gvp2%sp_rp(cb%ioz:cb%ioz+cb%mz-1,&
                       cb%iox:cb%iox+cb%mx-1,&
                       cb%ioy:cb%ioy+cb%my-1,:) = &
            gvp2%sp_rp(cb%ioz:cb%ioz+cb%mz-1,&
                       cb%iox:cb%iox+cb%mx-1,&
                       cb%ioy:cb%ioy+cb%my-1,:) + F1_star_lapE0%sp_rp(:,:,:,:) + F2_star_lapdE%sp_rp(:,:,:,:)
            
        endif

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    if(setup%get_str('JOB')=='forward modeling') then
        call mpiworld%final
        stop
    endif

    !collect PFEI misfit values
    !no need to collect WPI misfit values
    call mpi_allreduce(mpi_in_place, [fobj%misfit], 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked PFEI_misfit '//num2str(fobj%misfit))

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
    
    !collect global correlations
    !call mpi_allreduce(mpi_in_place, gtilD%sp_rp, m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, gvp2%sp_rp,  m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale
    call shls%scale(1,o_from_sampled=[fobj%misfit])
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)
    call shls%scale(m%n*3,o_from_sampled=m%correlate)


    !call alloc(m%gradient,m%nz,m%nx,m%ny,1)

    call alloc(correlation_gradient,m%nz,m%nx,m%ny,1)

    if(index(setup%get_str('JOB'),'build tilD')>0) then
        if(mpiworld%is_master) then
            call gtDt%write
        endif
        correlation_gradient=correlation_gradient +gtDt%sp_rp
        !correlation_gradient=correlation_gradient +gtDs%sp_rp
    endif

    if(index(setup%get_str('JOB'),'update velocity')>0) then
        if(mpiworld%is_master) then
            call F1_star_lapE0%write
            call F2_star_lapdE%write
            call gvp2%write
        endif
        correlation_gradient=correlation_gradient +gvp2%sp_rp
    endif

    ! if(ppg%if_compute_engy) then
    !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    ! endif
        
    call mpiworld%barrier

end subroutine
