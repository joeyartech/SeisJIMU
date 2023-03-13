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
    type(t_correlation) :: F2_star_E0, F1_star_E0, F2_star_dE
    real,dimension(:,:,:,:),allocatable :: gtilD, gvp2
    ! real,dimension(:,:),allocatable :: Wdres
    !real,dimension(:,:,:),allocatable :: Ddt2
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff

    character(:),allocatable :: dnorm
    
    call alloc(gtilD,m%nz,m%nx,m%ny,1)
    call alloc(gvp2 ,m%nz,m%nx,m%ny,1)

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

        if(index(setup%get_str('JOB'),'estimate wavelet')>0) then
            call hud('--------------------------------')
            call hud('        Estimate wavelet        ')
            call hud('--------------------------------')

            call wei%update

            call shot%update_wavelet(wei%weight) !call gradient_matchfilter_data
        
            !write synthetic data
            call shot%write('updated_RE0_',shot%dsyn)

        endif

        if(index(setup%get_str('JOB'),'build tilD')>0) then
            call hud('--------------------------------')
            call hud('        Build tilD model        ')
            call hud('--------------------------------')

            call F2_star_E0%init('F2_star_E0',shape123_from='model') !F₂★E₀ for gtDt
            call alloc(F2_star_E0%drp_dt_dsp_dt,m%nz,m%nx,m%ny,1)
            call alloc(F2_star_E0%nab_rp_nab_sp,m%nz,m%nx,m%ny,1)

            call wei%update!('_4IMAGING')
            fobj%misfit = fobj%misfit &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-(shot%dsyn+shot%dsyn_aux), shot%dt)
            
            call alloc(shot%dadj,shot%nt,shot%nrcv)
            call kernel_L2sq(shot%dadj)
            call shot%write('dadj_',shot%dadj)

            !adjoint modeling
            !AᴴF₂ = -RᴴN(E+δE)
            call ppg%init_field(fld_F2,name='fld_F2',ois_adjoint=.true.); call fld_F2%ignite
            call ppg%adjoint_tilD(fld_F2,fld_E0,o_F2_star_E0=F2_star_E0)
            call hud('---------------------------------')

            !call cb%project_back
            gtilD(cb%ioz:cb%ioz+cb%mz-1,&
                  cb%iox:cb%iox+cb%mx-1,&
                  cb%ioy:cb%ioy+cb%my-1,1) = &
            gtilD(cb%ioz:cb%ioz+cb%mz-1,&
                  cb%iox:cb%iox+cb%mx-1,&
                  cb%ioy:cb%ioy+cb%my-1,1) + &
                F2_star_E0%drp_dt_dsp_dt(:,:,:,1) - m%vp**2*F2_star_E0%nab_rp_nab_sp(:,:,:,1) !inverse scattering imaging condition
            ! gtDs%sp_rp(cb%ioz:cb%ioz+cb%mz-1,&
            !            cb%iox:cb%iox+cb%mx-1,&
            !            cb%ioy:cb%ioy+cb%my-1,shot%index) = F2_star_E0%sp_rp(:,:,:,1)

        endif
        
        if(index(setup%get_str('JOB'),'update velocity')>0) then
            call hud('---------------------------------')
            call hud('     Update velocity model       ')
            call hud('---------------------------------')

            call F1_star_E0%init('F1_star_E0',shape123_from='model') !F₁★∇²E₀ for gvp2 (one RE)
            call F2_star_dE%init('F2_star_dE',shape123_from='model') !F₂★∇²δE for gvp2 (the other RE)
            call F2_star_E0%init('F2_star_E0',shape123_from='model') !F₂★∇²E₀ for gvp2 (3rd RE?)
            call alloc(F1_star_E0%nab_rp_nab_sp,m%nz,m%nx,m%ny,1)
            call alloc(F2_star_dE%nab_rp_nab_sp,m%nz,m%nx,m%ny,1)
            call alloc(F2_star_E0%nab_rp_nab_sp,m%nz,m%nx,m%ny,1)

            call wei%update
            call alloc(shot%dadj,shot%nt,shot%nrcv)

            if(.not.allocated(dnorm)) dnorm=setup%get_str('DATA_NORM','DNORM',o_default='L2sq_sq')
            select case (dnorm)
                case ('L2sq')
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-(shot%dsyn+shot%dsyn_aux), shot%dt)
                call kernel_L2sq(shot%dadj)

                case ('L2sq_abs')
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, abs(shot%dobs)-abs(shot%dsyn+shot%dsyn_aux), shot%dt)
                call kernel_L2sq(shot%dadj)
                shot%dadj = shot%dadj*sgns(shot%dsyn+shot%dsyn_aux)!,shot%dobs)

		case default
		call error('No DNORM specified!')
                !case default
                !fobj%misfit = fobj%misfit &
                !    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, (shot%dobs)**2-(shot%dsyn+shot%dsyn_aux)**2, shot%dt)
                !call kernel_L2sq(shot%dadj)
                !shot%dadj = shot%dadj*2.*(shot%dsyn+shot%dsyn_aux)

            end select

            call shot%write('dadj_',shot%dadj)
            
            !adjoint modeling
            !AᴴF₂ =      +Rᴴ(d-(E₀+δE))
            !AᴴF₁ = -DF₂ +Rᴴ(d-(E₀+δE))
            call ppg%init_field(fld_F1,name='fld_F1',ois_adjoint=.true.); call fld_F1%ignite
            call ppg%init_field(fld_F2,name='fld_F2',ois_adjoint=.true.); call fld_F2%ignite
            call ppg%adjoint_vp2(fld_F1,fld_F2,fld_dE,fld_E0,o_F1_star_E0=F1_star_E0,o_F2_star_dE=F2_star_dE,o_F2_star_E0=F2_star_E0)
            call hud('---------------------------------')
    
            !project back    
            gvp2(cb%ioz:cb%ioz+cb%mz-1,&
                 cb%iox:cb%iox+cb%mx-1,&
                 cb%ioy:cb%ioy+cb%my-1,1) = &
            gvp2(cb%ioz:cb%ioz+cb%mz-1,&
                 cb%iox:cb%iox+cb%mx-1,&
                 cb%ioy:cb%ioy+cb%my-1,1) + &
                F1_star_E0%nab_rp_nab_sp(:,:,:,1) +F2_star_dE%nab_rp_nab_sp(:,:,:,1) -m%tilD*F2_star_E0%nab_rp_nab_sp(:,:,:,1)
            
        endif

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    if(setup%get_str('JOB')=='forward modeling') then
        call mpiworld%final
        stop
    endif
    if(setup%get_str('JOB')=='estimate wavelet') then
        call mpiworld%final
        stop
    endif

    !collect PFEI misfit values
    !no need to collect WPI misfit values
    call mpi_allreduce(mpi_in_place, [fobj%misfit], 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked PFEI_misfit '//num2str(fobj%misfit))

    fobj%dnorms=fobj%misfit

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
    
    !collect global correlations
    call mpi_allreduce(mpi_in_place, gtilD, m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place,  gvp2, m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale
    call shls%scale(1,o_from_sampled=[fobj%misfit])
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)
    call shls%scale(m%n*3,o_from_sampled=m%correlate)


    !call alloc(m%gradient,m%nz,m%nx,m%ny,1)

    !call alloc(correlation_gradient,m%nz,m%nx,m%ny,1)

    if(index(setup%get_str('JOB'),'build tilD')>0) then
        if(mpiworld%is_master) then
            call sysio_write('gtilD_dt2' ,F2_star_E0%drp_dt_dsp_dt,m%n)
            call sysio_write('gtilD_nab2',F2_star_E0%nab_rp_nab_sp,m%n)
            call sysio_write('gtilD',gtilD,m%n)
        endif
        correlation_gradient=gtilD
    endif

    if(index(setup%get_str('JOB'),'update velocity')>0) then
        if(mpiworld%is_master) then
            call sysio_write('gvp2_F1_star_E0',F1_star_E0%nab_rp_nab_sp,m%n)
            call sysio_write('gvp2_F2_star_dE',F2_star_dE%nab_rp_nab_sp,m%n)
            call sysio_write('gvp2_F2_star_E0',F2_star_E0%nab_rp_nab_sp,m%n)
            call sysio_write('gvp2',gvp2,m%n)
        endif
        correlation_gradient=gvp2
    endif

    ! if(ppg%if_compute_engy) then
    !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    ! endif
        
    call mpiworld%barrier


    contains
    pure function sgns(a) result(s)
!    use m_math, only: r_eps
        real,dimension(:,:),intent(in) :: a
        real,dimension(:,:),allocatable :: s
        s=a
        !where (a>r_eps)
        !    s=1.
        !elsewhere (a<-r_eps)
        !    s=-1.
        !elsewhere
        !    s=0.
        !endwhere
        where (a>0.)
            s=1.
        elsewhere
            s=-1.
        endwhere
    end function

end subroutine
