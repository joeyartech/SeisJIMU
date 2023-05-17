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
use m_hilbert

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_u,fld_E0, fld_dE, fld_F1, fld_F2
    type(t_correlate) :: F2_star_E0, F1_star_E0, F2_star_dE
    ! real,dimension(:,:),allocatable :: Wdres
    !real,dimension(:,:,:),allocatable :: Ddt2
    ! character(:),allocatable :: update_wavelet
    real,dimension(:,:),allocatable :: Eu, dfault

    character(:),allocatable :: dnorm
    

    !PFEI misfit
    fobj%misfit=0.

    call alloc(correlate_gradient,m%nz,m%nx,m%ny,ppg%ngrad)
    
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


        !forward modeling
        
        call hud('----  Solving Au=s  ----')
        call ppg%init_field(fld_u, name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire
        call shot%write('Ru_',shot%dsyn)
        ! Eu=shot%dsyn
    call alloc(Eu,shot%nt,shot%nrcv)
    call hilbert_transform(shot%dsyn,Eu,shot%nt,shot%nrcv)
    call shot%write('hilb',Eu)
    call alloc(Eu,shot%nt,shot%nrcv)
        call hilbert_envelope(shot%dsyn,Eu,shot%nt,shot%nrcv)
        call shot%write('E[u]_',Eu)

stop


        call hud('----  Computing dadj=dobs-RE[u]  ----')
        call wei%update
        call alloc(shot%dadj,shot%nt,shot%nrcv)
        fobj%misfit = fobj%misfit &
            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-Eu, shot%dt)
        call kernel_L2sq(shot%dadj)
        call shot%write('dadj_',shot%dadj)


        call hud('----  Solving AE₀=s  ----')
        call shot%read_wlenv
        call ppg%init_field(fld_E0,name='fld_E0');   call fld_E0%ignite
        call ppg%forward(fld_E0)
        call fld_E0%acquire
        call shot%write('RE0_',shot%dsyn)
        dfault=Eu-shot%dsyn
        call shot%write('E[u]-RE0_',dfault)

        
        call hud('----  Solving AδE=DE₀ and AᴴF₂=Rᴴ(d-RE₀-RδE)  ----')
        call hud('----  Forming F₂★δE  ----')
        call ppg%init_field(fld_dE,name='fld_dE');                     fld_dE%wavelet=transpose(dfault)
        call ppg%init_field(fld_F2,name='fld_F2',ois_adjoint=.true.);  call fld_F2%ignite
        call ppg%init_correlate(F2_star_dE,'F2_star_dE','gikpa_gbuo')
        call ppg%adjoint_F2_star_dE(fld_F2,fld_dE,fld_E0,F2_star_dE)
        call fld_dE%acquire; call shot%write('RdE_',shot%dsyn)


        call hud('----  Solving AE₀=s again ----')
        call ppg%init_field(fld_E0,name='fld_E0');   call fld_E0%ignite
        call ppg%forward(fld_E0)


        call hud('----  Solving AᴴF₁=DF₂+Rᴴ(d-RE₀-RδE)  ----')
        call hud('----  Forming F₁★E₀  ----')
        call ppg%init_field(fld_dE,name='fld_dE');                     fld_dE%wavelet=transpose(dfault)
        call ppg%init_field(fld_F2,name='fld_F2',ois_adjoint=.true.);  call fld_F2%ignite
        call ppg%init_field(fld_F1,name='fld_F1',ois_adjoint=.true.);  call fld_F1%ignite
        call ppg%init_correlate(F1_star_E0,'F1_star_E0','gikpa_gbuo')
        call ppg%adjoint_F1_star_E0(fld_F1,fld_F2,fld_dE,fld_E0,F1_star_E0)

        stop


        ! if(setup%get_str('JOB')=='forward modeling') cycle

        ! if(index(setup%get_str('JOB'),'estimate wavelet')>0) then
        !     call hud('--------------------------------')
        !     call hud('        Estimate wavelet        ')
        !     call hud('--------------------------------')

        !     call wei%update

        !     call shot%update_wavelet(wei%weight) !call gradient_matchfilter_data
        
        !     !write synthetic data
        !     call shot%write('updated_RE0_',shot%dsyn)

        ! endif

        ! if(index(setup%get_str('JOB'),'build tilD')>0) then
        !     call hud('--------------------------------')
        !     call hud('        Build tilD model        ')
        !     call hud('--------------------------------')

        !     call ppg%init_correlate(F2_star_E0,'F2_star_E0','gtilD') !F₂★E₀ for gtDt
            
        !     call wei%update!('_4IMAGING')
        !     fobj%misfit = fobj%misfit &
        !         + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-(shot%dsyn+shot%dsyn_aux), shot%dt)
            
        !     call alloc(shot%dadj,shot%nt,shot%nrcv)
        !     call kernel_L2sq(shot%dadj)
        !     call shot%write('dadj_',shot%dadj)

        !     !adjoint modeling
        !     !AᴴF₂ = -RᴴN(E+δE)
        !     call ppg%init_field(fld_F2,name='fld_F2',ois_adjoint=.true.); call fld_F2%ignite
        !     call ppg%adjoint_tilD(fld_F2,fld_E0,o_F2_star_E0=F2_star_E0)
        !     call hud('---------------------------------')

        ! endif

        ! if(index(setup%get_str('JOB'),'update velocity')>0) then
        !     call hud('---------------------------------')
        !     call hud('     Update velocity model       ')
        !     call hud('---------------------------------')

        !     call ppg%init_correlate(F1_star_E0,'F1_star_E0','gikpa_gbuo')    !F₁★∇²E₀
        !     call ppg%init_correlate(F2_star_dE,'F2_star_dE','gikpa_gbuo')    !F₂★∇²δE
        !     call ppg%init_correlate(F2_star_E0,'F2_star_E0','gikpa_gbuo_MI') !F₂★∇²E₀
            
        !     call wei%update
        !     call alloc(shot%dadj,shot%nt,shot%nrcv)

        !     if(.not.allocated(dnorm)) dnorm=setup%get_str('DATA_NORM','DNORM',o_default='L2sq')

        !     select case (dnorm)
        !         case ('L2sq')
        !         fobj%misfit = fobj%misfit &
        !             + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-(shot%dsyn+shot%dsyn_aux), shot%dt)
        !         call kernel_L2sq(shot%dadj)

        !         case ('L2sq_abs')
        !         fobj%misfit = fobj%misfit &
        !             + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, abs(shot%dobs)-abs(shot%dsyn+shot%dsyn_aux), shot%dt)
        !         call kernel_L2sq(shot%dadj)
        !         shot%dadj = shot%dadj*sgns(shot%dsyn+shot%dsyn_aux)

        !         case ('L2sq_modabs')
        !         fobj%misfit = fobj%misfit &
        !             + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, abs(shot%dobs)-abs(shot%dsyn+shot%dsyn_aux), shot%dt)
        !         call kernel_L2sq(shot%dadj)
        !         shot%dadj = shot%dadj*sgns2(shot%dsyn+shot%dsyn_aux,  shot%dobs-(shot%dsyn+shot%dsyn_aux)  )

        !         case default
        !         call error('No DNORM specified!')
        !         !fobj%misfit = fobj%misfit &
        !         !    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, (shot%dobs)**2-(shot%dsyn+shot%dsyn_aux)**2, shot%dt)
        !         !call kernel_L2sq(shot%dadj)
        !         !shot%dadj = shot%dadj*2.*(shot%dsyn+shot%dsyn_aux)

        !     end select

        !     call shot%write('dadj_',shot%dadj)
            
        !     !adjoint modeling
        !     !AᴴF₂ =      +Rᴴ(d-(E₀+δE))
        !     !AᴴF₁ = -DF₂ +Rᴴ(d-(E₀+δE))
        !     call ppg%init_field(fld_F1,name='fld_F1',ois_adjoint=.true.); call fld_F1%ignite
        !     call ppg%init_field(fld_F2,name='fld_F2',ois_adjoint=.true.); call fld_F2%ignite
        !     call ppg%adjoint(fld_F1,fld_F2,fld_dE,fld_E0,o_F1_star_E0=F1_star_E0,o_F2_star_dE=F2_star_dE,o_F2_star_E0=F2_star_E0)
        !     call hud('---------------------------------')
            
        ! endif

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    ! if(setup%get_str('JOB')=='forward modeling') then
    !     call mpiworld%final
    !     stop
    ! endif
    ! if(setup%get_str('JOB')=='estimate wavelet') then
    !     call mpiworld%final
    !     stop
    ! endif


    ! !allreduce PFEI misfit values
    ! call mpi_allreduce(mpi_in_place, [fobj%misfit], 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    ! call hud('Stacked PFEI_misfit '//num2str(fobj%misfit))

    ! fobj%dnorms=fobj%misfit

    ! call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
    
    ! !scale by shotlist
    ! call shls%scale(1,o_from_sampled=[fobj%misfit])
    ! call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)


    ! !write correlate
    ! if(mpiworld%is_master) then
    !     call F2_star_E0%write
    !     call F1_star_E0%write
    !     call F2_star_dE%write
    ! endif

    ! !allreduce energy, gradient
    ! ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    ! call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    ! !scale by shotlist
    ! call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)

    ! if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)

    call mpiworld%barrier


    contains
    pure function sgns(a) result(s)
    use m_math, only: r_eps
        real,dimension(:,:),intent(in) :: a
        real,dimension(:,:),allocatable :: s
        s=a
        where (a>r_eps)
            s=1.
        elsewhere (a<-r_eps)
            s=-1.
        elsewhere
            s=0.
        endwhere
        !where (a>0.)
        !    s=1.
        !elsewhere
        !    s=-1.
        !endwhere
    end function

    pure function sgns2(a,b) result(s)
        real,dimension(:,:),intent(in) :: a,b
        real,dimension(:,:),allocatable :: s
        s=a
        where (a>r_eps)
            s=1.
        elsewhere (a<-r_eps)
            s=-1.
        elsewhere
            where (b>r_eps)
                s=1.
            elsewhere (b<-r_eps)
                s=-1.
            elsewhere
                s=0.
            endwhere
        endwhere
    end function

end subroutine