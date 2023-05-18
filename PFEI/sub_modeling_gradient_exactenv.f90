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
use m_fracderi

    logical,save :: is_first_in=.true.

    character(:),allocatable :: update_wavelet
    type(t_weighter) :: wei_wl

    type(t_field) :: fld_u,fld_v, fld_F0, fld_dF
    type(t_correlate) :: F0_star_E, dF_star_E
    ! real,dimension(:,:),allocatable :: Wdres
    !real,dimension(:,:,:),allocatable :: Ddt2
    ! character(:),allocatable :: update_wavelet
    real,dimension(:,:),allocatable :: hilb

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


        call hud('----  Solving Au=s  ----')
        call ppg%init_field(fld_u, name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire; call shot%write('Ru_',shot%dsyn)


        if(index(setup%get_str('JOB',o_default='gradient'),'estimate wavelet')>0) then
            call hud('--------------------------------')
            call hud('        Estimate wavelet        ')
            call hud('--------------------------------')

            call wei_wl%update

            call shot%update_wavelet(wei_wl%weight) !call gradient_matchfilter_data
        
            !write synthetic data
            call shot%write('updated_RE0_',shot%dsyn)

            cycle

        endif

        call hud('----  Solving Av=H[s]  ----')
        call shot%read_wlhilb
        call ppg%init_field(fld_v, name='fld_v');    call fld_v%ignite
        call ppg%forward(fld_v)
        call fld_v%acquire; !call shot%write('Rv_',shot%dsyn)
        
        shot%dsyn_aux = shot%dsyn
        call fld_u%acquire
        shot%dsyn=sqrt(shot%dsyn**2+shot%dsyn_aux**2)
        call shot%write('RE_',shot%dsyn)

        if(setup%get_str('JOB')=='forward modeling') cycle

        call hud('----  Computing dadj=dobs-E[u]=d-RE₀-RδE  ----')
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

        call wei%update
        call alloc(shot%dadj,shot%nt,shot%nrcv)
        fobj%misfit = fobj%misfit &
            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
            ! + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn-shot%dsyn_aux, shot%dt)
        call kernel_L2sq(shot%dadj)
        call shot%write('dadj_',shot%dadj)

        
        call hud('----  Solving AᴴF₀=Rᴴdadj and AᴴδF=DF₀  ----')
        call hud('----  Forming (F₀+δF)★E  ----')
        call ppg%init_field(fld_dF,name='fld_dF',ois_adjoint=.true.)
        call ppg%init_field(fld_F0,name='fld_F0',ois_adjoint=.true.);  call fld_F0%ignite
        call ppg%init_correlate(dF_star_E,'dF_star_E','gradient')
        call ppg%init_correlate(F0_star_E,'F0_star_E','gradient')
        call ppg%adjoint(fld_dF,fld_F0, fld_v,fld_u, dF_star_E,F0_star_E)

        call hud('---------------------------------')

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    
    !allreduce PFEI misfit values
    call mpi_allreduce(mpi_in_place, [fobj%misfit], 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked PFEI_misfit '//num2str(fobj%misfit))

    fobj%dnorms=fobj%misfit

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
    
    !scale by shotlist
    call shls%scale(1,o_from_sampled=[fobj%misfit])
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)


    !write correlate
    if(mpiworld%is_master) then
        call dF_star_E%write
        call F0_star_E%write
    endif


call mpi_allreduce(mpi_in_place, dF_star_E%rp_ddsp, m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
call mpi_allreduce(mpi_in_place, F0_star_E%rp_ddsp, m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
if(mpiworld%is_master) then
    call dF_star_E%write(o_suffix='_stacked')
    call F0_star_E%write(o_suffix='_stacked')


    den1=sqrt(sum(F0_star_E%rp_ddsp**2,.not.m%is_freeze_zone))
    den2=sqrt(sum(dF_star_E%rp_ddsp**2,.not.m%is_freeze_zone))
    costh=sum(F0_star_E%rp_ddsp/den1*dF_star_E%rp_ddsp/den2,.not.m%is_freeze_zone)

    call hud('Angle between F₀★E & δF★E (w/ mask): '//num2str(acosd(costh))//'°')
    
endif


    !allreduce energy, gradient
    ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale by shotlist
    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)

    if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)


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