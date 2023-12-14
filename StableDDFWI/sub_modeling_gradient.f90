! subroutine modeling_gradient_ip
! use mpi
! use m_System
! use m_Modeling
! use m_separator
! use m_weighter
! use m_Lpnorm
! use m_fobjective
! use m_matchfilter
! use m_smoother_laplacian_sparse
! use m_resampler
!
!     logical,save :: is_first_in=.true.
!
!     type(t_field) :: fld_u, fld_a
!     type(t_correlate) :: a_star_u
!
!     character(:),allocatable :: s_job
!     ! real,dimension(:,:),allocatable :: Wdres
!     !real,dimension(:,:,:),allocatable :: Ddt2
!     ! character(:),allocatable :: update_wavelet
!     ! real,dimension(:,:),allocatable :: Rmu, Rdiff
!
!     !RWI misfit
!     fobj%misfit=0.
!     fobj%reflection=0.
!     fobj%diving=0.
!
!     call alloc(correlate_gradient,m%nz,m%nx,m%ny,ppg%ngrad)
!
!     call hud('===== START LOOP OVER SHOTS =====')
!
!     do i=1,shls%nshots_per_processor
!
!         call shot%init(shls%yield(i))
!         call shot%read_from_data
!         call shot%set_var_time
!         call shot%set_var_space(index(ppg%info,'FDSG')>0)
!
!         call hud('Modeling Shot# '//shot%sindex)
!
!         call cb%init(ppg%nbndlayer)
!         call cb%project
!         call ppg%check_discretization
!         call ppg%init
!         call ppg%init_abslayer
!
!
!         s_job=setup%get_str('JOB',o_mandatory=1)
!
!         !forward modeling
!         !A(m )u = s
!         call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
!         call ppg%forward(fld_u)
!         call fld_u%acquire;  call shot%write('Ru_',shot%dsyn)
!
!         !update wavelet
!         if(setup%get_str('UPDATE_WAVELET')/='no') then
! !            call wei%update
!             call shot%update_wavelet!(wei%weight)
!             call matchfilter_apply_to_data(shot%dsyn)
!
!             !write synthetic data
!             call shot%write('updated_Ru_',shot%dsyn)
!         endif
!
!         call hud('---------------------------------')
!         if(index(s_job,'modeling')>0) cycle
!
!
!         call sepa%update
!         call wei%update!('_4IMAGING')
!         call alloc(shot%dadj,shot%nt,shot%nrcv)
!
!         call hud('------------------------')
!         call hud('     Build Ip model     ')
!         call hud('------------------------')
!
!         fobj%misfit = fobj%misfit &
!             + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*sepa%nearoffset*sepa%reflection,shot%dobs-shot%dsyn, shot%dt)
! !            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*sepa%nearoffset, shot%dobs-shot%dsyn, shot%dt)
!
!         call kernel_L2sq(shot%dadj)
!         call shot%write('dadj_',shot%dadj)
!
!         !adjoint modeling
!         !A(m)ᴴa = Rʳ(d-u)
!         call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.); call fld_a%ignite
!         call ppg%init_correlate(a_star_u,'a_star_u','gradient') !a★u
!         call ppg%adjoint(fld_a,fld_u,o_a_star_u=a_star_u)
!         call hud('---------------------------------')
!
!     enddo
!
!     call hud('        END LOOP OVER SHOTS        ')
!
!     if(index(s_job,'modeling')>0) then
!         call mpiworld%final
!         stop
!     endif
!
!
!     !allreduce RWI misfit values
!     call mpi_allreduce(mpi_in_place, [fobj%reflection,fobj%diving,fobj%misfit], 3, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
!     call hud('Stacked RWI reflection/diving/total misfit = '//num2str(fobj%reflection)//'/'//num2str(fobj%diving)//'/'//num2str(fobj%misfit))
!
!     fobj%dnorms=fobj%misfit
!
!     call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
!
!     !scale by shotlist
!     call shls%scale(1,o_from_sampled=[fobj%misfit])
!     call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)
!
!     !write correlate
!     if(mpiworld%is_master) then
!         call a_star_u%write
!     endif
!
!
!     !allreduce energy, gradient
!     ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
!     call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
!
!     !scale by shotlist
!     call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)
!
!     if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)
!
!     call mpiworld%barrier
!
! end subroutine


subroutine modeling_gradient_m0
use mpi
use m_System
use m_Modeling
! use m_separator
use m_weighter
use m_Lpnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_ub, fld_um, fld_ab, fld_am
    type(t_correlate) :: am_star_um, ab_star_ub

    character(:),allocatable :: s_job
    ! real,dimension(:,:),allocatable :: Wdres
    !real,dimension(:,:,:),allocatable :: Ddt2
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff
    
    !misfit
    fobj%misfit=0.
    fobj%fourD=0.
    fobj%baseline=0.

    call alloc(correlate_gradient,m%nz,m%nx,m%ny,ppg%ngrad)
    
    call hud('===== START LOOP OVER SHOTS =====')
        
    do i=1,shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_4Ddata
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project(ois_monitor=.true.)
        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer


        s_job=setup%get_str('JOB',o_default='forward modeling')

        !forward modeling
        !A(m )u = s
        call ppg%init_field(fld_um,name='fld_um');    call fld_um%ignite
        call ppg%forward(fld_um)
        call fld_um%acquire(o_seismo=shot%dsynm);  call shot%write('Rum_',shot%dsynm)

        !A(m₀)u₀= s
        call cb%project
        call ppg%init
        call ppg%init_field(fld_ub,name='fld_ub');    call fld_ub%ignite
        call ppg%forward(fld_ub)
        call fld_ub%acquire(o_seismo=shot%dsynb);  call shot%write('Rub_',shot%dsynb)

!         !update wavelet
!         if(setup%get_str('UPDATE_WAVELET')/='no') then
!             call wei%update
!             call shot%update_wavelet!(wei%weight)
!             call matchfilter_apply_to_data(shot%dsyn)
!             call matchfilter_apply_to_data(shot%dsyn_aux)
!
!             !write synthetic data
!             call shot%write('updated_Ru_',shot%dsyn)
!             call shot%write('updated_Ru0_',shot%dsyn_aux)
!         endif
!
!         call hud('---------------------------------')
!         if(index(s_job,'modeling')>0) cycle


!         call sepa%update
        call wei%update!('_4IMAGING')
        call alloc(shot%dadj,shot%nt,shot%nrcv)

        call hud('-------------------------')
        call hud('     Update m0           ')
        call hud('-------------------------')

        !4D data
        fobj%fourD = fobj%fourD &
            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, (shot%dobsm-shot%dobsb)-(shot%dsynm-shot%dsynb), shot%dt)
!             + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*(1.-sepa%nearoffset)*sepa%reflection, shot%dobs-shot%dsyn, shot%dt)

        call kernel_L2sq(shot%dadj)
        call shot%write('dadj4D_',shot%dadj)
                    
        !adjoint modeling
        !A(m )ᴴa  = Rʳ(d-u)
        call cb%project(ois_monitor=.true.)
        call ppg%init
        call ppg%init_field(fld_am,name='fld_am' ,ois_adjoint=.true.); call fld_am%ignite
        call ppg%init_correlate(am_star_um,'am^4D_star_um','gradient') !a★u
        call ppg%adjoint(fld_am,fld_um,o_a_star_u=am_star_um)

        !baseline data
        fobj%baseline = fobj%baseline &
            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobsb-shot%dsynb, shot%dt)
!             + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*(1.-sepa%nearoffset)*sepa%diving, shot%dobs-shot%dsyn, shot%dt)
        
!         shot%dadj=-shot%dadj; call kernel_L2sq(shot%dadj,oif_stack=.true.) !diving minus reflection residuals
        call shot%write('dadjb_',shot%dadj)

        !adjoint modeling
        !A(m₀)ᴴa₀ = Rᵈ(d-u)-Rʳ(d-u)
        call cb%project
        call ppg%init
        call ppg%init_field(fld_ab,name='fld_ab',ois_adjoint=.true.); call fld_ab%ignite
        call ppg%init_correlate(ab_star_ub,'ab^bas_star_ub','gradient') !a₀★u₀
        call ppg%adjoint(fld_ab,fld_ub,o_a_star_u=ab_star_ub)

        !produce final misfit
        fobj%misfit = fobj%fourD+fobj%baseline

        call hud('---------------------------------')
            
    enddo

    call hud('        END LOOP OVER SHOTS        ')


    !allreduce RWI misfit values
    call mpi_allreduce(mpi_in_place, [fobj%fourD,fobj%baseline,fobj%misfit], 3, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked RWI fourD/baseline/total misfit = '//num2str(fobj%fourD)//'/'//num2str(fobj%baseline)//'/'//num2str(fobj%misfit))

    fobj%dnorms=fobj%misfit

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
    
    !scale by shotlist
    call shls%scale(1,o_from_sampled=[fobj%misfit])
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)

    !write correlate
    if(mpiworld%is_master) then
        call am_star_um%write
        call ab_star_ub%write
    endif
    

    !rebuild by summation
    !gkpa
    correlate_gradient(:,:,:,1) = -ppg%inv_kpa(1:m%nz,1:m%nx,1:m%ny)*(am_star_um%rp_div_sv+ab_star_ub%rp_div_sv)
    !grho0
    correlate_gradient(:,:,:,2) = (am_star_um%rv_grad_sp+ab_star_ub%rv_grad_sp) / m%rho
    


    !allreduce energy, gradient
    ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale by shotlist
    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)

    if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)

    call mpiworld%barrier
    
    if(index(s_job,'modeling')>0) then
        call mpiworld%final
        stop
    endif

end subroutine
