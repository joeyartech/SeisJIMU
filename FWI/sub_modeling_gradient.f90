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
! use m_hilbert_nofft
use m_butterworth

    logical,save :: is_first_in=.true.

    character(:),allocatable :: update_wavelet, dnorm
    type(t_weighter) :: wei_wl

    type(t_field) :: fld_reU,fld_imU, fld_reA, fld_imA
    type(t_correlate) :: A_star_U    
    
    !misfit
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


        call hud('----  Solving AU=S  ----')
        call ppg%init_field(fld_reU, name='fld_reU')
        call ppg%init_field(fld_imU, name='fld_imU')
        
        call fld_reU%ignite
        call fld_imU%ignite(o_wavelet=reshape(hilbert(shot%wavelet),[shot%nt,1]))

        call ppg%forward(fld_reU,fld_imU)
        call fld_imU%acquire; call shot%write('imRU_',shot%dsyn)
        call fld_reU%acquire; call shot%write('Ru_',shot%dsyn)       
                
        if(setup%get_str('JOB')=='forward modeling') cycle


        call ppg%init_field(fld_reA,name='fld_reA',ois_adjoint=.true.)
        call ppg%init_field(fld_imA,name='fld_imA',ois_adjoint=.true.)

        call hud('----  Computing obj func & dadj  ----')
            call wei%update
            call alloc(shot%dadj,shot%nt,shot%nrcv)


call hud('butterworth filtering')
fpasshi=setup%get_real('FILT_FPEAK',o_default=num2str(shot%fpeak))
fstophi=setup%get_real('FILT_FMAX', o_default=num2str(shot%fmax ))
call butterworth(shot%dobs,shot%nt,shot%nrcv,shot%dt, oif_locut=.false., o_fpasshi=fpasshi,o_fstophi=fstophi)
call butterworth(shot%dsyn,shot%nt,shot%nrcv,shot%dt, oif_locut=.false., o_fpasshi=fpasshi,o_fstophi=fstophi)
call shot%write('filt_dobs_',shot%dobs)
call shot%write('filt_dsyn_',shot%dsyn)

            if(.not.allocated(dnorm)) dnorm=setup%get_str('DATA_NORM','DNORM',o_default='L2')
            select case (dnorm)
                case ('L2')
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
                call kernel_L2sq(shot%dadj)

! call shot%write('raw_dadj_',shot%dadj)
! call hud('butterworth filtering')
! call butterworth(shot%dadj,shot%nt,shot%nrcv,shot%dt, &
!     oif_locut=.false., &
!     o_fpasshi=shot%fpeak,o_fstophi=shot%fmax)


                call fld_reA%ignite(o_wavelet=shot%dadj)
                call shot%write('dadj_',shot%dadj)
                call shot%write('Hdadj_',hilbert(shot%dadj))

                call fld_imA%ignite(o_wavelet=hilbert(shot%dadj))
                                
                case default
                call error('No DNORM specified!')

            end select

        
        call hud('----  Solving adjoint eqn & xcorrelate  ----')
        call ppg%init_correlate(A_star_U,'A_star_U')
        call ppg%adjoint(fld_reA,fld_imA, fld_reU,fld_imU, A_star_U)
        ! call ppg%adjoint(fld_reA, fld_reU, A_star_U)


        call hud('----  Assemble  ----')
        call ppg%assemble(A_star_U)

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
        call A_star_U%write
    endif

    !allreduce energy, gradient
    ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale by shotlist
    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)

    if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)


    call mpiworld%barrier


    ! contains
    ! function deri(a) result(d)
    !     real,dimension(:,:),allocatable :: a,d
    !     intent(in) :: a
    !     real inv_2dt

    !     inv_2dt = 1./2/shot%dt
        
    !     d=a
    !     do ir=1,shot%nrcv
    !         do it=2,shot%nt-1
    !             d(it,ir) = (a(it+1,ir)-a(it-1,ir))*inv_2dt
    !         enddo
    !             d(1,ir) = d(2,ir)
    !             d(shot%nt,ir) = d(shot%nt-1,ir)
    !     enddo

    ! end function


end subroutine
