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

    character(:),allocatable :: update_wavelet
    type(t_weighter) :: wei_wl

    type(t_field) :: fld_u,fld_a
    type(t_correlate) :: a_star_u
    real,dimension(:,:),allocatable :: tmp
    real,dimension(3) :: grad_term_weights
    character(:),allocatable :: dnorm
    
    type :: t_S
        real,dimension(:),allocatable :: scale
    end type
    type(t_S),dimension(:),allocatable,save :: S
    real,dimension(:,:),allocatable :: tmp_dsyn, gmwindow

    if(is_first_in) allocate(S(shls%nshots_per_processor)) !then can NOT randomly sample shots..


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


        call hud('----  Solving Au=s  ----')
        call ppg%init_field(fld_u, name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire; call shot%write('Ru_',shot%dsyn)


        ! if(index(setup%get_str('JOB',o_default='gradient'),'estimate wavelet')>0) then
        !     call hud('--------------------------------')
        !     call hud('        Estimate wavelet        ')
        !     call hud('--------------------------------')

        !     call wei_wl%update

        !     call shot%update_wavelet(wei_wl%weight) !call gradient_matchfilter_data
        
        !     !write synthetic data
        !     call shot%write('updated_Ru_',shot%dsyn)

        !     cycle

        ! endif
        
        if(setup%get_str('JOB')=='forward modeling') cycle

        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)

        call hud('----  Computing obj func & dadj  ----')
            call wei%update
            call alloc(shot%dadj,shot%nt,shot%nrcv)

            if(.not.allocated(dnorm)) dnorm=setup%get_str('DATA_NORM','DNORM',o_default='L2')
            select case (dnorm)
                case ('L2')
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
                call kernel_L2sq(shot%dadj)
                call fld_a%ignite(o_wavelet=shot%dadj)
                call shot%write('dadj_',shot%dadj)

                case ('L2averaged')
                length=setup%get_int('MOVING_AVERAGE_LENGTH','MA_LEN',o_mandatory=1)
                ! scaler=setup%get_real('MOVING_AVERAGE_SCALER','MA_SCALER',o_default='1.')
                call moving_average(shot%dsyn,length)!,scaler)
                call shot%write('avg_dsyn_',shot%dsyn)
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
                call kernel_L2sq(shot%dadj)
                call moving_average(shot%dadj,length)!,scaler)
                call fld_a%ignite(o_wavelet=shot%dadj)
                call shot%write('dadj_',shot%dadj)


                case ('strain_L2averaged')
                !first goto strain
                call derivative_x(shot%dsyn,shot)
                call integrate_t(shot%dsyn,shot%dt)
                !then add gauge length
                length=setup%get_int('MOVING_AVERAGE_LENGTH','MA_LEN',o_mandatory=1)
                ! scaler=setup%get_real('MOVING_AVERAGE_SCALER','MA_SCALER',o_default='1.')
                call moving_average(shot%dsyn,length)!,scaler)
                call shot%write('DAS_dsyn_',shot%dsyn)
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
                call kernel_L2sq(shot%dadj)
                !goback to velocity adjoint src
                call moving_average(shot%dadj,length)!,scaler)
                call rev_integrate_t(shot%dadj,shot%dt)
                call derivative_x(-shot%dadj,shot)

                call fld_a%ignite(o_wavelet=shot%dadj)
                call shot%write('dadj_',shot%dadj)


                case('L2_scaled')
                    if(is_first_in) call alloc(S(i)%scale,shot%nrcv) !then can NOT randomly sample shots..
                    call alloc(tmp_dsyn,shot%nt,shot%nrcv)
                    do j=1,shot%nrcv
                        if(is_first_in) S(i)%scale(j) = either(0., maxval(abs(shot%dobs(:,j))) / maxval(abs(shot%dsyn(:,j))) , shot%rcv(j)%is_badtrace)
                        tmp_dsyn(:,j)=shot%dsyn(:,j)*S(i)%scale(j)
                    enddo
                    if(is_first_in) then
                        open(12,file=dir_out//'dobs_dsyn_max_ratio',access='direct',recl=4*shot%nrcv)
                        write(12,rec=shot%index) S(i)%scale
                        close(12)
                    endif

                    !check if S is changing..
        !            if(shot%index==1)   print*, 'on '//shot%sindex,i,S(i)%scale
        !            if(shot%index==112) print*, 'on '//shot%sindex,i,S(i)%scale

                    fobj%misfit = fobj%misfit &
                        + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-tmp_dsyn, shot%dt)

                    call alloc(shot%dadj,shot%nt,shot%nrcv)
                    call kernel_L2sq(shot%dadj)
                    call fld_a%ignite(o_wavelet=shot%dadj)
                    call shot%write('dadj_',shot%dadj)

                case default
                call error('No DNORM specified!')

            end select

        
        call hud('----  Solving adjoint eqn & xcorrelate  ----')

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

    if(setup%get_str('JOB')=='gradient') then
        call mpiworld%final
        stop
    endif

    is_first_in=.false.

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


    contains

    subroutine moving_average(data,length)
        real,dimension(:,:) :: data
        real,dimension(:,:),allocatable :: tmp
        call alloc(tmp,shot%nt,shot%nrcv)
        tmp=data

        L = either( (length-1)/2 , length/2 , mod(length,2)/=0 )

        do it=1,shot%nt
            do itr=1,L
                denom = 1./(itr+L)
                data(it,itr) = sum(tmp(it,1:itr+L)) * denom
            enddo
        enddo

        denom = 1./(2*L+1)
        do it=1,shot%nt
            do itr=L+1,shot%nrcv-L
                data(it,itr) = sum(tmp(it,itr-L:itr+L)) * denom
            enddo
        enddo

        do it=1,shot%nt
            do itr=shot%nrcv-L+1,shot%nrcv
                denom = 1./(shot%nrcv-itr+L+1)
                data(it,itr) = sum(tmp(it,itr-L:shot%nrcv)) * denom
            enddo
        enddo

        deallocate(tmp)

    end subroutine

    subroutine derivative_x(data,shot) !gradient by cdiff
        real,dimension(:,:) :: data !nt x nrcv
        type(t_shot) :: shot 

        real,dimension(:,:), allocatable :: dout
        call alloc(dout,shot%nt,shot%nrcv)

        do ir=2,shot%nrcv-1
            dout(:,ir) = (data(:,ir+1)-data(:,ir-1)) / (shot%rcv(ir+1)%x-shot%rcv(ir-1)%x)
        enddo

        !padding
        dout(:,1)=dout(:,2)
        dout(:,shot%nrcv)=dout(:,shot%nrcv-1)

        data=dout

    end subroutine

    subroutine integrate_t(data,dt) !integrate
        real,dimension(:,:) :: data !nt x nrcv
        real :: dt

        real,dimension(:,:), allocatable :: dout
        call alloc(dout,shot%nt,shot%nrcv)

        do ir=1,shot%nrcv
                dout(1,ir)=data(1,ir)*dt
            do it=2,shot%nt
                dout(it,ir) = dout(it-1,ir) + data(it,ir)*dt
            enddo
        enddo

        data=dout

    end subroutine

    subroutine rev_integrate_t(data,dt) !integrate
        real,dimension(:,:) :: data !nt x nrcv
        real :: dt

        real,dimension(:,:), allocatable :: dout
        call alloc(dout,shot%nt,shot%nrcv)

        do ir=1,shot%nrcv
                dout(shot%nt,ir)=data(shot%nt,ir)*dt
            do it=shot%nt-1,1,-1
                dout(it,ir) = dout(it+1,ir) + data(it,ir)*dt
            enddo
        enddo

        data=dout

    end subroutine

end subroutine
