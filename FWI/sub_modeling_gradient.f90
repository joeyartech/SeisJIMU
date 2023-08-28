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
!use m_hilbert
use m_hilbert_nofft

    logical,save :: is_first_in=.true.

    character(:),allocatable :: update_wavelet
    type(t_weighter) :: wei_wl

    type(t_field) :: fld_u,fld_v, fld_p, fld_q
    type(t_correlate) :: U_star_D,D_star_D,U_star_U, D_star_U
    type(t_correlate) :: a_star_u
    real,dimension(:,:),allocatable :: tmp
    real,dimension(3) :: grad_term_weights
    character(:),allocatable :: dnorm

    real,dimension(:,:,:),allocatable :: term1,term2
    

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


        if(index(setup%get_str('JOB',o_default='gradient'),'estimate wavelet')>0) then
            call hud('--------------------------------')
            call hud('        Estimate wavelet        ')
            call hud('--------------------------------')

            call wei_wl%update

            call shot%update_wavelet(wei_wl%weight) !call gradient_matchfilter_data
        
            !write synthetic data
            call shot%write('updated_Ru_',shot%dsyn)

            cycle

        endif

        call hud('----  Solving Av=H[s]  ----')
        call shot%read_wlhilb
        call ppg%init_field(fld_v, name='fld_v');    call fld_v%ignite
        call ppg%forward(fld_v)
        call fld_v%acquire; call shot%write('Rv_',shot%dsyn); shot%dsyn_aux=shot%dsyn
        call fld_u%acquire
        
        if(setup%get_str('JOB')=='forward modeling') cycle


        call ppg%init_field(fld_p,name='fld_p',ois_adjoint=.true.)
        call ppg%init_field(fld_q,name='fld_q',ois_adjoint=.true.)

        call hud('----  Computing obj func & dadj  ----')
            call wei%update
            call alloc(shot%dadj,shot%nt,shot%nrcv)

            if(.not.allocated(dnorm)) dnorm=setup%get_str('DATA_NORM','DNORM',o_default='L2')
            select case (dnorm)
                case ('L2')
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
                call kernel_L2sq(shot%dadj)
                call fld_p%ignite(o_wavelet=shot%dadj)
                call shot%write('dadj_',shot%dadj)
                tmp=shot%dadj
                !call hilbert_transform(shot%dadj,tmp,shot%nt,shot%nrcv)
                call hilbert_nofft('generic',shot%dadj,tmp,shot%nt,shot%nrcv)
                call fld_q%ignite(o_wavelet=tmp)

                ! case default
                ! call error('No DNORM specified!')

            end select

        
        call hud('----  Solving adjoint eqn & xcorrelate  ----')

        call ppg%init_correlate(a_star_u,'a_star_u')
        call ppg%adjoint(fld_p,fld_u,a_star_u)

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
