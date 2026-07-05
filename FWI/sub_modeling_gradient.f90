subroutine modeling_gradient
use mpi
use m_System
use m_Modeling
use m_weighter
use m_Hilbert
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
    real,dimension(:,:),allocatable :: tmp, Eobs
    real,dimension(3) :: grad_term_weights
    
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


        if(setup%get_str('JOB')=='forward') cycle

        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)

        s_update_wavelet=setup%get_str('UPDATE_WAVELET')
        if(s_update_wavelet/='') call wei_wl%update(o_suffix='_4WAVELET')

        call hud('----  Computing obj func & dadj  ----')

            call wei%update
            call alloc(shot%dadj,shot%nt,shot%nrcv)

            select case (setup%get_str('DATA_NORM','DNORM',o_default='L2'))

            case ('L2'); call hud('0.5|| W(f*u - d)||² => adjsrc = f★ W W(f*u - d)')

if(s_update_wavelet/='') then
call hud('update wavelet')
call shot%update_wavelet(wei_wl%weight) !call matchfilter_apply_to_data(shot%dsyn)
call shot%write('updated_Ru_',shot%dsyn)
call suformat_write('updated_wavelet_'//shot%sindex,shot%wavelet,shot%nt,1,shot%dt)
endif

                call wei%update
                call alloc(shot%dadj,shot%nt,shot%nrcv)

                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
                call kernel_L2sq(shot%dadj)
          
if(s_update_wavelet/='') then
call hud('update adjoint source')
call shot%update_adjsource
endif

            case('L2_scaled'); call hud('0.5|| W(Su - d)||² => adjsrc = S W W(Su - d)')
                if(is_first_in) call alloc(S(i)%scale,shot%nrcv) !then can NOT randomly sample shots..
                do j=1,shot%nrcv
                    if(is_first_in) S(i)%scale(j) = either(0., maxval(abs(shot%dobs(:,j))) / maxval(abs(shot%dsyn(:,j))) , shot%rcv(j)%is_badtrace)
                    shot%dsyn(:,j)=shot%dsyn(:,j)*S(i)%scale(j)
                enddo
                if(is_first_in) then
                    open(12,file=dir_out//'dobs_dsyn_max_ratio',access='direct',recl=4*shot%nrcv)
                    write(12,rec=shot%index) S(i)%scale
                    close(12)
                endif

                !check if S is changing..
    !            if(shot%index==1)   print*, 'on '//shot%sindex,i,S(i)%scale
    !            if(shot%index==112) print*, 'on '//shot%sindex,i,S(i)%scale

if(s_update_wavelet/='') then
call hud('update wavelet')
call shot%update_wavelet(wei_wl%weight) !call matchfilter_apply_to_data(shot%dsyn)
call hud('premute')
shot%dsyn=shot%dsyn*wei%weight / (wei%weight+r_epsilon)
call shot%write('updated_Ru_',shot%dsyn)
call suformat_write('updated_wavelet_'//shot%sindex,shot%wavelet,shot%nt,1,shot%dt)
endif

                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
                call kernel_L2sq(shot%dadj)

if(s_update_wavelet/='') then
call hud('update adjoint source')
call shot%update_adjsource
endif

if(setup%get_bool('L2SCALED_SCALE_ADJSRC',o_default='T')) then
do j=1,shot%nrcv
shot%dadj(:,j)=shot%dadj(:,j)*S(i)%scale(j)
enddo
endif

            case('L2_scaled_filtered'); call hud('0.5|| W(f*Su - d)||² => adjsrc = S f★W W(f*Su - d)')
                if(is_first_in) call alloc(S(i)%scale,shot%nrcv) !then can NOT randomly sample shots..
                do j=1,shot%nrcv
                    if(is_first_in) S(i)%scale(j) = either(0., maxval(abs(shot%dobs(:,j))) / maxval(abs(shot%dsyn(:,j))) , shot%rcv(j)%is_badtrace)
                    shot%dsyn(:,j)=shot%dsyn(:,j)*S(i)%scale(j)
                enddo
                if(is_first_in) then
                    open(12,file=dir_out//'dobs_dsyn_max_ratio',access='direct',recl=4*shot%nrcv)
                    write(12,rec=shot%index) S(i)%scale
                    close(12)
                endif



if(s_update_wavelet/='') then
call hud('update wavelet')
call shot%update_wavelet(wei_wl%weight) !call matchfilter_apply_to_data(shot%dsyn)
call shot%write('updated_Ru_',shot%dsyn)
call suformat_write('updated_wavelet_'//shot%sindex,shot%wavelet,shot%nt,1,shot%dt)
endif

                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
                call kernel_L2sq(shot%dadj)

if(s_update_wavelet/='') then
call hud('update adjoint source')
call shot%update_adjsource
endif

if(setup%get_bool('L2SCALED_SCALE_ADJSRC',o_default='T')) then
do j=1,shot%nrcv
shot%dadj(:,j)=shot%dadj(:,j)*S(i)%scale(j)
enddo
endif


            case('Envsq'); call hud('0.5|| W(E[f*u] - E[d])||² => adjsrc = ...')
if(s_update_wavelet/='') then
call hud('update wavelet')
call shot%update_wavelet(wei_wl%weight) !call matchfilter_apply_to_data(shot%dsyn)
call hud('premute') !0 should be still 0 in the adjoint source, esp below first arrivals
shot%dsyn=shot%dsyn*wei%weight / (wei%weight+r_epsilon)
call shot%write('updated_Ru_',shot%dsyn)
call suformat_write('updated_wavelet_'//shot%sindex,shot%wavelet,shot%nt,1,shot%dt)
endif

                call alloc(Eobs,shot%nt,shot%nrcv)
                call hilbert_envelope(shot%dobs,Eobs,shot%nt,shot%nrcv)
                call shot%write('Eobs_',Eobs)

                fobj%misfit = fobj%misfit &
                    + Envsq(0.5, shot%nt, shot%nrcv, wei%weight, shot%dsyn, Eobs, shot%dt)
                call kernel_Envsq(shot%dadj,shot%nt,shot%nrcv)

if(s_update_wavelet/='') then
call hud('update adjoint source')
call shot%update_adjsource
endif


            case default
                call error('No DNORM specified!')

            end select
        


        call hud('remute') !0 should be still 0 in the adjoint source, esp below first arrivals
        shot%dadj=shot%dadj*wei%weight / (wei%weight+r_epsilon)

        call fld_a%ignite(o_wavelet=shot%dadj)
        call shot%write('dadj_',shot%dadj)


        call hud('----  Solving adjoint eqn & xcorrelate  ----')

        call ppg%init_correlate(a_star_u,'a_star_u')
        call ppg%adjoint(fld_a,fld_u,a_star_u)

if(allocated(a_star_u%gepsr)) call sysio_write('gepsr_'//shot%sindex,a_star_u%gepsr,cb%mz*cb%mx)
if(allocated(a_star_u%gmur )) call sysio_write('gmur_' //shot%sindex,a_star_u%gmur ,cb%mz*cb%mx)
if(allocated(a_star_u%gsgma)) call sysio_write('gsgma_'//shot%sindex,a_star_u%gsgma,cb%mz*cb%mx)

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
