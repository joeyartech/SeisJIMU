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
    character(:),allocatable :: s_conversion
    real :: dtr
    
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
dtr=setup%get_real('CHANNEL_SPACING',o_alias='DTR',o_default=num2str(shot%rcv(2)%x-shot%rcv(1)%x)) !assume constant DAS channel spacing

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

            case ('L2averaged'); call hud('0.5|| W(f*MA*u - d)||² => adjsrc = MA★f★W W(f*u - d)')

                length=nint( setup%get_real('MOVING_AVERAGE_LENGTH','MA_LEN',o_mandatory=1)/dtr )
                if(length>0) call moving_average(shot%dsyn,length)
                call shot%write('avg_dsyn_',shot%dsyn)

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

                if(length>0) call moving_average(shot%dadj,length)
                call fld_a%ignite(o_wavelet=shot%dadj)
                call shot%write('dadj_',shot%dadj)


            case ('strain_L2averaged'); call hud('0.5|| W(f*MA*conv*u - d)||² => adjsrc = conv★MA★f★W W(f*u - d)')

                s_conversion=setup%get_str('DATA_CONVERSION_METHOD',o_default='fk')

                !first goto strain
                if(index(ppg%info,'Momemtum-Strain')>0) then !m_propagator_DAS
                
                elseif(index(ppg%info,'Velocity-Stress')>0) then !m_propagator_PSV
                    if(s_conversion=='fk') then
                        call convert_in_fk(shot%dsyn,'k/w',dtr) !call shot%write('k_w_dsyn',shot%dsyn)

                    else ! s_conversion=='tx'
                        call hud('integrate_t(differentiate_x(shot%dsyn))')
                        call differentiate_x(shot%dsyn,dtr)
                        call integrate_t(shot%dsyn) !call shot%write('diffx_intt_dsyn',shot%dsyn)

                    endif

                elseif(index(ppg%info,'Displacement-Strain formulation')>0) then !m_propagator_PSV2ndUE
                
                elseif(index(ppg%info,'Displacement formulation')>0) then !m_propagator_PSV2nd
                    if(s_conversion=='fk') then
                        call convert_in_fk(shot%dsyn,'ik',dtr) !call shot%write('ik_dsyn',shot%dsyn)

                    else ! s_conversion=='tx'
                        call hud('differentiate_x(shot%dsyn)')
                        call differentiate_x(shot%dsyn,dtr) !call shot%write('diffx_dsyn',shot%dsyn)
                    endif

                endif


                !then add gauge length
                length=nint( setup%get_real('MOVING_AVERAGE_LENGTH','MA_LEN',o_mandatory=1)/dtr )
                if(length>0) call moving_average(shot%dsyn,length)
                call shot%write('DAS_dsyn_',shot%dsyn)


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

                if(length>0) call moving_average(shot%dadj,length)


                !finally convert back
                if(index(ppg%info,'Momemtum-Strain')>0) then !m_propagator_DAS
                
                elseif(index(ppg%info,'Velocity-Stress')>0) then !m_propagator_PSV
                    if(s_conversion=='fk') then
                        call convert_in_fk(shot%dadj,'k/w',dtr) !call shot%write('k_w_dadj',shot%dadj)

                    else ! s_conversion=='tx'
                        ! call hud('rev_integrate_t(differentiate_x(shot%dadj))')
                        call differentiate_x(shot%dadj,dtr)
                        ! call rev_integrate_t(shot%dadj)
                        call integrate_t(shot%dadj) !why? 
                        !call shot%write('diffx_intt_dadj',shot%dadj)

                    endif

                elseif(index(ppg%info,'Displacement-Strain formulation')>0) then !m_propagator_PSV2ndUE
                
                elseif(index(ppg%info,'Displacement formulation')>0) then !m_propagator_PSV2nd
                    if(s_conversion=='fk') then
                        shot%dadj=-shot%dadj
                        call convert_in_fk(shot%dadj,'ik',dtr) !call shot%write('ik_ndadj',shot%dadj)

                    else ! s_conversion=='tx'
                        call hud('differentiate_x(-shot%dadj)')
                        !call differentiate_x(-shot%dadj) !this is wrong in fortran..
                        shot%dadj=-shot%dadj
                        call differentiate_x(shot%dadj,dtr) !so let's use functions instead of subroutines..
                        !call shot%write('diffx_ndadj',shot%dadj)

                    endif

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

            case('L2averaged_filtered'); call hud('0.5|| W(f*MA*u - d)||² => adjsrc = MA★f★W W(f*MA*u - d)')
                length=nint( setup%get_real('MOVING_AVERAGE_LENGTH','MA_LEN',o_mandatory=1)/dtr )
                if(length>0) call moving_average(shot%dsyn,length)
                call shot%write('avg_dsyn_',shot%dsyn)

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

                if(length>0) call moving_average(shot%dadj,length)
                call fld_a%ignite(o_wavelet=shot%dadj)
                call shot%write('dadj_',shot%dadj)

            case('L2averaged_scaled_filtered'); call hud('0.5|| W(f*S MA*u - d)||² => adjsrc = MA★S f★W W(f*S MA*u - d)')
                length=nint( setup%get_real('MOVING_AVERAGE_LENGTH','MA_LEN',o_mandatory=1)/dtr )
                if(length>0) call moving_average(shot%dsyn,length)
                call shot%write('avg_dsyn_',shot%dsyn)

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

                if(length>0) call moving_average(shot%dadj,length)
                call fld_a%ignite(o_wavelet=shot%dadj)
                call shot%write('dadj_',shot%dadj)
                


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

if(allocated(a_star_u%gepsr )) call sysio_write('gepsr_' //shot%sindex,a_star_u%gepsr ,m%n)
if(allocated(a_star_u%gmur  )) call sysio_write('gmur_'  //shot%sindex,a_star_u%gmur  ,m%n)
if(allocated(a_star_u%gsgmar)) call sysio_write('gsgmar_'//shot%sindex,a_star_u%gsgmar,m%n)

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

    subroutine differentiate_x(data,dtr) !by central diff
        real,dimension(:,:) :: data !nt x nrcv
        real :: dtr

        real,dimension(:,:), allocatable :: dout
        call alloc(dout,shot%nt,shot%nrcv)

        do ir=2,shot%nrcv-1
            dout(:,ir) = (data(:,ir+1)-data(:,ir-1))
        enddo

        !padding
        dout(:,1)=dout(:,2)
        dout(:,shot%nrcv)=dout(:,shot%nrcv-1)

        data=dout /2./dtr  !assuming const interval..

        deallocate(dout)

    end subroutine

    subroutine integrate_t(data) !integrate
        real,dimension(:,:) :: data !nt x nrcv

        real,dimension(:,:), allocatable :: dout
        call alloc(dout,shot%nt,shot%nrcv)

        do ir=1,shot%nrcv
                dout(1,ir)=data(1,ir)
            do it=2,shot%nt
                dout(it,ir) = dout(it-1,ir) + data(it,ir)
            enddo
        enddo

        data=dout *shot%dt

        deallocate(dout)

    end subroutine

    subroutine rev_integrate_t(data) !reverse-time integrate
        real,dimension(:,:) :: data !nt x nrcv

        real,dimension(:,:), allocatable :: dout
        call alloc(dout,shot%nt,shot%nrcv)

        do ir=1,shot%nrcv
                dout(shot%nt,ir)=data(shot%nt,ir)
            do it=shot%nt-1,1,-1
                dout(it,ir) = dout(it+1,ir) + data(it,ir)
            enddo
        enddo

        data=dout *shot%dt

        deallocate(dout)

    end subroutine

    subroutine convert_in_fk(data,op,dtr)
    use m_math
    use singleton
        real,dimension(:,:) :: data !nt x nrcv
        character(*) :: op
        real :: dtr

        real :: w(shot%nt), k(shot%nrcv)

        complex,dimension(:,:),allocatable :: filter
        complex(fftkind),dimension(:,:),allocatable :: data_fft

        allocate(filter(shot%nt,shot%nrcv))

        if(allocated(data_fft)) deallocate(data_fft)
        allocate(data_fft(shot%nt,shot%nrcv))

        n=shot%nt
        if(mod(n,2)==0) then !if n is even, 1 is DC; 2:n/2 are +f; n/2+1:n are -f
            w(1:n/2    )= [(i,i=1,n/2)]-1
            w(  n/2+1:n)= -w(n/2:1:-1)
        else !if n is odd, 1 is DC, 2:(n+1)/2 are +f, (n+1)/2+1:n are -f
            w(1:(n+1)/2    )= [(i,i=1,(n+1)/2)]-1
            w(  (n+1)/2+1:n)= -w((n+1)/2:2:-1)
        endif

        w=w*2*r_pi/(n-1)/shot%dt

        n=shot%nrcv
        if(mod(n,2)==0) then !if n is even, 1 is DC; 2:n/2 are +f; n/2+1:n are -f
            k(1:n/2    )= [(i,i=1,n/2)]-1
            k(  n/2+1:n)= -k(n/2:1:-1)
        else !if n is odd, 1 is DC, 2:(n+1)/2 are +f, (n+1)/2+1:n are -f
            k(1:(n+1)/2    )= [(i,i=1,(n+1)/2)]-1
            k(  (n+1)/2+1:n)= -k((n+1)/2:2:-1)
        endif

        k=k*2*r_pi/(n-1)/dtr  !assuming const interval..

        call hud('convert_in_fk op: '//op)

        select case (op)
        case ('k/w') 
            ! eps=maxval(w)*1e-5
            eps=w(2) !smallest w

            do ik=1,shot%nrcv; do iw=1,shot%nt
                filter(iw,ik) = k(ik)*w(iw) / (w(iw)*w(iw)+eps)
            enddo; enddo

        case ('ik')
            ! eps=maxval(w)*1e-5
            eps=w(2) !smallest w

            do ik=1,shot%nrcv
                filter(:,ik) = c_i*k(ik)
            enddo

        case ('w/k')
            ! eps=maxval(k)*1e-5
            eps=k(2) !smallest k

            do ik=1,shot%nrcv; do iw=1,shot%nt
                filter(iw,ik) = w(iw)*k(ik) / (k(ik)*k(ik)+eps)
            enddo; enddo

        case ('1/ik')
            ! eps=maxval(k)*1e-5
            eps=k(2) !smallest k

            do ik=1,shot%nrcv
                filter(:,ik) = 1/(c_i*k(ik)+eps)
            enddo
            
        endselect

        data_fft = fft2d(dcmplx(data))
        
        data=real(fft2d(filter*data_fft,inv=.true.),kind=4)

        deallocate(filter)

    end subroutine

end subroutine
