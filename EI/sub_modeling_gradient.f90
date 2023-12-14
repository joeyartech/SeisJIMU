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

    type(t_field) :: fld_u,fld_v, fld_p, fld_q
    type(t_correlate) :: p_star_u,   q_star_v
    ! real,dimension(:,:),allocatable :: Wdres
    !real,dimension(:,:,:),allocatable :: Ddt2
    ! character(:),allocatable :: update_wavelet
    real,dimension(:,:),allocatable :: Esyn,phsyn, Eobs, phobs
    real,dimension(:,:,:),allocatable :: term1, term2

    character(:),allocatable :: dnorm

type :: t_S
    real,dimension(:),allocatable :: scale
end type
type(t_S),dimension(:),allocatable,save :: S
real,dimension(:,:),allocatable :: tmp_Esyn, gmwindow
if(is_first_in) allocate(S(shls%nshots_per_processor)) !then can NOT randomly sample shots..
    

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
            call shot%write('updated_Ru_',shot%dsyn)

            cycle

        endif

        call hud('----  Solving Av=H[s]  ----')
        call shot%read_wlhilb
        call ppg%init_field(fld_v, name='fld_v');    call fld_v%ignite
        call ppg%forward(fld_v)
        call fld_v%acquire; call shot%write('Rv_',shot%dsyn)
        
        shot%dsyn_aux = shot%dsyn
        call fld_u%acquire
        Esyn =sqrt(shot%dsyn**2+shot%dsyn_aux**2); !call shot%write('RE_', Esyn)
        ! ph = (atan2(v,u) -r_pi/2.) *E/(E+1e-4*maxval(E)) +r_pi/2.
        phsyn=atan(shot%dsyn_aux,shot%dsyn);       !call shot%write('Rph_',phsyn)

        if(setup%get_str('JOB')=='forward modeling') cycle


        if (setup%get_bool('IS_DOBS_WAVE',o_default='T')) then
            call hud('Converting observed seismogram (dobs) to envelope data via Hilbert transform')
            Eobs=shot%dobs
            phobs=shot%dobs
            call hilbert_envelope(shot%dobs,Eobs, shot%nt,shot%nrcv)
            call hilbert_phase(   shot%dobs,phobs,shot%nt,shot%nrcv)
            if(mpiworld%is_master) then
                call shot%write('denv_',Eobs)
                call shot%write('dph_',phobs)
            endif
        else
            Eobs=shot%dobs
            phobs=Eobs; phobs=0.
        endif


        call ppg%init_field(fld_p,name='fld_p',ois_adjoint=.true.)
        call ppg%init_field(fld_q,name='fld_q',ois_adjoint=.true.)  

        call hud('----  Computing obj func & dadj  ----')
            call wei%update
            call alloc(shot%dadj,shot%nt,shot%nrcv)

            if(.not.allocated(dnorm)) dnorm=setup%get_str('DATA_NORM','DNORM',o_default='DeltaE')
            select case (dnorm)
                case ('DeltaE')
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, Eobs-Esyn, shot%dt)
                call kernel_L2sq(shot%dadj)
                call fld_p%ignite(o_wavelet=shot%dadj*cos(phsyn))
                call fld_q%ignite(o_wavelet=shot%dadj*sin(phsyn))


case('DeltaE_scaled')
if(is_first_in) call alloc(S(i)%scale,shot%nrcv) !then can NOT randomly sample shots..
call alloc(tmp_Esyn,shot%nt,shot%nrcv)
do j=1,shot%nrcv
    if(is_first_in) S(i)%scale(j) = either(0., maxval(abs(Eobs(:,j))) / maxval(abs(Esyn(:,j))) , shot%rcv(j)%is_badtrace)
    tmp_Esyn(:,j)=Esyn(:,j)*S(i)%scale(j)
enddo
if(is_first_in) then
	open(12,file=dir_out//'Eobs_Esyn_max_ratio',access='direct',recl=4*shot%nrcv)
	write(12,rec=shot%index) S(i)%scale
	close(12)
endif
fobj%misfit = fobj%misfit &
    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, Eobs-tmp_Esyn, shot%dt)
call kernel_L2sq(shot%dadj)
call fld_p%ignite(o_wavelet=shot%dadj*cos(phsyn))
call fld_q%ignite(o_wavelet=shot%dadj*sin(phsyn))
                    

                case ('DeltaE*cos(phobs)')
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*cos(phobs), Eobs-Esyn, shot%dt)
                call kernel_L2sq(shot%dadj)
                call fld_p%ignite(o_wavelet=shot%dadj*cos(phsyn))
                call fld_q%ignite(o_wavelet=shot%dadj*sin(phsyn))

                case ('Deltau')
                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
                call kernel_L2sq(shot%dadj)
                call fld_p%ignite(o_wavelet=shot%dadj)
                call fld_q%ignite(o_wavelet=0.*shot%dadj)

                case default
                call error('No DNORM specified!')

            end select

        call shot%write('dres_',shot%dadj)
        
        ! call hud('----  Solving Aᴴp=Rᴴdadj*cosϕ and Aᴴq=Rᴴdadj*sinϕ  ----')
        ! call hud('----  and crosscorrelate  ----')
        call hud('----  Solving adjoint eqn & xcorrelate  ----')
        call ppg%init_correlate(p_star_u,'p_star_u','gradient')
        call ppg%init_correlate(q_star_v,'q_star_v','gradient')
        call ppg%adjoint(fld_q,fld_p, fld_v,fld_u, q_star_v,p_star_u)

        call hud('---------------------------------')

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

is_first_in=.false.

    !allreduce PFEI misfit values
    call mpi_allreduce(mpi_in_place, [fobj%misfit], 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked EI_misfit '//num2str(fobj%misfit))

    fobj%dnorms=fobj%misfit

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
    
    !scale by shotlist
    call shls%scale(1,o_from_sampled=[fobj%misfit])
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)


    !write correlate
    if(mpiworld%is_master) then
        call q_star_v%write
        call p_star_u%write
    endif


if(dnorm/='Deltau') then
call mpi_allreduce(mpi_in_place, q_star_v%rp_ddsp, m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
call mpi_allreduce(mpi_in_place, p_star_u%rp_ddsp, m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
if(mpiworld%is_master) then
    call q_star_v%write(o_suffix='_stacked')
    call p_star_u%write(o_suffix='_stacked')
    
    term1=q_star_v%rp_ddsp
    term2=p_star_u%rp_ddsp
    den1=sqrt(sum(term1**2,.not.m%is_freeze_zone))
    den2=sqrt(sum(term2**2,.not.m%is_freeze_zone))
    costh=sum(term1/den1*term2/den2,.not.m%is_freeze_zone)

    call hud('Angle between p★u & q★v (w/ mask): '//num2str(acosd(costh))//'°')
endif
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
