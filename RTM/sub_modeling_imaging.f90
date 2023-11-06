subroutine modeling_imaging
use mpi
use m_System
use m_Modeling
use m_hilbert
!use m_weighter

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_v, fld_u, fld_q, fld_p
    type(t_correlate) :: a_star_u
    real,dimension(:,:),allocatable :: tmp

    call alloc(correlate_image,m%nz,m%nx,m%ny,ppg%nimag)
    
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
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire; call shot%write('Ru_',shot%dsyn)

        call hud('----  Solving Av=H[s]  ----')
        call shot%read_wlhilb
        call ppg%init_field(fld_v,name='fld_v');    call fld_v%ignite
        call ppg%forward(fld_v)
        call fld_v%acquire; call shot%write('Rv_',shot%dsyn); shot%dsyn_aux=shot%dsyn
        call fld_u%acquire

        if(setup%get_str('JOB')=='forward modeling') cycle

        if(setup%get_str('JOB')=='imaging') then
            call hud('-----------------------')
            call hud('        Imaging        ')
            call hud('-----------------------')

            call ppg%init_field(fld_p,name='fld_p',ois_adjoint=.true.)
            call ppg%init_field(fld_q,name='fld_q',ois_adjoint=.true.)

            if(setup%get_str('RTM_ADJSRC',o_default='dobs')=='dobs') then
                shot%dadj=shot%dobs
            else
                shot%dadj=shot%dobs-shot%dsyn
            endif

            call fld_p%ignite(o_wavelet=shot%dadj)
                
            call alloc(tmp,shot%nt,shot%nrcv)
            call hilbert_transform(shot%dadj,tmp,shot%nt,shot%nrcv)
            call fld_q%ignite(o_wavelet=tmp)

            call ppg%init_correlate(a_star_u,'a_star_u')

            call hud('----  Solving adjoint eqn & xcorrelate  ----')
            !Aᴴa = -Rᴴd
            call ppg%adjoint_poynting(fld_q,fld_p,fld_v,fld_u,a_star_u)

            call hud('----  Assemble  ----')
            call ppg%assemble(a_star_u)

            call hud('---------------------------------')

        endif
        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    if(setup%get_str('JOB')=='forward modeling') then
        call mpiworld%final
        stop
    endif


    if(mpiworld%is_master) call a_star_u%write

    !collect global correlations
    call mpi_allreduce(mpi_in_place,correlate_image, m%n*ppg%nimag, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    if(mpiworld%is_master) call sysio_write('correlate_image',correlate_image,m%n*ppg%nimag)

    !write correlate
    if(mpiworld%is_master) then
        call a_star_u%write

    endif

    ! if(ppg%if_compute_engy) then
    !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    ! endif
        
    call mpiworld%barrier

end subroutine