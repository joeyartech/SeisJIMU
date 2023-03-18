subroutine modeling_imaging
use mpi
use m_System
use m_Modeling
!use m_weighter

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_u, fld_a
    type(t_correlation) :: a_star_u
    
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

        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite

        !forward modeling
        !A u= s
        call ppg%forward(fld_u)
        call fld_u%acquire; call shot%write('Ru',shot%dsyn)

        if(setup%get_str('JOB')=='forward modeling') cycle

        if(setup%get_str('JOB')=='imaging') then
            call hud('-----------------------')
            call hud('        Imaging        ')
            call hud('-----------------------')

            call a_star_u%init('a_star_u',shape123_from='model') !F₂★E₀ for gtDt
            call alloc(a_star_u%drp_dt_dsp_dt,m%nz,m%nx,m%ny,1)
            call alloc(a_star_u%nab_rp_nab_sp,m%nz,m%nx,m%ny,1)

            !adjoint modeling
            !Aᴴa = -Rᴴd
            call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.); call fld_a%ignite(shot%dobs)
            call ppg%adjoint(fld_a,fld_u,o_a_star_u=a_star_u)
            call hud('---------------------------------')

            !call cb%project_back
            correlation_image(cb%ioz:cb%ioz+cb%mz-1,&
                  cb%iox:cb%iox+cb%mx-1,&
                  cb%ioy:cb%ioy+cb%my-1,1) = &
            correlation_image(cb%ioz:cb%ioz+cb%mz-1,&
                  cb%iox:cb%iox+cb%mx-1,&
                  cb%ioy:cb%ioy+cb%my-1,1) + a_star_u%rp_sp(:,:,:,1)

        endif
        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    if(setup%get_str('JOB')=='forward modeling') then
        call mpiworld%final
        stop
    endif

    !collect

    !collect global correlations
    call mpi_allreduce(mpi_in_place,correlation_image, m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    if(mpiworld%is_master) then
        call sysio_write('correlation_image' ,correlation_image,m%n)
    endif

    ! if(ppg%if_compute_engy) then
    !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    ! endif
        
    call mpiworld%barrier

end subroutine