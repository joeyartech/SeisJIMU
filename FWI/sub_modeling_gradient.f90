subroutine modeling_gradient
use mpi
use m_System
use m_Modeling
use m_weighter
use m_Lpnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse

    logical,save :: is_first_in=.true.
    type(t_field) :: fld_u, fld_a
    type(t_correlation) :: a_star_u
    character(:),allocatable :: update_wavelet

    real,dimension(:,:,:),allocatable :: gikpa, gbuo
    
    fobj%misfit=0.
    
    call alloc(gikpa,m%nz,m%nx,m%ny)
    call alloc(gbuo ,m%nz,m%nx,m%ny)
    
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor

        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        call hud('Modeling '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project

        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer

        call ppg%init_field(fld_u,name='fld_u');  call fld_u%ignite
        
        !forward modeling
        call ppg%forward(fld_u);   call fld_u%acquire

        if(mpiworld%is_master) call shot%write('draw_',shot%dsyn)

        update_wavelet=setup%get_str('UPDATE_WAVELET',o_default='per shot')

        call wei%update
        
        if(update_wavelet/='no') call shot%update_wavelet(wei%weight) !call gradient_matchfilter_data
        
        !write synthetic data
        call shot%write('dsyn_',shot%dsyn)

        !data weighting
        call wei%update
        
        !Adjoint state method with Lagrangian formulation
        !to compute the gradient (L2 norm example)
        !C = ½║u║² = ½∫ (u-d)² δ(x-xr) dtdx³
        !KᵤC = (u-d)δ(x-xr)
        !L = C + <a|Au-s> ≐ C + <Aᴴa|u>
        !0 = KᵤL = KᵤC + Aᴴa => Aᴴa = -KᵤC = (d-u)δ(x-xr)
        !KₘL = aᴴ KₘA u =: a★Du

        call a_star_u%init('a_star_u',shape123_from='model')
        call alloc(a_star_u%rp_lap_sp,m%nz,m%nx,m%ny,1)
        call alloc(a_star_u%nab_rp_nab_sp,m%nz,m%nx,m%ny,1)

        ! !objective function and adjoint source
        ! call fobj%stack_dnorms
        ! shot%dadj=-shot%dadj
        call wei%update!('_4IMAGING')
        fobj%misfit = fobj%misfit &
            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)

        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call kernel_L2sq(shot%dadj)
        call shot%write('dadj_',shot%dadj)

        
        if(mpiworld%is_master) call fobj%print_dnorms('Shotloop-stacked','upto '//shot%sindex)
        ! if(either(oif_gradient,.true.,present(oif_gradient))) then
        
            !adjoint source
            if(update_wavelet/='no') call shot%update_adjsource
            
            call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)
            
            call fld_a%ignite
            
            !adjoint modeling
            call ppg%adjoint(fld_a,fld_u,o_a_star_u=a_star_u)
            
            !call cb%project_back
            gikpa(cb%ioz:cb%ioz+cb%mz-1,&
                 cb%iox:cb%iox+cb%mx-1,&
                 cb%ioy:cb%ioy+cb%my-1) = &
            gikpa(cb%ioz:cb%ioz+cb%mz-1,&
                 cb%iox:cb%iox+cb%mx-1,&
                 cb%ioy:cb%ioy+cb%my-1) + m%vp**2*m%rho*a_star_u%rp_lap_sp(:,:,:,1)

            gbuo(cb%ioz:cb%ioz+cb%mz-1,&
                 cb%iox:cb%iox+cb%mx-1,&
                 cb%ioy:cb%ioy+cb%my-1) = &
            gbuo(cb%ioz:cb%ioz+cb%mz-1,&
                 cb%iox:cb%iox+cb%mx-1,&
                 cb%ioy:cb%ioy+cb%my-1) + a_star_u%nab_rp_nab_sp(:,:,:,1)
            
        ! endif

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    !collect global objective function value
    call mpi_allreduce(mpi_in_place, [fobj%misfit], 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked FWI_misfit '//num2str(fobj%misfit))

    fobj%dnorms=fobj%misfit

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')

        ! if(either(oif_gradient,.true.,present(oif_gradient))) then
        !collect global gradient
        call mpi_allreduce(mpi_in_place, gikpa, m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        call mpi_allreduce(mpi_in_place, gbuo,  m%n, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        
        !scale by shotlist
        call shls%scale(1,o_from_sampled=[fobj%misfit])
        call shls%scale(m%n,o_from_sampled=gikpa)
        call shls%scale(m%n,o_from_sampled=gbuo)

        ! if(ppg%if_compute_engy) then
        !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        ! endif
    ! endif

    if(mpiworld%is_master) then
        call sysio_write('gikpa',gikpa,m%n)
        call sysio_write('gbuo', gbuo, m%n)
    endif

    call alloc(correlation_gradient,m%nz,m%nx,m%ny,ppg%ngrad)
    correlation_gradient(:,:,:,1)=gikpa
    correlation_gradient(:,:,:,2)=gbuo

    call mpiworld%barrier

end subroutine
