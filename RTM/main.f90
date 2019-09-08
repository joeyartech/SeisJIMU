program main
use m_model
use m_computebox
use m_shotlist
use m_shot
use m_field
use m_propagator
    
    real,dimension(:,:,:,:),allocatable :: image
    
    call init_mpiworld

    call hud('===================================')
    call hud('     WELCOME TO LEGO RTM CODE      ')
    call hud('===================================')
    
    call init_setup(istat)
    
    if(istat==0) then !print manual
        !call fwd_print_manual
        stop
    endif
    
    
    !read initial model
    call init_model
    call field_print_info
    
    nshots=get_setup_int('NSHOTS',default=1)
    
    !assign shots to processors
    call build_shotlist(nshots)
    
    call alloc(image,m%nz,m%nx,m%ny,2)
    
    
    call hud('      START LOOP OVER SHOTS          ')

    do i=1,nshot_per_processor
        
        call init_shot(i,'data')
        
        call hud('Modeling shot# '//shot%cindex)
        
        call build_computebox
        
        call check_model
        call check_discretization !CFL condition, dispersion etc.
        
        !propagation and wavefield
        call init_propagator(if_will_do_rfield=.true.)
        
        !*******************************
        !intensive computation
        call propagator_forward(if_will_backpropagate=.true.)
        !*******************************
        
!         !update shot%src%wavelet
!         !while wavelet in m_propagator is un-touched
!         update_wavelet=get_setup_char('UPDATE_WAVELET',default='per shot')
!         if(update_wavelet/='no') then
!             call gradient_matchfilter_data
!             
!             !write wavelet updates for QC
!             if(mpiworld%is_master) then
!                 open(12,file='wavelet_update',access='stream',position='append')
!                 write(12) shot%src%wavelet
!                 close(12)
!             endif
!         endif
        
        !write synthetic data
        open(12,file='synth_data_'//shot%cindex,access='stream')
        write(12) dsyn
        close(12)
        
        !adjoint source
        call alloc(dres,shot%rcv(1)%nt,shot%nrcv)
        dres=dobs-dsyn
!         if(update_wavelet/='no') then
!             call matchfilter_correlate_filter_residual(shot%src%nt,shot%nrcv,dres)
!         endif
        
        call alloc(cb%image,cb%mz,cb%mx,cb%my,2)
        !*******************************
        !intensive computation
        call propagator_adjoint(image=cb%image)
        !*******************************
        
        !put cb%gradient into global gradient
        image(cb%ioz:cb%ioz+cb%mz-1,&
              cb%iox:cb%iox+cb%mx-1,&
              cb%ioy:cb%ioy+cb%my-1,:) = &
        image(cb%ioz:cb%ioz+cb%mz-1,&
              cb%iox:cb%iox+cb%mx-1,&
              cb%ioy:cb%ioy+cb%my-1,:) + cb%image(:,:,:,:)
        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')
    
    !collect global gradient
    call mpi_barrier(mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(MPI_IN_PLACE, image, m%n*2, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    open(12,file='corr_sr',action='write',access='direct',recl=4*m%n)
    write(12,rec=1) image(:,:,:,1)
    close(12)
    open(12,file='corr_ss',action='write',access='direct',recl=4*m%n)
    write(12,rec=1) image(:,:,:,2)
    close(12)
    open(12,file='corr_div',action='write',access='direct',recl=4*m%n)
    write(12,rec=1) image(:,:,:,1) / (image(:,:,:,2)+1e-8)
    close(12)
    
    call mpiworld_finalize
    
end