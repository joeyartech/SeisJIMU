module m_gradient
use m_model
use m_shotlist
use m_shot
use m_computebox
use m_propagator
use m_field
use m_objectivefunc
use m_matchfilter
use m_smoother_laplacian_sparse

    
    public
    
    character(:),allocatable :: update_wavelet
    
    real fobjective !objective function value
    real,dimension(:,:,:,:),allocatable :: gradient
    
    contains
    
    subroutine gradient_modeling(if_gradient)
        logical,optional :: if_gradient
        

        !assign shots to processors
        if(.not. allocated(shotlist)) call build_shotlist
        
        fobjective=0.
        
        if(present(if_gradient)) then
        if(if_gradient) then
            call alloc(gradient,m%nz,m%nx,m%ny,ncorr)
        endif
        endif
       
        
        call hud('      START LOOP OVER SHOTS          ')
        
        do i=1,nshot_per_processor
        
            call init_shot(i,'data')
            
            call hud('Modeling shot# '//shot%cindex)
            
            call build_computebox
            
            call check_model
            call check_discretization !CFL condition, dispersion etc.
            
            !init propagator and wavefield
            call init_propagator(if_will_do_rfield=.true.)
            
            !*******************************
            !intensive computation
            call propagator_forward(if_will_backpropagate=.true.)
            !*******************************

!write synthetic data
if(mpiworld%is_master) then
open(12,file='synth_raw_'//shot%cindex,access='stream')
write(12) dsyn
close(12)
endif
            
            !update shot%src%wavelet
            !while wavelet in m_propagator is un-touched
            update_wavelet=setup_get_char('UPDATE_WAVELET',default='per shot')
            if(update_wavelet/='no') then
                call gradient_matchfilter_data
                
                !write wavelet updates for QC
                if(mpiworld%is_master) then
                    open(12,file='wavelet_update',access='stream',position='append')
                    write(12) shot%src%wavelet /shot%src%dt*m%cell_volume !scale to be dt, dx independent, see m_shot.f90 ..
                    close(12)
                endif
            endif

!!mitigate some noise due to matchfilter
!where(abs(dobs)<1e-6)
!    dsyn=0.
!endwhere
            
            !write synthetic data
            open(12,file='synth_data_'//shot%cindex,access='stream')
            write(12) dsyn
            close(12)
            
            !fobjective and data residual
            call alloc(dres,shot%rcv(1)%nt,shot%nrcv)
            call objectivefunc_data_norm_residual
            
            if(mpiworld%is_master) write(*,*) 'Shot# 0001: Data misfit norm', dnorm
            
            fobjective=fobjective+dnorm
            
            if(present(if_gradient)) then
            if(if_gradient) then
                
                !adjoint source
                if(update_wavelet/='no') then
                    call matchfilter_correlate_filter_residual(shot%src%nt,shot%nrcv,dres)
                endif
                
                call alloc(cb%gradient,cb%mz,cb%mx,cb%my,ncorr) !(:,:,:,1) is glda, (:,:,:,2) is gmu, (:,:,:,3) is grho0
                !*******************************
                !intensive computation
                call propagator_adjoint(gradient=cb%gradient)
                !*******************************
                
                !put cb%gradient into global gradient
                gradient(cb%ioz:cb%ioz+cb%mz-1,&
                         cb%iox:cb%iox+cb%mx-1,&
                         cb%ioy:cb%ioy+cb%my-1,:) = &
                gradient(cb%ioz:cb%ioz+cb%mz-1,&
                         cb%iox:cb%iox+cb%mx-1,&
                         cb%ioy:cb%ioy+cb%my-1,:) + cb%gradient(:,:,:,:)
            
            endif
            endif
            
        enddo
        
        call hud('        END LOOP OVER SHOTS        ')
        
        !call objectivefunc_model_norm ...
        fobjective=fobjective!+lambda*mnorm
        
        !collect global objective function value
        call mpi_barrier(mpiworld%communicator, mpiworld%ierr)
        call mpi_allreduce(MPI_IN_PLACE, fobjective, 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        
        if(present(if_gradient)) then
        if(if_gradient) then
        
            !collect global gradient
            call mpi_barrier(mpiworld%communicator, mpiworld%ierr)
            call mpi_allreduce(MPI_IN_PLACE, gradient, m%n*ncorr, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)

            !smoothing
            if(setup_get_logical('IF_SMOOTHING',default=.true.)) then
                call hud('Initialize Laplacian smoothing')
                call init_smoother_laplacian([m%nz,m%nx,m%ny],[m%dz,m%dx,m%dy],shot%src%fpeak)
                do icorr=1,ncorr
                    call smoother_laplacian_extend_mirror(gradient(:,:,:,icorr),m%itopo)
                    call smoother_laplacian_pseudo_nonstationary(gradient(:,:,:,icorr),m%vp)
                enddo
            endif
            
            !mask gradient
            do iy=1,m%ny
            do ix=1,m%nx
                gradient(1:m%itopo(ix,iy)-1,ix,iy,:) =0.
            enddo
            enddo
            
        endif
        endif
        
    end subroutine
    
    subroutine gradient_matchfilter_data
    
        if(update_wavelet=='stack') then
            !average wavelet across all processors
            !note: if more shots than processors, non-assigned shots will not contribute to this averaging
            call matchfilter_estimate(shot%src%nt,shot%nrcv,dsyn,dobs,shot%index,if_stack=.true.)
        else
            call matchfilter_estimate(shot%src%nt,shot%nrcv,dsyn,dobs,shot%index,if_stack=.false.)
        endif
        
        call matchfilter_apply_to_wavelet(shot%src%nt,shot%src%wavelet)
        
        call matchfilter_apply_to_data(shot%src%nt,shot%nrcv,dsyn)
        
    end subroutine
    
end
