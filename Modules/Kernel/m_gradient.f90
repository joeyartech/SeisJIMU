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
    
    subroutine gradient

!if nshot_per_processor is not same for each processor,
!update_wavelet='stack' mode will be stuck due to collective communication in m_matchfilter.f90
if(nshot_per_processor * mpiworld%nproc /= nshots) then
    fatal('Unequal shot numbers on processors. If you are using UPDATE_WAVELET=''stack'', the code will be stuck due to collective communication in m_matchfilter')
endif

        
        call alloc(gradient,m%nz,m%nx,m%ny,ncorr)
        
        call fobj%init
        
        call hud('===== START LOOP OVER SHOTS =====')
        
        do ishot=1,shotlist%nshot_per_processor        

            call shot%init

            call hud('Modeling shot# '//shot%sindex)
            
            call cb%init
            call cb%project

            call sfield%add_source

            call sfield%check_model
            call sfield%check_discretization
            
            call sfield%init

            call sfield%forward
            
            call shot%acquire(sfield)

            if(mpiworld%is_master) call suformat_write(file='synth_raw_'//shot%sindex,shot%dsyn)
            
            update_wavelet=setup_get_char('UPDATE_WAVELET',default='per shot')
            if(update_wavelet/='no') then
                call gradient_matchfilter_data
                call suformat_write(iproc=0,file='wavelet_update',append=.true.)
            endif
            
            !write synthetic data
            call suformat_write(file='synth_data_'//shot%sindex,shot%dsyn)
                        
            !objective function
            call fobj%compute_dnorm(shot)
        
            !adjoint source
            if(update_wavelet/='no') then
                call matchfilter_correlate_filter_residual(shot%src%nt,shot%nrcv,shot%dres)
            endif
            
            call alloc(cb%image,cb%mz,cb%mx,cb%my,ncorr) !(:,:,:,1) is glda, (:,:,:,2) is gmu, (:,:,:,3) is grho0

            call rfield%backward(sfield=sfield,gradient=cb%image)

            !put cb%gradient into global gradient
            call cb%project_back(gradient,cb%gradient)
            
        enddo
        
        call hud('        END LOOP OVER SHOTS        ')
        
        !collect global objective function value
        call mpiworld%barrier
        call mpi_allreduce(MPI_IN_PLACE, fobj%dnorm, fobj%n_dnorm, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        call fobj%print_dnorm

        !collect global gradient
        call mpiworld%barrier
        call mpi_allreduce(MPI_IN_PLACE, gradient, m%n*ncorr, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)

        !call objectivefunc_model_norm ...
        call fobj%compute_mnorm(m)

        !smoothing
        if(setup%get_bool('IF_SMOOTHING',default=.true.)) then
            call hud('Initialize Laplacian smoothing')
            call init_smoother_laplacian([m%nz,m%nx,m%ny],[m%dz,m%dx,m%dy],shot%src%fpeak)
            do icorr=1,ncorr
                call smoother_laplacian_extend_mirror(gradient(:,:,:,icorr),m%itopo)
                call smoother_laplacian_pseudo_nonstationary(gradient(:,:,:,icorr),m%vp)
            enddo
        endif
        
        !mask gradient
        call mask_gradient(gradient)
        
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

    subroutine mask_gradient

        !bathymetry as hard mask
        do iy=1,m%ny
        do ix=1,m%nx
            gradient(1:m%itopo(ix,iy)-1,ix,iy,:) =0.
        enddo
        enddo

        !soft mask
        inquire(file='grad_mask', exist=alive)
        if(alive) then
            call alloc(tmp_grad_mask,m%nz,m%nx,m%ny)
            open(12,file='grad_mask',access='direct',recl=4*m%n,action='read',status='old')
            read(12,rec=1) tmp_grad_mask
            close(12)
            call hud('grad_mask is read. gradient is masked in addition to topo.')
            
            do i=1,size(gradient,4)
                gradient(:,:,:,i)=gradient(:,:,:,i)*tmp_grad_mask
            enddo
        endif

    end subroutine
    
end
