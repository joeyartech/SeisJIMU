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

    contains
    
    subroutine gradient_modeling(fobj)
        
        call alloc(fobj%gradient,m%nz,m%nx,m%ny,ncorr)
        call alloc(G1,m%nz,m%nx,m%ny,ncorr)
                
        call hud('===== START LOOP OVER SHOTS =====')
        
        do ishot=1,shotlist%nshot_per_processor        

            call shot%init

            call hud('Modeling shot# '//shot%sindex)

            !G1 term            
            call cb%init
            call cb%project

            call sfield%add_source(shot%source)

            call sfield%check_model
            call sfield%check_discretization
            
            call sfield%init

            call sfield%forward
            
            call shot%acquire(sfield)

            if(mpiworld%is_master) call suformat_write(file='synth_raw_'//shot%sindex,shot%dsyn)
            
            update_wavelet=setup_get_char('UPDATE_WAVELET',default='per shot')
            if(update_wavelet/='no') then
                call shot%source_estim !call gradient_matchfilter_data
            endif
            
            !write synthetic data
            call suformat_write(file='synth_data_'//shot%sindex,shot%dsyn)
            
            !apply separation
            call data_sepa(shot,'rfl')

            !objective function
            call fobj%compute_dnorms(shot); F1=fobj%dnorms
        
            !adjoint source
            if(update_wavelet/='no') then
                call shot%adjsource_update
            endif

            call rfield%add_source(shot%adjsource)

            call alloc(cb%kernel,cb%mz,cb%mx,cb%my,ncorr) !(:,:,:,1) is glda, (:,:,:,2) is gmu, (:,:,:,3) is grho0

            call rfield%backward(sfield=sfield,xcorr=G1,dt_Nyquist=shot%dt_Nyquist)
            
            !G2 term
            call cb%init
            call cb%project

            call sfield%add_source(shot%source)

            call sfield%init
            call sfield%forward

            call shot%acquire(sfield)

            !apply  separation
            call data_sepa(shot,'old_rfl-div')
                        
            !objective function
            call fobj%compute_dnorms(shot); F2=fobj%dnorms

            call rfield%add_source(shot%adjsource)
            
            call alloc(cb%kernel,cb%mz,cb%mx,cb%my,ncorr) !(:,:,:,1) is glda, (:,:,:,2) is gmu, (:,:,:,3) is grho0

            call rfield%backward(sfield=sfield,xcorr=G2,dt_Nyquist=shot%dt_Nyquist)

            !Gradient=G1+G2
            cb%kernel=G1+G2

            !put cb%kernel into global gradient
            call cb%project_back(fobj%gradient,cb%kernel)
            
        enddo
        
        call hud('        END LOOP OVER SHOTS        ')
        
        !collect global objective function value
        call mpiworld%barrier
        call mpi_allreduce(MPI_IN_PLACE, fobj%dnorm, fobj%n_dnorm, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        call fobj%print_dnorm

        !collect global gradient
        call mpiworld%barrier
        call mpi_allreduce(MPI_IN_PLACE, gradient, m%n*ncorr, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)

    end subroutine

end
