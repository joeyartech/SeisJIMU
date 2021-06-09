module m_gradient_modeling
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
                
        call hud('===== START LOOP OVER SHOTS =====')
        
        do ishot=1,shotlist%nshot_per_processor        

            call shot%init

            call hud('Modeling shot# '//shot%sindex)
            
            call cb%init
            call cb%project

            call propagator%check_model
            call propagator%check_discretization
            
            call propagator%init

            call propagator%init_field(sfield)

            call sfield%add_RHS(shot%source)
            
            call propagator%forward(sfield,dsyn=shot%dsyn)
            !call shot%acquire(sfield)

            if(mpiworld%is_master) call suformat_write(file='synth_raw_'//shot%sindex,shot%dsyn)
            
            update_wavelet=setup_get_char('UPDATE_WAVELET',default='per shot')
            if(update_wavelet/='no') then
                call shot%source_estim !call gradient_matchfilter_data
            endif
            
            !write synthetic data
            call suformat_write(file='synth_data_'//shot%sindex,shot%dsyn)
                        
            !objective function
            call fobj%compute_dnorms(shot)
        
            !adjoint source
            if(update_wavelet/='no') then
                call shot%adjsource_update
            endif

            call propagator%init_field(rfield)

            call rfield%add_RHS(shot%adjsource)
            
            call alloc(cb%kernel,cb%mz,cb%mx,cb%my,ncorr) !(:,:,:,1) is glda, (:,:,:,2) is gmu, (:,:,:,3) is grho0

            call propagator%reverse(rfield,o_sfield=sfield,o_grad=cb%kernel,o_rdt=shot%dt_Nyquist)

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
        call mpi_allreduce(MPI_IN_PLACE, fobj%gradient, m%n*ncorr, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)

    end subroutine


    subroutine gradient_modeling_approximate
    end subroutine

end
