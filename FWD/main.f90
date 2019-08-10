program main
use m_model
use m_computebox
use m_gen_acquisition
use m_shotlist
use m_shot
use m_propagator
use m_field

    call init_mpiworld

    call hud('========================================')
    call hud('     WELCOME TO LEGO MODELING CODE      ')
    call hud('========================================')
    
    call init_setup(istat)
    
    if(istat==0) then !print manual
        !call fwd_print_manual
        stop
    endif
    
    
    !read initial model
    call init_model
    call field_print_info
    
    !generate acquisition and source wavelet
    call gen_acquisition
    
    !assign shots to processors
    call build_shotlist(acqui%nsrc)
    
    call hud('      START LOOP OVER SHOTS          ')

    do i=1,nshot_per_processor
        
        call init_shot(i,'setup')
        
        call hud('Modeling shot# '//shot%cindex)
        
        call build_computebox
        
        call check_model
        call check_discretization !CFL condition, dispersion etc.
    
        !propagator and wavefield
        call init_propagator
        
        !*******************************
        !intensive computation
        call propagator_forward
        !*******************************
        
        !write synthetic data
        open(12,file='synth_data_'//shot%cindex,access='stream')
        write(12) dsyn
        close(12)
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')
    
    call mpiworld_finalize
    
end