program main
use m_gradient
use m_preconditioner
use m_parameterization
use m_optimizer

    character(:),allocatable :: job
    
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
    
    call gradient_modeling(if_gradient=.true.)
    

    open(12,file='gradient',action='write',access='stream')
    write(12,rec=1) gradient
    close(12)
        
    
!     open(12,file='precond_gradient',action='write',access='direct',recl=4*m%n*2)
!     write(12,rec=1) precond_gradient
!     close(12)
    
    ! job=get_setup_char('JOB',default='optimization')
    ! !if just estimate the wavelet or compute the gradient then this is it.
    ! if(job/='optimization') then
    !     call mpiworld_finalize
    !     stop
    ! endif
    
    
    !initialize parameterization
    call init_parameterization
    
    !initialize optimizer
    call init_optimizer(m%n*npar)
    
    open(12,file='initial%pg',action='write',access='direct',recl=4*m%n*npar)
    write(12,rec=1) current%pg
    close(12)

    job=get_setup_char('JOB',default='optimization')
    !if just estimate the wavelet or compute the gradient then this is it.
    if(job/='optimization') then
        call mpiworld_finalize
        stop
    endif
    
    call optimizer
    
    
    call mpiworld_finalize
    
end