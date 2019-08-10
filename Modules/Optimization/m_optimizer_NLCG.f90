module m_optimizer
use m_arrayop
use m_sysio
use m_gradient
use m_parameterization
use m_preconditioner
use m_linesearcher

    private
    public init_optimizer, optimizer, current

    !Method: Nonlinear Conjugate Gradient
    !                                                     !
    ! x_{0}=x0                                            !
    ! x_{k+1}=x_k+\alpha_k d_k                            !
    !                                                     !
    ! where the descent direction d_k is                  !
    !                                                     !
    ! d_k=-Q_k \nabla f_k + \beta_k d_{k-1}               !
    !                                                     !
    ! where Q_k       : preconditioner at iteration k     !
    !      \nabla f_k : gradient of f in x_k              !
    !      \beta_k    : scalar computed through the       !
    !                   Dai and Yuan formula              !
    !                                                     !
    ! and alpha_k is the steplength computed through the  !
    ! common linesearch algorithm of the TOOLBOX          !
    !                                                     !
    ! The first call to the algorithm must be done with   !
    ! FLAG='INIT'. For this first call, the initial point !
    ! x0 is given through the variable x. The input       !
    ! variables fcost and grad, grad_preco must correspond!
    ! respectively to the misfit, gradient and            !
    ! preconditioned gradient at x0.                      !
    ! Computation of the descent direction given by the 
    ! preconditioned nonlinear conjugate gradient algorithm 
    ! of Dai and Yuan 
    ! Y. DAI AND Y. YUAN, A nonlinear conjugate gradient  !
    ! method with a strong global convergence property,   !
    ! SIAM Journal on Optimization, 10 (1999), pp. 177-182!
    !                                                     !
    ! See also Nocedal, Numerical optimization,           !
    ! 2nd edition p.132                                   !
    
    
    real f0 !initial fobjective
    real g0norm !initial gradient norm
    
    integer,parameter :: l=1 !save one previous gradient in history
    
    !counter
    integer :: iterate=0, max_iterate
    
    real min_update !minimum par update allowed
    
    type(t_forwardmap) :: current
    type(t_forwardmap) :: perturb
    type(t_forwardmap) :: previous
        
    real,dimension(:,:),allocatable :: gradient_history
    
    logical :: debug
    
    contains
    
    subroutine init_optimizer(nproblem)
        
        !problem size
        n=nproblem
        
        !read setup
        min_update=get_setup_real('MIN_UPDATE',default=1e-8)
        max_iterate=get_setup_int('MAX_ITERATE',default=30)
        max_modeling=get_setup_int('MAX_MODELING',default=60)
        
        !initialize current point
        call alloc(current%x, n,initialize=.false.) !quiry point
        current%f=fobjective
        call alloc(current%g, n,initialize=.false.) !gradient
        call alloc(current%pg,n,initialize=.false.) !preconditioned gradient
        call alloc(current%d, n,initialize=.false.) !descent direction
        
        !transform by parameterization
        call parameterization_transform('m2x',current%x,current%g)
        
        !scale problem
        call linesearch_scaling(current)
        f0=current%f
        
        !initialize preconditioner and apply
        call init_preconditioner
        call preconditioner_apply(current%g,current%pg)
        g0norm=sqrt(sum(current%g*current%g))
        
        !current descent direction
        current%d=-current%pg
        current%gdotd=sum(current%g*current%d)
        
        !gradient history
        call alloc(gradient_history,n,l,initialize=.false.)
        gradient_history(:,1)=current%g
        
        !initialize perturb point
        call alloc(perturb%x, n,initialize=.false.)
        call alloc(perturb%g, n,initialize=.false.)
        call alloc(perturb%pg,n,initialize=.false.)
        call alloc(perturb%d, n,initialize=.false.)
        
        !initialize previous point
        call alloc(previous%d,n,initialize=.false.)
        previous%d=current%d
        
    end subroutine
    
    subroutine optimizer
        character(7) :: result
        
        real nom, denom
        
        
        call optimizer_print_info('start')
        call hud('============ START OPTIMIZATON ============')
        
        !iteration loop
        loop: do iterate=1,max_iterate
        
            if (sqrt(sum(current%d*current%d)) < min_update) then
                call hud('Maximum iteration number reached or convergence criteria met')
                call optimizer_print_info('criteria')
                exit loop
            endif
            
            !linesearcher finds steplength
            call linesearcher(iterate,current,perturb,l,gradient_history,result)
            
            select case(result)
                case('success')
                !update point
                current%x=perturb%x
                current%f=perturb%f
                current%g=perturb%g
                
                previous%d=current%d
                
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!
                !! NLCG, Dai-Yuan's beta !!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!
                !compute beta
                nom=sum(perturb%g*perturb%pg)
                denom=sum((perturb%g-gradient_history(:,1))*previous%d)
                beta=nom/denom
                
                !Safeguard 
                if((beta.ge.1e5).or.(beta.le.-1e5)) then     
                    beta=0.
                endif
                
                !new descent direction
                current%d=-perturb%pg+beta*previous%d
                
                current%gdotd=sum(current%g*current%d)
                
                call optimizer_print_info('update')
                
                case('failure')
                call optimizer_print_info('failure')
                exit loop
                
                case('maximum')
                call optimizer_print_info('maximum')
                exit loop
                
            end select
            
        end do loop
        
        call hud('-------- FINALIZE OPTIMIZATON --------')
        call optimizer_print_info('finalize')
        
    end subroutine

    subroutine optimizer_print_info(task)
        character(*) :: task
        
        if(mpiworld%is_master) then
        
            select case (task)
                case('start')
                    open(16,file='iterate.log')
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '                 NONLINEAR CONJUGATE GRADIENT ALGORITHM                '
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a,es8.2)') '     Min update allowed        : ',  min_update
                    write(16,'(a,i5)'   ) '     Max iteration allowed     : ',  max_iterate
                    write(16,'(a,i5)'   ) '     Max linesearch allowed    : ',  max_search
                    write(16,'(a,i5)'   ) '     Max modeling allowed      : ',  max_modeling
                    write(16,'(a,es8.2)') '     Initial fobjective (f0)         : ',  f0
                    write(16,'(a,es8.2)') '     Initial gradient norm (||g0||)  : ',  g0norm
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '  Iter#      f         f/f0    ||g||/||g0||    alpha     nls  Modeling#'
                                   !e.g.  !    0    1.00E+00    1.00E+00    1.00E+00    1.00E+00      0       1
                    write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, current%f, current%f/f0, sqrt(sum(current%g*current%g))/g0norm, alpha, isearch, imodeling
                    open(18,file='model_update',access='stream')
                case('update')
                    write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, current%f, current%f/f0, sqrt(sum(current%g*current%g))/g0norm, alpha, isearch, imodeling
                    
                    call parameterization_transform('x2m',current%x)
                    if(if_par_vp) write(18) m%vp
                    if(if_par_rho) write(18) m%rho
                    if(if_par_ip) write(18) m%vp*m%rho
                case('maximum')
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '     STOP: MAXIMUM MODELING NUMBER REACHED                             '
                    write(16,'(a)'      ) ' **********************************************************************'
                case('criteria')
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '     STOP: CONVERGENCE CRITERIA SATISFIED                              '
                    write(16,'(a)'      ) ' **********************************************************************'
                case('failure')
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '     STOP: LINE SEARCH FAILURE                                         '
                    write(16,'(a)'      ) ' **********************************************************************'
                case('finalize')
                    close(16)
                    
                    call parameterization_transform('x2m',perturb%x)
                    if(if_par_vp) write(18) m%vp
                    if(if_par_rho) write(18) m%rho
                    if(if_par_ip) write(18) m%vp*m%rho
                    close(18)
                    
                    open(18,file='model_final',access='stream')
                    call parameterization_transform('x2m',current%x)
                    if(if_par_vp) write(18) m%vp
                    if(if_par_rho) write(18) m%rho
                    if(if_par_ip) write(18) m%vp*m%rho
                    close(18)
                    
                    write(*,'(a,i0.4)') 'ximage < model_final n1=',m%nz
                    write(*,'(a,i0.4,a,i0.4)') 'xmovie < model_update loop=1 title=%g n1=',m%nz,' n2=',m%nx
                end select
        endif
    end subroutine
end
