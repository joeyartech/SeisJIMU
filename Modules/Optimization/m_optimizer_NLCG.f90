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
    ! and alpha_k is the steplength computed by linesearch!
    !                                                     !
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
    
    !counter
    integer :: iterate=0, max_iterate
    
    real min_update !minimum par update allowed
    
    type t_forwardmap
        integer :: n !model dimension
        real,dimension(:),allocatable :: x !point
        real,dimension(:),allocatable :: g !gradient
        real,dimension(:),allocatable :: pg !preconditioned gradient
        real,dimension(:),allocatable :: d !descent direction
        real :: f
        real :: gdotd  !dot product of g and d
    end type

    type(t_forwardmap) :: current
    type(t_forwardmap) :: perturb
    type(t_forwardmap) :: previous
    
    !history of computed gradients
    real,dimension(:,:),allocatable :: gradient_history
    integer,parameter :: l=1 !save one previous gradient in history
    
    logical :: debug
    
    contains
    
    subroutine init(nproblem)
    
        !problem size
        n=nproblem
        
        !read setup
        min_update=setup_get_real('MIN_UPDATE',default=1e-8)
        max_iterate=setup_get_int('MAX_ITERATE',default=30)
        
        !initialize current point
        call alloc(current%x, n,initialize=.false.) !quiry point
        current%f=fobjective
        call alloc(current%g, n,initialize=.false.) !gradient
        call alloc(current%pg,n,initialize=.false.) !preconditioned gradient
        call alloc(current%d, n,initialize=.false.) !descent direction
        
        !transform by parameterization
        call parameterization%transform('m2x',current%x,current%g)
        
        !scale problem
        call ls%scaling(current)
        f0=current%f
        
        !initialize preconditioner and apply
        call preconditioner%init
        call preconditioner%apply(current%g,current%pg)
        g0norm=norm2(current%g)
        
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
        
            if (norm2(current%d) < min_update) then
                call hud('Maximum iteration number reached or convergence criteria met')
                call optimizer_print_info('criteria')
                exit loop
            endif
            
            !linesearcher finds steplength
            call ls%search(iterate,current,perturb,gradient_history,result=result)
            
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
                
                !inner product of g and d for Wolfe condition
                current%gdotd=sum(current%g*current%d)
                
                call print_info('update')
                
                case('failure')
                call print_info('failure')
                exit loop
                
                case('maximum')
                call print_info('maximum')
                exit loop
                
            end select
            
        end do loop
        
        call hud('-------- FINALIZE OPTIMIZATON --------')
        call print_info('finalize')
        
    end subroutine

    subroutine print_info(task)
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
                    write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, current%f, current%f/f0, norm2(current%g)/g0norm, alpha, isearch, imodeling
                    
                case('update')
                    write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, current%f, current%f/f0, norm2(current%g)/g0norm, alpha, isearch, imodeling
                    
                    call parameterization%transform('x2m',current%x)
                    call model%write(suffix='Iter'//int2str(iterate))
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
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '     STOP: LINE SEARCH FINISHED                                         '
                    write(16,'(a)'      ) ' **********************************************************************'
                    close(16)
                    
                    call parameterization%transform('x2m',perturb%x)
                    call model%write(suffix='Iter'//int2str(iterate))
                    
                    call parameterization%transform('x2m',current%x)
                    call model%write(suffix='final')
                    
                    write(*,'(a,i0.4)') 'ximage < model_final n1=',m%nz
                    write(*,'(a,i0.4,a,i0.4)') 'xmovie < model_update loop=1 title=%g n1=',m%nz,' n2=',m%nx
                end select
        endif
    end subroutine
end
