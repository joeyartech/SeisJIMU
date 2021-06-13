module m_optimizer
use m_global
use m_string
use m_arrayop
use m_setup
use m_mpienv
use m_model
use m_propagator
use m_fobjective
use m_nabla
use m_parameterization
use m_preconditioner
use m_linesearcher

    private

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
    
    real min_update !minimum parameter updates allowed
 
    type,public :: t_optimizer
        integer :: n !dimension of model manifold
        type(t_forwardmap) :: curr !current point
        type(t_forwardmap) :: pert !point perturbed from current point
        type(t_forwardmap) :: prev !previous point

        contains
        procedure :: init => init
        procedure :: optimize => optimize
        ! procedure :: write => write
        ! final :: fin
    end type
    
    !history of computed gradients
    real,dimension(:,:),allocatable :: gradient_history
    integer,parameter :: l=1 !save one previous gradient in history
    
    logical :: debug
    
    contains

    subroutine init(self,n)
        class(t_optimizer) :: self

        self%n=n
    
        !read setup
        min_update=setup%get_real('MIN_UPDATE',o_default='1e-8')
        call ls%init

        associate(curr=>self%curr, pert=>self%pert, prev=>self%prev)
        
        !initialize current point
        call curr%init(self%n,o_fobj=fobj%total_misfit)
        
            f0=curr%f

            !transform by parameterization
            call param%transform('m2x',curr%x,curr%g)
            
            !scale problem
            call ls%scale(curr)
            
            !initialize preconditioner and apply
            call preco%init
            call preco%apply(curr%g,curr%pg)
            g0norm=norm2(curr%g)
            
            !current descent direction
            curr%d=-curr%pg
            curr%gdotd=sum(curr%g*curr%d)
            
            !gradient history
            call alloc(gradient_history,n,l)
            gradient_history(:,1)=curr%g
            
            !initialize perturb point
            call pert%init(self%n)
            prev%d=curr%d

        end associate
        
    end subroutine
    
    subroutine optimize(self)
        class(t_optimizer) :: self

        character(7) :: result
        real nom, denom
        
        ! call self%write('start')
        call hud('============ START OPTIMIZATON ============')
        
        associate(curr=>self%curr,pert=>self%pert,prev=>self%prev)

            !iteration loop
            loop: do iterate=1,max_iterate
            
                if (norm2(curr%d) < min_update) then
                    call hud('Maximum iteration number reached or convergence criteria met')
                    ! call self%write('criteria')
                    exit loop
                endif
                
                !linesearcher finds steplength
                call ls%search(iterate,curr,pert,gradient_history,result=result)
                
                select case(result)
                    case('success')
                    !update point
                    curr%x=pert%x
                    curr%f=pert%f
                    curr%g=pert%g
                    
                    prev%d=curr%d
                    
                    
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !! NLCG, Dai-Yuan's beta !!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !compute beta
                    nom=sum(pert%g*pert%pg)
                    denom=sum((pert%g-gradient_history(:,1))*prev%d)
                    beta=nom/denom
                    
                    !Safeguard 
                    if((beta.ge.1e5).or.(beta.le.-1e5)) then     
                        beta=0.
                    endif
                    
                    !new descent direction
                    curr%d=-pert%pg+beta*prev%d
                    
                    !inner product of g and d for Wolfe condition
                    curr%gdotd=sum(curr%g*curr%d)
                    
                    ! call self%write('update')
                    
                    case('failure')
                    ! call self%write('failure')
                    exit loop
                    
                    case('maximum')
                    ! call self%write('maximum')
                    exit loop
                    
                end select
                
            end do loop

        end associate
        
        call hud('-------- FINALIZE OPTIMIZATON --------')
        ! call self%write('finalize')
        
    end subroutine

    ! subroutine write(self,task)
    !     class(t_optimizer) :: self
    !     character(*) :: task
        
    !     if(mpiworld%is_master) then
        
    !         select case (task)
    !             case('start')
    !                 open(16,file='iterate.log')
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !                 write(16,'(a)'      ) '                 NONLINEAR CONJUGATE GRADIENT ALGORITHM                '
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !                 write(16,'(a,es8.2)') '     Min update allowed        : ',  min_update
    !                 write(16,'(a,i5)'   ) '     Max iteration allowed     : ',  max_iterate
    !                 write(16,'(a,i5)'   ) '     Max linesearch allowed    : ',  max_search
    !                 write(16,'(a,i5)'   ) '     Max modeling allowed      : ',  max_modeling
    !                 write(16,'(a,es8.2)') '     Initial fobjective (f0)         : ',  f0
    !                 write(16,'(a,es8.2)') '     Initial gradient norm (||g0||)  : ',  g0norm
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !                 write(16,'(a)'      ) '  Iter#      f         f/f0    ||g||/||g0||    alpha     nls  Modeling#'
    !                                !e.g.  !    0    1.00E+00    1.00E+00    1.00E+00    1.00E+00      0       1
    !                 write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, current%f, current%f/f0, norm2(current%g)/g0norm, alpha, isearch, imodeling
                    
    !             case('update')
    !                 write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, current%f, current%f/f0, norm2(current%g)/g0norm, alpha, isearch, imodeling
                    
    !                 call parameterization%transform('x2m',current%x)
    !                 call model%write(suffix='Iter'//int2str(iterate))
    !             case('maximum')
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !                 write(16,'(a)'      ) '     STOP: MAXIMUM MODELING NUMBER REACHED                             '
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !             case('criteria')
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !                 write(16,'(a)'      ) '     STOP: CONVERGENCE CRITERIA SATISFIED                              '
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !             case('failure')
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !                 write(16,'(a)'      ) '     STOP: LINE SEARCH FAILURE                                         '
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !             case('finalize')
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !                 write(16,'(a)'      ) '     STOP: LINE SEARCH FINISHED                                         '
    !                 write(16,'(a)'      ) ' **********************************************************************'
    !                 close(16)
                    
    !                 call parameterization%transform('x2m',perturb%x)
    !                 call model%write(suffix='Iter'//int2str(iterate))
                    
    !                 call parameterization%transform('x2m',current%x)
    !                 call model%write(suffix='final')
                    
    !                 write(*,'(a,i0.4)') 'ximage < model_final n1=',m%nz
    !                 write(*,'(a,i0.4,a,i0.4)') 'xmovie < model_update loop=1 title=%g n1=',m%nz,' n2=',m%nx
    !             end select
    !     endif
    
    ! end subroutine

end