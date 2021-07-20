module m_optimizer
use m_string
use m_mpienv
use m_arrayop
use m_setup
use m_checkpoint
use m_sysio
use m_model
use m_shotlist
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
    
    real min_update !minimum parameter updates allowed
 
    type,public :: t_optimizer
        integer :: n !dimension of model manifold
        type(t_forwardmap) :: curr !current point
        type(t_forwardmap) :: pert !point perturbed from current point
        type(t_forwardmap) :: prev !previous point

        real f0 !initial fobjective
        real g0norm !initial gradient norm
 
        !counter
        integer :: iterate=0 !number of iteration performed 
        integer :: max_iterate  !max number of iteration allowed

        contains
        procedure :: init => init
        procedure :: optimize => optimize
        procedure :: write => write
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
        max_iterate=setup%get_int('MAX_ITERATE',o_default='30')
        call ls%init

        associate(curr=>self%curr, pert=>self%pert, prev=>self%prev)
        
            !initialize current point
            call curr%init(self%n)

            self%f0=curr%f
            
            !scale problem
            call ls%scale(curr)
            
            !initialize preconditioner and apply
            call preco%init
            call preco%apply(curr%g,curr%pg)
            self%g0norm=norm2(curr%g)
            
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

        type(t_checkpoint),save :: chp
        real nom, denom
       
        ! call self%write('start')
        call hud('============ START OPTIMIZATON ============')
        call chp%init('FWI_shotlist_optim','Iter#',oif_fuse=.true.)
        
        associate(curr=>self%curr,pert=>self%pert,prev=>self%prev)

            !iteration loop
            loop: do iterate=1,self%max_iterate
                call chp%count

                if(.not.shls%is_registered(chp,'sampled_shots')) then
                    call shls%sample
                    call shls%register(chp,'sampled_shots')
                endif
                call shls%write
                call shls%assign

                self%iterate=iterate
            
                if (norm2(curr%d) < min_update) then
                    call hud('Maximum iteration number reached or convergence criteria met')
                    call self%write('criteria')
                    exit loop
                endif
                
                !linesearcher finds steplength
                call ls%search(self%iterate,curr,pert,gradient_history)
                
                select case(ls%result)
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
                    
                    call self%write('update')

                    case('failure')
                    call self%write('failure')
                    exit loop
                    
                    case('maximum')
                    call self%write('maximum')
                    exit loop
                    
                end select
                
            end do loop

        end associate
        
        call hud('-------- FINALIZE OPTIMIZATON --------')
        call self%write('finalize')
        
    end subroutine

    subroutine write(self,message)
        class(t_optimizer) :: self
        character(*) :: message
        
        if(mpiworld%is_master) then
        
            select case (message)
                case('start')
                    open(16,file=dir_out//'iterate.log',position='append',action='write')
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '                 NONLINEAR CONJUGATE GRADIENT ALGORITHM                '
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a,es8.2)') '     Min update allowed          : ',  min_update
                    write(16,'(a,i5)'   ) '     Max iterates allowed        : ',  self%max_iterate
                    write(16,'(a,i5)'   ) '     Max linesearches allowed    : ',  ls%max_search
                    write(16,'(a,i5)'   ) '     Max gradients allowed       : ',  ls%max_gradient
                    write(16,'(a,es8.2)') '     Initial fobjective (f0)         : ',  self%f0
                    write(16,'(a,es8.2)') '     Initial gradient norm (||g0||)  : ',  self%g0norm
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '  Iter#      f         f/f0    ||g||/||g0||    alpha     nls  Modeling#'
                                   !e.g.  !    0    1.00E+00    1.00E+00    1.00E+00    1.00E+00      0       1
                    write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  self%iterate, self%curr%f, self%curr%f/self%f0, norm2(self%curr%g)/self%g0norm, ls%alpha, ls%isearch, ls%igradient
                    close(16)
                    
                case('update')
                    open(16,file=dir_out//'iterate.log',position='append',action='write')
                    write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  self%iterate, self%curr%f, self%curr%f/f0, norm2(self%curr%g)/g0norm, ls%alpha, ls%isearch, ls%igradient
                    close(16)
                    
                    call param%transform_model('x->m',self%curr%x)
                    call m%write(o_suffix='Iter'//int2str(self%iterate))

                case('maximum')
                    open(16,file=dir_out//'iterate.log',position='append',action='write')
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '     STOP: MAXIMUM MODELING NUMBER REACHED                             '
                    write(16,'(a)'      ) ' **********************************************************************'
                    close(16)

                case('criteria')
                    open(16,file=dir_out//'iterate.log',position='append',action='write')
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '     STOP: CONVERGENCE CRITERIA SATISFIED                              '
                    write(16,'(a)'      ) ' **********************************************************************'
                    close(16)

                case('failure')
                    open(16,file=dir_out//'iterate.log',position='append',action='write')
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '     STOP: LINE SEARCH FAILURE                                         '
                    write(16,'(a)'      ) ' **********************************************************************'
                    close(16)

                case('finalize')
                    open(16,file=dir_out//'iterate.log',position='append',action='write')
                    write(16,'(a)'      ) ' **********************************************************************'
                    write(16,'(a)'      ) '     FINALIZE                                         '
                    write(16,'(a)'      ) ' **********************************************************************'
                    close(16)
                    
                    call param%transform_model('x->m',self%pert%x)
                    call m%write(o_suffix='Iter'//int2str(self%iterate))
                    
                    call param%transform_model('x->m',self%curr%x)
                    call m%write(o_suffix='final')
                    
                    write(*,'(a,i0.4)') 'ximage < model_final n1=',m%nz
                    write(*,'(a,i0.4,a,i0.4)') 'xmovie < model_update loop=1 title=%g n1=',m%nz,' n2=',m%nx

                end select

        endif
    
    end subroutine

end