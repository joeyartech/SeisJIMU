module m_optimizer
use m_arrayop
use m_sysio
use m_gradient
use m_parameterization
use m_preconditioner
use m_linesearcher

    private
    public init_optimizer, optimizer, current

    !Method: limited-memory BFGS
    !                                                     !
    ! x_{0}=x0                                            !
    ! x_{k+1}=x_k+\alpha_k d_k                            !
    !                                                     !
    ! where the descent direction d_k is                  !
    !                                                     !
    ! d_k=-Q_k \nabla f_k                                 !
    !                                                     !
    ! where Q_k       : l-BFGS approximation of the       !
    !                   inverse Hessian at iteration k    !
    !                   including preconditioning info    !
    !      \nabla f_k : gradient of f in x_k              !
    !                                                     !
    ! and alpha_k is the steplength computed by linesearch!
    !                                                     !
    ! See also Numerical Optimization, Nocedal, 2nd edition, 2006 !
    ! Algorithm 7.4 & 7.5 p. 178-179 !
    
    
    real f0 !initial fobjective
    real g0norm !initial gradient norm
    
    !counter
    integer :: iterate=0, max_iterate
    
    real min_update !minimum par update allowed
    
    type(t_forwardmap) :: current
    type(t_forwardmap) :: perturb
    type(t_forwardmap) :: previous
    
    !history of vector pairs
    real,dimension(:,:),allocatable :: sk,yk
    integer :: l !use l vector pairs to approx inv Hessian
    
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
        if(if_scaling) call linesearch_scaling(current)
        f0=current%f
        
        !reinitialize alpha in each iterate
        if_reinitialize_alpha=get_setup_logical('IF_REINITIALIZE_ALPHA',default=.false.)
        
        !initialize preconditioner and apply
        call init_preconditioner
        call preconditioner_apply(current%g,current%pg)
        g0norm=norm2(current%g)
        
        !current descent direction
        current%d=-current%pg
        current%gdotd=sum(current%g*current%d)
        
        !vector pairs history
        l=get_setup_int('NPAIRS',default=5); if(l<1) l=1
        call alloc(sk,n,l)
        call alloc(yk,n,l)
        sk(:,1)=current%x(:)
        yk(:,1)=current%g(:)

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
        
        integer m !declaration here can avoid conflict with m in m_model.f90
        real,dimension(:),allocatable :: q, alpha, rho, r      
        
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
            call linesearcher(iterate,current,perturb,result=result)
            
            select case(result)
                case('success')
                !update point
                current%x=perturb%x
                current%f=perturb%f
                current%g=perturb%g
                
                previous%d=current%d
                

                !!!!!!!!!!!!!!!!!!!
                !! Precond LBFGS !!
                !!!!!!!!!!!!!!!!!!!
                !update sk,yk pairs
                if(iterate<=l) then
                    sk(:,iterate)=current%x(:)-sk(:,iterate)
                    yk(:,iterate)=current%g(:)-yk(:,iterate)
                else
                    sk(:,l)=current%x(:)-sk(:,l)
                    yk(:,l)=current%g(:)-yk(:,l)
                endif
                
                m=min(iterate,l)
                call alloc(q,    n,initialize=.false.)
                call alloc(alpha,m,initialize=.false.)
                call alloc(rho,  m,initialize=.false.)

                q = current%g

                !1st part of two-loop recursion
                do i=m,1,-1
                    rho(i) = 1./dot_product(yk(:,i),sk(:,i))
                    alpha(i) = rho(i)*dot_product(sk(:,i),q(:))
                    q(:) = q(:) - alpha(i)*yk(:,i)
                enddo
                
! print*,yk
! print*,sk
! print*,rho
! print*,alpha
! print*,q
                
                !preconditioner provides seed for initial guess of inv Hessian
                call preconditioner_apply(q,q)

                gamma_num   = dot_product(sk(:,m),yk(:,m))
                gamma_denom = norm2(yk(:,m))
                gamma_denom = gamma_denom*gamma_denom

                gamma=gamma_num/gamma_denom

                r = gamma*q

                !2nd part of two-loop recursion
                do i=1,m
                    beta = rho(i)*dot_product(yk(:,i),r(:))
                    r(:) = r(:) + (alpha(i)-beta)*sk(:,i)
                enddo
                current%d = -r
                
! print*,gamma,gamma_num,gamma_denom
! print*,beta
! print*,current%d
! pause
                
                !save sk,yk pairs
                if(iterate+1.le.l) then !fill in history
                    sk(:,iterate+1)=current%x(:)
                    yk(:,iterate+1)=current%g(:)
                else !shift and replace oldest pairs
                    sk(:,1:l-1)=sk(:,2:l);  sk(:,l)=current%x(:)
                    yk(:,1:l-1)=yk(:,2:l);  yk(:,l)=current%g(:)
                endif
                
                !inner product of g and d for Wolfe condition
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
                    write(16,'(a)'      ) '                     LIMITED MEMORY BFGS ALGORITHM                     '
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
                    open(18,file='model_update',access='stream')
                case('update')
                    write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, current%f, current%f/f0, norm2(current%g)/g0norm, alpha, isearch, imodeling
                    
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
