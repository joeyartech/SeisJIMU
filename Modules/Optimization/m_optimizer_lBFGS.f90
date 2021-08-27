module m_optimizer
use m_System
use m_Modeling
use m_Kernel
use m_linesearcher

    private
    public :: optimizer_init, optimizer_loop, optimizer_write

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
    
    character(*),parameter :: info='Optimization: limited-memory BFGS'//s_NL// &
                                   'with preconditioned gradients'
    
    type(t_querypoint),pointer :: curr, pert
    type(t_querypoint),target  :: qp1
    real,dimension(:,:,:,:),allocatable :: prev_d  !descent direction of previous point (no need for another t_querypoint instance)

    real f0 !initial fobjective
    real g0norm !initial gradient norm
    
    !counter
    integer :: iterate=0 !number of iteration performed 
    integer :: max_iterate  !max number of iteration allowed

    real min_descent !minimum descent allowed

    !history of vector pairs
    real,dimension(:,:),allocatable :: sk,yk
    integer :: l !use l vector pairs to approx inv Hessian
    
    contains
    
    subroutine optimizer_init(qp0)        
        type(t_querypoint),target :: qp0 !initial (model) parameters

        call hud(info)

        !current point
        curr=>qp0
        curr%d=-curr%pg
        curr%g_dot_d = sum(curr%g*curr%d)

        !initial values
        f0=curr%f
        g0norm=norm2(curr%g)
        
        !perturbed point
        pert=>qp1
        call pert%init        

        !descent direction of previous point
        prev_d=curr%d

        !vector pairs history
        l=setup%get_int('NPAIRS',o_default='5'); if(l<1) l=1
        call alloc(sk,param%n,l)
        call alloc(yk,param%n,l)
        sk(:,1)=pack(curr%x,.true.)
        yk(:,1)=pack(curr%g,.true.)

        !read setup
        min_descent=setup%get_real('MIN_DESCENT',o_default='1e-8')
        max_iterate=setup%get_int('MAX_ITERATE',o_default='30')
        ls%max_gradient=setup%get_int('MAX_GRADIENT',o_default=num2str(max_iterate+30))        
        
    end subroutine
    
    subroutine optimizer_loop
        type(t_checkpoint),save :: chp
        type(t_querypoint),pointer :: tmp
        ! real nom, denom

        integer m !declaration here can avoid conflict with m in m_model.f90
        real,dimension(:),allocatable :: q, alpha, rho, r             

        call hud('============ START OPTIMIZATON ============')
        call chp%init('FWI_shotlist_optim','Iter#',oif_fuse=.true.)

        call optimizer_write('start')

        !iteration loop
        loop: do iterate=1,max_iterate
            call chp%count

            if(.not.shls%is_registered(chp,'sampled_shots')) then
                call shls%sample
                call shls%register(chp,'sampled_shots')
            endif
            call shls%assign

            if (maxval(abs(curr%d)) < min_descent) then
                call hud('Minimum descent is met')
                call optimizer_write('criteria')
                exit loop
            endif

            !reinitialize linesearch
            ! ls%alphaL=0.
            ! ls%alphaR=huge(1.)
            !ls%alpha=ls%alpha0 !reinitialize alpha in each iterate may help convergence for LBFGS method..
            !linesearcher finds steplength
            call ls%search(.true.,iterate,curr,pert)
            
            select case(ls%result)
                case('success')

                !previous descent direction
                prev_d=curr%d

                !switch curr & pert
                tmp=>curr
                curr=>pert
                pert=>tmp

                !!!!!!!!!!!!!!!!!!!
                !! Precond LBFGS !!
                !!!!!!!!!!!!!!!!!!!
                !update sk,yk pairs
                m=min(iterate,l)
                sk(:,m)=reshape(curr%x,[param%n])-sk(:,m)
                yk(:,m)=reshape(curr%g,[param%n])-yk(:,m)

                call alloc(q,    param%n)
                call alloc(alpha,m)
                call alloc(rho,  m)
                q = reshape(curr%g,[param%n])

                !1st part of two-loop recursion
                do i=m,1,-1
                    rho(i) = 1./dot_product(yk(:,i),sk(:,i))
                    alpha(i) = rho(i)*dot_product(sk(:,i),q)
                    q = q - alpha(i)*yk(:,i)
                enddo

! print*,yk
! print*,sk
! print*,rho
! print*,alpha
! print*,q
                
                !preconditioner provides seed for initial guess of inv Hessian
                call preco%apply(q,q)
                gamma_num   = dot_product(sk(:,m),yk(:,m))
                gamma_denom = norm2(yk(:,m))
                gamma_denom = gamma_denom*gamma_denom
                gamma=gamma_num/gamma_denom
                r = gamma*q

                !2nd part of two-loop recursion
                do i=1,m
                    beta = rho(i)*dot_product(yk(:,i),r(:))
                    r    = r + (alpha(i)-beta)*sk(:,i)
                enddo

                curr%d = reshape(-r,[param%n1,param%n2,param%n3,param%npars])

! print*,gamma,gamma_num,gamma_denom
! print*,beta
! print*,current%d
! pause
                
                !save sk,yk pairs
                if(iterate+1<=l) then !fill in history
                    sk(:,iterate+1)=reshape(curr%x,[param%n])
                    yk(:,iterate+1)=reshape(curr%g,[param%n])
                else !shift and replace oldest pairs
                    ! sk=eoshift(sk,dim=2,shift=1,boundary=reshape(curr%x,[param%n]))
                    ! yk=eoshift(sk,dim=2,shift=1,boundary=reshape(curr%g,[param%n]))
                    !use eoshift may require "ulimit -s unlimited"
                    !so use vectorized expressions instead
                    sk(:,1:l-1)=sk(:,2:l);  sk(:,l)=reshape(curr%x,[param%n])
                    yk(:,1:l-1)=yk(:,2:l);  yk(:,l)=reshape(curr%g,[param%n])
                endif

                !inner product of g and d for Wolfe condition
                curr%g_dot_d=sum(curr%g*curr%d)

                call optimizer_write('update')

                case('failure')
                call optimizer_write('failure')
                exit loop
                
                case('maximum')
                call optimizer_write('maximum')
                exit loop
                
            end select
            
        end do loop
        
        call hud('-------- FINALIZE OPTIMIZATON --------')
        call optimizer_write('finalize')
        
    end subroutine

    subroutine optimizer_write(message)
        character(*) :: message
        
        if(mpiworld%is_master) then
        
            select case (message)
            case('start')
                ! call execute_command_line('rm '//dir_out//'iterate.log', wait=.true.)
                open(16,file=dir_out//'iterate.log',position='append',action='write')
                write(16,'(a)'      ) ' **********************************************************************'
                write(16,'(a)'      ) '                     LIMITED MEMORY BFGS ALGORITHM                     '
                write(16,'(a)'      ) ' **********************************************************************'
                write(16,'(a,es8.2)') '     Min descent allowed      =',  min_descent
                write(16,'(a,i5)'   ) '     Max iterates allowed     =',  max_iterate
                write(16,'(a,i5)'   ) '     Max linesearches allowed =',  ls%max_search
                write(16,'(a,i5)'   ) '     Max gradients allowed    =',  ls%max_gradient
                write(16,'(a,es8.2)') '     Initial fobjective (f0)        =',  f0
                write(16,'(a,es8.2)') '     Initial gradient norm (||g0||) =',  g0norm
                write(16,'(a)'      ) ' **********************************************************************'
                write(16,'(a)'      ) '  Iter#      f         f/f0    ||g||/||g0||    alpha     nls  Gradient#'
                               !e.g.  !    0    1.00E+00    1.00E+00    1.00E+00    1.00E+00      0       1
                write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, curr%f, curr%f/f0, norm2(curr%g)/g0norm, ls%alpha, ls%isearch, ls%igradient
                close(16)
                
            case('update')
                open(16,file=dir_out//'iterate.log',position='append',action='write')
                write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, curr%f, curr%f/f0, norm2(curr%g)/g0norm, ls%alpha, ls%isearch, ls%igradient
                close(16)
                
                call param%transform('x->m',o_x=curr%x)
                call m%write(o_suffix='_Iter'//int2str(iterate))
                call sysio_write('pg_Iter'//int2str(iterate),curr%pg,size(curr%pg))
                call shot%write('dsyn_Iter'//int2str(iterate)//'_',shot%dsyn)

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
                
                call param%transform('x->m',o_x=pert%x)
                call m%write(o_suffix='_Iter'//int2str(iterate))
                call sysio_write('pg_Iter'//int2str(iterate),curr%pg,size(curr%pg))
                call shot%write('dsyn_Iter'//int2str(iterate)//'_',shot%dsyn)

                write(*,'(a,i0.4)') 'ximage < model_final n1=',m%nz
                write(*,'(a,i0.4,a,i0.4)') 'xmovie < model_update loop=1 title=%g n1=',m%nz,' n2=',m%nx

            end select

        endif
    
    end subroutine

end
