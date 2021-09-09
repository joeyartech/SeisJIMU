module m_optimizer
use m_optimizer_common

    private
    public :: optimizer_loop

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

    !history of vector pairs
    real,dimension(:,:),allocatable :: sk,yk
    integer :: l !use l vector pairs to approx inv Hessian
    
    contains
    
    subroutine optimizer_loop
        type(t_checkpoint),save :: chp
        type(t_querypoint),pointer :: tmp
        ! real nom, denom

        integer m !declaration here can avoid conflict with m in m_model.f90
        real,dimension(:),allocatable :: q, alpha, rho, r


        !vector pairs history
        l=setup%get_int('NPAIRS',o_default='5'); if(l<1) l=1
        call alloc(sk,param%n,l)
        call alloc(yk,param%n,l)
        sk(:,1)=reshape(curr%x,[param%n])
        yk(:,1)=reshape(curr%g,[param%n])


        call hud(info)


        call hud('============ START OPTIMIZATON ============')
        call chp%init('FWI_shotlist_optimizer','Iter#',oif_fuse=.true.)

        call optimizer_write('start','LIMITED MEMORY BFGS ALGORITHM')

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

                call switch_curr_pert


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
                call preco%apply_ext(q,q)

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

end