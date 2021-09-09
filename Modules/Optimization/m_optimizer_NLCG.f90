module m_optimizer
use m_optimizer_common

    private
    public :: optimizer_loop

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

    character(*),parameter :: info='Optimization: NonLinear Conjugate Gradient'//s_NL// &
                                   'with Dai-Yuan formula'

    !history of computed gradients
    real,dimension(:,:,:,:,:),allocatable :: gradient_history
    integer,parameter :: l=1 !save one previous gradient in history
        
    contains
    
    subroutine optimizer_loop
        type(t_checkpoint),save :: chp
        
        real nom, denom

        !gradient history
        call alloc(gradient_history,param%n1,param%n2,param%n3,param%npars,l)
        gradient_history(:,:,:,:,l)=curr%g


        call hud(info)

        
        call hud('============ START OPTIMIZATON ============')
        call chp%init('FWI_shotlist_optimizer','Iterate#',oif_fuse=.true.)

        call optimizer_write('start','NONLINEAR CONJUGATE GRADIENT ALGORITHM')

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
            call ls%search(.false.,iterate,curr,pert,gradient_history)
            
            select case(ls%result)
                case('success')

                !previous descent direction
                prev_d=curr%d
                
                call switch_curr_pert

                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!
                !! NLCG, Dai-Yuan's beta !!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!
                !compute beta
                nom=sum(curr%g*curr%pg)
                denom=sum((curr%g-gradient_history(:,:,:,:,l))*prev_d)
                beta=nom/denom
                
                !Safeguard 
                if((beta.ge.1e5).or.(beta.le.-1e5)) beta=0.
                
                !new descent direction
                curr%d=-curr%pg+beta*prev_d
                
                !inner product of g and d for Wolfe condition
                curr%g_dot_d = sum(curr%g*curr%d)

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