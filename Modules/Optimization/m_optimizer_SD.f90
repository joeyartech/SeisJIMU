module m_optimizer
use m_optimizer_common

    private
    public :: optimizer_loop

    !Method: Steepest Descent
    !                                                     !
    ! x_{0}=x0                                            !
    ! x_{k+1}=x_k+\alpha_k d_k                            !
    !                                                     !
    ! where the descent direction d_k is                  !
    !                                                     !
    ! d_k=-Q_k \nabla f_k                                 !
    !                                                     !
    ! where Q_k       : preconditioner at iteration k     !
    !      \nabla f_k : gradient of f in x_k              !
    !                                                     !
    ! and alpha_k is the steplength computed by linesearch!
    !                                                     !
    
    character(*),parameter :: info='Optimization: Steepest Descent'

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

        call optimizer_write('start','STEEPEST DESCENT ALGORITHM')

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

                !new descent direction
                curr%d=prev_d
                
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