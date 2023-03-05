!Find proper steplength alpha (α) based
!first on a bracketing strategy,
!then on a dichotomy algorithm. 
!Ref: Numerical Optimizationn Theoretical and Practical Aspects, J.F.Bonnans, J.C.Gilbert, C. Lemaréchal, C.A. Sagastizábal,  Springer-Verlag, Universitext

!This alpha should satisfy Wolfe's conditions 
!1st_cond: sufficient descent
!2nd_cond: sufficient curvature
!Ref: Nocedal, Numerical Optimization, 2nd Ed, P. 33

module m_linesearcher
use m_System
use m_Modeling
use m_Kernel

    private
    
    character(*),parameter :: info='Linesearch based on a bracketing - dichotomy algoritm'//s_NL// &
                                   'Steplength (alpha) judged by Wolfe conditions'

    !Wolfe conditions parameters (Nocedal value)
    real,parameter :: c1=1e-4, c2=0.9 !c2=0.9 for (quasi-)Newton method, 0.1 for NLCG
    !Bracketting parameter (Gilbert value)
    real,parameter :: multiplier=10.
    
    !thresholding
    real,parameter :: thres=0.

    !priors for steplength
    real,parameter :: alpha0=1.
    real,parameter :: alphaL0=1e-4 !~=2^-10
    real,parameter :: alphaR0=1e4  !=10^4

    !Wolfe conditions
    logical :: if_1st_cond, if_2nd_cond

    type,public :: t_linesearcher
        real :: alpha  !steplength
        real :: alphaL, alphaR !search interval: alpha \in [alphaL, alphaR]
        real :: scaler
        character(7) :: result
        
        !counter
        integer :: isearch !number of linesearch performed in each iterate
        integer :: max_search !max number of linesearch allowed per iteration
        integer :: igradient=1 !total number of gradient computed
        integer :: max_gradient !max total number of gradient computation allowed
    
        contains
        procedure :: init
        procedure :: search
        procedure :: scale
        procedure :: write
    end type

    type(t_linesearcher),public :: ls

    
    contains
    
    subroutine init(self)
        class(t_linesearcher) :: self

        call hud(info)
        call hud('Wolfe condition parameters: c1='//num2str(c1)//', c2='//num2str(c2))

        self%alpha=alpha0

        !read setup        
        self%max_search=setup%get_int('MAX_SEARCH',o_default='12')

    end subroutine
    
    subroutine search(self,if_reinit_alpha,iterate,curr,pert,o_gradient_history)
        use mpi
        class(t_linesearcher) :: self
        logical :: if_reinit_alpha
        type(t_querypoint) :: curr,pert
        real,dimension(:,:,:,:,:),optional :: o_gradient_history

        type(t_checkpoint),save :: chp

        !alpha probably can't be outside [alphaL0,alphaR0]
        !a prudent range as we've appropriately scaled the optimization problem
        !by the norm of the gradient (linesearcher scaler)
        !if alpha < alphaL0 then probably no significant model updates
        !if alpha > alphaR0 then probably too-large model updates
        self%alphaL=alphaL0  !0.
        self%alphaR=alphaR0  !huge(1.)

        !reset alpha = alpha0 if requested by the optimization algorithm
        if(if_reinit_alpha) self%alpha=alpha0
        
        ! if(mpiworld%is_master) write(*,'(a,3(2x,es8.2))') ' Initial alphaL/alpha/alphaR =',alphaL,alpha,alphaR
        call hud('Initial alpha = '//num2str(self%alpha))
        call hud('Current qp%f, ║g║₁ = '//num2str(curr%f)//', '//num2str(sum(abs(curr%g))))

        !save gradients
        if(present(o_gradient_history)) then
            l=size(o_gradient_history,5) !number of gradient in history
            i=1
            o_gradient_history(:,:,:,:,i)=curr%g !?? and why not eoshift(o_gradient_history,1,curr%g)
            i=i+1; if(i>l) i=1
        endif
        
        call hud('------------ START LINESEARCH ------------')
        !linesearch loop
        loop: do isearch=1,self%max_search
            self%isearch=isearch !gfortran requires so

            !perturb current point
            pert%x = curr%x + self%alpha*curr%d
            !if(mpiworld%is_master) call sysio_write('pert%x',pert%x,size(pert%x),o_mode='append')
            call threshold(pert%x,size(pert%x))
            call hud('Modeling with perturbed models')
            call chp%init('FWI_querypoint_linesearcher','Gradient#','per_init')
            if(.not.pert%is_registered(chp)) then
                call fobj%eval(pert)
                call pert%register(chp)
            endif

            ! if(.not. curr%is_fitting_data) then
            !     call hud('Negate the sign of pert due to curr')
            !     call pert%set_sign(o_sign='-')
            ! endif
            
            call self%scale(pert)

            pert%g_dot_d = sum(pert%g*curr%d)

            self%igradient=self%igradient+1

            !save gradients
            if(present(o_gradient_history)) then
                o_gradient_history(:,:,:,:,i)=pert%g
                i=i+1; if(i>l) i=1
            endif

            call hud('Iterate.LineSearch.Gradient# = '//num2str(iterate)//'.'//num2str(self%isearch)//'.'//num2str(self%igradient))

            !if(mpiworld%is_master) write(*,'(a,3(2x,es8.2))') ' Linesearch alphaL/alpha/alphaR =',alphaL,alpha,alphaR
            call hud('alpha (α) = '//num2str(self%alpha)//' in ['//num2str(self%alphaL)//','//num2str(self%alphaR)//']')
            call hud('Perturb qp%f, ║g║₁ = '//num2str(pert%f)//', '//num2str(sum(abs(pert%g))))
            
            !Wolfe conditions
            if_1st_cond = (pert%f <= curr%f+c1*self%alpha*curr%g_dot_d) !sufficient descent condition
            !if_1st_cond = (pert%f <= curr%f)
            if_2nd_cond = (abs(pert%g_dot_d) <= c2*abs(curr%g_dot_d)) !strong curvature condition
            !if_2nd_cond = (pert%g_dot_d >= c2*curr%g_dot_d) !weak curvature condition (note the diff of inequal sign..)

            !occasionally optimizers on processors don't have same behavior
            !try to avoid this by broadcast controlling logicals.
            call mpi_bcast(if_1st_cond, 1, mpi_logical, 0, mpiworld%communicator, mpiworld%ierr)
            call mpi_bcast(if_2nd_cond, 1, mpi_logical, 0, mpiworld%communicator, mpiworld%ierr)

            call self%write(iterate,pert)

            !1st condition OK, 2nd condition OK => use alpha
            if(if_1st_cond .and. if_2nd_cond) then
                call hud('Wolfe conditions are satisfied')
                call hud('Enter new iterate')
                self%result='success'
                exit loop
            endif
            
            !bracketing-dichotomy strategy to find alpha
            if(self%isearch<self%max_search-1) then
            
                !1st condition BAD => [ <-]
                if(.not. if_1st_cond) then
                    call hud("Sufficient descent condition NOT satified. Now try a smaller alpha (α)")
                    self%result='perturb'
                    self%alphaR=self%alpha
                    self%alpha=0.5*(self%alphaL+self%alphaR) !shrink the search interval
                endif
                
                !2nd condition BAD => [-> ]
                if(if_1st_cond .and. .not. if_2nd_cond) then
                    call hud("Curvature condition is NOT satified. Now try a larger alpha (α)")
                    self%result='perturb'
                    self%alphaL=self%alpha
                    if(self%alphaR < alphaR0) then
                        self%alpha=0.5*(self%alphaL+self%alphaR) !shrink the search interval
                    else
                        self%alpha=self%alpha*multiplier !extend search interval
                    endif
                endif
                
            endif

            !exit loop, 1st condition has priority over 2nd condition
            !unless re-use alphaL if necessary
            if(self%isearch==self%max_search-1) then
                if(.not. if_1st_cond) then
                    call hud("Sufficient decrease condition NOT satified && LinS# = max_search-1.")
                    if(self%alphaL>alphaL0) then
                        call hud("Redo linesearch with alpha = alphaL")
                        self%alpha=self%alphaL
                    else
                        call hud("alphaL==alphaL0. No need to do the last search.")
                        call hud("Linesearch failure: can't find good alpha.")
                        self%result='failure'
                        exit loop
                    endif
                
                else
                    if(.not. if_2nd_cond) call hud("LinS# = max_search-1. Use alpha although curvature condition is not satisfied")
                    call hud('Enter new iterate')
                    self%result='success'
                    exit loop

                endif

            endif

            !last
            if(self%isearch==self%max_search) then
                if(.not. if_2nd_cond) call hud("Linesearch max_search reached. Use alpha although curvature condition is not satisfied")
                call hud('Enter new iterate')
                self%result='success'
                exit loop
                
            endif

        enddo loop


        ! if(pert%is_fitting_data) then
        !     call hud('Set positive sign to pert')
        !     call pert%set_sign(o_sign='+')
        ! else
        !     call hud('Set negative sign to pert')
        !     call pert%set_sign(o_sign='-')
        ! endif

        
        if(self%igradient>=self%max_gradient) then
            call hud('Maximum number of gradients reached. Finalize program now..')
            self%result='maximum'
        endif       
    
    end subroutine

    subroutine threshold(x,n)
        real,dimension(n) :: x

        logical,dimension(n) :: is_reached

        is_reached=.false.

        !x should reside in [0,1]
        where (x<0.) 
            is_reached=.true.
            x=thres
        endwhere

        where (x>1.)
            is_reached=.true.
            x=1.-thres
        endwhere

        if( size(pack(is_reached,.false.))>n/10 ) then
            call warn('x significantly reaches {0,1}. Something might be wrong.')
        endif

    end subroutine
    
    !scale the problem s.t. 
    !qp%g, qp%pg, to be negated as qp%d, have a similar scale (or unit) as qp%x, 
    !and alpha0 can be simply 1 (unitless).
    !update=alpha*qp%d=-alpha*qp%pg
    subroutine scale(self,qp)
        class(t_linesearcher) :: self
        type(t_querypoint) :: qp

        logical,save :: is_first_in=.true.
        character(:),allocatable :: str
        
        if(is_first_in) then

            !first do sanity check
            if(sum(abs(qp%pg))==0.) then
                call error('Gradient becomes absolutely zero!')
            endif

            str=setup%get_str('LINESEARCHER_SCALING','LS_SCALING',o_default='by total_volume/||pg(1)||1')

            if(str=='by total_volume/||pg(1)||1') then
                !total_volume = ∫   1     dy³ = n1*n2*n3  *d1*d2*d3
                !║pg(1)║₁     = ∫ |pg(1)| dy³ = Σ |pg(1)| *d1*d2*d3
                eps=param%n1*param%n2*param%n3/sum(abs(qp%pg(:,:,:,1)))
                
            elseif(len(str)>0) then
                eps=str2real(str)
                eps=eps/maxval(abs(qp%pg(:,:,:,1))) !eg. =0.05/║pg(1)║∞, ie. maximum 50m/s perturbation on velocity
            else
                call error('LINESEARCHER_SCALING input is zero!')
            endif
            
            self%scaler=eps/param%pars(1)%range
            
            call hud('Linesearch scaler= '//num2str(self%scaler))

            is_first_in=.false.

        endif
        
        qp%f=qp%f*self%scaler
        qp%g=qp%g*self%scaler
        qp%pg=qp%pg*self%scaler
        
    end subroutine

    subroutine write(self,iterate,pert)
        class(t_linesearcher) :: self
        integer :: iterate
        type(t_querypoint) :: pert

        character(*),parameter :: fmt='(x,5x,2x,i5,2x,i5,2x,es8.2,6x,es11.4,22x,es9.2,3x,l,x,l)'

        if(mpiworld%is_master) then
            
            open(16,file=dir_out//'optimization.log',position='append',action='write')
            write(16,fmt)  ls%isearch, ls%igradient, ls%alpha, pert%f, pert%g_dot_d, if_1st_cond, if_2nd_cond
            !write(16,'(a)') 'pert%f (pert%f-curr%f)/α pert%g·d curr%g·d'
            !write(16,*) pert%f, (pert%f-curr%f)/self%alpha, pert%g_dot_d, curr%g_dot_d
            close(16)

        endif

        if(setup%get_bool('IF_LINESEARCH_WRITE','IF_LS_WRITE',o_default='F')) then
            call sysio_write('model_Iter'//num2str(iterate)//'.LinS'//num2str(ls%isearch),m%vp,m%n)
            call sysio_write('image_Iter'//num2str(iterate)//'.LinS'//num2str(ls%isearch),m%image,m%n)
            call sysio_mv( 'Ru_Shot0001.su', 'Ru_Shot0001_Iter'//num2str(iterate)//'.LinS'//num2str(ls%isearch)//'.su')
            call sysio_mv('Rdu_Shot0001.su','Rdu_Shot0001_Iter'//num2str(iterate)//'.LinS'//num2str(ls%isearch)//'.su')
        endif

    end subroutine
    

end
