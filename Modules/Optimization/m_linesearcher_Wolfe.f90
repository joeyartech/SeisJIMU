module m_linesearcher
use m_System
use m_Modeling
use m_Kernel

    private
    
    ! This linesearch enforces the Wolfe          !
    ! conditions: sufficient decrease             !
    !             sufficient curvature            !
    ! The Wolfe conditions can be found in        !
    ! Nocedal, Numerical Optimization, 2nd        !
    ! edition, p.33                               !
    !                                             !
    ! The linesearch method implemented here is   !
    ! based first on a bracketing strategy, then  !
    ! on a dichotomy algorithm. A full description!
    ! of this strategy can be found in            !
    ! Numerical Optimizationn Theoretical and     !
    ! Practical Aspects, J.F.Bonnans, J.C.Gilbert,!
    ! C. Lemaréchal, C.A. Sagastizábal,           !
    ! Springer-Verlag, Universitext               !
    
    character(*),parameter :: info='Linesearch based on a bracketing - dichotomy algoritm'//s_NL// &
                                   'Steplength (alpha) judged by Wolfe conditions'

    !Wolfe conditions parameters (Nocedal value)
    real,parameter :: c1=1e-4, c2=0.9 !for (quasi-)Newton method, 0.1 for NLCG
    !Bracketting parameter (Gilbert value)
    integer,parameter :: multiplier=10
    
    !thresholding
    real,parameter :: thres=0.

    type,public :: t_linesearcher
        real :: alpha0 !initial steplength
        real :: alpha  !steplength
        real :: alphaL, alphaR !search interval: [ ]
        real :: scaler
        logical :: if_has_scaler=.false.
        character(7) :: result
        
        !counter
        integer :: igradient=1 !total number of gradient computed
        integer :: isearch=0 !number of linesearch performed in each iterate
        
        integer :: max_gradient !max total number of gradient computation allowed
        integer :: max_search   !max number of linesearch allowed per iteration
        
        contains
        procedure :: init
        procedure :: search
        procedure :: scale
    end type

    type(t_linesearcher),public :: ls

    contains
    
    subroutine init(self,oif_reinit_alpha)
        class(t_linesearcher) :: self
        logical,optional :: oif_reinit_alpha

        call hud(info)

        !read setup
        self%max_search=setup%get_int('MAX_SEARCH',o_default='12')
        self%max_gradient=setup%get_int('MAX_GRADIENT',o_default=num2str(max_iterate+30))
        
        self%alpha0=setup%get_real('ALPHA0',o_default='1.')
        self%alpha=self%alpha0
        self%alphaL=0.
        self%alphaR=huge(1.)

    end subroutine
    
    subroutine search(self,iterate,curr,pert,o_gradient_history)
        use mpi
        class(t_linesearcher) :: self
        type(t_querypoint) :: curr,pert
        real,dimension(:,:,:,:,:),optional :: o_gradient_history

        logical :: if_1st_cond, if_2nd_cond
        
        ! if(mpiworld%is_master) write(*,'(a,3(2x,es8.2))') ' Initial alphaL/alpha/alphaR =',alphaL,alpha,alphaR
        call hud('Initial alpha = '//num2str(self%alpha))
        call hud('Initial f, ||g|| = '//num2str(curr%f)//', '//num2str(norm2(curr%g)))

        !perturb current point
        pert%x = curr%x + self%alpha*curr%d
        call threshold(pert%x,size(pert%x))
        
        !save gradients
        if(present(o_gradient_history)) then
            l=size(o_gradient_history,5) !number of gradient in history
            i=1
            o_gradient_history(:,:,:,:,i)=curr%g
            i=i+1; if(i>l) i=1
        endif
        
        !linesearch loop
        loop: do isearch=1,self%max_search

            self%isearch=isearch
        
            if(mpiworld%is_master) write(*,'(a,3(2x,i5))') '  Iterate.Linesearch.Gradient#',iterate,self%isearch,self%igradient
            
            
            call hud('Modeling with perturbed models')
            call nabla%act(fobj,pert,oif_update_m=.true.,oif_approx=.true.)
            call self%scale(pert)
            self%igradient=self%igradient+1
            
            call hud('Judge alpha by Wolfe conditions')
            
            pert%gdotd=sum(pert%g*curr%d)
            
            !Wolfe conditions
            if_1st_cond = (pert%f <= (curr%f+c1*self%alpha*curr%gdotd))
            !if_2nd_cond = (abs(pert%gdotd) >= c2*abs(curr%gdotd)) !strong Wolfe condition
            if_2nd_cond = (pert%gdotd >= c2*curr%gdotd) !Wolfe condition
            
            ! print*,'1st cond',pert%f,curr%f,curr%gdotd, if_1st_cond
            ! print*,'2nd cond',pert%gdotd, if_2nd_cond
            ! print*,'alpha(3)',alphaL,alpha,alpha_R

            !occasionally optimizers on processors don't have same behavior
            !try to avoid this by broadcast controlling logicals.
            call mpi_bcast(if_1st_cond, 1, mpi_logical, 0, mpiworld%communicator, mpiworld%ierr)
            call mpi_bcast(if_2nd_cond, 1, mpi_logical, 0, mpiworld%communicator, mpiworld%ierr)
            ! call mpi_bcast(if_1st_cond, 1, mpi_logical, 0, mpi_comm_world, mpiworld%ierr)
            ! call mpi_bcast(if_2nd_cond, 1, mpi_logical, 0, mpi_comm_world, mpiworld%ierr)

            !1st condition OK, 2nd condition OK => use alpha
            if(if_1st_cond .and. if_2nd_cond) then
                call hud('Wolfe conditions are satisfied')
                call hud('Enter new iterate')
                self%result='success'
                exit loop
            endif
            
            if(self%isearch<self%max_search) then
            
                !1st condition BAD => [ <-]
                if(.not. if_1st_cond) then
                    call hud("Sufficient decrease condition is not satified. Now try a smaller alpha")
                    self%result='perturb'
                    self%alphaR=self%alpha
                    self%alpha=0.5*(self%alphaL+self%alphaR) !shrink the search interval
                endif
                
                !2nd condition BAD => [-> ]
                if(if_1st_cond .and. .not. if_2nd_cond) then
                    call hud("Curvature condition is not satified. Now try a larger alpha")
                    self%result='perturb'
                    self%alphaL=self%alpha
                    if(self%alphaR < huge(1.)) then
                        self%alpha=0.5*(self%alphaL+self%alphaR) !shrink the search interval
                    else
                        self%alpha=self%alpha*multiplier !extend search interval
                    endif
                endif

                !perturb current point
                pert%x = curr%x + self%alpha*curr%d
                call threshold(pert%x,size(pert%x))

                !save gradients
                if(present(o_gradient_history)) then
                    o_gradient_history(:,:,:,:,i)=pert%g
                    i=i+1; if(i>l) i=1
                endif
                
            endif

            if(self%isearch==self%max_search) then
            
                !1st condition BAD => failure
                if(.not. if_1st_cond) then
                    call hud("Linesearch failure: can't find good alpha.")
                    self%result='failure'
                    exit loop
                endif
                
                !2nd condition BAD => use alpha
                if(.not. if_2nd_cond) then
                    call hud("Maximum linesearch number reached. Use alpha although curvature condition is not satisfied")
                    call hud('Enter new iterate')
                    self%result='success'
                    exit loop
                endif
                
            endif

            !if(mpiworld%is_master) write(*,'(a,3(2x,es8.2))') ' Linesearch alphaL/alpha/alphaR =',alphaL,alpha,alphaR
            if(mpiworld%is_master) write(*,'(a, 2x,es8.2)')   ' Linesearch alpha =',self%alpha
            if(mpiworld%is_master) write(*,*) 'f,||g|| =',pert%f,norm2(pert%g)
        
        enddo loop
        
        if (self%igradient>=self%max_gradient) then
            call hud('Maximum number of gradients reached. Finalize program now..')
            self%result='maximum'
        endif       
    
    end subroutine

    subroutine threshold(x,n)
        real,dimension(n) :: x

        !x should be inside [0,1] due to scaling
        where (x<0.) x=thres
        where (x>1.) x=1.-thres

    end subroutine
        
    subroutine scale(self,qp)
        class(t_linesearcher) :: self
        type(t_querypoint) :: qp
        
        if(.not.self%if_has_scaler) then
            !self%scaler=1e3* 1e-2*m%n/ (sum(abs(qp%g(1:m%n)))) / (par_vp_max -par_vp_min)  !=1e3* |gp1|_L1 / |gvp|_L1 / (vpmax-vpmin)
            self%scaler=1e3* 1e-2*m%n/ (sum(abs(qp%g(:,:,:,1)))) / (param%pars(1)%range)  !=1e3* |gpar1|_L1 / |gpar|_L1 / par1_range
            if(mpiworld%is_master) write(*,*) 'Linesearch scaler:', self%scaler
            self%if_has_scaler=.true.
        endif

        qp%f=qp%f*self%scaler
        qp%g=qp%g*self%scaler

    end subroutine
    

end
