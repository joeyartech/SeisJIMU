module m_linesearcher
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
    
    !Wolfe conditions parameters (Nocedal value)
    real,parameter :: c1=1e-4, c2=0.9 !for (quasi-)Newton method, 0.1 for NLCG
    !Bracketting parameter (Gilbert value)
    integer,parameter :: multiplier=10
    
    !default
    real,parameter :: alpha0=1. !steplength
    logical :: if_reinit_alpha=.false. !reset alpha=1 in each new iterate

    !thresholding
    real,parameter :: thres=0.
    
    type,public :: t_forwardmap
        real,dimension(:),allocatable :: x !point in model manifold
        real,dimension(:),allocatable :: g !gradient
        real,dimension(:),allocatable :: pg !preconditioned gradient
        real,dimension(:),allocatable :: d !descent direction
        real :: f !value of objective function
        real :: gdotd  !g.d

        contains
        procedure :: init => forwardmap_init
        final :: forwardmap_fin
    end type

    type,public :: t_linesearcher
        real :: alpha=alpha0 !steplength initializd to alpha0
        real :: alphaL=0., alphaR=0.
        real :: scaler=1.
        logical :: is_scaled=.false.
        
        !counter
        integer :: imodeling=1 !number of modeling performed
        integer :: isearch=0 !number of linesearch performed
        integer :: iterate=0 !number of iteration performed    
        
        integer max_modeling !max total number of forward modeling allowed
        integer max_search   !max number of linesearch allowed per iteration
        integer max_iterate  !max number of iteration allowed

        contains
        procedure :: init => init
        procedure :: search => search
        procedure :: scale => scale
    end type

    type(t_linesearcher),public :: ls

    contains

    subroutine forwardmap_init(self,n,o_fobj)
        class(t_forwardmap) :: self

        real,optional :: o_fobj

        call alloc(self%x, n)
        call alloc(self%g, n)
        call alloc(self%pg,n)
        call alloc(self%d, n)

        if(present(o_fobj)) self%f=o_fobj
        
    end subroutine

    subroutine forwardmap_fin(self)
        type(t_forwardmap) :: self

        if(allocated(self%x )) deallocate(self%x )
        if(allocated(self%g )) deallocate(self%g )
        if(allocated(self%pg)) deallocate(self%pg)
        if(allocated(self%d )) deallocate(self%d )

    end subroutine

    subroutine init(self)
        class(t_linesearcher) :: self

        !read setup
        max_iterate=setup%get_int('MAX_ITERATE',o_default='30')
        max_search=setup%get_int('MAX_SEARCH',o_default='12')
        max_modeling=setup%get_int('MAX_MODELING',o_default=num2str(max_iterate+30))
        
        if_reinit_alpha=setup%get_bool('IF_REINIT_ALPHA',o_default='F')

    end subroutine
    
    subroutine search(self,iterate,curr,pert,gradient_history,result)
        class(t_linesearcher) :: self
        integer,intent(in) :: iterate
        type(t_forwardmap),intent(inout) :: curr, pert
        real,dimension(:,:),intent(inout),optional :: gradient_history
        character(7) :: result
        
        logical :: first_condition, second_condition
        
        !initialize
        alphaL=0.
        alphaR=0.
        if(if_reinit_alpha) alpha=alpha0 !reinitialize alpha in each iterate may help convergence for LBFGS method..
        ! if(mpiworld%is_master) write(*,'(a,3(2x,es8.2))') ' Initial alphaL/alpha/alphaR =',alphaL,alpha,alphaR
        if(mpiworld%is_master) write(*,'(a, 2x,es8.2)')   ' Initial alpha =',alpha
        if(mpiworld%is_master) write(*,*) 'Initial f,||g|| =',curr%f,norm2(curr%g)
        pert%x=curr%x+alpha*curr%d
        
        !save gradients
        if(present(gradient_history)) then
            l=size(gradient_history,2) !number of gradient in history
            i=1
            gradient_history(:,i)=curr%g
            i=i+1; if(i>l) i=1;
        endif
        
        !linesearch loop
        loop: do isearch=1,max_search
        
            if(mpiworld%is_master) write(*,'(a,3(2x,i5))') '  Iteration.Linesearch.Modeling#',iterate,isearch,imodeling
            
            call hud('Modeling with perturbed parameters')
            call threshold(pert%x,size(pert%x))
            call param%transform('x2m',pert%x)
            call nabla%act(fobj)
            pert%f=fobj%total_misfit
            call param%transform('m2x',pert%x,pert%g)
            call self%scale(pert)
            call preco%apply(pert%g,pert%pg)
            imodeling=imodeling+1
            
            call hud('Judge alpha by Wolfe conditions')
            
            pert%gdotd=sum(pert%g*curr%d)
            
            !Wolfe conditions
            first_condition = (pert%f <= (curr%f+c1*self%alpha*curr%gdotd))
            !second_condition= (abs(pert%gdotd) >= c2*abs(curr%gdotd)) !strong Wolfe condition
            second_condition= (pert%gdotd >= c2*curr%gdotd) !strong Wolfe condition
            
            ! print*,'1st cond',pert%f,curr%f,curr%gdotd, first_condition
            ! print*,'2nd cond',pert%gdotd, second_condition
            ! print*,'alpha(3)',alphaL,alpha,alpha_R

            !occasionally optimizers on processors don't have same behavior
            !try to avoid this by broadcast controlling logicals.
            call mpi_bcast(first_condition,  1, mpi_logical, 0, mpiworld%communicator, mpiworld%ierr)
            call mpi_bcast(second_condition, 1, mpi_logical, 0, mpiworld%communicator, mpiworld%ierr)

            !1st condition OK, 2nd condition OK => use alpha
            if(first_condition .and. second_condition) then
                call hud('Wolfe conditions are satisfied')
                call hud('Enter new iterate')
                result='success'
                exit loop
            endif
            
            if(isearch==max_search) then
            
                !1st condition BAD => failure
                if(.not. first_condition) then
                    call hud("Linesearch failure: can't find good alpha.")
                    result='failure'
                    exit loop
                endif
                
                !2nd condition BAD => use alpha
                if(.not. second_condition) then
                    call hud("Maximum linesearch number reached. Use alpha although curvature condition is not satisfied")
                    call hud('Enter new iterate')
                    result='success'
                    exit loop
                endif
                
            else
            
                !1st condition BAD => shrink alpha
                if(.not. first_condition) then
                    call hud("Armijo condition is not satified. Now try a smaller alpha")
                    result='pert'
                    self%alphaR=self%alpha
                    self%alpha=0.5*(self%alphaL+self%alphaR)
                    pert%x=curr%x+self%alpha*curr%d
                    !save gradient
                    if(present(gradient_history)) then
                        gradient_history(:,i)=pert%g
                        i=i+1; if(i>l) i=1;
                    endif
                endif
                
                !2nd condition BAD => increase alpha unless alphaR=0.
                if(first_condition .and. .not. second_condition) then
                    call hud("Curvature condition is not satified. Now try a larger alpha")
                    result='perturb'
                    self%alphaL=self%alpha
                    if(abs(self%alphaR) > tiny(1.e0)) then
                        self%alpha=0.5*(self%alphaL+self%alphaR)
                    else
                        self%alpha=self%alpha*multiplier
                    endif
                    pert%x=curr%x+alpha*curr%d
                    !save gradient
                    if(present(gradient_history)) then
                        gradient_history(:,i)=pert%g
                        i=i+1; if(i>l) i=1;
                    endif
                endif
                
            endif
            
            !if(mpiworld%is_master) write(*,'(a,3(2x,es8.2))') ' Linesearch alphaL/alpha/alphaR =',alphaL,alpha,alphaR
            if(mpiworld%is_master) write(*,'(a, 2x,es8.2)')   ' Linesearch alpha =',self%alpha
            if(mpiworld%is_master) write(*,*) 'f,||g|| =',pert%f,norm2(pert%g)
        
        enddo loop
        
        if (imodeling>=max_modeling) then
            call hud('Maximum modeling number reached. Finalize program now..')
            result='maximum'
        endif       
    
    end subroutine

    subroutine threshold(x,n)
        real,dimension(n) :: x

        !x should be inside [0,1] due to scaling
        where (x<0.) x=thres
        where (x>1.) x=1.-thres

    end subroutine
    
    
    !scale problem
    subroutine scale(self,fm)
        class(t_linesearcher) :: self
        type(t_forwardmap) :: fm
        
        if(.not.self%is_scaled) then

            !self%scaler=1e3* 1e-2*m%n/ (sum(abs(fm%g(1:m%n)))) / (par_vp_max -par_vp_min)  !=1e3* |gp1|_L1 / |gvp|_L1 / (vpmax-vpmin)

            self%scaler=1e3* 1e-2*m%n/ (sum(abs(fm%g(1:m%n)))) / (pars_max(1)-pars_min(1))  !=1e3* |gpar1|_L1 / |gpar|_L1 / par1_range

            if(mpiworld%is_master) write(*,*) 'Linesearch Scaling Factor:', scaling

            self%is_scaled=.true.

        endif
        
        fm%f=fm%f*self%scaler
        fm%g=fm%g*self%scaler
        
    end subroutine
    

end
