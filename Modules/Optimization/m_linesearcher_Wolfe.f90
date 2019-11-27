module m_linesearcher
use m_mpienv
use m_sysio
use m_gradient
use m_parameterization
use m_preconditioner

    private
    public alpha, isearch, imodeling, max_search, max_modeling, if_project_x, if_reinitialize_alpha, if_scaling, n
    public t_forwardmap, linesearcher, linesearch_scaling
    
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
    integer,parameter :: multiplyfactor=10
    
    !steplength
    real,parameter :: alpha0=1.
    real :: alpha=alpha0  !very initial steplength
    real :: alphaL=0., alphaR=0.
    
    logical :: if_reinitialize_alpha=.false.
    
    !counter
    integer :: isearch=0 !number of linesearch performed
    integer :: imodeling=1 !number of modeling performed
    !integer,parameter :: max_search=20 !max number of linesearch allowed
integer,parameter :: max_search=12
    integer :: max_modeling  !max number of modeling allowed
    
    !projection
    logical :: if_project_x=.true.
    real,parameter :: threshold=0.
    
    !scaling
    logical :: if_scaling=.true.
    
    type t_forwardmap
        real,dimension(:),allocatable :: x !point
        real,dimension(:),allocatable :: g !gradient
        real,dimension(:),allocatable :: pg !preconditioned gradient
        real,dimension(:),allocatable :: d !descent direction
        real :: f
        real :: gdotd  !dot product of g and d
    end type
    
    integer n !problem size
    
    contains
    
    subroutine linesearcher(iterate,current,perturb,gradient_history,result)
        integer,intent(in) :: iterate
        type(t_forwardmap),intent(inout) :: current, perturb
        real,dimension(:,:),intent(inout),optional :: gradient_history
        character(7) :: result
        
        logical :: first_condition, second_condition
        
        !initialize
        alphaL=0.
        alphaR=0.
        if(if_reinitialize_alpha) alpha=alpha0 !reinitialize alpha in each iterate may help convergence for LBFGS method..
!         if(mpiworld%is_master) write(*,'(a,3(2x,es8.2))') ' Initial alphaL/alpha/alphaR =',alphaL,alpha,alphaR
        if(mpiworld%is_master) write(*,'(a, 2x,es8.2)')   ' Initial alpha =',alpha
        if(mpiworld%is_master) write(*,*) 'Initial f,||g|| =',current%f,norm2(current%g)
        perturb%x=current%x+alpha*current%d
        
        !save gradients
        if(present(gradient_history)) then
            l=size(gradient_history,2) !number of gradient in history
            i=1
            gradient_history(:,i)=current%g
            i=i+1; if(i>l) i=1;
        endif
        
        !linesearch loop
        loop: do isearch=1,max_search
        
            if(mpiworld%is_master) write(*,'(a,3(2x,i5))') '  Iteration.Linesearch.Modeling#',iterate,isearch,imodeling
            
            call hud('Modeling with perturbed parameters')
            if(if_project_x) call linesearcher_project(perturb%x)
            call parameterization_transform('x2m',perturb%x)
            call gradient_modeling(if_gradient=.true.)
            perturb%f=fobjective
            call parameterization_transform('m2x',perturb%x,perturb%g)
            if(if_scaling) call linesearch_scaling(perturb)
            call preconditioner_apply(perturb%g,perturb%pg)
            imodeling=imodeling+1
            
            call hud('Judge alpha by Wolfe conditions')
            
            perturb%gdotd=sum(perturb%g*current%d)
            
            !Wolfe conditions
            first_condition = (perturb%f <= (current%f+c1*alpha*current%gdotd))
            !second_condition= (abs(perturb%gdotd) >= c2*abs(current%gdotd)) !strong Wolfe condition
            second_condition= (perturb%gdotd >= c2*current%gdotd) !strong Wolfe condition
            
! print*,'1st cond',perturb%f,current%f,current%gdotd, first_condition
! print*,'2nd cond',perturb%gdotd, second_condition
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
                    result='perturb'
                    alphaR=alpha
                    alpha=0.5*(alphaL+alphaR)
                    perturb%x=current%x+alpha*current%d
                    !save gradient
                    if(present(gradient_history)) then
                        gradient_history(:,i)=perturb%g
                        i=i+1; if(i>l) i=1;
                    endif
                endif
                
                !2nd condition BAD => increase alpha unless alphaR=0.
                if(first_condition .and. .not. second_condition) then
                    call hud("Curvature condition is not satified. Now try a larger alpha")
                    result='perturb'
                    alphaL=alpha
                    if(abs(alphaR) > tiny(1.e0)) then
                        alpha=0.5*(alphaL+alphaR)
                    else
                        alpha=alpha*multiplyfactor
                    endif
                    perturb%x=current%x+alpha*current%d
                    !save gradient
                    if(present(gradient_history)) then
                        gradient_history(:,i)=perturb%g
                        i=i+1; if(i>l) i=1;
                    endif
                endif
                
            endif
            
            !if(mpiworld%is_master) write(*,'(a,3(2x,es8.2))') ' Linesearch alphaL/alpha/alphaR =',alphaL,alpha,alphaR
            if(mpiworld%is_master) write(*,'(a, 2x,es8.2)')   ' Linesearch alpha =',alpha
            if(mpiworld%is_master) write(*,*) 'f,||g|| =',perturb%f,norm2(perturb%g)
        
        enddo loop
        
        if (imodeling>=max_modeling) then
            call hud('Maximum modeling number reached. Finalize program now..')
            result='maximum'
        endif       
    
    end subroutine

    
    
    subroutine linesearcher_project(x)
        real,dimension(n) :: x

        !x should be inside [0,1] due to scaling
        where (x<0.) x=threshold
        where (x>1.) x=1.-threshold

        !x should be fixed in the mask area (e.g. water)
        call parameterization_applymask(x)

    end subroutine
    
    
    !scale problem
    subroutine linesearch_scaling(fm)
        type(t_forwardmap) :: fm
        
        real,save :: scaling=1.
        logical,save :: if_computed_scale=.false.
        
        if(.not.if_computed_scale) then
        
            !scaling=1e3* 1e-2*m%n/ (sum(abs(fm%g(1:m%n)))) / (par_vp_max -par_vp_min)  !=1e3* |gp1|_L1 / |gvp|_L1 / (vpmax-vpmin)

            scaling=1e3* 1e-2*m%n/ (sum(abs(fm%g(1:m%n)))) / (pars_max(1)-pars_min(1))  !=1e3* |gpar1|_L1 / |gpar|_L1 / par1_range

            if(mpiworld%is_master) write(*,*) 'Linesearch Scaling Factor:', scaling

            if_computed_scale=.true.
        endif
        
        fm%f=fm%f*scaling
        fm%g=fm%g*scaling
        
        
    end subroutine
    

end
