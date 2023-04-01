module m_linesearcher
use m_System
use m_Modeling
use m_Kernel

    private
    
    character(*),parameter :: info="Enumerate on alpha's"
    
    !thresholding
    logical :: if_threshold
    real,parameter :: thres=0.
    
    real,dimension(:),allocatable :: alphas
    integer :: nalpha

    !Wolfe conditions
    logical :: if_1st_cond, if_2nd_cond

    type,public :: t_linesearcher
        real :: alpha  !steplength
        real :: scaler
        
        !counter
        integer :: isearch !number of linesearch performed in each iterate
        integer :: max_search !max number of linesearch allowed per iteration
        integer :: igradient=1 !total number of gradient computed
        integer :: max_gradient !max total number of gradient computation allowed
    
        contains
        procedure :: init
        procedure :: search
        procedure :: scale
    end type

    type(t_linesearcher),public :: ls

    
    contains
    
    subroutine init(self)
        class(t_linesearcher) :: self

        call hud(info)
        
        alphas = setup%get_reals('TEST_ALPHAS','ALPHAS',o_default='0.001 0.01 0.1 1')
        nalpha = size(alphas)
        call hud('nalpha = '//num2str(nalpha))

        if_threshold = setup%get_bool('LINESEARCHER_THRESHOLD',o_default='T')
        
    end subroutine
    
    subroutine search(self,curr,pert)
        use mpi
        class(t_linesearcher) :: self
        type(t_querypoint) :: curr,pert
        
        type(t_checkpoint),save :: chp

        logical :: if_1st_cond, if_2nd_cond
        
        call hud('Initial alpha = '//num2str(self%alpha))
        call hud('Current qp%f, ║g║₁ = '//num2str(curr%f)//', '//num2str(sum(abs(curr%g))))
        
        call hud('------------ START LINESEARCH ------------')
        !linesearch loop
        loop: do isearch=1,nalpha
            self%isearch=isearch !gfortran requires so

            self%alpha=alphas(isearch)

            !perturb current point
            pert%x = curr%x + self%alpha*curr%d
            call threshold(pert%x,size(pert%x))
            if(mpiworld%is_master) call sysio_write('pert%x',pert%x,size(pert%x),o_mode='append')
            call hud('Modeling with perturbed models')

            ! call chp%init('FWI_querypoint_linesearcher','Gradient#','per_init')
            ! if(.not.pert%is_registered(chp)) then
                call fobj%eval(pert)
                ! call pert%register(chp)
            ! endif

            call self%scale(pert)

            pert%g_dot_d = sum(pert%g*curr%d)

            self%igradient=self%igradient+1

            !save gradients to disk
            if(mpiworld%is_master) call sysio_write('pert%g',pert%g,size(pert%g),o_mode='append')

            call hud('Linesearch.Gradient# = '//num2str(self%isearch)//'.'//num2str(self%igradient))
            
            call hud('alpha = '//num2str(self%alpha))
            call hud('Perturb qp%f, ║g║₁ = '//num2str(pert%f)//', '//num2str(sum(abs(pert%g))))
            
            !Wolfe conditions
            if_1st_cond = (pert%f <= curr%f+c1*self%alpha*curr%g_dot_d) !sufficient descent condition
            !if_1st_cond = (pert%f <= curr%f)
            !if_2nd_cond = (abs(pert%g_dot_d) >= c2*abs(curr%g_dot_d)) !strong curvature condition
            if_2nd_cond = (pert%g_dot_d >= c2*curr%g_dot_d) !weak curvature condition

            if(mpiworld%is_master) then
                print*,'1st cond',self%alpha,pert%f,curr%f,(pert%f-curr%f)/self%alpha, curr%g_dot_d,  curr%g_dot_d/((pert%f-curr%f)/self%alpha), if_1st_cond
                print*,'2nd cond',pert%g_dot_d, curr%g_dot_d, if_2nd_cond
            endif

            !occasionally optimizers on processors don't have same behavior
            !try to avoid this by broadcast controlling logicals.
            call mpi_bcast(if_1st_cond, 1, mpi_logical, 0, mpiworld%communicator, mpiworld%ierr)
            call mpi_bcast(if_2nd_cond, 1, mpi_logical, 0, mpiworld%communicator, mpiworld%ierr)

            !1st condition OK, 2nd condition OK => use alpha
            if(if_1st_cond .and. if_2nd_cond) then
                call hud('Wolfe conditions are satisfied')
                cycle loop
            endif
            
            !1st condition BAD
            if(.not. if_1st_cond) then
                call hud("Sufficient descent condition NOT satified -> should try a smaller alpha")
            endif
            
            !2nd condition BAD
            if(if_1st_cond .and. .not. if_2nd_cond) then
                call hud("Curvature condition NOT satified -> should try a larger alpha")
            endif

        enddo loop
    
    end subroutine

    subroutine threshold(x,n)
        real,dimension(n) :: x

        logical,dimension(n) :: is_reached

        if(.not.if_threshold) return
        
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
            
            elseif(str=='none') then
                eps=param%pars(1)%range !st self%scaler=1.

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
    

end
