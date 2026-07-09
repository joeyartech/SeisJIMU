subroutine optimizer_init_loop(qp0)        
use m_System
use m_Modeling
use m_Kernel
use m_linesearcher

    type(t_querypoint),target :: qp0 !initial (model) parameters
    type(t_querypoint),target :: qp1
    type(t_querypoint),pointer :: curr, pert
    character(:),allocatable :: str

    !subroutine optimizer_init:
    
        !current point
        curr=>qp0
        !choose descent direction
        str=setup%get_str('DESCENT_DIR',o_default='-curr%pg')
        
        if(str=='-curr%pg') then !steepest descent
            curr%d=-curr%pg

        elseif(str=='stack') then !stacking -curr%pg
            call alloc(curr%d,param%n1,param%n2,param%n3,param%npars)
            do ipar=1,param%npars
            do i1=1,param%n1
                tmp=-sum(curr%pg(i1,:,:,ipar))
                curr%d(i1,:,:,ipar)=tmp
            enddo
            enddo

            call sysio_write('curr%d',curr%d,size(curr%d))

            
        elseif(str=='random') then !random descent
            allocate(curr%d,source=curr%pg)
            call random_number(curr%d) ! ∈[0,1)
            curr%d=curr%d*2.-1. ! ∈[-1,1)
            
        elseif(str=='random_perturb') then !random perturbation to steepest descent
            allocate(curr%d,source=curr%pg)
            call random_number(curr%d) ! ∈[0,1)
            curr%d=curr%d*2.-1. ! ∈[-1,1)
            
            curr%d = sum(abs(curr%pg))/sum(abs(curr%d)) *curr%d
            curr%d = -curr%pg +0.5*curr%d
            
        elseif(str=='random_normal') then !random normal direction to steepest descent
            allocate(curr%d,source=curr%pg)
            call random_number(curr%d) !unlikely to // with curr%pg
            !   d·pg/‖pg‖ = ‖d‖cosθ = projection of d onto pg
            !d-(d·pg/‖pg‖)pg/‖pg‖ = normal direction to pg
            tmp=sum(curr%d*curr%pg)/norm2(curr%pg)
            curr%d = curr%d - tmp*curr%pg
            
        endif
        
        !ensure good magnitudes
        curr%d = sum(abs(curr%pg))/sum(abs(curr%d)) *curr%d
            
        curr%g_dot_d = sum(curr%g*curr%d) !inner product
        
        !perturbed point
        pert=>qp1
        !perturbed querypoint
        call pert%init('qp')
    
    !subroutine optimizer_loop:
        !linesearch
        call ls%search(curr,pert)
    
end subroutine
