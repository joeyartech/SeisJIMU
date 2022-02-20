module m_optimizer_common
use m_System
use m_Modeling
use m_Kernel
use m_linesearcher

    type(t_querypoint),pointer :: curr, pert
    type(t_querypoint),target  :: qp1
    real,dimension(:,:,:,:),allocatable :: prev_d  !descent direction of previous point (no need for another t_querypoint instance)

    real f0 !initial fobjective
    real g0norm !initial gradient norm
    
    !counter
    integer :: iterate=0 !number of iteration performed 
    integer :: max_iterate  !max number of iteration allowed

    real min_descent !minimum descent allowed
    
    contains
    
    subroutine optimizer_init(qp0)        
        type(t_querypoint),target :: qp0 !initial (model) parameters

        !current point
        curr=>qp0
        curr%d=-curr%pg
if(mpiworld%is_master) print*,'minmax curr%g',minval(curr%g),maxval(curr%g)
if(mpiworld%is_master) print*,'minmax curr%d',minval(curr%d),maxval(curr%d)
        curr%g_dot_d = sum(curr%g*curr%d)

        !initial values
        f0=curr%f
        g0norm=norm2(curr%g)
        
        !perturbed point
        pert=>qp1
        call pert%init('qp')

        !descent direction of previous point
        prev_d=curr%d

        !read setup
        min_descent=setup%get_real('MIN_DESCENT',o_default='1e-8')
        max_iterate=setup%get_int('MAX_ITERATE',o_default='120')
        ls%max_gradient=setup%get_int('MAX_GRADIENT',o_default=num2str(max_iterate+30))        
        
    end subroutine

    subroutine switch_curr_pert
        type(t_querypoint),pointer :: tmp

        tmp=>curr
        curr=>pert
        pert=>tmp
    end subroutine

    subroutine optimizer_write(message,o_title)
        character(*) :: message
        character(*),optional :: o_title
        
        if(mpiworld%is_master) then
        
            select case (message)
            case('start')
                ! call execute_command_line('rm '//dir_out//'iterate.log', wait=.true.)
                open(16,file=dir_out//'iterate.log',position='append',action='write')
                write(16,'(a)'      ) ' **********************************************************************'
                if(present(o_title)) then
                write(16,'(a)'      ) '    '//o_title
                endif
                write(16,'(a)'      ) ' **********************************************************************'
                write(16,'(a,es8.2)') '     Min descent allowed      =',  min_descent
                write(16,'(a,i5)'   ) '     Max iterates allowed     =',  max_iterate
                write(16,'(a,i5)'   ) '     Max linesearches allowed =',  ls%max_search
                write(16,'(a,i5)'   ) '     Max gradients allowed    =',  ls%max_gradient
                write(16,'(a,es8.2)') '     Initial fobjective (f0)        =',  f0
                write(16,'(a,es8.2)') '     Initial gradient norm (||g0||) =',  g0norm
                write(16,'(a)'      ) ' **********************************************************************'
                write(16,'(a)'      ) '  Iter#      f         f/f0    ||g||/||g0||    alpha     nls  Gradient#'
                               !e.g.  !    0    z1.00E+00    1.00E+00    1.00E+00    1.00E+00      0       1
                write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, curr%f, curr%f/f0, norm2(curr%g)/g0norm, ls%alpha, ls%isearch, ls%igradient
                close(16)
                
            case('update')
                open(16,file=dir_out//'iterate.log',position='append',action='write')
                write(16,'(i5,4(4x,es8.2),2x,i5,3x,i5)')  iterate, curr%f, curr%f/f0, norm2(curr%g)/g0norm, ls%alpha, ls%isearch, ls%igradient
                close(16)
                
                call param%transform('x->m',o_x=curr%x)
                call m%write(o_suffix='_Iter'//int2str(iterate))
                call sysio_write('descent_Iter'//int2str(iterate),curr%d,size(curr%d))
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
                call execute_command_line('(cd '//dir_out//' ; ln -sf model_Iter'//int2str(iterate)//' model_final )')
                call sysio_write('descent_Iter'//int2str(iterate),curr%d,size(curr%d))
                call shot%write('dsyn_Iter'//int2str(iterate)//'_',shot%dsyn)

                write(*,'(a,i0.4)') 'ximage < model_Iter* n1=',m%nz

            end select

        endif
    
    end subroutine

end
