module m_optimizer_common
use m_System
use m_Modeling
use m_Kernel
use m_linesearcher

    type(t_querypoint),pointer :: curr, pert
    type(t_querypoint),target  :: qp1
    real,dimension(:,:,:,:),allocatable :: prev_d  !descent direction of previous point (no need for another t_querypoint instance)

    real f0 !initial fobjective
    real g0norm1 !initial gradient norm2
    
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
        curr%g_dot_d = sum(curr%g*curr%d)

        !initial values
        f0=curr%f
        g0norm1=sum(abs(curr%g))
        
        !perturbed point
        pert=>qp1
        call pert%init('qp')

        !previous descent direction
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

    subroutine optimizer_write(task,o_title)
    use m_shot

        character(*) :: task
        character(*),optional :: o_title

        character(*),parameter :: fmt='(x,i5,2x,5x,2x,i5,10x,6x,es11.4,x,f7.3,2x,f7.3,5x,es9.2)'
        character(:),allocatable :: siterate,sindex
        type(t_string),dimension(:),allocatable :: tmp
        integer,dimension(:),allocatable :: list

        
        if(mpiworld%is_master) then

            siterate=int2str(iterate)
        
            select case (task)
            case('start')
                ! call execute_command_line('rm '//dir_out//'iterate.log', wait=.true.)
                open(16,file=dir_out//'optimization.log',position='append',action='write')
                write(16,'(a)'      ) ' **********************************************************************'
                if(present(o_title)) write(16,'(a)'      ) '    '//o_title
                write(16,'(a)'       ) ' **********************************************************************'
                write(16,'(a,es10.4)') '     Min descent allowed      =',  min_descent
                write(16,'(a,i5)'    ) '     Max iterates allowed     =',  max_iterate
                write(16,'(a,i5)'    ) '     Max linesearches allowed =',  ls%max_search
                write(16,'(a,i5)'    ) '     Max gradients allowed    =',  ls%max_gradient
                write(16,'(a,es14.8)') '     Linesearch scaler        =',  ls%scaler
                write(16,'(a,es10.4)') '     Initial gradient L1norm (║g0║₁)  =',  g0norm1
                write(16,'(a)'       ) ' **********************************************************************'
                write(16,'(a)'       ) '         LinS#  Grad#     α    pert%       f                            g·d   Wolfe_cond'
                write(16,'(a)'       ) '  Iter#         Grad#          curr%       f    f/f0(%) ║g║₁/║g0║₁(%)   g·d'
                write(16,'(a)'       ) ' ========================================================================================'
                write(16,fmt)  iterate, ls%igradient, curr%f, curr%f/f0*100., sum(abs(curr%g))/g0norm1*100., curr%g_dot_d
                close(16)
                
            case('update')
                open(16,file=dir_out//'optimization.log',position='append',action='write')
                write(16,fmt)  iterate, ls%igradient, curr%f, curr%f/f0*100., sum(abs(curr%g))/g0norm1*100., curr%g_dot_d
                write(16,'(a)'      ) ' ----------------------------------------------------------------------------------------'
                close(16)

                call param%transform('x->m',o_x=curr%x)
                call m%write(o_suffix='_Iter'//siterate)
                call sysio_write('descent_Iter'//siterate,curr%d,size(curr%d))
                
                call sysio_mv('updated_wavelet.su','updated_wavelet_Iter'//siterate//'.su')

                tmp=split(shls%all)
                list=setup%get_ints('ITERATE_KEEP_SHOT',&
                    o_default='1 '// either( tmp(size(tmp)/2+1)%s, '', size(tmp)>1 ))

                do i=1,size(list)
                    sindex=num2str(list(i),o_format='(i0.4)')
                    call sysio_mv('dsyn_Shot'//sindex//'.su','dsyn_Iter'//siterate//'_Shot'//sindex//'.su')
                    call sysio_mv( 'updated_RE0_Shot'//sindex//'.su', 'updated_RE0_Iter'//siterate//'_Shot'//sindex//'.su')
                    call sysio_mv( 'updated_RdE_Shot'//sindex//'.su', 'updated_RdE_Iter'//siterate//'_Shot'//sindex//'.su')
                    call sysio_mv('dadj_Shot'//sindex//'.su','dadj_Iter'//siterate//'_Shot'//sindex//'.su')
                enddo

            case('maximum')
                open(16,file=dir_out//'optimization.log',position='append',action='write')
                write(16,'(a)'      ) ' **********************************************************************'
                write(16,'(a)'      ) '     STOP: MAXIMUM MODELING NUMBER REACHED                             '
                write(16,'(a)'      ) ' **********************************************************************'
                close(16)

            case('criteria')
                open(16,file=dir_out//'optimization.log',position='append',action='write')
                write(16,'(a)'      ) ' **********************************************************************'
                write(16,'(a)'      ) '     STOP: CONVERGENCE CRITERIA SATISFIED                              '
                write(16,'(a)'      ) ' **********************************************************************'
                close(16)

            case('failure')
                open(16,file=dir_out//'optimization.log',position='append',action='write')
                write(16,'(a)'      ) ' **********************************************************************'
                write(16,'(a)'      ) '     STOP: LINE SEARCH FAILURE                                         '
                write(16,'(a)'      ) ' **********************************************************************'
                close(16)

            case('finalize')
                open(16,file=dir_out//'optimization.log',position='append',action='write')
                write(16,'(a)'      ) ' **********************************************************************'
                write(16,'(a)'      ) '     FINALIZE                                                          '
                write(16,'(a)'      ) ' **********************************************************************'
                close(16)

                call param%transform('x->m',o_x=pert%x)
                call m%write(o_suffix='_Iter'//int2str(iterate))
                call execute_command_line('(cd '//dir_out//' ; ln -sf model_Iter'//int2str(iterate)//' model_final )')
                call sysio_write('descent_Iter'//int2str(iterate),curr%d,size(curr%d))

                write(*,'(a,i0.4)') 'ximage < model_Iter* n1=',m%nz

                call sysio_mv('updated_wavelet.su','updated_wavelet_Iter'//siterate//'.su')
                call execute_command_line('(cd '//dir_out//' ; ln -sf updated_wavelet_Iter'//int2str(iterate)//'.su updated_wavelet_final.su )')

                tmp=split(shls%all)
                list=setup%get_ints('ITERATE_KEEP_SHOT',&
                    o_default='1 '// either( tmp(size(tmp)/2+1)%s, '', size(tmp)>1 ))

                do i=1,size(list)
                    sindex=num2str(list(i),o_format='(i0.4)')
                    call sysio_mv('dsyn_Shot'//sindex//'.su','dsyn_Iter'//siterate//'_Shot'//sindex//'.su')
                    call sysio_mv( 'RE0_Shot'//sindex//'.su', 'RE0_Iter'//siterate//'_Shot'//sindex//'.su')
                    call sysio_mv( 'RdE_Shot'//sindex//'.su', 'RdE_Iter'//siterate//'_Shot'//sindex//'.su')
                    call sysio_mv('dadj_Shot'//sindex//'.su','dadj_Iter'//siterate//'_Shot'//sindex//'.su')
                enddo

            end select

        endif
    
    end subroutine

end