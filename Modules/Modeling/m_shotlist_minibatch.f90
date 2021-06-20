module m_shotlist
use m_string
use m_mpienv
use m_message
use m_arrayop
use m_setup
use m_sysio

    private

    type,public :: t_shotlist
        integer,dimension(:),allocatable :: list, list_per_processor
        integer :: nshot, nshot_per_processor

        contains
        procedure :: read_from_setup => read_from_setup
        procedure :: read_from_data  => read_from_data
        procedure :: assign => assign
        procedure :: assign_random => assign_random
        procedure :: print => print
    end type

    type(t_shotlist),public :: sl
    
    contains

    subroutine read_from_setup(self)
        class(t_shotlist) :: self

        select case (setup%get_str('ACQUI_GEOMETRY',o_default='spread'))

        case ('spread_irregular')
            open(13,file=setup%get_file('FILE_SOURCE_POSITION','SPOS',o_mandatory=.true.),action='read')
            !count number of sources
            self%nshot=0
            do
                read (13,*,iostat=msg) z,x,y
                if(msg/=0) exit
                self%nshot=self%nshot+1
            enddo
            close(13)
            call hud('Will read '//num2str(self%nshot)//' source positions.')

        case default
            self%nshot=setup%get_int('NUMBER_SOURCE','NS',o_mandatory=.true.)
            
        end select        

        call alloc(self%list,self%nshot); self%list=[(i,i=1,self%nshot)]

        call hud('Will compute '//num2str(self%nshot)//' synthetic shots.')

    end subroutine
    
    subroutine read_from_data(self)
        class(t_shotlist) :: self

        type(t_string),dimension(:),allocatable :: shots, subshots
        character(:),allocatable :: file
        logical :: exist
        integer file_size
        
        shots=setup%get_strs('SHOT_INDEX','ISHOT')
        
        if(size(shots)==0) then !not given
            call hud('SHOT_INDEX is not given. Now count how many data files exist in the directory..')
            file=setup%get_file('FILE_DATA')

            i=0; exist=.true.
            do
                i=i+1
                inquire(file=file//num2str(i,'(i0.4)')//'.su', size=file_size, exist=exist)
                if(file_size==0) exist=.false.
                if(.not.exist) exit
            enddo

            self%nshot=i

            call alloc(self%list,self%nshot); self%list=[(i,i=1,self%nshot)]
        
            call hud('Found '//num2str(self%nshot)//' sequential shots.')

        endif

        if(size(shots)>0) then !given
            
            allocate(self%list(1))
            
            do i=1,size(shots)

                subshots=split(shots(i)%s,o_sep=':')

                if(size(subshots)==1) then !add/rm
                    if(subshots(1)%s(1:1)=='-'.or.subshots(1)%s(1:1)=='!') then !rm
                        call rm(self%list,size(self%list),str2int(subshots(1)%s(2:)))
                    else
                        call add(self%list,size(self%list),str2int(subshots(1)%s))
                    endif
                endif
                if(size(subshots)==2) then !first:last
                    do j=str2int(subshots(1)%s),str2int(subshots(2)%s)
                        call add(self%list,size(self%list),j)
                    enddo
                endif
                if(size(subshots)==3) then !first:increment:last
                    do j=str2int(subshots(1)%s),str2int(subshots(3)%s),str2int(subshots(2)%s)
                        call add(self%list,size(self%list),j)
                    enddo
                endif

            enddo

            self%nshot=size(self%list)

            call hud('Total number of shots: '//num2str(self%nshot)//'')
            
        endif

        deallocate(shots, subshots)

        !Random decimation of shots:
        !First, shots are grouped into N groups. In each group, one shot is randomly selected to compute the objective function and gradient.
        !When the optimizer signals a new selection, new shots are randomly selected from each group
        !
        !batch size = number of shots to compute one gradient = number of groups
        !batch size=1: stochastic gradient descent
        !1<batch size<all shots: mini-batch gradient descent
        !batch size=number of total shots: batch gradient descent (ie no random selection)
        
        batch_size=setup%get_str('SHOT_BATCH_SIZE','NBATCH',o_default=num2str(self%nshot)) !no random selection by default

        integer,dimension(:),allocatable :: list

        call hud('No. of processors: '//num2str(mpiworld%nproc)//s_return// &
                'No. of shots: '//num2str(self%nshot))

        if(self%nshot<mpiworld%nproc) call warn('Shot number < Processor number. Some processors will be always idle..')

        call alloc(list, ceiling(self%nshot*1./mpiworld%nproc) )

!deprecated.. can't do the job..
k=1 !shotlist's index
j=0 !modulo shot#/processor#

do i=1,self%nshot
    if (j==mpiworld%nproc) then
        j=j-mpiworld%nproc
        k=k+1
    endif
    if (mpiworld%iproc==j) then
        list(k)=self%list(i)
    endif
    j=j+1
enddo

self%nshot_per_processor=count(list>=0)

call alloc(self%list_per_processor,self%nshot_per_processor)

self%list_per_processor=list(1:self%nshot_per_processor)

deallocate(list)
        
    end subroutine

    subroutine assign(self)
        class(t_shotlist) :: self

        type(t_checkpoint) :: chp

        call chp%init('FWI_randomshotlist','shotlist')

        if(.not.chp%check('shotlist',o_filesize=?) then
            call hud('The random shotlist file is considered invalid as checkpoint '//chp%name//' returns F when checking this file. '//s_return//&
                'For reproducibility arguments, set if_use_checkpoint to F.')
            if_use_checkpoint=.false.
        endif

    end subroutine

    subroutine print(self)
        class(t_shotlist) :: self
                
        !message
        write(*,*) ' Proc# '//mpiworld%sproc//' has '//num2str(self%nshot_per_processor)//' assigned shots.'
        call hud('See file "shotlist" for details.')

        !write shotlist to disk
        call mpiworld%write(dir_out//'shotlist', 'Proc# '//mpiworld%sproc//' has '//num2str(self%nshot_per_processor)//' assigned shots:'// &
            strcat(ints2strs(self%list_per_processor,self%nshot_per_processor),self%nshot_per_processor))
        
    end subroutine

end
