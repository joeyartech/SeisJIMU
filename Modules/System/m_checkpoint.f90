module m_checkpoint
use m_string
use m_setup

    private
    public :: checkpoint_init

    logical,public :: if_checkpoint !if use checkpoint system in the 1st run
    logical,public :: if_use_checkpoint !if use checkpointed variables in the 2nd run 

    character(:),allocatable :: dir_checkpoint

    type,public :: t_checkpoint
        character(:),allocatable :: name

        logical :: if_first_init=.true.

        !counters
        integer n_cnts !number of counters
        integer,dimension(:),allocatable :: i_cnts !indices of each counter
        type(t_string),dimension(:),allocatable :: s_cnts !counters names
        type(t_string),dimension(:),allocatable :: s_cnts_cond !increment condition

        logical exist
        integer fp

        contains
        procedure :: init => init
        procedure :: count => count
        procedure :: check => check
        procedure :: open => open
        procedure :: read => read
        procedure :: write => write
        procedure :: close => close

    end type

    contains

    subroutine checkpoint_init

        if_checkpoint=setup%get_bool('IF_CHECKPOINT',o_default=checkpoint)

        if_use_checkpoint = either(&
            setup%get_bool('IF_USE_CHECKPOINT', o_default='T'),&
            .false.,&
            if_checkpoint)

        dir_checkpoint = setup%get_str('DIR_CHECKPOINT','DIR_CHP',o_default='./checkpoints/')
        
        if(if_checkpoint) call execute_command_line('mkdir -p '//dir_checkpoint, wait=.true.)

    end subroutine

    subroutine init(self,name,cnts_name,o_cnts_cond)
        class(t_checkpoint) :: self
        character(*) :: name, cnts_name
        character(*),optional :: o_cnts_cond

        type(t_string),dimension(:),allocatable :: list

        if(self%if_first_init) then
            self%name=name

            !counters name
            list=split(cnts_name); self%n_cnts=size(list)
            
            allocate(self%s_cnts(self%n_cnts))
            
            do i=1,self%n_cnts
                self%s_cnts(i)%s=list(i)%s
            enddo

            deallocate(list)

            !counters indices
            allocate(self%i_cnts(self%n_cnts),source=0)
            
            !counters increment condition
            allocate(self%s_cnts_cond(self%n_cnts))
            do i=1,self%n_cnts
                self%s_cnts_cond(i)%s='per_count' !default to increase every cycle of the loop
            enddo

            if(present(o_cnts_cond)) then
                list=split(o_cnts_cond)

                do i=1,size(list)
                    self%s_cnts_cond(i)%s=list(i)%s
                enddo

            endif

            self%if_first_init=.false.

        endif

        if(.not.self%if_first_init) then

            do i=1,self%n_cnts
                if (self%s_cnts_cond(i)%s=='per_init') then
                    self%i_cnts(i)=self%i_cnts(i)+1
                endif

            enddo

        endif

    end subroutine

    subroutine count(self,o_given)
        class(t_checkpoint) :: self
        integer,optional :: o_given

        do i=1,self%n_cnts
            select case (self%s_cnts_cond(i)%s)
            case ('per_cycle')
                self%i_cnts(i)=self%i_cnts(i)+1

            case ('given')
                self%i_cnts(i)=either(o_given,0,present(o_given))
            end select
        enddo

    end subroutine

    logical function check(self,var,o_filesize)
        class(t_checkpoint) :: self
        character(*) :: var
        integer,optional :: o_filesize
        
        character(:),allocatable :: file
        integer filesize

        file=self%name//'_'//&
            strcat(          self%s_cnts,             self%n_cnts,o_glue=',')//'_'//&
            strcat(nums2strs(self%i_cnts,self%n_cnts),self%n_cnts,o_glue=',')//'_'//&
            strcat(          self%s_cnts_cond,        self%n_cnts,o_glue=',')//'_'//&
            'Var:'//var

        inquire(file=dir_checkpoint//file,exist=check,size=filesize)
        
        !requires to match file size
        if(check .and. present(oif_filesize)) then
            filesize=filesize/4
            if(filesize<oif_filesize) then
                call hud('At checkpoint '//self%name//','//&
                    ' the size of the requested file ('//num2str(filesize)//') is smaller than the prescribed size ('//(num2str(o_filesize))//')'//s_return//&
                    'Will recompute the associated variables.')
                check=.false.
            endif
        endif

    end function
    
    subroutine open(self,var)
        class(t_checkpoint) :: self
        character(*) :: var

        character(:),allocatable :: file

        file=self%name//'_'//&
            strcat(          self%s_cnts,             self%n_cnts,o_glue=',')//'_'//&
            strcat(nums2strs(self%i_cnts,self%n_cnts),self%n_cnts,o_glue=',')//'_'//&
            strcat(          self%s_cnts_cond,        self%n_cnts,o_glue=',')//'_'//&
            'Var:'//var

        open(11,file=dir_checkpoint//file,access='stream')
        self%fp=11

    end subroutine

    subroutine read(self,arr,size)
        class(t_checkpoint) :: self
        integer size
        real :: arr(size)
        read(self%fp) arr
    end subroutine
    
    subroutine write(self,arr,size)
        class(t_checkpoint) :: self
        integer size
        real :: arr(size)
        write(self%fp) arr
    end subroutine

    subroutine close(self)
        class(t_checkpoint) :: self
        close(self%fp)
    end subroutine

end module