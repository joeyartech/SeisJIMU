module m_checkpoint
use m_string
use m_setup

    private
    public :: checkpoint_init

    logical :: if_checkpoint !if use checkpoint system in the 1st run
    logical :: if_use_checkpoint !if use checkpointed variables in the 2nd run 

    character(:),allocatable :: dir_chp

    type,public :: t_checkpoint
        character(:),allocatable :: name

        logical :: if_first_init=.true.
        logical :: if_fuse

        !counters
        integer n_cnts !number of counters
        integer,dimension(:),allocatable :: i_cnts !indices of each counter
        type(t_string),dimension(:),allocatable :: s_cnts !counters names
        type(t_string),dimension(:),allocatable :: s_cnts_cond !increment condition

        logical exist
        integer fp

        contains
        procedure :: init
        procedure :: count
        procedure :: check
        procedure :: open
        procedure :: read_real1,  read_real2,  read_real3,  read_real4
        procedure :: write_real1, write_real2, write_real3, write_real4
        generic :: read  =>  read_real1,  read_real2,  read_real3,  read_real4
        generic :: write => write_real1, write_real2, write_real3, write_real4
        procedure :: close

    end type

    contains

    subroutine checkpoint_init

        if_checkpoint=setup%get_bool('IF_CHECKPOINT',o_default='F')

        if_use_checkpoint = either(&
            setup%get_bool('IF_USE_CHECKPOINT', o_default='T'),&
            .false.,&
            if_checkpoint)

        dir_chp = setup%get_str('DIR_CHECKPOINT','DIR_CHP',o_default='./checkpoints/')
        
        if(if_checkpoint) call execute_command_line('mkdir -p '//dir_chp, wait=.true.)

    end subroutine

    subroutine init(self,name,o_cnts_name,o_cnts_cond,oif_fuse)
        class(t_checkpoint) :: self
        character(*) :: name
        character(*),optional :: o_cnts_name, o_cnts_cond
        logical,optional :: oif_fuse

        type(t_string),dimension(:),allocatable :: list

        if(.not.if_checkpoint) return

        if(self%if_first_init) then
            self%if_first_init=.false.

            self%name=name

            self%if_fuse=either(oif_fuse,.false.,present(oif_fuse))

            !no counters
            if(.not.present(o_cnts_name)) then
                self%n_cnts=0
                return
            endif

            !counters name
            list=split(o_cnts_name); self%n_cnts=size(list)
            
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

            return

        endif

        if(.not.self%if_first_init) then

            do i=1,self%n_cnts
                if (self%s_cnts_cond(i)%s=='per_init') then
                    self%i_cnts(i)=self%i_cnts(i)+1
                endif

            enddo

            return

        endif

    end subroutine

    subroutine count(self,o_given)
        class(t_checkpoint) :: self
        integer,optional :: o_given

        if(.not.if_checkpoint) return

        do i=1,self%n_cnts
            select case (self%s_cnts_cond(i)%s)
            case ('per_cycle')
                self%i_cnts(i)=self%i_cnts(i)+1

            case ('given')
                self%i_cnts(i)=either(o_given,0,present(o_given))
            end select
        enddo

    end subroutine

    logical function check(self,var)
        class(t_checkpoint) :: self
        character(*) :: var
        
        character(:),allocatable :: file
        integer filesize

        if(.not. if_use_checkpoint) then
            check=.false.
            return
        endif

        file=self%name//'_'//&
            strcat(          self%s_cnts,     o_glue=',')//'_'//&
            strcat(nums2strs(self%i_cnts),    o_glue=',')//'_'//&
            strcat(          self%s_cnts_cond,o_glue=',')//'_'//&
            'Var:'//var

        inquire(file=dir_chp//file,exist=check) !,size=filesize)
        
        !will fuse the whole checkpoint system when check is false
        if(self%if_fuse) then
            if(check.eqv..false.) if_use_checkpoint=.false.
        endif

        ! !requires to match file size
        !but is difficult to check
        ! if(check .and. present(oif_filesize)) then
        !     filesize=filesize/4
        !     if(filesize<oif_filesize) then
        !         call hud('At checkpoint '//self%name//','//&
        !             ' the size of the requested file ('//num2str(filesize)//') is smaller than the prescribed size ('//(num2str(o_filesize))//')'//s_NL//&
        !             'Will recompute the associated variables.')
        !         check=.false.
        !     endif
        ! endif

    end function
    
    subroutine open(self,var)
        class(t_checkpoint) :: self
        character(*) :: var

        character(:),allocatable :: file

        if(if_use_checkpoint.or.if_checkpoint) then

            file=self%name//'_'//&
                strcat(          self%s_cnts,     o_glue=',')//'_'//&
                strcat(nums2strs(self%i_cnts),    o_glue=',')//'_'//&
                strcat(          self%s_cnts_cond,o_glue=',')//'_'//&
                'Var:'//var

            open(11,file=dir_chp//file,access='stream')
            self%fp=11

        endif

    end subroutine

    subroutine read_real1(self,a,b,c,d,e,f,g)
        class(t_checkpoint) :: self
        real,dimension(:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g

        if(if_use_checkpoint) then
                if(allocated(a)) read(self%fp) a
            if(present(b)) then
                if(allocated(b)) read(self%fp) b
            endif
            if(present(c)) then
                if(allocated(c)) read(self%fp) c
            endif
            if(present(d)) then
                if(allocated(d)) read(self%fp) d
            endif
            if(present(e)) then
                if(allocated(e)) read(self%fp) e
            endif
            if(present(f)) then
                if(allocated(f)) read(self%fp) f
            endif
            if(present(g)) then
                if(allocated(g)) read(self%fp) g
            endif

        endif

    end subroutine

    subroutine read_real2(self,a,b,c,d,e,f,g)
        class(t_checkpoint) :: self
        real,dimension(:,:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g

        if(if_use_checkpoint) then
                if(allocated(a)) read(self%fp) a
            if(present(b)) then
                if(allocated(b)) read(self%fp) b
            endif
            if(present(c)) then
                if(allocated(c)) read(self%fp) c
            endif
            if(present(d)) then
                if(allocated(d)) read(self%fp) d
            endif
            if(present(e)) then
                if(allocated(e)) read(self%fp) e
            endif
            if(present(f)) then
                if(allocated(f)) read(self%fp) f
            endif
            if(present(g)) then
                if(allocated(g)) read(self%fp) g
            endif

        endif

    end subroutine

    subroutine read_real3(self,a,b,c,d,e,f,g)
        class(t_checkpoint) :: self
        real,dimension(:,:,:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g

        if(if_use_checkpoint) then
                if(allocated(a)) read(self%fp) a
            if(present(b)) then
                if(allocated(b)) read(self%fp) b
            endif
            if(present(c)) then
                if(allocated(c)) read(self%fp) c
            endif
            if(present(d)) then
                if(allocated(d)) read(self%fp) d
            endif
            if(present(e)) then
                if(allocated(e)) read(self%fp) e
            endif
            if(present(f)) then
                if(allocated(f)) read(self%fp) f
            endif
            if(present(g)) then
                if(allocated(g)) read(self%fp) g
            endif

        endif

    end subroutine

    subroutine read_real4(self,a,b,c,d,e,f,g)
        class(t_checkpoint) :: self
        real,dimension(:,:,:,:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g

        if(if_use_checkpoint) then
                if(allocated(a)) read(self%fp) a
            if(present(b)) then
                if(allocated(b)) read(self%fp) b
            endif
            if(present(c)) then
                if(allocated(c)) read(self%fp) c
            endif
            if(present(d)) then
                if(allocated(d)) read(self%fp) d
            endif
            if(present(e)) then
                if(allocated(e)) read(self%fp) e
            endif
            if(present(f)) then
                if(allocated(f)) read(self%fp) f
            endif
            if(present(g)) then
                if(allocated(g)) read(self%fp) g
            endif

        endif

    end subroutine

    subroutine write_real1(self,a,b,c,d,e,f,g)
        class(t_checkpoint) :: self
        real,dimension(:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g

        if(if_use_checkpoint) then
                if(allocated(a)) write(self%fp) a
            if(present(b)) then
                if(allocated(b)) write(self%fp) b
            endif
            if(present(c)) then
                if(allocated(c)) write(self%fp) c
            endif
            if(present(d)) then
                if(allocated(d)) write(self%fp) d
            endif
            if(present(e)) then
                if(allocated(e)) write(self%fp) e
            endif
            if(present(f)) then
                if(allocated(f)) write(self%fp) f
            endif
            if(present(g)) then
                if(allocated(g)) write(self%fp) g
            endif

        endif

    end subroutine

    subroutine write_real2(self,a,b,c,d,e,f,g)
        class(t_checkpoint) :: self
        real,dimension(:,:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g

        if(if_use_checkpoint) then
                if(allocated(a)) write(self%fp) a
            if(present(b)) then
                if(allocated(b)) write(self%fp) b
            endif
            if(present(c)) then
                if(allocated(c)) write(self%fp) c
            endif
            if(present(d)) then
                if(allocated(d)) write(self%fp) d
            endif
            if(present(e)) then
                if(allocated(e)) write(self%fp) e
            endif
            if(present(f)) then
                if(allocated(f)) write(self%fp) f
            endif
            if(present(g)) then
                if(allocated(g)) write(self%fp) g
            endif

        endif

    end subroutine

    subroutine write_real3(self,a,b,c,d,e,f,g)
        class(t_checkpoint) :: self
        real,dimension(:,:,:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g

        if(if_use_checkpoint) then
                if(allocated(a)) write(self%fp) a
            if(present(b)) then
                if(allocated(b)) write(self%fp) b
            endif
            if(present(c)) then
                if(allocated(c)) write(self%fp) c
            endif
            if(present(d)) then
                if(allocated(d)) write(self%fp) d
            endif
            if(present(e)) then
                if(allocated(e)) write(self%fp) e
            endif
            if(present(f)) then
                if(allocated(f)) write(self%fp) f
            endif
            if(present(g)) then
                if(allocated(g)) write(self%fp) g
            endif

        endif

    end subroutine

    subroutine write_real4(self,a,b,c,d,e,f,g)
        class(t_checkpoint) :: self
        real,dimension(:,:,:,:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g

        if(if_use_checkpoint) then
                if(allocated(a)) write(self%fp) a
            if(present(b)) then
                if(allocated(b)) write(self%fp) b
            endif
            if(present(c)) then
                if(allocated(c)) write(self%fp) c
            endif
            if(present(d)) then
                if(allocated(d)) write(self%fp) d
            endif
            if(present(e)) then
                if(allocated(e)) write(self%fp) e
            endif
            if(present(f)) then
                if(allocated(f)) write(self%fp) f
            endif
            if(present(g)) then
                if(allocated(g)) write(self%fp) g
            endif

        endif

    end subroutine

    subroutine close(self)
        class(t_checkpoint) :: self

        if(if_use_checkpoint.or.if_checkpoint) close(self%fp)
        
    end subroutine

end module