module m_setup
use m_string
use m_message

    private

    character(:),allocatable :: file

    character(:),allocatable :: previous

    type,public :: t_setup
        logical :: exist=.true.
        character(:),allocatable :: dir_in

        contains
        procedure :: init
        procedure :: check
        procedure :: get_int
        procedure :: get_ints
        procedure :: get_real
        procedure :: get_reals
        procedure :: get_str
        procedure :: get_strs
        procedure :: get_file
        !procedure :: get_files
        procedure :: get_bool
        !procedure :: get_bools

        procedure :: prev_int
        procedure :: prev_ints
        procedure :: prev_real
        procedure :: prev_reals
        procedure :: prev_str
        procedure :: prev_strs
        procedure :: prev_file
        procedure :: prev_bool

    end type
    
    type(t_setup),public :: setup

    contains
    
    subroutine init(self)
        class(t_setup) :: self

        character(i_str_xlen) :: tmp
        
        call getarg(1,tmp)
        file=lalign(tmp)
        
        if(file=='') then
            self%exist=.false.
            return
        endif

        inquire(file=file, exist=self%exist)
        if(.not.self%exist) then
            call hud('Setup file '//file//' does NOT exist!')
            self%exist=.false.
            return
        endif
        
        call hud('Setup file: '//file)

        self%dir_in='./'

        previous=''
        
    end subroutine

    function find(key,o_alias) result(res)
        character(*) :: key
        character(*),optional :: o_alias
        character(:),allocatable :: res

        !character(i_str_len) :: text
        !character(i_str_len+4) :: text2
        !character(i_str_slen) :: tmp_key, tmp_val, tmp_res
        character(i_str_xlen) :: text
        character(i_str_xlen+4) :: text2
        character(i_str_xlen) :: tmp_key, tmp_val, tmp_res
        tmp_res=''

        open(10,file=file,action='read')
        do
            read(10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            text2=text//'  !!'
            
            read(text2,*)  tmp_key, tmp_val
            if(tmp_val=='!!') cycle !missing value

            if(key==tmp_key) tmp_res=tmp_val

            if(present(o_alias)) then
                if(o_alias==tmp_key) tmp_res=tmp_val
            endif

        end do
        close(10)

        res=lalign(tmp_res)

        !update variable previous
        if(allocated(previous)) deallocate(previous)
        previous=res

    end function

    ! function demand(key) result(res)
    !     character(*) :: key
    !     character(:),allocatable :: res

    !     !only master processor can do this job
    !     if(mpiworld%is_master) then
        
    !         !call email()

    !         open(12,file='tentative')
    !         write(12,'(a)') '!!  '//key//'    '
    !         close(12)

    !         do
    !             sleep(60)
    !             inquire(file='tentative', exist=exist)
    !             if(.not.exist) then
    !                 open(12,file='tentative')
    !                 write(12,'(a)') '!!  '//key//'    '
    !                 close(12)
    !             endif

    !             open(12,file='tentative')
    !             read(12,"(a)",iostat=msg) text ! read line into character variable
    !             close(12)
    !             if (msg < 0) stop 'tentative file is empty. Something is wrong.'
    !             if (msg > 0) stop 'Check tentative file.  Something is wrong.'
                
    !             if (text(1:2)=='!!') then
    !                 cycle !continue to wait for user input
    !             else
    !                 read(text,*) tmp_val
    !                 res=lalign(tmp_val)
    !                 exit
    !             endif
    !         enddo

    !         !delete tentative file
    !         open(12,file='tentative',status='old')
    !         close(12,status='delete')

    !         !append to setup file
    !         open(12,file=file,status='unknown',access='append')
    !         write(12,'(a)') key//'  '//res
    !         close(12)

    !     endif

    ! end function

    function check(self,key,o_alias) result(exist)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        logical :: exist

        if(find(key,o_alias)/='') then
            call hud(key//' is found in '//file)
            exist=.true.
        else
            call hud(key//' is NOT found in '//file)
            exist=.false.
        endif
        
    end function

    function read(key,o_alias,o_default,o_mandatory) result(res)
        character(*) :: key
        character(*),optional :: o_alias, o_default
        integer,optional :: o_mandatory
        character(:),allocatable :: res

        character(:),allocatable :: keys
        integer :: mandatory
        
        keys=key//either(' ('//o_alias//') ',' ',present(o_alias))

        mandatory=either(o_mandatory,0,present(o_mandatory))

        res=find(key,o_alias)

        if(res/='') call hud(keys//': '//res)

        if(res=='') then
            if(present(o_default)) then
                call hud(keys//'is NOT found, take default: '//o_default)
                deallocate(res)
                res=o_default
            endif

            if(.not.present(o_default) .and. mandatory>0) then
                call error(keys//'is NOT found, but is MANDATORY.')
                ! call hud(keys//'is NOT found, but is MANDATORY.'//s_NL &
                !     'SeisJIMU has to ask for the value of this key before running further.'//s_NL &
                !     'Please open the text file "tentative" in the working directory,'//s_NL &
                !     'provide the missing value, remove "!!" in the line,'//s_NL &
                !     '(like in the setup file)'//s_NL &
                !     'and close the file.'//s_NL &
                !     'SeisJIMU checks this file every 1 min'//s_NL &
                !     'and will take the new value when no "!!" is leading the line.')
                ! res=demand(key)
            endif

            if(.not.present(o_default) .and. mandatory==0) then
                call hud(keys//"is NOT found, take 0, '' or F")
                res=''
            endif
        endif

        deallocate(keys)

    end function

    function get_int(self,key,o_alias,o_default,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        integer,optional :: o_mandatory       
        integer :: res

        res=str2int(read(key,o_alias,o_default,o_mandatory))

    end function

    function get_ints(self,key,o_alias,o_default,o_sep,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        character(1),optional :: o_sep
        integer,optional :: o_mandatory
        integer,dimension(:),allocatable :: res

        type(t_string),dimension(:),allocatable :: tmp, tmp2

        tmp=split(read(key,o_alias,o_default,o_mandatory),o_sep=o_sep)

        if(present(o_mandatory)) then
            if(size(tmp)<o_mandatory) then
                call error('NOT enough inputs for '//key//either(' ('//o_alias//') ',' ',present(o_alias))//', which requires >= '//num2str(o_mandatory)//' inputs.')
            endif
        endif

        res=strs2ints(tmp)

        deallocate(tmp)

    end function

    function get_real(self,key,o_alias,o_default,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        integer,optional :: o_mandatory
        real :: res

        res=str2real(read(key,o_alias,o_default,o_mandatory))

    end function

    function get_reals(self,key,o_alias,o_default,o_sep,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value        
        character(1),optional :: o_sep
        integer,optional :: o_mandatory
        real,dimension(:),allocatable :: res
        
        type(t_string),dimension(:),allocatable :: tmp

        tmp=split(read(key,o_alias,o_default,o_mandatory),o_sep=o_sep)

        if(present(o_mandatory)) then
            if(size(tmp)<o_mandatory) then
                call error('NOT enough inputs for '//key//either(' ('//o_alias//') ',' ',present(o_alias))//', which requires >= '//num2str(o_mandatory)//' inputs.')
            endif
        endif

        res=strs2reals(tmp)

        deallocate(tmp)

    end function

    function get_str(self,key,o_alias,o_default,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key       
        character(*),optional :: o_default !default output value        
        integer,optional :: o_mandatory
        character(:),allocatable :: res

        res=read(key,o_alias,o_default,o_mandatory)

    end function
    
    function get_strs(self,key,o_alias,o_default,o_sep,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        character(1),optional :: o_sep
        integer,optional :: o_mandatory
        type(t_string),dimension(:),allocatable :: res

        res=split(read(key,o_alias,o_default,o_mandatory),o_sep=o_sep)

        if(present(o_mandatory)) then
            if(size(res)<o_mandatory) then
                call error('NOT enough inputs for '//key//either(' ('//o_alias//') ',' ',present(o_alias))//', which requires >= '//num2str(o_mandatory)//' inputs.')
            endif
        endif

    end function

    function get_file(self,key,o_alias,o_default,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        integer,optional :: o_mandatory
        character(:),allocatable :: res

        logical :: exist
        integer :: size

        res=read(key,o_alias,o_default,o_mandatory)

        if(res/='') then
            inquire(file=self%dir_in//res,exist=exist,size=size)
            if(size==0) exist=.false.
            if (.not. exist) then
                call hud('File '//res//' has 0 size or does NOT exist, return empty filename.')
                res=''
            endif
        endif

    end function
        
    function get_bool(self,key,o_alias,o_default,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        integer,optional :: o_mandatory
        logical :: res

        res=str2bool(read(key,o_alias,o_default,o_mandatory))

    end function



    function prev_int(self) result(res)
        class(t_setup) :: self
        integer :: res
        
        res=str2int(previous)

    end function

    function prev_ints(self,o_sep) result(res)
        class(t_setup) :: self
        character(1),optional :: o_sep
        integer,dimension(:),allocatable :: res
        
        res=strs2ints(split(previous,o_sep=o_sep))

    end function

    function prev_real(self) result(res)
        class(t_setup) :: self
        real :: res

        res=str2real(previous)

    end function

    function prev_reals(self,o_sep) result(res)
        class(t_setup) :: self
        character(1),optional :: o_sep
        real,dimension(:),allocatable :: res

        res=strs2reals(split(previous,o_sep=o_sep))

    end function

    function prev_str(self) result(res)
        class(t_setup) :: self
        character(:),allocatable :: res

        res=previous

    end function

    function prev_strs(self,o_sep) result(res)
        class(t_setup) :: self
        character(1),optional :: o_sep
        type(t_string),dimension(:),allocatable :: res

        res=split(previous,o_sep=o_sep)

    end function

    function prev_file(self) result(res)
        class(t_setup) :: self
        character(:),allocatable :: res

        res=previous

    end function

    function prev_bool(self) result(res)
        class(t_setup) :: self
        logical :: res

        res=str2bool(previous)
        
    end function

end