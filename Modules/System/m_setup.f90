module m_setup
use m_message
use m_string

    private

    character(:),allocatable :: file

    type,public :: t_setup
        logical :: exist=.true.

        contains
        procedure :: init => init
        procedure,nopass :: check => check
        procedure :: get_int   => get_int
        procedure :: get_ints  => get_ints
        procedure :: get_real  => get_real
        procedure :: get_reals => get_reals
        procedure :: get_str   => get_str
        procedure :: get_strs  => get_strs
        procedure :: get_file  => get_file
        !procedure :: get_files => get_files
        procedure :: get_bool  => get_bool
        !procedure :: get_bools => get_bools

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
        
    end subroutine

    function find(key,o_alias) result(res)
        character(*) :: key
        character(*),optional :: o_alias
        character(:),allocatable :: res

        character(i_str_len) :: text
        character(i_str_len+4) :: text2
        character(i_str_slen) :: tmp_key, tmp_val, tmp_res
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

    function check(key,o_alias) result(exist)
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        logical :: exist

        exist=.false.

        if(find(key,o_alias)/='') exist=.true.
        
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
                call hud(keys//"is NOT found, take 0, '' or F for number, string/filename or boolean types.")
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

        type(t_string),dimension(:),allocatable :: tmp

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

        res=read(key,o_alias,o_default,o_mandatory)

        if(res/='') then
            inquire(file=res,exist=exist)
            if (.not. exist) then
                call hud('File '//res//' does NOT exist, return empty filename.')
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

end