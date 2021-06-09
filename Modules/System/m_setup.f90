module m_setup
use m_global
use m_message
use m_string

    private

    character(:),allocatable :: file

    type t_setup
        logical :: exist=.true.
        character(:),allocatable :: dir_working
        character(:),allocatable :: dir_scratch

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
        else
            inquire(file=file, exist=self%exist)
            if(.not.self%exist) then
                call hud('Setup file '//file//' does NOT exist!')
                self%exist=.false.
                return
            endif
            
            call hud('Setup file: '//file)

            self%dir_working=self%get_str('DIR','DIR_WORKING',o_default='./')
            self%dir_scratch=self%get_str('DIR_SCRATCH','DIR_TMP',o_default=self%dir_working)

            return

        endif

    end subroutine

    function find(key,alias) result(res)
        character(*) :: key,alias
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
            
            read(text,*)  tmp_key, tmp_val
            if(tmp_val=='!!') cycle !missing value

            if(key==lalign(tmp_key) .or. alias==lalign(tmp_key)) then
                tmp_res=tmp_val
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

    function check(key,alias) result(exist)
        character(*) :: key
        character(*),optional :: alias !alias of inquired key
        logical :: exist

        exist=.false.

        if(present(alias)) then
            if(find(key,alias)/='') exist=.true.
        else
            if(find(key,key)/='') exist=.true.
        endif
        
    end function

    function read(key,alias,default,mandatory) result(res)
        character(*) :: key,alias
        character(*) :: default !default output value
        logical :: mandatory
        character(:),allocatable :: res

        character(:),allocatable :: keys
        keys=key//' ('//alias//') '

        res=find(key,alias)

        if(res/='') then
            call hud(keys//':'//res)
        endif

        if(res=='') then
            if(default/='') then
                call hud(keys//'is NOT found, take default: '//default)
                res=default
            elseif(mandatory) then
                call error(keys//'is NOT found, but is MANDATORY.')
                ! call hud(keys//'is NOT found, but is MANDATORY.'//s_return &
                !     'SeisJIMU has to ask for the value of this key before running further.'//s_return &
                !     'Please open the text file "tentative" in the working directory,'//s_return &
                !     'provide the missing value, remove "!!" in the line,'//s_return &
                !     '(like in the setup file)'//s_return &
                !     'and close the file.'//s_return &
                !     'SeisJIMU checks this file every 1 min'//s_return &
                !     'and will take the new value when no "!!" is leading the line.')
                ! res=demand(key)
            else
                call hud(keys//"is NOT found, take 0 for number(s), '' for character(s) and filename, or .false. for logical type(s)")
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
        logical,optional :: o_mandatory       
        integer :: res

        character(:),allocatable :: alias, default
        logical :: mandatory

        if(present(o_alias)) then
            alias=o_alias
        else
            alias=key
        endif

        if(present(o_default)) then
            default=o_default
        else
            default=''
        endif

        if(present(o_mandatory)) then
            mandatory=o_mandatory
        else
            mandatory=.false.
        endif

        res=str2int(read(key,alias,default,mandatory))

        deallocate(alias, default)

    end function

    function get_ints(self,key,o_alias,o_default,o_sep,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        character(1),optional :: o_sep
        logical,optional :: o_mandatory
        integer,dimension(:),allocatable :: res

        character(:),allocatable :: alias, default
        character(1) :: sep
        logical :: mandatory
        type(t_string),dimension(:),allocatable :: tmp

        if(present(o_alias)) then
            alias=o_alias
        else
            alias=key
        endif

        if(present(o_sep)) then
            sep=o_sep
        else
            sep=' '
        endif

        if(present(o_default)) then
            default=o_default
        else
            default=''
        endif

        if(present(o_mandatory)) then
            mandatory=o_mandatory
        else
            mandatory=.false.
        endif

        tmp=split(read(key,alias,default,mandatory),o_sep=sep)
        res=strs2ints(tmp,size(tmp))

        deallocate(alias, default, tmp)

    end function

    function get_real(self,key,o_alias,o_default,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        logical,optional :: o_mandatory
        real :: res

        character(:),allocatable :: alias, default
        logical :: mandatory

        if(present(o_alias)) then
            alias=o_alias
        else
            alias=key
        endif

        if(present(o_default)) then
            default=o_default
        else
            default=''
        endif

        if(present(o_mandatory)) then
            mandatory=o_mandatory
        else
            mandatory=.false.
        endif

        res=str2real(read(key,alias,default,mandatory))

        deallocate(alias, default)

    end function

    function get_reals(self,key,o_alias,o_default,o_sep,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value        
        character(1),optional :: o_sep
        logical,optional :: o_mandatory
        real,dimension(:),allocatable :: res
        
        character(:),allocatable :: alias, default
        character(1) :: sep
        logical :: mandatory
        type(t_string),dimension(:),allocatable :: tmp

        if(present(o_alias)) then
            alias=o_alias
        else
            alias=key
        endif

        if(present(o_sep)) then
            sep=o_sep
        else
            sep=' '
        endif

        if(present(o_default)) then
            default=o_default
        else
            default=''
        endif

        if(present(o_mandatory)) then
            mandatory=o_mandatory
        else
            mandatory=.false.
        endif

        tmp=split(read(key,alias,default,mandatory),o_sep=sep)
        res=strs2reals(tmp,size(tmp))

        deallocate(alias, default, tmp)

    end function

    function get_str(self,key,o_alias,o_default,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key       
        character(*),optional :: o_default !default output value        
        logical,optional :: o_mandatory
        character(:),allocatable :: res

        character(:),allocatable :: alias, default
        logical :: mandatory

        if(present(o_alias)) then
            alias=o_alias
        else
            alias=key
        endif

        if(present(o_default)) then
            default=o_default
        else
            default=''
        endif

        if(present(o_mandatory)) then
            mandatory=o_mandatory
        else
            mandatory=.false.
        endif

        res=read(key,alias,default,mandatory)

        deallocate(alias, default)

    end function
    
    function get_strs(self,key,o_alias,o_default,o_sep,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        character(1),optional :: o_sep
        logical,optional :: o_mandatory
        type(t_string),dimension(:),allocatable :: res

        character(:),allocatable :: alias, default
        character(1) :: sep
        logical :: mandatory

        if(present(o_alias)) then
            alias=o_alias
        else
            alias=key
        endif

        if(present(o_sep)) then
            sep=o_sep
        else
            sep=' '
        endif

        if(present(o_default)) then
            default=o_default
        else
            default=''
        endif

        if(present(o_mandatory)) then
            mandatory=o_mandatory
        else
            mandatory=.false.
        endif

        res=split(read(key,alias,default,mandatory),o_sep=sep)

        deallocate(alias, default)

    end function

    function get_file(self,key,o_alias,o_default,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        logical,optional :: o_mandatory
        character(:),allocatable :: res

        character(:),allocatable :: alias, default
        logical :: mandatory, exist

        if(present(o_alias)) then
            alias=o_alias
        else
            alias=key
        endif

        if(present(o_default)) then
            default=o_default
        else
            default=''
        endif

        if(present(o_mandatory)) then
            mandatory=o_mandatory
        else
            mandatory=.false.
        endif

        res=read(key,alias,default,mandatory)

        if(res/='') then
            inquire(file=res,exist=exist)
            if (.not. exist) then
                call hud('File '//res//' does NOT exist, return empty filename.')
                res=''
            endif
        endif

        deallocate(alias, default)
        
    end function
        
    function get_bool(self,key,o_alias,o_default,o_mandatory) result(res)
        class(t_setup) :: self
        character(*) :: key
        character(*),optional :: o_alias !alias of inquired key
        character(*),optional :: o_default !default output value
        logical,optional :: o_mandatory
        logical :: res

        character(:),allocatable :: alias, default
        logical :: mandatory

        if(present(o_alias)) then
            alias=o_alias
        else
            alias=key
        endif

        if(present(o_default)) then
            default=o_default
        else
            default='F'
        endif

        if(present(o_mandatory)) then
            mandatory=o_mandatory
        else
            mandatory=.false.
        endif

        res=str2bool(read(key,alias,default,mandatory))

        deallocate(alias, default)

    end function

end