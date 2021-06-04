module m_setup
use m_global
use m_message
use m_string

    type t_setup
        logical :: exist=.true.
        character(:),allocatable :: file
        character(:),allocatable :: dir_working
        character(:),allocatable :: dir_scratch

        contains
        procedure :: init => init
        procedure,nopass :: check => check

        procedure,nopass :: get => get

    end type
    
    type(t_setup) :: setup


    contains
    
    subroutine init(self)
        type(t_setup) :: self

        character(i_str_xlen) :: tmp
        
        call getarg(1,tmp)
        self%file=trim(adjustl(tmp))
        
        if(self%file=='') then
            self%exist=.false.
            return
        else
            inquire(file=self%file, exist=exist)
            if(.not.exist) then
                call hud('Setup file '//self%file//' does NOT exist!')
                self%exist=.false.
                return
            endif
            
            call hud('Setup file: '//self%file)

            self%dir_working=self%read_str('DIR','DIR_WORKING',default='./')
            self%dir_scratch=self%read_str('DIR_SCRATCH','DIR_TMP',default=dir_working)

            return
        endif

    end subroutine


    function find(key,key2)
        character(*) :: key,key2
        character(*) :: find

        character(i_str_len) :: text
        character(i_str_len+4) :: text2
        character(i_str_slen) :: tmp_key, tmp_val, tmp_find
        tmp_find=''

        open(10,file=file_setup,action='read')
        do
            read(10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            text2=text//'  !!'
            
            read(text,*)  tmp_key, tmp_val
            if(tmp_val=='!!') cycle !missing value

            if(key==lalign(tmp_key)) then
                tmp_find=tmp_val
            endif

            if(key2==lalign(tmp_key)) then
                tmp_find=tmp_val
            endif

        end do
        close(10)

        find=lalign(tmp_find)

    end function

    function demand(key) result(res)
        character(*) :: key
        character(:),allocatable :: res

        !call email()

        if(mpiworld%is_master) then
            if(.not.exist) then
                open(12,file='tentative')
                write(12,'(a)') '!!  '//key//'    '
                close(12)
            endif

            do
                sleep(60)
                inquire(file='tentative', exist=exist)
                if(.not.exist) then
                    open(12,file='tentative')
                    write(12,'(a)') '!!  '//key//'    '
                    close(12)
                endif

                open(12,file='tentative')
                read(12,"(a)",iostat=msg) text ! read line into character variable
                close(12)
                if (msg < 0) stop 'tentative file is empty. Something is wrong.'
                if (msg > 0) stop 'Check tentative file.  Something is wrong.'
                
                if (text(1:2)=='!!') then
                    cycle !continue to wait for user input
                else
                    read(text,*) tmp_val
                    res=lalign(tmp_val)
                    exit
                endif
            enddo

            !delete tentative file
            open(12,file='tentative',status='old')
            close(12,status='delete')

            !append to setup file
            open(12,file=setup%file,status='unknown',access='append')
            write(12,'(a)') key//'  '//res
            close(12)

        endif

    end function

    function check(key,alias) result(exist)
        character(*) :: key
        character(*),optional :: alias !alias of inquired key

        exist=.false.

        if(present(alias)) then
            if(find(key,alias)/='') exist=.true.
        else
            if(find(key,key)/='') exist=.true.
        endif
        
    end function

    function get(key,key2,default,mandatory,iproc) result(res)
        character(*) :: key,key2
        character(*) :: default !default output value
        logical :: mandatory
        character(:),allocatable :: res

        character(:),allocatable :: keys
        keys=key//' ('//key2//') '

        res=find(key,key2)

        if(res/='') then
            call hud(keys//':'//res)
        endif

        if(res=='') then
            if(default/='') then
                call hud(keys//'is NOT found, take default: '//default)
                res=default
            elseif(mandatory==.true.) then
                call hud(keys//'is NOT found, but is MANDATORY.'//s_return &
                    'SeisJIMU has to ask for the value of this key before running further.'//s_return &
                    'Please open the text file "tentative" in the working directory,'//s_return &
                    'provide the missing value, remove "!!" in the line,'//s_return &
                    '(like in the setup file)'//s_return &
                    'and close the file.'//s_return &
                    'SeisJIMU checks this file every 1 min'//s_return &
                    'and will take the new value when no "!!" is leading the line.')
                res=demand(key)
            else
                call hud(keys//"is NOT found, take 0 for number(s), '' for character(s), or .false. for logical type(s)")
                res=''
            endif
        endif

        deallocate(keys)
        
    end function

    function get_int(key,alias,default,mandatory,iproc) result(res)
        character(*) :: key
        character(*),optional :: alias !alias of inquired key
        integer,optional :: default !default output value
        logical,optional :: mandatory
        integer :: res

        character(:),allocatable :: key2,def
        logical :: if_mandatory

        key2=key; if(present(alias)) key2=alias

        def=''; if(present(default)) def=int2str(default)

        if_mandatory=.false.; if(present(mandatory)) if_mandatory=mandatory

        res=str2int(get(key,key2,sdefault,if_mandatory))

    end function

    function get_ints(key,alias,default,mandatory,seperator,iproc) result(res)
        character(*) :: key
        character(*),optional :: alias !alias of inquired key
        integer,dimension(:),optional :: default !default output value
        logical,optional :: mandatory
        character(1),optional :: seperator
        integer,dimension(:) :: res

        character(:),allocatable :: key2,def
        logical :: if_mandatory

        key2=key; if(present(alias)) key2=alias

        def=''; if(present(default)) def=strcat(ints2strs(default))

        if_mandatory=.false.; if(present(mandatory)) if_mandatory=mandatory

        sep=' '; if(present(seperator)) sep=seperator
        res=strs2ints(split(get(key,key2,def,if_mandatory),o_sep=sep))

    end function

    function get_real(key,alias,default,mandatory,iproc) result(res)
        character(*) :: key
        character(*),optional :: alias !alias of inquired key
        integer,optional :: default !default output value
        logical,optional :: mandatory
        real :: res

    end function

    function get_reals(key,alias,default,mandatory,iproc) result(res)
        character(*) :: key
        character(*),optional :: alias !alias of inquired key
        integer,optional :: default !default output value
        logical,optional :: mandatory
        real,dimension(:),allocatable :: res

    end function

    function get_str(key,alias,default,mandatory,iproc) result(res)
        character(*) :: key
        character(*),optional :: alias !alias of inquired key
        integer,optional :: default !default output value
        logical,optional :: mandatory
        character(:),allocatable :: res

    end function
    
    function get_file(key,alias,default,mandatory,iproc) result(res)
        character(*) :: key
        character(*),optional :: alias !alias of inquired key
        integer,optional :: default !default output value
        logical,optional :: mandatory
        character(:),allocatable :: res

        if(present(alias)) then
            keys=key//' ('//alias//') '
            res=find(key,alias)
        else
            keys=key
            res=find(key,key)
        endif

        if(res/='') then
            inquire(file=res,exist=exist)
            if (.not. exist) res=''
        endif

        if(res=='') then
            call hud(keys//'is NOT found or do NOT exist, take empty filename.')
            read_file=''
        else
            call hud(keys//':'//res)
            read_file=res
        endif

        deallocate(keys,res)
        
    end function

    function get_files(key,alias,default,mandatory,iproc) result(res)

    end function
        
    function read_bool(key,alias,default,mandatory,iproc) result(res)
        character(*) :: key
        character(*),optional :: alias !alias of inquired key
        logical,optional :: default

        character(:),allocatable :: keys, res
        
        if(present(alias)) then
            keys=key//' ('//alias//') '
            res=find(key,alias)
        else
            keys=key
            res=find(key,key)
        endif

        if(res=='') then
            if(present(default)) then
                call hud(keys//'is NOT found, take default: '//bool2str(default)
                read_bool=default
            else
                call hud(keys//'is NOT found, take .false.')
                read_bool=.false.
            endif
        else
            call hud(keys//':'//res)
            read_bool=str2bool(res)
        endif

        deallocate(keys,res)

    end function

    function read_bools(key,alias,default,mandatory,iproc) result(res)

    end function

end