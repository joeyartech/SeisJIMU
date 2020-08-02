module m_setup
use m_const
use m_message

    character(:),allocatable :: file_setup

    character(:),allocatable :: work_dir
    character(:),allocatable :: scratch_dir
    
    contains
    
    subroutine setup_init(istat)
        integer :: istat
        logical :: alive
        character(i_str_xlen) :: tmp
        
        !get setup file
        istat=0
        
        call getarg(1,tmp)
        file_setup=trim(adjustl(tmp))
        
        if(file_setup=='') then
            istat=0
            call hud('No input setup file. Print manual..')
            return
        else
            inquire(file=file_setup, exist=alive)
            if(.not.alive) then
                call hud('Setup file '//file_setup//' does NOT exist!')
                istat=-1
                return
            endif
            
            call hud('Setup file: '//file_setup)
            istat=1

            work_dir=setup_get_char('DIR','WORK_DIR',default='./')
            scratch_dir=setup_get_char('SCRATCH_DIR','TMP_DIR',default=work_dir)

            return
        endif

    end subroutine
    
    function setup_ask(word,word2,default) result(found)
        character(*) :: word
        character(*),optional :: word2 !alias of inquired word
        logical,optional :: default
        
        character(i_str_len) :: text
        character(i_str_len) :: tmp,tmp2,tmp3
        
        logical :: found
        
        found=.false.
        if(present(default)) found=default
        
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp)) ) then
                found=.true.
            endif
        end do
        close(10)
    end function
    
    integer function setup_get_int(word,word2,default)
        character(*) :: word
        character(*),optional :: word2 !alias of inquired word
        integer,optional :: default
        
        character(i_str_len) :: text
        character(i_str_len) :: tmp,tmp2,tmp3
        
        logical :: found
        
        found=.false.
                
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    read(text,*) tmp2, setup_get_int
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp)) ) then
                read(text,*) tmp2, setup_get_int
                found=.true.
            endif
        end do
        close(10)
        
        if(found) then
            write(tmp,*) setup_get_int
            call hud(word//' '//trim(adjustl(tmp)))
        else
            if(present(default)) then
                setup_get_int=default
                write(tmp,*) default
                call hud(word//' is NOT found, take default: '//trim(adjustl(tmp)))
            else
                setup_get_int=0
                call hud(word//' is NOT found, take 0')
            endif
        endif
        
    end function
    
    real function setup_get_real(word,word2,default)
        character(*) :: word
        character(*),optional :: word2
        real, optional :: default
        
        character(i_str_len) :: text
        character(i_str_len) :: tmp,tmp2,tmp3
        
        logical :: found
        
        found=.false.
        
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    read(text,*) tmp2, setup_get_real
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp)) ) then
                read(text,*) tmp2, setup_get_real
                found=.true.
            endif
        end do
        close(10)
        
        if(found) then
            write(tmp,*) setup_get_real
            call hud(word//' '//trim(adjustl(tmp)))
        else
            if(present(default)) then
                setup_get_real=default
                write(tmp,*) default
                call hud(word//' is NOT found, take default: '//trim(adjustl(tmp)))
                
            else
                setup_get_real=0.
                call hud(word//' is NOT found, take 0.')
            endif
        endif
        
    end function
    
    function setup_get_char(word,word2,default)
        character(:),allocatable :: setup_get_char
        character(*) :: word
        character(*),optional :: word2
        character(*),optional :: default
        
        character(i_str_xlen) :: text
        character(i_str_xlen) :: tmp,tmp2,tmp3
        
        logical :: found
        
        found=.false.
        
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    read(text,*) tmp2, tmp3
                    setup_get_char=trim(adjustl(tmp3))
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp))) then
                read(text,*) tmp2, tmp3
                setup_get_char=trim(adjustl(tmp3))
                found=.true.
            endif
        end do
        close(10)
        
        if(found) then
            if(mpiworld%is_master) write(*,*) word//' '//setup_get_char
        else
            if(present(default)) then
                setup_get_char=default
                call hud(word//' is NOT found, take default: '//setup_get_char)
            else
                setup_get_char=''
                call hud(word//' is NOT found, take empty string')
            endif
        endif
        
    end function
    
    function setup_get_file(word,word2)
        character(:),allocatable :: setup_get_file
        character(*) :: word
        character(*),optional :: word2
        
        character(i_str_len) :: text
        character(i_str_len) :: tmp,tmp2,tmp3,tmp4
        
        logical :: found, alive
        
        found=.false.
        alive=.false.
        
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    read(text,*) tmp2, tmp3
                    tmp4=tmp3
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp))) then
                read(text,*) tmp2, tmp3
                tmp4=tmp3
                found=.true.
            endif
        end do
        close(10)
        
        if(found) inquire(file=trim(adjustl(tmp4)),exist=alive)
        
        found=found.and.alive
        
        
        if(found) then
            setup_get_file=trim(adjustl(tmp4))
            if(mpiworld%is_master) write(*,*) word//' '//setup_get_file
        else
            setup_get_file=''
            call hud(word//' is NOT found, take empty filename')
        endif
        
    end function
        
    logical function setup_get_logical(word,word2,default)
        character(*) :: word
        character(*),optional :: word2
        logical, optional :: default
        
        character(i_str_len) :: text
        character(i_str_len) :: tmp,tmp2,tmp3
        
        logical :: found
        
        found=.false.
        
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    read(text,*) tmp2, setup_get_logical
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp)) ) then
                read(text,*) tmp2, setup_get_logical
                found=.true.
            endif
        end do
        close(10)
        
        if(found) then
            write(tmp,*) setup_get_logical
            call hud(word//' '//trim(adjustl(tmp)))
        else
            if(present(default)) then
                setup_get_logical=default
                write(tmp,*) default
                call hud(word//' is NOT found, take default: '//trim(adjustl(tmp)))
                
            else
                setup_get_logical=.false.
                call hud(word//' is NOT found, take .false.')
            endif
        endif
        
    end function

end