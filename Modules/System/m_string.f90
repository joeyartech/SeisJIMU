!Usual operations on the Fortran character type
!for more flexible IO purposes (m_setup.f90)
! 
!https://www.star.le.ac.uk/%7ecgp/fortran.html
!https://en.wikibooks.org/wiki/Fortran/strings#CHARACTER_Collating_Sequence
!https://github.com/szaghi/StringiFor

module m_string
    use m_either

    private :: fin

    !default string length
    integer,parameter :: i_str_slen  = 64
    integer,parameter :: i_str_len  = 128
    integer,parameter :: i_str_xlen = 256
    integer,parameter :: i_str_xxlen = 512
    integer,parameter :: i_str_xxxlen = 1024
    character(*), parameter :: s_NL = achar(10) !new line (NL) \n or linefeed (LF)
    character(*), parameter :: s_CR = achar(13) !carriage return (CR) \r
        !https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/494946
        !https://www.petefreitag.com/item/863.cfm
        !char function: CHAR(10) is linefeed LF. CHAR(13) is carriage return CR. If you are a little paranoid, ACHAR(10) is better - this is a little more robust to the default character kind not being ascii.
        !The NEW_LINE standard intrinsic is even more robust. There's also the C_NEW_LINE and C_CARRIAGE_RETURN constants from the ISO_C_BINDING module for kind C_CHAR characters.

    integer,parameter,private :: A2a = ichar('a') - ichar('A')
    character(*),parameter,private :: digits='1234567890'

    interface num2str
        module procedure int2str
        module procedure real2str
    end interface

    interface nums2strs
        module procedure ints2strs
        module procedure reals2strs
    end interface

    !encapsulate type 'character(len=*)' with methods
    !enable structures of variable-length strings
    type t_string
        character(:), allocatable :: s
        contains
        !procedure :: as_real => as_real
        !procedure :: as_complex => as_complex
        final :: fin
    end type

    contains

    subroutine fin(str)
        type(t_string) :: str
        deallocate(str%s)
    end subroutine

    function lalign(str)
        character(*),intent(in) :: str
        character(:),allocatable :: lalign
        lalign=trim(adjustl(str))

    end function

    function upper(str)
        character(*),intent(in) :: str
        character(len(str)) :: upper
        
        do i = 1,len(str)
            if('a'<=str(i:i).and.str(i:i)<='z') then !indexing with e.g. (i) won't allow compile. why?
                upper(i:i) = achar(iachar(str(i:i)) - A2a)
            else
                upper(i:i) = str(i:i)
            endif
        enddo
        
    end function

    function lower(str)
        character(*),intent(in) :: str
        character(len(str)) :: lower
        
        do i = 1,len(str)
            if('A'<=str(i:i).and.str(i:i)<='Z') then
                lower(i:i) = achar(iachar(str(i:i)) - a2A)
            else
                lower(i:i) = str(i:i)
            endif
        enddo
        
    end function

    function caps(str)
        character(*),intent(in) :: str
        character(len(str)) :: caps

        caps=str
        caps(1:1) = upper(str(1:1))
        
    end function

    function int2str(num,o_format) result(str)
        integer :: num
        character(*),optional :: o_format
        character(:),allocatable :: str
        
        character(i_str_len) :: tmp
        
        if(present(o_format)) then
            write(tmp,o_format) num
        else
            write(tmp,*) num
        endif
        
        str=lalign(tmp)

    end function
    
    function ints2strs(num,n,o_format) result(str)
        integer,dimension(n) :: num
        character(*),optional :: o_format
        type(t_string),dimension(:),allocatable :: str

        character(i_str_len) :: tmp
                
        if(allocated(str)) deallocate(str)
        allocate(str(n))

        do i=1,n
            if(present(o_format)) then
                write(tmp,o_format) num(i)
            else
                write(tmp,*) num(i)
            endif
            
            str(i)%s=lalign(tmp)
        enddo

    end function

    function real2str(num,o_format) result(str)
        real :: num
        character(*),optional :: o_format
        character(:),allocatable :: str
        
        character(i_str_len) :: tmp
        
        if(present(o_format)) then
            write(tmp,o_format) num
        else
            write(tmp,*) num
        endif
        
        str=lalign(tmp)

    end function

    function reals2strs(num,n,o_format) result(str)
        real,dimension(n) :: num
        character(*),optional :: o_format
        type(t_string),dimension(:),allocatable :: str
        
        character(i_str_len) :: tmp

        if(allocated(str)) deallocate(str)
        allocate(str(n))
        
        do i=1,n
            if(present(o_format)) then
                write(tmp,o_format) num(i)
            else
                write(tmp,*) num(i)
            endif

            str(i)%s=lalign(tmp)
        enddo
        
    end function

    function bool2str(bool) result(str)
        logical :: bool
        character(1) :: str

        if(bool) then
            str='T'
        else
            str='F'
        endif
    end function
    
    function str2int(str,o_format) result(num)
        character(*) :: str
        character(*),optional :: o_format
        integer :: num
        
        !in case of no digits
        if(scan(str,digits)==0) then
            num=0
            return
        endif
        
        !convert
        if(present(o_format)) then
            read(str,o_format) num
        else
            read(str,*) num
        endif
        
    end function

    function strs2ints(str,n,o_format) result(num)
        type(t_string),dimension(n) :: str
        character(*),optional :: o_format
        integer,dimension(:),allocatable :: num
        
        if(allocated(num)) deallocate(num)
        allocate(num(n))

        do i=1,n
            !in case of no digits
            if(scan(str(i)%s,digits)==0) then
                num(i)=0
                cycle
            endif
            
            !convert
            if(present(o_format)) then
                read(str(i)%s,o_format) num(i)
            else
                read(str(i)%s,*) num(i)
            endif
        enddo

    end function
    
    function str2real(str,o_format) result(num)
        character(*) :: str
        character(*),optional :: o_format
        real :: num
        
        !in case of no digits
        if(scan(str,digits)==0) then
            num=0
            return
        endif
        
        !convert
        if(present(o_format)) then
            read(str,o_format) num
        else
            read(str,*) num
        endif
        
    end function

    function strs2reals(str,n,o_format) result(num)
        type(t_string),dimension(n) :: str
        character(*),optional :: o_format
        real,dimension(:),allocatable :: num
        
        if(allocated(num)) deallocate(num)
        allocate(num(n))

        do i=1,n
            !in case of no digits
            if(scan(str(i)%s,digits)==0) then
                num(i)=0
                cycle
            endif
            
            !convert
            if(present(o_format)) then
                read(str(i)%s,o_format) num(i)
            else
                read(str(i)%s,*) num(i)
            endif
        enddo
        
    end function

    function str2bool(str) result(bool)
        character(*) :: str
        logical :: bool

        if(str(1:1)=='T'.or.str(1:1)=='t') then
            bool=.true.
        else
            bool=.false.
        endif
    end function

    !remove continuously repeated char from str_in
    !char is searched in backward order
    function remove_repetition(str_in,char) result(str_out)
        character(*) :: str_in
        character(1) :: char
        character(:),allocatable :: str_out
        
        character(:),allocatable :: tmp
        
        str_out=str_in
        i=index(str_out,char,back=.true.) !backward search
        
        do while (i>1)
            if(str_out(i-1:i-1)==char) then !remove repeated char
                tmp=str_out(1:i-2)//str_out(i:)  !can't let str_out on LHS due to mem conflict..
                str_out=tmp
                i=index(str_out,char,back=.true.)
            else !next backward search
                i=index(str_out(1:i-2),char,back=.true.)
            endif
        enddo
        
    end function
    
    !split a long string into string arrays separated by o_sep
    function split(str,o_sep) result(strs)
        character(*) :: str
        character(1),optional :: o_sep !separator must be 1 character
        type(t_string),dimension(:),allocatable :: strs

        character(1) :: sep
        character(:),allocatable :: text
        
        sep=either(o_sep,' ',present(o_sep))
        
        !regularize input string
        text = remove_repetition(sep//str//sep , sep)
        
        !trivial case
        if(text==sep) then
            allocate(strs(1)); strs(1)%s=''
            deallocate(text)
            return
        endif
        
        !count how many substrings
        n=0
        do k=1,len(text)
            if(text(k:k)==sep) n=n+1
        enddo
        n=n-1
        allocate(strs(n))
        
        !split
        i=1
        j=index(text(i+1:),sep)+i
        
        do k=1,n
            strs(k)%s = text(i+1:j-1)
            i=j
            j=index(text(i+1:),sep)+i
        enddo

        deallocate(text)

    end function

    !flatten a string array into a single long string glued by o_glue
    function strcat(strs,n,o_glue) result(str)
        type(t_string),dimension(n) :: strs
        character(1),optional :: o_glue
        character(:),allocatable :: str

        character(1) :: glue
        character(:),allocatable :: text

        glue=either(o_glue,' ',present(o_glue))

        maxlen=len(strs(1)%s)+1
        do k=2,either(n,1,n>1)
            length=len(strs(k)%s)+1
            maxlen=either(maxlen,length,maxlen>=length)
        enddo
        text=repeat(glue,maxlen*n)

        do k=0,n-1
            text(k*maxlen+1:k*maxlen+len(strs(k+1)%s))=strs(k+1)%s
        enddo

        text=text(1:len(text)-1) !remove very last glue
        
        str=lalign(remove_repetition(text,glue))

        deallocate(text)

    end function

end