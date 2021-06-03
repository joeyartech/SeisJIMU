!Usual operations on the Fortran character type
!for more flexible IO purposes (m_setup.f90)
! 
!https://www.star.le.ac.uk/%7ecgp/fortran.html
!https://en.wikibooks.org/wiki/Fortran/strings#CHARACTER_Collating_Sequence
!https://github.com/szaghi/StringiFor

module m_string
use m_global, only : i_str_len

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
    end type

    contains

    function lalign(str)
        character(*),intent(in) :: str
        character(:),allocatable :: lalign
        lalign=trim(adjustl(s))

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
    
    function ints2strs(num,o_format) result(str)
        integer,dimension(*) :: num
        character(*),optional :: o_format
        type(t_string),dimension(:),allocatable :: str

        character(i_str_len) :: tmp
                
        n=size(num); allocate(str(n))

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

    function reals2strs(num,o_format) result(str)
        real,dimension(*) :: num
        character(*),optional :: o_format
        type(t_string),dimension(:),allocatable :: str
        
        character(i_str_len) :: tmp
        
        do i=1,n
            if(present(o_format)) then
                write(tmp,o_format) num(i)
            else
                write(tmp,*) num(i)
            endif

            str(i)%s=lalign(tmp)
        enddo
        
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

    function strs2ints(str,o_format) result(num)
        type(t_string),dimension(*) :: str
        character(*),optional :: o_format
        integer,dimension(:),allocatable :: num
        
        n=size(str)

        do i=1,n
            !in case of no digits
            if(scan(str(i),digits)==0) then
                num(i)=0
                cycle
            endif
            
            !convert
            if(present(o_format)) then
                read(str(i),o_format) num(i)
            else
                read(str(i),*) num(i)
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

    function strs2reals(str,o_format) result(num)
        type(t_string),dimension(*) :: str
        character(*),optional :: o_format
        real,dimension(:),allocatable :: num
        
        n=size(str)

        do i=1,n
            !in case of no digits
            if(scan(str(i),digits)==0) then
                num(i)=0
                cycle
            endif
            
            !convert
            if(present(o_format)) then
                read(str(i),o_format) num(i)
            else
                read(str(i),*) num(i)
            endif
        enddo
        
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
        
        sep=' '
        if(present(o_sep)) sep = o_sep
        
        !regularize input string
        text = remove_repetition(sep//str//sep , sep)
        
        !trivial case
        if(text==sep) then
            allocate(strs(1)); strs(1)%s=''
            deallocate(text)
            return
        endif
        
        !count how many substrings
        n=count(text==sep)-1
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
    function strcat(strs,o_glue) result(str)
        type(t_string),dimension(*) :: strs
        character(1),optional :: o_glue
        character(:),allocatable :: str

        character(1) :: glue

        glue=' '
        if(present(o_glue)) glue=o_glue

        n=size(strs)
        text=repeat(' ' , n*maxval(len(strs(:)%s)) )

        i=1
        do k=1,n
            j=len(strs(k)%s)+1
            text(i,j)=str(k)%s//glue
            i=j
        enddo

        str=lalign(text)

    end function
end