!More procedures for character type,
!in addition to Fortran instrinsic procedures.
!could be futher developed for more OOP methods or mimic Python
!
!https://www.star.le.ac.uk/%7ecgp/fortran.html
!https://en.wikibooks.org/wiki/Fortran/strings#CHARACTER_Collating_Sequence
!https://github.com/szaghi/StringiFor

module m_string
!use m_const, only : i_str_len
i_str_len=80

    public
    private :: A2a

    integer,parameter :: A2a = ichar('a') - ichar('A')
    
    !particularly for partition()
    type t_string
        character(:), allocatable :: string
    end type

    interface num2str
        module procedure int2str
        module procedure real2str
    end interface

!     interface str2num   !can't differentiate str2int & str2real..
!         module procedure str2int
!         module procedure str2real
!     end interface
    
    contains

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


    function int2str(int,o_format) result(str)
        integer :: int
        character(*),optional :: o_format
        character(:),allocatable :: str
        
        character(i_str_len) :: tmp
        
        if(present(o_format)) then
            write(tmp,o_format) int
        else
            write(tmp,*) int
        endif
        
        str=trim(adjustl(tmp))

    end function
    
    function real2str(real,o_format) result(str)
        real :: real
        character(*),optional :: o_format
        character(:),allocatable :: str
        
        character(i_str_len) :: tmp
        
        if(present(o_format)) then
            write(tmp,o_format) real
        else
            write(tmp,*) real
        endif
        
        str=trim(adjustl(tmp))

    end function
    
    function str2int(str,o_format) result(int)
        character(*) :: str
        character(*),optional :: o_format
        integer :: int
        
        character(len(str)) :: tmp
        
        n=len(str)
        tmp=repeat(' ',n)
        
        !remove invalid characters in str
        j=n
        loop: do i = n,1,-1
            if('0'<=str(i:i).and.str(i:i)<='9') then
                tmp(j:j) = str(i:i)
                j=j-1
            elseif(str(i:i)=='+'.or.str(i:i)=='-') then
                tmp(j:j) = str(i:i)
                exit loop
            endif
        enddo loop
        
        !convert
        if(present(o_format)) then
            read(tmp,o_format) int
        else
            read(tmp,*) int
        endif
        
    end function
    
    function str2real(str,o_format) result(real)
        character(*) :: str
        character(*),optional :: o_format
        real :: real
        
        character(len(str)) :: tmp
        
        n=len(str)
        tmp=repeat(' ',n)
        
        !remove invalid characters in str
        j=n
        loop: do i = n,1,-1
            if('0'<=str(i:i).and.str(i:i)<='9') then
                tmp(j:j) = str(i:i)
                j=j-1
            elseif(str(i:i)=='.') then
                tmp(j:j) = str(i:i)
                j=j-1
            elseif(str(i:i)=='+'.or.str(i:i)=='-') then
                tmp(j:j) = str(i:i)
                exit loop
            endif
        enddo loop
        
        !convert
        if(present(o_format)) then
            read(tmp,o_format) real
        else
            read(tmp,*) real
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
    
    !count how many times (n) a char occurs in string
    function count_occurence(string,char) result(n)
        character(*) :: string
        character(1) :: char
        
        integer :: n
        
        n=0
        i=index(string,char,back=.true.)
        
        do while (i>0)
            n=n+1
            i=index(string(1:i-1),char,back=.true.)
        enddo
        
    end function
    
    function partition(string,o_separator) result(arr)
        character(*) :: string
        character(1),optional :: o_separator !separator must be 1 character
        type(t_string),dimension(:),allocatable :: arr

        character(:),allocatable :: tmp,str
        character(1) :: sep

        character(:),allocatable :: text
        
        sep=' '

        if(present(o_separator)) sep = o_separator
        
        !regularize input string
        str = remove_repetition(sep//string//sep , sep)
        
        !trivial case
        if(str==sep) then
            allocate(arr(1))
            arr(1)%string=''
            return
        endif
        
        !count how many elements
        n=count_occurence(str,sep)-1
        
        allocate(arr(n))
        
        !partition
        i=1
        j=index(str(i+1:),sep)+i
        
        do k=1,n
            arr(k)%string = str(i+1:j-1)
            i=j
            j=index(str(i+1:),sep)+i
        enddo

    end function

end