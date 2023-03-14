module m_either

    interface either
        module procedure :: either2_int
        module procedure :: either2_ints
        module procedure :: either2_real
        module procedure :: either2_reals
        module procedure :: either2_str
        module procedure :: either2_bool
    end interface

    contains

    function either2_int(i1,i2,switch) result(res)
        integer,optional   :: i1
        integer,intent(in) :: i2
        logical,intent(in) :: switch
        integer :: res
        if(switch) then
            res=i1
        else
            res=i2
        endif
    end function

    function either2_ints(i1,i2,switch) result(res)
        integer,dimension(:),optional   :: i1
        integer,dimension(:),intent(in) :: i2
        logical,intent(in) :: switch
        integer,dimension(:),allocatable :: res
        if(switch) then
            allocate(res,source=i1)
        else
            allocate(res,source=i2)
        endif
    end function

    function either2_real(r1,r2,switch) result(res)
        real,optional   :: r1
        real,intent(in) :: r2
        logical,intent(in) :: switch
        if(switch) then
            res=r1
        else
            res=r2
        endif
    end function

    function either2_reals(r1,r2,switch) result(res)
        real,dimension(:),optional   :: r1
        real,dimension(:),intent(in) :: r2
        logical,intent(in) :: switch
        real,dimension(:),allocatable :: res
        if(switch) then
            allocate(res,source=r1)
        else
            allocate(res,source=r2)
        endif
    end function
    
    function either2_str(s1,s2,switch) result(res)
        character(*),optional   :: s1
        character(*),intent(in) :: s2
        logical,intent(in) :: switch
        character(:),allocatable :: res
        if(switch) then
            res=s1
        else
            res=s2
        endif
    end function

    function either2_bool(b1,b2,switch) result(res)
        logical,optional   :: b1
        logical,intent(in) :: b2, switch
        logical :: res
        if(switch) then
            res=b1
        else
            res=b2
        endif
    end function

end
