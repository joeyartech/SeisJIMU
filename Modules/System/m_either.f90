module m_either

    interface either
        module procedure :: either2_int
        module procedure :: either2_ints
        module procedure :: either2_real
        module procedure :: either2_reals
        module procedure :: either2_str
    end interface

    contains

    function either2_int(i1,i2,switch) result(res)
        logical :: switch
        integer :: res
        if(switch) then
            res=i1
        else
            res=i2
        endif
    end function

    function either2_ints(i1,i2,switch) result(res)
        logical :: switch
        integer,dimension(:) :: i1,i2
        integer,dimension(:),allocatable :: res
        if(switch) then
            allocate(res,source=i1)
        else
            allocate(res,source=i2)
        endif
    end function

    function either2_real(r1,r2,switch) result(res)
        logical :: switch
        if(switch) then
            res=r1
        else
            res=r2
        endif
    end function

    function either2_reals(r1,r2,switch) result(res)
        logical :: switch
        real,dimension(:) :: r1,r2
        real,dimension(:),allocatable :: res
        if(switch) then
            allocate(res,source=r1)
        else
            allocate(res,source=r2)
        endif
    end function

    function either2_str(s1,s2,switch)
        logical :: switch
        character(*) :: s1,s2
        character(:),allocatable :: either2_str
        if(switch) then
            either2_str=s1
        else
            either2_str=s2
        endif
    end function

end