module m_integral
use m_arrayop

    private
    public :: integral, nabla_integral

    integer :: m
    real,dimension(:),allocatable :: Kxm
    real dt_scalar

    contains

    !compute int K*x^n dt
    real function integral(size,K,x,n,dt,scalar)
        integer size
        real,dimension(*) :: K,x
        real :: dt, scalar

        m=n-1

        call alloc(Kxm,size)

        select case(n)
        case (1)
            Kxm=K(1:size)
        case (2)
            Kxm=K(1:size)*x(1:size)
        case default
            Kxm=K(1:size)*x(1:size)**m
        end select

        dt_scalar = dt*scalar

        integral = sum(Kxm(1:size)*x(1:size)) * dt_scalar

    end function

    !compute nabla_x int K*x^n dt
    !only works for n=2
    function nabla_integral()
        real,dimension(:),allocatable :: nabla_integral

        nabla_integral =  (m+1)*Kxm*dt_scalar

        deallocate(Kxm)

    end function

end