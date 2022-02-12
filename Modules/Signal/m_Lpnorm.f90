module m_Lpnorm
use m_arrayop

    public

    private :: size, Wpu, dt_scalar

    integer :: size
    real,dimension(:),allocatable :: Wpu
    real dt_scalar

    contains

    !||u||1 = int |W*u| *dt *scalar
    real function L1norm(n,W,u,dt,scalar)
        real,dimension(*) :: W,u
        real :: dt, scalar

        size=n

        call alloc(Wpu,n)

        Wpu=W(1:n)*u(1:n)

        dt_scalar = dt*scalar

        L1norm = sum(abs(Wpu))*dt_scalar

    end function

    !nabla_u ||u||1 = sign(W*u) *dt *scalar
    subroutine nable_L1norm(nabla_u)
        real,dimension(*) :: nabla_u

        do i=1,size
            if (Wpu(i)>0.) then
                nabla_u(i) = dt_scalar
            else
                nabla_u(i) =-dt_scalar
            endif
        enddo

        deallocate(Wpu)

    end subroutine

    !||u||2^2 = int (W*u)^2 *dt *scalar
    real function L2norm_sq(n,W,u,dt,scalar)
        real,dimension(*) :: W,u
        real :: dt, scalar

        size=n

        call alloc(Wpu,n)

        Wpu=W(1:n)*W(1:n)*u(1:n)

        dt_scalar = dt*scalar

        L2norm_sq = sum(Wpu*u(1:n))*dt_scalar

    end function

    !nabla_u ||u||2^2 = 2* W*W*u *dt*scalar
    subroutine nabla_L2norm_sq(nabla_u)
        real,dimension(*) :: nabla_u

        nabla_u(1:size) = 2.*Wpu*dt_scalar

        deallocate(Wpu)

    end subroutine

end