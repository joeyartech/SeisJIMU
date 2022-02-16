module m_Lpnorm
use m_arrayop

    public
    private :: n, Wpu, dt_scalar

    integer :: n
    real,dimension(:),allocatable :: Wpu
    real dt_scalar

    contains

    !||u||_1 = int |W*u| *dt *scalar
    real function L1norm(size,W,u,dt,scalar)
        integer size
        real,dimension(*) :: W,u
        real :: dt, scalar

        n=size

        call alloc(Wpu,n)

        Wpu=W(1:n)*u(1:n)

        dt_scalar = dt*scalar

        L1norm = sum(abs(Wpu))*dt_scalar

    end function

    !nabla_u ||u||_1 = sign(W*u) *scalar
    !gradient = nabla*dt
    !adjsource= -gradient
    subroutine adjsrc_L1norm(adjsrc)
        real,dimension(*) :: adjsrc

        do i=1,n
            if (Wpu(i)>0.) then
                adjsrc(i) =-dt_scalar
            else
                adjsrc(i) = dt_scalar
            endif
        enddo

        deallocate(Wpu)

    end subroutine

    !||u||2^2 = int (W*u)^2 *dt *scalar
    real function L2norm_sq(size,W,u,dt,scalar)
        integer size
        real,dimension(*) :: W,u
        real :: dt, scalar
        
        n=size
        
        call alloc(Wpu,n)
        Wpu=W(1:n)*W(1:n)*u(1:n)
        dt_scalar = dt*scalar
        
        L2norm_sq = sum(Wpu*u(1:n))*dt_scalar
        
    end function

    !nabla_u ||u||_2^2 = 2* W*W*u *scalar
    !gradient = nabla*dt
    !adjsource= -gradient
    subroutine adjsrc_L2norm_sq(adjsrc)
        real,dimension(*) :: adjsrc
        
        adjsrc(1:n) = -2.*Wpu*dt_scalar

        deallocate(Wpu)

    end subroutine

end
