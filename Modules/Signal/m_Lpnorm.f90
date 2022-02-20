module m_Lpnorm
use m_arrayop

    public
    private :: n, Wpu, d
    
    real,dimension(:),allocatable :: Wpu

    contains

    !||u||_1 = int |W*u| *dt
    real function L1norm(size,W,u,sampling)
        integer size
        real,dimension(*) :: W,u

        n=size
        d=sampling

        call alloc(Wpu,n)
        Wpu=W(1:n)*u(1:n)

        L1norm = sum(abs(Wpu))*d

    end function

    !nabla_u ||u||_1 = sign(Wu)
    !gradient = nabla*dt
    !adjsource= -gradient
    subroutine adjsrc_L1norm(adjsrc,scaler)
        real,dimension(*) :: adjsrc
        real scaler !external scaler

        do i=1,n
            if (Wpu(i)>0.) then
                adjsrc(i) = adjsrc(i) -d *scaler
            else
                adjsrc(i) = adjsrc(i) +d *scaler
            endif
        enddo

        deallocate(Wpu)

    end subroutine

    !||u||2^2 = int (W*u)^2 *dt
    real function L2norm_sq(size,W,u,sampling)
        integer size
        real,dimension(*) :: W,u
        
        n=size
        d=sampling

        call alloc(Wpu,n)
        Wpu=W(1:n)*W(1:n)*u(1:n)
        
        L2norm_sq = sum(Wpu*u(1:n))*d
        
    end function

    !nabla_u ||u||_2^2 = 2* W*W*u
    !gradient = nabla*dt
    !adjsource= -gradient
    subroutine adjsrc_L2norm_sq(adjsrc,scaler)
        real,dimension(*) :: adjsrc
        real scaler !external scaler
        
        adjsrc(1:n) = adjsrc(i) -2.*Wpu*d *scaler

        deallocate(Wpu)

    end subroutine

end
