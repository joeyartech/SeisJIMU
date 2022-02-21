module m_Lpnorm
use m_arrayop

    public
    private :: n, Wpu, d_scaler
    
    real,dimension(:),allocatable :: Wpu

    contains

    !||u||_1 = int |W*u| *dt
    real function L1norm(size,W,u,sampling,scaler)
        integer size
        real,dimension(*) :: W,u

        n=size
        d_scaler=sampling*scaler

        call alloc(Wpu,n)
        Wpu=W(1:n)*u(1:n)

        L1norm = sum(abs(Wpu))*d_scaler

    end function

    !nabla_u ||u||_1 = sign(Wu)
    !gradient = nabla*dt
    !adjsource= -gradient
    subroutine adjsrc_L1norm(adjsrc)
        real,dimension(*) :: adjsrc

        do i=1,n
            if (Wpu(i)>0.) then
                adjsrc(i) = adjsrc(i) -d_scaler
            else
                adjsrc(i) = adjsrc(i) +d_scaler
            endif
        enddo

        deallocate(Wpu)

    end subroutine

    !||u||2^2 = int (W*u)^2 *dt
    real function L2norm_sq(size,W,u,sampling,scaler)
        integer size
        real,dimension(*) :: W,u
        
        n=size
        d_scaler=sampling*scaler

        call alloc(Wpu,n)
        Wpu=W(1:n)*W(1:n)*u(1:n)
        
        L2norm_sq = sum(Wpu*u(1:n))*d_scaler
        
    end function

    !nabla_u ||u||_2^2 = 2* W*W*u
    !gradient = nabla*dt
    !adjsource= -gradient
    subroutine adjsrc_L2norm_sq(adjsrc)
        real,dimension(*) :: adjsrc
        
        adjsrc(1:n) = adjsrc(1:n) -2.*Wpu*d_scaler

        deallocate(Wpu)

    end subroutine

end
