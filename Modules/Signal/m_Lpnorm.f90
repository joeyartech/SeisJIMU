module m_Lpnorm
use m_arrayop

    public
    private :: n, Wpu, d, a
    
    real,dimension(:),allocatable :: Wpu

    contains

    ! function s_L1(str)
    !     character(*) :: str
    !     character(:),allocatable :: s_L1
    !     s_L1norm = '║'//str//'║₁'
    ! end function

    ! function s_L2sq(str)
    !     character(*) :: str
    !     character(:),allocatable :: s_L2sq
    !     s_L2norm_sq = '║'//str//'║₂²'
    ! end function


    !║u║₁ = a*∫ |Wu| dt
    !∇║u║₁ = a*sign(Wu)
    !adjsource= -∇║u║₁
    !gradient =  ∇║u║₁*dt
    real function L1(scaler,size,W,u,sampling)
        integer size
        real,dimension(*) :: W,u

        a=scaler
        n=size
        d=sampling

        call alloc(Wpu,n)
        Wpu=W(1:n)*u(1:n)

        L1 = a*sum(abs(Wpu))*d

    end function

    subroutine adjsrc_L1(adjsrc)
        real,dimension(*) :: adjsrc

        do i=1,n
            if (Wpu(i)>0.) then
                adjsrc(i) = adjsrc(i) -a
            else
                adjsrc(i) = adjsrc(i) +a
            endif
        enddo

        deallocate(Wpu)

    end subroutine

    !║u║₂² = a*∫ (Wu)² dt
    !∇║u║₂² = 2a*W²u
    !adjsource= -∇║u║₂²
    !gradient =  ∇║u║₂²*dt
    real function L2sq(scaler,size,W,u,sampling)
        integer size
        real,dimension(*) :: W,u

        a=scaler
        n=size
        d=sampling
        
        call alloc(Wpu,n)
        Wpu=W(1:n)*W(1:n)*u(1:n)
        
        L2sq = a*sum(Wpu*u(1:n))*d
        
    end function

    subroutine adjsrc_L2sq(adjsrc)
        real,dimension(*) :: adjsrc
        
        adjsrc(1:n) = adjsrc(1:n) -2.*a*Wpu

        deallocate(Wpu)

    end subroutine

end
