module m_Lpnorm
use m_arrayop
use m_math

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
    !kernel = a*sgn(Wu)
    !gradient = kernel*dt
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

    subroutine kernel_L1(kernel,oif_stack)
        real,dimension(*) :: kernel
        logical,optional :: oif_stack

        if(either(oif_stack,.false.,present(oif_stack))) then
            do i=1,n
                if (Wpu(i)>0.) then
                    kernel(i) = kernel(i) +a
                else
                    kernel(i) = kernel(i) -a
                endif
            enddo
        else
            do i=1,n
                if (Wpu(i)>0.) then
                    kernel(i) =           a
                else
                    kernel(i) =          -a
                endif
            enddo
        endif

        deallocate(Wpu)

    end subroutine

    !║u║₂² = a*∫ (Wu)² dt
    !kernel = 2a*W²u
    !gradient = kernel*dt
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

    subroutine kernel_L2sq(kernel,oif_stack)
        real,dimension(*) :: kernel
        logical,optional :: oif_stack
        
        if(either(oif_stack,.false.,present(oif_stack))) then
            kernel(1:n) = kernel(1:n) +2.*a*Wpu
        else
            kernel(1:n) =             2.*a*Wpu
        endif

        deallocate(Wpu)

    end subroutine

end
