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
    !nabla = a*sgn(Wu)
    !gradient = nabla*dt
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

    subroutine nabla_L1(nabla,oif_stack)
        real,dimension(*) :: nabla
        logical,optional :: oif_stack

        if(either(oif_stack,.false.,present(oif_stack))) then
            do i=1,n
                if (Wpu(i)>0.) then
                    nabla(i) = nabla(i) +a
                else
                    nabla(i) = nabla(i) -a
                endif
            enddo
        else
            do i=1,n
                if (Wpu(i)>0.) then
                    nabla(i) =           a
                else
                    nabla(i) =          -a
                endif
            enddo
        endif

        deallocate(Wpu)

    end subroutine

    !║u║₂² = a*∫ (Wu)² dt
    !nabla = 2a*W²u
    !gradient = nabla*dt
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

    subroutine nabla_L2sq(nabla,oif_stack)
        real,dimension(*) :: nabla
        logical,optional :: oif_stack
        
        if(either(oif_stack,.false.,present(oif_stack))) then
            nabla(1:n) = nabla(1:n) +2.*a*Wpu
        else
            nabla(1:n) =             2.*a*Wpu
        endif

        deallocate(Wpu)

    end subroutine

end
