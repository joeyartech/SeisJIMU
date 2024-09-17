module m_Envnorm
use m_arrayop
use m_math
use m_Hilbert

    public
    private :: a, W2, v, E, DeltaE, water
    
    real,dimension(:,:),allocatable :: W2, v, E, DeltaE

    contains

    !E[u] = √(u² + v²)
    !Envsq[u] = a*∫ W²ΔE² dt = a*∫ W²(E-Eobs)² dt
    !KᵤEnvsq = 2a*W²*ΔE*KᵤE
    !KᵤE = E⁻¹(u + v*Kᵤv) ~= (u + v*Kᵤv)/(E+ε)
    !Since v=H[u], Kᵤv = -H[*]
    !where H denotes the Hilbert transform
    !So
    !KᵤEnvsq = 2a*  W²*ΔE*u/(E+ε)
    !         -2a*H[W²*ΔE*v/(E+ε)]

    !gradient = kernel*dt
    real function Envsq(scaler,W,u,Eobs,dt,nt,ntr)
        real,dimension(nt,ntr) :: W,u,Eobs

        a=scaler
        W2=W**2

        call alloc(v,nt,ntr)
        call hilbert_transform(u,v,nt,ntr)

        E = sqrt(u**2+v**2)
        DeltaE = E-Eobs

        water=maxval(E)*1e-5

        Envsq = a*sum( W2*DeltaE**2 )*dt
        
    end function

    subroutine kernel_Envsq(kernel,nt,ntr,oif_stack)
        real,dimension(nt,ntr) :: kernel
        logical,optional :: oif_stack
        
        real,dimension(:,:),allocatable :: tmp

        call alloc(tmp,nt,ntr)
        call hilbert_transform(W2*DeltaE*v/(E+water), tmp, nt,ntr)

        if(either(oif_stack,.false.,present(oif_stack))) then
            kernel = kernel - 2*a*( W2*DeltaE*u/(E+water) -tmp ) !why a minus?

        else
            kernel =        - 2*a*( W2*DeltaE*u/(E+water) -tmp )

        endif

        deallocate(W2,v,E,DeltaE,tmp)

    end subroutine

end
