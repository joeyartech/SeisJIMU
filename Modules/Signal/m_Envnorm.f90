module m_Envnorm
use m_arrayop
use m_math
use m_Hilbert

    public
    private :: v, E, DeltaE, water
    
    real,dimension(:,:),allocatable :: v, E, DeltaE

    contains

    !E[u] = √(u² + v²)
    !Envsq[u] = a*∫ W²ΔE² dt = a*∫ W²(E-Eobs)² dt
    !K_u Envsq = a*W²2ΔE*K_u E
    !K_u E = E⁻¹(u + v*K_u v) ~= (u + v*K_u v)/(E+ε)
    !Since v=H[u], where H is the Hilbert transform
    !K_u v = -H[*]
    !So
    !K_u Envsq =  a*W²   2ΔE*u/(E+ε)
    !            -a*W²*H[2ΔE*v/(E+ε)]

    !gradient = kernel*dt
    real function Envsq(a,W,u,Eobs,dt,nt,ntr)
        real,dimension(nt,ntr) :: W,u,Eobs

        call alloc(v,nt,ntr)
        call hilbert_transform(u,v,nt,ntr)

        E = sqrt(u**2+v**2)
        DeltaE = E-Eobs

        water=maxval(E)*1e-6

        Envsq = a*sum( (W*DeltaE)**2 )*dt
        
    end function

    subroutine kernel_Envsq(kernel,nt,ntr,oif_stack)
        real,dimension(nt,ntr) :: kernel
        logical,optional :: oif_stack
        
        if(either(oif_stack,.false.,present(oif_stack))) then
            kernel = kernel + a*W**2*2*(  DeltaE*u/(E+water) &
                                         -DeltaE*v/(E+water)  )

        else
            kernel =          a*W**2*2*(  DeltaE*u/(E+water) &
                                         -DeltaE*v/(E+water)  )

        endif

        deallocate(v,E,DeltaE)

    end subroutine

end
