module m_Envnorm
use m_arrayop
use m_math
use m_Hilbert

    public
    private :: a, W2, u, v, E, DeltaE, water

    real,dimension(:,:),allocatable :: W2, u, v, E, DeltaE

    real,private,dimension(:,:),allocatable :: DeltaQ

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
    real function Envsq(scaler,nt,ntr,W,dsyn,Eobs,dt)
        real,dimension(nt,ntr) :: W,dsyn,Eobs

        a=scaler
        W2=W**2

        u=dsyn

        call alloc(v,nt,ntr)
        call hilbert_transform(u,v,nt,ntr)

        E = sqrt(u**2+v**2)
        DeltaE = Eobs-E !as defined, tobe consistent w/ L2sq

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
            kernel = kernel + 2*a*( W2*DeltaE*u/(E+water) -tmp )

        else
            kernel =        + 2*a*( W2*DeltaE*u/(E+water) -tmp )

        endif

        deallocate(W2,u,v,E,DeltaE,tmp)

    end subroutine


    !Q[u] = u*∂ₜv - ∂ₜu*v
    !Qsq[u] = a*∫ W²ΔQ² dt = a*∫ W²(Q-Qobs)² dt
    !KᵤQsq = 2a*W²*ΔQ*KᵤQ
    !W²*ΔQ*Kᵤ(u*∂ₜv) = W²*ΔQ*∂ₜv - ∂ₜ(W²*ΔQ*u)*Kᵤv = W²*ΔQ*∂ₜv + H[∂ₜ(W²*ΔQ*u)]
    !W²*ΔQ*Kᵤ(∂ₜu*v) =-∂ₜ(W²*ΔQ*v) + W²*ΔQ*∂ₜu*Kᵤv =-∂ₜ(W²*ΔQ*v) - H[W²*ΔQ*∂ₜu]
    !So
    !KᵤQsq = 2a*W²*ΔQ
    !      = 2a*( W²*ΔQ*∂ₜv + H[∂ₜ(W²*ΔQ*u)] - ∂ₜ(W²*ΔQ*v) - H[W²*ΔQ*∂ₜu] )
    real function Qsq(scaler,nt,ntr,W,dsyn,dobs,dt)
        real,dimension(nt,ntr) :: W,dsyn,dobs

        real,dimension(:,:),allocatable :: du, Qobs, Q

        a=scaler
        W2=W**2

        u=dobs
        call alloc(v,nt,ntr)
        call hilbert_transform(u,v,nt,ntr)
        Qobs = u*time_deri(v,dt) - time_deri(u,dt)*v

        u=dsyn
        call alloc(v,nt,ntr)
        call hilbert_transform(u,v,nt,ntr)
        Q = u*time_deri(v,dt) - time_deri(u,dt)*v

        DeltaQ = Qobs-Q !as defined, tobe consistent w/ L2sq

        Qsq = a*sum( W2*DeltaQ**2 )*dt

    end function

    subroutine kernel_Qsq(kernel,nt,ntr,dt,oif_stack) !tobe modified
        real,dimension(nt,ntr) :: kernel
        logical,optional :: oif_stack

        real,dimension(:,:),allocatable :: tmp1, tmp2

        call alloc(tmp1,nt,ntr)
        call alloc(tmp2,nt,ntr)
        call hilbert_transform(time_deri(W2*DeltaQ*u,dt), tmp1, nt,ntr)
        call hilbert_transform(W2*DeltaQ*time_deri(u,dt), tmp2, nt,ntr)
        
        if(either(oif_stack,.false.,present(oif_stack))) then
            kernel = kernel + 2*a*( W2*DeltaQ*time_deri(v,dt) + tmp1 &
                                   -time_deri(W2*DeltaQ*v,dt) - tmp2 )

        else
            kernel =          2*a*( W2*DeltaQ*time_deri(v,dt) + tmp1 &
                                   -time_deri(W2*DeltaQ*v,dt) - tmp2 )


        endif

        deallocate(W2,u,v,DeltaQ,tmp1,tmp2)

    end subroutine




    function time_deri(data,dt) !gradient by cdiff
        real,dimension(:,:) :: data !nt x nrcv
        real,dimension(:,:),allocatable :: time_deri

        nt=size(data,1)
        nr=size(data,2)

        call alloc(time_deri,nt,nr)

            time_deri(1,:) = (data(2,:)-data(1,:))/dt

        do i=2,nt-1
            time_deri(i,:) = (data(i+1,:)-data(i-1,:))/2./dt
        enddo

            time_deri(nt,:) = (data(nt,:)-data(nt-1,:))/dt

    end function

end
