!Psuedotime-depth conversion
!Ref: R.E. Plessix, A pseudo-time formulation for acoustic full waveform inversion, Geophysical Journal International 192, No. 2 (2013): 613-630.
!
!   Given a background velocity model
!       regularly sampled in the depth domain v(z) or in the pseudotime domain v(t).
!
!   Domain conversion (coordinate transformaion):
!       depth(z) -> pseudotime(t):      dt = dz/v(z) ; t = ∫₀^∞ dz/v(z)
!       pseudotime(t) -> depth(z):      dz = v(t)*dt ; z = ∫₀^∞ v(t)*dt
!
!   Translating the model m (include the velocity) from one domain to another
!       Eg. z->t: dz is const but dt is not. However, we like regularly sampled v(t),
!       so we resample v(z) at const Δt rate,
!       where Δt is chosen as dz/vmax to be safe & simple.
!   Large model contrasts like seafloor should be enforced after conversion.
!
!   Translating the kernel KₘC from depth to pseudotime domain
!     - is a resampling problem as above if m is NOT the background velocity model; or
!     - if m IS the background velocity, the kernels have the following relation:
!       K_{v(t)}(z) = K_{v(z)}(z)*v(z) - ∫_z^∞ K_{v(z)}(z') dv(z')
!           where K_{v(t)}(z) is the kernel wrt v(t), as a function of z=∫₀^t v(t')dt'
!           K_{v(z)}(z) is the kernel wrt v(z), as a function of z
!           dv(z) = v'(z)dz is the differential of v(z)
!       then we can resample K_{v(t)}(z) in the pseudotime domain (same as 1st case)
!
!Brief proof of the relation:
!   Objective function C=C[v]. With perturbation δv(z) in depth domain, we have
!   δC = ∫₀^∞ K_{v(z)}(z) δv(z) dz
!        ∫₀^∞ K_{v(z)}(z) δv(z) v(t)dt
!   The rest is to convert δv(z) to δv(t). Two paths:
!       1) (t꜀≠tₚ) δv(t꜀) --δzₚ=0-> δtₚ --v'(tₚ)-> δv(tₚ) -> δv(zₚ)
!           ie. δv @t꜀ => stretch/squeeze time coord @tₚ => differential of v @zₚ
!       2) (t꜀=tₚ) δv(t꜀)                       -> δv(tₚ) -> δv(zₚ)
!   1st path: 0 =δzₚ =v(tₚ)δtₚ +∫₀^tₚ δv(t꜀) dt꜀
!       => δtₚ =-v(tₚ)⁻¹∫₀^tₚ δv(t꜀) dt꜀.   Since  δv(zₚ)=δv(tₚ)=v'(tₚ)δtₚ,
!       δv(zₚ) =-v'(tₚ)*v(tₚ)⁻¹∫₀^tₚ δv(t꜀) dt꜀
!   2nd path: δv(zₚ) = δv(t꜀)
!   Together, we have
!       δv(zₚ) = ∫₀^∞ [δ(t꜀-tₚ) -v'(tₚ)*v(tₚ)⁻¹B(t꜀,tₚ)] δv(t꜀)dt꜀
!       where B(t꜀,tₚ)=H(t-t꜀)-H(t-tₚ)
!   Inserting to the above integration, and swapping the seq of two integrals, we have
!   δC = ∫₀^∞ [K_{v(z)}(z꜀)v(t꜀) -∫_t꜀^∞ K_{v(z)}(zₚ)dv(tₚ)] δv(t꜀)dt꜀
!   that is
!   K_{v(t)}(z) = K_{v(z)}(z)v(t) - ∫ₜ^∞  K_{v(z)}(z')dv(t')
!               = K_{v(z)}(z)v(z) - ∫_z^∞ K_{v(z)}(z')dv(t')
!   Caveat:
!   - In the 2nd term, we don't wanna move the differential (d) to K_{v(z)}(z'),
!       because the kernel function can have small-scale artificial variations.
!   - IMPORTANT for linesearch, optimization etc:
!       In practice, we have a fixed max depth of the modeling domain (Z), 
!       K_{v(z)}(z) = G(z)/dz*B(0,Z), where G(z) is the input computed gradient.
!       Converting to the pseudotime domain, the corresponding max pseudotime
!           T=∫₀^Z dz/v(z) also changes with v(z)
!       that means ∫₀^T K_{v(t)} δv(t)dt /= δC  that we expect.
!       Implementation: zero the elements of K_{v(t)} from index nt_real+1 to nt
!           and nt_real is updated each time.

module m_pseudotime
use m_arrayop

    private
    public :: pseudotime_init, pseudotime_reset, pseudotime_convert, pseudotime_convert_gradient

    integer :: nx,ny
    integer :: nz,nt
    real :: Dz,Dt !Δz, Δt

    contains

    subroutine pseudotime_init(dir,vmin,vmax,nx_,ny_,nz_,Dz_,nt_,Dt_)
        character(4) :: dir

        nx=nx_; ny=ny_

        if(dir=='z->t') then
            Dz=Dz_;  nz=nz_

            !Dt = dz/v(z) >= Dz/vmax
            Dt=Dz/vmax

            !t = ∫₀ᶻ dz'/v(z') <= (nz-1)*Dz /vmin
            tmax=(nz-1)*Dz/vmin

            !nt <= tmax/Dt
            nt=ceiling(tmax/Dt)+1

            !write(*,*) 'Dt,nt = ',Dt,nt
            dt_=Dt;  nt_=nt
        endif

        if(dir=='t->z') then
            Dt=Dt;  nt=nt_

            !dz = dt*v(z) >= Dt*vmin
            dz=Dt*vmin
            
            !z = ∫₀ᵗ v(t')*dt' <= (nt-1)*Dt *vmax
            zmax=(nt-1)*Dt*vmax

            !nz = zmax/Dz
            nz=ceiling(zmax/Dz)+1

            !write(*,*) 'Dz,nz = ',Dz,nz
            dz_=Dz;  nz_=nz
        endif
      
    end subroutine

    subroutine pseudotime_reset(nz_,Dz_,nt_,Dt_)
        nz=nz_; Dz=Dz_
        nt=nt_; Dt=Dt_
    end subroutine
    
    subroutine pseudotime_convert(dir,din,dout,o_v,o_nreal)
        character(4) :: dir
        real,dimension(:,:,:),target :: din
        real,dimension(:,:,:),allocatable :: dout
        real,dimension(:,:,:),target,optional :: o_v
        integer,dimension(nx,ny),optional :: o_nreal

        real,dimension(:,:,:),pointer :: v
        integer,dimension(nx,ny) :: nreal

        !background velocity model
        if(present(o_v)) then
            v=>o_v
        else
            v=>din
        endif

        if(dir=='z->t') then !from model to parameter

            call alloc(dout,nt,nx,ny)

            !real nt needed
            nreal = ceiling(sum(Dz/v,dim=1)/Dt)
            nreal = min(nreal,nt)

            do iy=1,ny; do ix=1,nx

                !loop
                z=0.
                do it=1,nreal(ix,iy)
                    dout(it,ix,iy)=interp(din(:,ix,iy),z,nz,Dz)
                    z=z+interp(v(:,ix,iy),z,nz,Dz)*Dt
                enddo

                dout(nreal(ix,iy)+1:nt,ix,iy)=dout(nreal(ix,iy),ix,iy)

            enddo; enddo

        endif

        if(dir=='t->z') then !from parameter to model

            call alloc(dout,nz,nx,ny)

            !real nz needed
            nreal = ceiling(sum(v*Dt,dim=1)/Dz)
            nreal = min(nreal,nz)

            do iy=1,ny; do ix=1,nx

                !loop
                t=0.
                do iz=1,nreal(ix,iy)
                   dout(iz,ix,iy)=interp(din(:,ix,iy),t,nt,Dt)
                   t=t+Dz/interp(v(:,ix,iy),t,nt,Dt)
                enddo

                dout(nreal(ix,iy)+1:nz,ix,iy)=dout(nreal(ix,iy),ix,iy)

            enddo; enddo

        endif

        if(present(o_nreal)) o_nreal=nreal

    end subroutine

    subroutine pseudotime_convert_gradient(gin,v_z,gout)
        real,dimension(nz,nx,ny) :: gin, v_z
        real,dimension(:,:,:),allocatable :: gout

        integer,dimension(nx,ny) :: ntreal
        real,dimension(:,:,:),allocatable :: K, dv
        
        call alloc(dv,nz,nx,ny)

        call alloc(gout,nt,nx,ny)
        
        K = gin/Dz

        do iy=1,ny; do ix=1,nx
            !compute dv
            !using central differences at interior nodes
            do iz=2,nz-1
               dv(iz,ix,iy)=(v_z(iz+1,ix,iy)-v_z(iz-1,ix,iy))/2.
            enddo
            !using one-side differences at boundary nodes
            dv( 1,ix,iy)=v_z( 2,ix,iy)-v_z(   1,ix,iy)
            dv(nz,ix,iy)=v_z(nz,ix,iy)-v_z(nz-1,ix,iy)
        enddo; enddo

        do iy=1,ny; do ix=1,nx
            do iz=1,nz
                K(iz,ix,iy) =K(iz,ix,iy)*v_z(iz,ix,iy) -sum(K(iz:nz,ix,iy)*dv(iz:nz,ix,iy))
            enddo
        enddo; enddo

        ! open(12,file='Kt(z)',access='direct',recl=4*size(K))
        ! write(12,rec=1) K
        ! close(12)

        !go to time domain
        call pseudotime_convert('z->t',K,gout,o_v=v_z,o_nreal=ntreal)
        gout=gout*Dt

        do iy=1,ny; do ix=1,nx
            gout(ntreal(ix,iy)+1:nt,ix,iy)=0.
        enddo; enddo

    end subroutine
    
    
    ! selecting one velocity value out of two adjacent values
    function interp(vin,z,n,d) result(v)
        real,dimension(n) :: vin

        i1=min(  floor(z/d)+1,   n) !maximum lower index
        i2=min(ceiling(z/d)+1,   n) !minimum upper index

        !distance inside interval
        z1=z-(i1-1)*d
        z2=(i2-1)*d-z
          
        v1=vin(i1)
        v2=vin(i2)

        ! !Nearest neighbour interpolation
        ! v=either(v1, v2, z1<z2)
        
        !Linear interpolation
        if(i1==i2) then
            v = v1
        else
            v = (v2*z1 + v1*z2)/d
        endif

    end function

end module
