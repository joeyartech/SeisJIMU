!Psuedotime-depth conversion
!Ref: R.E. Plessix, A pseudo-time formulation for acoustic full waveform inversion, Geophysical Journal International 192, No. 2 (2013): 613-630.
!
!   Given a background velocity model
!       sampled in the depth domain v(z) or in the pseudotime domain v(t).
!
!   Domain conversion:
!       depth(z) -> pseudotime(t):      t = ∫₀ᶻ v(z')⁻¹ dz'
!       pseudotime(t) -> depth(z):      z = ∫₀ᵗ v(t') dt'
!
!   Gradient conversion t->z (chain rule):
!       ∇ₓC(tₙ) = ∇ₘC(zₙ) - ∫_zₙ^zmax ∇ₘC(z') v(z')⁻¹ dv(z')/dz' dz'
!         where zₙ= ∫₀^tₙ v(t') dt'
!         Note that v(tₙ) canbe different from v(zₙ) or v(z')

module m_pseudotime
use m_arrayop

    private
    public :: pseudotime_init, pseudotime_z2t, pseudotime_z2t_seafloor, pseudotime_t2z, pseudotime_z2t_gradient

    integer :: nz,nt,nx,ny
    real :: dz,dt

    contains

    subroutine pseudotime_init(dir,vmin,vmax,nx_,ny_,nz_,dz_,nt_,dt_)
        character(4) :: dir

        nx=nx_
        ny=ny_

        if(dir=='t->z') then
            dt=dt_;  nt=nt_

            !dz = dt*v(z) >= dt*vmin
            dz=dt*vmin
            
            !z = ∫₀ᵗ v(t')*dt' <= (nt-1)*dt vmax 
            zmax=(nt-1)*dt*vmax

            !nz = zmax/dz
            nz=ceiling(zmax/dz)+1

            !write(*,*) 'dz,nz = ',dz,nz
            dz_=dz;  nz_=nz
        endif

        if(dir=='z->t') then
            dz=dz_;  nz=nz_

            !dt = dz/v(z) >= dz/vmax
            dt=dz/vmax

            !t = ∫₀ᶻ dz'/v(z') <= (nz-1)*dz /vmin
            tmax=(nz-1)*dz/vmin

            !tmax=(nz-1)*dz/(     vmin*1.5   )
            ! tmax=(nz-1)*dz/(0.5*(vmin+vmax)) ! Mean velocity instead of max
            !   t_max=0.
            !   do ix=1,nx
            !      t= 2*sum(1./v(:,ix)) - 1./v(1,ix) - 1./v(nz,ix)  !the integral
            !      t=t*dz/2                                         !the integral
            !      if (t>t_max) t_max=t
            !   enddo

            !nt <= tmax/dt
            nt=ceiling(tmax/dt)+1

            !write(*,*) 'dt,nt = ',dt,nt
            dt_=dt;  nt_=nt
        endif
      
    end subroutine

    !converting model to parameter z->t
    subroutine pseudotime_z2t(min,xout,o_v_z)
        real,dimension(nz,nx,ny) :: min
        real,dimension(:,:,:),allocatable :: xout

        real,dimension(nz,nx,ny),optional :: o_v_z

        call alloc(xout,nt,nx,ny)
     
        if(present(o_v_z)) then !converting other model
            do iy=1,ny; do ix=1,nx
                z=0.
                do it=1,nt
                    xout(it,ix,iy)=interp(min(:,ix,iy),z,nz,dz)
                    z=z+interp(o_v_z(:,ix,iy),z,nz,dz)*dt
                enddo
            enddo; enddo

        else !converting vp model itself
            do iy=1,ny; do ix=1,nx
                z=0.
                do it=1,nt
                    v=interp(min(:,ix,iy),z,nz,dz)
                    xout(it,ix,iy)=v
                    z=z+v*dt
                enddo
            enddo; enddo

        endif

    end subroutine


    ! calculates the seafloor in the time domain
    ! enables to mute the time-domain gradient in the bathymetry
    function pseudotime_z2t_seafloor(vw,iz) result(it)
        it= nint(   (iz-1)*dz/vw   /dt) +1
    end function


    !converting parameter to model t->z
    subroutine pseudotime_t2z(xin,mout,o_v_t)
        real,dimension(nt,nx,ny) :: xin
        real,dimension(:,:,:),allocatable :: mout

        real,dimension(nt,nx,ny),optional :: o_v_t

        call alloc(mout,nz,nx,ny)
     
        if(present(o_v_t)) then !converting other parameter
            do iy=1,ny; do ix=1,nx
                t=0.
                do iz=1,nz
                   mout(iz,ix,iy)=interp(xin(:,ix,iy),t,nt,dt)
                   t=t+dz/interp(o_v_t(:,ix,iy),t,nt,dt)
                enddo
            enddo; enddo

        else !converting vp model itself
            do iy=1,ny; do ix=1,nx
                t=0.
                do iz=1,nz
                    v=interp(xin(:,ix,iy),t,nt,dt)
                    mout(iz,ix,iy)=v
                    t=t+dz/v
                enddo
            enddo; enddo

        endif

    end subroutine

    ! convert the gradient from depth domain to time domain
    subroutine pseudotime_z2t_gradient(gin,v_z,v_t,gout)
        real,dimension(nz,nx,ny) :: gin, v_z
        real,dimension(nt,nx,ny) :: v_t
        real,dimension(:,:,:),allocatable :: gout
        real,dimension(nz) :: dvdz

        call alloc(gout,nt,nx,ny)

        do iy=1,ny; do ix=1,nx
            !compute dv/dz
            !using central finite different at interior nodes
            do iz=2,nz-1
               dvdz(iz)=(v_z(iz+1,ix,iy)-v_z(iz-1,ix,iy))/2/dz
            enddo
            !using one-side finite different at boundary nodes
            dvdz( 1)=(v_z( 2,ix,iy)-v_z(   1,ix,iy))/dz
            dvdz(nz)=(v_z(nz,ix,iy)-v_z(nz-1,ix,iy))/dz

            zn=0.
            do it=1,nt

                iz1=min(  floor(zn/dz)+1,   nz) !maximum lower index
                iz2=min(ceiling(zn/dz)+1,   nz) !minimum upper index

                z1=zn-(iz1-1)*dz !distance inside interval
                z2=(iz2-1)*dz-zn !distance inside interval

                !init value, assuming the gradient is linear inside interval
                if (iz1==iz2) then
                    g_interp = gin(iz2,ix,iy)
                else
                    g_interp = (gin(iz2,ix,iy)*z1 + gin(iz1,ix,iy)*z2)/dz
                endif

                ! floor v ceiling selecting, too hash..
                ! g_interp=either(gin(iz1,ix), gin(iz2,ix), z1<z2)
                
                gout(it,ix,iy) = g_interp
              
              
                !integral from zn to (iz2-1)*dz
                !no computation as the integral should be very small
                gout(it,ix,iy) = gout(it,ix,iy) - 0.


                !integral from (iz2-1)*dz to (nz-1)*dz       
                ! rectangular approx to the integral
                s = sum(gin(iz2:nz,ix,iy)/v_z(iz2:nz,ix,iy)*dvdz(iz2:nz))*dz

                ! ! trapezoidal approx to the integral, not used as it is no better than above.
                ! s = sum(gin(iz2:nz,ix,iy)*dvdz(iz2:nz)/v_z(iz2:nz,ix,iy))
                ! s = 2*s - gin(iz2,ix,iy)*dvdz(iz2)/v_z(iz2,ix,iy) -  gin(nz,ix,iy)*dvdz(nz)/v_z(nz,ix,iy)
                ! s = s * (nz-iz2)*dz/2

                gout(it,ix,iy) = gout(it,ix,iy) - s

                zn=zn+v_t(it,ix,iy)*dt

            enddo

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
