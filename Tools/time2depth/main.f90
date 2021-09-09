! Conversion :
!   Axes:
!     Depth -> pseudo-time:      t = int_0^z dz' / v(z')
!     Pseudo-time -> depth:      z = int_0^t v(t') dt'
!   Gradient:
!     Depth -> pseudo-time:
!       dC/dv(tn) = dC/dv(zn) - int_zn^zmax dC/dv(zp) * 1/v(zp) dv(zp)/zp dzp
!         where zn= int_0^tn v(t') dt'
!         Note that v(tn) may be different from v(zn) or v(zp)

module m_pseudotime

    integer :: nz,nt,nx
    real :: dz,dt

    contains

    subroutine pseudotime_init(dir,vmin,vmax,nx_,nz_,dz_,nt_,dt_)
        character(4) :: dir

        nx=nx_

        !for time->depth
        if(dir=='t->z') then
            dt=dt_;  nt=nt_

            !dz = dt*v(z) >= dt*vmin
            !z = int_0^T v(t')*dt' <= (nt-1)*dt vmax 
            !nz = zmax/dz
            dz=dt*vmin
            zmax=(nt-1)*dt*vmax
            nz=ceiling(zmax/dz)+1

            write(*,*) 'dz,nz = ',dz,nz
            dz_=dz;  nz_=nz
        endif

        !for depth->time
        if(dir=='z->t') then
            dz=dz_;  nz=nz_

            !dt = dz/v(z) >= dz/vmax
            dt=dz/vmax
            !t = int_0^Z dz'/v(z') <= (nz-1)*dz /vmin
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

            write(*,*) 'dt,nt = ',dt,nt
            dt_=dt;  nt_=nt
        endif
      
    end subroutine

    ! convert the model from depth domain to time domain
    subroutine pseudotime_z2t(vin,vout)
        real,dimension(nz,nx) :: vin
        real,dimension(:,:),allocatable :: vout
    
        allocate(vout(nt,nx))

        do ix=1,nx
            z=0.
            do it=1,nt
                v=interp(vin(:,ix),z,nz,dz)
                vout(it,ix)=v
                z=z+v*dt
            enddo
        enddo

    end subroutine

    ! convert the model from depth domain to time domain
    subroutine pseudotime_z2t_reflectivity(vin,rin,rout)
        real,dimension(nz,nx) :: vin, rin
        real,dimension(:,:),allocatable :: rout

        allocate(rout(nt,nx))
     
        do ix=1,nx
            z=0.
            do it=1,nt
                rout(it,ix)=interp(rin(:,ix),z,nz,dz)
                z=z+interp(vin(:,ix),z,nz,dz)*dt
            enddo
        enddo

    end subroutine


    ! calculates the water bottom in the time domain
    ! enables to mute the time-domain gradient in the bathymetry
    function pseudotime_z2t_water(vw,iz) result(it)
        it= nint(   (iz-1)*dz/vw   /dt) +1
    end function


    ! convert the model from time domain to depth domain
    subroutine pseudotime_t2z(vin,vout)
        real,dimension(nt,nx) :: vin
        real,dimension(:,:),allocatable :: vout
         
        allocate(vout(nz,nx))

        do ix=1,nx
            t=0.
            do iz=1,nz
                v=interp(vin(:,ix),t,nt,dt)
                vout(iz,ix)=v
                t=t+dz/v
            enddo
        enddo

    end subroutine

    ! convert the model from time domain to depth domain
    subroutine pseudotime_t2z_reflectivity(vin,rin,rout)
        real,dimension(nt,nx) :: vin, rin
        real,dimension(:,:),allocatable :: rout

        allocate(rout(nz,nx))
         
        do ix=1,nx
            t=0.
            do iz=1,nz
               rout(iz,ix)=interp(rin(:,ix),t,nt,dt)
               t=t+dz/interp(vin(:,ix),t,nt,dt)
            enddo
        enddo

    end subroutine

    ! convert the gradient from depth domain to time domain
    subroutine pseudotime_z2t_gradient(gin,vdepth,vtime,gout)
        real,dimension(nz,nx) :: gin,  vdepth
        real,dimension(nt,nx) :: vtime
        real,dimension(:,:),allocatable :: gout
        real,dimension(nz) :: dcdz

        allocate(gout(nt,nx))

        do ix=1,nx
            !compute dc/dz
            !using central finite different at interior nodes
            do iz=2,nz-1
               dcdz(iz)=(vdepth(iz+1,ix)-vdepth(iz-1,ix))/2/dz
            enddo
            !using one-side finite different at boundary nodes
            dcdz(1) =(vdepth(2,ix) -vdepth(1,ix))   /dz
            dcdz(nz)=(vdepth(nz,ix)-vdepth(nz-1,ix))/dz

            zn=0.
            do it=1,nt

                iz1=min(  floor(zn/dz)+1,   nz) !maximum lower index
                iz2=min(ceiling(zn/dz)+1,   nz) !minimum upper index

                z1=zn-(iz1-1)*dz !distance inside interval
                z2=(iz2-1)*dz-zn !distance inside interval

                !init value, assuming the gradient is linear inside interval
                if (iz1==iz2) then
                    g_interp = gin(iz2,ix)
                else
                    g_interp = gin(iz2,ix)*z1 + gin(iz1,ix)*z2
                    g_interp = g_interp/dz
                endif

                ! floor v ceiling selecting, too hash..
                ! g_interp=either(gin(iz1,ix), gin(iz2,ix), z1<z2)
                
                gout(it,ix) = g_interp
              
              
                !integral from zn to (iz2-1)*dz
                !no computation as the integral should be very small
                gout(it,ix) = gout(it,ix) - 0.


                !integral from (iz2-1)*dz to (nz-1)*dz       
                ! rectangular approx to the integral
                s = sum(gin(iz2:nz,ix)/vdepth(iz2:nz,ix)*dcdz(iz2:nz))*dz

                ! ! trapezoidal approx to the integral, not used as it is no better than above.
                ! s = sum(gin(iz2:nz,ix)*dcdz(iz2:nz)/vdepth(iz2:nz,ix))
                ! s = 2*s - gin(iz2,ix)*dcdz(iz2)/vdepth(iz2,ix) -  gin(nz,ix)*dcdz(nz)/vdepth(nz,ix)
                ! s = s * (nz-iz2)*dz/2

                gout(it,ix) = gout(it,ix) - s

                zn=zn+vtime(it,ix)*dt

            enddo

        enddo

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

program main
use m_pseudotime

    character(120) :: filename
    integer :: file_size
    logical :: exist
    real,dimension(:,:),allocatable :: vtime
    real,dimension(:,:),allocatable :: vdepth

    write(*,*) 'Enter filename of the input time-domain model:'
    read(*,*) filename

    inquire(file=filename,size=file_size,exist=exist)
    if(file_size==0) exist=.false.
    if(.not.exist) then
        write(*,*) 'ERROR: file '//trim(adjustl(filename))//' does NOT exist or has 0 length.'
        stop
    endif

    write(*,*) '1st dim of the model (nt, dt):'
    read(*,*) nt,dt

    nx=file_size/4/nt
    write(*,*) 'Shape of the model (nt x nx):'
    write(*,*) nt,nx

    allocate(vtime(nt,nx))
    open(11,file=filename,access='direct',recl=4*nt*nx)
    read(11,rec=1)vtime
    close(11)

    call pseudotime_init('t->z',vmin=minval(vtime),vmax=maxval(vtime),nx_=nx,&
        nz_=nz, dz_=dz, &
        nt_=nt, dt_=dt)

    call pseudotime_t2z(vtime,vdepth)

    deallocate(vtime)

    write(*,*) 'Shape of the output depth-domain model (nz x nx):'
    write(*,*) nz,nx
    
    write(*,*) 'Sampling on the depth axis (dz):'
    write(*,*) dz
    
    write(*,*) 'Enter filename of the output depth-domain model:'
    read(*,*) filename

    open(11,file=filename,access='direct',recl=4*size(vdepth))
    write(11,rec=1) vdepth
    close(11)

end program
