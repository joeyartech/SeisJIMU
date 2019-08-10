module m_laplacian_smoothing_sparse
use m_sysio

    private pi,nz,nx,ny,iaddmirror,dz,dx,dy,freq,frac_z,frac_x,frac_y,is_cubic,scaling
    public
    
    real,parameter :: pi=3.1415927
    
    integer :: nz,nx,ny, iaddmirror
    real    :: dz,dx,dy
    real    :: freq
    real    :: frac_z,frac_x,frac_y
    logical :: is_cubic
    
    real,parameter :: scaling=0.2 !scaling factor to convert fraction of half wavelength to standard deviation. This magic number assures the gradient is multiplified by ~0.001 by the laplacian or gaussian functions half frac_wavelength away from the center (note that these two functions do not have max=1 due to normalization)

    contains
    
    !analytical version, preserves vector's energy after smoothing but can amplify values at the end (e.g. could be dangerous at the seafloor)
!     subroutine smooth_1D_stationary(n, h, vector, b)
!         integer n  !vector length
!         real    h  !sampling rate
!         real,dimension(n) :: vector !array to be smoothed
!         
!         real    b !scale, variance = 2b^2
!         real,dimension(n) :: d !diagonal terms
!         real,dimension(n) :: du,dl  !super- and sub-diagonal terms
!         
!         du=-b/h**3
!         d=1./b/h+2*b/h**3
!         dl=-b/h**3
! 
!         !Neumann boundary conditions on two ends to avoid leakage after the smoothing
!         d(1)=d(1)-b/h**3
!         d(n)=d(n)-b/h**3
!         
!         call sgtsv( n, 1, dl(2:n), d(:), du(1:n-1), vector, n, info )
!         
!         renormalization
!         vector=vector*2/h**2
!         
!         if(info /= 0 ) then
!             !call hud('Laplacian smoothing has a problem')
!             if(info<0) then
!                 write(*,*) '-i=',info, ': The i-th argument had an illegal value'
!             else
!                 write(*,*) 'i=',info, ': U(i,i) is exactly zero, and the solution has not been computed. The factorization has not been completed unless i = N.'
!             endif
!             stop 'stopped now'
!         endif
!         
!         !renormalization
!         vector=vector*0.5*h/b
!         
!     end subroutine
    
    !`damped version', no preservation of energy but quite precise to the discretized Laplacian function and no value amplications at the end (e.g. safe at the seafloor)
    subroutine smooth_1D_stationary(n, h, vector, b)
        integer n  !vector length
        real    h  !sampling rate
        real,dimension(n) :: vector !array to be smoothed
        real    b !scale, variance = 2b^2
        
        real,dimension(n) :: d !diagonal terms
        real,dimension(n) :: du,dl  !super- and sub-diagonal terms
        
        tmp=1./(exp(h/b)-exp(-h/b))
        
        du=-1.*tmp
        d = 1.+2.*exp(-h/b)*tmp
        dl=-1.*tmp
        
        !Neumann boundary conditions on two ends to avoid leakage after the smoothing
        d(1)=d(1)-b/h**3
        d(n)=d(n)-b/h**3
        
        call sgtsv( n, 1, dl(2:n), d(:), du(1:n-1), vector, n, info )
        
        if(info /= 0 ) then
            !call hud('Laplacian smoothing has a problem')
            if(info<0) then
                write(*,*) '-i=',info, ': The i-th argument had an illegal value'
            else
                write(*,*) 'i=',info, ': U(i,i) is exactly zero, and the solution has not been computed. The factorization has not been completed unless i = N.'
            endif
            stop 'stopped now'
        endif
        
        !normalization
        vector=vector*0.5*h/b
        
    end subroutine
    
    subroutine smooth_1D_nonstationary(n, h, vector, frac_lambda)
        integer n  !vector length
        real    h  !sampling rate
        real,dimension(n) :: vector, frac_lambda  !array to be smoothed and fraction of wavelength
        
        real :: b !Laplacian function parameter
        real,dimension(n) :: d !diagonal terms
        real,dimension(n) :: du,dl  !super- and sub-diagonal terms
        
        do i=1,n
            b=scaling*frac_lambda(i)/sqrt(2.)  !b=std/sqrt(2)
            
            tmp=1./(exp(h/b)-exp(-h/b))
            
            du(i)=-1.*tmp
            d(i) = 1.+2.*exp(-h/b)*tmp
            dl(i)=-1.*tmp
        enddo
        
        !Neumann boundary conditions on two ends to avoid edge effect after the smoothing
        d(1)=d(1)-dl(1)
        d(n)=d(n)-du(n)
        
        call sgtsv( n, 1, dl(2:n), d(:), du(1:n-1), vector, n, info )
        
        if(info /= 0 ) then
            !call hud('Gaussian smoothing problem along z direction')
            if(info<0) then
                write(*,*) '-i=',info, ': The i-th argument had an illegal value'
            else
                write(*,*) 'i=',info, ': U(i,i) is exactly zero, and the solution has not been computed. The factorization has not been completed unless i = N.'
            endif
            stop 'stopped now'
        endif
        
        !normalization
        do i=1,n
           b=scaling*frac_lambda(i)/sqrt(2.)  !b=std/sqrt(2)
           vector(i)=vector(i)*0.5*h/b
        enddo
        
    end subroutine
    
    
    subroutine init_laplacian_smoothing(n,d,cubic,frequency)
        integer,dimension(3) :: n
        real,dimension(3) :: d
        logical :: cubic
        
        character(:),allocatable :: tmp
        
        nz=n(1); dz=d(1)
        nx=n(2); dx=d(2)
        ny=n(3); dy=d(3)
        
        is_cubic=cubic
        
        freq=frequency
        
        call hud('Initialize Laplacian smoothing')
        iaddmirror=nint( get_setup_real('SMOOTHING_ADDMIRROR',default=dz) /dz ) +1
        
        tmp=get_setup_char('SMOOTH_GRADIENT_WAVELENGTH_FRACTION',default='1 1 1')
        read(tmp,*)frac_z,frac_x,frac_y
        
    end subroutine
    
    subroutine laplacian_smoothing_extend_mirror(gradient,itopo)
    !to avoid leakage into the mask zone after smoothing
    !this subroutine mirrors the gradient wrt midway between points itopo+iaddmirror-1 & itopo+iaddmirror-1 in depth
    !then another subroutine cleans the mask zone later
        real,dimension(nz,nx,ny) :: gradient
        integer,dimension(nx,ny) :: itopo
    
        do iy=1,ny
        do ix=1,nx
            gradient(itopo(ix,iy)+iaddmirror-1 :1                            :-1, ix,iy)=&
            gradient(itopo(ix,iy)+iaddmirror   :2*(itopo(ix,iy)+iaddmirror-1)   , ix,iy)
        enddo
        enddo
    
    end subroutine
    
    !Use convolutions of 1D Laplacian functions on x,y,z axes to roughly approx the 3D Laplacian function
    !Note that high-dimensional Laplacian function can't be decomposed as convolutions of 1D Laplacian functions, unlike Gaussian functions..
    !True 2D Laplacian function has a circular shape whereas this approx function has a diamond shape.   
    subroutine laplacian_smoothing_pseudo_stationary(gradient,b)
        real,dimension(nz,nx,ny) :: gradient
        real,dimension(3) :: b
        
        !real norm
        real,dimension(:),allocatable :: vector
        
        !norm=sum(abs(gradient))
        
        !z axis
        allocate(vector(nz))
        do iy=1,ny
        do ix=1,nx
            vector=gradient(:,ix,iy)
            call smooth_1D_stationary(nz,dz,vector,b(1))
            gradient(:,ix,iy)=vector
        enddo
        enddo
        deallocate(vector)
        
        !x axis
        allocate(vector(nx))
        do iy=1,ny
        do iz=1,nz
            vector=gradient(iz,:,iy)
            call smooth_1D_stationary(nx,dx,vector,b(2))
            gradient(iz,:,iy)=vector
        enddo
        enddo
        deallocate(vector)
        
        !y axis
        if(is_cubic) then
            allocate(vector(ny))
            do ix=1,nx
            do iz=1,nz
                vector=gradient(iz,ix,:)
                call smooth_1D_stationary(ny,dy,vector,b(3))
                gradient(iz,ix,:)=vector
            enddo
            enddo
            deallocate(vector)
        endif
        
        !renormalization
        if(is_cubic) then
            gradient=gradient/pi !due to approximation this theoretical scaling factor is no more valid..
            !gradient=gradient*norm/sum(abs(gradient))
        else
            gradient=gradient*2./pi !due to approximation this theoretical scaling factor is no more valid..
            !gradient=gradient*norm/sum(abs(gradient))
        endif
        
    end subroutine
    
    subroutine laplacian_smoothing_pseudo_nonstationary(gradient,velocity)
        real,dimension(nz,nx,ny) :: gradient,velocity
        
        real,dimension(:),allocatable :: vector,frac_lambda
        
        !z axis
        allocate(vector(nz),frac_lambda(nz))
        do iy=1,ny
        do ix=1,nx
            vector=gradient(:,ix,iy)
            frac_lambda=velocity(:,ix,iy)/freq*frac_z
            call smooth_1D_nonstationary(nz,dz,vector,frac_lambda)
            gradient(:,ix,iy)=vector
        enddo
        enddo
        deallocate(vector,frac_lambda)
        
        !x axis
        allocate(vector(nx),frac_lambda(nx))
        do iy=1,ny
        do iz=1,nz
            vector=gradient(iz,:,iy)
            frac_lambda=velocity(iz,:,iy)/freq*frac_x
            call smooth_1D_nonstationary(nx,dx,vector,frac_lambda)
            gradient(iz,:,iy)=vector
        enddo
        enddo
        deallocate(vector,frac_lambda)
        
        !y axis
        if(is_cubic) then
            allocate(vector(ny),frac_lambda(ny))
            do ix=1,nx
            do iz=1,nz
                vector=gradient(iz,ix,:)
                frac_lambda=velocity(iz,ix,:)/freq*frac_y
                call smooth_1D_nonstationary(ny,dy,vector,frac_lambda)
                gradient(iz,ix,:)=vector
            enddo
            enddo
            deallocate(vector,frac_lambda)
        endif
        
        !renormalization
        if(is_cubic) then
            gradient=gradient/pi
        else
            gradient=gradient*pi
        endif
        
    end subroutine
    
end
