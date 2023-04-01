module m_smoother_Laplacian_sparse
use m_System

    private
    public :: smoother_Laplacian_init, smoother_Laplacian_extend_mirror, smoother_Laplacian_pseudo_stationary, smoother_Laplacian_pseudo_nonstationary
    
    integer :: nz,nx,ny, iaddmirror
    real    :: dz,dx,dy
    real    :: freq
    real    :: frac_z,frac_x,frac_y
    logical :: is_cubic
    character(6) :: preserve

    contains
    
    
    subroutine smoother_1D_stationary(n, h, vector, b)
        integer n  !vector length
        real    h  !sampling rate
        real,dimension(n) :: vector !array to be smoothed
        real    b !scale, variance = 2b^2
        
        real,dimension(n) :: d !diagonal terms
        real,dimension(n) :: du,dl  !super- and sub-diagonal terms
        
        tmp=1./(exp(h/b)-exp(-h/b))
        
        du=   -tmp
        d = 1.+2.*exp(-h/b)*tmp
        dl=   -tmp
        
        !Neumann boundary conditions on two ends to avoid leakage after the smoothing
        d(1)=d(1)+dl(1)
        d(n)=d(n)+du(n)
        
        call sgtsv( n, 1, dl(2:n), d(:), du(1:n-1), vector, n, info )
        
        if(info /= 0 ) then
            if(info<0) then
                call error('Laplacian smoothing has a problem: -i='//num2str(info)//', The i-th argument had an illegal value')
            else
                call error('Laplacian smoothing has a problem: i='//num2str(info)//', U(i,i) is exactly zero, and the solution has not been computed. The factorization has not been completed unless i = N.')
            endif
        endif
        
        !normalization
        vector=vector/b*0.39894228
        
    end subroutine
    
    subroutine smoother_1D_nonstationary(n, h, vector, frac_lambda)
        integer n  !vector length
        real    h  !sampling rate
        real,dimension(n) :: vector, frac_lambda  !array to be smoothed and fraction of wavelength
        
        real :: b !Laplacian function parameter
        real,dimension(n) :: d !diagonal terms
        real,dimension(n) :: du,dl  !super- and sub-diagonal terms
        
        do i=1,n
            b=frac_lambda(i)/2.
            
            tmp=1./(exp(h/b)-exp(-h/b))
            
            du(i)=   -tmp
            d(i) = 1.+2.*exp(-h/b)*tmp
            dl(i)=   -tmp
        enddo
        
        !Neumann boundary conditions on two ends to avoid edge effect after the smoothing
        d(1)=d(1)+dl(1)
        d(n)=d(n)+du(n)
        
        call sgtsv( n, 1, dl(2:n), d(:), du(1:n-1), vector, n, info )
        
        if(info /= 0 ) then
            if(info<0) then
                call error('Laplacian smoothing has a problem: -i='//num2str(info)//', The i-th argument had an illegal value')
            else
                call error('Laplacian smoothing has a problem: i='//num2str(info)//', U(i,i) is exactly zero, and the solution has not been computed. The factorization has not been completed unless i = N.')
            endif
        endif
        
        !normalization
        do i=1,n
           b=frac_lambda(i)/2.
           vector(i)=vector(i)/b*0.39894228
        enddo
        
    end subroutine
    
    
    subroutine smoother_Laplacian_init(n,d,freq_,o_addmirror,o_frac,o_preserve)
        integer,dimension(3) :: n
        real,dimension(3) :: d
        
        real,optional :: o_addmirror
        real,dimension(3),optional :: o_frac
        character(6),optional :: o_preserve
        
        real,dimension(3) :: tmp
        
        nz=n(1); dz=d(1)
        nx=n(2); dx=d(2)
        ny=n(3); dy=d(3)
        
        is_cubic=either(.true.,.false.,ny>1)
        
        freq=freq_

        iaddmirror = nint( &
            either(o_addmirror, &
                setup%get_real('SMTH_ADDMIRROR',o_default=num2str(dz)), &
                present(o_addmirror)) &
            /dz)+1
    
        if(present(o_frac)) then
            frac_z=o_frac(1)
            frac_x=o_frac(2)
            frac_y=o_frac(3)
        else
            tmp=setup%get_reals('LAPLACIAN_SMTH_FRACTION_WAVELENGTH','LAP_FRAC',o_default='1 1 1')
            frac_z=tmp(1)
            frac_x=tmp(2)
            frac_y=tmp(3)
        endif
        
        preserve = either(o_preserve, &
            setup%get_str('LAPLACIAN_SMTH_PRESERVE_GRAD_NORM','LAP_PRESERVE',o_default='nopreserve'), &
            present(o_preserve))

    end subroutine
    
    subroutine smoother_Laplacian_extend_mirror(grad,itopo)
    !to avoid leakage into the mask zone after smoothing
    !this subroutine mirrors the gradient wrt midway between grid points #1 & #2(itopo+iaddmirror-1) in depth
    !then another subroutine cleans the mask zone later
        real,dimension(nz,nx,ny) :: grad
        integer,dimension(nx,ny) :: itopo
    
        do iy=1,ny
        do ix=1,nx
            grad(itopo(ix,iy)+iaddmirror-1 :1                            :-1, ix,iy)=&
            grad(itopo(ix,iy)+iaddmirror   :2*(itopo(ix,iy)+iaddmirror-1)   , ix,iy)
        enddo
        enddo
    
    end subroutine
    
    !Use convolutions of 1D Laplacian functions on x,y,z axes to roughly approx the 3D Laplacian function
    !Note that high-dimensional Laplacian function can't be decomposed as convolutions of 1D Laplacian functions, unlike Gaussian functions..
    !True 2D Laplacian function has a circular shape whereas this approx function has a diamond shape.   
    subroutine smoother_Laplacian_pseudo_stationary(grad,b)
        real,dimension(nz,nx,ny) :: grad
        real,dimension(3) :: b
        
        !real norm
        real,dimension(:),allocatable :: vector
        
        !norm=sum(abs(grad))
        
        !z axis
        allocate(vector(nz))
        do iy=1,ny; do ix=1,nx
            vector=grad(:,ix,iy)
            call smoother_1D_stationary(nz,dz,vector,b(1))
            grad(:,ix,iy)=vector
        enddo; enddo
        deallocate(vector)
        
        !x axis
        allocate(vector(nx))
        do iy=1,ny; do iz=1,nz
            vector=grad(iz,:,iy)
            call smoother_1D_stationary(nx,dx,vector,b(2))
            grad(iz,:,iy)=vector
        enddo; enddo
        deallocate(vector)
        
        !y axis
        if(is_cubic) then
            allocate(vector(ny))
            do ix=1,nx; do iz=1,nz
                vector=grad(iz,ix,:)
                call smoother_1D_stationary(ny,dy,vector,b(3))
                grad(iz,ix,:)=vector
            enddo; enddo
            deallocate(vector)
        endif
        
    end subroutine
    
    subroutine smoother_Laplacian_pseudo_nonstationary(grad,velocity)
        real,dimension(nz,nx,ny) :: grad,velocity
        
        real,dimension(:),allocatable :: vector,frac_lambda
        real :: old_norm
        real,dimension(:,:,:),allocatable :: table_one
        
        select case (preserve)
        case('L1norm')
            old_norm = sum(abs(grad))
        case('L2norm')
            old_norm = norm2(grad)
        case('scale1')
            allocate(table_one(nz,nx,ny))
            table_one=1.
        end select
        
        !z axis
        allocate(vector(nz),frac_lambda(nz))
        do iy=1,ny
        do ix=1,nx
            vector=grad(:,ix,iy)
            frac_lambda=velocity(:,ix,iy)/freq*frac_z
            call smoother_1D_nonstationary(nz,dz,vector,frac_lambda)
            grad(:,ix,iy)=vector
        enddo
        enddo
        deallocate(vector,frac_lambda)
        
        !x axis
        allocate(vector(nx),frac_lambda(nx))
        do iy=1,ny
        do iz=1,nz
            vector=grad(iz,:,iy)
            frac_lambda=velocity(iz,:,iy)/freq*frac_x
            call smoother_1D_nonstationary(nx,dx,vector,frac_lambda)
            grad(iz,:,iy)=vector
        enddo
        enddo
        deallocate(vector,frac_lambda)
        
        !y axis
        if(is_cubic) then
            allocate(vector(ny),frac_lambda(ny))
            do ix=1,nx
            do iz=1,nz
                vector=grad(iz,ix,:)
                frac_lambda=velocity(iz,ix,:)/freq*frac_y
                call smoother_1D_nonstationary(ny,dy,vector,frac_lambda)
                grad(iz,ix,:)=vector
            enddo
            enddo
            deallocate(vector,frac_lambda)
        endif
        
        !redo smoothing on table_one for normalization
        if(preserve=='scale1') then
            !z axis
            allocate(vector(nz),frac_lambda(nz))
            do iy=1,ny
            do ix=1,nx
                vector=table_one(:,ix,iy)
                frac_lambda=velocity(:,ix,iy)/freq*frac_z
                call smoother_1D_nonstationary(nz,dz,vector,frac_lambda)
                table_one(:,ix,iy)=vector
            enddo
            enddo
            deallocate(vector,frac_lambda)
            
            !x axis
            allocate(vector(nx),frac_lambda(nx))
            do iy=1,ny
            do iz=1,nz
                vector=table_one(iz,:,iy)
                frac_lambda=velocity(iz,:,iy)/freq*frac_x
                call smoother_1D_nonstationary(nx,dx,vector,frac_lambda)
                table_one(iz,:,iy)=vector
            enddo
            enddo
            deallocate(vector,frac_lambda)
            
            !y axis
            if(is_cubic) then
                allocate(vector(ny),frac_lambda(ny))
                do ix=1,nx
                do iz=1,nz
                    vector=table_one(iz,ix,:)
                    frac_lambda=velocity(iz,ix,:)/freq*frac_y
                    call smoother_1D_nonstationary(ny,dy,vector,frac_lambda)
                    table_one(iz,ix,:)=vector
                enddo
                enddo
                deallocate(vector,frac_lambda)
            endif
        endif
        
        !normalization
        select case (preserve)
        case('L1norm')
            grad = grad * old_norm / sum(abs(grad))
        case('L2norm')
            grad = grad * old_norm / norm2(grad)
        case('scale1')
            grad = grad / table_one
            deallocate(table_one)
        end select
        
    end subroutine
    
end
