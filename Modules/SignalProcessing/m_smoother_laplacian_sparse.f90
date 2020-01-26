module m_laplacian_smoothing_sparse
use m_sysio

    private pi,nz,nx,ny,iaddmirror,dz,dx,dy,freq,frac_z,frac_x,frac_y,is_cubic,preserve
    public
    
    real,parameter :: pi=3.1415927
    
    integer :: nz,nx,ny, iaddmirror
    real    :: dz,dx,dy
    real    :: freq
    real    :: frac_z,frac_x,frac_y
    logical :: is_cubic
    character(6) :: preserve

    contains
    
    
    subroutine smooth_1D_stationary(n, h, vector, b)
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
            !call hud('Laplacian smoothing has a problem')
            if(info<0) then
                write(*,*) '-i=',info, ': The i-th argument had an illegal value'
            else
                write(*,*) 'i=',info, ': U(i,i) is exactly zero, and the solution has not been computed. The factorization has not been completed unless i = N.'
            endif
            stop 'stopped now'
        endif
        
        !normalization
        vector=vector/b*0.39894228
        
    end subroutine
    
    subroutine smooth_1D_nonstationary(n, h, vector, frac_lambda)
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
                write(*,*) '-i=',info, ': The i-th argument had an illegal value'
            else
                write(*,*) 'i=',info, ': U(i,i) is exactly zero, and the solution has not been computed. The factorization has not been completed unless i = N.'
            endif
            stop 'stopped now'
        endif
        
        !normalization
        do i=1,n
           b=frac_lambda(i)/2.
           vector(i)=vector(i)/b*0.39894228
        enddo
        
    end subroutine
    
    
    subroutine init_laplacian_smoothing(n,d,frequency,o_addmirror,o_frac,o_preserve)
        integer,dimension(3) :: n
        real,dimension(3) :: d
        
        real,optional :: o_addmirror
        real,dimension(3),optional :: o_frac
        character(6),optional :: o_preserve
        
        character(:),allocatable :: tmp
        
        nz=n(1); dz=d(1)
        nx=n(2); dx=d(2)
        ny=n(3); dy=d(3)
        
        if(ny==1) then
            is_cubic=.false.
        else
            is_cubic=.true.
        endif
        
        freq=frequency
        
        if(present(o_addmirror)) then
            iaddmirror=nint(o_addmirror/dz)+1
        else
            iaddmirror=nint( get_setup_real('SMOOTHING_ADDMIRROR',default=dz) /dz ) +1
        endif
        
        if(present(o_frac)) then
            frac_z=o_frac(1)
            frac_x=o_frac(2)
            frac_y=o_frac(3)
        else
            tmp=get_setup_char('SMOOTH_GRADIENT_WAVELENGTH_FRACTION',default='1 1 1')
            read(tmp,*)frac_z,frac_x,frac_y
        endif
        
        if(present(o_preserve)) then
            preserve=o_preserve
        else
            preserve=get_setup_char('SMOOTH_GRADIENT_PRESERVE_MAGNITUDE',default='nopreserve')
        endif
        
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
        
    end subroutine
    
    subroutine laplacian_smoothing_pseudo_nonstationary(gradient,velocity)
        real,dimension(nz,nx,ny) :: gradient,velocity
        
        real,dimension(:),allocatable :: vector,frac_lambda
        real :: old_norm
        real,dimension(:,:,:),allocatable :: table_one
        
        if(preserve=='L1norm') then
            old_norm = sum(abs(gradient))
        elseif(preserve=='L2norm') then
            old_norm = norm2(gradient)
        elseif(preserve=='scale1') then
            allocate(table_one(nz,nx,ny))
            table_one=1.
        endif
        
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
        
        !redo smoothing on table_one for normalization
        if(preserve=='scale1') then
            !z axis
            allocate(vector(nz),frac_lambda(nz))
            do iy=1,ny
            do ix=1,nx
                vector=table_one(:,ix,iy)
                frac_lambda=velocity(:,ix,iy)/freq*frac_z
                call smooth_1D_nonstationary(nz,dz,vector,frac_lambda)
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
                call smooth_1D_nonstationary(nx,dx,vector,frac_lambda)
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
                    call smooth_1D_nonstationary(ny,dy,vector,frac_lambda)
                    table_one(iz,ix,:)=vector
                enddo
                enddo
                deallocate(vector,frac_lambda)
            endif
        endif
        
        !normalization
        if(preserve=='L1norm') then
            gradient = gradient * old_norm / sum(abs(gradient))
        elseif(preserve=='L2norm') then
            gradient = gradient * old_norm / norm2(gradient)
        elseif(preserve=='scale1') then
            gradient = gradient / table_one
            deallocate(table_one)
        endif
        
    end subroutine
    
end
