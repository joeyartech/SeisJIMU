module m_hicks
use m_arrayop

    private b_monopole,err_monopole,b_dipole,err_dipole
    private r,b,epsilon
    
    !Optimum Kaiser windowing parameter b for monopole point source
    !half window width r=1~10, kmax=pi/2 or 2pi/3
    real,dimension(2,10),parameter :: b_monopole = reshape([ &
        1.24,    2.94,    4.53,    6.31,    7.91,    9.42,   10.95,   12.53,   14.09,   14.18, & !kmax=pi/2
        0.00,    1.84,    3.04,    4.14,    5.26,    6.40,    7.51,    8.56,    9.56,   10.64  & !kmax=2pi/3
        ],[2,10])
    !Corresponding spectral error
    real,dimension(2,10),parameter :: err_monopole = reshape([ &
        1.72e-1, 2.32e-2, 4.79e-3, 1.34e-3, 2.43e-4, 4.88e-5, 9.79e-6, 2.19e-6, 8.23e-7, 7.99e-7, & !kmax=pi/2
        3.52e-1, 8.46e-2, 2.50e-2, 7.97e-3, 2.83e-3, 1.04e-3, 3.49e-4, 1.03e-4, 3.53e-5, 1.14e-5  & !kmax=2pi/3
        ],[2,10])
        
    !Optimum Kaiser windowing parameter b for dipole point source
    !which is not needed in the velocity/stress formulation
    !half window width r=1~10, kmax=pi/2 or 2pi/3
    real,dimension(2,10),parameter :: b_dipole = reshape([ &
        0.00,    3.33,    4.96,    6.42,    7.77,    9.52,    11.11,   12.52,   14.25,   16.09, & !kmax=pi/2
        0.00,    1.48,    3.25,    4.40,    5.44,    6.46,     7.41,    8.50,    9.66,   10.77  & !kmax=2pi/3
        ],[2,10])
    real,dimension(2,10),parameter :: err_dipole = reshape([ &
        1.57e-0, 1.97e-1, 3.39e-2, 6.81e-3, 1.52e-3, 2.46e-4, 5.29e-5, 1.40e-5, 1.06e-5, 1.02e-5, & !kmax=pi/2
        2.07e-0, 5.07e-1, 1.46e-1, 4.36e-2, 1.35e-2, 4.40e-3, 1.55e-3, 4.95e-4, 1.63e-4, 6.00e-5  & !kmax=2pi/3
        ],[2,10])
    
    
    !choose r, kmax, b
    integer,parameter :: r=4
    real,parameter :: b=b_monopole(1,r) !kmax=pi/2.
    
    !threshold for small number
    real,parameter :: epsilon=1e-3
    
    !variable carrier
    type t_hicks
        !intent(in)
        real :: x,y,z  !source point location
        real :: dx,dy,dz
        logical :: is_cubic,if_freesurface
        !intent(out)
        !range of interpolation points
        integer :: ifx,ify,ifz
        integer :: ix,iy,iz
        integer :: ilx,ily,ilz
    end type
    
    type(t_hicks) :: hicks
    
    contains

    subroutine build_hicks(h,fold_type,interp_coeff)
        type(t_hicks) :: h
        character(*) :: fold_type
        real,dimension(:,:,:),allocatable :: interp_coeff
        
        real,dimension(-r:r) :: x_interp_coeff,y_interp_coeff,z_interp_coeff
        
        x_interp_coeff=0.
        y_interp_coeff=0.
        z_interp_coeff=0.
        
        !x axis
        call build_1d_hicks(h%x,h%dx,h%ifx,h%ix,h%ilx,x_interp_coeff)
        !y axis
        call build_1d_hicks(h%y,h%dy,h%ify,h%iy,h%ily,y_interp_coeff)
        !z axis
        call build_1d_hicks(h%z,h%dz,h%ifz,h%iz,h%ilz,z_interp_coeff)
        
        !index starts from 1 outside scope of build_1d_hicks
        h%ifx=h%ifx+1; h%ix=h%ix+1; h%ilx=h%ilx+1
        h%ify=h%ify+1; h%iy=h%iy+1; h%ily=h%ily+1
        if(.not.h%is_cubic) then
            h%ify=1; h%iy=1; h%ily=1
        endif
        h%ifz=h%ifz+1; h%iz=h%iz+1; h%ilz=h%ilz+1
        
        !For free surface condition,
        !fold z_interp_coeff wrt inquiry point, which is
        !* freesurface point (at q(1)) in pressure case, or
        !* grid point that has global index=1 in velocity case (as freesurface is midway between vz(1) and vz(2))
        if(h%if_freesurface .and. h%ifz<1) then
            n=1-h%ifz   !compute the number of layers above the inquiry point
            iquiry=-r+n !compute the projected index of the inquiry point
            
            if(fold_type=='antisym') then
                !antisymmetric folding (pressure case): subtract points at and below inquiry point by points at and above inquiry point
                !this includes the inquiry point such that coeff=0. at the inquiry point
                z_interp_coeff(iquiry:iquiry+n) = z_interp_coeff(iquiry:iquiry+n) - z_interp_coeff(iquiry:-r:-1)
                
            elseif(fold_type=='symmetric') then
                !symmetric folding (velocity case): add points above and at and 1point below inquiry point to points below and at inquiry point
                !such that coeff(iquiry)=coeff(iquiry+1)
                z_interp_coeff(iquiry:iquiry+n+1) = z_interp_coeff(iquiry:iquiry+n+1) + z_interp_coeff(iquiry+1:-r:-1)
                
            endif
            
            !clean above inquiry point
            !!note that no antisym or symmetric condition required for vx & vz
            !!so just truncate the coeff above inquiry point (ie. freesurface as they are at same depth levels of q)
            z_interp_coeff(iquiry-1:-r:-1)=0.
            
        endif
        
        !tensor product to generate 3D tensor
        if(h%is_cubic) then
            call alloc(interp_coeff,[-r,r],[-r,r],[-r,r])
            do k=-r,r
            do j=-r,r
            do i=-r,r
                interp_coeff(i,j,k)=z_interp_coeff(i)*x_interp_coeff(j)*y_interp_coeff(k)
            enddo
            enddo
            enddo
        else
            call alloc(interp_coeff,[-r,r],[-r,r],[1,1])
            do j=-r,r
            do i=-r,r
                interp_coeff(i,j,1)=z_interp_coeff(i)*x_interp_coeff(j)
            enddo
            enddo
        endif
        
    end subroutine
    
    
    subroutine build_1d_hicks(x,dx,ifx,ix,ilx,x_interp_coeff)
        real,dimension(-r:r) :: x_interp_coeff, xx
        
        !---------------------------------------->
        !-4   -3   -2   -1   0   1   2   3   4     source local index (r=4)
        !                    *                     source point always on 0
        !                    |-|                   =alpha (-0.5<alpha<=0.5)
        !   +    +    +    +   +   +   +   +   +   interpolation grid points (2*r pts)
        !   9   10   11   12  13  14  15  16  17   grid global index
        !  -4   -3   -2   -1   0   1   2   3   4   grid projected index onto local indexation
        !  ifx                ix              ilx
        
        !find alpha such that
        ![-r:r]+alpha = global index of interpolation points
        tmp=x/dx
        ix = nint(tmp)
        alpha = ix - tmp
        if(alpha==-0.5) then
            ix=ix+1
            alpha=0.5
        endif
        
        ifx=ix-r  !must > me%ixmin
        ilx=ix+r  !must < me%ixmax
        
        if(abs(alpha)<epsilon) then !no need to interpolate
            x_interp_coeff(0)=1.
        else
            xx=[(i,i=-r,r)]+alpha
            x_interp_coeff=W(xx)*sinc(xx)
        endif
    end subroutine
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! elemental functions !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    elemental real function W(x)
        real, intent(in) :: x
        real y
         if (x<-r .or. x>r) then
             W=0.
         else
            y=x/r
            W=bessi0(b*sqrt(1.-y*y)) / bessi0(b)
         endif
    end function

    !Zero-order Modified Bessel function of the first kind
    !http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tbessi_f90.txt

    elemental real function bessi0(x)
        real, intent(in) :: x

        real ax
        double precision :: y
        
        double precision,parameter :: p1=1.0d0,p2=3.5156229d0,p3=3.0899424d0,p4=1.2067492d0,p5=0.2659732d0,p6=0.360768d-1,p7=0.45813d-2
        double precision,parameter :: q1=0.39894228d0,q2=0.1328592d-1,q3=0.225319d-2,q4=-0.157565d-2,q5=0.916281d-2,q6=-0.2057706d-1,q7=0.2635537d-1,q8=-0.1647633d-1,q9=0.392377d-2

        if (abs(x)<3.75) then
            y=x*x/14.0625  !=(x/3.75)^2
            bessi0=real( p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))) )
        else
            ax=abs(x)
            y=3.75/ax
            bessi0=real( (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9)))))))) )
        endif
    end function


    elemental real function sinc(x)
        real,intent(in) :: x
        
        real,parameter :: pi=3.141593
        
!         if (abs(x) < epsilon) then
!             sinc = 1.
!         else
            sinc = sin(pi*x) / (pi*x)
!         end if
    end function


end
