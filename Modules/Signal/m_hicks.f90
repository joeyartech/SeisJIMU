!Hicks interpolation
!To interpolate source & receiver off-grid positions
!
!Hicks, 2002, Arbitrary source and receiver positioning in finite‐difference schemes using Kaiser windowed sinc function, _Geophysics_
!https://doi.org/10.1190/1.1451454
!
!Author comment:
!We made interactive illustration of this paper's theory (e.g. Fig 9) for the free surface case:
! - hicks_freesurf_symmetric.ggb
! - hicks_freesurf_antisymmetric.ggb
!with the free GeoGebra tool.
!
module m_hicks
use m_arrayop
use m_math

    private
    public :: hicks_init, hicks_put_position, hicks_get_coef, hicks_get_position

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

    !intent(in)
    real :: z,x,y  !point location
    real :: dz,dx,dy
    logical :: is_cubic,is_freesurface

    !intent(out)
    !range of interpolation points
    integer :: ifz,ifx,ify
    integer :: iz,ix,iy
    integer :: ilz,ilx,ily

    contains

    subroutine hicks_init(dz_,dx_,dy_,is_cubic_,is_freesurface_)
        logical :: is_cubic_,is_freesurface_
        dz=dz_; dx=dx_; dy=dy_
        is_cubic=is_cubic_; is_freesurface=is_freesurface_
    end subroutine

    subroutine hicks_put_position(z_,x_,y_)
        z=z_; x=x_; y=y_
    end subroutine

    subroutine hicks_get_coef(fold_type,interp_coef)
        character(*),intent(in) :: fold_type
        real,dimension(:,:,:),allocatable,intent(out) :: interp_coef
        
        real,dimension(-r:r) :: z_interp_coef,x_interp_coef,y_interp_coef

        z_interp_coef=0.
        x_interp_coef=0.
        y_interp_coef=0.
        
        !z axis
        call build_1d_coef(z,dz,ifz,iz,ilz,z_interp_coef)
        !x axis
        call build_1d_coef(x,dx,ifx,ix,ilx,x_interp_coef)
        !y axis
        call build_1d_coef(y,dy,ify,iy,ily,y_interp_coef)
        
        !For free surface condition,
        !fold z_interp_coeff wrt inquiry point, which is
        !* freesurface point (at q(1)) in pressure case, or
        !* grid point that has global index=1 in velocity case (as freesurface is midway between vz(1) and vz(2))
        if(is_freesurface .and. ifz<1) then
            n=1-ifz   !compute the number of layers above the inquiry point
            iquiry=-r+n !compute the projected index of the inquiry point
            
            if(fold_type=='antisymm') then
                !antisymmetric folding (pressure case): subtract points at and below inquiry point by points at and above inquiry point
                !this includes the inquiry point such that coeff=0. at the inquiry point
                z_interp_coef(iquiry:iquiry+n) = z_interp_coef(iquiry:iquiry+n) - z_interp_coef(iquiry:-r:-1)
                
            elseif(fold_type=='symmetric') then
                !symmetric folding (velocity case): add points above and at and 1point below inquiry point to points below and at inquiry point
                !such that coeff(iquiry)=coeff(iquiry+1)
                z_interp_coef(iquiry:iquiry+n+1) = z_interp_coef(iquiry:iquiry+n+1) + z_interp_coef(iquiry+1:-r:-1)
                
            endif
            
            !clean above inquiry point
            !!note that no antisym or symmetric condition required for vx & vz
            !!so just truncate the coeff above inquiry point (ie. freesurface as they are at same depth levels of q)
            z_interp_coef(iquiry-1:-r:-1)=0.
            
        endif
        
        !tensor product to generate 3D tensor
        if(is_cubic) then
            call alloc(interp_coef,[-r,r],[-r,r],[-r,r])
            do k=-r,r
            do j=-r,r
            do i=-r,r
                interp_coef(i,j,k)=z_interp_coef(i)*x_interp_coef(j)*y_interp_coef(k)
            enddo
            enddo
            enddo
        else
            call alloc(interp_coef,[-r,r],[-r,r],[1,1])
            do j=-r,r
            do i=-r,r
                interp_coef(i,j,1)=z_interp_coef(i)*x_interp_coef(j)
            enddo
            enddo
        endif
        
    end subroutine

    subroutine hicks_get_position(ifz_,ifx_,ify_,iz_,ix_,iy_,ilz_,ilx_,ily_)
        !index starts from 1 outside scope of build_1d
        ifz_=ifz+1; iz_=iz+1; ilz_=ilz+1
        ifx_=ifx+1; ix_=ix+1; ilx_=ilx+1
        if(is_cubic) then
            ify_=ify+1; iy_=iy+1; ily_=ily+1
        else
            ify_=1;     iy_=1;    ily_=1
        endif

    end subroutine

    subroutine build_1d_coef(x,dx,ifx,ix,ilx,x_interp_coef)
        real,dimension(-r:r) :: x_interp_coef, xx
        
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
        
        if(abs(alpha)<r_epsilon) then !no need to interpolate
            x_interp_coef(0)=1.
        else
            xx=[(i,i=-r,r)]+alpha
            x_interp_coef=W(xx)*sinc(xx)
        endif

    end subroutine

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

end
