module m_absorbinglayer
use m_sysio
use m_arrayop
use m_model
use m_shot

    private

!geometry of model & computebox
!                                      ifx<---- nx --->ilx
!                                    ifz+~~~~~~~~~~~~~~~+
!       iox     iox+mx-1              ^ l   1      mx   l
!  ++====+-------+===========++       | l  1+-------+   l
!  ||    | *     |           ||       | l   |       |   l
!  || A  |   C   |     B     ||         l D |   C   |   l
!  ||    |       |           ||      nz l   |       |   l
!  ||    |       |           ||         l   |       |   l
!  ++====+-------+===========++       | l mz+-------+   l
!                                     v l boundarylayer l
!       Model: A+B+C                 ilz+~~~~~~~~~~~~~~~+
!       * source point
!                                        Computebox: C+D

    !CPML parameter
    real,parameter :: k_x = 1.
    real,parameter :: k_y = 1.
    real,parameter :: k_z = 1.
    real,parameter :: npower = 2.
    real,parameter :: Rcoeff = 0.001  
    real,parameter :: kappa_max=7.!increase this number if absorbing is not satisfactory at grazing incident angle
                                    !(make in/outside PML reflection more separate..)
                                    !decrease this number if grid dispersion is not satisfactory

    type,public :: t_absorbinglayer
        real,dimension(:),allocatable :: b_x,b_x_half,a_x,a_x_half
        real,dimension(:),allocatable :: b_y,b_y_half,a_y,a_y_half
        real,dimension(:),allocatable :: b_z,b_z_half,a_z,a_z_half

        real,dimension(:),allocatable :: kappa_x,kappa_x_half
        real,dimension(:),allocatable :: kappa_y,kappa_y_half
        real,dimension(:),allocatable :: kappa_z,kappa_z_half

        integer :: nlayer

        contains
        procedure :: init => init
    end type
        
    type(t_absorbinglayer),public :: abl

    contains
    
    !CPML implementation from Komatitsch & Martin, 2007, Geophysics
    !code from Geodynamics or Komatitsch's GitHub site:
    !https://github.com/geodynamics/seismic_cpml
    subroutine init(self)
        class(t_absorbinglayer) :: abl

        self%nlayer=setup%get_int('CPML_SIZE','NCPML',o_default='20')+2  !for FDTDo4 and hicks interpolation of radius=2
        !self%nlayer=setup_get_int('CPML_SIZE','NCPML')+4  !for FDTDo8
        
        real,dimension(:),allocatable   :: alpha_x,alpha_x_half
        real,dimension(:),allocatable   :: alpha_y,alpha_y_half
        real,dimension(:),allocatable   :: alpha_z,alpha_z_half
        real,dimension(:),allocatable   :: d_x,d_x_half
        real,dimension(:),allocatable   :: d_y,d_y_half
        real,dimension(:),allocatable   :: d_z,d_z_half
                
        alpha_max = r_pi*shot%src%fpeak

        thickness_pml_x = self%nlayer * m%dx
        thickness_pml_y = self%nlayer * m%dy
        thickness_pml_z = self%nlayer * m%dz
        
        d0_x = -(npower+1)/(2.*thickness_pml_x) * cb%velmax * log(Rcoeff)
        d0_y = -(npower+1)/(2.*thickness_pml_y) * cb%velmax * log(Rcoeff) 
        d0_z = -(npower+1)/(2.*thickness_pml_z) * cb%velmax * log(Rcoeff)

        call alloc(alpha_x,[cb%ifx,cb%ilx]); call alloc(alpha_x_half,[cb%ifx,cb%ilx])
        call alloc(alpha_y,[cb%ify,cb%ily]); call alloc(alpha_y_half,[cb%ify,cb%ily])
        call alloc(alpha_z,[cb%ifz,cb%ilz]); call alloc(alpha_z_half,[cb%ifz,cb%ilz])

        call alloc(d_x,[cb%ifx,cb%ilx]); call alloc(d_x_half,[cb%ifx,cb%ilx])
        call alloc(d_y,[cb%ify,cb%ily]); call alloc(d_y_half,[cb%ify,cb%ily])
        call alloc(d_z,[cb%ifz,cb%ilz]); call alloc(d_z_half,[cb%ifz,cb%ilz])

        call alloc(self%b_x,[cb%ifx,cb%ilx]); call alloc(self%b_x_half,[cb%ifx,cb%ilx])
        call alloc(self%a_x,[cb%ifx,cb%ilx]); call alloc(self%a_x_half,[cb%ifx,cb%ilx])
        call alloc(self%b_y,[cb%ify,cb%ily]); call alloc(self%b_y_half,[cb%ify,cb%ily])
        call alloc(self%a_y,[cb%ify,cb%ily]); call alloc(self%a_y_half,[cb%ify,cb%ily])
        call alloc(self%b_z,[cb%ifz,cb%ilz]); call alloc(self%b_z_half,[cb%ifz,cb%ilz])
        call alloc(self%a_z,[cb%ifz,cb%ilz]); call alloc(self%a_z_half,[cb%ifz,cb%ilz])

        call alloc(self%kappa_x,[cb%ifx,cb%ilx],initialize=.false.); call alloc(self%kappa_x_half,[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(self%kappa_y,[cb%ify,cb%ily],initialize=.false.); call alloc(self%kappa_y_half,[cb%ify,cb%ily],initialize=.false.)
        call alloc(self%kappa_z,[cb%ifz,cb%ilz],initialize=.false.); call alloc(self%kappa_z_half,[cb%ifz,cb%ilz],initialize=.false.)

        self%kappa_x=1.; self%kappa_x_half=1.
        self%kappa_y=1.; self%kappa_y_half=1.
        self%kappa_z=1.; self%kappa_z_half=1.


        !x dir
        do i = cb%ifx,cb%ilx

            !!left edge
            abscissa_in_pml = -(i-1)*m%dx
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x(i)          = d0_x * abscissa_normalized**npower
                self%kappa_x(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_x(i)      = alpha_max * (1. - abscissa_normalized)
            endif

            !!left edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dx + m%dx/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x_half(i)     = d0_x * abscissa_normalized**npower
                self%kappa_x_half(i) = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_x_half(i) = alpha_max * (1. - abscissa_normalized)
            endif

            !!right edge
            abscissa_in_pml = (i-cb%mx)*m%dx
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x(i)          = d0_x * abscissa_normalized**npower
                self%kappa_x(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_x(i)      = alpha_max * (1.- abscissa_normalized)
            endif

            !!right edge half gridpoint
            abscissa_in_pml = (i-cb%mx)*m%dx - m%dx/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x_half(i)     = d0_x * abscissa_normalized**npower
                self%kappa_x_half(i) = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_x_half(i) = alpha_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alpha_x(i)     <0.) alpha_x(i)=0.
            if(alpha_x_half(i)<0.) alpha_x_half(i)=0.

            self%b_x(i)      = exp(-(d_x(i)     /self%kappa_x(i)     +alpha_x(i)     )*shot%src%dt)
            self%b_x_half(i) = exp(-(d_x_half(i)/self%kappa_x_half(i)+alpha_x_half(i))*shot%src%dt)

            if(abs(d_x(i)     )>1e-6) self%a_x(i)     =d_x(i)     *(self%b_x(i)     -1.)/(self%kappa_x(i)     *(d_x(i)     +self%kappa_x(i)     *alpha_x(i)     ))
            if(abs(d_x_half(i))>1e-6) self%a_x_half(i)=d_x_half(i)*(self%b_x_half(i)-1.)/(self%kappa_x_half(i)*(d_x_half(i)+self%kappa_x_half(i)*alpha_x_half(i)))

        enddo

        self%kappa_x     =1./self%kappa_x
        self%kappa_x_half=1./self%kappa_x_half


        !y dir
        do i = cb%ify,cb%ily

            !!front edge
            abscissa_in_pml = -(i-1)*m%dy
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y(i)          = d0_y * abscissa_normalized**npower
                self%kappa_y(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_y(i)      = alpha_max * (1. - abscissa_normalized)
            endif

            !!front edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dy + m%dy/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y_half(i)     = d0_y * abscissa_normalized**npower
                self%kappa_y_half(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_y_half(i) = alpha_max * (1. - abscissa_normalized)
            endif

            !!rear edge
            abscissa_in_pml = (i-cb%my)*m%dy
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y(i)          = d0_y * abscissa_normalized**npower
                self%kappa_y(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_y(i)      = alpha_max * (1.- abscissa_normalized)
            endif

            !!rear edge half gridpoint
            abscissa_in_pml = (i-cb%my)*m%dy - m%dy/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y_half(i)     = d0_y * abscissa_normalized**npower
                self%kappa_y_half(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_y_half(i) = alpha_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alpha_y(i)     <0.) alpha_y(i)=0.
            if(alpha_y_half(i)<0.) alpha_y_half(i)=0.

            self%b_y(i)      = exp(-(d_y(i)     /self%kappa_y(i)     +alpha_y(i)     )*shot%src%dt)
            self%b_y_half(i) = exp(-(d_y_half(i)/self%kappa_y_half(i)+alpha_y_half(i))*shot%src%dt)

            if(abs(d_y(i)     )>1e-6) self%a_y(i)     =d_y(i)     *(self%b_y(i)     -1.)/(self%kappa_y(i)     *(d_y(i)     +self%kappa_y(i)     *alpha_y(i)     ))
            if(abs(d_y_half(i))>1e-6) self%a_y_half(i)=d_y_half(i)*(self%b_y_half(i)-1.)/(self%kappa_y_half(i)*(d_y_half(i)+self%kappa_y_half(i)*alpha_y_half(i)))

        enddo

        self%kappa_y     =1./self%kappa_y
        self%kappa_y_half=1./self%kappa_y_half

        !z dir
        do i = cb%ifz,cb%ilz

            !!top edge
            abscissa_in_pml = -(i-1)*m%dz
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z(i)          = d0_z * abscissa_normalized**npower
                self%kappa_z(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_z(i)      = alpha_max * (1. - abscissa_normalized)
            endif

            !!top edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dz + m%dz/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z_half(i)     = d0_z * abscissa_normalized**npower
                self%kappa_z_half(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_z_half(i) = alpha_max * (1. - abscissa_normalized)
            endif

            !!bottom edge
            abscissa_in_pml = (i-cb%mz)*m%dz
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z(i)          = d0_z * abscissa_normalized**npower
                self%kappa_z(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_z(i)      = alpha_max * (1.- abscissa_normalized)
            endif

            !!bottom edge half gridpoint
            abscissa_in_pml = (i-cb%mz)*m%dz - m%dz/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z_half(i)     = d0_z * abscissa_normalized**npower
                self%kappa_z_half(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_z_half(i) = alpha_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alpha_z(i)     <0.) alpha_z(i)=0.
            if(alpha_z_half(i)<0.) alpha_z_half(i)=0.

            self%b_z(i)      = exp(-(d_z(i)     /self%kappa_z(i)     +alpha_z(i)     )*shot%src%dt)
            self%b_z_half(i) = exp(-(d_z_half(i)/self%kappa_z_half(i)+alpha_z_half(i))*shot%src%dt)

            if(abs(d_z(i)     )>1e-6) self%a_z(i)     =d_z(i)     *(self%b_z(i)     -1.)/(self%kappa_z(i)     *(d_z(i)     +self%kappa_z(i)     *alpha_z(i)     ))
            if(abs(d_z_half(i))>1e-6) self%a_z_half(i)=d_z_half(i)*(self%b_z_half(i)-1.)/(self%kappa_z_half(i)*(d_z_half(i)+self%kappa_z_half(i)*alpha_z_half(i)))

        enddo

        self%kappa_z     =1./self%kappa_z
        self%kappa_z_half=1./self%kappa_z_half

        
        !no more need these variables
        deallocate(alpha_x,alpha_x_half)
        deallocate(alpha_y,alpha_y_half)
        deallocate(alpha_z,alpha_z_half)
        deallocate(d_x,d_x_half)
        deallocate(d_y,d_y_half)
        deallocate(d_z,d_z_half)
        
    end subroutine
    
end
