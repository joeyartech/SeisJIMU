module m_cpml
use m_System
use m_math, only: r_pi
use m_resampler
use m_model
use m_shot
use m_computebox

    private

    !CPML parameter
    real,parameter :: npower = 2.
    real,parameter :: logR = log(0.001)
    real :: kpa_max
    !empirically kpa_max=1 for acoustic modeling,
    !and kpa_max=7 for elastic (P-SV) modeling
    !Comments by Komatitsch:
    !increase kpa_max if absorbing is not satisfactory at grazing incident angle
    !(making in/outside PML reflection more separate..)
    !decrease this number if grid dispersion is not satisfactory
    real :: kpa_max_m1

    type,public :: t_cpml
        real,dimension(:),allocatable :: b_z,b_z_half, a_z,a_z_half
        real,dimension(:),allocatable :: b_x,b_x_half, a_x,a_x_half
        real,dimension(:),allocatable :: b_y,b_y_half, a_y,a_y_half

        real,dimension(:),allocatable :: kpa_z,kpa_z_half
        real,dimension(:),allocatable :: kpa_x,kpa_x_half
        real,dimension(:),allocatable :: kpa_y,kpa_y_half

        integer :: nlayer

        contains
        procedure :: init

    end type
        
    type(t_cpml),public :: cpml

    contains
    
    !CPML implementation from Komatitsch & Martin, 2007, Geophysics
    !code from Geodynamics or Komatitsch's GitHub site:
    !https://github.com/geodynamics/seismic_cpml
    subroutine init(self,o_kpa_max)
        class(t_cpml) :: self
        real,optional :: o_kpa_max
        
        real,dimension(:),allocatable :: alfa_z,alfa_z_half
        real,dimension(:),allocatable :: alfa_x,alfa_x_half
        real,dimension(:),allocatable :: alfa_y,alfa_y_half
        real,dimension(:),allocatable :: d_z,d_z_half
        real,dimension(:),allocatable :: d_x,d_x_half
        real,dimension(:),allocatable :: d_y,d_y_half       

        kpa_max=either(o_kpa_max,7.,present(o_kpa_max))
        kpa_max_m1=kpa_max-1

        !self%nlayer=setup%get_int('CPML_SIZE','NCPML',o_default='20')+add_thickness
        self%nlayer=cb%nabslayer

        alfa_max = r_pi*shot%fpeak

        thickness_pml_z = self%nlayer * m%dz
        thickness_pml_x = self%nlayer * m%dx
        thickness_pml_y = self%nlayer * m%dy
        
        d0_z = -(npower+1)/(2.*thickness_pml_z) * cb%velmax * logR
        d0_x = -(npower+1)/(2.*thickness_pml_x) * cb%velmax * logR
        d0_y = -(npower+1)/(2.*thickness_pml_y) * cb%velmax * logR
        
        call alloc(alfa_z,[cb%ifz,cb%ilz]); call alloc(alfa_z_half,[cb%ifz,cb%ilz])
        call alloc(alfa_x,[cb%ifx,cb%ilx]); call alloc(alfa_x_half,[cb%ifx,cb%ilx])
        call alloc(alfa_y,[cb%ify,cb%ily]); call alloc(alfa_y_half,[cb%ify,cb%ily])

        call alloc(d_z,[cb%ifz,cb%ilz]); call alloc(d_z_half,[cb%ifz,cb%ilz])       
        call alloc(d_x,[cb%ifx,cb%ilx]); call alloc(d_x_half,[cb%ifx,cb%ilx])
        call alloc(d_y,[cb%ify,cb%ily]); call alloc(d_y_half,[cb%ify,cb%ily])

        call alloc(self%a_z,[cb%ifz,cb%ilz]); call alloc(self%a_z_half,[cb%ifz,cb%ilz])
        call alloc(self%a_x,[cb%ifx,cb%ilx]); call alloc(self%a_x_half,[cb%ifx,cb%ilx])
        call alloc(self%a_y,[cb%ify,cb%ily]); call alloc(self%a_y_half,[cb%ify,cb%ily])
        call alloc(self%b_z,[cb%ifz,cb%ilz]); call alloc(self%b_z_half,[cb%ifz,cb%ilz])
        call alloc(self%b_x,[cb%ifx,cb%ilx]); call alloc(self%b_x_half,[cb%ifx,cb%ilx])
        call alloc(self%b_y,[cb%ify,cb%ily]); call alloc(self%b_y_half,[cb%ify,cb%ily])
                
        call alloc(self%kpa_z,[cb%ifz,cb%ilz],o_init=1.); call alloc(self%kpa_z_half,[cb%ifz,cb%ilz],o_init=1.)
        call alloc(self%kpa_x,[cb%ifx,cb%ilx],o_init=1.); call alloc(self%kpa_x_half,[cb%ifx,cb%ilx],o_init=1.)
        call alloc(self%kpa_y,[cb%ify,cb%ily],o_init=1.); call alloc(self%kpa_y_half,[cb%ify,cb%ily],o_init=1.)

        !z dir
        do i = cb%ifz,cb%ilz

            !!top edge
            abscissa_in_pml = -(i-1)*m%dz
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z(i)          = d0_z * abscissa_normalized**npower
                self%kpa_z(i)   = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_z(i)        = alfa_max * (1. - abscissa_normalized)
            endif

            !!top edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dz + m%dz/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z_half(i)     = d0_z * abscissa_normalized**npower
                self%kpa_z_half(i)   = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_z_half(i)   = alfa_max * (1. - abscissa_normalized)
            endif

            !!bottom edge
            abscissa_in_pml = (i-cb%mz)*m%dz
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z(i)          = d0_z * abscissa_normalized**npower
                self%kpa_z(i)   = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_z(i)        = alfa_max * (1.- abscissa_normalized)
            endif

            !!bottom edge half gridpoint
            abscissa_in_pml = (i-cb%mz)*m%dz - m%dz/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z_half(i)        = d0_z * abscissa_normalized**npower
                self%kpa_z_half(i) = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_z_half(i)      = alfa_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alfa_z(i)     <0.) alfa_z(i)=0.
            if(alfa_z_half(i)<0.) alfa_z_half(i)=0.

            self%b_z(i)      = exp(-(d_z(i)     /self%kpa_z(i)     +alfa_z(i)     )*shot%dt)
            self%b_z_half(i) = exp(-(d_z_half(i)/self%kpa_z_half(i)+alfa_z_half(i))*shot%dt)

            if(abs(d_z(i)     )>1e-6) self%a_z(i)     =d_z(i)     *(self%b_z(i)     -1.)/(self%kpa_z(i)     *(d_z(i)     +self%kpa_z(i)     *alfa_z(i)     ))
            if(abs(d_z_half(i))>1e-6) self%a_z_half(i)=d_z_half(i)*(self%b_z_half(i)-1.)/(self%kpa_z_half(i)*(d_z_half(i)+self%kpa_z_half(i)*alfa_z_half(i)))

        enddo

        self%kpa_z     =1./self%kpa_z
        self%kpa_z_half=1./self%kpa_z_half

        !x dir
        do i = cb%ifx,cb%ilx

            !!left edge
            abscissa_in_pml = -(i-1)*m%dx
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x(i)          = d0_x * abscissa_normalized**npower
                self%kpa_x(i)   = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_x(i)      = alfa_max * (1. - abscissa_normalized)
            endif

            !!left edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dx + m%dx/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x_half(i)     = d0_x * abscissa_normalized**npower
                self%kpa_x_half(i) = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_x_half(i) = alfa_max * (1. - abscissa_normalized)
            endif

            !!right edge
            abscissa_in_pml = (i-cb%mx)*m%dx
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x(i)          = d0_x * abscissa_normalized**npower
                self%kpa_x(i)   = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_x(i)      = alfa_max * (1.- abscissa_normalized)
            endif

            !!right edge half gridpoint
            abscissa_in_pml = (i-cb%mx)*m%dx - m%dx/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x_half(i)     = d0_x * abscissa_normalized**npower
                self%kpa_x_half(i) = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_x_half(i) = alfa_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alfa_x(i)     <0.) alfa_x(i)=0.
            if(alfa_x_half(i)<0.) alfa_x_half(i)=0.

            self%b_x(i)      = exp(-(d_x(i)     /self%kpa_x(i)     +alfa_x(i)     )*shot%dt)
            self%b_x_half(i) = exp(-(d_x_half(i)/self%kpa_x_half(i)+alfa_x_half(i))*shot%dt)

            if(abs(d_x(i)     )>1e-6) self%a_x(i)     =d_x(i)     *(self%b_x(i)     -1.)/(self%kpa_x(i)     *(d_x(i)     +self%kpa_x(i)     *alfa_x(i)     ))
            if(abs(d_x_half(i))>1e-6) self%a_x_half(i)=d_x_half(i)*(self%b_x_half(i)-1.)/(self%kpa_x_half(i)*(d_x_half(i)+self%kpa_x_half(i)*alfa_x_half(i)))

        enddo

        self%kpa_x     =1./self%kpa_x
        self%kpa_x_half=1./self%kpa_x_half


        !y dir
        do i = cb%ify,cb%ily

            !!front edge
            abscissa_in_pml = -(i-1)*m%dy
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y(i)          = d0_y * abscissa_normalized**npower
                self%kpa_y(i)   = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_y(i)      = alfa_max * (1. - abscissa_normalized)
            endif

            !!front edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dy + m%dy/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y_half(i)     = d0_y * abscissa_normalized**npower
                self%kpa_y_half(i)   = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_y_half(i) = alfa_max * (1. - abscissa_normalized)
            endif

            !!rear edge
            abscissa_in_pml = (i-cb%my)*m%dy
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y(i)          = d0_y * abscissa_normalized**npower
                self%kpa_y(i)   = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_y(i)      = alfa_max * (1.- abscissa_normalized)
            endif

            !!rear edge half gridpoint
            abscissa_in_pml = (i-cb%my)*m%dy - m%dy/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y_half(i)     = d0_y * abscissa_normalized**npower
                self%kpa_y_half(i)   = 1. + kpa_max_m1*abscissa_normalized**npower
                alfa_y_half(i) = alfa_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alfa_y(i)     <0.) alfa_y(i)=0.
            if(alfa_y_half(i)<0.) alfa_y_half(i)=0.

            self%b_y(i)      = exp(-(d_y(i)     /self%kpa_y(i)     +alfa_y(i)     )*shot%dt)
            self%b_y_half(i) = exp(-(d_y_half(i)/self%kpa_y_half(i)+alfa_y_half(i))*shot%dt)

            if(abs(d_y(i)     )>1e-6) self%a_y(i)     =d_y(i)     *(self%b_y(i)     -1.)/(self%kpa_y(i)     *(d_y(i)     +self%kpa_y(i)     *alfa_y(i)     ))
            if(abs(d_y_half(i))>1e-6) self%a_y_half(i)=d_y_half(i)*(self%b_y_half(i)-1.)/(self%kpa_y_half(i)*(d_y_half(i)+self%kpa_y_half(i)*alfa_y_half(i)))

        enddo

        self%kpa_y     =1./self%kpa_y
        self%kpa_y_half=1./self%kpa_y_half
        
        !no more need these variables
        deallocate(alfa_z,alfa_z_half)
        deallocate(alfa_x,alfa_x_half)
        deallocate(alfa_y,alfa_y_half)
        deallocate(d_z,d_z_half)
        deallocate(d_x,d_x_half)
        deallocate(d_y,d_y_half)
        
    end subroutine
    
end
