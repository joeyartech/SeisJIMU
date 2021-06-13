module m_cpml
use m_global
use m_string
use m_mpienv
use m_arrayop
use m_setup
use m_resampler
use m_model
use m_shot
use m_computebox

    private

    !CPML parameter
    real,parameter :: k_z = 1.
    real,parameter :: k_x = 1.
    real,parameter :: k_y = 1.
    real,parameter :: npower = 2.
    real,parameter :: logRcoeff = log(0.001)
    real,parameter :: kpa_max=7.  !increase this number if absorbing is not satisfactory at grazing incident angle
                                    !(make in/outside PML reflection more separate..)
                                    !decrease this number if grid dispersion is not satisfactory
    real,parameter :: kpa_max_=kpa_max-1.

    type,public :: t_cpml
        real,dimension(:),allocatable :: b_z,b_z_half, a_z,a_z_half
        real,dimension(:),allocatable :: b_x,b_x_half, a_x,a_x_half
        real,dimension(:),allocatable :: b_y,b_y_half, a_y,a_y_half

        real,dimension(:),allocatable :: kpa_z,kpa_z_half
        real,dimension(:),allocatable :: kpa_x,kpa_x_half
        real,dimension(:),allocatable :: kpa_y,kpa_y_half

        integer :: nlayer

        contains
        procedure :: init => init

    end type
        
    type(t_cpml),public :: cpml

    contains
    
    !CPML implementation from Komatitsch & Martin, 2007, Geophysics
    !code from Geodynamics or Komatitsch's GitHub site:
    !https://github.com/geodynamics/seismic_cpml
    subroutine init(self,addwidth)
        class(t_cpml) :: self
        integer :: addwidth !=2 for FDTD-O(dx4) and Hicks interpolation of radius=2
                            !=4 for FDTD-O(dx8)

        real,dimension(:),allocatable :: alp_z,alp_z_half
        real,dimension(:),allocatable :: alp_x,alp_x_half
        real,dimension(:),allocatable :: alp_y,alp_y_half
        real,dimension(:),allocatable :: d_z,d_z_half
        real,dimension(:),allocatable :: d_x,d_x_half
        real,dimension(:),allocatable :: d_y,d_y_half       

        self%nlayer=setup%get_int('CPML_SIZE','NCPML',o_default='20')+addwidth  

        alp_max = r_pi*shot%fpeak

        thickness_pml_z = self%nlayer * m%dz
        thickness_pml_x = self%nlayer * m%dx
        thickness_pml_y = self%nlayer * m%dy
        
        d0_z = -(npower+1)/(2.*thickness_pml_z) * cb%velmax * logRcoeff
        d0_x = -(npower+1)/(2.*thickness_pml_x) * cb%velmax * logRcoeff
        d0_y = -(npower+1)/(2.*thickness_pml_y) * cb%velmax * logRcoeff
        
        call alloc(alp_z,[cb%ifz,cb%ilz]); call alloc(alp_z_half,[cb%ifz,cb%ilz])
        call alloc(alp_x,[cb%ifx,cb%ilx]); call alloc(alp_x_half,[cb%ifx,cb%ilx])
        call alloc(alp_y,[cb%ify,cb%ily]); call alloc(alp_y_half,[cb%ify,cb%ily])

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
                self%kpa_z(i)   = 1. + kpa_max_*abscissa_normalized**npower
                alp_z(i)        = alp_max * (1. - abscissa_normalized)
            endif

            !!top edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dz + m%dz/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z_half(i)     = d0_z * abscissa_normalized**npower
                self%kpa_z_half(i)   = 1. + kpa_max_*abscissa_normalized**npower
                alp_z_half(i)   = alp_max * (1. - abscissa_normalized)
            endif

            !!bottom edge
            abscissa_in_pml = (i-cb%mz)*m%dz
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z(i)          = d0_z * abscissa_normalized**npower
                self%kpa_z(i)   = 1. + kpa_max_*abscissa_normalized**npower
                alp_z(i)        = alp_max * (1.- abscissa_normalized)
            endif

            !!bottom edge half gridpoint
            abscissa_in_pml = (i-cb%mz)*m%dz - m%dz/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z_half(i)        = d0_z * abscissa_normalized**npower
                self%kpa_z_half(i) = 1. + kpa_max_*abscissa_normalized**npower
                alp_z_half(i)      = alp_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alp_z(i)     <0.) alp_z(i)=0.
            if(alp_z_half(i)<0.) alp_z_half(i)=0.

            self%b_z(i)      = exp(-(d_z(i)     /self%kpa_z(i)     +alp_z(i)     )*shot%dt)
            self%b_z_half(i) = exp(-(d_z_half(i)/self%kpa_z_half(i)+alp_z_half(i))*shot%dt)

            if(abs(d_z(i)     )>1e-6) self%a_z(i)     =d_z(i)     *(self%b_z(i)     -1.)/(self%kpa_z(i)     *(d_z(i)     +self%kpa_z(i)     *alp_z(i)     ))
            if(abs(d_z_half(i))>1e-6) self%a_z_half(i)=d_z_half(i)*(self%b_z_half(i)-1.)/(self%kpa_z_half(i)*(d_z_half(i)+self%kpa_z_half(i)*alp_z_half(i)))

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
                self%kpa_x(i)   = 1. + kpa_max_*abscissa_normalized**npower
                alp_x(i)      = alp_max * (1. - abscissa_normalized)
            endif

            !!left edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dx + m%dx/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x_half(i)     = d0_x * abscissa_normalized**npower
                self%kpa_x_half(i) = 1. + kpa_max_*abscissa_normalized**npower
                alp_x_half(i) = alp_max * (1. - abscissa_normalized)
            endif

            !!right edge
            abscissa_in_pml = (i-cb%mx)*m%dx
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x(i)          = d0_x * abscissa_normalized**npower
                self%kpa_x(i)   = 1. + kpa_max_*abscissa_normalized**npower
                alp_x(i)      = alp_max * (1.- abscissa_normalized)
            endif

            !!right edge half gridpoint
            abscissa_in_pml = (i-cb%mx)*m%dx - m%dx/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x_half(i)     = d0_x * abscissa_normalized**npower
                self%kpa_x_half(i) = 1. + kpa_max_*abscissa_normalized**npower
                alp_x_half(i) = alp_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alp_x(i)     <0.) alp_x(i)=0.
            if(alp_x_half(i)<0.) alp_x_half(i)=0.

            self%b_x(i)      = exp(-(d_x(i)     /self%kpa_x(i)     +alp_x(i)     )*shot%dt)
            self%b_x_half(i) = exp(-(d_x_half(i)/self%kpa_x_half(i)+alp_x_half(i))*shot%dt)

            if(abs(d_x(i)     )>1e-6) self%a_x(i)     =d_x(i)     *(self%b_x(i)     -1.)/(self%kpa_x(i)     *(d_x(i)     +self%kpa_x(i)     *alp_x(i)     ))
            if(abs(d_x_half(i))>1e-6) self%a_x_half(i)=d_x_half(i)*(self%b_x_half(i)-1.)/(self%kpa_x_half(i)*(d_x_half(i)+self%kpa_x_half(i)*alp_x_half(i)))

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
                self%kpa_y(i)   = 1. + kpa_max_*abscissa_normalized**npower
                alp_y(i)      = alp_max * (1. - abscissa_normalized)
            endif

            !!front edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dy + m%dy/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y_half(i)     = d0_y * abscissa_normalized**npower
                self%kpa_y_half(i)   = 1. + kpa_max_*abscissa_normalized**npower
                alp_y_half(i) = alp_max * (1. - abscissa_normalized)
            endif

            !!rear edge
            abscissa_in_pml = (i-cb%my)*m%dy
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y(i)          = d0_y * abscissa_normalized**npower
                self%kpa_y(i)   = 1. + kpa_max_*abscissa_normalized**npower
                alp_y(i)      = alp_max * (1.- abscissa_normalized)
            endif

            !!rear edge half gridpoint
            abscissa_in_pml = (i-cb%my)*m%dy - m%dy/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y_half(i)     = d0_y * abscissa_normalized**npower
                self%kpa_y_half(i)   = 1. + kpa_max_*abscissa_normalized**npower
                alp_y_half(i) = alp_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alp_y(i)     <0.) alp_y(i)=0.
            if(alp_y_half(i)<0.) alp_y_half(i)=0.

            self%b_y(i)      = exp(-(d_y(i)     /self%kpa_y(i)     +alp_y(i)     )*shot%dt)
            self%b_y_half(i) = exp(-(d_y_half(i)/self%kpa_y_half(i)+alp_y_half(i))*shot%dt)

            if(abs(d_y(i)     )>1e-6) self%a_y(i)     =d_y(i)     *(self%b_y(i)     -1.)/(self%kpa_y(i)     *(d_y(i)     +self%kpa_y(i)     *alp_y(i)     ))
            if(abs(d_y_half(i))>1e-6) self%a_y_half(i)=d_y_half(i)*(self%b_y_half(i)-1.)/(self%kpa_y_half(i)*(d_y_half(i)+self%kpa_y_half(i)*alp_y_half(i)))

        enddo

        self%kpa_y     =1./self%kpa_y
        self%kpa_y_half=1./self%kpa_y_half
        
        !no more need these variables
        deallocate(alp_z,alp_z_half)
        deallocate(alp_x,alp_x_half)
        deallocate(alp_y,alp_y_half)
        deallocate(d_z,d_z_half)
        deallocate(d_x,d_x_half)
        deallocate(d_y,d_y_half)
        
    end subroutine
    
end
