!CPML implementation from Komatitsch & Martin, 2007, Geophysics
!code from Geodynamics or Komatitsch's GitHub site:
!https://github.com/geodynamics/seismic_cpml
!!
! Dimitri Komatitsch, CNRS, Marseille, July 2018.
!
! The C-PML implementation is based in part on formulas given in Roden and Gedney (2000).
! If you use this code for your own research, please cite some (or all) of these
! articles:
!
! @ARTICLE{MaKoEz08,
! author = {Roland Martin and Dimitri Komatitsch and Abdela\^aziz Ezziani},
! title = {An unsplit convolutional perfectly matched layer improved at grazing
! incidence for seismic wave equation in poroelastic media},
! journal = {Geophysics},
! year = {2008},
! volume = {73},
! pages = {T51-T61},
! number = {4},
! doi = {10.1190/1.2939484}}
!
! @ARTICLE{MaKo09,
! author = {Roland Martin and Dimitri Komatitsch},
! title = {An unsplit convolutional perfectly matched layer technique improved
! at grazing incidence for the viscoelastic wave equation},
! journal = {Geophysical Journal International},
! year = {2009},
! volume = {179},
! pages = {333-344},
! number = {1},
! doi = {10.1111/j.1365-246X.2009.04278.x}}
!
! @ARTICLE{MaKoGe08,
! author = {Roland Martin and Dimitri Komatitsch and Stephen D. Gedney},
! title = {A variational formulation of a stabilized unsplit convolutional perfectly
! matched layer for the isotropic or anisotropic seismic wave equation},
! journal = {Computer Modeling in Engineering and Sciences},
! year = {2008},
! volume = {37},
! pages = {274-304},
! number = {3}}
!
! @ARTICLE{KoMa07,
! author = {Dimitri Komatitsch and Roland Martin},
! title = {An unsplit convolutional {P}erfectly {M}atched {L}ayer improved
!          at grazing incidence for the seismic wave equation},
! journal = {Geophysics},
! year = {2007},
! volume = {72},
! number = {5},
! pages = {SM155-SM167},
! doi = {10.1190/1.2757586}}
!
! The original CPML technique for Maxwell's equations is described in:
!
! @ARTICLE{RoGe00,
! author = {J. A. Roden and S. D. Gedney},
! title = {Convolution {PML} ({CPML}): {A}n Efficient {FDTD} Implementation
!          of the {CFS}-{PML} for Arbitrary Media},
! journal = {Microwave and Optical Technology Letters},
! year = {2000},
! volume = {27},
! number = {5},
! pages = {334-339},
! doi = {10.1002/1098-2760(20001205)27:5 < 334::AID-MOP14>3.0.CO;2-A}}

module m_cpml
use m_System
use m_math, only: r_pi
use m_resampler
use m_model
use m_shot
use m_computebox

    private

    !CPML parameter
    real,parameter :: npower = 2. !must be >1
    real,parameter :: logR = log(0.001)
    real,parameter :: kpa_max=1.
    !Comments by Komatitsch:
    !from stephen gedney's unpublished class notes for class ee699, lecture 8, slide 8-11
    !and increase kpa_max (eg. =7) if absorbing is not satisfactory at grazing incident angle
    !(making in/outside PML reflection more separate..)
    !decrease this number if grid dispersion is not satisfactory
    real,parameter :: kpa_max_m1=kpa_max-1.

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
    
    subroutine init(self)
        class(t_cpml) :: self
        
        real,dimension(:),allocatable :: alfa_z,alfa_z_half
        real,dimension(:),allocatable :: alfa_x,alfa_x_half
        real,dimension(:),allocatable :: alfa_y,alfa_y_half
        real,dimension(:),allocatable :: d_z,d_z_half
        real,dimension(:),allocatable :: d_x,d_x_half
        real,dimension(:),allocatable :: d_y,d_y_half       

        !self%nlayer=setup%get_int('CPML_SIZE','NCPML',o_default='20')+add_thickness
        self%nlayer=cb%nabslayer

        alfa_max = r_pi*shot%fpeak  !Komatitsch: from festa and vilotte

        thickness_pml_z = self%nlayer * m%dz
        thickness_pml_x = self%nlayer * m%dx
        thickness_pml_y = self%nlayer * m%dy
        
        !Komatitsch: compute d0 from inria report section 6.1 http://hal.inria.fr/docs/00/07/32/19/pdf/rr-3471.pdf
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
                self%kpa_z(i)   = 1. + kpa_max_m1*abscissa_normalized**npower ! from stephen gedney's unpublished class notes for class ee699, lecture 8, slide 8-2
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
