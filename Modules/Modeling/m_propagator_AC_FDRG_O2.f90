!code from Geodynamics or Komatitsch's GitHub site:
!https://github.com/geodynamics/seismic_cpml/blob/master/seismic_CPML_2D_pressure_second_order.f90


! 2D acoustic finite-difference code in pressure formulation
! with Convolutional-PML (C-PML) absorbing conditions for an heterogeneous isotropic acoustic medium

! Dimitri Komatitsch, CNRS, Marseille, July 2018.

! The pressure wave equation in an inviscid heterogeneous fluid is:
!
! 1/Kappa d2p / dt2 = div(grad(p) / rho) = d(1/rho dp/dx)/dx + d(1/rho dp/dy)/dy
!
! (see for instance Komatitsch and Tromp, Geophysical Journal International, vol. 149, p. 390-412 (2002), equations (19) and (21))
!
! The second-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used:
!
!            ^ y
!            |
!            |
!
!            +-------------------+
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!      dp/dy +---------+         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            +---------+---------+  ---> x
!            p       dp/dx
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

    !---
    !--- program starts here
    !---


!module m_propagator
program main

    ! private

    ! type,public :: t_propagator
    !   !info
    !   character(i_str_xxlen) :: info = &
    !       'Time-domain ISOtropic 2D ACoustic propagation'//s_NL// &
    !       '2nd-order Pressure formulation'//s_NL// &
    !       'Regular (non-staggered) Grid Finite-Difference (FDRG) method'//s_NL// &
    !       'Cartesian O(x²,t²) stencil'//s_NL// &
    !       ! 'CFL = Σ|coef| *Vmax *dt /rev_cell_diagonal'//s_NL// &
    !       ! '   -> dt ≤ 0.606(for 2D) or 0.494(3D) *Vmax/dx'//s_NL// &
    !       'Required model attributes: vp, rho'//s_NL// &
    !       'Required field components: p'//s_NL// &
    !       'Required boundary layer thickness: 2'//s_NL// &
    !       'Imaging conditions: P-Pxcorr'//s_NL// &
    !       'Energy terms: Σ_shot ∫ sfield%p² dt'//s_NL// &
    !       'Basic gradients: grho gkpa'

    !   integer :: nbndlayer=max(2,hicks_r) !minimum absorbing layer thickness
    !   integer :: ngrad=2 !number of basic gradients
    !   integer :: nimag=1 !number of basic images
    !   integer :: nengy=1 !number of energy terms

    !   logical :: if_compute_engy=.false.

    !   !local models shared between fields
    !   real,dimension(:,:,:),allocatable :: buoz, buox, buoy, kpa, inv_kpa

    !   !time frames
    !   integer :: nt
    !   real :: dt

    !   contains
    !   procedure :: print_info
    !   procedure :: estim_RAM
    !   procedure :: check_model
    !   procedure :: check_discretization
    !   procedure :: init
    !   procedure :: init_field
    !   procedure :: init_abslayer

    !   procedure :: forward
    !   procedure :: adjoint
      
    !   procedure :: inject_velocities
    !   procedure :: inject_stresses
    !   procedure :: update_velocities
    !   procedure :: update_stresses
    !   procedure :: extract


    !   final :: final

    ! end type


    ! flags to add PML layers to the edges of the grid
    logical, parameter :: USE_PML_XMIN = .true.
    logical, parameter :: USE_PML_XMAX = .true.
    logical, parameter :: USE_PML_YMIN = .true.
    logical, parameter :: USE_PML_YMAX = .true.

    ! total number of grid points in each direction of the grid
    integer, parameter :: NX = 201
    integer, parameter :: NY = 201

    ! size of a grid cell
    real, parameter :: DELTAX = 1.5
    real, parameter :: DELTAY = DELTAX

    ! thickness of the PML layer in grid points
    integer, parameter :: NPOINTS_PML = 10

    ! P-velocity and density
    ! the unrelaxed value is the value at frequency = 0 (the relaxed value would be the value at frequency = +infinity)
    real, parameter :: cp_unrelaxed = 2000
    real, parameter :: density = 2000

    ! total number of time steps
    integer, parameter :: NSTEP = 500

    ! time step in seconds
    real, parameter :: DELTAT = 5.2d-4

    ! parameters for the source
    real, parameter :: f0 = 35
    real, parameter :: t0 = 1.20d0 / f0
    real, parameter :: factor = 1

    ! source (in pressure)
    real, parameter :: xsource = 150
    real, parameter :: ysource = 150
    integer, parameter :: ISOURCE = xsource / DELTAX + 1
    integer, parameter :: JSOURCE = ysource / DELTAY + 1

    ! receivers
    integer, parameter :: NREC = 1
    !! DK DK I use 2301 here instead of 2300 in order to fall exactly on a grid point
    real, parameter :: xdeb = 7.5   ! first receiver x in meters
    real, parameter :: ydeb = 7.5   ! first receiver y in meters
    real, parameter :: xfin = 300   ! last receiver x in meters
    real, parameter :: yfin = 7.5   ! last receiver y in meters

    ! display information on the screen from time to time
    integer, parameter :: IT_DISPLAY = 100

    ! value of PI
    real, parameter :: PI = 3.1415927

    ! zero
    real, parameter :: ZERO = 0

    ! large value for maximum
    real, parameter :: HUGEVAL = huge(1.)

    ! threshold above which we consider that the code became unstable
    real, parameter :: STABILITY_THRESHOLD = huge(1.)

    ! main arrays
    real, dimension(NX,NY) :: pressure_past,pressure_present,pressure_future, &
      pressure_xx,pressure_yy,dpressurexx_dx,dpressureyy_dy,kappa_unrelaxed,rho,Kronecker_source

    ! to interpolate material parameters or velocity at the right location in the staggered grid cell
    real :: rho_half_x,rho_half_y

    ! power to compute d0 profile
    real, parameter :: NPOWER = 2

    ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
    real, parameter :: K_MAX_PML = 1
    real, parameter :: ALPHA_MAX_PML = 2*PI*(f0/2) ! from Festa and Vilotte

    ! arrays for the memory variables
    ! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
    real, dimension(NX,NY) :: &
      memory_dpressure_dx, &
      memory_dpressure_dy, &
      memory_dpressurexx_dx, &
      memory_dpressureyy_dy

    real :: &
      value_dpressure_dx, &
      value_dpressure_dy, &
      value_dpressurexx_dx, &
      value_dpressureyy_dy

    ! 1D arrays for the damping profiles
    real, dimension(NX) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half,K_x_half,alpha_x_half,a_x_half,b_x_half
    real, dimension(NY) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half,K_y_half,alpha_y_half,a_y_half,b_y_half

    real :: thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
    real :: Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized

    ! for the source
    real :: a,t,source_term

    ! for receivers
    real xspacerec,yspacerec,distval,dist
    integer, dimension(NREC) :: ix_rec,iy_rec
    real, dimension(NREC) :: xrec,yrec
    integer :: myNREC

    ! for seismograms
    real, dimension(NSTEP,NREC) :: sispressure

    integer :: i,j,it,irec

    real :: Courant_number,pressurenorm



    print *
    print *,'2D acoustic finite-difference code in pressure formulation with C-PML'
    print *

    ! display size of the model
    print *
    print *,'NX = ',NX
    print *,'NY = ',NY
    print *
    print *,'size of the model along X = ',(NX - 1) * DELTAX
    print *,'size of the model along Y = ',(NY - 1) * DELTAY
    print *
    print *,'Total number of grid points = ',NX * NY
    print *

    !--- define profile of absorption in PML region

    ! thickness of the PML layer in meters
    thickness_PML_x = NPOINTS_PML * DELTAX
    thickness_PML_y = NPOINTS_PML * DELTAY

    ! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    Rcoef = 0.001d0

    ! check that NPOWER is okay
    if (NPOWER < 1) stop 'NPOWER must be greater than 1'

    ! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
    d0_x = - (NPOWER + 1) * cp_unrelaxed * log(Rcoef) / (2 * thickness_PML_x)
    d0_y = - (NPOWER + 1) * cp_unrelaxed * log(Rcoef) / (2 * thickness_PML_y)

    print *,'d0_x = ',d0_x
    print *,'d0_y = ',d0_y
    print *

    d_x(:) = ZERO
    d_x_half(:) = ZERO
    K_x(:) = 1
    K_x_half(:) = 1
    alpha_x(:) = ZERO
    alpha_x_half(:) = ZERO
    a_x(:) = ZERO
    a_x_half(:) = ZERO

    d_y(:) = ZERO
    d_y_half(:) = ZERO
    K_y(:) = 1
    K_y_half(:) = 1
    alpha_y(:) = ZERO
    alpha_y_half(:) = ZERO
    a_y(:) = ZERO
    a_y_half(:) = ZERO

    ! damping in the X direction

    ! origin of the PML layer (position of right edge minus thickness, in meters)
    xoriginleft = thickness_PML_x
    xoriginright = (NX-1)*DELTAX - thickness_PML_x

    do i = 1,NX

    ! abscissa of current grid point along the damping profile
    xval = DELTAX * (i-1)

    !---------- left edge
    if (USE_PML_XMIN) then

    ! define damping profile at the grid points
      abscissa_in_PML = xoriginleft - xval
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
    ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x(i) = 1 + (K_MAX_PML - 1) * abscissa_normalized**NPOWER
        alpha_x(i) = ALPHA_MAX_PML * (1 - abscissa_normalized)
      endif

    ! define damping profile at half the grid points
      abscissa_in_PML = xoriginleft - (xval + DELTAX/2)
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
    ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_half(i) = 1 + (K_MAX_PML - 1) * abscissa_normalized**NPOWER
        alpha_x_half(i) = ALPHA_MAX_PML * (1 - abscissa_normalized)
      endif

    endif

    !---------- right edge
    if (USE_PML_XMAX) then

    ! define damping profile at the grid points
      abscissa_in_PML = xval - xoriginright
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
    ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x(i) = 1 + (K_MAX_PML - 1) * abscissa_normalized**NPOWER
        alpha_x(i) = ALPHA_MAX_PML * (1 - abscissa_normalized)
      endif

    ! define damping profile at half the grid points
      abscissa_in_PML = xval + DELTAX/2 - xoriginright
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
    ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_half(i) = 1 + (K_MAX_PML - 1) * abscissa_normalized**NPOWER
        alpha_x_half(i) = ALPHA_MAX_PML * (1 - abscissa_normalized)
      endif

    endif

    ! just in case, for -5 at the end
    if (alpha_x(i) < ZERO) alpha_x(i) = ZERO
    if (alpha_x_half(i) < ZERO) alpha_x_half(i) = ZERO

    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DELTAT)
    b_x_half(i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * DELTAT)

    ! this to avoid division by zero outside the PML
    if (abs(d_x(i)) > 1.e-6) a_x(i) = d_x(i) * (b_x(i) - 1) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)))
    if (abs(d_x_half(i)) > 1.e-6) a_x_half(i) = d_x_half(i) * &
      (b_x_half(i) - 1) / (K_x_half(i) * (d_x_half(i) + K_x_half(i) * alpha_x_half(i)))

    enddo

    ! damping in the Y direction

    ! origin of the PML layer (position of right edge minus thickness, in meters)
    yoriginbottom = thickness_PML_y
    yorigintop = (NY-1)*DELTAY - thickness_PML_y

    do j = 1,NY

    ! abscissa of current grid point along the damping profile
    yval = DELTAY * (j-1)

    !---------- bottom edge
    if (USE_PML_YMIN) then

    ! define damping profile at the grid points
      abscissa_in_PML = yoriginbottom - yval
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y(j) = d0_y * abscissa_normalized**NPOWER
    ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y(j) = 1 + (K_MAX_PML - 1) * abscissa_normalized**NPOWER
        alpha_y(j) = ALPHA_MAX_PML * (1 - abscissa_normalized)
      endif

    ! define damping profile at half the grid points
      abscissa_in_PML = yoriginbottom - (yval + DELTAY/2)
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half(j) = d0_y * abscissa_normalized**NPOWER
    ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y_half(j) = 1 + (K_MAX_PML - 1) * abscissa_normalized**NPOWER
        alpha_y_half(j) = ALPHA_MAX_PML * (1 - abscissa_normalized)
      endif

    endif

    !---------- top edge
    if (USE_PML_YMAX) then

    ! define damping profile at the grid points
      abscissa_in_PML = yval - yorigintop
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y(j) = d0_y * abscissa_normalized**NPOWER
    ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y(j) = 1 + (K_MAX_PML - 1) * abscissa_normalized**NPOWER
        alpha_y(j) = ALPHA_MAX_PML * (1 - abscissa_normalized)
      endif

    ! define damping profile at half the grid points
      abscissa_in_PML = yval + DELTAY/2 - yorigintop
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half(j) = d0_y * abscissa_normalized**NPOWER
    ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y_half(j) = 1 + (K_MAX_PML - 1) * abscissa_normalized**NPOWER
        alpha_y_half(j) = ALPHA_MAX_PML * (1 - abscissa_normalized)
      endif

    endif

    b_y(j) = exp(- (d_y(j) / K_y(j) + alpha_y(j)) * DELTAT)
    b_y_half(j) = exp(- (d_y_half(j) / K_y_half(j) + alpha_y_half(j)) * DELTAT)

    ! this to avoid division by zero outside the PML
    if (abs(d_y(j)) > 1.e-6) a_y(j) = d_y(j) * (b_y(j) - 1) / (K_y(j) * (d_y(j) + K_y(j) * alpha_y(j)))
    if (abs(d_y_half(j)) > 1.e-6) a_y_half(j) = d_y_half(j) * &
      (b_y_half(j) - 1) / (K_y_half(j) * (d_y_half(j) + K_y_half(j) * alpha_y_half(j)))

    enddo

    ! compute the Lame parameter and density
    do j = 1,NY
    do i = 1,NX
      rho(i,j) = density
      kappa_unrelaxed(i,j) = density*cp_unrelaxed*cp_unrelaxed
    enddo
    enddo

    ! print position of the source
    print *,'Position of the source:'
    print *
    print *,'x = ',xsource
    print *,'y = ',ysource
    print *

    ! define location of the source
    Kronecker_source(:,:) = 0
    Kronecker_source(ISOURCE,JSOURCE) = 1

    ! define location of receivers
    print *,'There are ',nrec,' receivers'
    print *
    if (NREC > 1) then
    ! this is to avoid a warning with GNU gfortran at compile time about division by zero when NREC = 1
    myNREC = NREC
    xspacerec = (xfin-xdeb) / (myNREC-1)
    yspacerec = (yfin-ydeb) / (myNREC-1)
    else
    xspacerec = 0
    yspacerec = 0
    endif
    do irec=1,nrec
    xrec(irec) = xdeb + (irec-1)*xspacerec
    yrec(irec) = ydeb + (irec-1)*yspacerec
    enddo

    ! find closest grid point for each receiver
    do irec=1,nrec
    dist = HUGEVAL
    do j = 1,NY
    do i = 1,NX
      distval = sqrt((DELTAX*(i-1) - xrec(irec))**2 + (DELTAY*(j-1) - yrec(irec))**2)
      if (distval < dist) then
        dist = distval
        ix_rec(irec) = i
        iy_rec(irec) = j
      endif
    enddo
    enddo
    print *,'receiver ',irec,' x_target,y_target = ',xrec(irec),yrec(irec)
    print *,'closest grid point found at distance ',dist,' in i,j = ',ix_rec(irec),iy_rec(irec)
    print *
    enddo

    ! check the Courant stability condition for the explicit time scheme
    ! R. Courant et K. O. Friedrichs et H. Lewy (1928)
    Courant_number = cp_unrelaxed * DELTAT * sqrt(1/DELTAX**2 + 1/DELTAY**2)
    print *,'Courant number is ',Courant_number
    print *
    if (Courant_number > 1) stop 'time step is too large, simulation will be unstable'

    ! initialize arrays
    pressure_present(:,:) = ZERO
    pressure_past(:,:) = ZERO

    ! PML
    memory_dpressure_dx(:,:) = ZERO
    memory_dpressure_dy(:,:) = ZERO
    memory_dpressurexx_dx(:,:) = ZERO
    memory_dpressureyy_dy(:,:) = ZERO

    ! initialize seismograms
    sispressure(:,:) = ZERO

    !---
    !---  beginning of time loop
    !---
    open(11,file='pressure_present',access='stream')

    do it = 1,NSTEP

        ! output information
        if (mod(it,IT_DISPLAY) == 0 .or. it == 5) then

            ! print maximum of pressure and of norm of velocity
            pressurenorm = maxval(abs(pressure_future))
            print *,'Time step # ',it,' out of ',NSTEP
            print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
            print *,'Max absolute value of pressure = ',pressurenorm
            print *
            ! check stability of the code, exit if unstable
            if (pressurenorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

        endif



        ! compute the first spatial derivatives divided by density

        do j = 1,NY
            do i = 1,NX-1
            value_dpressure_dx = (pressure_present(i+1,j) - pressure_present(i,j)) / DELTAX

            memory_dpressure_dx(i,j) = b_x_half(i) * memory_dpressure_dx(i,j) + a_x_half(i) * value_dpressure_dx

            value_dpressure_dx = value_dpressure_dx / K_x_half(i) + memory_dpressure_dx(i,j)

            rho_half_x = 0.5d0 * (rho(i+1,j) + rho(i,j))
            pressure_xx(i,j) = value_dpressure_dx / rho_half_x
            enddo
        enddo

        do j = 1,NY-1
            do i = 1,NX
            value_dpressure_dy = (pressure_present(i,j+1) - pressure_present(i,j)) / DELTAY

            memory_dpressure_dy(i,j) = b_y_half(j) * memory_dpressure_dy(i,j) + a_y_half(j) * value_dpressure_dy

            value_dpressure_dy = value_dpressure_dy / K_y_half(j) + memory_dpressure_dy(i,j)

            rho_half_y = 0.5d0 * (rho(i,j+1) + rho(i,j))
            pressure_yy(i,j) = value_dpressure_dy / rho_half_y
            enddo
        enddo

        ! compute the second spatial derivatives

        do j = 1,NY
            do i = 2,NX
            value_dpressurexx_dx = (pressure_xx(i,j) - pressure_xx(i-1,j)) / DELTAX

            memory_dpressurexx_dx(i,j) = b_x(i) * memory_dpressurexx_dx(i,j) + a_x(i) * value_dpressurexx_dx

            value_dpressurexx_dx = value_dpressurexx_dx / K_x(i) + memory_dpressurexx_dx(i,j)

            dpressurexx_dx(i,j) = value_dpressurexx_dx
            enddo
        enddo

        do j = 2,NY
            do i = 1,NX
            value_dpressureyy_dy = (pressure_yy(i,j) - pressure_yy(i,j-1)) / DELTAY

            memory_dpressureyy_dy(i,j) = b_y(j) * memory_dpressureyy_dy(i,j) + a_y(j) * value_dpressureyy_dy

            value_dpressureyy_dy = value_dpressureyy_dy / K_y(j) + memory_dpressureyy_dy(i,j)

            dpressureyy_dy(i,j) = value_dpressureyy_dy
            enddo
        enddo

        ! add the source (pressure located at a given grid point)
        a = pi*pi*f0*f0
        t = (it-1)*DELTAT

        ! Gaussian
        ! source_term = - factor * exp(-a*(t-t0)**2) / (2 * a)

        ! first derivative of a Gaussian
        ! source_term = factor * (t-t0)*exp(-a*(t-t0)**2)

        ! Ricker source time function (second derivative of a Gaussian)
        source_term = factor * (1 - 2*a*(t-t0)**2)*exp(-a*(t-t0)**2)

        ! apply the time evolution scheme
        ! we apply it everywhere, including at some points on the edges of the domain that have not be calculated above,
        ! which is of course wrong (or more precisely undefined), but this does not matter because these values
        ! will be erased by the Dirichlet conditions set on these edges below
        pressure_future(:,:) = - pressure_past(:,:) + 2 * pressure_present(:,:) + &
                                      DELTAT*DELTAT * ((dpressurexx_dx(:,:) + dpressureyy_dy(:,:)) * kappa_unrelaxed(:,:) + &
                                      4 * PI * cp_unrelaxed**2 * source_term * Kronecker_source(:,:))

        ! apply Dirichlet conditions at the bottom of the C-PML layers,
        ! which is the right condition to implement in order for C-PML to remain stable at long times

        ! Dirichlet condition for pressure on the left boundary
        pressure_future(1,:) = ZERO

        ! Dirichlet condition for pressure on the right boundary
        pressure_future(NX,:) = ZERO

        ! Dirichlet condition for pressure on the bottom boundary
        pressure_future(:,1) = ZERO

        ! Dirichlet condition for pressure on the top boundary
        pressure_future(:,NY) = ZERO

        ! store seismograms
        do irec = 1,NREC
            sispressure(it,irec) = pressure_future(ix_rec(irec),iy_rec(irec))
        enddo

        
        write(11) pressure_present

        ! move new values to old values (the present becomes the past, the future becomes the present)
        pressure_past = pressure_present
        pressure_present = pressure_future

    enddo   ! end of the time loop
    close(11)

    open(13,file='sispressure',access='stream')
    write(13) sispressure
    close(13)

    print *
    print *,'End of the simulation'
    print *

end