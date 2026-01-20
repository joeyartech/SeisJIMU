! common variables
  !MPI
  INTEGER :: mype,nproc,infompi
  COMMON mype,nproc,infompi

  !OpenMP
  INTEGER :: nrang
  COMMON nrang

  !others
  REAL :: epsilon,pi
  PARAMETER (epsilon =0.001) 
  PARAMETER (pi = 3.14159265)

  ! FD weights
  REAL :: a1,a2
  PARAMETER(  a1=1.125,  a2=-1./24.)
  REAL :: b1,b2,b3,b4
  PARAMETER(b1=1225./1024,b2=-245./3072.,b3=49./5120.,b4=-5./7168.)

  !---C-PML
  REAL,PARAMETER :: K_x = 1.
  REAL,PARAMETER :: K_y = 1
  REAL,PARAMETER :: K_z = 1.
  REAL,PARAMETER :: K_x1 = 1.
  REAL,PARAMETER :: K_y1 = 1
  REAL,PARAMETER :: K_z1 = 1.
  REAL,PARAMETER ::NPOWER=2.
  REAL,PARAMETER :: Rcoef=0.001
