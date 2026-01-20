TYPE optimaltransport
                                                                
SEQUENCE

  !RESIDUAL
  INTEGER :: nt_resamp

  !SDMM
  INTEGER :: ichoice_dmm
  INTEGER :: niter_max
  REAL :: conv
  INTEGER :: ifpost_sdmm_weight

  !WASSERSTEIN DISTANCE
  REAL :: scale_residual,scal_WAS
  REAL :: bound_adjoint
  REAL :: const
  REAL :: taper_rfl

!  REAL,DIMENSION(:),POINTER :: D
!  REAL,DIMENSION(:,:),POINTER :: ATA,ATA_sparse,Z,Zt
!  REAL,DIMENSION(:,:),POINTER :: residuals_L2save,residuals_Wsave,data_cal_save
!  DOUBLE PRECISION,DIMENSION(:),POINTER :: dD
!  DOUBLE PRECISION,DIMENSION(:,:),POINTER :: dZ,dZt

!  !3D WASSERSTEIN DISTANCE
!  integer :: nrec_x,nrec_y
!  integer,dimension(:),POINTER :: A3D_row,A3D_row_CSR,A3D_col
!  double precision,dimension(:),POINTER :: A3D

END TYPE optimaltransport
