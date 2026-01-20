TYPE fullwaveinv
                                                                
SEQUENCE

!INVERSION/OPTIMIZATION
   INTEGER :: misfit !misfit choice
   REAL :: fcost_div,fcost_refl,fcost_full,fcost,fcost_model
   INTEGER :: npar,npar0,invfamily,invpar(10),opt_meth, ibound,iprecond
   REAL :: lb(10),ub(10)
   REAL :: zpower
   REAL :: lambda, lambda_x,lambda_z,lambda_y 

   INTEGER,DIMENSION(:,:),POINTER     :: ibathy
   REAL,DIMENSION(:,:),POINTER     :: bathy
   INTEGER :: ideadzone !zone bellow bathy where gradient is not computed, but is allowed to change through regularization or preconditionning

!  integer,dimension(5)::iprint !print this parameter or not !  !Pengliang Yang, 06/2016:


   INTEGER :: firstgrad
   REAL :: scalingfactor
   INTEGER :: n1
   REAL :: h
   REAL,DIMENSION(:,:,:),allocatable :: gradient,gradient_tmp
   REAL,DIMENSION(:,:,:),allocatable :: grho,gkappa,gqp   !Pengliang Yang, 07/2016:basic parameter family for chain rule 
   REAL,DIMENSION(:,:,:),allocatable :: model, model_update, model_ip, model_ip_smth

!PSEUDO TIME INVERSION
   INTEGER :: iftime_domain !if go to pseudo time
   
!ENVELOPE OF DATA
   REAL,DIMENSION(:,:),POINTER :: data_env
   REAL,DIMENSION(:),POINTER :: data_mass
   
!PARAMETRIZATION
 REAL,DIMENSION(10) :: refvalue



!DATA
   INTEGER :: idata_p,idata_vx,idata_vy,idata_vz
   INTEGER :: nt_data_weight,  nx_data_weight
   INTEGER :: nt_data_weight2, nx_data_weight2
   REAL    :: dt_data_weight,  dx_data_weight
   REAL    :: dt_data_weight2, dx_data_weight2

   CHARACTER(LEN=80) :: file_data_p, file_data_vx, file_data_vy, file_data_vz
   CHARACTER(LEN=80) :: data_weight_file, data_weight_file2
   CHARACTER(LEN=80) :: mute_file, mute_file2

   REAL,DIMENSION(:,:),POINTER :: data_p ,data_vx ,data_vy ,data_vz
   REAL,DIMENSION(:,:),POINTER :: residual_p,residual_vx,residual_vy,residual_vz

   REAL,DIMENSION(:,:),POINTER :: residual_L2, residual_KR

   REAL    :: dt_boundary,bkwtime_stop

!SEPARATION

  CHARACTER(LEN=80)  :: twind_boundary_file

  REAL    ::  ratiod2r

!I/O
  INTEGER :: nshot_out, ishot_out(10)


END TYPE fullwaveinv
