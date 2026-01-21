!Copyright 2013-2015 SEISCOPEII project, All rights reserved.

!==============================================================
subroutine read_toyxdac_time_input(acqui,pbdir)
  implicit none
  include 'acqui.h'
  include 'pbdirect.h'
  !acquisition
  type (acquisition) :: acqui
  !pb direct 
  type (pbdirect) :: pbdir

  open(10,file='toyxdac_time_input')
  read(10,*)pbdir%mode !0--5,11,14
  read(10,*)pbdir%tool !1--4
  read(10,*)acqui%type
  if(acqui%type==1) read(10,*)acqui%file_acqui  !read acquisition file if acqui%type=1. if acqui%type=2, acquisition is embeded in the su file.
  if(acqui%type==2)acqui%first_su=1
  close(10)
end subroutine read_toyxdac_time_input

!============================================================================
subroutine read_fdtd_input(pbdir,file_vp,file_rho,file_epsilon,file_delta,file_qp,file_vp_smth,file_rho_smth,file_source,file_topo)
  implicit none
  include 'pbdirect.h'
  include 'common.h'
 
  !pb direct 
  type (pbdirect) :: pbdir

  character(len=80) :: file_vp,file_rho,file_epsilon,file_delta,file_qp,file_source,file_topo
  CHARACTER(LEN=80) :: file_vp_smth,file_rho_smth


  open(11,file='fdtd_input')
  read(11,*)pbdir%n1,pbdir%n2,pbdir%n3
  read(11,*)pbdir%h
  if(pbdir%tool==1) read(11,*)file_vp,file_rho                                 !isotropic model
  if(pbdir%tool==2) read(11,*)file_vp,file_rho,file_epsilon,file_delta         !VTI model
  if(pbdir%tool==3) read(11,*)file_vp,file_rho,                        file_qp !viscous isotropic model
  if(pbdir%tool==4) read(11,*)file_vp,file_rho,file_epsilon,file_delta,file_qp !viscous VTI model
  IF(pbdir%mode==11.or.pbdir%mode==14) READ(11,*)file_vp_smth,file_rho_smth!the second model to implement JFWI
  read(11,*)pbdir%ibnd,pbdir%nlayer
  read(11,*)pbdir%oop
  read(11,*)file_source
  read(11,*)pbdir%ihicks
  read(11,*)pbdir%fm
  read(11,*)pbdir%nts,pbdir%dt,pbdir%nt
  read(11,*)pbdir%itypesource
  read(11,*)pbdir%ifreesurf
  read(11,*)file_topo
  read(11,*)pbdir%idata_p,pbdir%idata_vx,pbdir%idata_vy,pbdir%idata_vz
  if(pbdir%mode==5) read(11,*)pbdir%family,pbdir%kt
  close(11)

  if(pbdir%tool==3 .or. pbdir%tool==4) then
     open(12,file='attenuation_management')
     read(12,*)pbdir%nl !number of attenuation mechanisms
     read(12,*)pbdir%Q0 !reference Q0 for least-squares fitting of constant Q
     read(12,*)pbdir%fmin,pbdir%fmax !minimum and maximum reference frequencies
     if(pbdir%mode==1.or.pbdir%mode==2.or.pbdir%mode==4.or.pbdir%mode==6.or.pbdir%mode==11.or.pbdir%mode==14) read(12,*) pbdir%nc!number of checkpoints if inversion
     close(12)
  endif
end subroutine read_fdtd_input

!====================================================================
subroutine read_fwi_input(inv,ot,optim,name_bathy,depth_deadzone)
  implicit none
  include 'inversion.h'
  include 'optimaltransport.h'
  include 'optim_type.h'

  !inversion
  type(fullwaveinv) :: inv
  !optimaltransport
  type(optimaltransport) :: ot
  !optimization
  type (optim_type ) :: optim

  !data
  character(len=80) :: name_bathy
  real depth_deadzone
  integer ::i

  open(15,file='fwi_input')
  read(15,*)inv%idata_p,inv%idata_vx,inv%idata_vy,inv%idata_vz
  read(15,*)inv%file_data_p,inv%file_data_vx,inv%file_data_vy,inv%file_data_vz
  read(15,*)inv%invfamily
  read(15,*)inv%npar,(inv%invpar(i),i=1,inv%npar)
  read(15,*)inv%misfit
  READ(15,*)inv%iftime_domain
  read(15,*)inv%opt_meth, optim%l, optim%conv,optim%niter_max, inv%iprecond, inv%zpower
  read(15,*)inv%ibound,(inv%lb(i),inv%ub(i),i=1,inv%npar) 
  read(15,*)inv%lambda, inv%lambda_x, inv%lambda_y,inv%lambda_z
  read(15,*)name_bathy
  read(15,*)depth_deadzone
  READ(15,*)inv%dt_boundary
  READ(15,*)inv%bkwtime_stop
  READ(15,*)inv%nshot_out, (inv%ishot_out(i),i=1,inv%nshot_out)
  READ(15,*)
  READ(15,*)inv%data_weight_file
  READ(15,*)inv%nt_data_weight,inv%dt_data_weight,inv%nx_data_weight,inv%dx_data_weight
  READ(15,*)inv%mute_file
  READ(15,*)inv%data_weight_file2 !post sdmm weighting
  READ(15,*)inv%nt_data_weight2, inv%dt_data_weight2,inv%nx_data_weight2,inv%dx_data_weight2
  READ(15,*)inv%mute_file2
  close(15)
  
  IF (inv%misfit>0) CALL read_ot_management(ot)
    
end subroutine read_fwi_input

!====================================================================
subroutine read_ot_management(ot)
  implicit none
  include 'optimaltransport.h'
  !optimaltransport
  type(optimaltransport) :: ot
  open(16,file='ot_management')
  read(16,*)ot%ichoice_dmm
  read(16,*)ot%conv, ot%niter_max
  read(16,*)ot%scale_residual, ot%bound_adjoint
  read(16,*)ot%taper_rfl
  read(16,*)ot%ifpost_sdmm_weight
  close(16)
end subroutine read_ot_management

!====================================================================
subroutine read_ipwi_input(inv,optim,name_bathy,depth_deadzone)
  implicit none
  include 'inversion.h'
  include 'optim_type.h'

  !acquisition
  TYPE(fullwaveinv) :: inv
  !optimization
  type (optim_type ) :: optim

  !DATA
  CHARACTER(LEN=80) :: name_bathy
  REAL depth_deadzone
  integer :: i

  OPEN(15,file='ipwi_input')
  read(15,*)inv%idata_p,inv%idata_vx,inv%idata_vy,inv%idata_vz
  read(15,*)inv%file_data_p,inv%file_data_vx,inv%file_data_vy,inv%file_data_vz
  read(15,*)inv%invfamily
  read(15,*)inv%npar,(inv%invpar(i),i=1,inv%npar)
  read(15,*)inv%opt_meth, optim%l, optim%conv,optim%niter_max, inv%iprecond, inv%zpower
  read(15,*)inv%ibound,(inv%lb(i),inv%ub(i),i=1,inv%npar) 
  read(15,*)inv%lambda, inv%lambda_x, inv%lambda_y,inv%lambda_z
  read(15,*)name_bathy
  read(15,*)depth_deadzone
  READ(15,*)inv%dt_boundary
  READ(15,*)inv%bkwtime_stop
  READ(15,*)inv%nshot_out, (inv%ishot_out(i),i=1,inv%nshot_out)
  READ(15,*)
  READ(15,*)inv%data_weight_file !weight only on reflection data
  READ(15,*)inv%nt_data_weight,inv%dt_data_weight,inv%nx_data_weight,inv%dx_data_weight
  READ(15,*)inv%mute_file
  CLOSE(15)

end subroutine read_ipwi_input

!====================================================================
subroutine read_rwi_input(inv,ot,optim,name_bathy,depth_deadzone)
  implicit none
  include 'inversion.h'
  include 'optimaltransport.h'
  include 'optim_type.h'

  !acquisition
  TYPE(fullwaveinv) :: inv
  !optimaltransport
  type(optimaltransport) :: ot
  !optimization
  type (optim_type ) :: optim

  !DATA
  CHARACTER(LEN=80) :: name_bathy
  REAL depth_deadzone
  integer :: i

  OPEN(15,file='rwi_input')
  READ(15,*)inv%idata_p,inv%idata_vx,inv%idata_vy,inv%idata_vz
  READ(15,*)inv%file_data_p,inv%file_data_vx,inv%file_data_vy,inv%file_data_vz
  READ(15,*)inv%invfamily
  READ(15,*)inv%npar,(inv%invpar(i),i=1,inv%npar)
  READ(15,*)inv%misfit
  READ(15,*)inv%iftime_domain
  READ(15,*)inv%opt_meth, optim%l, optim%conv,optim%niter_max, inv%iprecond, inv%zpower
  READ(15,*)inv%ibound,(inv%lb(i),inv%ub(i),i=1,inv%npar)
  READ(15,*)inv%lambda, inv%lambda_x, inv%lambda_y,inv%lambda_z
  READ(15,*)name_bathy
  READ(15,*)depth_deadzone
  READ(15,*)inv%dt_boundary
  READ(15,*)inv%bkwtime_stop
  READ(15,*)inv%nshot_out, (inv%ishot_out(i),i=1,inv%nshot_out)
  READ(15,*)
  READ(15,*)inv%data_weight_file  !weight for reflection data
  READ(15,*)inv%nt_data_weight,  inv%dt_data_weight,  inv%nx_data_weight,  inv%dx_data_weight
  READ(15,*)inv%mute_file
!   READ(15,*)inv%data_weight_file2 !weight for diving data
!   READ(15,*)inv%nt_data_weight2, inv%dt_data_weight2, inv%nx_data_weight2, inv%dx_data_weight2
!   READ(15,*)inv%mute_file2
  CLOSE(15)
  
  IF (inv%misfit>0) CALL read_ot_management(ot)

end subroutine read_rwi_input

!====================================================================
subroutine read_datasep_management(inv,acqui)
  implicit none
  include 'acqui.h'
  include 'inversion.h'

  !acquisition
  TYPE (acquisition) :: acqui
  !acquisition
  TYPE(fullwaveinv) :: inv

  !DATA
  integer :: i

  OPEN(17,file='datasep_management')
  READ(17,*)inv%twind_boundary_file
  READ(17,*)inv%ratiod2r  !to be multiplied on reflected waves
  READ(17,*)acqui%offset_cut
  CLOSE(17)

end subroutine read_datasep_management

