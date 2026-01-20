!Copyright 2013-2015 SEISCOPEII project, All rights reserved.

SUBROUTINE fcost_source_adjoint_ot(pbdir,acqui,inv,ot,FLAG)

  implicit none
  INCLUDE 'mpif.h'
  include 'pbdirect.h'
  include 'acqui.h'
  include 'inversion.h'
  include 'optimaltransport.h'
  include 'common.h'


  !PBDIRECT
  TYPE(pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE(fullwaveinv) :: inv
  !OPTIMAL TRANSPORT
  TYPE(optimaltransport) :: ot
  
  CHARACTER(4) :: FLAG

  INTEGER :: irec,it,i, itwind
  REAL :: dist
  integer :: ntaper
  real,dimension(:),allocatable :: taper
  REAL :: tmp
  REAL :: fcost_tmp
  REAL :: const
  REAL :: fcost_comm(3),tmp_comm(3)

  character(4) :: str_mype
  
  double precision :: dot_prod(2), dot_prod_tmp(2), angle

#ifdef TIME_PROFILING
  REAL(kind=8) :: time_fcost,tmp_time,time_comm
#endif

#ifdef TIME_PROFILING
  tmp_time=MPI_WTIME()
#endif

  fcost_tmp=0.

  const = 0.5* pbdir%dt

  
  !COMPUTE L2 RESIDUAL
  
	inv%residual_p(it,irec)=(inv%data_p(it,irec)-pbdir%data_p(it,irec))*acqui%weight_data_time(it,irec)
	!remultiply by the weight for the source of the adjoint !
	inv%residual_p(it,irec)=inv%residual_p(it,irec)*acqui%weight_data_time(it,irec)
	  
  !WRITE L2 RESIDUAL & FCOST
  inv%residual_L2=inv%residual_p
  CALL write_residu(inv,pbdir%nt,acqui%nrec,inv%residual_L2,'residu_L2_')
  CALL write_fcost(inv,1,[fcost_tmp],'fcost_L2_')
  
  !ZW check mass conservation condition
  if( sum(inv%residual_p) > 1e-2*maxval(inv%residual_p) ) then
    write(*,*)'************************************************************ '
    write(*,*)'*** WARNING : MASS CONSERVATION MAY NOT SATISFACTED !!! **** '
    write(*,*)'*** SHOT #',mype,' HAS NON-NEGLIGIBLE SUM OF DATA RESIDUALS: ', sum(inv%residual_p)/maxval(inv%residual_p)
    write(*,*)'************************************************************ '
  end if
  
  !COMPUTE KR RESIDUAL
  inv%residual_p(:,:)=inv%residual_p(:,:)*ot%scale_residual
  SELECT CASE (ot%ichoice_dmm)
    CASE (7);  CALL ot_sdmm_fishpack(ot,inv,pbdir%nt,acqui%nrec)    !ZW FFT solver for Poisson eq, recommended
  CASE (8);  CALL ot_sdmm_mudpack(ot,inv,pbdir%nt,acqui%nrec)     !ZW Multigraid solver for Poisson eq.
  CASE (11); CALL ot_sdmm_mudpack_3D(ot,inv,pbdir%nt,acqui%nrec)  !ZW Multigrid solver for Poisson eq
  END SELECT
  
  inv%residual_KR=inv%residual_p
    

  !ZW weight after SDMM
  IF(ot%ifpost_sdmm_weight==1) THEN
    inv%residual_p=inv%residual_p*acqui%weight_data_time
    ntaper=nint(ot%taper_rfl/pbdir%dt)
    allocate(taper(ntaper))
    do i=1,ntaper
	taper(i)=(cos(  (i-1)*3.14/ntaper + 3.14  )+1)/2.
    enddo

    DO irec=1,acqui%nrec
    dist=SQRT((acqui%s1(1)-acqui%r1(irec))**2+(acqui%s2(1)-acqui%r2(irec))**2+(acqui%s3(1)-acqui%r3(irec))**2)
    i=1  !ntaper
    IF (dist>=acqui%offset_cut) THEN
	itwind= int(acqui%twind_boundary(irec)/pbdir%dt)+1

	DO it=1, itwind-ntaper-1, 1 !before reflection & transition
	  inv%residual_p(it,irec)=0.
	END DO

	DO it=itwind-ntaper , itwind-1, 1 !before reflection,transition, need tapering to maintain DC=0 for OT
	  if (it<1 .or. it>pbdir%nt) then
	    i=i+1
	    cycle
	  end if
	  inv%residual_p(it,irec)=inv%residual_p(it,irec)*taper(i)
	  i=i+1
	END DO
    ELSE
      inv%residual_p(:,irec)=0.
    END IF
    END DO

    deallocate(taper)
  END IF
  !END ZW
  
  CALL MPI_BARRIER(MPI_COMM_WORLD, INFOMPI)

!     !Dot product of L2 residual and unscaled Wasserstein residual
!     dot_prod(1)=sum(dble(inv%residual_L2)*dble(inv%residual_KR))
!     CALL MPI_REDUCE (dot_prod(1), dot_prod_tmp(1), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,  infompi )
!     if(mype==0) print*, 'Dot product of L2 and unscaled W residuals',dot_prod_tmp(1)
! 
!     !Angle between L2 residual and unscaled Wasserstein residual
!     dot_prod(1)=sum(dble( inv%residual_L2 * inv%residual_L2 ))
!     dot_prod(2)=sum(dble( inv%residual_KR * inv%residual_KR ))
!     CALL MPI_ALLREDUCE (dot_prod, dot_prod_tmp, 2, MPI_DOUBLE_PRECISION,  MPI_SUM, MPI_COMM_WORLD,  infompi )
! 
!     dot_prod(1)=sum(dble( inv%residual_L2/dsqrt(dot_prod_tmp(1)) * inv%residual_KR/dsqrt(dot_prod_tmp(2)) ))
! 
!     CALL MPI_REDUCE (dot_prod(1), angle, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,  infompi )
! 
!     angle=acos(angle)*180/3.1415927
!     if(mype==0) print*,'angle between L2 and unscaled KR residuals',angle
    
  
  
  !COMMUNICATE FCOST
  fcost_comm(1)=inv%fcost
  
  !WRITE W RESIDUAL & FCOST   
  CALL write_residu(inv,pbdir%nt,acqui%nrec,inv%residual_KR,'residu_KR_')
  CALL write_fcost(inv,1,[fcost_comm(1)],'fcost_KR_')
  

  !SUM ON ALL PROCESSORS
  tmp_comm=0.
  CALL MPI_ALLREDUCE (  fcost_comm(1),  tmp_comm,   1,   MPI_REAL,  MPI_SUM,  MPI_COMM_WORLD,  infompi )
  
  inv%fcost=tmp_comm(1)


END SUBROUTINE FCOST_SOURCE_ADJOINT_OT
