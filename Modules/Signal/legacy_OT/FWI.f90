!Copyright 2013-2015 SEISCOPEII project, All rights reserved.
!============================================================
subroutine FWI(pbdir,acqui)

  implicit none
  INCLUDE 'mpif.h'
  include 'pbdirect.h'
  include 'acqui.h'
  include 'common.h'
  include 'inversion.h'
  include 'optimaltransport.h'
  include 'optim_type.h'


  !PBDIRECT
  TYPE(pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE(fullwaveinv) :: inv
  !OPTIMAL TRANSPORT
  TYPE(optimaltransport) :: ot
  !OPTIMIZATION
  TYPE(optim_type) :: optim

  
  !SPECIFIC TO OPTIMIZATION TOOL BOX
  real f
  real,dimension(:),allocatable:: x, g,precond_g
  integer :: n
  character (LEN=4) :: FLAG

  !INTERN
  INTEGER :: i,j,iipar,iiipar,cpt_iter,ifound
  !DATA
  CHARACTER(LEN=80) :: name_bathy
  REAL depth_deadzone


  !tuning of optimization
  optim%debug =.false.
  !print out optimization info only on mype=0
  if(mype==0) then  
     optim%print_flag=1
  else
     optim%print_flag=0
  end if


  !=================================================================
  !step 1 -> read input
  !=================================================================
  call  read_fwi_input(inv,ot,optim,name_bathy,depth_deadzone)

  ot%nt_resamp=pbdir%nt !to be considered in the future
  
  !READ DATA
  call read_sismo(pbdir,inv,acqui)

!   !UPDATE THE HICKS ACQUISITION TABLES FOR INVERSION 
!   !attention, the following routine assumes only one type of data can be 
!   !used : velocity or pressure
!   if(pbdir%ihicks==1) call hicks_constant_inv(pbdir,acqui,inv)

  !COMPUTE WEIGHT AND SEPARATION ON DATA
  CALL compute_weight_data(inv,acqui,pbdir)
  

  pbdir%idata_p=inv%idata_p
  pbdir%idata_vx=inv%idata_vx
  pbdir%idata_vy=inv%idata_vy
  pbdir%idata_vz=inv%idata_vz


  !read bathy
  allocate(inv%bathy(pbdir%n2,pbdir%n3))
  allocate(inv%ibathy(pbdir%n2,pbdir%n3))
  open(17,file=name_bathy,access='direct',recl=4*pbdir%n2*pbdir%n3)
  read(17,rec=1)inv%bathy(:,:)
  close(17)
  !compute ibathy
  do j=1,pbdir%n3
     do i=1,pbdir%n2
        inv%ibathy(i,j)=nint(inv%bathy(i,j)/pbdir%h)+1
     end do
  end do
  inv%ideadzone=nint(depth_deadzone/pbdir%h)

  !compute the it_boundary for memory on boundary
  pbdir%it_boundary=nint(inv%dt_boundary/pbdir%dt)
  pbdir%nt_boundary=pbdir%nt/pbdir%it_boundary
  if(pbdir%ibnd/=3) then
     if(mype==0) then
        WRITE(*,*)'****************************************'
        WRITE(*,*)' pbdir%it_boundary',pbdir%it_boundary
        WRITE(*,*)' pbdir%nt_boundary',pbdir%nt_boundary
        WRITE(*,*)'FREQUENCY UP TO',0.5/((pbdir%it_boundary)*pbdir%dt),'HZ'
        WRITE(*,*)'       CAN BE CONSIDERED  !!!! '
        WRITE(*,*)'****************************************'
        IF(0.5/((pbdir%it_boundary)*pbdir%dt) < pbdir%fm) THEN
           WRITE(*,*)'****************************************'
           WRITE(*,*)' DECIMATION OF DATA ON BOUNDARY IS NOT '
           WRITE(*,*)'COMPATIBLE WITH FREQUENCY OF SIMULATION'
           WRITE(*,*)'----------------------------------------'
           WRITE(*,*)'           CODE STOPPED NOW !!!'
           WRITE(*,*)'****************************************'
           CALL MPI_ABORT()
        end if
     end if
  end if

  !====================================================================
  ! STEP 2 -> define the model to be inverted
  !=====================================================================
  IF (inv%iftime_domain==0) THEN !INVERSION IN DEPTH DOMAIN
     inv%n1=pbdir%n1; inv%h=pbdir%h
  ELSE !INVERSION IN PSEUDO TIME DOMAIN
     call init_d2t(inv%lb, inv%ub, pbdir%n1, pbdir%h, pbdir%n2, inv%n1, inv%h)
     if(mype==0)write(*,*)  'model size in time domain:', inv%n1, inv%h
  END IF
    
  allocate(inv%model(inv%n1,pbdir%n2,pbdir%n3))
  allocate(inv%model_update(inv%n1,pbdir%n2,pbdir%n3))
  
  IF (inv%iftime_domain==0) THEN !INVERSION IN DEPTH DOMAIN
     CALL sub_setmodelpbdir2inv(pbdir,inv)
  ELSE !INVERSION IN PSEUDO TIME DOMAIN
     call depth2time(pbdir%vp, pbdir%n1, pbdir%h, inv%n1, inv%h, pbdir%n2, inv%model)
  END IF

  
  !=================================================================
  !STEP 3 -> FWI OPTIMIZATION
  !=================================================================
   
  
  if(mype==0) write(*,*)'============= START OPTIMIZATION ==================='
  FLAG='INIT'
  inv%firstgrad=0
  
  !compute the size of the inverse problem
  n=inv%n1*pbdir%n2*pbdir%n3
  !allocate tables for gradient and model
  allocate(x(n), g(n))
  IF (inv%opt_meth == 1.OR.inv%opt_meth == 2.OR.inv%opt_meth == 4) allocate(precond_g(n))

  !manage the bounds
  IF(inv%ibound==1) THEN
     allocate(optim%lb(n),optim%ub(n))
        optim%lb=inv%lb(1)
        optim%ub=inv%ub(1)
        optim%threshold=1.
        optim%bound=1
  ELSE
     optim%bound=0
  END IF


  !FIRST STEP IS TO COMPUTE COST FUNCTION AND GRADIENT
  CALL modeling_FWI(pbdir,acqui,inv,ot)
  
  IF(pbdir%mode==4) GOTO 100  !exit of the code if mode==4 -> single gradient computation for benchmarking

  CALL sub_setmodel2x(n,g,pbdir,inv,inv%model_update)
  CALL sub_setmodel2x(n,x,pbdir,inv,inv%model)


  !1=PSD,2=PCG,3=LBFGS,4=PLBFGS
  if (inv%opt_meth == 1.or.inv%opt_meth == 2.or.inv%opt_meth == 4)then
     precond_g(:)=g(:)
     if(mype==0) write(*,*) '------------APPLY PRECONDITIONER----------------'
     if(inv%iprecond==1)  call precond_optimization(inv,pbdir,precond_g,n)
  end if

  cpt_iter=0
  DO WHILE (FLAG .ne. 'CONV' .and. FLAG .ne. 'FAIL')
  if(mype==0)write(70,*)cpt_iter,FLAG,inv%fcost, inv%fcost/inv%scalingfactor
     IF (inv%opt_meth == 1) call PSTD(n,x,inv%fcost,g,precond_g,optim,FLAG)
     IF (inv%opt_meth == 2) call PNLCG(n,x,inv%fcost,g,precond_g,optim,FLAG)
     IF (inv%opt_meth == 3) call LBFGS(n,x,inv%fcost,g,optim,FLAG)
     IF (inv%opt_meth == 4) call PLBFGS(n,x,inv%fcost,g,precond_g,optim,FLAG)

     !! fluch buffer !!
     CALL flushc()

     IF(FLAG .eq. 'GRAD') THEN
        if(mype==0) write(*,*)'---------COMPUTE GRADIENT--------'
        !compute fcost and gradient for model x
        CALL sub_setx2model(n,x,pbdir,inv,inv%model)
        CALL sub_modelinv2pbdir(pbdir,inv)
        CALL modeling_FWI(pbdir,acqui,inv,ot)
        CALL sub_setmodel2x(n,g,pbdir,inv,inv%model_update)
        CALL sub_setmodel2x(n,x,pbdir,inv,inv%model)

        IF (inv%opt_meth == 1.OR.inv%opt_meth == 2)THEN!precondition the descent direction in case of P-XXXX
           precond_g(:)=g(:)
           IF(mype==0) WRITE(*,*) 'APPLY PRECONDITIONER'
           IF(inv%iprecond==1)  CALL precond_optimization(inv,pbdir,precond_g,n)
        END IF

     elseif(FLAG.eq.'NSTE') then
! 	if(cpt_iter==0.and.mype==0) then
! 	  allocate(inv%grad_base(inv%n1,pbdir%n2,pbdir%n3), inv%grad_monitor(inv%n1,pbdir%n2,pbdir%n3))
! 	  inv%grad_base = inv%gradient
! 	  amod_base=sqrt(sum(inv%grad_base*inv%grad_base))
! 	end if
	
        cpt_iter=cpt_iter+1
        CALL sub_setx2model(n,x,pbdir,inv,inv%model)
        
        IF(mype==0) THEN !write only on master
            CALL sub_modelinv2pbdirforwrinting(pbdir,inv,cpt_iter,0)

! 	    ! monitor the gradient
! 	    !!amplitude
! 	    inv%grad_monitor = inv%gradient
! 	    amod_monitor=sqrt(sum(inv%grad_monitor*inv%grad_monitor))
! 	    print*,'Relative amplitude of gradient:',real(amod_monitor/amod_base)
! 	    !!phase
! 	    dot_prod=sum(inv%grad_base*inv%grad_monitor)
! 	    angle= dot_prod / amod_base / amod_monitor
! 	    angle= acos(angle)*180/3.1415927
! 	    print*, 'Relative angle of gradient:', angle
! 		      
! ! 	    write(10,*) cpt_iter, angle, real(amod_monitor/amod_base)
        END IF
        
	IF(mype==0) WRITE(*,*) 'NEW ITERATE',cpt_iter
	
     elseif(FLAG.eq.'PREC') then !apply preconditionner
        IF(mype==0) WRITE(*,*) 'APPLY PRECONDITIONER'
        IF(inv%iprecond==1)   CALL precond_optimization(inv,pbdir,optim%q_plb,n)
     END IF
     
     
  END DO
  !FINISH, WE WRITE THE FINAL MODEL
  CALL sub_setx2model(n,x,pbdir,inv,inv%model)
  IF(mype==0) CALL sub_modelinv2pbdirforwrinting(pbdir,inv,cpt_iter,1)!write only on master

100 continue

end subroutine FWI
