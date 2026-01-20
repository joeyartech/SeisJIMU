!*********************************************************************
!  TEST CODE FOR IMPLEMENTING WASSERSTEIN LIKE DISTANCE IN FORTRAN   !    
! V0.0  - 03/2014  - L. Metivier                                     !
!--------------------------------------------------------------------!
!********************************************************************!

!****************************************************************!
!*    TRANSPOSE TO REDUCE COMPLEXITY IN OUR CASE (NREC< NTIME STEPS)  *!
!****************************************************************!

! 2017.1.14 W. Zhou : detected a pointer issue in fishpack when n1 equal certain numbers
!                     cannot fix this bug due to the old Fortran fashion.

subroutine ot_sdmm_fishpack(ot,inv,nt,nrec) 
  
  implicit none
  
  include 'inversion.h'
  include 'optimaltransport.h'
  include 'common.h'
  
  !IN
  integer :: nt,nrec
  TYPE(optimaltransport) :: ot
  !IN/OUT
  TYPE(fullwaveinv) :: inv     

  integer :: check_n1(20) = (/ 32 ,  64 , 96 , 128 , 192 , 255 , 256 , 384 , 511 , 512 , 640 , 768 , 1023 , 1024 , 1280 , 1535 , 1536 , 1792 , 1920 , 1984 /)
!  2016 , 2032 , 2040 , 2044 , 2045 , 2046 , 2047 , 2048 , 2560 , 3071 , 3072  !to be added...

  !Local variables
  integer :: NA,N,m,i,i1,i2,n1,n2,k
  real :: gamma,t2,t1,tol,fcost,fcost_prev
  real,dimension(:),allocatable :: L,P
  real,dimension(:),allocatable :: x,xt
  real,dimension(:),allocatable :: y1,y2,z1,z2,d1,d2,ATd2
  real,dimension(:),allocatable :: s1,s2,s1temp,s2temp  
  real,dimension(:),allocatable :: fcost_profile

  !Local fishpack variables
  integer :: Mfish,Nfish,MBDCND,NBDCND,IDIMF,IERROR  
  double precision :: A,B,C,D,ELMBDA,PERTRB
  double precision,dimension(:),allocatable :: BDA,BDB,BDC,BDD
  double precision,dimension(:,:),allocatable :: F
  
  if(any( (check_n1-nrec)==0 )) then
    write(*,*)'*********************************************'
    write(*,*)'*** WARNING : FISHPACK MAY CRASH !!! ********'
    write(*,*)'*** PROC #',mype,' HAS DANGEROUS NO. OF RECV:',nrec
    write(*,*)'*********************************************'
  end if

  !Dimensions  
  n1=nrec    
  n2=nt  
  N=n1*n2  
  NA=N+n1*(n2-1)+n2*(n1-1)
  
  !Allocation
  allocate(xt(N))
  xt(:)=0e0
  allocate(x(N))
  x(:)=0e0
  allocate(y1(N),z1(N))
  allocate(y2(NA),z2(NA))
  y1(:)=0e0
  z1(:)=0e0
  y2(:)=0e0
  z2(:)=0e0
  allocate(d1(N),d2(NA))
  d1(:)=0e0
  d2(:)=0e0
  allocate(ATd2(N))
  ATd2(:)=0e0
  allocate(s1(N),s1temp(N),s2(NA),s2temp(NA))
  s1(:)=0e0
  s2(:)=0e0
  s1temp(:)=0e0
  s2temp(:)=0e0
  

  !Allocation standard residuals: extraction of the residual each n_samp_t time steps
  allocate(L(N))
  do i2=1,n2
     do i1=1,n1                
        L(i1+(i2-1)*n1)=inv%residual_p(i2,i1)
     enddo
  enddo

!   open(10,file='sismo_L',access='direct',recl=4*N)
!   read(10,rec=1) L
!   close(10)


  gamma=0.9  !the magic parameter
  !gamma=1.5

  !CONFIG FISHPACK
  A=0d0
  B=1d0
  Mfish=n2-1 !nt-1
  Nfish=n1-1 !nrec-1
  MBDCND=3
  allocate(BDA(Nfish+1))
  BDA(:)=0d0
  allocate(BDB(Nfish+1))
  BDB(:)=0d0
  C=0d0
  D=1d0
  NBDCND=3
  allocate(BDC(Mfish+1))
  BDC(:)=0d0
  allocate(BDD(Mfish+1))
  BDD(:)=0d0
  ELMBDA=dble(-1./ot%bound_adjoint/ot%bound_adjoint-1.)
  allocate(F(Mfish+1,Nfish+1))
  F(:,:)=0d0
  IDIMF=Mfish+1
  PERTRB=0d0

  fcost=0.
  fcost_prev=0.
  
  if(mype==0) open(10,file='F_fishpack',access='direct',recl=4*n1*n2)
  
  !allocate(fcost_profile(inv%niter_max))
  !Main Loop 
  do i=1,ot%niter_max
     
     !Set to 0
     xt(:)=0e0
     
     !PROX1  !ZW update RHS
     d1(:)=y1(:)-z1(:)     
     xt(:)=xt(:)+d1(:)     
     
     !PROX2  !ZW update RHS
     d2(:)=y2(:)-z2(:)
     call sub_compute_ATx2(NA,N,n1,n2,d2,ATd2)    
     xt(:)=xt(:)+ATd2(:)
     
     !CALL FISHPACK HWSCRT     
     do i2=1,nrec
        do i1=1,nt
           F(i1,i2)=-xt(i2+(i1-1)*nrec)           !note the minus sign: we are solving min (-F) instead of max F
        enddo
     enddo
     
     if(i==1)  call CPU_time(t1)        
     call HWSCRT (A,B,Mfish,MBDCND,BDA,BDB,C,D,Nfish,NBDCND,BDC,BDD,&
          ELMBDA,F,IDIMF,PERTRB,IERROR)
     if(i==1)  call CPU_time(t2)
     if(mype==0) then
        write(*,'(A,I3)')  ' ITER # :',i
        if(i==1)     write(*,*) '   - Elapsed time for FISHPACK :', t2-t1, ' s'   
        if(i==1)     write(*,*) '   - PERTRB, IERROR :', sngl(PERTRB), IERROR   
        write(*,*) '   - Inner misfit :', fcost
     end if
     
     do i2=1,nrec
        do i1=1,nt
           x(i2+(i1-1)*nrec)=F(i1,i2)
        enddo
     enddo
     
     if(mype==0) then
        write(10,rec=i) sngl(-F)
        write(*,*) '   - Norm of F on master :',maxval(abs(F))
     end if
     
     !fcost_profile(i)=0.
     !do k=1,n1*n2
     !   fcost_profile(i)=fcost_profile(i)-x(k)*L(k)
     !enddo
     fcost_prev=fcost
     fcost=-sum(x*L)
     if (ot%conv>0..and.i>10) then
        if (abs( (fcost-fcost_prev)/fcost_prev ).le.ot%conv) then
                       write(*,'(A,I4,A,I3,A,E10.5,A,E10.5)') ' Proc #',mype, ' CONVerged at iter ',i,'. Inner misfit :', fcost, '; Norm of F :',maxval(abs(F))
           exit
        end if
     end if
     !if(i>1) then
     !   diff=fcost_profile(i)-fcost_profile(i-1)
     !   write(*,*) 'diff is : ',diff
     !   write(*,*) 'diff/fcost :',diff/fcost_profile(i)
     !endif

     !PROX 1   !ZW APPLY
     s1(:)=x(:)
     s1temp(:)=s1(:)+z1(:)     
     call sub_prox_linear(N,L,s1temp,gamma,y1)     
     z1(:)=z1(:)+s1(:)-y1(:)
     
     !PROX 2   !ZW APPLY
     call sub_compute_Ax2(NA,N,n1,n2,x,s2)     
     s2temp(:)=s2(:)+z2(:)     
     call sub_prox_cube(NA,s2temp,y2)     
     z2(:)=z2(:)+s2(:)-y2(:)  


  enddo
  if (i==ot%niter_max+1) write(*,'(A,I4,A,E10.5,A,E10.5)')     ' Proc #',mype, ' MAX ITERation reached.  Inner misfit :', fcost, '; Norm of F :',maxval(abs(F))
  
  if(mype==0) close(10)
  
  !open(50,file='fcost_profile',access='direct',recl=4*ot%niter_max)
  !write(50,rec=1) fcost_profile(:)
  !close(50)
  
  !COMPUTE FCOST AS \int(x,residual)
  inv%fcost=-sum(x*L)
  
  inv%fcost=inv%fcost*ot%nt_resamp
!   write(*,*) 'inv%fcost in sdmm is : ',inv%fcost
!   open(10,file='res_L2',access='direct',recl=4*n1*n2)
!   write(10,rec=1) inv%residual_p(:,:)
!   close(10)
  
  
  !Extrapolate residuals  
  do i2=1,n2
     do i1=1,n1        
        inv%residual_p(i2,i1)=-x(i1+(i2-1)*n1)
     enddo
  enddo
  
  !Deallocation
  deallocate(xt,x)
  deallocate(y1,z1)
  deallocate(y2,z2)
  deallocate(d1,d2,ATd2)
  deallocate(s1,s1temp,s2,s2temp)
  deallocate(BDA,BDB,BDC,BDD,F)
  
end subroutine ot_sdmm_fishpack








































! !STRC VERSION
! subroutine sub_sdmm_fishpack_stcr(inv,nt,nrec) 
!   
!   implicit none
!   
!   include 'inversion.h'
!   include 'common.h'
!   
!   !IN
!   integer :: nt,nrec
!   !IN/OUT
!   TYPE(fullwaveinv) :: inv     
!   !Local variables
!   integer :: NA,N,m,i,i1,i2,n1,n2,k
!   integer :: i1_ex,i1_samp,stcr,cpt_iter
!   real :: gamma,t2,t1,tol,fcost,fcost_prev,diff
!   real,dimension(:),allocatable :: L,P,x0
!   real,dimension(:),allocatable :: x,xt,xm1
!   real,dimension(:),allocatable :: y1,y2,z1,z2,d1,d2,ATd2
!   real,dimension(:),allocatable :: s1,s2,s1temp,s2temp  
!   real,dimension(:),allocatable :: fcost_profile
! 
!   !Local fishpack variables
!   integer :: Mfish,Nfish,MBDCND,NBDCND,IDIMF,IERROR  
!   double precision :: A,B,C,D,ELMBDA,PERTRB
!   double precision,dimension(:),allocatable :: BDA,BDB,BDC,BDD!,W
!   double precision,dimension(:,:),allocatable :: F
!   
!   !Dimensions  
!   n1=nrec    
!   n2=nt  
!   N=n1*n2  
!   NA=N+n1*(n2-1)+n2*(n1-1)
!   
!   !Allocation
!   allocate(xt(N))
!   xt(:)=0e0
!   allocate(x(N))
!   x(:)=0e0
!   allocate(y1(N),z1(N))
!   allocate(y2(NA),z2(NA))
!   y1(:)=0e0
!   z1(:)=0e0
!   y2(:)=0e0
!   z2(:)=0e0
!   allocate(d1(N),d2(NA))
!   d1(:)=0e0
!   d2(:)=0e0
!   allocate(ATd2(N))
!   ATd2(:)=0e0
!   allocate(s1(N),s1temp(N),s2(NA),s2temp(NA))
!   s1(:)=0e0
!   s2(:)=0e0
!   s1temp(:)=0e0
!   s2temp(:)=0e0
!   allocate(xm1(N))
!   xm1(:)=0e0
!   
!   !Allocation standard residuals: extraction of the residual each n_samp_t time steps
!   allocate(L(N))
!   do i2=1,n2
!      do i1=1,n1                
!         L(i1+(i2-1)*n1)=inv%residual_p(i2,i1)
!      enddo
!   enddo
!   !Write Wasserstein residuals after extrapolation
! 
!   gamma=0.9
!   !gamma=1.5
! 
!   !CONFIG FISHPACK
!   A=0d0
!   B=1d0
!   Mfish=n2-1
!   Nfish=n1-1
!   MBDCND=3
!   allocate(BDA(Nfish+1))
!   BDA(:)=0d0
!   allocate(BDB(Nfish+1))
!   BDB(:)=0d0
!   C=0d0
!   D=1d0
!   NBDCND=3
!   allocate(BDC(Mfish+1))
!   BDC(:)=0d0
!   allocate(BDD(Mfish+1))
!   BDD(:)=0d0
!   ELMBDA=-2d0
!   allocate(F(Mfish+1,Nfish+1))
!   F(:,:)=0d0
!   IDIMF=Mfish+1
!   PERTRB=0d0
! 
!   fcost=0.
!   stcr=0
!   cpt_iter=0
!   !Main Loop 
!   do while (stcr.eq.0) 
!      
!      !Set to 0
!      xt(:)=0e0
!      
!      !PROX1
!      d1(:)=y1(:)-z1(:)     
!      xt(:)=xt(:)+d1(:)     
!      
!      !PROX2
!      d2(:)=y2(:)-z2(:)
!      call sub_compute_ATx2(NA,N,n1,n2,d2,ATd2)    
!      xt(:)=xt(:)+ATd2(:)
!      
!      !CALL FISHPACK HWSCRT     
!      F(:,:)=0d0
!      do i2=1,nrec
!         do i1=1,nt
!            F(i1,i2)=-xt(i2+(i1-1)*nrec)           
!         enddo
!      enddo
!      
!      BDA(:)=0d0
!      BDB(:)=0d0
!      BDC(:)=0d0
!      BDD(:)=0d0
!      call CPU_time(t1)        
!      call HWSCRT (A,B,Mfish,MBDCND,BDA,BDB,C,D,Nfish,NBDCND,BDC,BDD,&
!           ELMBDA,F,IDIMF,PERTRB,IERROR)     
!      call CPU_time(t2)
!      if(mype.eq.0) then
!         write(*,*) 't2-t1 : ',t2-t1, ' s'   
!      endif
!      do i2=1,nrec
!         do i1=1,nt
!            x(i2+(i1-1)*nrec)=F(i1,i2)
!         enddo
!      enddo
!      
!      !PROX 1   
!      s1(:)=x(:)
!      s1temp(:)=s1(:)+z1(:)     
!      call sub_prox_linear(N,L,s1temp,gamma,y1)     
!      z1(:)=z1(:)+s1(:)-y1(:)
!      
!      !PROX 2
!      call sub_compute_Ax2(NA,N,n1,n2,x,s2)     
!      s2temp(:)=s2(:)+z2(:)     
!      call sub_prox_cube(NA,s2temp,y2)     
!      z2(:)=z2(:)+s2(:)-y2(:)  
!      
!      cpt_iter=cpt_iter+1
! 
!      fcost_prev=fcost
!      fcost=0.
!      do k=1,n1*n2
!         fcost=fcost-x(k)*L(k)
!      enddo
!      diff=fcost-fcost_prev
!      if( ((diff/fcost)<1e-3).or.(cpt_iter.ge.ot%niter_max) ) then
!         stcr=1
!      endif
!      
!   enddo
!   if(mype.eq.0) then
!      write(*,*) 'niter_max sdmm : ',cpt_iter
!   endif
! 
!   !COMPUTE FCOST AS (x,residual)
!   inv%fcost=0.
!   do i=1,n1*n2
!      inv%fcost=inv%fcost-x(i)*L(i)
!   enddo
!   inv%fcost=inv%fcost*inv%n1_samp
!   
!   !Save L2 res
!   inv%residuals_L2save(:,:)=inv%residual_p(:,:)
!   
!   !Extrapolate residuals  
!   do i2=1,n2
!      do i1=1,n1        
!         inv%residual_p(i2,i1)=-x(i1+(i2-1)*n1)
!      enddo
!   enddo
!   
!   !Save W res
!   inv%residuals_Wsave(:,:)=inv%residual_p(:,:)
!   
!   
! end subroutine ot_sdmm_fishpack_stcr


