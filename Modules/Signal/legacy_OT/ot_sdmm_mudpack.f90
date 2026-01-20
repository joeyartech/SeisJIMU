!***************************************************************
!  TEST CODE FOR IMPLEMENTING 3D WASSERSTEIN LIKE DISTANCE     !    
!             V1.0  - 11/2015  - L. Metivier                   !
!--------------------------------------------------------------!
!**************************************************************!

subroutine ot_sdmm_mudpack(ot,inv,nt,nrec) !ZW this subroutine is not ready so far 
  
  implicit none
  
  include 'inversion.h'
  include 'optimaltransport.h'
  include 'common.h'
  
  !IN
  integer :: nt,nrec
  TYPE(optimaltransport) :: ot
  !IN/OUT
  TYPE(fullwaveinv) :: inv     
  !Local variables
  integer :: NA,N,m,niter,i,i1,i2,n1,n2,k
  real :: gamma,t2,t1,tol,fcost,fcost_prev
  real,dimension(:),allocatable :: L,P
  real,dimension(:),allocatable :: x,xt
  real,dimension(:),allocatable :: y1,y2,z1,z2,d1,d2,ATd2
  real,dimension(:),allocatable :: s1,s2,s1temp,s2temp  
  real,dimension(:),allocatable :: fcost_profile

  !Local mudpack variables
  integer :: ierror
  integer :: iparm(17), mgopt(4)   !input parameters, multigrid options
  double precision :: fparm(8) !float input parameters
  double precision,dimension(:),allocatable :: work
  double precision,dimension(:,:),allocatable :: rhs,F
  external :: bndyc,sigt,sigp,xlmbda
    
  !Dimensions  
  n1=1537  !nt    
  n2=129   !nrec
  
  N=n1*n2
  NA=(n1-1)*n2+n1*(n2-1)+N
  
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
  k=1
  do i2=1,n2
     do i1=1,n1                
        L(i1+(i2-1)*n1)=inv%residual_p(i1,k)           
     enddo
     k=k+1
  enddo
  
!   open(10,file='sismo_L',access='direct',recl=4*N)
!   read(10,rec=1) L
!   close(10)


  gamma=0.9  !the magic parameter
  !gamma=1.5

  !CONFIG MUDPACK 3D  
  iparm(1)=0     !call for initialization and descritization 
  iparm(2)=2     !Neumann boundary conditions
  iparm(3)=2     !Neumann boundary conditions
  iparm(4)=2     !Neumann boundary conditions
  iparm(5)=2     !Neumann boundary conditions
  iparm(6)=3  !5     !ixp, Coarsest grid in n1
  iparm(7)=2  !5     !jyq, Coarsest grid in n2
  iparm(8)=10   !8     !iex, Exponent for finest grid in n1
  iparm(9)=7   !5     !jey, Exponent for finest grid in n2
  iparm(10)=3*(2**(10-1))+1   !n1   !n1=nx=ixp*(2^(iex-1))+1
  iparm(11)=2*(2**(7-1))+1  !n2   !n2=ny=jyq*(2^(jey-1))+1
print*,iparm(10),iparm(11)
  iparm(12)=0    !No initial guess
  iparm(13)=2    !maximum number of cycles
  iparm(14)=1    !Relaxation on n1 (if n1 >> n2)
  iparm(15)=ceiling( (n1*n2*(10+5+0)+8*(n1+n2+2)) *1.3333334 ) ! memory size
print*,iparm(15)
  fparm(1)=0e0   !x1_min 
  fparm(2)=1e0   !x1_max 
  fparm(3)=0e0   !x2_min 
  fparm(4)=1e0   !x2_max 
  fparm(5)=0e0   !No error control on cycling. (Maximum cycle stopping criterion iparm(18))
  allocate(work(iparm(16))) !working array  

  !Default options for multigrid solver
  mgopt(1)=2 ! w cycle by defaut
  mgopt(2)=2 ! pre-relaxation   (default 2)
  mgopt(3)=1 ! post-relaxaction (default 1)
  mgopt(4)=3 ! interpolation (1 = linear, 3=cubic)
  !mgopt(1)=1 ! w cycle by defaut
  !mgopt(2)=2 ! pre-relaxation   (default 2)
  !mgopt(3)=1 ! post-relaxaction (default 1)
  !mgopt(4)=1 ! interpolation (1 = linear, 3=cubic)

  allocate(rhs(n1,n2))  !right hand side
  allocate(F(n1,n2))  !solution and initial guess
  
  !DISCRETIZATION   
  call mud2sa(iparm,fparm,work,sigt,sigp,xlmbda,bndyc,rhs,F,mgopt,ierror)        
  if(mype==0)write(*,*) 'MUDPACK: ierror = ',ierror
  if(mype==0)write(*,*) 'MUDPACK: memory required : ',iparm(15), iparm(16)

  IF(ierror.ne.0) THEN
    if(mype==0) then
       write(*,*) '********************************************'
       write(*,*) '**** FAILURE OF MUDPACK INTITIALIZATION ****'
       write(*,*) '****        CODE STOP NOW !!!           ****'
       write(*,*) '********************************************'
       CALL MPI_ABORT()
    end if
  END IF

  iparm(1)=1 !call for real work

  fcost=0.
  fcost_prev=0.

  if(mype==0) open(10,file='F_mudpack',access='direct',recl=4*n1*n2)

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
     
     if(i.eq.1) then
        x(:)=xt(:)
     else
        !CALL MUDPACK
        F(:,:)=0d0
        rhs(:,:)=0d0
        do i2=1,n2
           do i1=1,n1     
              rhs(i1,i2)=-xt(i1+(i2-1)*n1)
           enddo
        enddo

        call CPU_time(t1)  
        call mud3sa(iparm,fparm,work,sigt,sigp,xlmbda,bndyc,rhs,F,mgopt,ierror)        
        call CPU_time(t2)

        if(mype==0) then
          write(*,'(A,I3)')  ' ITER # :',i
          write(*,*) '   - Elapsed time for MUDPACK :', t2-t1, ' s'
          write(*,*) '   - ierror :', ierror
          write(*,*) '   - Inner misfit :', fcost
        end if
        
        do i2=1,n2
           do i1=1,n1           
              x(i1+(i2-1)*n1)=F(i1,i2)
           enddo
        enddo
     endif

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
                       write(*,'(A,I4,A,I3,A,E10.4,A,E10.4)') ' Proc #',mype, ' CONVerged at iter ',i,'. Inner misfit :', fcost, '; Norm of F :',maxval(abs(F))
           exit
        end if
     end if
     !if(i>1) then
     !   diff=fcost_profile(i)-fcost_profile(i-1)
     !   write(*,*) 'diff is : ',diff
     !   write(*,*) 'diff/fcost :',diff/fcost_profile(i)
     !endif

     !PROX 1  !ZW APPLY 
     s1(:)=x(:)
     s1temp(:)=s1(:)+z1(:)     
     call sub_prox_linear(N,L,s1temp,gamma,y1)     
     z1(:)=z1(:)+s1(:)-y1(:)
     
     !PROX 2  !ZW APPLY
     call sub_compute_Ax2(NA,N,n1,n2,x,s2)     
     s2temp(:)=s2(:)+z2(:)     
     call sub_prox_cube(NA,s2temp,y2)     
     z2(:)=z2(:)+s2(:)-y2(:)  
     
  enddo
  if (i==ot%niter_max+1) write(*,'(A,I4,A,E10.4,A,E10.4)')     ' Proc #',mype, ' MAX ITERation reached.  Inner misfit :', fcost, '; Norm of F :',maxval(abs(F))

  if(mype==0) close(10)


  !COMPUTE FCOST AS \int(x,residual)
  inv%fcost=-sum(x*L)
    
  inv%fcost=inv%fcost*ot%nt_resamp

  !Extrapolate residuals 
  k=1 
  do i2=1,n2
     do i1=1,n1        
        inv%residual_p(i1,k)=-x(i1+(i2-1)*n1) 
     enddo
     k=k+1
  enddo
  
  !Deallocation
  deallocate(xt,x)
  deallocate(y1,z1)
  deallocate(y2,z2)
  deallocate(d1,d2,ATd2)
  deallocate(s1,s1temp,s2,s2temp)
  deallocate(work)  
  deallocate(rhs)
  deallocate(F)
  
end subroutine ot_sdmm_mudpack














function sigt(t,p) !coefficient for theta derivative (mud3sa will call sigt off grid)
  implicit none
  double precision sigt,t,p
  sigt = 1d0
  return
end function sigt


function sigp(t,p) !coefficient for phi derivative (mud3sa will call sigp off grid)
  implicit none
  double precision sigp,t,p
  sigp = 1d0
  return
end function sigp


function xlmbda(t,p) 
   implicit none
  include 'optimaltransport.h'
  TYPE(optimaltransport) :: ot
  double precision xlmbda,t,p
  
  xlmbda = dble(1./ot%bound_adjoint/ot%bound_adjoint+1.)
  return
  
end function xlmbda


subroutine bndyc(kbdy,torp,alfa,gbdy)
  implicit none
  !IN
  integer :: kbdy
  !IN/OUT
  double precision :: torp,alfa,gbdy
  alfa=0d0
  gbdy=0d0
end subroutine bndyc

