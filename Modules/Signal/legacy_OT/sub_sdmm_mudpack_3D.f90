!***************************************************************
!  TEST CODE FOR IMPLEMENTING 3D WASSERSTEIN LIKE DISTANCE     !    
!             V1.0  - 11/2015  - L. Metivier                   !
!--------------------------------------------------------------!
!**************************************************************!

subroutine sub_sdmm_mudpack_3D(inv,nt,nrec_x,nrec_y) 
  
  implicit none
  
  include 'inversion.h'
  include 'common.h'
  
  !IN
  integer :: nt,nrec_x,nrec_y
  !IN/OUT
  TYPE(fullwaveinv) :: inv     
  !Local variables
  integer :: NA,N,m,niter,i,i1,i2,i3,n1,n2,n3,k
  integer :: i1_ex,i1_samp
  real :: gamma,t2,t1,tol,diff,tmp
  real,dimension(:),allocatable :: L,P,x0
  real,dimension(:),allocatable :: x,xt
  real,dimension(:),allocatable :: y1,y2,z1,z2,d1,d2,ATd2
  real,dimension(:),allocatable :: s1,s2,s1temp,s2temp  
  real,dimension(:),allocatable :: fcost_profile

  !Local mudpack variables
  integer :: ierror
  integer,dimension(:),allocatable :: iparm,mgopt
  double precision,dimension(:),allocatable :: fparm,work
  double precision,dimension(:,:,:),allocatable :: rhs,phi
  external :: bndyc,sigr,sigt,sigp,lam
    
  !Dimensions  
  n1=nt    
  n2=nrec_x
  n3=nrec_y
  
  N=n1*n2*n3
  NA=(n1-1)*n2*n3+n1*(n2-1)*n3+n1*n2*(n3-1)+N
  
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
  
  !################
  !## TEMP #######
  !################
  deallocate(inv%residual_p)
  allocate(inv%residual_p(nt,nrec_x*nrec_y))
  open(10,file='sismo_1D_cube_small',access='direct',recl=4*nt*nrec_x*nrec_y)
  read(10,rec=1) inv%residual_p(:,:)
  close(10)

  !Allocation standard residuals: extraction of the residual each n_samp_t time steps
  allocate(L(N))
  k=1
  do i3=1,n3         
     do i2=1,n2
        do i1=1,n1                
           L(i1+(i2-1)*n1+(i3-1)*n1*n2)=inv%residual_p(i1,k)           
        enddo
        k=k+1
     enddo
  enddo
  
  !Niter max  
  niter=inv%niter_WAS  
  gamma=0.9

  !CONFIG MUDPACK 3D  
  allocate(iparm(23))
  iparm(1)=0     !int1=0 for init call   
  iparm(2)=2     !Neumann boundary conditions
  iparm(3)=2     !Neumann boundary conditions
  iparm(4)=2     !Neumann boundary conditions
  iparm(5)=2     !Neumann boundary conditions
  iparm(6)=2     !Neumann boundary conditions
  iparm(7)=2     !Neumann boundary conditions
  iparm(8)=5     !Coarsest grid in n1
  iparm(9)=5     !Coarsest grid in n2
  iparm(10)=5    !Coarsest grid in n3
  iparm(11)=8    !Exponent for finest grid in n1
  iparm(12)=5    !Exponent for finest grid in n2
  iparm(13)=5    !Exponent for finest grid in n3
  iparm(14)=n1   !n1
  iparm(15)=n2   !n2
  iparm(16)=n3   !n3
  iparm(17)=0    !No initial guess
  iparm(18)=2    !maximum number of cycles
  iparm(19)=1    !relaxation strategy (if n1 >> n2, n1 >> n3 and n2~n3)
  iparm(20)=0    !special choice for relax strategy 7,8,9
  iparm(21)=(n1+2)*(n2+2)*(n3+2)*(10+5+0+0) ! memory size
  allocate(fparm(8))
  fparm(1)=0e0   !x1_min 
  fparm(2)=1e0   !x1_max 
  fparm(3)=0e0   !x2_min 
  fparm(4)=1e0   !x2_max 
  fparm(5)=0e0   !x3_min 
  fparm(6)=1e0   !x3_max 
  fparm(7)=0e0   !error control on cycling
  allocate(work(iparm(21))) !working array  
  allocate(mgopt(4))       !multigrid options => default options are chosen
  mgopt(1)=2 ! w cycle by defaut
  mgopt(2)=2 ! pre-relaxation   (default 2)
  mgopt(3)=1 ! post-relaxaction (default 1)
  mgopt(4)=3 ! interpolation (1 = linear, 3=cubic)
  !mgopt(1)=1 ! w cycle by defaut
  !mgopt(2)=2 ! pre-relaxation   (default 2)
  !mgopt(3)=1 ! post-relaxaction (default 1)
  !mgopt(4)=1 ! interpolation (1 = linear, 3=cubic)


  
  allocate(rhs(n1,n2,n3))  !right hand side
  allocate(phi(n1,n2,n3))  !solution and initial guess
  !allocate(rhs(n2,n3,n1))  !right hand side
  !allocate(phi(n2,n3,n1))  !solution and initial guess
  
  !DISCRETIZATION   
  call mud3sa(iparm,fparm,work,sigr,sigt,sigp,lam,bndyc,rhs,phi,&
       mgopt,ierror)        
  write(*,*) 'ierror is : ',ierror
  write(*,*) 'iparm(22) memory required : ',iparm(22)
  write(*,*) 'iparm(21) memory allocated: ',iparm(21)
  iparm(1)=1

  if(ierror.ne.0) goto 101

  write(*,*) 'INTO THE LOOP ! '

  !Main Loop 
  do i=1,niter
     
     !Set to 0
     xt(:)=0e0
     
     !PROX1
     d1(:)=y1(:)-z1(:)     
     xt(:)=xt(:)+d1(:)     
     
     !PROX2
     d2(:)=y2(:)-z2(:)
     write(*,*) 'BEFORE ATx2_3D '
     call sub_compute_ATx2_3D(NA,N,n1,n2,n3,d2,ATd2)    
     write(*,*) 'AFTER ATx2_3D '
     xt(:)=xt(:)+ATd2(:)
     
     if(i.ge.2) then
        !CALL MUDPACK
        phi(:,:,:)=0d0
        rhs(:,:,:)=0d0
        do i3=1,n3
           do i2=1,n2
              do i1=1,n1     
                 rhs(i1,i2,i3)=-xt(i1+(i2-1)*n1+(i3-1)*n1*n2)                         
                 !rhs(i2,i3,i1)=-xt(i1+(i2-1)*n1+(i3-1)*n1*n2)                         
              enddo
           enddo
        enddo
        call CPU_time(t1)  
        call mud3sa(iparm,fparm,work,sigr,sigt,sigp,lam,bndyc,rhs,phi,&
             mgopt,ierror)        
        call CPU_time(t2)
        write(*,*) 'IERROR : ',IERROR
        write(*,*) 't2-t1 : ',t2-t1, ' s'   
        do i3=1,n3
           do i2=1,n2
              do i1=1,n1           
                 x(i1+(i2-1)*n1+(i3-1)*n1*n2)=phi(i1,i2,i3)
                 !x(i1+(i2-1)*n1+(i3-1)*n1*n2)=phi(i2,i3,i1)
              enddo
           enddo
        enddo
     else
        x(:)=xt(:)
     endif

     !PROX 1   
     s1(:)=x(:)
     s1temp(:)=s1(:)+z1(:)     
     call sub_prox_linear(N,L,s1temp,gamma,y1)     
     z1(:)=z1(:)+s1(:)-y1(:)
     
     !PROX 2
     write(*,*) 'BEFORE Ax2_3D '
     call sub_compute_Ax2_3D(NA,N,n1,n2,n3,x,s2)     
     s2temp(:)=s2(:)+z2(:)     
     call sub_prox_cube(NA,s2temp,y2)     
     z2(:)=z2(:)+s2(:)-y2(:)  
     
  enddo
  !COMPUTE FCOST AS (x,residual)
  inv%fcost=0.
  do i=1,N
     inv%fcost=inv%fcost-x(i)*L(i)
  enddo
    
  !Save L2 res
  inv%residuals_L2save(:,:)=inv%residual_p(:,:)
  open(10,file='res_L2',access='direct',recl=4*n1*n2*n3)
  write(10,rec=1) inv%residual_p(:,:)
  close(10)  
  
  
  !Optimal transport residuals  
  k=1
  do i3=1,n3  
     do i2=1,n2
        do i1=1,n1        
           inv%residual_p(i1,k)=-x(i1+(i2-1)*n1+(i3-1)*n1*n2) 
        enddo
        k=k+1
     enddo
  enddo
  
  !Save W res
  inv%residuals_Wsave(:,:)=inv%residual_p(:,:)
  open(10,file='res_W_mudpack_3D',access='direct',recl=4*n1*n2*n3)
  write(10,rec=1) inv%residual_p(:,:)
  close(10)

101 return
  
end subroutine sub_sdmm_mudpack_3D

function sigr(r,t,p)
  
  !coefficient for r derivative (mud3sa will call sigr off grid)
  
  implicit none
  double precision sigr,r,t,p
  sigr = 1d0
  return
end function sigr

function sigt(r,t,p)
  
  !coefficient for theta derivative (mud3sa will call sigt off grid)
  
  implicit none
  double precision sigt,r,t,p
  sigt = 1d0
  return
end function sigt

function sigp(r,t,p)
  
  !coefficient for phi derivative (mud3sa will call sigp off grid)
  
  implicit none
  double precision sigp,r,t,p
  sigp = 1d0
  return
end function sigp

function lam(r,t,p)

  implicit none
  double precision lam,r,t,p
  lam = 2d0
  return
  
end function lam

subroutine bndyc(kbdy,xory,yorz,alfa,gbdy)
  
  implicit none
  !IN
  integer :: kbdy
  double precision :: xory,yorz
  !IN/OUT
  double precision :: alfa,gbdy
  
  alfa=0d0
  gbdy=0d0
  
end subroutine bndyc

