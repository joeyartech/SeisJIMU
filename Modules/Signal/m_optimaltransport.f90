module m_KROT
!Kantorovic-Robinstein Wasserstein-1 norm 
!An optimal transport based distance measurement
!from legacy Seiscope code
!written by L. Métivier, modified by W. Zhou 
use m_proxy

	private

	public :: krot

	type,public :: t_KROT
		!residual
		integer :: nt_resamp

		!SDMM
		integer :: ichoice_dmm
		integer :: niter_max
		real :: conv
		integer :: ifpost_sdmm_weight

		integer :: ierr

		!Wasserstein distance
		real :: scale_residual,scal_was
		real :: bound_adjoint
		real :: const
		real :: taper_rfl

		!  real,dimension(:),pointer :: d
		!  real,dimension(:,:),pointer :: ata,ata_sparse,z,zt
		!  real,dimension(:,:),pointer :: residuals_l2save,residuals_wsave,data_cal_save
		!  double precision,dimension(:),pointer :: dd
		!  double precision,dimension(:,:),pointer :: dz,dzt

		!  !3D Wasserstein distance
		!  integer :: nrec_x,nrec_y
		!  integer,dimension(:),pointer :: a3d_row,a3d_row_csr,a3d_col
		!  double precision,dimension(:),pointer :: a3d

	end type

	type(t_KROT),public :: ot

	subroutine krot(shot%dadj)
		call SDMM_fishpack(ot,inv,pbdir%nt,acqui%nrec)    !ZW FFT solver for poisson eq, recommended
        ! call SDMM_mudpack(ot,inv,pbdir%nt,acqui%nrec)     !ZW multigraid solver for poisson eq.
		! call SDMM_mudpack_3d(ot,inv,pbdir%nt,acqui%nrec)  !ZW multigrid solver for poisson eq
	end subroutine


	subroutine SDMM_fishpack(ot,inv,nt,nrec)
	!V0 - 03/2014 - L. Métivier:
	!	Test code. Transpose the seismogram to reduce complexity in our case (nrcv < nt)
	!V3 - Jan14/2017 - W. Zhou: 
	!	Detected a pointer issue in fishpack when n1 equal certain numbers
	!   cannot fix this bug due to the old Fortran style.
	!V4 - Jan20/2026 - W. Zhou:
	!   Imported to SeisJIMU

		integer :: check_n1(20) = [ 32 , 64 , 96 , 128 , 192 , 255 , 256 , 384 , 511 , 512 , 640 , 768 , 1023 , 1024 , 1280 , 1535 , 1536 , 1792 , 1920 , 1984 ]
		!  2016 , 2032 , 2040 , 2044 , 2045 , 2046 , 2047 , 2048 , 2560 , 3071 , 3072  !to be added...

		real,dimension(:),allocatable :: L,P
		real,dimension(:),allocatable :: x,xt
		real,dimension(:),allocatable :: y1,y2,z1,z2,d1,d2,ATd2
		real,dimension(:),allocatable :: s1,s2,s1temp,s2temp  
		real,dimension(:),allocatable :: fcost_profile

		!config FISHPACK
		double precision,parameter :: A=0d0, B=1d0, C=0d0, D=1d0, PERTRB=0d0
		integer,parameter :: MBDCND=3, NBDCND=3
		double precision :: ELMBDA

		double precision,dimension(:),allocatable :: BDA,BDB,BDC,BDD
		double precision,dimension(:,:),allocatable :: F

		if(any( (check_n1-shot%nrcv)==0 )) then
			call warn(shot%sindex//' has a dangerous shot%nrcv ('//num2str(shot%nrcv)'). FISHPACK may crash.',mpiworld%iproc)
		end if

		!Dimensions
		n1=nrec
		n2=nt
		N=n1*n2
		NA=N+n1*(n2-1)+n2*(n1-1)

		Mfish=n2-1 !nt-1
		Nfish=n1-1 !nrec-1
		IDIMF=Mfish+1

		!Allocation
		call alloc(xt,N);  call alloc(x,N)
		call alloc(y1,N ); call alloc(z1,N )
		call alloc(y2,NA); call alloc(z2,NA)
		call alloc(d1,N);  call alloc(d2,NA)

		call alloc(ATd2,N)

		call alloc(s1,N); 		call alloc(s1temp,N)
		call alloc(s2,NA);		call alloc(s2temp,NA)

		call alloc(L,N)
		do i2=1,n2
		do i1=1,n1                
			L(i1+(i2-1)*n1)=shot%dadj(i2,i1)
		enddo
		enddo

		!   open(10,file='sismo_L',access='direct',recl=4*N)
		!   read(10,rec=1) L
		!   close(10)

		gamma=0.9  !magic parameter
		!gamma=1.5

		call alloc(BDA,Nfish+1)
		call alloc(BDB,Nfish+1)

		call alloc(BDC,Mfish+1)
		call alloc(BDD,Mfish+1)

		ELMBDA=dble(-1./ot%bound_adjoint/ot%bound_adjoint-1.)

		call alloc(F,Mfish+1,Nfish+1)

		fcost=0.
		fcost_prev=0.

		!allocate(fcost_profile(inv%niter_max))
		!Main Loop 
		do i=1,ot%niter_max

			!Set to 0
			xt=0.

			!PROX1 - update RHS
			d1 = y1 - z1
			xt = xt + d1

			!PROX2 - update RHS
			d2 = y2 - z2
			ATd2 = ATx(d2,NA,N,n1,n2)
			xt = xt + ATd2

			!CALL FISHPACK HWSCRT
			do i2=1,nrec
			do i1=1,nt
				F(i1,i2)=-xt(i2+(i1-1)*nrec) !note the minus sign: we are solving min (-F) instead of max F
			enddo
			enddo

			call cpu_time(tic)
			call HWSCRT(A,B,Mfish,MBDCND,BDA,BDB,C,D,Nfish,NBDCND,BDC,BDD,&
						ELMBDA,F,IDIMF,PERTRB,ot%ierr)
			call cpu_time(toc)

			if(mype==0) then
				write(*,'(A,I3)')  ' ITER # :',i
				if(i==1)     write(*,*) '   - Elapsed time for FISHPACK :', toc-tic, ' s'   
				if(i==1)     write(*,*) '   - PERTRB, ot%ierr :', sngl(PERTRB), ot%ierr   
				write(*,*) '   - Inner misfit :', fcost
			end if

			do i2=1,nrec
			do i1=1,nt
				x(i2+(i1-1)*nrec)=F(i1,i2)
			enddo
			enddo

			if(mype==0) open(10,file='F_fishpack',access='direct',recl=4*n1*n2)
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

			!PROX 1 - apply
			s1 = x
			s1temp = s1 + z1
			y1 = proxy_linear(s1temp,gamma,L,N)
			z1 = z1 + s1 - y1

			!PROX 2 - apply
			s2 = Ax(x,NA,N,n1,n2)
			s2temp = s2 + z2
			y2 = proxy_cube(s2temp,NA)
			z2 = z2 + s2 - y2

		enddo

		if (i==ot%niter_max+1) write(*,'(A,I4,A,E10.5,A,E10.5)')     ' Proc #',mype, ' MAX ITERation reached.  Inner misfit :', fcost, '; Norm of F :',maxval(abs(F))

		!open(50,file='fcost_profile',access='direct',recl=4*ot%niter_max)
		!write(50,rec=1) fcost_profile(:)
		!close(50)

		!compute fcost = \int x residual 
		fobj%misfit = -sum(x*L) *ot%nt_resamp
		!   write(*,*) 'inv%fcost in sdmm is : ',inv%fcost
		!   open(10,file='res_L2',access='direct',recl=4*n1*n2)
		!   write(10,rec=1) inv%residual_p(:,:)
		!   close(10)
	  
	  
		!Extrapolate residuals  
		do i2=1,n2
		do i1=1,n1        
			shot%dadj(i2,i1)=-x(i1+(i2-1)*n1)
		enddo
		enddo

	end subroutine

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


	! subroutine ot_sdmm_mudpack(ot,inv,nt,nrec)
	! !V1 - 11/2015 - L. Métivier:
	! !	Test code to implement 3D Wasserstein-like distance
	! !Vx - Jan14/2017 - W. Zhou: 
	! !	this subroutine is not ready so far 
	! !V4 - Jan20/2026 - W. Zhou:
	! !   Imported to SeisJIMU

	! 	implicit none

	! 	include 'inversion.h'
	! 	include 'optimaltransport.h'
	! 	include 'common.h'

	! 	!IN
	! 	integer :: nt,nrec
	! 	TYPE(optimaltransport) :: ot
	! 	!IN/OUT
	! 	TYPE(fullwaveinv) :: inv     
	! 	!Local variables
	! 	integer :: NA,N,m,niter,i,i1,i2,n1,n2,k
	! 	real :: gamma,t2,t1,tol,fcost,fcost_prev
	! 	real,dimension(:),allocatable :: L,P
	! 	real,dimension(:),allocatable :: x,xt
	! 	real,dimension(:),allocatable :: y1,y2,z1,z2,d1,d2,ATd2
	! 	real,dimension(:),allocatable :: s1,s2,s1temp,s2temp  
	! 	real,dimension(:),allocatable :: fcost_profile

	! 	!Local mudpack variables
	! 	integer :: ierror
	! 	integer :: iparm(17), mgopt(4)   !input parameters, multigrid options
	! 	double precision :: fparm(8) !float input parameters
	! 	double precision,dimension(:),allocatable :: work
	! 	double precision,dimension(:,:),allocatable :: rhs,F
	! 	external :: bndyc,sigt,sigp,xlmbda

	! 	!Dimensions  
	! 	n1=1537  !nt    
	! 	n2=129   !nrec

	! 	N=n1*n2
	! 	NA=(n1-1)*n2+n1*(n2-1)+N

	! 	!Allocation
	! 	allocate(xt(N))
	! 	xt(:)=0e0
	! 	allocate(x(N))
	! 	x(:)=0e0
	! 	allocate(y1(N),z1(N))
	! 	allocate(y2(NA),z2(NA))
	! 	y1(:)=0e0
	! 	z1(:)=0e0
	! 	y2(:)=0e0
	! 	z2(:)=0e0
	! 	allocate(d1(N),d2(NA))
	! 	d1(:)=0e0
	! 	d2(:)=0e0
	! 	allocate(ATd2(N))
	! 	ATd2(:)=0e0
	! 	allocate(s1(N),s1temp(N),s2(NA),s2temp(NA))
	! 	s1(:)=0e0
	! 	s2(:)=0e0
	! 	s1temp(:)=0e0
	! 	s2temp(:)=0e0


	! 	!Allocation standard residuals: extraction of the residual each n_samp_t time steps
	! 	allocate(L(N))
	! 	k=1
	! 	do i2=1,n2
	! 	 do i1=1,n1                
	! 	    L(i1+(i2-1)*n1)=inv%residual_p(i1,k)           
	! 	 enddo
	! 	 k=k+1
	! 	enddo

	! 	!   open(10,file='sismo_L',access='direct',recl=4*N)
	! 	!   read(10,rec=1) L
	! 	!   close(10)


	! 	gamma=0.9  !the magic parameter
	! 	!gamma=1.5

	! 	!CONFIG MUDPACK 3D  
	! 	iparm(1)=0     !call for initialization and descritization 
	! 	iparm(2)=2     !Neumann boundary conditions
	! 	iparm(3)=2     !Neumann boundary conditions
	! 	iparm(4)=2     !Neumann boundary conditions
	! 	iparm(5)=2     !Neumann boundary conditions
	! 	iparm(6)=3  !5     !ixp, Coarsest grid in n1
	! 	iparm(7)=2  !5     !jyq, Coarsest grid in n2
	! 	iparm(8)=10   !8     !iex, Exponent for finest grid in n1
	! 	iparm(9)=7   !5     !jey, Exponent for finest grid in n2
	! 	iparm(10)=3*(2**(10-1))+1   !n1   !n1=nx=ixp*(2^(iex-1))+1
	! 	iparm(11)=2*(2**(7-1))+1  !n2   !n2=ny=jyq*(2^(jey-1))+1
	! 	print*,iparm(10),iparm(11)
	! 	iparm(12)=0    !No initial guess
	! 	iparm(13)=2    !maximum number of cycles
	! 	iparm(14)=1    !Relaxation on n1 (if n1 >> n2)
	! 	iparm(15)=ceiling( (n1*n2*(10+5+0)+8*(n1+n2+2)) *1.3333334 ) ! memory size
	! 	print*,iparm(15)
	! 	fparm(1)=0e0   !x1_min 
	! 	fparm(2)=1e0   !x1_max 
	! 	fparm(3)=0e0   !x2_min 
	! 	fparm(4)=1e0   !x2_max 
	! 	fparm(5)=0e0   !No error control on cycling. (Maximum cycle stopping criterion iparm(18))
	! 	allocate(work(iparm(16))) !working array  

	! 	!Default options for multigrid solver
	! 	mgopt(1)=2 ! w cycle by defaut
	! 	mgopt(2)=2 ! pre-relaxation   (default 2)
	! 	mgopt(3)=1 ! post-relaxaction (default 1)
	! 	mgopt(4)=3 ! interpolation (1 = linear, 3=cubic)
	! 	!mgopt(1)=1 ! w cycle by defaut
	! 	!mgopt(2)=2 ! pre-relaxation   (default 2)
	! 	!mgopt(3)=1 ! post-relaxaction (default 1)
	! 	!mgopt(4)=1 ! interpolation (1 = linear, 3=cubic)

	! 	allocate(rhs(n1,n2))  !right hand side
	! 	allocate(F(n1,n2))  !solution and initial guess

	! 	!DISCRETIZATION   
	! 	call mud2sa(iparm,fparm,work,sigt,sigp,xlmbda,bndyc,rhs,F,mgopt,ierror)        
	! 	if(mype==0)write(*,*) 'MUDPACK: ierror = ',ierror
	! 	if(mype==0)write(*,*) 'MUDPACK: memory required : ',iparm(15), iparm(16)

	! 	IF(ierror.ne.0) THEN
	! 	if(mype==0) then
	! 	   write(*,*) '********************************************'
	! 	   write(*,*) '**** FAILURE OF MUDPACK INTITIALIZATION ****'
	! 	   write(*,*) '****        CODE STOP NOW !!!           ****'
	! 	   write(*,*) '********************************************'
	! 	   CALL MPI_ABORT()
	! 	end if
	! 	END IF

	! 	iparm(1)=1 !call for real work

	! 	fcost=0.
	! 	fcost_prev=0.

	! 	if(mype==0) open(10,file='F_mudpack',access='direct',recl=4*n1*n2)

	! 	!Main Loop 
	! 	do i=1,ot%niter_max
		 
	! 	 !Set to 0
	! 	 xt(:)=0e0
		 
	! 	 !PROX1  !ZW update RHS
	! 	 d1(:)=y1(:)-z1(:)     
	! 	 xt(:)=xt(:)+d1(:)     
		 
	! 	 !PROX2  !ZW update RHS
	! 	 d2(:)=y2(:)-z2(:)
	! 	 call sub_compute_ATx2(NA,N,n1,n2,d2,ATd2)    
	! 	 xt(:)=xt(:)+ATd2(:)
		 
	! 	 if(i.eq.1) then
	! 	    x(:)=xt(:)
	! 	 else
	! 	    !CALL MUDPACK
	! 	    F(:,:)=0d0
	! 	    rhs(:,:)=0d0
	! 	    do i2=1,n2
	! 	       do i1=1,n1     
	! 	          rhs(i1,i2)=-xt(i1+(i2-1)*n1)
	! 	       enddo
	! 	    enddo

	! 	    call CPU_time(t1)  
	! 	    call mud3sa(iparm,fparm,work,sigt,sigp,xlmbda,bndyc,rhs,F,mgopt,ierror)        
	! 	    call CPU_time(t2)

	! 	    if(mype==0) then
	! 	      write(*,'(A,I3)')  ' ITER # :',i
	! 	      write(*,*) '   - Elapsed time for MUDPACK :', t2-t1, ' s'
	! 	      write(*,*) '   - ierror :', ierror
	! 	      write(*,*) '   - Inner misfit :', fcost
	! 	    end if
		    
	! 	    do i2=1,n2
	! 	       do i1=1,n1           
	! 	          x(i1+(i2-1)*n1)=F(i1,i2)
	! 	       enddo
	! 	    enddo
	! 	 endif

	! 	 if(mype==0) then
	! 	    write(10,rec=i) sngl(-F)
	! 	    write(*,*) '   - Norm of F on master :',maxval(abs(F))
	! 	 end if

	! 	 !fcost_profile(i)=0.
	! 	 !do k=1,n1*n2
	! 	 !   fcost_profile(i)=fcost_profile(i)-x(k)*L(k)
	! 	 !enddo
	! 	 fcost_prev=fcost
	! 	 fcost=-sum(x*L)
	! 	 if (ot%conv>0..and.i>10) then
	! 	    if (abs( (fcost-fcost_prev)/fcost_prev ).le.ot%conv) then
	! 	                   write(*,'(A,I4,A,I3,A,E10.4,A,E10.4)') ' Proc #',mype, ' CONVerged at iter ',i,'. Inner misfit :', fcost, '; Norm of F :',maxval(abs(F))
	! 	       exit
	! 	    end if
	! 	 end if
	! 	 !if(i>1) then
	! 	 !   diff=fcost_profile(i)-fcost_profile(i-1)
	! 	 !   write(*,*) 'diff is : ',diff
	! 	 !   write(*,*) 'diff/fcost :',diff/fcost_profile(i)
	! 	 !endif

	! 	 !PROX 1  !ZW APPLY 
	! 	 s1(:)=x(:)
	! 	 s1temp(:)=s1(:)+z1(:)     
	! 	 call sub_prox_linear(N,L,s1temp,gamma,y1)     
	! 	 z1(:)=z1(:)+s1(:)-y1(:)
		 
	! 	 !PROX 2  !ZW APPLY
	! 	 call sub_compute_Ax2(NA,N,n1,n2,x,s2)     
	! 	 s2temp(:)=s2(:)+z2(:)     
	! 	 call sub_prox_cube(NA,s2temp,y2)     
	! 	 z2(:)=z2(:)+s2(:)-y2(:)  
		 
	! 	enddo
	! 	if (i==ot%niter_max+1) write(*,'(A,I4,A,E10.4,A,E10.4)')     ' Proc #',mype, ' MAX ITERation reached.  Inner misfit :', fcost, '; Norm of F :',maxval(abs(F))

	! 	if(mype==0) close(10)


	! 	!COMPUTE FCOST AS \int(x,residual)
	! 	inv%fcost=-sum(x*L)

	! 	inv%fcost=inv%fcost*ot%nt_resamp

	! 	!Extrapolate residuals 
	! 	k=1 
	! 	do i2=1,n2
	! 	 do i1=1,n1        
	! 	    inv%residual_p(i1,k)=-x(i1+(i2-1)*n1) 
	! 	 enddo
	! 	 k=k+1
	! 	enddo

	! 	!Deallocation
	! 	deallocate(xt,x)
	! 	deallocate(y1,z1)
	! 	deallocate(y2,z2)
	! 	deallocate(d1,d2,ATd2)
	! 	deallocate(s1,s1temp,s2,s2temp)
	! 	deallocate(work)  
	! 	deallocate(rhs)
	! 	deallocate(F)

	! 	contains

	! 	!coef for theta derivative (mud3sa will call sigt off grid)
	! 	function sigt(t,p)
	! 		implicit none
	! 		double precision sigt,t,p
	! 		sigt = 1d0
	! 		return
	! 	end function

	! 	!coef for phi derivative (mud3sa will call sigp off grid)
	! 	function sigp(t,p)
	! 		implicit none
	! 		double precision sigp,t,p
	! 		sigp = 1d0
	! 		return
	! 	end function

	! 	function xlmbda(t,p) 
	! 		implicit none
	! 		include 'optimaltransport.h'
	! 		TYPE(optimaltransport) :: ot
	! 		double precision xlmbda,t,p

	! 		xlmbda = dble(1./ot%bound_adjoint/ot%bound_adjoint+1.)

	! 		return
	! 	end function

	! 	subroutine bndyc(kbdy,torp,alfa,gbdy)
	! 		implicit none
	! 		integer(in) :: kbdy
	! 		double precision(inout) :: torp,alfa,gbdy
	! 		alfa=0d0
	! 		gbdy=0d0

	! 	end subroutine

	! end subroutine

end