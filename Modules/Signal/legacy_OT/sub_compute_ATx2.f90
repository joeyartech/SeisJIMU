!*********************************************************************
!  TEST CODE FOR IMPLEMENTING WASSERSTEIN LIKE DISTANCE IN FORTRAN   !    
! V0.0  - 03/2014  - L. Metivier                                     !
!--------------------------------------------------------------------!
!********************************************************************!

subroutine sub_compute_ATx2(NA,N,n1,n2,x,ATx)
  
  implicit none
  
  !IN
  integer :: NA,N,n1,n2
  real,dimension(NA) :: x
  !IN/OUT
  real,dimension(N) :: ATx
  
  !Local variables
  integer :: i1,i2,i,k,dec1,dec2
  real :: h1,h2
  
  !Initialize
  ATx(:)=0e0
  
  h1=real(n1)
  h2=real(n2)
  
  k=1
  
  dec1=n1*(n2-1)
  dec2=n2*(n1-1)

  !FIRST BLOCK i2=1
  ATx(1)=&
       -h2*x(1)&
       -h1*x(1+dec1)&
       +x(1+dec1+dec2)
  do i1=2,n1-1
     ATx(i1)=&
          -h2*x(i1)&
          +h1*x(i1-1+dec1)&
          -h1*x(i1+dec1)&
          +x(i1+dec1+dec2)
  enddo
  ATx(n1)=&
       -h2*x(n1)&
       +h1*x(n1-1+dec1)&          
       +x(i1+dec1+dec2)
  
  k=n1
  !INTERNAL BLOCKS
  do i2=2,n2-1
     ATx(1+(i2-1)*n1)=&
          +h2*x(1+(i2-2)*n1)&
          -h2*x(1+(i2-1)*n1)+&          
          -h1*x(k+dec1)&
          +x(1+(i2-1)*n1+dec1+dec2)
     k=k+1
     do i1=2,n1-1
        ATx(i1+(i2-1)*n1)=&
             +h2*x(i1+(i2-2)*n1)&
             -h2*x(i1+(i2-1)*n1)+&
             +h1*x(k-1+dec1)&
             -h1*x(k+dec1)&
             +x(i1+(i2-1)*n1+dec1+dec2)
        k=k+1
     enddo
     ATx(n1+(i2-1)*n1)=&
          +h2*x(n1+(i2-2)*n1)&
          -h2*x(n1+(i2-1)*n1)+&
          +h1*x(k-1+dec1)&          
          +x(n1+(i2-1)*n1+dec1+dec2)
  enddo
  
  !LAST BLOCK
  ATx(1+(n2-1)*n1)=&
       +h2*x(1+(n2-2)*n1)&       
       -h1*x(k+dec1)&
       +x(1+(n2-1)*n1+dec1+dec2)
  k=k+1
  do i1=2,n1-1
     ATx(i1+(i2-1)*n1)=&
          +h2*x(i1+(n2-2)*n1)&          
          +h1*x(k-1+dec1)&
          -h1*x(k+dec1)&
          +x(i1+(n2-1)*n1+dec1+dec2)
     k=k+1
  enddo
  ATx(n1+(i2-1)*n1)=&
       +h2*x(n1+(n2-2)*n1)&       
       +h1*x(k-1+dec1)&          
       +x(n1+(n2-1)*n1+dec1+dec2)
  
end subroutine sub_compute_ATx2


!DEBUG
subroutine sub_compute_ATx4(NA,N,n1,n2,scal,x,ATx)
  
  implicit none
  
  !IN
  integer :: NA,N,n1,n2
  real :: scal
  real,dimension(NA) :: x
  !IN/OUT
  real,dimension(N) :: ATx
  
  !Local variables
  integer :: i1,i2,i,k,dec1,dec2
  real :: h1,h2
  
  !Initialize
  ATx(:)=0e0
  
  h1=real(n1)
  h2=real(n2)
  
  k=1
  
  dec1=n1*(n2-1)
  dec2=n2*(n1-1)

  !FIRST BLOCK i2=1
  ATx(1)=&
       -h2*x(1)&
       -h1*x(1+dec1)&
       +scal*x(1+dec1+dec2)
  do i1=2,n1-1
     ATx(i1)=&
          -h2*x(i1)&
          +h1*x(i1-1+dec1)&
          -h1*x(i1+dec1)&
          +scal*x(i1+dec1+dec2)
  enddo
  ATx(n1)=&
       -h2*x(n1)&
       +h1*x(n1-1+dec1)&          
       !+x(i1+dec1+dec2)
       +scal*x(n1+dec1+dec2)
  
  k=n1
  !INTERNAL BLOCKS
  do i2=2,n2-1
     ATx(1+(i2-1)*n1)=&
          +h2*x(1+(i2-2)*n1)&
          -h2*x(1+(i2-1)*n1)+&          
          -h1*x(k+dec1)&
          +scal*x(1+(i2-1)*n1+dec1+dec2)
     k=k+1
     do i1=2,n1-1
        ATx(i1+(i2-1)*n1)=&
             +h2*x(i1+(i2-2)*n1)&
             -h2*x(i1+(i2-1)*n1)+&
             +h1*x(k-1+dec1)&
             -h1*x(k+dec1)&
             +scal*x(i1+(i2-1)*n1+dec1+dec2)
        k=k+1
     enddo
     ATx(n1+(i2-1)*n1)=&
          +h2*x(n1+(i2-2)*n1)&
          -h2*x(n1+(i2-1)*n1)+&
          +h1*x(k-1+dec1)&          
          +scal*x(n1+(i2-1)*n1+dec1+dec2)          
  enddo
  
  !LAST BLOCK
  ATx(1+(n2-1)*n1)=&
       +h2*x(1+(n2-2)*n1)&       
       -h1*x(k+dec1)&
       +scal*x(1+(n2-1)*n1+dec1+dec2)
  k=k+1
  do i1=2,n1-1
     !ATx(i1+(i2-1)*n1)=&
     ATx(i1+(n2-1)*n1)=&
          +h2*x(i1+(n2-2)*n1)&          
          +h1*x(k-1+dec1)&
          -h1*x(k+dec1)&
          +scal*x(i1+(n2-1)*n1+dec1+dec2)
     k=k+1
  enddo
  !ATx(n1+(i2-1)*n1)=&
  ATx(n1+(n2-1)*n1)=&
       +h2*x(n1+(n2-2)*n1)&       
       +h1*x(k-1+dec1)&          
       +scal*x(n1+(n2-1)*n1+dec1+dec2)
  
end subroutine sub_compute_ATx4
























subroutine sub_compute_ATx2_bis(NA,N,n1,n2,x,ATx)
  
  implicit none
  
  !IN
  integer :: NA,N,n1,n2
  real,dimension(NA) :: x
  !IN/OUT
  real,dimension(N) :: ATx
  
  !Local variables
  integer :: i1,i2,i,k,dec1,dec2
  real :: h1,h2
  
  !Initialize
  ATx(:)=0e0
  
  h1=real(n1)
  h2=real(n2)
  
  k=1
  
  dec1=n1*(n2-1)
  dec2=n2*(n1-1)
  
  !Upper block (horizontal constraints)
  i2=1
  do i1=1,n1
     ATx(i1+(i2-1)*n1)=&
          ATx(i1+(i2-1)*n1)&
          -h2*x(i1+(i2-1)*n1)
  enddo
  do i2=2,n2-1
     do i1=1,n1
        ATx(i1+(i2-1)*n1)=&
             ATx(i1+(i2-1)*n1)&
          -h2*x(i1+(i2-1)*n1)&
          +h2*x(i1+(i2-2)*n1)
     enddo
  enddo
  i2=n2
  do i1=1,n1
     ATx(i1+(i2-1)*n1)=&
          ATx(i1+(i2-1)*n1)&          
          +h2*x(i1+(i2-2)*n1)
  enddo
  
  !Middle block (vertical constraints)
  do i2=1,n2
     i1=1
     ATx(i1+(i2-1)*n1)=&
          ATx(i1+(i2-1)*n1)&
          -h1*x(i1+(i2-1)*(n1-1)+dec1)
     do i1=2,n1-1
        ATx(i1+(i2-1)*n1)=&
             ATx(i1+(i2-1)*n1)&
             -h1*x(i1+(i2-1)*(n1-1)+dec1)&
             +h1*x(i1-1+(i2-1)*(n1-1)+dec1)
     enddo
     i1=n1
     ATx(i1+(i2-1)*n1)=&
          ATx(i1+(i2-1)*n1)&
          +h1*x(i1-1+(i2-1)*(n1-1)+dec1)
  enddo
  
  !Lower block (identity)
  do i=1,N
     ATx(i)=ATx(i)+x(i+dec1+dec2)
  enddo
  
end subroutine sub_compute_ATx2_bis






























!*********************************************************************
!  TEST CODE FOR IMPLEMENTING WASSERSTEIN LIKE DISTANCE IN FORTRAN   !    
! V0.1  - 04/2015  - L. Metivier                                     !
!--------------------------------------------------------------------!
!********************************************************************!

subroutine sub_compute_ATx3(NA,N,n1,n2,x,ATx)
  
  implicit none
  
  !IN
  integer :: NA,N,n1,n2
  real,dimension(NA) :: x
  !IN/OUT
  real,dimension(N) :: ATx
  
  !Local variables
  integer :: i1,i2,i
  real :: h1,h2
  
  !Initialize
  ATx(:)=0e0
    
  h1=real(n1)
  h2=real(n2)
  
  !h1=1.
  !h2=1.

  !Upper block (vertical constraints)  
  do i2=1,n2
     !i1=1
     i1=1
     ATx(i1+(i2-1)*n1)=-h1*x(i1+(i2-1)*n1)
     do i1=2,n1-1
        ATx(i1+(i2-1)*n1)=-h1*x(i1+(i2-1)*n1)+h1*x(i1-1+(i2-1)*n1)        
     enddo
     i1=n1
     ATx(i1+(i2-1)*n1)=h1*x(i1-1+(i2-1)*n1)
  enddo
  !Middle block (horizontal constraints)
  i2=1
  do i1=1,n1        
     ATx(i1+(i2-1)*n1)=&
          ATx(i1+(i2-1)*n1)-h2*x(i1+(i2-1)*n1)
  enddo
  do i2=2,n2-1
     do i1=1,n1        
        ATx(i1+(i2-1)*n1)=&
             ATx(i1+(i2-1)*n1)-h2*x(i1+(i2-1)*n1)+h2*x(i1+(i2-2)*n1)
     enddo
  enddo
  i2=n2
  do i1=1,n1        
     ATx(i1+(i2-1)*n1)=&
          ATx(i1+(i2-1)*n1)+h2*x(i1+(i2-2)*n1)
  enddo
  !Lower block (identity)  
  do i=1,n1*n2
     ATx(i)=ATx(i)+x(i)
  enddo
  

  
end subroutine sub_compute_ATx3
