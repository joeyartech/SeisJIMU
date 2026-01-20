!*********************************************************************
!  TEST CODE FOR IMPLEMENTING WASSERSTEIN LIKE DISTANCE IN FORTRAN   !    
! V0.0  - 03/2014  - L. Metivier                                     !
!--------------------------------------------------------------------!
!********************************************************************!

subroutine sub_compute_ATAx(N,n1,ATA_sparse,x,ATAx)
  
  implicit none
  
  !IN
  integer :: N,n1
  real,dimension(N) :: x
  real,dimension(3,N) :: ATA_sparse
  !IN/OUT
  real,dimension(N) :: ATAx
  
  !Local variables
  integer :: i1,i2,i,inc,k
  real :: h1,h2,c1,c2,alpha
  
  ATAx(:)=0e0  
  do i=1,N
     ATAx(i)=ATA_sparse(1,i)*x(i)
  enddo
  
  do i=1,N-1
     ATAx(i)=ATAx(i)+ATA_sparse(2,i)*x(i+1)
  enddo
  
  do i=2,N
     ATAx(i)=ATAx(i)+ATA_sparse(2,i-1)*x(i-1)
  enddo
  
  do i=1,N-n1
     !ATAx(i)=ATA(n1+1,i)*x(i+n1)
     ATAx(i)=ATAx(i)+ATA_sparse(3,i)*x(i+n1)
  enddo
  
  do i=n1+1,N
     !ATAx(i)=ATA(n1+1,i)*x(i-n1)
     ATAx(i)=ATAx(i)+ATA_sparse(3,i-n1)*x(i-n1)
  enddo
   
end subroutine sub_compute_ATAx


!BETTER IMPLEMENTATION
subroutine sub_compute_ATAx2(N,n1,ATA_sparse,x,ATAx)
  
  implicit none
  
  !IN
  integer :: N,n1
  real,dimension(N) :: x
  real,dimension(N,3) :: ATA_sparse
  !IN/OUT
  real,dimension(N) :: ATAx
  
  !Local variables
  integer :: i1,i2,i,inc,k
  real :: h1,h2,c1,c2,alpha
  
  ATAx(:)=0e0  
  do i=1,N
     ATAx(i)=ATA_sparse(i,1)*x(i)
  enddo
  
  do i=1,N-1
     ATAx(i)=ATAx(i)+ATA_sparse(i,2)*x(i+1)
  enddo
  
  do i=2,N
     ATAx(i)=ATAx(i)+ATA_sparse(i-1,2)*x(i-1)
  enddo
  
  do i=1,N-n1
     ATAx(i)=ATAx(i)+ATA_sparse(i,3)*x(i+n1)
  enddo
  
  do i=n1+1,N
     ATAx(i)=ATAx(i)+ATA_sparse(i-n1,3)*x(i-n1)
  enddo
   
end subroutine sub_compute_ATAx2

