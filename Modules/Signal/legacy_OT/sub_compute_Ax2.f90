!*********************************************************************
!  TEST CODE FOR IMPLEMENTING WASSERSTEIN LIKE DISTANCE IN FORTRAN   !    
! V0.0  - 03/2014  - L. Metivier                                     !
!--------------------------------------------------------------------!
!********************************************************************!

subroutine sub_compute_Ax2(NA,N,n1,n2,x,Ax)
  
  implicit none
  
  !IN
  integer :: NA,N,n1,n2
  real,dimension(N) :: x
  !IN/OUT
  real,dimension(NA) :: Ax
  
  !Local variables
  integer :: i1,i2,i,inc,k
  real :: h1,h2
  
  !Initialize
  Ax(:)=0e0
  
  h1=real(n1)
  h2=real(n2)
  
  k=1
  
  !Upper block (horizontal constraints)  
  do i2=1,n2-1         
     do i1=1,n1                             
        Ax(k)=-h2*x(i1+(i2-1)*n1)+h2*x(i1+i2*n1)        
        k=k+1
     enddo
  enddo
  
  !Middle block (vertical constraints)
  do i2=1,n2
     do i1=1,n1-1        
        Ax(k)=-h1*x(i1+(i2-1)*n1)+h1*x(i1+1+(i2-1)*n1)
        k=k+1
     enddo
  enddo
  
  !Lower block (identity)  
  do i=1,N     
     Ax(k)=x(i)
     k=k+1
  enddo
  
end subroutine sub_compute_Ax2

!**************************************!
! TEST
!**************************************!
subroutine sub_compute_Ax4(NA,N,n1,n2,scal,x,Ax)
  
  implicit none
  
  !IN
  integer :: NA,N,n1,n2
  real :: scal
  real,dimension(N) :: x
  !IN/OUT
  real,dimension(NA) :: Ax
  
  !Local variables
  integer :: i1,i2,i,inc,k
  real :: h1,h2
  
  !Initialize
  Ax(:)=0e0
  
  h1=real(n1)
  h2=real(n2)
  
  k=1
  
  !Upper block (horizontal constraints)  
  do i2=1,n2-1         
     do i1=1,n1                             
        Ax(k)=-h2*x(i1+(i2-1)*n1)+h2*x(i1+i2*n1)        
        k=k+1
     enddo
  enddo
  
  !Middle block (vertical constraints)
  do i2=1,n2
     do i1=1,n1-1        
        Ax(k)=-h1*x(i1+(i2-1)*n1)+h1*x(i1+1+(i2-1)*n1)
        k=k+1
     enddo
  enddo
  
  !Lower block (identity)  
  do i=1,N     
     Ax(k)=x(i)*scal
     k=k+1
  enddo
  
end subroutine sub_compute_Ax4

!!$!*********************************************************************
!!$!  TEST CODE FOR IMPLEMENTING WASSERSTEIN LIKE DISTANCE IN FORTRAN   !    
!!$! V0.1  - 03/2014  - L. Metivier                                     !
!!$!--------------------------------------------------------------------!
!!$!********************************************************************!
!!$
!!$subroutine sub_compute_Ax3(NA,N,n1,n2,x,Ax)
!!$  
!!$  implicit none
!!$  
!!$  !IN
!!$  integer :: NA,N,n1,n2
!!$  real,dimension(N) :: x
!!$  !IN/OUT
!!$  real,dimension(NA) :: Ax
!!$  
!!$  !Local variables
!!$  integer :: i1,i2,i
!!$  real :: h1,h2
!!$  
!!$  !Initialize
!!$  Ax(:)=0e0
!!$  
!!$  h1=real(n1)
!!$  h2=real(n2)
!!$  !h1=1.
!!$  !h2=1.
!!$  
!!$  !Upper block (vertical constraints)  
!!$  do i2=1,n2         
!!$     do i1=1,n1-1
!!$        Ax(i1+(i2-1)*n1)=-h1*x(i1+(i2-1)*n1)+h1*x(i1+1+(i2-1)*n1)        
!!$     enddo
!!$  enddo
!!$  
!!$  !Middle block (horizontal constraints)
!!$  do i2=1,n2-1
!!$     do i1=1,n1        
!!$        Ax(i1+(i2-1)*n1+n1*n2)=-h2*x(i1+(i2-1)*n1)+h2*x(i1+i2*n1)        
!!$     enddo
!!$  enddo
!!$  
!!$  !Lower block (identity)  
!!$  do i=1,n1*n2
!!$     Ax(i+2*n1*n2)=x(i)
!!$  enddo
!!$  
!!$end subroutine sub_compute_Ax3
