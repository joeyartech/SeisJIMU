!*********************************************************************
!  TEST CODE FOR IMPLEMENTING WASSERSTEIN LIKE DISTANCE IN FORTRAN   !    
! V0.0  - 03/2014  - L. Metivier                                     !
!--------------------------------------------------------------------!
!********************************************************************!

subroutine sub_prox_linear(N,L,x,gamma,y)
  
  implicit none
  
  !IN
  integer :: N
  real :: gamma
  real,dimension(N) :: L,x
  !IN/OUT
  real,dimension(N) :: y
  
  y(:)=x(:)-gamma*L(:)

  !TEST?
  !y(:)=gamma*x(:)-gamma*L(:)
  
end subroutine sub_prox_linear

subroutine sub_prox_linear_dble(N,L,x,gamma,y)
  
  implicit none
  
  !IN
  integer :: N
  double precision :: gamma
  double precision,dimension(N) :: L,x
  !IN/OUT
  double precision,dimension(N) :: y
  
  y(:)=x(:)-gamma*L(:)

end subroutine sub_prox_linear_dble




subroutine sub_prox_cube(NA,x,y)
  
  implicit none
  
  !IN
  integer :: NA
  real,dimension(NA) :: x
  !IN/OUT
  real,dimension(NA) :: y
  
  !Local variables
  integer :: i

  y(:)=x(:)
  do i=1,NA
     if(x(i)>1e0) then
        y(i)=1e0
     elseif(x(i)<-1e0) then
        y(i)=-1e0
     endif
  enddo
  
end subroutine sub_prox_cube


subroutine sub_prox_cube_dble(NA,x,y)
  
  implicit none
  
  !IN
  integer :: NA
  double precision,dimension(NA) :: x
  !IN/OUT
  double precision,dimension(NA) :: y
  
  !Local variables
  integer :: i

  y(:)=x(:)
  do i=1,NA
     if(x(i)>1d0) then
        y(i)=1d0
     elseif(x(i)<-1d0) then
        y(i)=-1d0
     endif
  enddo
  
end subroutine sub_prox_cube_dble
