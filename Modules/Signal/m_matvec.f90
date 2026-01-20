!*********************************************************************
! TEST CODE FOR IMPLEMENTING WASSERSTEIN LIKE DISTANCE IN FORTRAN   !    
!********************************************************************!
module matvec

    subroutine sub_compute_Ax2(NA,N,n1,n2,x,Ax)
    ! V0.0  - 03/2014  - L. Metivier 
        integer,intent(in) :: NA,N,n1,n2
        real,dimension(N),intent(in) :: x
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
      
    end subroutine

    ! test
    subroutine sub_compute_Ax4(NA,N,n1,n2,scal,x,Ax)

        integer,intent(in) :: NA,N,n1,n2
        real,intent(in) :: scal
        real,dimension(N),intent(in) :: x
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


    subroutine sub_compute_ATx2(NA,N,n1,n2,x,ATx)
    ! V0.0  - 03/2014  - L. Metivier
        integer,intent(in) :: NA,N,n1,n2
        real,dimension(NA),intent(in) :: x
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
      
    end subroutine


    !devug
    subroutine sub_compute_ATx4(NA,N,n1,n2,scal,x,ATx)
        integer,intent(in) :: NA,N,n1,n2
        real,intent(in) :: scal
        real,dimension(NA),intent(in) :: x
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
        integer,intent(in):: NA,N,n1,n2
        real,dimension(NA),intent(in) :: x
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
            ATx(i1+(i2-1)*n1) &
          -h2*x(i1+(i2-1)*n1)
        enddo
        do i2=2,n2-1
        do i1=1,n1
            ATx(i1+(i2-1)*n1)=&
            ATx(i1+(i2-1)*n1) &
          -h2*x(i1+(i2-1)*n1) &
          +h2*x(i1+(i2-2)*n1)
        enddo
        enddo

        i2=n2
        do i1=1,n1
            ATx(i1+(i2-1)*n1)=&
            ATx(i1+(i2-1)*n1) & 
          +h2*x(i1+(i2-2)*n1)
        enddo
      
        !Middle block (vertical constraints)
        do i2=1,n2
            i1=1
            ATx(i1+(i2-1)*n1)=&
            ATx(i1+(i2-1)*n1) &
          -h1*x(i1+(i2-1)*(n1-1)+dec1)
        do i1=2,n1-1
            ATx(i1+(i2-1)*n1)=&
            ATx(i1+(i2-1)*n1) &
          -h1*x(i1+(i2-1)*(n1-1)+dec1)&
          +h1*x(i1-1+(i2-1)*(n1-1)+dec1)
        enddo

        i1=n1
            ATx(i1+(i2-1)*n1)=&
            ATx(i1+(i2-1)*n1) &
          +h1*x(i1-1+(i2-1)*(n1-1)+dec1)
        enddo

        !Lower block (identity)
        do i=1,N
            ATx(i)=ATx(i)+x(i+dec1+dec2)
        enddo
      
    end subroutine sub_compute_ATx2_bis



    subroutine sub_compute_ATx3(NA,N,n1,n2,x,ATx)
    ! V0.1  - 04/2015  - L. Metivier
        integer,intent(in) :: NA,N,n1,n2
        real,dimension(NA),intent(in) :: x
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
      
    end subroutine

    subroutine sub_compute_ATAx(N,n1,ATA_sparse,x,ATAx)
    ! V0.0  - 03/2014  - L. Metivier
        integer,intent(in) :: N,n1
        real,dimension(N),intent(in) :: x
        real,dimension(3,N),intent(in) :: ATA_sparse
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

    end subroutine

    subroutine sub_compute_ATAx2(N,n1,ATA_sparse,x,ATAx)
    !BETTER IMPLEMENTATION
        integer,intent(in) :: N,n1
        real,dimension(N),intent(in) :: x
        real,dimension(N,3),intent(in) :: ATA_sparse
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
       
    end subroutine

end
