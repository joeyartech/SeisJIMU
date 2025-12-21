! ----------------------------------------------------------------------------------------------------
! PROGRAMM SMOOTHGAUSS2D: smooth a 2D grid with a Gaussian filter
!           smoothing is applied below topography
!
! itypsmo : 0 ::: smooth the parameter
!           1 ::: smooth the inverse of the parameter
!           2 ::: smooth the square of the inverse of the parameter 
! (good for bulk modulus which should be averaged using the harmonic mean .... and the bulk modulus
!  behaves as the square of the velocity and therefore its inverse as the square of the slowness
!
! ----------------------------------------------------------------------------------------------------
program smooth_kinematic
  IMPLICIT NONE
  CHARACTER(LEN=132) :: name_in,name_out,name_topo
  !------------------- arrays to be smoothed 
  REAL(kind=4),ALLOCATABLE :: x(:,:),x1(:,:),xcp(:,:)
  !------------------- arrays for topography manipulation
  REAL(kind=4), ALLOCATABLE :: topo(:)
  INTEGER(kind=4),ALLOCATABLE :: itopo(:)
  !------------------- arrays for gaussian distribution
  REAL(kind=4), ALLOCATABLE :: beta2(:),beta1(:)

  INTEGER(kind=4) :: n1,n2,itypsmo,i1,i2,il2,il1,k2,k1,j1,j2,jj1,jj2
  REAL(kind=4) :: t
  REAL(kind=4) :: d1,d2,sig1,sig2,xl2,xl1,sig12,sig22,d,betatot,tmp,betatot1
  ! -----------------------------------------------------------------------------------

  write(*,*) "Input / Output file names (name_in name_out)"
  read(*,*) name_in,name_out
  write(*,*) 'Grid dimensions (n1 n2)'
  read(*,*) n1,n2
  write(*,*) 'Grid interval (d1 d2)'
  read(*,*) d1,d2
  write(*,*) 'Standard deviation of the Gaussian function (sig1, sig2)'
  write(*,*) 'Wavenumbers (without factor 2*pi) higher than 0.48/sig will be suppressed'
  write(*,*) "c'est à dire que les nombres d'onde qui sont plus large que 0.48/sig ont l'amplitude sous 1% du crête. "
  read(*,*) sig1,sig2
  write(*,*) 'Smooth parameter (0) or its inverse (1) or its inverse square (2)'
  read(*,*) itypsmo
  write(*,*) 'Constant topography (>=0) or variant (<0)'
  read(*,*) t
  if (t < 0.) then
     write(*,*) ' enter the variant topography file name'
     read(*,'(a)') name_topo
  end if

  ! -----------------------------------------------------------------------------------
  open(1,file=name_in,access='direct',recl=n1*n2*4)
  open(2,file=name_out,access='direct',recl=n1*n2*4)
  ! ----------------------------------------------------------------
  ! DEBUG
  ! open(3,file='fdebug',access='direct',recl=n1*n2*4)
  ! ----------------------------------------------------------------
  if (t < 0.) open(10,file=name_topo,access='direct',recl=n2*4)
  ! -----------------------------------------------------------------------------------
  ALLOCATE (x(n1,n2),xcp(n1,n2))
  ALLOCATE (x1(n1,n2))

  ALLOCATE (topo(n2))
  ALLOCATE (itopo(n2))

  read(1,rec=1) ((x(i1,i2),i1=1,n1),i2=1,n2)
  xcp=x
  
  if (t < 0.) then
     read(10,rec=1) (topo(i2),i2=1,n2)
  else
     topo(:)=t
  end if

  do i2=1,n2
     itopo(i2)=nint(topo(i2)/d1)+1
     !write(*,*) i2,itopo(i2)
  end do

  if(itypsmo == 1) then ! we consider the inverse
    write(*,*) ' go to the inverse'
     do i2=1,n2
        do i1=itopo(i2),n1
           x(i1,i2)=1./x(i1,i2)
        end do
     end do
  elseif(itypsmo == 2) then ! we consider the square of the inverse
     write(*,*) ' go to the square inverse'
     do i2=1,n2
        do i1=itopo(i2),n1
           x(i1,i2)=1./x(i1,i2)/x(i1,i2)
        end do
     end do
  end if
  
  !extrapolate to water column
  do i2=1,n2
    do i1=1,itopo(i2)-1
        x(i1,i2)=x(itopo(i2),i2)
    end do
  end do  
    
!========================== we need to make the gaussian smoothing using two 1D operators

  xl2=4.*sig2     ! not necessary to go beyond
  xl1=4.*sig1     ! idem
  il2=int(xl2/d2)+1
  il1=int(xl1/d1)+1
  write(*,*) " grids points involved il2 il1 = ",il2,il1

  ALLOCATE (beta2(2*il2+1))
  ALLOCATE (beta1(2*il1+1))

  sig12=2*sig1*sig1
  sig22=2*sig2*sig2

  if(sig22>1e-8) then
    k2=0
    do i2=-il2,il2 ! reduction of the sampling
       k2=k2+1
       d=float(i2)*d2
       beta2(k2)=exp(-d*d/sig22)
    end do
  else
    print*,'sig2 is smaller than 7.1e-5. NO smoothing in the 2nd dim.'
  end if
  
  if(sig12>1e-8) then
    k1=0
    do i1=-il1,il1
       k1=k1+1
       d=float(i1)*d1
       beta1(k1)=exp(-d*d/sig12)
    end do
  else
    print*,'sig1 is smaller than 7.1e-5. NO smoothing in the 1st dim.'
  end if

  
  !==================================   1D loop
  if(sig22>1e-8) then
    do i2=1,n2
    do i1=itopo(i2),n1
       x1(i1,i2)=0.
       betatot=0.
       jj2=0
       do j2=-il2,il2
	  jj2=jj2+1
	  k2=j2+i2
	  if(k2 >= 1 .and. k2 <= n2) then
	    x1(i1,i2)=x1(i1,i2)+beta2(jj2)*x(i1,k2)
	    betatot=betatot+beta2(jj2)
	  endif
       end do
       x1(i1,i2)=x1(i1,i2)/betatot   ! get the weighted new value
    end do
    end do
  else
    x1=x
  end if

  ! -------------------------------------------------------------------------
  ! DEBUG
  !         if (itypsmo.eq.1) then
  !        do i2=1,n2
  !        do i1=itopo(i2),n1
  !        x1(i1,i2)=1./x1(i1,i2)
  !        end do
  !        end do
  !        end if
  !
  !         write(3,rec=1) ((x1(i1,i2),i1=1,n1),i2=1,n2)
  ! ------------------------------------------------------------------------
  !================================= The other loop

  if(sig12>1e-8) then
    do i2=1,n2
    do i1=itopo(i2),n1
       x(i1,i2)=0.
       betatot=0.
       jj1=0
       do j1=-il1,il1
	  jj1=jj1+1
	  k1=j1+i1
	  if(k1 >= itopo(i2) .and. k1 <= n1) then ! only values below the topography
	    x(i1,i2)=x(i1,i2)+beta1(jj1)*x1(k1,i2)
	    betatot=betatot+beta1(jj1)
	  endif
       end do
       x(i1,i2)=x(i1,i2)/betatot ! get the average value
    end do
    end do
  else
    x=x1
  end if

  if(itypsmo == 1) then
    write(*,*) ' go back from the inverse'
     do i2=1,n2
        do i1=1,n1
           x(i1,i2)=1./x(i1,i2)
        end do
     end do
  elseif(itypsmo == 2) then
    write(*,*) ' go back from the square inverse'
     do i2=1,n2
        do i1=1,n1
           x(i1,i2)=1./dsqrt(dble(x(i1,i2)))
        enddo
     enddo
  endif
  
  !put back water column
  do i2=1,n2
    do i1=1,itopo(i2)-1
        x(i1,i2)=xcp(i1,i2)
    end do
  end do  

  write(2,rec=1) ((x(i1,i2),i1=1,n1),i2=1,n2)
  print*,sum(sum(abs(x(:,:)),1))

  close(2)
  close(1)
  !============= close the topography in case
  if (t<0.) close(10)

  DEALLOCATE (x)
  DEALLOCATE (x1)

  DEALLOCATE (topo)
  DEALLOCATE (itopo)
  DEALLOCATE (beta2)
  DEALLOCATE (beta1)

  stop
end program smooth_kinematic



