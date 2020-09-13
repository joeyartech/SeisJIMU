module m_objectivefunc
use m_shot
use m_matchfilter
use m_weighter_polygon
use m_weighter_table

    real dnorm, mnorm !norm of data and model residuals

    real,dimension(:,:),allocatable :: weight

    contains

    subroutine objectivefunc_data_norm_residual
        
        real,save :: ref_modulus

character(:),allocatable :: tmp_char
real :: multi_t, multi_x

        ref_modulus=m%ref_vp**2*m%ref_rho

!        if(.not. allocated(weight)) then
            call alloc(weight,shot%rcv(1)%nt,shot%nrcv,initialize=.false.);  weight=1

           call build_weight_polygon(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,weight) !so far the weighting is for mono component data only
            call build_weight_table(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,weight) !so far the weighting is for mono component data only

            
tmp_char=get_setup_char('WEIGHT_MINMAXFILT',default='1e-3   2  1')
read(tmp_char,*) thres, multi_t, multi_x

!threshold weighting
tmp=maxval(abs(dobs)); where (abs(dobs)<tmp*thres) weight=0.

!minmax filt weighting
call minmaxfilt(weight,shot%rcv(1)%nt,shot%nrcv, &
nint(multi_t/shot%src%fpeak/shot%rcv(1)%dt),&
nint(multi_x*1500./shot%src%fpeak/12.5)) !ref_vp=1500, dx=12.5

if(mpiworld%is_master) print*,'minmaxfilt lent lenx =',&
nint(multi_t/shot%src%fpeak/shot%rcv(1)%dt),&
nint(multi_x*1500./shot%src%fpeak/12.5)

!mute bad traces
do ir=1,shot%nrcv
if(shot%rcv(ir)%icomp==2) then
weight(:,ir)=0.
endif
enddo
open(33,file='weight_'//shot%cindex,access='stream')
write(33) weight
close(33)
            
!mute bad traces
do ir=1,shot%nrcv
if(shot%rcv(ir)%icomp==2) then
dsyn(:,ir)=0.; dobs(:,ir)=0.
endif
enddo

!        endif


        dres = (dsyn-dobs)*weight
        
        dnorm= 0.
        
        !set unit of dnorm to be [Nm], same as Lagrangian
        !this also help balance contributions from different component data
        do ir=1,shot%nrcv
            if(shot%rcv(ir)%icomp==1) then !for pressure data
                dnorm = dnorm + sum(dres(:,ir)*dres(:,ir))*0.5/ref_modulus*m%cell_volume
            else !for velocities data
                dnorm = dnorm + sum(dres(:,ir)*dres(:,ir))*0.5*m%ref_rho  *m%cell_volume
            endif
        enddo

! !write balanced residuals
! if(mpiworld%is_master) then
! open(12,file='residu_unit_'//shot%cindex,access='stream')
! do ir=1,shot%nrcv
!   if(shot%rcv(ir)%icomp==1) then !for pressure data
!     write(12) dres(:,ir)/sqrt(ref_modulus)
!   else !for velocities data
!     write(12) dres(:,ir)*sqrt(m%ref_rho)
!   endif
! enddo
! close(12)
! endif

        
        !compute adjoint source and set proper units
        dres = - dres*weight
        do ir=1,shot%nrcv
            if(shot%rcv(ir)%icomp==1) then !for pressure data
                dres(:,ir) = dres(:,ir) / ref_modulus/shot%rcv(1)%dt
            else !for velocities data
                dres(:,ir) = dres(:,ir) * m%ref_rho  /shot%rcv(1)%dt
            endif
        enddo

! !write balanced adjoint source
! if(mpiworld%is_master) then
! open(12,file='adjsrc_unit_'//shot%cindex,access='stream')
! do ir=1,shot%nrcv
!   if(shot%rcv(ir)%icomp==1) then !for pressure data
!     write(12) dres(:,ir)/ref_modulus
!   else !for velocities data
!     write(12) dres(:,ir)*m%ref_rho
!   endif
! enddo
! close(12)
! endif

        
        !multi-valued objectivefunc ..
        
    end subroutine
    
    subroutine objectivefunc_model_norm_residual
        
        mnorm=0. !to be developed...

    end subroutine

    
subroutine minmaxfilt(weight,nt,nx,lent,lenx)
real,dimension(nt,nx) :: weight
real,dimension(:,:),allocatable :: copy
real,dimension(:),allocatable :: tmp
integer :: h

!lent, lenx must be odd integer
if(mod(lent,2)==0) lent=lent+1
if(lent==1) lent=3

if(mod(lenx,2)==0) lenx=lenx+1
if(lenx==1) lenx=3

!time dir
h=(lent-1)/2
allocate(copy(1-h:nt+h,nx))
copy(1:nt,:)=weight
do ix=1,nx
copy(1 -h:1,ix)=copy(1 ,ix)
copy(nt+h:1,ix)=copy(nt,ix)
enddo
allocate(tmp(lent))

do ix=1,nx
do it=1,nt
    tmp=copy(it-h:it+h,ix)
    n_ones=sum(tmp)
    n_zeros= lent-n_ones
    
    if(n_ones>n_zeros) then
        weight(it,ix) = 1.
    else
        weight(it,ix) = 0.
    endif
enddo
enddo

deallocate(copy,tmp)


!space dir
h=(lenx-1)/2
allocate(copy(nt,1-h:nx+h))
copy(:,1:nx)=weight
do it=1,nt
copy(it,1 -h:1)=copy(it,1 )
copy(it,nx+h:1)=copy(it,nx)
enddo
allocate(tmp(lenx))

do it=1,nt
do ix=1,nx
    tmp=copy(it,ix-h:ix+h)
    n_ones=sum(tmp)
    n_zeros= lenx-n_ones
    
    if(n_ones>n_zeros) then
        weight(it,ix) = 1.
    else
        weight(it,ix) = 0.
    endif
enddo
enddo

deallocate(copy,tmp)

end subroutine

end
