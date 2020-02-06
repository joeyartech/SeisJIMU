module m_objectivefunc
use m_shot
use m_matchfilter
use m_weighter_polygon
use m_weighter_table

    real dnorm, mnorm !norm of data and model residuals

    real,dimension(:,:),allocatable :: weight

real,dimension(:,:),allocatable :: dsyn2, dobs2, tmp_dobs2

    contains

    subroutine objectivefunc_data_norm_residual

        real,save :: ref_modulus

        ref_modulus=m%ref_vp**2*m%ref_rho

        if(.not. allocated(weight)) then
            call alloc(weight,shot%rcv(1)%nt,shot%nrcv,initialize=.false.);  weight=1
            call build_weight_polygon(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,weight) !so far the weighting is for mono component data only
            call build_weight_table(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,weight) !so far the weighting is for mono component data only
open(33,file='weight_'//shot%cindex,access='stream')
write(33) weight
close(33)
        endif


if(.not. allocated(dobs2)) then
allocate(    dobs2(shot%rcv(1)%nt,    shot%nrcv))
allocate(tmp_dobs2(shot%rcv(1)%nt+60, shot%nrcv))
open(34,file=get_setup_char('FILE_DATA2_OBS')//shot%cindex//'.su',access='direct',recl=4*(shot%rcv(1)%nt+60)*shot%nrcv)
read(34,rec=1) tmp_dobs2
close(34)
dobs2=tmp_dobs2(61:shot%rcv(1)%nt+60,:)
deallocate(tmp_dobs2)
endif

if(.not. allocated(dsyn2)) then
allocate(    dsyn2(shot%rcv(1)%nt,     shot%nrcv))
open(35,file=get_setup_char('FILE_DATA2_SYN')//shot%cindex,access='direct',recl=4* shot%rcv(1)%nt    *shot%nrcv)
read(35,rec=1) dsyn2
close(35)
endif

        !dres = (dsyn-dobs)*weight
dres = ( (dsyn-dsyn2) - (dobs-dobs2) )*weight

!if(mpiworld%is_master) then
!open(12,file='dsyns',access='stream')
!write(12) dsyn
!write(12) dsyn2
!write(12) dsyn-dsyn2
!close(12)
!open(12,file='dobss',access='stream')
!write(12) dobs
!write(12) dobs2
!write(12) dobs-dobs2
!close(12)
!endif

        
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

!write balanced residuals
if(mpiworld%is_master) then
open(12,file='residu_unit_'//shot%cindex,access='stream')
do ir=1,shot%nrcv
  if(shot%rcv(ir)%icomp==1) then !for pressure data
    write(12) dres(:,ir)/sqrt(ref_modulus)
  else !for velocities data
    write(12) dres(:,ir)*sqrt(m%ref_rho)
  endif
enddo
close(12)
endif

        
        !compute adjoint source and set proper units
        dres = - dres*weight
        do ir=1,shot%nrcv
            if(shot%rcv(ir)%icomp==1) then !for pressure data
                dres(:,ir) = dres(:,ir) / ref_modulus/shot%rcv(1)%dt
            else !for velocities data
                dres(:,ir) = dres(:,ir) * m%ref_rho  /shot%rcv(1)%dt
            endif
        enddo

!write balanced adjoint source
if(mpiworld%is_master) then
open(12,file='adjsrc_unit_'//shot%cindex,access='stream')
do ir=1,shot%nrcv
  if(shot%rcv(ir)%icomp==1) then !for pressure data
    write(12) dres(:,ir)/ref_modulus
  else !for velocities data
    write(12) dres(:,ir)*m%ref_rho
  endif
enddo
close(12)
endif

        
        !multi-valued objectivefunc ..
        
    end subroutine
    
    subroutine objectivefunc_model_norm_residual
        
        mnorm=0. !to be developed...

    end subroutine

end
