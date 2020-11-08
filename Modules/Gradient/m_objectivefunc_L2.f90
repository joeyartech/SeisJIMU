module m_objectivefunc
use m_field, only:survey
use m_shot
use m_matchfilter
use m_weighter_polygon
use m_weighter_table

    real dnorm, mnorm !norm of data and model residuals

    real,dimension(:,:),allocatable :: weight

    contains

    subroutine objectivefunc_data_norm_residual

        real,save :: ref_modulus

        ref_modulus=m%ref_vp**2*m%ref_rho

        if(.not. allocated(weight)) then
            call alloc(weight,shot%rcv(1)%nt,shot%nrcv,initialize=.false.);  weight=1
!             call build_weight_polygon(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,weight) !so far the weighting is for mono component data only
!             call build_weight_table(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,weight) !so far the weighting is for mono component data only
! open(33,file='weight_'//shot%cindex,access='stream')
! write(33) weight
! close(33)
        endif

        if(survey=='base') then !survey is public in m_field_AC.f90
            dres  = (dsyn -dobs )*weight
        else
            dres2 = (dsyn2-dobs2)*weight
        endif
        
        dnorm= 0.
        
        !set unit of dnorm to be [Nm], same as Lagrangian
        !this also help balance contributions from different component data
        if(survey=='base') then
            do ir=1,shot%nrcv
                if(shot%rcv(ir)%icomp==1) then !for pressure data
                    dnorm = dnorm + sum(dres(:,ir)*dres(:,ir))*0.5/ref_modulus*m%cell_volume
                else !for velocities data
                    dnorm = dnorm + sum(dres(:,ir)*dres(:,ir))*0.5*m%ref_rho  *m%cell_volume
                endif
            enddo
        else
            do ir=1,shot2%nrcv
                if(shot2%rcv(ir)%icomp==1) then !for pressure data
                    dnorm = dnorm + sum(dres2(:,ir)*dres2(:,ir))*0.5/ref_modulus*m%cell_volume
                else !for velocities data
                    dnorm = dnorm + sum(dres2(:,ir)*dres2(:,ir))*0.5*m%ref_rho  *m%cell_volume
                endif
            enddo
        endif

!write balanced residuals
if(survey=='base') then
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
else
    if(mpiworld%is_master) then
    open(12,file='residu_unit_2_'//shot%cindex,access='stream')
    do ir=1,shot2%nrcv
      if(shot2%rcv(ir)%icomp==1) then !for pressure data
        write(12) dres2(:,ir)/sqrt(ref_modulus)
      else !for velocities data
        write(12) dres2(:,ir)*sqrt(m%ref_rho)
      endif
    enddo
    close(12)
    endif
endif
        
        !compute adjoint source and set proper units
        if(survey=='base') then
            dres = - dres*weight
            do ir=1,shot%nrcv
                if(shot%rcv(ir)%icomp==1) then !for pressure data
                    dres(:,ir) = dres(:,ir) / ref_modulus/shot%rcv(1)%dt
                else !for velocities data
                    dres(:,ir) = dres(:,ir) * m%ref_rho  /shot%rcv(1)%dt
                endif
            enddo
        else
            dres2 = - dres2*weight
            do ir=1,shot2%nrcv
                if(shot2%rcv(ir)%icomp==1) then !for pressure data
                    dres2(:,ir) = dres2(:,ir) / ref_modulus/shot2%rcv(1)%dt
                else !for velocities data
                    dres2(:,ir) = dres2(:,ir) * m%ref_rho  /shot2%rcv(1)%dt
                endif
            enddo
        endif

!write balanced adjoint source
if(survey=='base') then
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
else
    if(mpiworld%is_master) then
    open(12,file='adjsrc_unit_2_'//shot2%cindex,access='stream')
    do ir=1,shot2%nrcv
      if(shot2%rcv(ir)%icomp==1) then !for pressure data
        write(12) dres2(:,ir)/ref_modulus
      else !for velocities data
        write(12) dres2(:,ir)*m%ref_rho
      endif
    enddo
    close(12)
    endif
endif
        
        !multi-valued objectivefunc ..
        
    end subroutine
    
    subroutine objectivefunc_model_norm_residual
        
        mnorm=0. !to be developed...

    end subroutine

end
