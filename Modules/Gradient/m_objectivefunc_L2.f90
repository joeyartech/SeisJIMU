module m_objectivefunc
use m_shot
use m_matchfilter
use m_weighting

    real dnorm, mnorm !norm of data and model residuals

    contains

    subroutine objectivefunc_data_norm_residual
        
        real,save :: ref_modulus

        ref_modulus=m%ref_vp**2*m%ref_rho

        call build_weighting(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset) !so far the weighting is for mono component data only
        
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
