module m_objectivefunc
use m_shot
use m_matchfilter
use m_weighting

    real dnorm, mnorm !norm of data and model residuals

    contains

    subroutine objectivefunc_data_norm_residual
        
        call build_weighting(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset) !so far the weighting is for mono component data only
        
        dres = (dsyn-dobs)*weight
        
        dnorm= 0.
        
        !set unit of dnorm to be [Nm], same as Lagrangian
        !this also help balance contributions from different component data
        do ir=1,shot%nrcv
            if(shot%rcv(ir)%icomp==1) then !for pressure data
                dnorm = dnorm + sum(dres(:,ir)*dres(:,ir))*0.5/m%ref_kpa*m%cell_size
            else !for velocities data
                dnorm = dnorm + sum(dres(:,ir)*dres(:,ir))*0.5*m%ref_rho*m%cell_size
            endif
        enddo
        
        !compute adjoint source and set proper units
        dres = - dres*weight
        do ir=1,shot%nrcv
            if(shot%rcv(ir)%icomp==1) then !for pressure data
                dres(:,ir) = dres(:,ir) / m%ref_kpa/shot%rcv(1)%dt
            else !for velocities data
                dres(:,ir) = dres(:,ir) * m%ref_rho/shot%rcv(1)%dt
            endif
        enddo
        
        !multi-valued objectivefunc ..
        
    end subroutine
    
    subroutine objectivefunc_model_norm_residual
        
        mnorm=0. !to be developed...

    end subroutine

end
