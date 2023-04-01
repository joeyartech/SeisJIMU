module m_objectivefunc
use m_shot
use m_matchfilter
use m_weighter_polygon
use m_weighter_table
use m_separater

    real dnorm, mnorm !norm of data and model residuals

    real,dimension(:,:),allocatable :: weight

    real,dimension(:,:),allocatable :: sepa_div, sepa_rfl, tmp_dres_rfl

    contains

    subroutine objectivefunc_data_norm_residual(o_residual)
        character(*),optional :: o_residual
        
        real,save :: ref_modulus

        ref_modulus=m%ref_vp**2*m%ref_rho

        ! if(.not. allocated(weight)) then
            call alloc(weight,shot%rcv(1)%nt,shot%nrcv,initialize=.false.);  weight=1
            call build_weight_polygon(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,weight) !so far the weighting is for mono component data only
            call build_weight_table(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,weight) !so far the weighting is for mono component data only
! open(33,file='weight',access='stream')
! write(33) weight
! close(33)
        ! endif

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


        if(present(o_residual)) then
            ! if(.not. allocated(sepa_div) .or. .not. allocated(sepa_rfl)) then
                call alloc(sepa_div,shot%rcv(1)%nt,shot%nrcv)
                call alloc(sepa_rfl,shot%rcv(1)%nt,shot%nrcv)
                call build_separater(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,sepa_div,sepa_rfl)
            ! endif

            if(o_residual=='rfl') then
                dres =       sepa_rfl*dres

                call alloc(tmp_dres_rfl,shot%rcv(1)%nt,shot%nrcv,initialize=.false.)
                tmp_dres_rfl=dres

            elseif(o_residual=='div-rfl') then
                dres =       sepa_div*dres
                dres = dres -sepa_rfl*tmp_dres_rfl

            endif

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

end
