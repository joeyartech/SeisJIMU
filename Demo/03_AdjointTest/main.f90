program main
use m_System
use m_Modeling
use m_hilbert

    !mpiworld lives in t_mpienv
    call mpiworld%init(name='MPIWorld')

    call hud('======================================'//s_NL// &
             '       WELCOME TO SeisJIMU FWI        '//s_NL// &
             '======================================')
    
    call setup%init
    call sysio_init
    
    if(.not. setup%exist) then    
        call hud('No input setup file given. Stop.')
        call mpiworld%final
        stop
    endif

    !print propagator info
    call ppg%print_info

    !model
    call m%init
    call m%read
    call ppg%check_model

    !shotlist
    call shls%read_from_setup
    call shls%build(o_batchsize=shls%nshot)
    call shls%sample
    call shls%assign
    
    call modeling_gradient

    call mpiworld%final

    ! stop

end

subroutine modeling_gradient
use m_System
use m_Modeling
use m_hilbert

    double precision :: LHS=0., RHS=0.

    real,dimension(:,:),allocatable :: u,v,Lu,Lv !forward
    real,dimension(:,:),allocatable :: q,p,Ladj_q,Ladj_p, Ladj_p_hilb !adjoint
    type(t_field) :: sfield_u, sfield_v, rfield_p, rfield_q
    type(t_correlate) :: a_star_u

    logical :: if_use_random

!     call alloc(m%gradient,m%nz,m%nx,m%ny,ppg%ngrad)
            
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,1 !shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_setup
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project

        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer
        
        call ppg%init_field(sfield_u,name='sfield_u')
        call ppg%init_field(sfield_v,name='sfield_v')

        call ppg%init_correlate(a_star_u,'a_star_u')
        
        if_use_random=setup%get_bool('IF_USE_RANDOM',o_default='T')
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !variables for dotproduct test
        call alloc(u     ,ppg%nt,1        )
        call alloc(v     ,ppg%nt,1        )
        call alloc(q     ,ppg%nt,1        )
        call alloc(p     ,ppg%nt,1        )
        call alloc(Lu    ,ppg%nt,shot%nrcv)
        call alloc(Lv    ,ppg%nt,shot%nrcv)


        if(if_use_random) then
            call random_number(u)
        else
            u(:,1)=shot%wavelet
            ! v(:,1)=shot%wavelet
            call hilbert_transform(u,v,ppg%nt,1)
        endif
        call suformat_write('u',u,ppg%nt,1        ,o_dt=ppg%dt)
        call suformat_write('v',v,ppg%nt,1        ,o_dt=ppg%dt)
        ! call suformat_write('v',v,ppg%nt,shot%nrcv,o_dt=ppg%dt)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call sfield_u%ignite(o_wavelet=u)
        call sfield_v%ignite(o_wavelet=v)
        
        !forward modeling
        call ppg%forward(sfield_u,sfield_v)
        
        call sfield_u%acquire(o_seismo=Lu)
        call sfield_v%acquire(o_seismo=Lv)

        !call shot%write('dsyn_')
        call suformat_write('Lu',Lu,ppg%nt,shot%nrcv,o_dt=ppg%dt)
        call suformat_write('Lv',Lv,ppg%nt,shot%nrcv,o_dt=ppg%dt)

        call ppg%init_field(rfield_q,name='rfield_q',ois_adjoint=.true.)
        call ppg%init_field(rfield_p,name='rfield_p',ois_adjoint=.true.)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !variables for dotproduct test
        call alloc(q     ,ppg%nt,shot%nrcv)
        call alloc(p     ,ppg%nt,shot%nrcv)
        call alloc(Ladj_q,ppg%nt,1        )
        call alloc(Ladj_p,ppg%nt,1        )
        call alloc(Ladj_p_hilb,ppg%nt,1        )
        if(if_use_random) then
            call random_number(q)
        else
            q=Lu
            p=Lv
        endif
        call suformat_write('q',q,ppg%nt,shot%nrcv,o_dt=ppg%dt)
        call suformat_write('p',p,ppg%nt,shot%nrcv,o_dt=ppg%dt)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call rfield_q%ignite(o_wavelet=q)
        call rfield_p%ignite(o_wavelet=p)
        
        !adjoint modeling
        call ppg%adjoint(rfield_q,rfield_p,sfield_u,sfield_v,a_star_u)

        call rfield_q%acquire(o_seismo=Ladj_q)
        call rfield_p%acquire(o_seismo=Ladj_p)
        Ladj_p=-Ladj_p
        ! call hilbert_transform(Ladj_p,Ladj_p_hilb,ppg%nt,1)
        ! call hilbert_transform(Ladj_p_hilb,Ladj_p,,ppg%nt,1)
        ! call hilbert_transform(Ladj_p,Ladj_p_hilb,ppg%nt,1)
        call suformat_write('Ladj_q',Ladj_q,ppg%nt,1,o_dt=ppg%dt)
        call suformat_write('Ladj_p',Ladj_p,ppg%nt,1,o_dt=ppg%dt)
        
!         call cb%project_back
        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

!     call sysio_write('gradient',m%gradient,size(m%gradient))

    print*,'shape(u)=',     shape(u),    ', ║u║₂=',      norm2(u)*sqrt(ppg%dt)
    print*,'shape(v)=',     shape(v),    ', ║v║₂=',      norm2(v)*sqrt(ppg%dt)
    print*,'shape(q)=',     shape(q),    ', ║q║₂=',      norm2(q)*sqrt(ppg%dt)
    print*,'shape(p)=',     shape(p),    ', ║v║₂=',      norm2(p)*sqrt(ppg%dt)
    print*,'shape(Lu)=',    shape(Lu),   ', ║Lu║₂=',     norm2(Lu)*sqrt(ppg%dt)
    print*,'shape(Lv)=',    shape(Lv),   ', ║Lv║₂=',     norm2(Lv)*sqrt(ppg%dt)
    print*,'shape(Lᴴq)=',   shape(Ladj_q),', ║Lᴴq║₂=', norm2(Ladj_q)*sqrt(ppg%dt)
    print*,'shape(Lᴴp)=',   shape(Ladj_p),', ║Lᴴp║₂=', norm2(Ladj_p)*sqrt(ppg%dt)
    print*,'Remind that ║u║₂ := √ (∫ u² dt) = norm2(u)*sqrt(dt)'

    !<v|Lu> =?= <L^Tv|u>
    !<v|Lu>=int v*Lu*dt = sum(v*Lu)*dt
    ! LHS=sum(dprod(cmplx(q,-p)    ,cmplx(Lu,lv)))*ppg%dt
    ! RHS=sum(dprod(cmplx(Ladj_q,-Ladj_p),cmplx(u,v)))*ppg%dt

    ! LHS=sum(dprod(q,Lu))*ppg%dt
    ! RHS=sum(dprod(Ladj_q,u))*ppg%dt

    ! LHS=sum(dprod(p,Lv))*ppg%dt
    ! RHS=sum(dprod(Ladj_p,v))*ppg%dt

    
    LHS=sum(dprod(q,Lu))*ppg%dt + sum(dprod(p,Lv))*ppg%dt
    RHS=sum(dprod(Ladj_q,u))*ppg%dt + sum(dprod(Ladj_p,v))*ppg%dt

    print*,'LHS = <  v|Lu> = ', LHS
    print*,'RHS = <Lᴴv| u> = ', RHS
    print*,'relative difference = ', (LHS-RHS)/LHS


    LHS=sum(dprod(p,Lu))*ppg%dt 
    RHS=sum(dprod(Ladj_p,u))*ppg%dt 

    print*,'LHS = <  v|Lu> = ', LHS
    print*,'RHS = <Lᴴv| u> = ', RHS
    print*,'relative difference = ', (LHS-RHS)/LHS

    LHS=sum(dprod(q,Lv))*ppg%dt
    RHS=sum(dprod(Ladj_q,v))*ppg%dt

    print*,'LHS = <  v|Lu> = ', LHS
    print*,'RHS = <Lᴴv| u> = ', RHS
    print*,'relative difference = ', (LHS-RHS)/LHS

    LHS=sum(dprod(p,Lu))*ppg%dt - sum(dprod(q,Lv))*ppg%dt
    RHS=sum(dprod(Ladj_p,u))*ppg%dt - sum(dprod(Ladj_q,v))*ppg%dt


    print*,'LHS = <  v|Lu> = ', LHS
    print*,'RHS = <Lᴴv| u> = ', RHS
    print*,'relative difference = ', (LHS-RHS)/LHS

    call mpiworld%barrier

end subroutine
