program main
use m_System
use m_Modeling

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

    call sysio_write('gradient',m%gradient,size(m%gradient))

    call mpiworld%final

    ! stop

end

subroutine modeling_gradient
use m_System
use m_Modeling

    double precision :: LHS=0., RHS=0.

    real,dimension(:,:),allocatable :: u,v,Lu,Ladj_v
    type(t_field) :: sfield, rfield

    logical :: if_use_random

    call alloc(m%gradient,m%nz,m%nx,m%ny,ppg%ngrad)
            
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
        call ppg%init_field(sfield,name='sfield',origin='src',oif_will_reconstruct=.true.)
        call ppg%init_abslayer

        
        if_use_random=setup%get_bool('IF_USE_RANDOM',o_default='T')
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !variables for dotproduct test
        call alloc(u     ,ppg%nt,1        )
        call alloc(Lu    ,ppg%nt,shot%nrcv)

        if(if_use_random) then
            call random_number(u)
        else
            u(:,1)=shot%wavelet
        endif
        call suformat_write('u',u,ppg%nt,1        ,o_dt=ppg%dt)
        ! call suformat_write('v',v,ppg%nt,shot%nrcv,o_dt=ppg%dt)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call sfield%ignite(o_wavelet=u)
        
        !forward modeling
        call ppg%forward(sfield)
        
        call sfield%acquire(o_seismo=Lu)

        !call shot%write('dsyn_')
        call suformat_write('Lu',Lu,ppg%nt,shot%nrcv,o_dt=ppg%dt)

        call ppg%init_field(rfield,name='rfield',origin='rcv')

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !variables for dotproduct test
        call alloc(v     ,ppg%nt,shot%nrcv)
        call alloc(Ladj_v,ppg%nt,1        )
        if(if_use_random) then
            call random_number(v)
        else
            v=Lu
        endif
        call suformat_write('v',v,ppg%nt,shot%nrcv,o_dt=ppg%dt)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call rfield%ignite(o_wavelet=v,ois_adjoint=.true.)

        !adjoint modeling
        call ppg%adjoint(rfield,oif_record_adjseismo=.true.,o_sf=sfield,o_grad=cb%grad)

        call rfield%acquire(o_seismo=Ladj_v)
        call suformat_write('Ladj_v',Ladj_v,ppg%nt,1,o_dt=ppg%dt)
        
        call cb%project_back(m%gradient,cb%grad,ppg%ngrad)
        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    print*,'shape(u)=',     shape(u),    '||u||=',      norm2(u)*sqrt(ppg%dt)
    !||u||=sqrt(int u^2*dt) = sqrt(sum(u^2)*dt) = norm2(u)*sqrt(dt)
    print*,'shape(v)=',     shape(v),    '||v||=',      norm2(v)*sqrt(ppg%dt)
    print*,'shape(Lu)=',    shape(Lu),   '||Lu||=',     norm2(Lu)*sqrt(ppg%dt)
    print*,'shape(Ladj_v)=',shape(Ladj_v),'||Ladj_v||=',norm2(Ladj_v)*sqrt(ppg%dt)

    !<v|Lu> =?= <L^Tv|u>
    !<v|Lu>=int v*Lu*dt = sum(v*Lu)*dt
    LHS=sum(dprod(v    ,Lu))*ppg%dt
    RHS=sum(dprod(Ladj_v,u))*ppg%dt

    print*,'LHS = <   v|Lu> = ', LHS
    print*,'RHS = <L^Tv| u> = ', RHS
    print*,'relative difference = ', (LHS-RHS)/LHS

    call mpiworld%barrier

end subroutine
