program main
use m_System
use m_Modeling

    !mpiworld lives in t_mpienv
    call mpiworld%init(name='MPIWorld')

    call hud('======================================'//s_NL// &
             '       WELCOME TO SeisJIMU FWI        '//s_NL// &
             '======================================')
    
    call setup%init
    
    if(.not. setup%exist) then    
        call hud('No input setup file given. Stop.')
        call mpiworld%fin
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
    call shls%build
    call shls%sample
    call shls%assign
    
    call gradient_modeling

    call sysio_write('gradient',m%gradient,size(m%gradient))

    call mpiworld%fin
    
    stop

end

subroutine gradient_modeling
use m_System
use m_Modeling

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

        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !variables for dotproduct test
        call alloc(u     ,1,        ppg%nt)
        call alloc(v     ,shot%nrcv,ppg%nt)
        call alloc(Lu    ,shot%nrcv,ppg%nt)
        call alloc(Ladj_v,1        ,ppg%nt)

        if_use_random=setup%get_bool('IF_USE_RANDOM',o_default='T')

        if(if_use_random) then
            call random_number(u)
            call random_number(v)
        else
            u(1,:)=shot%wavelet
            do ir=1,shot%nrcv; v(ir,:)=shot%wavelet; enddo
        endif
        call suformat_write('u',u,size(u,1),size(u,2))
        call suformat_write('v',v,size(v,1),size(v,2))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call sfield%ignite(o_wavelet=u)
        
        !forward modeling
        call ppg%forward(sfield)
        
        ! call sfield%acquire
        Lu=sfield%seismo

        !call shot%write('dsyn_')
        call suformat_write('Lu',Lu,size(Lu,1),size(Lu,2))

        call ppg%init_field(rfield,name='rfield',origin='rcv')

        call rfield%ignite(o_wavelet=v,ois_adjoint=.true.)

        call alloc(cb%grad,cb%mz,cb%mx,cb%my,ppg%ngrad)

        !adjoint modeling
        call ppg%adjoint(rfield,oif_record_adjseismo=.true.,o_sf=sfield,o_grad=cb%grad)
        Ladj_v=rfield%seismo
        call suformat_write('Ladj_v',Ladj_v,size(Ladj_v,1),size(Ladj_v,2))
        
        call cb%project_back(m%gradient,cb%grad,ppg%ngrad)
        deallocate(cb%grad)
        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')
    
    print*,'shape(u),||u|| = ', shape(u), norm2(u)*sqrt(ppg%dt) !sqrt(int u^2*delta(x)*dx3*dt) = sqrt(u^2*dt) = norm2(u)*sqrt(dt)
    print*,'shape(v),||v|| = ', shape(v), norm2(v)*sqrt(ppg%dt)
    print*,'shape(Lu),    ||Lu|| = ',     shape(Lu),     norm2(Lu)*sqrt(ppg%dt)
    print*,'shape(Ladj_v),||Ladj_v|| = ', shape(Ladj_v), norm2(Ladj_v)*sqrt(ppg%dt)

    !<v|Lu> =?= <L^Tv|u>
    LHS=sum(dprod(v,Lu))*ppg%dt  !int v*Lu*delta(x)*dx3*dt = sum(v*Lu)*dt
    RHS=sum(dprod(Ladj_v,u))*ppg%dt
    print*,'LHS = <   v|Lu> = ', LHS
    print*,'RHS = <L^Tv| u> = ', RHS

    call mpiworld%barrier

end subroutine

