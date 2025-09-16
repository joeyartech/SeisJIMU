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

    double precision :: LHS1, LHS2, LHS3
    double precision :: RHS1, RHS2, RHS3

    real,dimension(:,:),allocatable :: reS,imS, reR,imR, reLS,imLS, reLadjR, imLadjR
    type(t_field) :: reU, imU, reA, imA
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
        
        call ppg%init_field(reU,name='reU')
        call ppg%init_field(imU,name='imU')

        call ppg%init_correlate(a_star_u,'a_star_u')
        
        if_use_random=setup%get_bool('IF_USE_RANDOM',o_default='T')
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !variables for dotproduct test
        call alloc(reS     ,ppg%nt,1        )
        call alloc(imS     ,ppg%nt,1        )
        call alloc(reR     ,ppg%nt,shot%nrcv)
        call alloc(imR     ,ppg%nt,shot%nrcv)
        call alloc(reLS    ,ppg%nt,shot%nrcv)
        call alloc(imLS    ,ppg%nt,shot%nrcv)
        call alloc(reLadjR ,ppg%nt,1        )
        call alloc(imLadjR ,ppg%nt,1        )


        if(if_use_random) then
            call random_number(u)
        else
            reS(:,1)=shot%wavelet
            call hilbert_transform(reS,imS,ppg%nt,1)
        endif
        call suformat_write('reS',reS,ppg%nt,1        ,o_dt=ppg%dt)
        call suformat_write('imS',imS,ppg%nt,1        ,o_dt=ppg%dt)
        ! call suformat_write('v',v,ppg%nt,shot%nrcv,o_dt=ppg%dt)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call reU%ignite(o_wavelet=reS)
        call imU%ignite(o_wavelet=imS)
        
        !forward modeling
        call ppg%forward(reU,imU)
        
        call reU%acquire(o_seismo=reLS)
        call imU%acquire(o_seismo=imLS)

        !call shot%write('dsyn_')
        call suformat_write('reLS',reLS,ppg%nt,shot%nrcv,o_dt=ppg%dt)
        call suformat_write('imLS',imLS,ppg%nt,shot%nrcv,o_dt=ppg%dt)


        call ppg%init_field(reA,name='reA',ois_adjoint=.true.)
        call ppg%init_field(imA,name='imA',ois_adjoint=.true.)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(if_use_random) then
            call random_number(reR)
        else
            !reR(:,1)=shot%wavelet
            reR=reLS
!             imR=-imLS
            call hilbert_transform(reR,imR,ppg%nt,shot%nrcv)
        endif
        call suformat_write('reR',reR,ppg%nt,shot%nrcv,o_dt=ppg%dt)
        call suformat_write('imR',imR,ppg%nt,shot%nrcv,o_dt=ppg%dt)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call reA%ignite(o_wavelet=reR)
        call imA%ignite(o_wavelet=imR)
        
        !adjoint modeling
        call ppg%adjoint(reA,imA,reU,imU,a_star_u)

        call reA%acquire(o_seismo=reLadjR)
        call imA%acquire(o_seismo=imLadjR)
!         Ladj_p=-Ladj_p
        call suformat_write('reLadjR',reLadjR,ppg%nt,1,o_dt=ppg%dt)
        call suformat_write('imLadjR',imLadjR,ppg%nt,1,o_dt=ppg%dt)
        
!         call cb%project_back
        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

!     call sysio_write('gradient',m%gradient,size(m%gradient))

    print*,'Vector   |       shape       |   ║*║₂'
    print*,'reS  '   ,     shape(reS),     norm2(reS)*sqrt(ppg%dt)
    print*,'imS  '   ,     shape(imS),     norm2(imS)*sqrt(ppg%dt)
    print*,'reLS '   ,     shape(reLS),    norm2(reLS)*sqrt(ppg%dt)
    print*,'imLS '   ,     shape(imLS),    norm2(imLS)*sqrt(ppg%dt)
    print*,'reR  '   ,     shape(reR),     norm2(reR)*sqrt(ppg%dt)
    print*,'imR  '   ,     shape(imR),     norm2(imR)*sqrt(ppg%dt)
    print*,'reLᴴR'   ,     shape(reLadjR), norm2(reLadjR)*sqrt(ppg%dt)
    print*,'imLᴴR'   ,     shape(imLadjR), norm2(imLadjR)*sqrt(ppg%dt)
!     print*,'Remind that ║u║₂ := √ (∫ u² dt) = norm2(u)*sqrt(dt)'
    print*,''
    print*,'<R|LS> =?= <LᴴR|S>'
    print*,'LHS= (re  R-iim  R)·(reLS+iimLS) = re  R·reLS + im  R·imLS +i(re  R·imLS - im  R·reLS)'
    print*,'RHS= (reLᴴR-iimLᴴR)·(re S+iim S) = reLᴴR·re S + imLᴴR·im S +i(reLᴴR·im S - imLᴴR·re S)'
    print*,''

    print*,'Real parts:'
    LHS1=sum(dprod(    reR,reLS)); LHS2=sum(dprod(    imR,imLS)); LHS3=LHS1+LHS2
    RHS1=sum(dprod(reLadjR, reS)); RHS2=sum(dprod(imLadjR, imS)); RHS3=RHS1+RHS2
    print*,'LHS=', LHS1, LHS2, LHS3
    print*,'RHS=', RHS1, RHS2, RHS3
    print*,'Relative diff=', (LHS1-RHS1)/LHS1, (LHS2-RHS2)/LHS2, (LHS3-RHS3)/LHS3
    print*,''

    print*,'Imag parts:'
    LHS1=sum(dprod(    reR,imLS)); LHS2=sum(dprod(    imR,reLS)); LHS3=LHS1-LHS2
    RHS1=sum(dprod(reLadjR, imS)); RHS2=sum(dprod(imLadjR, reS)); RHS3=RHS1-RHS2
    print*,'LHS=', LHS1, LHS2, LHS3
    print*,'RHS=', RHS1, RHS2, RHS3
    print*,'Relative diff=', (LHS1-RHS1)/LHS1, (LHS2-RHS2)/LHS2, (LHS3-RHS3)/LHS3


    !<v|Lu> =?= <L^Tv|u>
    !<v|Lu>=int v*Lu*dt = sum(v*Lu)*dt
    ! LHS=sum(dprod(cmplx(q,-p)    ,cmplx(Lu,lv)))*ppg%dt
    ! RHS=sum(dprod(cmplx(Ladj_q,-Ladj_p),cmplx(u,v)))*ppg%dt

    ! LHS=sum(dprod(q,Lu))*ppg%dt
    ! RHS=sum(dprod(Ladj_q,u))*ppg%dt

    ! LHS=sum(dprod(p,Lv))*ppg%dt
    ! RHS=sum(dprod(Ladj_p,v))*ppg%dt
!
!
!     LHS=sum(dprod(q,Lu))*ppg%dt + sum(dprod(p,Lv))*ppg%dt
!     RHS=sum(dprod(Ladj_q,u))*ppg%dt + sum(dprod(Ladj_p,v))*ppg%dt
!
!     print*,'LHS = <  v|Lu> = ', LHS
!     print*,'RHS = <Lᴴv| u> = ', RHS
!     print*,'relative difference = ', (LHS-RHS)/LHS
!
!
!     LHS=sum(dprod(p,Lu))*ppg%dt
!     RHS=sum(dprod(Ladj_p,u))*ppg%dt
!
!     print*,'LHS = <  v|Lu> = ', LHS
!     print*,'RHS = <Lᴴv| u> = ', RHS
!     print*,'relative difference = ', (LHS-RHS)/LHS
!
!     LHS=sum(dprod(q,Lv))*ppg%dt
!     RHS=sum(dprod(Ladj_q,v))*ppg%dt
!
!     print*,'LHS = <  v|Lu> = ', LHS
!     print*,'RHS = <Lᴴv| u> = ', RHS
!     print*,'relative difference = ', (LHS-RHS)/LHS
!
!     LHS=sum(dprod(p,Lu))*ppg%dt - sum(dprod(q,Lv))*ppg%dt
!     RHS=sum(dprod(Ladj_p,u))*ppg%dt - sum(dprod(Ladj_q,v))*ppg%dt
!
!
!     print*,'LHS = <  v|Lu> = ', LHS
!     print*,'RHS = <Lᴴv| u> = ', RHS
!     print*,'relative difference = ', (LHS-RHS)/LHS

    call mpiworld%barrier

end subroutine
