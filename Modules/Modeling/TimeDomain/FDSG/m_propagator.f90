module m_propagator
use m_sysio
use m_arrayop
use m_shot
use m_field
use m_boundarystore

    private
    public sbloom,rbloom
    public init_propagator,propagator_forward,propagator_adjoint
    
    integer :: it      !index of time instant
    integer :: ift,ilt !first & last time instant
    integer :: nt
    real dt
    
    real,dimension(:),allocatable :: wavelet
    
    type(t_field) :: sfield, rfield
    
    !blooming
    integer,dimension(:,:),allocatable :: sbloom, rbloom
    
    !snapshot
    logical :: if_snapshot
    integer :: it_delta_snapshot
    
    contains
        
    subroutine init_propagator(if_will_do_rfield)
        logical,optional :: if_will_do_rfield
        
        integer,parameter :: initial_half_bloom=5 !half bloom width at ift, should involve hicks points
        
        !time window
        ift=1
        ilt=shot%src%nt
        nt=ilt-ift+1
        dt=shot%src%dt
        
        !scaling source wavelet
        call alloc(wavelet,nt)
        wavelet = shot%src%wavelet !*dt/m%cell_size
        
        !field
        call init_field_localmodel
        call init_field(sfield,if_save_previous=.true.) !source field, allocate prev arrays no matter of use them or not
        if(present(if_will_do_rfield)) then
        if(if_will_do_rfield) then
            call init_field(rfield) !receiver field
        endif
        endif
         
        !blooming
        call alloc(sbloom,6,nt) !sfield blooming
        call alloc(rbloom,6,nt) !rfield blooming
        
        !directional maximum propagation distance per time step
        distz = cb%velmax * dt / m%dz
        distx = cb%velmax * dt / m%dx
        disty = cb%velmax * dt / m%dy
        
        sbloom(1,1)=shot%src%iz -initial_half_bloom
        sbloom(2,1)=shot%src%iz +initial_half_bloom
        sbloom(3,1)=shot%src%ix -initial_half_bloom
        sbloom(4,1)=shot%src%ix +initial_half_bloom
        sbloom(5,1)=shot%src%iy -initial_half_bloom
        sbloom(6,1)=shot%src%iy +initial_half_bloom
        do it=2,nt
            sbloom(1,it)=max(nint(sbloom(1,1)-it*distz),cb%ifz) !bloombox ifz
            sbloom(2,it)=min(nint(sbloom(2,1)+it*distz),cb%ilz) !bloombox ilz
            sbloom(3,it)=max(nint(sbloom(3,1)-it*distx),cb%ifx) !bloombox ifx
            sbloom(4,it)=min(nint(sbloom(4,1)+it*distx),cb%ilx) !bloombox ilx
            sbloom(5,it)=max(nint(sbloom(5,1)-it*disty),cb%ify) !bloombox ify
            sbloom(6,it)=min(nint(sbloom(6,1)+it*disty),cb%ily) !bloombox ily
        enddo
        
        if(.not.m%is_cubic) sbloom(5:6,:)=1
        
        if(present(if_will_do_rfield)) then
        if(if_will_do_rfield) then
            rbloom(1,nt)=minval(shot%rcv(:)%iz) -initial_half_bloom
            rbloom(2,nt)=maxval(shot%rcv(:)%iz) +initial_half_bloom
            rbloom(3,nt)=minval(shot%rcv(:)%ix) -initial_half_bloom
            rbloom(4,nt)=maxval(shot%rcv(:)%ix) +initial_half_bloom
            rbloom(5,nt)=minval(shot%rcv(:)%iy) -initial_half_bloom
            rbloom(6,nt)=maxval(shot%rcv(:)%iy) +initial_half_bloom
            do it=nt-1,1,-1
                it_fwd=nt-it+1
                rbloom(1,it)=max(nint(rbloom(1,nt)-it_fwd*distz),cb%ifz) !bloombox ifz
                rbloom(2,it)=min(nint(rbloom(2,nt)+it_fwd*distz),cb%ilz) !bloombox ilz
                rbloom(3,it)=max(nint(rbloom(3,nt)-it_fwd*distx),cb%ifx) !bloombox ifx
                rbloom(4,it)=min(nint(rbloom(4,nt)+it_fwd*distx),cb%ilx) !bloombox ilx
                rbloom(5,it)=max(nint(rbloom(5,nt)-it_fwd*disty),cb%ify) !bloombox ify
                rbloom(6,it)=min(nint(rbloom(6,nt)+it_fwd*disty),cb%ily) !bloombox ily
            enddo
            
            if(.not.m%is_cubic) rbloom(5:6,:)=1
        endif
        endif
        
        !boundarystore
        call init_boundarystore
        
        !snapshot
        if_snapshot=get_setup_logical('IF_SNAPSHOT',default=.false.)
        if_snapshot=if_snapshot.and.mpiworld%is_master
        if(if_snapshot) it_delta_snapshot=get_setup_int('SNAPSHOT_DELTA_IT',default=50)
        
    end subroutine
    
    subroutine propagator_forward(if_will_backpropagate)
        logical,optional :: if_will_backpropagate
        
        integer,parameter :: time_dir=1 !time direction
        
        real,dimension(:,:),allocatable :: seismo
        
        call alloc(seismo,shot%nrcv,nt,initialize=.false.)
        
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.

        if(if_snapshot) open(16,file='snapshot_sfield%vz',access='stream')
        
        do it=ift,ilt
            if(mod(it,100)==0 .and. mpiworld%is_master) write(*,*) 'it----',it
            if(mod(it,100)==0) call check_field(sfield,'sfield')
            
            call cpu_time(t0)
            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            !step 1: add forces to v^it
            call put_velocities(time_dir,it,wavelet(it),sfield)
            call cpu_time(t1)
            tt1=tt1+t1-t0
            
            !step 2: from v^it to v^it+1 by differences of s^it+0.5
            call update_velocities(time_dir,it,sfield,sbloom(:,it))
            call cpu_time(t2)
            tt2=tt2+t2-t1
            
            !step 3: add pressure to s^it+0.5
            call put_stresses(time_dir,it,wavelet(it),sfield)
            call cpu_time(t3)
            tt3=tt3+t3-t2
            
            !step 4: from s^it+0.5 to s^it+1.5 by differences of v^it+1
            call update_stresses(time_dir,it,sfield,sbloom(:,it))
            call cpu_time(t4)
            tt4=tt4+t4-t3
            
            !step 5: sample v^it+1 or s^it+1.5 at receivers
            call get_field(sfield,seismo(:,it))
            
            !snapshot
            if(if_snapshot) then
            if(mod(it-1,it_delta_snapshot)==0 .or. it==ilt) then
                call write_field(16,sfield)
            endif
            endif
            
            !step 6: save v^it+1 in boundary layers
            if(present(if_will_backpropagate)) then
            if(if_will_backpropagate) then
                call cpu_time(t5)
                tt5=tt5+t5-t4
                call boundarystore_transport('save',it,sfield)
                call cpu_time(t6)
                tt6=tt6+t6-t5
            endif
            endif
            
        enddo
        
        if(if_snapshot) then
            close(16)
            write(*,'(a,i0.5,a)') 'ximage < snapshot_sfield%vz  n1=',cb%nz,' perc=99'
            write(*,'(a,i0.5,a,i0.5,a)') 'xmovie < snapshot_sfield%vz  n1=',cb%nz,' n2=',cb%nx,' clip=?e-?? loop=2 title=%g'
        endif
        
        if(mpiworld%is_master) then
            write(*,*) 'time add source velocities',tt1
            write(*,*) 'time add source velocities',tt1
            write(*,*) 'time update velocities    ',tt2
            write(*,*) 'time add source stresses  ',tt3
            write(*,*) 'time update stresses      ',tt4
            write(*,*) 'time extract&write fields ',tt5
            write(*,*) 'time save boundary        ',tt6
        endif
        
        !synthetic data
        call alloc(dsyn,shot%rcv(1)%nt,shot%nrcv)
        dsyn=transpose(seismo)
        
        deallocate(seismo)

    end subroutine
    
    subroutine propagator_adjoint(dout,gradient)
        real,dimension(nt),optional :: dout
        real,dimension(cb%mz,cb%mx,cb%my,2),optional :: gradient
        
        real,dimension(:,:),allocatable :: seismo
        
        integer,parameter :: time_dir=-1 !time direction
        
        
        !reinitialize memory for incident wavefield reconstruction
        call field_cpml_reinitialize(sfield)
        
        !save previous time-step wavefield
        if(present(gradient)) then
            call field_save_previous(sfield,sbloom(:,it))
        endif
        
        call alloc(seismo,shot%nrcv,nt,initialize=.false.)
        seismo=transpose(dres) !to save mem space, seismo could be just time slices..
        
        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.

        if(if_snapshot) then
            open(20,file='snapshot_sfield%vz_back',access='stream')
!             open(22,file='snapshot_sfield%vz_deri',access='stream')
            open(24,file='snapshot_rfield%vz',access='stream')
            open(26,file='snapshot_scorr',access='stream')
            open(28,file='snapshot_vcorr',access='stream')
        endif
        
        
        do it=ilt,ift,time_dir
            
            if(mod(it,100)==0 .and. mpiworld%is_master) write(*,*) 'it----',it
            if(mod(it,100)==0) call check_field(sfield,'sfield')
            if(mod(it,100)==0) call check_field(rfield,'rfield')
            
            call cpu_time(t0)
            !do backward time stepping to reconstruct the source (incident) wavefield
            !step# conforms with forward time stepping
            !step 6: retrieve v^it+1 at boundary layers (BC)
            call boundarystore_transport('load',it,sfield)
            call cpu_time(t1)
            tt1=tt1+t1-t0
            
            !step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
            call update_stresses(time_dir,it,sfield,sbloom(:,it))
            call cpu_time(t2)
            tt2=tt2+t2-t1
            
            !step 3: rm pressure from s^it+0.5
            call put_stresses(time_dir,it,wavelet(it),sfield)
            call cpu_time(t3)
            tt3=tt3+t3-t2
            
            !step 2: v^it+1 -> v^it by FD of s^it+0.5
            call update_velocities(time_dir,it,sfield,sbloom(:,it))
            call cpu_time(t4)
            tt4=tt4+t4-t3
            
            !step 1: rm forces from v^it
            call put_velocities(time_dir,it,wavelet(it),sfield)
            call cpu_time(t5)
            tt5=tt5+t5-t4

            !do adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping
            !step 5: fill s^it+1.5 at receivers
            call put_stresses_adjoint(time_dir,it,seismo(:,it),rfield)
            call cpu_time(t6)
            tt6=tt6+t6-t5
            
            !step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call update_stresses_adjoint(time_dir,it,rfield,rbloom(:,it))
            call cpu_time(t7)
            tt7=tt7+t7-t6
            
            !step 5: fill v^it+1 at receivers
            call put_velocities_adjoint(time_dir,it,seismo(:,it),rfield)
            call cpu_time(t8)
            tt8=tt8+t8-t7
            
            !step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call update_velocities_adjoint(time_dir,it,rfield,rbloom(:,it))
            call cpu_time(t9)
            tt9=tt9+t9-t8
            
            !step 1: sample v^it or s^it+0.5 at source position
            if(present(dout)) call get_field_adjoint(rfield,dout(it))
            
            !snapshot
            if(if_snapshot) then
            if(mod(it-1,it_delta_snapshot)==0 .or. it==ift) then
                call write_field(20,sfield)
!                 call write_field(22,sfield,if_difference=.true.) !s^it+1.5 - s^it+0.5
                call write_field(24,rfield)
                write(26)gradient(:,:,:,1)
                write(28)gradient(:,:,:,2)
            endif
            endif
            
            call cpu_time(t10)
            tt10=tt10+t10-t9

            !crosscorrelation between
            !  sfield%s^it+0.5, v^it, prev_s^it+1.5, prev_v^it+1, and
            !  rfield%s^it+0.5, v^it, prev_s^it+1.5, prev_v^it+1
            if(present(gradient)) then
                call field_correlation_stresses(it,sfield,rfield,sbloom(:,it),rbloom(:,it),gradient(:,:,:,1))
                call field_correlation_velocities(it,sfield,rfield,sbloom(:,it),rbloom(:,it),gradient(:,:,:,2))
            endif
            
            call cpu_time(t11)
            tt11=tt11+t11-t10
            
            !save previous time-step wavefield
            !!actually we don't need to do so because the updates are available during timestepping, thus save mem space
            !!however this would lead to updating sfield and rfield simultaneously which makes the code harder to manage..
            if(present(gradient)) then
                call field_save_previous(sfield,sbloom(:,it))
            endif
            
            call cpu_time(t12)
            tt12=tt12+t12-t11
            
        enddo
        
        !scale gradient
        if(present(gradient)) then
            call field_correlation_scaling(gradient)
        endif
        
        if(if_snapshot) then
            close(20)
!             close(22)
            close(24)
            close(26)
            close(28)
            write(*,'(a,i0.5,a)') 'ximage < snapshot_rfield%vz  n1=',cb%nz,' perc=99'
            write(*,'(a,i0.5,a,i0.5,a)') 'xmovie < snapshot_rfield%vz  n1=',cb%nz,' n2=',cb%nx,' clip=?e-?? loop=2 title=%g'
            
            write(*,'(a,i0.5,a,i0.5,a)') 'xmovie < snapshot_scorr  n1=',cb%mz,' n2=',cb%mx,' clip=?e-?? loop=2 title=%g'
            write(*,'(a,i0.5,a,i0.5,a)') 'xmovie < snapshot_scorr  n1=',cb%mz,' n2=',cb%mx,' clip=?e-?? loop=2 title=%g'
        endif
        
        
        if(mpiworld%is_master) then
            write(*,*) 'time load boundary            ',tt1
            write(*,*) 'time update stresses          ',tt2
            write(*,*) 'time rm source stresses       ',tt3
            write(*,*) 'time update velocities        ',tt4
            write(*,*) 'time rm source velocities     ',tt5
            write(*,*) 'time -------------------------'
            write(*,*) 'time add adjsource stresses   ',tt6
            write(*,*) 'time update adj stresses      ',tt7
            write(*,*) 'time add adjsource velocities ',tt8
            write(*,*) 'time update adj velocities    ',tt9
            write(*,*) 'time extract&write fields     ',tt10
            write(*,*) 'time correlation              ',tt11
            write(*,*) 'time save previous sfields    ',tt12
        endif
        
        deallocate(seismo)

    end subroutine
    
end
