module m_gen_acquisition
use m_sysio
use m_arrayop
use m_model, only:m

    private t_receiver, t_source, acqui_type, source_line, receiver_line
    
    type t_receiver
        real    :: x,y,z
        !integer :: icomp        !component: 1=P, 2=vx, 3=vy, 4=vz
        !integer :: nt
    end type
    
    type t_source
        real    :: x,y,z
        !integer :: icomp        !component: 1=P, 2=vx, 3=vy, 4=vz
        !integer :: nt
        !real,dimension(:),allocatable :: wavelet
        
        integer :: nrcv
        type(t_receiver),dimension(:),allocatable :: rcv  !should be in size of nrcv
    end type
    
    type t_acquisition
        !character(*) :: type
        integer :: nsrc
        type(t_source),dimension(:),allocatable :: src  !should be in size of nsrc
        
        integer :: iscomp      !source component: 1=P, 2=vx, 3=vy, 4=vz
        integer :: ircomp      !receiver component: 1=P, 2=vx, 3=vy, 4=vz
        
        integer :: nt
        real :: dt
    end type
    
    type(t_acquisition) :: acqui
    
    character(:),allocatable :: acqui_type, source_line, receiver_line
    
    contains
    
    subroutine gen_acquisition
        acqui_type=get_setup_char('ACQUI_TYPE')
        source_line=get_setup_char('SOURCE_LINE')
        receiver_line=get_setup_char('RECEIVER_LINE')
        
        if(acqui_type=='spread') then
            call gen_acqui_spread
        elseif(acqui_type=='streamer') then
            call gen_acqui_streamer
        else
            call hud('Sorry, only Spread or Streamer acquisitions are implemented now')
            stop
        endif
        
        call check_acqui
        
        acqui%iscomp=get_setup_int('SOURCE_COMPONENT')
        acqui%ircomp=get_setup_int('RECEIVER_COMPONENT')
        
        acqui%nt=get_setup_int('TIME_STEP')
        acqui%dt=get_setup_real('TIME_INTERVAL')
        
    end subroutine
    
    subroutine gen_acqui_spread
        real :: fsx,fsy,fsz, lsx,lsy,lsz, dsx,dsy,dsz
        real :: frx,fry,frz, lrx,lry,lrz, drx,dry,drz
        integer :: ns, nr
        
        read(source_line,*)   fsz,fsx,fsy, lsz,lsx,lsy, ns
        read(receiver_line,*) frz,frx,fry, lrz,lrx,lry, nr
        
        if(.not.m%is_cubic) then
            fsy=0; lsy=0;
            fry=0; lry=0;
        endif
        
        dsx=(lsx-fsx)/(ns-1)
        dsy=(lsy-fsy)/(ns-1)
        dsz=(lsz-fsz)/(ns-1)
        drx=(lrx-frx)/(nr-1)
        dry=(lry-fry)/(nr-1)
        drz=(lrz-frz)/(nr-1)
        
        ox=m%ox
        oy=m%oy
        oz=m%oz
        
        acqui%nsrc=ns
        if(allocated(acqui%src))deallocate(acqui%src)
        allocate(acqui%src(ns))
        
        
        sx=fsx; sy=fsy; sz=fsz
        do i=1,ns
            acqui%src(i)%x=sx-ox
            acqui%src(i)%y=sy-oy
            acqui%src(i)%z=sz-oz
            sx=sx+dsx
            sy=sy+dsy
            sz=sz+dsz
            
            !acqui%src(i)%icomp=iscomp
            
            acqui%src(i)%nrcv=nr
            if(allocated(acqui%src(i)%rcv))deallocate(acqui%src(i)%rcv)
            allocate(acqui%src(i)%rcv(nr))
            
            rx=frx; ry=fry; rz=frz;
            do j=1,nr
                acqui%src(i)%rcv(j)%x=rx-ox
                acqui%src(i)%rcv(j)%y=ry-oy
                acqui%src(i)%rcv(j)%z=rz-oz
                rx=rx+drx
                ry=ry+dry
                rz=rz+drz
                
                !acqui%src(i)%rcv(j)%icomp=ircomp
            enddo
        enddo
        
    end subroutine
    
    subroutine gen_acqui_streamer
        real :: fsx,fsy,fsz, lsx,lsy,lsz, dsx,dsy,dsz
        real :: foffx,foffy,fz, loffx,loffy,lz, doffx,doffy,dz
        integer :: ns, noff
        
        read(source_line,*)   fsz,fsx,fsy, lsz,lsx,lsy, ns
        read(receiver_line,*) fz,foffx,foffy, lz,loffx,loffy, noff
        
        if(.not.m%is_cubic) then
            fsy=0; lsy=0;
            foffy=0; loffy=0;
        endif
        
        dsx=(lsx-fsx)/(ns-1)
        dsy=(lsy-fsy)/(ns-1)
        dsz=(lsz-fsz)/(ns-1)
        doffx=(loffx-foffx)/(noff-1)
        doffy=(loffy-foffy)/(noff-1)
        dz=(lz-fz)/(noff-1)
        
        ox=m%ox
        oy=m%oy
        oz=m%oz
        
        acqui%nsrc=ns
        if(allocated(acqui%src))deallocate(acqui%src)
        allocate(acqui%src(ns))
        
        
        sx=fsx; sy=fsy; sz=fsz
        do i=1,ns
            acqui%src(i)%x=sx-ox
            acqui%src(i)%y=sy-oy
            acqui%src(i)%z=sz-oz
            
            acqui%src(i)%nrcv=noff
            if(allocated(acqui%src(i)%rcv))deallocate(acqui%src(i)%rcv)
            allocate(acqui%src(i)%rcv(noff))
            
            rx=sx+foffx; ry=sy+foffy; rz=fz;
            do j=1,noff
                acqui%src(i)%rcv(j)%x=rx-ox
                acqui%src(i)%rcv(j)%y=ry-oy
                acqui%src(i)%rcv(j)%z=rz-oz
                rx=rx+doffx
                ry=ry+doffy
                rz=rz+dz
            enddo
            
            sx=sx+dsx
            sy=sy+dsy
            sz=sz+dsz
            
        enddo
        
    end subroutine
    
    subroutine check_acqui
        !source position
        if(any(acqui%src(:)%x<0.)) then
            call hud('ERROR: some source x position is outside left bound of the model')
            stop
        endif
        if(any(acqui%src(:)%x>(m%nx-1)*m%dx)) then
            call hud('ERROR: some source x position is outside right bound of the model')
            stop
        endif
        if(any(acqui%src(:)%y<0.)) then
            call hud('ERROR: some source y position is outside front bound of the model')
            stop
        endif
        if(any(acqui%src(:)%y>(m%ny-1)*m%dy)) then
            call hud('ERROR: some source y position is outside rear bound of the model')
            stop
        endif
        if(any(acqui%src(:)%z<0.)) then
            call hud('ERROR: some source z position is outside top bound of the model')
            stop
        endif
        if(any(acqui%src(:)%z>(m%nz-1)*m%dz)) then
            call hud('ERROR: some source z position is outside bottom bound of the model')
            stop
        endif
        
        !receiver position
        do i=1,acqui%nsrc
            if(any(acqui%src(i)%rcv(:)%x<0.)) then
                if(mpiworld%is_master) write(*,*) 'ERROR: source#',i,'''s receiver x position is outside left bound of the model'
                stop
            endif
            if(any(acqui%src(i)%rcv(:)%x>(m%nx-1)*m%dx)) then
                if(mpiworld%is_master) write(*,*) 'ERROR: source#',i,'''s receiver x position is outside right bound of the model'
                stop
            endif
            if(any(acqui%src(i)%rcv(:)%y<0.)) then
                if(mpiworld%is_master) write(*,*) 'ERROR: source#',i,'''s receiver y position is outside front bound of the model'
                stop
            endif
            if(any(acqui%src(i)%rcv(:)%y>(m%ny-1)*m%dy)) then
                if(mpiworld%is_master) write(*,*) 'ERROR: source#',i,'''s receiver y position is outside rear bound of the model'
                stop
            endif
            if(any(acqui%src(i)%rcv(:)%z<0.)) then
                if(mpiworld%is_master) write(*,*) 'ERROR: source#',i,'''s receiver z position is outside top bound of the model'
                stop
            endif
            if(any(acqui%src(i)%rcv(:)%z>(m%nz-1)*m%dz)) then
                if(mpiworld%is_master) write(*,*) 'ERROR: source#',i,'''s receiver z position is outside bottom bound of the model'
                stop
            endif
        enddo
    end subroutine
    
end