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
        source_line=get_setup_char('SOURCE_LINE','SOURCE_AREA')
        receiver_line=get_setup_char('RECEIVER_LINE','RECEIVER_AREA')
        
        if(acqui_type=='spread') then
            call gen_acqui_spread
        elseif(acqui_type=='streamer') then
            call gen_acqui_streamer
        elseif(acqui_type=='irregularOBN') then
            call gen_acqui_irregularOBN
        elseif(acqui_type=='spread3D') then
            call gen_acqui_spread3D
        else
            call hud('Sorry, only spread, streamer, or irregular OBN acquisition geometries are implemented now')
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
        
        acqui%nsrc=ns
        if(allocated(acqui%src))deallocate(acqui%src)
        allocate(acqui%src(ns))
        
        sx=fsx; sy=fsy; sz=fsz
        do i=1,ns
            acqui%src(i)%x=sx
            acqui%src(i)%y=sy
            acqui%src(i)%z=sz
            sx=sx+dsx
            sy=sy+dsy
            sz=sz+dsz
            
            !acqui%src(i)%icomp=iscomp
            
            acqui%src(i)%nrcv=nr
            if(allocated(acqui%src(i)%rcv))deallocate(acqui%src(i)%rcv)
            allocate(acqui%src(i)%rcv(nr))
            
            rx=frx; ry=fry; rz=frz;
            do j=1,nr
                acqui%src(i)%rcv(j)%x=rx
                acqui%src(i)%rcv(j)%y=ry
                acqui%src(i)%rcv(j)%z=rz
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
        
        acqui%nsrc=ns
        if(allocated(acqui%src))deallocate(acqui%src)
        allocate(acqui%src(ns))
        
        sx=fsx; sy=fsy; sz=fsz
        do i=1,ns
            acqui%src(i)%x=sx
            acqui%src(i)%y=sy
            acqui%src(i)%z=sz
            
            acqui%src(i)%nrcv=noff
            if(allocated(acqui%src(i)%rcv))deallocate(acqui%src(i)%rcv)
            allocate(acqui%src(i)%rcv(noff))
            
            rx=sx+foffx; ry=sy+foffy; rz=fz;  !NOTE no need foffz for rz !
            do j=1,noff
                acqui%src(i)%rcv(j)%x=rx
                acqui%src(i)%rcv(j)%y=ry
                acqui%src(i)%rcv(j)%z=rz
                rx=rx+doffx
                ry=ry+doffy
                rz=rz+dz
            enddo
            
            sx=sx+dsx
            sy=sy+dsy
            sz=sz+dsz
            
        enddo
        
    end subroutine

    subroutine gen_acqui_irregularOBN
        character(80) :: text

        !read sources
        open(13,file=source_line,action='read')

            !count number of sources
            n=0
            do
                read (13,*,iostat=msg) z,x,y
                if(msg/=0) exit
                n=n+1
            end do
            if(mpiworld%is_master) write(*,*) 'Will read',n,'sources.'
            if(allocated(acqui%src))deallocate(acqui%src)
            allocate(acqui%src(n))
            acqui%nsrc=n

        !read source positions
        rewind(13)
            i=1
            do
                read (13,*,iostat=msg) z,x,y
                if(msg/=0) exit
                acqui%src(i)%z=z
                acqui%src(i)%x=x
                acqui%src(i)%y=y
                i=i+1
            end do

        close(13)


        !read receivers
        open(15,file=receiver_line,action='read')

            !count number of receivers
            n=0
            do
                read (15,*,iostat=msg) z,x,y
                if(msg/=0) exit
                n=n+1
            end do
            if(mpiworld%is_master) write(*,*) 'Will read',n,'receivers.'
            do i=1,acqui%nsrc
                if(allocated(acqui%src(i)%rcv))deallocate(acqui%src(i)%rcv)
                allocate(acqui%src(i)%rcv(n))
            enddo
            acqui%src(:)%nrcv=n

        !read receiver positions
        rewind(15)
            i=1
            do
                read (15,*,iostat=msg) z,x,y
                if(msg/=0) exit
                do j=1,acqui%nsrc
                    acqui%src(j)%rcv(i)%z=z
                    acqui%src(j)%rcv(i)%x=x
                    acqui%src(j)%rcv(i)%y=y
                enddo
                i=i+1
            end do
        
        close(15)

    end subroutine

    subroutine gen_acqui_spread3D
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Quadrilateral Geometry  !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !  IL directio ->
        ! XL  P1---------------P2
        ! dir -------------------
        ! |   -------------------
        ! v   P3---------------P4
        !
        !4 anchor points, P1-4, and number of points in IL and XL directions
        !should be read from setup (SOURCE_AREA, RECEIVER_AREA)
        !other points are interpolated from P1-4 and nil, nxl

        real :: frx, lrx, fry, lry, frz, lrz
        
        read(source_line,*)    sz1,sx1,sy1, sz2,sx2,sy2, sz3,sx3,sy3, sz4,sx4,sy4, nsil, nsxl
        read(receiver_line,*)  rz1,rx1,ry1, rz2,rx2,ry2, rz3,rx3,ry3, rz4,rx4,ry4, nril, nrxl
        
        if(.not.m%is_cubic) then
            sx3=sx1; sy3=sy1; sz3=sz1;
            sx4=sx2; sy4=sy2; sz4=sz2;  nsxl=1

            rx3=rx1; ry3=ry1; rz3=rz1;
            rx4=rx2; ry4=ry2; rz4=rz2;  nrxl=1
        endif

        ns=nsil*nsxl
        nr=nril*nrxl

        acqui%nsrc=ns
        if(allocated(acqui%src))deallocate(acqui%src)
        allocate(acqui%src(ns))

        do ixl=1,nsxl

            fsx = ( sx1*(nsxl-ixl) + sx3*(ixl-1) )/(nsxl-1)
            fsy = ( sy1*(nsxl-ixl) + sy3*(ixl-1) )/(nsxl-1)
            fsz = ( sz1*(nsxl-ixl) + sz3*(ixl-1) )/(nsxl-1)

            lsx = ( sx2*(nsxl-ixl) + sx4*(ixl-1) )/(nsxl-1)
            lsy = ( sy2*(nsxl-ixl) + sy4*(ixl-1) )/(nsxl-1)
            lsz = ( sz2*(nsxl-ixl) + sz4*(ixl-1) )/(nsxl-1)

        do iil=1,nsil

            sx = ( fsx*(nsil-iil) + lsx*(iil-1) )/(nsil-1)
            sy = ( fsy*(nsil-iil) + lsy*(iil-1) )/(nsil-1)
            sz = ( fsz*(nsil-iil) + lsz*(iil-1) )/(nsil-1)

            i  = iil + (ixl-1)*nsil

            acqui%src(i)%x=sx
            acqui%src(i)%y=sy
            acqui%src(i)%z=sz
            
            !acqui%src(i)%icomp=iscomp

        
            acqui%src(i)%nrcv=nr
            if(allocated(acqui%src(i)%rcv))deallocate(acqui%src(i)%rcv)
            allocate(acqui%src(i)%rcv(nr))
            
            do jxl=1,nrxl

                frx = ( rx1*(nrxl-jxl) + rx3*(jxl-1) )/(nrxl-1)
                fry = ( ry1*(nrxl-jxl) + ry3*(jxl-1) )/(nrxl-1)
                frz = ( rz1*(nrxl-jxl) + rz3*(jxl-1) )/(nrxl-1)

                lrx = ( rx2*(nrxl-jxl) + rx4*(jxl-1) )/(nrxl-1)
                lry = ( ry2*(nrxl-jxl) + ry4*(jxl-1) )/(nrxl-1)
                lrz = ( rz2*(nrxl-jxl) + rz4*(jxl-1) )/(nrxl-1)

            do jil=1,nril

                rx = ( frx*(nril-jil) + lrx*(jil-1) )/(nril-1)
                ry = ( fry*(nril-jil) + lry*(jil-1) )/(nril-1)
                rz = ( frz*(nril-jil) + lrz*(jil-1) )/(nril-1)

                j  = jil + (jxl-1)*nril

                acqui%src(i)%rcv(j)%x=rx
                acqui%src(i)%rcv(j)%y=ry
                acqui%src(i)%rcv(j)%z=rz
                
                !acqui%src(i)%rcv(j)%icomp=ircomp

            enddo

            enddo

        enddo

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