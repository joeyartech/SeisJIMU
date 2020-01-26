program main
use m_mpienv
use m_suformat
use m_shot

    integer :: scomp, rcomp
    integer :: scalel, scalco
    character(:),allocatable :: data_file, sline, rline
    type(t_suformat),dimension(:),allocatable :: sudata
    logical alive
    integer :: file_size

    !character(:),allocatable :: job
    
    call init_mpiworld

    call hud('=================')
    call hud('     WELCOME     ')
    call hud('=================')
    
    call init_setup(istat)
    
    if(istat==0) then !print manual
        call print_manual
        call mpiworld_finalize
        stop
    endif
    
    
    !read setup
    nt=get_setup_int('TIME_STEP')
    dt=get_setup_real('TIME_INTERVAL')
    
    scomp=get_setup_int('SOURCE_COMPONENT',default=1)
    rcomp=get_setup_int('RECEIVER_COMPONENT',default=1)
    
    scalel=get_setup_int('SCALE_ELEVATION',default=0)
    scalco=get_setup_int('SCALE_COORDINATE',default=0)
    
    if(scalel> 31072) scalel= 31072
    if(scalel<-31072) scalel=-31072
    if(scalco> 31072) scalco= 31072
    if(scalco<-31072) scalco=-31072
    
    
    nshots=get_setup_int('NSHOTS',default=1)
    
    !assign shots to processors
    call build_shotlist(nshots)
    
    call hud('      START LOOP OVER SHOTS          ')
    
    do k=1,nshot_per_processor
        
        
        !read data
        shot%index=shotlist(k)
        write(shot%cindex,'(i0.4)') shot%index
        
        data_file=get_setup_char('FILE_DATA')//shot%cindex
        
        inquire(file=data_file, size=file_size, exist=alive) !file_size in bytes
        file_size=file_size/4 !file_size in sizeof(float)
        
        ntr=file_size/nt
        
        if(.not.alive) then
            write(*,*) 'ERROR: data file '//data_file//'does NOT exist!'
            stop
        endif
        
        if( ntr*nt /= file_size ) then
            write(*,*) 'ERROR: size of '//data_file//' is NOT a multiple of TIME_STEPS (nt)!'
            stop
        endif
        
        call alloc(dobs, nt,ntr)
        open(11,file=data_file,action='read',access='direct',recl=4*nt*ntr)
        read(11,rec=1) dobs
        close(11)
        
        
        !make header
        sline=get_setup_file('SOURCE_LINE')
        rline=get_setup_file('RECEIVER_LINE')
        
        open(11,file=sline,status='old')
        do i=1,shot%index-1; read(11,*); enddo
        read(11,*) shot%src%z, shot%src%x, shot%src%y
        close(11)
        
        if(allocated(shot%rcv)) deallocate(shot%rcv)
        allocate(shot%rcv(ntr))
        
        open(11,file=rline,status='old')
        do itr=1,ntr
            read(11,*) shot%rcv(itr)%z, shot%rcv(itr)%x, shot%rcv(itr)%y
        enddo
        close(11)
        
        
        !write su data
        if(allocated(sudata)) deallocate(sudata)
        allocate(sudata(ntr))
        call init_suheader(sudata,ntr)
        
        sudata(:)%hdr%tracl = [(itr,itr=1,ntr)]
        sudata(:)%hdr%fldr = shot%index
        sudata(:)%hdr%ns = nt
        sudata(:)%hdr%dt = dt*1e6
        sudata(:)%hdr%ntr = ntr
        
        !sudata(:)%hdr%year=
        !sudata(:)%hdr%day
        !sudata(:)%hdr%hour
        !sudata(:)%hdr%minute
        !sudata(:)%hdr%sec
        
        select case (rcomp)
            case (1) !pressure
            sudata(:)%hdr%trid = 1   !pressure
            case (2) !vx
            sudata(:)%hdr%trid = 14  !in-line component
            case (3) !vy
            sudata(:)%hdr%trid = 13  !croxx-line component
            case (4) !vz
            sudata(:)%hdr%trid = 12  !vertical component
        end select
        
        if(scalel>0) then
            sudata(:)%hdr%sdepth =  nint(shot%src%z    / scalel)
            sudata(:)%hdr%gelev  = -nint(shot%rcv(:)%z / scalel)
        elseif(scalel<0) then
            sudata(:)%hdr%sdepth =  nint(shot%src%z    * (-scalel))
            sudata(:)%hdr%gelev  = -nint(shot%rcv(:)%z * (-scalel))
        else
            sudata(:)%hdr%sdepth =  nint(shot%src%z)
            sudata(:)%hdr%gelev  = -nint(shot%rcv(:)%z)
        endif
        
        if(scalco>0) then
            sudata(:)%hdr%sx = nint(shot%src%x     / scalco)
            sudata(:)%hdr%sy = nint(shot%src%y     / scalco)
            sudata(:)%hdr%gx = nint(shot%rcv(:)%x  / scalco)
            sudata(:)%hdr%gy = nint(shot%rcv(:)%y  / scalco)
        elseif(scalco<0) then
            sudata(:)%hdr%sx = nint(shot%src%x     * (-scalco))
            sudata(:)%hdr%sy = nint(shot%src%y     * (-scalco))
            sudata(:)%hdr%gx = nint(shot%rcv(:)%x  * (-scalco))
            sudata(:)%hdr%gy = nint(shot%rcv(:)%y  * (-scalco))
        else
            sudata(:)%hdr%sx = nint(shot%src%x)
            sudata(:)%hdr%sy = nint(shot%src%y)
            sudata(:)%hdr%gx = nint(shot%rcv(:)%x)
            sudata(:)%hdr%gy = nint(shot%rcv(:)%y)
        endif
        
        sudata(:)%hdr%offset = sudata(:)%hdr%gx - sudata(:)%hdr%sx
        
        sudata(:)%hdr%scalco = scalco
        sudata(:)%hdr%scalel = scalel
        
!         sudata(:)%hdr%offset = sqrt( (shot%src%x-shot%rcv(:)%x)**2 &
!                                     +(shot%src%y-shot%rcv(:)%y)**2 )
        
        do itr=1,ntr
            call alloc(sudata(itr)%trace, nt)
            sudata(itr)%trace(:) = dobs(:,itr)
        enddo
        
        open(12,file=data_file//'.su',action='write',access='stream') !access='direct',recl=4*(ns+60))
        write(12) (sudata(itr)%hdr,sudata(itr)%trace, itr=1,ntr)
        close(12)
        
        if(mpiworld%is_master) write(*,*) 'Shot# '//shot%cindex//' has written ',ntr,' traces, each trace has',nt,'samples.'
        
        call hud('Data write success.')

        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')
    
    
    call mpiworld_finalize
    
end


subroutine print_manual
use m_mpienv
use m_field

    if(mpiworld%is_master) then
        write(*,'(a)') ""
        write(*,'(a)') "----------------------------"
        write(*,'(a)') "To launch the program, do:"
        write(*,'(a)') "----------------------------"
        write(*,'(a)') ""
        write(*,'(a)') "bash $ mpirun -np $np ADD_SUHEADER setup.in"
        write(*,'(a)') ""
        write(*,'(a)') "----------------------------"
        write(*,'(a)') "Mandatory items in setup.in:"
        write(*,'(a)') "----------------------------"
        write(*,'(a)') ""
        write(*,'(a)') "TIME_STEP               500                #Total number of time step (nt)"
        write(*,'(a)') "TIME_INTERVAL           0.006              #Time step interval (dt)"
        write(*,'(a)') ""
        write(*,'(a)') "SOURCE_COMPONENT        1"
        write(*,'(a)') "RECEIVER_COMPONENT      1"
        write(*,'(a)') "                        #source & receiver component: 1=Isotropic Pressure, 2=Particle velocity in x (vx), 3=in y (vy), 4=in z (vz)"
        write(*,'(a)') "                        #Multi-component data will be considered and developed in future."
        write(*,'(a)') ""
        write(*,'(a)') "SCALE_ELEVATION         1"
        write(*,'(a)') "SCALE_COORDINATE        1"       
        write(*,'(a)') ""
        write(*,'(a)') "NSHOTS                  1                  #Number of shots to add header"
        write(*,'(a)') ""
        write(*,'(a)') "FILE_DATA               'synth_data_'      #Prefix of input data filename, which has the form synth_data_????"
        write(*,'(a)') "                                           #where ???? represents 4 digits starting from 0001"
        write(*,'(a)') ""
        write(*,'(a)') "SOURCE_LINE             'sline.txt'"
        write(*,'(a)') "RECEIVER_LINE           'rline.txt'"
        write(*,'(a)') "                        #Source & receiver positions:"
        write(*,'(a)') "                        #ACQUI_TYPE=='irregularOBN': Needing file name in SOURCE_LINE, which contains a Nx3 matrix (columns separated by white spaces),"
        write(*,'(a)') "                        #                  where N=total number of sources, and 3 columns specify (z,x,y) coordinate of sources"
        write(*,'(a)') "                        #                  Same for RECEIVER_LINE"
        write(*,'(a)') "#Other acqui types or source & receiver lines will be considered and developed in future."
        write(*,'(a)') ""
        write(*,'(a)') "---------------------------"
        write(*,'(a)') "Optional items in setup.in:"
        write(*,'(a)') "---------------------------"      
        write(*,'(a)') ""
        write(*,'(a)') "---------------------------"
        write(*,'(a)') "Notes:"
        write(*,'(a)') "---------------------------"
        write(*,'(a)') ""

    endif

end subroutine