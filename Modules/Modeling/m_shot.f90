module m_shot
use m_System
use m_hicks
use m_wavelet
use m_resampler
use m_matchfilter
use m_model

    private

    !abbreviation:
    ! s,src : source
    ! r,rcv : receiver
    ! o : origin
    ! d : increment, spacing
    ! n : number of
    ! z : depth
    ! x : inline
    ! y : crossline
    ! pos : position
    ! comp : component
    
    type,public :: t_source
        real    :: z,x,y
        integer :: iz,ix,iy
        integer :: ifz,ilz,ifx,ilx,ify,ily
        character(4) :: comp
        real,dimension(:,:,:),allocatable :: interp_coef
    end type

    type,public :: t_receiver
        real    :: z,x,y, aoffset
        integer :: iz,ix,iy
        integer :: ifz,ilz,ifx,ilx,ify,ily
        logical :: is_badtrace=.false.
        character(4) :: comp
        real,dimension(:,:,:),allocatable :: interp_coef
    end type
    
    type,public :: t_shot
        integer :: index
        character(4) :: sindex
        
        real,dimension(:),allocatable :: wavelet
        integer :: nt
        real :: dt, fmin, fmax, fpeak

        logical :: if_hicks=.true.

        type(t_source) :: src
        type(t_receiver),dimension(:),allocatable :: rcv
        integer :: nrcv   !=size(rcv)

        real,dimension(:,:),allocatable :: dobs !observed seismogram
        real,dimension(:,:),allocatable :: dsyn !synthetic seismogram
        !real,dimension(:,:),allocatable :: dres !residual seismogram
        real,dimension(:,:),allocatable :: dadj !adjoint source seismogram
        
        contains
        procedure :: init
        procedure :: read_from_setup
        procedure :: read_from_data
        procedure :: set_var_time
        procedure :: set_var_space
        procedure :: update_wavelet
        procedure :: update_adjsource
        procedure :: write
        
    end type
    
    type(t_shot),public :: shot
    
    contains
    
    subroutine init(self,index)
        class(t_shot) :: self
        !character(*) :: sindex

        !self%sindex=sindex
        !self%index=str2int(sindex)

        self%index=index
        self%sindex=num2str(index,'(i0.4)') !add leading zeros
    
    end subroutine

    subroutine read_from_setup(self)
        class(t_shot) :: self

        logical,save :: is_first_in=.true.
        type(t_string),dimension(:),allocatable :: scomp

        if(is_first_in) then
            self%nt=setup%get_int('TIME_STEP','NT',o_mandatory=1)
            self%dt=setup%get_real('TIME_INTERVAL','DT',o_mandatory=1)
            
            scomp=setup%get_strs('SOURCE_COMPONENT','SCOMP',o_default='p')
            if(size(scomp)>1) call hud('SeisJIMU only considers the 1st component from SOURCE_COMPONENT: '//scomp(1)%s)
            self%src%comp=scomp(1)%s

            is_first_in=.false.

        endif

        select case (setup%get_str('ACQUI_GEOMETRY',o_default='spread'))
        case ('spread');    call geometry_spread(self)
        case ('streamer');  call geometry_streamer(self)
        case ('emptybox');  call geometry_emptybox(self)    
        case ('spread_irregular')
            call geometry_spread_irregular(self)
        case ('spread3D')
        end select        

        !shift positions to be 0-based
        !then m%oz,ox,oy is no longer of use before writing shots
        self%src%z   =self%src%z    - m%oz
        self%src%x   =self%src%x    - m%ox
        self%src%y   =self%src%y    - m%oy
        self%rcv(:)%z=self%rcv(:)%z - m%oz
        self%rcv(:)%x=self%rcv(:)%x - m%ox
        self%rcv(:)%y=self%rcv(:)%y - m%oy

        call check_range(self)

    end subroutine

    subroutine read_from_data(self)
        class(t_shot) :: self

        type(t_suformat) :: sudata
        character(:),allocatable :: str, zerophase, locut, hicut
        character(:),allocatable :: fstoplo,fpasslo,fpasshi,fstophi
        character(:),allocatable :: astoplo,apasslo,apasshi,astophi

        call sudata%read(setup%get_str('FILE_DATA_PREFIX')//self%sindex//'.su',self%sindex)

        shot%nrcv=sudata%ntr
        if(allocated(self%rcv))deallocate(self%rcv)
        allocate(self%rcv(shot%nrcv))

        scalel=sudata%hdrs(1)%scalel !assume same scalel for all traces
        scalco=sudata%hdrs(1)%scalco !assume same scalco for all traces

        if(scalel==0.) scalel=1.
        if(scalco==0.) scalco=1.

        if(scalel<0.) scalel=-1./scalel
        if(scalco<0.) scalco=-1./scalco

        !source & receiver geometry
        self%src%z=sudata%hdrs(1)%sdepth*scalel !assume single shot
        self%src%x=sudata%hdrs(1)%sx    *scalco
        self%src%y=sudata%hdrs(1)%sy    *scalco
        
        self%src%comp='p' !don't know which su header tells this info..
        
        do i=1,shot%nrcv
            self%rcv(i)%z=-sudata%hdrs(i)%gelev*scalel
            self%rcv(i)%x= sudata%hdrs(i)%gx   *scalco
            self%rcv(i)%y= sudata%hdrs(i)%gy   *scalco

            self%rcv(i)%is_badtrace = sudata%hdrs(i)%trid==2 .or. sudata%hdrs(i)%trid==3  !dead or dummy trace
            
            select case (sudata%hdrs(i)%trid)
            case (11)
                self%rcv(i)%comp='p'  !pressure
            case (12)
                self%rcv(i)%comp='vz' !vertical velocity
            case (13)
                self%rcv(i)%comp='vy' !horizontal velocity
            case (14)
                self%rcv(i)%comp='vx'
            case default
                self%rcv(i)%comp='p'
            end select
            
        enddo

        !shift positions to be 0-based
        !then m%oz,ox,oy is no longer of use before writing shots
        self%src%z   =self%src%z    - m%oz
        self%src%x   =self%src%x    - m%ox
        self%src%y   =self%src%y    - m%oy
        self%rcv(:)%z=self%rcv(:)%z - m%oz
        self%rcv(:)%x=self%rcv(:)%x - m%ox
        self%rcv(:)%y=self%rcv(:)%y - m%oy

        call check_range(self)


        !if want to filter the traces
        str=setup%get_str('FILTER_DATA','BFILT',o_default='none')

        if(str/='none') then
            
            !default
            zerophase=read_key(str,'zerophase',o_default='1')
            locut=read_key(str,'locut',o_default='1')
            hicut=read_key(str,'hicut',o_default='1')
            rdt = 0.5/sudata%dt
            fstoplo=read_key(str,'fstoplo',o_default=num2str(.10*rdt))
            fpasslo=read_key(str,'fpasslo',o_default=num2str(.15*rdt))
            fpasshi=read_key(str,'fpasshi',o_default=num2str(.40*rdt))
            fstophi=read_key(str,'fstophi',o_default=num2str(.55*rdt))

            astoplo=read_key(str,'astoplo',o_default=num2str(.05))
            apasslo=read_key(str,'apasslo',o_default=num2str(.95))
            apasshi=read_key(str,'apasshi',o_default=num2str(.95))
            astophi=read_key(str,'astophi',o_default=num2str(.05))

            !butterworth filtering
            call hud('subfilt < data > data_filt '//&
                ' zerophase='//zerophase//' locut='//locut//' hicut='//hicut// &
                ' fstoplo='//fstoplo//' fpasslo='//fpasslo//' fpasshi='//fpasshi//' fstophi='//fstophi// &
                ' astoplo='//astoplo//' apasslo='//apasslo//' apasshi='//apasshi//' astophi='//astophi   )
            
            call butterworth(sudata%trs,sudata%ns,sudata%ntr,sudata%dt, &
                ois_zerophase=either(.true.,.false.,zerophase=='1'),oif_locut=either(.true.,.false.,locut=='1'),oif_hicut=either(.true.,.false.,hicut=='1'),       &
                o_fstoplo=str2real(fstoplo),o_fpasslo=str2real(fpasslo),o_fpasshi=str2real(fpasshi),o_fstophi=str2real(fstophi), &
                o_astoplo=str2real(astoplo),o_apasslo=str2real(apasslo),o_apasshi=str2real(apasshi),o_astophi=str2real(astophi)  )

            sudata%hdrs(:)%lcf=str2int(fstoplo)
            sudata%hdrs(:)%hcf=str2int(fstophi)
            
            call sudata%write('dobs_filt_'//self%sindex)

        endif

        !load traces
        self%nt=sudata%ns !assume all traces have same ns
        self%dt=sudata%dt !assume all traces have same dt

        call alloc(self%dobs,self%nt,self%nrcv)
        self%dobs=sudata%trs

    end subroutine

    subroutine set_var_time(self)
        class(t_shot) :: self

        character(:),allocatable :: file, str
        type(t_suformat) :: wavelet

        self%fpeak=setup%get_real('PEAK_FREQUENCY','FPEAK')
        self%fmax=setup%get_real('MAX_FREQUENCY','FMAX',o_default=num2str(self%fpeak*2.))
        !self%fmin=setup%get_real('MIN_FREQUENCY','FMIN',o_default='1')

        !source time function, which should not have dt, dx info
        file=setup%get_file('FILE_SRC_TIME_FUNC','FILE_WAVELET')

        if(file=='') then !not given
            if(setup%get_str('WAVELET_TYPE',o_default='sinexp')=='sinexp') then
                call hud('Use filtered sinexp wavelet')
                self%wavelet=sinexp(self%nt,self%dt,self%fpeak)
            else
                call hud('Use Ricker wavelet')
                self%wavelet=ricker(self%nt,self%dt,self%fpeak)
            endif

        else !wavelet file exists
            call alloc(self%wavelet,self%nt)
            call wavelet%read(file)
            call resampler(wavelet%trs(:,1),self%wavelet,1,&
                            din=wavelet%dt,nin=wavelet%ns, &
                            dout=self%dt,  nout=self%nt)

        endif

        str=setup%get_str('WAVELET_SCALING')
        if(str=='') then
        elseif(str=='by dx3dt' .or. str=='by dtdx3') then
            self%wavelet=self%wavelet/self%dt*m%cell_volume
        else
            self%wavelet=self%wavelet*str2real(str)
        endif

        if(mpiworld%is_master) call suformat_write('wavelet',self%wavelet,self%nt,ntr=1,o_dt=self%dt)

    end subroutine

    subroutine set_var_space(self,is_fdsg)
        class(t_shot) :: self
        logical :: is_fdsg

        logical,save :: is_first_in=.true.

        !absolute offset
        do ir=1,self%nrcv
            self%rcv(ir)%aoffset=sqrt( (self%src%z-self%rcv(ir)%z)**2 &
                                      +(self%src%x-self%rcv(ir)%x)**2 &
                                      +(self%src%y-self%rcv(ir)%y)**2 )
        enddo

        if(is_first_in) then
            !Hicks interpolation
            self%if_hicks=setup%get_bool('IF_HICKS',o_default='T')

            call hicks_init(m%dz,m%dx,m%dy,m%is_cubic,m%is_freesurface)
            
            !velocity grids are shift from normal stress grids by half grid length (staggered-grid FD)
            !e.g. v(1) is actually v[0.5], v(2) is actually v[1.5],
            !thus to correctly Hicks-interpolate the recorded s/r positions of forces,
            !we need to artificially add 0.5 to the positions.
            halfz=either(m%dz/2.,0.,is_fdsg)
            halfx=either(m%dx/2.,0.,is_fdsg)
            halfy=either(m%dy/2.,0.,is_fdsg)

            is_first_in=.false.

        endif

        !source
        select case (self%src%comp)
        case('p') !explosive source or non-vertical force
            call hicks_put_position(self%src%z,       self%src%x,       self%src%y)
        case('vz') !vertical force
            call hicks_put_position(self%src%z+halfz, self%src%x,       self%src%y)
        case('vx')
            call hicks_put_position(self%src%z,       self%src%x+halfx, self%src%y)
        case('vy')
            call hicks_put_position(self%src%z,       self%src%x,       self%src%y+halfy)
        end select

        call hicks_get_position(self%src%ifz, self%src%ifx, self%src%ify,&
                                self%src%iz , self%src%ix , self%src%iy ,&
                                self%src%ilz, self%src%ilx, self%src%ily )

        select case (self%src%comp)
        case('p') !explosive source or non-vertical force
            call hicks_get_coefficient('antisymm', self%src%interp_coef)
        case('vz') !vertical force
            call hicks_get_coefficient('symmetric',self%src%interp_coef)
        case default
            call hicks_get_coefficient('truncate', self%src%interp_coef)
        end select

        !receivers
        do i=1,self%nrcv

            select case (self%rcv(i)%comp)
            case('p') !explosive source or non-vertical force
                call hicks_put_position(self%rcv(i)%z,       self%rcv(i)%x,       self%rcv(i)%y)
            case('vz') !vertical force
                call hicks_put_position(self%rcv(i)%z+halfz, self%rcv(i)%x,       self%rcv(i)%y)
            case('vx')
                call hicks_put_position(self%rcv(i)%z,       self%rcv(i)%x+halfx, self%rcv(i)%y)
            case('vy')
                call hicks_put_position(self%rcv(i)%z,       self%rcv(i)%x,       self%rcv(i)%y+halfy)
            end select

            call hicks_get_position(self%rcv(i)%ifz, self%rcv(i)%ifx, self%rcv(i)%ify,&
                                    self%rcv(i)%iz , self%rcv(i)%ix , self%rcv(i)%iy ,&
                                    self%rcv(i)%ilz, self%rcv(i)%ilx, self%rcv(i)%ily )

            select case (self%rcv(i)%comp)
            case('p') !explosive source or non-vertical force
                call hicks_get_coefficient('antisymm', self%rcv(i)%interp_coef)
            case('vz') !vertical force
                call hicks_get_coefficient('symmetric',self%rcv(i)%interp_coef)
            case default
                call hicks_get_coefficient('truncate', self%rcv(i)%interp_coef)
            end select

        enddo

        if(.not.self%if_hicks) then
            self%src%ifz=0; self%src%ilz=0
            self%src%ifx=0; self%src%ilx=0
            self%src%ify=0; self%src%ily=0
            self%rcv(:)%ifz=0; self%rcv(:)%ilz=0
            self%rcv(:)%ifx=0; self%rcv(:)%ilx=0
            self%rcv(:)%ify=0; self%rcv(:)%ily=0
        endif

        if(mpiworld%is_master) then
            write(*,*)'================================='
            write(*,*)'Shot# '//self%sindex//' info:'
            write(*,*)'================================='
            write(*,*)'  nt,dt:',self%nt,self%dt
            write(*,*)'---------------------------------'
            write(*,*)'S/R positions after removing m%oz,ox,iy,'
            write(*,*)'  sz,iz:',self%src%z,self%src%iz
            write(*,*)'  sx,ix:',self%src%x,self%src%ix
            write(*,*)'  sy,iy:',self%src%y,self%src%iy
            write(*,*)'  ifz,ilz:',self%src%ifz,self%src%ilz
            write(*,*)'  ifx,ilx:',self%src%ifx,self%src%ilx
            write(*,*)'  ify,ily:',self%src%ify,self%src%ily
            write(*,*)'---------------------------------'
            write(*,*)'  minmax rz,iz:',minval(self%rcv(:)%z),maxval(self%rcv(:)%z),minval(self%rcv(:)%iz),maxval(self%rcv(:)%iz)
            write(*,*)'  minmax rx,ix:',minval(self%rcv(:)%x),maxval(self%rcv(:)%x),minval(self%rcv(:)%ix),maxval(self%rcv(:)%ix)
            write(*,*)'  minmax ry,iy:',minval(self%rcv(:)%y),maxval(self%rcv(:)%y),minval(self%rcv(:)%iy),maxval(self%rcv(:)%iy)
            write(*,*)'  minmax ifz,ilz:',minval(self%rcv(:)%ifz),maxval(self%rcv(:)%ifz),minval(self%rcv(:)%ilz),maxval(self%rcv(:)%ilz)
            write(*,*)'  minmax ifx,ilx:',minval(self%rcv(:)%ifx),maxval(self%rcv(:)%ifx),minval(self%rcv(:)%ilx),maxval(self%rcv(:)%ilx)
            write(*,*)'  minmax ify,ily:',minval(self%rcv(:)%ify),maxval(self%rcv(:)%ify),minval(self%rcv(:)%ily),maxval(self%rcv(:)%ily)
            write(*,*)'  nrcv:',self%nrcv
            write(*,*)'---------------------------------'
        endif

    end subroutine

    subroutine update_wavelet(self)
        class(t_shot) :: self

        call matchfilter_estimate(self%dsyn,self%dobs,self%nt,self%nrcv)!,self%index)
        
        call matchfilter_apply_to_wavelet(self%wavelet)
        
        call matchfilter_apply_to_data(self%dsyn)

        if(mpiworld%is_master) call suformat_write('wavelet',self%wavelet,self%nt,1,self%dt,o_mode='append')

    end subroutine
    
    subroutine update_adjsource(self)
        class(t_shot) :: self

        call matchfilter_correlate_filter_residual(self%dadj)

    end subroutine

    subroutine geometry_spread(shot)
        type(t_shot) :: shot

        logical,save :: is_first_in=.true.
        real,dimension(3),save :: fs=0.,fr=0.
        real,dimension(3),save :: ds=0.,dr=0.
        integer,save :: nr=1
        type(t_string),dimension(:),allocatable,save :: rcomp

        if(is_first_in) then
            fs=setup%get_reals('SOURCE_FIRST',  'FS',o_default='0. 0. 0.')
            fr=setup%get_reals('RECEIVER_FIRST','FR',o_default='0. 0. 0.')

            ds=setup%get_reals('SOURCE_SPACING', 'DS',o_default='0. 0. 0.')
            dr=setup%get_reals('RECEIVER_SPACING','DR',o_default='0. 0. 0.')
            
            nr=setup%get_int('NUMBER_RECEIVER','NR',o_default='1')
            rcomp=setup%get_strs('RECEIVER_COMPONENT','RCOMP',o_default='p')

            if(.not.m%is_cubic) then
                fs(3)=0.; fr(3)=0.
                ds(3)=0.; dr(3)=0.
            endif

            is_first_in=.false.
        
        endif

        !source side
        shot%src%z=fs(1)+(shot%index-1)*ds(1)
        shot%src%x=fs(2)+(shot%index-1)*ds(2)
        shot%src%y=fs(3)+(shot%index-1)*ds(3)

        !receiver side
        shot%nrcv=nr*size(rcomp)

        if(allocated(shot%rcv)) deallocate(shot%rcv)
        allocate(shot%rcv(shot%nrcv))

        do j=0,size(rcomp)-1
            do i=0,nr-1
                shot%rcv(1+i+j*nr)%comp = rcomp(j+1)%s
                shot%rcv(1+i+j*nr)%z = fr(1) + i*dr(1)
                shot%rcv(1+i+j*nr)%x = fr(2) + i*dr(2)
                shot%rcv(1+i+j*nr)%y = fr(3) + i*dr(3)
            enddo
        enddo

    end subroutine

    subroutine geometry_streamer(shot)
        type(t_shot) :: shot

        logical,save :: is_first_in=.true.
        real,dimension(3),save :: fs=0.,foff=0.
        real,dimension(3),save :: ds=0.,doff=0.
        integer,save :: noff=1
        type(t_string),dimension(:),allocatable,save :: rcomp

        if(is_first_in) then
            fs  =setup%get_reals('SOURCE_FIRST','FS',o_default='0. 0. 0.')
            foff=setup%get_reals('OFFSET_FIRST','FOFF',o_default='0. 0. 0.')

            ds  =setup%get_reals('SOURCE_SPACING', 'DS',o_default='0. 0. 0.')
            doff=setup%get_reals('OFFSET_SPACING', 'DOFF',o_default='0. 0. 0.')

            noff=setup%get_int('NUMBER_OFFSET','NOFF',o_default='1')
            rcomp=setup%get_strs('RECEIVER_COMPONENT','RCOMP',o_default='p')

            if(.not.m%is_cubic) then
                fs(3)=0.; foff(3)=0.
                ds(3)=0.; doff(3)=0.
            endif

            is_first_in=.false.

        endif

        !source side
        shot%src%z=fs(1)+(shot%index-1)*ds(1)
        shot%src%x=fs(2)+(shot%index-1)*ds(2)
        shot%src%y=fs(3)+(shot%index-1)*ds(3)

        !receiver side
        shot%nrcv=noff*size(rcomp)

        if(allocated(shot%rcv)) deallocate(shot%rcv)
        allocate(shot%rcv(shot%nrcv))

        do j=0,size(rcomp)-1
            do i=0,noff-1
                shot%rcv(1+i+j*noff)%comp = rcomp(j+1)%s
                shot%rcv(1+i+j*noff)%z =              foff(1) + i*doff(1)
                shot%rcv(1+i+j*noff)%x = shot%src%x + foff(2) + i*doff(2)
                shot%rcv(1+i+j*noff)%y = shot%src%y + foff(3) + i*doff(3)
            enddo
        enddo

    end subroutine

    subroutine geometry_spread_irregular(shot)
        type(t_shot) :: shot

        logical,save :: is_first_in=.true.
        character(:),allocatable,save :: spos, rpos
        integer, save :: nr
        type(t_string),dimension(:),allocatable,save :: rcomp

        if(is_first_in) then            
            spos=setup%get_file('FILE_SOURCE_POSITION','SPOS',o_mandatory=1)
            rpos=setup%get_file('FILE_RECEIVER_POSITION','RPOS',o_mandatory=1)

            rcomp=setup%get_strs('RECEIVER_COMPONENT','RCOMP',o_default='p')

            is_first_in=.false.

        endif

        !read source
        open(13,file=spos,action='read')

        k=0
        do
            read (13,*,iostat=msg) shot%src%z,shot%src%x,shot%src%y
            if(msg/=0) exit
            k=k+1
            if(k==shot%index) exit
        enddo

        close(13)


        !read receivers
        open(15,file=rpos,action='read')
        !count number of receivers
        nr=0
        do
            read (15,*,iostat=msg) z,x,y
            if(msg/=0) exit
            nr=nr+1
        end do
        if(mpiworld%is_master) write(*,*) 'Will read',nr,'receiver positions.'

        shot%nrcv=nr*size(rcomp)

        if(allocated(shot%rcv)) deallocate(shot%rcv)
        allocate(shot%rcv(shot%nrcv))
        
        i=0

        do l=1,size(rcomp)
            rewind(15)

            do
                i=i+1
                read (15,*,iostat=msg) shot%rcv(i)%z,shot%rcv(i)%x,shot%rcv(i)%y
                shot%rcv(i)%comp=rcomp(l)%s
                if(msg/=0) exit
            end do

        enddo

        close(15)

    end subroutine

    subroutine geometry_emptybox(shot)
        type(t_shot) :: shot

        logical,save :: is_first_in=.true.
        real,dimension(3),save :: spos
        real,dimension(4),save :: tline, rline, bline, lline
        real,dimension(4) :: line
        integer,save :: nr=1
        !type(t_string),dimension(:),allocatable,save :: rcomp
        character(:),allocatable,save :: rcomp

        if(is_first_in) then
            spos=setup%get_reals('SOURCE_POSITION',  'SPOS',o_default='0. 0. 0.') !for now just 1 shot
            
            lline=setup%get_reals('LEFT_RCV_LINE',o_mandatory=4)
            tline=setup%get_reals('TOP_RCV_LINE',o_mandatory=4) !const_z, fx, dx, nx
            rline=setup%get_reals('RIGHT_RCV_LINE',o_mandatory=4) !const_x, fz, dz, nz
            bline=setup%get_reals('BOTTOM_RCV_LINE',o_mandatory=4)

            nr=tline(4)+rline(4)+bline(4)+lline(4)

            !rcomp=setup%get_strs('RECEIVER_COMPONENT','RCOMP',o_default='p')
            rcomp=setup%get_str('RECEIVER_COMPONENT','RCOMP',o_default='p')

            ! if(.not.m%is_cubic) then
            !     fs(3)=0.; fr(3)=0.
            !     ds(3)=0.; dr(3)=0.
            ! endif

            is_first_in=.false.
        
        endif

        !source side
        shot%src%z=spos(1)
        shot%src%x=spos(2)
        shot%src%y=spos(3)

        !receiver side
        shot%nrcv=nr*1  !size(rcomp)
        if(allocated(shot%rcv)) deallocate(shot%rcv)
        allocate(shot%rcv(shot%nrcv))
        shot%rcv(:)%comp = rcomp
        shot%rcv(:)%y = 0.

        k=0

        !left receivers
        line=lline
        do i=k+1,k+int(line(4))
            shot%rcv(i)%x = line(1)
            shot%rcv(i)%z = line(2) + (i-k-1)*line(3)
        enddo
        k=k+line(4)

        !top receiver line
        line=tline
        do i=k+1,k+int(line(4))
            shot%rcv(i)%z = line(1)
            shot%rcv(i)%x = line(2) + (i-k-1)*line(3)
        enddo
        k=k+line(4)

        !right receiver line
        line=rline
        do i=k+1,k+int(line(4))
            shot%rcv(i)%x = line(1)
            shot%rcv(i)%z = line(2) + (i-k-1)*line(3)
        enddo
        k=k+line(4)

        !bottom receivers
        line=bline
        do i=k+1,k+int(line(4))
            shot%rcv(i)%z = line(1)
            shot%rcv(i)%x = line(2) + (i-k-1)*line(3)
        enddo
        k=k+line(4)

    end subroutine

    ! subroutine geometry_3Dspread
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     !  Quadrilateral Geometry  !
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     !
    !     !  IL direction ->
    !     ! XL  P1---------------P2
    !     ! dir -------------------
    !     ! |   -------------------
    !     ! v   P3---------------P4
    !     !
    !     !4 anchor points, P1-4, and number of points in IL and XL directions
    !     !should be read from setup (SOURCE_POSITION, RECEIVER_POSITION)
    !     !other points are interpolated from P1-4 and nil, nxl

    !     real :: fsx,fsy,fsz, frx,fry,frz
    !     real :: lsx,lsy,lsz, lrx,lry,lrz
        
    !     read(src_pos,*)  sz1,sx1,sy1, sz2,sx2,sy2, sz3,sx3,sy3, sz4,sx4,sy4, nsil, nsxl
    !     read(rcv_pos,*)  rz1,rx1,ry1, rz2,rx2,ry2, rz3,rx3,ry3, rz4,rx4,ry4, nril, nrxl
        
    !     if(.not.m%is_cubic) then
    !         sx3=sx1; sy3=sy1; sz3=sz1;
    !         sx4=sx2; sy4=sy2; sz4=sz2;  nsxl=1

    !         rx3=rx1; ry3=ry1; rz3=rz1;
    !         rx4=rx2; ry4=ry2; rz4=rz2;  nrxl=1
    !     endif

    !     ns=nsil*nsxl
    !     nr=nril*nrxl

    !     acqui%nsrc=ns*nscomp
    !     if(allocated(acqui%src))deallocate(acqui%src)
    !     allocate(acqui%src(acqui%nsrc))

    !     do ixl=1,nsxl

    !         fsx = ( sx1*(nsxl-ixl) + sx3*(ixl-1) )/(nsxl-1)
    !         fsy = ( sy1*(nsxl-ixl) + sy3*(ixl-1) )/(nsxl-1)
    !         fsz = ( sz1*(nsxl-ixl) + sz3*(ixl-1) )/(nsxl-1)

    !         lsx = ( sx2*(nsxl-ixl) + sx4*(ixl-1) )/(nsxl-1)
    !         lsy = ( sy2*(nsxl-ixl) + sy4*(ixl-1) )/(nsxl-1)
    !         lsz = ( sz2*(nsxl-ixl) + sz4*(ixl-1) )/(nsxl-1)

    !     do iil=1,nsil

    !         sx  = ( fsx*(nsil-iil) + lsx*(iil-1) )/(nsil-1)
    !         sy  = ( fsy*(nsil-iil) + lsy*(iil-1) )/(nsil-1)
    !         sz  = ( fsz*(nsil-iil) + lsz*(iil-1) )/(nsil-1)

    !         i   = iil + (ixl-1)*nsil

    !         acqui%src(i)%x=sx
    !         acqui%src(i)%y=sy
    !         acqui%src(i)%z=sz
    !         acqui%src(i)%comp=src_comp(1)
            

    !         acqui%src(i)%nrcv=nr*nrcomp
    !         if(allocated(acqui%src(i)%rcv))deallocate(acqui%src(i)%rcv)
    !         allocate(acqui%src(i)%rcv(acqui%src(i)%nrcv))
            
    !         do jxl=1,nrxl

    !             frx = ( rx1*(nrxl-jxl) + rx3*(jxl-1) )/(nrxl-1)
    !             fry = ( ry1*(nrxl-jxl) + ry3*(jxl-1) )/(nrxl-1)
    !             frz = ( rz1*(nrxl-jxl) + rz3*(jxl-1) )/(nrxl-1)

    !             lrx = ( rx2*(nrxl-jxl) + rx4*(jxl-1) )/(nrxl-1)
    !             lry = ( ry2*(nrxl-jxl) + ry4*(jxl-1) )/(nrxl-1)
    !             lrz = ( rz2*(nrxl-jxl) + rz4*(jxl-1) )/(nrxl-1)

    !         do jil=1,nril

    !             rx  = ( frx*(nril-jil) + lrx*(jil-1) )/(nril-1)
    !             ry  = ( fry*(nril-jil) + lry*(jil-1) )/(nril-1)
    !             rz  = ( frz*(nril-jil) + lrz*(jil-1) )/(nril-1)

    !             j   = jil + (jxl-1)*nril

    !             acqui%src(i)%rcv(j)%x=rx
    !             acqui%src(i)%rcv(j)%y=ry
    !             acqui%src(i)%rcv(j)%z=rz
    !             acqui%src(i)%rcv(j)%z=rcv_comp(1)
                
    !         enddo

    !         enddo

    !     enddo

    !     enddo

    !     !duplicate positions for other components
    !     do l=2,nrcomp
    !         do j=1,nr
    !             jj=j+(l-1)*nr
    !             acqui%src(1)%rcv(jj)%z=acqui%src(1)%rcv(j)%z
    !             acqui%src(1)%rcv(jj)%x=acqui%src(1)%rcv(j)%x
    !             acqui%src(1)%rcv(jj)%y=acqui%src(1)%rcv(j)%y
    !         enddo
    !     enddo

    !     do k=2,nscomp
    !         do i=1,ns
    !             ii=i+(k-1)*ns
    !             acqui%src(ii)%rcv(:)%z=acqui%src(i)%rcv(:)%z
    !             acqui%src(ii)%rcv(:)%x=acqui%src(i)%rcv(:)%x
    !             acqui%src(ii)%rcv(:)%y=acqui%src(i)%rcv(:)%y
    !         enddo
    !     enddo

    ! end subroutine

    subroutine write(self,file,data)
        class(t_shot) :: self
        character(*) :: file
        real,dimension(shot%nt,shot%nrcv) :: data
        
        character(:),allocatable,save :: format
        type(t_suformat) :: sudata
        logical,save :: is_first_in=.true.
        real,save :: scalel, scalco

        if(is_first_in) format=setup%get_str('DATA_FORMAT',o_default='su')

        if(format=='su') then

            call sudata%init(self%nt,self%nrcv,o_dt=self%dt,o_data=data)

            !assume all shots share same same scalel & scalco
            !or given by user
            if(is_first_in) then
                scalel=setup%get_real('SU_SCALEL',o_default=num2str(find_proper_scalel()))
                scalco=setup%get_real('SU_SCALCO',o_default=num2str(find_proper_scalco()))
            endif

            sudata%hdrs(:)%scalel=int(-scalel)
            sudata%hdrs(:)%scalco=int(-scalco)

            do i=1,self%nrcv
                sudata%hdrs(i)%tracl=i

                sudata%hdrs(i)%sdepth= (self%src%z+m%oz)    *scalel 
                sudata%hdrs(i)%gelev =-(self%rcv(i)%z+m%oz) *scalel

                sudata%hdrs(i)%sx=(self%src%x+m%ox)         *scalco
                sudata%hdrs(i)%sy=(self%src%y+m%oy)         *scalco
                sudata%hdrs(i)%gx=(self%rcv(i)%x+m%ox)      *scalco
                sudata%hdrs(i)%gy=(self%rcv(i)%y+m%oy)      *scalco

                select case (self%rcv(i)%comp)
                case ('p') !pressure
                    sudata%hdrs(i)%trid=11
                case ('vz') !vertical velocity
                    sudata%hdrs(i)%trid=12
                case ('vx') !horizontal velocity
                    sudata%hdrs(i)%trid=13
                case ('vy')
                    sudata%hdrs(i)%trid=14
                end select

                if(self%rcv(i)%is_badtrace) sudata%hdrs(i)%trid=3 !dummy trace

            enddo
            
            call sudata%write(file//self%sindex)

        endif

        if(format/='su') then !binary format
            call sysio_write(file//self%sindex,data,size(data))
            
        endif
        
        is_first_in=.false.

    end subroutine

    function find_proper_scalel() result(scalel)

        scalel=1
        
        do
            test=(shot%src%z+m%oz)*scalel
            if(floor(test) == ceiling(test)) exit
            scalel=scalel*10.

            if(scalel>1e4) then
                call warn('scalel > 1e4! Will have to lose some precision in the SU header.')
                return
            endif
        enddo

        do i=1,shot%nrcv
            do
                test=(shot%rcv(i)%z+m%oz)*scalel
                if(floor(test) == ceiling(test)) exit
                scalel=scalel*10.

                if(scalel>1e4) then
                    call warn('scalel > 1e4! Will have to lose some precision in the SU header.')
                    return
                endif
            enddo
        enddo

        call hud('Found proper SU scalel = '//num2str(scalel))

    end function

    function find_proper_scalco() result(scalco)

        scalco=1
        
        do
            test=(shot%src%x+m%ox)*scalco
            if(floor(test) == ceiling(test)) exit
            scalco=scalco*10.

            if(scalco>1e4) then
                call warn('scalco > 1e4! Will have to lose some precision in the SU header.')
                return
            endif
        enddo

        do
            test=(shot%src%y+m%oy)*scalco
            if(floor(test) == ceiling(test)) exit
            scalco=scalco*10.

            if(scalco>1e4) then
                call warn('scalco > 1e4! Will have to lose some precision in the SU header.')
                return
            endif
        enddo

        do i=1,shot%nrcv
            do
                test=(shot%rcv(i)%x+m%ox)*scalco
                if(floor(test) == ceiling(test)) exit
                scalco=scalco*10.

                if(scalco>1e4) then
                    call warn('scalco > 1e4! Will have to lose some precision in the SU header.')
                    return
                endif
            enddo
        enddo

        do i=1,shot%nrcv
            do
                test=(shot%rcv(i)%y+m%oy)*scalco
                if(floor(test) == ceiling(test)) exit
                scalco=scalco*10.

                if(scalco>1e4) then
                    call warn('scalco > 1e4! Will have to lose some precision in the SU header.')
                    return
                endif
            enddo
        enddo

        call hud('Found proper SU scalco = '//num2str(scalco))

    end function

    subroutine check_range(shot)
        type(t_shot) :: shot

        character(:),allocatable :: str

        !check source positions
        str='Shot# '//shot%sindex//' : '
        if(shot%src%z<0.) then
            shot%src%z = m%dz
            call hud(str//'sz above top bnd of model! Set sz = dz.',mpiworld%iproc)
        endif
        if(shot%src%z>(m%nz-1)*m%dz) then
            shot%src%z = (m%nz-2)*m%dz
            call hud(str//'sz below bottom bnd of model! Set sz = (nz-2)*dz.',mpiworld%iproc)
        endif
        if(shot%src%x<0.) then
            shot%src%x = m%dx
            call hud(str//'sx beyond left bnd of model! Set sx = dx.',mpiworld%iproc)
        endif
        if(shot%src%x>(m%nx-1)*m%dx) then
            shot%src%x = (m%nx-2)*m%dx
            call hud(str//'sx beyond right bnd of model! Set sx = (nx-2)*dx.',mpiworld%iproc)
        endif
        if(shot%src%y<0.) then
            shot%src%y = m%dy
            call hud(str//'sy beyond front bnd of model! Set sy = dy.',mpiworld%iproc)
        endif
        if(shot%src%y>(m%ny-1)*m%dy) then
            shot%src%y = (m%ny-2)*m%dy
            call hud(str//'sy beyond rear bnd of model! Set sy = (ny-2)*dy.',mpiworld%iproc)
        endif
        deallocate(str)

        !check receiver positions
        str='Shot# '//shot%sindex//' has receiver(s) '
        if(any(shot%rcv(:)%z<0.)) then
            where(shot%rcv(:)%z<0.) shot%rcv(:)%z=m%dz
            call hud(str//'rz above top of model! Set rz = dz.',mpiworld%iproc)
        endif
        if(any(shot%rcv(:)%z>(m%nz-1)*m%dz)) then
            where(shot%rcv(:)%z>(m%nz-1)*m%dz) shot%rcv(:)%z=(m%nz-2)*m%dz
            call hud(str//'rz below bottom of model! Set rz = (nz-2)*dz.',mpiworld%iproc)
        endif
        if(any(shot%rcv(:)%x<0.)) then
            where(shot%rcv(:)%x<0.) shot%rcv(:)%x=m%dx
            call hud(str//'rx beyond left bnd of model! Set rx = dx.',mpiworld%iproc)
        endif
        if(any(shot%rcv(:)%x>(m%nx-1)*m%dx)) then
            where(shot%rcv(:)%x>(m%nx-1)*m%dx) shot%rcv(:)%x=(m%nx-2)*m%dx
            call hud(str//'rx beyond right bnd of model! Set rx = (nx-2)*dx.',mpiworld%iproc)
        endif
        if(any(shot%rcv(:)%y<0.)) then
            where(shot%rcv(:)%y<0.) shot%rcv(:)%y=m%dy
            call hud(str//'ry outside front end of model! Set ry = dy.',mpiworld%iproc)
        endif
        if(any(shot%rcv(:)%y>(m%ny-1)*m%dy)) then
            where(shot%rcv(:)%y>(m%ny-1)*m%dy) shot%rcv(:)%y=(m%ny-2)*m%dy
            call hud(str//'ry beyond rear end of model! Set sy = (ny-2)*dy.',mpiworld%iproc)
        endif
    end subroutine

end

