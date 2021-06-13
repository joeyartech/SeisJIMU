module m_shot
use m_string
use m_mpienv
use m_setup
use m_format_su
use m_hicks
use m_resampler
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
        procedure :: init => init
        procedure :: read_from_setup => read_from_setup
        procedure :: read_from_data  => read_from_data
        procedure :: set_var_time  => set_var_time
        procedure :: set_var_space => set_var_space
        procedure :: update_wavelet => update_wavelet
        procedure :: update_adjsource => update_adjsource
        
    end type
    
    type(t_shot),public :: shot
    
    contains
    
    subroutine init(self,index)
        class(t_shot) :: self
        integer index

        self%index=index
        self%sindex=num2str(index,'(i0.4)')
    
    end subroutine

    subroutine read_from_setup(self)
        class(t_shot) :: self

        logical,save :: first_in=.true.
        type(t_string),dimension(:),allocatable :: scomp

        if(first_in) then
            self%nt=setup%get_int('TIME_STEP','NT')
            self%dt=setup%get_real('TIME_INTERVAL','DT')
            
            scomp=setup%get_strs('SOURCE_COMPONENT','SCOMP',o_default='p')
            if(size(scomp)>1) call hud('SeisJIMU only considers the 1st component from SOURCE_COMPONENT.')
            self%src%comp=scomp(1)%s
        endif

        first_in=.false.

        select case (setup%get_str('ACQUI_GEOMETRY',o_default='spread'))
        case ('spread')
            call geometry_spread(self)
        case ('streamer')
            call geometry_streamer(self)
        case ('spread_irregular')
            call geometry_spread_irregular(self)
        case ('spread3D')
        end select

        call check_range(self)

    end subroutine

    subroutine read_from_data(self)
        class(t_shot) :: self

        type(t_format_su) :: data
        
        call data%read(setup%get_file('FILE_DATA',o_mandatory=.true.),self%sindex)

        !source & receiver geometry
        self%src%x=data%hdrs(1)%sx
        self%src%y=data%hdrs(1)%sy
        self%src%z=data%hdrs(1)%sdepth
        
        self%src%comp='p' !don't know which su header tells this info..
        
        self%nt=data%ns  !assume all traces have same ns
        self%dt=data%dt  !assume all traces have same dt
        
        self%nrcv=data%ntr
        if(allocated(self%rcv))deallocate(self%rcv)
        allocate(self%rcv(self%nrcv))

        do ir=1,self%nrcv
            self%rcv(ir)%x= data%hdrs(ir)%gx
            self%rcv(ir)%y= data%hdrs(ir)%gy
            self%rcv(ir)%z=-data%hdrs(ir)%gelev
            
            select case (data%hdrs(ir)%trid)
            case (11)
                self%rcv(ir)%comp='p'  !pressure
            case (12)
                self%rcv(ir)%comp='vz' !vertical velocity
            case (13)
                self%rcv(ir)%comp='vy' !horizontal velocity
            case (14)
                self%rcv(ir)%comp='vx'
                
            case default
                self%rcv(ir)%comp='p'
            end select
            
        enddo
        

        !scale back elevation
        if(data%hdrs(1)%scalel > 0) then
            self%src%z    = self%src%z    * data%hdrs(1)%scalel
            self%rcv(:)%z = self%rcv(:)%z * data%hdrs(1)%scalel
        elseif(data%hdrs(1)%scalel < 0) then
            self%src%z    = self%src%z    / (-data%hdrs(1)%scalel)
            self%rcv(:)%z = self%rcv(:)%z / (-data%hdrs(1)%scalel)
        endif
        
        !scale back coordinates
        if(data%hdrs(1)%scalco > 0) then
            self%src%x    = self%src%x    * data%hdrs(1)%scalco
            self%src%y    = self%src%y    * data%hdrs(1)%scalco
            self%rcv(:)%x = self%rcv(:)%x * data%hdrs(1)%scalco
            self%rcv(:)%y = self%rcv(:)%y * data%hdrs(1)%scalco
        elseif(data%hdrs(1)%scalco < 0) then
            self%src%x    = self%src%x    / (-data%hdrs(1)%scalco)
            self%src%y    = self%src%y    / (-data%hdrs(1)%scalco)
            self%rcv(:)%x = self%rcv(:)%x / (-data%hdrs(1)%scalco)
            self%rcv(:)%y = self%rcv(:)%y / (-data%hdrs(1)%scalco)
        endif
        
        !load obs traces
        call alloc(self%dobs,self%nt,self%nrcv)
        self%dobs=data%trs
        
        call check_setup_positions

    end subroutine

    subroutine set_var_time(self)
        class(t_shot) :: self

        character(:),allocatable :: file
        type(t_format_su)::wavelet

        self%fpeak=setup%get_real('PEAK_FREQUENCY','FPEAK')
        self%fmax=setup%get_real('MAX_FREQUENCY','FMAX',o_default=num2str(self%fpeak*2.5))
        !self%fmin=setup%get_real('MIN_FREQUENCY','FMIN',o_default='1')

        !wavelet
        file=setup%get_file('FILE_SOURCE_WAVELET')

        if(file=='') then !not given
            if(setup%get_str('WAVELET_TYPE',o_default='sinexp')=='sinexp') then
                call hud('Use filtered sinexp wavelet')
                self%wavelet=sinexp(self%nt,self%dt,self%fpeak)
            else
                call hud('Use Ricker wavelet')
                self%wavelet=ricker(self%nt,self%dt,self%fpeak)
            endif

        else !wavelet file exists
            call wavelet%read(file)
            call resampler(wavelet%trs(:,1),self%wavelet,1,&
                            din=wavelet%dt,nin=wavelet%ns, &
                            dout=self%dt,  nout=self%nt)

        endif

        if(mpiworld%is_master) then
            open(12,file='source_wavelet',access='direct',recl=4*self%nt)
            write(12,rec=1) self%wavelet
            close(12)
        endif

    end subroutine

    subroutine set_var_space(self,is_fdsg)
        class(t_shot) :: self
        logical :: is_fdsg

        !shift positions to be 0-based
        !then m%oz,ox,oy is no longer of use
        self%src%z   =self%src%z    - m%oz
        self%src%x   =self%src%x    - m%ox
        self%src%y   =self%src%y    - m%oy
        self%rcv(:)%z=self%rcv(:)%z - m%oz
        self%rcv(:)%x=self%rcv(:)%x - m%ox
        self%rcv(:)%y=self%rcv(:)%y - m%oy

        !absolute offset
        do ir=1,self%nrcv
            self%rcv(ir)%aoffset=sqrt( (self%src%z-self%rcv(ir)%z)**2 &
                                      +(self%src%x-self%rcv(ir)%x)**2 &
                                      +(self%src%y-self%rcv(ir)%y)**2 )
        enddo

        self%if_hicks=setup%get_bool('IF_HICKS',o_default='T')
        
        if(self%if_hicks) then

            !Hicks interpolation
            call hicks_init(m%dz,m%dx,m%dy,m%is_cubic,m%is_freesurface)
            
            !hicks coeff for source point
            if(is_fdsg) then
                !for vx,vy,vz components, hicks%x,y,z should be added by m%dx/2,dy/2,dz/2, respectively,
                !because v(1) is actually v[0.5], v(2) is v[1.5] etc.
                if(self%src%comp=='vz') call hicks_init_position(self%src%z+m%dz/2.,self%src%x        ,self%src%y        )
                if(self%src%comp=='vx') call hicks_init_position(self%src%z        ,self%src%x+m%dx/2.,self%src%y        )
                if(self%src%comp=='vy') call hicks_init_position(self%src%z        ,self%src%x        ,self%src%y+m%dy/2.)
            else
                call hicks_init_position(self%src%z,self%src%x,self%src%y)
            endif

            if(self%src%comp=='p') then !explosive source or non-vertical force
                call hicks_get_coef('antisymm', self%src%interp_coef)
            elseif(self%src%comp=='vz') then !vertical force
                call hicks_get_coef('symmetric',self%src%interp_coef)
            else
                call hicks_get_coef('truncate', self%src%interp_coef)
            endif

            call hicks_get_position(self%src%ifz, self%src%ifx, self%src%ify, &
                                    self%src%iz,  self%src%ix,  self%src%iy , &
                                    self%src%ilz, self%src%ilx, self%src%ily  )

            !hicks coeff for receivers
            do i=1,self%nrcv

                if(is_fdsg) then
                    !for vx,vy,vz components, hicks%x,y,z should be added by m%dx/2,dy/2,dz/2, respectively,
                    !because v(1) is actually v[0.5], v(2) is v[1.5] etc.
                    if(self%rcv(i)%comp=='vz') call hicks_init_position(self%rcv(i)%z+m%dz/2.,self%rcv(i)%x        ,self%rcv(i)%y        )
                    if(self%rcv(i)%comp=='vx') call hicks_init_position(self%rcv(i)%z        ,self%rcv(i)%x+m%dx/2.,self%rcv(i)%y        )
                    if(self%rcv(i)%comp=='vy') call hicks_init_position(self%rcv(i)%z        ,self%rcv(i)%x        ,self%rcv(i)%y+m%dy/2.)
                else
                    call hicks_init_position(self%rcv(i)%z,self%rcv(i)%x,self%rcv(i)%y)
                endif

                if(self%rcv(i)%comp=='p') then !explosive source or non-vertical force
                    call hicks_get_coef('antisymm', self%rcv(i)%interp_coef)
                elseif(self%rcv(i)%comp=='vz') then !vertical force
                    call hicks_get_coef('symmetric',self%rcv(i)%interp_coef)
                else
                    call hicks_get_coef('truncate', self%rcv(i)%interp_coef)
                endif

                call hicks_get_position(self%rcv(i)%ifz, self%rcv(i)%ifx, self%rcv(i)%ify, &
                                        self%rcv(i)%iz , self%rcv(i)%ix , self%rcv(i)%iy , &
                                        self%rcv(i)%ilz, self%rcv(i)%ilx, self%rcv(i)%ily  )

            enddo

        endif
                
        if(mpiworld%is_master) then
            write(*,*)'================================='
            write(*,*)'Shot# '//self%sindex//' info:'
            write(*,*)'================================='
            write(*,*)'  nt,dt:',self%nt,self%dt
            write(*,*)'---------------------------------'
            write(*,*)'  sz,isz:',self%src%z,self%src%iz
            write(*,*)'  sx,isx:',self%src%x,self%src%ix
            write(*,*)'  sy,isy:',self%src%y,self%src%iy
            write(*,*)'  ifz,ilz:',self%src%ifz,self%src%ilz
            write(*,*)'  ifx,ilx:',self%src%ifx,self%src%ilx
            write(*,*)'  ify,ily:',self%src%ify,self%src%ily
            write(*,*)'---------------------------------'
            write(*,*)'  minmax rz,irz:',minval(self%rcv(:)%z),maxval(self%rcv(:)%z),minval(self%rcv(:)%iz),maxval(self%rcv(:)%iz)
            write(*,*)'  minmax rx,irx:',minval(self%rcv(:)%x),maxval(self%rcv(:)%x),minval(self%rcv(:)%ix),maxval(self%rcv(:)%ix)
            write(*,*)'  minmax ry,iry:',minval(self%rcv(:)%y),maxval(self%rcv(:)%y),minval(self%rcv(:)%iy),maxval(self%rcv(:)%iy)
            write(*,*)'  minmax ifz,ilz:',minval(self%rcv(:)%ifz),maxval(self%rcv(:)%ifz),minval(self%rcv(:)%ilz),maxval(self%rcv(:)%ilz)
            write(*,*)'  minmax ifx,ilx:',minval(self%rcv(:)%ifx),maxval(self%rcv(:)%ifx),minval(self%rcv(:)%ilx),maxval(self%rcv(:)%ilx)
            write(*,*)'  minmax ify,ily:',minval(self%rcv(:)%ify),maxval(self%rcv(:)%ify),minval(self%rcv(:)%ily),maxval(self%rcv(:)%ily)
            write(*,*)'  nrcv:',self%nrcv
            write(*,*)'---------------------------------'
        endif

    end subroutine

    subroutine update_wavelet(self)
        class(t_shot) :: self

        call matchfilter_estimate(self%nt,self%nrcv,self%dsyn,self%dobs,self%sindex)
        
        call matchfilter_apply_to_wavelet(self%nt,self%wavelet)
        
        call matchfilter_apply_to_data(self%nt,self%nrcv,self%dsyn)

        !call suformat_write(iproc=0,file='wavelet_update',append=.true.)        
        if(mpiworld%is_master) then
            open(12,file='source_wavelet',access='stream',position='append')
            write(12) self%wavelet
            close(12)
        endif

    end subroutine
    
    subroutine update_adjsource(self)
        class(t_shot) :: self

        call matchfilter_correlate_filter_residual(self%nt,self%nrcv,self%dadj)

    end subroutine

    subroutine geometry_spread(shot)
        type(t_shot) :: shot
        logical :: first_in=.true.
        real,dimension(3),save :: os,or, ds,dr
        integer,save :: nr
        type(t_string),dimension(:),allocatable,save :: rcomp

        if(first_in) then
            os=setup%get_reals('SOURCE_ORIGIN',  'OS')
            or=setup%get_reals('RECEIVER_ORIGIN','OR')

            ds=setup%get_reals('SOURCE_SPACING', 'DS')
            dr=setup%get_reals('RECEIVER_SPACING','DR')

            nr=setup%get_int('NUMBER_RECEIVER','NR')
            rcomp=setup%get_strs('RECEIVER_COMPONENT','RCOMP',o_default='p')
        endif

        first_in=.false.

        if(.not.m%is_cubic) then
            os(3)=0.; or(3)=0.
            ds(3)=0.; dr(3)=0.
        endif

        !source side
        shot%src%z=os(1)+(shot%index-1)*ds(1)
        shot%src%x=os(2)+(shot%index-1)*ds(2)
        shot%src%y=os(3)+(shot%index-1)*ds(3)

        !receiver side
        shot%nrcv=nr*size(rcomp)

        if(allocated(shot%rcv)) deallocate(shot%rcv)
        allocate(shot%rcv(shot%nrcv))

        do j=1,size(rcomp)
            do i=0,nr-1
                shot%rcv(1+i+j*nr)%comp = rcomp(j)%s
                shot%rcv(1+i+j*nr)%z = or(1) + (i-1)*dr(1)
                shot%rcv(1+i+j*nr)%x = or(2) + (i-1)*dr(2)
                shot%rcv(1+i+j*nr)%y = or(3) + (i-1)*dr(3)
            enddo
        enddo

    end subroutine

    subroutine geometry_streamer(shot)
        type(t_shot) :: shot
        logical :: first_in=.true.
        real,dimension(3),save :: os,ooff, ds,doff
        integer,save :: noff
        type(t_string),dimension(:),allocatable,save :: rcomp

        if(first_in) then
            os  =setup%get_reals('SOURCE_ORIGIN','OS')
            ooff=setup%get_reals('OFFSET_ORIGIN','OOFF')

            ds  =setup%get_reals('SOURCE_SPACING', 'DS')
            doff=setup%get_reals('OFFSET_SPACING', 'DOFF')

            noff=setup%get_int('NUMBER_OFFSET','NOFF')
            rcomp=setup%get_strs('RECEIVER_COMPONENT','RCOMP',o_default='p')
        endif

        first_in=.false.

        if(.not.m%is_cubic) then
            os(3)=0.; ooff(3)=0.
            ds(3)=0.; doff(3)=0.
        endif

        !source side
        shot%src%z=os(1)+(shot%index-1)*ds(1)
        shot%src%x=os(2)+(shot%index-1)*ds(2)
        shot%src%y=os(3)+(shot%index-1)*ds(3)

        !receiver side
        shot%nrcv=nr*size(rcomp)

        if(allocated(shot%rcv)) deallocate(shot%rcv)
        allocate(shot%rcv(shot%nrcv))

        do j=1,size(rcomp)
            do i=0,nr-1
                shot%rcv(1+i+j*nr)%comp = rcomp(j)%s
                shot%rcv(1+i+j*nr)%z = shot%src%z + off(1) + (i-1)*doff(1)
                shot%rcv(1+i+j*nr)%x = shot%src%x + off(2) + (i-1)*doff(2)
                shot%rcv(1+i+j*nr)%y = shot%src%y + off(3) + (i-1)*doff(3)
            enddo
        enddo

    end subroutine

    subroutine geometry_spread_irregular(shot)
        type(t_shot) :: shot
        logical :: first_in=.true.
        character(:),allocatable,save :: spos, rpos
        integer, save :: nr
        type(t_string),dimension(:),allocatable,save :: rcomp

        if(first_in) then            
            spos=setup%get_file('FILE_SOURCE_POSITION','SPOS',o_mandatory=.true.)
            rpos=setup%get_file('FILE_RECEIVER_POSITION','RPOS',o_mandatory=.true.)

            rcomp=setup%get_strs('RECEIVER_COMPONENT','RCOMP',o_default='p')
        endif

        first_in=.false.

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

    ! end subroutin

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

        !check receiver positions
        str='Shot# '//shot%sindex//' has receiver(s) '
        if(any(shot%rcv(:)%z<0.)) then
            where(shot%rcv(:)%z<0.) shot%rcv(:)%z=m%dz
            call hud(str//'rz above top of model! Set rz = dz.',mpiworld%iproc)
        endif
        if(any(shot%rcv(:)%z>(m%nz-1)*m%dz)) then
            where(shot%rcv(:)%z>(m%nz-1)*m%dz) shot%rcv(:)%z=(m%nz-2)*dz
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
