module m_shot
use m_shotlist
use m_gen_acqui
use m_suformat
use m_hicks
use m_butterworth

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
    
    type t_source
        real    :: x,y,z
        integer :: ix,iy,iz
        integer :: ifz,ilz,ifx,ilx,ify,ily
        character(4) :: comp
        real,dimension(:,:,:),allocatable :: interp_coeff
    end type

    type t_receiver
        real    :: x,y,z, aoffset
        integer :: ix,iy,iz
        integer :: ifz,ilz,ifx,ilx,ify,ily
        character(4) :: comp
        real,dimension(:,:,:),allocatable :: interp_coeff
    end type
    
    type t_shot
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
        real,dimension(:,:),allocatable :: dres !residual seismogram
        !real,dimension(:,:),allocatable :: dadj !adjoint seismogram
        
    end type
    
    type(t_shot) :: shot
    
    contains
    
    subroutine init(ishot)
        integer ishot

        self%index=shotlist(ishot)
        self%sindex=num2str(shot%index,'(i0.4)')
    
    end subroutine

    subroutine read_from_setup
        logical,save :: first_in=.true.

        if(first_in) then
            shot%nt=setup%get_int('TIME_STEP','NT')
            shot%dt=setup%get_int('TIME_INTERVAL','DT')
            shot%src%comp=setup%get_str( 'SOURCE_COMPONENT',  'SCOMP',default='p'))
        endif

        first_in=.false.

        select case (setup%get_char('ACQUI_GEOMETRY',default='spread'))
        case ('spread')
            call geometry_spread(self)
        case ('streamer')
            call geometry_streamer(self)
        case ('irregularOBN')
        case ('spread3D')
        end select

        call check_setup_positions

    end subroutine

    subroutine read_from_data
        type(t_format_su) :: data
        
        file=setup%get_file('FILE_DATA',mandatory=.true.)

        call data%read(file//shot%sindex,shot%sindex)
        
        !source & receiver geometry
        shot%src%x=data%hdr(1)%sx
        shot%src%y=data%hdr(1)%sy
        shot%src%z=data%hdr(1)%sdepth
        
        shot%src%icomp=1 !don't know which su header tells this info..
        
        shot%nt=data%ns  !assume all traces have same ns
        shot%dt=data%hdr(1)%dt*1e-6 !assume all traces have same dt
        
        shot%nrcv=data%ntr
        if(allocated(shot%rcv))deallocate(shot%rcv)
        allocate(shot%rcv(shot%nrcv))

        do ir=1,shot%nrcv
            shot%rcv(ir)%x= data%hdr(ir)%gx
            shot%rcv(ir)%y= data%hdr(ir)%gy
            shot%rcv(ir)%z=-data%hdr(ir)%gelev
            
            select case (data%hdr(ir)%trid)
                case (11)
                shot%rcv(ir)%comp='p'  !pressure
                case (12)
                shot%rcv(ir)%comp='vz' !vertical velocity
                case (13)
                shot%rcv(ir)%comp='vy' !horizontal velocity
                case (14)
                shot%rcv(ir)%comp='vx'
                
                case default
                shot%rcv(ir)%comp='p'
            end select
            
        enddo
        

        !scale back elevation
        if(data%hdr(1)%scalel > 0) then
            shot%src%z    = shot%src%z    * data(1)%hdr%scalel
            shot%rcv(:)%z = shot%rcv(:)%z * data(1)%hdr%scalel
        elseif(data%hdr(1)%scalel < 0) then
            shot%src%z    = shot%src%z    / (-data(1)%hdr%scalel)
            shot%rcv(:)%z = shot%rcv(:)%z / (-data(1)%hdr%scalel)
        endif
        
        !scale back coordinates
        if(data%hdr(1)%scalco > 0) then
            shot%src%x    = shot%src%x    * data(1)%hdr%scalco
            shot%src%y    = shot%src%y    * data(1)%hdr%scalco
            shot%rcv(:)%x = shot%rcv(:)%x * data(1)%hdr%scalco
            shot%rcv(:)%y = shot%rcv(:)%y * data(1)%hdr%scalco
        elseif(data%hdr(1)%scalco < 0) then
            shot%src%x    = shot%src%x    / (-data(1)%hdr%scalco)
            shot%src%y    = shot%src%y    / (-data(1)%hdr%scalco)
            shot%rcv(:)%x = shot%rcv(:)%x / (-data(1)%hdr%scalco)
            shot%rcv(:)%y = shot%rcv(:)%y / (-data(1)%hdr%scalco)
        endif
        
        !load obs traces
        call alloc(self%dobs,shot%nt,shot%nrcv)
        self%dobs=data%trace
        
        call check_setup_positions

    end subroutine

    subroutine set_time

        type(t_format_su)::wavelet

        shot%fpeak=setup%get_real('PEAK_FREQUENCY','FPEAK')
        shot%fmax=setup%get_real('MAX_FREQUENCY','FMAX',default=shot%fpeak*2.5)

        !wavelet
        file=setup%get_file('FILE_SOURCE_WAVELET')

        if(file=='') then !not given
            if(setup%get_str('WAVELET_TYPE',default='sinexp')=='sinexp') then
                call hud('Use filtered sinexp wavelet')
                shot%wavelet=sinexp(shot%nt,shot%dt,shot%fpeak)
            else
                call hud('Use Ricker wavelet')
                shot%wavelet=ricker(shot%nt,shot%dt,shot%fpeak)
            endif
            shot%wavelet=shot%wavelet *shot%dt/m%cell_volume

        else !wavelet file exists
            call wavelet%read(file)
            call resample(wavelet%trs(:,1),shot%wavelet,shot%nt,shot%dt)

        endif

        if(mpiworld%is_master) then
            open(12,file='source_wavelet',access='direct',recl=4*shot%nt)
            write(12,rec=1) shot%wavelet
            close(12)
        endif

    end subroutine

    subroutine update_wavelet

        call matchfilter_estimate(shot%nt,shot%nrcv,shot%dsyn,shot%dobs,shot%sindex)
        
        call matchfilter_apply_to_wavelet(shot%nt,shot%wavelet)
        
        call matchfilter_apply_to_data(shot%nt,shot%nrcv,shot%dsyn)

        !call suformat_write(iproc=0,file='wavelet_update',append=.true.)        
        if(mpiworld%is_master) then
            open(12,file='source_wavelet',access='stream',recl=4*shot%nt,position='append')
            write(12) shot%wavelet
            close(12)
        endif

    end subroutine
    
    subroutine update_residuals

        call matchfilter_correlate_filter_residual(shot%nt,shot%nrcv,shot%dres)

    end subroutine

    subroutine acquire(self)
        type(t_shot) :: self

        call alloc(self%dsyn,self%nt,self%nrcv)
        do i=1,self%nrcv
            select case (self%rcv(i)%comp)
            case ('p')
                call resample(f%seismo%p(i,:), self%dsyn(:,i))
            case ('vx')
                call resample(f%seismo%vx(i,:),self%dsyn(:,i))
            case ('vy')
                call resample(f%seismo%vy(i,:),self%dsyn(:,i))
            case ('vz')
                call resample(f%seismo%vz(i,:),self%dsyn(:,i))
            end select
        enddo

    end subroutine


    subroutine set_space

        !shift positions to be 0-based
        shot%src%x   =shot%src%x    - m%ox
        shot%src%y   =shot%src%y    - m%oy
        shot%src%z   =shot%src%z    - m%oz
        shot%rcv(:)%x=shot%rcv(:)%x - m%ox
        shot%rcv(:)%y=shot%rcv(:)%y - m%oy
        shot%rcv(:)%z=shot%rcv(:)%z - m%oz

        !absolute offset
        do ir=1,shot%nrcv
            shot%rcv(ir)%aoffset=sqrt( (shot%src%x-shot%rcv(ir)%x)**2 &
                                      +(shot%src%y-shot%rcv(ir)%y)**2 &
                                      +(shot%src%z-shot%rcv(ir)%z)**2 )
        enddo

        !shift velocity components in case of staggered grid
        if_fdsg=index(projector%info,'FDSG')>0

        if(if_fdsg) then
            !source side
            if(self%src%comp(1:1)=='v') then
                shot%src%x=shot%src%x-(1-1)*m%dx !?
                shot%src%y=shot%src%y-(1-1)*m%dy
                shot%src%z=shot%src%z-(1-1)*m%dz
            endif
            
            !receiver side
            do i=1,shot%nrcv
                if(self%rcv%comp(1:1)=='v') then
                    shot%rcv(i)%x=shot%rcv(i)%x-(1-1)*m%dx !?
                    shot%rcv(i)%y=shot%rcv(i)%y-(1-1)*m%dy
                    shot%rcv(i)%z=shot%rcv(i)%z-(1-1)*m%dz
                endif
            enddo
            
            call hud('For velocity components, src & rcv positions have been shifted one grid point due to staggered grid stencils.')
        endif

        self%if_hicks=setup_get_logical('IF_HICKS',default=.true.)
        
        if(self%if_hicks) then

            !Hicks interpolation
            call hicks_init([m%dx,m%dy,m%dz],m%is_cubic, m%if_freesurface)
            
            !hicks coeff for source point
            if(if_fdsg) then
                !for vx,vy,vz components, hicks%x,y,z should be added by m%dx/2,dy/2,dz/2, respectively,
                !because v(1) is actually v[0.5], v(2) is v[1.5] etc.
                if(shot%src%comp=='vx') call hicks_init_position([shot%src%x+m%dx/2.,shot%src%y,        shot%src%z        ])
                if(shot%src%comp=='vy') call hicks_init_position([shot%src%x        ,shot%src%y+m%dy/2.,shot%src%z        ])
                if(shot%src%comp=='vz') call hicks_init_position([shot%src%x        ,shot%src%y,        shot%src%z+m%dz/2.])
            else
                call hicks_init_position([shot%src%x,shot%src%y,shot%src%z])
            endif

            if(shot%src%comp(1)=='p') then !explosive source or non-vertical force
                call hicks_build_coeff('antisymm', shot%src%interp_coeff)
            elseif(shot%src%comp=='vz') then !vertical force
                shot%src%interp_coeff=hicks%build_coeff('symmetric')
                call hicks_build_coeff('symmetric',shot%src%interp_coeff)
            else
                shot%src%interp_coeff=hicks%build_coeff('truncate')
                call hicks_build_coeff('truncate', shot%src%interp_coeff)
            endif

            call hicks_get_position([shot%src%ifx, shot%src%ify, shot%src%ifz], &
                                    [shot%src%ix,  shot%src%iy,  shot%src%iz ], &
                                    [shot%src%ilx, shot%src%ily, shot%src%ilz]  )

            !hicks coeff for receivers
            do i=1,shot%nrcv

                if(if_fdsg) then
                    !for vx,vy,vz components, hicks%x,y,z should be added by m%dx/2,dy/2,dz/2, respectively,
                    !because v(1) is actually v[0.5], v(2) is v[1.5] etc.
                    if(shot%rcv(i)%comp=='vx') call hicks_init_position([shot%rcv(i)%x+m%dx/2.,shot%rcv(i)%y,        shot%rcv(i)%        ])
                    if(shot%rcv(i)%comp=='vy') call hicks_init_position([shot%rcv(i)%x        ,shot%rcv(i)%y+m%dy/2.,shot%rcv(i)%        ])
                    if(shot%rcv(i)%comp=='vz') call hicks_init_position([shot%rcv(i)%x        ,shot%rcv(i)%y,        shot%rcv(i)%+m%dz/2.])
                else
                    call hicks_init_position([shot%rcv(i)%,shot%rcv(i)%y,shot%rcv(i)%z])
                endif

                if(shot%rcv(i)%comp(1)=='p') then !explosive source or non-vertical force
                    call hicks_build_coeff('antisymm', shot%rcv(i)%interp_coeff)
                elseif(shot%src%comp=='vz') then !vertical force
                    call hicks_build_coeff('symmetric',shot%rcv(i)%interp_coeff)
                else
                    call hicks_build_coeff('truncate', shot%rcv(i)%interp_coeff)
                endif

                call hicks_get_position([shot%rcv(i)%ifx, shot%rcv(i)%ify, shot%rcv(i)%ifz], &
                                        [shot%rcv(i)%ix,  shot%rcv(i)%iy,  shot%rcv(i)%iz ], &
                                        [shot%rcv(i)%ilx, shot%rcv(i)%ily, shot%rcv(i)%ilz]  )

            enddo

        endif
                
        if(mpiworld%is_master) then
            write(*,*)'================================='
            write(*,*)'Shot# '//shot%sindex//' info:'
            write(*,*)'================================='
            write(*,*)'  nt,dt:',shot%src%nt,shot%src%dt
            write(*,*)'---------------------------------'
            write(*,*)'  sz,isz:',shot%src%z,shot%src%iz
            write(*,*)'  sx,isx:',shot%src%x,shot%src%ix
            write(*,*)'  sy,isy:',shot%src%y,shot%src%iy
            write(*,*)'  ifz,ilz:',shot%src%ifz,shot%src%ilz
            write(*,*)'  ifx,ilx:',shot%src%ifx,shot%src%ilx
            write(*,*)'  ify,ily:',shot%src%ify,shot%src%ily
            write(*,*)'---------------------------------'
            write(*,*)'  minmax rz,irz:',minval(shot%rcv(:)%z),maxval(shot%rcv(:)%z),minval(shot%rcv(:)%iz),maxval(shot%rcv(:)%iz)
            write(*,*)'  minmax rx,irx:',minval(shot%rcv(:)%x),maxval(shot%rcv(:)%x),minval(shot%rcv(:)%ix),maxval(shot%rcv(:)%ix)
            write(*,*)'  minmax ry,iry:',minval(shot%rcv(:)%y),maxval(shot%rcv(:)%y),minval(shot%rcv(:)%iy),maxval(shot%rcv(:)%iy)
            write(*,*)'  minmax ifz,ilz:',minval(shot%rcv(:)%ifz),maxval(shot%rcv(:)%ifz),minval(shot%rcv(:)%ilz),maxval(shot%rcv(:)%ilz)
            write(*,*)'  minmax ifx,ilx:',minval(shot%rcv(:)%ifx),maxval(shot%rcv(:)%ifx),minval(shot%rcv(:)%ilx),maxval(shot%rcv(:)%ilx)
            write(*,*)'  minmax ify,ily:',minval(shot%rcv(:)%ify),maxval(shot%rcv(:)%ify),minval(shot%rcv(:)%ily),maxval(shot%rcv(:)%ily)
            write(*,*)'  nrcv:',shot%nrcv
            write(*,*)'---------------------------------'
        endif

    end subroutine
        
    subroutine geometry_spread(shot)
        type(t_shot) :: shot
        logical :: first_in=.true.
        real,dimension(3),save :: os,or, ds,dr
        real,save :: nr
        type(t_string),dimension(:),allocatable,save :: rcomp

        if(first_in) then
            os=setup%get_reals('SOURCE_ORIGIN',  'OS')
            or=setup%get_reals('RECEIVER_ORIGIN','OR')

            ds=setup%get_reals('SOURCE_SPACING', 'DS')
            dr=setup%get_reals('RECEIVER_SPACING','DR')

            nr=setup%get_int('NUMBER_RECEIVER','NR')
            rcomp=setup%get_strs('RECEIVER_COMPONENT','RCOMP',default='p'))
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

        if(allocated(self%rcv)) deallocate(self%rcv)
        allocate(self%rcv(shot%nrcv))

        do j=1,size(rcomp)
            do i=0,nr-1
                self%rcv(1+i+j*nr)%comp = rcomp(j)%s
                self%rcv(1+i+j*nr)%rz = or(1) + (i-1)*dr(1)
                self%rcv(1+i+j*nr)%rx = or(2) + (i-1)*dr(2)
                self%rcv(1+i+j*nr)%ry = or(3) + (i-1)*dr(3)
            enddo
        enddo

    end subroutine

    subroutine geometry_streamer(shot)
        type(t_shot) :: shot
        logical :: first_in=.true.
        real,dimension(3),save :: os,ooff, ds,doff
        real,save :: noff
        type(t_string),dimension(:),allocatable,save :: rcomp

        if(first_in) then
            os  =setup%get_reals('SOURCE_ORIGIN','OS')
            ooff=setup%get_reals('OFFSET_ORIGIN','OOFF')

            ds  =setup%get_reals('SOURCE_SPACING', 'DS')
            doff=setup%get_reals('OFFSET_SPACING', 'DOFF')

            noff=setup%get_int('NUMBER_OFFSET','NOFF')
            rcomp=setup%get_strs('RECEIVER_COMPONENT','RCOMP',default='p'))
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

        if(allocated(self%rcv)) deallocate(self%rcv)
        allocate(self%rcv(shot%nrcv))

        do j=1,size(rcomp)
            do i=0,nr-1
                self%rcv(1+i+j*nr)%comp = rcomp(j)%s
                self%rcv(1+i+j*nr)%rz = self%src%z + off(1) + (i-1)*doff(1)
                self%rcv(1+i+j*nr)%rx = self%src%x + off(2) + (i-1)*doff(2)
                self%rcv(1+i+j*nr)%ry = self%src%y + off(3) + (i-1)*doff(3)
            enddo
        enddo

    end subroutine

    subroutine geometry_irregularOBN(shot)
        type(t_shot) :: shot
        logical :: first_in=.true.
        character(:),allocatable,save :: spos, rpos
        integer, save :: nr
        type(t_string),dimension(:),allocatable,save :: rcomp

        if(first_in) then            
            spos=setup%get_file('FILE_SOURCE_POSITION','SPOS',mandatory=.true.)
            rpos=setup%get_file('FILE_RECEIVER_POSITION','RPOS',mandatory=.true.)

            rcomp=setup%get_strs('RECEIVER_COMPONENT','RCOMP',default='p'))
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
        open(15,file=rcv_pos,action='read')
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

        first_in=.false.

    end subroutine

    ! subroutine gen_acqui_3Dspread
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

    subroutine check_setup_positions(shot)
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
            where(shot%rcv(:)<0.) shot%rcv(:)=m%dz
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
