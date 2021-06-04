module m_shot
use m_shotlist
use m_gen_acqui
use m_suformat
use m_hicks
use m_butterworth

    private t_source, t_receiver, t_shot, file_wavelet
    
    type t_source
        real    :: x,y,z
        integer :: ix,iy,iz
        integer :: ifz,ilz,ifx,ilx,ify,ily
        character(4) :: comp
        real,dimension(:,:,:),allocatable :: interp_coeff
        integer :: nt
        real :: dt, fpeak
    end type

    type t_receiver
        real    :: x,y,z, aoffset
        integer :: ix,iy,iz
        integer :: ifz,ilz,ifx,ilx,ify,ily
        character(4) :: comp
        real,dimension(:,:,:),allocatable :: interp_coeff
        integer :: nt
        real :: dt
    end type
    
    type t_shot
        integer :: index
        character(4) :: sindex
        type(t_source) :: src
        type(t_receiver),dimension(:),allocatable :: rcv
        integer :: nrcv   !=size(rcv)

        real,dimension(:),allocatable :: wavelet
        character(:),allocatable :: file_wavelet

        real,dimension(:,:),allocatable :: dobs !observed seismogram
        real,dimension(:,:),allocatable :: dsyn !synthetic seismogram
        real,dimension(:,:),allocatable :: dres !residual seismogram
        !real,dimension(:,:),allocatable :: dadj !adjoint seismogram
        
    end type
    
    type(t_shot) :: shot

    logical :: if_staggered_grid=.false.
    
    contains
    
    subroutine shot_init(ishot,from)
        integer ishot
        character(*) :: from
        
        character(:),allocatable :: scale_wavelet

        shot%index=shotlist(ishot)
        shot%sindex=num2str(shot%index,'(i0.4)')
        
        !read geometry, nt and dt
        if(from=='gen_acqui') then
            call from_gen_acqui
        elseif(from=='data') then
            call from_data
        endif

        
        !! time !!

        shot%src%fpeak=setup_get_real('PEAK_FREQUENCY','FPEAK')
        
        !read wavelet
        shot%file_wavelet=setup_get_file('FILE_WAVELET')
        if(file_wavelet=='') then
            if(setup_get_char('WAVELET_TYPE',default='sinexp')=='sinexp') then
                call hud('Use filtered sinexp wavelet')
                call source_wavelet_sinexp
            else
                call hud('Use Ricker wavelet')
                call source_wavelet_ricker
            endif
        else !wavelet file exists
            call alloc(shot%src%wavelet,shot%src%nt)
            open(11,file=file_wavelet,access='direct',recl=4*shot%src%nt)
            read(11,rec=1) shot%src%wavelet
            close(11)
        endif

        scale_wavelet=setup_get_char('SCALE_WAVELET',default='by dtdx')
        
        if(scale_wavelet/='no') then
            if(scale_wavelet=='by dtdx' .or. scale_wavelet=='by dxdt') then
                !scale wavelet to be dt, dx independent
                shot%src%wavelet = shot%src%wavelet* shot%src%dt/m%cell_volume
            else
                !user defined scaler
                read(scale_wavelet,*) scaler
                shot%src%wavelet = shot%src%wavelet* scaler
            endif
        endif

        open(12,file='source_wavelet',access='direct',recl=4*shot%src%nt)
        write(12,rec=1) shot%src%wavelet
        close(12)


        !! space !!

        !shift position to be 0-based
        shot%src%x=shot%src%x - m%ox
        shot%src%y=shot%src%y - m%oy
        shot%src%z=shot%src%z - m%oz
        shot%rcv(:)%x=shot%rcv(:)%x - m%ox
        shot%rcv(:)%y=shot%rcv(:)%y - m%oy
        shot%rcv(:)%z=shot%rcv(:)%z - m%oz

        !Hicks interpolation
        call hicks_init([m%dx,m%dy,m%dz],m%is_cubic, m%if_freesurface)
        
        !hicks coeff for source point
        if(if_staggered_grid) then
            !for vx,vy,vz components, hicks%x,y,z should be added by m%dx/2,dy/2,dz/2, respectively,
            !because v(1) is actually v[0.5], v(2) is v[1.5] etc.
            if(shot%src%comp=='vx') call hicks_init_position([shot%src%x+m%dx/2.,shot%src%y,        shot%src%z        ])
            if(shot%src%comp=='vy') call hicks_init_position([shot%src%x        ,shot%src%y+m%dy/2.,shot%src%z        ])
            if(shot%src%comp=='vz') call hicks_init_position([shot%src%x        ,shot%src%y,        shot%src%z+m%dz/2.])
        else
            call hicks_init_position([shot%src%x,shot%src%y,shot%src%z])
        endif

        if(shot%src%comp(1)=='P') then !explosive source or non-vertical force
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

            if(if_staggered_grid) then
                !for vx,vy,vz components, hicks%x,y,z should be added by m%dx/2,dy/2,dz/2, respectively,
                !because v(1) is actually v[0.5], v(2) is v[1.5] etc.
                if(shot%rcv(i)%comp=='vx') call hicks_init_position([shot%rcv(i)%x+m%dx/2.,shot%rcv(i)%y,        shot%rcv(i)%        ])
                if(shot%rcv(i)%comp=='vy') call hicks_init_position([shot%rcv(i)%x        ,shot%rcv(i)%y+m%dy/2.,shot%rcv(i)%        ])
                if(shot%rcv(i)%comp=='vz') call hicks_init_position([shot%rcv(i)%x        ,shot%rcv(i)%y,        shot%rcv(i)%+m%dz/2.])
            else
                call hicks_init_position([shot%rcv(i)%,shot%rcv(i)%y,shot%rcv(i)%z])
            endif

            if(shot%rcv(i)%comp(1)=='P') then !explosive source or non-vertical force
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
    
    subroutine shot_shift_by_computebox(iox,ioy,ioz)
        !source side
        shot%src%ix=shot%src%ix-iox+1
        shot%src%iy=shot%src%iy-ioy+1
        shot%src%iz=shot%src%iz-ioz+1
        shot%src%ifx=shot%src%ifx-iox+1; shot%src%ilx=shot%src%ilx-iox+1
        shot%src%ify=shot%src%ify-ioy+1; shot%src%ily=shot%src%ily-ioy+1
        shot%src%ifz=shot%src%ifz-ioz+1; shot%src%ilz=shot%src%ilz-ioz+1
        shot%src%x=shot%src%x-(iox-1)*m%dx
        shot%src%y=shot%src%y-(ioy-1)*m%dy
        shot%src%z=shot%src%z-(ioz-1)*m%dz
        
        !receiver side
        do ir=1,shot%nrcv
            shot%rcv(ir)%ix=shot%rcv(ir)%ix-iox+1
            shot%rcv(ir)%iy=shot%rcv(ir)%iy-ioy+1
            shot%rcv(ir)%iz=shot%rcv(ir)%iz-ioz+1
            shot%rcv(ir)%ifx=shot%rcv(ir)%ifx-iox+1; shot%rcv(ir)%ilx=shot%rcv(ir)%ilx-iox+1
            shot%rcv(ir)%ify=shot%rcv(ir)%ify-ioy+1; shot%rcv(ir)%ily=shot%rcv(ir)%ily-ioy+1
            shot%rcv(ir)%ifz=shot%rcv(ir)%ifz-ioz+1; shot%rcv(ir)%ilz=shot%rcv(ir)%ilz-ioz+1
            shot%rcv(ir)%x=shot%rcv(ir)%x-(iox-1)*m%dx
            shot%rcv(ir)%y=shot%rcv(ir)%y-(ioy-1)*m%dy
            shot%rcv(ir)%z=shot%rcv(ir)%z-(ioz-1)*m%dz
        enddo
        
        call hud('If forces, shot''s source & receiver positions have been shifted according to computebox ox,oy,oz')
    end subroutine

    !======= private procedures =======

    subroutine from_gen_acqui
        character(4),dimension(gen_acqui_ntr)   :: rc !rcv comp
        real,dimension(gen_acqui_ntr)           :: rz,rx,ry !rcv pos


        call gen_acqui_shotgather(shot%index, &&
                                shot%src%comp,shot%src%z,shot%src%x,shot%src%y, &&
                                rc,rz,rx,ry, &&
                                shot%nt, shot%dt)
        
        shot%nrcv=gen_acqui_ntr
        if(allocated(shot%rcv))deallocate(shot%rcv)
        allocate(shot%rcv(gen_acqui_ntr))
        
        do i=1,shot%nrcv
            shot%rcv(i)%comp=rc(i)
            shot%rcv(i)%z=acqui%src(shot%index)%rcv(ir)%z
            shot%rcv(i)%x=acqui%src(shot%index)%rcv(ir)%x
            shot%rcv(i)%y=acqui%src(shot%index)%rcv(ir)%y
            shot%rcv(i)%aoffset=sqrt( (shot%src%x-shot%rcv(ir)%x)**2 &
                                      +(shot%src%y-shot%rcv(ir)%y)**2 )
        enddo

        shot%rcv(:)%nt=shot%nt
        shot%rcv(:)%dt=shot%dt
        
    end subroutine
    
    subroutine from_data
        type(t_suformat),dimension(:),allocatable :: sudata
        
        call read_sudata(shot%cindex,sudata)
        
        !source & receiver geometry
        shot%src%x=sudata(1)%hdr%sx
        shot%src%y=sudata(1)%hdr%sy
        shot%src%z=sudata(1)%hdr%sdepth
        
        shot%src%icomp=1 !I don't know which su header tells this info..
        
        shot%src%nt=sudata(1)%hdr%ns     !let's assume source's nt is same as receiver's nt
        shot%src%dt=sudata(1)%hdr%dt*1e-6 !let's assume source's dt is same as receiver's dt
        
        shot%nrcv=sudata(1)%hdr%ntr
        if(allocated(shot%rcv))deallocate(shot%rcv)
        allocate(shot%rcv(shot%nrcv))
        do ir=1,shot%nrcv
            shot%rcv(ir)%x= sudata(ir)%hdr%gx
            shot%rcv(ir)%y= sudata(ir)%hdr%gy
            shot%rcv(ir)%z=-sudata(ir)%hdr%gelev
            
            select case (sudata(ir)%hdr%trid)
                case (11)
                shot%rcv(ir)%comp='P'  !pressure
                case (12)
                shot%rcv(ir)%comp='vz' !vertical velocity
                case (13)
                shot%rcv(ir)%comp='vy' !horizontal velocity
                case (14)
                shot%rcv(ir)%comp='vx'
                
                case (1)
                shot%rcv(ir)%comp='P'
                case default
                shot%rcv(ir)%comp='P'
            end select
            
            shot%rcv(ir)%nt=sudata(ir)%hdr%ns
            shot%rcv(ir)%dt=sudata(ir)%hdr%dt*1e-6
        enddo
        
        
        !scaling elevation
        if(sudata(1)%hdr%scalel > 0) then
            shot%src%z    = shot%src%z    * sudata(1)%hdr%scalel
            shot%rcv(:)%z = shot%rcv(:)%z * sudata(1)%hdr%scalel
        elseif(sudata(1)%hdr%scalel < 0) then
            shot%src%z    = shot%src%z    / (-sudata(1)%hdr%scalel)
            shot%rcv(:)%z = shot%rcv(:)%z / (-sudata(1)%hdr%scalel)
        endif
        
        !scaling coordinates
        if(sudata(1)%hdr%scalco > 0) then
            shot%src%x    = shot%src%x    * sudata(1)%hdr%scalco
            shot%src%y    = shot%src%y    * sudata(1)%hdr%scalco
            shot%rcv(:)%x = shot%rcv(:)%x * sudata(1)%hdr%scalco
            shot%rcv(:)%y = shot%rcv(:)%y * sudata(1)%hdr%scalco
        elseif(sudata(1)%hdr%scalco < 0) then
            shot%src%x    = shot%src%x    / (-sudata(1)%hdr%scalco)
            shot%src%y    = shot%src%y    / (-sudata(1)%hdr%scalco)
            shot%rcv(:)%x = shot%rcv(:)%x / (-sudata(1)%hdr%scalco)
            shot%rcv(:)%y = shot%rcv(:)%y / (-sudata(1)%hdr%scalco)
        endif
        
        do ir=1,shot%nrcv
!            shot%rcv(ir)%aoffset=sqrt( (shot%src%x-shot%rcv(ir)%x)**2 &
!                                      +(shot%src%y-shot%rcv(ir)%y)**2 &
!                                      +(shot%src%z-shot%rcv(ir)%z)**2 )
            shot%rcv(ir)%aoffset=sqrt( (shot%src%x-shot%rcv(ir)%x)**2 &
                                      +(shot%src%y-shot%rcv(ir)%y)**2 )
        enddo
        
        !load obs traces
        call alloc(dobs,shot%rcv(1)%nt,shot%nrcv)
        do ir=1,shot%nrcv
            dobs(:,ir)=sudata(ir)%trace
        enddo
        
        !clean su data
        do ir=1,shot%nrcv
            deallocate(sudata(ir)%trace)
        enddo
        deallocate(sudata)
        
    end subroutine
    
    subroutine source_wavelet_sinexp

        a=-3.3333333*shot%src%fpeak
        
        call alloc(shot%src%wavelet,shot%src%nt)
        
        do it=1,shot%src%nt
            t=(it-1)*shot%src%dt
            
            shot%src%wavelet(it)=sin(2.*r_pi*shot%src%fpeak*t)*exp(a*t)
        enddo
        
        !butterworth filtering to mitigate spectrum high-end tail
        call butterworth(1,shot%src%nt, shot%src%dt, shot%src_wavelet,&
        o_zerophase=.false.,o_locut=.false.,&
        o_fpasshi=fpeak,o_fstophi=2.*fpeak,&
                        o_astophi=0.1)
        
    end subroutine
    
    subroutine source_wavelet_ricker

        t0=setup_get_real('RICKER_DELAYTIME','T0',default=1./shot%src%fpeak)
        
        if (shot%src%fpeak*2.5 > 1./shot%src%dt) then
            call hud('Ricker wavelet peak frequency too high (fpeak*2.5 > 1/dt). Reduce it.')
        endif

        call alloc(shot%src%wavelet,shot%src%nt)

        do it=1,shot%src%nt

            t=(it-1)*shot%src%dt-t0

            x=r_pi*shot%src%fpeak*t
            x=-x*x
            shot%src%wl(it)=(1.+2.*x)*exp(x)

        enddo

    end subroutine


    subroutine update_wavelet
            
        character(:),allocatable :: update_wavelet
        
        !if nshot_per_processor is not same for each processor,
        !update_wavelet='stack' mode will be stuck due to collective communication in m_matchfilter.f90
        if(nshot_per_processor * mpiworld%nproc /= nshots) then
            fatal('Unequal shot numbers on processors. If you are using UPDATE_WAVELET=''stack'', the code will be stuck due to collective communication in m_matchfilter')
        endif

        if(update_wavelet=='stack') then
            !average wavelet across all processors
            !note: if more shots than processors, non-assigned shots will not contribute to this averaging
            call matchfilter_estimate(shot%src%nt,shot%nrcv,dsyn,dobs,shot%index,if_stack=.true.)
        else
            call matchfilter_estimate(shot%src%nt,shot%nrcv,dsyn,dobs,shot%index,if_stack=.false.)
        endif
        
        call matchfilter_apply_to_wavelet(shot%src%nt,shot%src%wavelet)
        
        call matchfilter_apply_to_data(shot%src%nt,shot%nrcv,dsyn)
        
        call suformat_write(iproc=0,file='wavelet_update',append=.true.)
    end subroutine
    
    subroutine update_adjsrc
        call matchfilter_correlate_filter_residual(shot%src%nt,shot%nrcv,shot%dres)
    end subroutine

!     subroutine shot_check_rcv_ranges(mx,my,mz)
!         integer :: problem
!         problem=0
!         
!         do ir=1,shot%nrcv
!             if(shot%rcv(ir)%ix>mx .or. shot%rcv(ir)%iy>my .or. shot%rcv(ir)%iz>mz) then
!                 problem=problem+1
!             endif
!         enddo
!         
!         if(problem>0) then
!             call hud('WARNING: some receivers are outside the computebox! They will be truncated.')
!         endif
!     end subroutine
    
!     subroutine read_shot_wavelet
!         logical alive
!         inquire(file=setup%file_wavelet, exist=alive)
!         if(alive) then
!             open(12,file=setup%file_wavelet,access='direct',recl=4*setup%nt,action='read',status='old')
!             call alloc(shot%src%wavelet(setup%nt))
!                 read(12,rec=1) shot%src%wavelet
!             enddo
!             close(12)
!         else
!             stop 'FILE_WAVELET doesn''t exist.'
!         endif
!        
!    end subroutine

    
end
