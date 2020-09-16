module m_shot
use m_sysio
use m_shotlist
use m_gen_acquisition, only: acqui
! use m_gen_wavelet
use m_suformat
use m_model, only: m
use m_hicks

    private t_receiver, t_source, t_shot!, file_wavelet
    
    type t_receiver
        real    :: x,y,z, aoffset
        integer :: ix,iy,iz
        integer :: ifz,ilz,ifx,ilx,ify,ily
        integer :: icomp        !component: 1=P, 2=vx, 3=vy, 4=vz
        real,dimension(:,:,:),allocatable :: interp_coeff
        integer :: nt
        real :: dt
    end type
    
    type t_source
        real    :: x,y,z
        integer :: ix,iy,iz
        integer :: ifz,ilz,ifx,ilx,ify,ily
        integer :: icomp        !component: 1=P, 2=vx, 3=vy, 4=vz
        real,dimension(:,:,:),allocatable :: interp_coeff
        real,dimension(:),allocatable :: wavelet
        integer :: nt
        real :: dt, fpeak
        
    end type
    
    type t_shot
        integer :: index
        character(4) :: cindex
        type(t_source) :: src
        integer :: nrcv
        type(t_receiver),dimension(:),allocatable :: rcv  !should be in size of nrcv
        logical :: if_hicks
    end type
    
    type(t_shot) :: shot
    
    real,dimension(:,:),allocatable :: dobs !observed seismogram
    real,dimension(:,:),allocatable :: dsyn !synthetic seismogram
    real,dimension(:,:),allocatable :: dres !residual as well as adjoint source seismogram
    
    ! character(:),allocatable :: file_wavelet
    
    contains
    
    subroutine init_shot(k,from)
        integer k
        character(*) :: from

        logical :: if_staggered_grid
        
        character(:),allocatable :: scale_wavelet

        shot%index=shotlist(k)
        write(shot%cindex,'(i0.4)') shot%index
        
        !read geometry, nt and dt
        if(from=='setup') then
            call init_shot_from_setup
        elseif(from=='data') then
            call init_shot_from_data
        endif
        
        !shift position to be 0-based
        shot%src%x=shot%src%x - m%ox
        shot%src%y=shot%src%y - m%oy
        shot%src%z=shot%src%z - m%oz
        shot%rcv(:)%x=shot%rcv(:)%x - m%ox
        shot%rcv(:)%y=shot%rcv(:)%y - m%oy
        shot%rcv(:)%z=shot%rcv(:)%z - m%oz
        
        ! !read wavelet
        ! shot%src%fpeak=get_setup_real('PEAK_FREQUENCY')
        ! file_wavelet=get_setup_file('FILE_WAVELET')
        ! if(file_wavelet=='') then
        !     if(get_setup_char('WAVELET_TYPE',default='sinexp')=='sinexp') then
        !         call hud('Use filtered sinexp wavelet')
        !         shot%src%wavelet=gen_wavelet_sinexp(shot%src%nt,shot%src%dt,shot%src%fpeak)
        !     else
        !         call hud('Use Ricker wavelet')
        !         shot%src%wavelet=gen_wavelet_ricker(shot%src%nt,shot%src%dt,shot%src%fpeak)
        !     endif
        ! else !wavelet file exists
        !     call alloc(shot%src%wavelet,shot%src%nt)
        !     open(11,file=file_wavelet,access='direct',recl=4*shot%src%nt)
        !     read(11,rec=1) shot%src%wavelet
        !     close(11)
        ! endif

        ! scale_wavelet=get_setup_char('SCALE_WAVELET',default='no')

        ! if(trim(adjustl(scale_wavelet))/='no') then
        !     if(trim(adjustl(scale_wavelet))=='by dtdx') then
        !         !scale wavelet to be dt, dx independent
        !         shot%src%wavelet = shot%src%wavelet* shot%src%dt/m%cell_volume
        !     else
        !         !user defined scaler
        !         read(scale_wavelet,*) scaler
        !         shot%src%wavelet = shot%src%wavelet* scaler
        !     endif
        ! endif

        
        !hicks coeff for source point
        hicks%x=shot%src%x; hicks%dx=m%dx
        hicks%y=shot%src%y; hicks%dy=m%dy
        hicks%z=shot%src%z; hicks%dz=m%dz
        hicks%is_cubic=m%is_cubic
        hicks%if_freesurface=m%if_freesurface
        
        !for vx,vy,vz, hicks%x,y,z should be added by m%dx/2,dy/2,dz/2, respectively, as v(1) is actually v[0.5], v(2) is v[1.5] etc.
        if_staggered_grid=get_setup_logical('IF_STAGGERED_GRID',default=.true.)
        if(if_staggered_grid) then
            if(shot%src%icomp==2) hicks%x=shot%src%x+m%dx/2.
            if(shot%src%icomp==3) hicks%y=shot%src%y+m%dy/2.
            if(shot%src%icomp==4) hicks%z=shot%src%z+m%dz/2.
        endif
        
        if(shot%src%icomp==1) then !explosive source or non-vertical force
            call build_hicks(hicks,'antisym',  shot%src%interp_coeff)
        elseif(shot%src%icomp==4) then !vertical force
            call build_hicks(hicks,'symmetric',shot%src%interp_coeff)
        else
            call build_hicks(hicks,'truncate', shot%src%interp_coeff)
        endif
        
        shot%src%ifz=hicks%ifz; shot%src%iz=hicks%iz; shot%src%ilz=hicks%ilz
        shot%src%ifx=hicks%ifx; shot%src%ix=hicks%ix; shot%src%ilx=hicks%ilx
        shot%src%ify=hicks%ify; shot%src%iy=hicks%iy; shot%src%ily=hicks%ily
        
        !hicks coeff for receivers
        do i=1,shot%nrcv
            hicks%x=shot%rcv(i)%x
            hicks%y=shot%rcv(i)%y
            hicks%z=shot%rcv(i)%z
            
            !for vx,vy,vz, hicks%x,y,z should be added by m%dx/2,dy/2,dz/2, respectively, as v(1) is actually v[0.5], v(2) is v[1.5] etc.
            if(if_staggered_grid) then
                if(shot%rcv(i)%icomp==2) hicks%x=shot%rcv(i)%x+m%dx/2.
                if(shot%rcv(i)%icomp==3) hicks%y=shot%rcv(i)%y+m%dy/2.
                if(shot%rcv(i)%icomp==4) hicks%z=shot%rcv(i)%z+m%dz/2.
            endif
            
            if(shot%rcv(i)%icomp==1) then !pressure component
                call build_hicks(hicks,'antisym',  shot%rcv(i)%interp_coeff)
            elseif(shot%rcv(i)%icomp==4) then !vz component
                call build_hicks(hicks,'symmetric',shot%rcv(i)%interp_coeff)
            else
                call build_hicks(hicks,'truncate', shot%rcv(i)%interp_coeff)
            endif
            
            shot%rcv(i)%ifz=hicks%ifz; shot%rcv(i)%iz=hicks%iz; shot%rcv(i)%ilz=hicks%ilz
            shot%rcv(i)%ifx=hicks%ifx; shot%rcv(i)%ix=hicks%ix; shot%rcv(i)%ilx=hicks%ilx
            shot%rcv(i)%ify=hicks%ify; shot%rcv(i)%iy=hicks%iy; shot%rcv(i)%ily=hicks%ily
        
        enddo
                
        if(mpiworld%is_master) then
            write(*,*)'================================='
            write(*,*)'Shot# '//shot%cindex//' info:'
            write(*,*)'================================='
            ! write(*,*)'  nt,dt:',shot%src%nt,shot%src%dt
            ! write(*,*)'---------------------------------'
            ! write(*,*)'  sz,isz:',shot%src%z,shot%src%iz
            ! write(*,*)'  sx,isx:',shot%src%x,shot%src%ix
            ! write(*,*)'  sy,isy:',shot%src%y,shot%src%iy
            ! write(*,*)'  ifz,ilz:',shot%src%ifz,shot%src%ilz
            ! write(*,*)'  ifx,ilx:',shot%src%ifx,shot%src%ilx
            ! write(*,*)'  ify,ily:',shot%src%ify,shot%src%ily
            ! write(*,*)'---------------------------------'
            ! write(*,*)'  minmax rz,irz:',minval(shot%rcv(:)%z),maxval(shot%rcv(:)%z),minval(shot%rcv(:)%iz),maxval(shot%rcv(:)%iz)
            ! write(*,*)'  minmax rx,irx:',minval(shot%rcv(:)%x),maxval(shot%rcv(:)%x),minval(shot%rcv(:)%ix),maxval(shot%rcv(:)%ix)
            ! write(*,*)'  minmax ry,iry:',minval(shot%rcv(:)%y),maxval(shot%rcv(:)%y),minval(shot%rcv(:)%iy),maxval(shot%rcv(:)%iy)
            ! write(*,*)'  minmax ifz,ilz:',minval(shot%rcv(:)%ifz),maxval(shot%rcv(:)%ifz),minval(shot%rcv(:)%ilz),maxval(shot%rcv(:)%ilz)
            ! write(*,*)'  minmax ifx,ilx:',minval(shot%rcv(:)%ifx),maxval(shot%rcv(:)%ifx),minval(shot%rcv(:)%ilx),maxval(shot%rcv(:)%ilx)
            ! write(*,*)'  minmax ify,ily:',minval(shot%rcv(:)%ify),maxval(shot%rcv(:)%ify),minval(shot%rcv(:)%ily),maxval(shot%rcv(:)%ily)
            ! write(*,*)'  nrcv:',shot%nrcv
            ! write(*,*)'---------------------------------'
            write(*,*)'  sz,sx:',shot%src%z,shot%src%x
            write(*,*)'  ifz:ilz:',shot%src%ifz,shot%src%ilz
            write(*,*)'  ifx:ilx:',shot%src%ifx,shot%src%ilx
            write(*,*)'---------------------------------'
            write(*,*)'  minmax rz:',minval(shot%rcv(:)%z),maxval(shot%rcv(:)%z)
            write(*,*)'  minmax rx:',minval(shot%rcv(:)%x),maxval(shot%rcv(:)%x)
            write(*,*)'  nrcv:',shot%nrcv
            write(*,*)'---------------------------------'
        endif
        
    end subroutine
    
    subroutine init_shot_from_setup
        
        !source side
        shot%src%x=acqui%src(shot%index)%x
        shot%src%y=acqui%src(shot%index)%y
        shot%src%z=acqui%src(shot%index)%z
        shot%src%icomp=acqui%iscomp
        ! shot%src%nt=acqui%nt
        ! shot%src%dt=acqui%dt
        
        !receiver side
        shot%nrcv=acqui%src(shot%index)%nrcv
        if(allocated(shot%rcv))deallocate(shot%rcv)
        allocate(shot%rcv(shot%nrcv))
        do ir=1,shot%nrcv
            shot%rcv(ir)%x=acqui%src(shot%index)%rcv(ir)%x
            shot%rcv(ir)%y=acqui%src(shot%index)%rcv(ir)%y
            shot%rcv(ir)%z=acqui%src(shot%index)%rcv(ir)%z
!            shot%rcv(ir)%aoffset=sqrt( (shot%src%x-shot%rcv(ir)%x)**2 &
!                                      +(shot%src%y-shot%rcv(ir)%y)**2 &
!                                      +(shot%src%z-shot%rcv(ir)%z)**2 )
            shot%rcv(ir)%aoffset=sqrt( (shot%src%x-shot%rcv(ir)%x)**2 &
                                      +(shot%src%y-shot%rcv(ir)%y)**2 )
            shot%rcv(ir)%icomp=acqui%ircomp
            ! shot%rcv(ir)%nt=acqui%nt
            ! shot%rcv(ir)%dt=acqui%dt
        enddo
        
    end subroutine
    
    
    subroutine init_shot_from_data
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
                case (11) !pressure
                shot%rcv(ir)%icomp=1
                case (12) !vertical component
                shot%rcv(ir)%icomp=4  !vz
                case (13) !cross-line component
                shot%rcv(ir)%icomp=3  !vy
                case (14) !in-line component
                shot%rcv(ir)%icomp=2  !vx
                
                case (1)  !pressure
                shot%rcv(ir)%icomp=1
                case default !pressure
                shot%rcv(ir)%icomp=1
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
