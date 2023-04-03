module m_field
use m_System
use m_resampler
use m_model
use m_shot
use m_computebox

use, intrinsic :: ieee_arithmetic

    private
    
    public :: field_init

    !boundary components for wavefield recontruction
    type t_boundary

        real,dimension(:,:),allocatable :: vz_top,  vz_bot
        real,dimension(:,:),allocatable :: vx_left, vx_right
        real,dimension(:,:),allocatable :: vy_front,vy_rear
        real,dimension(:,:),allocatable :: vx_top,vx_bot,vz_left,vz_right

        real,dimension(:,:),allocatable :: p_top,  p_bot
        real,dimension(:,:),allocatable :: p_left, p_right
        real,dimension(:,:),allocatable :: p_front,p_rear

    end type

    !fields
    type,public :: t_field

        character(:),allocatable :: name

        !adjointness
        logical :: is_adjoint
        
        !wavefield components in computation domain
        real,dimension(:,:,:),allocatable :: vz,vx,vy !velocities
        real,dimension(:,:,:),allocatable :: szz,szx,szy !stress tensor
        real,dimension(:,:,:),allocatable ::     sxx,sxy
        real,dimension(:,:,:),allocatable ::         syy
        real,dimension(:,:,:),allocatable :: shh !szz, sxx or syy
        real,dimension(:,:,:),allocatable :: p, p_prev, p_next !negated pressure
        !N.B. pressure is defined >0 for inward stress, but here tobe compatible with szz etc, p is defined >0 for outward stress

        !boundary components for wavefield recontruction
        logical :: if_will_reconstruct=.false.
        type(t_boundary) :: bnd

        !cpml components for absorbing boundary wavefield
        real,dimension(:,:,:),allocatable :: dvz_dz,dvz_dx,dvz_dy
        real,dimension(:,:,:),allocatable :: dvx_dz,dvx_dx,dvx_dy
        real,dimension(:,:,:),allocatable :: dvy_dz,dvy_dx,dvy_dy
        
        real,dimension(:,:,:),allocatable :: dszz_dz,dszx_dz
        real,dimension(:,:,:),allocatable :: dszx_dx,dsxx_dx        
        real,dimension(:,:,:),allocatable :: dshh_dz,dshh_dx,dshh_dy
        real,dimension(:,:,:),allocatable :: dp_dz,dp_dx,dp_dy
        real,dimension(:,:,:),allocatable :: dpzz_dz,dpxx_dx,dpyy_dy
        real,dimension(:,:,:),allocatable :: lap
        real,dimension(:,:,:),allocatable :: lapz,lapx,lapy

        !source time function
        ! real,dimension(:,:),allocatable :: fz,fx,fy !forces
        ! real,dimension(:,:),allocatable :: mzz,mzx,mzy !moments
        ! real,dimension(:,:),allocatable ::     mxx,mxy
        ! real,dimension(:,:),allocatable ::         myy
        ! real,dimension(:,:),allocatable :: mhh,p
        real,dimension(:,:),allocatable :: wavelet

        !synthetic seismograms (receiver time function)
        real,dimension(:,:),allocatable :: seismo
        
        !blooming
        integer,dimension(:,:),allocatable :: bloom

        !snapshot

        contains
        ! procedure :: init
        procedure :: init_bloom
        procedure :: init_boundary
        procedure :: init_boundary_pressure
        procedure :: reinit
        procedure :: check_value
        procedure :: ignite
        procedure :: acquire
        procedure :: write
        procedure :: write_ext
        procedure :: boundary_transport
        procedure :: boundary_transport_pressure
        final :: final
        
        procedure :: is_registered
        procedure :: register

    end type

    !propagator's nt, dt
    integer :: nt
    real :: dt

    !shear components
    logical :: if_shear

    !bloombox
    logical :: if_bloom
    integer,parameter :: initial_half_bloomwidth=5 !half bloom width at ift, should involve hicks points

    !snapshot
    logical :: if_snapshot
    type(t_string),dimension(:),allocatable :: snapshot
    integer :: i_snapshot, n_snapshot
    
    contains
    
    subroutine field_init(if_shear_in,nt_in,dt_in)
        logical if_shear_in

        logical,save :: is_first_in=.true.

        if_shear=if_shear_in
        nt=nt_in
        dt=dt_in

        !snapshot
        snapshot=setup%get_strs('SNAPSHOT')
        if_snapshot=size(snapshot)>0 .and. mpiworld%is_master
        if(if_snapshot) then
            n_snapshot=setup%get_int('REF_NUMBER_SNAPSHOT','NSNAPSHOT',o_default='50')
            if(n_snapshot==0) n_snapshot=50
            i_snapshot=ceiling(nt*1./n_snapshot)

            !rm existing snap files
            if(is_first_in) then
                call sysio_rm('snap*')
                is_first_in=.false.
            endif

        endif

    end subroutine

    ! subroutine init(self,name)
    !     class(t_field) :: self
    !     character(*) :: name

    !     self%name=name

    ! end subroutine

    subroutine init_bloom(self)
        class(t_field) :: self

        !blooming
        call alloc(self%bloom,6,nt)

        if_bloom=setup%get_bool('IF_BLOOM',o_default='T')      

        if(if_bloom) then
            !directional maximum propagation distance per time step
            distz = cb%velmax * dt / m%dz
            distx = cb%velmax * dt / m%dx
            disty = cb%velmax * dt / m%dy

            if(.not.self%is_adjoint) then
                self%bloom(1,1)=max(shot%src%iz -initial_half_bloomwidth, cb%ifz)
                self%bloom(2,1)=min(shot%src%iz +initial_half_bloomwidth, cb%ilz)
                self%bloom(3,1)=max(shot%src%ix -initial_half_bloomwidth, cb%ifx)
                self%bloom(4,1)=min(shot%src%ix +initial_half_bloomwidth, cb%ilx)
                self%bloom(5,1)=max(shot%src%iy -initial_half_bloomwidth, cb%ify)
                self%bloom(6,1)=min(shot%src%iy +initial_half_bloomwidth, cb%ily)
                do it=2,nt
                    self%bloom(1,it)=max(nint(self%bloom(1,1)-it*distz),cb%ifz) !bloombox ifz
                    self%bloom(2,it)=min(nint(self%bloom(2,1)+it*distz),cb%ilz) !bloombox ilz
                    self%bloom(3,it)=max(nint(self%bloom(3,1)-it*distx),cb%ifx) !bloombox ifx
                    self%bloom(4,it)=min(nint(self%bloom(4,1)+it*distx),cb%ilx) !bloombox ilx
                    self%bloom(5,it)=max(nint(self%bloom(5,1)-it*disty),cb%ify) !bloombox ify
                    self%bloom(6,it)=min(nint(self%bloom(6,1)+it*disty),cb%ily) !bloombox ily
                enddo
            
            else
                self%bloom(1,nt)=max(minval(shot%rcv(:)%iz) -initial_half_bloomwidth, cb%ifz)
                self%bloom(2,nt)=min(maxval(shot%rcv(:)%iz) +initial_half_bloomwidth, cb%ilz)
                self%bloom(3,nt)=max(minval(shot%rcv(:)%ix) -initial_half_bloomwidth, cb%ifx)
                self%bloom(4,nt)=min(maxval(shot%rcv(:)%ix) +initial_half_bloomwidth, cb%ilx)
                self%bloom(5,nt)=max(minval(shot%rcv(:)%iy) -initial_half_bloomwidth, cb%ify)
                self%bloom(6,nt)=min(maxval(shot%rcv(:)%iy) +initial_half_bloomwidth, cb%ily)
                do it=nt-1,1,-1
                    it_fwd=nt-it+1
                    self%bloom(1,it)=max(nint(self%bloom(1,nt)-it_fwd*distz),cb%ifz) !bloombox ifz
                    self%bloom(2,it)=min(nint(self%bloom(2,nt)+it_fwd*distz),cb%ilz) !bloombox ilz
                    self%bloom(3,it)=max(nint(self%bloom(3,nt)-it_fwd*distx),cb%ifx) !bloombox ifx
                    self%bloom(4,it)=min(nint(self%bloom(4,nt)+it_fwd*distx),cb%ilx) !bloombox ilx
                    self%bloom(5,it)=max(nint(self%bloom(5,nt)-it_fwd*disty),cb%ify) !bloombox ify
                    self%bloom(6,it)=min(nint(self%bloom(6,nt)+it_fwd*disty),cb%ily) !bloombox ily
                enddo

            endif

        else
            self%bloom(1,:)=cb%ifz
            self%bloom(2,:)=cb%ilz
            self%bloom(3,:)=cb%ifx
            self%bloom(4,:)=cb%ilx
            self%bloom(5,:)=cb%ify
            self%bloom(6,:)=cb%ily

        endif

        if(.not.m%is_cubic) self%bloom(5:6,:)=1

    end subroutine

    subroutine init_boundary(self)
        class(t_field) :: self
        !save 3 grid points, for 4th order FD only
        !different indexing
        n=3*cb%mx*cb%my
        call alloc(self%bnd%vz_top,n,nt)
        call alloc(self%bnd%vz_bot,n,nt)
        if(if_shear) then
            call alloc(self%bnd%vx_top,n,nt)
            call alloc(self%bnd%vx_bot,n,nt)
        endif
        
        n=cb%mz*3*cb%my
        call alloc(self%bnd%vx_left, n,nt)
        call alloc(self%bnd%vx_right,n,nt)
        if(if_shear) then
            call alloc(self%bnd%vz_left, n,nt)
            call alloc(self%bnd%vz_right,n,nt)
        endif

        if(m%is_cubic) then
            n=cb%mz*cb%mx*3
            call alloc(self%bnd%vy_front,n,nt)
            call alloc(self%bnd%vy_rear, n,nt)
        endif

    end subroutine

    subroutine init_boundary_pressure(self)
        class(t_field) :: self
        !save 3 grid points, for 4th order FD only
        !different indexing
        n=3*cb%mx*cb%my
        call alloc(self%bnd%p_top,n,nt)
        call alloc(self%bnd%p_bot,n,nt)
        ! if(if_shear) then
        !     call alloc(self%bnd%p_top,n,nt)
        !     call alloc(self%bnd%p_bot,n,nt)
        ! endif
        
        n=cb%mz*3*cb%my
        call alloc(self%bnd%p_left, n,nt)
        call alloc(self%bnd%p_right,n,nt)
        ! if(if_shear) then
        !     call alloc(self%bnd%p_left, n,nt)
        !     call alloc(self%bnd%p_right,n,nt)
        ! endif

        ! if(m%is_cubic) then
        !     n=cb%mz*cb%mx*3
        !     call alloc(self%bnd%vy_front,n,nt)
        !     call alloc(self%bnd%vy_rear, n,nt)
        ! endif

    end subroutine

    subroutine reinit(self)
        class(t_field) :: self

        if(allocated(self%dp_dz))   self%dp_dz=0.
        if(allocated(self%dp_dx))   self%dp_dx=0.
        if(allocated(self%dp_dy))   self%dp_dy=0.

        if(allocated(self%dpzz_dz))   self%dpzz_dz=0.
        if(allocated(self%dpxx_dx))   self%dpxx_dx=0.
        if(allocated(self%dpyy_dy))   self%dpyy_dy=0.

        if(allocated(self%lap)) self%lap=0.

        if(allocated(self%lapz)) self%lapz=0.
        if(allocated(self%lapx)) self%lapx=0.
        if(allocated(self%lapy)) self%lapy=0.

    end subroutine
    
    subroutine check_value(self)
        class(t_field) :: self
        
        if(mpiworld%is_master) write(*,*) self%name//' minmax values:',minval(self%p),maxval(self%p)
        
        if(any(.not. ieee_is_finite(self%p))) then
            call error(self%name//' values become Infinity on '//shot%sindex//' !!')
        endif
        if(any(ieee_is_nan(self%p))) then
            call error(self%name//' values become NaN on '//shot%sindex//' !!')
        endif
        
    end subroutine
        
    subroutine write(self,it,o_suffix)
        class(t_field) :: self
        character(*),optional :: o_suffix

        character(:),allocatable :: suf

        suf=either(o_suffix,'',present(o_suffix))

        if(if_snapshot) then

            if(it==1 .or. mod(it,i_snapshot)==0 .or. it==nt) then
                do i=1,size(snapshot)
                    select case (snapshot(i)%s)
                    ! case ('vz')
                    !     call sysio_write('snap_'//self%name//'%vz'//suf,self%vz,size(self%vz),o_mode='append')
                    ! case ('vx')
                    !     call sysio_write('snap_'//self%name//'%vx'//suf,self%vx,size(self%vx),o_mode='append')
                    ! case ('vy')
                    !     call sysio_write('snap_'//self%name//'%vx'//suf,self%vy,size(self%vy),o_mode='append')
                        
                    case ('p')
                        call sysio_write('snap_'//self%name//'%p'//suf,self%p,size(self%p),o_mode='append')
                    ! case ('szz')
                    !     call sysio_write('snap_'//self%name//'%szz'//suf,self%szz,size(self%szz),o_mode='append')
                    ! case ('sxx')
                    !     call sysio_write('snap_'//self%name//'%sxx'//suf,self%sxx,size(self%sxx),o_mode='append')
                    ! case ('szx')
                    !     call sysio_write('snap_'//self%name//'%szx'//suf,self%szx,size(self%szx),o_mode='append')
                    case ('dp_dz')
                        call sysio_write('snap_'//self%name//'%dp_dz'//suf,self%dp_dz,size(self%dp_dz),o_mode='append')
                    case ('lap')
                        call sysio_write('snap_'//self%name//'%lap'//suf,self%lap,size(self%lapz),o_mode='append')
                    end select
                enddo
            endif

        endif

    end subroutine

    subroutine write_ext(self,it,name,data,n)
        class(t_field) :: self
        character(*),optional :: name
        real,dimension(n) :: data

        if(if_snapshot) then

            if(it==1 .or. mod(it,i_snapshot)==0 .or. it==nt) then
                call sysio_write('snap_'//name,data,n,o_mode='append')
            endif

        endif

    end subroutine
 
    subroutine ignite(self,o_wavelet)
        class(t_field) :: self
        real,dimension(:,:),optional :: o_wavelet !use external wavelet instead of shot%wavelet or %dadj

        !add adjoint source
        if(self%is_adjoint) then

            if(present(o_wavelet)) then
                if(all(shape(o_wavelet)==[nt,shot%nrcv])) then
                    self%wavelet=transpose(o_wavelet)
                else
                    call hud('shape(o_wavelet) = '    //strcat(nums2strs(shape(o_wavelet))))
                    call hud('required shape = '//num2str(nt)//' , '//num2str(shot%nrcv))
                    call error('External o_wavelet do NOT have the required shape!')
                endif

            else
                call alloc(self%wavelet,shot%nrcv,nt)
                do i=1,shot%nrcv !implicit transpose
                    call resampler(shot%dadj(:,i),self%wavelet(i,:),1,&
                                    din=shot%dt,nin=shot%nt,&
                                    dout=dt,nout=nt)
                enddo
                
            endif

            return

        endif

        !add source
            if(present(o_wavelet)) then
                if(all(shape(o_wavelet)==[nt,1])) then
                    self%wavelet=transpose(o_wavelet)
                else
                    call hud('shape(o_wavelet) = '//strcat(nums2strs(shape(o_wavelet))))
                    call hud('required shape = '//num2str(nt)//', 1')
                    call error('External o_wavelet do NOT have the required shape!')
                endif

            else
                call alloc(self%wavelet,1,nt)
                call resampler(shot%wavelet,self%wavelet(1,:),1,&
                                din=shot%dt,nin=shot%nt,&
                                dout=dt,nout=nt)
            endif

    end subroutine

    subroutine acquire(self,o_seismo)
        class(t_field) :: self
        real,dimension(:,:),allocatable,optional :: o_seismo

        if(present(o_seismo)) then
            o_seismo=transpose(self%seismo)

        else

            call alloc(shot%dsyn,shot%nt,shot%nrcv)
            do i=1,shot%nrcv
                call resampler(self%seismo(i,:),shot%dsyn(:,i),1,&
                                din=dt,nin=nt,&
                                dout=shot%dt,nout=shot%nt)
            enddo

        endif

    end subroutine
        
    subroutine boundary_transport(self,action,it)
        class(t_field) :: self
        character(4) :: action
        integer :: it
        
        nz=cb%mz
        nx=cb%mx
        ny=cb%my
        if(m%is_cubic) then
            !top
            call copy(action,self%vz,self%bnd%vz_top(:,it),  [1,3],    [1,nx],[1,ny])  !old version: [0,2],[1,nx],[1,nx]
            !bottom
            call copy(action,self%vz,self%bnd%vz_bot(:,it),  [nz-1,nz+1],[1,nx],[1,ny])  !old version: [nz,nz+2],[1,nx],[1,nx]
            !left
            call copy(action,self%vx,self%bnd%vx_left(:,it), [1,nz],[1,3],    [1,ny])
            !right
            call copy(action,self%vx,self%bnd%vx_right(:,it),[1,nz],[nx-1,nx+1],[1,ny])
            !front
            call copy(action,self%vy,self%bnd%vy_front(:,it),[1,nz],[1,nx],[1,3])
            !rear
            call copy(action,self%vy,self%bnd%vy_rear(:,it), [1,nz],[1,nx],[ny-1,ny+1])
        else
            !top
            call copy(action,self%vz,self%bnd%vz_top(:,it),  [1,3],    [1,nx],[1,1])
            !bottom
            call copy(action,self%vz,self%bnd%vz_bot(:,it),  [nz-1,nz+1],[1,nx],[1,1])
            !left
            call copy(action,self%vx,self%bnd%vx_left(:,it), [1,nz],[1,3],    [1,1])
            !right
            call copy(action,self%vx,self%bnd%vx_right(:,it),[1,nz],[nx-1,nx+1],[1,1])
        endif

        !shear part
        if(if_shear) then
            if(m%is_cubic) then
            else
                !top
                call copy(action,self%vx,self%bnd%vx_top(:,it),  [1,3],    [1,nx],[1,1])
                !bottom
                call copy(action,self%vx,self%bnd%vx_bot(:,it),  [nz-2,nz  ],[1,nx],[1,1])
                !left
                call copy(action,self%vz,self%bnd%vz_left(:,it), [1,nz],[1,3],    [1,1])
                !right
                call copy(action,self%vz,self%bnd%vz_right(:,it),[1,nz],[nx-2,nx  ],[1,1])
            endif
        endif
        
    end subroutine
    
    subroutine boundary_transport_pressure(self,action,it)
        class(t_field) :: self
        character(4) :: action
        integer :: it
        
        nz=cb%mz
        nx=cb%mx
        ny=cb%my
        if(m%is_cubic) then
            !top
            ! call copy(action,self%p,self%bnd%vz_top(:,it),  [1,3],    [1,nx],[1,ny])  !old version: [0,2],[1,nx],[1,nx]
            ! !bottom
            ! call copy(action,self%p,self%bnd%vz_bot(:,it),  [nz-1,nz+1],[1,nx],[1,ny])  !old version: [nz,nz+2],[1,nx],[1,nx]
            ! !left
            ! call copy(action,self%p,self%bnd%vx_left(:,it), [1,nz],[1,3],    [1,ny])
            ! !right
            ! call copy(action,self%p,self%bnd%vx_right(:,it),[1,nz],[nx-1,nx+1],[1,ny])
            ! !front
            ! call copy(action,self%p,self%bnd%vy_front(:,it),[1,nz],[1,nx],[1,3])
            ! !rear
            ! call copy(action,self%p,self%bnd%vy_rear(:,it), [1,nz],[1,nx],[ny-1,ny+1])
        else
            !top
            call copy(action,self%p,self%bnd%p_top(:,it),  [1,3],    [1,nx],[1,1])
            !bottom
            call copy(action,self%p,self%bnd%p_bot(:,it),  [nz-1,nz+1],[1,nx],[1,1])
            !left
            call copy(action,self%p,self%bnd%p_left(:,it), [1,nz],[1,3],    [1,1])
            !right
            call copy(action,self%p,self%bnd%p_right(:,it),[1,nz],[nx-1,nx+1],[1,1])
        endif

        ! !shear part
        ! if(if_shear) then
        !     if(m%is_cubic) then
        !     else
        !         !top
        !         call copy(action,self%vx,self%bnd%vx_top(:,it),  [1,3],    [1,nx],[1,1])
        !         !bottom
        !         call copy(action,self%vx,self%bnd%vx_bot(:,it),  [nz-2,nz  ],[1,nx],[1,1])
        !         !left
        !         call copy(action,self%vz,self%bnd%vz_left(:,it), [1,nz],[1,3],    [1,1])
        !         !right
        !         call copy(action,self%vz,self%bnd%vz_right(:,it),[1,nz],[nx-2,nx  ],[1,1])
        !     endif
        ! endif
        
    end subroutine
    
    subroutine copy(action,v,bv,iiz,iix,iiy)
        character(4) :: action
        real,dimension(cb%n) :: v
        real,dimension(*) :: bv
        integer,dimension(2),intent(in) :: iiz,iix,iiy
        
        ifz=iiz(1); ilz=iiz(2)
        ifx=iix(1); ilx=iix(2)
        ify=iiy(1); ily=iiy(2)
        
        nnz=iiz(2)-iiz(1)+1
        nnx=iix(2)-iix(1)+1
        nny=iiy(2)-iiy(1)+1
        
        if(action=='save') then !save
            do iy=ify,ily
            do ix=ifx,ilx
            do iz=ifz,ilz
                i = (iz-cb%ifz) + (ix-cb%ifx)*cb%nz + (iy-cb%ify)*cb%nz*cb%nx +1 !field indexing
                k = (iz-ifz)    + (ix-ifx)*nnz      + (iy-ify)*nnz*nnx +1 !boundary_field indexing
                
                bv(k) = v(i)
            enddo
            enddo
            enddo
            
        else !load
            do iy=ify,ily
            do ix=ifx,ilx
            do iz=ifz,ilz
                i = (iz-cb%ifz) + (ix-cb%ifx)*cb%nz + (iy-cb%ify)*cb%nz*cb%nx +1
                k = (iz-ifz)    + (ix-ifx)*nnz      + (iy-ify)*nnz*nnx +1
                
                v(i) = bv(k)
            enddo
            enddo
            enddo
        endif
        
    end subroutine

    subroutine final(self)
        type(t_field) :: self
        
        !deallocate(self%name)

        call dealloc(self%vz,self%vx,self%vy)
        call dealloc(self%szz,self%szx,self%szy,self%sxx,self%sxy,self%shh,self%p)

        call dealloc(self%bnd%vz_top,  self%bnd%vz_bot)
        call dealloc(self%bnd%vx_left, self%bnd%vx_right)
        call dealloc(self%bnd%vy_front,self%bnd%vy_rear)
        call dealloc(self%bnd%vx_top,  self%bnd%vx_bot)
        call dealloc(self%bnd%vz_left, self%bnd%vz_right)

        call dealloc(self%dvz_dz,self%dvz_dx,self%dvz_dy)
        call dealloc(self%dvx_dz,self%dvx_dx,self%dvx_dy)
        call dealloc(self%dvy_dz,self%dvy_dx,self%dvy_dy)
        call dealloc(self%dszz_dz,self%dszx_dz)
        call dealloc(self%dszx_dx,self%dsxx_dx)
        call dealloc(self%dshh_dz,self%dshh_dx,self%dshh_dy)
        call dealloc(self%dp_dz,self%dp_dx,self%dp_dy)

        call dealloc(self%lap)

        call dealloc(self%wavelet)

        call dealloc(self%seismo)

        if(allocated(self%bloom)) deallocate(self%bloom)

        !snapshot
    end subroutine

    
    logical function is_registered(self,chp,str)
        class(t_field) :: self
        type(t_checkpoint) :: chp
        character(*) :: str

        type(t_string),dimension(:),allocatable :: list

        list=split(str)

        do i=1,size(list)
            is_registered=chp%check(self%name//'%'//list(i)%s)
            if(.not.is_registered) return
        enddo

        do i=1,size(list)
            select case (list(i)%s)
            case ('seismo')
                call chp%open(self%name//'%seismo')
                call chp%read(self%seismo)
                call chp%close
                call hud('Read '//self%name//'%seismo from '//chp%name//', size='//num2str(size(self%seismo)))
            case ('comp')
                call chp%open(self%name//'%comp')
                call chp%read(self%vz, self%vx, self%vy )
                call chp%read(self%szz,self%szx,self%szy)
                call chp%read(self%sxx,self%sxy,self%syy)
                call chp%read(self%shh,self%p)
                call chp%close
                call hud('Read '//self%name//'%vz,vx,vy from '//chp%name//', size='//num2str(total_size(self%vz,self%vx,self%vy)))
                call hud('Read '//self%name//'%szz,szx,szy from '//chp%name//', size='//num2str(total_size(self%szz,self%szx,self%szy)))
                call hud('Read '//self%name//'%sxx,sxy,syy from '//chp%name//', size='//num2str(total_size(self%sxx,self%sxy,self%syy)))
                call hud('Read '//self%name//'%shh,p from '//chp%name//', size='//num2str(total_size(self%shh,self%p)))
            case ('boundary')
                call chp%open(self%name//'%boundary')
                call chp%read(self%bnd%vz_top,  self%bnd%vz_bot  )
                call chp%read(self%bnd%vx_left, self%bnd%vx_right)
                call chp%read(self%bnd%vy_front,self%bnd%vy_rear )
                call chp%read(self%bnd%vx_top  ,self%bnd%vx_bot  )
                call chp%read(self%bnd%vz_left ,self%bnd%vz_right)
                call chp%close
                call hud('Read '//self%name//'%boundary vz_top, vz_bot from '//chp%name//', size='//num2str(total_size(self%bnd%vz_top,self%bnd%vz_bot)))
            end select

        enddo

    end function

    subroutine register(self,chp,str)
        class(t_field) :: self
        type(t_checkpoint) :: chp
        character(*) :: str

        type(t_string),dimension(:),allocatable :: list

        list=split(str)

        do i=1,size(list)
            select case (list(i)%s)
            case ('seismo')
                call chp%open(self%name//'%seismo')
                call chp%write(self%seismo)
                call chp%close
            case ('comp')
                call chp%open(self%name//'%comp')
                call chp%write(self%vz, self%vx, self%vy )
                call chp%write(self%szz,self%szx,self%szy)
                call chp%write(self%sxx,self%sxy,self%syy)
                call chp%write(self%shh,self%p)
                call chp%close
            case ('boundary')
                call chp%open(self%name//'%boundary')
                call chp%write(self%bnd%vz_top  ,self%bnd%vz_bot  )
                call chp%write(self%bnd%vx_left ,self%bnd%vx_right)
                call chp%write(self%bnd%vy_front,self%bnd%vy_rear )
                call chp%write(self%bnd%vx_top  ,self%bnd%vx_bot  )
                call chp%write(self%bnd%vz_left ,self%bnd%vz_right)
                call chp%close
            end select

        enddo

    end subroutine

end
