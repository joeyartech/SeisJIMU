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

        real,dimension(:,:),allocatable :: top_z, top_x, top_y, top_curr, top_next
        real,dimension(:,:),allocatable :: bot_z, bot_x, bot_y, bot_curr, bot_next
        real,dimension(:,:),allocatable :: left_z, left_x, left_y, left_curr, left_next
        real,dimension(:,:),allocatable :: rite_z, rite_x, rite_y, rite_curr, rite_next
        real,dimension(:,:),allocatable :: frnt_z, frnt_x, frnt_y
        real,dimension(:,:),allocatable :: rear_z, rear_x, rear_y

    end type

    !fields
    type,public :: t_field

        character(:),allocatable :: name

        !adjointness
        logical :: is_adjoint
        
        !following the staggered-grid scheme

        !components on the shifted grid
        real,dimension(:,:,:),allocatable :: vz,vx,vy !velocity
        real,dimension(:,:,:),allocatable :: pz,px,py !momenta
        real,dimension(:,:,:),allocatable :: az,ax,ay !acceleration

        real,dimension(:,:,:),pointer :: uz=>null(),      ux=>null()      !displacement
        real,dimension(:,:,:),pointer :: uz_prev=>null(), ux_prev=>null() !previous displacement
        real,dimension(:,:,:),pointer :: uz_next=>null(), ux_next=>null() !future displacement

        !components on the shifted grid
        real,dimension(:,:,:),allocatable :: szx,szy,sxy !shear stress
        real,dimension(:,:,:),allocatable :: ss          !shear stress in 2D
        real,dimension(:,:,:),allocatable :: es          !shear strains in 2D

        !components on the reference grid
        real,dimension(:,:,:),allocatable :: szz,sxx,syy !normal stress
        real,dimension(:,:,:),allocatable :: sz, sx      !normal stress in 2D
        real,dimension(:,:,:),pointer :: p=> null(), p_prev=> null(), p_next=> null() !pressure
        
        real,dimension(:,:,:),allocatable :: ez,ex       !normal strains in 2D

        !laplacian
        real,dimension(:,:,:),allocatable :: lap,lapz,lapx

        ! !Poynting vector
        ! real,dimension(:,:,:),allocatable :: poynz,poynx

        !boundary components for wavefield recontruction
        logical :: if_will_reconstruct=.false.
        logical :: if_boundary_vector=.true.
        type(t_boundary) :: bnd

        !derivatives for absorbing boundary wavefield
        !1st-order
        real,dimension(:,:,:),allocatable :: dz_z, dz_x, dz_y, dz_p, dz_zz, dz_zx
        real,dimension(:,:,:),allocatable :: dx_z, dx_x, dx_y, dx_p, dx_xx, dx_zx
        real,dimension(:,:,:),allocatable :: dy_z, dy_x, dy_y, dy_p
        !2nd-order
        real,dimension(:,:,:),allocatable :: dz_dz_z, dz_dz_x, dz_dz_y, dz_dz_p
        real,dimension(:,:,:),allocatable :: dx_dx_z, dx_dx_x, dx_dx_y, dx_dx_p

        !source wavelet
        real,dimension(:,:),allocatable :: wavelet

        !synthetic seismograms (receiver time function)
        real,dimension(:,:),allocatable :: seismo
        
        !blooming
        integer,dimension(:,:),allocatable :: bloom

        
        contains
        ! procedure :: init
        procedure :: init_bloom
        procedure :: init_boundary
        procedure :: boundary_transport
        procedure :: reinit
        procedure :: check_value
        procedure :: ignite
        procedure :: acquire
        procedure :: write
        procedure :: write_ext
        final :: final
        
        ! procedure :: is_registered
        ! procedure :: register

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
    integer :: i_snapshot, n_snapshot
    
    contains
    
    subroutine field_init(if_shear_in,nt_in,dt_in)
        logical if_shear_in

        logical,save :: is_first_in=.true.

        if_shear=if_shear_in
        nt=nt_in
        dt=dt_in

        !snapshot
        if_snapshot=setup%get_bool('IF_SNAPSHOT') .and. mpiworld%is_master
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
        
        if(self%if_boundary_vector) then
            !save 3 grid points, for 4th order FD only
            !different indexing
            n=3*cb%mx*cb%my
            call alloc(self%bnd%top_z,n,nt)
            call alloc(self%bnd%bot_z,n,nt)
            if(if_shear) then
                call alloc(self%bnd%top_x,n,nt)
                call alloc(self%bnd%bot_x,n,nt)
            endif
            
            n=cb%mz*3*cb%my
            call alloc(self%bnd%left_x, n,nt)
            call alloc(self%bnd%rite_x,n,nt)
            if(if_shear) then
                call alloc(self%bnd%left_z, n,nt)
                call alloc(self%bnd%rite_z,n,nt)
            endif

        return; endif

        !boundary values are scalar
            !save 3 grid points, for 4th order FD only
            !different indexing
            n=3*cb%mx*cb%my
            call alloc(self%bnd%top_curr,n,nt)
            call alloc(self%bnd%top_next,n,nt)
            call alloc(self%bnd%bot_curr,n,nt)
            call alloc(self%bnd%bot_next,n,nt)
            ! if(if_shear) then
            !     call alloc(self%bnd%p_top,n,nt)
            !     call alloc(self%bnd%p_bot,n,nt)
            ! endif

            n=cb%mz*3*cb%my
            call alloc(self%bnd%left_curr, n,nt)
            call alloc(self%bnd%left_next, n,nt)
            call alloc(self%bnd%rite_curr, n,nt)
            call alloc(self%bnd%rite_next, n,nt)
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

!         if(allocated(self%dp_dz))   self%dp_dz=0.
!         if(allocated(self%dp_dx))   self%dp_dx=0.
!         if(allocated(self%dp_dy))   self%dp_dy=0.
!
!         if(allocated(self%dpzz_dz))   self%dpzz_dz=0.
!         if(allocated(self%dpxx_dx))   self%dpxx_dx=0.
!         if(allocated(self%dpyy_dy))   self%dpyy_dy=0.
!
!         if(allocated(self%lap)) self%lap=0.
!
!         if(allocated(self%lapz)) self%lapz=0.
!         if(allocated(self%lapx)) self%lapx=0.
!         if(allocated(self%lapy)) self%lapy=0.

    end subroutine
    
    subroutine check_value(self)
        class(t_field) :: self
        real,dimension(:,:,:),pointer :: a

        if(allocated(self%vz)) then
            a=self%vz
        elseif(associated(self%uz)) then
            a=self%uz
        elseif(allocated(self%pz)) then
            a=self%pz
        elseif(allocated(self%az)) then
            a=self%az
        elseif(associated(self%p)) then
            a=self%p
        endif
                
        if(mpiworld%is_master) write(*,*) self%name//' minmax values:',minval(a),maxval(a)
                
        if(any(.not. ieee_is_finite(a))) then
            call error(self%name//' values become Infinity on '//shot%sindex//' !!')
        endif
        if(any(ieee_is_nan(a))) then
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

                call write_snaps_pointer(self%name,suf,self%uz,'uz', self%ux,'ux')!, self%uy,'uy')
                call write_snaps_pointer(self%name,suf,self%p,'p')

                call write_snaps(self%name,suf,self%vz,'vz', self%vx,'vx', self%vy,'vy')
                call write_snaps(self%name,suf,self%pz,'pz', self%px,'px', self%py,'py')
                call write_snaps(self%name,suf,self%az,'az', self%ax,'ax', self%ay,'ay')

                call write_snaps(self%name,suf,self%szz,'szz', self%sxx,'sxx', self%szx,'szx')
                call write_snaps(self%name,suf,self%sz,'sz', self%sx,'sx', self%ss,'ss')

            endif

        endif

    end subroutine

    subroutine write_snaps(name,suf, a1,c1, o_a2,o_c2, o_a3,o_c3)
        character(*) :: name,suf
        real,dimension(:,:,:),allocatable :: a1, o_a2, o_a3
        character(*) :: c1, o_c2, o_c3

        optional :: o_a2, o_c2, o_a3, o_c3

        if(allocated(a1)) then
            call sysio_write('snap_'//name//'%'//  c1//suf,  a1,cb%n,o_mode='append')
        endif

        if(present(o_a2).and.allocated(o_a2)) then
            call sysio_write('snap_'//name//'%'//o_c2//suf,o_a2,cb%n,o_mode='append')
        endif
        if(present(o_a3).and.allocated(o_a3)) then
            call sysio_write('snap_'//name//'%'//o_c3//suf,o_a3,cb%n,o_mode='append')
        endif

    end subroutine

    subroutine write_snaps_pointer(name,suf, a1,c1, o_a2,o_c2, o_a3,o_c3)
        character(*) :: name,suf
        real,dimension(:,:,:),pointer :: a1, o_a2, o_a3
        character(*) :: c1, o_c2, o_c3

        optional :: o_a2, o_c2, o_a3, o_c3

        if(associated(a1)) then
            call sysio_write('snap_'//name//'%'//  c1//suf,  a1,cb%n,o_mode='append')
        endif

        if(present(o_a2).and.associated(o_a2)) then
            call sysio_write('snap_'//name//'%'//o_c2//suf,o_a2,cb%n,o_mode='append')
        endif
        if(present(o_a3).and.associated(o_a3)) then
            call sysio_write('snap_'//name//'%'//o_c3//suf,o_a3,cb%n,o_mode='append')
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

        nr=either(shot%nrcv,1,self%is_adjoint)
        call alloc(self%wavelet,nr,nt)

        if(present(o_wavelet)) then !w/ wavelet
            if(size(o_wavelet,2)/=nr) then
                call hud('size(o_wavelet,2) vs required size = '//num2str(size(o_wavelet,2))//' vs '//num2str(nr))
                call error('field%ignite: External o_wavelet do NOT have the required shape!')
            endif

            if(size(o_wavelet,1)/=nt) then
                call hud('size(o_wavelet,1) vs required size = '//num2str(size(o_wavelet,1))//' vs '//num2str(nt))
                call warn('field%ignite: Resample o_wavelet to have the required size!')
                 do i=1,nr !implicit transpose
                    call resampler(o_wavelet(:,i),self%wavelet(i,:),1,&
                                    din=shot%dt,nin=size(o_wavelet,1),&
                                    dout=dt,nout=nt)
                enddo

            else
                self%wavelet=transpose(o_wavelet)
            
            endif
            
        return; endif

        !w/o o_wavelet
            if(self%is_adjoint) then
                do i=1,shot%nrcv !implicit transpose
                    call resampler(shot%dadj(:,i),self%wavelet(i,:),1,&
                                    din=shot%dt,nin=shot%nt,&
                                    dout=dt,nout=nt)
                enddo
                
            else
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
        
        !aliases
        real,dimension(:,:,:),pointer :: vz, vx, vy

        nz=cb%mz
        nx=cb%mx
        ny=cb%my

        if(self%if_boundary_vector) then

            if(allocated(self%vz)) then
                                vz=self%vz
                                vx=self%vx
                if(m%is_cubic)  vy=self%vy
            elseif(associated(self%uz)) then
                                vz=self%uz
                                vx=self%ux
                !if(m%is_cubic)  vy=self%uy
            elseif(allocated(self%pz)) then
                                vz=self%pz
                                vx=self%px
                if(m%is_cubic)  vy=self%py
            elseif(allocated(self%az)) then
                                vz=self%az
                                vx=self%ax
                if(m%is_cubic)  vy=self%ay
            endif
        
                call copy(action,vz,self%bnd%top_z (:,it), [1,3],      [1,nx],[1,ny])  !old version: [0,2],[1,nx],[1,nx]
                call copy(action,vz,self%bnd%bot_z (:,it), [nz-1,nz+1],[1,nx],[1,ny])  !old version: [nz,nz+2],[1,nx],[1,nx]
                call copy(action,vx,self%bnd%left_x(:,it), [1,nz],[1,3],      [1,ny])
                call copy(action,vx,self%bnd%rite_x(:,it), [1,nz],[nx-1,nx+1],[1,ny])

            if(m%is_cubic) then
                call copy(action,vy,self%bnd%frnt_y(:,it), [1,nz],[1,nx], [1,3])
                call copy(action,vy,self%bnd%rear_y(:,it), [1,nz],[1,nx], [ny-1,ny+1])
            endif

            !shear part
            if(if_shear) then
                call copy(action,vx,self%bnd%top_x (:,it), [1,3],     [1,nx],[1,ny])
                call copy(action,vx,self%bnd%bot_x (:,it), [nz-2,nz], [1,nx],[1,ny])
                call copy(action,vz,self%bnd%left_z(:,it), [1,nz],[1,3],     [1,ny])
                call copy(action,vz,self%bnd%rite_z(:,it), [1,nz],[nx-2,nx], [1,ny])

            if(m%is_cubic) then
            endif

            endif

        return; endif


        !boundary values are scalar
                call copy(action,self%p,self%bnd%top_curr (:,it), [1,3],      [1,nx],[1,ny])
                call copy(action,self%p,self%bnd%bot_curr (:,it), [nz-1,nz+1],[1,nx],[1,ny])
                call copy(action,self%p,self%bnd%left_curr(:,it), [1,nz],[1,3],      [1,ny])
                call copy(action,self%p,self%bnd%rite_curr(:,it), [1,nz],[nx-1,nx+1],[1,ny])

            if(m%is_cubic) then

            endif
        
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

        !components on the shifted grid
        !call dealloc(self%uz,self%ux,self%uy) !displacement
        call dealloc(self%uz,self%ux) !displacement
        call dealloc(self%vz,self%vx,self%vy) !velocity
        call dealloc(self%pz,self%px,self%py) !momenta
        call dealloc(self%az,self%ax,self%ay) !acceleration

        call dealloc(self%uz,self%ux)           !displacement
        call dealloc(self%uz_prev,self%ux_prev) !previous displacement
        call dealloc(self%uz_next,self%ux_next) !future displacement

        !components on the shifted grid
        call dealloc(self%szx,self%szy,self%sxy) !shear stress
        call dealloc(self%ss) !shear stress in 2D
        call dealloc(self%es) !shear strains in 2D

        !components on the reference grid
        call dealloc(self%szz,self%sxx,self%syy) !normal stress
        call dealloc(self%sz, self%sx, self%ss ) !normal stress in 2D
        call dealloc(self%p, self%p_prev, self%p_next) !pressure
        
        call dealloc(self%ez,self%ex)       !normal strains in 2D

        !laplacian
        call dealloc(self%lap,self%lapz,self%lapx)


        !boundary components
        call dealloc(self%bnd%top_z,  self%bnd%top_x,  self%bnd%top_y,  self%bnd%top_curr,  self%bnd%top_next)
        call dealloc(self%bnd%bot_z,  self%bnd%bot_x,  self%bnd%bot_y,  self%bnd%bot_curr,  self%bnd%bot_next)
        call dealloc(self%bnd%left_z, self%bnd%left_x, self%bnd%left_y, self%bnd%left_curr, self%bnd%left_next)
        call dealloc(self%bnd%rite_z, self%bnd%rite_x, self%bnd%rite_y, self%bnd%rite_curr, self%bnd%rite_next)
        call dealloc(self%bnd%frnt_z, self%bnd%frnt_x, self%bnd%frnt_y)
        call dealloc(self%bnd%rear_z, self%bnd%rear_x, self%bnd%rear_y)

        !derivatives
        call dealloc(self%dz_z, self%dz_x, self%dz_y, self%dz_p)
        call dealloc(self%dx_z, self%dx_x, self%dx_y, self%dx_p)
        call dealloc(self%dy_z, self%dy_x, self%dy_y, self%dy_p)
        call dealloc(self%dz_zz, self%dz_zx)
        call dealloc(self%dx_xx, self%dx_zx)

        !etc
        call dealloc(self%wavelet)

        call dealloc(self%seismo)

        if(allocated(self%bloom)) deallocate(self%bloom)

        !snapshot
    end subroutine

    
    ! logical function is_registered(self,chp,str)
    !     class(t_field) :: self
    !     type(t_checkpoint) :: chp
    !     character(*) :: str

    !     type(t_string),dimension(:),allocatable :: list

    !     list=split(str)

    !     do i=1,size(list)
    !         is_registered=chp%check(self%name//'%'//list(i)%s)
    !         if(.not.is_registered) return
    !     enddo

    !     do i=1,size(list)
    !         select case (list(i)%s)
    !         case ('seismo')
    !             call chp%open(self%name//'%seismo')
    !             call chp%read(self%seismo)
    !             call chp%close
    !             call hud('Read '//self%name//'%seismo from '//chp%name//', size='//num2str(size(self%seismo)))
    !         case ('comp')
    !             call chp%open(self%name//'%comp')
    !             call chp%read(self%vz, self%vx, self%vy )
    !             call chp%read(self%szz,self%szx,self%szy)
    !             call chp%read(self%sxx,self%sxy,self%syy)
    !             call chp%read(self%shh)!,self%p)
    !             call chp%close
    !             call hud('Read '//self%name//'%vz,vx,vy from '//chp%name//', size='//num2str(total_size(self%vz,self%vx,self%vy)))
    !             call hud('Read '//self%name//'%szz,szx,szy from '//chp%name//', size='//num2str(total_size(self%szz,self%szx,self%szy)))
    !             call hud('Read '//self%name//'%sxx,sxy,syy from '//chp%name//', size='//num2str(total_size(self%sxx,self%sxy,self%syy)))
    !             !call hud('Read '//self%name//'%shh,p from '//chp%name//', size='//num2str(total_size(self%shh,self%p)))
    !         case ('boundary')
    !             call chp%open(self%name//'%boundary')
    !             call chp%read(self%bnd%vz_top,  self%bnd%vz_bot  )
    !             call chp%read(self%bnd%vx_left, self%bnd%vx_right)
    !             call chp%read(self%bnd%vy_front,self%bnd%vy_rear )
    !             call chp%read(self%bnd%vx_top  ,self%bnd%vx_bot  )
    !             call chp%read(self%bnd%vz_left ,self%bnd%vz_right)
    !             call chp%close
    !             call hud('Read '//self%name//'%boundary vz_top, vz_bot from '//chp%name//', size='//num2str(total_size(self%bnd%vz_top,self%bnd%vz_bot)))
    !         end select

    !     enddo

    ! end function

    ! subroutine register(self,chp,str)
    !     class(t_field) :: self
    !     type(t_checkpoint) :: chp
    !     character(*) :: str

    !     type(t_string),dimension(:),allocatable :: list

    !     list=split(str)

    !     do i=1,size(list)
    !         select case (list(i)%s)
    !         case ('seismo')
    !             call chp%open(self%name//'%seismo')
    !             call chp%write(self%seismo)
    !             call chp%close
    !         case ('comp')
    !             call chp%open(self%name//'%comp')
    !             call chp%write(self%vz, self%vx, self%vy )
    !             call chp%write(self%szz,self%szx,self%szy)
    !             call chp%write(self%sxx,self%sxy,self%syy)
    !             call chp%write(self%shh)!,self%p)
    !             call chp%close
    !         case ('boundary')
    !             call chp%open(self%name//'%boundary')
    !             call chp%write(self%bnd%vz_top  ,self%bnd%vz_bot  )
    !             call chp%write(self%bnd%vx_left ,self%bnd%vx_right)
    !             call chp%write(self%bnd%vy_front,self%bnd%vy_rear )
    !             call chp%write(self%bnd%vx_top  ,self%bnd%vx_bot  )
    !             call chp%write(self%bnd%vz_left ,self%bnd%vz_right)
    !             call chp%close
    !         end select

    !     enddo

    ! end subroutine

end
