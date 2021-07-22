module m_field
use m_string
use m_mpienv
use m_arrayop
use m_setup
use m_sysio
use m_checkpoint
use m_resampler
use m_model
use m_shot
use m_computebox

use, intrinsic :: ieee_arithmetic

    private
    
    !boundary components for wavefield recontruction
    type t_boundary

        real,dimension(:,:),allocatable :: vz_top,  vz_bot
        real,dimension(:,:),allocatable :: vx_left, vx_right
        real,dimension(:,:),allocatable :: vy_front,vy_rear
        real,dimension(:,:),allocatable :: vx_top,vx_bot,vz_left,vz_right

    end type

    !fields
    type,public :: t_field

        character(:),allocatable :: name

        !wavefield components in computation domain
        real,dimension(:,:,:),allocatable :: vz,vx,vy !velocities
        real,dimension(:,:,:),allocatable :: szz,szx,szy !stress tensor
        real,dimension(:,:,:),allocatable ::     sxx,sxy
        real,dimension(:,:,:),allocatable ::         syy
        real,dimension(:,:,:),allocatable :: shh,p

        !boundary components for wavefield recontruction
        type(t_boundary) :: bnd

        !cpml components for absorbing boundary wavefield
        real,dimension(:,:,:),allocatable :: dvz_dz,dvz_dx,dvz_dy
        real,dimension(:,:,:),allocatable :: dvx_dz,dvx_dx,dvx_dy
        real,dimension(:,:,:),allocatable :: dvy_dz,dvy_dx,dvy_dy
        
        real,dimension(:,:,:),allocatable :: dszz_dz,dszx_dz
        real,dimension(:,:,:),allocatable :: dszx_dx,dsxx_dx        
        real,dimension(:,:,:),allocatable :: dshh_dz,dshh_dx,dshh_dy
        real,dimension(:,:,:),allocatable :: dp_dz,dp_dx,dp_dy

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
        procedure :: init  => init
        procedure :: init_bloom    => init_bloom
        procedure :: init_boundary => init_boundary
        procedure :: check_value => check_value
        procedure :: add_src => add_src
        procedure :: acquire => acquire
        procedure :: write => write
        procedure :: boundary_transport => boundary_transport
        final :: fin

        procedure :: is_registered => is_registered
        procedure :: register => register

    end type

    !info shared by all t_field instances
    integer :: nt
    real :: dt

    logical :: if_shear

    logical :: if_bloom
    integer,parameter :: initial_half_bloomwidth=5 !half bloom width at ift, should involve hicks points

    logical :: if_snapshot
    type(t_string),dimension(:),allocatable :: snapshot
    integer :: i_snapshot
    
    contains
        
    subroutine init(self,name,nt_in,dt_in,if_shear_in)
        class(t_field) :: self
        character(*) :: name
        logical if_shear_in

        self%name=name

        nt=nt_in
        dt=dt_in
        
        if_shear=if_shear_in

        !snapshot
        snapshot=setup%get_strs('SNAPSHOT')
        if_snapshot=size(snapshot)>0 .and. mpiworld%is_master
        if(if_snapshot) then
            n_snapshot=setup%get_int('REF_NUMBER_SNAPSHOT','NSNAPSHOT',o_default='50')
            if(n_snapshot==0) n_snapshot=50
            i_snapshot=ceiling(nt*1./n_snapshot)

            !rm existing snap files
            call execute_command_line('rm '//dir_out//'snap*', wait=.true.)

        endif

        !seismo
        call alloc(self%seismo,shot%nrcv,nt)

    end subroutine

    subroutine init_bloom(self,origin)
        class(t_field) :: self
        character(3) :: origin

        !blooming
        call alloc(self%bloom,6,nt)

        if_bloom=setup%get_bool('IF_BLOOM',o_default='T')      

        if(if_bloom) then
            !directional maximum propagation distance per time step
            distz = cb%velmax * dt / m%dz
            distx = cb%velmax * dt / m%dx
            disty = cb%velmax * dt / m%dy

            if(origin=='src') then
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
    
    subroutine check_value(self)
        class(t_field) :: self
        
        if(mpiworld%is_master) write(*,*) self%name//' sample values:',minval(self%vz),maxval(self%vz)
        
        if(any(.not. ieee_is_finite(self%vz))) then
            call error(self%name//' values become Infinity on Shot# '//shot%sindex//' !!')
        endif
        if(any(ieee_is_nan(self%vz))) then
            call error(self%name//' values become NaN on Shot# '//shot%sindex//' !!')
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
                    case ('vz')
                        call sysio_write('snap_'//self%name//'%vz'//suf,self%vz,size(self%vz),o_mode='append')
                    case ('p')
                        call sysio_write('snap_'//self%name//'%p'//suf,self%p,size(self%p),o_mode='append')
                    end select
                enddo
            endif

        endif

    end subroutine
 
    subroutine add_src(self,ois_adjoint)
        class(t_field) :: self
        logical,optional :: ois_adjoint

        !add adjoint source
        if(either(ois_adjoint,.false.,present(ois_adjoint))) then
            call alloc(self%wavelet,shot%nrcv,nt)
            do i=1,shot%nrcv
                call resampler(shot%dadj(:,i),self%wavelet(i,:),1,din=shot%dt,nin=shot%nt,dout=dt,nout=nt)
            enddo

            return

        endif

        !add source
        call alloc(self%wavelet,1,nt)
        call resampler(shot%wavelet,self%wavelet(1,:),1,din=shot%dt,nin=shot%nt,dout=dt,nout=nt)
        self%wavelet(1,:)=self%wavelet(1,:)*dt/m%cell_volume
            
    end subroutine

    subroutine acquire(self)
        class(t_field) :: self

        call alloc(shot%dsyn,shot%nt,shot%nrcv)
        do i=1,shot%nrcv
            call resampler(self%seismo(i,:), shot%dsyn(:,i),1,din=dt,nin=nt,dout=shot%dt,nout=shot%nt)
        enddo

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
            call copy(action,self%vz,self%bnd%vz_top(:,it),[1,3],    [1,nx],[1,1])
            !bottom
            call copy(action,self%vz,self%bnd%vz_bot(:,it),[nz-1,nz+1],[1,nx],[1,1])
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
                call copy(action,self%vx,self%bnd%vx_top(:,it),[1,3],    [1,nx],[1,1])
                !bottom
                call copy(action,self%vx,self%bnd%vx_bot(:,it),[nz-2,nz  ],[1,nx],[1,1])
                !left
                call copy(action,self%vz,self%bnd%vz_left(:,it), [1,nz],[1,3],    [1,1])
                !right
                call copy(action,self%vz,self%bnd%vz_right(:,it),[1,nz],[nx-2,nx  ],[1,1])
            endif
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
        
        nz=iiz(2)-iiz(1)+1
        nx=iix(2)-iix(1)+1
        ny=iiy(2)-iiy(1)+1
        
        if(action=='save') then !save
            do iy=ify,ily
            do ix=ifx,ilx
            do iz=ifz,ilz
                i = (iz-cb%ifz) + (ix-cb%ifx)*cb%nz + (iy-cb%ify)*cb%nz*cb%nx +1 !field indexing
                k = (iz-ifz)    + (ix-ifx)*nz       + (iy-ify)*nz*nx +1 !boundary_field indexing
                
                bv(k) = v(i)
            enddo
            enddo
            enddo
            
        else !load
            do iy=ify,ily
            do ix=ifx,ilx
            do iz=ifz,ilz
                i = (iz-cb%ifz) + (ix-cb%ifx)*cb%nz + (iy-cb%ify)*cb%nz*cb%nx +1
                k = (iz-ifz)    + (ix-ifx)*nz      + (iy-ify)*nz*nx +1
                
                v(i) = bv(k)
            enddo
            enddo
            enddo
        endif
        
    end subroutine

    subroutine fin(s)
        type(t_field) :: s
        
        deallocate(s%name)

        call dealloc(s%vz,s%vx,s%vy)
        call dealloc(s%szz,s%szx,s%szy,s%sxx,s%sxy,s%shh,s%p)

        call dealloc(s%bnd%vz_top,  s%bnd%vz_bot)
        call dealloc(s%bnd%vx_left, s%bnd%vx_right)
        call dealloc(s%bnd%vy_front,s%bnd%vy_rear)
        call dealloc(s%bnd%vx_top,  s%bnd%vx_bot)
        call dealloc(s%bnd%vz_left, s%bnd%vz_right)
        
        call dealloc(s%dvz_dz,s%dvz_dx,s%dvz_dy)
        call dealloc(s%dvx_dz,s%dvx_dx,s%dvx_dy)
        call dealloc(s%dvy_dz,s%dvy_dx,s%dvy_dy)
        call dealloc(s%dszz_dz,s%dszx_dz)
        call dealloc(s%dszx_dx,s%dsxx_dx)
        call dealloc(s%dshh_dz,s%dshh_dx,s%dshh_dy)
        call dealloc(s%dp_dz,s%dp_dx,s%dp_dy)
        
        call dealloc(s%wavelet,s%seismo)
        deallocate(s%bloom)

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
                if(allocated(self%seismo)) call chp%read(self%seismo,size(self%seismo))
                call chp%close
            case ('comp')
                call chp%open(self%name//'%comp')
                call chp%read(self%vz,size(self%vz),self%vx,size(self%vx),self%vy,size(self%vy))
                call chp%read(self%szz,size(self%szz),self%szx,size(self%szx),self%szy,size(self%szy))
                call chp%read(self%sxx,size(self%sxx),self%sxy,size(self%sxy),self%syy,size(self%syy))
                call chp%read(self%shh,size(self%shh),self%p  ,size(self%p  ))
                call chp%close
            case ('boundary')
                call chp%open(self%name//'%boundary')
                call chp%read(self%bnd%vz_top,  size(self%bnd%vz_top),  self%bnd%vz_bot,  size(self%bnd%vz_bot  ))
                call chp%read(self%bnd%vx_left, size(self%bnd%vx_left), self%bnd%vx_right,size(self%bnd%vx_right))
                call chp%read(self%bnd%vy_front,size(self%bnd%vy_front),self%bnd%vy_rear ,size(self%bnd%vy_rear ))
                call chp%read(self%bnd%vx_top  ,size(self%bnd%vx_top  ),self%bnd%vx_bot  ,size(self%bnd%vx_bot  ))
                call chp%read(self%bnd%vz_left ,size(self%bnd%vz_left ),self%bnd%vz_right,size(self%bnd%vz_right))
                call chp%close
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
                if(allocated(self%seismo)) call chp%write(self%seismo,size(self%seismo))
                call chp%close
            case ('comp')
                call chp%open(self%name//'%comp')
                call chp%write(self%vz,size(self%vz),self%vx,size(self%vx),self%vy,size(self%vy))
                call chp%write(self%szz,size(self%szz),self%szx,size(self%szx),self%szy,size(self%szy))
                call chp%write(self%sxx,size(self%sxx),self%sxy,size(self%sxy),self%syy,size(self%syy))
                call chp%write(self%shh,size(self%shh),self%p  ,size(self%p  ))
                call chp%close
            case ('boundary')
                call chp%open(self%name//'%boundary')
                call chp%write(self%bnd%vz_top  ,size(self%bnd%vz_top  ),self%bnd%vz_bot  ,size(self%bnd%vz_bot  ))
                call chp%write(self%bnd%vx_left ,size(self%bnd%vx_left ),self%bnd%vx_right,size(self%bnd%vx_right))
                call chp%write(self%bnd%vy_front,size(self%bnd%vy_front),self%bnd%vy_rear ,size(self%bnd%vy_rear ))
                call chp%write(self%bnd%vx_top  ,size(self%bnd%vx_top  ),self%bnd%vx_bot  ,size(self%bnd%vx_bot  ))
                call chp%write(self%bnd%vz_left ,size(self%bnd%vz_left ),self%bnd%vz_right,size(self%bnd%vz_right))
                call chp%close
            end select

        enddo

    end subroutine

end