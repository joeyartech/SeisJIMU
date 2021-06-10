module m_field
use m_mpienv
use m_arrayop
use m_sysio
use m_model, only:m
use m_shot, only:shot
use m_computebox, only: cb
use m_FDSG

use, intrinsic :: ieee_arithmetic

    
    !boundary components for wavefield recontruction
    type t_boundary
        real,dimension(:,:),allocatable :: vz_top,vz_bot
        real,dimension(:,:),allocatable :: vx_left,vx_right
        real,dimension(:,:),allocatable :: vy_front,vy_rear
        real,dimension(:,:),allocatable :: vx_top,vx_bot,vz_left,vz_right

    end type

    !cpml components for absorbing boundary wavefield
    type t_cpml
        real,dimension(:,:,:),allocatable :: dvx_dx,dvx_dy,dvx_dz
        real,dimension(:,:,:),allocatable :: dvy_dx,dvy_dy,dvy_dz
        real,dimension(:,:,:),allocatable :: dvz_dx,dvz_dy,dvz_dz
        
        real,dimension(:,:,:),allocatable :: dsxx_dx
        real,dimension(:,:,:),allocatable :: dsxz_dx,        dsxz_dz
        real,dimension(:,:,:),allocatable ::                 dszz_dz
        real,dimension(:,:,:),allocatable :: dshh_dx,dshh_dy,dshh_dz
        real,dimension(:,:,:),allocatable :: dp_dx,  dp_dy,  dp_dz

    end type

    !seismogram components
    type t_seismo
        real,dimension(:,:,:),allocatable :: vx,vy,vz !velocities
        real,dimension(:,:,:),allocatable :: sxx,sxy,sxz !stress tensor
        real,dimension(:,:,:),allocatable :: syx,syy,syz
        real,dimension(:,:,:),allocatable :: szx,szy,szz
        real,dimension(:,:,:),allocatable :: shh,p
    end type

    !fields
    type t_field

        character(:),allocatable :: name

        !wavefield components in computation domain
        real,dimension(:,:,:),allocatable :: vx,vy,vz !velocities
        real,dimension(:,:,:),allocatable :: sxx,sxy,sxz !stress tensor
        real,dimension(:,:,:),allocatable :: syx,syy,syz
        real,dimension(:,:,:),allocatable :: szx,szy,szz
        real,dimension(:,:,:),allocatable :: shh,p
        
        !RHS (source terms)
        real,dimension(:,:),allocatable :: fx,fy,fx !forces
        real,dimension(:,:),allocatable :: mxx,mxy,mxz !moments
        real,dimension(:,:),allocatable :: myx,myy,myz
        real,dimension(:,:),allocatable :: mzx,mzy,mzz
        real,dimension(:,:),allocatable :: mhh,p

        !boundary components for wavefield recontruction
        type(t_boundary) :: bnd

        !cpml components for absorbing boundary wavefield
        type(t_cpml) :: cpml

        !synthetic seismogram components
        type(t_seismo) :: seismo
        
        !source function
        real,dimension(:),allocatable :: source

        !blooming
        integer,dimension(:,:),allocatable :: bloom
        integer,parameter :: initial_half_bloomwidth=5 !half bloom width at ift, should involve hicks points
        
        !snapshot
        type(t_string),dimension(:),allocatable :: snapshot
        integer :: it_delta_snapshot

        contains
        procedure :: init  => init
        procedure :: init_bloom    => init_bloom
        procedure :: init_boundary => init_boundary
        procedure :: check_value => check_value
        procedure :: reinit_boundary => reinit_boundary
        procedure :: write => write
        procedure :: add_source => add_source

        procedure :: propagate => propagate
        procedure :: propagate_reverse => propagate_reverse

    end type
    

    contains
        
    subroutine init(name)
        logical :: if_snapshot

        self%name=name
        
        !snapshot
        self%snapshot=setup%get_strs('SNAPSHOT',default='')
        self%if_snapshot=size(snapshot)>0 .and. mpiworld%is_master
        if(self%if_snapshot) it_snapshot=setup%get_int('SNAPSHOT_DELTA_IT',default=50)

    end subroutine

    subroutine init_bloom(nt)

        !blooming
        call alloc(self%bloom,6,nt)
        
        !directional maximum propagation distance per time step
        distz = cb%velmax * dt / m%dz
        distx = cb%velmax * dt / m%dx
        disty = cb%velmax * dt / m%dy

        if(field%name=='sfield') then
            bloom(1,1)=max(shot%src%iz -initial_half_bloomwidth, cb%ifz)
            bloom(2,1)=min(shot%src%iz +initial_half_bloomwidth, cb%ilz)
            bloom(3,1)=max(shot%src%ix -initial_half_bloomwidth, cb%ifx)
            bloom(4,1)=min(shot%src%ix +initial_half_bloomwidth, cb%ilx)
            bloom(5,1)=max(shot%src%iy -initial_half_bloomwidth, cb%ify)
            bloom(6,1)=min(shot%src%iy +initial_half_bloomwidth, cb%ily)
            do it=2,nt
                bloom(1,it)=max(nint(bloom(1,1)-it*distz),cb%ifz) !bloombox ifz
                bloom(2,it)=min(nint(bloom(2,1)+it*distz),cb%ilz) !bloombox ilz
                bloom(3,it)=max(nint(bloom(3,1)-it*distx),cb%ifx) !bloombox ifx
                bloom(4,it)=min(nint(bloom(4,1)+it*distx),cb%ilx) !bloombox ilx
                bloom(5,it)=max(nint(bloom(5,1)-it*disty),cb%ify) !bloombox ify
                bloom(6,it)=min(nint(bloom(6,1)+it*disty),cb%ily) !bloombox ily
            enddo
        
        else
            bloom(1,nt)=max(minval(shot%rcv(:)%iz) -initial_half_bloomwidth, cb%ifz)
            bloom(2,nt)=min(maxval(shot%rcv(:)%iz) +initial_half_bloomwidth, cb%ilz)
            bloom(3,nt)=max(minval(shot%rcv(:)%ix) -initial_half_bloomwidth, cb%ifx)
            bloom(4,nt)=min(maxval(shot%rcv(:)%ix) +initial_half_bloomwidth, cb%ilx)
            bloom(5,nt)=max(minval(shot%rcv(:)%iy) -initial_half_bloomwidth, cb%ify)
            bloom(6,nt)=min(maxval(shot%rcv(:)%iy) +initial_half_bloomwidth, cb%ily)
            do it=nt-1,1,-1
                it_fwd=nt-it+1
                bloom(1,it)=max(nint(bloom(1,nt)-it_fwd*distz),cb%ifz) !bloombox ifz
                bloom(2,it)=min(nint(bloom(2,nt)+it_fwd*distz),cb%ilz) !bloombox ilz
                bloom(3,it)=max(nint(bloom(3,nt)-it_fwd*distx),cb%ifx) !bloombox ifx
                bloom(4,it)=min(nint(bloom(4,nt)+it_fwd*distx),cb%ilx) !bloombox ilx
                bloom(5,it)=max(nint(bloom(5,nt)-it_fwd*disty),cb%ify) !bloombox ify
                bloom(6,it)=min(nint(bloom(6,nt)+it_fwd*disty),cb%ily) !bloombox ily
            enddo

        endif

        if(.not.m%is_cubic) bloom(5:6,:)=1

    end subroutine

    subroutine init_boundary
        !save 3 grid points, for 4th order FD only
        !different indexing
        n=3*cb%mx*cb%my
        call alloc(bnd%vz_top,n,nt)
        call alloc(bnd%vz_bot,n,nt)
        call alloc(bnd%vx_top,n,nt)
        call alloc(bnd%vx_bot,n,nt)
        
        n=cb%mz*3*cb%my
        call alloc(bnd%vx_left, n,nt)
        call alloc(bnd%vx_right,n,nt)
        call alloc(bnd%vz_left, n,nt)
        call alloc(bnd%vz_right,n,nt)

        if(m%is_cubic) then
            n=cb%mz*cb%mx*3
            call alloc(bnd%vy_front,n,nt)
            call alloc(bnd%vy_rear, n,nt)
        endif

    end subroutine
    
    subroutine check_value(f,name)
        class(t_field), intent(in) :: f
        character(*) :: name
        
        if(mpiworld%is_master) write(*,*) name//' sample values:',minval(f%p),maxval(f%p)
        
        if(any(.not. ieee_is_finite(f%vz))) then
            error(name//' values become Infinity on Shot# '//shot%sindex//' !!')
        endif
        if(any(ieee_is_nan(f%vz))) then
            error(name//' values become NaN on Shot# '//shot%sindex//' !!')
        endif
        
    end subroutine
        
    subroutine snapshot(self)
        do i=1,size(self%snapshot)
            select case (self%snapshot(i)%s)
            case ('vz')
                open(12,file='snap_'//self%name//'_vz',access='stream',position='append')
                write(12) self%vz
                close(12)
            case ('p')
                open(12,file='snap_'//self%name//'_p',access='stream',position='append')
                write(12) self%p
                close(12)
            end select
        enddo
    end subroutine
 
    subroutine add_RHS(wavelet)
        self%wavelet=>wavelet
    end subroutine
    

    subroutine acquire(self)
        type(t_shot) :: self

        call alloc(self%dsyn,self%nt,self%nrcv)
        do i=1,self%nrcv
            select case (self%rcv(i)%comp)
            case ('p')
                call resampler(f%seismo%p(i,:), self%dsyn(:,i),1,din=propagator%dt,nin=propagator%nt,dout=shot%dt,nout=shot%nt)
            case ('vx')
                call resampler(f%seismo%vx(i,:),self%dsyn(:,i),1,din=propagator%dt,nin=propagator%nt,dout=shot%dt,nout=shot%nt)
            case ('vy')
                call resampler(f%seismo%vy(i,:),self%dsyn(:,i),1,din=propagator%dt,nin=propagator%nt,dout=shot%dt,nout=shot%nt)
            case ('vz')
                call resampler(f%seismo%vz(i,:),self%dsyn(:,i),1,din=propagator%dt,nin=propagator%nt,dout=shot%dt,nout=shot%nt)
            end select
        enddo

    end subroutine


    subroutine acquire(f,seismo)
        class(t_field), intent(in) :: f
        real,dimension(*) :: seismo
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            ify=shot%rcv(ircv)%ify; iy=shot%rcv(ircv)%iy; ily=shot%rcv(ircv)%ily
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (1)
                    seismo(ircv)=sum(  f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                    case (2)
                    seismo(ircv)=sum( f%vx(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                    case (3)
                    seismo(ircv)=sum( f%vy(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                    case (4)
                    seismo(ircv)=sum( f%vz(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (1) !p[iz,ix,iy]
                    seismo(ircv)= f%p(iz,ix,iy)
                    case (2) !vx[iz,ix-0.5,iy]
                    seismo(ircv)=f%vx(iz,ix,iy)
                    case (3) !vy[iz,ix,iy-0.5]
                    seismo(ircv)=f%vy(iz,ix,iy)
                    case (4) !vz[iz-0.5,ix,iy]
                    seismo(ircv)=f%vz(iz,ix,iy)
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    subroutine acquire_reverse(f,w)
        class(t_field), intent(in) :: f
        real w
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        ify=shot%src%ify; iy=shot%src%iy; ily=shot%src%ily
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (1)
                w = sum(lm%invkpa(ifz:ilz,ifx:ilx,ify:ily)*f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
                case (2)
                w = sum(     f%vx(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
                case (3)
                w = sum(     f%vy(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
                case (4)
                w = sum(     f%vz(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
            end select
            
        else
            select case (shot%src%icomp)
                case (1) !p[iz,ix,iy]
                w=lm%invkpa(iz,ix,iy)*f%p(iz,ix,iy)
                
                case (2) !vx[iz,ix-0.5,iy]
                w=f%vx(iz,ix,iy)
                
                case (3) !vy[iz,ix,iy-0.5]
                w=f%vy(iz,ix,iy)
                
                case (4) !vz[iz-0.5,ix,iy]
                w=f%vz(iz,ix,iy)
            end select
            
        endif
        
    end subroutine
    
    subroutine boundary_transport(action,it,f)
        character(4) :: action
        integer :: it
        type(t_field) :: f
        
        nz=cb%mz
        nx=cb%mx
        ny=cb%my
        if(m%is_cubic) then
            !top
            call bndcpy(action,f%vz,bnd%vz_top(:,it),[1,3],    [1,nx],[1,ny])  !old version: [0,2],[1,nx],[1,nx]
            !bottom
            call bndcpy(action,f%vz,bnd%vz_bot(:,it),[nz-1,nz+1],[1,nx],[1,ny])  !old version: [nz,nz+2],[1,nx],[1,nx]
            !left
            call bndcpy(action,f%vx,bnd%vx_left(:,it), [1,nz],[1,3],    [1,ny])
            !right
            call bndcpy(action,f%vx,bnd%vx_right(:,it),[1,nz],[nx-1,nx+1],[1,ny])
            !front
            call bndcpy(action,f%vy,bnd%vy_front(:,it),[1,nz],[1,nx],[1,3])
            !rear
            call bndcpy(action,f%vy,bnd%vy_rear(:,it), [1,nz],[1,nx],[ny-1,ny+1])
        else
            !top
            call bndcpy(action,f%vz,bnd%vz_top(:,it),[1,3],    [1,nx],[1,1])
            !bottom
            call bndcpy(action,f%vz,bnd%vz_bot(:,it),[nz-1,nz+1],[1,nx],[1,1])
            !left
            call bndcpy(action,f%vx,bnd%vx_left(:,it), [1,nz],[1,3],    [1,1])
            !right
            call bndcpy(action,f%vx,bnd%vx_right(:,it),[1,nz],[nx-1,nx+1],[1,1])
        endif

        !shear part
        if(present(if_shear)) then
        if(if_shear) then
            if(m%is_cubic) then
            else
                !top
                call bndcpy(action,f%vx,bnd%vx_top(:,it),[1,3],    [1,nx],[1,1])
                !bottom
                call bndcpy(action,f%vx,bnd%vx_bot(:,it),[nz-2,nz  ],[1,nx],[1,1])
                !left
                call bndcpy(action,f%vz,bnd%vz_left(:,it), [1,nz],[1,3],    [1,1])
                !right
                call bndcpy(action,f%vz,bnd%vz_right(:,it),[1,nz],[nx-2,nx  ],[1,1])
            endif
        endif
        endif
        
    end subroutine
    
    subroutine bndcpy(action,v,bv,iiz,iix,iiy)
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
    
end
