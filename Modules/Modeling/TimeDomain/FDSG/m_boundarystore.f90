module m_boundarystore
use m_arrayop
use m_model, only: m
use m_shot, only: shot
use m_computebox, only: cb
use m_field, only: t_field

    private
    public init_boundarystore, boundarystore_transport

    type t_boundarystore
        real,dimension(:,:),allocatable :: vz_top,vz_bot,vx_left,vx_right,vy_front,vy_rear
        real,dimension(:,:),allocatable :: vx_top,vx_bot,vz_left,vz_right
    end type
    
    type(t_boundarystore) :: bnd

    contains
    
    subroutine init_boundarystore
    
        integer::n

        !nt=setup%nt_boundary
        nz=cb%mz
        nx=cb%mx
        ny=cb%my
        nt=shot%rcv(1)%nt
        
        !save 3 grid points, for 4th order FD only
        !different indexing in bfield
        
        n=3*nx*ny
        call alloc(bnd%vz_top,n,nt)
        call alloc(bnd%vx_top,n,nt)
        call alloc(bnd%vz_bot,n,nt)
        call alloc(bnd%vx_bot,n,nt)
        n=nz*3*ny
        call alloc(bnd%vx_left, n,nt)
        call alloc(bnd%vz_left, n,nt)
        call alloc(bnd%vx_right,n,nt)
        call alloc(bnd%vz_right, n,nt)
        if(m%is_cubic) then
            n=nz*nx*3
            call alloc(bnd%vy_front,n,nt)
            call alloc(bnd%vy_rear, n,nt)
        endif
        
    end subroutine
    
    subroutine boundarystore_transport(action,it,f)
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
            call bndcpy(action,f%vz,bnd%vz_bot(:,it),[nz-2,nz],[1,nx],[1,ny])  !old version: [nz,nz+2],[1,nx],[1,nx]
            !left
            call bndcpy(action,f%vx,bnd%vx_left(:,it), [1,nz],[1,3],    [1,ny])
            !right
            call bndcpy(action,f%vx,bnd%vx_right(:,it),[1,nz],[nx-2,nx],[1,ny])
            !front
            call bndcpy(action,f%vy,bnd%vy_front(:,it),[1,nz],[1,nx],[1,3])
            !rear
            call bndcpy(action,f%vy,bnd%vy_rear(:,it), [1,nz],[1,nx],[ny-2,ny])
        else
            !top
            call bndcpy(action,f%vz,bnd%vz_top(:,it),[1,3],    [1,nx],[1,1])
            call bndcpy(action,f%vx,bnd%vx_top(:,it),[1,3],    [1,nx],[1,1])
            !bottom
            call bndcpy(action,f%vz,bnd%vz_bot(:,it),[nz-2,nz],[1,nx],[1,1])
            call bndcpy(action,f%vx,bnd%vx_bot(:,it),[nz-2,nz],[1,nx],[1,1])
            !left
            call bndcpy(action,f%vx,bnd%vx_left(:,it), [1,nz],[1,3],    [1,1])
            call bndcpy(action,f%vz,bnd%vz_left(:,it), [1,nz],[1,3],    [1,1])
            !right
            call bndcpy(action,f%vx,bnd%vx_right(:,it),[1,nz],[nx-2,nx],[1,1])
            call bndcpy(action,f%vz,bnd%vz_right(:,it),[1,nz],[nx-2,nx],[1,1])
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