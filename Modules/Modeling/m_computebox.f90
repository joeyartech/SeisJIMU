module m_computebox
use m_sysio
use m_arrayop
use m_model
use m_shot


!geometry of model & computebox
!                                      ifx<---- nx --->ilx
!                                    ifz+~~~~~~~~~~~~~~~+
!       iox     iox+mx-1              ^ l   1      mx   l
!  ++====+-------+===========++       | l  1+-------+   l
!  ||    | *     |           ||       | l   |       |   l
!  || A  |   C   |     B     ||         l D |   C   |   l
!  ||    |       |           ||      nz l   |       |   l
!  ||    |       |           ||         l   |       |   l
!  ++====+-------+===========++       | l mz+-------+   l
!                                     v l boundarylayer l
!       Model: A+B+C                 ilz+~~~~~~~~~~~~~~~+
!       * source point
!                                        Computebox: C+D

    type t_computebox
        
        !C's index in Model
        integer :: iox,ioy,ioz
        
        !C's index in Computebox
        integer :: mx,my,mz
        
        !C+D's index in Computebox
        integer :: ifx,ify,ifz
        integer :: ilx,ily,ilz
        integer :: nx,ny,nz,n
        
        real velmin, velmax


        real,dimension(:,:,:),allocatable :: vp,vs,rho
        real,dimension(:,:,:),allocatable :: eps,del,eta
        real,dimension(:,:,:),allocatable :: qp,qs
        
        real,dimension(:,:,:,:),allocatable :: kernel
        
    end type

    type(t_computebox) :: cb

    contains

    subroutine init
        character(:),allocatable :: tmp
        real :: aperture(4)=0.
        
        !shot aperture, default is whole model
        tmp=setup_get_char('APERTURE',default='-99999 99999 -99999 99999')
        read(tmp,*) aperture
        
    end subroutine
    
    subroutine project

        !C's origin index in model
        cb%ioz=1 !always from top of model
        x=minval(shot%rcv(:)%x)
        x=min(x,shot%src%x,shot%src%x+aperture(1))
        y=minval(shot%rcv(:)%y)
        y=min(y,shot%src%y,shot%src%y+aperture(3))
        cb%iox=max(1,   nint(x/m%dx)+1) !can't exceed size of model
        cb%ioy=max(1,   nint(y/m%dy)+1) !can't exceed size of model
        
        !C's size
        cb%mz=m%nz !always down to bottom of model
        x=maxval(shot%rcv(:)%x)
        x=max(x,shot%src%x,shot%src%x+aperture(2))
        y=maxval(shot%rcv(:)%y)
        y=max(y,shot%src%y,shot%src%y+aperture(4))
        cb%mx=min(m%nx,nint(x/m%dx)+1) +1 -cb%iox
        cb%my=min(m%ny,nint(y/m%dy)+1) +1 -cb%ioy
        
        !C+D's index
        cb%ifx = 1     - cpml%nlayer
        cb%ilx = cb%mx + cpml%nlayer
        cb%ify = 1     - cpml%nlayer
        cb%ily = cb%my + cpml%nlayer
        cb%ifz = 1     - cpml%nlayer
        cb%ilz = cb%mz + cpml%nlayer
        
        !take care of y
        if(.not.m%is_cubic) then
            cb%ioy=1
            cb%ify=1
            cb%ily=1
            cb%my=1
            cb%ny=1
        endif
        
        cb%nz=cb%ilz-cb%ifz+1
        cb%nx=cb%ilx-cb%ifx+1
        cb%ny=cb%ily-cb%ify+1
        cb%n=cb%nz*cb%nx*cb%ny
        if(mpiworld%is_master) then
            write(*,*)'Computebox Size = [ifz,ilz] x [ifx,ilx] x [ify,ily] = ',cb%n
            write(*,*)'  [ifz,ilz],nz:',cb%ifz,cb%ilz,cb%nz
            write(*,*)'  [ifx,ilx],nx:',cb%ifx,cb%ilx,cb%nx
            write(*,*)'  [ify,ilx],ny:',cb%ify,cb%ily,cb%ny
            write(*,*)'Inner area of Computebox:'
            write(*,*)'  ioz,mz:',cb%ioz,cb%mz
            write(*,*)'  iox,mx:',cb%iox,cb%mx
            write(*,*)'  ioy,my:',cb%ioy,cb%my
        endif
        
        !shift source and receiver positions by computebox origin
        call shot_shift_by_computebox(cb%iox,cb%ioy,cb%ioz)
        
        !models in computebox
        if(allocated(m%vp)  call m2cb(m%vp ,cb%vp)
        if(allocated(m%vs)  call m2cb(m%vs ,cb%vs)
        if(allocated(m%rho) call m2cb(m%rho,cb%rho)
        if(allocated(m%eps) call m2cb(m%eps,cb%eps)
        if(allocated(m%del) call m2cb(m%del,cb%del)
        if(allocated(m%eta) call m2cb(m%eta,cb%eta)
        if(allocated(m%qp)  call m2cb(m%qp ,cb%qp)
        if(allocated(m%qs)  call m2cb(m%qs ,cb%qs)
        
        
        cb%velmin=min(minval(cb%vp), minval(cb%vs,cb%vs>0.))
        cb%velmax=maxval(cb%vp*sqrt(1.+2*cb%eps)) !I don't think negative epsilon value can play a role here..
        
        call hud('Computebox value ranges:')
        if(mpiworld%iproc==0) then
                                 write(*,*)'vp' ,minval(cb%vp),maxval(cb%vp)
            if(allocated(cb%vs)  write(*,*)'vs' ,minval(cb%vs),maxval(cb%vs)
                                 write(*,*)'rho',minval(cb%rho),maxval(cb%rho)
                                 write(*,*)'ip' ,minval(cb%vp*cb%rho),maxval(cb%vp*cb%rho)
            if(allocated(cb%eps) write(*,*)'eps',minval(cb%eps),maxval(cb%eps)
            if(allocated(cb%del) write(*,*)'del',minval(cb%del),maxval(cb%del)
        end if
        
        call build_cpml
        
    end subroutine
    
    subroutine m2cb(big,small)
        real,dimension(*),intent(in) :: big
        real,dimension(:,:,:),allocatable,intent(out) :: small
       
        call alloc(small, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        !values inside
        small(1:cb%mz,1:cb%mx,1:cb%my)=big(cb%ioz:cb%ioz+cb%mz-1,cb%iox:cb%iox+cb%mx-1,cb%ioy:cb%ioy+cb%my-1)

        !values in boundary layers
        !!top
        do iz=cb%ifz,0 ;  small(iz,:,:)=small(1,:,:)
        !!bottom
        do iz=cb%mz+1,cb%ilz ; small(iz,:,:)=small(cb%mz,:,:)
        !!left
        do ix=cb%ifx,0 ; small(:,ix,:)=small(:,1,:)
        !!right
        do ix=cb%mx+1,cb%ilx ; small(:,ix,:)=small(:,cb%mx,:)
        !!front
        do iy=cb%ify,0 ; small(:,:,iy)=small(:,:,1)
        !!rear
        do iy=cb%my+1,cb%ily ; small(:,:,iy)=small(:,:,cb%my)

    end subroutine

    subroutine project_back
            gradient(cb%ioz:cb%ioz+cb%mz-1,&
                     cb%iox:cb%iox+cb%mx-1,&
                     cb%ioy:cb%ioy+cb%my-1,:) = &
            gradient(cb%ioz:cb%ioz+cb%mz-1,&
                     cb%iox:cb%iox+cb%mx-1,&
                     cb%ioy:cb%ioy+cb%my-1,:) + cb%kernel(:,:,:,:)
    end subroutine
    
end
