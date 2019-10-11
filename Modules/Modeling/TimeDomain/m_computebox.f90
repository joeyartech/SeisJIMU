module m_computebox
use m_sysio
use m_arrayop
use m_model, only: m
use m_shot

    private
    public cb, build_computebox
    
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
        integer :: nboundarylayer
        
        !C's index in Model
        integer :: iox,ioy,ioz
        
        !C's index in Computebox
        integer :: mx,my,mz
        
        !C+D's index in Computebox
        integer :: ifx,ify,ifz
        integer :: ilx,ily,ilz
        integer :: nx,ny,nz,n
        
        real velmin, velmax
        
        real,dimension(:,:,:),allocatable :: vp,vs,rho,eps,del
        real,dimension(:),allocatable :: b_x,b_y,b_z,a_x,a_y,a_z
        real,dimension(:,:,:,:),allocatable :: gradient,image
        
    end type
    
    type(t_computebox) :: cb
    
    contains
    
    subroutine build_computebox
        character(:),allocatable :: c_aperture
        real :: aperture(4)=0.
        
        nboundarylayer=get_setup_int('NBOUNDARYLAYER')
        
        cb%nboundarylayer=nboundarylayer+2  !for FDTDo4 and hicks interpolation of radius=2
        !cb%nboundarylayer=setup%nboundarylayer+4  !for FDTDo8
        
        !shot aperture, default is whole model
        c_aperture=get_setup_char('APERTURE',default='-99999 99999 -99999 99999')
        read(c_aperture,*) aperture
        
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
        cb%ifx = 1     - cb%nboundarylayer
        cb%ilx = cb%mx + cb%nboundarylayer
        cb%ify = 1     - cb%nboundarylayer
        cb%ily = cb%my + cb%nboundarylayer
        cb%ifz = 1     - cb%nboundarylayer
        cb%ilz = cb%mz + cb%nboundarylayer
        
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
        call alloc(cb%vp, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(cb%vs, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(cb%rho,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(cb%eps,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        call alloc(cb%del,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        
        !model -> computebox
        call m2cb(shape(m%vp), m%vp, cb%vp)
        call m2cb(shape(m%vs), m%vs, cb%vs)
        call m2cb(shape(m%rho),m%rho,cb%rho)
        call m2cb(shape(m%eps),m%eps,cb%eps)
        call m2cb(shape(m%del),m%del,cb%del)
        
        !values in boundary layers
        !!top
        do iz=cb%ifz,0
            cb%vp (iz,:,:)=cb%vp (1,:,:)
            cb%vs (iz,:,:)=cb%vs (1,:,:)
            cb%rho(iz,:,:)=cb%rho(1,:,:)
            cb%eps(iz,:,:)=cb%eps(1,:,:)
            cb%del(iz,:,:)=cb%del(1,:,:)
        enddo
        !!bottom
        do iz=cb%mz+1,cb%ilz
            cb%vp (iz,:,:)=cb%vp (cb%mz,:,:)
            cb%vs (iz,:,:)=cb%vs (cb%mz,:,:)
            cb%rho(iz,:,:)=cb%rho(cb%mz,:,:)
            cb%eps(iz,:,:)=cb%eps(cb%mz,:,:)
            cb%del(iz,:,:)=cb%del(cb%mz,:,:)
        enddo
        !!left
        do ix=cb%ifx,0
            cb%vp (:,ix,:)=cb%vp (:,1,:)
            cb%vs (:,ix,:)=cb%vs (:,1,:)
            cb%rho(:,ix,:)=cb%rho(:,1,:)
            cb%eps(:,ix,:)=cb%eps(:,1,:)
            cb%del(:,ix,:)=cb%del(:,1,:)
        enddo
        !!right
        do ix=cb%mx+1,cb%ilx
            cb%vp (:,ix,:)=cb%vp (:,cb%mx,:)
            cb%vs (:,ix,:)=cb%vs (:,cb%mx,:)
            cb%rho(:,ix,:)=cb%rho(:,cb%mx,:)
            cb%eps(:,ix,:)=cb%eps(:,cb%mx,:)
            cb%del(:,ix,:)=cb%del(:,cb%mx,:)
        enddo
        !!front
        do iy=cb%ify,0
            cb%vp (:,:,iy)=cb%vp (:,:,1)
            cb%vs (:,:,iy)=cb%vs (:,:,1)
            cb%rho(:,:,iy)=cb%rho(:,:,1)
            cb%eps(:,:,iy)=cb%eps(:,:,1)
            cb%del(:,:,iy)=cb%del(:,:,1)
        enddo
        !!rear
        do iy=cb%my+1,cb%ily
            cb%vp (:,:,iy)=cb%vp (:,:,cb%my)
            cb%vs (:,:,iy)=cb%vs (:,:,cb%my)
            cb%rho(:,:,iy)=cb%rho(:,:,cb%my)
            cb%eps(:,:,iy)=cb%eps(:,:,cb%my)
            cb%del(:,:,iy)=cb%del(:,:,cb%my)
        enddo
        
        cb%velmin=min(minval(cb%vp), minval(cb%vs,cb%vs>0.))
        cb%velmax=maxval(cb%vp*sqrt(1.+2*cb%eps)) !I don't think negative epsilon value can play a role here..
        
        call hud('Computebox value ranges:')
        if(mpiworld%iproc==0) then
            write(*,*)'vp' ,minval(cb%vp),maxval(cb%vp)
            write(*,*)'vs' ,minval(cb%vs),maxval(cb%vs)
            write(*,*)'rho',minval(cb%rho),maxval(cb%rho)
            write(*,*)'ip' ,minval(cb%vp*cb%rho),maxval(cb%vp*cb%rho)
            write(*,*)'eps',minval(cb%eps),maxval(cb%eps)
            write(*,*)'del',minval(cb%del),maxval(cb%del)
        end if
        
        call build_cpml
        
    end subroutine
    
    subroutine m2cb(n,large,little)
        integer :: n(3)
        real :: large(n(1),n(2),n(3))
        real :: little(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily)
        if(n(1)*n(2)*n(3)==1) then
            little(1:cb%mz,1:cb%mx,1:cb%my)=large(1,1,1)
        else
            little(1:cb%mz,1:cb%mx,1:cb%my)=large(cb%ioz:cb%ioz+cb%mz-1,cb%iox:cb%iox+cb%mx-1,cb%ioy:cb%ioy+cb%my-1)
        endif
    end subroutine
    
    subroutine build_cpml
        real,parameter :: k_x = 1.
        real,parameter :: k_y = 1.
        real,parameter :: k_z = 1.
        real,parameter :: npower=2.
        real,parameter :: Rcoef=0.0001
    
        real                            :: thickness_pml_x,thickness_pml_y,thickness_pml_z
        real                            :: d0_x,d0_y,d0_z
        real,dimension(:),allocatable   :: alpha_x,alpha_y,alpha_z
        real,dimension(:),allocatable   :: d_x,d_y,d_z
        real                            :: xoriginleft,xoriginright
        real                            :: yoriginleft,yoriginright
        real                            :: zoriginbottom,zorigintop
        real                            :: abscissa_in_pml,abscissa_normalized  
        real                            :: xval,yval,zval
        real                            :: alpha_max_pml
        
        alpha_max_pml = 3.1415927*shot%src%fpeak
!         if(mpiworld%is_master)  write(*,*) 'Coeff alpha_max in CPML',alpha_max_pml
        !---thickness of the pml layer in meters
        thickness_pml_x = cb%nboundarylayer * m%dx
        thickness_pml_y = cb%nboundarylayer * m%dy
        thickness_pml_z = cb%nboundarylayer * m%dz
        
        d0_x = -(npower+1)/(2.*thickness_pml_x) * cb%velmax * log(Rcoef)  ! ref or max vp?
        d0_y = -(npower+1)/(2.*thickness_pml_y) * cb%velmax * log(Rcoef) 
        d0_z = -(npower+1)/(2.*thickness_pml_z) * cb%velmax * log(Rcoef)    
!         if(mpiworld%is_master)  write(*,*) 'Coeff d0 in CPML',d0_z,d0_y,d0_x

        call alloc(alpha_x,[cb%ifx,cb%ilx])
        call alloc(alpha_y,[cb%ify,cb%ily])
        call alloc(alpha_z,[cb%ifz,cb%ilz])
        
        call alloc(d_x,[cb%ifx,cb%ilx])
        call alloc(d_y,[cb%ify,cb%ily])
        call alloc(d_z,[cb%ifz,cb%ilz])

        call alloc(cb%b_x,[cb%ifx,cb%ilx])
        call alloc(cb%a_x,[cb%ifx,cb%ilx])
        call alloc(cb%b_y,[cb%ify,cb%ily])
        call alloc(cb%a_y,[cb%ify,cb%ily])
        call alloc(cb%b_z,[cb%ifz,cb%ilz])
        call alloc(cb%a_z,[cb%ifz,cb%ilz])

        !---damping in the x direction
        ! origin of the pml layer (position of right edge minus thickness, in meters)
        xoriginleft = 0.
        xoriginright = (cb%mx-1)*m%dx

        do i = cb%ifx,cb%ilx
            ! abscissa of current grid point along the damping profile
            xval = m%dx*(i-1)

            !---left edge
            !define damping profile at the grid points
            abscissa_in_pml = xoriginleft - xval
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x(i) = d0_x * abscissa_normalized**npower
                alpha_x(i) = alpha_max_pml * (1. - abscissa_normalized)
            endif

            !---right edge
            !define damping profile at the grid points
            abscissa_in_pml = xval - xoriginright
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x(i) = d0_x * abscissa_normalized**npower
                alpha_x(i) = alpha_max_pml * (1.- abscissa_normalized)
            endif

            !just in case ?
            if(alpha_x(i) < 0.) alpha_x(i) = 0.

            !c-pml parameters
            cb%b_x(i) = exp(- (d_x(i) / k_x + alpha_x(i)) * shot%src%dt)

            !this to avoid division by zero outside the pml
            if(abs(d_x(i)) > 1.e-6) cb%a_x(i) = d_x(i) * (cb%b_x(i) - 1.) / (k_x * (d_x(i) + k_x * alpha_x(i)))

        enddo

        !---damping in the y direction
        ! origin of the pml layer (position of right edge minus thickness, in meters)
        yoriginleft = 0.
        yoriginright = (cb%my-1)*m%dy

        do i = cb%ify,cb%ily
            ! abscissa of current grid point along the damping profile
            yval =  m%dy*(i-1)

            !---left edge
            !define damping profile at the grid points
            abscissa_in_pml = yoriginleft - yval
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y(i) = d0_y * abscissa_normalized**npower
                alpha_y(i) = alpha_max_pml * (1. - abscissa_normalized)
            !  cb%icpml(:,:,i)=.true.
            endif

            !---right edge
            !define damping profile at the grid points
            abscissa_in_pml = yval - yoriginright
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y(i) = d0_y * abscissa_normalized**npower
                alpha_y(i) = alpha_max_pml * (1.- abscissa_normalized)
                ! cb%icpml(:,:,i)=.true.
            endif

            !just in case ?
            if(alpha_y(i) < 0.) alpha_y(i) = 0.

            !c-pml parameters
            cb%b_y(i) = exp(- (d_y(i) / k_y + alpha_y(i)) * shot%src%dt)

            !this to avoid division by zero outside the pml
            if(abs(d_y(i)) > 1.e-6) cb%a_y(i) = d_y(i) * (cb%b_y(i) - 1.) / (k_y * (d_y(i) + k_y * alpha_y(i)))

        enddo

        !---damping in the z direction
        ! origin of the pml layer (position of right edge minus thickness, in meters)
        zorigintop = 0.
        zoriginbottom = (cb%mz-1)*m%dz

        do i = cb%ifz,cb%ilz
            ! abscissa of current grid point along the damping profile
            zval = m%dz*(i-1)
            
            !---top edge
            !define damping profile at the grid points
            abscissa_in_pml = zorigintop - zval
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z(i) = d0_z * abscissa_normalized**npower
                alpha_z(i) = alpha_max_pml * (1. - abscissa_normalized)
                ! cb%icpml(i,:,:)=.true.
            endif
            
            !---bottom edge
            !define damping profile at the grid points
            abscissa_in_pml = zval - zoriginbottom
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z(i) = d0_z * abscissa_normalized**npower
                alpha_z(i) = alpha_max_pml * (1. - abscissa_normalized)
            !  cb%icpml(i,:,:)=.true.
            endif
            
            !just in case ?
            if(alpha_z(i) < 0.) alpha_z(i) = 0.
            
            !c-pml parameters
            cb%b_z(i) = exp(- (d_z(i) / k_z + alpha_z(i)) * shot%src%dt)
            
            !this to avoid division by zero outside the pml
            if(abs(d_z(i)) > 1.e-6) cb%a_z(i) = d_z(i) * (cb%b_z(i) - 1.) / (k_z * (d_z(i) + k_z * alpha_z(i)))
            
        enddo
        
        !no more need these variables
        deallocate(alpha_x)
        deallocate(alpha_y)
        deallocate(alpha_z)
        deallocate(d_x)
        deallocate(d_y)
        deallocate(d_z)
        
    end subroutine
    
end
