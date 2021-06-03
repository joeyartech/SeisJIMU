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

    !CPML parameter
    real,parameter :: k_x = 1.
    real,parameter :: k_y = 1.
    real,parameter :: k_z = 1.
    real,parameter :: npower = 2.
    real,parameter :: Rcoeff = 0.001  
    real,parameter :: kappa_max=7.!increase this number if absorbing is not satisfactory at grazing incident angle
                                    !(make in/outside PML reflection more separate..)
                                    !decrease this number if grid dispersion is not satisfactory

    type t_computebox
        integer :: ncpml
        
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
        real,dimension(:),allocatable :: b_x,b_x_half,a_x,a_x_half
        real,dimension(:),allocatable :: b_y,b_y_half,a_y,a_y_half
        real,dimension(:),allocatable :: b_z,b_z_half,a_z,a_z_half

        real,dimension(:),allocatable :: kappa_x,kappa_x_half
        real,dimension(:),allocatable :: kappa_y,kappa_y_half
        real,dimension(:),allocatable :: kappa_z,kappa_z_half

        real,dimension(:,:,:,:),allocatable :: gradient,image
        
    end type
    
    type(t_computebox) :: cb   

    contains
    
    subroutine project_onto
        character(:),allocatable :: tmp
        real :: aperture(4)=0.
        
        cb%ncpml=setup_get_int('CPML_SIZE','NCPML')+2  !for FDTDo4 and hicks interpolation of radius=2
        !cb%ncpml=setup_get_int('CPML_SIZE','NCPML')+4  !for FDTDo8
        
        !shot aperture, default is whole model
        tmp=setup_get_char('APERTURE',default='-99999 99999 -99999 99999')
        read(tmp,*) aperture
        
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
        cb%ifx = 1     - cb%ncpml
        cb%ilx = cb%mx + cb%ncpml
        cb%ify = 1     - cb%ncpml
        cb%ily = cb%my + cb%ncpml
        cb%ifz = 1     - cb%ncpml
        cb%ilz = cb%mz + cb%ncpml
        
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

    subroutine project_back
            gradient(cb%ioz:cb%ioz+cb%mz-1,&
                     cb%iox:cb%iox+cb%mx-1,&
                     cb%ioy:cb%ioy+cb%my-1,:) = &
            gradient(cb%ioz:cb%ioz+cb%mz-1,&
                     cb%iox:cb%iox+cb%mx-1,&
                     cb%ioy:cb%ioy+cb%my-1,:) + cb%gradient(:,:,:,:)
    end subroutine
    
    !CPML implementation from Komatitsch & Martin, 2007, Geophysics
    !code from Geodynamics or Komatitsch's GitHub site:
    !https://github.com/geodynamics/seismic_cpml
    subroutine build_cpml

        real,dimension(:),allocatable   :: alpha_x,alpha_x_half
        real,dimension(:),allocatable   :: alpha_y,alpha_y_half
        real,dimension(:),allocatable   :: alpha_z,alpha_z_half
        real,dimension(:),allocatable   :: d_x,d_x_half
        real,dimension(:),allocatable   :: d_y,d_y_half
        real,dimension(:),allocatable   :: d_z,d_z_half
                
        alpha_max = r_pi*shot%src%fpeak

        thickness_pml_x = cb%ncpml * m%dx
        thickness_pml_y = cb%ncpml * m%dy
        thickness_pml_z = cb%ncpml * m%dz
        
        d0_x = -(npower+1)/(2.*thickness_pml_x) * cb%velmax * log(Rcoeff)
        d0_y = -(npower+1)/(2.*thickness_pml_y) * cb%velmax * log(Rcoeff) 
        d0_z = -(npower+1)/(2.*thickness_pml_z) * cb%velmax * log(Rcoeff)

        call alloc(alpha_x,[cb%ifx,cb%ilx]); call alloc(alpha_x_half,[cb%ifx,cb%ilx])
        call alloc(alpha_y,[cb%ify,cb%ily]); call alloc(alpha_y_half,[cb%ify,cb%ily])
        call alloc(alpha_z,[cb%ifz,cb%ilz]); call alloc(alpha_z_half,[cb%ifz,cb%ilz])

        call alloc(d_x,[cb%ifx,cb%ilx]); call alloc(d_x_half,[cb%ifx,cb%ilx])
        call alloc(d_y,[cb%ify,cb%ily]); call alloc(d_y_half,[cb%ify,cb%ily])
        call alloc(d_z,[cb%ifz,cb%ilz]); call alloc(d_z_half,[cb%ifz,cb%ilz])

        call alloc(cb%b_x,[cb%ifx,cb%ilx]); call alloc(cb%b_x_half,[cb%ifx,cb%ilx])
        call alloc(cb%a_x,[cb%ifx,cb%ilx]); call alloc(cb%a_x_half,[cb%ifx,cb%ilx])
        call alloc(cb%b_y,[cb%ify,cb%ily]); call alloc(cb%b_y_half,[cb%ify,cb%ily])
        call alloc(cb%a_y,[cb%ify,cb%ily]); call alloc(cb%a_y_half,[cb%ify,cb%ily])
        call alloc(cb%b_z,[cb%ifz,cb%ilz]); call alloc(cb%b_z_half,[cb%ifz,cb%ilz])
        call alloc(cb%a_z,[cb%ifz,cb%ilz]); call alloc(cb%a_z_half,[cb%ifz,cb%ilz])

        call alloc(cb%kappa_x,[cb%ifx,cb%ilx],initialize=.false.); call alloc(cb%kappa_x_half,[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(cb%kappa_y,[cb%ify,cb%ily],initialize=.false.); call alloc(cb%kappa_y_half,[cb%ify,cb%ily],initialize=.false.)
        call alloc(cb%kappa_z,[cb%ifz,cb%ilz],initialize=.false.); call alloc(cb%kappa_z_half,[cb%ifz,cb%ilz],initialize=.false.)

        cb%kappa_x=1.; cb%kappa_x_half=1.
        cb%kappa_y=1.; cb%kappa_y_half=1.
        cb%kappa_z=1.; cb%kappa_z_half=1.


        !x dir
        do i = cb%ifx,cb%ilx

            !!left edge
            abscissa_in_pml = -(i-1)*m%dx
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x(i)          = d0_x * abscissa_normalized**npower
                cb%kappa_x(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_x(i)      = alpha_max * (1. - abscissa_normalized)
            endif

            !!left edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dx + m%dx/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x_half(i)     = d0_x * abscissa_normalized**npower
                cb%kappa_x_half(i) = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_x_half(i) = alpha_max * (1. - abscissa_normalized)
            endif

            !!right edge
            abscissa_in_pml = (i-cb%mx)*m%dx
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x(i)          = d0_x * abscissa_normalized**npower
                cb%kappa_x(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_x(i)      = alpha_max * (1.- abscissa_normalized)
            endif

            !!right edge half gridpoint
            abscissa_in_pml = (i-cb%mx)*m%dx - m%dx/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_x
                d_x_half(i)     = d0_x * abscissa_normalized**npower
                cb%kappa_x_half(i) = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_x_half(i) = alpha_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alpha_x(i)     <0.) alpha_x(i)=0.
            if(alpha_x_half(i)<0.) alpha_x_half(i)=0.

            cb%b_x(i)      = exp(-(d_x(i)     /cb%kappa_x(i)     +alpha_x(i)     )*shot%src%dt)
            cb%b_x_half(i) = exp(-(d_x_half(i)/cb%kappa_x_half(i)+alpha_x_half(i))*shot%src%dt)

            if(abs(d_x(i)     )>1e-6) cb%a_x(i)     =d_x(i)     *(cb%b_x(i)     -1.)/(cb%kappa_x(i)     *(d_x(i)     +cb%kappa_x(i)     *alpha_x(i)     ))
            if(abs(d_x_half(i))>1e-6) cb%a_x_half(i)=d_x_half(i)*(cb%b_x_half(i)-1.)/(cb%kappa_x_half(i)*(d_x_half(i)+cb%kappa_x_half(i)*alpha_x_half(i)))

        enddo

        cb%kappa_x     =1./cb%kappa_x
        cb%kappa_x_half=1./cb%kappa_x_half


        !y dir
        do i = cb%ify,cb%ily

            !!front edge
            abscissa_in_pml = -(i-1)*m%dy
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y(i)          = d0_y * abscissa_normalized**npower
                cb%kappa_y(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_y(i)      = alpha_max * (1. - abscissa_normalized)
            endif

            !!front edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dy + m%dy/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y_half(i)     = d0_y * abscissa_normalized**npower
                cb%kappa_y_half(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_y_half(i) = alpha_max * (1. - abscissa_normalized)
            endif

            !!rear edge
            abscissa_in_pml = (i-cb%my)*m%dy
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y(i)          = d0_y * abscissa_normalized**npower
                cb%kappa_y(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_y(i)      = alpha_max * (1.- abscissa_normalized)
            endif

            !!rear edge half gridpoint
            abscissa_in_pml = (i-cb%my)*m%dy - m%dy/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_y
                d_y_half(i)     = d0_y * abscissa_normalized**npower
                cb%kappa_y_half(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_y_half(i) = alpha_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alpha_y(i)     <0.) alpha_y(i)=0.
            if(alpha_y_half(i)<0.) alpha_y_half(i)=0.

            cb%b_y(i)      = exp(-(d_y(i)     /cb%kappa_y(i)     +alpha_y(i)     )*shot%src%dt)
            cb%b_y_half(i) = exp(-(d_y_half(i)/cb%kappa_y_half(i)+alpha_y_half(i))*shot%src%dt)

            if(abs(d_y(i)     )>1e-6) cb%a_y(i)     =d_y(i)     *(cb%b_y(i)     -1.)/(cb%kappa_y(i)     *(d_y(i)     +cb%kappa_y(i)     *alpha_y(i)     ))
            if(abs(d_y_half(i))>1e-6) cb%a_y_half(i)=d_y_half(i)*(cb%b_y_half(i)-1.)/(cb%kappa_y_half(i)*(d_y_half(i)+cb%kappa_y_half(i)*alpha_y_half(i)))

        enddo

        cb%kappa_y     =1./cb%kappa_y
        cb%kappa_y_half=1./cb%kappa_y_half

        !z dir
        do i = cb%ifz,cb%ilz

            !!top edge
            abscissa_in_pml = -(i-1)*m%dz
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z(i)          = d0_z * abscissa_normalized**npower
                cb%kappa_z(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_z(i)      = alpha_max * (1. - abscissa_normalized)
            endif

            !!top edge half gridpoint
            abscissa_in_pml = -(i-1)*m%dz + m%dz/2.
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z_half(i)     = d0_z * abscissa_normalized**npower
                cb%kappa_z_half(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_z_half(i) = alpha_max * (1. - abscissa_normalized)
            endif

            !!bottom edge
            abscissa_in_pml = (i-cb%mz)*m%dz
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z(i)          = d0_z * abscissa_normalized**npower
                cb%kappa_z(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_z(i)      = alpha_max * (1.- abscissa_normalized)
            endif

            !!bottom edge half gridpoint
            abscissa_in_pml = (i-cb%mz)*m%dz - m%dz/2
            if(abscissa_in_pml >= 0.)then
                abscissa_normalized = abscissa_in_pml / thickness_pml_z
                d_z_half(i)     = d0_z * abscissa_normalized**npower
                cb%kappa_z_half(i)   = 1. + (kappa_max-1.)*abscissa_normalized**npower
                alpha_z_half(i) = alpha_max * (1.- abscissa_normalized)
            endif

            !CPML parameters
            if(alpha_z(i)     <0.) alpha_z(i)=0.
            if(alpha_z_half(i)<0.) alpha_z_half(i)=0.

            cb%b_z(i)      = exp(-(d_z(i)     /cb%kappa_z(i)     +alpha_z(i)     )*shot%src%dt)
            cb%b_z_half(i) = exp(-(d_z_half(i)/cb%kappa_z_half(i)+alpha_z_half(i))*shot%src%dt)

            if(abs(d_z(i)     )>1e-6) cb%a_z(i)     =d_z(i)     *(cb%b_z(i)     -1.)/(cb%kappa_z(i)     *(d_z(i)     +cb%kappa_z(i)     *alpha_z(i)     ))
            if(abs(d_z_half(i))>1e-6) cb%a_z_half(i)=d_z_half(i)*(cb%b_z_half(i)-1.)/(cb%kappa_z_half(i)*(d_z_half(i)+cb%kappa_z_half(i)*alpha_z_half(i)))

        enddo

        cb%kappa_z     =1./cb%kappa_z
        cb%kappa_z_half=1./cb%kappa_z_half

        
        !no more need these variables
        deallocate(alpha_x,alpha_x_half)
        deallocate(alpha_y,alpha_y_half)
        deallocate(alpha_z,alpha_z_half)
        deallocate(d_x,d_x_half)
        deallocate(d_y,d_y_half)
        deallocate(d_z,d_z_half)
        
    end subroutine
    
end
