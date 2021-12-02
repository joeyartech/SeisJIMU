module m_computebox
use m_System
use m_model
use m_shot

    private

!geometry of model & computebox
!                                      ifx<---- nx ---->ilx
!                                    ifz+~~~~~~~~~~~~~~~~+
!       iox     iox+mx-1              ^ l   1      mx    l
!  ++====+-------+===========++       | l  1+-------+    l
!  ||    | *     |           ||       | l   |       |    l
!  || A  |   C   |     B     ||         l D |   C   |    l
!  ||    |       |           ||      nz l   |       |    l
!  ||    |       |           ||         l   |       |    l
!  ++====+-------+===========++       | l mz+-------+    l
!                                     v l absorbinglayer l
!       Model: A+B+C                 ilz+~~~~~~~~~~~~~~~~+
!       * source point
!                                        Computebox: C+D

    real,dimension(4) :: add_aperture
    
    type,public :: t_computebox

        integer :: nabslayer !thickness of D

        !C's index in Model
        integer :: ioz,iox,ioy
        
        !C's index in Computebox
        integer :: mz,mx,my !=m%nz,nx,ny
        
        !C+D's index in Computebox
        integer :: ifz,ifx,ify
        integer :: ilz,ilx,ily
        integer :: nz,nx,ny,n !>m%nz,nx,ny,n

        ! real :: cell_volume, cell_diagonal, cell_inv_diagonal

        real velmin, velmax

        real,dimension(:,:,:),allocatable :: vp,vs,rho
        real,dimension(:,:,:),allocatable :: eps,del,eta
        real,dimension(:,:,:),allocatable :: qp,qs
        
        real,dimension(:,:,:,:),allocatable :: grad, imag, engy, corr
        
        contains
        procedure :: init
        procedure :: project
        procedure :: project_back
        final :: final

        procedure :: is_registered
        procedure :: register
    end type

    type(t_computebox),public :: cb

    contains

    subroutine init(self,add_abslayer)
        class(t_computebox) :: self
        integer :: add_abslayer

        !add aperture; aperture default to whole model
        add_aperture=setup%get_reals('ADD_APERTURE',o_default='-99999 99999 -99999 99999')

        ! self%cell_volume = m%dz*m%dx*m%dy

        ! if(m%is_cubic) then
        !     self%cell_diagonal=sqrt(m%dz**2+m%dx**2+m%dy**2)
        !     self%cell_inv_diagonal=sqrt(m%dz**(-2) + m%dx**(-2) + m%dy**(-2))
        ! else
        !     self%cell_diagonal=sqrt(m%dz**2+m%dx**2)
        !     self%cell_inv_diagonal=sqrt(m%dz**(-2) + m%dx**(-2))
        ! endif

        !thickness of D = 
        !thickness of user given + thickness required by propagator
        self%nabslayer=setup%get_int('REF_BNDLAYER_THICKNESS','NCPML',o_default='20')+add_abslayer 

    end subroutine
    
    subroutine project(self)
        class(t_computebox) :: self

        !C's origin index in model
        self%ioz=1 !always from top of model

        x=min( shot%src%x, minval(shot%rcv(:)%x)) +add_aperture(1)
        y=min( shot%src%y, minval(shot%rcv(:)%y)) +add_aperture(3)

        self%iox=max(1,   nint(x/m%dx)+1) !can't exceed size of model
        self%ioy=max(1,   nint(y/m%dy)+1) !can't exceed size of model
        
        !C's size
        self%mz=m%nz !always down to bottom of model

        x=max( shot%src%x, maxval(shot%rcv(:)%x)) +add_aperture(2)
        y=max( shot%src%y, maxval(shot%rcv(:)%y)) +add_aperture(4)
        
        self%mx=min(m%nx,nint(x/m%dx)+1) -self%iox +1
        self%my=min(m%ny,nint(y/m%dy)+1) -self%ioy +1
        
        !C+D's index
        self%ifx = 1       - self%nabslayer
        self%ilx = self%mx + self%nabslayer
        self%ify = 1       - self%nabslayer
        self%ily = self%my + self%nabslayer
        self%ifz = 1       - self%nabslayer
        self%ilz = self%mz + self%nabslayer
        
        !take care of y
        if(.not.m%is_cubic) then
            self%ioy=1
            self%ify=1
            self%ily=1
            self%my=1
            self%ny=1
        endif
        
        self%nz=self%ilz-self%ifz+1
        self%nx=self%ilx-self%ifx+1
        self%ny=self%ily-self%ify+1
        self%n=self%nz*self%nx*self%ny
        if(mpiworld%is_master) then
            write(*,*)'Computebox Size = [ifz,ilz] x [ifx,ilx] x [ify,ily] = ',self%n
            write(*,*)'  [ifz,ilz],nz:',self%ifz,self%ilz,self%nz
            write(*,*)'  [ifx,ilx],nx:',self%ifx,self%ilx,self%nx
            write(*,*)'  [ify,ilx],ny:',self%ify,self%ily,self%ny
            write(*,*)'Inner area of Computebox:'
            write(*,*)'  ioz,mz:',self%ioz,self%mz
            write(*,*)'  iox,mx:',self%iox,self%mx
            write(*,*)'  ioy,my:',self%ioy,self%my
        endif

        !move this piece of code to m_propagator.f90
        !as we need to write dsyn during modelings in WPI
        ! !shift source-receiver positions by computebox origin
        ! !then positions are 0-based inside computebox
        ! !source side
        ! associate(s=>shot%src)
        !     s%iz=s%iz-self%ioz+1
        !     s%ix=s%ix-self%iox+1
        !     s%iy=s%iy-self%iox+1

        !     s%ifz=s%ifz-self%ioz+1; s%ilz=s%ilz-self%ioz+1
        !     s%ifx=s%ifx-self%iox+1; s%ilx=s%ilx-self%iox+1
        !     s%ify=s%ify-self%ioy+1; s%ily=s%ily-self%ioy+1
            
        !     s%z=s%z-(self%ioz-1)*m%dz
        !     s%x=s%x-(self%iox-1)*m%dx
        !     s%y=s%y-(self%ioy-1)*m%dy
        ! end associate

        ! !receiver side
        ! associate(r=>shot%rcv)
        !     do i=1,shot%nrcv
        !         r(i)%iz=r(i)%iz-self%ioz+1
        !         r(i)%ix=r(i)%ix-self%iox+1
        !         r(i)%iy=r(i)%iy-self%ioy+1

        !         r(i)%ifz=r(i)%ifz-self%ioz+1; r(i)%ilz=r(i)%ilz-self%ioz+1
        !         r(i)%ifx=r(i)%ifx-self%iox+1; r(i)%ilx=r(i)%ilx-self%iox+1
        !         r(i)%ify=r(i)%ify-self%ioy+1; r(i)%ily=r(i)%ily-self%ioy+1

        !         r(i)%z=r(i)%z-(self%ioz-1)*m%dz
        !         r(i)%x=r(i)%x-(self%iox-1)*m%dx
        !         r(i)%y=r(i)%y-(self%ioy-1)*m%dy
        !     enddo
        ! end associate

        !models in computebox
        call m2cb(m%vp ,self%vp )
        call m2cb(m%vs ,self%vs )
        call m2cb(m%rho,self%rho)
        call m2cb(m%eps,self%eps)
        call m2cb(m%del,self%del)
        call m2cb(m%eta,self%eta)
        call m2cb(m%qp ,self%qp )
        call m2cb(m%qs ,self%qs )

        self%velmin=minval(self%vp)
        if(allocated(self%vs)) then
            self%velmin=min( self%velmin, minval(self%vs,self%vs>0.) )
        endif

        self%velmax=maxval(self%vp)
        if(allocated(self%eps)) then
            self%velmax=max( self%velmax, maxval(self%vp*sqrt(1.+2*self%eps)) )  !I don't think negative epsilon value can play a role here..
        endif

        call hud('Computebox value ranges:')
        if(mpiworld%is_master) then
                                    write(*,*)'vp' ,minval(self%vp),maxval(self%vp)
            if(allocated(self%vs )) write(*,*)'vs' ,minval(self%vs),maxval(self%vs)
                                    write(*,*)'rho',minval(self%rho),maxval(self%rho)
                                    write(*,*)'ip' ,minval(self%vp*self%rho),maxval(self%vp*self%rho)
            if(allocated(self%eps)) write(*,*)'eps',minval(self%eps),maxval(self%eps)
            if(allocated(self%del)) write(*,*)'del',minval(self%del),maxval(self%del)
        end if

    end subroutine
    
    subroutine m2cb(big,small)
        real,dimension(:,:,:),allocatable :: big, small
       
        if(.not.allocated(big)) return

        call alloc(small, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        !values inside..
        !resample z axis? Maybe interesting for elastic modeling..
        small(1:cb%mz,1:cb%mx,1:cb%my) = &
            big(cb%ioz:cb%ioz+cb%mz-1, &
                cb%iox:cb%iox+cb%mx-1, &
                cb%ioy:cb%ioy+cb%my-1)

        !values in boundary layers
        !!top
        do iz=cb%ifz ,0      ; small(iz,:,:)=small(1    ,:,:) ; enddo
        !!bottom
        do iz=cb%mz+1,cb%ilz ; small(iz,:,:)=small(cb%mz,:,:) ; enddo
        !!left
        do ix=cb%ifx ,0      ; small(:,ix,:)=small(:,1    ,:) ; enddo
        !!right
        do ix=cb%mx+1,cb%ilx ; small(:,ix,:)=small(:,cb%mx,:) ; enddo
        !!front
        do iy=cb%ify ,0      ; small(:,:,iy)=small(:,:,1    ) ; enddo
        !!rear
        do iy=cb%my+1,cb%ily ; small(:,:,iy)=small(:,:,cb%my) ; enddo

    end subroutine

    subroutine project_back(self)
        class(t_computebox) :: self
        
        call cb2m(m%gradient,cb%grad)
        call cb2m(m%image   ,cb%imag)
        call cb2m(m%energy  ,cb%engy)
        call cb2m(m%correlate,cb%corr)

        call final(self)

    end subroutine

    subroutine cb2m(big,small)
        real,dimension(:,:,:,:),allocatable :: big, small

        if(.not. allocated(small)) return

        call alloc(big,m%nz,m%nx,m%ny,size(small,4),oif_protect=.true.)

        big(cb%ioz:cb%ioz+cb%mz-1,&
            cb%iox:cb%iox+cb%mx-1,&
            cb%ioy:cb%ioy+cb%my-1,:) = &
        big(cb%ioz:cb%ioz+cb%mz-1,&
            cb%iox:cb%iox+cb%mx-1,&
            cb%ioy:cb%ioy+cb%my-1,:) + small(:,:,:,:)

    end subroutine
    
    subroutine final(self)
        type(t_computebox) :: self

        call dealloc(self%vp, self%vs, self%rho)
        call dealloc(self%eps,self%del,self%eta)
        call dealloc(self%qp,self%qs)
        call dealloc(self%grad,self%imag,self%engy)

    end subroutine

    logical function is_registered(self,chp,str)
        class(t_computebox) :: self
        type(t_checkpoint) :: chp
        character(*) :: str

        type(t_string),dimension(:),allocatable :: list

        list=split(str)

        do i=1,size(list)
            is_registered=chp%check('computebox%'//list(i)%s)
            if(.not.is_registered) return
        enddo

        do i=1,size(list)
            select case (list(i)%s)
            case ('corr')
                call chp%open('computebox%corr')
                call chp%read(self%imag,self%grad,self%engy)
                call chp%close
                call hud('Read computebox%corr from '//chp%name//', size='//num2str(total_size(self%imag,self%grad,self%engy)))
            end select

        enddo

    end function

    subroutine register(self,chp,str)
        class(t_computebox) :: self
        type(t_checkpoint) :: chp
        character(*) :: str

        type(t_string),dimension(:),allocatable :: list

        list=split(str)

        do i=1,size(list)
            select case (list(i)%s)
            case ('corr')
                call chp%open('computebox%corr')
                call chp%write(self%imag,self%grad,self%engy)
                call chp%close
            end select

        enddo

    end subroutine

end