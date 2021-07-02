module m_computebox
use m_string
use m_arrayop
use m_setup
use m_checkpoint
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

    real,dimension(4) :: aperture

    type,public :: t_computebox
        
        !C's index in Model
        integer :: ioz,iox,ioy
        
        !C's index in Computebox
        integer :: mz,mx,my
        
        !C+D's index in Computebox
        integer :: ifz,ifx,ify
        integer :: ilz,ilx,ily
        integer :: nz,nx,ny,n
        
        ! real :: cell_volume, cell_diagonal, cell_inv_diagonal

        real velmin, velmax

        real,dimension(:,:,:),allocatable :: vp,vs,rho
        real,dimension(:,:,:),allocatable :: eps,del,eta
        real,dimension(:,:,:),allocatable :: qp,qs
        
        real,dimension(:,:,:,:),allocatable :: kernel
        
        contains
        procedure :: init => init
        procedure :: project => project
        procedure :: project_back => project_back
        final :: fin

        procedure :: is_registered => is_registered
        procedure :: register => register
    end type

    type(t_computebox),public :: cb

    contains

    subroutine init(self)
        class(t_computebox) :: self
        
        !shot aperture, default to whole model
        aperture=setup%get_reals('APERTURE',o_default='-99999 99999 -99999 99999')

        ! self%cell_volume = m%dz*m%dx*m%dy

        ! if(m%is_cubic) then
        !     self%cell_diagonal=sqrt(m%dz**2+m%dx**2+m%dy**2)
        !     self%cell_inv_diagonal=sqrt(m%dz**(-2) + m%dx**(-2) + m%dy**(-2))
        ! else
        !     self%cell_diagonal=sqrt(m%dz**2+m%dx**2)
        !     self%cell_inv_diagonal=sqrt(m%dz**(-2) + m%dx**(-2))
        ! endif

    end subroutine
    
    subroutine project(self,is_fdsg,abslayer_width)
        class(t_computebox) :: self
        logical :: is_fdsg
        integer :: abslayer_width

        !C's origin index in model
        self%ioz=1 !always from top of model
        x=minval(shot%rcv(:)%x)
        x=min(x,shot%src%x,shot%src%x+aperture(1))
        y=minval(shot%rcv(:)%y)
        y=min(y,shot%src%y,shot%src%y+aperture(3))
        self%iox=max(1,   nint(x/m%dx)+1) !can't exceed size of model
        self%ioy=max(1,   nint(y/m%dy)+1) !can't exceed size of model
        
        !C's size
        self%mz=m%nz !always down to bottom of model
        x=maxval(shot%rcv(:)%x)
        x=max(x,shot%src%x,shot%src%x+aperture(2))
        y=maxval(shot%rcv(:)%y)
        y=max(y,shot%src%y,shot%src%y+aperture(4))
        self%mx=min(m%nx,nint(x/m%dx)+1) +1 -self%iox
        self%my=min(m%ny,nint(y/m%dy)+1) +1 -self%ioy
        
        !C+D's index
        self%ifx = 1       - abslayer_width
        self%ilx = self%mx + abslayer_width
        self%ify = 1       - abslayer_width
        self%ily = self%my + abslayer_width
        self%ifz = 1       - abslayer_width
        self%ilz = self%mz + abslayer_width
        
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
        
        !shift source-receiver positions by computebox origin
        !then positions are 0-based inside computebox
        !source side
        shot%src%iz=shot%src%iz-self%ioz+1
        shot%src%ix=shot%src%ix-self%iox+1
        shot%src%iy=shot%src%iy-self%ioy+1

        shot%src%ifz=shot%src%ifz-self%ioz+1; shot%src%ilz=shot%src%ilz-self%ioz+1
        shot%src%ifx=shot%src%ifx-self%iox+1; shot%src%ilx=shot%src%ilx-self%iox+1
        shot%src%ify=shot%src%ify-self%ioy+1; shot%src%ily=shot%src%ily-self%ioy+1
        
        shot%src%z=shot%src%z-(self%ioz-1)*m%dz
        shot%src%x=shot%src%x-(self%iox-1)*m%dx
        shot%src%y=shot%src%y-(self%ioy-1)*m%dy
                    
        !receiver side
        do i=1,shot%nrcv
            shot%rcv(i)%iz=shot%rcv(i)%iz-self%ioz+1
            shot%rcv(i)%ix=shot%rcv(i)%ix-self%iox+1
            shot%rcv(i)%iy=shot%rcv(i)%iy-self%ioy+1
            shot%rcv(i)%ifz=shot%rcv(i)%ifz-self%ioz+1; shot%rcv(i)%ilz=shot%rcv(ir)%ilz-self%ioz+1
            shot%rcv(i)%ifx=shot%rcv(i)%ifx-self%iox+1; shot%rcv(i)%ilx=shot%rcv(ir)%ilx-self%iox+1
            shot%rcv(i)%ify=shot%rcv(i)%ify-self%ioy+1; shot%rcv(i)%ily=shot%rcv(ir)%ily-self%ioy+1
            shot%rcv(i)%z=shot%rcv(i)%z-(self%ioz-1)*m%dz
            shot%rcv(i)%x=shot%rcv(i)%x-(self%iox-1)*m%dx
            shot%rcv(i)%y=shot%rcv(i)%y-(self%ioy-1)*m%dy
        enddo
        

        !models in computebox
        if(allocated(m%vp )) call m2cb(m%vp ,self%vp )
        if(allocated(m%vs )) call m2cb(m%vs ,self%vs )
        if(allocated(m%rho)) call m2cb(m%rho,self%rho)
        if(allocated(m%eps)) call m2cb(m%eps,self%eps)
        if(allocated(m%del)) call m2cb(m%del,self%del)
        if(allocated(m%eta)) call m2cb(m%eta,self%eta)
        if(allocated(m%qp )) call m2cb(m%qp ,self%qp )
        if(allocated(m%qs )) call m2cb(m%qs ,self%qs )
                
        self%velmin=min(minval(self%vp), minval(self%vs,self%vs>0.))
        self%velmax=maxval(self%vp*sqrt(1.+2*self%eps)) !I don't think negative epsilon value can play a role here..
        
        call hud('Computebox value ranges:')
        if(mpiworld%iproc==0) then
                                    write(*,*)'vp' ,minval(self%vp),maxval(self%vp)
            if(allocated(self%vs )) write(*,*)'vs' ,minval(self%vs),maxval(self%vs)
                                    write(*,*)'rho',minval(self%rho),maxval(self%rho)
                                    write(*,*)'ip' ,minval(self%vp*self%rho),maxval(self%vp*self%rho)
            if(allocated(self%eps)) write(*,*)'eps',minval(self%eps),maxval(self%eps)
            if(allocated(self%del)) write(*,*)'del',minval(self%del),maxval(self%del)
        end if
        
    end subroutine
    
    subroutine m2cb(big,small)
        real,dimension(m%nz,m%nx,m%ny),intent(in) :: big
        real,dimension(:,:,:),allocatable,intent(out) :: small
       
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

    subroutine project_back(self,big,small,n)
        class(t_computebox) :: self
        real,dimension(m%nz,m%nx,m%ny,n) :: big
        real,dimension(cb%mz,cb%mx,cb%my, n) :: small

            big(self%ioz:self%ioz+self%mz-1,&
                self%iox:self%iox+self%mx-1,&
                self%ioy:self%ioy+self%my-1,:) = &
            big(self%ioz:self%ioz+self%mz-1,&
                self%iox:self%iox+self%mx-1,&
                self%ioy:self%ioy+self%my-1,:) + small(:,:,:,:)

    end subroutine
    
    subroutine fin(s)
        type(t_computebox) :: s

        call dealloc(s%vp,s%vs,s%rho)
        call dealloc(s%eps,s%del,s%eta)
        call dealloc(s%qp,s%qs)
        call dealloc(s%kernel)

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
            case ('kernel')
                call chp%open('computebox%kernel')
                call chp%read(self%kernel,size(self%kernel))
                call chp%close
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
            case ('kernel')
                call chp%open('computebox%kernel')
                call chp%write(self%kernel,size(self%kernel))
                call chp%close
            end select

        enddo

    end subroutine

end