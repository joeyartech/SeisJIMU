module m_computebox
use m_System
use m_math
use m_model
use m_shot

    private

!geometry of model & computebox
!                                      ifx<---- nx ---->ilx
!                                    ifz+~~~~~~~~~~~~~~~~+
!                                     ^ l   1      mx    l
!                                     | l   +-------+    l
!       iox     iox+mx-1              | l   |   A   |    l
!  ++====+-------+===========++       | l  1+-------+    l
!  ||    | *     |           ||         l   |       |    l
!  || B  |   C   |     B     ||      nz l D |   C   |    l
!  ||    |       |           ||         l   |       |    l
!  ||    |       |           ||       | l   |       |    l
!  ++====+-------+===========++       | l mz+-------+    l
!                                     v l absorbinglayer l
!       Model: B+C+B                 ilz+~~~~~~~~~~~~~~~~+
!       * source point
!                                        Computebox: A+C+D

    real,dimension(4) :: add_aperture
    
    type,public :: t_computebox

        integer :: nairlayer !thickness of A
        integer :: nabslayer !thickness of D

        !C's index in Model
        integer :: ioz,iox,ioy
        
        !C's index in Computebox
        integer :: mz,mx,my !=m%nz,nx,ny
        
        !A+C+D's index in Computebox
        integer :: ifz,ifx,ify
        integer :: ilz,ilx,ily
        integer :: nz,nx,ny,n !>m%nz,nx,ny,n

        ! real :: cell_volume, cell_diagonal, cell_inv_diagonal

        real :: velmin=huge(1.), velmax=0.

        ! real,dimension(:,:,:),allocatable :: vp,vs,rho
        ! real,dimension(:,:,:),allocatable :: eps,del,eta
        ! real,dimension(:,:,:),allocatable :: qp,qs

        real,dimension(:,:,:),allocatable :: mur,epsr,sgma
        real,dimension(:,:,:),allocatable :: celerity
        
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

        !thickness of A
        self%nairlayer=setup%get_int('AIRLAYER_THICKNESS','NAIRL',o_default='10')

        !thickness of D = 
        !thickness of user given + thickness required by propagator
        self%nabslayer=setup%get_int('REF_BNDLAYER_THICKNESS','NCPML',o_default='20')+add_abslayer 

    end subroutine
    
    subroutine project(self,ois_background)
        class(t_computebox) :: self
        logical,optional :: ois_background

        real,dimension(:,:,:),allocatable :: tmp_imp

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
        self%ifz = 1       - self%nabslayer - self%nairlayer
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

        !models in computebox
        if(either(ois_background,.false.,present(ois_background))) then
            call m2cb(m%epsr0,self%epsr,1.)
            call m2cb(m%mur0, self%mur, 1.)
        else
            call m2cb(m%epsr,self%epsr,1.)
            call m2cb(m%mur, self%mur, 1.)
        endif
        call m2cb(m%sgma,self%sgma,0.)

        self%celerity = r_c0/sqrt(self%epsr*self%mur)
        
        self%velmin=minval(self%celerity)
        self%velmax=maxval(self%celerity)

        call hud('Computebox value ranges:')
        if(mpiworld%is_master) then
            write(*,*)'epsr' ,minval(self%epsr),maxval(self%epsr)
            write(*,*)'mur'  ,minval(self%mur ),maxval(self%mur )
            write(*,*)'celerity', minval(self%celerity),maxval(self%celerity)
            tmp_imp = sqrt((r_mu0/r_eps0)*(self%mur/self%epsr))
            write(*,*)'impedance',minval(tmp_imp),maxval(tmp_imp)
            write(*,*)'sgma',minval(self%sgma),maxval(self%sgma)
            
        end if

    end subroutine
    
    subroutine m2cb(big,small,air_value)
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
        !!air
        do iz=0,-cb%nairlayer+1,-1         ; small(iz,:,:)=air_value    ; enddo
        !!top
        do iz=cb%ifz,cb%ifz+cb%nabslayer-1 ; small(iz,:,:)=small(cb%ifz+cb%nabslayer,:,:) ; enddo
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
        
        ! call cb2m(correlation_gradient,cb%grad)
        ! !call cb2m(m%image   ,cb%imag)
        ! call cb2m(m%energy  ,cb%engy)
        ! call cb2m(m%correlate,cb%corr)

        call final(self)

    end subroutine

    subroutine cb2m(big,small)
        real,dimension(:,:,:,:),allocatable :: big, small

        if(.not. allocated(small)) return
        if(.not. allocated(big  )) return

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

        call dealloc(self%epsr,self%mur, self%sgma)
        call dealloc(self%celerity)
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
