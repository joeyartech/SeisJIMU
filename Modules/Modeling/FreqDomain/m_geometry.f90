module m_geometry
use m_shotlist
use m_shot

    private t_receiver, t_source
    
    type t_receiver
        integer :: ix,iy,iz
        integer :: ifz,ilz,ifx,ilx,ify,ily
        real,dimension(:,:),allocatable :: interp_coeff
    end type
    
    type t_source
        integer :: ix,iy,iz
        integer :: ifz,ilz,ifx,ilx,ify,ily
        real,dimension(:,:),allocatable :: interp_coeff

        type(t_receiver),dimension(:),allocatable :: rcv
        integer :: nrcv
    end type
    
    type t_geometry
        type(t_source),dimension(:),allocatable :: src
        integer :: nsrc
        integer :: ntr !total number of traces
    end type

    type(t_geometry) :: geo

    contains

    subroutine build_geom_acqui
        integer,parameter :: r=4 !r chosen to be 4 in m_hicks

        call build_shotlist(acqui%nsrc)

        !allocate
        allocate(geo%src(nshots))
        do i=1,nshots
            allocate(geo%src(i)%interp_coeff(-r:r,-r:r))
            allocate(geo%src(i)%rcv(shot%nrcv))

            do j=1,shot%nrcv
                allocate(geo%src(i)%rcv(j)%interp_coeff(-r:r,-r:r))
            enddo
        enddo

        geo%ntr=0

        !sources
        geo%nsrc=nshots
        do i=1,nshots

            call init_shot(i,'setup')
            geo%src(i)%ifz=shot%src%ifz; geo%src(i)%iz=shot%src%iz; geo%src(i)%ilz=shot%src%ilz
            geo%src(i)%ifx=shot%src%ifx; geo%src(i)%ix=shot%src%ix; geo%src(i)%ilx=shot%src%ilx

            geo%src(i)%interp_coeff=shot%src%interp_coeff(:,:,1)

            !receivers
            geo%src(i)%nrcv=shot%nrcv

            geo%src(i)%rcv(:)%ifz=shot%rcv(:)%ifz; geo%src(i)%rcv(:)%iz=shot%rcv(:)%iz; geo%src(i)%rcv(:)%ilz=shot%rcv(:)%ilz
            geo%src(i)%rcv(:)%ifx=shot%rcv(:)%ifx; geo%src(i)%rcv(:)%ix=shot%rcv(:)%ix; geo%src(i)%rcv(:)%ilx=shot%rcv(:)%ilx

            geo%src(i)%interp_coeff=shot%src%interp_coeff(:,:,1)

            do j=1,shot%nrcv
                geo%src(i)%rcv(j)%interp_coeff=shot%rcv(j)%interp_coeff(:,:,1)
            enddo

            geo%ntr=geo%ntr+geo%src(i)%nrcv

        enddo

    end subroutine

    ! !reassign indices considering extended model
    ! subroutine geom_reassign
    !     npml=get_setup_int('NBOUNDARYLAYER','NPML',default=10)

    !     !sources
    !     geo%src(:)%ifz=geo%src(:)%ifz+npml;  geo%src(:)%iz=geo%src(:)%iz+npml; geo%src(:)%ilz=geo%src(:)%ilz+npml
    !     geo%src(:)%ifx=geo%src(:)%ifx+npml;  geo%src(:)%iz=geo%src(:)%ix+npml; geo%src(:)%ilx=geo%src(:)%ilx+npml

    !     !receivers
    !     do i=1,geo%nsrc
    !         geo%src(i)%rcv(:)%ifz=geo%src(i)%rcv(:)%ifz+npml; geo%src(i)%rcv(:)%iz=geo%src(i)%rcv(:)%iz+npml; geo%src(i)%rcv(:)%ilz=geo%src(i)%rcv(:)%ilz+npml
    !         geo%src(i)%rcv(:)%ifx=geo%src(i)%rcv(:)%ifx+npml; geo%src(i)%rcv(:)%ix=geo%src(i)%rcv(:)%ix+npml; geo%src(i)%rcv(:)%ilz=geo%src(i)%rcv(:)%ilx+npml
    !     enddo

    ! end subroutine

end