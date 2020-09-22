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

    type(t_geometry) :: geom

    contains

    subroutine build_geom_acqui
        integer,parameter :: r=4 !r chosen to be 4 in m_hicks

        call build_shotlist(acqui%nsrc)
        
        !sources
        geom%nsrc=nshots
        allocate(geom%src(nshots))
        
        geom%ntr=0

        do i=1,nshots

            call init_shot(i,'setup')
            
            !sz,sx
            geom%src(i)%ifz=shot%src%ifz; geom%src(i)%iz=shot%src%iz; geom%src(i)%ilz=shot%src%ilz
            geom%src(i)%ifx=shot%src%ifx; geom%src(i)%ix=shot%src%ix; geom%src(i)%ilx=shot%src%ilx

            !hicks coeff
            allocate(geom%src(i)%interp_coeff(-r:r,-r:r))
            geom%src(i)%interp_coeff=shot%src%interp_coeff(:,:,1)

            !receivers
            geom%src(i)%nrcv=shot%nrcv
            allocate(geom%src(i)%rcv(shot%nrcv))
            
            !total number of traces
            geom%ntr=geom%ntr+shot%nrcv
            
            !rz,rx
            geom%src(i)%rcv(:)%ifz=shot%rcv(:)%ifz; geom%src(i)%rcv(:)%iz=shot%rcv(:)%iz; geom%src(i)%rcv(:)%ilz=shot%rcv(:)%ilz
            geom%src(i)%rcv(:)%ifx=shot%rcv(:)%ifx; geom%src(i)%rcv(:)%ix=shot%rcv(:)%ix; geom%src(i)%rcv(:)%ilx=shot%rcv(:)%ilx

            !hicks coeff
            do j=1,shot%nrcv
                allocate(geom%src(i)%rcv(j)%interp_coeff(-r:r,-r:r))
                geom%src(i)%rcv(j)%interp_coeff=shot%rcv(j)%interp_coeff(:,:,1)
            enddo

        enddo

    end subroutine

    ! !reassign indices considering extended model
    ! subroutine geom_reassign
    !     npml=get_setup_int('NBOUNDARYLAYER','NPML',default=10)

    !     !sources
    !     geom%src(:)%ifz=geom%src(:)%ifz+npml;  geom%src(:)%iz=geom%src(:)%iz+npml; geom%src(:)%ilz=geom%src(:)%ilz+npml
    !     geom%src(:)%ifx=geom%src(:)%ifx+npml;  geom%src(:)%iz=geom%src(:)%ix+npml; geom%src(:)%ilx=geom%src(:)%ilx+npml

    !     !receivers
    !     do i=1,geom%nsrc
    !         geom%src(i)%rcv(:)%ifz=geom%src(i)%rcv(:)%ifz+npml; geom%src(i)%rcv(:)%iz=geom%src(i)%rcv(:)%iz+npml; geom%src(i)%rcv(:)%ilz=geom%src(i)%rcv(:)%ilz+npml
    !         geom%src(i)%rcv(:)%ifx=geom%src(i)%rcv(:)%ifx+npml; geom%src(i)%rcv(:)%ix=geom%src(i)%rcv(:)%ix+npml; geom%src(i)%rcv(:)%ilz=geom%src(i)%rcv(:)%ilx+npml
    !     enddo

    ! end subroutine

end