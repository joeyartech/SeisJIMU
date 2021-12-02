module m_preconditioner
use m_sysio
use m_arrayop
use m_model, only: m
use m_parameterization, only: npar

    real zpower
    real,dimension(:),allocatable :: precond

    real,dimension(:,:),allocatable :: topo_zpow_stop
    integer,dimension(:,:),allocatable :: itopo_zpow_stop

    contains

    subroutine init_preconditioner
        logical :: alive

        zpower=get_setup_real('ZPOWER',default=1.)

        call alloc(precond,m%n*2)
        do iz=1,m%nz
            precond(iz) = ((iz-1)*m%dz)**zpower
        enddo
!        precond=[(((iz-1)*m%dz)**zpower,iz=1,m%nz)]


        call alloc( topo_zpow_stop,m%nx,m%ny);  topo_zpow_stop=(m%nz-1)*m%dz
        call alloc(itopo_zpow_stop,m%nx,m%ny); itopo_zpow_stop=m%nz

        inquire(file='topo_zpow_stop', exist=alive)
        if(alive) then
            open(12,file='topo_zpow_stop',access='direct',recl=4*m%nx*m%ny,action='read',status='old')
            read(12,rec=1) topo_zpow_stop
            close(12)
            call hud('topo_zpow_stop is read. zpow reset to 1 below this level.')
            itopo_zpow_stop=nint(topo_zpow_stop/m%dz)+1
        endif


    end subroutine

    subroutine preconditioner_apply(g,pg)
        real,dimension(m%nz,m%nx,m%ny,npar) :: g,pg

        real :: old_norm, new_norm

! !precond per parameter        
! do i=1,npar
!    old_norm = norm2(g(:,:,:,i))
!    do iy=1,m%ny
!    do ix=1,m%nx
!        pg(:,ix,iy,i)=g(:,ix,iy,i)*precond
!    enddo
!    enddo
!    pg(:,:,:,i) = pg(:,:,:,i) * old_norm / norm2(pg(:,:,:,i))
! enddo

        !precond whole grediant
        old_norm = norm2(g)
        do i=1,npar
        do iy=1,m%ny
        do ix=1,m%nx
            do iz=1,m%nz
                pg(iz,ix,iy,i)=g(iz,ix,iy,i)*precond(iz)
            enddo
        enddo
        enddo
        enddo
        new_norm = norm2(pg)

        do i=1,npar
        do iy=1,m%ny
        do ix=1,m%nx
            do iz=itopo_zpow_stop(ix,iy)+1,m%nz
                pg(iz,ix,iy,i)=g(iz,ix,iy,i)*(iz-1)*m%dz
            enddo
        enddo
        enddo
        enddo

        pg = pg * old_norm / new_norm

    end subroutine

end
