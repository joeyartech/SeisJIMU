module m_preconditioner
use m_sysio
use m_arrayop
use m_model, only: m
use m_parameterization, only: npar

    real zpower
    real,dimension(:),allocatable :: precond

    real,dimension(:,:,:),allocatable :: zpow_user
    
    contains
    
    subroutine init_preconditioner
        logical :: alive=.false.

        zpower=get_setup_real('ZPOWER',default=1.)
        
        call alloc(precond,m%n*2)
        precond=[(((iz-1)*m%dz)**zpower,iz=1,m%nz)]

        inquire(file='zpow_user', exist=alive)
        if(alive) then
            call alloc(zpow_user,m%nz,m%nx,m%ny)
            open(12,file='zpow_user',access='direct',recl=4*m%n,action='read',status='old')
            read(12,rec=1) zpow_user
            close(12)
            call hud('zpow_user is read. ZPOWER input is not used.')
        endif
        
    end subroutine
    
    subroutine preconditioner_apply(g,pg)
        real,dimension(m%nz,m%nx,m%ny,npar) :: g,pg
        
        real :: old_norm

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
                pg(:,ix,iy,i)=g(:,ix,iy,i)*precond
            enddo
            enddo
            if(allocated(zpow_user)) then
               pg(:,:,:,i)=g(:,:,:,i)*zpow_user
            endif
        enddo
        pg = pg * old_norm / norm2(pg)

    end subroutine

end
