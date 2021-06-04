module m_preconditioner
use m_sysio
use m_arrayop
use m_model, only: m
use m_parameterization, only: npar

    real,dimension(:,:,:),allocatable :: precond

    contains
    
    subroutine init_preconditioner
        file=setup%get_file('FILE_PRECONDITIONER')
        if(file/='') then
            call alloc(precond,m%nz,m%nx,m%ny)
            open(12,file=file,access='direct',recl=4*m%n,action='read',status='old')
            read(12,rec=1) precond
            close(12)
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
            pg(:,:,:,i)=g(:,:,:,i)*precond
        enddo
        pg = pg * old_norm / norm2(pg)

    end subroutine

end
