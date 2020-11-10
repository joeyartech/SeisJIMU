module m_preconditioner
use m_sysio
use m_arrayop
use m_model, only: m
use m_parameterization, only: npar

    character(:),allocatable :: tl_precond_strategy

    real zpower
    real,dimension(:),allocatable :: precond
    
    contains
    
    subroutine init_preconditioner
        tl_precond_strategy=get_setup_char('4D_PRECOND',default='simul')

        zpower=get_setup_real('ZPOWER',default=1.)
        
        call alloc(precond,m%n*2)
        precond=[(((iz-1)*m%dz)**zpower,iz=1,m%nz)]
        
    end subroutine
    
    subroutine preconditioner_apply(g,pg)
        real,dimension(m%nz,m%nx,m%ny,npar,2) :: g,pg
        
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

        if(tl_precond_strategy(1:5)=='simul') then !rescale base & monitor gradients simultaneously
            old_norm = norm2(g)
            do j=1,2
            do i=1,npar
            do iy=1,m%ny
            do ix=1,m%nx
                pg(:,ix,iy,i,j)=g(:,ix,iy,i,j)*precond
            enddo
            enddo
            enddo
            enddo
            pg = pg * old_norm / norm2(pg)

        else  !rescale base & monitor gradients separately
            !precond baseline grediant
            old_norm = norm2(g(:,:,:,:,1))
            do i=1,npar
            do iy=1,m%ny
            do ix=1,m%nx
                pg(:,ix,iy,i,1)=g(:,ix,iy,i,1)*precond
            enddo
            enddo
            enddo
            pg(:,:,:,:,1) = pg(:,:,:,:,1) * old_norm / norm2(pg(:,:,:,:,1))

            !precond monitor grediant
            old_norm = norm2(g(:,:,:,:,2))
            do i=1,npar
            do iy=1,m%ny
            do ix=1,m%nx
                pg(:,ix,iy,i,2)=g(:,ix,iy,i,2)*precond
            enddo
            enddo
            enddo
            pg(:,:,:,:,2) = pg(:,:,:,:,2) * old_norm / norm2(pg(:,:,:,:,2))

        endif

    end subroutine

end
