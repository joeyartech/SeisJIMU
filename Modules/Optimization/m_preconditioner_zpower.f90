module m_preconditioner
use m_sysio
use m_arrayop
use m_model, only: m
use m_parameterization, only: npar

    real zpower
    real,dimension(:),allocatable :: precond
    
    contains
    
    subroutine init_preconditioner
        zpower=get_setup_real('ZPOWER',default=1.)
        
        call alloc(precond,m%n*2)
        precond=[(((iz-1)*m%dz)**zpower,iz=1,m%nz)]
        
    end subroutine
    
    subroutine preconditioner_apply(g,pg)
        real,dimension(m%nz,m%nx,m%ny,npar) :: g,pg
        
        real(kind=8) :: norm1,norm2
        
        do i=1,npar
            norm1=sqrt(sum(dprod(g(:,:,:,i),g(:,:,:,i))))
            do iy=1,m%ny
            do ix=1,m%nx
                pg(:,ix,iy,i)=g(:,ix,iy,i)*precond
            enddo
            enddo
            norm2=sqrt(sum(dprod(pg(:,:,:,i),pg(:,:,:,i))))
            pg(:,:,:,i)=pg(:,:,:,i)*norm1/norm2
        enddo
    end subroutine

end