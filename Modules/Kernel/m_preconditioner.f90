module m_preconditioner
use m_System
use m_Modeling
use m_parametrizer

    private

    logical :: if_will_update=.false.
    real :: depth

    type,public :: t_preconditioner

        real,dimension(:,:,:),allocatable :: preco

        contains

        procedure :: init
        procedure :: update
        procedure :: apply

    end type

    type(t_preconditioner),public :: preco
    
    contains
    
    subroutine init(self)
        class(t_preconditioner) :: self

        type(t_string),dimension(:),allocatable :: list, sublist
        character(:),allocatable :: file

        call alloc(self%preco,m%nz,m%nx,m%ny,o_init=1.)

        list=setup%get_strs('PRECONDITIONING','PRECO',o_default='z^1')

        do i=1,size(list)
            select case (list(i)%s)
            case ('z^')
                sublist=split(list(i)%s,o_sep='^')
                call hud('Will precondition the gradient by z^'//sublist(2)%s)
                call by_depth(self%preco,o_power=str2real(sublist(2)%s))

            case ('z*')
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will precondition the gradient by z*'//sublist(2)%s)
                call by_depth(self%preco,o_factor=str2real(sublist(2)%s))

            case ('sfield')
                sublist=split(list(i)%s,o_sep='/')
                call hud('Will precondition the gradient by autocorrelation of sfield')
                if_will_update=.true.
                depth=str2real(sublist(2)%s)
                call by_sfield(self%preco)

            case ('custom')
                sublist=split(list(i)%s,o_sep=':')
                if(size(sublist)==1) then !filename is not attached, ask for it
                    file=setup%get_file('FILE_PRECONDITION_CUSTOM',o_mandatory=1)
                else !filename is attached
                    file=sublist(2)%s
                endif

                call hud('Will precondition the gradient in a custom way defined in '//file)
                call by_custom(self%preco,file)

            end select

        enddo

    end subroutine

    subroutine update(self)
        class(t_preconditioner) :: self
        call by_sfield(self%preco)
    end subroutine

    subroutine by_depth(preco,o_power,o_factor)
        real,dimension(m%nz,m%nx,m%ny) :: preco
        real,optional :: o_power, o_factor

        if(present(o_power)) then
            preco(:,1,1)=[(((iz-1)*m%dz)**o_power,iz=1,m%nz)]
            do iy=1,m%ny; do ix=1,m%nx
                preco(:,ix,iy)=preco(:,1,1)
            enddo; enddo
        endif

        if(present(o_factor)) then
            preco(:,1,1)=[(((iz-1)*m%dz)*o_factor,iz=1,m%nz)]
            do iy=1,m%ny; do ix=1,m%nx
                preco(:,ix,iy)=preco(:,1,1)
            enddo; enddo
        endif
        
    end subroutine

    subroutine by_custom(preco,file)
        real,dimension(m%nz,m%nx,m%ny) :: preco
        character(*) :: file

        integer :: file_size
        real,dimension(:,:,:),allocatable :: tmp

        inquire(file=file,size=file_size)
        if(file_size < 4*m%n) call warn('FILE_PRECONDITION_CUSTOM has size '//num2str(file_size/4)//' < nz*nx*ny. Some where of the gradient will not be preconditioned.')
        
        call alloc(tmp,m%nz,m%nx,m%ny,o_init=1.)

        call sysio_read(file,tmp,size(tmp))

        preco=preco*tmp

        deallocate(tmp)

    end subroutine

    subroutine by_sfield(preco)
        real,dimension(m%nz,m%nx,m%ny) :: preco

        ! !convert gradient wrt model to gradient wrt parameters
        ! call param%transform_gradient('m->x',fobj%gradient)
        
        ! !start_depth

        ! !H11 = H_[kappainv,kappainv] = (wadj dA/dVp w)^adj (wadj dA/dVp w) = autocorr of pressure field
        ! !tmp1=pbdir%p*(pbdir%p_bwd-pbdir%p)*pbdir%kappa
        ! !inv%h11+=tmp1*tmp1 /pbdir%dt
        ! !inv%h11=inv%h11/pbdir%rho**4/pbdir%vp**8 !kappa is the model parameter instead of kappainv
        
        ! idepth=nint(depth/m%dz)

        ! do i=1,ppg%ngrad
        !     do iy=1,m%ny; do ix=1,m%nx
        !     self%preco(1:idepth,ix,iy)=[(i-1)*m%dz,i=1,idepth)] &
        !         /(depth*sfield%amp(1:idepth,ix,iy))

        !     self%preco(idepth+1:m%nz,ix,iy)=1./sfield%amp(idepth+1:m%nz,:,:)
        !     enddo; enddo
        ! enddo

        ! where(ieee_is_nan(self%preco) .or. .not. ieee_is_finite(self%preco))
        !     self%preco=0.
        ! endwhere

    end subroutine
    
    subroutine apply(self,g,pg)
        class(t_preconditioner) :: self
        real,dimension(m%nz,m%nx,m%ny,ppg%ngrad) :: g,pg
        
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
        do i=1,ppg%ngrad
            pg(:,:,:,i)=g(:,:,:,i)*self%preco(:,:,:)
        enddo
        pg = pg * old_norm / norm2(pg)
        
    end subroutine

end
