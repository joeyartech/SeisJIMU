module m_preconditioner
use m_string
use m_arrayop
use m_setup
use m_mpienv
use m_model
use m_propagator
use m_parameterization

    private

    logical :: if_will_update=.false.
    real :: depth

    type,public :: t_preconditioner

        real,dimension(:,:,:),allocatable :: preco

        contains

        procedure :: init => init
        procedure :: by_depth => by_depth
        procedure :: by_sfield => by_sfield
        procedure :: by_custom => by_custom
        procedure :: update => update
        procedure :: apply => apply

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
            if(list(i)%s(1:2)=='z^') then
                sublist=split(list(i)%s,o_sep='^')
                call hud('Will precondition the gradient by z^'//sublist(2)%s)
                call self%by_depth(o_power=str2real(sublist(2)%s))
            endif

            if(list(i)%s(1:2)=='z*') then
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will precondition the gradient by z*'//sublist(2)%s)
                call self%by_depth(o_factor=str2real(sublist(2)%s))
            endif

            if(list(i)%s(1:2)=='sfield') then
                sublist=split(list(i)%s,o_sep='/')
                call hud('Will precondition the gradient by autocorrelation of sfield')
                if_will_update=.true.
                depth=str2real(sublist(2)%s)
                call self%by_sfield
            endif

            if(list(i)%s(1:6)=='custom') then
                sublist=split(list(i)%s,o_sep=':')
                if(size(sublist)==1) then !filename is not attached, ask for it
                    file=setup%get_file('FILE_PRECONDITION_CUSTOM',o_mandatory=.true.)
                else !filename is attached
                    file=sublist(2)%s
                endif

                call hud('Will precondition the gradient in a custom way defined in '//file)
                call self%by_custom(file)

            endif

        enddo

    end subroutine

    subroutine update(self)
        class(t_preconditioner) :: self
        call self%by_sfield
    end subroutine

    subroutine by_depth(self,o_power,o_factor)
        class(t_preconditioner) :: self
        real,optional :: o_power, o_factor

        if(present(o_power)) then
            self%preco(:,1,1)=[(((iz-1)*m%dz)**o_power,iz=1,m%nz)]
            do iy=1,m%ny; do ix=1,nx
                self%preco(:,ix,iy)=self%preco(:,1,1)
            enddo; enddo
        endif

        if(present(o_factor)) then
            self%preco(:,1,1)=[(((iz-1)*m%dz)*o_factor,iz=1,m%nz)]
            do iy=1,m%ny; do ix=1,nx
                self%preco(:,ix,iy)=self%preco(:,1,1)
            enddo; enddo
        endif
        
    end subroutine

    subroutine by_custom(self,file)
        class(t_preconditioner) :: self
        character(*) :: file

        integer :: file_size
        real,dimension(:,:,:),allocatable :: tmp

        inquire(file=file,size=file_size)
        if(file_size < 4*m%n) call warn('FILE_PRECONDITION_CUSTOM has size '//num2str(file_size/4)//' < nz*nx*ny. Some where of the gradient will not be preconditioned.')
        
        call alloc(tmp,m%nz,m%nx,m%ny,o_init=1.)

        open(10,file=file,access='stream',action='read')
        read(10) tmp
        close(10)

        self%preco=self%preco*tmp

        deallocate(tmp)

    end subroutine

    subroutine by_sfield(self)
        class(t_preconditioner) :: self

        ! !convert gradient wrt model to gradient wrt parameters
        ! call param%transform_gradient('m->x',fobj%gradient)
        
        ! !start_depth

        ! !H11 = H_[kappainv,kappainv] = (wadj dA/dVp w)^adj (wadj dA/dVp w) = autocorr of pressure field
        ! !tmp1=pbdir%p*(pbdir%p_bwd-pbdir%p)*pbdir%kappa
        ! !inv%h11+=tmp1*tmp1 /pbdir%dt
        ! !inv%h11=inv%h11/pbdir%rho**4/pbdir%vp**8 !kappa is the model parameter instead of kappainv
        
        ! idepth=nint(depth/m%dz)

        ! do i=1,param%npars
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
        real,dimension(m%nz,m%nx,m%ny,param%npars) :: g,pg
        
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
        do i=1,param%npars
            pg(:,:,:,i)=g(:,:,:,i)*self%preco
        enddo
        pg = pg * old_norm / norm2(pg)

    end subroutine

end
