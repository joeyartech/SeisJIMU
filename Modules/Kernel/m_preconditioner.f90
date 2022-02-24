module m_preconditioner
use m_System
use m_Modeling
use m_parametrizer

    private

    type,public :: t_preconditioner

        contains

        procedure :: update
        procedure :: apply
        procedure :: apply_ext

    end type

    type(t_preconditioner),public :: preco

    real,dimension(:,:,:),allocatable :: preco_in_x
    real,dimension(:,:,:),allocatable :: preco_in_m
    
    contains
    
    subroutine update(self)
        class(t_preconditioner) :: self

        type(t_string),dimension(:),allocatable :: list, sublist
        character(:),allocatable :: file

        call alloc(preco_in_m,m%nz,m%nx,m%ny,o_init=1.)

        list=setup%get_strs('PRECONDITIONING','PRECO',o_default='z^1')

        do i=1,size(list)
            if(index(list(i)%s,'z^')>0) then
                sublist=split(list(i)%s,o_sep='^')
                call hud('Will precondition the gradient by z^'//sublist(2)%s)
                call by_depth(o_power=str2real(sublist(2)%s))
            endif

            if(index(list(i)%s,'z*')>0) then
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will precondition the gradient by z*'//sublist(2)%s)
                call by_depth(o_factor=str2real(sublist(2)%s))
            endif
            
            if(index(list(i)%s,'energy')>0) then
                sublist=split(list(i)%s,o_sep=':')
                iengy=either(str2int(sublist(2)%s),1,len(sublist(2)%s)>0)
                call hud('Will precondition the gradient by '//num2str(iengy)//"'th energy term")
                if((iengy) > size(m%energy,4)) then
                    call warn('The chosen energy term does NOT exist! Use instead the 1st term provided from propagator.')
                    iengy=1
                endif
                call by_energy(iengy)
            endif

            if(index(list(i)%s,'custom')>0) then
                sublist=split(list(i)%s,o_sep=':')
                if(size(sublist)==1) then !filename is not attached, ask for it
                    file=setup%get_file('FILE_PRECONDITION_CUSTOM',o_mandatory=1)
                else !filename is attached
                    file=sublist(2)%s
                endif

                call hud('Will precondition the gradient in a custom way defined in '//file)
                call by_custom(file)
            endif

        enddo

        call alloc(preco_in_x,param%n1,param%n2,param%n3)
        call param%transform_preconditioner(preco_in_m,preco_in_x)

    end subroutine

    subroutine by_depth(o_power,o_factor)
        real,optional :: o_power, o_factor

        if(present(o_power)) then
            preco_in_m(:,1,1)=[(((iz-1)*m%dz)**o_power,iz=1,m%nz)]
            do iy=1,m%ny; do ix=1,m%nx
                preco_in_m(:,ix,iy)=preco_in_m(:,1,1)
            enddo; enddo
        endif

        if(present(o_factor)) then
            preco_in_m(:,1,1)=[(((iz-1)*m%dz)*o_factor,iz=1,m%nz)]
            do iy=1,m%ny; do ix=1,m%nx
                preco_in_m(:,ix,iy)=preco_in_m(:,1,1)
            enddo; enddo
        endif
        
    end subroutine

    subroutine by_energy(iengy)

        preco_in_m=1./m%energy(:,:,:,iengy)

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

    subroutine by_custom(file)
        character(*) :: file

        integer :: file_size
        real,dimension(:,:,:),allocatable :: tmp

        inquire(file=file,size=file_size)
        if(file_size < 4*m%n) call warn('FILE_PRECONDITION_CUSTOM has size '//num2str(file_size/4)//' < nz*nx*ny. Some part of the gradient will not be preconditioned!')
        
        call alloc(tmp,m%nz,m%nx,m%ny,o_init=1.)

        call sysio_read(file,tmp,size(tmp))

        preco_in_m=preco_in_m*tmp

        deallocate(tmp)

    end subroutine
    
    subroutine apply(self,g,pg)
        class(t_preconditioner) :: self
        real,dimension(:,:,:,:),allocatable :: g,pg
        
        real :: old_norm

        call alloc(pg,param%n1,param%n2,param%n3,param%npars,oif_protect=.true.)

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

        !precond whole grediant to keep the gradient norm
        old_norm = norm2(g)
        do i=1,param%npars
            pg(:,:,:,i)=g(:,:,:,i)*preco_in_x(:,:,:)
        enddo
        pg = sqrt(old_norm/norm2(pg)) * pg
        
        !save some RAM
        call dealloc(preco_in_x,preco_in_m)

    end subroutine

    subroutine apply_ext(self,g,pg)
        class(t_preconditioner) :: self
        real,dimension(param%n1,param%n2,param%n3,param%npars) :: g,pg
        real :: old_norm

        old_norm = norm2(g)
        do i=1,param%npars
            pg(:,:,:,i)=g(:,:,:,i)*preco_in_x(:,:,:)
        enddo
        pg = sqrt(old_norm/norm2(pg)) * pg
    end subroutine

end
