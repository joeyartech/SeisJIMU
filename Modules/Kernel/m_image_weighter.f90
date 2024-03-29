module m_image_weighter
use m_System
use m_Modeling

    private

    type,public :: t_image_weighter

        integer nz,nx,ny
        real,dimension(:,:,:),allocatable :: weight

        contains

        procedure :: update

        procedure :: by_depth
        procedure :: by_custom

    end type

    type(t_image_weighter),public :: iwei
    
    contains

    subroutine update(self,o_nz,o_nx,o_ny)
        class(t_image_weighter) :: self
        integer,optional :: o_nz, o_nx, o_ny

        type(t_string),dimension(:),allocatable :: list,sublist
        character(:),allocatable :: file

        self%nz = either(o_nz, m%nz, present(o_nz))
        self%nx = either(o_nx, m%nx, present(o_nx))
        self%ny = either(o_ny, m%ny, present(o_ny))
        call alloc(self%weight,self%nz,self%nx,self%ny,o_init=1.)

        list=setup%get_strs('IMAGE_WEIGHTING','IWEI',o_default='one')

        do i=1,size(list)

            if (index(list(i)%s,'topo')>0) then
                where(m%is_freeze_zone) self%weight=0.
            endif

            if(index(list(i)%s,'z^')>0) then
                sublist=split(list(i)%s,o_sep='^')
                call hud('Will weight the image by z^'//sublist(2)%s)
                call self%by_depth(o_power=str2real(sublist(2)%s))
            endif
            
            if (index(list(i)%s,'custom')>0) then
                sublist=split(list(i)%s,o_sep=':')
                if(size(sublist)==1) then !filename is not attached, ask for it
                    file=setup%get_file('FILE_IMAGE_WEIGHT_CUSTOM',o_mandatory=1)
                else !filename is attached
                    file=sublist(2)%s
                endif

                call hud('Will weight the image in a custom way defined in '//file)
                call self%by_custom(file)
            endif

        enddo

        if(mpiworld%is_master) call sysio_write('Imag_weights',self%weight,size(self%weight))
        
    end subroutine

    subroutine by_depth(self,o_power,o_factor)
        class(t_image_weighter) :: self
        real,optional :: o_power, o_factor

        real,dimension(:),allocatable :: tmp

        if(present(o_power)) then
            call alloc(tmp,m%nz,o_init=1.)

            tmp=[(((iz-1)*m%dz)**o_power,iz=1,m%nz)]
            do iy=1,m%ny; do ix=1,m%nx
                self%weight(:,ix,iy)=tmp(:)
            enddo; enddo

            deallocate(tmp)

        endif

        ! if(present(o_factor)) then
        !     preco_in_m(:,1,1)=[(((iz-1)*m%dz)*o_factor,iz=1,m%nz)]
        !     do iy=1,m%ny; do ix=1,m%nx
        !         preco_in_m(:,ix,iy)=preco_in_m(:,1,1)
        !     enddo; enddo
        ! endif
        
    end subroutine

    subroutine by_custom(self,file)
        class(t_image_weighter) :: self
        character(*) :: file

        integer :: file_size
        real,dimension(:,:,:),allocatable :: tmp

        inquire(file=file,size=file_size)
        if(file_size<4*self%nz*self%nx*self%ny) call warn('FILE_IMAGE_WEIGHT_CUSTOM has size '//num2str(file_size/4)//' < nz*nx*ny. Some part of the image will not be weighted.')

        call alloc(tmp,self%nz,self%nx,self%ny,o_init=1.)

        call sysio_read(file,tmp,size(tmp))

        !self%weight=self%weight*tmp
        self%weight=tmp

        deallocate(tmp)

    end subroutine

end
