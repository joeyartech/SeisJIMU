module m_correlation
use m_System
use m_model
use m_computebox

    private

    !correlation
    type,public :: t_correlation

        character(:),allocatable :: name

        !dimension
        integer :: nz,nx,ny,nt

        !index
        integer :: ifz,ilz,ifx,ilx,ify,ily,ift,ilt

        !wavefield components in computation domain
        real,dimension(:,:,:,:),allocatable :: sp_rp
        !real,dimension(:,:,:,:),allocatable :: svz_rstress

        real,dimension(:,:,:,:),allocatable :: drp_dt_dsp_dt, nab_rp_nab_sp
        real,dimension(:,:,:,:),allocatable :: rp_lap_sp

        contains

        procedure :: init
        ! procedure :: init_bounds
        ! procedure :: check_value
        ! procedure :: change_dim
        procedure :: write
        final :: final
        
        ! procedure :: is_registered
        ! procedure :: register

    end type

    real,dimension(:,:,:,:),allocatable,public :: correlation_image, correlation_gradient

    ! !image, as a child of correlation
    ! type,public :: t_image
    ! end type

    ! !gradient, as a child of correlation
    ! type,public :: t_gradient
    ! end type

    contains

    ! subroutine init_bounds(self,name,nz,nx,ny,nt,o_init)
    !     class(t_correlation) :: self
    !     character(*) :: name
    !     integer,dimension(2) :: nz,nx,ny,nt
    !     real,optional :: o_init

    !     self%name=name

    !     call alloc(self%sp_rp,nz,nx,ny,nt,o_init=o_init)
    
    ! 	self%ifz=nz(1); self%ilz=nz(2); self%nz=nz(2)-nz(1)+1
    ! 	self%ifx=nx(1); self%ilx=nx(2); self%nx=nx(2)-nx(1)+1
    ! 	self%ify=ny(1); self%ily=ny(2); self%ny=ny(2)-ny(1)+1
    ! 	self%ift=nt(1); self%ilt=nt(2); self%nt=nt(2)-nt(1)+1

    ! end subroutine

    !subroutine init_shapes(self,name,nz,nx,ny,nt,o_init)
    subroutine init(self,name,shape123_from,o_dim4,o_init)
        class(t_correlation) :: self
        character(*) :: name
        character(*) :: shape123_from
        integer,dimension(2),optional :: o_dim4
        real,optional :: o_init

        self%name=name

        select case (shape123_from)
            case ('model')
            self%ifz=1; self%ilz=m%nz; self%nz=m%nz
            self%ifx=1; self%ilx=m%nx; self%nx=m%nx
            self%ify=1; self%ily=m%ny; self%ny=m%ny

            case ('computebox')
            self%ifz=cb%ifz; self%ilz=cb%ilz; self%nz=cb%nz
            self%ifx=cb%ifx; self%ilx=cb%ilx; self%nx=cb%nx
            self%ify=cb%ify; self%ily=cb%ily; self%ny=cb%ny

        end select

        if(present(o_dim4)) then
            self%ift=o_dim4(1); self%ilt=o_dim4(2); self%nt=o_dim4(2)-o_dim4(1)+1
        else
            self%ift=1;         self%ilt=1;         self%nt=1
        endif

        call alloc(self%sp_rp,[self%ifz,self%ilz],[self%ifx,self%ilx],[self%ify,self%ily],[self%ift,self%ilt],o_init=o_init)

    end subroutine

    ! subroutine change_dim(self,new_nz,new_nx,new_ny,new_nt)
    ! 	class(t_correlation) :: self
    ! 	integer,dimension(2) :: new_nz,new_nx,new_ny,new_nt

    ! 	real,dimension(:,:,:,:),allocatable :: tmp

    ! 	call alloc(tmp,[self%ifz,self%ilz],[self%ifx,self%ilx],[self%ify,self%ily],[self%ift,self%ilt],o_init=self%sp_rp)

    ! 	call alloc(self%sp_rp,new_nz,new_nx,new_ny,new_nt)
    
    ! 	self%sp_rp(self%ifz:self%ilz,self%ifx:self%ilx,self%ify:self%ily,self%ift:self%ilt) = tmp

    ! 	deallocate(tmp)

    ! 	self%ifz=new_nz(1); self%ilz=new_nz(2); self%nz=new_nz(2)-new_nz(1)+1
    ! 	self%ifx=new_nx(1); self%ilx=new_nx(2); self%nx=new_nx(2)-new_nx(1)+1
    ! 	self%ify=new_ny(1); self%ily=new_ny(2); self%ny=new_ny(2)-new_ny(1)+1
    ! 	self%ift=new_nt(1); self%ilt=new_nt(2); self%nt=new_nt(2)-new_nt(1)+1

    ! end subroutine

    ! subroutine project_to(self,big)
    ! 	class(t_correlation) :: self
    ! 	type(t_correlation) :: big

    ! 	big%sp_rp(self%ifz:self%ilz,self%ifx:self%ilx,self%ify:self%ily,self%ift:self%ilt) = self%sp_rp

    ! end subroutine

    subroutine write(self,o_prefix,o_suffix,o_mode)
        class(t_correlation) :: self
        character(*),optional :: o_prefix,o_suffix,o_mode

        character(:),allocatable :: prf,suf

        prf=either(o_prefix,'',present(o_prefix))
        suf=either(o_suffix,'',present(o_suffix))

        if(present(o_mode)) then
            call sysio_write(prf//self%name//suf,self%sp_rp,size(self%sp_rp),o_mode=o_mode)
        else
            call sysio_write(prf//self%name//suf,self%sp_rp,size(self%sp_rp))
        endif

    end subroutine

    subroutine final(self)
        type(t_correlation) :: self

        call dealloc(self%sp_rp)

    end subroutine

end