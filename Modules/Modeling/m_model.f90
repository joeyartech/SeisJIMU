module m_model
use m_System

    private

    type,public :: t_model
        integer :: nx,ny,nz, n
        real    :: dx,dy,dz
        real    :: ox,oy,oz

        real :: cell_volume, cell_diagonal, cell_inv_diagonal
                
        character(:),allocatable :: file
        type(t_string),dimension(:),allocatable :: attributes_read, attributes_write

        logical :: is_cubic, is_freesurface, if_has_prior
        
        real,dimension(:,:,:),allocatable :: vp,vs,rho
        real,dimension(:,:,:),allocatable :: eps,del,eta
        real,dimension(:,:,:),allocatable :: qp,qs

        real,dimension(:,:,:),allocatable :: vp_prior,vs_prior,rho_prior

        real :: ref_vel=1500., ref_rho=1000., ref_kpa=1000.*1500.**2

        logical,dimension(:,:,:),allocatable :: is_freeze_zone

        real,dimension(:,:,:,:),allocatable :: gradient

        contains
        procedure :: init => init
        procedure :: estim_RAM => estim_RAM
        procedure :: read => read
        procedure :: read_prior => read_prior
        procedure :: write => write
        procedure :: apply_freeze_zone => apply_freeze_zone

    end type
    
    type(t_model),public :: m

    contains

    subroutine estim_RAM(self)
        class(t_model) :: self

    end subroutine

    subroutine init(self)
        class(t_model) :: self

        integer,dimension(3) :: itmp=[1,1,1]
        real,dimension(3) :: rtmp=[1.,1.,1.]

        itmp=setup%get_ints('MODEL_SIZE','NZXY',o_mandatory=.true.)
        
        self%nz=itmp(1); self%nx=itmp(2); self%ny=itmp(3)

        self%n=self%nx*self%ny*self%nz

        if(self%ny==1) then
            call hud('2D geometry')
            self%is_cubic=.false.
            self%dy=1.
        else
            call hud('3D geometry')
            self%is_cubic=.true.
        endif

        rtmp=setup%get_reals('MODEL_SPACING','DZXY',o_mandatory=.true.)
        self%dz=rtmp(1); self%dx=rtmp(2); self%dy=rtmp(3)

        rtmp=setup%get_reals('MODEL_ORIGIN','OZXY',o_default='0. 0. 0.')
        self%oz=rtmp(1); self%ox=rtmp(2); self%oy=rtmp(3)

        self%cell_volume = self%dz*self%dx*self%dy
        if(self%is_cubic) then
            self%cell_diagonal=sqrt(self%dz**2+self%dx**2+self%dy**2)
            self%cell_inv_diagonal=sqrt(self%dz**(-2) + self%dx**(-2) + self%dy**(-2))
        else
            self%cell_diagonal=sqrt(self%dz**2+self%dx**2)
            self%cell_inv_diagonal=sqrt(self%dz**(-2) + self%dx**(-2))
        endif

        self%file=setup%get_file('FILE_MODEL',o_mandatory=.true.) 

        self%attributes_read =setup%get_strs('MODEL_ATTRIBUTES',o_default='vp rho')
        self%attributes_write=setup%get_strs('MODEL_ATTRIBUTES_WRITE',o_default='vp')
        
    end subroutine

    subroutine read(self)
        class(t_model) :: self

        real,dimension(:,:,:),allocatable :: tmp
        character(:),allocatable :: file

        open(12,file=self%file,access='direct',recl=4*self%n,action='read',status='old')
        
        do i=1,size(self%attributes_read)
            select case(self%attributes_read(i)%s)
            case ('vp')
                call alloc(self%vp,self%nz,self%nx,self%ny)
                read(12,rec=i) self%vp
                call hud('vp model is read.')
                self%ref_vel=self%vp(1,1,1)

            case ('vs')
                call alloc(self%vs,self%nz,self%nx,self%ny)
                read(12,rec=i) self%vs
                call hud('vs model is read.')

            case ('rho')
                call alloc(self%rho,self%nz,self%nx,self%ny)
                read(12,rec=i) self%rho
                call hud('rho model is read.')
                self%ref_rho=self%rho(1,1,1)

            case ('eps')
                call alloc(self%eps,self%nz,self%nx,self%ny)
                read(12,rec=i) self%eps
                call hud('eps model is read.')

            case ('del')
                call alloc(self%del,self%nz,self%nx,self%ny)
                read(12,rec=i) self%del
                call hud('del model is read.')

            case ('eta')
                call alloc(self%eta,self%nz,self%nx,self%ny)
                read(12,rec=i) self%eta
                call hud('eta model is read.')

            case ('qp')
                call alloc(self%qp,self%nz,self%nx,self%ny)
                read(12,rec=i) self%qp
                call hud('qp model is read.')

            case ('qs')
                call alloc(self%qs,self%nz,self%nx,self%ny)
                read(12,rec=i) self%qs
                call hud('qs model is read.')

            end select

        enddo

        close(12)

        self%ref_kpa=self%ref_rho*self%ref_vel**2

        !freesurface
        self%is_freesurface=setup%get_bool('IS_FREESURFACE',o_default='T')
                
        !freeze zone
        allocate(self%is_freeze_zone(self%nz,self%nx,self%ny),source=.false.)

        !check file topo
        file=setup%get_file('FILE_TOPO')
        if(file/='') then
            call alloc(tmp,self%nx,self%ny,1)
            open(12,file=file,access='direct',recl=4*self%n,action='read')
            read(12,rec=1) tmp
            close(12)
            call hud('topo minmax value:'//num2str(minval(tmp))//num2str(maxval(tmp)))

            do iy=1,self%ny; do ix=1,self%nx
                self%is_freeze_zone(1:nint(tmp(ix,iy,1)/self%dz)+1,ix,iy)=.true.
            enddo; enddo
            call hud('Freeze zone is set from FILE_TOPO.')
        endif

        !2nd check vs model
        if(allocated(self%vs)) then
            if(setup%get_bool('IF_TOPO_FROM_VS',o_default='T')) then
                !self%itopo = maxloc(self%vs, dim=1, mask=(self%vs<10), back=.true.)+1 !the "back" argument isn't implemented in gfortran until version 9 ..
                do iy=1,self%ny; do ix=1,self%nx
                    loopz: do iz=1,self%nz
                        if(self%vs(iz,ix,iy)<10.) then
                            self%is_freeze_zone(iz,ix,iy) = .true.
                        else
                            exit loopz
                        endif
                    enddo loopz
                enddo; enddo
                call hud('Freeze zone is additionally set from vs model.')
            endif
        endif

        !3rd check if file freeze is given
        file=setup%get_file('FILE_FREEZE_ZONE')
        if(file/='') then
            call alloc(tmp,self%nz,self%nx,self%ny)
            open(12,file=file,access='direct',recl=4*self%n,action='read')
            read(12,rec=1) tmp
            close(12)
            
            where(tmp==0.) self%is_freeze_zone=.true.
            call hud('Freeze zone is additionally set from FILE_FREEZE_ZONE.')
        endif

        if(allocated(tmp)) deallocate(tmp)

    end subroutine

    subroutine read_prior(self)
        class(t_model) :: self

        character(:),allocatable :: file
        type(t_string),dimension(:),allocatable :: attributes_read

        file=setup%get_file('FILE_PRIOR_MODEL')

        self%if_has_prior=either(.false.,.true.,file=='')

        if(.not.self%if_has_prior) return

        attributes_read =setup%get_strs('PRIOR_MODEL_ATTRIBUTES',o_default='vp rho')

        open(12,file=file,access='direct',recl=4*self%n,action='read',status='old')

        do i=1,size(attributes_read)
            select case(attributes_read(i)%s)
            case ('vp')
                call alloc(self%vp_prior,self%nz,self%nx,self%ny)
                read(12,rec=i) self%vp_prior
                call hud('vp prior model is read.')
                
            case ('vs')
                call alloc(self%vs_prior,self%nz,self%nx,self%ny)
                read(12,rec=i) self%vs_prior
                call hud('vs prior model is read.')

            case ('rho')
                call alloc(self%rho_prior,self%nz,self%nx,self%ny)
                read(12,rec=i) self%rho_prior
                call hud('rho prior model is read.')
                
            end select

        enddo

        close(12)

    end subroutine

    subroutine apply_freeze_zone(self)
        class(t_model) :: self

        real,dimension(:,:,:),allocatable :: tmp

        call alloc(tmp,self%nz,self%nx,self%ny)

        open(12,file=self%file,access='direct',recl=4*self%n,action='read',status='old')

        do i=1,size(self%attributes_read)
            select case(self%attributes_read(i)%s)
            case ('vp')        
                read(12,rec=i) tmp
                where(self%is_freeze_zone) self%vp=tmp

            case ('vs')
                read(12,rec=i) tmp
                where(self%is_freeze_zone) self%vs=tmp

            case ('rho')
                read(12,rec=i) tmp
                where(self%is_freeze_zone) self%rho=tmp

            case ('eps')
                read(12,rec=i) tmp
                where(self%is_freeze_zone) self%eps=tmp

            case ('del')
                read(12,rec=i) tmp
                where(self%is_freeze_zone) self%del=tmp

            case ('eta')
                read(12,rec=i) tmp
                where(self%is_freeze_zone) self%eta=tmp

            case ('qp')
                read(12,rec=i) tmp
                where(self%is_freeze_zone) self%qp=tmp

            case ('qs')
                read(12,rec=i) tmp
                where(self%is_freeze_zone) self%qs=tmp

            end select

        enddo

        close(12)

        deallocate(tmp)

    end subroutine

    subroutine write(self,o_suffix)
        class(t_model) :: self
        character(*),optional :: o_suffix

        character(:),allocatable :: suf

        suf=either(o_suffix,'',present(o_suffix))

        open(13,file=dir_out//'model'//suf,access='direct',recl=4*self%n,action='write')
        
        do i=1,size(self%attributes_write)
            select case(self%attributes_write(i)%s)
            case ('vp')
                if(allocated(self%vp)) then
                    write(13,rec=i) self%vp
                    ! call hud('vp model is written.')
                endif

            case ('vs')
                if(allocated(self%vp)) then
                    write(13,rec=i) self%vs
                    ! call hud('vs model is written.')
                endif

            case ('rho')
                if(allocated(self%rho)) then
                    write(13,rec=i) self%rho
                    ! call hud('rho model is written.')
                endif

            case ('ip')
                if(allocated(self%vp).and.allocated(self%rho)) then
                    write(13,rec=i) self%vp*self%rho
                    ! call hud('ip model is written.')
                endif

            case ('is')
                if(allocated(self%vs).and.allocated(self%rho)) then
                    write(13,rec=i) self%vs*self%rho
                    ! call hud('is model is written.')
                endif

            case ('eps')
                if(allocated(self%eps)) then
                    write(13,rec=i) self%eps
                    ! call hud('eps model is written.')
                endif

            case ('del')
                if(allocated(self%del)) then
                    write(13,rec=i) self%del
                    ! call hud('del model is written.')
                endif

            case ('eta')
                if(allocated(self%eta)) then
                    write(13,rec=i) self%eta
                    ! call hud('eta model is written.')
                endif

            case ('qp')
                if(allocated(self%qp)) then
                    write(13,rec=i) self%qp
                    ! call hud('qp model is written.')
                endif

            case ('qs')
                if(allocated(self%qs)) then
                    write(13,rec=i) self%qs
                    ! call hud('qs model is written.')
                endif

            end select

        enddo

        close(13)

    end subroutine

end
