module m_model
use m_System
use m_math
use m_smoother_laplacian_sparse

    private

    type,public :: t_model
        integer :: nx,ny,nz, n
        real    :: dx,dy,dz, dmin
        real    :: ox,oy,oz

        real :: cell_volume, cell_diagonal, rev_cell_diagonal
                
        character(:),allocatable :: file
        type(t_string),dimension(:),allocatable :: attributes_read, attributes_write

        logical :: is_cubic, is_freesurface, if_has_prior
        
        real,dimension(:,:,:),allocatable :: vp,vs,rho
        ! real,dimension(:,:,:),allocatable :: eps,del,eta
        real,dimension(:,:,:),allocatable :: qp,qs

        real,dimension(:,:,:),allocatable :: epsr,mur,sgma


        !prior models
        real,dimension(:,:,:),allocatable :: vp_prior,vs_prior,rho_prior

        !reference values
        real :: ref_inv_vel, ref_rho

        integer,dimension(:,:),allocatable :: ibathy
        logical,dimension(:,:,:),allocatable :: is_freeze_zone

        contains
        procedure :: init
        procedure :: estim_RAM
        procedure :: read
        procedure :: read_prior
        procedure :: set_reference
        procedure :: write
        procedure :: apply_freeze_zone
        procedure :: apply_relativity
        ! procedure :: apply_elastic_continuum

    end type
    
    type(t_model),public :: m

    contains

    subroutine estim_RAM(self)
        class(t_model) :: self

        call hud('stack size:')
        if(mpiworld%is_master) call execute_command_line('ulimit -s', wait=.true.)

    end subroutine

    subroutine init(self)
        class(t_model) :: self

        integer,dimension(3) :: itmp=1
        real,dimension(3) :: rtmp=1.

        itmp=setup%get_ints('MODEL_SIZE','NZXY',o_mandatory=3)
        
        self%nz=itmp(1); self%nx=itmp(2); self%ny=itmp(3)

        self%n=self%nz*self%nx*self%ny

        rtmp=setup%get_reals('MODEL_SPACING','DZXY',o_mandatory=3)
        self%dz=rtmp(1); self%dx=rtmp(2); self%dy=rtmp(3)

        if(self%ny==1) then
            call hud('2D geometry')
            self%is_cubic=.false.
            self%dy=1.
        else
            call hud('3D geometry')
            self%is_cubic=.true.
        endif

        self%dmin=min(self%dz,self%dx)
        if(self%is_cubic) self%dmin=min(self%dmin,self%dy)

        rtmp=setup%get_reals('MODEL_ORIGIN','OZXY',o_default='0. 0. 0.')
        rtmp=0.
        self%oz=rtmp(1); self%ox=rtmp(2); self%oy=rtmp(3)

        !discretization of the model
        self%cell_volume = self%dz*self%dx*self%dy
        if(self%is_cubic) then
            self%cell_diagonal=sqrt(self%dz**2+self%dx**2+self%dy**2)
            self%rev_cell_diagonal=sqrt(self%dz**(-2) + self%dx**(-2) + self%dy**(-2))
        else
            self%cell_diagonal=sqrt(self%dz**2+self%dx**2)
            self%rev_cell_diagonal=sqrt(self%dz**(-2) + self%dx**(-2))
        endif

        self%file=setup%get_file('FILE_MODEL') 

        self%attributes_read =setup%get_strs('MODEL_ATTRIBUTES',o_default='eps mu sgma')
        self%attributes_write=setup%get_strs('MODEL_ATTRIBUTES_WRITE',o_default=strcat(self%attributes_read))
        
        ! !convert from second to microsecond: F/m=A·s/V/m to A·µs/V/m, H/m=V·s/A/m to V·µs/A/m
        ! !but no need to convert S/m=A²·s³/kg/m² to A²·µs³/kg/m²
        ! self%time_scale=setup%get_real('MODEL_TIME_UNIT',o_default='1e6') 
        ! call hud('If absolute eps, mu, cel models are read, apply temporal scaling: '//num2str(self%time_scale))

    end subroutine

    subroutine read(self)
        class(t_model) :: self

        real,dimension(:,:,:),allocatable :: tmp
        character(:),allocatable :: file, str

        if(self%file=='') then !no read
            ! !freesurface
            ! self%is_freesurface=setup%get_bool('IS_FREESURFACE',o_default='T')
            !freeze zone
            allocate(self%is_freeze_zone(self%nz,self%nx,self%ny),source=.false.)
            return
        endif

        open(12,file=self%file,access='direct',recl=4*self%n,action='read',status='old')
        
        do i=1,size(self%attributes_read)
            select case(self%attributes_read(i)%s)
            ! case ('vp')
            !     call alloc(self%vp,self%nz,self%nx,self%ny)
            !     read(12,rec=i) self%vp
            !     call hud('vp model is read.')

            ! case ('vs')
            !     call alloc(self%vs,self%nz,self%nx,self%ny)
            !     read(12,rec=i) self%vs
            !     call hud('vs model is read.')

            ! case ('rho')
            !     call alloc(self%rho,self%nz,self%nx,self%ny)
            !     read(12,rec=i) self%rho
            !     call hud('rho model is read.')

            ! case ('eps')
            !     call alloc(self%eps,self%nz,self%nx,self%ny)
            !     read(12,rec=i) self%eps
            !     call hud('eps model is read.')

            ! case ('del')
            !     call alloc(self%del,self%nz,self%nx,self%ny)
            !     read(12,rec=i) self%del
            !     call hud('del model is read.')

            ! case ('eta')
            !     call alloc(self%eta,self%nz,self%nx,self%ny)
            !     read(12,rec=i) self%eta
            !     call hud('eta model is read.')

            ! case ('qp')
            !     call alloc(self%qp,self%nz,self%nx,self%ny)
            !     read(12,rec=i) self%qp
            !     call hud('qp model is read.')

            ! case ('qs')
            !     call alloc(self%qs,self%nz,self%nx,self%ny)
            !     read(12,rec=i) self%qs
            !     call hud('qs model is read.')

            case ('eps')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                self%epsr = tmp/r_eps0 !convert from absolute to relative permittivity
                call hud('eps model is read and is converted to epsr.')

            case ('epsr')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                self%epsr = tmp
                call hud('epsr model is read.')

            case ('mu')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                self%mur = tmp/r_mu0 !convert from absolute to relative permeability
                call hud('mu model is read and is converted to mur.')

            case ('mur')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                self%mur = tmp
                call hud('mur model is read.')

            case ('cel')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                !cel**-2 = eps0mu0 * epsr*mur
                self%epsr = tmp**(-2) / r_eps0mu0/self%mur
                call hud('cel model is read and is converted to epsr based on mur.')

            case ('sgma')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                self%sgma = tmp
                call hud('sgma model is read.')

            end select

        enddo

        close(12)

        !freesurface
        self%is_freesurface=setup%get_bool('IS_FREESURFACE',o_default='T')
        
        !bathymetry or topography
        call alloc(self%ibathy,self%nx,self%ny)

        !freeze zone
        allocate(self%is_freeze_zone(self%nz,self%nx,self%ny),source=.false.)

        !check file bathymetry or topography
        file=setup%get_file('FILE_BATHYMETRY','FILE_BATHY')
        if(file=='') file=setup%get_file('FILE_TOPOGRAPHY','FILE_TOPO')
        if(file/='') then
            call alloc(tmp,self%nx,self%ny,1)
            call sysio_read(file,tmp,self%nx*self%ny)
            call hud('bathy minmax value: '//num2str(minval(tmp))//' , '//num2str(maxval(tmp)))
            call hud('water or air layer is from #1 to #(floor(bathy/dz)+1) grid points in depth')

            do iy=1,self%ny; do ix=1,self%nx    
                self%ibathy(ix,iy)=floor(tmp(ix,iy,1)/self%dz)+1
            enddo; enddo

            do iy=1,self%ny; do ix=1,self%nx    
                self%is_freeze_zone(1:self%ibathy(ix,iy),ix,iy)=.true.
            enddo; enddo
            
            call hud('Freeze zone is set from FILE_BATHYMETRY or FILE_TOPOGRAPHY.')
        endif

        !2nd check eps_r model
        if(allocated(self%vs)) then
            if(setup%get_bool('IF_TOPO_FROM_EPS',o_default='T')) then
                !self%ibathy = maxloc(self%vs, dim=1, mask=(self%vs<10), back=.true.)+1 !the "back" argument isn't implemented in gfortran until version 9 ..
                do iy=1,self%ny; do ix=1,self%nx
                    loopz: do iz=1,self%nz
                        if(self%epsr(iz,ix,iy)<1.1) then !air
                            self%is_freeze_zone(iz,ix,iy) = .true.
                        else
                            exit loopz
                        endif
                    enddo loopz
                enddo; enddo
                call hud('Freeze zone is additionally set from epsr model.')
            endif
        endif

        !3rd check if file freeze is given
        file=setup%get_file('FILE_FREEZE_ZONE')
        if(file/='') then
            call alloc(tmp,self%nz,self%nx,self%ny)
            call sysio_read(file,tmp,self%n)
            
            where(tmp==0.) self%is_freeze_zone=.true.
            call hud('Freeze zone is additionally set from FILE_FREEZE_ZONE.')
        endif

        call dealloc(tmp)

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

    subroutine set_reference(self,iz,ix,iy)
    use mpi
        class(t_model) :: self
        real :: tmp(2)
        
        if(mpiworld%is_master) then

            if(setup%check('MODEL_REFERENCE','MREF')) then
                tmp=setup%get_reals('MODEL_REFERENCE','MREF',o_mandatory=2)
                self%ref_inv_vel=1./tmp(1)
                self%ref_rho   =tmp(2)

            else
                self%ref_inv_vel= 1./sqrt(r_eps0mu0*self%epsr(iz,ix,iy)*self%mur(iz,ix,iy))
                self%ref_rho    = self%mur(iz,ix,iy)

            endif

            write(*,*) 'Reference velocity value =',1./self%ref_inv_vel
            write(*,*) 'Reference rho value =',self%ref_rho
        
        endif
        
        call mpi_bcast(self%ref_inv_vel,1,mpi_real,0,mpiworld%communicator,mpiworld%ierr)
        call mpi_bcast(self%ref_rho    ,1,mpi_real,0,mpiworld%communicator,mpiworld%ierr)

        !ref_inv_vel & _rho should also be checkpointed
    end subroutine

    subroutine apply_freeze_zone(self)
        class(t_model) :: self

        real,dimension(:,:,:),allocatable :: tmp

        call alloc(tmp,self%nz,self%nx,self%ny)

        open(12,file=self%file,access='direct',recl=4*self%n,action='read',status='old')

        do i=1,size(self%attributes_read)
            select case(self%attributes_read(i)%s)
            case ('eps')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                where(self%is_freeze_zone) self%epsr=tmp/r_eps0 !convert from absolute to relative permittivity

            case ('epsr')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                where(self%is_freeze_zone) self%epsr=tmp
                
            case ('mu')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                where(self%is_freeze_zone) self%mur = tmp/r_mu0 !convert from absolute to relative permeability

            case ('mur')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                where(self%is_freeze_zone) self%mur=tmp

            case ('sgma')
                call alloc(tmp,self%nz,self%nx,self%ny); read(12,rec=i) tmp
                where(self%is_freeze_zone) self%sgma=tmp

            end select

        enddo

        close(12)

        deallocate(tmp)

    end subroutine

    subroutine apply_elastic_continuum(self)
        class(t_model) :: self

    !     !vp²=(K+4/3*G)/ρ; vs²=G/ρ
    !     !lowest possible vp²=K/ρ=500 m/s
    !     !so vp²  500² + 4/3vs²
    !     !0.75(vp²-500²) >= vs²
        
    !     real,dimension(:,:,:),allocatable :: tmp
        
    !     tmp = sqrt(0.75* (self%vp**2 - 500**2))

    !     where (tmp<self%vs) self%vs=tmp

    !     deallocate(tmp)

    end subroutine

    subroutine apply_relativity(self)
        class(t_model) :: self

        real,dimension(:,:,:),allocatable :: tmp

        tmp = self%epsr*self%mur

        !celerity cannot surpass the speed of light
        where(tmp < 1.) tmp = 1.

        !epsr < 1 is more un-usual than mur < 1 (diamagnetic)
        self%epsr = tmp / self%mur

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
            case ('eps')
                write(13,rec=i) self%epsr*r_eps0

            case ('epsr')
                write(13,rec=i) self%epsr

            case ('mu')
                write(13,rec=i) self%mur*r_mu0

            case ('mur')
                write(13,rec=i) self%mur

            case ('sgma')
                write(13,rec=i) self%sgma

            case ('cel')
                write(13,rec=i) r_c0/sqrt(self%epsr*self%mur)

            end select

        enddo

        close(13)

    end subroutine

end
