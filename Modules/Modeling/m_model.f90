module m_model
use m_mpienv
use m_setup
use m_arrayop

    real,dimension(:,:),allocatable :: topo

    type t_model
        integer :: nx,ny,nz, n
        real    :: dx,dy,dz
        real    :: ox,oy,oz
        
        real cell_volume, cell_diagonal, cell_inv_diagonal
        
        character(:),allocatable :: file
        type(t_string),dimension(:),allocatable :: attributes_read, attributes_write

        logical :: is_cubic, is_freesurface
        
        real,dimension(:,:,:),allocatable :: vp,vs,rho
        real,dimension(:,:,:),allocatable :: eps,del,eta
        real,dimension(:,:,:),allocatable :: qp,qs

        real :: ref_vel, ref_rho, ref_kpa

        logical,dimension(:,:,:),allocatable :: is_freeze_zone

    end type
    
    type(t_model) :: m

    contains

    subroutine estim_RAM
        character(:),allocatable :: tmp

        ! tmp=setup_get_char('MODEL_SIZE','NZXY')
        ! read(tmp1,*) m%nz, m%nx, m%ny
        ! m%n=m%nx*m%ny*m%nz

    end subroutine

    subroutine init        
        itmp=setup%get_ints('MODEL_SIZE','NZXY',mandatory=.true.)
        m%nz=itmp(1); m%nx=itmp(2); m%ny=itmp(3)
        m%n=m%nx*m%ny*m%nz

        if(m%ny==1) then
            call hud('2D geometry')
            m%is_cubic=.false.
            m%dy=1.
        else
            call hud('3D geometry')
            m%is_cubic=.true.
        endif

        rtmp=setup%get_reals('MODEL_SPACING','DZXY',mandatory=.true.)
        m%dz=rtmp(1); m%dx=rtmp(2); m%dy=rtmp(3)

        m%cell_volume = m%dx*m%dy*m%dz
        m%cell_diagonal=sqrt(m%dx**2+m%dz**2)
        m%cell_inv_diagonal=sqrt(m%dx**(-2) + m%dz**(-2))

        if(m%is_cubic) then
            m%cell_diagonal=sqrt(m%dx**2+m%dy**2+m%dz**2)
            m%cell_inv_diagonal=sqrt(m%dx**(-2) + m%dy**(-2) + m%dz**(-2))
        endif

        rtmp=setup%get_reals('MODEL_ORIGIN','OZXY',default=[0.,0.,0.])
        m%oz=rtmp(1); m%ox=rtmp(2); m%oy=rtmp(3)
        
        m%file=setup%get_file('FILE_MODEL',mandatory=.true.) 

        self%attributes_read =setup%get_strs('MODEL_ATTRIBUTES',default='vp rho')
        self%attributes_write=setup%get_strs('WRITE_MODEL_ATTRIBUTES',default='vp')
        
    end subroutine

    subroutine read
        real,dimension(:,:,:),allocatable :: tmp

        open(12,file=m%file,access='direct',recl=4*m%n,action='read',status='old')
        
        do i=1,size(self%attributes_read)
            select case(self%attributes_read%s)
            case ('vp')
                call alloc(m%vp,m%nz,m%nx,m%ny)
                read(12,rec=i) m%vp
                call hud('vp model is read.')
                m%ref_vel=m%vp(1,1,1)

            case ('vs')
                call alloc(m%vs,m%nz,m%nx,m%ny)
                read(12,rec=i) m%vs
                call hud('vs model is read.')

            case ('rho')
                call alloc(m%rho,m%nz,m%nx,m%ny)
                read(12,rec=i) m%rho
                call hud('rho model is read.')
                m%ref_rho=m%rho(1,1,1)
                m%ref_kpa=m%ref_rho*m%ref_vel**2

            case ('eps')
                call alloc(m%eps,m%nz,m%nx,m%ny)
                read(12,rec=i) m%eps
                call hud('eps model is read.')

            case ('del')
                call alloc(m%del,m%nz,m%nx,m%ny)
                read(12,rec=i) m%del
                call hud('del model is read.')

            case ('eta')
                call alloc(m%eta,m%nz,m%nx,m%ny)
                read(12,rec=i) m%eta
                call hud('eta model is read.')

            case ('qp')
                call alloc(m%qp,m%nz,m%nx,m%ny)
                read(12,rec=i) m%qp
                call hud('qp model is read.')

            case ('qs')
                call alloc(m%qs,m%nz,m%nx,m%ny)
                read(12,rec=i) m%qs
                call hud('qs model is read.')
            end select

        enddo

        close(12)


        !freesurface
        m%is_freesurface=setup_get_logical('IF_FREESURFACE',default=.true.)
                
        !freeze zone
        call alloc(m%is_freeze_zone(m%nz,m%nx,m%ny),initialize=.false.)
        m%is_freeze_zone=.false.

        !check file topo
        associate(file=>setup%get_file('FILE_TOPO'))
            if(file/='') then
                call alloc(tmp,m%nx,m%ny,1)
                open(12,file=file,access='direct',recl=4*m%n,action='read')
                read(12,rec=1) tmp
                close(12)
                if(mpiworld%is_master) write(*,*) 'topo minmax value:', minval(topo), maxval(topo)

                do iy=1,m%ny; do ix=1,m%nx
                    m%is_freeze_zone(1:nint(tmp(ix,iy,1)/dz)+1,ix,iy)=.true.
                enddo; enddo
                call hud('Freeze zone is set from FILE_TOPO.')
            endif
        end associate

        !2nd check vs model
        if(allocated(m%vs)>1) then
            if(setup_get_logical('IF_TOPO_FROM_VS',default=.true.)) then
                !m%itopo = maxloc(m%vs, dim=1, mask=(m%vs<10), back=.true.)+1 !the "back" argument isn't implemented in gfortran until version 9 ..
                do iy=1,m%ny; do ix=1,m%nx
                    loopz: do iz=1,m%nz
                        if(m%vs(iz,ix,iy)<10.) then
                            m%is_freeze_zone(iz,ix,iy) = .true.
                        else
                            exit loopz
                        endif
                    enddo
                enddo; enddo
                call hud('Freeze zone is additionally set from vs model.')
            endif
        endif

        !3rd check if file freeze is given
        associate(file=>setup%get_file('FILE_FREEZE_ZONE'))
            if(file/='') then
                call alloc(tmp,m%nz,m%nx,m%ny)
                open(12,file=file,access='direct',recl=4*m%n,action='read')
                read(12,rec=1) tmp
                close(12)
                
                where(tmp==0.) m%is_freeze_zone=.true.
                call hud('Freeze zone is additionally set from FILE_FREEZE_ZONE.')
            endif
        end associate

        if(allocated(deallocated(tmp)))
        
        
    end subroutine

    subroutine apply_freeze_zone
        real,dimension(:,:,:),allocatable :: tmp

        call alloc(tmp(m%nz,m%nx,m%ny))

        open(12,file=m%file,access='direct',recl=4*m%n,action='read',status='old')

        do i=1,size(self%attributes_read)
            select case(self%attributes_read%s)
            case ('vp')        
                read(12,rec=i) tmp
                where(m%is_freeze_zone) m%vp=tmp

            case ('vs')
                read(12,rec=i) tmp
                where(m%is_freeze_zone) m%vs=tmp

            case ('rho')
                read(12,rec=i) tmp
                where(m%is_freeze_zone) m%rho=tmp

            case ('eps')
                read(12,rec=i) tmp
                where(m%is_freeze_zone) m%eps=tmp

            case ('del')
                read(12,rec=i) tmp
                where(m%is_freeze_zone) m%del=tmp

            case ('eta')
                read(12,rec=i) tmp
                where(m%is_freeze_zone) m%eta=tmp

            case ('qp')
                read(12,rec=i) tmp
                where(m%is_freeze_zone) m%qp=tmp

            case ('qs')
                read(12,rec=i) tmp
                where(m%is_freeze_zone) m%qs=tmp

            end select

        enddo

        close(12)

        deallocate(tmp)

    end subroutine

    subroutine write(o_suffix)
        character(*),optional :: o_suffix
        character(:),allocatable :: suffix

        if(present(o_suffix)) then
            suffix=o_suffix
        else
            suffix=''
        endif

        open(13,file=setup%dir_working//'model'//suffix,access='direct',recl=4*m%n,action='write')
        
        do i=1,size(self%attributes_write)
            select case(self%attributes_write%s)
            case ('vp')
                if(allocated(m%vp)) then
                    write(13,rec=i) m%vp
                    call hud('vp model is written.')
                endif

            case ('vs')
                if(allocated(m%vp)) then
                    write(13,rec=i) m%vs
                    call hud('vs model is written.')
                endif

            case ('rho')
                if(allocated(m%rho)) then
                    write(13,rec=i) m%rho
                    call hud('rho model is written.')
                endif

            case ('ip')
                if(allocated(m%vp).and.allocated(m%rho)) then
                    write(13,rec=i) m%vp*m%rho
                    call hud('ip model is written.')
                endif

            case ('is')
                if(allocated(m%vs).and.allocated(m%rho)) then
                    write(13,rec=i) m%vs*m%rho
                    call hud('is model is written.')
                endif

            case ('eps')
                if(allocated(m%eps)) then
                    write(13,rec=i) m%eps
                    call hud('eps model is written.')
                endif

            case ('del')
                if(allocated(m%del)) then
                    write(13,rec=i) m%del
                    call hud('del model is written.')
                endif

            case ('eta')
                if(allocated(m%eta)) then
                    write(13,rec=i) m%eta
                    call hud('eta model is written.')
                endif

            case ('qp')
                if(allocated(m%qp)) then
                    write(13,rec=i) m%qp
                    call hud('qp model is written.')
                endif

            case ('qs')
                if(allocated(m%qs)) then
                    write(13,rec=i) m%qs
                    call hud('qs model is written.')
                endif

            end select

        enddo

        close(13)

    end subroutine

end
