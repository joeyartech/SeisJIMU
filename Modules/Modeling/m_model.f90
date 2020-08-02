module m_model
use m_mpienv
use m_setup
use m_arrayop

    type t_model
        integer :: nx,ny,nz, n
        real    :: dx,dy,dz
        real    :: ox,oy,oz
        logical :: is_cubic, is_isotropic, if_freesurface
        real,dimension(:,:,:),allocatable :: vp,vs,rho,eps,del,eta
        real,dimension(:,:),allocatable :: topo
        integer,dimension(:,:),allocatable :: itopo
        real,dimension(:,:,:),allocatable :: mask_vp,mask_vs,mask_rho
                
        real ref_vp,ref_vs,ref_rho

        real velmin, velmax
        
        real cell_volume, cell_diagonal, cell_inv_diagonal
        
    end type
    
    type(t_model) :: m

    contains

    subroutine model_estim_RAM
        character(:),allocatable :: tmp

        ! tmp=setup_get_char('MODEL_SIZE','NZXY')
        ! read(tmp1,*) m%nz, m%nx, m%ny
        ! m%n=m%nx*m%ny*m%nz

    end subroutine

    
    subroutine model_read
        character(:),allocatable :: tmp1,tmp2,tmp3,tmp4
        logical alive
        
        tmp1=setup_get_char('MODEL_SIZE','NZXY')
        read(tmp1,*) m%nz, m%nx, m%ny
        m%n=m%nx*m%ny*m%nz
        
        tmp2=setup_get_char('MODEL_SPACING','DZXY')
        read(tmp2,*) m%dz, m%dx, m%dy
        
        tmp3=setup_get_char('MODEL_ORIGIN','OZXY')
        read(tmp3,*) m%oz, m%ox, m%oy
        
        if(m%ny==1) then
            call hud('2D geometry')
            m%is_cubic=.false.
            m%dy=1.
        else
            call hud('3D geometry')
            m%is_cubic=.true.
        endif
        
        m%cell_volume = m%dx*m%dy*m%dz
        
        m%cell_diagonal=sqrt(m%dx**2+m%dz**2)
        m%cell_inv_diagonal=sqrt(m%dx**(-2) + m%dz**(-2))

        if(m%is_cubic) then
            m%cell_diagonal=sqrt(m%dx**2+m%dy**2+m%dz**2)
            m%cell_inv_diagonal=sqrt(m%dx**(-2) + m%dy**(-2) + m%dz**(-2))
        endif


        n=4*m%nx*m%ny*m%nz
        
        !read models
        tmp4=setup_get_char('FILE_MODELS')
        
        !vp
        inquire(file=tmp4//'%vp', exist=alive)
        if(.not.alive) call error('Unable to find vp model!')
        call alloc(m%vp,m%nz,m%nx,m%ny)
        open(12,file=tmp4//'%vp',access='direct',recl=n,action='read',status='old')
        read(12,rec=1) m%vp 
        close(12)
        call hud('vp model is read.')
        
        !vs
        inquire(file=tmp4//'%vs', exist=alive)
        if(alive) then
            call alloc(m%vs,m%nz,m%nx,m%ny)
            open(12,file=tmp4//'%vs',access='direct',recl=n,action='read',status='old')
            read(12,rec=1) m%vs
            close(12)
            call hud('vs model is read.')
        else
            call alloc(m%vs,1,1,1)
        endif
        
        !rho
        inquire(file=tmp4//'%rho', exist=alive)
        if(alive) then
            call alloc(m%rho,m%nz,m%nx,m%ny)
            open(12,file=tmp4//'%rho',access='direct',recl=n,action='read',status='old')
            read(12,rec=1) m%rho
            close(12)
            call hud('rho model is read.')
        else
            call alloc(m%rho,1,1,1); m%rho=1000. !in [kg/m3]
        endif
        
        !epsilon
        inquire(file=tmp4//'%eps', exist=alive)
        if(alive) then
            m%is_isotropic=.false.
            call alloc(m%eps,m%nz,m%nx,m%ny)
            open(12,file=tmp4//'%eps',access='direct',recl=n,action='read',status='old')
            read(12,rec=1) m%eps
            close(12)
            call hud('eps model is read.')
        else
            call alloc(m%eps,1,1,1)
        endif
        
        !delta
        inquire(file=tmp4//'%del', exist=alive)
        if(alive) then
            m%is_isotropic=.false.
            call alloc(m%del,m%nz,m%nx,m%ny)
            open(12,file=tmp4//'%del',access='direct',recl=n,action='read',status='old')
            read(12,rec=1) m%del
            close(12)
            call hud('del model is read.')
        else
            call alloc(m%del,1,1,1)
        endif
        
        m%ref_vp=m%vp(1,1,1)
        m%ref_vs=m%vs(1,1,1)
        m%ref_rho=m%rho(1,1,1)
        
        !refresh isotropicness if user requires
        if(setup_ask('IF_ISOTROPIC')) then
            m%is_isotropic=setup_get_logical('IF_ISOTROPIC')
        endif
        
        !topography
        call alloc(m%topo,m%nx,m%ny)
        call alloc(m%itopo,m%nx,m%ny,initialize=.false.); m%itopo=1

        if(size(m%vs)>1) then
        if(setup_get_logical('IF_TOPO_FROM_VS',default=.true.)) then
            !m%itopo = maxloc(m%vs, dim=1, mask=(m%vs<10), back=.true.)+1 !the "back" argument isn't implemented in gfortran until version 9 ..
            !m%itopo = minloc(m%vs, dim=1, mask=(m%vs>=10.)); where(m%itopo==0) m%itopo=m%nz+1  !still not correct
            do i3=1,m%ny; do i2=1,m%nx
            loop: do i1=1,m%nz
                if(m%vs(i1,i2,i3)>=10) then
                    m%itopo(i2,i3) = i1
                    exit loop
                endif
            enddo loop
            enddo; enddo
            m%topo = (m%itopo-1)*m%dz
        endif
        endif
        
        inquire(file=tmp4//'%topo', exist=alive)
        if(alive) then
            open(12,file=tmp4//'%topo',access='direct',recl=4*m%nx*m%ny,action='read',status='old')
            read(12,rec=1) m%topo
            close(12)
            call hud('topo model is read.')
            m%itopo=nint(m%topo/m%dz)+1
        endif

        if(mpiworld%is_master) then
            write(*,*) 'm%itopo minmax value:', minval(m%itopo), maxval(m%itopo)
        endif
        ! if(mpiworld%is_master) then
        !     open(12,file='itopo_from_vs')
        !     write(12,*) m%itopo
        !     close(12)
        ! endif

            call alloc( m%mask_vp,  maxval(m%itopo),m%nx,m%ny ); m%mask_vp(:,:,:)  = m%vp(1:maxval(m%itopo),:,:)
        if(size(m%vs)>1) then
            call alloc( m%mask_vs,  maxval(m%itopo),m%nx,m%ny ); m%mask_vs(:,:,:)  = m%vs(1:maxval(m%itopo),:,:)
        endif
            call alloc( m%mask_rho, maxval(m%itopo),m%nx,m%ny ); m%mask_rho(:,:,:) = m%rho(1:maxval(m%itopo),:,:)
        
        !freesurface
        m%if_freesurface=setup_get_logical('IF_FREESURFACE',default=.true.)
        
    end subroutine
    
    subroutine model_realloc(cmodel)
        character(*) :: cmodel
        
        select case (cmodel)
            case ('rho')
                if(size(m%rho)==1) then
                    deallocate(m%rho)
                    call alloc(m%rho,m%nz,m%nx,m%ny,initialize=.false.)
                    m%rho=1000.
                endif
        end select
        
    end subroutine
end
