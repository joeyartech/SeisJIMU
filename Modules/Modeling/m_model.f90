module m_model
use m_mpienv
use m_sysio
use m_arrayop

    type t_model
        integer :: nx,ny,nz, n
        real    :: dx,dy,dz
        real    :: ox,oy,oz
        logical :: is_cubic, is_isotropic, if_freesurface
        real,dimension(:,:,:),allocatable :: vp,vs,rho,eps,del,eta
        real,dimension(:,:,:),allocatable :: vp2,vs2,rho2 !,eps,del,eta
        real,dimension(:,:),allocatable :: topo
        integer,dimension(:,:),allocatable :: itopo
        real,dimension(:,:,:),allocatable :: vp_mask,vs_mask,rho_mask
                
        real ref_vp,ref_vs,ref_rho

        real velmin, velmax
        
        real cell_volume, cell_diagonal, cell_inv_diagonal
        
    end type
    
    type(t_model) :: m

    contains
    
    subroutine init_model
        character(:),allocatable :: tmp1,tmp2,tmp3,tmp4,tmp5
        logical alive
        
        tmp1=get_setup_char('MODEL_DIMENSION')
        read(tmp1,*) m%nz, m%nx, m%ny
        m%n=m%nx*m%ny*m%nz
        
        tmp2=get_setup_char('MODEL_SPACING')
        read(tmp2,*) m%dz, m%dx, m%dy
        
        tmp3=get_setup_char('MODEL_ORIGIN')
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
        tmp4=get_setup_char('FILE_MODELS')
        tmp5=get_setup_char('FILE_MODELS_2')
        
        !vp
        inquire(file=tmp4//'_vp', exist=alive)
        if(.not.alive) then
            call hud('ERROR: Unable to find vp model!')
            stop
        endif
        call alloc(m%vp,m%nz,m%nx,m%ny)
        open(12,file=tmp4//'_vp',access='direct',recl=n,action='read',status='old')
        read(12,rec=1) m%vp 
        close(12)
        call hud('vp model is read.')

        inquire(file=tmp5//'_vp', exist=alive)
        if(.not.alive) then
            call hud('ERROR: Unable to find vp2 model!')
            stop
        endif
        call alloc(m%vp2,m%nz,m%nx,m%ny)
        open(12,file=tmp5//'_vp',access='direct',recl=n,action='read',status='old')
        read(12,rec=1) m%vp2
        close(12)
        call hud('vp2 model is read.')
        
        !vs
        inquire(file=tmp4//'_vs', exist=alive)
        if(alive) then
            call alloc(m%vs,m%nz,m%nx,m%ny)
            open(12,file=tmp4//'_vs',access='direct',recl=n,action='read',status='old')
            read(12,rec=1) m%vs
            close(12)
            call hud('vs model is read.')
        else
            call alloc(m%vs,1,1,1)
        endif
        
        inquire(file=tmp5//'_vs', exist=alive)
        if(alive) then
            call alloc(m%vs2,m%nz,m%nx,m%ny)
            open(12,file=tmp5//'_vs',access='direct',recl=n,action='read',status='old')
            read(12,rec=1) m%vs2
            close(12)
            call hud('vs 2 model is read.')
        else
            call alloc(m%vs2,1,1,1)
        endif

        !rho
        inquire(file=tmp4//'_rho', exist=alive)
        if(alive) then
            call alloc(m%rho,m%nz,m%nx,m%ny)
            open(12,file=tmp4//'_rho',access='direct',recl=n,action='read',status='old')
            read(12,rec=1) m%rho
            close(12)
            call hud('rho model is read.')
        else
            call alloc(m%rho,1,1,1); m%rho=1000. !in [kg/m3]
        endif
        
        inquire(file=tmp5//'_rho', exist=alive)
        if(alive) then
            call alloc(m%rho2,m%nz,m%nx,m%ny)
            open(12,file=tmp5//'_rho',access='direct',recl=n,action='read',status='old')
            read(12,rec=1) m%rho2
            close(12)
            call hud('rho2 model is read.')
        else
            call alloc(m%rho2,1,1,1); m%rho2=1000. !in [kg/m3]
        endif

        !epsilon
        inquire(file=tmp4//'_eps', exist=alive)
        if(alive) then
            m%is_isotropic=.false.
            call alloc(m%eps,m%nz,m%nx,m%ny)
            open(12,file=tmp4//'_eps',access='direct',recl=n,action='read',status='old')
            read(12,rec=1) m%eps
            close(12)
            call hud('eps model is read.')
        else
            call alloc(m%eps,1,1,1)
        endif
        
        !delta
        inquire(file=tmp4//'_del', exist=alive)
        if(alive) then
            m%is_isotropic=.false.
            call alloc(m%del,m%nz,m%nx,m%ny)
            open(12,file=tmp4//'_del',access='direct',recl=n,action='read',status='old')
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
        if(ask_setup('IF_ISOTROPIC')) then
            m%is_isotropic=get_setup_logical('IF_ISOTROPIC')
        endif
        
        !topography
        call alloc(m%topo,m%nx,m%ny)
        call alloc(m%itopo,m%nx,m%ny,initialize=.false.); m%itopo=1

        if(size(m%vs)>1) then
        if(get_setup_logical('IF_TOPO_FROM_VS',default=.true.)) then
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
        
        inquire(file=tmp4//'_topo', exist=alive)
        if(alive) then
            open(12,file=tmp4//'_topo',access='direct',recl=4*m%nx*m%ny,action='read',status='old')
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

            call alloc( m%vp_mask,  maxval(m%itopo),m%nx,m%ny ); m%vp_mask(:,:,:)  = m%vp(1:maxval(m%itopo),:,:)
        if(size(m%vs)>1) then
            call alloc( m%vs_mask,  maxval(m%itopo),m%nx,m%ny ); m%vs_mask(:,:,:)  = m%vs(1:maxval(m%itopo),:,:)
        endif
            call alloc( m%rho_mask, maxval(m%itopo),m%nx,m%ny ); m%rho_mask(:,:,:) = m%rho(1:maxval(m%itopo),:,:)
        
        !freesurface
        m%if_freesurface=get_setup_logical('IF_FREESURFACE',default=.true.)
        
    end subroutine
    
    subroutine model_reallocate(cmodel)
        character(*) :: cmodel
        
        select case (cmodel)
            case ('rho')
                if(size(m%rho)==1) then
                    deallocate(m%rho,m%rho2)
                    call alloc(m%rho ,m%nz,m%nx,m%ny,initialize=.false.)
                    call alloc(m%rho2,m%nz,m%nx,m%ny,initialize=.false.)
                    m%rho=1000.
                    m%rho2=1000.
                endif
        end select
        
    end subroutine
end
