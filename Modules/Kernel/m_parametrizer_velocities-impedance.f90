module m_parametrizer
use m_System
use m_Modeling
use m_empirical

    !PARAMETERIZATION     -- ALLOWED PARAMETERS
    !velocities-impedance -- vp vs ip

    !acoustic:
    !kpa = rho*vp^2 = vp*ip
    !rho0= rho      = ip/vp
    !gvp = (gkpa*vp - grho0/vp)*rho
    !gip =  gkpa*vp + grho0/vp

    !P-SV:
    !lda = rho(vp^2-2vs^2) = vp*ip - 2vs^2*ip/vp
    !mu  = rho*vs^2        = vs^2*ip/vp
    !rho0= rho             = ip/vp
    !gvp = (glda*vp^2 + (2glda-gmu)vs^2 - grho0)*rho/vp
    !gvs = (-2glda + gmu)*2vs*rho
    !gip = (glda*vp^2 + (-2glda+gmu)*vs^2 +grho0) /vp

    private

    type t_parameter
        character(:),allocatable :: name
        real :: min, max, range
    end type

    type,public :: t_parametrizer
        !info
        character(i_str_xxlen) :: info = &
            'Parameterization: velocities-impedance'//s_NL// &
            'Allowed pars: vp, vs, ip'

        type(t_parameter),dimension(:),allocatable :: pars
        integer :: npars

        integer :: n1,n2,n3,n
        real :: d1,d2,d3
        
        contains
        procedure :: init
        procedure :: transform
        procedure :: transform_preconditioner
    end type

    type(t_parametrizer),public :: param

    logical :: is_AC=.false., is_EL=.false.
    integer :: i_vp=0, i_vs=0, i_ip=0

    contains
    
    subroutine init(self)
        class(t_parametrizer) :: self

        type(t_string),dimension(:),allocatable :: list,sublist

        call hud('Invoked parametrizer module info : '//s_NL//self%info)

        if(allocated(   list)) deallocate(   list)
        if(allocated(sublist)) deallocate(sublist)
        
        
        !PDE info
        is_AC = index(ppg%info,'AC')>0
        is_EL = index(ppg%info,'EL')>0

        !read in active parameters and their allowed ranges
        list=setup%get_strs('PARAMETER',o_default='vp:1500:3400')
        
        self%npars=size(list)
        allocate(self%pars(self%npars))

        !re-count to remove illegal parameters
        self%npars=0
        loop: do i=1,size(list)
            sublist=split(list(i)%s,o_sep=':') !=[name, min, max]

            select case (sublist(1)%s)
            case ('vp' )
                i_vp=i
                self%pars(i)%name='vp'
                self%npars=self%npars+1

            case ('vs' )
                if(is_AC) then
                    call hud('vs in PARAMETER is neglected as the PDE is ACoustic.')
                    cycle loop
                endif
                i_vs=i
                self%pars(i)%name='vs'
                self%npars=self%npars+1

            case ('ip')
                i_ip=i
                self%pars(i)%name='ip'
                self%npars=self%npars+1
                
            end select

            self%pars(i)%min=str2real(sublist(2)%s)
            self%pars(i)%max=str2real(sublist(3)%s)
            self%pars(i)%range=self%pars(i)%max-self%pars(i)%min

        enddo loop
        
        deallocate(list,sublist)

        !check vp,vs,rho [min,max] is in the same range as m%ref_vp,ref_vs,ref_rho
        !

        call empirical_init

        self%n1=m%nz
        self%n2=m%nx
        self%n3=m%ny
        self%n=self%n1*self%n2*self%n3*self%npars

        self%d1=m%dz
        self%d2=m%dx
        self%d3=m%dy

    end subroutine
    

    subroutine transform(self,o_dir,o_x,o_xprior,o_g)
        class(t_parametrizer) :: self
        character(4),optional :: o_dir
        real,dimension(:,:,:,:),allocatable,optional :: o_x,o_xprior,o_g

        real,dimension(:,:,:),allocatable :: tmp_vp

        if(present(o_x)) then
            call alloc(o_x,self%n1,self%n2,self%n3,self%npars,oif_protect=.true.)

            if(either(o_dir,'m->x',present(o_dir))=='m->x') then
                if(i_vp >0) o_x(:,:,:,i_vp ) = (m%vp      -self%pars(i_vp)%min)/self%pars(i_vp)%range 
                if(i_vs >0) o_x(:,:,:,i_vs ) = (m%vs      -self%pars(i_vs)%min)/self%pars(i_vs)%range 
                if(i_ip >0) o_x(:,:,:,i_ip ) = (m%vp*m%rho-self%pars(i_ip)%min)/self%pars(i_ip)%range 

                call empirical_m2x('velocities-impedance')

            else !x->m

                if(i_vp >0) then
                    tmp_vp = o_x(:,:,:,i_vp)*self%pars(i_vp)%range +self%pars(i_vp)%min  !implicit allocation
                    m%rho = m%vp*m%rho / tmp_vp
                    m%rho0= m%vp*m%rho0/ tmp_vp !rho0 in m%rho0 has diff meaning from grho0
                    m%vp  = tmp_vp
                    deallocate(tmp_vp)
                endif
                if(i_vs >0) m%vs  =  o_x(:,:,:,i_vs)*self%pars(i_vs)%range +self%pars(i_vs)%min
                if(i_ip >0) m%rho = (o_x(:,:,:,i_ip)*self%pars(i_ip)%range +self%pars(i_ip)%min)/m%vp

                call empirical_x2m('velocities-impedance')
                
                call m%apply_elastic_continuum
                call m%apply_freeze_zone

            endif

        endif

        if(present(o_xprior)) then
            call alloc(o_xprior,self%n1,self%n2,self%n3,self%npars)

            if(i_vp >0) o_x(:,:,:,i_vp ) = (m%vp_prior            -self%pars(i_vp)%min)/self%pars(i_vp)%range
            if(i_vs >0) o_x(:,:,:,i_vs ) = (m%vs_prior            -self%pars(i_vs)%min)/self%pars(i_vs)%range
            if(i_ip >0) o_x(:,:,:,i_ip ) = (m%vp_prior*m%rho_prior-self%pars(i_ip)%min)/self%pars(i_ip)%range

            call empirical_m2x('velocities-impedance')

        endif

        if(present(o_g)) then
            call alloc(o_g,self%n1,self%n2,self%n3,self%npars)

            if(is_AC) then
                !correlate_gradient(:,:,:,1) = grho0
                !correlate_gradient(:,:,:,2) = gkpa
                if(i_vp >0) o_g(:,:,:,i_vp ) = (correlate_gradient(:,:,:,2)*m%vp - correlate_gradient(:,:,:,1)/m%vp)*m%rho
                if(i_ip >0) o_g(:,:,:,i_ip ) =  correlate_gradient(:,:,:,2)*m%vp + correlate_gradient(:,:,:,1)/m%vp

                call empirical_gradient('velocities-impedance',o_gvp=o_g(:,:,:,i_vp),o_gip=o_g(:,:,:,i_ip))
            endif

            if(is_EL) then
                !correlate_gradient(:,:,:,1) = grho0
                !correlate_gradient(:,:,:,2) = glda
                !correlate_gradient(:,:,:,2) = gmu
                if(i_vp >0) o_g(:,:,:,i_vp ) =(correlate_gradient(:,:,:,1)*m%vp**2 + (2*correlate_gradient(:,:,:,1)-correlate_gradient(:,:,:,2))*m%vs**2 - correlate_gradient(:,:,:,3))*m%rho/m%vp
                if(i_vs >0) o_g(:,:,:,i_vs ) =(-2*correlate_gradient(:,:,:,1) + correlate_gradient(:,:,:,2))*2*m%rho*m%vs
                if(i_ip >0) o_g(:,:,:,i_ip ) =(correlate_gradient(:,:,:,1)*m%vp**2 + (-2*correlate_gradient(:,:,:,1)+correlate_gradient(:,:,:,2))*m%vs**2 + correlate_gradient(:,:,:,3))/m%vp

                call empirical_gradient('velocities-impedance',o_gvp=o_g(:,:,:,i_vp),o_gvs=o_g(:,:,:,i_vs),o_gip=o_g(:,:,:,i_ip))
            endif

            !normaliz g by allowed parameter range
            !s.t. g is in unit [Nm]
            do i=1,param%npars
                o_g(:,:,:,i)=o_g(:,:,:,i)*param%pars(i)%range
            enddo
            
        endif

    end subroutine

    subroutine transform_preconditioner(self,preco_in_m,preco_in_x)
        class(t_parametrizer) :: self
        real,dimension(m%nz,m%nx,m%ny) :: preco_in_m
        real,dimension(self%n1,self%n2,self%n3) :: preco_in_x

        preco_in_x=preco_in_m
    end subroutine

end module
