module m_parametrizer
use m_System
use m_Modeling
use m_empirical

    private

    type t_parameter
        character(:),allocatable :: name
        real :: min, max, range
    end type

    type,public :: t_parametrizer
        !info
        character(i_str_xxlen) :: info = &
            'Parameterization: velocities-density'//s_NL// &
            'Allowed pars: vp, vs, rho'

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

    logical :: is_grho,is_gbuo,is_gkpa,is_gikpa,is_glda,is_gmu, is_gqp
    integer :: i_vp=0, i_vs=0, i_rho=0, i_qp=0

    contains
    
    subroutine init(self)
        class(t_parametrizer) :: self

        type(t_string),dimension(:),allocatable :: list,sublist

        call hud('Invoked parametrizer module info : '//s_NL//self%info)

        if(allocated(   list)) deallocate(   list)
        if(allocated(sublist)) deallocate(sublist)
        
        
        !check basic gradients provided from propagator
        is_grho = index(ppg%info,'grho')>0
        is_gbuo = index(ppg%info,'gbuo')>0
        is_gkpa = index(ppg%info,'gkpa')>0 
        is_gikpa= index(ppg%info,'gikpa')>0
        is_glda = index(ppg%info,'glda')>0
        is_gmu  = index(ppg%info,'gmu')>0
        is_gqp  = index(ppg%info,'gqp')>0

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
                if(index(ppg%info,'AC')>0) then
                    call hud('vs in PARAMETER is neglected as the PDE is ACoustic.')
                    cycle loop
                endif
                i_vs=i
                self%pars(i)%name='vs'
                self%npars=self%npars+1

            case ('rho')
                i_rho=i
                self%pars(i)%name='rho'
                self%npars=self%npars+1

            case ('qp')
                i_qp=i
                self%pars(i)%name='qp'
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

        if(present(o_x)) then
            call alloc(o_x,self%n1,self%n2,self%n3,self%npars,oif_protect=.true.)

            if(either(o_dir,'m->x',present(o_dir))=='m->x') then
                if(i_vp >0) o_x(:,:,:,i_vp ) = (m%vp -self%pars(i_vp )%min)/self%pars(i_vp )%range
                if(i_vs >0) o_x(:,:,:,i_vs ) = (m%vs -self%pars(i_vs )%min)/self%pars(i_vs )%range
                if(i_rho>0) o_x(:,:,:,i_rho) = (m%rho-self%pars(i_rho)%min)/self%pars(i_rho)%range
                if(i_qp>0)  o_x(:,:,:,i_qp)  = (m%qp - self%pars(i_qp)%min)/self%pars(i_qp)%range

                call empirical_m2x('velocities-density')

            else !x->m
                if(i_vp >0) m%vp = o_x(:,:,:,i_vp )*self%pars(i_vp )%range +self%pars(i_vp )%min
                if(i_vs >0) m%vs = o_x(:,:,:,i_vs )*self%pars(i_vs )%range +self%pars(i_vs )%min
                if(i_rho>0) m%rho= o_x(:,:,:,i_rho)*self%pars(i_rho)%range +self%pars(i_rho)%min
                if(i_qp>0)  m%qp = o_x(:,:,:,i_qp) *self%pars(i_qp)%range  +self%pars(i_qp)%min

                call empirical_x2m('velocities-density')
                
                if(index(ppg%info,'EL')>0) call m%apply_elastic_continuum
                call m%apply_freeze_zone

            endif

        endif

        if(present(o_xprior)) then
            call alloc(o_xprior,self%n1,self%n2,self%n3,self%npars)

            if(i_vp >0) o_x(:,:,:,i_vp ) = (m%vp_prior -self%pars(i_vp )%min)/self%pars(i_vp )%range
            if(i_vs >0) o_x(:,:,:,i_vs ) = (m%vs_prior -self%pars(i_vs )%min)/self%pars(i_vs )%range
            if(i_rho>0) o_x(:,:,:,i_rho) = (m%rho_prior-self%pars(i_rho)%min)/self%pars(i_rho)%range
            ! if(i_qp>0)  o_x(:,:,:,i_qp)  = (m%qp_prior-self%pars(i_qp)%min)/self%pars(i_qp)%range
                
            call empirical_m2x('velocities-density')

        endif

        if(present(o_g)) then
            call alloc(o_g,self%n1,self%n2,self%n3,self%npars)

            n_entry=0

            if(is_grho.and.is_gkpa) then
                call hud('Parametrizer finds grho & gkpa')
                !correlate_gradient(:,:,:,1) = grho
                !correlate_gradient(:,:,:,2) = gkpa
                !
                !kpa = rho*vp² = vp*ip
                !rho0= rho      = ip/vp
                !So,
                !gvp = gkpa*2*rho*vp
                !grho= gkpa*vp² + grho0
                if(i_vp >0) o_g(:,:,:,i_vp ) = correlate_gradient(:,:,:,2)*2*m%rho*m%vp
                if(i_rho>0) o_g(:,:,:,i_rho) = correlate_gradient(:,:,:,2)*m%vp**2 + correlate_gradient(:,:,:,1)
                ! if(i_qp>0)  o_g(:,:,:,i_qp)  = correlate_gradient(:,:,:,3)

                call empirical_gradient('velocities-density',o_gvp=o_g(:,:,:,i_vp),o_grho=o_g(:,:,:,i_rho))

                n_entry=n_entry+1
            endif

            if(is_gbuo.and.is_gikpa.and.is_gqp) then 
                call hud('Parametrizer finds gbuo & gikpa & gqp')    
                !correlate_gradient(:,:,:,1) = gbuo
                !correlate_gradient(:,:,:,2) = gikpa
                !
                !ikpa= kpa⁻¹ = rho⁻¹ vp⁻²
                !buo = rho⁻¹
                !So,
                !gvp = gikpa* rho⁻¹*(-2)vp⁻³
                !grho= -rho⁻²*( gbuo + gikpa*vp⁻² )
                if(i_vp >0) o_g(:,:,:,i_vp ) = correlate_gradient(:,:,:,2)*(-2.)/m%rho/(m%vp**3)
                if(i_rho>0) o_g(:,:,:,i_rho) = -m%rho**(-2)*( &
                    correlate_gradient(:,:,:,1) + correlate_gradient(:,:,:,2)/(m%vp**2) )
                ! if(i_qp>0)  o_g(:,:,:,i_qp)  = correlate_gradient(:,:,:,3)

                call empirical_gradient('velocities-density',o_gvp=o_g(:,:,:,i_vp),o_grho=o_g(:,:,:,i_rho))

                n_entry=n_entry+1
            endif

            if(is_grho.and.is_glda.and.is_gmu) then
                call hud('Parametrizer finds grho glda & gmu')
                !correlate_gradient(:,:,:,1) = grho
                !correlate_gradient(:,:,:,2) = glda
                !correlate_gradient(:,:,:,2) = gmu 
                !
                !lda = rho(vp²-2vs²)
                !mu  = rho*vs²
                !rho0= rho
                !So,
                !gvp = glda*2rho*vp
                !gvs = (glda*-2 + gmu)*2rho*vs
                !grho= glda*vp² + (-2glda+gmu)*vs² + grho0

                if(i_vp >0) o_g(:,:,:,i_vp ) = correlate_gradient(:,:,:,2)*2*m%rho*m%vp
                if(i_vs >0) o_g(:,:,:,i_vs ) =(correlate_gradient(:,:,:,2)*(-2) + correlate_gradient(:,:,:,3))*2*m%rho*m%vs
                if(i_rho>0) o_g(:,:,:,i_rho) = correlate_gradient(:,:,:,2)*m%vp**2 + (-2*correlate_gradient(:,:,:,2)+correlate_gradient(:,:,:,3))*m%vs**2 + correlate_gradient(:,:,:,1)

                call empirical_gradient('velocities-density',o_gvp=o_g(:,:,:,i_vp),o_gvs=o_g(:,:,:,i_vs),o_grho=o_g(:,:,:,i_rho))

                n_entry=n_entry+1
            endif

            if(n_entry/=1) then
                call error('Parametrizer has n_entry='//num2str(n_entry))
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