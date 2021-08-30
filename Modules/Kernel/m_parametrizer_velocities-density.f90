module m_parametrizer
use m_System
use m_Modeling

    !PARAMETERIZATION     -- ALLOWED PARAMETERS
    !velocities-density   -- vp vs ip

    !acoustic:
    !kpa = rho*vp^2 = vp*ip
    !rho0= rho      = ip/vp
    !gvp = gkpa*2rho*vp
    !grho= gkpa*vp^2 + grho0

    !acoustic + gardner:
    !kpa = a*vp^(b+2)
    !rho0= a*vp^b
    !gvp = (gkpa*(b+2)/b*vp^2 + grho0)*ab*vp^(b-1)

    !P-SV:
    !lda = rho(vp^2-2vs^2)
    !mu  = rho*vs^2
    !rho0= rho
    !gvp = glda*2rho*vp
    !gvs = (glda*-2 + gmu)*2rho*vs
    !grho= glda*vp^2 + (-2glda+gmu)*vs^2 + grho0

    !P-SV + gardner:
    !lda = a*vp^(b+2) - 2a*vp^b*vs^2
    !mu  = a*vp^b*vs^2
    !rho0= a*vp^b
    !gvp = (glda*(b+2)/b*vp + (-2glda + gmu)vs^2 + grho0)*ab*vp^(b-1)
    !gvs = (glda*-2 + gmu)*2a*vp^b*vs

    private

    type t_parameter
        character(:),allocatable :: name
        real :: min, max, range
    end type

    type,public :: t_parametrizer
        !info
        character(i_str_xxlen) :: info = &
            'Parameterization: velocities-density'//s_NL// &
            'Allowed pars: vp, vs, ip'//s_NL// &
            'Available empirical laws: Gardner, Castagna (will implemented in future)'

        type(t_parameter),dimension(:),allocatable :: pars
        integer :: npars
        type(t_string),dimension(:),allocatable :: empirical

        logical,dimension(:,:,:,:),allocatable :: is_freeze_zone

        integer :: n1,n2,n3,n
        real :: d1,d2,d3,cell_volume_in_Pa
        
        contains
        procedure :: init
        procedure :: transform
        procedure :: transform_preconditioner
    end type

    type(t_parametrizer),public :: param

    logical :: is_empirical, is_gardner !, is_castagna
    logical :: is_AC, is_EL

    real :: a,b

    contains
    
    subroutine init(self)
        class(t_parametrizer) :: self

        type(t_string),dimension(:),allocatable :: list,sublist

        !read in empirical law
        list=setup%get_strs('EMPIRICAL_LAW')

        is_empirical=size(list)>0

        if(is_empirical) then
            do i=1,size(list)
                if(list(i)%s(1:7)=='Gardner') then
                    !Gardner law rho=a*vp^b
                    !passive rho will be updated according to vp
                    !https://wiki.seg.org/wiki/Dictionary:Gardner%E2%80%99s_equation
                    !http://www.subsurfwiki.org/wiki/Gardner%27s_equation
                    !https://en.wikipedia.org/wiki/Gardner%27s_relation
                    is_gardner=.true.
                    if(len(list(1)%s)<=7) then
                        a=310; b=0.25
                        if(m%ref_rho<1000.) a=0.31
                        
                    else
                        sublist=split(list(1)%s,o_sep=',')
                        a=str2real(sublist(2)%s)
                        b=str2real(sublist(3)%s)
                    endif

                    call hud('Gardner law is enabled: a='//num2str(a)//', b='//num2str(b)//s_NL// &
                        'Parameter rho will not be active in the inversion.')

                elseif(list(i)%s(1:8)=='Castagna') then
                    ! !Castagna mudrock line vp=a*vp+b
                    ! !passive vs will be updated according to vp
                    ! !https://en.wikipedia.org/wiki/Mudrock_line
                    ! if(ia_castagna) then
                    !     c=1/1.16; d=-1.36/1.16
                    ! endif

                endif

            enddo

        endif

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
                self%pars(i)%name='vp'
                self%npars=self%npars+1

            case ('vs' )
                if(is_AC) then
                    call hud('vs parameter from PARAMETER is neglected as the PDE is ACoustic.')
                    cycle loop
                endif
                self%pars(i)%name='vs'
                self%npars=self%npars+1

            case ('rho')
                if(is_empirical) then
                    call hud('rho parameter from PARAMETER is neglected as EMPIRICAL_LAW is read (rho becomes a passive parameter).')
                    cycle loop
                endif
                self%pars(i)%name='rho'
                self%npars=self%npars+1
                
            end select

            self%pars(i)%min=str2real(sublist(2)%s)
            self%pars(i)%max=str2real(sublist(3)%s)
            self%pars(i)%range=self%pars(i)%max-self%pars(i)%min

        enddo loop
        
        deallocate(list,sublist)

        !check vp,vs,rho [min,max] is in the same range as m%ref_vp,ref_vs,ref_rho
        !

        self%n1=m%nz
        self%n2=m%nx
        self%n3=m%ny
        self%n=self%n1*self%n2*self%n3*self%npars

        self%d1=m%dz
        self%d2=m%dx
        self%d3=m%dy
        self%cell_volume_in_Pa=self%d1*self%d2*self%d3*m%ref_kpa

    end subroutine
    
    ! subroutine init_forwardmap(self,fm,oif_g,oif_)
    !     class(t_parametrizer) :: self
    !     type(t_forwardmap) :: fm

    !     call alloc(fm%x, n1,n2,n3,self%npars)
    !     call alloc(fm%g, n1,n2,n3,self%npars)
    !     call alloc(fm%pg,n1,n2,n3,self%npars)
    !     call alloc(fm%d, n1,n2,n3,self%npars)

    ! end subroutine

    subroutine transform(self,o_dir,o_x,o_xprior,o_g)
        class(t_parametrizer) :: self
        character(4),optional :: o_dir
        real,dimension(:,:,:,:),allocatable,optional :: o_x,o_xprior,o_g

        if(present(o_x)) then
            call alloc(o_x,self%n1,self%n2,self%n3,self%npars,oif_protect=.true.)

            if(either(o_dir,'m->x',present(o_dir))=='m->x') then
                do i=1,self%npars
                        select case (self%pars(i)%name)
                        case ('vp' ); o_x(:,:,:,i) = (m%vp -self%pars(i)%min)/self%pars(i)%range
                        case ('vs' ); o_x(:,:,:,i) = (m%vs -self%pars(i)%min)/self%pars(i)%range
                        case ('rho'); o_x(:,:,:,i) = (m%rho-self%pars(i)%min)/self%pars(i)%range
                        end select
                enddo

            else
                do i=1,self%npars
                    select case (self%pars(i)%name)
                    case ('vp' ); m%vp = o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                    case ('vs' ); m%vs = o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                    case ('rho'); m%rho= o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                    end select
                enddo
                ! + gardner
                if(is_gardner) m%rho = a*m%vp**b
                call m%apply_freeze_zone

            endif

        endif

        if(present(o_xprior)) then
            call alloc(o_xprior,self%n1,self%n2,self%n3,self%npars,oif_protect=.true.)
            do i=1,self%npars
                select case (self%pars(i)%name)
                case ('vp' ); o_x(:,:,:,i) = (m%vp_prior -self%pars(i)%min)/self%pars(i)%range
                case ('vs' ); o_x(:,:,:,i) = (m%vs_prior -self%pars(i)%min)/self%pars(i)%range
                case ('rho'); o_x(:,:,:,i) = (m%rho_prior-self%pars(i)%min)/self%pars(i)%range
                end select
            enddo
        endif

        if(present(o_g)) then
            call alloc(o_g,self%n1,self%n2,self%n3,self%npars,oif_protect=.true.)
            call transform_gradient(m%gradient,o_g)
        endif
        
        ! if(present(o_pg)) then
        !     call alloc(o_pg,param%n1,param%n2,param%n3,param%npars,oif_protect=.true.)
        !     call transform_gradient(preco%preco,o_preco)
        ! endif

    end subroutine

    subroutine transform_preconditioner(self,preco_in_m,preco_in_x)
        class(t_parametrizer) :: self
        real,dimension(m%nz,m%nx,m%ny) :: preco_in_m
        real,dimension(self%n1,self%n2,self%n3) :: preco_in_x

        preco_in_x=preco_in_m
    end subroutine

    subroutine transform_gradient(gradient,g)
        real,dimension(m%nz,m%nx,m%ny,ppg%ngrad) :: gradient
        real,dimension(param%n1,param%n2,param%n3,param%npars) :: g
        !gradient(:,:,:,1) = grho
        !gradient(:,:,:,2) = gkpa or glda
        !gradient(:,:,:,3) = gmu

        !acoustic
        if(is_AC .and. .not. is_empirical) then
            do i=1,param%npars
                select case (param%pars(i)%name)
                case ('vp' ); g(:,:,:,i) = gradient(:,:,:,2)*2*m%rho*m%vp
                case ('rho'); g(:,:,:,i) = gradient(:,:,:,2)*m%vp**2 + gradient(:,:,:,1)
                end select
            enddo
        endif

        !acoustic + gardner
        if(is_AC .and. is_gardner) then
            g(:,:,:,1) =(gradient(:,:,:,2)*(b+2)/b*m%vp**2 + gradient(:,:,:,1))*a*b*m%vp**(b-1)
        endif

        !elastic
        if(is_EL .and. .not. is_empirical) then
            do i=1,param%npars
                select case (param%pars(i)%name)
                case ('vp' ); g(:,:,:,i) = gradient(:,:,:,2)*2*m%rho*m%vp
                case ('vs' ); g(:,:,:,i) =(gradient(:,:,:,2)*(-2) + gradient(:,:,:,3))*2*m%rho*m%vs
                case ('rho'); g(:,:,:,i) = gradient(:,:,:,2)*m%vp**2 + (-2*gradient(:,:,:,2)+gradient(:,:,:,3))*m%vs**2 + gradient(:,:,:,1)
                end select
            enddo
        endif

        !elastic + gardner
        if(is_EL .and. is_gardner) then
            do i=1,param%npars
                select case (param%pars(i)%name)
                case ('vp' ); g(:,:,:,i) =(gradient(:,:,:,2)*(b+2)/b*m%vp + (-2*gradient(:,:,:,2)+gradient(:,:,:,3))*m%vs**2 + gradient(:,:,:,1))*a*b*m%vp**(b-1)
                case ('vs' ); g(:,:,:,i) =(gradient(:,:,:,2)*(-2) + gradient(:,:,:,3))*2*a*m%vp**b*m%vs
                end select
            enddo
        endif

        !normaliz g by allowed parameter range
        !s.t. g is in unit [Nm]
        do i=1,param%npars
            g(:,:,:,i)=g(:,:,:,i)*param%pars(i)%range
        enddo
        
        ! where (param%is_freeze_zone) g=0.

    end subroutine

end module
