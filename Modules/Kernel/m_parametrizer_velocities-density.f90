module m_parametrizer
use m_System
use m_Modeling

    !PARAMETERIZATION     -- ALLOWED PARAMETERS
    !velocities-density   -- vp vs rho

    !acoustic:
    !kpa = rho*vp^2 = vp*ip
    !rho0= rho      = ip/vp
    !gvp = gkpa*2rho*vp
    !grho= gkpa*vp^2 + grho0

    !acoustic + gardner:
    !kpa = a*vp^(b+2)
    !rho0= a*vp^b
    !gvp = (gkpa*(b+2)/b*vp^2 + grho0)*b*rho/vp

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
    !gvp = (glda*(b+2)/b*vp^2 + (-2glda + gmu)vs^2 + grho0)*b*rho/vp
    !gvs = (glda*-2 + gmu)*2rho*vs

    !P-SV + castagna:
    !lda = rho( (1-2a^2)vp^2 -4ab*vp -2b^2)
    !mu  = rho( a^2*vp^2 +2ab*vp +b^2)
    !rho0= rho
    !gvp = 2rho( glda*(vp-2a*vs) + gmu*a*vs
    !grho= glda*(vp^2-2vs^2) +gmu*vs^2 + grho0

    !P-SV + poisson:
    !since Vp/Vs=const., gvs is scaled of gvp
    !gvp = glda*2rho*vp

    private

    type t_parameter
        character(:),allocatable :: name
        real :: min, max, range
    end type

    type,public :: t_parametrizer
        !info
        character(i_str_xxlen) :: info = &
            'Parameterization: velocities-density'//s_NL// &
            'Allowed pars: vp, vs, rho'//s_NL// &
            'Available empirical law: Gardner, Castagna, Poisson'

        type(t_parameter),dimension(:),allocatable :: pars
        integer :: npars
        type(t_string),dimension(:),allocatable :: empirical

        integer :: n1,n2,n3,n
        real :: d1,d2,d3
        
        contains
        procedure :: init
        procedure :: transform
        procedure :: transform_preconditioner
    end type

    type(t_parametrizer),public :: param

    logical :: is_empirical=.false., is_gardner=.false., is_castagna=.false., is_poisson=.false.
    logical :: is_AC=.false., is_EL=.false.

    real :: a,b

    contains
    
    subroutine init(self)
        class(t_parametrizer) :: self

        type(t_string),dimension(:),allocatable :: list,sublist

        call hud('Invoked parametrizer module info : '//s_NL//self%info)

        !read in empirical law
        list=setup%get_strs('EMPIRICAL_LAW')

        is_empirical=size(list)>0

        if(is_empirical) then
            do i=1,size(list)
                if(index(list(i)%s,'Gardner')>0) then
                    !Gardner law rho=a*vp^b
                    !passive rho will be updated according to vp
                    !https://wiki.seg.org/wiki/Dictionary:Gardner%E2%80%99s_equation
                    !http://www.subsurfwiki.org/wiki/Gardner%27s_equation
                    !https://en.wikipedia.org/wiki/Gardner%27s_relation
                    is_gardner=.true.
!                     if(len(list(i)%s)<=7) then
                        a=either(0.31,310.,m%rho(1,1,1)<1000.) !g/cm³ or kg/m³
                        b=0.25
                        
!                     else
!                         sublist=split(list(i)%s,o_sep=',') !ifort generates 'catastropic error ... internal compiler error' and lets me report...
!                         a=str2real(sublist(2)%s)
!                         b=str2real(sublist(3)%s)
                        
!                     endif
! 
                    call hud('Gardner law is enabled: a='//num2str(a)//', b='//num2str(b)//s_NL// &
                        'Parameter rho will be passive in the inversion.')

                elseif(list(i)%s(1:8)=='Castagna') then
                    !Castagna mudrock line vs=a*vp+b
                    !passive vs will be updated according to vp
                    !https://en.wikipedia.org/wiki/Mudrock_line
                    is_castagna=.true.
                        a=1/1.16
                        b=-1360./1.16 !m/s

                    call hud('Castagna law is enabled: a='//num2str(a)//', b='//num2str(b)//s_NL// &
                        'Parameter vs will be passive in the inversion.')

                    if(is_gardner) call error("Sorry. Hasn't implemented both Gardner & Castagna's laws")

                elseif(list(i)%s(1:7)=='Poisson') then
                    !Poisson solid vs=a*vp
                    !passive vs will be updated according to vp
                    is_poisson=.true.
                        a=1./sqrt(3.)

                    call hud('Poisson solid is enabled: a='//num2str(a)//s_NL// &
                        'Parameter vs will be passive in the inversion.')

                    if(is_gardner) call error("Sorry. Hasn't implemented both Gardner & Poisson's laws")

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
                    call hud('vs in PARAMETER is neglected as the PDE is ACoustic.')
                    cycle loop
                endif
                if(is_castagna.or.is_poisson) then
                    call hud('vs in PARAMETER is neglected and becomes passive as Castagna/Poisson law is used.')
                    cycle loop
                endif
                self%pars(i)%name='vs'
                self%npars=self%npars+1

            case ('rho')
                if(is_gardner) then
                    call hud('rho in PARAMETER is neglected and becomes passive as Gardner/Poisson law is used')
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

    end subroutine
    
!deprecated
!     subroutine init_forwardmap(self,fm,oif_g,oif_)
!         class(t_parametrizer) :: self
!         type(t_forwardmap) :: fm
! 
!         call alloc(fm%x, n1,n2,n3,self%npars)
!         call alloc(fm%g, n1,n2,n3,self%npars)
!         call alloc(fm%pg,n1,n2,n3,self%npars)
!         call alloc(fm%d, n1,n2,n3,self%npars)
! 
!     end subroutine

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

            else !x->m
                do i=1,self%npars
                    select case (self%pars(i)%name)
                    case ('vp' ); m%vp = o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                    case ('vs' ); m%vs = o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                    case ('rho'); m%rho= o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                    end select
                enddo
                ! + gardner
                if(is_gardner) m%rho = a*m%vp**b
                ! + castagna
                if(is_castagna) m%vs = a*m%vp + b
                ! + poisson
                if(is_poisson) m%vs = a*m%vp
                !
                call m%apply_elastic_continuum
                call m%apply_freeze_zone

            endif

        endif

        if(present(o_xprior)) then
            call alloc(o_xprior,self%n1,self%n2,self%n3,self%npars)

            do i=1,self%npars
                select case (self%pars(i)%name)
                case ('vp' ); o_x(:,:,:,i) = (m%vp_prior -self%pars(i)%min)/self%pars(i)%range
                case ('vs' ); o_x(:,:,:,i) = (m%vs_prior -self%pars(i)%min)/self%pars(i)%range
                case ('rho'); o_x(:,:,:,i) = (m%rho_prior-self%pars(i)%min)/self%pars(i)%range
                end select
            enddo
        endif

        if(present(o_g)) then
            call alloc(o_g,self%n1,self%n2,self%n3,self%npars)
            !correlate_gradient(:,:,:,1) = gkpa
            !correlate_gradient(:,:,:,2) = grho0

            !acoustic
            if(is_AC .and. .not. is_empirical) then
                do i=1,param%npars
                    select case (param%pars(i)%name)
                    case ('vp' ); o_g(:,:,:,i) = correlate_gradient(:,:,:,1)*2*m%rho*m%vp
                    case ('rho'); o_g(:,:,:,i) = correlate_gradient(:,:,:,1)*m%vp**2 + correlate_gradient(:,:,:,2)
                    end select
                enddo
            endif

!            !acoustic + gardner
!            if(is_AC .and. is_gardner) then
!                o_g(:,:,:,1) =(correlate_gradient(:,:,:,2)*(b+2)/b*m%vp**2 + correlate_gradient(:,:,:,1))*b*m%rho/m%vp  !m%rho has tobe updated in prior
!            endif
!
!            !elastic
!            if(is_EL .and. .not. is_empirical) then
!                do i=1,param%npars
!                    select case (param%pars(i)%name)
!                    case ('vp' ); o_g(:,:,:,i) = correlate_gradient(:,:,:,2)*2*m%rho*m%vp
!                    case ('vs' ); o_g(:,:,:,i) =(correlate_gradient(:,:,:,2)*(-2) + correlate_gradient(:,:,:,3))*2*m%rho*m%vs
!                    case ('rho'); o_g(:,:,:,i) = correlate_gradient(:,:,:,2)*m%vp**2 + (-2*correlate_gradient(:,:,:,2)+correlate_gradient(:,:,:,3))*m%vs**2 + correlate_gradient(:,:,:,1)
!                    end select
!                enddo
!            endif
!
!            !elastic + gardner
!            if(is_EL .and. is_gardner) then
!                do i=1,param%npars
!                    select case (param%pars(i)%name)
!                    case ('vp' ); o_g(:,:,:,i) =(correlate_gradient(:,:,:,2)*(b+2)/b*m%vp**2 + (-2*correlate_gradient(:,:,:,2)+correlate_gradient(:,:,:,3))*m%vs**2 + correlate_gradient(:,:,:,1))*b*m%rho/m%vp !m%rho has tobe updated in prior
!                    case ('vs' ); o_g(:,:,:,i) =(correlate_gradient(:,:,:,2)*(-2) + correlate_gradient(:,:,:,3))*2*m%rho*m%vs !m%rho has tobe updated in prior
!                    end select
!                enddo
!            endif
!
!            !elastic + castagna
!            if(is_EL .and. is_castagna) then
!                do i=1,param%npars
!                    select case (param%pars(i)%name)
!                    case ('vp' )
!                        o_g(:,:,:,i) =2*m%rho*( correlate_gradient(:,:,:,2)*(m%vp-2*a*m%vs) + correlate_gradient(:,:,:,3)*a*m%vs )
!                    case ('rho')
!                        o_g(:,:,:,i) =correlate_gradient(:,:,:,2)*(m%vp**2-2*m%vs**2) + correlate_gradient(:,:,:,3)*m%vs**2 + correlate_gradient(:,:,:,1)
!                    end select
!                enddo
!            endif

            !normaliz g by allowed parameter range
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
