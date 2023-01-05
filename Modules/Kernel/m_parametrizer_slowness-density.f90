module m_parametrizer
use m_System
use m_Modeling

    !PARAMETERIZATION   -- ALLOWED PARAMETERS
    !slowness-density   -- sp spr rho
    !sp: P-wave slowness, spr: vs/vp ratio = vs*sp

    !acoustic:
    !kpa = rho*vp^2 = rho*sp^-2
    !rho0= rho
    !gsp = gkpa*(-2)rho*vp^3
    !grho= gkpa*vp^2 + grho0

    !acoustic + gardner:
    !kpa = a*vp^(b+2) = a*sp^(-b-2)
    !rho0= a*vp^b     = a*sp^-b
    !gsp = (gkpa*(b+2)/b + grho0*vp^-2)*(-ab)*vp^(b+3)

    !P-SV:
    !lda = rho(vp^2-2vs^2) = rho*sp^-2*(1-2spr^2)
    !mu  = rho*vs^2        = rho*(spr/sp)^2
    !rho0= rho
    !gsp = (glda*vp^3 + (-2glda+gmu)vp*vs^2)*(-2)rho
    !gspr= (-2glda + gmu)*2rho*vp*vs
    !grho= glda*vp^2 + (-2glda+gmu)*vp*vs + grho0

    !P-SV + gardner:
    !lda = a*vp^(b+2) - 2a*vp^b*vs^2 = a*sp^(-b-2)*(1-2spr^2)
    !mu  = a*vp^b*vs^2               = a*sp^(-b-2)*spr^2
    !rho0= a*vp^b                    = a*sp^-b
    !gsp = (glda*vp^2 + (-2glda + gmu)vs^2 + grho0*b(b+2))*-a/(b+2)*vp^(b+1)
    !gspr= (-2glda + gmu)*2a*vp^(b+1)*vs

    private 
    
    type t_parameter
        character(:),allocatable :: name
        real :: min, max, range
    end type

    type,public :: t_parametrizer
        !info
        character(i_str_xxlen) :: info = &
            'Parameterization: slowness-density'//s_NL// &
            'Allowed pars: sp, spr, rho'//s_NL// &
            'Available empirical laws: Gardner'

        type(t_parameter),dimension(:),allocatable :: pars
        integer :: npars
        type(t_string),dimension(:),allocatable :: empirical

        logical,dimension(:,:,:,:),allocatable :: is_freeze_zone

        integer :: n1,n2,n3,n
        real :: d1,d2,d3
        
        contains
        procedure :: init
        procedure :: transform
        procedure :: transform_preconditioner
    end type

    type(t_parametrizer),public :: param

    logical :: is_empirical=.false., is_gardner=.false.
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
                        a=either(310.,0.31,m%ref_rho<1000.)
                        b=0.25
                        
!                     else
!                         sublist=split(list(i)%s,o_sep=',') !ifort generates 'catastropic error ... internal compiler error' and lets me report...
!                         a=str2real(sublist(2)%s)
!                         b=str2real(sublist(3)%s)
                        
!                     endif
! 
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
        list=setup%get_strs('PARAMETER',o_default='sp:0.3:0.67')
        
        self%npars=size(list)
        allocate(self%pars(self%npars))

        !re-count to remove illegal parameters
        self%npars=0
        loop: do i=1,size(list)
            sublist=split(list(i)%s,o_sep=':') !=[name, min, max]

            select case (sublist(1)%s)
            case ('sp' )
                self%pars(i)%name='sp'
                self%npars=self%npars+1

            case ('spr')
                if(is_AC) then
                    call hud('spr in PARAMETER is neglected as the PDE is ACoustic.')
                    cycle loop
                endif
                self%pars(i)%name='spr'
                self%npars=self%npars+1

            case ('rho')
                if(is_empirical) then
                    call hud('rho in PARAMETER is neglected as EMPIRICAL_LAW is read (ip becomes a passive parameter).')
                    cycle loop
                endif
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
                do i=1,self%npars
                    select case (self%pars(i)%name)
                    case ('sp' ); o_x(:,:,:,i) = (  1./m%vp-self%pars(i)%min)/self%pars(i)%range
                    case ('spr'); o_x(:,:,:,i) = (m%vs/m%vp-self%pars(i)%min)/self%pars(i)%range
                    case ('rho'); o_x(:,:,:,i) = (m%rho    -self%pars(i)%min)/self%pars(i)%range
                    end select
                enddo

            else !x->m
                do i=1,self%npars
                    select case (self%pars(i)%name)
                    case ('sp')
                        tmp_vp = 1./(o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min)   !implicit allocation
                        if(allocated(m%vs)) m%vs = m%vs/m%vp * tmp_vp
                        m%vp  = tmp_vp
                        deallocate(tmp_vp)
                    case ('spr'); m%vs =  (o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min)*m%vp
                    case ('rho'); m%rho = o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                    end select
                enddo

                ! + gardner
                if(is_gardner) m%rho = a*m%vp**b
                call m%apply_freeze_zone

            endif

        endif

        if(present(o_xprior)) then
            call alloc(o_xprior,self%n1,self%n2,self%n3,self%npars)

            do i=1,self%npars
                select case (self%pars(i)%name)
                case ('sp' ); o_x(:,:,:,i) = (1./m%vp_prior        -self%pars(i)%min)/self%pars(i)%range
                case ('spr'); o_x(:,:,:,i) = (m%vs_prior/m%vp_prior-self%pars(i)%min)/self%pars(i)%range
                case ('rho'); o_x(:,:,:,i) = (m%rho_prior          -self%pars(i)%min)/self%pars(i)%range
                end select
            enddo
        endif

        if(present(o_g)) then
            call alloc(o_g,self%n1,self%n2,self%n3,self%npars)
            !m%gradient(:,:,:,1) = grho
            !m%gradient(:,:,:,2) = gkpa or glda
            !m%gradient(:,:,:,3) = gmu

            !acoustic
            if(is_AC .and. .not. is_empirical) then
                do i=1,param%npars
                    select case (param%pars(i)%name)
                    case ('sp' ); o_g(:,:,:,i) = m%gradient(:,:,:,2)*(-2)*m%rho*m%vp**3
                    case ('rho'); o_g(:,:,:,i) = m%gradient(:,:,:,2)*m%vp**2 + m%gradient(:,:,:,1)
                    end select
                enddo
            endif

            !acoustic + gardner
            if(is_AC .and. is_gardner) then
                o_g(:,:,:,1) =(m%gradient(:,:,:,2)*(b+2)/b + m%gradient(:,:,:,1)/m%vp**2)*(-a)*b*m%vp**(b+3)
            endif

            !elastic
            if(is_EL .and. .not. is_empirical) then
                do i=1,param%npars
                    select case (param%pars(i)%name)
                    case ('sp' ); o_g(:,:,:,i) =(m%gradient(:,:,:,2)*m%vp**3 + (-2*m%gradient(:,:,:,2)+m%gradient(:,:,:,3)*m%vp*m%vs**2))*(-2)*m%rho
                    case ('spr'); o_g(:,:,:,i) =(-2*m%gradient(:,:,:,2) + m%gradient(:,:,:,3))*2*m%rho*m%vp*m%vs
                    case ('rho'); o_g(:,:,:,i) = m%gradient(:,:,:,2)*m%vp**2 + (-2*m%gradient(:,:,:,2)+m%gradient(:,:,:,3))*m%vp*m%vs + m%gradient(:,:,:,1)
                    end select
                enddo
            endif

            !elastic + gardner
            if(is_EL .and. .not. is_gardner) then
                do i=1,param%npars
                    select case (param%pars(i)%name)
                    case ('sp' ); o_g(:,:,:,i) =(m%gradient(:,:,:,2)*m%vp**2 + (-2*m%gradient(:,:,:,2)+m%gradient(:,:,:,3))*m%vs**2 + m%gradient(:,:,:,1)*b*(b+2))*(-a)/(b+2)*m%vp**(b+1)
                    case ('spr'); o_g(:,:,:,i) =(-2*m%gradient(:,:,:,2) + m%gradient(:,:,:,3))*2*a*m%vp**(b+1)*m%vs
                    end select
                enddo
            endif

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
