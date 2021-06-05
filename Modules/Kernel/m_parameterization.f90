module m_parameterization
use m_sysio
use m_arrayop
use m_model
use m_field, only:waveeq_info
use m_gradient, only: gradient


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

    type t_parameter
        character(:),allocatable :: name
        real :: min, max, range
    end type

    type t_paramerization

        character(:),allocatable :: parameterization
        type(t_parameter),dimension(:),allocatable :: pars
        integer :: npar
        type(t_string),dimension(:),allocatable :: empirical        
        
        contains
        procedure :: init => init
    end type

    logical :: is_ac, is_el
    logical :: is_empirical,is_gardner,is_castagna

    contains
    
    subroutine init
        type(t_string),dimension(:),allocatable :: list,tmp

        !read in empirical law
        is_empirical=setup%check('EMPIRICAL_LAW')

        if(is_empirical) then
            if(empirical(1)%s=='Gardner') then
                !Gardner law rho=a*vp^b
                !passive rho will be updated according to vp
                !https://wiki.seg.org/wiki/Dictionary:Gardner%E2%80%99s_equation
                !http://www.subsurfwiki.org/wiki/Gardner%27s_equation
                !https://en.wikipedia.org/wiki/Gardner%27s_relation
                is_gardner=.true.
                if(size(empirical)>1) then
                    a=str2real(empirical(2)%s)
                    b=str2real(empirical(3)%s)
                else
                    a=310; b=0.25
                    if(m%ref_rho<1000.) a=0.31
                endif

                call hud('Gardner law is enabled. Parameter rho will not be active in the inversion.')

            elseif(empirical(1)%s=='Castagna') then
                ! !Castagna mudrock line vp=a*vp+b
                ! !passive vs will be updated according to vp
                ! !https://en.wikipedia.org/wiki/Mudrock_line
                ! if(if_castagna) then
                !     c=1/1.16; d=-1.36/1.16
                ! endif

            endif

        endif
        
        !pde info
        is_ac = index(field%pde_info,'acoustic')>0
        is_el = index(field%pde_info,'elastic')>0

        !read in active parameters and their allowed ranges
        list=setup%get_strs('ACTIVE_PARAMETER',default='vp:1500:3400')

        self%npar=size(list)
        call alloc(self%pars(npar))

        do i=1,npar
            tmp=split(list(i),o_sep=':') !=[name, min, max]

            select case (tmp(1)%s)
            case ('vp' )
                self%pars(i)%name='vp'
            case ('vs' )
                if(is_el) then
                    self%pars(i)%name='vs'
                else
                    call hud('vs parameter from ACTIVE_PARAMETER is neglected since the PDE is acoustic')
                endif
            case ('rho')
                if(.not. is_empirical) then
                    self%pars(i)%name='rho'
                else
                    call hud('rho parameter from ACTIVE_PARAMETER is neglected since EMPIRICAL_LAW is read')
                endif
            end select

            self%pars(i)%min=str2real(tmp(2)%s)
            self%pars(i)%max=str2real(tmp(3)%s)
            self%pars(i)%range=self%pars(i)%max-self%pars(i)%min

        enddo

        !check vp,vs,rho [min,max] is in the same range as m%ref_vp,ref_vs,ref_rho
        !

    end subroutine
        
    subroutine transform(dir,x,g)
        character(3) :: dir
        real,dimension(m%nz,m%nx,m%ny,npar) :: x
        real,dimension(m%nz,m%nx,m%ny,npar),optional :: g

        !model
        if(dir=='m2x') then
            do i=1,npar
                select case (self%pars(i)%name)
                case ('vp' ); x(:,:,:,i) = (m%vp -self%pars(i)%min))/self%pars(i)%range
                case ('vs' ); x(:,:,:,i) = (m%vs -self%pars(i)%min))/self%pars(i)%range
                case ('rho'); x(:,:,:,i) = (m%rho-self%pars(i)%min))/self%pars(i)%range
                end select
            enddo

        else !x2m
            do i=1,npar
                select case (self%pars(i)%name)
                case ('vp' ); m%vp = x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                case ('vs' ); m%vs = x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                case ('rho'); m%rho= x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                end select
            enddo

            ! + gardner
            if(if_gardner) m%rho = a*m%vp**b

            call m%apply_freeze_zone

        endif

        !gradient
        !!for units of gradient and g, see m_field*.f90
        if(present(g)) then
        if(dir=='m2x') then
            
            !acoustic
            if(is_ac .and. .not. is_empirical) then
                do i=1,self%npar
                    select case (self%pars(i))
                    case ('vp' ); g(:,:,:,ipar) = fobj%gradient(:,:,:,1)*2*m%rho*m%vp
                    case ('rho'); g(:,:,:,ipar) = fobj%gradient(:,:,:,1)*m%vp**2 + fobj%gradient(:,:,:,2)
                    end select
                enddo
            endif

            !acoustic + gardner
            if(is_ac .and. is_gardner) then
                g(:,:,:,1) =(fobj%gradient(:,:,:,1)*(b+2)/b*m%vp**2 + fobj%gradient(:,:,:,2))*a*b*m%vp**(b-1)
            endif

            !elastic
            if(is_el .and. .not. is_empirical) then
                do i=1,self%npar
                    select case (self%pars(i))
                    case ('vp' ); g(:,:,:,ipar) = fobj%gradient(:,:,:,1)*2*m%rho*m%vp
                    case ('vs' ); g(:,:,:,ipar) =(fobj%gradient(:,:,:,1)*(-2) + fobj%gradient(:,:,:,2))*2*m%rho*m%vs
                    case ('rho'); g(:,:,:,ipar) = fobj%gradient(:,:,:,1)*m%vp**2 + (-2*fobj%gradient(:,:,:,1)+fobj%gradient(:,:,:,2))*m%vs**2 + fobj%gradient(:,:,:,3)
                    end select
                enddo
            endif

            !elastic + gardner
            if(is_el .and. is_gardner) then
                do i=1,self%npar
                    select case (self%pars(i))
                    case ('vp' ); g(:,:,:,i) =(fobj%gradient(:,:,:,1)*(b+2)/b*m%vp + (-2*fobj%gradient(:,:,:,1)+fobj%gradient(:,:,:,2))*m%vs**2 + fobj%gradient(:,:,:,3))*a*b*m%vp**(b-1)
                    case ('vs' ); g(:,:,:,i) =(fobj%gradient(:,:,:,1)*(-2) + fobj%gradient(:,:,:,2))*2*a*m%vp**b*m%vs
                    end select
                enddo
            endif


            !normaliz g to be unitless
            do i=1,self%npar
                g(:,:,:,i)=g(:,:,:,i)*self%pars(i)%range
            enddo

        endif
        endif
        
    end subroutine

end module
