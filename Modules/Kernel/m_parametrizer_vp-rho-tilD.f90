module m_parametrizer
use m_System
use m_Modeling

    !PARAMETERIZATION         -- ALLOWED PARAMETERS
    !velocities-rho-tilde D   -- vp rho tilD tilrho

    !acoustic:
    !ikpa= 1/kpa = 1/rho/vp^2
    !buo = 1/rho
    !gvp = gikpa*(-2/rho/vp^3)
    !grho= -1/rho^2*( gikpa/vp^2 + gbuo )


    private

    type t_parameter
        character(:),allocatable :: name
        real :: min, max, range
    end type

    type,public :: t_parametrizer
        !info
        character(i_str_xxlen) :: info = &
            'Parameterization: vp-rho-tilD'//s_NL// &
            'Allowed pars: vp, rho, tilD, tilrho'!//s_NL// &
            !'Available empirical law: Gardner, Castagna'

        type(t_parameter),dimension(:),allocatable :: pars
        integer :: npars
        !type(t_string),dimension(:),allocatable :: empirical

        integer :: n1,n2,n3,n
        real :: d1,d2,d3
        
        contains
        procedure :: init
        procedure :: transform
        procedure :: transform_preconditioner
    end type

    type(t_parametrizer),public :: param

    ! logical :: is_empirical=.false., is_gardner=.false., is_castagna=.false.
    ! logical :: is_AC=.false., is_EL=.false.

    ! real :: a,b

    contains
    
    subroutine init(self)
        class(t_parametrizer) :: self

        type(t_string),dimension(:),allocatable :: list,sublist

        call hud('Invoked parametrizer module info : '//s_NL//self%info)

        !!read in empirical law
        !list=setup%get_strs('EMPIRICAL_LAW')

!         is_empirical=size(list)>0

!         if(is_empirical) then
!             do i=1,size(list)
!                 if(index(list(i)%s,'Gardner')>0) then
!                     !Gardner law rho=a*vp^b
!                     !passive rho will be updated according to vp
!                     !https://wiki.seg.org/wiki/Dictionary:Gardner%E2%80%99s_equation
!                     !http://www.subsurfwiki.org/wiki/Gardner%27s_equation
!                     !https://en.wikipedia.org/wiki/Gardner%27s_relation
!                     is_gardner=.true.
! !                     if(len(list(i)%s)<=7) then
!                         a=either(0.31,310.,m%rho(1,1,1)<1000.) !g/cm³ or kg/m³
!                         b=0.25
                        
! !                     else
! !                         sublist=split(list(i)%s,o_sep=',') !ifort generates 'catastropic error ... internal compiler error' and lets me report...
! !                         a=str2real(sublist(2)%s)
! !                         b=str2real(sublist(3)%s)
                        
! !                     endif
! ! 
!                     call hud('Gardner law is enabled: a='//num2str(a)//', b='//num2str(b)//s_NL// &
!                         'Parameter rho will passive in the inversion.')

!                 elseif(list(i)%s(1:8)=='Castagna') then
!                     !Castagna mudrock line vs=a*vs+b
!                     !passive vs will be updated according to vp
!                     !https://en.wikipedia.org/wiki/Mudrock_line
!                     is_castagna=.true.
!                         a=1/1.16
!                         b=-1360./1.16 !m/s

!                     call hud('Castagna law is enabled: a='//num2str(a)//', b='//num2str(b)//s_NL// &
!                         'Parameter vs will be passive in the inversion.')

!                     if(is_gardner) then
!                         call error("Sorry. Hasn't implemented both Gardner & Castagna's laws")
!                     endif

!                 endif

!             enddo

!         endif

        if(allocated(   list)) deallocate(   list)
        if(allocated(sublist)) deallocate(sublist)
        
        
        ! !PDE info
        ! is_AC = index(ppg%info,'AC')>0
        ! is_EL = index(ppg%info,'EL')>0

        !read in active parameters and their allowed ranges
        ! list=setup%get_strs('PARAMETER',o_default='vp2:2250000:11560000')
        list=setup%get_strs('PARAMETER',o_default='vp:1500:3400')
        !list=setup%get_strs('PARAMETER',o_default='tilD:-0.001:0.001')
        !list=setup%get_strs('PARAMETER',o_default='tilrho:-1000:1000')
                
        self%npars=size(list)
        allocate(self%pars(self%npars))

        !re-count to remove illegal parameters
        self%npars=0
        loop: do i=1,size(list)
            sublist=split(list(i)%s,o_sep=':') !=[name, min, max]

            select case (sublist(1)%s)
            ! case ('vp2' )
            !     self%pars(i)%name='vp2'
            !     self%npars=self%npars+1

            case ('vp' )
                self%pars(i)%name='vp'
                self%npars=self%npars+1

            case ('rho')
                ! if(is_gardner) then
                !     call hud('rho in PARAMETER is neglected and becomes passive as Gardner law is used')
                !     cycle loop
                ! endif
                self%pars(i)%name='rho'
                self%npars=self%npars+1

            case ('tilD')
                self%pars(i)%name='tilD'
                self%npars=self%npars+1

            case ('tilrho')
                self%pars(i)%name='tilrho' !tilrho:=tilD/buo^2=rho^2*tilD, ranging -1000:1000, in [kg/m3]
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
                        ! case ('vp2' ); o_x(:,:,:,i) = (m%vp**2 -self%pars(i)%min)/self%pars(i)%range
                        case ('vp'  ); o_x(:,:,:,i) = (m%vp  -self%pars(i)%min)/self%pars(i)%range
                        case ('rho');  o_x(:,:,:,i) = (m%rho -self%pars(i)%min)/self%pars(i)%range
                        case ('tilD'); o_x(:,:,:,i) = (m%tilD-self%pars(i)%min)/self%pars(i)%range
                        case ('tilrho'); o_x(:,:,:,i) = (m%tilD*m%rho**2-self%pars(i)%min)/self%pars(i)%range
                        end select
                enddo

            else !x->m
                do i=1,self%npars
                    select case (self%pars(i)%name)
                    ! case ('vp2' ); m%vp = sqrt(o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min)
                    case ('vp'  ); m%vp  = o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                    case ('rho');  m%rho = o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                    case ('tilD'); m%tilD= o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min
                    case ('tilrho'); m%tilD= o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min; m%tilD=m%tilD/m%rho**2
                    end select
                enddo
                ! ! + gardner
                ! if(is_gardner) m%rho = a*m%vp**b
                ! ! + castagna
                ! if(is_castagna) m%vs = a*m%vp + b
                ! call m%apply_freeze_zone

            endif

        endif

        ! if(present(o_xprior)) then
        !     call alloc(o_xprior,self%n1,self%n2,self%n3,self%npars)

        !     do i=1,self%npars
        !         select case (self%pars(i)%name)
        !         case ('vp' ); o_x(:,:,:,i) = (m%vp_prior -self%pars(i)%min)/self%pars(i)%range
        !         case ('vs' ); o_x(:,:,:,i) = (m%vs_prior -self%pars(i)%min)/self%pars(i)%range
        !         case ('rho'); o_x(:,:,:,i) = (m%rho_prior-self%pars(i)%min)/self%pars(i)%range
        !         end select
        !     enddo
        ! endif

        if(present(o_g)) then
            call alloc(o_g,self%n1,self%n2,self%n3,self%npars)

                do i=1,param%npars
                    select case (param%pars(i)%name)
                    !correlate_gradient(:,:,:,1) = gikpa
                    !correlate_gradient(:,:,:,2) = gbuo
                    !correlate_gradient(:,:,:,2) = gtilD
                    case ('vp' ); o_g(:,:,:,i) = correlate_gradient(:,:,:,1)*(-2.)/m%rho/(m%vp**3)
                    case ('rho'); o_g(:,:,:,i) = -m%rho**(-2)*( &
                        correlate_gradient(:,:,:,1)/(m%vp**2) + correlate_gradient(:,:,:,2) )
                    case ('tilD'); o_g(:,:,:,i) = correlate_gradient(:,:,:,3)
                    case ('tilrho'); o_g(:,:,:,i) = correlate_gradient(:,:,:,3)/m%rho**2
                    end select
                enddo


            !m%gradient(:,:,:,1) = grho
            !m%gradient(:,:,:,2) = gkpa or glda
            !m%gradient(:,:,:,3) = gmu

            ! !acoustic
            ! if(is_AC .and. .not. is_empirical) then
            !     do i=1,param%npars
            !         select case (param%pars(i)%name)
            !         case ('vp' ); o_g(:,:,:,i) = m%gradient(:,:,:,2)*2*m%rho*m%vp
            !         case ('rho'); o_g(:,:,:,i) = m%gradient(:,:,:,2)*m%vp**2 + m%gradient(:,:,:,1)
            !         end select
            !     enddo
            ! endif

            ! !acoustic + gardner
            ! if(is_AC .and. is_gardner) then
            !     o_g(:,:,:,1) =(m%gradient(:,:,:,2)*(b+2)/b*m%vp**2 + m%gradient(:,:,:,1))*b*m%rho/m%vp  !m%rho has tobe updated in prior
            ! endif

            ! !elastic
            ! if(is_EL .and. .not. is_empirical) then
            !     do i=1,param%npars
            !         select case (param%pars(i)%name)
            !         case ('vp' ); o_g(:,:,:,i) = m%gradient(:,:,:,2)*2*m%rho*m%vp
            !         case ('vs' ); o_g(:,:,:,i) =(m%gradient(:,:,:,2)*(-2) + m%gradient(:,:,:,3))*2*m%rho*m%vs
            !         case ('rho'); o_g(:,:,:,i) = m%gradient(:,:,:,2)*m%vp**2 + (-2*m%gradient(:,:,:,2)+m%gradient(:,:,:,3))*m%vs**2 + m%gradient(:,:,:,1)
            !         end select
            !     enddo
            ! endif

            ! !elastic + gardner
            ! if(is_EL .and. is_gardner) then
            !     do i=1,param%npars
            !         select case (param%pars(i)%name)
            !         case ('vp' ); o_g(:,:,:,i) =(m%gradient(:,:,:,2)*(b+2)/b*m%vp**2 + (-2*m%gradient(:,:,:,2)+m%gradient(:,:,:,3))*m%vs**2 + m%gradient(:,:,:,1))*b*m%rho/m%vp !m%rho has tobe updated in prior
            !         case ('vs' ); o_g(:,:,:,i) =(m%gradient(:,:,:,2)*(-2) + m%gradient(:,:,:,3))*2*m%rho*m%vs !m%rho has tobe updated in prior
            !         end select
            !     enddo
            ! endif

            ! !elastic + castagna
            ! if(is_EL .and. is_castagna) then
            !     do i=1,param%npars
            !         select case (param%pars(i)%name)
            !         case ('vp' )
            !             o_g(:,:,:,i) =2*m%rho*( m%gradient(:,:,:,2)*(m%vp-2*a*m%vs) + m%gradient(:,:,:,3)*a*m%vs )
            !         case ('rho')
            !             o_g(:,:,:,i) =m%gradient(:,:,:,2)*(m%vp**2-2*m%vs**2) + m%gradient(:,:,:,3)*m%vs**2 + m%gradient(:,:,:,1)
            !         end select
            !     enddo
            ! endif

            !normaliz g by allowed parameter range
            do i=1,param%npars
                o_g(:,:,:,i)=o_g(:,:,:,i)*param%pars(i)%range
            enddo
            
        endif
        
!deprecated
!         if(present(o_pg)) then
!             call alloc(o_pg,param%n1,param%n2,param%n3,param%npars,oif_protect=.true.)
!             call transform_gradient(preco%preco,o_preco)
!         endif

    end subroutine

    subroutine transform_preconditioner(self,preco_in_m,preco_in_x)
        class(t_parametrizer) :: self
        real,dimension(m%nz,m%nx,m%ny) :: preco_in_m
        real,dimension(self%n1,self%n2,self%n3) :: preco_in_x

        preco_in_x=preco_in_m
    end subroutine

end module
