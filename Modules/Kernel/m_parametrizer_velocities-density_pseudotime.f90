module m_parametrizer
use m_System
use m_pseudotime
use m_Modeling

    !PARAMETERIZATION     -- ALLOWED PARAMETERS
    !velocities-density   -- vp rho

    !models (m) in depth (z) domain
    !parameters (x) in pseudotime (t) domain

    !acoustic:
    !kpa = rho*vp^2 = vp*ip
    !rho0= rho      = ip/vp
    !gvp = gkpa*2rho*vp
    !grho= gkpa*vp^2 + grho0

    !acoustic + gardner:
    !kpa = a*vp^(b+2)
    !rho0= a*vp^b
    !gvp = (gkpa*(b+2)/b*vp^2 + grho0)*ab*vp^(b-1)

    private

    type t_parameter
        character(:),allocatable :: name
        real :: min, max, range
    end type

    type,public :: t_parametrizer
        !info
        character(i_str_xxlen) :: info = &
            'Parameterization: velocities-density in pseudotime'//s_NL// &
            'Allowed pars: vp, rho'//s_NL// &
            'Available empirical laws: Gardner'

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

    integer :: index_vp !index of vp in param list

    logical :: is_empirical=.false., is_gardner=.false. !, is_castagna

    real :: a,b

    real,dimension(:,:,:),allocatable :: freeze_zone_in_m, freeze_zone_in_x

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
                        a=either(310.,0.31,m%rho(1,1,1)<1000.)
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
        
        
        !check PDE
        if(index(ppg%info,'EL')>0) call error('Current implementation of pseudotime does not consider elastic inversion! Sorry.')

        !read in active parameters and their allowed ranges
        list=setup%get_strs('PARAMETER',o_default='vp:1500:3400'); index_vp=1
        
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
                index_vp=i
                vmin=str2real(sublist(2)%s)
                vmax=str2real(sublist(3)%s)

            case ('rho')
                if(is_empirical) then
                    call hud('rho in PARAMETER is neglected as EMPIRICAL_LAW is read (rho becomes a passive parameter).')
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
        
        call pseudotime_init('z->t',vmin,vmax,m%nx,m%ny,&
            nz_=m%nz,   Dz_=m%dz, &
            nt_=self%n1,Dt_=self%d1)

        call hud('pseudotime dimension nt, dt = '//num2str(self%n1)//' , '//num2str(self%d1))

        self%n2=m%nx
        self%n3=m%ny
        self%n=self%n1*self%n2*self%n3*self%npars

        self%d2=m%dx
        self%d3=m%dy

        call alloc(freeze_zone_in_m,m%nz,m%nx,m%ny,o_init=1.)
        where(m%is_freeze_zone) freeze_zone_in_m=0.
        call pseudotime_convert('z->t',freeze_zone_in_m,freeze_zone_in_x,o_v=m%vp)
        deallocate(freeze_zone_in_m)

    end subroutine
    
    subroutine transform(self,o_dir,o_x,o_xprior,o_g)
        class(t_parametrizer) :: self
        character(4),optional :: o_dir
        real,dimension(:,:,:,:),allocatable,optional :: o_x,o_xprior,o_g
        
        real,dimension(:,:,:),allocatable :: v_t !velocity model in pseudotime domain
        real,dimension(:,:,:),allocatable :: tmp

        if(present(o_x)) then
            call alloc(o_x,self%n1,self%n2,self%n3,self%npars,oif_protect=.true.)

            if(either(o_dir,'m->x',present(o_dir))=='m->x') then !z->t
                do i=1,self%npars
                    select case (self%pars(i)%name)
                    case ('vp' )
                        call pseudotime_convert('z->t',m%vp,v_t)
                        o_x(:,:,:,i) = (v_t-self%pars(i)%min)/self%pars(i)%range
                    case ('rho')
                        call pseudotime_convert('z->t',m%rho,tmp,o_v=m%vp)
                        o_x(:,:,:,i) = (tmp-self%pars(i)%min)/self%pars(i)%range
                    end select
                enddo

            else !x->m, t->z
                !first convert velocity
                v_t=o_x(:,:,:,index_vp)*self%pars(index_vp)%range +self%pars(index_vp)%min
                call pseudotime_convert('t->z',v_t,m%vp)
                
                !then convert other parameters
                do i=1,self%npars
                    select case (self%pars(i)%name)
                    case ('rho')
                        call pseudotime_convert('t->z',o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min,m%rho, o_v=v_t)
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
                case ('vp' )
                    call pseudotime_convert('z->t',m%vp_prior,tmp,o_v=m%vp)
                    o_x(:,:,:,i) = (tmp -self%pars(i)%min)/self%pars(i)%range
                case ('rho')
                    call pseudotime_convert('z->t',m%rho_prior,tmp,o_v=m%vp)
                    o_x(:,:,:,i) = (m%rho_prior-self%pars(i)%min)/self%pars(i)%range
                end select
            enddo
        endif

        if(present(o_g)) then
            call alloc(o_g,self%n1,self%n2,self%n3,self%npars)
            !m%gradient(:,:,:,1) = grho
            !m%gradient(:,:,:,2) = gkpa

            !acoustic
            ! if(is_AC .and. .not. is_empirical) then
            if(.not. is_empirical) then
                do i=1,self%npars
                    select case (self%pars(i)%name)
                    case ('vp' )
                        call pseudotime_convert_gradient(m%gradient(:,:,:,2)*2*m%rho*m%vp, &
                            m%vp,tmp)
                        o_g(:,:,:,i) = tmp
                    case ('rho')
                        call pseudotime_convert_gradient(m%gradient(:,:,:,2)*m%vp**2 + m%gradient(:,:,:,1), &
                            m%vp,tmp)
                        o_g(:,:,:,i) = tmp
                    end select
                enddo
            endif

            !acoustic + gardner
            ! if(is_AC .and. is_gardner) then
            if(is_gardner) then
                call pseudotime_convert_gradient((m%gradient(:,:,:,2)*(b+2)/b*m%vp**2 + m%gradient(:,:,:,1))*a*b*m%vp**(b-1), &
                                            m%vp,tmp)
                o_g(:,:,:,1) = tmp
            endif

            !normaliz g by allowed parameter range
            do i=1,self%npars
                o_g(:,:,:,i)=o_g(:,:,:,i)*self%pars(i)%range
            enddo
            
            !apply bathymetry
            do i=1,self%npars
                o_g(:,:,:,i)=o_g(:,:,:,i)*freeze_zone_in_x
            enddo

        endif
        
        call dealloc(tmp,v_t)

    end subroutine

    subroutine transform_preconditioner(self,preco_in_m,preco_in_x)
        class(t_parametrizer) :: self
        real,dimension(m%nz,m%nx,m%ny) :: preco_in_m
        real,dimension(:,:,:),allocatable :: preco_in_x

        call pseudotime_convert('z->t',preco_in_m,preco_in_x,o_v=m%vp)
        
    end subroutine

end module
