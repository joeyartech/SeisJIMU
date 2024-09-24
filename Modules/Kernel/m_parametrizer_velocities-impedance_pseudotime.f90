module m_parametrizer
use m_System
use m_pseudotime
use m_Modeling

    !PARAMETERIZATION     -- ALLOWED PARAMETERS
    !velocities-impedance -- vp vs ip

    !acoustic:
    !kpa = rho*vp^2 = vp*ip
    !rho0= rho      = ip/vp
    !gvp = (gkpa*vp - grho0/vp)*rho
    !gip =  gkpa*vp + grho0/vp

    private

    type t_parameter
        character(:),allocatable :: name
        real :: min, max, range
    end type

    type,public :: t_parametrizer
        !info
        character(i_str_xxlen) :: info = &
            'Parameterization: velocities-impedance in pseudotime'//s_NL// &
            'Allowed pars: vp, ip'

        type(t_parameter),dimension(:),allocatable :: pars
        integer :: npars
        !type(t_string),dimension(:),allocatable :: empirical

        logical,dimension(:,:,:,:),allocatable :: is_freeze_zone

        integer :: n1,n2,n3,n
        real :: d1,d2,d3
        
        contains
        procedure :: init
        procedure :: transform
        procedure :: transform_preconditioner
    end type

    type(t_parametrizer),public :: param

    integer :: index_vp !index of vp in param list

    !logical :: is_empirical=.false., is_gardner=.false. !, is_castagna

    real :: a,b

    contains
    
    subroutine init(self)
        class(t_parametrizer) :: self

        type(t_string),dimension(:),allocatable :: list,sublist

        call hud('Invoked parametrizer module info : '//s_NL//self%info)

        if(allocated(   list)) deallocate(   list)
        if(allocated(sublist)) deallocate(sublist)
        
        
        !check PDE
        if(index(ppg%info,'EL')>0) call error('Current implementation of pseudotime does not consider elastic inversion! Sorry.')


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
                index_vp=i
                vmin=str2real(sublist(2)%s)
                vmax=str2real(sublist(3)%s)

            case ('ip')
                self%pars(i)%name='ip'
                self%npars=self%npars+1
                
            end select

            self%pars(i)%min=str2real(sublist(2)%s)
            self%pars(i)%max=str2real(sublist(3)%s)
            self%pars(i)%range=self%pars(i)%max-self%pars(i)%min

        enddo loop
                
        deallocate(list,sublist)

        
        call pseudotime_init('z->t',vmin,vmax,m%nx,m%ny,&
            nz_=m%nz,   Dz_=m%dz, &
            nt_=self%n1,Dt_=self%d1)

        call hud('pseudotime dimension nt, dt = '//num2str(self%n1)//' , '//num2str(self%d1))

        self%n2=m%nx
        self%n3=m%ny
        self%n=self%n1*self%n2*self%n3*self%npars

        self%d2=m%dx
        self%d3=m%dy

        !call alloc(freeze_zone_in_m,m%nz,m%nx,m%ny,o_init=1.)
        !where(m%is_freeze_zone) freeze_zone_in_m=0.
        !call pseudotime_convert('z->t',freeze_zone_in_m,freeze_zone_in_x,o_v=m%vp)
        !deallocate(freeze_zone_in_m)

    end subroutine
    
    subroutine transform(self,o_dir,o_x,o_xprior,o_g)
        class(t_parametrizer) :: self
        character(4),optional :: o_dir
        real,dimension(:,:,:,:),allocatable,optional :: o_x,o_xprior,o_g
        
        real,dimension(:,:,:),allocatable :: v_t !velocity model in pseudotime domain
        real,dimension(:,:,:),allocatable :: tmp
        real,dimension(:,:,:),allocatable :: freeze_zone_in_t !gradient in pseudotime domain

        if(present(o_x)) then
            call alloc(o_x,self%n1,self%n2,self%n3,self%npars,oif_protect=.true.)

            if(either(o_dir,'m->x',present(o_dir))=='m->x') then !z->t
                do i=1,self%npars
                    select case (self%pars(i)%name)
                    case ('vp')
                        call pseudotime_convert('z->t',m%vp,v_t)
                        o_x(:,:,:,i) = (v_t-self%pars(i)%min)/self%pars(i)%range
                    case ('ip')
                        call pseudotime_convert('z->t',m%vp*m%rho,tmp,o_v=m%vp)
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
                    case ('ip')
                        call pseudotime_convert('t->z',o_x(:,:,:,i)*self%pars(i)%range +self%pars(i)%min,tmp, o_v=v_t)
                        m%rho = tmp/m%vp
                    end select
                enddo
                
                call m%apply_freeze_zone

            endif

        endif

        ! if(present(o_xprior)) then
        !     call alloc(o_xprior,self%n1,self%n2,self%n3,self%npars)

        !     do i=1,self%npars
        !         select case (self%pars(i)%name)
        !         case ('vp' )
        !             call pseudotime_convert('z->t',m%vp_prior,tmp,o_v=m%vp)
        !             o_x(:,:,:,i) = (tmp -self%pars(i)%min)/self%pars(i)%range
        !         case ('rho')
        !             call pseudotime_convert('z->t',m%rho_prior,tmp,o_v=m%vp)
        !             o_x(:,:,:,i) = (m%rho_prior-self%pars(i)%min)/self%pars(i)%range
        !         end select
        !     enddo
        ! endif

        if(present(o_g)) then
            call alloc(o_g,self%n1,self%n2,self%n3,self%npars)
            !correlate_gradient(:,:,:,1) = grho0
            !correlate_gradient(:,:,:,2) = gkpa

            !acoustic
            do i=1,self%npars
                select case (self%pars(i)%name)
                case ('vp' )
                    call pseudotime_convert_gradient((correlate_gradient(:,:,:,2)*m%vp - correlate_gradient(:,:,:,1)/m%vp)*m%rho, &
                        m%vp,tmp)

                case ('ip')
                    call pseudotime_convert_gradient(correlate_gradient(:,:,:,2)*m%vp + correlate_gradient(:,:,:,1)/m%vp, &
                        m%vp,tmp)
                end select
                o_g(:,:,:,i) = tmp
            enddo

            !normaliz g by allowed parameter range
            do i=1,self%npars
                o_g(:,:,:,i)=o_g(:,:,:,i)*self%pars(i)%range
            enddo


            !apply bathymetry
            call pseudotime_convert('z->t',bools2reals(m%is_freeze_zone),freeze_zone_in_t,o_v=m%vp)
            do i=1,self%npars
               o_g(:,:,:,i)=o_g(:,:,:,i)*(1.-freeze_zone_in_t)
            enddo

            
        endif
        
        call dealloc(tmp,v_t,freeze_zone_in_t)

    end subroutine

    subroutine transform_preconditioner(self,preco_in_m,preco_in_x)
        class(t_parametrizer) :: self
        real,dimension(m%nz,m%nx,m%ny) :: preco_in_m
        real,dimension(:,:,:),allocatable :: preco_in_x

        call pseudotime_convert('z->t',preco_in_m,preco_in_x,o_v=m%vp)
        
    end subroutine

end module
