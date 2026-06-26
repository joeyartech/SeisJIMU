module m_parametrizer
use m_System
use m_math
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
            'Parameterization: celerity-permeability'//s_NL// &
            'Allowed pars: cel, mu, sgma'

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

    integer :: i_cel=0, i_mu=0, i_sgma=0

    contains
    
    subroutine init(self)
        class(t_parametrizer) :: self

        type(t_string),dimension(:),allocatable :: list,sublist

        call hud('Invoked parametrizer module info : '//s_NL//self%info)

        if(allocated(   list)) deallocate(   list)
        if(allocated(sublist)) deallocate(sublist)
        
        
        ! !check basic gradients provided from propagator
        ! is_grho = index(ppg%info,'grho')>0
        ! is_gbuo = index(ppg%info,'gbuo')>0
        ! is_gkpa = index(ppg%info,'gkpa')>0 
        ! is_gikpa= index(ppg%info,'gikpa')>0
        ! is_glda = index(ppg%info,'glda')>0
        ! is_gmu  = index(ppg%info,'gmu')>0
        ! is_gimu = index(ppg%info,'gimu')>0

        !read in active parameters and their allowed ranges
        list=setup%get_strs('PARAMETER',o_default='cel:55:299')
        
        self%npars=size(list)
        allocate(self%pars(self%npars))

        !re-count to remove illegal parameters
        self%npars=0
        loop: do i=1,size(list)
            sublist=split(list(i)%s,o_sep=':') !=[name, min, max]

            select case (sublist(1)%s)
            case ('cel' )
                i_cel=i
                self%pars(i)%name='cel'
                self%npars=self%npars+1

            case ('mu')
                i_mu=i
                self%pars(i)%name='mu'
                self%npars=self%npars+1

            case ('sgma')
                i_sgma=i
                self%pars(i)%name='sgma'
                self%npars=self%npars+1
                
            end select

            self%pars(i)%min=str2real(sublist(2)%s)
            self%pars(i)%max=str2real(sublist(3)%s)
            self%pars(i)%range=self%pars(i)%max-self%pars(i)%min

        enddo loop
        
        deallocate(list,sublist)

        !check vp,vs,rho [min,max] is in the same range as m%ref_vp,ref_vs,ref_rho
        !

        ! call empirical_init

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

        real,dimension(:,:,:),allocatable :: tmp_cel

        !c⁻² = εμ = ε₀εᵣμ  ==>  εᵣ = c⁻²/(ε₀μ) = c⁻²/(ε₀μ₀μᵣ) ; μᵣ = μ/μ₀
        
        call alloc(tmp_cel, m%nz,m%nx,m%ny)

        if(present(o_x)) then
            call alloc(o_x,self%n1,self%n2,self%n3,self%npars,oif_protect=.true.)

            if(either(o_dir,'m->x',present(o_dir))=='m->x') then
                if(i_cel>0) then
                    tmp_cel = r_c0/sqrt(m%epsr*m%mur)
                    o_x(:,:,:,i_cel) = (tmp_cel     -self%pars(i_cel )%min)/self%pars(i_cel )%range
                endif
                if(i_mu >0) then
                    o_x(:,:,:,i_mu ) = (r_mu0*m%mur -self%pars(i_mu  )%min)/self%pars(i_mu  )%range
                endif
                if(i_sgma>0) then
                    o_x(:,:,:,i_sgma) = (     m%sgma-self%pars(i_sgma)%min)/self%pars(i_sgma)%range
                endif

                ! call empirical_m2x('velocities-density')

            else !x->m
                if(i_mu >0) then
                    tmp_cel = r_c0/sqrt(m%epsr*m%mur)
                    m%mur  = (o_x(:,:,:,i_mu )*self%pars(i_mu )%range +self%pars(i_mu )%min) / r_mu0
                    m%epsr = (r_c0/tmp_cel)**2/m%mur
                endif
                if(i_cel>0) then
                    tmp_cel=  o_x(:,:,:,i_cel)*self%pars(i_cel)%range +self%pars(i_cel)%min
                    m%epsr =  1./tmp_cel/tmp_cel/r_eps0mu0/m%mur
                endif
                if(i_sgma>0) then
                    m%sgma = (o_x(:,:,:,i_sgma)*self%pars(i_sgma)%range +self%pars(i_sgma)%min)
                endif

                ! call empirical_x2m('velocities-density')
                
                call m%apply_relativity
                call m%apply_freeze_zone

            endif

        endif

        ! if(present(o_xprior)) then
        !     call alloc(o_xprior,self%n1,self%n2,self%n3,self%npars)

        !     if(i_vp >0) o_x(:,:,:,i_vp ) = (m%vp_prior -self%pars(i_vp )%min)/self%pars(i_vp )%range
        !     if(i_vs >0) o_x(:,:,:,i_vs ) = (m%vs_prior -self%pars(i_vs )%min)/self%pars(i_vs )%range
        !     if(i_rho>0) o_x(:,:,:,i_rho) = (m%rho_prior-self%pars(i_rho)%min)/self%pars(i_rho)%range
                
        !     call empirical_m2x('velocities-density')

        ! endif

        if(present(o_g)) then
            call alloc(o_g,self%n1,self%n2,self%n3,self%npars)

            ! call hud('Parametrizer finds grho & gkpa')
            !correlate_gradient(:,:,:,1) = gmur
            !correlate_gradient(:,:,:,2) = gepsr
            !correlate_gradient(:,:,:,3) = gsgma
            
            !c⁻² = εμ = ε₀εᵣμ  ==>  εᵣ = c⁻²μ⁻¹/ε₀; μᵣ = μ/μ₀
            !So,
            !gcel = gepsr ∂εᵣ/∂c + gmur ∂μᵣ/∂c = gepsr (-2)(c⁻³μ⁻¹/ε₀)
            !     = gepsr (-2)(εᵣ/c)
            !gmu  = gepsr ∂εᵣ/∂μ + gmur ∂μᵣ/∂μ = gepsr (-1)(c⁻²μ⁻²/ε₀) + gmur 1/μ₀
            !     = gepsr (-1)(εᵣμ⁻¹) + gmur μ₀⁻¹

            tmp_cel = r_c0/sqrt(m%epsr*m%mur)

            if(i_cel >0) then
                o_g(:,:,:,i_cel ) = correlate_gradient(:,:,:,2)*(-2)*m%epsr/tmp_cel
            endif

            if(i_mu  >0) then
                o_g(:,:,:,i_mu  ) = correlate_gradient(:,:,:,2)*(-1)*m%epsr/m%mur/r_mu0 &
                                   +correlate_gradient(:,:,:,1)/r_mu0
            endif

            if(i_sgma>0) then
                o_g(:,:,:,i_sgma) = correlate_gradient(:,:,:,3)
            endif
            
            ! if(i_sgma>0) o_g(:,:,:,i_sgma) = correlate_gradient(:,:,:,3)

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