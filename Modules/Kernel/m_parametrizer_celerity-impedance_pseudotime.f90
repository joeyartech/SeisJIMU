module m_parametrizer
use m_System
use m_math
use m_pseudotime
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
            'Parameterization: celerity-impedance in pseudotime'//s_NL// &
            'Allowed pars: cel, imp'
!In electrical engineering, the electrical impedance Z is the total opposition to an alternating current (AC) in a circuit.
!In electromagnetics, the intrinsic impedance η describes how a medium interacts with an open wave.

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

    integer :: i_cel=0, i_imp=0, i_sgma=0

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
                vmin=str2real(sublist(2)%s)
                vmax=str2real(sublist(3)%s)

            case ('imp')
                i_imp=i
                self%pars(i)%name='imp'
                self%npars=self%npars+1

            ! case ('sgma')
            !     i_sgma=i
            !     self%pars(i)%name='sgma'
            !     self%npars=self%npars+1
                
            end select

            self%pars(i)%min=str2real(sublist(2)%s)
            self%pars(i)%max=str2real(sublist(3)%s)
            self%pars(i)%range=self%pars(i)%max-self%pars(i)%min

        enddo loop
        
        deallocate(list,sublist)

        !check vp,vs,rho [min,max] is in the same range as m%ref_vp,ref_vs,ref_rho
        !

        ! call empirical_init
        
        call pseudotime_init('z->t',vmin,vmax,m%nx,m%ny,&
            nz_=m%nz,   Dz_=m%dz, &
            nt_=self%n1,Dt_=self%d1)

        call hud('pseudotime dimension nt, dt = '//num2str(self%n1)//' , '//num2str(self%d1))
        
        self%n2=m%nx
        self%n3=m%ny
        self%n=self%n1*self%n2*self%n3*self%npars

        self%d2=m%dx
        self%d3=m%dy

    end subroutine
    

    subroutine transform(self,o_dir,o_x,o_xprior,o_g)
        class(t_parametrizer) :: self
        character(4),optional :: o_dir
        real,dimension(:,:,:,:),allocatable,optional :: o_x,o_xprior,o_g

        real,dimension(:,:,:),allocatable :: v_t !velocity model in pseudotime domain
        real,dimension(:,:,:),allocatable :: tmp_imp, tmp_cel, tmp_gimp, tmp_gcel, tmp

        !c⁻² = εμ ; η² = (cμ)² = μ/ε ==>  
        !ε = 1/(ηc) ; μ = η/c


        if(present(o_x)) then
            call alloc(o_x,self%n1,self%n2,self%n3,self%npars,oif_protect=.true.)

            tmp_cel = r_c0/sqrt(m%epsr*m%mur)

            if(either(o_dir,'m->x',present(o_dir))=='m->x') then !z->t
                if(i_cel>0) then
                    call pseudotime_convert('z->t',tmp_cel,v_t)
                    o_x(:,:,:,i_cel) =( v_t-self%pars(i_cel)%min )/self%pars(i_cel)%range
                endif
                if(i_imp>0) then
                    call pseudotime_convert('z->t',sqrt( (r_mu0*m%mur)/(r_eps0*m%epsr) ), &
                                                    tmp_imp,o_v=tmp_cel)
                    o_x(:,:,:,i_imp) =( tmp_imp-self%pars(i_imp)%min )/self%pars(i_imp)%range
                endif
                ! if(i_sgma>0) then
                !     o_x(:,:,:,i_sgma)=(m%sgma-self%pars(i_sgma)%min)/self%pars(i_sgma)%range
                ! endif

                ! call empirical_m2x('velocities-density')

            else !x->m, t->z
                ! if(i_cel>0) then
                
                !first convert velocity
                v_t=o_x(:,:,:,i_cel)*self%pars(i_cel)%range +self%pars(i_cel)%min
                call pseudotime_convert('t->z',v_t,tmp_cel)
                
                tmp_imp =  sqrt( (r_mu0*m%mur)/(r_eps0*m%epsr) )
                m%epsr = 1/(tmp_imp*tmp_cel) /r_eps0
                m%mur  =    tmp_imp/tmp_cel  /r_mu0

                if(allocated(m%epsr0).and.allocated(m%mur0)) then
                    call hud('transform m%epsr0 m%mur0 too')
                    tmp_imp =  sqrt( (r_mu0*m%mur0)/(r_eps0*m%epsr0) )
                    m%epsr0 = 1/(tmp_imp*tmp_cel) /r_eps0
                    m%mur0  =    tmp_imp/tmp_cel  /r_mu0
                endif

                ! if(i_imp>0) then
                !     call pseudotime_convert('t->z',o_x(:,:,:,i_imp)*self%pars(i_imp)%range +self%pars(i_imp)%min, &
                !             tmp_imp, o_v=v_t)
                !     tmp_cel = r_c0/sqrt(m%epsr*m%mur)

                !     m%epsr = 1/(tmp_imp*tmp_cel) /r_eps0
                !     m%mur  =    tmp_imp/tmp_cel  /r_mu0
                !     if(allocated(m%epsr0).and.allocated(m%mur0)) then
                !         call hud('transform m%epsr0 m%mur0 too')
                !         tmp_cel = r_c0/sqrt(m%epsr0*m%mur0)
                !         m%epsr0 = 1/(tmp_imp*tmp_cel) /r_eps0
                !         m%mur0  =    tmp_imp/tmp_cel  /r_mu0
                !     endif
                ! endif

                ! call empirical_x2m('velocities-density')

                ! endif
                
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
            
            !c⁻² = (ε₀μ₀)εᵣμᵣ ; η² = (μ₀/ε₀)μᵣ/εᵣ ==>  
            !εᵣ = ε₀⁻¹η⁻¹c⁻¹ ; μᵣ = μ₀⁻¹ηc⁻¹
            !So,
            !gcel = gepsr ∂εᵣ/∂c + gmur ∂μᵣ/∂c = gepsr (-1)ε₀⁻¹η⁻¹c⁻² + gmur (-1)μ₀⁻¹ηc⁻²
            !     =-( (gepsr/ε₀)/η + (gmur/μ₀)η )c⁻²
            !gimp = gepsr ∂εᵣ/∂η + gmur ∂μᵣ/∂η = gepsr (-1)ε₀⁻¹η⁻²c⁻¹ + gmur μ₀⁻¹c⁻¹
            !     = (-(gepsr/ε₀)/η² + gmur/μ₀ )c⁻¹

            tmp_cel = r_c0/sqrt(m%epsr*m%mur)
            tmp_imp = sqrt( (r_mu0*m%mur)/(r_eps0*m%epsr) )

            tmp_gcel =-( correlate_gradient(:,:,:,2)/r_eps0/tmp_imp + correlate_gradient(:,:,:,1)/r_mu0*tmp_imp ) &
                        /tmp_cel/tmp_cel
            tmp_gimp = (-correlate_gradient(:,:,:,2)/r_eps0/tmp_imp/tmp_imp + correlate_gradient(:,:,:,1)/r_mu0) &
                        /tmp_cel

            if(i_cel>0) then
                call pseudotime_convert_gradient(tmp_gcel,tmp_cel,tmp)
                o_g(:,:,:,i_cel) = tmp
            endif

            if(i_imp>0) then
                call pseudotime_convert_gradient(tmp_gimp,tmp_cel,tmp)
                o_g(:,:,:,i_imp) = tmp
            endif


            !normaliz g by allowed parameter range
            !s.t. g is in unit [Nm]
            do i=1,param%npars
                o_g(:,:,:,i)=o_g(:,:,:,i)*param%pars(i)%range
            enddo

            ! !apply bathymetry
            ! call pseudotime_convert('z->t',bools2reals(m%is_freeze_zone),freeze_zone_in_t,o_v=m%vp)
            ! do i=1,self%npars
            !    o_g(:,:,:,i)=o_g(:,:,:,i)*(1.-freeze_zone_in_t)
            ! enddo
            
        endif

        call dealloc(tmp_cel,tmp_imp,tmp_gcel,tmp_gimp,tmp)
        call dealloc(v_t)

    end subroutine

    subroutine transform_preconditioner(self,preco_in_m,preco_in_x)
        class(t_parametrizer) :: self
        real,dimension(m%nz,m%nx,m%ny) :: preco_in_m
        real,dimension(:,:,:),allocatable :: preco_in_x

        real,dimension(:,:,:),allocatable :: tmp_cel

        tmp_cel = r_c0/sqrt(m%epsr*m%mur)
        call pseudotime_convert('z->t',preco_in_m,preco_in_x,o_v=tmp_cel)

    end subroutine

end module  