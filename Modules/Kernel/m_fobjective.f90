module m_fobjective
use m_System
use m_Modeling
use m_weighter
use m_parametrizer

    private

    type,public :: t_querypoint
        
        real,dimension(:,:,:,:),allocatable :: x  !point in parameter space X (unitless)
        real,dimension(:,:,:,:),allocatable :: g  !gradient, unit in [1/Nm]
        real,dimension(:,:,:,:),allocatable :: pg !preconditioned gradient
        real,dimension(:,:,:,:),allocatable :: d  !descent direction
        real :: f !objective function value, unit in [Nm]
        real :: gdotd  !g.d

        contains
        procedure :: init => querypoint_init
        final :: querypoint_fin

    end type

    real,dimension(:,:,:,:),allocatable :: xprior !prior parameter

    type,public :: t_fobjective

        type(t_string),dimension(:),allocatable :: s_dnorms !norms for data residuals
        type(t_string),dimension(:),allocatable :: s_xnorms !norms for model residuals

        real,dimension(:),allocatable :: dnorms !unit in [Nm]
        real,dimension(:),allocatable :: xnorms !unit in [Nm]

        integer :: n_dnorms, n_xnorms
        integer :: i_dnorm_4adjsource !which dnorm is used for the adjoint source (only 1 is allowed so far)

        real,dimension(:),allocatable :: xpenalizers

        contains

        procedure :: init => fobjective_init
        procedure :: compute_dnorms
        procedure :: compute_xnorms
        procedure :: total_loss
        procedure :: print_dnorms
        procedure :: print_xnorms

    end type

    type(t_fobjective),public :: fobj

    real,dimension(:,:),allocatable :: Wdres
    real,dimension(:,:,:,:),allocatable :: Wxres
    real :: xnorms_weights(4) !weights on 1st, 2nd, 3rd & 4th dimensions

    contains

    subroutine querypoint_init(self,oif_g,oif_pg,oif_d)
        class(t_querypoint) :: self
        logical,optional :: oif_g,oif_pg,oif_d

        logical,save :: is_first_in=.true.

                    call alloc(self%x, param%n1,param%n2,param%n3,param%npars)
        if(oif_g)   call alloc(self%g, param%n1,param%n2,param%n3,param%npars)
        if(oif_pg)  call alloc(self%pg,param%n1,param%n2,param%n3,param%npars)
        if(oif_d)   call alloc(self%d, param%n1,param%n2,param%n3,param%npars)
        
        self%f=0.
        self%gdotd=0.

                    call param%transform('m->x',o_x=self%x)

        if(oif_g)   call param%transform(o_g=self%g)

        if(is_first_in) then
            !prior models
            call m%read_prior
            if(m%if_has_prior) then
                call alloc(xprior,param%n1,param%n2,param%n3,param%npars)
                !transform prior models
                call param%transform('m->x',o_xprior=xprior)
                !save some RAM
                call dealloc(m%vp_prior,m%vs_prior,m%rho_prior)
            endif

            is_first_in=.false.

        endif

    end subroutine

    subroutine querypoint_fin(self)
        type(t_querypoint) :: self

        call dealloc(self%x,self%g,self%pg,self%d)

    end subroutine


    subroutine fobjective_init(self)
        class(t_fobjective) :: self
        
        !data norms
        self%s_dnorms=setup%get_strs('DATA_NORMS','DNORMS',o_default='L1 =>L2 Linf')
        self%n_dnorms=size(self%s_dnorms)
        call alloc(self%dnorms,self%n_dnorms)

        do i=1,self%n_dnorms
            if(self%s_dnorms(i)%s(1:2)=='=>') then
                self%i_dnorm_4adjsource=i
                self%s_dnorms(i)%s(1:2)='  '; self%s_dnorms(i)%s=lalign(self%s_dnorms(i)%s) !remove the indicator '=>'
                exit
            endif
        enddo

        !parameter norms
        self%s_xnorms=setup%get_strs('PARAMETER_NORMS','XNORMS',o_default='none')
        self%n_xnorms=size(self%s_xnorms)
        call alloc(self%xnorms,self%n_xnorms)
        call alloc(self%xpenalizers,self%n_xnorms)

        if(self%n_xnorms>0) then
            xnorms_weights=setup%get_reals('PARAMETER_NORMS_WEIGHTS','XNORMS_WEI',o_default='1. 1. 1. 1.')
        endif
        
    end subroutine
    
    subroutine compute_dnorms(self)
        class(t_fobjective) :: self
        
        call alloc(Wdres,shot%nt,shot%nrcv)

        self%dnorms=0.

        do i=1,self%n_dnorms
            select case (self%s_dnorms(i)%s)
            case ('L2')
                Wdres = (shot%dsyn-shot%dobs)*wei%weight

                if(mpiworld%is_master) call shot%write('Wdres')

                do ir=1,shot%nrcv
                    if(shot%rcv(ir)%is_badtrace) cycle

                    self%dnorms(i) = self%dnorms(i) + 0.5*m%cell_volume*sum(Wdres(:,ir)**2)

                    !set unit of dnorm to be [Nm], same as Lagrangian
                    !this also help balance contributions from different component data
                    if(shot%rcv(ir)%comp=='p') then !for pressure data
                        self%dnorms(i) = self%dnorms(i) /m%ref_kpa
                    else !for velocities data
                        self%dnorms(i) = self%dnorms(i) *m%ref_rho
                    endif

                enddo

                !compute adjoint source and set proper units
                if(i==self%i_dnorm_4adjsource) then
                    shot%dadj = -wei%weight*Wdres*m%cell_volume

                    do ir=1,shot%nrcv
                        if(shot%rcv(ir)%comp=='p') then !for pressure data
                            shot%dadj(:,ir) = shot%dadj(:,ir) /m%ref_kpa
                        else !for velocities data
                            shot%dadj(:,ir) = shot%dadj(:,ir) *m%ref_rho
                        endif

                    enddo

                    call shot%write('dadj')

                endif
                
            case ('L1')
                !call Lp(dres,shot,wei,1)

            end select

        enddo

    end subroutine

    subroutine compute_xnorms(self,x,g)
        class(t_fobjective) :: self
        real,dimension(param%n1,param%n2,param%n3,param%npars) :: x,g

        real,dimension(:,:,:,:),allocatable :: spatgrad !spatial gradient of x

        call alloc(spatgrad,param%n1,param%n2,param%n3,param%npars,o_init=0.)

        self%xnorms=0.

        do i=1,self%n_xnorms
            select case (self%s_xnorms(i)%s)
            case ('Ridge','L2')
                !fcost
                call alloc(Wxres,param%n1,param%n2,param%n3,param%npars)
                if(m%if_has_prior) then
                    Wxres=x-xprior
                else
                    Wxres=x
                endif

                self%xnorms(i) = sum(Wxres*Wxres, .not.param%is_freeze_zone)

                !spatial gradient
                spatgrad = spatgrad + Wxres

            case ('Lasso','L1')

            case ('1st-ord_Ridge')

                !fcost
                do j=1,param%n1
                    self%xnorms(i) = self%xnorms(i) + xnorms_weights(1)*sum((x(j+1,:,:,:)-x(j,:,:,:))**2,.not.param%is_freeze_zone(j,:,:,:))
                enddo

                do j=1,param%n2
                    self%xnorms(i) = self%xnorms(i) + xnorms_weights(2)*sum((x(:,j+1,:,:)-x(:,j,:,:))**2,.not.param%is_freeze_zone(:,j,:,:))
                enddo

                if(m%is_cubic) then
                do j=1,param%n3
                    self%xnorms(i) = self%xnorms(i) + xnorms_weights(3)*sum((x(:,:,j+1,:)-x(:,:,j,:))**2,.not.param%is_freeze_zone(:,:,j,:))
                enddo
                endif

                !spatial gradient
                do i4=1,param%npars
                do i3=1,param%n3
                do i2=1,param%n2
                do i1=1,param%n1
                    spatgrad(i1,i2,i3,i4) = spatgrad(i1,i2,i3,i4) &
                        -xnorms_weights(1)*(x(i1-1,i2,i3,i4)-2.*x(i1,i2,i3,i4)+x(i1+1,i2,i3,i4)) &
                        -xnorms_weights(2)*(x(i1,i2-1,i3,i4)-2.*x(i1,i2,i3,i4)+x(i1,i2+1,i3,i4)) &
                        -xnorms_weights(3)*(x(i1,i2,i3-1,i4)-2.*x(i1,i2,i3,i4)+x(i1,i2,i3+1,i4))
                enddo
                enddo
                enddo
                enddo

            case ('TV')

            end select

        enddo

        self%xnorms = 0.5*self%xnorms*param%cell_volume_in_Pa

        where(param%is_freeze_zone) spatgrad=0.

        !add to dnorm's gradient wrt parameters
        g = g + spatgrad*param%cell_volume_in_Pa

    end subroutine

    real function total_loss(self)
        class(t_fobjective) :: self

        total_loss = self%dnorms(self%i_dnorm_4adjsource) + sum(self%xnorms)

    end function

    subroutine print_dnorms(self)
        class(t_fobjective) :: self
    end subroutine

    subroutine print_xnorms(self)
        class(t_fobjective) :: self
    end subroutine

end