module m_fobjective
use m_string
use m_arrayop
use m_setup
use m_format_su
use m_model
use m_shot
use m_weighter
use m_parametrizer

    private

    type,public :: t_fobjective
        real,dimension(:,:,:,:),allocatable :: x, xprior !unitless
        real :: h3
        real :: total_loss !unit in [Nm]
        real,dimension(:,:,:,:),allocatable :: gradient !unit in [1/Nm]

        type(t_string),dimension(:),allocatable :: s_dnorms !norms for data residuals
        type(t_string),dimension(:),allocatable :: s_xnorms !norms for model residuals

        real,dimension(:),allocatable :: dnorms !unit in [Nm]
        real,dimension(:),allocatable :: xnorms !unit in [Nm]

        integer :: n_dnorms, n_xnorms
        integer :: i_dnorm_4adjsource !which dnorm is used for the adjoint source (only 1 is allowed so far)

        real,dimension(:),allocatable :: xpenalizers

        contains

        procedure :: init => init
        procedure :: update => update
        procedure :: compute_dnorms => compute_dnorms
        procedure :: compute_xnorms => compute_xnorms
        procedure :: print_dnorms => print_dnorms
        procedure :: print_xnorms => print_xnorms

    end type

    type(t_fobjective),public :: fobj

    real,dimension(:,:),allocatable :: Wdres
    real,dimension(:,:,:,:),allocatable :: Wxres
    real :: xnorms_weights(4) !weights on 1st, 2nd, 3rd & 4th dimensions

    contains

    subroutine init(self)
        class(t_fobjective) :: self
        
        call alloc(self%x,     param%n1,param%n2,param%n3,param%npars)
        call param%transform_model('m->x',self%x)

        !prior models
        call m%read_prior
        if(m%if_exist_prior) then
            !convert prior models
            call param%transform_model('m->x',self%xprior,oif_prior=.true.)
            !save some RAM
            call dealloc(m%vp,m%vs,m%rho)
        else
            call alloc(self%xprior,param%n1,param%n2,param%n3,param%npars)
        endif

        self%h3=param%d1*param%d2*param%d3

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

        self%total_loss=0.

    end subroutine
    
    subroutine update(self,x)
        class(t_fobjective) :: self
        real,dimension(param%n1,param%n2,param%n3,param%npars) :: x
        self%x=x
    end subroutine

    subroutine compute_dnorms(self)
        class(t_fobjective) :: self
        
        call alloc(Wdres,shot%nt,shot%nrcv)

        self%dnorms=0.

        do i=1,self%n_dnorms
            select case (self%s_dnorms(i)%s)
            case ('L2')
                Wdres = (shot%dsyn-shot%dobs)*wei%weight

                call write('Wdres')

                do ir=1,shot%nrcv
                    if(shot%rcv(ir)%if_badtrace) cycle

                    self%dnorms(i) = self%dnorms(i) + self%h3/shot%dt*0.5*sum(Wdres(:,ir)**2)

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
                    shot%dadj = -Wdres*wei%weight*param%h3/shot%dt

                    do ir=1,shot%nrcv
                        if(shot%rcv(ir)%comp=='p') then !for pressure data
                            shot%dadj(:,ir) = shot%dadj(:,ir) /m%ref_kpa
                        else !for velocities data
                            shot%dadj(:,ir) = shot%dadj(:,ir) *m%ref_rho
                        endif

                    enddo

                    call write('shot%dadj')

                endif
                
            case ('L1')
                !call Lp(dres,shot,wei,1)

            end select

        enddo

    end subroutine

    subroutine compute_xnorms(self)
        class(t_fobjective) :: self

        real,dimension(:,:,:,:),allocatable :: grad

        call alloc(grad,param%n1,param%n2,param%n3,param%npars,o_init=0.)

        self%xnorms=0.

        do i=1,self%n_xnorms
            select case (self%s_xnorms(i)%s)
            case ('Ridge','L2')
                !fcost
                call alloc(Wxres,param%n1,param%n2,param%n3,param%npars)
                Wxres=self%x-self%xprior

                self%xnorms(i) = sum(Wxres*Wxres, .not.param%is_freeze_zone)

                !gradient
                grad = grad + Wxres

            case ('Lasso','L1')

            case ('1st-ord_Ridge')

                !fcost
                do j=1,param%n1
                    self%xnorms(i) = self%xnorms(i) + xnorms_weights(1)*sum((self%x(j+1,:,:,:)-self%x(j,:,:,:))**2,.not.param%is_freeze_zone(j,:,:,:))
                enddo

                do j=1,param%n2
                    self%xnorms(i) = self%xnorms(i) + xnorms_weights(2)*sum((self%x(:,j+1,:,:)-self%x(:,j,:,:))**2,.not.param%is_freeze_zone(:,j,:,:))
                enddo

                if(m%is_cubic) then
                do j=1,param%n3
                    self%xnorms(i) = self%xnorms(i) + xnorms_weights(3)*sum((self%x(:,:,j+1,:)-self%x(:,:,j,:))**2,.not.param%is_freeze_zone(:,:,j,:))
                enddo
                endif

                !gradient
                do i4=1,param%npars
                do i3=1,param%n3
                do i2=1,param%n2
                do i1=1,param%n1
                    grad(i1,i2,i3,i4) = grad(i1,i2,i3,i4) &
                        -xnorms_weights(1)*(self%x(i1-1,i2,i3,i4)-2.*self%x(i1,i2,i3,i4)+self%x(i1+1,i2,i3,i4)) &
                        -xnorms_weights(2)*(self%x(i1,i2-1,i3,i4)-2.*self%x(i1,i2,i3,i4)+self%x(i1,i2+1,i3,i4)) &
                        -xnorms_weights(3)*(self%x(i1,i2,i3-1,i4)-2.*self%x(i1,i2,i3,i4)+self%x(i1,i2,i3+1,i4))
                enddo
                enddo
                enddo
                enddo

            case ('TV')

            end select

        enddo

        self%xnorms = 0.5*m%ref_kpa*self%xnorms*param%h3

        where(param%is_freeze_zone) grad=0.

        self%gradient = self%gradient + m%ref_kpa*param%h3*grad

    end subroutine

    subroutine compute_total_loss(self)
        class(t_fobjective) :: self

        self%total_loss = self%dnorms(self%i_dnorm_4adjsource) + sum(self%xnorms)

    end subroutine

    subroutine print_dnorms(self)
        class(t_fobjective) :: self
    end subroutine

    subroutine print_xnorms(self)
        class(t_fobjective) :: self
    end subroutine

end