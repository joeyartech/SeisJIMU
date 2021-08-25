module m_fobjective
use m_System
use m_Modeling
use m_smoother_laplacian_sparse
use m_parametrizer
use m_querypoint
use m_weighter
use m_preconditioner

    private

    type,public :: t_fobjective

        type(t_string),dimension(:),allocatable :: s_dnorms !norms for data residuals
        type(t_string),dimension(:),allocatable :: s_xnorms !norms for model residuals

        real,dimension(:),allocatable :: dnorms !unit in [Nm]
        real,dimension(:),allocatable :: xnorms !unit in [Nm]

        integer :: n_dnorms, n_xnorms
        integer :: i_dnorm_4adjsource !which dnorm is used for the adjoint source (only 1 is allowed so far)

        real,dimension(:),allocatable :: xpenalizers

        contains

        procedure :: init
        procedure :: compute_dnorms
        procedure :: compute_xnorms
        procedure :: print_dnorms
        procedure :: print_xnorms
        procedure :: eval

    end type

    type(t_fobjective),public :: fobj

    real,dimension(:,:),allocatable :: Wdres
    real,dimension(:,:,:,:),allocatable :: Wxres
    real,dimension(:,:,:,:),allocatable :: xprior !prior parameter
    real :: xnorms_weights(4) !weights on 1st, 2nd, 3rd & 4th dimensions

    contains

    subroutine init(self)
        class(t_fobjective) :: self
        
        !data norms
        self%s_dnorms=setup%get_strs('DATA_NORMS','DNORMS',o_default='L1 =>L2')
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

        !prior models
        call m%read_prior

        if(m%if_has_prior) then
            call alloc(xprior,param%n1,param%n2,param%n3,param%npars)
            !transform prior models
            call param%transform('m->x',o_xprior=xprior)
            !save some RAM
            call dealloc(m%vp_prior,m%vs_prior,m%rho_prior)
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

                if(mpiworld%is_master) call shot%write('Wdres_',Wdres)

                do ir=1,shot%nrcv
                    if(shot%rcv(ir)%is_badtrace) cycle

                    !set unit of dnorm to be [Nm], same as Lagrangian
                    !this also help balance contributions from different component data
                    if(shot%rcv(ir)%comp=='p') then !for pressure data
                        self%dnorms(i) = self%dnorms(i) + 0.5*m%cell_volume*norm2(Wdres(:,ir)) /m%ref_kpa
                    else !for velocities data
                        self%dnorms(i) = self%dnorms(i) + 0.5*m%cell_volume*norm2(Wdres(:,ir)) *m%ref_rho
                    endif

                enddo

                !compute adjoint source and set proper units
                if(i==self%i_dnorm_4adjsource) then
                    shot%dadj = -wei%weight*Wdres*m%cell_volume/shot%dt

                    do ir=1,shot%nrcv
                        if(shot%rcv(ir)%comp=='p') then !for pressure data
                            shot%dadj(:,ir) = shot%dadj(:,ir) /m%ref_kpa
                        else !for velocities data
                            shot%dadj(:,ir) = shot%dadj(:,ir) *m%ref_rho
                        endif

                    enddo

                    call shot%write('dadj_',shot%dadj)

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

    subroutine print_dnorms(self)
        class(t_fobjective) :: self
    end subroutine

    subroutine print_xnorms(self)
        class(t_fobjective) :: self
    end subroutine

    subroutine eval(self,qp,oif_update_m,oif_approx,oif_gradient)
        class(t_fobjective) :: self
        type(t_querypoint) :: qp
        logical,optional :: oif_update_m,oif_approx,oif_gradient

        type(t_string),dimension(:),allocatable :: smoothings
        character(:),allocatable :: smask
        real,dimension(:,:,:),allocatable :: mask

        !update model
        if(either(oif_update_m,.true.,present(oif_update_m))) then
            call param%transform(o_dir='x->m',o_x=qp%x)
        endif
        
        !compute dnorm's gradient by adjoint-state method
        ! if(either(oif_approx,.false.,present(oif_approx))) then
        !     call modeling_gradient_approximate(fobj)
        ! else
            call modeling_gradient(oif_gradient)
        ! endif

        qp%f = self%dnorms(self%i_dnorm_4adjsource) + sum(self%xnorms)

        if(.not. either(oif_gradient,.true.,present(oif_gradient))) return

        !post-modeling smoothing in physical (model) domain
        smoothings=setup%get_strs('SMOOTHING','SMTH',o_default='Laplacian')

        do i=1,size(smoothings)
            !Laplacian smoothing
            if(smoothings(i)%s=='Laplacian') then
                call hud('Laplacian smoothing')
                call smoother_laplacian_init([m%nz,m%nx,m%ny],[m%dz,m%dx,m%dy],shot%fpeak)
                do j=1,ppg%ngrad
                    !call smoother_laplacian_extend_mirror(m%gradient(:,:,:,i),m%itopo)
                    call smoother_laplacian_pseudo_nonstationary(m%gradient(:,:,:,j),m%vp)
                enddo    
            endif
        enddo

        !freeze_zone as hard mask
        do i=1,ppg%ngrad
            where(m%is_freeze_zone) m%gradient(:,:,:,i)=0.
        enddo

        !soft mask
        smask=setup%get_file('GRADIENT_SOFT_MASK','MASK')
        if(smask/='') then
            call alloc(mask,m%nz,m%nx,m%ny)
            ! call sysio_read(smask,mask,size(mask))
            
            do i=1,ppg%ngrad
                m%gradient(:,:,:,i)=m%gradient(:,:,:,i)*mask
            enddo
        endif

        ! call sysio_write('grho',m%gradient(:,:,:,1),size(m%gradient(:,:,:,1)))
        ! call sysio_write('gkpa',m%gradient(:,:,:,2),size(m%gradient(:,:,:,2)))

        !preconditioning
        call alloc(m%pgradient,m%nz,m%nx,m%ny,ppg%ngrad)
        call preco%update
        call preco%apply(m%gradient,m%pgradient)

        !transform
        call alloc(qp%g ,param%n1,param%n2,param%n3,param%npars)
        call alloc(qp%pg,param%n1,param%n2,param%n3,param%npars)
        call param%transform(o_g=qp%g,o_pg=qp%pg)
        !call param%transform_model('m->x',sfield%autocorr)

        !save some RAM
        call dealloc(m%gradient,m%pgradient)

        ! !Tikhonov regularization
        ! if(either(oif_approx,.false.,present(oif_approx))) then
        !     call regularize_approximate(fobj,qp)
        ! else
        !     call regularize(fobj,qp)
        ! endif

    end subroutine

    ! subroutine regularize(fobj,qp)
    !     type(t_fobjective) :: fobj
    !     type(t_querypoint) :: qp

    !     call fobj%compute_xnorms(qp%x,qp%g)

    ! end subroutine

end