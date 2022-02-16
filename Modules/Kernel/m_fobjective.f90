module m_fobjective
use m_System
use m_Modeling
use m_Lpnorm
use m_matchfilter
use m_smoother_laplacian_sparse
use m_parametrizer
use m_querypoint
use m_weighter
use m_preconditioner

    private

    type,public :: t_fobjective

        type(t_string),dimension(:),allocatable :: s_dnorms !norms for data residuals
        type(t_string),dimension(:),allocatable :: s_xnorms !norms for (unknown) parameter residuals

        real,dimension(:),allocatable :: dnorms !unit in [Nm]
        real,dimension(:),allocatable :: xnorms !unit in [Nm]

        real FWI_misfit
        
        integer :: n_dnorms, n_xnorms

        real,dimension(:),allocatable :: dnorm_weights !scalar weights before each dnorm in the objective function
        real,dimension(:),allocatable :: xnorm_weights !scalar weights before each xnorm in the objective function

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
    ! real :: xnorms_weights(4) !weights on 1st, 2nd, 3rd & 4th dimensions

    contains

    subroutine init(self)
        class(t_fobjective) :: self
        
        !data norms
        self%s_dnorms=setup%get_strs('DATA_NORMS','DNORMS',o_default='L1 L2')
        self%n_dnorms=size(self%s_dnorms)
        call alloc(self%dnorms,self%n_dnorms)

        self%dnorm_weights=setup%get_reals('DATA_NORM_WEIGHTS','DWEIGHTS',o_default='0. 1.')

        !parameter norms
        self%s_xnorms=setup%get_strs('PARAMETER_NORMS','XNORMS',o_default='none')
        self%n_xnorms=size(self%s_xnorms)
        call alloc(self%xnorms,self%n_xnorms)

        self%xnorm_weights=setup%get_reals('PARAMETER_NORM_WEIGHTS','XWEIGHTS',o_default='0.')

        ! if(self%n_xnorms>0) then
        !     xnorms_weights=setup%get_reals('PARAMETER_NORMS_WEIGHTS','XNORMS_WEI',o_default='1. 1. 1. 1.')
        ! endif

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

        real,dimension(:,:),allocatable :: Wdres
        real :: numer, denom, time(shot%nt), v(shot%nt)
        logical :: is_adjsrc
        
        !Data residual of 1st shot as an example
        if(mpiworld%is_master) then
            Wdres=wei%weight*(shot%dsyn-shot%dobs)
            call shot%write('Wdres_',Wdres)
            deallocate(Wdres)
        endif

        self%dnorms=0.

        do i=1,self%n_dnorms
            select case (self%s_dnorms(i)%s)

            case ('L2')

                is_adjsrc=abs(self%dnorm_weights(i))>0.

                do ir=1,shot%nrcv
                    if(shot%rcv(ir)%is_badtrace) cycle

                    !unit of dnorms = [Nm], same as Lagrangian
                    !this also helps balance contributions from different component data
                    if(shot%rcv(ir)%comp=='p') then !for pressure data
                        self%dnorms(i) = self%dnorms(i) + L2norm_sq(shot%nt,wei%weight(:,ir),&
                            shot%dsyn(:,ir)-shot%dobs(:,ir),shot%dt,0.5/m%unit_pressure*m%unit_volume/m%unit_time)
                            !   [Pa]^2                      [s]      [1/Pa m^3/s]=[Nm]

                        if(is_adjsrc) call adjsrc_L2norm_sq(shot%dadj(:,ir))

                    else !for velocities data
                        self%dnorms(i) = self%dnorms(i) + L2norm_sq(shot%nt,wei%weight(:,ir),&
                            shot%dsyn(:,ir)-shot%dobs(:,ir),shot%dt,0.5*m%unit_mass/m%unit_time)
                            !   [m/s]^2                     [s]        [kg/s]=[kg m^2/s^2]=[Nm]

                        if(is_adjsrc) call adjsrc_L2norm_sq(shot%dadj(:,ir))

                    endif

                enddo

            case ('L1')

                is_adjsrc=abs(self%dnorm_weights(i))>0.

                do ir=1,shot%nrcv
                    if(shot%rcv(ir)%is_badtrace) cycle

                    !unit of dnorms = [Nm], same as Lagrangian
                    !this also helps balance contributions from different component data
                    if(shot%rcv(ir)%comp=='p') then !for pressure data
                        self%dnorms(i) = self%dnorms(i) + L1norm(shot%nt,wei%weight(:,ir),&
                            shot%dsyn(:,ir)-shot%dobs(:,ir),shot%dt,m%unit_volume/m%unit_time)
                            !   [Pa]                        [s]     [m^3/s]=[Nm]

                        if(is_adjsrc) call adjsrc_L1norm(shot%dadj(:,ir))

                    else !for velocities data
                        self%dnorms(i) = self%dnorms(i) + L1norm(shot%nt,wei%weight(:,ir),&
                            shot%dsyn(:,ir)-shot%dobs(:,ir),shot%dt,m%unit_mass*m%unit_length/m%unit_time/m%unit_time)
                            !   [m/s]                       [s]     [kg m/s^2]=[kg m^2/s^2]=[Nm]

                        if(is_adjsrc) call adjsrc_L1norm(shot%dadj(:,ir))

                    endif

                enddo

            case ('ndecon') !normalized deconvolution
                !aka Reverse AWI (Warner & Guasch, 2016, Geophysics, 81-6)
                !Eq D-4. Note this eq has a minor type; check the paper for correct one.

                ! if(is_first_in) then
                !     ref_p = 0.
                !     do ir=1,shot%nrcv
                !         if(shot%rcv(ir)%is_badtrace) cycle
                !         if(shot%rcv(ir)%comp=='p') then !for pressure data
                !             tmp=maxval(abs(Wdres(:,ir)))
                !             ref_p=either(ref_p,tmp,ref_p>tmp)
                !         endif
                !     enddo
                ! endif

                is_adjsrc=abs(self%dnorm_weights(i))>0.

                numer=0.
                denom=0.
                time=[(it-1.,it=1,shot%nt)]*shot%dt !why shot%nt?

                do ir=1,shot%nrcv
                    if(shot%rcv(ir)%is_badtrace) cycle
                    if(shot%rcv(ir)%comp=='p') then
                        call matchfilter_estimate(shot%dobs(:,ir),shot%dsyn(:,ir),shot%nt,1,o_filter_time=v)
                        numer=numer+norm2(time*v) 
                        denom=denom+norm2(v*v)

                        self%dnorms(i) = self%dnorms(i) + 0.5*numer/denom !/shot%dt*m%cell_volume*ref_p*ref_p/m%ref_kpa
                    
                        if(is_adjsrc) then
                            call matchfilter_adjointsrc(shot%dobs(:,ir),shot%dadj(:,ir), &
                                -(time*time-2.*self%dnorms(i))/denom*v)

                        endif
                
                    endif
                enddo               

            end select

        enddo

        !scale dnorms by shotlist
        call shls%scale(self%dnorms)

        call shot%write('dadj_',shot%dadj)

    end subroutine

    subroutine compute_xnorms(self,x,g)
        class(t_fobjective) :: self
        real,dimension(param%n1,param%n2,param%n3,param%npars) :: x,g

        real,dimension(:,:,:,:),allocatable :: xres !weighted residuals of x
        real,dimension(:,:,:,:),allocatable :: spatgrad !spatial gradient of x

        call alloc(spatgrad,param%n1,param%n2,param%n3,param%npars,o_init=0.)

        self%xnorms=0.

        do i=1,self%n_xnorms
            select case (self%s_xnorms(i)%s)
            case ('Ridge','L2')
                !fcost
                call alloc(xres,param%n1,param%n2,param%n3,param%npars)
                if(m%if_has_prior) then
                    xres=x-xprior
                else
                    xres=x
                endif

                self%xnorms(i) = self%xnorms(i) + sum(x*x)
                !self%xnorms(i) = sum(Wxres*Wxres, .not.param%is_freeze_zone)

                !spatial gradient
                spatgrad = spatgrad + xres

                deallocate(xres)

            case ('Lasso','L1')
                !to be implemented..

            case ('1st_Ridge','1-L2')

                do j=1,param%n1
                    !self%xnorms(i) = self%xnorms(i) + xnorms_weights(1)*sum((x(j+1,:,:,:)-x(j,:,:,:))**2,.not.param%is_freeze_zone(j,:,:,:))
                    self%xnorms(i) = self%xnorms(i) + sum((x(j+1,:,:,:)-x(j,:,:,:))**2)
                enddo

                do j=1,param%n2
                    self%xnorms(i) = self%xnorms(i) + sum((x(:,j+1,:,:)-x(:,j,:,:))**2)
                enddo

                if(m%is_cubic) then
                do j=1,param%n3
                    self%xnorms(i) = self%xnorms(i) + sum((x(:,:,j+1,:)-x(:,:,j,:))**2)
                enddo
                endif

                !spatial gradient
                do i4=1,param%npars
                do i3=1,param%n3
                do i2=1,param%n2
                do i1=1,param%n1
                    ! spatgrad(i1,i2,i3,i4) = spatgrad(i1,i2,i3,i4) &
                    !     -xnorms_weights(1)*(x(i1-1,i2,i3,i4)-2.*x(i1,i2,i3,i4)+x(i1+1,i2,i3,i4)) &
                    !     -xnorms_weights(2)*(x(i1,i2-1,i3,i4)-2.*x(i1,i2,i3,i4)+x(i1,i2+1,i3,i4)) &
                    !     -xnorms_weights(3)*(x(i1,i2,i3-1,i4)-2.*x(i1,i2,i3,i4)+x(i1,i2,i3+1,i4))
                    spatgrad(i1,i2,i3,i4) = spatgrad(i1,i2,i3,i4) &
                        -x(i1-1,i2,i3,i4)-2.*x(i1,i2,i3,i4)+x(i1+1,i2,i3,i4) &
                        -x(i1,i2-1,i3,i4)-2.*x(i1,i2,i3,i4)+x(i1,i2+1,i3,i4) &
                        -x(i1,i2,i3-1,i4)-2.*x(i1,i2,i3,i4)+x(i1,i2,i3+1,i4)
                enddo
                enddo
                enddo
                enddo

            case ('TV')
                !to be implemented

            end select

        enddo

        self%xnorms = 0.5*self%xnorms *m%cell_volume*m%unit_pressure
                                      ![m^3 Pa]=[Nm]

        !where(param%is_freeze_zone) spatgrad=0.

        !add to dnorm's gradient wrt parameters
        g = g + spatgrad *m%cell_volume*m%unit_pressure

    end subroutine

    subroutine print_dnorms(self,prefix,suffix)
        class(t_fobjective) :: self
        character(*) :: prefix,suffix

        if(mpiworld%is_master) then
            write(*,*) prefix//' fobj%dnorms '//suffix//': ', (' ',self%s_dnorms(i)%s,' ',self%dnorms(i),i=1,self%n_dnorms)
        endif

    end subroutine

    subroutine print_xnorms(self)
        class(t_fobjective) :: self

        if(mpiworld%is_master) then
            write(*,*) 'fobj%xnorms:', (' ',self%s_xnorms(i)%s,' ',self%xnorms(i),i=1,self%n_xnorms)
        endif

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
            call modeling_gradient(qp%is_fitting_data)
        ! endif
        
        qp%f = sum(self%dnorm_weights*self%dnorms) ! + sum(self%xnorm_weights*self%xnorms)

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

        !transform to x-space
        call param%transform(o_g=qp%g)

        !Regularization in x-space
        ! if(either(oif_approx,.false.,present(oif_approx))) then
        !     call regularize_approximate(fobj,qp)
        ! else
        !    call regularize(self,qp)
        ! endif

        !preconditioner
        call preco%update
        call preco%apply(qp%g,qp%pg)

        !save some RAM
        call dealloc(m%gradient,m%energy)

    end subroutine

    subroutine regularize(fobj,qp)
        type(t_fobjective) :: fobj
        type(t_querypoint) :: qp

        if(mpiworld%is_master) call fobj%compute_xnorms(qp%x,qp%g)

    end subroutine

end
