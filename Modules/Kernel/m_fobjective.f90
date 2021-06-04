module m_fobjective
use m_shot
use m_matchfilter
use m_weighter_polygon
use m_weighter_table

    type t_fobj
        type(t_string),dimension(:),allocatable :: s_dnorms !norms for data residuals
        type(t_string),dimension(:),allocatable :: s_mnorms !norms for model residuals

        real,dimension(:),allocatable :: dnorms
        real,dimension(:),allocatable :: mnorms

        integer :: n_dnorms, n_mnorms

        real,dimension(:,:),allocatable :: dweight
        real,dimension(:),allocatable :: mpenalizer

        integer :: idnorm_for_adjsrc !which dnorm for the adjoint source (only 1 is allowed)

        real :: total_misfit

        real,dimension(:,:,:,:),allocatable :: gradient

        contains
        procedure :: compute_dnorms => compute_dnorms
        procedure :: compute_mnorms => compute_mnorms
        procedure :: print_dnorms => print_dnorms
        procedure :: print_mnorms => print_mnorms

    end type

    contains

    subroutine init(self)
        type(t_fobj) :: self
        
        !data norms
        self%s_dnorms=setup%get_strs('DATA_NORMS',default='L1 =>L2 Linf')
        self%n_dnorms=size(self%s_dnorms)
        call alloc(self%dnorms(self%n_dnorms))

        do i=1,self%n_dnorms
            if(self%s_dnorms(i)%s(1:2)=='=>') then
                self%idnorm_for_adjsrc=i
                self%s_dnorms(i)%s(1:2)=='  '; self%s_dnorms(i)%s=lalign(self%s_dnorms(i)%s) !remove the indicator '=>'
                exit
            endif
        enddo

        !model norms
        self%s_mnorms=setup%get_strs('MODEL_NORMS',default='none')
        call alloc(self%mnorms(size(self%_mnorms)))

        self%n_mnorms=size(self%s_mnorms)

        !data weights
        if(.not. allocated(weight)) then
                    call alloc(weight,shot%rcv(1)%nt,shot%nrcv,initialize=.false.);  weight=1
                    call build_weight_polygon(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,weight) !so far the weighting is for mono component data only
                    call build_weight_table(shot%rcv(1)%nt,shot%rcv(1)%dt,shot%nrcv,shot%rcv(:)%aoffset,weight) !so far the weighting is for mono component data only
        ! open(33,file='weight',access='stream')
        ! write(33) weight
        ! close(33)
                endif

    end subroutine

    subroutine compute_dnorms(shot)
        type(t_shot) :: shot        
        
        shot%dres = (shot%dsyn-shot%dobs)*shot%weight
        
        dnorm= 0.
        
        !set unit of dnorm to be [Nm], same as Lagrangian
        !this also help balance contributions from different component data
        do ir=1,shot%nrcv
            if(shot%rcv(ir)%icomp==1) then !for pressure data
                dnorm = dnorm + sum(dres(:,ir)*dres(:,ir))*0.5/ref_modulus*m%cell_volume
            else !for velocities data
                dnorm = dnorm + sum(dres(:,ir)*dres(:,ir))*0.5*m%ref_rho  *m%cell_volume
            endif
        enddo
        
        !compute adjoint source and set proper units
        dres = - dres*weight
        do ir=1,shot%nrcv
            if(shot%rcv(ir)%icomp==1) then !for pressure data
                dres(:,ir) = dres(:,ir) / ref_modulus/shot%rcv(1)%dt
            else !for velocities data
                dres(:,ir) = dres(:,ir) * m%ref_rho  /shot%rcv(1)%dt
            endif
        enddo

    end subroutine

    subroutine compute_total_misfit

        self%total_misfit=self%total_misfit+self%dnorms(self%idnorm_for_adjsrc)

    end subroutine

    subroutine compute_gradient

        !postprocessing

        !smoothing
        if(setup%get_bool('IF_SMOOTHING',default=.true.)) then
            call hud('Initialize Laplacian smoothing')
            call init_smoother_laplacian([m%nz,m%nx,m%ny],[m%dz,m%dx,m%dy],shot%fpeak)
            do icorr=1,ncorr
                call smoother_laplacian_extend_mirror(gradient(:,:,:,icorr),m%itopo)
                call smoother_laplacian_pseudo_nonstationary(gradient(:,:,:,icorr),m%vp)
            enddo
        endif
        
        !bathymetry as hard mask
        do iy=1,m%ny
        do ix=1,m%nx
            gradient(1:m%itopo(ix,iy)-1,ix,iy,:) =0.
        enddo
        enddo

        !soft mask
        if(setup%get_file('GRADIENT_SOFT_MASK','SOFT_MASK')/='') then
            call alloc(mask,m%nz,m%nx,m%ny)
            open(12,file='grad_mask',access='direct',recl=4*m%n,action='read',status='old')
            read(12,rec=1) mask
            close(12)
            
            do i=1,size(gradient,4)
                gradient(:,:,:,i)=gradient(:,:,:,i)*mask
            enddo
        endif

    end subroutine
end
