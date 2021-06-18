module m_fobjective
use m_string
use m_arrayop
use m_setup
use m_format_su
use m_model
use m_shot
use m_weighter

    private

    type,public :: t_fobjective
        type(t_string),dimension(:),allocatable :: s_dnorms !norms for data residuals
        type(t_string),dimension(:),allocatable :: s_mnorms !norms for model residuals

        real,dimension(:),allocatable :: dnorms
        real,dimension(:),allocatable :: mnorms

        integer :: n_dnorms, n_mnorms

        real,dimension(:),allocatable :: mpenalizers

        integer :: i_dnorm_4adjsource !which dnorm is used for the adjoint source (only 1 is allowed)

        real :: total_misfit

        real,dimension(:,:,:,:),allocatable :: gradient

        contains

        procedure :: init => init
        procedure :: compute_dnorms => compute_dnorms
        procedure :: compute_mnorms => compute_mnorms
        procedure :: compute_adjsource => compute_adjsource
        procedure :: print_dnorms => print_dnorms
        procedure :: print_mnorms => print_mnorms

    end type

    type(t_fobjective),public :: fobj

    real,dimension(:,:),allocatable :: dres

    contains

    subroutine init(self)
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

        !model norms
        self%s_mnorms=setup%get_strs('MODEL_NORMS','MNORMS',o_default='none')
        self%n_mnorms=size(self%s_mnorms)
        call alloc(self%mnorms,self%n_mnorms)
        call alloc(self%mpenalizers,self%n_mnorms)

        self%total_misfit=0.

    end subroutine
    
    subroutine compute_dnorms(self)
        class(t_fobjective) :: self
        
        call alloc(dres,shot%nt,shot%nrcv)
        self%dnorms=0.

        do i=1,self%n_dnorms
            select case (self%s_dnorms(i)%s)
            case ('L2')
                dres = (shot%dsyn-shot%dobs)*wei%weight
                
                do ir=1,shot%nrcv
                    if(shot%rcv(ir)%if_badtrace) cycle

                    self%dnorms(i) = self%dnorms(i) + sum(dres(:,ir)**2)

                    !set unit of dnorm to be [Nm], same as Lagrangian
                    !this also help balance contributions from different component data
                    if(shot%rcv(ir)%comp=='p') then !for pressure data
                        self%dnorms(i) = self%dnorms(i) *0.5/m%ref_kpa*m%cell_volume
                    else !for velocities data
                        self%dnorms(i) = self%dnorms(i) *0.5*m%ref_rho*m%cell_volume
                    endif

                enddo
                
            case ('L1')
                dres = (shot%dsyn-shot%dobs)*wei%weight

                do ir=1,shot%nrcv
                    if(shot%rcv(ir)%if_badtrace) cycle
                    
                    self%dnorms(i) = self%dnorms(i) + sum(abs(dres(:,ir)))

                    if(shot%rcv(ir)%comp=='p') then !for pressure data
                        self%dnorms(i) = self%dnorms(i) *0.5/m%ref_kpa*m%cell_volume
                    else !for velocities data
                        self%dnorms(i) = self%dnorms(i) *0.5*m%ref_rho*m%cell_volume
                    endif

                enddo

            end select

        enddo

    end subroutine

    subroutine compute_adjsource(self)
        class(t_fobjective) :: self

        !compute adjoint source and set proper units
        select case (self%s_dnorms(self%i_dnorm_4adjsource)%s)
        case ('L2')
            dres = (shot%dsyn-shot%dobs)*wei%weight
            shot%dadj = - dres*wei%weight
            do ir=1,shot%nrcv
                if(shot%rcv(ir)%comp=='p') then !for pressure data
                    shot%dadj(:,ir) = shot%dadj(:,ir) /m%ref_kpa/shot%dt
                else !for velocities data
                    shot%dadj(:,ir) = shot%dadj(:,ir) *m%ref_rho/shot%dt
                endif
            enddo

        end select

    end subroutine

    subroutine compute_mnorms(self)
        class(t_fobjective) :: self
    end subroutine

    subroutine compute_total_misfit(self)
        class(t_fobjective) :: self

        self%total_misfit = self%total_misfit +self%dnorms(self%i_dnorm_4adjsource)

    end subroutine

    subroutine print_dnorms(self)
        class(t_fobjective) :: self
    end subroutine

    subroutine print_mnorms(self)
        class(t_fobjective) :: self
    end subroutine

end
