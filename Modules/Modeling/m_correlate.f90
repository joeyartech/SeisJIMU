module m_correlate
use m_System
use m_model
use m_computebox

    private

    public :: correlate_init, correlate_assemble

    real,dimension(:,:,:,:),allocatable,public :: correlate_energy, correlate_image, correlate_gradient, correlate_pgradient

    !correlate
    type,public :: t_correlate

        character(:),allocatable :: name

        ! !dimension
        ! integer :: nz,nx,ny,nt

        ! !index
        ! integer :: ifz,ilz,ifx,ilx,ify,ily,ift,ilt

        !gradient components
        real,dimension(:,:,:),allocatable :: grho, gbuo
        real,dimension(:,:,:),allocatable :: gkpa, gikpa
        real,dimension(:,:,:),allocatable :: glda, gmu

        !image components
        real,dimension(:,:,:),allocatable :: ipp
        real,dimension(:,:,:),allocatable :: ibksc, ifwsc !backward & forward scatters

        !energy components
        real,dimension(:,:,:),allocatable :: epp

        contains

        ! procedure :: init
        ! procedure :: init_bounds
        ! procedure :: check_value
        ! procedure :: change_dim
        procedure :: scale
        procedure :: stack
        procedure :: write
        final :: final
        
        ! procedure :: is_registered
        ! procedure :: register

    end type

    !propagator's nt, dt
    integer :: nt
    real :: dt

    !snapshot
    logical :: if_snapshot
    integer :: i_snapshot, n_snapshot

    contains

    subroutine correlate_init(nt_in,dt_in)
        logical,save :: is_first_in=.true.

        nt=nt_in
        dt=dt_in

        !snapshot
        if_snapshot=setup%get_bool('IF_SNAPSHOT') .and. mpiworld%is_master
        if(if_snapshot) then
            n_snapshot=setup%get_int('REF_NUMBER_SNAPSHOT','NSNAPSHOT',o_default='50')
            if(n_snapshot==0) n_snapshot=50
            i_snapshot=ceiling(nt*1./n_snapshot)

            !rm existing snap files
            if(is_first_in) then
                call sysio_rm('snap*')
                is_first_in=.false.
            endif

        endif

    end subroutine

    subroutine correlate_assemble(small,big)
        real,dimension(:,:,:) :: small, big

        big(cb%ioz:cb%ioz+cb%mz-1,&
            cb%iox:cb%iox+cb%mx-1,&
            cb%ioy:cb%ioy+cb%my-1  ) = &
        big(cb%ioz:cb%ioz+cb%mz-1,&
            cb%iox:cb%iox+cb%mx-1,&
            cb%ioy:cb%ioy+cb%my-1  ) + small

    end subroutine

    subroutine scale(self,scaler)
        class(t_correlate) :: self

        if(allocated(self%gkpa))       call scale_copy(self%gkpa,scaler)
        if(allocated(self%gikpa))      call scale_copy(self%gikpa,scaler)
        if(allocated(self%grho))       call scale_copy(self%grho,scaler)
        if(allocated(self%gbuo))       call scale_copy(self%gbuo,scaler)
        
    end subroutine

    subroutine scale_copy(array,scaler)
        real,dimension(:,:,:) :: array
        array=array*scaler
        array(1,:,:)=array(2,:,:)
        array(:,1,:)=array(:,2,:)
    end subroutine

    subroutine stack(self)
        class(t_correlate) :: self

        if(allocated(self%gkpa))  call mpi_reduce(mpi_in_place, self%gkpa,  m%n, mpi_real, mpi_sum, 0, mpiworld%communicator, mpiworld%ierr)
        if(allocated(self%gikpa)) call mpi_reduce(mpi_in_place, self%gikpa, m%n, mpi_real, mpi_sum, 0, mpiworld%communicator, mpiworld%ierr)
        if(allocated(self%grho))  call mpi_reduce(mpi_in_place, self%grho,  m%n, mpi_real, mpi_sum, 0, mpiworld%communicator, mpiworld%ierr)
        if(allocated(self%gbuo))  call mpi_reduce(mpi_in_place, self%gbuo,  m%n, mpi_real, mpi_sum, 0, mpiworld%communicator, mpiworld%ierr)
        
    end subroutine

    subroutine write(self,o_it,o_suffix)
        class(t_correlate) :: self
        integer,optional :: o_it
        character(*),optional :: o_suffix

        character(:),allocatable :: suf

        suf=either(o_suffix,'',present(o_suffix))

        if(.not.present(o_it)) then !just write
            call write_snaps(self%name,suf, self%grho,'grho', self%gbuo,'gbuo')
            call write_snaps(self%name,suf, self%gkpa,'gkpa', self%gikpa,'gikpa')
            call write_snaps(self%name,suf, self%glda,'glda', self%gmu,'gmu')

            call write_snaps(self%name,suf, self%ipp,'ipp')
            ! call write_snaps(self%name,suf, self%ibksc,'ibksc', self%ifwsc,'ifwsc', o_mode='append')

        return; endif

        if(if_snapshot) then !write snapshots

            if(o_it==1 .or. mod(o_it,i_snapshot)==0 .or. o_it==nt) then
                call write_snaps(self%name,suf, self%grho,'grho', self%gbuo,'gbuo', o_mode='append')
                call write_snaps(self%name,suf, self%gkpa,'gkpa', self%gikpa,'gikpa', o_mode='append')
                call write_snaps(self%name,suf, self%glda,'glda', self%gmu,'gmu', o_mode='append')
                
                call write_snaps(self%name,suf, self%ipp,'ipp', o_mode='append')
                ! call write_snaps(self%name,suf, self%ibksc,'ibksc', self%ifwsc,'ifwsc', o_mode='append')
                
            endif

        endif

    end subroutine

    subroutine write_snaps(name,suf, a1,c1, o_a2,o_c2, o_a3,o_c3, o_mode)
        character(*) :: name,suf        
        real,dimension(:,:,:),allocatable :: a1, o_a2, o_a3
        character(*) :: c1, o_c2, o_c3
        optional :: o_a2, o_c2, o_a3, o_c3
        character(*),optional :: o_mode

        if(allocated(a1)) then
            call sysio_write('snap_'//name//'%'//  c1//suf,  a1,m%n,o_mode=o_mode)
        endif

        if(present(o_a2)) then; if(allocated(o_a2)) then
            call sysio_write('snap_'//name//'%'//o_c2//suf,o_a2,m%n,o_mode=o_mode)
        endif; endif
        if(present(o_a3)) then; if(allocated(o_a3)) then
            call sysio_write('snap_'//name//'%'//o_c3//suf,o_a3,m%n,o_mode=o_mode)
        endif; endif

    end subroutine

    subroutine final(self)
        type(t_correlate) :: self

        call dealloc(self%grho, self%gbuo)
        call dealloc(self%gkpa, self%gikpa)

        call dealloc(self%ipp,self%ibksc,self%ifwsc)

    end subroutine

end
