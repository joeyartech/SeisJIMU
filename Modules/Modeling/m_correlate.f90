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
        real,dimension(:,:,:),allocatable :: gkpa, gikpa
        real,dimension(:,:,:),allocatable :: grho, gbuo

        !image components
        real,dimension(:,:,:),allocatable :: ipp
        real,dimension(:,:,:),allocatable :: ibksc, ifwsc !backward & forward scatters

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
    type(t_string),dimension(:),allocatable :: snapshot
    integer :: i_snapshot, n_snapshot

    contains

    subroutine correlate_init(nt_in,dt_in)
        logical,save :: is_first_in=.true.

        nt=nt_in
        dt=dt_in

        !snapshot
        snapshot=setup%get_strs('SNAPSHOT')
        if_snapshot=size(snapshot)>0 .and. mpiworld%is_master
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
            if(allocated(self%grho))   call sysio_write(self%name//'%grho'//suf  ,self%grho,  m%n)
            if(allocated(self%gbuo))   call sysio_write(self%name//'%gbuo'//suf  ,self%gbuo,  m%n)
            if(allocated(self%gkpa))   call sysio_write(self%name//'%gkpa'//suf  ,self%gkpa,  m%n)
            if(allocated(self%gikpa))  call sysio_write(self%name//'%gikpa'//suf ,self%gikpa, m%n)

            if(allocated(self%ipp))    call sysio_write(self%name//'%ipp'//suf   ,self%ipp,   m%n)
            if(allocated(self%ibksc))  call sysio_write(self%name//'%ibksc'//suf ,self%ibksc, m%n)
            if(allocated(self%ifwsc))  call sysio_write(self%name//'%ifwsc'//suf ,self%ifwsc, m%n)
            return

        endif

        if(if_snapshot) then !write snapshots

            if(o_it==1 .or. mod(o_it,i_snapshot)==0 .or. o_it==nt) then
                if(allocated(self%grho))  call sysio_write('snap_'//self%name//'%grho'//suf, self%grho, m%n,o_mode='append')
                if(allocated(self%gbuo))  call sysio_write('snap_'//self%name//'%gbuo'//suf, self%gbuo, m%n,o_mode='append')
                if(allocated(self%gkpa))  call sysio_write('snap_'//self%name//'%gkpa'//suf, self%gkpa, m%n,o_mode='append')
                if(allocated(self%gikpa)) call sysio_write('snap_'//self%name//'%gikpa'//suf,self%gikpa,m%n,o_mode='append')
  
                if(allocated(self%ipp))    call sysio_write('snap_'//self%name//'%ipp'//suf   ,self%ipp,   m%n,o_mode='append')
                if(allocated(self%ibksc))  call sysio_write('snap_'//self%name//'%ibksc'//suf ,self%ibksc, m%n,o_mode='append')
                if(allocated(self%ifwsc))  call sysio_write('snap_'//self%name//'%ifwsc'//suf ,self%ifwsc, m%n,o_mode='append')

            endif

        endif

    end subroutine

    subroutine final(self)
        type(t_correlate) :: self

        call dealloc(self%grho, self%gbuo)
        call dealloc(self%gkpa, self%gikpa)

        call dealloc(self%ipp,self%ibksc,self%ifwsc)

    end subroutine

end