module m_correlate
use m_System
use m_model
use m_computebox

    private

    public :: correlate_init, correlate_stack

    real,dimension(:,:,:,:),allocatable,public :: correlate_energy, correlate_image, correlate_gradient

    !correlate
    type,public :: t_correlate

        character(:),allocatable :: name,purpose

        ! !dimension
        ! integer :: nz,nx,ny,nt

        ! !index
        ! integer :: ifz,ilz,ifx,ilx,ify,ily,ift,ilt

        !wavefield components in computation domain
        real,dimension(:,:,:),allocatable :: sp_sp
        real,dimension(:,:,:),allocatable :: rp_sp
        real,dimension(:,:,:),allocatable :: drp_dt_dsp_dt, div_rp_div_sp
        real,dimension(:,:,:),allocatable :: rp_lap_sp, lap_rp_sp

        contains

        ! procedure :: init
        ! procedure :: init_bounds
        ! procedure :: check_value
        ! procedure :: change_dim
        procedure :: scale
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

    ! subroutine init_bounds(self,name,nz,nx,ny,nt,o_init)
    !     class(t_correlate) :: self
    !     character(*) :: name
    !     integer,dimension(2) :: nz,nx,ny,nt
    !     real,optional :: o_init

    !     self%name=name

    !     call alloc(self%sp_rp,nz,nx,ny,nt,o_init=o_init)
    
    ! 	self%ifz=nz(1); self%ilz=nz(2); self%nz=nz(2)-nz(1)+1
    ! 	self%ifx=nx(1); self%ilx=nx(2); self%nx=nx(2)-nx(1)+1
    ! 	self%ify=ny(1); self%ily=ny(2); self%ny=ny(2)-ny(1)+1
    ! 	self%ift=nt(1); self%ilt=nt(2); self%nt=nt(2)-nt(1)+1

    ! end subroutine

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

    ! subroutine init(self,name,shape123_from,o_comp)
    !     class(t_correlate) :: self
    !     character(*) :: name, shape123_from
    !     character(*),optional :: o_comp

    !     type(t_string),dimension(:),allocatable :: comp

    !     self%name=name

    !     select case (shape123_from)
    !         case ('model')
    !         self%ifz=1; self%ilz=m%nz; self%nz=m%nz
    !         self%ifx=1; self%ilx=m%nx; self%nx=m%nx
    !         self%ify=1; self%ily=m%ny; self%ny=m%ny

    !         case ('computebox')
    !         self%ifz=cb%ifz; self%ilz=cb%ilz; self%nz=cb%nz
    !         self%ifx=cb%ifx; self%ilx=cb%ilx; self%nx=cb%nx
    !         self%ify=cb%ify; self%ily=cb%ily; self%ny=cb%ny

    !     end select

    !     if(present(o_comp)) then
    !         comp=split(o_comp,o_sep=',')
    !         do i=1,size(comp)
    !             select case (comp(i)%s)
    !             case ('rp_sp')
    !                 call alloc(self%rp_sp,        [self%ifz,self%ilz],[self%ifx,self%ilx],[self%ify,self%ily])
    !             case ('drp_dt_dsp_dt')
    !                 call alloc(self%drp_dt_dsp_dt,[self%ifz,self%ilz],[self%ifx,self%ilx],[self%ify,self%ily])
    !             case ('nab_rp_nab_sp')
    !                 call alloc(self%nab_rp_nab_sp,[self%ifz,self%ilz],[self%ifx,self%ilx],[self%ify,self%ily])
    !             case ('rp_lap_sp')
    !                 call alloc(self%rp_lap_sp,    [self%ifz,self%ilz],[self%ifx,self%ilx],[self%ify,self%ily])

    !             end select

    !         enddo
    !     endif

    ! end subroutine

    subroutine scale(self,scaler)
        class(t_correlate) :: self

        if(allocated(self%sp_sp)) call scale_copy(self%sp_sp,scaler)
        if(allocated(self%rp_sp)) call scale_copy(self%rp_sp,scaler)

        if(allocated(self%drp_dt_dsp_dt)) call scale_copy(self%drp_dt_dsp_dt,scaler)
        if(allocated(self%div_rp_div_sp)) call scale_copy(self%div_rp_div_sp,scaler)
        if(allocated(self%rp_lap_sp))     call scale_copy(self%rp_lap_sp,    scaler)
        if(allocated(self%lap_rp_sp))     call scale_copy(self%lap_rp_sp,    scaler)

        ! if(allocated(self%rp_sp)) then
        !     self%rp_sp        = self%rp_sp*scaler
        !     self%rp_sp(1,:,:) = self%rp_sp(2,:,:)
        ! endif
        ! if(allocated(self%drp_dt_dsp_dt)) then
        !     self%drp_dt_dsp_dt        = self%drp_dt_dsp_dt*scaler
        !     self%drp_dt_dsp_dt(1,:,:) = self%drp_dt_dsp_dt(2,:,:)
        ! endif
        ! if(allocated(self%nab_rp_nab_sp)) then
        !     self%nab_rp_nab_sp        = self%nab_rp_nab_sp*scaler
        !     self%nab_rp_nab_sp(1,:,:) = self%nab_rp_nab_sp(2,:,:)
        ! endif
        ! if(allocated(self%rp_lap_sp)) then
        !     self%rp_lap_sp        = self%rp_lap_sp*scaler
        !     self%rp_lap_sp(1,:,:) = self%rp_lap_sp(2,:,:)
        ! endif

    end subroutine

    subroutine scale_copy(array,scaler)
        real,dimension(:,:,:) :: array
        array=array*scaler
        array(1,:,:)=array(2,:,:)
        array(:,1,:)=array(:,2,:)
    end subroutine


    subroutine correlate_stack(small,big)
        real,dimension(:,:,:) :: small, big

        big(cb%ioz:cb%ioz+cb%mz-1,&
            cb%iox:cb%iox+cb%mx-1,&
            cb%ioy:cb%ioy+cb%my-1  ) = &
        big(cb%ioz:cb%ioz+cb%mz-1,&
            cb%iox:cb%iox+cb%mx-1,&
            cb%ioy:cb%ioy+cb%my-1  ) + small

    end subroutine

    ! subroutine change_dim(self,new_nz,new_nx,new_ny,new_nt)
    ! 	class(t_correlate) :: self
    ! 	integer,dimension(2) :: new_nz,new_nx,new_ny,new_nt

    ! 	real,dimension(:,:,:,:),allocatable :: tmp

    ! 	call alloc(tmp,[self%ifz,self%ilz],[self%ifx,self%ilx],[self%ify,self%ily],[self%ift,self%ilt],o_init=self%sp_rp)

    ! 	call alloc(self%sp_rp,new_nz,new_nx,new_ny,new_nt)
    
    ! 	self%sp_rp(self%ifz:self%ilz,self%ifx:self%ilx,self%ify:self%ily,self%ift:self%ilt) = tmp

    ! 	deallocate(tmp)

    ! 	self%ifz=new_nz(1); self%ilz=new_nz(2); self%nz=new_nz(2)-new_nz(1)+1
    ! 	self%ifx=new_nx(1); self%ilx=new_nx(2); self%nx=new_nx(2)-new_nx(1)+1
    ! 	self%ify=new_ny(1); self%ily=new_ny(2); self%ny=new_ny(2)-new_ny(1)+1
    ! 	self%ift=new_nt(1); self%ilt=new_nt(2); self%nt=new_nt(2)-new_nt(1)+1

    ! end subroutine

    ! subroutine project_to(self,big)
    ! 	class(t_correlate) :: self
    ! 	type(t_correlate) :: big

    ! 	big%sp_rp(self%ifz:self%ilz,self%ifx:self%ilx,self%ify:self%ily,self%ift:self%ilt) = self%sp_rp

    ! end subroutine

    subroutine write(self,o_it,o_suffix)
        class(t_correlate) :: self
        integer,optional :: o_it
        character(*),optional :: o_suffix

        character(:),allocatable :: suf

        if(.not.present(o_it)) then !just write
            if(allocated(self%sp_sp))         call sysio_write(self%name//'%sp_sp'        ,self%sp_sp,        m%n)
            if(allocated(self%rp_sp))         call sysio_write(self%name//'%rp_sp'        ,self%rp_sp,        m%n)
            if(allocated(self%drp_dt_dsp_dt)) call sysio_write(self%name//'%drp_dt_dsp_dt',self%drp_dt_dsp_dt,m%n)
            if(allocated(self%div_rp_div_sp)) call sysio_write(self%name//'%div_rp_div_sp',self%div_rp_div_sp,m%n)
            if(allocated(self%rp_lap_sp))     call sysio_write(self%name//'%rp_lap_sp'    ,self%rp_lap_sp,    m%n)
            if(allocated(self%lap_rp_sp))     call sysio_write(self%name//'%lap_rp_sp'    ,self%lap_rp_sp,    m%n)

            return

        endif

        if(if_snapshot) then !write snapshots

            suf=either(o_suffix,'',present(o_suffix))

            if(o_it==1 .or. mod(o_it,i_snapshot)==0 .or. o_it==nt) then
                if(allocated(self%sp_sp))         call sysio_write('snap_'//self%name//'%sp_sp'//suf,        self%sp_sp,        m%n,o_mode='append')
                if(allocated(self%rp_sp))         call sysio_write('snap_'//self%name//'%rp_sp'//suf,        self%rp_sp,        m%n,o_mode='append')
                if(allocated(self%drp_dt_dsp_dt)) call sysio_write('snap_'//self%name//'%drp_dt_dsp_dt'//suf,self%drp_dt_dsp_dt,m%n,o_mode='append')
                if(allocated(self%div_rp_div_sp)) call sysio_write('snap_'//self%name//'%div_rp_div_sp'//suf,self%div_rp_div_sp,m%n,o_mode='append')
                if(allocated(self%rp_lap_sp))     call sysio_write('snap_'//self%name//'%rp_lap_sp'//suf,    self%rp_lap_sp,    m%n,o_mode='append')
                if(allocated(self%lap_rp_sp))     call sysio_write('snap_'//self%name//'%lap_rp_sp'//suf,    self%lap_rp_sp,    m%n,o_mode='append')

            endif

            ! if(it==1 .or. mod(it,i_snapshot)==0 .or. it==nt) then
            !     do i=1,size(snapshot)
            !         select case (snapshot(i)%s)
            !         case ('rp_sp')
            !             call sysio_write('snap_'//self%name//'%rp_sp'//suf,self%rp_sp,m%n,o_mode='append')
            !         case ('drp_dt_dsp_dt')
            !             call sysio_write('snap_'//self%name//'%drp_dt_dsp_dt'//suf,self%drp_dt_dsp_dt,m%n,o_mode='append')
            !         case ('nab_rp_nab_sp')
            !             call sysio_write('snap_'//self%name//'%nab_rp_nab_sp'//suf,self%nab_rp_nab_sp,m%n,o_mode='append')
            !         case ('rp_lapsp')
            !             call sysio_write('snap_'//self%name//'%rp_lapsp'//suf,self%rp_lapsp,m%n,o_mode='append')
            !         case default

            !         end select

            !     enddo

            ! endif

        endif

    end subroutine

    subroutine final(self)
        type(t_correlate) :: self

        call dealloc(self%sp_sp,self%rp_sp)
        call dealloc(self%drp_dt_dsp_dt)
        call dealloc(self%div_rp_div_sp)
        call dealloc(self%rp_lap_sp,self%lap_rp_sp)

    end subroutine

end