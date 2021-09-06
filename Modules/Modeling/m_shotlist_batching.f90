module m_shotlist
use m_System

!Procedure read: find all shots with their indices. E.g.
!{1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 18} (16 shots in total)
!
!Procedure build: build lists with batch size and offsets.
!E.g. with batch size=5, offsets=[1 3 5 9 15] 
!1st list: { 1  2  3}
!2nd list: { 3  4  5}
!3rd list: { 5  6  7  8}
!4th list: { 9 10 11 12}
!5th list: {15 16 17 18}
!batch size = number of lists. Since one shot will be selected from each list
!batch size is also the number of shots for one objective function and gradient computation.
!offsets = first shot index in each list.
!
!Two methods:
!   a) Read batch size from setup (default to number of MPI processors). Offsets are set such that first list starts with first shot and no overlapping of shots between lists.
!      Three typical cases:
!      1. Stochastic gradient descent (batch size=1):
!         {1}
!         {2}
!         {3}
!         ...
!         {18}
!
!      2. Minibatch gradient descent (1<batch size<all shots). E.g.
!         When batch size=3:
!         { 1  2  3  4  5  6} (length=ceiling(16/3))
!         { 7  8  9 10 11 12}
!         {15 16 17 18}       (last list has fewer shots, implying each shot has higher probability to be chosen than shots in other lists)
!
!         When batch size=4:
!         { 1  2  3  4}
!         { 5  6  7  8}
!         { 9 10 11 12}
!         {15 16 17 18}
!
!         When batch size=5:
!         { 1  2  3  4}
!         { 5  6  7  8}
!         { 9 10 11 12}
!         {15 16 17 18}
!         {18}                (Shot #18 is always chosen)
! 
!         When batch size=6:
!         { 1  2  3}
!         { 4  5  6}
!         { 7  8  9}
!         {10 11 12}
!         {15 16 17}
!         {18}
! 
!         When batch size=7:
!         { 1  2  3}
!         { 4  5  6}
!         { 7  8  9}
!         {10 11 12}
!         {15 16 17}
!         {18}
!         {18}                (Shot #18 will be always computed at least twice)
!
!       3. Batch gradient descent (batch size=number of all shots)
!          {1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 18}
!          In this case, no random shot selections.
! 
!   b) User provides a shotlists file (via setup)
!      for arbitrary lengths, arbitrary offsets etc.
!
!Procedure sample: select one shot from each list. Two methods:
!   a) Cyclic sampling. For the very beginning example,
!      lists         -> 1st sampled shots:
!      { 1  2  3}    -> 1
!      { 3  4  5}    -> 3
!      { 5  6  7  8} -> 5
!      { 9 10 11 12} -> 9
!      {15 16 17 18} -> 15
!      which gives 5 shots (=batch size) for one objective function and gradient computation.
!      
!      lists         -> 2nd sampled shots:
!      { 1  2  3}    -> 2
!      { 3  4  5}    -> 3
!      { 5  6  7  8} -> 6
!      { 9 10 11 12} -> 9
!      {15 16 17 18} -> 16
!      ...
!      lists         -> 4th sampled shots:
!      { 1  2  3}    -> 1
!      { 3  4  5}    -> 3
!      { 5  6  7  8} -> 8
!      { 9 10 11 12} -> 12
!      {15 16 17 18} -> 18
!                  
!   b) Random sampling.
!      1st sampled shots:
!      { 1  2  3}    -> 2
!      { 3  4  5}    -> 3
!      { 5  6  7  8} -> 5
!      { 9 10 11 12} -> 12
!      {15 16 17 18} -> 17
!         
!      2nd sampled shots:
!      { 1  2  3}    -> 1
!      { 3  4  5}    -> 5
!      { 5  6  7  8} -> 5
!      { 9 10 11 12} -> 9
!      {15 16 17 18} -> 16
!
!Procedure assign: assign the sampled shots to MPI processors. Each process is assigned with 1 or more shots to each process.
!E.g. sampled_shots={1 3 5 9 15}, if we have 2 processors,
!Proc #0 (master) is assigned with {1 5 15}
!Proc #1 is assigned with {3 9}
!
!Procedure yield: return one shot from sampled_shots inside shot loop.
!
!Procedure scale: scale the data misfit computed by using a portion of shots to approximate the misfit by all shots

    private

    type,public :: t_shotlist
        character(:),allocatable :: all
        integer :: nshot

        type(t_string),dimension(:),allocatable :: lists
        integer :: nlists
        
        integer,dimension(:),allocatable :: icycle
        type(t_string),dimension(:),allocatable :: sampled_shots
        
        type(t_string),dimension(:),allocatable :: shots_per_processor
        integer :: nshots_per_processor

        contains
        procedure :: read_from_setup
        procedure :: read_from_data
        procedure :: build
        procedure :: sample
        procedure :: assign
        procedure :: yield
        procedure :: scale

        procedure :: is_registered
        procedure :: register
    end type

    type(t_shotlist),public :: shls

    character(:),allocatable :: sample_method
    
    contains

    subroutine read_from_setup(self)
        class(t_shotlist) :: self

        select case (setup%get_str('ACQUI_GEOMETRY',o_default='spread'))

        case ('spread_irregular')
            open(13,file=setup%get_file('FILE_SOURCE_POSITION','SPOS',o_mandatory=1),action='read')
            !count number of sources
            self%nshot=0
            do
                read (13,*,iostat=msg) z,x,y
                if(msg/=0) exit
                self%nshot=self%nshot+1
            enddo
            close(13)
            call hud('Will read '//num2str(self%nshot)//' source positions.')

        case default
            self%nshot=setup%get_int('NUMBER_SOURCE','NS',o_default='1')
            
        end select        

        self%all=strcat(nums2strs([(i,i=1,self%nshot)]))

        call hud('Will compute '//num2str(self%nshot)//' synthetic shots.')

    end subroutine
    
    subroutine read_from_data(self)
        class(t_shotlist) :: self

        type(t_string),dimension(:),allocatable :: list, sublist
        integer,dimension(:),allocatable :: ishots
        character(:),allocatable :: file
        logical :: exist
        integer file_size
        
        file=setup%get_str('FILE_DATA_PREFIX')
        list=setup%get_strs('SHOT_INDEX','ISHOT')
        
        if(size(list)==0) then !not given
            call hud('SHOT_INDEX is not given. Now count how many data files exist in the directory..')

            i=0; exist=.true.
            do
                inquire(file=file//num2str(i+1,'(i0.4)')//'.su', size=file_size, exist=exist)
                if(file_size==0) exist=.false.
                if(.not.exist) exit
                i=i+1
            enddo

            self%nshot=i

            call hud('Found '//num2str(self%nshot)//' sequential shots.')

            if(self%nshot==0) call error('No shots are found!')

            self%all=strcat(nums2strs([(i,i=1,self%nshot)]))

        endif

        if(size(list)>0) then !given
            
            allocate(ishots(0))
            
            do i=1,size(list)

                sublist=split(list(i)%s,o_sep=':')

                if(size(sublist)==1) then !add/rm
                    if(sublist(1)%s(1:1)=='-'.or.sublist(1)%s(1:1)=='!') then !rm
                        call rm(ishots,str2int(sublist(1)%s(2:)))
                    else
                        call add(ishots,str2int(sublist(1)%s))
                    endif
                endif
                if(size(sublist)==2) then !first:last
                    do j=str2int(sublist(1)%s),str2int(sublist(2)%s)
                        call add(ishots,j)
                    enddo
                endif
                if(size(sublist)==3) then !first:increment:last
                    do j=str2int(sublist(1)%s),str2int(sublist(3)%s),str2int(sublist(2)%s)
                        call add(ishots,j)
                    enddo
                endif

            enddo

            !check sanity of each shotgather file
            do i=1,size(ishots)
                inquire(file=file//num2str(ishots(i),'(i0.4)')//'.su', size=file_size, exist=exist)
                if(file_size==0) exist=.false.
                
                if(.not.exist) then
                    call hud('Neglect invalid shot file: Shot #'//int2str(ishots(i)))
                    call rm(ishots,ishots(i))
                endif
            enddo

            self%all=strcat(nums2strs(ishots))

            self%nshot=size(ishots)

            call hud('Total number of shots: '//num2str(self%nshot)//'')
            
        endif

        if(allocated(   list)) deallocate(   list)
        if(allocated(sublist)) deallocate(sublist)
        if(allocated( ishots)) deallocate( ishots)

        call hud('No. of processors: '//num2str(mpiworld%nproc)//s_NL// &
                'No. of shots: '//num2str(self%nshot))

        if(self%nshot<mpiworld%nproc) call warn('Shot number < Processor number. Some processors will be always idle..')
                
    end subroutine

    subroutine build(self,o_batchsize)
        class(t_shotlist) :: self
        integer,optional :: o_batchsize
        
        character(:),allocatable :: file
        type(t_string),dimension(:),allocatable :: sall
        
        file=setup%get_file('FILE_SHOT_BATCH')
        
        if(file/='') then
            !read ascii to self%lists
            return
        endif
        
        sall=split(self%all)

        self%nlists=either(o_batchsize,&
            setup%get_int('SHOT_BATCH_SIZE','NBATCH',o_default=num2str(mpiworld%nproc)),&
            present(o_batchsize))

        if(self%nlists>self%nshot) then
            call warn('shotlist%nlists > shotlist%nshot! Reset shotlist%nlists = shotlist%nshot')
            self%nlists=self%nshot
        endif

        allocate(self%lists(self%nlists))

        maxlength=ceiling(self%nshot*1./self%nlists)

        !initialization
        do k=1,self%nlists
            self%lists(k)%s=''
        enddo

        k=1
        do i=1,self%nshot
            self%lists(k)%s=self%lists(k)%s//' '//sall(i)%s
            if(mod(i,maxlength)==0) k=k+1
        enddo

        !some lists could be empty
        do k=1,self%nlists
            if(self%lists(k)%s=='') then
                self%lists(k)%s=sall(self%nshot)%s
            endif
        enddo

        !write
        call sysio_rm('shotlist')
        if(mpiworld%is_master) then
            open(10,file=dir_out//'shotlist')
            write(10,'(a)') 'Build lists:'
            do k=1,self%nlists
                write(10,'(a)') '{'//self%lists(k)%s//'}'
            enddo
            write(10,'(a)') ' '
            close(10)
        endif
        
    end subroutine

    subroutine sample(self)
        class(t_shotlist) :: self

        logical,save :: is_first_in=.true.
        type(t_string),dimension(:),allocatable :: sublist
        
        if(is_first_in) then
            sample_method=setup%get_str('SHOT_SAMPLE_METHOD',o_default='cyclic')
            allocate(self%icycle(self%nlists),source=1)
            allocate(self%sampled_shots(self%nlists))

            is_first_in=.false.

        endif

        if(sample_method=='cyclic') then
            do i=1,self%nlists
                sublist=split(self%lists(i)%s)
                self%sampled_shots(i)%s=sublist(self%icycle(i))%s
                
                self%icycle(i)=self%icycle(i)+1
                if(self%icycle(i)>size(sublist)) self%icycle(i)=1

            enddo
        endif

        if(sample_method=='random') then
            do i=1,self%nlists
                sublist=split(self%lists(i)%s)
                
                call random_number(r) !gives random real number 0<=r<1
                !convert to random integer from n to m: i=n+floor((m+1-n)*r)            
                self%sampled_shots(i)%s=sublist(1+floor(size(sublist)*r))%s
                
            enddo
        endif

        !write
        if(mpiworld%is_master) then
            open(10,file=dir_out//'shotlist',access='append')
            write(10,'(a)') 'Sampled shots:'
            do k=1,self%nlists
                write(10,'(a)') '{'//self%sampled_shots(k)%s//'}'
            enddo
            write(10,'(a)') ' '
            close(10)
        endif
        
    end subroutine
    
    subroutine assign(self)
        class(t_shotlist) :: self

        logical,save :: is_first_in=.true.
        
        if(is_first_in) then
            !count how many shots per processors
            self%nshots_per_processor=0
            k=1 !list index
            j=0 !modulo list#/processor#
            do i=1,self%nlists
                if (j==mpiworld%nproc) then
                    j=j-mpiworld%nproc
                    k=k+1
                endif
                if (mpiworld%iproc==j) then !iproc starts from 0
                    self%nshots_per_processor=self%nshots_per_processor+1
                endif
                j=j+1
            enddo
            
            allocate(self%shots_per_processor(self%nshots_per_processor))
            
            is_first_in=.false.

        endif

        k=1 !list index
        j=0 !modulo list#/processor#
        do i=1,self%nlists
            if (j==mpiworld%nproc) then
                j=j-mpiworld%nproc
                k=k+1
            endif
            if (mpiworld%iproc==j) then
                self%shots_per_processor(k)%s=self%sampled_shots(i)%s
            endif
            j=j+1
        enddo

        !message
        write(*,*) ' Proc# '//mpiworld%sproc//' has '//num2str(self%nshots_per_processor)//' assigned shots.'
        call hud('See file "shotlist" for details.')

        !write shotlist to disk
        call mpiworld%write('shotlist', 'Proc# '//mpiworld%sproc//' has '//num2str(self%nshots_per_processor)//' assigned shots:'// &
            strcat(self%shots_per_processor))

    end subroutine

    function yield(self,i)
        class(t_shotlist) :: self
        !character(:),allocatable :: yield
        integer :: yield
        
        yield=str2int(self%shots_per_processor(i)%s)

    end function

    subroutine scale(self,dnorms)
        class(t_shotlist) :: self
        real,dimension(:) :: dnorms

        logical,save :: is_first_in=.true.
        integer,save :: n=0, m=1

        if(is_first_in) then
            m=size(self%lists) !=sampled shots

            !total number of elements in lists
            !/= number of all provided shots, since some shots might be repeated
            do i=1,size(self%lists)
                n = n + size(split(self%lists(i)%s))
            enddo

            is_first_in=.false.
            
        endif

        dnorms = dnorms * n/m

    end subroutine


    logical function is_registered(self,chp,str)
        class(t_shotlist) :: self
        type(t_checkpoint) :: chp
        character(*) :: str

        type(t_string),dimension(:),allocatable :: list
        real,dimension(:),allocatable :: tmp

        list=split(str)

        do i=1,size(list)
            is_registered=chp%check('shotlist%'//list(i)%s)
            if(.not.is_registered) return
        enddo

        do i=1,size(list)
            select case (list(i)%s)
            case ('sampled_shots')
                call chp%open('shotlist%sampled_shots')
                call alloc(tmp,self%nlists)
                call chp%read(tmp)
                self%sampled_shots=nums2strs(tmp)
                call chp%close
            end select

        enddo

    end function

    subroutine register(self,chp,str)
        class(t_shotlist) :: self
        type(t_checkpoint) :: chp
        character(*) :: str

        type(t_string),dimension(:),allocatable :: list
        real,dimension(:),allocatable :: tmp

        list=split(str)

        do i=1,size(list)
            select case (list(i)%s)
            case ('sampled_shots')
                call chp%open('shotlist%sampled_shots')
                call alloc(tmp,self%nlists)
                tmp=strs2reals(self%sampled_shots)
                call chp%write(tmp)
                call chp%close
            end select

        enddo

    end subroutine

end
