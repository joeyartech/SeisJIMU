module m_shotlist
use m_string
use m_mpienv
use m_message
use m_arrayop
use m_setup
use m_sysio

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
!batch size = number of lists. Since one shot will be randomly selected from each list
!batch size is also the number of shots for one objective function and gradient computation.
!offsets = first shot index in each list.
!
!Two methods to build lists:
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
!Procedure assign: Assign lists to MPI processors. Each process is assigned with 1 or more lists.
!In the first example, if we have 2 processors,
!Proc #0 (master):
!1st list: { 1  2  3}
!3rd list: { 5  6  7  8}
!5th list: {15 16 17 18}
!
!Proc #1:
!2nd list: { 3  4  5}
!4th list: { 9 10 11 12}
!
!Procedure yield: During shot loop, randomly select shots from each list per processor. E.g.
!1st iteration:
!   Proc #0: { 1  2  3}    -> 2
!   Proc #1: { 3  4  5}    -> 3
!
!2nd iteration:
!   Proc #0: { 5  6  7  8} -> 5
!   Proc #1: { 9 10 11 12} -> 12
!
!3rd iteration:
!   Proc #1: {15 16 17 18} -> 17
!
!which gives 5 shots (=batch size) for one objective function and gradient computation.

    private

    type,public :: t_shotlist
        character(:),allocatable :: wholelist
        integer :: nshot
        type(t_string),dimension(:),allocatable :: lists, lists_per_processor
        integer :: nlists, nlists_per_processor

        contains
        procedure :: read_from_setup => read_from_setup
        procedure :: read_from_data  => read_from_data
        procedure :: build => build
        procedure :: assign => assign
        procedure :: yield => yield
        procedure :: print => print
    end type

    type(t_shotlist),public :: shls
    
    contains

    subroutine read_from_setup(self)
        class(t_shotlist) :: self

        select case (setup%get_str('ACQUI_GEOMETRY',o_default='spread'))

        case ('spread_irregular')
            open(13,file=setup%get_file('FILE_SOURCE_POSITION','SPOS',o_mandatory=.true.),action='read')
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
            self%nshot=setup%get_int('NUMBER_SOURCE','NS',o_mandatory=.true.)
            
        end select        

        self%wholelist=strcat(nums2strs([(i,i=1,self%nshot)],self%nshot),self%nshot)

        call hud('Will compute '//num2str(self%nshot)//' synthetic shots.')

    end subroutine
    
    subroutine read_from_data(self)
        class(t_shotlist) :: self

        type(t_string),dimension(:),allocatable :: list, sublist
        integer,dimension(:),allocatable :: ishots
        character(:),allocatable :: file
        logical :: exist
        integer file_size
        
        list=setup%get_strs('SHOT_INDEX','ISHOT')
        
        if(size(list)==0) then !not given
            call hud('SHOT_INDEX is not given. Now count how many data files exist in the directory..')
            file=setup%get_file('FILE_DATA')

            i=0; exist=.true.
            do
                i=i+1
                inquire(file=file//num2str(i,'(i0.4)')//'.su', size=file_size, exist=exist)
                if(file_size==0) exist=.false.
                if(.not.exist) exit
            enddo

            self%nshot=i

            self%wholelist=strcat(nums2strs([(i,i=1,self%nshot)],self%nshot),self%nshot)
        
            call hud('Found '//num2str(self%nshot)//' sequential shots.')

        endif

        if(size(list)>0) then !given
            
            allocate(ishots(1))
            
            do i=1,size(list)

                sublist=split(list(i)%s,o_sep=':')

                if(size(sublist)==1) then !add/rm
                    if(sublist(1)%s(1:1)=='-'.or.sublist(1)%s(1:1)=='!') then !rm
                        call rm(ishots,size(ishots),str2int(sublist(1)%s(2:)))
                    else
                        call add(ishots,size(ishots),str2int(sublist(1)%s))
                    endif
                endif
                if(size(sublist)==2) then !first:last
                    do j=str2int(sublist(1)%s),str2int(sublist(2)%s)
                        call add(ishots,size(ishots),j)
                    enddo
                endif
                if(size(sublist)==3) then !first:increment:last
                    do j=str2int(sublist(1)%s),str2int(sublist(3)%s),str2int(sublist(2)%s)
                        call add(ishots,size(ishots),j)
                    enddo
                endif

            enddo

            !check sanity of each shotgather file
            do i=1,size(ishots)
                inquire(file=num2str(ishots(i))//'.su', size=file_size, exist=exist)
                if(file_size==0) exist=.false.
                
                if(.not.exist) then
                    call hud('Neglect invalid shot file: Shot #'//int2str(ishots(i)))
                    call rm(ishots,size(ishots),ishots(i))
                endif
            enddo

            self%wholelist=strcat(nums2strs(ishots,size(ishots)),size(ishots))

            self%nshot=size(ishots)

            call hud('Total number of shots: '//num2str(self%nshot)//'')
            
        endif

        deallocate(list, sublist, ishots)

        call hud('No. of processors: '//num2str(mpiworld%nproc)//s_return// &
                'No. of shots: '//num2str(self%nshot))

        if(self%nshot<mpiworld%nproc) call warn('Shot number < Processor number. Some processors will be always idle..')
                
    end subroutine

    subroutine build(self)
        class(t_shotlist) :: self
        
        character(:),allocatable :: file
        type(t_string),dimension(:),allocatable :: swholelist
        
        file=setup%get_file('FILE_SHOT_BATCH')
        
        if(file/='') then
            !read ascii to self%lists
            return
        endif
        
        swholelist=split(self%wholelist)
        
        if(NBATCH>self%nshot) then
            print*,'NBATCH>self%nshot! Reset NBATCH=NPROC'
            NBATCH=NPROC
        endif
        allocate(self%lists(NBATCH))
        self%nlists=NBATCH
        
        maxlength=ceiling(self%nshot*1./NBATCH)

        do k=1,NBATCH
            self%lists(k)%s='' !initialization
        enddo
        
        k=1
        do i=1,self%nshot
            self%lists(k)%s=self%lists(k)%s//' '//swholelist(i)%s
            if(mod(i,maxlength)==0) k=k+1
        enddo
        
        !handle empty lists
        do k=NBATCH/2,NBATCH
            if(self%lists(k)%s=='') then
                self%lists(k)%s=swholelist(self%nshot)%s
            endif
        enddo
        
    end subroutine

    subroutine assign(self)
        class(t_shotlist) :: self
                
        !count how many lists per processors
        self%nlists_per_processor=0
        k=1 !list index
        j=0 !modulo list#/processor#
        do i=1,NBATCH
            if (j==NPROC) then
                j=j-NPROC
                k=k+1
            endif
            if (IPROC==j) then !IPROC starts from 0
                self%nlists_per_processor=self%nlists_per_processor+1
            endif
            j=j+1
        enddo
        
        allocate(self%lists_per_processor(ceiling(self%nlists*1./NPROC)))
        
        !redo the loop
        k=1 !list index
        j=0 !modulo list#/processor#
        do i=1,NBATCH
            if (j==NPROC) then
                j=j-NPROC
                k=k+1
            endif
            if (IPROC==j) then
                self%lists_per_processor(k)=self%lists(i)
            endif
            j=j+1
        enddo
        
    end subroutine
    
    function yield(self,i)
        class(t_shotlist) :: self
        character(:),allocatable :: yield
        
        type(t_string),dimension(:),allocatable :: slists
        
        slists=split(self%lists_per_processor(i)%s)
        
        call random_number(r) !gives random real number 0<=r<1
        !convert to random integer from n to m: i=n+floor((m+1-n)*r) 
        
        yield=slists(1+floor(size(slists)*r))%s
        
    end function

    subroutine print(self)
        class(t_shotlist) :: self
                
        ! !message
        ! write(*,*) ' Proc# '//mpiworld%sproc//' has '//num2str(self%nlist_per_processor)//' assigned shots.'
        ! call hud('See file "shotlist" for details.')

        ! !write shotlist to disk
        ! call mpiworld%write(dir_out//'shotlist', 'Proc# '//mpiworld%sproc//' has '//num2str(self%nshot_per_processor)//' assigned shots:'// &
        !     strcat(ints2strs(self%list_per_processor,self%nshot_per_processor),self%nshot_per_processor))
        
    end subroutine

end
