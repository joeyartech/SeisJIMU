module m_mpienv
use omp_lib
use mpi
use m_string

    private

    !compilation info
    character(*),parameter :: commit = git_commit
    character(*),parameter :: branch = git_branch

#ifdef GNU
    character(*),parameter :: compiler = 'gfortran v' // __VERSION__
    character(*),parameter :: version = __VERSION__
    integer,parameter      :: endian = __BYTE_ORDER__
#endif

#ifdef INTEL
    character(*),parameter :: compiler = 'ifort'
    integer,parameter      :: version(2) = [__INTEL_COMPILER, __INTEL_COMPILER_UPDATE]
    character(*),parameter :: endian = 'unavailable'
#endif

    type,public :: t_mpienv
        integer iproc, nproc, communicator
        character(4) :: sproc
        logical :: is_master=.false.
        integer ierr
        character(MPI_MAX_PROCESSOR_NAME) :: node_name
        integer :: node_name_length
        
        integer :: thread_level
        integer :: max_threads

        contains
        procedure :: init   => init
        procedure :: barrier => barrier
        procedure :: write  => write
        procedure :: fin => fin
    end type
    
    type(t_mpienv),public :: mpiworld

    contains

    subroutine init(self,name,o_communicator,o_thread_level)
        class(t_mpienv) :: self
        character(*) :: name
        integer,optional :: o_communicator,o_thread_level
        
        self%communicator=either(o_communicator,MPI_COMM_WORLD,   present(o_communicator))
        
        call mpi_init_thread(either(o_thread_level,MPI_THREAD_SINGLE,present(o_thread_level)),self%thread_level,self%ierr)
        self%max_threads=OMP_GET_MAX_THREADS()

        call mpi_comm_rank(self%communicator,self%iproc,self%ierr)
        write(self%sproc,'(i0.4)') self%iproc
        if(self%iproc==0) self%is_master=.true.

        call mpi_comm_size(self%communicator,self%nproc,self%ierr)
        
        call mpi_get_processor_name(self%node_name, self%node_name_length, self%ierr)

        if(self%is_master) then
            write(*,*) name//' info:'
            write(*,*) ' MPI_INIT_THREAD level:', self%thread_level
            write(*,*) ' Number of MPI processors:',self%nproc
            write(*,*) ' Max number of OMP threads / processor:',self%max_threads
        endif

        !exe info
        if(self%is_master) call exe_info
        
        !time stamp
        if(self%is_master) print*, time_stamp()

        call self%barrier

    end subroutine

    subroutine barrier(self)
        class(t_mpienv) :: self
        call mpi_barrier(self%communicator,self%ierr)
    end subroutine

    subroutine write(self,filename,string)
        class(t_mpienv) :: self
        character(*) :: filename, string

        integer fhandle
        integer(kind=mpi_offset_kind) disp, offset
        character(:),allocatable :: stamp, str

        !time stamp
        stamp='========================'//s_NL// &
            ' MPI File IO Time Stamp '//s_NL// &
            '========================'//s_NL//time_stamp()

        if(self%is_master) then
            open(10,file=filename,access='append')
            write(10,'(a)') stamp
            close(10)
        endif

        call self%barrier
        
        ishift=len(stamp)
        
        !mpi write only + append
        str=string//s_NL
        disp=ishift+len(str)*self%iproc
        
        call mpi_file_open(self%communicator, filename, mpi_mode_append+mpi_mode_wronly, mpi_info_null, fhandle, self%ierr)
        call mpi_file_set_view(fhandle, disp, mpi_char, mpi_char, 'native', mpi_info_null, self%ierr)
        call mpi_file_write(fhandle, str, len(str), mpi_char, mpi_status_ignore, self%ierr)        
        call mpi_file_close(fhandle, self%ierr)

    end subroutine
    
    subroutine fin(self)
        class(t_mpienv) :: self

        if(self%is_master) print*, time_stamp()

        call self%barrier
        call mpi_abort(self%communicator,ierrcode,self%ierr)
        call mpi_finalize(self%ierr)
        
    end subroutine

    
    subroutine exe_info
        character(i_str_xlen) :: exe

        call getarg(0,exe)
        write(*,*) 'Working directory: (pwd)'
        call execute_command_line('pwd', wait=.true.)
        
        write(*,*) 'Using executable: (ls -l $exe)'
        call execute_command_line('ls -l '//trim(adjustl(exe)), wait=.true.)
        
        write(*,*) 'Git Commit: ', commit
        write(*,*) 'Git Branch: ', branch
        write(*,*) 'Compiler: ',   compiler
        write(*,*) 'Version: ',    version
        write(*,*) 'Endianness: ', endian

        !!allow unlimited stack memory
        !call execute_command_line('ulimit -s unlimited', wait=.true.)

    end subroutine

    function time_stamp() result(str)
        character(10)         :: date
        character(10)         :: time
        character(10)         :: zone
        character(:),allocatable :: str

        call date_and_time(time=time,date=date,zone=zone)
        str =   'System date: '//date(5:6)//'/'//date(7:8)//'/'//date(1:4)//s_NL// &
                'System time: '//time(1:2)//':'//time(3:4)//':'//time(5:6)//s_NL// &
                'System timezone: '//zone(1:3)//':'//zone(4:5)//s_NL// &
                '                        '

    end function


end
