module m_mpienv
use omp_lib
use mpi
use m_global

    private :: exe_info, time_stamp
    private :: init, barrier, write, fin

    type t_mpienv
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
    
    type(t_mpienv) :: mpiworld

    contains

    subroutine init(self,name,communicator,thread_level)
        type(t_mpienv) :: self
        character(*) :: name
        integer communicator, thread_level
        
        self%communicator=communicator

        call mpi_init_thread(thread_level,self%thread_level,self%ierr)
        self%max_threads=OMP_GET_MAX_THREADS()

        call mpi_comm_rank(self%communicator,self%iproc,self%ierr)
        write(self%sproc,'(i0.4)') self%iproc
        if(self%iproc==0) self%is_master=.true.

        call mpi_comm_size(self%communicator,self%nproc,self%ierr)
        
        call mpi_get_processor_name(self%node_name, self%node_name_length, self%ierr)

        if(self%is_master) then
            write(*,*) name//' info:'
            write(*,'(a,i2,a,i2)') ' MPI_INIT_THREAD, required level:',thread_level,', provided level:', self%thread_level
            write(*,'(a,i5)') ' Number of MPI processors:',self%nproc
            write(*,'(a,i5)') ' Max number of OMP threads / processor:',self%max_threads
        endif

        !exe info
        if(self%is_master) call exe_info
        
        !time stamp
        if(self%is_master) print*, time_stamp()

        call self%barrier

    end subroutine

    subroutine barrier(self)
        type(t_mpienv) :: self
        call mpi_barrier(self%communicator,self%ierr)
    end subroutine

    subroutine write(self,filename,string)
        type(t_mpienv) :: self
        character(*) :: filename, string

        integer file_handler
        integer(kind=mpi_offset_kind) disp, offset
        character(:),allocatable :: str

        !time stamp
        str='========================'//s_return// &
            ' MPI File IO Time Stamp '//s_return// &
            '========================'//s_return//time_stamp()
        
        if(self%is_master) then
            open(10,file=filename)
            write(10,'(a)') str
            close(10)
        endif

        call self%barrier
        
        ishift=len(str)
        
        !mpi write only + append
        str=string//s_return
        disp=ishift+len(str)*self%iproc
        
        call mpi_file_open(self%communicator, filename, mpi_mode_append+mpi_mode_wronly, mpi_info_null, file_handler, self%ierr)
        call mpi_file_set_view(file_handler, disp, mpi_char, mpi_char, 'native', mpi_info_null, self%ierr)
        call mpi_file_write(file_handler, str, len(str), mpi_char, mpi_status_ignore, self%ierr)        
        call mpi_file_close(file_handler, self%ierr)

    end subroutine
    
    subroutine fin(self)
        type(t_mpienv) :: self

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
        
        write(*,*) 'Git Commit: ', s_commit
        write(*,*) 'Git Branch: ', s_branch
        write(*,*) 'Compiler: ',   s_compiler
        write(*,*) 'Version: ',    s_version
        write(*,*) 'Endianness: ', i_endian

        !!allow unlimited stack memory
        !call execute_command_line('ulimit -s unlimited', wait=.true.)

    end subroutine

    function time_stamp() result(str)
        character(10)         :: date
        character(10)         :: time
        character(10)         :: zone
        character(:),allocatable :: str

        call date_and_time(time=time,date=date,zone=zone)
        str =   'System date: '//date(5:6)//'/'//date(7:8)//'/'//date(1:4)//s_return// &
                'System time: '//time(1:2)//':'//time(3:4)//':'//time(5:6)//s_return// &
                'System timezone: '//zone(1:3)//':'//zone(4:5)//s_return// &
                '                        '

    end function


end
