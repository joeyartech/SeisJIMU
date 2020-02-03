module m_mpienv
use omp_lib
use mpi

    type t_mpienv
        integer iproc, nproc, communicator
        character(4) :: cproc
        logical :: is_master=.false.
        integer ierr
        character(MPI_MAX_PROCESSOR_NAME) :: node_name
        integer :: node_name_length
        
        integer :: thread_level
        integer :: max_threads
    end type
    
    type(t_mpienv) :: mpiworld

    contains
    
    subroutine init_mpiworld(o_thread_level)
        integer,optional :: o_thread_level
        integer thread_level
        
        character(256) :: exe
        character(10)  :: date
        character(10)  :: time
        character(10)  :: zone
        
        !mpiworld
        mpiworld%communicator=MPI_COMM_WORLD !by default
        if(present(o_thread_level)) then
            thread_level=o_thread_level
        else
            thread_level=MPI_THREAD_SINGLE
        endif
        call mpi_init_thread(thread_level,mpiworld%thread_level,mpiworld%ierr)
        mpiworld%max_threads=OMP_GET_MAX_THREADS()
        call mpi_comm_rank(mpiworld%communicator,mpiworld%iproc,mpiworld%ierr)
        write(mpiworld%cproc,'(i0.4)') mpiworld%iproc
        if(mpiworld%iproc==0) mpiworld%is_master=.true.
        call mpi_comm_size(mpiworld%communicator,mpiworld%nproc,mpiworld%ierr)
        call mpi_get_processor_name(mpiworld%node_name, mpiworld%node_name_length, mpiworld%ierr)
        
        if(mpiworld%is_master) then
            write(*,'(a,i2,a,i2)') ' MPI_INIT_THREAD, required level:',thread_level,', provided level:', mpiworld%thread_level
            write(*,'(a,i5)') ' Number of MPI processors:',mpiworld%nproc
            write(*,'(a,i5)') ' Max number of OMP threads / processor:',mpiworld%max_threads
        endif
        
        !exe info
        call getarg(0,exe)
        if(mpiworld%is_master) write(*,*) 'Working directory: (pwd)'
        if(mpiworld%is_master) call execute_command_line('pwd', wait=.true.)
        
        if(mpiworld%is_master) write(*,*) 'Using executable: (ls -l $exe)'
        if(mpiworld%is_master) call execute_command_line('ls -l '//trim(adjustl(exe)), wait=.true.)
        
!         !ulimit -s unlimited
!         call execute_command_line('ulimit -s unlimited', wait=.true.)
        
        !time info
        if(mpiworld%is_master) then
            call date_and_time(time=time,date=date,zone=zone)
            write(*,*) 'System date: '//date(5:6)//'/'//date(7:8)//'/'//date(1:4)
            write(*,*) 'System time: '//time(1:2)//':'//time(3:4)//':'//time(5:6)
            write(*,*) 'System timezone: '//zone(1:3)//':'//zone(4:5)
        endif
        
        call mpi_barrier(mpiworld%communicator,mpiworld%ierr)
        
    end subroutine

    subroutine mpiworld_file_write(filename,string)
        character(*) :: filename, string      

        character(10)  :: date
        character(10)  :: time
        character(10)  :: zone

        integer file_handler
        integer(kind=mpi_offset_kind) disp, offset
        character(:),allocatable :: str

        !time info
        call date_and_time(time=time,date=date,zone=zone)
        str =   '========================'//achar(10)// &
                ' MPI File IO Time Stamp '//achar(10)// &
                '========================'//achar(10)// &
                'System date: '//date(5:6)//'/'//date(7:8)//'/'//date(1:4)//achar(10)// &
                'System time: '//time(1:2)//':'//time(3:4)//':'//time(5:6)//achar(10)// &
                'System timezone: '//zone(1:3)//':'//zone(4:5)//achar(10)// &
                '                        '//achar(10)
        !https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/494946
        !char function: CHAR(10) is linefeed LF. CHAR(13) is carriage return CR. If you are a little paranoid, ACHAR(10) is better - this is a little more robust to the default character kind not being ascii.
        !The NEW_LINE standard intrinsic is even more robust. There's also the C_NEW_LINE and C_CARRIAGE_RETURN constants from the ISO_C_BINDING module for kind C_CHAR characters.
        
        if(mpiworld%is_master) then
            open(10,file=filename)
            write(10,'(a)') str
            close(10)
        endif
        
        call mpi_barrier(mpiworld%communicator,mpiworld%ierr)
        
        ishift=len(str)
        
        !mpi write only + append
        str=string//achar(10)
        disp=ishift+len(str)*mpiworld%iproc
        
        call mpi_file_open(mpiworld%communicator, filename, mpi_mode_append+mpi_mode_wronly, mpi_info_null, file_handler, mpiworld%ierr)
        call mpi_file_set_view(file_handler, disp, mpi_char, mpi_char, 'native', mpi_info_null, mpiworld%ierr)
        call mpi_file_write(file_handler, str, len(str), mpi_char, mpi_status_ignore, mpiworld%ierr)        
        call mpi_file_close(file_handler, mpiworld%ierr)

    end subroutine
    
    subroutine mpiworld_finalize
        character(10)  :: date
        character(10)  :: time
        character(10)  :: zone
        
        !time info
        if(mpiworld%is_master) then
            call date_and_time(time=time,date=date,zone=zone)
            write(*,*) 'System date: '//date(5:6)//'/'//date(7:8)//'/'//date(1:4)
            write(*,*) 'System time: '//time(1:2)//':'//time(3:4)//':'//time(5:6)
            write(*,*) 'System timezone: '//zone(1:3)//':'//zone(4:5)
        endif
        
        !finalize mpi
        if(mpiworld%is_master) write(*,*) 'Waiting other processes before finish ...'
        call mpi_barrier(mpiworld%communicator,mpiworld%ierr)
        call mpi_finalize(mpiworld%ierr)
        
    end subroutine
    
    !abort an mpienv
    subroutine mpienv_abort(mpienv)
        type(t_mpienv) :: mpienv
        call mpi_abort(mpienv%communicator,ierrcode,mpienv%ierr)
    end subroutine
    
!     subroutine mpi_time(flag)
!         integer, parameter              :: p = selected_real_kind(12)
!         real(kind=p)                    :: eoverhead, coverhead
!         real(kind=p), dimension(2)      :: etime, ctime
!         integer, dimension(8)           :: values
!         character(len=8), dimension(2)  :: date
!         character(len=10), dimension(2) :: time
!         character(len=5)                :: zone
! 
!         !... input dummy parameter
!         integer, intent(in) :: flag
! 
!         !... local variables
!         integer                                 :: rank, nb_procs, i, code
!         integer, allocatable, dimension(:)      :: all_rank
!         real(kind=p), allocatable, dimension(:) :: all_etime, all_ctime, all_ratio
!         real(kind=p)                            :: total_etime,total_ctime,total_ratio,&
!             max_etime, max_ctime, max_ratio, &
!             min_etime, min_ctime, min_ratio, &
!             avg_etime, avg_ctime, avg_ratio, &
!             dummy
!         character(len=128), dimension(8) :: lignes
!         character(len=128)               :: hline, start_date, final_date
!         character(len=2048)              :: fmt
! 
!         select case(flag)
!         case(0)
! 
!             !... compute clock overhead
!             eoverhead = mpi_wtime()
!             eoverhead = mpi_wtime() - eoverhead
!             call cpu_time(dummy)
!             call cpu_time(coverhead)
!             if (dummy < 0.0_p) &
!                 print *,"warning, mpi_time: cpu user time is not available on this machine."
!             coverhead = coverhead - dummy
!             call mpi_comm_rank(mpi_comm_world, rank, code)
!             !... start of timings on "date & time"
!             if ( rank == 0 ) &
!                 call date_and_time(date(1),time(1),zone,values)
!             !... start elapsed and cpu time counters
!             etime(1) = mpi_wtime()
!             call cpu_time(ctime(1))
! 
!         case(1)
!             !... final cpu and elapsed times
!             call cpu_time(ctime(2))
!             etime(2) = mpi_wtime() - etime(1) - eoverhead - coverhead
!             ctime(2) = ctime(2) - ctime(1) - coverhead
!             !... gather all times
!             call mpi_comm_rank(mpi_comm_world, rank, code)
!             call mpi_comm_size(mpi_comm_world, nb_procs, code)
!             !modif for intel compiler r.brossier 03/09
!             !   if ( rank == 0) allocate(all_etime(nb_procs), &
!             allocate(all_etime(nb_procs), &
!                 all_ctime(nb_procs), &
!                 all_ratio(nb_procs), &
!                 all_rank(nb_procs) )
!             call mpi_gather(etime(2), 1, mpi_double_precision,  &
!                 all_etime, 1, mpi_double_precision, &
!                 0, mpi_comm_world, code)
!             call mpi_gather(ctime(2), 1, mpi_double_precision,  &
!                 all_ctime, 1, mpi_double_precision, &
!                 0, mpi_comm_world, code)
!             if ( rank == 0) then
!             all_rank(:) = (/ (i,i=0,nb_procs-1) /)
! 
!             !... compute elapse user time
!             total_etime = sum(all_etime(:))
!             avg_etime   = total_etime/real(nb_procs,kind=p)
!             max_etime   = maxval(all_etime(:))
!             min_etime   = minval(all_etime(:))
!             if( min_etime <= 0.0_p ) then
!                 print *,"warning, mpi_time: measured elapsed user time seems to be too short"
!                 print *,"compared to the clock precision. timings could be erroneous."
!             end if
! 
!             !... compute cpu user time
!             total_ctime = sum(all_ctime(:))
!             avg_ctime   = total_ctime/real(nb_procs,kind=p)
!             max_ctime   = maxval(all_ctime(:))
!             min_ctime   = minval(all_ctime(:))
!             if( min_ctime <= 0.0_p ) then
!                 print *,"warning, mpi_time: measured cpu user time seems to be too short"
!                 print *,"compared to the clock precision. timings could be erroneous."
!             end if
! 
!             !... compute cpu/elapsed ratio
!             all_ratio(:) = all_ctime(:) / all_etime(:)
!             total_ratio  = sum(all_ratio(:))
!             avg_ratio    = total_ratio/real(nb_procs,kind=p)
!             max_ratio    = maxval(all_ratio(:))
!             min_ratio    = minval(all_ratio(:))
! 
!             !... end of timings on "date & time"
!             call date_and_time(date(2),time(2),zone,values)
! 
!             !... output format
!             hline    ='10x,13("-"),"|",18("-"),"|",14("-"),"|",18("-"),/,'
!             lignes(1)='(//,10x,"(c) January 2003, CNRS/IDRIS, France.",/,'
!             lignes(2)='10x,"MPI_Time (Release 3.4) Summary Report:",//,'
!             lignes(3)='10x,"Process Rank |"," Elapsed Time (s) |"," CPU Time (s) |"," Ratio CPU/Elapsed",/,'
!             lignes(4)='    (10x,i4,9(" "),"|",f12.3,6(" "),"|",f12.3,2(" "),"|",4(" "),f7.3,/),'
!             write(lignes(4)(1:4),'(i4)') nb_procs
!             lignes(5)='10x,"Total        |",f12.3,6(" "),"|",f12.3,2(" "),"|",4(" "),f7.3,/,'
!             lignes(6)='10x,"Minimum      |",f12.3,6(" "),"|",f12.3,2(" "),"|",4(" "),f7.3,/,'
!             lignes(7)='10x,"Maximum      |",f12.3,6(" "),"|",f12.3,2(" "),"|",4(" "),f7.3,/,'
!             lignes(8)='10x,"Average      |",f12.3,6(" "),"|",f12.3,2(" "),"|",4(" "),f7.3,/,'
!             start_date='/,10x,"MPI_Time started on ",2(a2,"/"),a4," at ",2(a2,":"),a2," met ",a3,":",a2," from GMT",/,'
!             final_date='10x,  "MPI_Time   ended on ",2(a2,"/"),a4," at ",2(a2,":"),a2," met ",a3,":",a2," from GMT",//)'
!             fmt=trim(lignes(1))//trim(lignes(2))//trim(lignes(3))//           &
!                     & trim(hline)//trim(lignes(4))//trim(hline)//trim(lignes(5))//  &
!                     & trim(hline)//trim(lignes(6))//trim(hline)//trim(lignes(7))//  &
!                     & trim(hline)//trim(lignes(8))//trim(hline)//trim(start_date)// &
!                     & trim(final_date)
!             write(*, trim(fmt)) &
!                     (all_rank(i),all_etime(i),all_ctime(i),all_ratio(i),i=1, nb_procs), &
!                     total_etime,  total_ctime,  total_ratio,  &
!                     min_etime,    min_ctime,    min_ratio,    &
!                     max_etime,    max_ctime,    max_ratio,    &
!                     avg_etime,    avg_ctime,    avg_ratio,    &
!                     date(1)(7:8), date(1)(5:6), date(1)(1:4), &
!                     time(1)(1:2), time(1)(3:4), time(1)(5:6), &
!                     zone(1:3),    zone(4:5),                  &
!                     date(2)(7:8), date(2)(5:6), date(2)(1:4), &
!                     time(2)(1:2), time(2)(3:4), time(2)(5:6), &
!                     zone(1:3),    zone(4:5)
!     !modif for intel compiler r.brossier 03/09
!             ! deallocate(all_etime, all_ctime, all_ratio, all_rank)
!             end if
!             deallocate(all_etime, all_ctime, all_ratio, all_rank)
!         case default
!             print *,"Error, MPI_Time: Invalid input parameter"
! 
!         end select
!     end subroutine mpi_time
! 
!     subroutine sub_mpi_abort()
!     !####################################################################
!     ! suborutine to kill the code through a mpi_abort call in case of 
!     ! problem detected on a single process
!     !####################################################################
! 
!     integer :: errorcode,ierror
!     write(*,*)'######################################################'
!     write(*,*)' '
!     write(*,*)'code is asked to be killed now by process ',mype
!     write(*,*)' '
!     write(*,*)'******************************************************'
!     write(*,*)' '
!     write(*,*)'a mpi_abort call is sent now to kill all process'
!     write(*,*)' '
!     write(*,*)'               bye bye'
!     write(*,*)' '
!     write(*,*)'######################################################'
!     !     call flushc()
! 
!     call mpi_abort(mpi_comm_world,errorcode,ierror)
! 
!   end subroutine sub_mpi_abort

end
