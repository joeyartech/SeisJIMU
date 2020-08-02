module m_shotlist
use m_string
use m_mpienv
use m_message
use m_arrayop

    private
    public shotlist_build, shotlist, shotlist_nshot_per_processor

    integer,dimension(:),allocatable :: shotlist   

    integer :: nshot, shotlist_nshot_per_processor
        
    contains
    
    subroutine shotlist_build(o_nshot)
        integer,optional :: o_nshot

        integer :: fshot, dshot, lshot
        integer :: file_size
        logical :: alive
        character(:),allocatable :: data_file
        type(t_string),dimension(:),allocatable :: shots
        character(i_str_xxxlen) :: tmp

        fshot=1; dshot=1; lshot=1

        if(present(o_nshot)) then
            lshot=o_nshot
            nshot=o_nshot

        else
            shots=partition(setup_get_char('SHOT_INDEX','ISHOT'),o_separator=':')

            if(shots(1)%string=='') then !if not given, check disk
                call hud('SHOT_INDEX is not given. Now count how many data files exist in the directory..')
                data_file=setup_get_char('FILE_DATA')

                i=0
                alive=.true.
                do while (alive)
                    i=i+1
                    inquire(file=data_file//num2str(i,'(i0.4)')//'.su', size=file_size, exist=alive)
                    if(file_size==0) alive=.false.
                enddo

                lshot=i-1 !upto continuous indexing
                nshot=i-1
                
                call hud('Found '//num2str(lshot)//' sequential shots.')

            else
                if     (size(shots)==1) then
                    lshot=str2int(shots(1)%string)
                elseif (size(shots)==2) then
                    fshot=str2int(shots(1)%string)
                    lshot=str2int(shots(2)%string)
                else
                    fshot=str2int(shots(1)%string)
                    dshot=str2int(shots(2)%string)
                    lshot=str2int(shots(3)%string)
                endif

                lshot = lshot-mod((lshot-fshot),dshot)  !reduce lshot in case (lshot-fshot) is not a multiple of dshot
                nshot = (lshot-fshot)/dshot +1

                if(mpiworld%is_master) then
                    write(*,'(4(a,i0.4),a)') 'Will use (',lshot,'-',fshot,')/',dshot,'+1 = ',nshot, ' shots.'
                endif

            endif

        endif
        
        !message
        if(nshot<mpiworld%nproc) then
            call warn('No. of processors: '//num2str(mpiworld%nproc)//s_return// &
                    'No. of shots: '//num2str(nshot)//s_return// &
                    'Shot number < Processor number. Some processors will be always idle..')
        endif
        
        !build shotlist
        call alloc(shotlist, ceiling(nshot*1./mpiworld%nproc) )


        k=1 !shotlist's index
        j=0 !modulo shot#/processor#
        !loop over shots
        do i=fshot,lshot,dshot !shot#
            if (j==mpiworld%nproc) then
                j=j-mpiworld%nproc
                k=k+1
            endif
            if (mpiworld%iproc==j) then
                shotlist(k)=i
            endif
            j=j+1
        enddo
        
        !actual number of shots for each processor
        shotlist_nshot_per_processor=count(shotlist>0)       
        
        !message
        write(*,'(a,i2,a)') ' Proc# '//mpiworld%sproc//' has ',shotlist_nshot_per_processor,' assigned shots.'
        call hud('See file "shotlist" for details.')

        !write shotlist to disk
        write(tmp, *)  shotlist
        call mpiworld_file_write('shotlist', 'Proc# '//mpiworld%sproc//' has '//int2str(shotlist_nshot_per_processor,'(i4)')//' assigned shots:'//trim(tmp))
        
    end subroutine

end
