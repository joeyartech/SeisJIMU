module m_shotlist
use m_mpienv
use m_message
use m_arrayop

    private
    public build_shotlist, shotlist, nshots!, nshot_per_processor

    integer,dimension(:),allocatable :: shotlist   

    integer :: nshots!, nshot_per_processor
        
    contains
    
    subroutine build_shotlist(o_nshots)
        integer,optional :: o_nshots

        integer :: fshot, dshot, lshot
        character(4) :: cindex, cnshot_per_processor
        integer :: file_size
        logical :: alive
        character(:),allocatable :: cshotno, data_file, text
        ! character(240)::string

        fshot=1; dshot=1; lshot=1

        if(present(o_nshots)) then !number of shots is internally given
            nshots=o_nshots
            lshot =o_nshots

        else !use other info
            cshotno=get_setup_char('SHOTNO')

            if(cshotno=='') then !if SHOTNO not given, check dist
                call hud('SHOTNO is not given. Now count how many data file exists in the directory..')
                data_file=get_setup_char('FILE_DATA')

                i=0
                alive=.true.
                do while (alive)
                    i=i+1
                    write(cindex,'(i0.4)') i
                    inquire(file=data_file//cindex//'.su', size=file_size, exist=alive)
                    if(file_size==0) alive=.false.
                enddo

                lshot=i-1 !upto continuous indexing
                nshots=i-1
                
                if(mpiworld%is_master) then
                    write(*,*) 'Found', lshot, 'sequential shots.'
                endif

            else !extract the shot range and increment from SHOTNO
                k1=index(cshotno,':')
                text=cshotno(1:k1-1)
                if(len(text)>0) then
                    read(text,*) fshot
                else
                    fshot=-99999  !if fshot not given, fshot will be same as lshot
                endif

                k2=index(cshotno(k1+1:),':')+k1
                text=cshotno(k1+1:k2-1)
                if(len(text)>0) read(text,*) dshot

                !k3=index(string(k2+1:),':')+k2
                text=cshotno(k2+1:)
                if(len(text)>0) read(text,*) lshot

                if(fshot==-99999) fshot=lshot

                lshot = lshot-mod((lshot-fshot),dshot)  !reduce lshot in case (lshot-fshot) is not a multiple of dshot
                nshots = (lshot-fshot)/dshot +1

                if(mpiworld%is_master) then
                    write(*,'(4(a,i0.4),a)') 'Will use (',lshot,'-',fshot,')/',dshot,'+1 = ',nshots, ' shots.'
                endif

            endif

        endif
        
        !message
        ! if(mpiworld%is_master) then
        !     write(*,*) 'No. of processors:', mpiworld%nproc
        !     write(*,*) 'No. of shots:',nshots
        !     if(nshots<mpiworld%nproc) then
        !         call hud('ERROR: Shot number < Processor number. Some processors will be always idle..')
        !         stop
        !     endif
        ! endif
        
        !build shotlist
        ! l=ceiling(nshots*1./mpiworld%nproc)
        ! 
        ! call alloc(shotlist,l)
        call alloc(shotlist,nshots)
        shotlist=[(i,i=fshot,dshot,lshot)]

    end subroutine

end
