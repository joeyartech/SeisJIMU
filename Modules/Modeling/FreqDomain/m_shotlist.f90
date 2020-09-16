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


        ! k=1 !shotlist's index
        ! j=0 !modulo shot#/processor#
        ! !loop over shots
        ! do i=fshot,lshot,dshot !shot#
        !     if (j==mpiworld%nproc) then
        !         j=j-mpiworld%nproc
        !         k=k+1
        !     endif
        !     if (mpiworld%iproc==j) then
        !         shotlist(k)=i
        !     endif
        !     j=j+1
        ! enddo
        
        ! !real number of shots for each processor
        ! nshot_per_processor=count(shotlist>0)       
        
        ! !message
        ! write(*,'(a,i2,a)') ' Proc# '//mpiworld%cproc//' has ',nshot_per_processor,' assigned shots.'
        ! call hud('See file "shotlist" for details.')

        ! !write shotlist to disk
        ! write(cnshot_per_processor,'(i4)') nshot_per_processor
        ! write(string, *)  shotlist
        ! call mpiworld_file_write('shotlist', 'Proc# '//mpiworld%cproc//' has '//cnshot_per_processor//' assigned shots:'//trim(string))

        ! !if nshot_per_processor is not same for each processor,
        ! !update_wavelet='stack' mode will be stuck due to collective communication in m_matchfilter.f90
        ! if(mpiworld%is_master) then
        !     if(nshot_per_processor * mpiworld%nproc /= nshots) then
        !         write(*,*) 'WARNING: unequal shot numbers on processors. If you are using UPDATE_WAVELET=''stack'', the code will be stuck due to collective communication in m_matchfilter'
        !         write(*,*) 'WARNING: Therefore you should STOP right now!'
        !     endif
        ! endif
        
    end subroutine

end
