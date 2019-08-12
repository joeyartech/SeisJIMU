module m_shotlist
use m_message
use m_arrayop

    integer,dimension(:),allocatable :: shotlist
    
    integer nshot_per_processor

    contains
    
    subroutine build_shotlist(nshots)
        
        if(mpiworld%is_master) then
            write(*,*) 'No. of processors:', mpiworld%nproc
            write(*,*) 'No. of shots:',nshots
            if(nshots<mpiworld%nproc) then
                call hud('ERROR: Shot number < Processor number. Some processors will be always idle..')
                stop
            endif
        endif
        
        l=ceiling(nshots*1./mpiworld%nproc)
        
        call alloc(shotlist,l)
        
        k=1 !shotlist's index
        j=0 !modulo shot#/processor#
        !loop over shots
        do i=1,nshots !shot# starts from 1
            if (j==mpiworld%nproc) then
                j=j-mpiworld%nproc
                k=k+1
            endif
            if (mpiworld%iproc==j) then
                shotlist(k)=i
            endif
            j=j+1
        enddo
        
        !real number of shots for each processor
        nshot_per_processor=count(shotlist>0)
        
        write(*,'(a,i2,a)') ' Proc# '//mpiworld%cproc//' has ',nshot_per_processor,' assigned shots.'
        
        !if nshot_per_processor is not same for each processor,
        !update_wavelet='stack' mode will be stuck due to collective communication in m_matchfilter.f90
        if(mpiworld%is_master) then
            if(nshot_per_processor * mpiworld%nproc /= nshots) then
                write(*,*) 'WARNING: unequal shot numbers on processors. If you are using UPDATE_WAVELET=''stack'', the code will be stuck due to collective communication in m_matchfilter'
                write(*,*) 'WARNING: Therefore you should STOP right now!'
            endif
        endif
        
    end subroutine

end