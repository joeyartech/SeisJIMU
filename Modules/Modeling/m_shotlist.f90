module m_shotlist
use m_string
use m_mpienv
use m_message
use m_arrayop

    type t_shotlist
        integer,dimension(:),allocatable :: list, list_per_processor
        integer :: nshot, nshot_per_processor
    end type

    contains

    subroutine read_from_setup

        select case (setup%get_char('ACQUI_GEOMETRY',default='spread'))

        case ('irregularOBN')
            open(13,file=setup%get_file('FILE_SOURCE_POSITION','SPOS',mandatory=.true.),action='read')
            !count number of sources
            self%nshot=0
            do
                read (13,*,iostat=msg) z,x,y
                if(msg/=0) exit
                self%nshot=self%nshot+1
            enddo
            close(13)
            call hud('Will read'//num2str(self%nshot)//'source positions.')

        case default
            self%nshot=setup%get_int('NUMBER_SOURCE','NS',mandatory=.true.)
            
        end select        

        call alloc(self%list(self%nshot)); self%list=[(i,1,self%nshot)]

        call hud('Will compute '//num2str(self%nshot)//' synthetic shots.')

    end subroutine
    
    subroutine read_from_data
        integer file_size
        type(t_string),dimension(:),allocatable :: shots, subshots
        
        shots=setup%get_strs('SHOT_INDEX','ISHOT'))
        
        if(size(shots)==0) then !not given
            call hud('SHOT_INDEX is not given. Now count how many data files exist in the directory..')
            file=setup%get_char('FILE_DATA')

            i=0; exist=.true.
            do
                i=i+1
                inquire(file=self%file//num2str(i,'(i0.4)')//'.su', size=file_size, exist=exist)
                if(file_size==0) exist=.false.
                if(.not.exist) exit
            enddo

            self%nshot=i

            call alloc(self%list(self%nshot)); self%list=[(i,1,self%nshot)]
        
            call hud('Found '//num2str(self%nshot)//' sequential shots.')

        endif

        if(size(shots)>0) then !given
            
            call alloc(self%list(1))
            
            do i=1,size(shots)

                subshots=partition(shots(i)%s,o_sep=':')
            
                if(size(subshots)==1) then !add/rm
                    if(subshots(1)%s(1:1)='-') !rm
                        call rm(self%list,size(self%list),str2int(subshots(1)%s(2:)))
                    else
                        call add(self%list,size(self%list),str2int(subshots(1)%s))
                    endif
                endif
                if(size(subshots)==2) then !first:last
                    do j=str2int(subshots(1)%s),str2int(subshots(2)%s)
                        call add(self%list,size(self%list),j)
                    enddo
                endif
                if(size(subshots)==3) then !first:increment:last
                    do j=str2int(subshots(1)%s),str2int(subshots(3)%s),str2int(subshots(2)%s)
                        call add(self%list,size(self%list),j)
                    enddo
                endif               

            enddo

            self%nshot=size(self%list)

            call hud('Total number of shots: '//num2str(nshot)//'')
            
        endif

        deallocate(shots, subshots)

    end subroutine

    subroutine assign

        if(self%nshot<mpiworld%nproc) then
            call warn('No. of processors: '//num2str(mpiworld%nproc)//s_return// &
                    'No. of shots: '//num2str(self%nshot)//s_return// &
                    'Shot number < Processor number. Some processors will be always idle..')
        endif

        call alloc(self%list_per_processor, ceiling(self%nshot*1./mpiworld%nproc) )


        k=1 !shotlist's index
        j=0 !modulo shot#/processor#
        
        do i=1,self%nshot
            if (j==mpiworld%nproc) then
                j=j-mpiworld%nproc
                k=k+1
            endif
            if (mpiworld%iproc==j) then
                self%list_per_processor(k)=i
            endif
            j=j+1
        enddo

        self%nshot_per_processor=count(self%list_per_processor>=0)
        
    end subroutine

    subroutine assign_random
    end subroutine

    subroutine print
                
        !message
        write(*,'(a,i2,a)') ' Proc# '//mpiworld%sproc//' has ',self%nshot_per_processor,' assigned shots.'
        call hud('See file "shotlist" for details.')

        !write shotlist to disk
        write(tmp, *)  shotlist
        call mpiworld_file_write('shotlist', 'Proc# '//mpiworld%sproc//' has '//int2str(self%nshot_per_processor,'(i4)')//' assigned shots:'//trim(tmp))
        
    end subroutine

end
