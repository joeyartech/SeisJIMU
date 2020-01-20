module m_sysio
use m_message

    character(:),allocatable :: file_setup
    
    contains
    
    subroutine init_setup(istat)
        integer :: istat
        logical :: alive
        character(256) :: tmp
        
        !get setup file
        istat=0
        
        call getarg(1,tmp)
        file_setup=trim(adjustl(tmp))
        
        if(file_setup=='') then
            istat=0
            call hud('No input setup file. Print manual..')
            return
        else
            inquire(file=file_setup, exist=alive)
            if(.not.alive) then
                call hud('Setup file '//file_setup//' does NOT exist!')
                istat=-1
                return
            endif
            
            call hud('Setup file: '//file_setup)
            istat=1
            return
            
        endif
    end subroutine
    
    function ask_setup(word,word2,default) result(found)
        character(*) :: word
        character(*),optional :: word2 !alias of inquired word
        logical,optional :: default
        
        character(80) :: text
        character(80) :: tmp,tmp2,tmp3
        
        logical :: found
        
        found=.false.
        if(present(default)) found=default
        
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp)) ) then
                found=.true.
            endif
        end do
        close(10)
    end function
    
    integer function get_setup_int(word,word2,default)
        character(*) :: word
        character(*),optional :: word2 !alias of inquired word
        integer,optional :: default
        
        character(80) :: text
        character(80) :: tmp,tmp2,tmp3
        
        logical :: found
        
        found=.false.
                
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    read(text,*) tmp2, get_setup_int
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp)) ) then
                read(text,*) tmp2, get_setup_int
                found=.true.
            endif
        end do
        close(10)
        
        if(found) then
            write(tmp,*) get_setup_int
            call hud(word//' '//trim(adjustl(tmp)))
        else
            if(present(default)) then
                get_setup_int=default
                write(tmp,*) default
                call hud(word//' is NOT found, take default: '//trim(adjustl(tmp)))
            else
                get_setup_int=0
                call hud(word//' is NOT found, take 0')
            endif
        endif
        
    end function
    
    real function get_setup_real(word,word2,default)
        character(*) :: word
        character(*),optional :: word2
        real, optional :: default
        
        character(80) :: text
        character(80) :: tmp,tmp2,tmp3
        
        logical :: found
        
        found=.false.
        
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    read(text,*) tmp2, get_setup_real
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp)) ) then
                read(text,*) tmp2, get_setup_real
                found=.true.
            endif
        end do
        close(10)
        
        if(found) then
            write(tmp,*) get_setup_real
            call hud(word//' '//trim(adjustl(tmp)))
        else
            if(present(default)) then
                get_setup_real=default
                write(tmp,*) default
                call hud(word//' is NOT found, take default: '//trim(adjustl(tmp)))
                
            else
                get_setup_real=0.
                call hud(word//' is NOT found, take 0.')
            endif
        endif
        
    end function
    
    function get_setup_char(word,word2,default)
        character(:),allocatable :: get_setup_char
        character(*) :: word
        character(*),optional :: word2
        character(*), optional :: default
        
        character(160) :: text
        character(160) :: tmp,tmp2,tmp3
        
        logical :: found
        
        found=.false.
        
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    read(text,*) tmp2, tmp3
                    get_setup_char=trim(adjustl(tmp3))
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp))) then
                read(text,*) tmp2, tmp3
                get_setup_char=trim(adjustl(tmp3))
                found=.true.
            endif
        end do
        close(10)
        
        if(found) then
            if(mpiworld%is_master) write(*,*) word//' '//get_setup_char
        else
            if(present(default)) then
                get_setup_char=default
                call hud(word//' is NOT found, take default: '//get_setup_char)
            else
                get_setup_char=''
                call hud(word//' is NOT found, take empty string')
            endif
        endif
        
    end function
    
    function get_setup_file(word,word2)
        character(:),allocatable :: get_setup_file
        character(*) :: word
        character(*),optional :: word2
        
        character(80) :: text
        character(80) :: tmp,tmp2,tmp3,tmp4
        
        logical :: found, alive
        
        found=.false.
        alive=.false.
        
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    read(text,*) tmp2, tmp3
                    tmp4=tmp3
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp))) then
                read(text,*) tmp2, tmp3
                tmp4=tmp3
                found=.true.
            endif
        end do
        close(10)
        
        if(found) inquire(file=trim(adjustl(tmp4)),exist=alive)
        
        found=found.and.alive
        
        
        if(found) then
            get_setup_file=trim(adjustl(tmp4))
            if(mpiworld%is_master) write(*,*) word//' '//get_setup_file
        else
            get_setup_file=''
            call hud(word//' is NOT found, take empty filename')
        endif
        
    end function
        
    logical function get_setup_logical(word,word2,default)
        character(*) :: word
        character(*),optional :: word2
        logical, optional :: default
        
        character(80) :: text
        character(80) :: tmp,tmp2,tmp3
        
        logical :: found
        
        found=.false.
        
        open(10,file=file_setup,action='read')
        do
            read (10,"(a)",iostat=msg) text ! read line into character variable
            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check setup file.  Something is wrong.'
            if (text=='') cycle !blank line
            
            read (text,*) tmp
            if(present(word2)) then
                if(word2==trim(adjustl(tmp))) then
                    read(text,*) tmp2, get_setup_logical
                    found=.true.
                endif
            endif
            if(word==trim(adjustl(tmp)) ) then
                read(text,*) tmp2, get_setup_logical
                found=.true.
            endif
        end do
        close(10)
        
        if(found) then
            write(tmp,*) get_setup_logical
            call hud(word//' '//trim(adjustl(tmp)))
        else
            if(present(default)) then
                get_setup_logical=default
                write(tmp,*) default
                call hud(word//' is NOT found, take default: '//trim(adjustl(tmp)))
                
            else
                get_setup_logical=.false.
                call hud(word//' is NOT found, take .false.')
            endif
        endif
        
    end function

!     subroutine write_master(array,filename)
!         real,dimension(*) :: array
!         character(*) :: filename
!         if(mpiworld%iproc==0) then
!             open(66,file=trim(adjustl(filename)),access='direct',recl=4*size(array))
!             write(66,rec=1)array
!             close(66)
!         endif
!     end
!     
!     
!     subroutine write_mpi_ascii(array,prefix)
!         real,dimension(*) :: array
!         character(*) :: prefix
!         integer :: i
!             
!         if(setup%nshot_out<0) then !all processors write 
!             open(66,file=trim(adjustl(prefix//mpiworld%cproc)),action="write",position='append')
!             write(66,*) "Iter, fcost's",inv%iter, array
!             close(66)
!         else !selected processors write
!             do i=1,setup%nshot_out
!                 if(mpiworld%iproc==setup%ishot_out(i)) then
!                 open(66,file=trim(adjustl(prefix//mpiworld%cproc)),action="write",position='append')
!                 write(66,*) "Iter, fcost's", inv%iter, array
!                 close(66)
!                 end if
!             end do
!         end if
!     end subroutine
!   
!     subroutine write_mpi_binary(array,prefix)
!         real,dimension(*) :: array
!         character(*) :: prefix
!         integer :: i
!         
!         if(setup%nshoto<0) then !all processors write 
!             open(66,file=trim(adjustl(prefix//mpiworld%cproc)),access='direct',recl=4*size(array))
!             write(66,rec=1)array
!             close(66)
!         else !selected processors write
!             do i=1,setup%nshoto
!                 if(mpiworld%iproc==setup%ishoto(i)) then
!                 open(66,file=trim(adjustl(prefix//mpiworld%cproc)),access='direct',recl=4*size(array))
!                 write(66,rec=1)array
!                 close(66)
!                 end if
!             enddo
!         end if
!     end subroutine
!     
!
!     subroutine write_parameters(mode)
!         integer mode, n
!         integer,save :: k=1
!         n=4*m%nx*m%ny*m%nz
!         
!         if(mpi_world%iproc==0) then
!         
!             select case (mode)
!                 case (1) !write during line search
!                     !just for fun
!                     open(42,file='velocity_history',access='direct',recl=n)
!                     write(42,rec=k)m%vp
!                     close(42)
!                     k=k+1
!                     
!                 case (2) !write for new iteration
!                     if(if_par_vp) then
!                         open(42,file='param_vp_inter',access='direct',recl=n)
!                         write(42,rec=iter)m%vp
!                         close(42)
!                     endif
!                     if(if_par_rho) then
!                         open(42,file='param_rho_inter',access='direct',recl=n)
!                         write(42,rec=iter)m%rho
!                         close(42)
!                     endif
!                     if(if_par_ip) then
!                         open(42,file='param_ip_inter',access='direct',recl=n)
!                         write(42,rec=iter)m%vp*m%rho
!                         close(42)
!                     endif
!                     if(if_par_kappa) then
!                         open(42,file='param_kappa_inter',access='direct',recl=n)
!                         write(42,rec=iter)m%vp*m%vp*m%rho
!                         close(42)
!                     endif
!                 
!                 case (3) ! final writing
!                     open(42,file='param_vp_final',access='direct',recl=n)
!                     write(42,rec=1)m%vp
!                     close(42)
!                     open(42,file='param_rho_final',access='direct',recl=n)
!                     write(42,rec=1)m%rho
!                     close(42)
!                     open(42,file='param_ip_final',access='direct',recl=n)
!                     write(42,rec=1)m%vp * m%rho
!                     close(42)
!                 
!             end select
! 
!         end if
!     end subroutine

end