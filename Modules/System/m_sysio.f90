module m_sysio
use m_mpienv
use m_message
use m_setup

    character(:),allocatable :: dir_in, dir_out

    contains

    subroutine sysio_init
        logical exist

        !input directory
        dir_in=setup%get_str('DIR_IN',o_default='./') 

#ifdef GNU
        inquire(file=dir_in,exist=exist)
#endif
#ifdef INTEL
        inquire(directory=dir_in,exist=exist)
#endif
        
        if(.not.exist) then !input directory not exist
            call hud('Input directory: '//dir_in//' does NOT exist, reset to ./')
            dir_in='./'
        endif

        call hud('Input directory: '//dir_in)

        !output directory
        dir_out=setup%get_str('DIR_OUT',o_default='./results')

        call execute_command_line('mkdir -p '//dir_out, wait=.true.)
        
        call hud('Output directory:'//dir_out)

    end subroutine

end
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
