module m_message
use m_mpienv, only: mpiworld

!ANSI colors
!http://fortranwiki.org/fortran/show/ansi_colors
!https://rosettacode.org/wiki/Terminal_control/Coloured_text#Fortran
!http://pueblo.sourceforge.net/doc/manual/ansi_color_codes.html

!foreground
character(*), parameter :: red     =achar(27)//'[31m'
character(*), parameter :: green   =achar(27)//'[32m'
character(*), parameter :: yellow  =achar(27)//'[33m'
character(*), parameter :: blue    =achar(27)//'[34m'
character(*), parameter :: magenta =achar(27)//'[35m'
character(*), parameter :: cyan    =achar(27)//'[36m'
character(*), parameter :: grey    =achar(27)//'[90m' !Bright-Black 
!one background color
character(*), parameter :: bg_black =achar(27)//'[40m'
character(*), parameter :: bg_default =achar(27)//'[49m'
!others
character(*), parameter :: bold = achar(27)//'[1m'
character(*), parameter :: bold_blink = achar(27)//'[1;5m'
character(*), parameter :: reset = achar(27)//'[0m' ! Terminates an ANSI code.

    contains
    
    subroutine hud(msg,o_iproc)
        character(*) :: msg
        integer,optional :: o_iproc
        
        if(present(o_iproc)) then
            if(mpiworld%iproc==o_iproc) then
                write(*,*) 'Proc# '//mpiworld%sproc//' : '//msg
            endif
        else
            if(mpiworld%is_master) then
                write(*,*) msg
            endif
        endif
        
    end subroutine

    subroutine warn(msg,o_iproc)
        character(*) :: msg
        integer,optional :: o_iproc

        if(present(o_iproc)) then
            if(mpiworld%iproc==o_iproc) then
                write(*,'(a,x,a)') 'Proc# '//mpiworld%sproc//' : '//bg_black//yellow//bold//'WARNING:'//reset, msg
            endif
        else
            if(mpiworld%is_master) then
                write(*,'(a,x,a)') bg_black//yellow//bold//'WARNING:'//reset, msg
            endif
        endif
        
    end subroutine

    subroutine error(msg,solution,o_iproc)
        character(*) :: msg
        character(*),optional :: solution
        integer,optional :: o_iproc

        if(present(o_iproc)) then
            if(mpiworld%iproc==o_iproc) then
                write(*,'(a,x,a)') 'Proc# '//mpiworld%sproc//' : '//bg_black//red//bold_blink//'ERROR:'//reset, msg
                if(present(solution)) then
                    write(*,'(a)') 'Possible solutions:'
                    write(*,'(a)') solution
                endif
            endif
        else
            if(mpiworld%is_master) then
                write(*,'(a,x,a)') bg_black//red//bold_blink//'ERROR:'//reset, msg
                if(present(solution)) then
                    write(*,'(a)') 'Possible solutions:'
                    write(*,'(a)') solution
                endif
            endif
        endif
        
        call mpiworld%fin
        error stop
        
    end subroutine

    subroutine fatal(msg)
        character(*) :: msg
        call error(msg)
    end subroutine

end
