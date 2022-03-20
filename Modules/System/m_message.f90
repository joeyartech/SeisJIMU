module m_message
use m_either
use m_string
use m_mpienv

    private
    public :: hud, warn, error, fatal

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
        
        if(mpiworld%iproc==either(o_iproc,0,present(o_iproc))) then
            write(*,*) either('Proc# '//mpiworld%sproc//' : ','',present(o_iproc))//msg
        endif
        
    end subroutine

    subroutine warn(msg,o_iproc)
        character(*) :: msg
        integer,optional :: o_iproc

        if(mpiworld%iproc==either(o_iproc,0,present(o_iproc))) then
            write(*,*) either('Proc# '//mpiworld%sproc//' : ','',present(o_iproc))//&
                bg_black//yellow//bold//'WARNING:'//reset//' '//msg
        endif
        
    end subroutine

    subroutine error(msg,o_solution,o_iproc)
        character(*) :: msg
        character(*),optional :: o_solution
        integer,optional :: o_iproc

        if(mpiworld%iproc==either(o_iproc,0,present(o_iproc))) then
            write(*,*) either('Proc# '//mpiworld%sproc//' : ','',present(o_iproc))//&
                bg_black//red//bold_blink//'ERROR:'//reset//' '//msg//&
                either(s_NL//'Possible solutions:'//o_solution,'',present(o_solution))
        endif

        call mpiworld%final
        error stop
        
    end subroutine

    subroutine fatal(msg)
        character(*) :: msg
        call error(msg)
    end subroutine

end
