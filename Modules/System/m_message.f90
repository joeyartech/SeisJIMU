module m_message
use m_mpienv, only: mpiworld

    contains
    
    subroutine hud(str,iproc)
        character(*) :: str
        integer,optional :: iproc
        
        if(present(iproc)) then
            if(mpiworld%iproc==iproc) then
                write(*,*) str
            endif
        else
            if(mpiworld%is_master) then
                write(*,*) str
            endif
        endif
        
    end subroutine

end