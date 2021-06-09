! Modified from sumute.c in Seismic Unix
! dir: SeisUnix/src/su/main/windowing_sorting_muting

module m_separater
use m_sysio
use m_arrayop

    contains
    
    subroutine build_separater(nt,dt,ntr,aoffset,sepa_div,sepa_rfl)
        integer nt,ntr
        real    dt,aoffset(ntr)
        real,dimension(nt,ntr) :: sepa_div,sepa_rfl
        
        character(:),allocatable :: file_sepa
        character(80) :: text
        integer,parameter :: max_sepa_length=20 !maximum number of points
        real,dimension(max_sepa_length) :: xsepa,tsepa
        

        file_sepa=get_setup_file('FILE_SEPARATER')
        if(file_sepa=='') then
            if(mpiworld%is_master) write(*,*) 'ERROR: FILE_SEPARATER does NOT exist!'
            error stop
        endif


        open(10,file=file_sepa,action='read')

        sepa_div=0.; sepa_rfl=0.

        !initialize
        xsepa=-99999.
        tsepa=-99999.

        !read xsepa vector as string
        loopxline: do
            read(10,"(a)",iostat=msg) text
            text=trim(adjustl(text))
            i=index(text,"#"); if(i>0) text=text(1:i-1) !remove comments after #
            i=index(text,"!"); if(i>0) text=text(1:i-1) !remove comments after !

            if (msg < 0) then !end of file
                call hud('ERROR: FILE_SEPARATER is missing xsepa points.')
                call hud('Code stop now!')
                stop
            endif
            if (msg > 0) stop 'Check FILE_SEPARATER.  Something is wrong..'
            if (text=='') then !blank line
                cycle
            else
                exit loopxline
            endif
        enddo loopxline

        !convert string to real numbers
        read(text,*,iostat=msg) xsepa
        mx=count(xsepa>=0.) !number of input xsepa

        !xsepa should be increasing
        do i=1,mx-1
            if(xsepa(i+1)<xsepa(i)) then
                call hud('ERROR: xsepa from FILE_SEPARATER is NOT increasing!')
                call hud('Code stop now!')
                stop
            end if
        end do


        !read tsepa vector as string
        looptline: do
            read(10,"(a)",iostat=msg) text
            text=trim(adjustl(text))
            i=index(text,"#"); if(i>0) text=text(1:i-1) !remove comments after #
            i=index(text,"!"); if(i>0) text=text(1:i-1) !remove comments after !

            if (msg < 0) then !end of file
                call hud('ERROR: FILE_SEPARATER is missing tsepa points.')
                call hud('Code stop now!')
                stop
            endif
            if (msg > 0) stop 'Check FILE_SEPARATER.  Something is wrong..'
            if (text=='') then !blank line
                cycle
            else 
                exit looptline
            endif
        enddo looptline

        !convert string to real numbers
        read(text,*,iostat=msg) tsepa(1:mx)

        if (any(tsepa(1:mx)<0.)) then
            call hud('ERROR: FILE_SEPARATER has unequal tsepa & xsepa pairs.')
            call hud('Code stop now!')
            stop
        endif


        !loop of offset
        do itr=1,ntr
            
            !find the two xsepa's closest to aoffset
            if(aoffset(itr)<=xsepa(1)) then
                jx1=1
                jx2=1
                dx1=0.
                dx2=1.
            elseif(aoffset(itr)>=xsepa(mx)) then
                jx1=mx
                jx2=mx
                dx1=1.
                dx2=0.
            else
                jx1=maxloc(xsepa(1:mx), dim=1, mask=xsepa(1:mx)<aoffset(itr))
                jx2=jx1+1
                dx1=(aoffset(itr)-xsepa(jx1)) / (xsepa(jx2)-xsepa(jx1))  !reduced distance from aoffset to xsepa(jx1)
                dx2=(xsepa(jx2)-aoffset(itr)) / (xsepa(jx2)-xsepa(jx1))  !reduced distance from aoffset to xsepa(jx2)
            endif
            
            !interp time by the two tsepa & xsepa pairs
            time = tsepa(jx1)*dx2 + tsepa(jx2)*dx1

            itime = nint(time/dt)+1  !assume all receivers have same dt ...
            itime = min(max(itime,1),nt)

            sepa_div(1:itime,itr)=1.
            sepa_rfl(itime+1:nt,itr)=1.

        end do
            
        close(10)

    end subroutine

end
