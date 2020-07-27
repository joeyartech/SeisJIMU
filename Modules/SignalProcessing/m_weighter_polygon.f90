! Modified from sumute.c in Seismic Unix
! dir: SeisUnix/src/su/main/windowing_sorting_muting

module m_weighter_polygon
use m_sysio
use m_arrayop

    contains
    
    subroutine build_weight_polygon(nt,dt,ntr,aoffset,weight)
        integer nt,ntr
        real    dt,aoffset(ntr)
        real    weight(nt,ntr)
        
        character(:),allocatable :: file_weight
        character(80) :: text
        integer,parameter :: max_gain_length=20 !maximum number of points
        real,dimension(max_gain_length) :: xgain,tgain
        

        file_weight=setup_get_file('FILE_WEIGHT_POLYGON')
        if(file_weight=='') then
            call hud('FILE_WEIGHT_POLYGON does NOT exist; weight unchanged.')
            return
        endif


        open(10,file=file_weight,action='read')

        gain_old=1.

        loopfile: do
            !initialize
            xgain=-99999.
            tgain=-99999.

            !read xgain vector as string
            loopxline: do
                read(10,"(a)",iostat=msg) text
                text=trim(adjustl(text))
                i=index(text,"#"); if(i>0) text=text(1:i-1) !remove comments after #
                i=index(text,"!"); if(i>0) text=text(1:i-1) !remove comments after !

                if (msg < 0) then !end of file
                    call hud('ERROR: FILE_WEIGHT_POLYGON is missing xgain points.')
                    call hud('Code stop now!')
                    stop
                endif
                if (msg > 0) stop 'Check FILE_WEIGHT_POLYGON.  Something is wrong..'
                if (text=='') then !blank line
                    cycle
                else
                    exit loopxline
                endif
            enddo loopxline

            if (text=='end' .or. text=='END') exit loopfile !end of file

            !convert string to real numbers
            read(text,*,iostat=msg) xgain
            mx=count(xgain>=0.) !number of input xgain

            !xgain should be increasing
            do i=1,mx-1
                if(xgain(i+1)<xgain(i)) then
                    call hud('ERROR: xgain from FILE_WEIGHT_POLYGON is NOT increasing!')
                    call hud('Code stop now!')
                    stop
                end if
            end do

            !read tgain vector as string
            looptline: do
                read(10,"(a)",iostat=msg) text
                text=trim(adjustl(text))
                i=index(text,"#"); if(i>0) text=text(1:i-1) !remove comments after #
                i=index(text,"!"); if(i>0) text=text(1:i-1) !remove comments after !

                if (msg < 0) then !end of file
                    call hud('ERROR: FILE_WEIGHT_POLYGON is missing tgain points.')
                    call hud('Code stop now!')
                    stop
                endif
                if (msg > 0) stop 'Check FILE_WEIGHT_POLYGON.  Something is wrong..'
                if (text=='') then !blank line
                    cycle
                else 
                    exit looptline
                endif
            enddo looptline

            !convert string to real numbers
            read(text,*,iostat=msg) tgain(1:mx)

            if (any(tgain(1:mx)<0.)) then
                call hud('ERROR: FILE_WEIGHT_POLYGON has unequal tgain & xgain pairs.')
                call hud('Code stop now!')
                stop
            endif


            !read in gain value and taper
            loopgline: do
                read(10,"(a)",iostat=msg) text
                text=trim(adjustl(text))
                i=index(text,"#"); if(i>0) text=text(1:i-1) !remove comments after #
                i=index(text,"!"); if(i>0) text=text(1:i-1) !remove comments after !

                if (msg < 0) then !end of file
                    call hud('ERROR: FILE_WEIGHT_POLYGON is missing tgain points.')
                    call hud('Code stop now!')
                    stop
                endif
                if (msg > 0) stop 'Check FILE_WEIGHT_POLYGON.  Something is wrong..'
                if (text=='') then !blank line
                    cycle
                else 
                    exit loopgline
                endif
            enddo loopgline

            !convert string to real numbers
            read(text,*,iostat=msg) gain, taper
            ntaper=nint(taper/dt)

            !loop of offset
            do itr=1,ntr
                
                !find the two xgain's closest to aoffset
                if(aoffset(itr)<=xgain(1)) then
                    jx1=1
                    jx2=1
                    dx1=0.
                    dx2=1.
                elseif(aoffset(itr)>=xgain(mx)) then
                    jx1=mx
                    jx2=mx
                    dx1=1.
                    dx2=0.
                else
                    jx1=maxloc(xgain(1:mx), dim=1, mask=xgain(1:mx)<aoffset(itr))
                    jx2=jx1+1
                    dx1=(aoffset(itr)-xgain(jx1)) / (xgain(jx2)-xgain(jx1))  !reduced distance from aoffset to xgain(jx1)
                    dx2=(xgain(jx2)-aoffset(itr)) / (xgain(jx2)-xgain(jx1))  !reduced distance from aoffset to xgain(jx2)
                endif
                
                !interp time by the two tgain & xgain pairs
                time = tgain(jx1)*dx2 + tgain(jx2)*dx1

                itime = nint(time/dt)+1  !assume all receivers have same dt ...
                itime_start = itime-ntaper  !1 <= itime_start <= itime <= nt

                !apply linear taper from itime_start+1 to itime
                do it=min(max(itime_start+1,1),nt) , &
                      min(max(itime        ,1),nt)
                    
                    k=it-itime_start-1

                    weight(it,itr) = ( gain_old*(ntaper-k) + gain*k ) / ntaper
                end do               

                !apply constant gain from itime+1 to nt
                it = min(max(itime+1,1),nt)
                weight(it:nt,itr) = gain
                
            end do
            
            gain_old=gain

        end do loopfile

        close(10)

    end subroutine

end
