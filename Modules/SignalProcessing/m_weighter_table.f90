module m_weighter_table
use m_sysio
use m_arrayop

use, intrinsic :: ieee_arithmetic

    contains
    
    subroutine build_weight_table(nt,dt,ntr,aoffset,weight)
        integer nt,ntr
        real    dt,aoffset(ntr)
        real    weight(nt,ntr)
        
        character(:),allocatable :: file_weight
        character(80) :: text
        character(1)  :: delim
        integer,parameter :: max_gain_length=20 !maximum number of points per line
        real,dimension(max_gain_length) :: xgain, tgain, gain
        real,dimension(:,:),allocatable :: table, tmp_table


        file_weight=get_setup_file('FILE_WEIGHT_TABLE')       
        if(file_weight=='') then
            call hud('FILE_WEIGHT_TABLE does NOT exist; weight unchanged.')
            return
        endif


        open(10,file=file_weight,action='read')

        xgain=ieee_value(1.0,ieee_quiet_nan)

        !loop to read xgain vector as string
        do
            read(10,"(a)",iostat=msg) text
            text=trim(adjustl(text))
            if(index(text,"#")) text=text(1:index(text,"#")-1) !remove comments after #
            if(index(text,"!")) text=text(1:index(text,"!")-1) !remove comments after !

            if (msg < 0) then !end of file
                call hud('ERROR: FILE_WEIGHT_TABLE is missing xgain points.')
                call hud('Code stop now!')
                stop
            endif
            if (msg > 0) stop 'Check FILE_WEIGHT_TABLE.  Something is wrong..'
            if (text=='') cycle !blank line
            if (msg == 0) exit
        enddo

        !convert string to real numbers
        read(text,*,iostat=msg) xgain
        mx=count(.not. ieee_is_nan(xgain)) !number of input (non-nan) xgain

        !xgain should be increasing
        do i=1,mx-1
            if(xgain(i+1)<xgain(i)) then
                call hud('ERROR: xgain from FILE_WEIGHT_TABLE is NOT increasing!')
                call hud('Code stop now!')
                stop
            end if
        end do

        allocate(tmp_table(max_gain_length,mx))
        
        gain=ieee_value(1.0,ieee_quiet_nan)
        tmp_table=ieee_value(1.0,ieee_quiet_nan)
        
        !loop to read tgain and gain vector as string
        k=1
        do           
            read(10,"(a)",iostat=msg) text
            text=trim(adjustl(text))
            if(index(text,"#")) text=text(1:index(text,"#")-1) !remove comments after #
            if(index(text,"!")) text=text(1:index(text,"!")-1) !remove comments after !

            if (msg < 0) exit !end of file
            if (msg > 0) stop 'Check FILE_WEIGHT_TABLE.  Something is wrong..'
            if (text=='') cycle !blank line
            if (text=='end' .or. text=='END') exit

            !convert string to real numbers
            read(text,*,iostat=msg) tgain(k), delim , gain

            if (k>1) then; if(tgain(k)<tgain(k-1)) then
                call hud('Please make sure tgain in FILE_WEIGHT_TABLE is increasing!')
                call hud('Code stop now!')
                stop
            endif; endif
            if ( count(.not. ieee_is_nan(gain)) /= mx ) then
                call hud('ERROR: FILE_WEIGHT_TABLE has gain vector length < xgain vector length.')
                call hud('Code stop now!')
                stop
            endif

            tmp_table(k,:) = gain
            k=k+1
        enddo

        !reshape table without nan elements
        mt=count(.not. ieee_is_nan(tmp_table(:,1)))  !number of input (non-nan) tgain
        allocate(table(mt,mx))
        table=tmp_table(1:mt,1:mx)
        deallocate(tmp_table)


        !map table to weight with interpolation
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


            do it=1,nt
                t=(it-1)*dt

                !find the two tgain's closest to aoffset
                if(t<=tgain(1)) then
                    jt1=1
                    jt2=1
                    dt1=0.
                    dt2=1.
                elseif(t>=tgain(mt)) then
                    jt1=mt
                    jt2=mt
                    dt1=1.
                    dt2=0.
                else
                    jt1=maxloc(tgain(1:mt), dim=1, mask=tgain(1:mt)<t)
                    jt2=jt1+1
                    dt1=(t-tgain(jt1)) / (tgain(jt2)-tgain(jt1))
                    dt2=(tgain(jt2)-t) / (tgain(jt2)-tgain(jt1))
                endif

                
                weight(it,itr) = weight(it,itr) * ( &
                      ( table(jt1,jx1)*dt2 + table(jt2,jx1)*dt1 ) *dx2  &
                    + ( table(jt1,jx2)*dt2 + table(jt2,jx2)*dt1 ) *dx1  )

            enddo

        enddo

        deallocate(table)

    end subroutine

end
