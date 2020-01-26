module m_weighting
use m_sysio
use m_arrayop

    real,dimension(:,:),allocatable :: weight

    contains
    
    subroutine build_weighting(nt,dt,nx,aoffset)
        integer nt,nx
        real    dt,aoffset(nx)
        
        character(:),allocatable :: file_weight
        character(80) :: file_weight_table
        logical :: alive
        real,dimension(:),allocatable :: xgain,tgain
        real,dimension(:,:),allocatable :: table
        real :: local_dist
        
        call alloc(weight,nt,nx)
        
        weight=1.
        gain_old=0
        
        file_weight=get_setup_file('FILE_WEIGHT')
        
        if(file_weight=='') then
            call hud('FILE_WEIGHT does NOT exist. Use weight=1.')
            return
        endif
        
        open(11,file=file_weight)
        
        read(11,*) i_fgain
        
        !compute gain now
        do while(i_fgain>0) 
            !manage the front gain info
            read(11,*)ngain, gain, taper
            allocate(xgain(ngain),tgain(ngain))
            read(11,*)xgain
            read(11,*)tgain
            
            !check the increase of distance in the gain function
            do i=1,ngain-1
            if(xgain(i+1)<xgain(i)) then
                call hud('Please make sure xgain in FILE_WEIGHT is monotically increase')
                stop
            end if
            end do

            ntaper=nint(taper/dt)
            
            !now apply the gain function in the table
            do ix=1,nx
                !find the xgain interval
                idist=ngain
                do i=1,ngain
                    if(aoffset(ix)<xgain(i))then
                        idist=i-1
                        exit
                    end if
                end do
                if(idist>ngain) then
                    call hud('Some problem for xgain in FILE_WEIGHT')
                    stop
                elseif(idist==ngain) then
                    time_gain=tgain(ngain)
                else
                    local_dist= (aoffset(ix)-xgain(idist)) &
                            /(xgain(idist+1)-xgain(idist))
                    time_gain=tgain(idist)*(1-local_dist)+tgain(idist+1)*local_dist
                end if

                itime_gain=nint(time_gain/dt)+1  !assume all receivers have same dt ...
                itime_gain_start=itime_gain-ntaper  !1 <= itime_gain_start <= itime_gain <= nt
                ntaper_true=itime_gain-itime_gain_start

                !apply taper, linear function from itime_gain_start to itime_gain
                k=1
                do it=itime_gain_start+1,itime_gain
                    if(it<1.or.it>nt) cycle  !assume all receivers have same nt ...
                    weight(it,ix) = (gain-gain_old)*k/ntaper_true + gain_old
                    k=k+1
                end do
                
                !apply sepa from itime_sepa+1 to nt
                do it= itime_gain+1, nt
                    if(it<1.or.it>nt) cycle 
                    weight(it,ix)=gain
                end do
                
            end do

            deallocate(xgain,tgain)
            
            gain_old=gain
            i_fgain=i_fgain-1

        enddo
        

        read(11,*)file_weight_table
        inquire(file=file_weight_table, exist=alive)
        
        if (alive) then
            read(11,*)ns,ds,noff,doff
            allocate(table(ns,noff))
            open(13,file=file_weight_table,access='direct',recl=4*ns*noff)
            read(13,rec=1)table
            close(13)

            !compute weight now
            do ix=1,nx
                idist = int(aoffset(ix)/doff)+1
                red_dist=aoffset(ix)/doff-idist+1 !reduced coordinate between 0 and 1
                do it=1,nt
                    !compute time
                    time=(it-1)*dt
                    itime = int(time/ds)+1
                    red_time=time/ds-itime+1 !reduced coordinate between 0 and 1
                    !linear interpolation +extrapolation of the last value if the table is not big enough
                    temp1 = table( min(itime  ,ns) , min(idist  ,noff) ) * (1.-red_dist) &
                          + table( min(itime  ,ns) , min(idist+1,noff) ) * red_dist
                    temp2 = table( min(itime+1,ns) , min(idist  ,noff) ) * (1.-red_dist) &
                          + table( min(itime+1,ns) , min(idist+1,noff) ) * red_dist
                    weight(it,ix)=weight(it,ix)* (temp1*(1.-red_time)+temp2*(red_time))
                end do
            end do
            deallocate(table)

        end if
        
        close(13)
        
        close(11)

    end subroutine

end
