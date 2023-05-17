module m_weighter
use m_System
use m_Modeling

    private

    type,public :: t_weighter

        real,dimension(:,:),allocatable :: weight

        contains

        procedure :: update

    end type

    type(t_weighter),public :: wei

    integer nt,ntr
    real dt
    
    contains

    subroutine update(self,o_suffix,o_nt,o_dt,o_ntr)
        class(t_weighter) :: self
        character(*),optional :: o_suffix
        integer,optional :: o_nt, o_ntr
        real,optional :: o_dt

        real,dimension(:,:),allocatable :: weight
        character(:),allocatable :: suf,op
        type(t_string),dimension(:),allocatable :: list,sublist

        nt = either(o_nt,  shot%nt  , present(o_nt ))
        ntr= either(o_ntr, shot%nrcv, present(o_ntr))
        dt = either(o_dt,  shot%dt  , present(o_dt ))
        call alloc(self%weight,nt,ntr,o_init=1.)

        suf = either(o_suffix,'',present(o_suffix))

        list=setup%get_strs('WEIGHTING'//suf,'WEI'//suf,o_default=either('aoffset^1','aoffset^0.5',m%is_cubic))
        op=setup%get_str('WEIGHTING_OPERATIONS'//suf,'WEI_OP'//suf,o_default='multiply')

        do i=1,size(list)
            call alloc(weight,nt,ntr,o_init=1.)

            if(index(list(i)%s,'p*')>0) then !weight pressure components
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will multiply pressure component by '//sublist(2)%s)
                call by_components(weight,o_p=str2real(sublist(2)%s))
            endif

            if (index(list(i)%s,'vz*')>0) then
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will multiply vz component by '//sublist(2)%s)
                call by_components(weight,o_vz=str2real(sublist(2)%s))
            endif

            if (index(list(i)%s,'vx*')>0) then
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will multiply vx component by '//sublist(2)%s)
                call by_components(weight,o_vx=str2real(sublist(2)%s))
            endif

            if (index(list(i)%s,'vy*')>0) then
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will multiply vy component by '//sublist(2)%s)
                call by_components(weight,o_vy=str2real(sublist(2)%s))
            endif

            if (index(list(i)%s,'aoffset^')>0) then !weight traces based on power of aoffset
                sublist=split(list(i)%s,o_sep='^')
                call hud('Will weight traces by aoffset^'//sublist(2)%s)
                call by_aoffset(weight,o_power=str2real(sublist(2)%s))
            endif
            
            if (index(list(i)%s,'aoffset*')>0) then !use aoffset to multiply on traces
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will weight traces by aoffset*'//sublist(2)%s)
                call by_aoffset(weight,o_factor=str2real(sublist(2)%s))
            endif

            if (index(list(i)%s,'aoffset_range')>0) then !use window defined by aoffset
                sublist=split(list(i)%s,o_sep=':')
                call hud('Will weight traces by aoffset range:'//sublist(2)%s//':'//sublist(3)%s)
                call by_aoffset_range(weight,str2real(sublist(2)%s),str2real(sublist(3)%s))
            endif

            if (index(list(i)%s,'time^')>0) then !weight traces by power of time
                sublist=split(list(i)%s,o_sep='^')
                call hud('Will weight traces by time^'//sublist(2)%s)
                call by_time(weight,o_power=str2real(sublist(2)%s))
            endif
            
            if (index(list(i)%s,'time_window')>0) then !use window defined by time
                sublist=split(list(i)%s,o_sep=':')
                call hud('Will weight traces by time window:'//sublist(2)%s//':'//sublist(3)%s)
                call by_time_window(weight,str2real(sublist(2)%s),str2real(sublist(3)%s))
            endif

            if (index(list(i)%s,'polygon')>0) then !polygon defined weights to multiply on traces
                sublist=split(list(i)%s,o_sep=':')
                call hud('Will weight traces with polygons defined in '//sublist(2)%s)
                call by_polygon(weight,sublist(2)%s)
            endif
            
            if (index(list(i)%s,'table')>0) then !table defined weights to multiply on traces
                sublist=split(list(i)%s,o_sep=':')
                call hud('Will weight traces with tables defined in '//sublist(2)%s)
                call by_table(weight,sublist(2)%s)
            endif
                
            if (index(list(i)%s,'custom')>0) then
                sublist=split(list(i)%s,o_sep=':')
                call hud('Will weight traces in a custom way defined in '//sublist(2)%s)
                call by_custom(weight,sublist(2)%s)
            endif

            if (index(list(i)%s,'per_shot')>0) then
                sublist=split(list(i)%s,o_sep=':')
                call hud('Will weight traces per shot defined in '//sublist(2)%s)
                call by_pershot(weight,sublist(2)%s)
            endif


	    select case(op)
                case ('multiply')
                self%weight = self%weight * weight

                case ('add')
                self%weight = self%weight + weight

                case ('overlay')
                self%weight = weight

            end select

        enddo

        deallocate(weight)

        do ir=1,shot%nrcv
            if (shot%rcv(ir)%is_badtrace) self%weight(:,ir)=0. !bad trace
        enddo

        if(mpiworld%is_master) call suformat_write('weights'//suf,self%weight,nt,ntr,o_dt=dt)
        
    end subroutine


    subroutine by_components(weight,o_p,o_vz,o_vx,o_vy)
        real,dimension(:,:) :: weight
        real,optional :: o_p, o_vz, o_vx, o_vy

        if(present(o_p)) then
            do i=1,shot%nrcv
                if(shot%rcv(i)%comp=='p') weight(:,i)=weight(:,i)*o_p
            enddo
        endif

        if(present(o_vz)) then
            do i=1,shot%nrcv
                if(shot%rcv(i)%comp=='vz') weight(:,i)=weight(:,i)*o_vz
            enddo
        endif

        if(present(o_vx)) then
            do i=1,shot%nrcv
                if(shot%rcv(i)%comp=='vx') weight(:,i)=weight(:,i)*o_vx
            enddo
        endif

        if(present(o_vy)) then
            do i=1,shot%nrcv
                if(shot%rcv(i)%comp=='vy') weight(:,i)=weight(:,i)*o_vy
            enddo
        endif

    end subroutine

    subroutine by_aoffset(weight,o_power,o_factor)
        real,dimension(:,:) :: weight
        real,optional :: o_power, o_factor

        if(present(o_power)) then
            do i=1,shot%nrcv    
                weight(:,i)=weight(:,i)*shot%rcv(i)%aoffset**o_power
            enddo
        endif

        if(present(o_factor)) then
            do i=1,shot%nrcv    
                weight(:,i)=weight(:,i)*shot%rcv(i)%aoffset*o_factor
            enddo
        endif

    end subroutine

    subroutine by_time(weight,o_power,o_factor)
        real,dimension(:,:) :: weight
        real,optional :: o_power, o_factor

        if(present(o_power)) then
            do i=1,shot%nt
                weight(i,:)=weight(i,:)*((i+0.5-1)*shot%dt)**o_power
            enddo
        endif

!        if(present(o_factor)) then
!            do i=1,shot%nrcv    
!                weight(:,i)=weight(:,i)*shot%rcv(i)%aoffset*o_factor
!            enddo
!        endif

    end subroutine

    subroutine by_aoffset_range(weight,aoff_min,aoff_max)
        real,dimension(:,:) :: weight

        do i=1,shot%nrcv    
            if (shot%rcv(i)%aoffset<aoff_min) weight(:,i)=0.
            if (shot%rcv(i)%aoffset>aoff_max) weight(:,i)=0.
        enddo

    end subroutine

    subroutine by_time_window(weight,tmin,tmax)
        real,dimension(:,:) :: weight

        do i=1,shot%nt
            t=(i-1)*shot%dt
            if (t<tmin) weight(i,:)=0.
            if (t>tmax) weight(i,:)=0.
        enddo

    end subroutine

    ! Modified from sumute.c in Seismic Unix
    ! dir: SeisUnix/src/su/main/windowing_sorting_muting
    subroutine by_polygon(weight,file)
        real,dimension(:,:) :: weight
        character(*) :: file

        character(i_str_slen) :: text
        integer,parameter :: max_gain_length=20 !maximum number of points
        real,dimension(max_gain_length) :: xgain,tgain

        real,dimension(:,:),allocatable :: weight

        call alloc(weight,self%nt,self%ntr)
        
        open(10,file=file,action='read')

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
                    call error(file//' is missing xgain points.')
                endif
                if (msg > 0) then
                    call error('Check '//file//'.  Something is wrong..')
                endif
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
                    call error('xgain from '//file//' is NOT increasing!')
                end if
            end do

            !read tgain vector as string
            looptline: do
                read(10,"(a)",iostat=msg) text
                text=trim(adjustl(text))
                i=index(text,"#"); if(i>0) text=text(1:i-1) !remove comments after #
                i=index(text,"!"); if(i>0) text=text(1:i-1) !remove comments after !

                if (msg < 0) then !end of file
                    call error(file//' is missing tgain points.')
                endif
                if (msg > 0) then
                    call error('Check '//file//'.  Something is wrong..')
                endif
                if (text=='') then !blank line
                    cycle
                else 
                    exit looptline
                endif
            enddo looptline

            !convert string to real numbers
            read(text,*,iostat=msg) tgain(1:mx)

            if (any(tgain(1:mx)<0.)) then
                call error(file//' has unequal tgain & xgain pairs.')
            endif


            !read in gain value and taper
            loopgline: do
                read(10,"(a)",iostat=msg) text
                text=trim(adjustl(text))
                i=index(text,"#"); if(i>0) text=text(1:i-1) !remove comments after #
                i=index(text,"!"); if(i>0) text=text(1:i-1) !remove comments after !

                if (msg < 0) then !end of file
                    call error(file//' is missing tgain points.')
                endif
                if (msg > 0) then
                    call error('Check '//file//'.  Something is wrong..')
                endif
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
                if(shot%rcv(itr)%aoffset<=xgain(1)) then
                    jx1=1
                    jx2=1
                    dx1=0.
                    dx2=1.
                elseif(shot%rcv(itr)%aoffset>=xgain(mx)) then
                    jx1=mx
                    jx2=mx
                    dx1=1.
                    dx2=0.
                else
                    jx1=maxloc(xgain(1:mx), dim=1, mask=xgain(1:mx)<shot%rcv(itr)%aoffset)
                    jx2=jx1+1
                    dx1=(shot%rcv(itr)%aoffset-xgain(jx1)) / (xgain(jx2)-xgain(jx1))  !reduced distance from aoffset to xgain(jx1)
                    dx2=(xgain(jx2)-shot%rcv(itr)%aoffset) / (xgain(jx2)-xgain(jx1))  !reduced distance from aoffset to xgain(jx2)
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

    subroutine by_table(weight,file)
        real,dimension(:,:) :: weight
        character(*) :: file
        
        character(i_str_slen) :: text
        character(1)  :: delim
        integer,parameter :: max_gain_length=20 !maximum number of points per line
        real,dimension(max_gain_length) :: xgain, tgain, gain
        real,dimension(:,:),allocatable :: table, tmp_table

        real,dimension(:,:),allocatable :: weight

        call alloc(weight,self%nt,self%ntr)

        open(10,file=file,action='read')

        xgain=-99999.

        !loop to read xgain vector as string
        loopxline: do
            read(10,"(a)",iostat=msg) text
            text=trim(adjustl(text))
            i=index(text,"#"); if(i>0) text=text(1:i-1) !remove comments after #
            i=index(text,"!"); if(i>0) text=text(1:i-1) !remove comments after !

            if (msg < 0) then !end of file
                call error(file//' is missing xgain points.')
            endif
            if (msg > 0) then
               call error('Check file '//file//'.  Something is wrong..')
            endif
            if (text=='') then !blank line
                cycle
            else
                exit loopxline
            endif
        enddo loopxline

        !convert string to real numbers
        read(text,*,iostat=msg) xgain
        mx=count(xgain>=0.) !number of input xgain

        !xgain should be increasing
        do i=1,mx-1
            if(xgain(i+1)<xgain(i)) then
                call error('xgain from '//file//' is NOT increasing!')
            endif
        end do

        allocate(tmp_table(max_gain_length,mx))

        !initialize
        tgain=-99999.
        tmp_table=-99999.

        !loop to read tgain and gain vector as string
        k=1
        loopfile: do
            
            !initialize
            gain=-99999.

            looptline: do           
                read(10,"(a)",iostat=msg) text
                text=trim(adjustl(text))
                i=index(text,"#"); if(i>0) text=text(1:i-1) !remove comments after #
                i=index(text,"!"); if(i>0) text=text(1:i-1) !remove comments after !

                if (msg < 0) then !end of file
                    call error(file//' is missing tgain points.')
                endif
                if (msg > 0) then
                    call error('Check file '//file//'. Something is wrong..')
                endif
                if (text=='') then !blank line
                    cycle
                else
                    exit looptline
                endif
            enddo looptline

            if (text=='end' .or. text=='END') exit loopfile

            !convert string to real numbers
            read(text,*,iostat=msg) tgain(k), delim , gain

            if (k>1) then; if(tgain(k)<tgain(k-1)) then
                call error('Please make sure tgain in '//file//' is increasing!')
            endif; endif
            if ( count(gain>=0.) < mx ) then
                call error(file//' has gain vector length < xgain vector length.')
            endif

            tmp_table(k,:) = gain(1:mx)
            k=k+1

        enddo loopfile

        close(10)

        !reshape table without -99999 elements
        mt=count(tmp_table(:,1)>=0.)  !number of input tgain
        allocate(table(mt,mx))
        table=tmp_table(1:mt,1:mx)
        deallocate(tmp_table)

        !map table to weight with interpolation
        do itr=1,ntr

            !find the two xgain's closest to aoffset
            if(shot%rcv(itr)%aoffset<=xgain(1)) then
                jx1=1
                jx2=1
                dx1=0.
                dx2=1.
            elseif(shot%rcv(itr)%aoffset>=xgain(mx)) then
                jx1=mx
                jx2=mx
                dx1=1.
                dx2=0.
            else
                jx1=maxloc(xgain(1:mx), dim=1, mask=xgain(1:mx)<shot%rcv(itr)%aoffset)
                jx2=jx1+1
                dx1=(shot%rcv(itr)%aoffset-xgain(jx1)) / (xgain(jx2)-xgain(jx1))  !reduced distance from aoffset to xgain(jx1)
                dx2=(xgain(jx2)-shot%rcv(itr)%aoffset) / (xgain(jx2)-xgain(jx1))  !reduced distance from aoffset to xgain(jx2)
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

		        weight(it,itr) = 1.* ( &
                      ( table(jt1,jx1)*dt2 + table(jt2,jx1)*dt1 ) *dx2  &
                    + ( table(jt1,jx2)*dt2 + table(jt2,jx2)*dt1 ) *dx1  )

            enddo

        enddo

        deallocate(table)

    end subroutine

    subroutine by_custom(weight,file)
        real,dimension(:,:) :: weight
        character(*) :: file

        integer :: file_size

        inquire(file=file,size=file_size)
        if(file_size<4*nt*ntr) call warn('weight cumstom has size '//num2str(file_size/4)//' < nt*ntr. Some traces will not be weighted.')

        call sysio_read(file,weight,size(weight))

    end subroutine

    subroutine by_pershot(weight,file)
        real,dimension(:,:) :: weight
        character(*) :: file

        integer :: file_size

        inquire(file=file//shot%sindex(5:8),size=file_size)
        if(file_size<4*nt*ntr) call warn('weight per shot file '//file//shot%sindex(5:8)//' has size '//num2str(file_size/4)//' < nt*ntr. Some traces will not be weighted.',mpiworld%iproc)

        call sysio_read(file//shot%sindex(5:8),weight,size(weight))

    end subroutine


end
