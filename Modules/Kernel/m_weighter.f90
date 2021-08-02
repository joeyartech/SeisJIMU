module m_weighter
use m_System
use m_Modeling

    private

    type,public :: t_weighter

        integer nt,ntr
        real dt
        real,dimension(:,:),allocatable :: weight

        contains

        procedure :: init
        procedure :: by_components
        procedure :: by_aoffset
        procedure :: by_polygon
        procedure :: by_table
        procedure :: by_custom

    end type

    type(t_weighter),public :: wei
    
    contains

    subroutine init(self,o_nt,o_dt,o_ntr)
        class(t_weighter) :: self
        integer,optional :: o_nt, o_ntr
        real,optional :: o_dt

        type(t_string),dimension(:),allocatable :: list,sublist
        character(:),allocatable :: file

        self%nt = either(o_nt,  shot%nt  , present(o_nt ))
        self%ntr= either(o_ntr, shot%nrcv, present(o_ntr))
        self%dt = either(o_dt,  shot%dt  , present(o_dt ))

        call alloc(self%weight,self%nt,self%ntr,o_init=1.)

        list=setup%get_strs('WEIGHTING','WEI',o_default=either('aoffset^1','aoffset^0.5',m%is_cubic))

        do i=1,size(list)
            !weight pressure components
            if(list(i)%s(1:2)=='p*') then
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will multiply pressure component by '//sublist(2)%s)
                call self%by_components(o_p=str2real(sublist(2)%s))
            endif

            if(list(i)%s(1:3)=='vz*') then
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will multiply vz component by '//sublist(2)%s)
                call self%by_components(o_vz=str2real(sublist(2)%s))
            endif

            if(list(i)%s(1:3)=='vx*') then
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will multiply vx component by '//sublist(2)%s)
                call self%by_components(o_vx=str2real(sublist(2)%s))
            endif

            if(list(i)%s(1:3)=='vy*') then
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will multiply vy component by '//sublist(2)%s)
                call self%by_components(o_vy=str2real(sublist(2)%s))
            endif

            !weight traces based on power of aoffset
            if(list(i)%s(1:8)=='aoffset^') then
                sublist=split(list(i)%s,o_sep='^')
                call hud('Will weight traces by aoffset^'//sublist(2)%s)
                call self%by_aoffset(o_power=str2real(sublist(2)%s))
            endif

            !use aoffset to multiply on traces
            if(list(i)%s(1:8)=='aoffset*') then
                sublist=split(list(i)%s,o_sep='*')
                call hud('Will weight traces by aoffset*'//sublist(2)%s)
                call self%by_aoffset(o_factor=str2real(sublist(2)%s))
            endif

            !polygon defined weights to multiply on traces
            if(list(i)%s(1:7)=='polygon') then
                sublist=split(list(i)%s,o_sep=':')
                if(size(sublist)==1) then !filename is not attached, ask for it
                    file=setup%get_file('FILE_WEIGHT_POLYGON',o_mandatory=1)
                else !filename is attached
                    file=sublist(2)%s
                endif

                call hud('Will weight traces with polygons defined in '//file)
                call self%by_polygon(file)
            endif

            !table defined weights to multiply on traces
            if(list(i)%s(1:5)=='table') then
                sublist=split(list(i)%s,o_sep=':')
                if(size(sublist)==1) then !filename is not attached, ask for it
                    file=setup%get_file('FILE_WEIGHT_TABLE',o_mandatory=1)
                else !filename is attached
                    file=sublist(2)%s
                endif

                call hud('Will weight traces with tables defined in '//file)
                call self%by_table(file)

            endif

            if(list(i)%s(1:6)=='custom') then
                sublist=split(list(i)%s,o_sep=':')
                if(size(sublist)==1) then !filename is not attached, ask for it
                    file=setup%get_file('FILE_WEIGHT_CUSTOM',o_mandatory=1)
                else !filename is attached
                    file=sublist(2)%s
                endif

                call hud('Will weight traces in a custom way defined in '//file)
                call self%by_custom(file)

            endif

        enddo

        if(mpiworld%is_master) call suformat_write(dir_out//'weights',self%weight,self%nt,self%ntr,o_dt=self%dt)
        
    end subroutine

    subroutine by_components(self,o_p,o_vz,o_vx,o_vy)
        class(t_weighter) :: self
        real,optional :: o_p, o_vz, o_vx, o_vy

        if(present(o_p)) then
            do i=1,shot%nrcv
                if(shot%rcv(i)%comp=='p') self%weight(:,i)=o_p
            enddo
        endif

        if(present(o_vz)) then
            do i=1,shot%nrcv
                if(shot%rcv(i)%comp=='vz') self%weight(:,i)=o_vz
            enddo
        endif

        if(present(o_vx)) then
            do i=1,shot%nrcv
                if(shot%rcv(i)%comp=='vx') self%weight(:,i)=o_vx
            enddo
        endif

        if(present(o_vy)) then
            do i=1,shot%nrcv
                if(shot%rcv(i)%comp=='vy') self%weight(:,i)=o_vy
            enddo
        endif

    end subroutine

    subroutine by_aoffset(self,o_power,o_factor)
        class(t_weighter) :: self
        real,optional :: o_power, o_factor

        if(present(o_power)) then
            do i=1,shot%nrcv    
                self%weight(:,i)=self%weight(:,i)*shot%rcv(i)%aoffset**o_power
            enddo
        endif

        if(present(o_factor)) then
            do i=1,shot%nrcv    
                self%weight(:,i)=self%weight(:,i)*shot%rcv(i)%aoffset*o_factor
            enddo
        endif

    end subroutine

    ! Modified from sumute.c in Seismic Unix
    ! dir: SeisUnix/src/su/main/windowing_sorting_muting
    subroutine by_polygon(self,file)
        class(t_weighter) :: self
        character(*) :: file

        character(i_str_slen) :: text
        integer,parameter :: max_gain_length=20 !maximum number of points
        real,dimension(max_gain_length) :: xgain,tgain
        
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
                if (msg > 0) stop 'Check '//file//'.  Something is wrong..'
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
                if (msg > 0) stop 'Check '//file//'.  Something is wrong..'
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
                if (msg > 0) stop 'Check '//file//'.  Something is wrong..'
                if (text=='') then !blank line
                    cycle
                else 
                    exit loopgline
                endif
            enddo loopgline

            !convert string to real numbers
            read(text,*,iostat=msg) gain, taper
            ntaper=nint(taper/self%dt)

            !loop of offset
            do itr=1,self%ntr
                
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

                itime = nint(time/self%dt)+1  !assume all receivers have same dt ...
                itime_start = itime-ntaper  !1 <= itime_start <= itime <= nt

                !apply linear taper from itime_start+1 to itime
                do it=min(max(itime_start+1,1),self%nt) , &
                      min(max(itime        ,1),self%nt)
                    
                    k=it-itime_start-1

                    self%weight(it,itr) = self%weight(it,itr)* ( gain_old*(ntaper-k) + gain*k ) / ntaper
                end do               

                !apply constant gain from itime+1 to nt
                it = min(max(itime+1,1),self%nt)
                self%weight(it:self%nt,itr) = self%weight(it:self%nt,itr)* gain
                
            end do
            
            gain_old=gain

        end do loopfile

        close(10)

    end subroutine

    subroutine by_table(self,file)
        class(t_weighter) :: self
        character(*) :: file
        
        character(i_str_slen) :: text
        character(1)  :: delim
        integer,parameter :: max_gain_length=20 !maximum number of points per line
        real,dimension(max_gain_length) :: xgain, tgain, gain
        real,dimension(:,:),allocatable :: table, tmp_table

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
            if (msg > 0) stop 'Check file '//file//'.  Something is wrong..'
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
                if (msg > 0) stop 'Check file '//file//'. Something is wrong..'
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
        do itr=1,self%ntr

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


            do it=1,self%nt
                t=(it-1)*self%dt

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

                
                self%weight(it,itr) = self%weight(it,itr) * ( &
                      ( table(jt1,jx1)*dt2 + table(jt2,jx1)*dt1 ) *dx2  &
                    + ( table(jt1,jx2)*dt2 + table(jt2,jx2)*dt1 ) *dx1  )

            enddo

        enddo

        deallocate(table)

    end subroutine

    subroutine by_custom(self,file)
        class(t_weighter) :: self
        character(*) :: file

        integer :: file_size
        real,dimension(:,:),allocatable :: tmp

        inquire(file=file,size=file_size)
        if(file_size<4*self%nt*self%ntr) call warn('FILE_WEIGHT_CUSTOM has size '//num2str(file_size/4)//' < nt*ntr. Some traces will not be weighted.')

        call alloc(tmp,self%nt,self%ntr,o_init=1.)

        call sysio_read(file,tmp,size(tmp))

        self%weight=self%weight*tmp

        deallocate(tmp)

    end subroutine

end