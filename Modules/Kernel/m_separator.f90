! Modified from sumute.c in Seismic Unix
! dir: SeisUnix/src/su/main/windowing_sorting_muting

module m_separator
use m_System
use m_shot

    private

    type,public :: t_separator

        real,dimension(:,:),allocatable :: nearoffset, reflection, diving

        contains

        procedure :: update
        procedure :: by_aoffset
        procedure :: by_polygon

    end type

    type(t_separator),public :: sepa

    contains
    
    subroutine update(self)
        class(t_separator) :: self

        character(:),allocatable :: suf
        type(t_string),dimension(:),allocatable :: list,sublist,subsublist
        character(:),allocatable :: file

        call alloc(self%nearoffset,shot%nt,shot%nrcv,o_init=0.)
        call alloc(self%reflection,shot%nt,shot%nrcv,o_init=0.)
        call alloc(self%diving    ,shot%nt,shot%nrcv,o_init=0.)

        list=setup%get_strs('SEPARATING','SEPA',o_mandatory=2)

        do i=1,size(list)
            
            if (index(list(i)%s,'aoffset')>0) then !use aoffset for either Ip or Vp inversion, but not for both
                sublist=split(list(i)%s,o_sep=':')
                call hud('Will separate data by aoffset:'//sublist(2)%s)
                call self%by_aoffset(str2real(sublist(2)%s))
            endif
            
            if (index(list(i)%s,'polygon')>0) then !use polygon to separate diving vs reflections
                sublist=split(list(i)%s,o_sep=':')
                if(size(sublist)==1) then !filename is not attached, ask for it
                    file=setup%get_file('FILE_SEPARATOR_POLYGON','FILE_SEPA',o_mandatory=1)
                else !filename is attached
                    file=sublist(2)%s
                endif
                
                call hud('Will separate data with polygons defined in '//file)
                call self%by_polygon(file)
            endif

            if (index(list(i)%s,'d/r')>0) then !weighting on diving vs reflections
                sublist=split(list(i)%s,o_sep=':')
                call hud('Will separate data with weight diving/reflection= '//sublist(2)%s)
                subsublist = split(sublist(2)%s,o_sep='/')
                self%diving    =self%diving    *str2real(subsublist(1)%s)
                self%reflection=self%reflection*str2real(subsublist(2)%s)
            endif

        enddo

        if(mpiworld%is_master) then
            call suformat_write('sepa%nearoffset',self%nearoffset,shot%nt,shot%nrcv,o_dt=shot%dt)
            call suformat_write('sepa%reflection',self%reflection,shot%nt,shot%nrcv,o_dt=shot%dt)
            call suformat_write('sepa%diving'    ,self%diving    ,shot%nt,shot%nrcv,o_dt=shot%dt)
        endif
        
    end subroutine

    subroutine by_aoffset(self,aoffset)
        class(t_separator) :: self
        real :: aoffset

        do i=1,shot%nrcv    
            if (shot%rcv(i)%aoffset<=aoffset) self%nearoffset(:,i)=1.
        enddo

    end subroutine

    subroutine by_polygon(self,file)
        class(t_separator) :: self
        character(*) :: file

        character(i_str_slen) :: text
        integer,parameter :: max_sepa_length=20 !maximum number of points
        real,dimension(max_sepa_length) :: xsepa,tsepa
        
        open(10,file=file,action='read')

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
                call error(file//' is missing tsepa points.')
            endif
            if (msg > 0) call error('Check '//file//'.  Something is wrong..')
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
                call error('xsepa from '//file//' is NOT increasing!')
            end if
        end do


        !read tsepa vector as string
        looptline: do
            read(10,"(a)",iostat=msg) text
            text=trim(adjustl(text))
            i=index(text,"#"); if(i>0) text=text(1:i-1) !remove comments after #
            i=index(text,"!"); if(i>0) text=text(1:i-1) !remove comments after !

            if (msg < 0) then !end of file
                call error(file//' is missing tsepa points.')
            endif
            if (msg > 0) call error('Check '//file//'.  Something is wrong..')
            if (text=='') then !blank line
                cycle
            else 
                exit looptline
            endif
        enddo looptline

        !convert string to real numbers
        read(text,*,iostat=msg) tsepa(1:mx)

        if (any(tsepa(1:mx)<0.)) then
            call error(file//' has unequal tsepa & xsepa pairs.')
        endif


        !loop of offset
        do itr=1,shot%nrcv
            
            !find the two xsepa's closest to aoffset
            if(shot%rcv(itr)%aoffset<=xsepa(1)) then
                jx1=1
                jx2=1
                dx1=0.
                dx2=1.
            elseif(shot%rcv(itr)%aoffset>=xsepa(mx)) then
                jx1=mx
                jx2=mx
                dx1=1.
                dx2=0.
            else
                jx1=maxloc(xsepa(1:mx), dim=1, mask=xsepa(1:mx)<shot%rcv(itr)%aoffset)
                jx2=jx1+1
                dx1=(shot%rcv(itr)%aoffset-xsepa(jx1)) / (xsepa(jx2)-xsepa(jx1))  !reduced distance from aoffset to xsepa(jx1)
                dx2=(xsepa(jx2)-shot%rcv(itr)%aoffset) / (xsepa(jx2)-xsepa(jx1))  !reduced distance from aoffset to xsepa(jx2)
            endif
            
            !interp time by the two tsepa & xsepa pairs
            time = tsepa(jx1)*dx2 + tsepa(jx2)*dx1

            itime = nint(time/shot%dt)+1  !assume all receivers have same shot%dt ...
            itime = min(max(itime,1),shot%nt)

            self%diving(1:itime,itr)=1.
            self%reflection(itime+1:shot%nt,itr)=1.

        end do
            
        close(10)

    end subroutine

end
