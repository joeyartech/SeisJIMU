module m_empirical
use m_System
use m_pseudotime
use m_Modeling
    
    private
    public :: empirical_init, empirical_x2m, empirical_m2x, empirical_gradient

    logical,public :: is_empirical=.false. !needed by m_parametrizer
    logical :: is_gardner=.false., is_castagna=.false., is_vpqp=.false.

    real :: a,b

    contains
    
    subroutine empirical_init()

        type(t_string),dimension(:),allocatable :: list,sublist

        !print info
        call hud('Empirical laws module is invoked.')
        call hud('Available empirical law: Gardner, Castagna, VpQp.')
        call hud('  - If Gardner law is used, rho from PARAMETER setup (if loaded) will be neglected and become passive.')
        call hud('  - If Castagna law is used, vs from PARAMETER setup (if loaded) will be neglected and become passive.')
        call hud('  - If VpQp law is used, vs, qp from PARAMETER setup (if loaded) will be neglected and become passive.')

        !read in empirical law
        list=setup%get_strs('EMPIRICAL_LAW')

        is_empirical=size(list)>0

        if(is_empirical) then
            do i=1,size(list)
                if(list(i)%s(1:7)=='Gardner') then
                    !Gardner law rho=a*vp^b
                    !passive rho will be updated according to vp
                    !https://wiki.seg.org/wiki/Dictionary:Gardner%E2%80%99s_equation
                    !http://www.subsurfwiki.org/wiki/Gardner%27s_equation
                    !https://en.wikipedia.org/wiki/Gardner%27s_relation
                    is_gardner=.true.

                    if(len(list(i)%s)<=7) then
                        a=either(0.31,310.,m%rho(1,1,1)<1000.) !g/cm³ or kg/m³
                        b=0.25
                    else
                        sublist=split(list(i)%s,o_sep=',') !ifort generates 'catastropic error ... internal compiler error' and lets me report...
                        a=str2real(sublist(2)%s)
                        b=str2real(sublist(3)%s)
                    endif

                    call hud('Gardner law is enabled: a='//num2str(a)//either(' g/cm³',' kg/m³',m%rho(1,1,1)<1000.)//', b='//num2str(b))

                elseif(list(i)%s(1:8)=='Castagna') then
                    !Castagna mudrock line vs=a*vp+b
                    !passive vs will be updated according to vp
                    !https://en.wikipedia.org/wiki/Mudrock_line
                    !Note that Poisson solid vs=a*vp
                    !can be realized by setting a=1./sqrt(3.) and b=0.
                    is_castagna=.true.
                    
                    if(len(list(i)%s)<=8) then
                        a=1/1.16
                        b=-1360./1.16 !m/s
                    else
                        sublist=split(list(i)%s,o_sep=',') !ifort generates 'catastropic error ... internal compiler error' and lets me report...
                        a=str2real(sublist(2)%s)
                        b=str2real(sublist(3)%s)
                    endif

                    call hud('Castagna law is enabled: a='//num2str(a)//', b='//num2str(b)//' m/s')

                elseif(list(i)%s(1:4)=='VpQp') then
                    is_vpqp=.true.
                    
                    call hud('Qp=sqrt(Vp) law is enabled')

                endif

            enddo

        endif

        if(allocated(   list)) deallocate(   list)
        if(allocated(sublist)) deallocate(sublist)
        
    end subroutine
    
    subroutine empirical_m2x(parametrization)
        character(*) :: parametrization

        if(.not.is_empirical) return

        if(parametrization=='velocities-density') then

        endif

    end subroutine

    subroutine empirical_x2m(parametrization)
        character(*) :: parametrization

        if(.not.is_empirical) return

        if(parametrization=='velocities-density') then
            if(is_gardner)  m%rho = a*m%vp**b
            if(is_castagna) m%vs = a*m%vp + b
            if(is_vpqp)     m%qp = sqrt(m%vp)
        endif

        if(parametrization=='velocities-impedance') then
            if(is_gardner)  m%rho = a*m%vp**b
            if(is_castagna) m%vs = a*m%vp + b
        endif

    end subroutine

    subroutine empirical_gradient(parametrization,o_gvp,o_gvs,o_grho,o_gip,o_gqp)
        character(*) :: parametrization
        real,dimension(:,:,:),optional :: o_gvp,o_gvs,o_grho,o_gip,o_gqp

        real,dimension(:,:,:),allocatable :: v_t !velocity model in pseudotime domain

        if(.not.is_empirical) return

        if(parametrization=='velocities-density') then

            !Gardner
            if(is_gardner) then
                o_gvp  = o_gvp + o_grho * a*b*m%vp**(b-1)
                o_grho = 0.
            endif

            !Castagna
            if(is_castagna) then
                o_gvp = o_gvp + o_gvs * a
                o_gvs = 0.
            endif

            if(is_vpqp) then
                !to be completed
            endif
            
        endif

        if(parametrization=='velocities-impedance') then

            !Gardner
            if(is_gardner) then
                o_gvp  = o_gvp + o_gip * a*(b+1)*m%vp**b
                o_gip = 0.
            endif

            !Castagna
            if(is_castagna) then
                o_gvp = o_gvp + o_gvs * a
                o_gvs = 0.
           endif
            
        endif

        if(parametrization=='velocities-impedance_pseudotime') then

            !Gardner
            if(is_gardner) then
                call pseudotime_convert('z->t',m%vp,v_t)
                o_gvp  = o_gvp + o_gip * a*(b+1)*v_t**b
                o_gip = 0.
                deallocate(v_t)
            endif

            !Castagna
            if(is_castagna) then
                o_gvp = o_gvp + o_gvs * a
                o_gvs = 0.
           endif
            
        endif

    end subroutine

end module
