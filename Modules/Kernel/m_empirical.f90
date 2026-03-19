module m_empirical
use m_System
use m_pseudotime
use m_Modeling
    
    private
    public :: empirical_init, empirical_x2m, empirical_m2x, empirical_gradient

    logical,public :: is_empirical=.false. !needed by m_parametrizer
    logical :: is_gardner=.false., is_castagna=.false., is_abdullah=.false., is_vpqp=.false., is_vpvs=.false.

    real :: const_a,const_b,const_k,const_p
    real :: a_vpqp

    contains
    
    subroutine empirical_init()

        type(t_string),dimension(:),allocatable :: list,sublist

        !print info
        call hud('Empirical laws module is invoked.')
        call hud('Available empirical law: Gardner, Castagna, VpQp.')
        call hud('  - If Gardner law is used, rho from PARAMETER setup (if loaded) will be neglected and become passive.')
        call hud('  - If Castagna law is used, vs from PARAMETER setup (if loaded) will be neglected and become passive.')
        call hud('  - If Abdullah law is used, vs from PARAMETER setup (if loaded) will be neglected and become passive.')
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
                        const_k=either(0.31,310.,m%rho(1,1,1)<1000.) !g/cm³ or kg/m³
                        const_p=0.25
                    else
                        sublist=split(list(i)%s,o_sep=',') !ifort generates 'catastropic error ... internal compiler error' and lets me report...
                        const_k=str2real(sublist(2)%s)
                        const_p=str2real(sublist(3)%s)
                    endif

                    call hud('Gardner law is enabled: k='//num2str(const_k)//either(' g/cm³',' kg/m³',m%rho(1,1,1)<1000.)//', p='//num2str(const_p))

                elseif(list(i)%s(1:8)=='Castagna') then
                    !Castagna mudrock line vs=a*vp+b
                    !passive vs will be updated according to vp
                    !https://en.wikipedia.org/wiki/Mudrock_line
                    !Note that Poisson solid vs=a*vp
                    !can be realized by setting a=1./sqrt(3.) and b=0.
                    is_castagna=.true.
                    
                    if(len(list(i)%s)<=8) then
                        const_a=1/1.16
                        const_b=-1360./1.16 !m/s
                    else
                        sublist=split(list(i)%s,o_sep=',') !ifort generates 'catastropic error ... internal compiler error' and lets me report...
                        const_a=str2real(sublist(2)%s)
                        const_b=str2real(sublist(3)%s)
                    endif

                    call hud('Castagna law is enabled: a='//num2str(const_a)//', b='//num2str(const_b)//' m/s')

                elseif(list(i)%s(1:8)=='Abdullah') then
                    !Abdullah = Castagna + Gardner
                    !vs=a*vp+b, rho=a*vp^b
                    !passive vs will be updated according to vp & ip
                    is_abdullah=.true.

                    if(len(list(i)%s)<=8) then
                        const_a=1/1.16
                        const_b=-1360./1.16 !m/s
                        const_k=either(0.31,310.,m%rho(1,1,1)<1000.) !g/cm³ or kg/m³
                        const_p=0.25
                    else
                        sublist=split(list(i)%s,o_sep=',') !ifort generates 'catastropic error ... internal compiler error' and lets me report...
                        const_a=str2real(sublist(2)%s)
                        const_b=str2real(sublist(3)%s)
                        const_k =str2real(sublist(4)%s)
                        const_p =str2real(sublist(5)%s)
                    endif

                    call hud('Abdullah law is enabled: '//num2str(const_a)//', '//num2str(const_b)//', ' &
                                                        //num2str(const_k)//', '//num2str(const_p))

                elseif(list(i)%s(1:4)=='VpQp') then
                    is_vpqp=.true.
                    
                    call hud('Qp=sqrt(Vp) law is enabled')

                elseif(list(i)%s(1:4)=='VpVs') then
                    is_vpvs=.true.
                    
                    sublist=split(list(i)%s,o_sep=':')
                    a_vpqp=str2real(sublist(2)%s)
                    
                    call hud('Vp/Vs ='//num2str(a_vpqp)//' is enforced')

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
            if(is_gardner)  m%rho = const_k*m%vp**const_p
            if(is_castagna) m%vs  = const_a*m%vp + const_b
            if(is_vpqp)     m%qp  = sqrt(m%vp)
        endif

        if(parametrization=='velocities-impedance') then
            if(is_gardner)  m%rho = const_k*m%vp**const_p
            if(is_castagna) m%vs  = const_a*m%vp + const_b
        endif

    end subroutine

    subroutine empirical_gradient(parametrization,o_gvp,o_gvs,o_grho,o_gip,o_gqp,  o_ip)
        character(*) :: parametrization
        real,dimension(:,:,:),optional :: o_gvp,o_gvs,o_grho,o_gip,o_gqp
        real,dimension(:,:,:),optional :: o_ip

        real,dimension(:,:,:),allocatable :: v_t !velocity model in pseudotime domain

        if(.not.is_empirical) return

        if(parametrization=='velocities-density') then

            !Gardner
            if(is_gardner) then
                o_gvp  = o_gvp + o_grho * const_k*const_p*m%vp**(const_p-1)
            endif

            !Castagna
            if(is_castagna) then
                o_gvp = o_gvp + o_gvs * const_a
            endif

            if(is_vpqp) then
                !to be completed
            endif
            
        endif

        if(parametrization=='velocities-impedance') then

            !Gardner
            if(is_gardner) then
                o_gvp  = o_gvp + o_gip * a*(const_p+1)*m%vp**const_p
            endif

            !Castagna
            if(is_castagna) then
                o_gvp = o_gvp + o_gvs * const_a
            endif

            !Abdullah
            if(is_abdullah) then
                o_gip =         o_gvs* (const_a*m%vp+const_b)*const_k*m%vp**(const_p+1)/(-o_ip**2) + o_grho/m%vp
                o_gvp = o_gvp + o_gvs* (const_a*const_k*(const_p+2)  *m%vp**(const_p+1) + const_b*const_k*(const_p+1)*m%vp**const_p)/o_ip &
                              + o_grho* o_ip/(-m%vp**2)
            endif

            !Vp/Vs=a
            if(is_vpvs) then
                o_gvp = o_gvp + o_gvs / a_vpqp
            endif
            
        endif

        if(parametrization=='velocities-impedance_pseudotime') then

            !Gardner
            if(is_gardner) then
                call pseudotime_convert('z->t',m%vp,v_t)
                o_gvp  = o_gvp + o_gip * const_k*(const_p+1)*v_t**const_p
                deallocate(v_t)
            endif

            !Castagna
            if(is_castagna) then
                o_gvp = o_gvp + o_gvs * const_a
            endif


            !Abdullah
            if(is_abdullah) then
                o_gip =         o_gvs* (const_a*m%vp+const_b)*const_k*m%vp**(const_p+1)/(-o_ip**2) + o_grho/m%vp
                o_gvp = o_gvp + o_gvs* (const_a*const_k*(const_p+2)  *m%vp**(const_p+1) + const_b*const_k*(const_p+1)*m%vp**const_p)/o_ip &
                              + o_grho* o_ip/(-m%vp**2)
            endif

            !Vp/Vs=a
            if(is_vpvs) then
                o_gvp = o_gvp + o_gvs / a_vpqp
            endif
            
        endif

    end subroutine

end module
