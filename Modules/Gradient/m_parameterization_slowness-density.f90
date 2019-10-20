module m_parameterization
use m_sysio
use m_arrayop
use m_model
use m_field, only:waveeq_info
use m_gradient, only: gradient


    !PARAMETERIZATION   -- ALLOWED PARAMETERS
    !slowness-impedance -- sp sps rho

    !acoustic:
    !kpa = rho*vp^2 = rho*sp^-2
    !rho0= rho
    !gsp = gkpa*(-2)rho*vp^3
    !grho= gkpa*vp^2 + grho0

    !acoustic + gardner:
    !kpa = a*vp^(b+2) = a*sp^(-b-2)
    !rho0= a*vp^b     = a*sp^-b
    !gsp = (gkpa*(b+2)/b + grho0*vp^-2)*(-ab)*vp^(b+3)

    !P-SV:
    !lda = rho(vp^2-2vs^2) = rho*sp^-2*(1-2sps^2)
    !mu  = rho*vs^2        = rho*(sps/sp)^2
    !rho0= rho
    !gsp = (glda*vp^3 + (-2glda+gmu)vp*vs^2)*(-2)rho
    !gsps= (-2glda + gmu)*2rho*vp*vs
    !grho= glda*vp^2 + (-2glda+gmu)*vp*vs + grho0

    !P-SV + gardner:
    !lda = a*vp^(b+2) - 2a*vp^b*vs^2 = a*sp^(-b-2)*(1-2sps^2)
    !mu  = a*vp^b*vs^2               = a*sp^(-b-2)*sps^2
    !rho0= a*vp^b                    = a*sp^-b
    !gsp = (glda*vp^2 + (-2glda + gmu)vs^2 + grho0*b(b+2))*-a/(b+2)*vp^(b+1)
    !gsps= (-2glda + gmu)*2a*vp^(b+1)*vs


    public  parameterization, npar
    private 
    
    character(:),parameter :: parameterization='slowness-density'

    character(:),allocatable :: empirical, parlist
    character(3),dimension(3) :: pars !max 3 active parameters, max 3 letters for each
    real,dimension(3) :: pars_min, pars_max
    integer :: npar=0

    !real,dimension(3) :: hyper=0. !max 3 hyper-parameters for empirical law
    real :: a,b

    logical :: if_empirical=.false.
    logical :: if_gardner=.false., if_castagna=.false.

    contains
    
    subroutine init_parameterization
        
        !read in empirical law
        empirical=get_setup_char('EMPIRICAL_LAW')
        if_empirical = len(empirical)>0
        if(if_empirical) call read_empirical

        !read in active parameters and their allowed ranges
        parlist=get_setup_char('ACTIVE_PARAMETER',default='sp2.94e-4:6.67e-4')
        call read_parlist

    end subroutine

    subroutine read_empirical
        
        !Gardner law rho=a*vp^b
        !passive rho will be updated according to vp
        !https://wiki.seg.org/wiki/Dictionary:Gardner%E2%80%99s_equation
        !http://www.subsurfwiki.org/wiki/Gardner%27s_equation
        !https://en.wikipedia.org/wiki/Gardner%27s_relation
        if_gardner = index(empirical,'gardner')>0
        if(if_gardner) then
            a=310; b=0.25 !modify if you want
        endif

        ! !Castagna mudrock line vp=a*vp+b
        ! !passive vs will be updated according to vp
        ! !https://en.wikipedia.org/wiki/Mudrock_line
        ! if_castagna = index(empirical,'castagna')>0
        ! if(if_castagna) then
        !     c=1/1.16; d=-1.36/1.16 !modify if you want
        ! endif

    end subroutine
    
    subroutine read_parlist
        character(:),allocatable :: text
        !e.g. text='vp1500:3400', will output:
        !pars='vp', pars_min=1500, pars_max=3400
        !  k1=2,    k2=7,   k3=11

        npar=0

        loop: do while (npar<=3) !maximum 3 parameters
            
            if(len(parlist)==1) exit loop !no more parameters to read

            k3=index(parlist,' ') !parameters are separated by spaces
            if(k3==0) k3=len(parlist) !k3 can be 0 when looking at last parameter

            text=parlist(1:k3)
            if(len(text)==1) exit loop !no more parameters to read

            !ignore parameters depending on parameterization and empirical law
            if(.not. check_par(text) ) then
                !truncate parlist
                parlist=trim(adjustl(parlist(k3:len(parlist))))

                cycle loop
            endif

            !counter of valid active parameters
            npar=npar+1

            !read in parameter
            pars(npar)=text(1:2); k1=2
            if (text(1:3)=='sps' .or. text(1:3)=='rho') then
                pars(npar)=text(1:3); k1=3
            endif
        
            !read in parameter range
            k2=index(text,':')
            read(text(k1+1:k2-1),*) pars_min(npar)
            read(text(k2+1:k3  ),*) pars_max(npar)

            ! if(pars(i)=='vp') then
            !         par_vp_min=pars_min(i)
            !         par_vp_max=pars_max(i)
            ! endif

            !truncate parlist
            parlist=trim(adjustl(parlist(k3:len(parlist))))

            if(len(parlist)==1) exit loop !fin du parlist

        enddo loop


        if(mpiworld%is_master) then
            write(*,*) 'Number of inversion parameters:', npar
            write(*,*) 'Read parameters:',pars(1),' ',pars(2),' ',pars(3)
        endif

        if(npar==0) then
            stop 'ERROR: npar=0'
        endif

    end subroutine

    function check_par(text)
        character(*) :: text
        logical :: check_par
        
        check_par=.false.
        
        if(text(1:2)=='sp' ) check_par = .true.

        if(text(1:2)=='sps') check_par = index(waveeq_info,'elastic')>0
        
        if(text(1:3)=='rho') check_par = .not. if_gardner

        if(.not. check_par) then !ce paramètre n'est plus actif
            if(mpiworld%is_master) write(*,*) 'Parameter ',text(1:3),' will not be active in the inversion..'
        endif
        
    end function
    
    subroutine parameterization_transform(dir,x,g)
        character(3) :: dir
        real,dimension(m%nz,m%nx,m%ny,npar) :: x
        real,dimension(m%nz,m%nx,m%ny,npar),optional :: g

        logical :: if_acoustic, if_elastic
        
        !model
        if(dir=='m2x') then
            do ipar=1,npar
                select case (pars(ipar))
                case ('sp' ); x(:,:,:,ipar) = 1/m%vp
                case ('sps'); x(:,:,:,ipar) = m%vs/m%vp
                case ('rho'); x(:,:,:,ipar) = m%rho
                end select

                !normalize x to be unitless
                x(:,:,:,ipar)=(x(:,:,:,ipar)-pars_min(ipar))/(pars_max(ipar)-pars_min(ipar))
            enddo

        else !x2m
            !first run
            do ipar=1,npar
                select case (pars(ipar))
                case ('sp' ); m%vp = 1./(x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar))
                case ('rho'); m%rho=     x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar)
                end select
            enddo

            !second run
            do ipar=1,npar
                if(pars(ipar)=='sps') then
                    m%vs = (x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar)) *m%vp  !m%vp should have been updated in the first run
                endif
            enddo

        endif


        !gradient
        !!for units of gradient and g, see m_field*.f90
        if(present(g)) then
        if(dir=='m2x') then

            if_acoustic = index(waveeq_info,'acoustic')>0
            if_elastic  = index(waveeq_info,'elastic')>0

            !acoustic
            if(if_acoustic .and. .not. if_empirical) then
                do ipar=1,npar
                    select case (pars(ipar))
                    case ('sp' ); g(:,:,:,ipar) = gradient(:,:,:,1)*(-2)*m%rho*m%vp**3
                    case ('rho'); g(:,:,:,ipar) = gradient(:,:,:,1)*m%vp**2 + gradient(:,:,:,2)
                    end select
                enddo
            endif

            !acoustic + gardner
            if(if_acoustic .and. if_gardner) then
                g(:,:,:,1) =(gradient(:,:,:,1)*(b+2)/b + gradient(:,:,:,2)/m%vp**2)*(-a)*b*m%vp**(b+3)
            endif

            !elastic
            if(if_elastic .and. .not. if_empirical) then
                do ipar=1,npar
                    select case (pars(ipar))
                    case ('sp' ); g(:,:,:,ipar) =(gradient(:,:,:,1)*m%vp**3 + (-2*gradient(:,:,:,1)+gradient(:,:,:,2)*m%vp*m%vs**2))*(-2)*m%rho
                    case ('sps'); g(:,:,:,ipar) =(-2*gradient(:,:,:,1) + gradient(:,:,:,2))*2*m%rho*m%vp*m%vs
                    case ('rho'); g(:,:,:,ipar) = gradient(:,:,:,1)*m%vp**2 + (-2*gradient(:,:,:,1)+gradient(:,:,:,2))*m%vp*m%vs + gradient(:,:,:,3)
                    end select
                enddo
            endif

            !elastic + gardner
            if(if_elastic .and. if_gardner) then
                do ipar=1,npar
                    select case (pars(ipar))
                    case ('sp' ); g(:,:,:,ipar) =(gradient(:,:,:,1)*m%vp**2 + (-2*gradient(:,:,:,1)+gradient(:,:,:,2))*m%vs**2 + gradient(:,:,:,3)*b*(b+2))*(-a)/(b+2)*m%vp**(b+1)
                    case ('sps'); g(:,:,:,ipar) =(-2*gradient(:,:,:,1) + gradient(:,:,:,2))*2*a*m%vp**(b+1)*m%vs
                    end select
                enddo
            endif


            !normaliz g to be unitless
            do ipar=1,npar
                g(:,:,:,ipar)=g(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))
            enddo

        endif
        endif
        
    end subroutine
    
end module
