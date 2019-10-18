module m_parameterization
use m_sysio
use m_arrayop
use m_model
use m_field, only:waveeq_info
use m_gradient, only: gradient


    !SITUATION    -- PARAMETERIZATION     -- ALLOWED PARAMETERS
    !WaveEquation -- moduli-density       -- kpa rho (acoustic) or lda mu rho (elastic)
    !I/O models   -- velocities-density   -- vp vs rho
    !seismic      -- velocities-impedance -- vp vs ip
    !tomography   -- slowness-density     -- sp ss rho
    !optimization -- can be any of above

    !Physical meaning of symbols uniformly used in this code:
    !vp                   -- P-wave velocity
    !vs                   -- S-wave velocity
    !sp =1/vp             -- P-wave slowness
    !sps=vs/vp            -- inv Vp-Vs ratio
    !ip =vp*rho           -- P-wave (acoustic) impedance
    !rho                  -- density
    !kpa=rho*vp^2         -- bulk modulus
    !lda=rho*(vp^2-2vs^2) -- 1st Lamé parameter
    !mu =rho*vs^2         -- 2nd Lamé parameter, shear modulus

    !source/destination trilogy:
    !model m =FWD=> WaveEq lm =PARAMETERIZATION=> optim x,g =FWD=> model m

    !active parameters will be converted to optim x via feature scaling
    !(e.g. (par-par_min)/(par_max-par_min))
    !
    !passive parameters will NOT be converted to optim x, but will be updated
    !according to user-specified empirical law (e.g. Gardner)
    !(user has to modified the code where necessary to insert such a law
    ! as symbolic computation is not straightforward in fortran..)
    !
    !both has influence on the conversion of gradient from WaveEq lm to optim x,g


    public
    private :: passive, read_passive
    private ::  active, read_active, check_par
    
    character(:),allocatable :: parameterization, passive, active

    character(3),dimension(3) :: pars !max 3 active parameters, max 3 letters for each
    real,dimension(3) :: pars_min, pars_max
    integer :: npar=0

    real,dimension(3) :: hyper=0. !max 3 hyper-parameters for passive parameters

    logical :: if_gardner=.false.


    contains
    
    subroutine init_parameterization
                
        parameterization=get_setup_char('PARAMETERIZATION',default='velocities-density')

        !read in passive parameters which depends on active parameters
        passive=trim(adjustl(get_setup_char('PASSIVE_PARAMETER')))
        if(len(passive)>0) call read_passive

        !read in active parameters and their allowed ranges
        active=get_setup_char('ACTIVE_PARAMETER',default='vp1500:3400')
        active=trim(adjustl(active)) !no head or tail spaces
        call read_active

        !expect rational users..
        ! !resolve conflict of parameters e.g. vp & sp
        !
        ! !reallocate model if initial model is not provided
        ! if(if_par_rho.or.if_par_kpa.or.if_par_ip) then
        !     call model_reallocate('rho')
        ! endif

    end subroutine

    subroutine read_passive
        !gardner law
        !passive rho will be updated according to vp
        if_gardner = index(passive,'gardner')>0
        if(if_gardner) then
            hyper(1)=310; hyper(2)=0.25 !http://www.subsurfwiki.org/wiki/Gardner%27s_equation
            !modify if you want
        endif

    end subroutine
    
    subroutine read_active
        character(:),allocatable :: text
        !e.g. text='vp1500:3400', will output:
        !pars='vp', pars_min=1500, pars_max=3400
        !  k1=2,    k2=7,   k3=11

        npar=0

        loop: do while (npar<=3) !maximum 3 parameters
            
            if(len(active)==1) exit loop !no more parameters to read

            k3=index(active,' ') !parameters are separated by spaces
            if(k3==0) k3=len(active) !k3 can be 0 when looking at last parameter

            text=active(1:k3)
            if(len(text)==1) exit loop !no more parameters to read

            !ignore parameters depending on parameterization and passive
            if(.not. check_par(text) ) then
                !truncate active
                active=trim(adjustl(active(k3:len(active))))

                cycle loop
            endif

            !counter of valid active parameters
            npar=npar+1

            !read in parameter
            pars(npar)=text(1:2); k1=2
            if (text(1:3)=='sps'.or.text(1:3)=='rho'.or.text(1:3)=='kpa'.or.text(1:3)=='lda') then
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

            !truncate active
            active=trim(adjustl(active(k3:len(active))))

            if(len(active)==1) exit loop !fin du active

        enddo loop


        if(mpiworld%is_master) then
            write(*,*) 'Number of inversion parameters:', npar
            write(*,*) 'Read parameters:',pars(1),' ',pars(2),' ',pars(3)
        endif

    end subroutine

    function check_par(text)
        character(*) :: text
        logical :: check_par
        
        check_par=.false.
        
        if    (text(1:2)=='vp' ) then !P-wave velocity
            check_par = index(parameterization,'velocities')>0
            
        elseif(text(1:2)=='vs' ) then !S-wave velocity
            check_par = index(parameterization,'velocities')>0 .and. index(waveeq_info,'elastic')>0
            
        elseif(text(1:2)=='sp' ) then !P-wave slowness
            check_par = index(parameterization,'slowness')>0

        elseif(text(1:3)=='sps') then !inv Vp-Vs ratio
            check_par = index(parameterization,'slowness')>0  .and. index(waveeq_info,'elastic')>0
            
        elseif(text(1:2)=='ip' ) then !P-wave (acoustic) impedance
            check_par = index(parameterization,'impedance')>0
            
        elseif(text(1:3)=='rho') then !density
            check_par = index(parameterization,'density')>0 .and. .not. if_gardner
            
        elseif(text(1:3)=='kpa') then !bulk modulus
            check_par = index(parameterization,'moduli')>0 .and. index(waveeq_info,'acoustic')>0
            
        elseif(text(1:3)=='lda') then !1st Lamé parameter
            check_par = index(parameterization,'moduli')>0 .and. index(waveeq_info,'elastic')>0
            
        elseif(text(1:2)=='mu' ) then !2nd Lamé parameter, shear modulus
            check_par = index(parameterization,'moduli')>0 .and. index(waveeq_info,'elastic')>0
        
        endif

        if(.not. check_par) then !ce paramètre n'est plus actif
            if(mpiworld%is_master) write(*,*) 'Parameter ',text(1:3),' will not be active in the inversion..'
        endif
        
    end function
    
    subroutine parameterization_transform(dir,x)
        character(3) :: dir
        real,dimension(m%nz,m%nx,m%ny,npar) :: x
        
        if(dir=='m2x') then
            do ipar=1,npar
                select case (pars(ipar))
                case ('vp' ); x(:,:,:,ipar) = m%vp
                case ('vs' ); x(:,:,:,ipar) = m%vs
                case ('sp' ); x(:,:,:,ipar) = 1./m%vp
                case ('sps'); x(:,:,:,ipar) = m%vs/m%vp
                case ('ip' ); x(:,:,:,ipar) = m%vp*m%rho
                case ('rho'); x(:,:,:,ipar) = m%rho
                case ('kpa'); x(:,:,:,ipar) = m%rho* m%vp*m%vp
                case ('lda'); x(:,:,:,ipar) = m%rho*(m%vp*m%vp-2.*m%vs*m%vs)
                case ('mu' ); x(:,:,:,ipar) = m%rho*(             m%vs*m%vs)
                end select

                !normalize x to be unitless
                x(:,:,:,ipar)=(x(:,:,:,ipar)-pars_min(ipar))/(pars_max(ipar)-pars_min(ipar))
            enddo

        else !x2m
            !first run
            do ipar=1,npar
                select case (pars(ipar))
                case ('vp' ); m%vp = (x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar))
                case ('vs' ); m%vs = (x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar))
                case ('sp' ); m%vp = 1. / (x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar))
                case ('rho'); m%rho= (x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar))
                end select
            enddo

            !second run
            do ipar=1,npar
                select case (pars(ipar))
                case ('sps'); m%vs = (x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar)) * m%vp !m%vp should have been updated after first run
                case ('ip' ); m%rho= (x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar)) / m%vp
                case ('kpa'); m%vp = sqrt( (x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar)) /m%rho )
                case ('lda'); m%vp = sqrt( (x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar)) /m%rho +2*m%vs*m%vs )
                case ('mu' ); m%vs = sqrt( (x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))+pars_min(ipar)) /m%rho )
                end select
            enddo

        endif

    end subroutine

    subroutine parameterization_transform_gradient(dir,g)
        character(3) :: dir
        real,dimension(m%nz,m%nx,m%ny,npar) :: g  !for units of gradient and g, see m_field*.f90
        
        if(dir=='m2x') then

            select case (parameterization)

            case ('moduli-density')
                do ipar=1,npar
                    select case (pars(ipar))

                    case ('kpa')
                        g(:,:,:,ipar)=gradient(:,:,:,1)

                    case ('lda')
                        g(:,:,:,ipar)=gradient(:,:,:,1)

                    case ('mu')
                        g(:,:,:,ipar)=gradient(:,:,:,2)

                    case ('rho')
                        
                        if(index(waveeq_info,'acoustic')>0) then
                            g(:,:,:,ipar)=gradient(:,:,:,2)
                        else !elastic
                            g(:,:,:,ipar)=gradient(:,:,:,3)
                        endif

                    end select
                enddo

            case ('velocities-density')
                do ipar=1,npar
                    select case (pars(ipar))

                    case ('vp')
                        if(index(waveeq_info,'acoustic')>0) then
                            if(if_gardner) then
                                g(:,:,:,ipar)=(gradient(:,:,:,1)*2.*m%vp*m%rho)
                            else
                            endif
                        endif
                        if(index(waveeq_info,'elastic')>0) then
                            if(if_gardner) then
                            else
                            endif
                        endif

                    case ('vs')

                    case ('rho')

                    end select
                enddo

            case ('velocities-impedance')
                do ipar=1,npar
                    select case (pars(ipar))

                    case ('vp')
                        if(index(waveeq_info,'acoustic')>0) then
                            g(:,:,:,1)=(gradient(:,:,:,1)*2.*m%vp*m%rho)
                        else !elastic
                        endif

                    case ('vs')

                    case ('ip')

                    end select
                enddo

            case ('slowness-density')
                do ipar=1,npar
                    select case (pars(ipar))

                    case ('sp')
                        if(index(waveeq_info,'acoustic')>0) then
                            g(:,:,:,1)=(gradient(:,:,:,1)*2.*m%vp*m%rho)
                        else !elastic
                        endif

                    case ('sps')

                    case ('rho')

                    end select
                enddo

            end select

            !normaliz g to be unitless
            do ipar=1,npar
                g(:,:,:,ipar)=g(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar))
            enddo


    !         select case (parameterization)
    !             case ('velocities-density') !unit of g = [kg/m2/s2]
                    
    !                     g(:,:,:,1)=(gradient(:,:,:,1)*2.*m%vp*m%rho)               *(par_vp_max-par_vp_min)
    !                 if(if_par_rho) g(:,:,:,2)=(gradient(:,:,:,1)*m%vp*m%vp+gradient(:,:,:,2)) *(par_rho_max-par_rho_min)
    ! !             case ('vprho-gardner')
    !             case ('vp-ip')
    !                 if(if_par_vp) g(:,:,:,1)=(gradient(:,:,:,1)*m%vp*m%rho - gradient(:,:,:,2)*m%rho/m%vp) *(par_vp_max-par_vp_min)
    !                 if(if_par_ip) g(:,:,:,2)=(gradient(:,:,:,1)*m%vp       + gradient(:,:,:,2)/m%vp      ) *(par_ip_max-par_ip_min)
    !         end select
        endif
        
    end subroutine
    
end module
