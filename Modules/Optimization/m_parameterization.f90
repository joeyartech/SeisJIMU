module m_parameterization
use m_sysio
use m_arrayop
use m_model
use m_field, only:waveeq_info
use m_gradient, only: gradient

    !PARAMETERIZATION                              !allowed ACTIVE_PARAMETER
    !WaveEquation:  moduli-density:                 kpa-rho or lda-mu-rho
    !I/O models:    velocities-density:             vp-rho or vp-vs-rho
    !Empirical:     velocities-density-empirical
    !seismic:       velocities-impedances           vp-ip or vp-vs-ip
    !tomography:    slowness-density                sp-rho
    !optimization:  can be any of above

    !optim x =FWD=> model m =FWD=> WaveEq lm =PARAMETERIZATION=> optim x,g

    !Symbol meaning:
    !vp    P-wave velocity
    !vs    S-wave velocity
    !sp    P-wave slowness
    !ip    P-wave (acoustic) impedance
    !rho   density
    !kpa   bulk modulus
    !lda   1st Lamé parameter
    !mu    2nd Lamé parameter, shear modulus

    public
    private :: active, read_active, check_par
    
    character(:),allocatable :: parameterization, active
    
    character(3),dimension(3) :: pars !max 3 letters for each par, max 3 parameters
    real,dimension(3) :: pars_min
    real,dimension(3) :: pars_max
    
    integer :: npar

    contains
    
    subroutine init_parameterization
                
        parameterization=get_setup_char('PARAMETERIZATION',default='velocities-density')

        active=get_setup_char('ACTIVE_PARAMETER',default='vp1500:3400')

        active=trim(adjustl(active)) !no head or tail spaces

        !read in parameters and their allowed ranges
        call read_active

        ! !adjust parameters according to parameterization
        ! call adjust_pars
        
        ! !reallocate model if initial model is not provided
        ! if(if_par_rho.or.if_par_kpa.or.if_par_ip) then
        !     call model_reallocate('rho')
        ! endif
        ! 
        ! npar=count([if_par_kpa,if_par_rho,if_par_vp,if_par_ip].eqv..true.)
        ! if(mpiworld%is_master) write(*,*) 'Number of inversion parameters:', npar

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
            if(len(text)==1) exit !no more parameters to read

            !ignore parameters depending on parameterization
            if(.not. check_par(text) ) then
                !truncate active
                active=trim(adjustl(active(k3:len(active))))

                cycle loop
            endif

            !counter of valid parameters
            npar=npar+1

            !read in parameter
            pars(npar)=text(1:2); k1=2
            if (text(1:3)=='rho'.or.text(1:3)=='kpa'.or.text(1:3)=='lda') then
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
        
        check_par=.true.
        
        if    (text(1:2)=='vp' ) then !P-wave velocity
            check_par = index(parameterization,'velocities')>0
            
        elseif(text(1:2)=='vs' ) then !S-wave velocity
            check_par = index(parameterization,'velocities')>0 .and. index(waveeq_info,'elastic')>0
            
        elseif(text(1:2)=='sp' ) then !P-wave slowness
            check_par = index(parameterization,'slowness')>0
            
        elseif(text(1:2)=='ip' ) then !P-wave (acoustic) impedance
            check_par = index(parameterization,'impedance')>0
            
        elseif(text(1:3)=='rho') then !density
            check_par = index(parameterization,'density')>0
            
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
    
    subroutine parameterization_transform(dir,x,g)
        character(3) :: dir
        real,dimension(m%nz,m%nx,m%ny,npar) :: x
        real,dimension(m%nz,m%nx,m%ny,npar),optional :: g
        
        !! Issue:
        !! what if I just want to invert for rho, ip ?
        !!
        
        !model
        if(dir=='m2x') then
            ! select case (parameterization)
            !     case ('vp-rho')
            !         if(if_par_vp)  x(:,:,:,1) = (m%vp -par_vp_min)  / (par_vp_max -par_vp_min)
            !         if(if_par_rho) x(:,:,:,2) = (m%rho-par_rho_min) / (par_rho_max-par_rho_min)
            !     case ('vp-ip')
            !         if(if_par_vp)  x(:,:,:,1) = (m%vp -par_vp_min)  / (par_vp_max -par_vp_min)
            !         if(if_par_ip)  x(:,:,:,2) = (m%vp*m%rho -par_ip_min)  / (par_ip_max -par_ip_min)
            !     end select
            do ipar=1,npar
                if(pars(ipar)=='vp' ) x(:,:,:,ipar)= (m%vp      -pars_min(ipar))  / (pars_max(ipar) -pars_min(ipar))
                if(pars(ipar)=='rho') x(:,:,:,ipar)= (m%rho     -pars_min(ipar))  / (pars_max(ipar) -pars_min(ipar))
                if(pars(ipar)=='ip' ) x(:,:,:,ipar)= (m%vp*m%rho-pars_min(ipar))  / (pars_max(ipar) -pars_min(ipar))
            enddo
        else 
!             select case (parameterization)
! !                 case ('kappa-rho')
!                 case ('vp-rho')
!                     if(if_par_vp)  m%vp = x(:,:,:,1)*(par_vp_max -par_vp_min )+par_vp_min
!                     if(if_par_rho) m%rho= x(:,:,:,2)*(par_rho_max-par_rho_min)+par_rho_min
! !                case ('vprho-gardner')
!                 case ('vp-ip')
!                     if(if_par_vp)  m%vp = x(:,:,:,1)*(par_vp_max -par_vp_min )+par_vp_min
!                     if(if_par_ip)  m%rho=(x(:,:,:,2)*(par_ip_max -par_ip_min )+par_ip_min) /m%vp
!             end select
            do ipar=1,npar
                if(pars(ipar)=='vp' ) m%vp = x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar)) + pars_min(ipar)
                if(pars(ipar)=='rho') m%rho= x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar)) + pars_min(ipar)
                if(pars(ipar)=='ip' ) m%rho=(x(:,:,:,ipar)*(pars_max(ipar)-pars_min(ipar)) + pars_min(ipar)) /m%vp
            enddo
        endif
        
        
        !gradient
        if(present(g)) then
        if(dir=='m2x') then
    !         select case (parameterization)
    ! !             case ('kappa-rho')
    !             case ('vp-rho') !unit of g = [kg/m/s2]
    !                 if(if_par_vp)  g(:,:,:,1)=(gradient(:,:,:,1)*2.*m%vp*m%rho)               *(par_vp_max-par_vp_min)
    !                 if(if_par_rho) g(:,:,:,2)=(gradient(:,:,:,1)*m%vp*m%vp+gradient(:,:,:,2)) *(par_rho_max-par_rho_min)
    ! !             case ('vprho-gardner')
    !             case ('vp-ip')
    !                 if(if_par_vp) g(:,:,:,1)=(gradient(:,:,:,1)*m%vp*m%rho - gradient(:,:,:,2)*m%rho/m%vp) *(par_vp_max-par_vp_min)
    !                 if(if_par_ip) g(:,:,:,2)=(gradient(:,:,:,1)*m%vp       + gradient(:,:,:,2)/m%vp      ) *(par_ip_max-par_ip_min)
    !         end select
        endif
        endif
        
    end subroutine
    
end module
