module m_parameterization
use m_sysio
use m_arrayop
use m_model
use m_gradient, only: gradient

    public
    
    character(:),allocatable :: parameterization, active
    
    logical :: if_par_vp=.false., if_par_rho=.false., if_par_ip=.false., if_par_kpa=.false.
    
    real :: par_vp_min,  par_vp_max
    real :: par_rho_min, par_rho_max
    real :: par_ip_min,  par_ip_max
    real :: par_kpa_min, par_kpa_max
    
    integer :: npar

    contains
    
    subroutine init_parameterization
        
        parameterization=get_setup_char('PARAMETERIZATION',default='vp-rho')
        
        active=get_setup_char('ACTIVE_PARAMETER',default='vp1500-3400')
        
        !turn on inv pars and read in value ranges
        call read_active_par('vp',   if_par_vp, par_vp_min, par_vp_max)
        call read_active_par('rho',  if_par_rho,par_rho_min,par_rho_max)
        call read_active_par('ip',   if_par_ip, par_ip_min, par_ip_max)
        call read_active_par('kappa',if_par_kpa,par_kpa_min,par_kpa_max)
        
        !turn off inv pars depending on parameterization
        select case (parameterization)
            case ('kappa-rho')
                if_par_vp=.false.
                if_par_ip=.false.
            case ('vp-rho')
                if_par_kpa=.false.
                if_par_ip=.false.
            case ('vprho-gardner')
                if_par_kpa=.false.
                if_par_rho=.false.
                if_par_ip=.false.
            case ('vp-ip')
                if_par_kpa=.false.
                if_par_rho=.false.
        end select
        
        !reallocate model if initial model is not provided
        if(if_par_rho.or.if_par_kpa.or.if_par_ip) then
            call model_reallocate('rho')
        endif
        
        npar=count([if_par_kpa,if_par_rho,if_par_vp,if_par_ip].eqv..true.)
        if(mpiworld%is_master) write(*,*) 'Number of inversion parameters:', npar
        
    end subroutine
    
    subroutine read_active_par(par,if_par,amin,amax)
        character(*),intent(in) :: par
        logical,intent(out) :: if_par
        real,intent(out) :: amin,amax

        i=index(active,par)
        !print*,'i=',i

        if(i>0) then
            if_par=.true.
            j=index(active(i:),':')
            !print*,'j=',j
            
            read(active(i+len(par):i+j-2),*) amin
            
            k=index(active(i:),' ')
            if(k==0) k=len(active)+1-i
            !print*,'k=',k
            
            read(active(i+j:i+k-1),*) amax
            
            !print*,par_min, par_max
        else
            if_par=.false.
        endif
    end subroutine
    
    subroutine parameterization_transform(dir,x,g)
        character(3) :: dir
        real,dimension(m%nz,m%nx,m%ny,npar) :: x
        real,dimension(m%nz,m%nx,m%ny,npar),optional :: g
        
        !! Issue:
        !! what if I just want to invert for rho, ip ?
        !!
        
        !model
        if(dir=='m2x') then
            select case (parameterization)
                case ('vp-rho')
                    if(if_par_vp)  x(:,:,:,1) = (m%vp -par_vp_min)  / (par_vp_max -par_vp_min)
                    if(if_par_rho) x(:,:,:,2) = (m%rho-par_rho_min) / (par_rho_max-par_rho_min)
                case ('vp-ip')
                    if(if_par_vp)  x(:,:,:,1) = (m%vp -par_vp_min)  / (par_vp_max -par_vp_min)
                    if(if_par_ip)  x(:,:,:,2) = (m%vp*m%rho -par_ip_min)  / (par_ip_max -par_ip_min)
                end select
        else 
            select case (parameterization)
!                 case ('kappa-rho')
                case ('vp-rho')
                    if(if_par_vp)  m%vp = x(:,:,:,1)*(par_vp_max -par_vp_min )+par_vp_min
                    if(if_par_rho) m%rho= x(:,:,:,2)*(par_rho_max-par_rho_min)+par_rho_min
!                case ('vprho-gardner')
                case ('vp-ip')
                    if(if_par_vp)  m%vp = x(:,:,:,1)*(par_vp_max -par_vp_min )+par_vp_min
                    if(if_par_ip)  m%rho=(x(:,:,:,2)*(par_ip_max -par_ip_min )+par_ip_min) /m%vp
            end select
        endif
        
        
        !gradient
        if(present(g)) then
        if(dir=='m2x') then
            select case (parameterization)
    !             case ('kappa-rho')
                case ('vp-rho') !unit of g = [kg/m/s2]
                    if(if_par_vp)  g(:,:,:,1)=(gradient(:,:,:,1)*2.*m%vp*m%rho)               *(par_vp_max-par_vp_min)
                    if(if_par_rho) g(:,:,:,2)=(gradient(:,:,:,1)*m%vp*m%vp+gradient(:,:,:,2)) *(par_rho_max-par_rho_min)
    !             case ('vprho-gardner')
                case ('vp-ip')
                    if(if_par_vp) g(:,:,:,1)=(gradient(:,:,:,1)*m%vp*m%rho - gradient(:,:,:,2)*m%rho/m%vp) *(par_vp_max-par_vp_min)
                    if(if_par_ip) g(:,:,:,2)=(gradient(:,:,:,1)*m%vp       + gradient(:,:,:,2)/m%vp      ) *(par_ip_max-par_ip_min)
            end select
        endif
        endif
        
    end subroutine
    
end module
