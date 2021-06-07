module m_field
use m_mpienv
use m_arrayop
use m_sysio
use m_model, only:m
use m_shot, only:shot
use m_computebox, only: cb
use m_FDSG

use, intrinsic :: ieee_arithmetic
    
    !fields
    type t_field

        !shorthand for greek letters
        !alp bta gma  del(dta) eps zta 
        !eta tht iota kpa lda mu
        !nu xi omi pi rho sgm
        !tau ups phi chi psi oga
        !
        !buo : buoyancy
        !vx,vy : horizontal velocities
        !leading _ : inverse
        !
        !physical components
        real,dimension(:,:,:),allocatable :: vx,vy,vz,p
        !local models in computebox
        real,dimension(:,:,:),allocatable :: buox, buoy, buoz, kpa, _kpa
        !cpml components
        real,dimension(:,:,:),allocatable :: cpml_dvx_dx,cpml_dvy_dy,cpml_dvz_dz
        real,dimension(:,:,:),allocatable :: cpml_dp_dz, cpml_dp_dx, cpml_dp_dy

        integer :: it      !index of time instant
        integer :: ift,ilt !first & last time instant
        integer :: nt
        real dt
        
        real,dimension(:),allocatable :: wavelet
            
        !snapshot
        logical :: if_snapshot
        integer :: it_delta_snapshot

        contains
        procedure :: print_info => print_info
        procedure :: estim_RAM => estim_RAM
        procedure :: check_model => check_model
        procedure :: check_discretization => check_discretization
        procedure :: init  => init
        procedure :: check_value => check_value
        procedure :: reinit_boundary => reinit_boundary
        procedure :: write => write
        procedure :: add_source => add_source

        procedure :: propagate => propagate
        procedure :: propagate_reverse => propagate_reverse

    end type
    

    contains
    
    !========= public procedures =================

    subroutine estim_RAM
    end subroutine
    
    subroutine check_model  !not required
        
    end
    
    subroutine check_discretization
        !grid dispersion condition
        if (5.*m%cell_diagonal > cb%velmin/shot%src%fpeak/2.) then  !FDTDo4 rule
            call warn('WARNING: Shot# '//shot%cindex//' can have grid dispersion!')
            if(mpiworld%is_master) write(*,*) 'Shot# '//shot%cindex//' 5*dx, velmin, fpeak:',5.*m%cell_diagonal, cb%velmin,shot%src%fpeak
        endif
        
        !CFL condition
        cfl = cb%velmax*shot%src%dt*m%cell_inv_diagonal*sum(abs(fdcoeff_o4))! ~0.494 (3D); ~0.606 (2D)
        if(mpiworld%is_master) write(*,*) 'CFL value:',CFL
        
        if(cfl>1.) then
            if(mpiworld%is_master) write(*,*) 'Shot# '//shot%cindex//' velmax, dt, 1/dx:',cb%velmax,shot%src%dt,m%cell_inv_diagonal
            call error('CFL > 1 on shot# '//shot%cindex//'!')
            stop
        endif
        
    end subroutine
    
    subroutine init
                       
        !propagation time frames
        ift=1
        ilt=shot%src%nt
        nt=ilt-ift+1
        dt=shot%src%dt       

        call alloc(wavelet,nt)
        wavelet = shot%src%wavelet
        
        !snapshot
        if_snapshot=setup_get_logical('IF_SNAPSHOT',default=.false.)
        if_snapshot=if_snapshot.and.mpiworld%is_master
        if(if_snapshot) it_delta_snapshot=setup_get_int('SNAPSHOT_DELTA_IT',default=50)
        
    end subroutine
    
    subroutine check_value(f,name)
        class(t_field), intent(in) :: f
        character(*) :: name
        
        if(mpiworld%is_master) write(*,*) name//' sample values:',minval(f%p),maxval(f%p)
        
        if(any(.not. ieee_is_finite(f%vz))) then
            error(name//' values become Infinity on Shot# '//shot%cindex//' !!')
        endif
        if(any(ieee_is_nan(f%vz))) then
            error(name//' values become NaN on Shot# '//shot%cindex//' !!')
        endif
        
    end subroutine
        
    subroutine write(f,iunit)
        class(t_field), intent(in) :: f
        write(iunit) f%vz
    end subroutine
 
    subroutine add_RHS(wavelet)
        self%wavelet=>wavelet
    end subroutine
    
    subroutine acquire(f,seismo)
        class(t_field), intent(in) :: f
        real,dimension(*) :: seismo
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            ify=shot%rcv(ircv)%ify; iy=shot%rcv(ircv)%iy; ily=shot%rcv(ircv)%ily
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (1)
                    seismo(ircv)=sum(  f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                    case (2)
                    seismo(ircv)=sum( f%vx(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                    case (3)
                    seismo(ircv)=sum( f%vy(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                    case (4)
                    seismo(ircv)=sum( f%vz(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (1) !p[iz,ix,iy]
                    seismo(ircv)= f%p(iz,ix,iy)
                    case (2) !vx[iz,ix-0.5,iy]
                    seismo(ircv)=f%vx(iz,ix,iy)
                    case (3) !vy[iz,ix,iy-0.5]
                    seismo(ircv)=f%vy(iz,ix,iy)
                    case (4) !vz[iz-0.5,ix,iy]
                    seismo(ircv)=f%vz(iz,ix,iy)
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    subroutine acquire_reverse(f,w)
        class(t_field), intent(in) :: f
        real w
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        ify=shot%src%ify; iy=shot%src%iy; ily=shot%src%ily
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (1)
                w = sum(lm%invkpa(ifz:ilz,ifx:ilx,ify:ily)*f%p(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
                case (2)
                w = sum(     f%vx(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
                case (3)
                w = sum(     f%vy(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
                case (4)
                w = sum(     f%vz(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
            end select
            
        else
            select case (shot%src%icomp)
                case (1) !p[iz,ix,iy]
                w=lm%invkpa(iz,ix,iy)*f%p(iz,ix,iy)
                
                case (2) !vx[iz,ix-0.5,iy]
                w=f%vx(iz,ix,iy)
                
                case (3) !vy[iz,ix,iy-0.5]
                w=f%vy(iz,ix,iy)
                
                case (4) !vz[iz-0.5,ix,iy]
                w=f%vz(iz,ix,iy)
            end select
            
        endif
        
    end subroutine
    
end
