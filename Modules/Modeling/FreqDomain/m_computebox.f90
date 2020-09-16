module m_computebox
use m_global
use m_sysio
use m_arrayop
use m_model, only: m

    private
    public cb, build_computebox, computebox_dealloc_model
  
    type t_computebox
        integer :: nz,nx, n
        integer :: npml
        real    :: apml, bpml

        real,dimension(:,:),allocatable :: vp, rho

        complex,dimension(:),allocatable :: damp_z,      damp_x
        complex,dimension(:),allocatable :: damp_z_half, damp_x_half
    end type
    
    type(t_computebox) :: cb
    
    contains
    
    subroutine build_computebox(freq)
        complex :: freq, omega

        cb%npml=get_setup_int('NBOUNDARYLAYER','NPML',default=10)
        cb%apml=get_setup_int('PML_COEFF','APML',default=90)
        cb%bpml=0.

        cb%nz=m%nz+2*cb%npml
        cb%nx=m%nx+2*cb%npml

        cb%n=cb%nz*cb%nx

        if (mpiworld%is_master) then
            !pml
            call alloc(cb%damp_z,     [0,cb%nz+1])
            call alloc(cb%damp_x,     [0,cb%nx+1])
            call alloc(cb%damp_z_half,[0,cb%nz+1])
            call alloc(cb%damp_x_half,[0,cb%nx+1])
            
            !models in computebox
            call alloc(cb%vp, cb%nz,cb%nx); call m2cb(m%vp, cb%vp )
            call alloc(cb%rho,cb%nz,cb%nx); call m2cb(m%rho,cb%rho)

            !build pml
            omega=freq*2.*r_pi
            call build_pml(cb%damp_z,cb%damp_z_half,cb%nz,omega)
            call build_pml(cb%damp_x,cb%damp_x_half,cb%nx,omega)

        endif

    end subroutine

    subroutine computebox_dealloc_model
        deallocate(cb%vp)
        deallocate(cb%rho)
    end subroutine

    subroutine m2cb(in,out)
        real :: in ( m%nz, m%nx)
        real :: out(cb%nz,cb%nx)

        integer,dimension(:),allocatable :: indz, indx, indb

        call alloc(indz,m%nz,initialize=.false.);    indz=[(i,i=1,m%nz)] !construct index array without semicolumns..
        call alloc(indx,m%nx,initialize=.false.);    indx=[(i,i=1,m%nx)]
        call alloc(indb,cb%npml,initialize=.false.); indb=[(i,i=1,cb%npml)]

        out(indz+cb%npml,indx+cb%npml) = in

        !copy values for PML zones
        !left & right zones
        do ix=1,cb%npml
            out(indz+cb%npml ,ix)               =in(indz+cb%npml ,1)
            out(indz+cb%npml ,ix+cb%npml+m%nx)  =in(indz+cb%npml ,m%nx)
        end do

        !top & bottom zones
        do iz=1,cb%npml
            out(iz,             indx+cb%npml) =in(1,   indx)
            out(iz+m%nz+cb%npml,indx+cb%npml) =in(m%nz,indx)
        end do

        !4 corners
        out(indb,indb)                          =in(1,1)
        out(indb,indb+cb%npml+m%nx)             =in(1,m%nx)
        out(indb+cb%npml+m%nz,indb)             =in(m%nz,1)
        out(indb+cb%npml+m%nz,indb+cb%npml+m%nx)=in(m%nz,m%nx)

    end subroutine

    subroutine build_pml(damp,damp_half,n,omega)
        complex,dimension(0:n+1) :: damp, damp_half
        complex :: omega
        ! damp=1+i*gamma/omega
        ! where gamma ia a smoothly decreasing function
        ! damp      for classical grid
        ! damp_half for staggered grid

        complex :: beta_w, betad, betad_half
        ! alphad: add a constant to the x function
        ! betad: add an imaginary part to frequency

        real,parameter :: alphad=1., beta_max=0.        
        real,parameter :: pi2=r_pi/2
        
        !nu=c_omega+cmplx(omega,betad*omega)

        ! values of damp and dampn outside PML (equal to 1)    
        damp     =cmplx(1.,0.)
        damp_half=cmplx(1.,0.)
    
        beta_w=beta_max*omega
        
        xpml=real(cb%npml)*m%h
        xmax=real(n-1)*m%h


        do i=1,cb%npml

            x     =real(i-1)*m%h
            x_half=real(i-1)*m%h+0.5*m%h

            betad     =x     *beta_w/xpml
            betad_half=x_half*beta_w/xpml

            !compute eps and alpha
            !!with polynomial function
            !eps      =cb%apml*exp(cb%bpml*log(((xpml-x     )/xpml)+1e-10))
            !eps_half =cb%apml*exp(cb%bpml*log(((xpml-x_half)/xpml)+1e-10))
            !alpha     =(alphad-1.)*exp(cb%bpml*log(((xpml-x     )/xpml)+1e-10))+1.
            !alpha_half=(alphad-1.)*exp(cb%bpml*log(((xpml-x_half)/xpml)+1e-10))+1.

            !with cosine function
            tmp     =1.-cos((xpml-x     )*pi2/xpml)
            tmp_half=1.-cos((xpml-x_half)*pi2/xpml)
            eps     =      cb%apml*tmp
            eps_half      =cb%apml*tmp_half
            alpha     =(alphad-1.)*tmp     +1.
            alpha_half=(alphad-1.)*tmp_half+1.

            if (abs(beta_max)<1e-3) then
                betad     =cmplx(0.,0.)
                betad_half=cmplx(0.,0.)
            end if
            
            damp(i)     =1./(alpha     +c_i*eps     /(omega+c_i*betad     ))
            damp_half(i)=1./(alpha_half+c_i*eps_half/(omega+c_i*betad_half))
            
            damp(n-i+1)=damp(i)

        end do

        damp(0)  =damp(1)
        damp(n+1)=damp(n)


        do i=1,cb%npml+1

            x_half = xmax+0.5*m%h-real(i-1)*m%h
            y_half=       real(i-1)*m%h-0.5*m%h

            betad_half=y_half*beta_w/xpml

            !compute eps and alpha
            !!with polynomial function
            !eps_half  =    cb%apml*exp(cb%bpml*log(((x_half-(xmax-xpml))/xpml)+1e-10))
            !alpha_half=(alphad-1.)*exp(cb%bpml*log(((x_half-(xmax-xpml))/xpml)+1e-10))+1.

            !with cosine function
            tmp_half=1.-cos((x_half-(xmax-xpml))*pi2/xpml)
            eps_half=      cb%apml*tmp_half
            alpha_half=(alphad-1.)*tmp_half+1.
        
            if (abs(beta_max)<1e-3) then
                betad     =cmplx(0.,0.)
                betad_half=cmplx(0.,0.)
            end if
            
            damp_half(n-i+1)=1./(alpha_half+c_i*eps_half/(omega+c_i*betad_half))

        end do

        damp_half(0)  =damp_half(1)
        damp_half(n+1)=damp_half(n)

    end subroutine
    
end
