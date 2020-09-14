module m_computebox
use m_sysio
use m_arrayop
use m_model, only: m
use m_shot

    private
    public cb, build_computebox
  
    type t_computebox
        integer :: nz,nx
        integer :: npml
        real    :: apml, bpml

        complex,dimension(:),allocatable :: damp_z, damp_x
        complex,dimension(:),allocatable :: damp_z_half, damp_x_half
    end type
    
    type(t_computebox) :: cb
    
    contains
    
    subroutine build_computebox(omega)

        cb%npml=get_setup_int('NBOUNDARYLAYER','NPML',default=10)
        cb%apml=get_setup_int('PML_COEFF','NPML',default=90)
        cb%bpml=0.

        cb%nz=m%nz+2*cb%npml
        cb%nx=m%nx+2*cb%npml

        if (mpiworld%is_master) then
            !pml
            call alloc(cb%damp_z,     [0,cb%nz+1])
            call alloc(cb%damp_x,     [0,cb%nx+1])
            call alloc(cb%damp_z_half,[0,cb%nz+1])
            call alloc(cb%damp_x_half,[0,cb%nx+1])
            
            ! !computebox
            ! call extend(m%vp, m%n2,m%n1, cb%vp, cb%npml,cb%npml)
            ! call extend(m%rho,m%n2,m%n1, cb%rho,cb%npml,cb%npml)

            !build pml
            call build_pml(cb%damp_z,cb%damp_z_half,cb%nz,omega)
            call build_pml(cb%damp_x,cb%damp_x_half,cb%nx,omega)

        endif

    end subroutine

    subroutine extend(m,nx,nz,me,nspx,nspz)
        real m(nz,nx),me(nz+2*nspz,nx+2*nspx)

        do ix=1,nx
        do iz=1,nz
            cb(iz+nspz,ix+nspx)     =m(iz,ix)
        end do
        end do

        !Extension of the model defined by m(nz,nx) into cb(nz+2*nsp,nz+2*nsp)
        !copying values at edges into PML zones
        
        !left & right
        do ix=1,nspx
        do iz=1,nz
            cb(iz+nspz,ix)          =m(iz,1)
            cb(iz+nspz,ix+nspx+nx)  =m(iz,nx)
        end do
        end do

        !top & bottom
        do ix=1,nx
        do iz=1,nspz
            cb(iz,ix+nspx)          =m(1,ix)
            cb(nz+iz+nspz,ix+nspx)  =m(nz,ix)
        end do
        end do

        !4 corners
        do ix=1,nspx
        do iz=1,nspz
            cb(iz,ix)                =m(1,1)
            cb(iz,ix+nspx+nx)        =m(1,nx)
            cb(iz+nz+nspz,ix)        =m(nz,1)
            cb(iz+nz+nspz,ix+nspx+nx)=m(nz,nx)
        end do
        end do

    end subroutine


    subroutine build_pml(damp,damp_half,n,omega)
        ! damp=1+i*gamma/omega
        ! where gamma ia a smoothly decreasing function
        ! damp      for classical grid
        ! damp_half for staggered grid

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
                betad     =0.
                betad_half=0.
            end if
            
            damp(i)     =1./(alpha     +c_i*eps     /(c_omega+c_i*betad     ))
            damp_half(i)=1./(alpha_half+c_i*eps_half/(c_omega+c_i*betad_half))
            
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
                betad     =0.
                betad_half=0.
            end if
            
            damp_half(n-i+1)=1./(alpha_half+c_i*eps_half/(c_omega+c_i*betad_half))

        end do

        damp_half(0)  =damp_half(1)
        damp_half(n+1)=damp_half(n)

    end subroutine
    
end
