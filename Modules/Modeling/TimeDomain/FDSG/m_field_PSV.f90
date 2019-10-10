module m_field
use m_mpienv
use m_arrayop
use m_sysio
use m_model, only:m
use m_shot, only:shot
use m_computebox, only: cb

    private fdcoeff_o4,fdcoeff_o8,c1x,c1z,c2x,c2z,c3x,c3z,c4x,c4z
    private k_x,k_z,npower,rcoef
    public
    
    !FD coeff
    real,dimension(2),parameter :: fdcoeff_o4 = [1.125,      -1./24.]
    real,dimension(4),parameter :: fdcoeff_o8 = [1225./1024, -245./3072., 49./5120., -5./7168]
    
    real :: c1x, c1z
    real :: c2x, c2z
    real :: c3x, c3z
    real :: c4x, c4z
    
    !CPML
    real,parameter :: k_x = 1.
    real,parameter :: k_z = 1.
    real,parameter :: npower=2.
    real,parameter :: rcoef=0.001
    
    !local models in computebox
    type t_localmodel
        sequence
        real,dimension(:,:),allocatable :: buox, buoz, ldap2mu, lda, mu
        real,dimension(:,:),allocatable :: inv_ldapmu_4mu
    end type
    
    type(t_localmodel) :: lm
    
    !fields
    type t_field
        sequence
        real,dimension(:,:),allocatable :: vx,vy,vz,sxx,szz,sxz  !vy is just for fun
        real,dimension(:,:),allocatable :: cpml_dvx_dx,cpml_dvz_dz,cpml_dvx_dz,cpml_dvz_dx !for cpml
        real,dimension(:,:),allocatable :: cpml_dsxx_dx,cpml_dszz_dz,cpml_dsxz_dx,cpml_dsxz_dz
    end type

    !sampling factor
    real,parameter :: factor_xx=0.6666667 , factor_zz=0.3333333
    
    !hicks interpolation
    logical :: if_hicks
    
    !gradient info
    character(*),parameter :: gradient_info='lda-mu-rho'
    
    contains
    
    !========= use before propagation =================
    subroutine field_print_info
        !modeling method
        call hud('WaveEq : Time-domain isotropic P-SV (2D elastic) system')
        call hud('1st-order Velocity-Stress formulation')
        call hud('Staggered-grid Finite-difference method')
        call hud('Cartesian O(x4,t2) stencil')
        
        !stencil constant
        if(mpiworld%is_master) then
            write(*,*) 'Coeff:',fdcoeff_o4
        endif
        c1x=fdcoeff_o4(1)/m%dx; c1z=fdcoeff_o4(1)/m%dz
        c2x=fdcoeff_o4(2)/m%dx; c2z=fdcoeff_o4(2)/m%dz
        
        !hicks point interpolation
        if_hicks=get_setup_logical('IF_HICKS',default=.true.)
        
    end subroutine
    
    subroutine check_model  !not required
        
    end
    
    subroutine check_discretization
        !grid dispersion condition
        if (5.*m%cell_diagonal > cb%velmin/shot%src%fpeak/2.) then  !FDTDo4 rule
            write(*,*) 'WARNING: Shot# '//shot%cindex//' can have grid dispersion!'
            write(*,*) 'Shot# '//shot%cindex//' 5*dx, velmin, fpeak:',5.*m%cell_diagonal, cb%velmin,shot%src%fpeak
        endif
        
        !CFL condition
        cfl = cb%velmax*shot%src%dt*m%cell_inv_diagonal*sum(abs(fdcoeff_o4))! ~0.606
        if(mpiworld%is_master) write(*,*) 'CFL value:',CFL
        
        if(cfl>1.) then
            write(*,*) 'ERROR: CFL > 1 on shot# '//shot%cindex//'!'
            write(*,*) 'Shot# '//shot%cindex//' velmax, dt, dx:',cb%velmax,shot%src%dt,m%cell_inv_diagonal
            stop
        endif
        
    end subroutine
    
    !========= use inside propagation =================
    
    subroutine init_field_localmodel
        real,dimension(:,:),allocatable :: temp_mu
        
        ! (Lamé par)  (Rock phy)  (Seismic wave)     (Voigt)
        ! lda+2mu   = K+4/3G    = rho*Vp^2         = c11=c33
        ! lda       = K-2/3G    = rho*(Vp^2-2Vs^2) = c13
        ! mu        = G         = rho*Vs^2         = c44
        !
        ! Poisson ratio nv = (3K-2G)/(6k+2G) = 0.5-3/(6sqrt(Vp/Vs)+2)

        call alloc(lm%buox,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(lm%buoz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(lm%ldap2mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(lm%lda, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        call alloc(lm%mu,  [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)

        call alloc(temp_mu,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],initialize=.false.)
        
        !Staggered grid:
        ! <---*---<  
        ! |   |   |  * integer point   - lda+2mu,lda 
        ! x---v---x  < half in x dir   - buox  
        ! |   |   |  v half in z dir   - buoz
        ! <---*---<  x half in x&z dir - mu
        !     (grid index)      (real index)
        ! lda+2mu,lda(iz,ix) := lda+2mu,lda[iz,ix]
        !        buox(iz,ix) := buox[iz,ix-0.5]
        !        buoz(iz,ix) := buoz[iz-0.5,ix]
        !          mu(iz,ix) :=   mu[iz-0.5,ix-0.5]

        lm%ldap2mu(:,:)=cb%rho(:,:,1)*cb%vp(:,:,1)**2
           temp_mu(:,:)=cb%rho(:,:,1)*cb%vs(:,:,1)**2

        lm%lda=lm%ldap2mu-2.*temp_mu

        temp_mu=1./temp_mu

        do ix=cb%ifx+1,cb%ilx
        do iz=cb%ifz+1,cb%ilz
                lm%mu(iz,ix)=4./( temp_mu(iz-1,ix-1) &
                                 +temp_mu(iz-1,ix  ) &
                                 +temp_mu(iz  ,ix-1) &
                                 +temp_mu(iz  ,ix  ))
        end do
        end do

        where( isnan(lm%mu) .or. lm%mu==lm%mu+1. )
            lm%mu=0.
        endwhere

        lm%mu(cb%ifz,:)=lm%mu(cb%ifz+1,:)
        lm%mu(:,cb%ifx)=lm%mu(:,cb%ifx+1)

        !check mu values
        if(mpiworld%is_master) then
            write(*,*) 'lm%mu sanity:', minval(lm%mu),maxval(lm%mu), any(isnan(lm%mu)), any(lm%mu==lm%mu+1.)
        endif
        
! if(mpiworld%is_master) then
! write(*,*) 'lm%ldap2mu sanity:', minval(lm%ldap2mu),maxval(lm%ldap2mu)
! write(*,*) 'lm%lda     sanity:', minval(lm%lda),maxval(lm%lda)
! endif

        do iz=cb%ifz+1,cb%ilz
            lm%buoz(iz,:)=0.5/cb%rho(iz,:,1)+0.5/cb%rho(iz-1,:,1)
        enddo
        
        do ix=cb%ifx+1,cb%ilx
            lm%buox(:,ix)=0.5/cb%rho(:,ix,1)+0.5/cb%rho(:,ix-1,1)
        enddo

        lm%buoz(cb%ifz,:)=1./cb%rho(cb%ifz,:,1)
        lm%buox(:,cb%ifx)=1./cb%rho(:,cb%ifx,1)

        deallocate(temp_mu)

    end subroutine
    
    subroutine init_field(f)
        type(t_field) :: f
        
        call alloc(f%vx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%vz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%sxx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%szz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%sxz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        
        call alloc(f%cpml_dvx_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dvz_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dvx_dz, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dvz_dx, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dsxx_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dszz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dsxz_dx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        call alloc(f%cpml_dsxz_dz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx])
        
    end subroutine
    
    subroutine check_field(f,name)
        type(t_field) :: f
        character(*) :: name
        
        if(mpiworld%is_master) write(*,*) name//' sample values:',minval(f%sxz),maxval(f%sxz)
        
        !if(any(f%p-1.==f%p)) then !this is not a good numerical judgement..
        !    write(*,*) 'ERROR: '//name//' values become +-Infinity on Shot# '//shot%cindex//' !!'
        !    stop
        !endif
        if(any(isnan(f%vz))) then
           write(*,*) 'ERROR: '//name//' values become NaN on Shot# '//shot%cindex//' !!'
           stop
        endif
        
    end subroutine
    
    subroutine field_cpml_reinitialize(f)
        type(t_field) :: f
        f%cpml_dvx_dx=0.
        f%cpml_dvz_dz=0.
        f%cpml_dvx_dz=0.
        f%cpml_dvz_dx=0.
        f%cpml_dsxx_dx=0.
        f%cpml_dszz_dz=0.
        f%cpml_dsxz_dx=0.
        f%cpml_dsxz_dz=0.
    end subroutine
    
    subroutine write_field(iunit,f)
        type(t_field) :: f
        write(iunit) f%vz
    end subroutine
    
    !========= forward propagation =================
    !WE: du_dt = MDu
    !u=[vx vy vz sxx szz sxz]^T, p=(2sxx+szz)/3
    !  [b                     ]    [      dx 0  dz]
    !  |  b                   |    |      0  dz dx|
    !M=|   lda+2mu   lda      |, D=|dx 0          |
    !  |      lda   lda+2mu   |    |0  dz         |
    !  [                    mu]    [dz dx         ]
    !  (b=1/rho)
    !
    !Discretization (staggered grid in space and time):
    !  (grid index)     (real index)
    ! sxx,szz(iz,ix) := sxx,szz[iz,ix]^it+0.5,it+1.5,...
    !      vx(iz,ix) := vx[iz,ix-0.5]^it,it+1,...
    !      vz(iz,ix) := vy[iz-0.5,ix]^it,it+1,...
    !     sxz(iz,ix) := sxz[iz-0.5,ix-0.5]^it,it+1,...
    ! in space:
    ! <---*---<  
    ! |   |   |  * integer point   - sxx,szz - lda+2mu,lda 
    ! x---v---x  < half in x dir   - vx      - buox  
    ! |   |   |  v half in z dir   - vz      - buoz
    ! <---*---<  x half in x&z dir - sxz     - mu
    
    !add RHS to v^it
    subroutine put_velocities(time_dir,it,w,f)
        integer :: time_dir,it
        real :: w
        type(t_field) :: f
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        
        source_term=time_dir*w
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (2)
                f%vx(ifz:ilz,ifx:ilx) = f%vx(ifz:ilz,ifx:ilx) + source_term *lm%buox(ifz:ilz,ifx:ilx) *shot%src%interp_coeff(:,:,1)
                
                case (4)
                f%vz(ifz:ilz,ifx:ilx) = f%vz(ifz:ilz,ifx:ilx) + source_term *lm%buoz(ifz:ilz,ifx:ilx) *shot%src%interp_coeff(:,:,1)
                
            end select
            
        else
            select case (shot%src%icomp)
                case (2) !horizontal x force on vx[iz,ix-0.5,iy]
                f%vx(iz,ix) = f%vx(iz,ix) + source_term*lm%buox(iz,ix)
                
                case (4) !vertical force     on vz[iz-0.5,ix,iy]
                f%vz(iz,ix) = f%vz(iz,ix) + source_term*lm%buoz(iz,ix)
                
            end select
            
        endif
        
    end subroutine
    
    !v^it -> v^it+1 by FD of s^it+0.5
    subroutine update_velocities(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl

        ifz=bl(1)+2
        ilz=bl(2)-2
        ifx=bl(3)+2
        ilx=bl(4)-2
        !ify=bl(5)+2
        !ily=bl(6)-2
                
        if(m%if_freesurface) ifz=max(ifz,1)
        
        call fd2d_flat_velocities(f%vx,f%vz,f%sxx,f%szz,f%sxz,                                &
                                  f%cpml_dsxx_dx,f%cpml_dszz_dz,f%cpml_dsxz_dx,f%cpml_dsxz_dz,&
                                  lm%buox,lm%buoz,                                            &
                                  ifz,ilz,ifx,ilx,time_dir*shot%src%dt)
        
        dz_dx = m%dz/m%dx

        !apply free surface boundary condition if needed
        !free surface is located at [1,ix,iy] level
        !dszz=0 -> -(lda+2mu)dvz_dz = lda dvx_dx
        !       -> -ldap2mu[1,ix]*(vz[1.5,ix]-vz[0.5,ix])/dz = lda[1,ix]*(vx[1,ix-0.5]-vx[1,ix+0.5])/dx
        !       -> -ldap2mu(1,ix)*(vz(2  ,ix)-vz(1  ,ix))/dz = lda(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
        !       ->  ldap2mu(1,ix)*(vz(1  ,ix)-vz(2  ,ix))/dz = lda(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
        !dsxz=0 -> -dvx_dz = dvz_dx
        !dsxz=0 -> -(vx[2,ix+0.5]-vx[0,ix+0.5])/2dz = ( (vz[0.5,ix]+vz[1.5,ix])/2 - (vz[0.5,ix+1]+vz[1.5,ix+1])/2 )/dx
        !       -> -(vx(2,ix+1  )-vx(0,ix+1  ))/2dz = ( (vz(1  ,ix)+vz(2  ,ix))/2 - (vz(1  ,ix+1)+vz(2  ,ix+1))/2 )/dx
        !       ->  (vx(0,ix+1  )-vx(2,ix+1  ))/ dz = ( (vz(1  ,ix)+vz(2  ,ix))   -  vz(1  ,ix+1)-vz(2  ,ix+1)    )/dx
        if(m%if_freesurface) then
            do ix=ifx,ilx
                f%vz(1,ix)= f%vz(2,ix) + lm%lda(1,ix)*(f%vx(1,ix)-f%vx(1,ix+1))*dz_dx/lm%ldap2mu(1,ix)
                !f%vx(0,ix)= f%vx(2,ix) + (f%vz(1,ix)+f%vz(2,ix)-f%vz(1,ix+1)-f%vz(2,ix+1))*dz_dx/lm%ldap2mu(1,ix) !toy2del says this condition is not needed
            enddo
        endif
        
    end subroutine
    
    !add RHS to s^it+0.5
    subroutine put_stresses(time_dir,it,w,f)
        integer :: time_dir
        real :: w
        type(t_field) :: f
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        
        source_term=time_dir*w
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (1)
                f%sxx(ifz:ilz,ifx:ilx) = f%sxx(ifz:ilz,ifx:ilx) + source_term *shot%src%interp_coeff(:,:,1)
                f%szz(ifz:ilz,ifx:ilx) = f%szz(ifz:ilz,ifx:ilx) + source_term *shot%src%interp_coeff(:,:,1)
            end select
            
        else
            select case (shot%src%icomp)
                case (1) !explosion on s[iz,ix,iy]
                f%sxx(iz,ix) = f%sxx(iz,ix) + source_term
                f%szz(iz,ix) = f%szz(iz,ix) + source_term
            end select
            
        endif
        
    end subroutine
    
    !s^it+0.5 -> s^it+1.5 by FD of v^it+1
    subroutine update_stresses(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl
        
        ifz=bl(1)+2
        ilz=bl(2)-2
        ifx=bl(3)+2
        ilx=bl(4)-2
        !ify=bl(5)+2
        !ily=bl(6)-2
        
        if(m%if_freesurface) ifz=max(ifz,1)
        
        call fd2d_flat_stresses(f%vx,f%vz,f%sxx,f%szz,f%sxz,                            &
                                f%cpml_dvx_dx,f%cpml_dvz_dz,f%cpml_dvx_dz,f%cpml_dvz_dx,&
                                lm%ldap2mu,lm%lda,lm%mu,                                &
                                ifz,ilz,ifx,ilx,time_dir*shot%src%dt)

        !apply free surface boundary condition if needed
        !free surface is located at [1,ix] level
        !so explicit boundary condition: szz(1,ix)=0
        !and antisymmetric mirroring: sxz[0.5,ix-0.5]=-sxz[1.5,ix-0.5] -> sxz(1,ix)=-sxz(2,ix)
        if(m%if_freesurface) then
            f%szz(1,:)=0.
            f%sxz(1,:)=-f%sxz(2,:)
        endif
        
    end subroutine
    
    !get v^it+1 or s^it+1.5
    subroutine get_field(f,seismo)
        type(t_field) :: f
        real,dimension(*) :: seismo
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (1)
                    seismo(ircv)=sum( (factor_xx*f%sxx(ifz:ilz,ifx:ilx) + factor_zz*f%szz(ifz:ilz,ifx:ilx)) *shot%rcv(ircv)%interp_coeff(:,:,1) )
                    case (2)
                    seismo(ircv)=sum( f%vx(ifz:ilz,ifx:ilx) *shot%rcv(ircv)%interp_coeff(:,:,1) )
                    case (4)
                    seismo(ircv)=sum( f%vz(ifz:ilz,ifx:ilx) *shot%rcv(ircv)%interp_coeff(:,:,1) )
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (1) !sxx,szz[iz,ix,iy]
                    seismo(ircv) = factor_xx * f%sxx(iz,ix) + factor_zz * f%szz(iz,ix)
                    case (2) !vx[iz,ix-0.5,iy]
                    seismo(ircv)=f%vx(iz,ix)
                    case (4) !vz[iz-0.5,ix,iy]
                    seismo(ircv)=f%vz(iz,ix)
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    
    !========= adjoint propagation =================
    !WE:        du_dt   = MDu   + src
    !->        Ndu_dt   = Du    + Nsrc, N=M^-1
    !Adjoint:  Ndv_dt^T = D^Tv  + adjsrc
    !->        -dv_dt   =-MDv   + Mr
    !     (reverse time) (negative derivatives)
    !
    !Discrete form (staggered grid in space and time, 2D as example):
    !    [ vx^it+1  ] [ vx^it    ]   [            dx(-)      dz(-)][ vx^it+1  ]
    ! N  | vz^it+1  | | vz^it    |   |                 dz(-) dx(-)|| vz^it+1  |
    !----|sxx^it+1.5|-|sxx^it+0.5| = |dx(+)                       ||sxx^it+0.5|  +Nsrc
    ! dt |szz^it+1.5| |szz^it+0.5|   |      dz(+)                 ||szz^it+0.5|
    !    [sxz^it+1.5] [sxz^it+0.5]   [dz(+) dx(+)                 ][sxz^it+0.5]
    !dx(-):=a1(s(ix)-s(ixm1))+a2(s(ixp1)-s(ixm2))  (O(x4))
    !dx(+):=a1(v(ixp1)-v(ix))+a2(v(ixp2)-v(ixm1))  (O(x4))
    !Adjoint:
    !    [ vx^it    ] [ vx^it+1  ]   [              dx(+)^T        dz(+)^T][ vx^it+1  ]
    ! N  | vz^it    |-| vz^it+1  |   |                     dz(+)^T dx(+)^T|| vz^it+1  |
    !----|sxx^it+0.5| [sxx^it+1.5] = |dx(-)^T                             ||sxx^it+0.5]  +Nsrc
    ! dt |szz^it+0.5| [szz^it+1.5]   |        dz(-)^T                     ||szz^it+0.5]
    !    [sxz^it+0.5] [sxz^it+1.5]   [dz(-)^T dx(-)^T                     ][sxz^it+0.5]
    !dx(+)^T=a1(s(ixm1)-s(ix))+a2(s(ixm2)-s(ixp1))=-dx(-)
    !dx(-)^T=a1(v(ix)-v(ixp1))+a2(v(ixm1)-v(ixp2))=-dx(+)
    !-dx,-dz can be regarded as -dt, thus allowing to use same code with flip of dt sign.
    !transposes of time derivatives centered on v^it+0.5 and s^it+1
    !so v^it+1   - v^it     => v^it - v^it+1
    !   s^it+1.5 - s^it+0.5 => s^it+0.5 - s^it+1.5
    !
    !In each time step:
    !d=RGAsrc, src:source term, A:put source, G:propagator, R:sampling at receivers
    !d=[vx vz p]^T, R=[I  0   0 ], A=[diag2(b) 0], src=[fx fz p]^T
    !                 |0 2/3 1/3|    [   0     1]
    !                 [   0    1]
    !Adjoint: src=A^T N G^T M R^T d
    !M R^T:put adjoint sources, 
    !G^T:  adjoint propagator,
    !A^T N:get adjoint fields
    
    
    !add RHS to s^it+1.5
    subroutine put_stresses_adjoint(time_dir,it,adjsource,f)
        integer :: time_dir
        real,dimension(*) :: adjsource
        type(t_field) :: f
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            
            tmp=adjsource(ircv)
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (1) !pressure adjsource
                    f%sxx(ifz:ilz,ifx:ilx) = f%sxx(ifz:ilz,ifx:ilx)  &
                        +tmp*( factor_xx*lm%ldap2mu(ifz:ilz,ifx:ilx) &
                              +factor_zz*lm%lda(ifz:ilz,ifx:ilx) )   &
                        *shot%rcv(ircv)%interp_coeff(:,:,1)
        
                    f%szz(ifz:ilz,ifx:ilx) = f%szz(ifz:ilz,ifx:ilx)    &
                        +tmp*( factor_xx*lm%lda(ifz:ilz,ifx:ilx)       &
                              +factor_zz*lm%ldap2mu(ifz:ilz,ifx:ilx) ) &
                        *shot%rcv(ircv)%interp_coeff(:,:,1)
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (1) !pressure adjsource
                    !sxx[iz,ix]
                    f%sxx(iz,ix) = f%sxx(iz,ix)            &
                        +tmp*( factor_xx*lm%ldap2mu(iz,ix) &
                              +factor_zz*lm%lda(iz,ix) )
                    
                    !szz[iz,ix]
                    f%szz(iz,ix) = f%szz(iz,ix)            &
                        +tmp*( factor_xx*lm%lda(iz,ix)     &
                              +factor_zz*lm%ldap2mu(iz,ix) )
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    !s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
    subroutine update_stresses_adjoint(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl
        
        !same code with -dt
        call update_stresses(time_dir,it,f,bl)
        
    end subroutine
    
    !add RHS to v^it+1
    subroutine put_velocities_adjoint(time_dir,it,adjsource,f)
        integer :: time_dir
        real,dimension(*) :: adjsource
        type(t_field) :: f
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            
            tmp=adjsource(ircv)
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (2) !horizontal x adjsource
                    f%vx(ifz:ilz,ifx:ilx) = f%vx(ifz:ilz,ifx:ilx) + tmp*lm%buox(ifz:ilz,ifx:ilx) *shot%rcv(ircv)%interp_coeff(:,:,1)
                    
                    case (4) !horizontal z adjsource
                    f%vz(ifz:ilz,ifx:ilx) = f%vz(ifz:ilz,ifx:ilx) + tmp*lm%buoz(ifz:ilz,ifx:ilx) *shot%rcv(ircv)%interp_coeff(:,:,1)
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (2) !horizontal x adjsource
                    !vx[ix-0.5,iy,iz]
                    f%vx(iz,ix) = f%vx(iz,ix) + tmp*lm%buox(iz,ix)
                    
                    case (4) !horizontal z adjsource
                    !vz[ix,iy,iz-0.5]
                    f%vz(iz,ix) = f%vz(iz,ix) + tmp*lm%buoz(iz,ix)
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    !v^it+1 -> v^it by FD^T of s^it+0.5
    subroutine update_velocities_adjoint(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl
        
        !same code with -dt
        call update_velocities(time_dir,it,f,bl)
        
    end subroutine
    
    !get v^it or s^it+0.5
    subroutine get_field_adjoint(f,w)
        type(t_field) :: f
        real w
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (1)
                w = sum( &
        ( lm%ldap2mu(ifz:ilz,ifx:ilx)*f%sxx(ifz:ilz,ifx:ilx) -lm%lda    (ifz:ilz,ifx:ilx)*f%szz(ifz:ilz,ifx:ilx)   &
         -lm%lda    (ifz:ilz,ifx:ilx)*f%sxx(ifz:ilz,ifx:ilx) +lm%ldap2mu(ifz:ilz,ifx:ilx)*f%szz(ifz:ilz,ifx:ilx) ) &
        *lm%inv_ldapmu_4mu(ifz:ilz,ifx:ilx)  *shot%src%interp_coeff(:,:,1) )
                
                case (2)
                w = sum( f%vx(ifz:ilz,ifx:ilx) *shot%src%interp_coeff(:,:,1) )
                
                case (4)
                w = sum( f%vz(ifz:ilz,ifx:ilx) *shot%src%interp_coeff(:,:,1) )
                
            end select
            
        else
            select case (shot%src%icomp)
                case (1) !s[iz,ix,iy]
                w = ( lm%ldap2mu(iz,ix)*f%sxx(iz,ix) -lm%lda    (iz,ix)*f%szz(iz,ix)   &
                     -lm%lda    (iz,ix)*f%sxx(iz,ix) +lm%ldap2mu(iz,ix)*f%szz(iz,ix) ) &
                    *lm%inv_ldapmu_4mu(iz,ix)
                
                case (2) !vx[iz,ix-0.5,iy]
                w=f%vx(iz,ix)
                
                case (4) !vz[iz-0.5,ix,iy]
                w=f%vz(iz,ix)
            end select
            
        endif
        
    end subroutine
    
    !========= for wavefield correlation ===================
    !dC/dm = \int a  dN/dm du/dt
    !      = \int a  dN/dm M Du
    !      = \int a -N dM/dm Du
    !                            [0                    ]
    !                            |  0                  |
    !-NdM/dlda = -1/(lda+mu)/4mu |    lda+2mu  -lda    |
    !                            |     -lda   lda+2mu  |
    !                            [                    0]
    ! => corr_lda := [adjsxx adjszz][lda+2mu  -lda  ][dvx_dx]
    !                               [ -lda   lda+2mu][dvz_dz]
    !                           [0                          ]
    !                           |  0                        |
    !-NdM/dmu = -2/(lda+mu)/4mu |    lda+2mu                |
    !                           |           lda+2mu         |
    !                           [                  2(lda+mu)]
    !                                     [lda+2mu                ][   dvx_dx    ]
    ! => corr_mu := [adjsxx adjszz adjsxz][       lda+2mu         ][   dvz_dz    ]
    !                                     [              2(lda+mu)][dvx_dz+dvz_dx]
    !            [b    ]
    !dN/drho M = |  b  |
    !            [    0]
    ! => corr_rho := [adjvx adjvz][sxx_dx+sxz_dz]
    !                             [szz_dz+sxz_dx]
    !
    !For gmu, we need to spatially interpolate adjsxz and dvx_dz+dvz_dx.
    !For grho, we need to spatially interpolate adjvx, adjvz and sxx_dx+sxz_dz, szz_dz+sxz_dx
    !
    !For time, we just use v^it+1, s^it+0.5, adjs^it+0.5
    
    subroutine field_correlation_moduli(it,sf,rf,sb,rb,corr)
        type(t_field) :: sf,rf
        integer,dimension(6) :: sb,rb
        real,dimension(cb%mz,cb%mx,3) :: corr
        
        !nonzero only when sf touches rf
        ifz=max(sb(1),rb(1),2)
        ilz=min(sb(2),rb(2),cb%mz-2)
        ifx=max(sb(3),rb(3),2)
        ilx=min(sb(4),rb(4),cb%mx-2)
        ify=max(sb(5),rb(5),2)
        ily=min(sb(6),rb(6),cb%my-2)
        
        call corr2d_flat_moduli(sf%vx,sf%vz,            &
                                rf%sxx,rf%szz,rf%sxz,   &
                                lm%ldap2mu,lm%lda,lm%mu,&
                                corr(:,:,1),corr(:,:,2),&
                                ifz,ilz,ifx,ilx)
        
    end subroutine
    
    subroutine field_correlation_density(it,sf,rf,sb,rb,corr)
        type(t_field) :: sf,rf
        integer,dimension(6) :: sb,rb
        real,dimension(cb%mz,cb%mx,3) :: corr
        
        !nonzero only when sf touches rf
        ifz=max(sb(1),rb(1),2)
        ilz=min(sb(2),rb(2),cb%mz-2)
        ifx=max(sb(3),rb(3),2)
        ilx=min(sb(4),rb(4),cb%mx-2)
        ify=max(sb(5),rb(5),2)
        ily=min(sb(6),rb(6),cb%my-2)
        
        call corr2d_flat_density(sf%sxx,sf%szz,sf%sxz,&
                                 rf%vx,rf%vz,         &
                                 corr(:,:,3),         &
                                 ifz,ilz,ifx,ilx)

    end subroutine
    
    subroutine field_correlation_scaling(corr)
        real,dimension(cb%mz,cb%mx,3) :: corr
        
        corr(:,:,1)=corr(:,:,1) * (-lm%inv_ldapmu_4mu(1:cb%mz,1:cb%mx))

        corr(:,:,2)=corr(:,:,2) * (-2.*lm%inv_ldapmu_4mu(1:cb%mz,1:cb%mx))

        corr(:,:,3)=corr(:,:,3) / cb%rho(1:cb%mz,1:cb%mx,1)


        corr(1,:,:) = corr(2,:,:)
        corr(cb%mz-1,:,:) = corr(cb%mz-2,:,:)
        corr(cb%mz,  :,:) = corr(cb%mz-2,:,:)

        corr(:,1,:) = corr(:,2,:)
        corr(:,cb%mx-1,:) = corr(:,cb%mx-2,:)
        corr(:,cb%mx  ,:) = corr(:,cb%mx-2,:)
        
        !set unit of gkpa to be [m3], grho to be [m5/s2]
        !such that after multiplied by (lda_max-lda_min) or (rho_max-rho_min) (will be done in m_parameterization.f90)
        !the unit of parameter update is [Nm], same as Lagrangian
        !and the unit of gradient scaling factor is [1/N/m] (in m_scaling.f90)
        !therefore parameters become unitless
        corr=corr*m%cell_volume*shot%src%dt
        
    end subroutine
    
    !========= Finite-Difference on flattened arrays ==================
    
    subroutine fd2d_flat_velocities(vx,vz,sxx,szz,sxz,                                  &
                                    cpml_dsxx_dx,cpml_dszz_dz,cpml_dsxz_dx,cpml_dsxz_dz,&
                                    buox,buoz,                                          &
                                    ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vx,vz,sxx,szz,sxz
        real,dimension(*) :: cpml_dsxx_dx,cpml_dszz_dz,cpml_dsxz_dx,cpml_dsxz_dz
        real,dimension(*) :: buox,buoz
        
        nz=cb%nz
        nx=cb%nx
        
        dsxx_dx=0.
        dsxz_dz=0.
        dsxz_dx=0.
        dszz_dz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dsxx_dx,dsxz_dz,dsxz_dx,dszz_dz)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1
                
                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2=i  -2*nz  !iz,ix-2
                iz_ixm1=i    -nz  !iz,ix-1
                iz_ixp1=i    +nz  !iz,ix+1
                iz_ixp2=i  +2*nz  !iz,ix+2
                

                dsxx_dx= c1x*(sxx(iz_ix)-sxx(iz_ixm1)) +c2x*(sxx(iz_ixp1)-sxx(iz_ixm2))
                dsxz_dz= c1z*(sxz(izp1_ix)-sxz(iz_ix)) +c2z*(sxz(izp2_ix)-sxz(izm1_ix))

                dsxz_dx= c1x*(sxz(iz_ixp1)-sxz(iz_ix)) +c2x*(sxz(iz_ixp2)-sxz(iz_ixm1))
                dszz_dz= c1z*(szz(iz_ix)-szz(izm1_ix)) +c2z*(szz(izp1_ix)-szz(izm2_ix))
                                
                !cpml
                cpml_dsxx_dx(i)= cb%b_x(ix)*cpml_dsxx_dx(i) + cb%a_x(ix)*dsxx_dx
                cpml_dsxz_dz(i)= cb%b_z(iz)*cpml_dsxz_dz(i) + cb%a_z(iz)*dsxz_dz 
                cpml_dsxz_dx(i)= cb%b_x(ix)*cpml_dsxz_dx(i) + cb%a_x(ix)*dsxz_dx
                cpml_dszz_dz(i)= cb%b_z(iz)*cpml_dszz_dz(i) + cb%a_z(iz)*dszz_dz

                dsxx_dx=dsxx_dx*k_x + cpml_dsxx_dx(i)
                dsxz_dz=dsxz_dz*k_z + cpml_dsxz_dz(i)
                dsxz_dx=dsxz_dx*k_x + cpml_dsxz_dx(i)
                dszz_dz=dszz_dz*k_z + cpml_dszz_dz(i)
                
                !velocity
                vx(i)=vx(i) + dt*buox(i)*(dsxx_dx+dsxz_dz)
                vz(i)=vz(i) + dt*buoz(i)*(dsxz_dx+dszz_dz)
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine fd2d_flat_stresses(vx,vz,sxx,szz,sxz,                              &
                                  cpml_dvx_dx,cpml_dvz_dz,cpml_dvx_dz,cpml_dvz_dx,&
                                  ldap2mu,lda,mu,                                 &
                                  ifz,ilz,ifx,ilx,dt)
        real,dimension(*) :: vx,vz,sxx,szz,sxz
        real,dimension(*) :: cpml_dvx_dx,cpml_dvz_dz,cpml_dvx_dz,cpml_dvz_dx
        real,dimension(*) :: ldap2mu,lda,mu
        
        nz=cb%nz
        nx=cb%nx
        
        dvx_dx=0.
        dvz_dz=0.
        dvx_dz=0.
        dvz_dx=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,iz_ixp1,iz_ixp2,&
        !$omp         dvx_dx,dvz_dz,&
        !$omp         dvx_dz,dvz_dx)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*nz+1
                
                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2=i  -2*nz !iz,ix-2
                iz_ixm1=i    -nz !iz,ix-1
                iz_ixp1=i    +nz !iz,ix+1
                iz_ixp2=i  +2*nz !iz,ix+2
                

                dvx_dx= c1x*(vx(iz_ixp1)-vx(iz_ix))  +c2x*(vx(iz_ixp2)-vx(iz_ixm1))
                dvz_dz= c1z*(vz(izp1_ix)-vz(iz_ix))  +c2z*(vz(izp2_ix)-vz(izm1_ix))
                
                !cpml
                cpml_dvx_dx(i)=cb%b_x(ix)*cpml_dvx_dx(i)+cb%a_x(ix)*dvx_dx
                cpml_dvz_dz(i)=cb%b_z(iz)*cpml_dvz_dz(i)+cb%a_z(iz)*dvz_dz

                dvx_dx=dvx_dx*k_x + cpml_dvx_dx(iz_ix)
                dvz_dz=dvz_dz*k_z + cpml_dvz_dz(iz_ix)
                
                !normal stresses
                sxx(i) = sxx(i) + dt * (ldap2mu(i)*dvx_dx+lda(i)    *dvz_dz)
                szz(i) = szz(i) + dt * (lda(i)    *dvx_dx+ldap2mu(i)*dvz_dz)


                dvx_dz= c1z*(vx(iz_ix)-vx(izm1_ix))  +c2z*(vx(izp1_ix)-vx(izm2_ix))
                dvz_dx= c1x*(vz(iz_ix)-vz(iz_ixm1))  +c2x*(vz(iz_ixp1)-vz(iz_ixm2))

                !cpml
                cpml_dvx_dz(i)=cb%b_z(iz)*cpml_dvx_dz(i)+cb%a_z(iz)*dvx_dz
                cpml_dvz_dx(i)=cb%b_x(ix)*cpml_dvz_dx(i)+cb%a_x(ix)*dvz_dx

                dvx_dz=dvx_dz*k_z + cpml_dvx_dz(i)
                dvz_dx=dvz_dx*k_x + cpml_dvz_dx(i)

                !shear stress
                sxz(i) = sxz(i) + dt * mu(i)*(dvx_dz+dvz_dx)
                
            enddo
            
        enddo
        !$omp enddo 
        !$omp end parallel
        
    end subroutine
    
    subroutine corr2d_flat_moduli(sf_vx,sf_vz,         &
                                  rf_sxx,rf_szz,rf_sxz,&
                                  ldap2mu,lda,mu,      &
                                  corr_lda,corr_mu,    &
                                  ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_vx,sf_vz
        real,dimension(*) :: rf_sxx,rf_szz,rf_sxz
        real,dimension(*) :: ldap2mu,lda,mu
        real,dimension(*) :: corr_lda,corr_mu
        
        nz=cb%nz
        
        sf_dvx_dx=0.
        sf_dvz_dz=0.
        sf_4_dvzdx_p_dvxdz=0.
        rf_4_sxz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,izp1_ixm2,iz_ixm1,izp1_ixm1,&
        !$omp         izm2_ixp1,izm1_ixp1,iz_ixp1,izp1_ixp1,izp2_ixp1,&
        !$omp         iz_ixp2,izp1_ixp2,&
        !$omp         sf_dvx_dx,sf_dvz_dz,&
        !$omp         sf_4_dvzdx_p_dvxdz,rf_4_sxz)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
                
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz !corr has no boundary layers

                izm2_ix= i-2  !iz-2,ix
                izm1_ix= i-1  !iz-1,ix
                iz_ix  = i    !iz  ,ix
                izp1_ix= i+1  !iz+1,ix
                izp2_ix= i+2  !iz+2,ix

                iz_ixm2  = i    -2*nz  !iz  ,ix-2
                izp1_ixm2= i+1  -2*nz  !iz+1,ix-2
                iz_ixm1  = i      -nz  !iz  ,ix-1
                izp1_ixm1= i+1    -nz  !iz+1,ix-1

                izm2_ixp1= i-2    +nz  !iz-2,ix+1
                izm1_ixp1= i-1    +nz  !iz-1,ix+1
                iz_ixp1  = i      +nz  !iz  ,ix+1
                izp1_ixp1= i+1    +nz  !iz+1,ix+1
                izp2_ixp1= i+2    +nz  !iz+2,ix+1

                iz_ixp2  = i    +2*nz  !iz  ,ix+2
                izp1_ixp2= i+1  +2*nz  !iz+1,ix+2
                

                sf_dvx_dx = c1x*(sf_vx(iz_ixp1)-sf_vx(iz_ix)) +c2x*(sf_vx(iz_ixp2)-sf_vx(iz_ixm1))
                sf_dvz_dz = c1z*(sf_vz(izp1_ix)-sf_vz(iz_ix)) +c2z*(sf_vz(izp2_ix)-sf_vz(izm1_ix))

                corr_lda(j)=corr_lda(j) + ldap2mu(i)*rf_sxx(i)*sf_dvx_dx -     lda(i)*rf_sxx(i)*sf_dvz_dz &
                                        -     lda(i)*rf_szz(i)*sf_dvx_dx + ldap2mu(i)*rf_szz(i)*sf_dvz_dz


                sf_4_dvzdx_p_dvxdz = &  !(dvz_dx+dvx_dz)(iz,ix)     [iz-0.5,ix-0.5]
                                         c1x*(sf_vz(iz_ix)-sf_vz(iz_ixm1)) +c2x*(sf_vz(iz_ixp1)-sf_vz(iz_ixm2)) &
                                       + c1z*(sf_vx(iz_ix)-sf_vx(izm1_ix)) +c2z*(sf_vx(izp1_ix)-sf_vx(izm2_ix)) &
                                     & &
                                     & & !(dvz_dx+dvx_dz)(iz,ix+1)   [iz-0.5,ix+0.5]
                                       + c1x*(sf_vz(iz_ixp1)-sf_vz(iz_ix    )) +c2x*(sf_vz(iz_ixp2  )-sf_vz(iz_ixm1  )) &
                                       + c1z*(sf_vx(iz_ixp1)-sf_vx(izm1_ixp1)) +c2z*(sf_vx(izp1_ixp1)-sf_vx(izm2_ixp1)) &
                                     & &
                                     & & !(dvz_dx+dvx_dz)(iz+1,ix)   [iz+0.5,ix-0.5]
                                       + c1x*(sf_vz(izp1_ix)-sf_vz(izp1_ixm1)) +c2x*(sf_vz(izp1_ixp1)-sf_vz(izp1_ixm2)) &
                                       + c1z*(sf_vx(izp1_ix)-sf_vx(iz_ix    )) +c2z*(sf_vx(izp2_ix  )-sf_vx(izm1_ix  )) &
                                     & &
                                     & & !(dvz_dx+dvx_dz)(iz+1,ix+1) [iz+0.5,ix+0.5]
                                       + c1x*(sf_vz(izp1_ixp1)-sf_vz(izp1_ix)) +c2x*(sf_vz(izp1_ixp2)-sf_vz(izp1_ixm1)) &
                                       + c1z*(sf_vx(izp1_ixp1)-sf_vx(iz_ixp1)) +c2z*(sf_vx(izp2_ixp1)-sf_vx(izm1_ixp1))

                rf_4_sxz = rf_sxz(iz_ix) + rf_sxz(izp1_ix) + rf_sxz(iz_ixp1) + rf_sxz(izp1_ixp1)

                corr_mu(j)=corr_mu(j) + ldap2mu(i)*rf_sxx(i)*sf_dvx_dx &
                                      + ldap2mu(i)*rf_szz(i)*sf_dvz_dz &
                                      + 2.*(lda(i)+mu(i))*0.0625*rf_4_sxz*sf_4_dvzdx_p_dvxdz
                
            end do
            
        end do
        !$omp end do
        !$omp end parallel

    end subroutine

    subroutine corr2d_flat_density(sf_sxx,sf_szz,sf_sxz,&
                                   rf_vx,rf_vz,         &
                                   corr,                &
                                   ifz,ilz,ifx,ilx)
        real,dimension(*) :: sf_sxx,sf_szz,sf_sxz
        real,dimension(*) :: rf_vx,rf_vz
        real,dimension(*) :: corr
        
        nz=cb%nz
        
        sf_2_dsxxdx_p_dsxzdz=0.
        sf_2_dsxzdx_p_dszzdz=0.
        rf_2vx=0.
        rf_2vz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,i,j,&
        !$omp         izm2_ix,izm1_ix,iz_ix,izp1_ix,izp2_ix,&
        !$omp         iz_ixm2,iz_ixm1,izp1_ixm1,&
        !$omp         izm1_ixp1,iz_ixp1,izp1_ixp1,izp2_ixp1,&
        !$omp         iz_ixp2,izp1_ixp2,&
        !$omp         sf_2_dsxxdx_p_dsxzdz,sf_2_dsxzdx_p_dszzdz,&
        !$omp         rf_2vx,rf_2vz)
        !$omp do schedule(dynamic)
        do ix=ifx,ilx
        
            !dir$ simd
            do iz=ifz,ilz
            
                i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+1 !field has boundary layers
                j=(iz-1)     +(ix-1)     *cb%mz+1 !corr has no boundary layers
                
                izm2_ix=i-2  !iz-2,ix
                izm1_ix=i-1  !iz-1,ix
                iz_ix  =i    !iz,ix
                izp1_ix=i+1  !iz+1,ix
                izp2_ix=i+2  !iz+2,ix
                
                iz_ixm2  = i    -2*nz  !iz  ,ix-2
                iz_ixm1  = i      -nz  !iz  ,ix-1
                izp1_ixm1= i+1    -nz  !iz-1,ix-1

                izm1_ixp1= i-1    +nz  !iz-1,ix+1
                iz_ixp1  = i      +nz  !iz  ,ix+1
                izp1_ixp1= i+1    +nz  !iz+1,ix+1
                izp2_ixp1= i+2    +nz  !iz+2,ix+1

                iz_ixp2  = i    +2*nz  !iz ,ix+2
                izp1_ixp2= i+1  +2*nz  !iz+1,ix+2
                
                sf_2_dsxxdx_p_dsxzdz = &   !(dsxx_dx+dsxz_dz)(iz,ix)   [iz,ix-0.5]
                                           c1x*(sf_sxx(iz_ix  )-sf_sxx(iz_ixm1)) +c2x*(sf_sxx(iz_ixp1)-sf_sxx(iz_ixm2)) &
                                         + c1z*(sf_sxz(izp1_ix)-sf_sxz(iz_ix  )) +c2z*(sf_sxz(izp2_ix)-sf_sxz(izm1_ix)) &
                                       & &
                                       & & !(dsxx_dx+dsxz_dz)(iz,ix+1) [iz,ix+0.5]
                                         + c1x*(sf_sxx(iz_ixp1  )-sf_sxx(iz_ix  )) +c2x*(sf_sxx(iz_ixp2  )-sf_sxx(iz_ixm1  )) &
                                         + c1z*(sf_sxz(izp1_ixp1)-sf_sxz(iz_ixp1)) +c2z*(sf_sxz(izp2_ixp1)-sf_sxz(izm1_ixp1))

                sf_2_dsxzdx_p_dszzdz = &  !(dsxz_dx+dszz_dz)(iz,ix)   [iz-0.5,ix]
                                           c1x*(sf_sxz(iz_ixp1)-sf_sxz(iz_ix  )) +c2x*(sf_sxz(iz_ixp2)-sf_sxz(iz_ixm1)) &
                                         + c1z*(sf_szz(iz_ix  )-sf_szz(izm1_ix)) +c2z*(sf_szz(izp1_ix)-sf_szz(izm1_ix)) &
                                       & &
                                       & & !(dsxz_dx+dszz_dz)(iz+1,ix) [iz+0.5,ix]
                                         + c1x*(sf_sxz(izp1_ixp1)-sf_sxz(izp1_ix)) +c2x*(sf_sxz(izp1_ixp2)-sf_sxz(izp1_ixm1)) &
                                         + c1z*(sf_szz(izp1_ix  )-sf_szz(iz_ix  )) +c2z*(sf_szz(izp2_ix  )-sf_szz(iz_ix    ))

                rf_2vx = rf_vx(iz_ix) + rf_vx(iz_ixp1)
                rf_2vz = rf_vz(iz_ix) + rf_vz(izp1_ix)
                
                corr(j)=corr(j) + 0.25*( rf_2vx*sf_2_dsxxdx_p_dsxzdz + rf_2vz*sf_2_dsxzdx_p_dszzdz )
                
            enddo
            
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    

end
