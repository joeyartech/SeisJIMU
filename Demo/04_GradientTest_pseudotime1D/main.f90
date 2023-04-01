program main
use m_pseudotime

    integer,parameter :: nz=373 !101 !373
    real,parameter :: dz=13.375 !20. !13.375
    integer,parameter :: iwaterbottom=20, ifreeze=iwaterbottom+1
    
    integer,dimension(1,1) :: nt_real, nz_real
    
    real,dimension(:,:,:),allocatable :: v_z, v_t !velocity
    real,dimension(:,:,:),allocatable :: v_z2, v_t2 !copies
    real,dimension(:,:,:),allocatable :: g_z, g_t !gradient
    real,dimension(:,:,:),allocatable :: eta !perturbation
    
    real,dimension(:,:,:),allocatable :: v_z_pert, v_t_pert
    
        
    allocate(v_z(nz,1,1))
    
    open(11,file='v_z',access='direct',recl=4*nz)
    read(11,rec=1) v_z
    close(11)
    vmin=minval(v_z); vmax=maxval(v_z)

! v_z=100.
! v_z(nz-60:nz-10,:,:)=100.
! v_z(nz,:,:)=100.
! ! do iz=50,nz
! ! v_z(iz,:,:)=  100+(iz-50)/(nz-50)*100.
! ! enddo
! vmin=50.
! vmax=500.
    
    call pseudotime_init('z->t',vmin,vmax,nx_=1,ny_=1,&
        nz_=nz, dz_=dz,&
        nt_=nt, dt_=dt)
    print*,'nz, dz =', nz, dz
    print*,'nt, dt =', nt, dt
    
    call pseudotime_convert('z->t',v_z,v_t,o_nreal=nt_real)
    open(11,file='v_t',access='direct',recl=4*nt)
    write(11,rec=1) v_t
    close(11)    
    print*,'real nt needed #1 =', nt_real
        
!     !Test #1: Reproducibility
!     open(12,file='v_z_reproduced',access='direct',recl=4*nz)
!     open(13,file='v_t_reproduced',access='direct',recl=4*nt)
!     write(12,rec=1) v_z
!     write(13,rec=1) v_t
!     
!     call pseudotime_convert('t->z',v_t,v_z2,o_nreal=nz_real)
!     v_z2(1:ifreeze,:,:)=v_z(1:ifreeze,:,:) !apply freeze zone
!     call pseudotime_convert('z->t',v_z2,v_t2,o_nreal=nt_real)
!     print*,'returned nz #2 =', nz_real
!     print*,'real nt needed #2 =', nt_real
!     write(12,rec=2) v_z2
!     write(13,rec=2) v_t2
!     
!     call pseudotime_convert('t->z',v_t2,v_z2,o_nreal=nz_real)
!     v_z2(1:ifreeze,:,:)=v_z(1:ifreeze,:,:) !apply freeze zone
!     call pseudotime_convert('z->t',v_z2,v_t2,o_nreal=nt_real)
!     print*,'returned nz #3 =', nz_real
!     print*,'real nt needed #3 =', nt_real
!     write(12,rec=3) v_z2
!     write(13,rec=3) v_t2
!     
!     call pseudotime_convert('t->z',v_t2,v_z2,o_nreal=nz_real)
!     v_z2(1:ifreeze,:,:)=v_z(1:ifreeze,:,:) !apply freeze zone
!     call pseudotime_convert('z->t',v_z2,v_t2,o_nreal=nt_real)
!     print*,'returned nz #4 =', nz_real
!     print*,'real nt needed #4 =', nt_real
!     write(12,rec=4) v_z2 !diverge
!     write(13,rec=4) v_t2 !diverge
!     
!     close(12)
!     close(13)
    

    !Task #2: Validation of gradient in time domain
    ! lim_α→0 (C[v+α*η]-C[v])/α =: ∂C/∂η = ∫KᵥC η dz
    ! Kernel: KᵥC; Gradient: ∇ᵥC:=KᵥC*dz
    ! Matrix-vector product: ∂C/∂η=∇ᵥC*η
    ! Validation:
    ! ∂C/∂η /dz ?= converted KᵥC
    !            = K_{v(z)}(z)v(z) - ∫_z^∞ K_{v(z)}(z')dv(t')

    call pseudotime_convert_gradient(gradient(v_z,dz),v_z,g_t)
    open(12,file='g_t_convert',access='direct',recl=4*nt)
    write(12,rec=1) g_t
    close(12)
    
    eta=g_t*dt !/50.  !g_t/dt
!     open(12,file='eta',access='direct',recl=4*nt)
!     write(12,rec=1) eta
!     close(12)
    
    print*,'∫₀^∞ K_vt η dt =', sum(g_t*eta)
    print*,'∂C/∂η = lim_α→0 (C[v(t)+α*eta]-C[v(t)])/α'
    print*,'      = 0.5*v(T)²∂T/∂η + ∫₀^T vη dt '
    print*,'     /= ∫₀^T vη dt =',sum(v_t(1:nt_real(1,1),1,1)*eta(1:nt_real(1,1),1,1))*dt

    C=objective(v_z,dz,nz)
    print*,'C=', C
    
    !limit
    open(12,file='v_z_pert',access='direct',recl=4*nz)
    write(12,rec=1) v_z
    
    print*,'∂C/∂η ≈'
    print*, 'alpha            C_pert           (C_pert-C)/alpha   ratio'
    do i=10,1,-1
        alpha=i
        
        call pseudotime_convert('t->z',v_t+alpha*eta,v_z_pert)
        write(12,rec=i+1) v_z_pert
                
        C_pert  =objective(v_z_pert,dz,nz)

        print*, alpha,C_pert, (C_pert - C)/alpha, (C_pert - C)/alpha/sum(g_t*eta)
        
    enddo    
    close(12)
    
    contains
    
    real function objective(v_z,dz,nz)
        real,dimension(:,:,:) :: v_z
        objective=0.5*sum(v_z(1:nz,:,:)**2)*dz
    end function
    
    function gradient(v_z,dz)
        real,dimension(:,:,:) :: v_z
        real,dimension(:,:,:),allocatable :: gradient
        gradient=v_z*dz
    end function
    
end
