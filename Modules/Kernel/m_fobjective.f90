module m_fobjective
use m_System
use m_Modeling
use m_Lpnorm
use m_matchfilter
use m_smoother_laplacian_sparse
use m_parametrizer
use m_querypoint
use m_weighter
use m_preconditioner

    private

    type,public :: t_fobjective

        type(t_string),dimension(:),allocatable :: s_dnorms !norms for data residuals
        type(t_string),dimension(:),allocatable :: s_xnorms !norms for (unknown) parameter residuals

        integer :: n_dnorms, n_xnorms !number of norms

        real,dimension(:),allocatable :: dnorms, xnorms !values of norms
        real,dimension(:),allocatable :: dnorm_normalizers, xnorm_normalizers
        real,dimension(:),allocatable :: dnorm_weights,     xnorm_weights

        real FWI_misfit

        contains

        procedure :: init
        procedure :: stack_dnorms
        procedure :: stack_xnorms
        procedure :: print_dnorms
        procedure :: print_xnorms
        procedure :: eval

    end type

    type(t_fobjective),public :: fobj

    logical :: if_has_dnorm_normalizers, if_has_xnorm_normalizers

    real,dimension(:,:),allocatable :: Wdres
    real,dimension(:,:,:,:),allocatable :: Wxres
    real,dimension(:,:,:,:),allocatable :: xprior !prior parameter
    ! real :: xnorms_weights(4) !weights on 1st, 2nd, 3rd & 4th dimensions

    contains

    subroutine init(self)
        class(t_fobjective) :: self
        
        !data norms
        self%s_dnorms=setup%get_strs('DATA_NORMS','DNORMS',o_default='L1 L2sq')
        self%n_dnorms=size(self%s_dnorms)
        
        call alloc(self%dnorms,           self%n_dnorms)
        
        self%dnorm_weights=setup%get_reals('DATA_NORM_WEIGHTS','DWEIGHTS',o_default='0. 1.')
        if(size(self%dnorm_weights)/=self%n_dnorms) then
            stop 'error reading DATA_NORM_WEIGHTS (DWEIGHTS)'
        endif

        if_has_dnorm_normalizers=.false.
        call alloc(self%dnorm_normalizers,self%n_dnorms,o_init=1.)
        if(setup%check('DATA_NORM_NORMALIZERS','DNORMALIZERS')) then
            self%dnorm_normalizers=setup%get_reals('DATA_NORM_NORMALIZERS','DNORMALIZERS',o_mandatory=self%n_dnorms)
            if_has_dnorm_normalizers=.true.
        endif

        !parameter norms
        self%s_xnorms=setup%get_strs('PARAMETER_NORMS','XNORMS',o_default='none')
        self%n_xnorms=size(self%s_xnorms)
        call alloc(self%xnorms,self%n_xnorms)

        self%xnorm_weights=setup%get_reals('PARAMETER_NORM_WEIGHTS','XWEIGHTS',o_default='0.')

        ! if(self%n_xnorms>0) then
        !     xnorms_weights=setup%get_reals('PARAMETER_NORMS_WEIGHTS','XNORMS_WEI',o_default='1. 1. 1. 1.')
        ! endif

        !prior models
        call m%read_prior

        if(m%if_has_prior) then
            call alloc(xprior,param%n1,param%n2,param%n3,param%npars)
            !transform prior models
            call param%transform('m->x',o_xprior=xprior)
            !save some RAM
            call dealloc(m%vp_prior,m%vs_prior,m%rho_prior)
        endif
        
    end subroutine
    
    recursive subroutine stack_dnorms(self)
    use mpi
        class(t_fobjective) :: self

        real,dimension(:,:),allocatable :: Wdres
        logical :: is_4adjsrc
        real :: numer, denom, time(shot%nt), v(shot%nt)
        
        !In theory, (L2 norm sq for example)
        !C(u) = ½║u║₂² = ½Σ_xr ∫ (u-d)² δ(x-xr) dx3 dt, and
        !∇C = Σ_xr ∫ (u-d)² δ(x-xr)
        !where δ(x-xr) is the Dirac delta function centered at receiver positions.
        !In the discretized world, δ(x-xr)={0,1}/dx3,
        !in the implementation we take
        !C(u) = 0.5*Σ_xr ∫ (u-d)² dt, and
        !the presence of δ(x-xr) in the adjoint source is accounted for
        !by the RHS of adjoint weq
        !(ie. when calling rfield%ignite, to be same as sfield%ignite)

        !reinitialize adjoint source
        if(if_has_dnorm_normalizers) call alloc(shot%dadj,shot%nt,shot%nrcv)

        do i=1,self%n_dnorms
            select case (self%s_dnorms(i)%s)

            case ('L1')

                is_4adjsrc = self%dnorm_weights(i)>0. .and. if_has_dnorm_normalizers

                do ir=1,shot%nrcv
                    if(shot%rcv(ir)%is_badtrace) cycle
                    
                    if(shot%rcv(ir)%comp=='p') then !pressure data, use 1/ref_vp to balance amplitudes versus velocities data
                        self%dnorms(i) = self%dnorms(i) + L1(self%dnorm_normalizers(i), shot%nt, &
                            wei%weight(:,ir)*m%ref_inv_vp, shot%dsyn(:,ir)-shot%dobs(:,ir), shot%dt)
                            
                        if(is_4adjsrc) call nabla_L1(shot%dadj(:,ir),oif_stack=.true.)
                        
                    else !velocities data, use ref_rho to balance amplitudes versus pressure data
                        self%dnorms(i) = self%dnorms(i) + L1(self%dnorm_normalizers(i), shot%nt, &
                            wei%weight(:,ir)*m%ref_rho,    shot%dsyn(:,ir)-shot%dobs(:,ir), shot%dt)

                        if(is_4adjsrc) call nabla_L1(shot%dadj(:,ir),oif_stack=.true.)

                    endif

                enddo
                
            case ('L2sq')
            
                is_4adjsrc = self%dnorm_weights(i)>0. .and. if_has_dnorm_normalizers

                do ir=1,shot%nrcv
                    if(shot%rcv(ir)%is_badtrace) cycle
                    
                    if(shot%rcv(ir)%comp=='p') then !pressure data
                        self%dnorms(i) = self%dnorms(i) + L2sq(0.5*self%dnorm_normalizers(i), shot%nt, &
                            wei%weight(:,ir)*m%ref_inv_vp, shot%dsyn(:,ir)-shot%dobs(:,ir), shot%dt)
                            
                        if(is_4adjsrc) call nabla_L2sq(shot%dadj(:,ir),oif_stack=.true.)
                        
                    else !velocity data
                        self%dnorms(i) = self%dnorms(i) + L2sq(0.5*self%dnorm_normalizers(i), shot%nt, &
                            wei%weight(:,ir)*m%ref_rho,    shot%dsyn(:,ir)-shot%dobs(:,ir), shot%dt)

                        if(is_4adjsrc) call nabla_L2sq(shot%dadj(:,ir),oif_stack=.true.)

                    endif

                enddo

            ! case ('ndecon') !normalized deconvolution
            !     ! !Data residual of 1st shot as an example
            !     ! if(mpiworld%is_master) then
            !     !     Wdres=wei%weight*(shot%dsyn-shot%dobs)
            !     !     call shot%write('Wdres_',Wdres)
            !     !     deallocate(Wdres)
            !     ! endif

            !     !aka Reverse AWI (Warner & Guasch, 2016, Geophysics, 81-6)
            !     !Note Eq D-4 in Appendix D has a minor typo, the other eq in the paper is correct.

            !     ! if(is_first_in) then
            !     !     ref_p = 0.
            !     !     do ir=1,shot%nrcv
            !     !         if(shot%rcv(ir)%is_badtrace) cycle
            !     !         if(shot%rcv(ir)%comp=='p') then !for pressure data
            !     !             tmp=maxval(abs(Wdres(:,ir)))
            !     !             ref_p=either(ref_p,tmp,ref_p>tmp)
            !     !         endif
            !     !     enddo
            !     ! endif

            !     is_4adjsrc=abs(self%dnorm_weights(i))>0. .and. if_compute_adjsrc
            !     if(is_4adjsrc) call alloc(shot%dadj,shot%nt,shot%nrcv,oif_protect=.true.)

            !     numer=0.
            !     denom=0.
            !     time=[(it-1.,it=1,shot%nt)]*shot%dt !why shot%nt?

            !     do ir=1,shot%nrcv
            !         if(shot%rcv(ir)%is_badtrace) cycle
            !         if(shot%rcv(ir)%comp=='p') then
            !             call matchfilter_estimate(shot%dobs(:,ir),shot%dsyn(:,ir),shot%nt,1,o_filter_time=v)
            !             numer=numer+norm2(time*v) 
            !             denom=denom+norm2(v*v)

            !             self%dnorms(i) = self%dnorms(i) + 0.5*numer/denom !/shot%dt*m%cell_volume*ref_p*ref_p/m%ref_kpa
                    
            !             if(is_4adjsrc) then
            !                 call matchfilter_adjointsrc(shot%dobs(:,ir),shot%dadj(:,ir), &
            !                     -(time*time-2.*self%dnorms(i))/denom*v)

            !             endif
                
            !         endif
            !     enddo

            end select

        enddo
        
        if(if_has_dnorm_normalizers) call shot%write('dadj_',shot%dadj)
        
        if(if_has_dnorm_normalizers) return
        
        
        !otherwise estimate normalizers
        if(mpiworld%is_master) then
            self%dnorm_normalizers=1./self%dnorms
            write(*,*) 'fobj%dnorm_normalizers before shotlist-scaling =',self%dnorm_normalizers
            call shls%scale(self%n_dnorms,o_from_1st_shot=self%dnorm_normalizers)
            write(*,*) 'fobj%dnorm_normalizers =',self%dnorm_normalizers
        endif
        
        call mpi_bcast(self%dnorm_normalizers, self%n_dnorms, mpi_real, 0, mpiworld%communicator, mpiworld%ierr)
        !dnorm_normalizers should also be checkpointed

        if_has_dnorm_normalizers=.true.
        
        !and redo the computations for adjoint sources
        self%dnorms=0.
        call self%stack_dnorms
        
    end subroutine

    subroutine stack_xnorms(self,x,g)
        class(t_fobjective) :: self
        real,dimension(param%n1,param%n2,param%n3,param%npars) :: x,g

        real,dimension(:,:,:,:),allocatable :: xres !weighted residuals of x
        real,dimension(:,:,:,:),allocatable :: spatgrad !spatial gradient of x

        call alloc(spatgrad,param%n1,param%n2,param%n3,param%npars,o_init=0.)

        self%xnorms=0.

        do i=1,self%n_xnorms
            select case (self%s_xnorms(i)%s)
            case ('Ridge','L2sq')
                !fcost
                call alloc(xres,param%n1,param%n2,param%n3,param%npars)
                if(m%if_has_prior) then
                    xres=x-xprior
                else
                    xres=x
                endif

                self%xnorms(i) = self%xnorms(i) + sum(x*x)
                !self%xnorms(i) = sum(Wxres*Wxres, .not.param%is_freeze_zone)

                !spatial gradient
                spatgrad = spatgrad + xres

                deallocate(xres)

            case ('Lasso','L1')
                !to be implemented..

            case ('1st_Ridge','1st_L2sq')

                do j=1,param%n1
                    !self%xnorms(i) = self%xnorms(i) + xnorms_weights(1)*sum((x(j+1,:,:,:)-x(j,:,:,:))**2,.not.param%is_freeze_zone(j,:,:,:))
                    self%xnorms(i) = self%xnorms(i) + sum((x(j+1,:,:,:)-x(j,:,:,:))**2)
                enddo

                do j=1,param%n2
                    self%xnorms(i) = self%xnorms(i) + sum((x(:,j+1,:,:)-x(:,j,:,:))**2)
                enddo

                if(m%is_cubic) then
                do j=1,param%n3
                    self%xnorms(i) = self%xnorms(i) + sum((x(:,:,j+1,:)-x(:,:,j,:))**2)
                enddo
                endif

                !spatial gradient
                do i4=1,param%npars
                do i3=1,param%n3
                do i2=1,param%n2
                do i1=1,param%n1
                    ! spatgrad(i1,i2,i3,i4) = spatgrad(i1,i2,i3,i4) &
                    !     -xnorms_weights(1)*(x(i1-1,i2,i3,i4)-2.*x(i1,i2,i3,i4)+x(i1+1,i2,i3,i4)) &
                    !     -xnorms_weights(2)*(x(i1,i2-1,i3,i4)-2.*x(i1,i2,i3,i4)+x(i1,i2+1,i3,i4)) &
                    !     -xnorms_weights(3)*(x(i1,i2,i3-1,i4)-2.*x(i1,i2,i3,i4)+x(i1,i2,i3+1,i4))
                    spatgrad(i1,i2,i3,i4) = spatgrad(i1,i2,i3,i4) &
                        -x(i1-1,i2,i3,i4)-2.*x(i1,i2,i3,i4)+x(i1+1,i2,i3,i4) &
                        -x(i1,i2-1,i3,i4)-2.*x(i1,i2,i3,i4)+x(i1,i2+1,i3,i4) &
                        -x(i1,i2,i3-1,i4)-2.*x(i1,i2,i3,i4)+x(i1,i2,i3+1,i4)
                enddo
                enddo
                enddo
                enddo

            case ('TV')
                !to be implemented

            end select

        enddo

        self%xnorms = 0.5*self%xnorms *m%cell_volume

        !where(param%is_freeze_zone) spatgrad=0.

        !add to dnorm's gradient wrt parameters
        g = g + spatgrad *m%cell_volume

    end subroutine

    subroutine print_dnorms(self,prefix,suffix)
        class(t_fobjective) :: self
        character(*) :: prefix,suffix

        if(mpiworld%is_master) then
            write(*,*) prefix//' fobj%dnorms '//suffix//': ', (' ',self%s_dnorms(i)%s//'= ',self%dnorms(i),i=1,self%n_dnorms)
        endif

    end subroutine

    subroutine print_xnorms(self)
        class(t_fobjective) :: self

        if(mpiworld%is_master) then
            write(*,*) 'fobj%xnorms:', (' ',self%s_xnorms(i)%s//'= ',self%xnorms(i),i=1,self%n_xnorms)
        endif

    end subroutine

    subroutine eval(self,qp,oif_update_m,oif_approx,oif_gradient)
        class(t_fobjective) :: self
        type(t_querypoint) :: qp
        logical,optional :: oif_update_m,oif_approx,oif_gradient

        type(t_string),dimension(:),allocatable :: smoothings
        character(:),allocatable :: smask
        real,dimension(:,:,:),allocatable :: mask

        real,dimension(:,:,:),allocatable ::  FWI, RE, DR
        !real,dimension(:),allocatable :: tmp
        real,dimension(:),allocatable :: precond

        !update model
        if(either(oif_update_m,.true.,present(oif_update_m))) then
            call param%transform(o_dir='x->m',o_x=qp%x)
        endif
        
        !compute dnorm's gradient by adjoint-state method
        ! if(either(oif_approx,.false.,present(oif_approx))) then
        !     call modeling_gradient_approximate(fobj)
        ! else
            call modeling_gradient(qp%is_fitting_data)
        ! endif
        
        qp%f = sum(self%dnorm_weights*self%dnorms) ! + sum(self%xnorm_weights*self%xnorms)

        if(.not. either(oif_gradient,.true.,present(oif_gradient))) return

        !post-modeling smoothing in physical (model) domain
        smoothings=setup%get_strs('SMOOTHING','SMTH',o_default='Laplacian')

        do i=1,size(smoothings)
            !Laplacian smoothing
            if(smoothings(i)%s=='Laplacian') then
                call hud('Laplacian smoothing')
                call smoother_laplacian_init([m%nz,m%nx,m%ny],[m%dz,m%dx,m%dy],shot%fpeak)
                do j=1,ppg%ngrad
                    !call smoother_laplacian_extend_mirror(m%gradient(:,:,:,i),m%itopo)
                    call smoother_laplacian_pseudo_nonstationary(m%gradient(:,:,:,j),m%vp)
                enddo    
            endif
        enddo

        !freeze_zone as hard mask
        do i=1,ppg%ngrad
            where(m%is_freeze_zone) m%gradient(:,:,:,i)=0.
        enddo

print*,'Im here'
if(setup%get_bool('IF_STICK_TO_RE',o_default='false')) then
call hud('Stick to RE: qp%pg -> RE')
call alloc(qp%pg,m%nz,m%nx,m%ny,1)
qp%pg(:,:,:,1) = (m%correlate(:,:,:,2)+m%correlate(:,:,:,3))
where(m%is_freeze_zone) qp%pg(:,:,:,1)=0.
endif
print*,'Im here2'

        !soft mask
        smask=setup%get_file('GRADIENT_SOFT_MASK','MASK')
        if(smask/='') then
            call alloc(mask,m%nz,m%nx,m%ny)
            ! call sysio_read(smask,mask,size(mask))
            
            do i=1,ppg%ngrad
                m%gradient(:,:,:,i)=m%gradient(:,:,:,i)*mask
            enddo
        endif

        ! call sysio_write('grho',m%gradient(:,:,:,1),size(m%gradient(:,:,:,1)))
        ! call sysio_write('gkpa',m%gradient(:,:,:,2),size(m%gradient(:,:,:,2)))

        !transform to x-space
        call param%transform(o_g=qp%g)

        !Regularization in x-space
        ! if(either(oif_approx,.false.,present(oif_approx))) then
        !     call regularize_approximate(fobj,qp)
        ! else
        !    call regularize(self,qp)
        ! endif

        !preconditioner
        call preco%update
        call preco%apply(qp%g,qp%pg)

if(setup%get_bool('IF_STICK_TO_RE',o_default='false')) then
!transform to parameter space
qp%pg(:,:,:,1) = qp%pg(:,:,:,1)*2*m%rho*m%vp
!normaliz by allowed vp range
qp%pg=qp%pg*param%pars(1)%range
!should be in good direction
if(sum(qp%pg*qp%g)<0.) then
    call hud('flip sign of qp%pg (RE) tobe in good direction')
    qp%pg=-qp%pg
endif
call hud('L1norm(qp%pg) before precond: '//num2str(sum(abs(qp%pg))))
!depth precond
call alloc(precond,m%nz)
precond=[(((iz-1)*m%dz)**1.,iz=1,m%nz)]
do iy=1,m%ny; do ix=1,m%nx
    qp%pg(:,ix,iy,1)=qp%pg(:,ix,iy,1)*precond
enddo; enddo
call hud('L1norm(qp%pg) after precond: '//num2str(sum(abs(qp%pg))))
!scale to gradient magnitude
qp%pg = sum(abs(qp%g))/sum(abs(qp%pg)) *qp%pg
endif

        !save some RAM
        call dealloc(m%gradient,m%energy)


        if(allocated(m%correlate)) then !in case WPI
            !correlate(:,:,:,1) = FWI gradient a★Du
            !correlate(:,:,:,2) = one RE -δa★Du
            !correlate(:,:,:,3) = one RE  δu★Da
            !correlate(:,:,:,5) = demig-remig (DR) Adj(Rᴴδu)★Du
            call alloc(FWI,m%nz,m%nx,m%ny); FWI=m%correlate(:,:,:,1)
            call alloc( RE,m%nz,m%nx,m%ny); RE =m%correlate(:,:,:,2)+m%correlate(:,:,:,3)
            call alloc( DR,m%nz,m%nx,m%ny); DR =m%correlate(:,:,:,5)
            where(m%is_freeze_zone)
                FWI=0.; RE=0.; DR=0.
            endwhere

            call hud('angle btw FWI & RE terms=    '//angle_terms(FWI,RE))
            call hud('angle btw FWI & DR terms=    '//angle_terms(FWI,DR))
            call hud('angle btw RE  & DR terms=    '//angle_terms(RE,DR))
            call hud('angle btw RE  & RE+DR terms= '//angle_terms(RE,RE+DR))

            ! !prcondition by z^1
            ! call alloc(tmp,m%nz,o_init=1.)
            ! tmp=[(((iz-1)*m%dz)**1,iz=1,m%nz)]
            ! do iy=1,m%ny; do ix=1,m%nx
            !     RE(:,ix,iy)=RE(:,ix,iy)*tmp(:)
            !     DR(:,ix,iy)=DR(:,ix,iy)*tmp(:)
            ! enddo; enddo

            ! norm1=sqrt(sum(RE*RE))
            ! norm2=sqrt(sum(DR*DR))
            ! call hud('cos(pRE,pDR) = '//num2str( sum(RE*DR)/norm1/norm2 ))

            deallocate(FWI,RE,DR)
            ! deallocate(tmp)

            !deallocate(m%correlate)
        endif

    end subroutine

    subroutine regularize(fobj,qp)
        type(t_fobjective) :: fobj
        type(t_querypoint) :: qp

        if(mpiworld%is_master) call fobj%stack_xnorms(qp%x,qp%g)

    end subroutine

    function angle_terms(vec1,vec2) result(res)
        real,dimension(m%n):: vec1,vec2
        character(:),allocatable :: res

        real :: n1, n2

        n1=norm2(vec1)
        n2=norm2(vec2)

        res = num2str(acos(sum(vec1*vec2)/n1/n2)*180./r_pi,o_format='(f6.2)')// &
        ' ( '//num2str(n1,o_format='(es10.4)')//' , '//num2str(n2,o_format='(es10.4)')//' )'

    end function

end
