!adjoint eqn from FWI: Aᴴa = Rᴴ(d-u)
!imaging condition I := a★u (~FWI gradient),
!leading to similar phases in Rᴴδu and Rᴴ(d-u) (data residuals)
!therefore, Rᴴδu is a migrated-then-demigrated version of Rᴴ(d-u) (simpler, smoother)
!on the other hand, this choice gives a minus signs in front of the RE terms of the WPI gradient

!L2 norm of image
!I(x) = ∫ a(x,t)u(x,t) dt
!C= ½║I║² = ½∫ I(x)² dx³ = ½∫ a(x,t₁)u(x,t₁) a(x,t₂)u(x,t₂) dt₁dt₂dx³
!δC = ½∫ a(t₁)a(t₂)(u(t₁)δu(t₂)+δu(t₁)u(t₂)) dt₁dt₂dx³
!   =  ∫ a(t₁)a(t₂)u(t₁)δu(t₂) dt₁dt₂dx³
!   =  ∫ I a δu dtdx³
!KᵤC = Ia
!similarly, KₐC = Iu
!
!Adjoint state method with Lagrangian formulation to compute the gradient
!A u=f
!Aᴴa =Rᴴ(d-u)
!L = C +<λ|Au-f> +<μ|Aᴴa-Rᴴ(d-u)>
!  ≐ C +<Aᴴλ|u>  +<Aμ|a> +<Rμ|u>
!0 = KᵤL = KᵤC +Aᴴλ +Rμ => Aᴴλ =-Ia -Rμ => λ=-δa +Adj(Rᴴδu)
!0 = KₐL = KₐC +Aμ      => Aμ  =-Iu     => μ=-δu
!where
!Aᴴδa =Ia, Aᴴ(Adj(Rᴴδu))=Adj(Rᴴδu)
!A δu =Iu
!KₘL = λᴴ KₘA u + μᴴ KₘAᴴ a
!    = -δa★Du +δu★Da +Adj(Rᴴδu)★Du
!
!For PDE:      A u = M ∂ₜ u - D u = f
!    Adjoint:  Aᵀa = M ∂ₜᵀa - Dᵀa = d
!    μᴴ KₘAᵀ a
!= -δuᴴ KₘM ∂ₜᵀa
!≐ -δuᴴ KₘM (M⁻¹Dᵀa)
!=  δuᴴ KₘM (M⁻¹D a)
!=  δu★Da  (more accurate)
!= -a★Dδu  (symmetric to -δa★Du, looks nicer)

subroutine modeling_imaging
use mpi
use m_System
use m_Modeling
use m_propagator_WPI
use m_weighter
use m_Lpnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse

    type(t_field) :: fld_u, fld_a
    ! real,dimension(:,:),allocatable :: Wdres
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff
    
    call alloc(m%image,m%nz,m%nx,m%ny,1)
    
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg_WPI%info,'FDSG')>0)
        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg_WPI%nbndlayer)
        call cb%project
        
        call ppg_WPI%check_discretization
        call ppg_WPI%init
        call ppg_WPI%init_abslayer
        
        call hud('---------------------------------')
        call hud('     Imaging (aka Wavepath)      ')
        call hud('---------------------------------')
        !forward modeling on u
        call ppg_WPI%init_field(fld_u,name='Imag_fld_u');    call fld_u%ignite
        call ppg_WPI%forward(fld_u);                         call fld_u%acquire
        call shot%write('Imag_Ru_',shot%dsyn)

        !weighting on the adjoint source for the image
        call wei%update!('_4IMAGING')
        
        if(index(setup%get_str('MODE',o_default='min I w/ data residual'),'residual')>0) then
            tmp = L2sq(0.5,shot%nrcv*shot%nt,&
                wei%weight, shot%dobs-shot%dsyn, shot%dt)
        else
            tmp = L2sq(0.5,shot%nrcv*shot%nt,&
                wei%weight, shot%dobs          , shot%dt)
        endif
        
        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call kernel_L2sq(shot%dadj)
        call shot%write('Imag_dadj_',shot%dadj)
        
        !adjoint modeling
        call ppg_WPI%init_field(fld_a,name='Imag_fld_a',ois_adjoint=.true.); call fld_a%ignite
        call ppg_WPI%adjoint(fld_a,fld_u,oif_compute_imag=.true.)
        call hud('---------------------------------')

        call cb%project_back

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    !collect
    call mpi_allreduce(mpi_in_place, m%image,  size(m%image), mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)

    call mpiworld%barrier

    if(setup%get_str('JOB')=='imaging') then
        call sysio_write('Image',m%image,size(m%image))
        call mpiworld%final
        stop
    endif

end subroutine


subroutine modeling_gradient
use mpi
use m_System
use m_Modeling
use m_propagator_WPI
use m_weighter
use m_image_weighter
use m_Lpnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    logical,save :: is_first_in=.true.

    character(:),allocatable :: corrs
    type(t_field) :: fld_u, fld_du, fld_Adj_du, fld_a, fld_da
    ! real,dimension(:,:),allocatable :: Wdres
    real,dimension(:,:,:),allocatable :: W2Idt
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff

    !info
    call hud('Entering modeling_gradient')

    !First do imaging
    ! call modeling_imaging
    if(setup%get_file('FILE_IMAGE')/='') then
        call alloc(m%image,m%nz,m%nx,m%ny,1)
        call sysio_read(setup%get_file('FILE_IMAGE'),m%image,size(m%image))
        call sysio_write('loaded_I',m%image,size(m%image))    
        
        !get shot%dt
        call shot%init(shls%yield(1))
        call shot%read_from_data
        call shot%set_var_time
    else
        call hud('Will do imaging')
        call modeling_imaging
    endif
        
    
    !FWI misfit
    fobj%FWI_misfit=0.

    call alloc(m%correlate, m%nz,m%nx,m%ny,5)

    corrs=setup%get_str('COMPUTE_CORRELATIONS','CORRS','RE+DR')
    ! corrs='RE+DR'
    
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg_WPI%info,'FDSG')>0)

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg_WPI%nbndlayer)
        call cb%project
        
        call ppg_WPI%check_discretization
        call ppg_WPI%init
        call ppg_WPI%init_abslayer

        call alloc(cb%corr,     cb%mz,cb%mx,cb%my,5)
        !corr(:,:,:,1) = FWI gradient a★Du
        !corr(:,:,:,2) = one RE -δa★Du
        !corr(:,:,:,3) = one RE  δu★Da
        !corr(:,:,:,4) = 2ndMI δa★Dδu
        !corr(:,:,:,5) = demig-remig (DR) Adj(Rᴴδu)★Du


        call ppg_WPI%init_field(fld_u, name='fld_u');    call fld_u%ignite
        call ppg_WPI%init_field(fld_du,name='fld_du')

        !forward modeling
        !A u = f
        !Aδu = Iu
        call ppg_WPI%forward_scattering(fld_du,fld_u,W2Idt)
        call fld_du%acquire; call shot%write('Rdu_',shot%dsyn)
        call fld_u%acquire;  call shot%write('Ru_', shot%dsyn)
        shot%dsyn_aux=shot%dsyn

        !compute FWI data misfit C_data=║Δd║²
        call wei%update
        fobj%FWI_misfit = fobj%FWI_misfit &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)

        if(index(setup%get_str('MODE',o_default='min I w/ data residual'),'residual')>0) then
            tmp = fobj%FWI_misfit 
        else
            tmp = L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs          , shot%dt)
        endif
        
        !adjoint modeling
        !Aᴴ a = Rᴴ(d-u)
        !Aᴴδa = Ia
        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call kernel_L2sq(shot%dadj)

        call ppg_WPI%init_field(fld_a, name='fld_a' ,ois_adjoint=.true.); call fld_a%ignite
        call ppg_WPI%init_field(fld_da,name='fld_da',ois_adjoint=.true.)
        
        call ppg_WPI%init_field(fld_Adj_du,name='fld_Adj_du',ois_adjoint=.true.)
        call fld_du%acquire
        shot%dadj=shot%dsyn*wei%weight*wei%weight
        if(mpiworld%is_master) call shot%write('Wei_Rdu_',shot%dadj)
        call fld_Adj_du%ignite
                
        call ppg_WPI%adjoint_WPI(fld_da,fld_a,fld_Adj_du,fld_du,fld_u,W2Idt,corrs)
        !if(index(corrs,'RE')>0) cb%corr(:,:,:,2)=-cb%corr(:,:,:,2)
        cb%corr(:,:,:,2)=-cb%corr(:,:,:,2)

        call cb%project_back

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')


    !collect FWI misfit values
    !no need to collect WPI misfit values
    call mpi_allreduce(mpi_in_place, fobj%FWI_misfit, 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked FWI_misfit '//num2str(fobj%FWI_misfit))

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
    
    !collect global correlations
    call mpi_allreduce(mpi_in_place, m%correlate,  m%n*5, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale
    call shls%scale(1,o_from_sampled=[fobj%FWI_misfit])
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)
    call shls%scale(m%n*5,o_from_sampled=m%correlate)

    if(mpiworld%is_master) then
            call sysio_write( 'a_star_Du',     m%correlate(:,:,:,1),m%n)
        if(index(corrs,'RE')>0) then
            call sysio_write('-da_star_Du',    m%correlate(:,:,:,2),m%n)
            call sysio_write( 'a_star_Ddu',    m%correlate(:,:,:,3),m%n)
            call sysio_write('RE',             m%correlate(:,:,:,2)+m%correlate(:,:,:,3),m%n)
        endif
        if(index(corrs,'DR')>0) then
            call sysio_write('Adj(du)_star_Du',m%correlate(:,:,:,5),m%n)
        endif
    endif

    call alloc(m%gradient,m%nz,m%nx,m%ny,2)
    ! if(index(corrs,'FWI')>0) then
    !     m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,1)
    ! endif
    ! if(index(corrs,'RE')>0) then
    !     m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,2)+m%correlate(:,:,:,3)
    ! endif
    ! if(index(corrs,'DR')>0) then
    !     m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,5)
    ! endif
    
    !RE
    m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,2)+m%correlate(:,:,:,3)

    !DR
    if(index(setup%get_str('MODE',o_default='min I w/ data residual'),'residual')>0) then
        m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,5)
    endif


    !deallocate(m%correlate)
    
    !required, otherwise cb%project_back will project back cb%corr
    !in the sequential calling of modeling_imaging
    deallocate(cb%corr)
        
    ! if(ppg_WPI%if_compute_engy) then
    !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg_WPI%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    ! endif
    
    
    call mpiworld%barrier

end subroutine
