!adjoint eqn from FWI: Aᴴa = Rᴴ(d-u)
!imaging condition I := a★u (~FWI gradient),
!which leads to Rδu w/ same reflection phases as in d
!but minus signs in front of the RE terms of the WPI gradient

!L2 norm of image
!I(x) = ∫ a(x,t)u(x,t) dt
!C= ½║I║² = ½∫ I(x)² dx³ = ½∫ a(x,t₁)u(x,t₁) a(x,t₂)u(x,t₂) dt₁dt₂dx³
!δC = ½∫ a(t₁)a(t₂)(u(t₁)δu(t₂)+δu(t₁)u(t₂)) dt₁dt₂dx³
!   =  ∫ a(t₁)a(t₂)u(t₁)δu(t₂) dt₁dt₂dx³
!   =  ∫ I a δu dtdx³
!∇ᵤC = Ia
!similarly, ∇ₐC = Iu
!
!Adjoint state method with Lagrangian formulation to compute the gradient
!A u=f
!Aᴴa =Rᴴ(d-u)
!L = C +<λ|Au-f> +<μ|Aᴴa-Rᴴ(d-u)>
!  ≐ C +<Aᴴλ|u>  +<Aμ|a> +<Rμ|u>
!0 = ∇ᵤL = ∇ᵤC +Aᴴλ +Rμ => Aᴴλ =-Ia -Rμ => λ=-δa +Adj(Rᴴδu)
!0 = ∇ₐL = ∇ₐC +Aμ      => Aμ  =-Iu     => μ=-δu
!where
!Aᴴδa =Ia, Aᴴ(Adj(Rᴴδu))=Adj(Rᴴδu)
!A δu =Iu
!∇ₘL = λᴴ ∇ₘA u + μᴴ ∇ₘAᴴ a
!    = -δa★Du +δu★Da +Adj(Rᴴδu)★Du
!
!For PDE:      A u = M ∂ₜ u - D u = f
!    Adjoint:  Aᵀa = M ∂ₜᵀa - Dᵀa = d
!    μᴴ ∇ₘAᵀ a
!= -δuᴴ ∇ₘM ∂ₜᵀa
!≐ -δuᴴ ∇ₘM (M⁻¹Dᵀa)
!=  δuᴴ ∇ₘM (M⁻¹D a)
!=  δu★Da  (more accurate)
!= -a★Dδu  (symmetric to -δa★Du, looks nicer)

subroutine modeling_imaging
use mpi
use m_System
use m_Modeling
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
        call shot%set_var_space(index(ppg%info,'FDSG')>0)
        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project
        
        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer
        
        call hud('---------------------------------')
        call hud('     Imaging (aka Wavepath)      ')
        call hud('---------------------------------')
        !forward modeling on u
        call ppg%init_field(fld_u,name='Imag_fld_u');    call fld_u%ignite
        call ppg%forward(fld_u);                         call fld_u%acquire
        call shot%write('Imag_Ru_',shot%dsyn)

        !weighting on the adjoint source for the image
        call wei%update!('_4IMAGING')
        
        tmp = L2sq(0.5,shot%nrcv*shot%nt,&
            wei%weight, shot%dobs-shot%dsyn, shot%dt)
        
        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call nabla_L2sq(shot%dadj)
        call shot%write('Imag_dadj_',shot%dadj)
        
        !adjoint modeling
        call ppg%init_field(fld_a,name='Imag_fld_a',ois_adjoint=.true.); call fld_a%ignite
        call ppg%adjoint(fld_a,fld_u,oif_compute_imag=.true.)
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


subroutine modeling_gradient_slow(is_fitting_data)
use mpi
use m_System
use m_Modeling
use m_weighter
use m_image_weighter
use m_Lpnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    logical :: is_fitting_data

    logical,save :: is_first_in=.true.

    character(:),allocatable :: corrs
    type(t_field) :: fld_u, fld_du, fld_Adj_du, fld_a, fld_da
    ! real,dimension(:,:),allocatable :: Wdres
    ! real,dimension(:,:,:),allocatable :: WImag
    real,dimension(:,:,:),allocatable :: W2Idt
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff

    !info
    call hud('Entering modeling_gradient_slow')

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

    !weights for the image
    call iwei%update

    !L1 norm of image
    fobj%dnorms(1) = L1  (1. ,m%n,iwei%weight,m%image,m%cell_volume)

    fobj%dnorms(2) = L2sq(0.5,m%n,iwei%weight,m%image,m%cell_volume)
    call alloc(W2Idt,m%nz,m%nx,m%ny)
    call nabla_L2sq(W2Idt)

    !times dt in the RHS for propagation
    W2Idt=W2Idt*shot%dt
    call sysio_write('W2Idt',W2Idt,size(W2Idt))
    
    
    !FWI misfit
    fobj%FWI_misfit=0.

    call alloc(m%correlate, m%nz,m%nx,m%ny,5)

    corrs=setup%get_str('CORRELATIONS','CORRS','RE+DR')
    
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project
        
        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer

        call alloc(cb%corr,     cb%mz,cb%mx,cb%my,5)

    ! if(.false.) then
        call hud('-------- FWI SK (a_star_Du) --------')
        !forward modeling on u
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite       
        call ppg%forward(fld_u);                    call fld_u%acquire

        call shot%write('Ru_',shot%dsyn)

        call wei%update!('_4FWI')

        !compute FWI data misfit C_data=║Δd║²
        !adjoint modeling Aᴴa = Rᴴ(d-u)
        !FWI gradient = a★Du
        fobj%FWI_misfit = fobj%FWI_misfit &
            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
            
        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call nabla_L2sq(shot%dadj)
        call shot%write('FWI_dadj_',shot%dadj)

        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)
        call fld_a%ignite
        call ppg%adjoint(fld_a,fld_u,oif_compute_grad=.true.)

        cb%corr(:,:,:,1)=cb%grad(:,:,:,2) !=gkpa, propto gvp under Vp-Rho
        call hud('---------------------------------')
    ! endif

        if(index(corrs,'RE')>0) then

            call hud('--- extd Rabbit Ear (-da_star_Du) ---')
            !re-model u
            call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
            call ppg%forward(fld_u); call fld_u%acquire
            
            !adjoint eqn Aᴴa = Rᴴ(d-u), Aδa = Ia
            call wei%update
            tmp = L2sq(0.5,shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
            call nabla_L2sq(shot%dadj)
            
            !adjoint modeling
            !grad = -δa★Du
            call ppg%init_field(fld_a, name='fld_a', ois_adjoint=.true.); call fld_a%ignite
            call ppg%init_field(fld_da,name='fld_da',ois_adjoint=.true.)
            call ppg%adjoint_da_star_Du(fld_da,fld_u,fld_a,W2Idt,oif_compute_grad=.true.)

            cb%corr(:,:,:,2)=-cb%grad(:,:,:,2) !=gkpa, propto gvp under Vp-Rho
            call hud('---------------------------------')

            !pause

            call hud('--- extd Rabbit Ear (du_star_Da) ---')
            call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
            call ppg%init_field(fld_du,name='fld_du')
            
            !forward scattering
            !Aδu = Iu
            call ppg%forward_scattering(fld_du,fld_u,W2Idt)
            call fld_du%acquire; call shot%write('Rdu_',shot%dsyn)
            call fld_u%acquire !shot%dsyn = Ru
            
            !adjoint eqn Aᴴa = Rᴴ(d-u)
            call wei%update
            tmp = L2sq(0.5,shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
            call nabla_L2sq(shot%dadj)
            
            !adjoint modeling
            !grad = δu★Da
            call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.); call fld_a%ignite
            call ppg%adjoint_du_star_Da(fld_a,fld_du,fld_u,W2Idt,oif_compute_grad=.true.) !,oif_compute_imag=.true.)

            cb%corr(:,:,:,3)=cb%grad(:,:,:,2) !=gkpa, propto gvp under Vp-Rho
            call hud('---------------------------------')

        endif

        if(index(corrs,'2ndMI')>0) then

            call hud('-------- 2ndMI da_star_Ddu --- (just for curiosity)')
            call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
            call ppg%init_field(fld_du,name='fld_du')
            call ppg%forward_scattering(fld_du,fld_u,W2Idt)
            call fld_u%acquire !shot%dsyn = Ru

            !adjoint eqn Aᴴa = Rᴴ(d-u), Aδa = Ia
            call wei%update
            tmp = L2sq(0.5,shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
            call nabla_L2sq(shot%dadj)
            
            !adjoint modeling
            !grad = δu★Dδa
            call ppg%init_field(fld_a, name='fld_a', ois_adjoint=.true.); call fld_a%ignite
            call ppg%init_field(fld_da,name='fld_da',ois_adjoint=.true.)
            call ppg%adjoint_da_star_Ddu(fld_da,fld_du,fld_a,fld_u,W2Idt,oif_compute_grad=.true.)

            cb%corr(:,:,:,4)=cb%grad(:,:,:,2) !=gkpa, propto gvp under Vp-Rho
            call hud('---------------------------------')

        endif


        if(index(corrs,'DR')>0) then

            call hud('--- demig-remig (Adj(du)_star_Du) ---')
            call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
            call ppg%init_field(fld_du,name='fld_du')
            call ppg%forward_scattering(fld_du,fld_u,W2Idt); call fld_du%acquire !shot%dsyn = Rδu
            
            !adjoint eqn Aᴴa = RᴴRδu
            call wei%update
            shot%dadj=shot%dsyn*wei%weight*wei%weight
            if(mpiworld%is_master) call shot%write('Wei_Rdu_',shot%dadj)
                        
            call ppg%init_field(fld_Adj_du,name='fld_Adj_du',ois_adjoint=.true.)
            call fld_Adj_du%ignite
            
            !adjoint modeling
            !grad = Adj(δu)★Du
            call ppg%adjoint(fld_Adj_du,fld_u,oif_compute_grad=.true.)
            
            cb%corr(:,:,:,5)=cb%grad(:,:,:,2) !=gkpa, propto gvp under Vp-Rho
            call hud('---------------------------------')

        endif
        
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
            call sysio_write('a_star_Du',      m%correlate(:,:,:,1),m%n)
        if(index(corrs,'RE')>0) then
            call sysio_write('-da_star_Du',    m%correlate(:,:,:,2),m%n)
            call sysio_write(' du_star_Da',    m%correlate(:,:,:,3),m%n)
            call sysio_write('RE',             m%correlate(:,:,:,2)+m%correlate(:,:,:,3),m%n)
        endif
        if(index(corrs,'2ndMI')>0) then
            call sysio_write('du_star_Dda',    m%correlate(:,:,:,4),m%n)
        endif
        if(index(corrs,'DR')>0) then
            call sysio_write('Adj(du)_star_Du',m%correlate(:,:,:,5),m%n)
        endif
    endif

    m%gradient=0.
    if(index(corrs,'FWI')>0) then
        m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,1)
    endif
    if(index(corrs,'RE')>0) then
        m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,2)+m%correlate(:,:,:,3)
    endif
    if(index(corrs,'DR')>0) then
        m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,5)
    endif

    !check if fitting the t-x domain data
    is_fitting_data = sum(m%correlate(:,:,:,1)*m%gradient(:,:,:,2))>0.

    !deallocate(m%correlate)
    
    !required, otherwise cb%project_back will project_back cb%corr
    !in the sequential calling of modeling_imaging
    deallocate(cb%corr)
    
    if(ppg%if_compute_engy) then
        call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    endif
    
    
    call mpiworld%barrier

end subroutine


subroutine modeling_gradient_costRAM(is_fitting_data)
use mpi
use m_System
use m_Modeling
use m_weighter
use m_image_weighter
use m_Lpnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    logical :: is_fitting_data

    logical,save :: is_first_in=.true.

    character(:),allocatable :: corrs
    type(t_field) :: fld_u, fld_du, fld_Adj_du, fld_a, fld_da
    ! real,dimension(:,:),allocatable :: Wdres
    real,dimension(:,:,:),allocatable :: W2Idt
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff

    !info
    call hud('Entering modeling_gradient_costRAM')

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

    !weights for the image
    call iwei%update

    !L1 norm of image
    fobj%dnorms(1) = L1  (1. ,m%n,iwei%weight,m%image,m%cell_volume)

    fobj%dnorms(2) = L2sq(0.5,m%n,iwei%weight,m%image,m%cell_volume)
    call alloc(W2Idt,m%nz,m%nx,m%ny)
    call nabla_L2sq(W2Idt)
    
    !times dt in the RHS for propagation
    W2Idt=W2Idt*shot%dt
    call sysio_write('W2Idt',W2Idt,size(W2Idt))
        
    
    !FWI misfit
    fobj%FWI_misfit=0.

    call alloc(m%correlate, m%nz,m%nx,m%ny,5)

    corrs=setup%get_str('CORRELATIONS','CORRS','RE+DR')
    
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project
        
        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer

        call alloc(cb%corr,     cb%mz,cb%mx,cb%my,5)
        !corr(:,:,:,1) = FWI gradient a★Du
        !corr(:,:,:,2) = one RE -δa★Du
        !corr(:,:,:,3) = one RE  δu★Da
        !corr(:,:,:,4) = 2ndMI δa★Dδu
        !corr(:,:,:,5) = demig-remig (DR) Adj(Rᴴδu)★Du


        call ppg%init_field(fld_u, name='fld_u');    call fld_u%ignite
        call ppg%init_field(fld_du,name='fld_du')

        !forward modeling
        !A u = f
        !Aδu = Iu
        call ppg%forward_scattering(fld_du,fld_u,W2Idt)
        call fld_du%acquire; call shot%write('Rdu_',shot%dsyn)
        call fld_u%acquire;  call shot%write('Ru_', shot%dsyn)

        !compute FWI data misfit C_data=║Δd║²
        call wei%update
        fobj%FWI_misfit = fobj%FWI_misfit &
            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)
        
        !adjoint modeling
        !Aᴴ a = Rᴴ(d-u)
        !Aᴴδa = Ia
        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call nabla_L2sq(shot%dadj)
        call shot%write('FWI_dadj_',shot%dadj)

        call ppg%init_field(fld_a, name='fld_a' ,ois_adjoint=.true.); call fld_a%ignite
        call ppg%init_field(fld_da,name='fld_da',ois_adjoint=.true.)
        
        call ppg%init_field(fld_Adj_du,name='fld_Adj_du',ois_adjoint=.true.)
        call fld_du%acquire
        shot%dadj=shot%dsyn*wei%weight*wei%weight
        if(mpiworld%is_master) call shot%write('Wei_Rdu_',shot%dadj)
        call fld_Adj_du%ignite
                
        call ppg%adjoint_WPI(fld_da,fld_a,fld_Adj_du,fld_du,fld_u,W2Idt,corrs)
        if(index(corrs,'RE')>0) cb%corr(:,:,:,2)=-cb%corr(:,:,:,2)

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
    if(index(corrs,'FWI')>0) then
        m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,1)
    endif
    if(index(corrs,'RE')>0) then
        m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,2)+m%correlate(:,:,:,3)
    endif
    if(index(corrs,'DR')>0) then
        m%gradient(:,:,:,2)=m%gradient(:,:,:,2) +m%correlate(:,:,:,5)
    endif

    !check if fitting the t-x domain data
    is_fitting_data = sum(m%correlate(:,:,:,1)*m%gradient(:,:,:,2))>0.

    !deallocate(m%correlate)
    
    !required, otherwise cb%project_back will project_back cb%corr
    !in the sequential calling of modeling_imaging
    deallocate(cb%corr)
        
    ! if(ppg%if_compute_engy) then
    !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    ! endif
    
    
    call mpiworld%barrier

end subroutine
