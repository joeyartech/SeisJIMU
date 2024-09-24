subroutine modeling_gradient_ip
use mpi
use m_System
use m_Modeling
use m_separator
use m_weighter
use m_Lpnorm
use m_Envnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_u, fld_a
    type(t_correlate) :: a_star_u


    ! real,dimension(:,:),allocatable :: Wdres
    !real,dimension(:,:,:),allocatable :: Ddt2
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff
    
    !RWI misfit
    fobj%misfit=0.
    fobj%reflection=0.
    fobj%diving=0.

    call alloc(correlate_gradient,m%nz,m%nx,m%ny,ppg%ngrad)
    
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


        call hud('----  Solving A(m)u=s  ----')
        call ppg%init_field(fld_u, name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire; call shot%write('Ru_',shot%dsyn)

        !update wavelet
        if(setup%get_str('UPDATE_WAVELET')/='no') then
!            call wei%update
            call shot%update_wavelet!(wei%weight)
            call matchfilter_apply_to_data(shot%dsyn)

            !write synthetic data
            call shot%write('updated_Ru_',shot%dsyn)
        endif


        if(setup%get_str('JOB')=='forward modeling') cycle


        call sepa%update
        call wei%update!('_4IMAGING')
        call alloc(shot%dadj,shot%nt,shot%nrcv)

        call hud('------------------------')
        call hud('     Build Ip model     ')
        call hud('------------------------')
        
        fobj%misfit = fobj%misfit &
            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*sepa%nearoffset*sepa%reflection,shot%dobs-shot%dsyn, shot%dt)
!            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*sepa%nearoffset, shot%dobs-shot%dsyn, shot%dt)

        call kernel_L2sq(shot%dadj)
        call shot%write('dadj_',shot%dadj)

        call hud('----  Solving A(m)ᴴa = RʳᴴΔdʳ and a★u  ----')
        call ppg%init_correlate(a_star_u,'a_star_u')
        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.); call fld_a%ignite
        call ppg%adjoint(fld_a,fld_u,a_star_u)

        call hud('----  Assemble a★u  ----')
        call ppg%assemble(a_star_u)
        
        call hud('---------------------------------')

    enddo

    call hud('        END LOOP OVER SHOTS        ')
    

    !allreduce RWI misfit values
    call mpi_allreduce(mpi_in_place, [fobj%reflection,fobj%diving,fobj%misfit], 3, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked RWI reflection/diving/total misfit = '//num2str(fobj%reflection)//'/'//num2str(fobj%diving)//'/'//num2str(fobj%misfit))

    fobj%dnorms=fobj%misfit

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
    
    !scale by shotlist
    call shls%scale(1,o_from_sampled=[fobj%misfit])
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)

    !write correlate
    if(mpiworld%is_master) then
        call a_star_u%write
    endif
    

    !allreduce energy, gradient
    ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale by shotlist
    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)

    if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)

    call mpiworld%barrier

    is_first_in=.false.

end subroutine


subroutine modeling_gradient_vp
use mpi
use m_System
use m_Modeling
use m_separator
use m_weighter
use m_Lpnorm
use m_Envnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_u0, fld_u, fld_a0, fld_a
    type(t_correlate) :: a0_star_u0, a_star_u

    character(:),allocatable :: s_dnorm
    real,dimension(:,:),allocatable :: Eobs
    ! character(:),allocatable :: update_wavelet
    
    !RWI misfit
    fobj%misfit=0.
    fobj%reflection=0.
    fobj%diving=0.

    s_dnorm=setup%get_str('DATA_NORM','DNORM',o_default='L2sq')
    if(s_dnorm/='L2sq'.and.s_dnorm/='Envsq') then
        call warn("Sorry, other DATA_NORM hasn't been implemented. Set DNORM to L2sq.")
        s_dnorm='L2sq'
    endif

    call alloc(correlate_gradient,m%nz,m%nx,m%ny,ppg%ngrad)
    
    call hud('===== START LOOP OVER SHOTS =====')
        
    do i=1,shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        if(s_dnorm=='Envsq') then
            call alloc(Eobs,shot%nt,shot%nrcv)
            call hilbert_envelope(shot%dobs,Eobs,shot%nt,shot%nrcv)
            call shot%write('Eobs_',Eobs)
        endif

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project
        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer


        call hud('----  Solving A(m)u=s  ----')
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire;  call shot%write('Ru_',shot%dsyn)

        call hud('----  Solving A(m₀)u₀=s  ----')
        call cb%project(ois_background=.true.)
        call ppg%init
        call ppg%init_field(fld_u0,name='fld_u0');    call fld_u0%ignite
        call ppg%forward(fld_u0)
        call fld_u0%acquire(o_seismo=shot%dsyn_aux);  call shot%write('Ru0_',shot%dsyn_aux)

        !update wavelet
        if(setup%get_str('UPDATE_WAVELET')/='no') then
            call wei%update
            call shot%update_wavelet!(wei%weight)
            call matchfilter_apply_to_data(shot%dsyn)
            call matchfilter_apply_to_data(shot%dsyn_aux)

            !write synthetic data
            call shot%write('updated_Ru_',shot%dsyn)
            call shot%write('updated_Ru0_',shot%dsyn_aux)
        endif

        if(setup%get_str('JOB')=='forward modeling') cycle

        call sepa%update
        call wei%update!('_4IMAGING')
        call alloc(shot%dadj,shot%nt,shot%nrcv)

        call hud('-------------------------')
        call hud('     Update Vp model     ')
        call hud('-------------------------')

        
        !reflection data
        if(s_dnorm=='L2sq') then
            fobj%reflection = fobj%reflection &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*(1.-sepa%nearoffset)*sepa%reflection, shot%dobs-shot%dsyn, shot%dt)
            call kernel_L2sq(shot%dadj)

        elseif(s_dnorm=='Envsq') then
            fobj%reflection = fobj%reflection &
                + Envsq(0.5, wei%weight*(1.-sepa%nearoffset)*sepa%reflection, shot%dsyn, Eobs, shot%dt, shot%nt, shot%nrcv)
            call kernel_Envsq(shot%dadj,shot%nt,shot%nrcv)

        endif

        call shot%write('dadj_rfl_',shot%dadj)
                    
        call hud('----  Solving A(m)ᴴa = RʳᴴΔdʳ and a★u  ----')
        call cb%project
        call ppg%init
        call ppg%init_field(fld_a ,name='fld_a' ,ois_adjoint=.true.); call fld_a%ignite
        call ppg%init_correlate(a_star_u,'a_star_u')
        call ppg%adjoint(fld_a,fld_u,a_star_u)

        !diving waves

        if(s_dnorm=='L2sq') then
            fobj%diving = fobj%diving &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*(1.-sepa%nearoffset)*sepa%diving, shot%dobs-shot%dsyn, shot%dt)
            shot%dadj=-shot%dadj; call kernel_L2sq(shot%dadj,oif_stack=.true.) !diving minus reflection residuals

        elseif(s_dnorm=='Envsq') then
            fobj%diving = fobj%diving &
                + Envsq(0.5, wei%weight*(1.-sepa%nearoffset)*sepa%diving, shot%dsyn, Eobs, shot%dt, shot%nt, shot%nrcv)
            shot%dadj=-shot%dadj; call kernel_Envsq(shot%dadj,shot%nt,shot%nrcv,oif_stack=.true.) !diving minus reflection residuals
            
        endif

        call shot%write('dadj_div-rfl_',shot%dadj)    

        call hud('----  Solving A(m₀)ᴴa₀ = RᵈᴴΔdᵈ-RʳᴴΔdʳ and a₀★u₀  ----')
        call cb%project(ois_background=.true.)
        call ppg%init
        call ppg%init_field(fld_a0,name='fld_a0',ois_adjoint=.true.); call fld_a0%ignite
        call ppg%init_correlate(a0_star_u0,'a0_star_u0')
        call ppg%adjoint(fld_a0,fld_u0,a0_star_u0)

        call hud('----  Assemble a★u+a₀★u₀  ----')
        call ppg%assemble(correlate_add(a_star_u,a0_star_u0))

        !produce total misfit
        fobj%misfit = fobj%reflection+fobj%diving

        call hud('---------------------------------')
            
    enddo

    call hud('        END LOOP OVER SHOTS        ')
    

    !allreduce RWI misfit values
    call mpi_allreduce(mpi_in_place, [fobj%reflection,fobj%diving,fobj%misfit], 3, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked RWI reflection/diving/total misfit = '//num2str(fobj%reflection)//'/'//num2str(fobj%diving)//'/'//num2str(fobj%misfit))

    fobj%dnorms=fobj%misfit

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')
    
    !scale by shotlist
    call shls%scale(1,o_from_sampled=[fobj%misfit])
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)

    !write correlate
    if(mpiworld%is_master) then
        call a_star_u%write
        call a0_star_u0%write
    endif
    
    !allreduce energy, gradient
    ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale by shotlist
    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)

    if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)

    call mpiworld%barrier

    is_first_in=.false.
    
end subroutine
