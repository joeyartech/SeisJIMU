subroutine modeling_gradient
use mpi
use m_System
use m_Modeling
use m_separator
use m_weighter
use m_Lpnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_u0, fld_u, fld_a0, fld_a
    type(t_correlate) :: a0_star_u0, a_star_u

    character(:),allocatable :: s_job
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


        s_job=setup%get_str('JOB',o_default='forward modeling')

        call hud('--------------------------------')
        call hud('        Forward modeling        ')
        call hud('--------------------------------')

        !forward modeling
        !A(m )u = s
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire;  call shot%write('Ru_',shot%dsyn)

        if(index(s_job,'build Ip')<=0) then
            !A(m₀)u₀= s
            call cb%project(ois_background=.true.)
            call ppg%init
            call ppg%init_field(fld_u0,name='fld_u0');    call fld_u0%ignite
            call ppg%forward(fld_u0)
            call fld_u0%acquire(o_seismo=shot%dsyn_aux);  call shot%write('Ru0_',shot%dsyn_aux)
        endif

        !estimate wavelet
        call wei%update
        call shot%update_wavelet(wei%weight)
        call matchfilter_apply_to_data(shot%dsyn)

        !write synthetic data
        call shot%write('updated_Ru',shot%dsyn)
        call shot%write('updated_Ru0',shot%dsyn_aux)

        call hud('---------------------------------')
        if(s_job=='forward modeling') cycle


        call sepa%update
        call wei%update!('_4IMAGING')
        call alloc(shot%dadj,shot%nt,shot%nrcv)

        if(index(s_job,'build Ip')>0) then
            call hud('------------------------')
            call hud('     Build Ip model     ')
            call hud('------------------------')
            
            fobj%misfit = fobj%misfit &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*sepa%nearoffset*sepa%reflection, shot%dobs-shot%dsyn, shot%dt)

            call kernel_L2sq(shot%dadj)
            call shot%write('dadj_',shot%dadj)

            !adjoint modeling
            !A(m)ᴴa = Rʳ(d-u)
            call cb%project
            call ppg%init
            call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.); call fld_a%ignite
            call ppg%init_correlate(a_star_u,'a_star_u','gradient') !a★u
            call ppg%adjoint(fld_a,fld_u,o_a_star_u=a_star_u)
            call hud('---------------------------------')

        endif

        if(index(s_job,'update Vp')>0) then
            call hud('-------------------------')
            call hud('     Update Vp model     ')
            call hud('-------------------------')

            !reflection data
            fobj%reflection = fobj%reflection &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*(1.-sepa%nearoffset)*sepa%reflection, shot%dobs-shot%dsyn, shot%dt)

            call kernel_L2sq(shot%dadj)
            call shot%write('dadj_refl_',shot%dadj)
                        
            !adjoint modeling
            !A(m )ᴴa  = Rʳ(d-u)
            call cb%project
            call ppg%init
            call ppg%init_field(fld_a ,name='fld_a' ,ois_adjoint=.true.); call fld_a%ignite
            call ppg%init_correlate(a_star_u,'a_star_u','gradient') !a★u
            call ppg%adjoint(fld_a,fld_u,o_a_star_u=a_star_u)

            !diving waves
            fobj%diving = fobj%diving &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight*(1.-sepa%nearoffset)*sepa%diving, shot%dobs-shot%dsyn, shot%dt)
            
            shot%dadj=-shot%dadj; call kernel_L2sq(shot%dadj,oif_stack=.true.) !diving minus reflection residuals
            call shot%write('dadj_div-refl_',shot%dadj)    

            !adjoint modeling
            !A(m₀)ᴴa₀ = Rᵈ(d-u)-Rʳ(d-u)
            call cb%project(ois_background=.true.)
            call ppg%init
            call ppg%init_field(fld_a0,name='fld_a0',ois_adjoint=.true.); call fld_a0%ignite
            call ppg%init_correlate(a0_star_u0,'a0_star_u0','gradient') !a₀★u₀                    
            call ppg%adjoint(fld_a0,fld_u0,o_a_star_u=a0_star_u0)

            !produce final misfit
            fobj%misfit = fobj%reflection+fobj%diving

            call hud('---------------------------------')
            
        endif

    enddo

    call hud('        END LOOP OVER SHOTS        ')

    if(s_job=='forward modeling'.or.s_job=='estimate wavelet') then
        call mpiworld%final
        stop
    endif
    

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
    
    if(index(s_job,'update Vp')>0) then
        !gkpa
        correlate_gradient(:,:,:,1) = -ppg%inv_kpa(1:m%nz,1:m%nx,1:m%ny)*(a_star_u%rp_div_sv+a0_star_u0%rp_div_sv)
        !grho0
        correlate_gradient(:,:,:,2) = (a_star_u%rv_grad_sp+a0_star_u0%rv_grad_sp) / m%rho
    endif


    !allreduce energy, gradient
    ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale by shotlist
    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)

    if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)

    call mpiworld%barrier
    
end subroutine