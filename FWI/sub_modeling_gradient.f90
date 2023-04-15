subroutine modeling_gradient_ip
end subroutine


subroutine modeling_gradient_vp
use mpi
use m_System
use m_Modeling
use m_weighter
use m_Lpnorm
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse

    logical,save :: is_first_in=.true.
    type(t_field) :: fld_u, fld_a
    type(t_correlate) :: u_star_u,a_star_u
    character(:),allocatable :: update_wavelet
    
    fobj%misfit=0.
    
    call alloc(correlate_gradient,m%nz,m%nx,m%ny,ppg%ngrad)
    
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor

        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        call hud('Modeling '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project

        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer

        call ppg%init_field(fld_u,name='fld_u');  call fld_u%ignite
        call ppg%init_correlate(u_star_u,'u_star_u','energy')
                
        !forward modeling
        call ppg%forward(fld_u); call fld_u%acquire

        if(mpiworld%is_master) call shot%write('draw_',shot%dsyn)

        update_wavelet=setup%get_str('UPDATE_WAVELET',o_default='per shot')

        call wei%update
        
        if(update_wavelet/='no') call shot%update_wavelet!(wei%weight) !call gradient_matchfilter_data
        
        !write synthetic data
        call shot%write('dsyn_',shot%dsyn)

        !data weighting
        call wei%update
        
        !Adjoint state method with Lagrangian formulation
        !to compute the gradient (L2 norm example)
        !C = ½║u║² = ½∫ (u-d)² δ(x-xr) dtdx³
        !KᵤC = (u-d)δ(x-xr)
        !L = C + <a|Au-s> ≐ C + <Aᴴa|u>
        !0 = KᵤL = KᵤC + Aᴴa => Aᴴa = -KᵤC = (d-u)δ(x-xr)
        !KₘL = aᴴ KₘA u =: a★Du

        ! !objective function and adjoint source
        ! call fobj%stack_dnorms
        ! shot%dadj=-shot%dadj
        call wei%update!('_4IMAGING')
        fobj%misfit = fobj%misfit &
            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)

        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call kernel_L2sq(shot%dadj)
        call shot%write('dadj_',shot%dadj)

        
        if(mpiworld%is_master) call fobj%print_dnorms('Shotloop-stacked','upto '//shot%sindex)
        ! if(either(oif_gradient,.true.,present(oif_gradient))) then
        
        !adjoint source
        if(update_wavelet/='no') call shot%update_adjsource
        
        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.); call fld_a%ignite
        call ppg%init_correlate(a_star_u,'a_star_u','gradient')
        
        !adjoint modeling
        call ppg%adjoint(fld_a,fld_u,o_a_star_u=a_star_u)
            
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')


    !allreduce objective function value
    call mpi_allreduce(mpi_in_place, [fobj%misfit], 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked misfit '//num2str(fobj%misfit))

    fobj%dnorms=fobj%misfit

    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')

    !scale by shotlist
    call shls%scale(1,o_from_sampled=[fobj%misfit])
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)


    !write correlate
    if(mpiworld%is_master) then
!        call u_star_u%write
        call a_star_u%write
    endif

    !allreduce energy, gradient
    ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale by shotlist
    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)

    if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)
    
    call mpiworld%barrier

end subroutine
