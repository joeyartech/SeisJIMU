subroutine modeling_imaging
use mpi
use m_System
use m_Modeling
use m_hilbert
!use m_weighter

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_a1, fld_a2
    type(t_correlate) :: a1_star_a2
    real,dimension(:,:),allocatable :: tmp

    call alloc(correlate_image,m%nz,m%nx,m%ny,ppg%nimag)
    
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%read_from_data2
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)
        call shot%set_var_space2(index(ppg%info,'FDSG')>0)

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project
        
        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer

        call hud('-----------------------')
        call hud('        Imaging        ')
        call hud('-----------------------')

        call ppg%init_field(fld_a1,name='fld_a1',ois_adjoint=.true.)
        call ppg%init_field(fld_a2,name='fld_a2',ois_adjoint=.true.)

        call fld_a1%ignite(o_wavelet=shot%dobs)
        call fld_a2%ignite(o_wavelet=shot%dobs2)
        
        call ppg%init_correlate(a1_star_a2,'a1_star_a2')

        call hud('----  Solving adjoint eqn & xcorrelate  ----')
        !Aᴴa = -Rᴴd
        call ppg%adjoint_hyperbola(fld_a1,fld_a2,a1_star_a2)

        call hud('----  Assemble  ----')
        call ppg%assemble(a1_star_a2)

        call hud('---------------------------------')
        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')


    !collect global correlations
    call mpi_allreduce(mpi_in_place,correlate_image, m%n*ppg%nimag, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    if(mpiworld%is_master) call sysio_write('correlate_image',correlate_image,m%n*ppg%nimag)

    !write correlate
    if(mpiworld%is_master) then
        call a1_star_a2%write

    endif

    ! if(ppg%if_compute_engy) then
    !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    ! endif
        
    call mpiworld%barrier

end subroutine