module m_propagator_WPI
use m_System
use m_hicks, only : hicks_r
use m_resampler
use m_model
use m_shot
use m_computebox
use m_field
use m_cpml

    private

    !FD coef
    real,dimension(2),parameter :: coef = [9./8.,-1./24.] !Fornberg, 1988, Generation of Finite-Difference Formulas on Arbitrary Spaced Grids.
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z

    type,extends(t_propagator),public :: t_propagator_WPI

        contains
        procedure :: forward_scattering
        
        procedure :: adjoint_da_star_Du
        procedure :: adjoint_du_star_Da
        procedure :: adjoint_da_star_Ddu
        procedure :: adjoint_WPI
        
        procedure :: inject_stresses_scattering
        
    end type

    type(t_propagator),public :: ppg

    logical :: if_hicks
    integer :: irdt
    real :: rdt

    ! real,dimension(:,:,:),allocatable :: sf_p_save

    contains
    
    ! => FD eqn:
    !   [     ] [ vz^n  ]   [ 0     0    ∂zᵇ               ][ vz^n+1]   
    !   | ∂ₜᶠ | | vx^n  |   | 0     0    ∂ₓᵇ       0       || vx^n+1|   
    ! M |     |⨂|  p^n+1| = |∂zᶠ   ∂ₓᶠ    0                ||  p^n+½| +f  
    !   |     | |vzᵃ^n  |   |                0     0    ∂zᵇ||vzᵃ^n+1|   
    !   | ∂ₜᵇ | |vxᵃ^n  |   |       0        0     0    ∂ₓᵇ||vxᵃ^n+1|   
    !   [     ] [ pᵃ^n+1]   [               ∂zᶠ   ∂ₓᶠ    0 ][ pᵃ^n+½]   
    ! where
    ! f=[fz fx fp ... -dz -dx -dp]ᵀ

    ! Time marching:
    ! [ vz^n+1 ]   [ vz^n  ]      [∂zᵇ p^n+½                ]
    ! | vx^n+1 |   | vx^n  |      |∂ₓᵇ p^n+½                |
    ! |  p^n+1½| = |  p^n+½| + M⁻¹[∂zᶠ vz^n+1 + ∂ₓᶠ vx^n+1  |dt +M⁻¹f*dt
    ! |vzᵃ^n+1 |   |vzᵃ^n  |      [∂zᵇ pᵃ^n+½               |
    ! |vxᵃ^n+1 |   |vxᵃ^n  |      |∂ₓᵇ pᵃ^n+½               |
    ! [ pᵃ^n+1½]   [ pᵃ^n+½]      [∂zᶠ vzᵃ^n+1 + ∂ₓᶠ vxᵃ^n+1]

    ! Reverse time marching:
    ! [ vz^n  ]   [ vz^n+1 ]      [∂zᵇ p^n+½                ]
    ! | vx^n  |   | vx^n+1 |      |∂ₓᵇ p^n+½                |
    ! |  p^n+½| = |  p^n+1½| - M⁻¹[∂zᶠ vz^n+1 + ∂ₓᶠ vx^n+1  |dt -M⁻¹f*dt
    ! |vzᵃ^n  |   |vzᵃ^n+1 |      [∂zᵇ pᵃ^n+½               |
    ! |vxᵃ^n  |   |vxᵃ^n+1 |      |∂ₓᵇ pᵃ^n+½               |
    ! [ pᵃ^n+½]   [ pᵃ^n+1½]      [∂zᶠ vzᵃ^n+1 + ∂ₓᶠ vxᵃ^n+1]

    ! <= FD eqn:
    !   [     ] [ vzλ^n  ]   [ 0     0    ∂zᵇ               ][ vzλ^n+1]   
    !   | ∂ₜᵇ | | vxλ^n  |   | 0     0    ∂ₓᵇ       0       || vxλ^n+1|   
    ! M |     |⨂|  pλ^n+1| = |∂zᶠ   ∂ₓᶠ    0                ||  pλ^n+½| -g
    !   |     | | vzμ^n  |   |                0     0    ∂zᵇ|| vzμ^n+1|   
    !   | ∂ₜᶠ | | vxμ^n  |   |       0        0     0    ∂ₓᵇ|| vxμ^n+1|   
    !   [     ] [  pμ^n+1]   [               ∂zᶠ   ∂ₓᶠ    0 ][  pμ^n+½]   

    ! Time marching:
    ! [ vzλ^n+1 ]   [ vzλ^n  ]      [∂zᵇ  pλ^n+½              ]
    ! | vxλ^n+1 |   | vxλ^n  |      |∂ₓᵇ  pλ^n+½              |
    ! |  pλ^n+1½| = |  pλ^n+½| + M⁻¹[∂zᶠ vzλ^n+1 + ∂ₓᶠ vxλ^n+1|dt -M⁻¹g*dt
    ! | vzμ^n+1 |   | vzμ^n  |      [∂zᵇ  pμ^n+½              |
    ! | vxμ^n+1 |   | vxμ^n  |      |∂ₓᵇ  pμ^n+½              |
    ! [  pμ^n+1½]   [  pμ^n+½]      [∂zᶠ vzμ^n+1 + ∂ₓᶠ vxμ^n+1]

    ! Reverse time marching:
    ! [ vzλ^n  ]   [ vzλ^n+1 ]      [∂zᵇ  pλ^n+½              ]
    ! | vxλ^n  |   | vxλ^n+1 |      |∂ₓᵇ  pλ^n+½              |
    ! |  pλ^n+½| = |  pλ^n+1½| - M⁻¹[∂zᶠ vzλ^n+1 + ∂ₓᶠ vxλ^n+1|dt +M⁻¹g*dt
    ! | vzμ^n  |   | vzμ^n+1 |      [∂zᵇ  pμ^n+½              |
    ! | vxμ^n  |   | vxμ^n+1 |      |∂ₓᵇ  pμ^n+½              |
    ! [  pμ^n+½]   [  pμ^n+1½]      [∂zᶠ vzμ^n+1 + ∂ₓᶠ vxμ^n+1]

    subroutine forward_scattering(self,fld_du,fld_u,W2Idt)
        class(t_propagator) :: self
        type(t_field) :: fld_du,fld_u
        real,dimension(cb%mz,cb%mx,cb%my) :: W2Idt

        real,parameter :: time_dir=1. !time direction

        !seismo
        call alloc(fld_u%seismo, shot%nrcv,self%nt)
        call alloc(fld_du%seismo,shot%nrcv,self%nt)
        
        tt1=0.; tt2=0.; tt3=0.; tt4=0.; tt5=0.; tt6=0.
        
        ift=1; ilt=self%nt

        do it=ift,ilt
            if(mod(it,100)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_du%check_value
                 call fld_u%check_value
            endif

            !do forward time stepping (step# conforms with backward & adjoint time stepping)
            ! !step 1: add forces to v^it
            ! call cpu_time(tic)
            ! call self%inject_velocities(fld_u,time_dir,it)
            ! call self%inject_velocities_scattering(fld_du,fld_u,W2Idt,time_dir)
            ! call cpu_time(toc)
            ! tt1=tt1+toc-tic

            !step 2: from v^it to v^it+1 by differences of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities(fld_u, time_dir,it)
            call self%update_velocities(fld_du,time_dir,it)
            call cpu_time(toc)
            tt2=tt2+toc-tic

            !step 3: add pressure to s^it+0.5
            call cpu_time(tic)
            call self%inject_stresses(fld_u,time_dir,it)
            call self%inject_stresses_scattering(fld_du,fld_u,W2Idt,time_dir)
            call cpu_time(toc)
            tt3=tt3+toc-tic

            !step 4: from s^it+0.5 to s^it+1.5 by differences of v^it+1
            call cpu_time(tic)
            call self%update_stresses(fld_u, time_dir,it)
            call self%update_stresses(fld_du,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !step 5: sample v^it+1 or s^it+1.5 at receivers
            call cpu_time(tic)
            call self%extract(fld_u,it)
            call self%extract(fld_du,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !snapshot
            ! call fld_u%write(it,o_suffix='_for2')
            call fld_du%write(it)

            !step 6: save v^it+1 in boundary layers
            ! if(fld_u%if_will_reconstruct) then
                call cpu_time(tic)
                call  fld_u%boundary_transport('save',it)
                call fld_du%boundary_transport('save',it)
                call cpu_time(toc)
                tt6=tt6+toc-tic
            ! endif

        enddo

        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to add source velocities',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update velocities    ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to add source stresses  ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update stresses      ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract field        ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to save boundary        ',tt6/mpiworld%max_threads
        endif

    end subroutine

    subroutine adjoint_du_star_Da(self,fld_a,fld_du,fld_u,W2Idt,oif_compute_imag,oif_compute_grad)
        class(t_propagator) :: self
        type(t_field) :: fld_du,fld_a,fld_u
        real,dimension(cb%mz,cb%mx,cb%my) :: W2Idt
        logical,optional :: oif_compute_imag,oif_compute_grad
        
        real,parameter :: time_dir=-1. !time direction
        logical :: if_compute_imag,if_compute_grad

        ! if_record_adjseismo =either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))
        if_compute_imag     =either(oif_compute_imag,    .false.,present(oif_compute_imag))
        if_compute_grad     =either(oif_compute_grad,    .false.,present(oif_compute_grad))
        ! self%if_compute_engy=self%if_compute_engy.and.(if_compute_imag.or.if_compute_grad)
        
        ! if(if_record_adjseismo)  call alloc(fld_a%seismo,1,self%nt)
        if(if_compute_imag)      call alloc(cb%imag,cb%mz,cb%mx,cb%my,self%nimag)
        if(if_compute_grad)      call alloc(cb%grad,cb%mz,cb%mx,cb%my,self%ngrad)
        ! if(self%if_compute_engy) call alloc(cb%engy,cb%mz,cb%mx,cb%my,self%nengy)

        !reinitialize absorbing boundary for incident wavefield reconstruction
        ! if(present(o_sf)) then
                                          fld_du%dvz_dz=0.
                                           fld_u%dvz_dz=0.
                                          fld_du%dvx_dx=0.
                                           fld_u%dvx_dx=0.
            ! if(allocated(fld_du%dvy_dy))  fld_du%dvy_dy=0.
            ! if(allocated( fld_u%dvy_dy))   fld_u%dvy_dy=0.
                                          fld_du%dp_dz=0.
                                           fld_u%dp_dz=0.
                                          fld_du%dp_dx=0.
                                           fld_u%dp_dx=0.
            ! if(allocated(fld_du%dp_dy))   fld_du%dp_dy=0.
            ! if(allocated( fld_u%dp_dy))    fld_u%dp_dy=0.
        ! endif
                    
        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.
        
        ift=1; ilt=self%nt

        ! call alloc(sf_p_save,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        
        do it=ilt,ift,int(time_dir)
            if(mod(it,100)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_a%check_value
                call fld_du%check_value
                call fld_u%check_value
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then

                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_du%boundary_transport('load',it)
                call  fld_u%boundary_transport('load',it)
                call cpu_time(toc)
                tt1=tt1+toc-tic
                
                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%update_stresses(fld_du,time_dir,it)
                call self%update_stresses( fld_u,time_dir,it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_stresses_scattering(fld_du,fld_u,W2Idt,time_dir)
                call self%inject_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_stresses(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_stresses(fld_a,time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: sfield%s_dt^it+0.5 \dot rfield%s^it+0.5
            !use sfield%v^it+1 to compute sfield%s_dt^it+0.5, as backward step 4
            if(if_compute_grad.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call gradient_moduli(fld_du,fld_a,it,cb%grad(:,:,:,2))
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif

            ! if(if_compute_imag.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call imaging(fld_du,fld_a,it,cb%imag)
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif

            ! !energy term of sfield
            ! if(self%if_compute_engy.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call energy(fld_u,it,cb%engy)
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
                
            !========================================================!

            ! if(present(o_sf)) then
                !backward step 2: v^it+1 -> v^it by FD of s^it+0.5
                call cpu_time(tic)
                call self%update_velocities(fld_du,time_dir,it)
                call self%update_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt7=tt7+toc-tic

                !backward step 1: rm forces from v^it
                call cpu_time(tic)
                ! call self%inject_velocities_ext(fld_du,fld_u,imag)
                call self%inject_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt8=tt8+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 3: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_velocities(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities(fld_a,time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            ! !adjoint step 1: sample v^it or s^it+0.5 at source position
            ! if(if_record_adjseismo) then
            !     call cpu_time(tic)
            !     call self%extract(fld_a,it)
            !     call cpu_time(toc)
            !     tt11=tt11+toc-tic
            ! endif
            
            ! !grho: sfield%v_dt^it \dot rfield%v^it
            ! !use sfield%s^it+0.5 to compute sfield%v_dt^it, as backward step 2
            ! if(if_compute_grad.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call gradient_density(o_sf,rf,it,cb%grad(:,:,:,1))
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
            
            !snapshot
            ! call fld_a%write(it,o_suffix='_back2')
            call fld_du%write(it,o_suffix='_rev')
            ! call fld_u%write(it,o_suffix='_back2')
            if(if_compute_imag) then
                call fld_a%write_ext(it,'imag_a_star_du' ,cb%imag,size(cb%imag))
            endif
            if(if_compute_grad) then
                call fld_a%write_ext(it,'grad_a_star_Ddu' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
            endif
            ! if(self%if_compute_engy) then
            !     call fld_a%write_ext(it,'engy',cb%engy(:,:,:,1),size(cb%engy(:,:,:,1)))
            ! endif

        enddo
        
        !postprocess gradient
        if(if_compute_grad) call gradient_postprocess
        
        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to load boundary            ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update stresses          ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source stresses       ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update velocities        ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source velocities     ',tt8/mpiworld%max_threads
            write(*,*) 'Elapsed time ----------------------------'
            write(*,*) 'Elapsed time to add adjsource stresses   ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj stresses      ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to add adjsource velocities ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj velocities    ',tt10/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract&write fields     ',tt11/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate                ',tt6/mpiworld%max_threads

        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_rfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_rfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')
        call hud('ximage < snap_*  n1='//num2str(cb%mz)//' perc=99')
        call hud('xmovie < snap_*  n1='//num2str(cb%mz)//' n2='//num2str(cb%mx)//' clip=?e-?? loop=2 title=%g')
        
    end subroutine

    subroutine adjoint_da_star_Du(self,fld_da,fld_u,fld_a,W2Idt,oif_record_adjseismo,oif_compute_grad)
        class(t_propagator) :: self
        type(t_field) :: fld_da,fld_u,fld_a
        real,dimension(cb%mz,cb%mx,cb%my) :: W2Idt
        logical,optional :: oif_record_adjseismo, oif_compute_grad
        
        real,parameter :: time_dir=-1. !time direction
        logical :: if_record_adjseismo, if_compute_imag, if_compute_grad

        if_record_adjseismo =either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))
        ! if_compute_imag     =either(oif_compute_imag,    .false.,present(oif_compute_imag))
        if_compute_grad     =either(oif_compute_grad,    .false.,present(oif_compute_grad))
        ! self%if_compute_engy=self%if_compute_engy.and.(if_compute_imag.or.if_compute_grad)
        
        if(if_record_adjseismo)  call alloc(fld_da%seismo,1,self%nt)
        ! if(if_compute_imag)      call alloc(cb%imag,cb%mz,cb%mx,cb%my,self%nimag)
        if(if_compute_grad)      call alloc(cb%grad,cb%mz,cb%mx,cb%my,self%ngrad)
        ! if(self%if_compute_engy) call alloc(cb%engy,cb%mz,cb%mx,cb%my,self%nengy)

        !reinitialize absorbing boundary for incident wavefield reconstruction
        ! if(present(o_sf)) then
                                         fld_u%dvz_dz=0.
                                         fld_u%dvx_dx=0.
            if(allocated(fld_u%dvy_dy))  fld_u%dvy_dy=0.
                                         fld_u%dp_dz=0.
                                         fld_u%dp_dx=0.
            if(allocated(fld_u%dp_dy))   fld_u%dp_dy=0.
        ! endif

        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.
        
        ift=1; ilt=self%nt
        
        ! call alloc(sf_p_save,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])

        do it=ilt,ift,int(time_dir)
            if(mod(it,100)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_da%check_value
                call fld_a%check_value
                call fld_u%check_value
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then
                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_u%boundary_transport('load',it)
                call cpu_time(toc)
                tt1=tt1+toc-tic
                
                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%update_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_stresses_scattering(fld_da,fld_a,W2Idt,time_dir)
            call self%inject_stresses(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_stresses(fld_da,time_dir,it)
            call self%update_stresses(fld_a, time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: sfield%s_dt^it+0.5 \dot rfield%s^it+0.5
            !use sfield%v^it+1 to compute sfield%s_dt^it+0.5, as backward step 4
            if(if_compute_grad.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call gradient_moduli(fld_da,fld_u,it,cb%grad(:,:,:,2))
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif

            ! if(if_compute_imag.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call imaging(fld_u,fld_a,it,cb%imag)
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif

            ! !energy term of sfield
            ! if(self%if_compute_engy.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call energy(fld_u,it,cb%engy)
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
                
            !========================================================!

            ! if(present(o_sf)) then
                !backward step 2: v^it+1 -> v^it by FD of s^it+0.5
                call cpu_time(tic)
                call self%update_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt7=tt7+toc-tic

                !backward step 1: rm forces from v^it
                call cpu_time(tic)
                call self%inject_velocities(fld_u,time_dir,it)
                call cpu_time(toc)
                tt8=tt8+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 3: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_velocities(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities(fld_da,time_dir,it)
            call self%update_velocities(fld_a, time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            !adjoint step 1: sample v^it or s^it+0.5 at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_da,it)
                call cpu_time(toc)
                tt11=tt11+toc-tic
            endif
            
            ! !grho: sfield%v_dt^it \dot rfield%v^it
            ! !use sfield%s^it+0.5 to compute sfield%v_dt^it, as backward step 2
            ! if(if_compute_grad.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call gradient_density(o_sf,rf,it,cb%grad(:,:,:,1))
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
            
            !snapshot
            call fld_da%write(it,o_suffix='_rev')
            ! call fld_a%write(it,o_suffix='_rev')
            ! call fld_u%write(it,o_suffix='_rev')
            ! if(if_compute_imag) then
                ! call fld_lda%write_ext(it,'imag' ,cb%imag,size(cb%imag))
            ! endif
            if(if_compute_grad) then
                call fld_da%write_ext(it,'grad_da_star_Du' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
            endif
            ! if(self%if_compute_engy) then
                ! call fld_a%write_ext(it,'engy',cb%engy(:,:,:,1),size(cb%engy(:,:,:,1)))
            ! endif

        enddo
        
        !postprocess gradient
        if(if_compute_grad) call gradient_postprocess
        
        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to load boundary            ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update stresses          ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source stresses       ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update velocities        ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source velocities     ',tt8/mpiworld%max_threads
            write(*,*) 'Elapsed time ----------------------------'
            write(*,*) 'Elapsed time to add adjsource stresses   ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj stresses      ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to add adjsource velocities ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj velocities    ',tt10/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract&write fields     ',tt11/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate                ',tt6/mpiworld%max_threads

        endif

    end subroutine

    subroutine adjoint_da_star_Ddu(self,fld_da,fld_du,fld_a,fld_u,W2Idt,oif_record_adjseismo,oif_compute_grad)
        class(t_propagator) :: self
        type(t_field) :: fld_da,fld_a
        type(t_field) :: fld_du,fld_u
        real,dimension(cb%mz,cb%mx,cb%my) :: W2Idt
        logical,optional :: oif_record_adjseismo, oif_compute_grad
        
        real,parameter :: time_dir=-1. !time direction
        logical :: if_record_adjseismo, if_compute_imag, if_compute_grad

        ! if_record_adjseismo =either(oif_record_adjseismo,.false.,present(oif_record_adjseismo))
        ! if_compute_imag     =either(oif_compute_imag,    .false.,present(oif_compute_imag))
        if_compute_grad     =either(oif_compute_grad,    .false.,present(oif_compute_grad))
        ! self%if_compute_engy=self%if_compute_engy.and.(if_compute_imag.or.if_compute_grad)
        
        ! if(if_record_adjseismo)  call alloc(fld_da%seismo,1,self%nt)
        ! if(if_compute_imag)      call alloc(cb%imag,cb%mz,cb%mx,cb%my,self%nimag)
        if(if_compute_grad)      call alloc(cb%grad,cb%mz,cb%mx,cb%my,self%ngrad)
        ! if(self%if_compute_engy) call alloc(cb%engy,cb%mz,cb%mx,cb%my,self%nengy)

        !reinitialize absorbing boundary for incident wavefield reconstruction
        ! if(present(o_sf)) then
                                          fld_du%dvz_dz=0.
                                           fld_u%dvz_dz=0.
                                          fld_du%dvx_dx=0.
                                           fld_u%dvx_dx=0.
            ! if(allocated(fld_du%dvy_dy))  fld_du%dvy_dy=0.
            ! if(allocated( fld_u%dvy_dy))   fld_u%dvy_dy=0.
                                          fld_du%dp_dz=0.
                                           fld_u%dp_dz=0.
                                          fld_du%dp_dx=0.
                                           fld_u%dp_dx=0.
        ! endif

        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.
        
        ift=1; ilt=self%nt

        ! call alloc(sf_p_save,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        
        do it=ilt,ift,int(time_dir)
            if(mod(it,100)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_da%check_value
                call fld_a%check_value
                call fld_du%check_value
                call fld_u%check_value
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then
                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_du%boundary_transport('load',it)
                call  fld_u%boundary_transport('load',it)
                call cpu_time(toc)
                tt1=tt1+toc-tic
                
                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%update_stresses(fld_du,time_dir,it)
                call self%update_stresses( fld_u,time_dir,it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_stresses_scattering(fld_du,fld_u,W2Idt,time_dir)
                call self%inject_stresses( fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_stresses_scattering(fld_da,fld_a,W2Idt,time_dir)
            call self%inject_stresses(fld_a,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_stresses(fld_da,time_dir,it)
            call self%update_stresses(fld_a, time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: sfield%s_dt^it+0.5 \dot rfield%s^it+0.5
            !use sfield%v^it+1 to compute sfield%s_dt^it+0.5, as backward step 4
            if(if_compute_grad.and.mod(it,irdt)==0) then
                call cpu_time(tic)
                call gradient_moduli(fld_da,fld_du,it,cb%grad(:,:,:,2))
                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif

            ! if(if_compute_imag.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call imaging(fld_a,fld_u,it,cb%imag)
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif

            ! !energy term of sfield
            ! if(self%if_compute_engy.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call energy(fld_u,it,cb%engy)
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
                
            !========================================================!

            ! if(present(o_sf)) then
                !backward step 2: v^it+1 -> v^it by FD of s^it+0.5
                call cpu_time(tic)
                call self%update_velocities(fld_du,time_dir,it)
                call self%update_velocities( fld_u,time_dir,it)
                call cpu_time(toc)
                tt7=tt7+toc-tic

                !backward step 1: rm forces from v^it
                call cpu_time(tic)
                ! call self%inject_velocities(fld_du,fld_u,it)
                call self%inject_velocities( fld_u,time_dir,it)
                call cpu_time(toc)
                tt8=tt8+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 3: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_velocities(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities(fld_da,time_dir,it)
            call self%update_velocities(fld_a, time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            !adjoint step 1: sample v^it or s^it+0.5 at source position
            if(if_record_adjseismo) then
                call cpu_time(tic)
                call self%extract(fld_da,it)
                call cpu_time(toc)
                tt11=tt11+toc-tic
            endif
            
            ! !grho: sfield%v_dt^it \dot rfield%v^it
            ! !use sfield%s^it+0.5 to compute sfield%v_dt^it, as backward step 2
            ! if(if_compute_grad.and.mod(it,irdt)==0) then
            !     call cpu_time(tic)
            !     call gradient_density(o_sf,rf,it,cb%grad(:,:,:,1))
            !     call cpu_time(toc)
            !     tt6=tt6+toc-tic
            ! endif
            
            !snapshot
            call fld_da%write(it,o_suffix='_rev')
            ! call fld_a%write(it,o_suffix='_back3')
            ! call fld_u%write(it,o_suffix='_back3')
            ! if(if_compute_imag) then
                ! call fld_lda%write_ext(it,'imag' ,cb%imag,size(cb%imag))
            ! endif
            if(if_compute_grad) then
                call fld_da%write_ext(it,'grad_da_star_Ddu' ,cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
            endif
            ! if(self%if_compute_engy) then
                ! call fld_a%write_ext(it,'engy',cb%engy(:,:,:,1),size(cb%engy(:,:,:,1)))
            ! endif

        enddo
        
        !postprocess gradient
        if(if_compute_grad) call gradient_postprocess
        
        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to load boundary            ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update stresses          ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source stresses       ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update velocities        ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source velocities     ',tt8/mpiworld%max_threads
            write(*,*) 'Elapsed time ----------------------------'
            write(*,*) 'Elapsed time to add adjsource stresses   ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj stresses      ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to add adjsource velocities ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj velocities    ',tt10/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract&write fields     ',tt11/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate                ',tt6/mpiworld%max_threads

        endif

    end subroutine

    subroutine adjoint_WPI(self,fld_da,fld_a,fld_Adj_du,fld_du,fld_u,W2Idt,corrs)
        class(t_propagator) :: self
        type(t_field) :: fld_da,fld_a, fld_Adj_du
        type(t_field) :: fld_du,fld_u
        real,dimension(cb%mz,cb%mx,cb%my) :: W2Idt
        character(*) :: corrs
        
        real,parameter :: time_dir=-1. !time direction
        
        !reinitialize absorbing boundary for incident wavefield reconstruction
        ! if(present(o_sf)) then
                                          fld_du%dvz_dz=0.
                                           fld_u%dvz_dz=0.
                                          fld_du%dvx_dx=0.
                                           fld_u%dvx_dx=0.
            ! if(allocated(fld_du%dvy_dy))  fld_du%dvy_dy=0.
            ! if(allocated( fld_u%dvy_dy))   fld_u%dvy_dy=0.
                                          fld_du%dp_dz=0.
                                           fld_u%dp_dz=0.
                                          fld_du%dp_dx=0.
                                           fld_u%dp_dx=0.
        ! endif

        !timing
        tt1=0.; tt2=0.; tt3=0.
        tt4=0.; tt5=0.; tt6=0.
        tt7=0.; tt8=0.; tt9=0.
        tt10=0.;tt11=0.; tt12=0.; tt13=0.
        
        ift=1; ilt=self%nt

        ! call alloc(sf_p_save,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily])
        
        do it=ilt,ift,int(time_dir)
            if(mod(it,100)==0 .and. mpiworld%is_master) then
                write(*,*) 'it----',it
                call fld_da%check_value
                call fld_a%check_value
                call fld_Adj_du%check_value
                call fld_du%check_value
                call fld_u%check_value
            endif            

            !do backward time stepping to reconstruct the source (incident) wavefield
            !and adjoint time stepping to compute the receiver (adjoint) field
            !step# conforms with forward time stepping

            ! if(present(o_sf)) then
                !backward step 6: retrieve v^it+1 at boundary layers (BC)
                call cpu_time(tic)
                call fld_du%boundary_transport('load',it)
                call  fld_u%boundary_transport('load',it)
                call cpu_time(toc)
                tt1=tt1+toc-tic
                
                !backward step 4: s^it+1.5 -> s^it+0.5 by FD of v^it+1
                call cpu_time(tic)
                call self%update_stresses(fld_du,time_dir,it)
                call self%update_stresses( fld_u,time_dir,it)
                call cpu_time(toc)
                tt2=tt2+toc-tic

                !backward step 3: rm pressure from s^it+0.5
                call cpu_time(tic)
                call self%inject_stresses_scattering(fld_du,fld_u,W2Idt,time_dir)
                call self%inject_stresses(fld_u,time_dir,it)
                call cpu_time(toc)
                tt3=tt3+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 5: inject to s^it+1.5 at receivers
            call cpu_time(tic)
            call self%inject_stresses_scattering(fld_da,fld_a,W2Idt,time_dir)
            call self%inject_stresses(fld_a,time_dir,it)
            call self%inject_stresses(fld_Adj_du,time_dir,it)
            call cpu_time(toc)
            tt4=tt4+toc-tic

            !adjoint step 4: s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
            call cpu_time(tic)
            call self%update_stresses(fld_da,time_dir,it)
            call self%update_stresses(fld_a, time_dir,it)
            call self%update_stresses(fld_Adj_du, time_dir,it)
            call cpu_time(toc)
            tt5=tt5+toc-tic

            !gkpa: sfield%s_dt^it+0.5 \dot rfield%s^it+0.5
            !use sfield%v^it+1 to compute sfield%s_dt^it+0.5, as backward step 4
            if(mod(it,irdt)==0) then
                call cpu_time(tic)
                
                !a★Du
                call gradient_moduli(fld_a,fld_u,it,cb%corr(:,:,:,1))

                if(index(corrs,'RE')>0) then
                    !δa★Du
                    call gradient_moduli(fld_da,fld_u,it,cb%corr(:,:,:,2))
                    !δu★Da
                    call gradient_moduli(fld_du,fld_a,it,cb%corr(:,:,:,3))
                endif
                
                ! if(index(corrs,'2ndMI')>0) then
                !     !corr(:,:,:,4) = da_star_Ddu
                !     call gradient_moduli(fld_da,fld_du,it,cb%corr(:,:,:,4))
                ! endif
                
                if(index(corrs,'DR')>0) then
                    !Adj(Rᴴδu)★Du
                    call gradient_moduli(fld_Adj_du,fld_u,it,cb%corr(:,:,:,5))
                    !will add the minus sign outside timestepping loop
                endif

                call cpu_time(toc)
                tt6=tt6+toc-tic
            endif
                
            !========================================================!

            ! if(present(o_sf)) then
                !backward step 2: v^it+1 -> v^it by FD of s^it+0.5
                call cpu_time(tic)
                call self%update_velocities(fld_du,time_dir,it)
                call self%update_velocities( fld_u,time_dir,it)
                call cpu_time(toc)
                tt7=tt7+toc-tic

                !backward step 1: rm forces from v^it
                call cpu_time(tic)
                ! call self%inject_velocities(fld_du,fld_u,it)
                call self%inject_velocities( fld_u,time_dir,it)
                call cpu_time(toc)
                tt8=tt8+toc-tic
            ! endif

            !--------------------------------------------------------!

            !adjoint step 3: inject to v^it+1 at receivers
            call cpu_time(tic)
            call self%inject_velocities(fld_a,time_dir,it)
            call cpu_time(toc)
            tt9=tt9+toc-tic

            !adjoint step 2: v^it+1 -> v^it by FD^T of s^it+0.5
            call cpu_time(tic)
            call self%update_velocities(fld_da,time_dir,it)
            call self%update_velocities(fld_a, time_dir,it)
            call self%update_velocities(fld_Adj_du, time_dir,it)
            call cpu_time(toc)
            tt10=tt10+toc-tic
            
            ! !adjoint step 1: sample v^it or s^it+0.5 at source position
            ! if(if_record_adjseismo) then
            !     call cpu_time(tic)
            !     call self%extract(fld_da,it)
            !     call cpu_time(toc)
            !     tt11=tt11+toc-tic
            ! endif
                        
            !snapshot
            call fld_u%write(it,o_suffix='_rev')
            call fld_du%write(it,o_suffix='_rev')
            call fld_a%write(it,o_suffix='_rev')
            call fld_da%write(it,o_suffix='_rev')
            
                call fld_u%write_ext(it,'grad_a_star_Du' ,cb%corr(:,:,:,1),size(cb%corr(:,:,:,1)))

            if(index(corrs,'RE')>0) then
                call fld_u%write_ext(it,'grad_du_star_Da' ,cb%corr(:,:,:,2),size(cb%corr(:,:,:,2)))
                call fld_u%write_ext(it,'grad_da_star_Du' ,cb%corr(:,:,:,3),size(cb%corr(:,:,:,3)))
            endif

            if(index(corrs,'DR')>0) then
                call fld_u%write_ext(it,'grad_Adj(du)_star_Du' ,cb%corr(:,:,:,5),size(cb%corr(:,:,:,5)))
            endif

        enddo
        
        !postprocess gradient
        
        !scale by m%cell_volume*rdt tobe a gradient in the discretized world
        cb%corr = cb%corr*m%cell_volume*rdt

        !gkpa
        do i=1,5
        cb%corr(:,:,:,i) = cb%corr(:,:,:,i) * (-ppg%inv_kpa(1:cb%mz,1:cb%mx,1:cb%my))
        enddo

        !preparing for cb%project_back
        cb%corr(1,:,:,:) = cb%corr(2,:,:,:)

        
        if(mpiworld%is_master) then
            write(*,*) 'Elapsed time to load boundary            ',tt1/mpiworld%max_threads
            write(*,*) 'Elapsed time to update stresses          ',tt2/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source stresses       ',tt3/mpiworld%max_threads
            write(*,*) 'Elapsed time to update velocities        ',tt7/mpiworld%max_threads
            write(*,*) 'Elapsed time to rm source velocities     ',tt8/mpiworld%max_threads
            write(*,*) 'Elapsed time ----------------------------'
            write(*,*) 'Elapsed time to add adjsource stresses   ',tt4/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj stresses      ',tt5/mpiworld%max_threads
            write(*,*) 'Elapsed time to add adjsource velocities ',tt9/mpiworld%max_threads
            write(*,*) 'Elapsed time to update adj velocities    ',tt10/mpiworld%max_threads
            write(*,*) 'Elapsed time to extract&write fields     ',tt11/mpiworld%max_threads
            write(*,*) 'Elapsed time to correlate                ',tt6/mpiworld%max_threads

        endif

        call hud('Viewing the snapshots (if written) with SU ximage/xmovie:')
        call hud('ximage < snap_rfield%*  n1='//num2str(cb%nz)//' perc=99')
        call hud('xmovie < snap_rfield%*  n1='//num2str(cb%nz)//' n2='//num2str(cb%nx)//' clip=?e-?? loop=2 title=%g')
        call hud('ximage < snap_*  n1='//num2str(cb%mz)//' perc=99')
        call hud('xmovie < snap_*  n1='//num2str(cb%mz)//' n2='//num2str(cb%mx)//' clip=?e-?? loop=2 title=%g')

    end subroutine

    
    subroutine inject_stresses_scattering(self,fld_d,fld,W2Idt,time_dir)
        class(t_propagator) :: self
        type(t_field) :: fld_d, fld
        real,dimension(cb%mz,cb%mx,cb%my) :: W2Idt
        real :: time_dir

        if(.not.fld%is_adjoint) then
            fld_d%p(1:cb%mz,1:cb%mx,1:cb%my) = fld_d%p(1:cb%mz,1:cb%mx,1:cb%my) + time_dir*W2Idt*fld%p(1:cb%mz,1:cb%mx,1:cb%my)*self%kpa(1:cb%mz,1:cb%mx,1:cb%my)
            return
        endif

            fld_d%p(1:cb%mz,1:cb%mx,1:cb%my) = fld_d%p(1:cb%mz,1:cb%mx,1:cb%my) +          W2Idt*fld%p(1:cb%mz,1:cb%mx,1:cb%my)*self%kpa(1:cb%mz,1:cb%mx,1:cb%my)

    end subroutine

end
