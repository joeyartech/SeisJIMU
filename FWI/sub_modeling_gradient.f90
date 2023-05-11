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
    character(:),allocatable :: update_wavelet
    type(t_weighter) :: wei_wl
    type(t_field) :: fld_u, fld_a
    type(t_correlate) :: u_star_u,a_star_u

    type :: t_S
        real,dimension(:),allocatable :: scale
    end type
    type(t_S),dimension(:),allocatable,save :: S
    real,dimension(:,:),allocatable :: tmp_dsyn, gmwindow
    
    if(is_first_in) allocate(S(shls%nshots_per_processor)) !then can NOT randomly sample shots..

    
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
        !call ppg%init_correlate(u_star_u,'u_star_u','energy')
                
        !forward modeling
        !call ppg%forward(fld_u,u_star_u); call fld_u%acquire
        call ppg%forward(fld_u); call fld_u%acquire

        if(mpiworld%is_master) call shot%write('draw_',shot%dsyn)

        !update wavelet
	update_wavelet=setup%get_str('UPDATE_WAVELET',o_default='per shot')
        if(update_wavelet/='no') then
            call wei_wl%update('_4WAVELET')
            call shot%update_wavelet(wei_wl%weight)
        endif
        !write synthetic data
        call shot%write('dsyn_',shot%dsyn)

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
        call wei%update
        
        select case(setup%get_str('DNORM',o_default='L2'))
        case('L2')
            fobj%misfit = fobj%misfit &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)

        case('L2_scaled')
            if(is_first_in) call alloc(S(i)%scale,shot%nrcv) !then can NOT randomly sample shots..
            call alloc(tmp_dsyn,shot%nt,shot%nrcv)
            do j=1,shot%nrcv
                if(is_first_in) S(i)%scale(j) = either(0., maxval(abs(shot%dobs(:,j))) / maxval(abs(shot%dsyn(:,j))) , shot%rcv(j)%is_badtrace)
                tmp_dsyn(:,j)=shot%dsyn(:,j)*S(i)%scale(j)
            enddo
            if(is_first_in) then
            	open(12,file=dir_out//'dobs_dsyn_max_ratio',access='direct',recl=4*shot%nrcv)
            	write(12,rec=shot%index) S(i)%scale
            	close(12)
            endif
            
            !check if S is changing..
!            if(shot%index==1)   print*, 'on '//shot%sindex,i,S(i)%scale
!            if(shot%index==112) print*, 'on '//shot%sindex,i,S(i)%scale

            fobj%misfit = fobj%misfit &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-tmp_dsyn, shot%dt)

        case('imaging')
            fobj%misfit = fobj%misfit &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs, shot%dt)
                
        end select

        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call kernel_L2sq(shot%dadj)
        call shot%write('dadj_',shot%dadj)

        
        ! if(mpiworld%is_master) call fobj%print_dnorms('Shotloop-stacked','upto '//shot%sindex)
                
        !adjoint source
        if(update_wavelet/='no') call shot%update_adjsource
        
        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.); call fld_a%ignite
        call ppg%init_correlate(a_star_u,'a_star_u','gradient')
        
        !adjoint modeling
        call ppg%adjoint(fld_a,fld_u,o_a_star_u=a_star_u)


!output grad for each shot
if(setup%get_bool('IF_SAVE_GRADIENT_PER_SHOT')) then
open(12,file=dir_out//'a_star_u_'//shot%sindex,action='write',access='direct',recl=4*m%n)
write(12,rec=1) a_star_u%rp_div_sv
write(12,rec=2) a_star_u%rv_grad_sp
close(12)
endif
!gaussian masking on source singularities
if(setup%get_bool('IF_MASK_SINGULARITY')) then
    nmute=nint(1.2*1900/shot%fpeak/m%dz) !=19
    !nmute=nint(1.2*m%vp(shot%src%iz,shot%src%ix)/shot%fpeak/m%dz)
    jsigma=3
    call alloc(gmwindow,[-nmute,nmute],[-nmute,nmute])
    do jx=-nmute,nmute
    do jz=-nmute,nmute
        R=sqrt(1.0*(jz**2+jx**2))
        if(R>nmute) R=nmute
        gmwindow(jz,jx) = exp(-0.5*((nmute-R)/jsigma)**2)
    enddo
    enddo

    do jx=-nmute,+nmute
        ix=shot%src%ix+jx; if(ix<1.or.ix>m%nx) cycle
    do jz=-nmute,+nmute
        iz=shot%src%iz+jz; if(iz<1.or.iz>m%nz) cycle
        correlate_gradient(iz,ix,1,:)=correlate_gradient(iz,ix,1,:)*gmwindow(jz,jx)
    enddo
    enddo
endif
          
  
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    is_first_in=.false.

    if(mpiworld%is_master) then
        call execute_command_line('cat '//dir_out//'updated_wavelet_Shot????.su > '//dir_out//'updated_wavelet.su')
    endif

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
