module m_nabla
use m_System
use m_Modeling
use m_smoother_laplacian_sparse
use m_parametrizer
use m_fobjective

    private

    type,public :: t_nabla
        contains
        procedure :: act => act
    end type

    type(t_nabla),public :: nabla

    real,dimension(:,:,:),allocatable :: mask

    contains

    !nabla f: gradient of f
    subroutine act(self,fobj,qp,oif_update_m,oif_approx)
        class(t_nabla) :: self
        type(t_fobjective) :: fobj
        type(t_querypoint) :: qp
        logical,optional :: oif_update_m,oif_approx

        type(t_string),dimension(:),allocatable :: smoothings
        character(:),allocatable :: smask

        !update model
        if(either(oif_update_m,.false.,present(oif_update_m))) then
            call param%transform(o_dir='x->m',o_x=qp%x)
        endif

        !compute dnorm's gradient by adjoint-state method
        if(either(oif_approx,.false.,present(oif_approx))) then
            call gradient_modeling_approximate(fobj)
        else
            call gradient_modeling(fobj)
        endif

        !post-modeling smoothing in physical (model) domain
        smoothings=setup%get_strs('SMOOTHING','SMTH',o_default='Laplacian')

        do i=1,size(smoothings)
            !Laplacian smoothing
            if(smoothings(i)%s=='Laplacian') then
                call hud('Laplacian smoothing')
                call smoother_laplacian_init([m%nz,m%nx,m%ny],[m%dz,m%dx,m%dy],shot%fpeak)
                do j=1,ppg%ngrad
                    !call smoother_laplacian_extend_mirror(cb%gradient(:,:,:,i),m%itopo)
                    call smoother_laplacian_pseudo_nonstationary(cb%grad(:,:,:,j),m%vp)
                enddo    
            endif
        enddo

        !freeze_zone as hard mask
        !deprecated, hard mask is now applied on the model side
        !where(m%is_freeze_zone) fobj%gradient=0.

        !soft mask
        smask=setup%get_file('GRADIENT_SOFT_MASK')
        if(smask/='') then
            call alloc(mask,m%nz,m%nx,m%ny)
            ! call sysio_read(smask,mask,size(mask))
            
            do i=1,ppg%ngrad
                cb%grad(:,:,:,i)=cb%grad(:,:,:,i)*mask
            enddo
        endif

        !reparameterization
        call param%transform(o_g=qp%g)
        !call param%transform_model('m->x',sfield%autocorr)

        !Tikhonov regularization
        if(either(oif_approx,.false.,present(oif_approx))) then
            call regularize_approximate(fobj,qp)
        else
            call regularize(fobj,qp)
        endif

        qp%f=fobj%total_loss()

    end subroutine
    
    !nabla^2 f: Hessian of f
    subroutine act2
    end subroutine

    !nabla dot v: divergence of vector v
    subroutine dot
    end subroutine

    !nabla curl v: curl (rotation) of vector v
    subroutine curl
    end subroutine

    ! subroutine gradient_modeling(fobj)
    !     type(t_fobjective) :: fobj

    !     type(t_field) :: sfield, rfield
    !     character(:),allocatable :: update_wavelet

    !     type(t_checkpoint),save :: chp

    !     call alloc(fobj%gradient,m%nz,m%nx,m%ny,ppg%ngrad)
                
    !     call hud('===== START LOOP OVER SHOTS =====')
        
    !     call chp%init('FWI_shotloop','Init# Shot#','per_init given')
    !     do i=1,shls%nshots_per_processor
    !         call chp%count(str2int(shls%yield(i)))

    !         call shot%init(shls%yield(i))
    !         call shot%read_from_data
    !         call shot%set_var_time
    !         call shot%set_var_space(index(ppg%info,'FDSG')>0)

    !         call hud('Modeling shot# '//shot%sindex)
            
    !         call cb%init
    !         call cb%project(is_fdsg=index(ppg%info,'FDSG')>0,abslayer_width=cpml%nlayer)

    !         call ppg%check_discretization
            
    !         call ppg%init

    !         call ppg%init_field(sfield,name='sfield',origin='src',oif_will_reconstruct=.true.)

    !         call sfield%add_RHS
            
    !         !forward modeling
    !         if(.not.sfield%is_registered(chp,'seismo comp boundary')) then
    !             call ppg%forward(sfield)
    !             call sfield%register(chp,'seismo comp boundary')
    !         endif

    !         call sfield%acquire

    !         if(mpiworld%is_master) call suformat_write('draw_'//shot%sindex,shot%dsyn,shot%nt,shot%nrcv,o_dt=shot%dt)
            
    !         update_wavelet=setup%get_str('UPDATE_WAVELET',o_default='per shot')

    !         if(update_wavelet/='no') call shot%update_wavelet !call gradient_matchfilter_data
            
    !         !write synthetic data
    !         call suformat_write('dsyn_'//shot%sindex,shot%dsyn,shot%nt,shot%nrcv,o_dt=shot%dt,o_sindex=shot%sindex)
                        
    !         !objective function
    !         call fobj%compute_dnorms
    !         ! call fobj%compute_adjsource
        
    !         !adjoint source
    !         if(update_wavelet/='no') call shot%update_adjsource

    !         call ppg%init_field(rfield,'rfield','rcv')

    !         call rfield%add_RHS(ois_adjoint=.true.)
            
    !         call alloc(cb%kernel,cb%mz,cb%mx,cb%my,ppg%ngrad)

    !         !adjoint modeling
    !         if(.not.cb%is_registered(chp,'kernel')) then
    !             call ppg%adjoint(rfield,o_sf=sfield,o_grad=cb%kernel)
    !             call cb%register(chp,'kernel')
    !         endif

    !         call cb%project_back(fobj%gradient,cb%kernel,ppg%ngrad)
            
    !     enddo
        
    !     call hud('        END LOOP OVER SHOTS        ')
        
    !     !collect global objective function value
    !     call mpi_allreduce(MPI_IN_PLACE, fobj%dnorms, fobj%n_dnorms, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    !     call fobj%print_dnorms

    !     !collect global gradient
    !     call mpi_allreduce(MPI_IN_PLACE, fobj%gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)

    !     call mpiworld%barrier

    ! end subroutine

    subroutine gradient_modeling_approximate(fobj)
        type(t_fobjective) :: fobj

        call gradient_modeling(fobj)

    end subroutine

    subroutine regularize(fobj,qp)
        type(t_fobjective) :: fobj
        type(t_querypoint) :: qp

        call fobj%compute_xnorms(qp%x,qp%g)

    end subroutine

    subroutine regularize_approximate(fobj,qp)
        type(t_fobjective) :: fobj
        type(t_querypoint) :: qp

        call regularize(fobj,qp)

    end subroutine

end