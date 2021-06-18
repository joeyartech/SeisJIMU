module m_nabla
use m_string
use m_setup
use m_mpienv
use m_format_su
use m_checkpoint
use m_model
use m_shotlist
use m_shot
use m_computebox
use m_field
use m_cpml
use m_propagator
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse

    private

    type,public :: t_nabla
        contains
        procedure :: act => act
    end type

    type(t_nabla),public :: nabla

    contains

    !nabla f: gradient of f
    subroutine act(self,fobj,oif_approx)
        class(t_nabla) :: self
        type(t_fobjective) :: fobj
        logical,optional :: oif_approx

        type(t_string),dimension(:),allocatable :: smoothings
        character(:),allocatable :: smask

        call alloc(fobj%gradient,m%nz,m%nx,m%ny,propagator%ngrad)

        !compute gradient by adjoint-state method
        if(either(oif_approx,.false.,present(oif_approx))) then
            call gradient_modeling_approximate(fobj)
        else
            call gradient_modeling(fobj)
        endif

        !postprocessing

        !smoothing
        smoothings=setup%get_strs('SMOOTHING','SMTH',o_default='Laplacian')

        do i=1,size(smoothings)
            !Laplcian smoothing
            if(smoothings(i)%s=='Laplacian') then
                call hud('Laplacian smoothing')
                call smoother_laplacian_init([m%nz,m%nx,m%ny],[m%dz,m%dx,m%dy],shot%fpeak)
                do j=1,propagator%ngrad
                    !call smoother_laplacian_extend_mirror(fobj%gradient(:,:,:,i),m%itopo)
                    call smoother_laplacian_pseudo_nonstationary(fobj%gradient(:,:,:,j),m%vp)
                enddo    
            endif
        enddo

        !freeze_zone as hard mask
        !deprecated
        !where(m%is_freeze_zone) fobj%gradient=0.

        !soft mask
        smask=setup%get_file('GRADIENT_SOFT_MASK')
        if(smask/='') then
            call alloc(mask,m%nz,m%nx,m%ny)
            open(12,file=smask,access='direct',recl=4*m%n,action='read',status='old')
            read(12,rec=1) mask
            close(12)
            
            do i=1,propagator%ngrad
                fobj%gradient(:,:,:,i)=fobj%gradient(:,:,:,i)*mask
            enddo
        endif

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

    subroutine gradient_modeling(fobj)
        type(t_fobjective) :: fobj

        type(t_field) :: sfield, rfield
        character(:),allocatable :: update_wavelet

        type(t_checkpoint),save :: chp
                
        call hud('===== START LOOP OVER SHOTS =====')
        
        call chp%init('FWI_shotlist','Init# Shot#','per_init given')
        do i=1,sl%nshot_per_processor
            call chp%count(sl%list_per_processor(i))

            call shot%init(sl%list_per_processor(i))
            call shot%read_from_data
            call shot%set_var_time
            call shot%set_var_space(index(propagator%info,'FDSG')>0)

            call hud('Modeling shot# '//shot%sindex)
            
            call cb%init
            call cb%project(is_fdsg=index(propagator%info,'FDSG')>0,abslayer_width=cpml%nlayer)

            call propagator%check_discretization
            
            call propagator%init

            call propagator%init_field(sfield,name='sfield',origin='src',oif_will_reconstruct=.true.)

            call sfield%add_RHS
            
            !forward modeling
            if(.not.sfield%is_registered(chp,'seismo comp boundary')) then
                call propagator%forward(sfield)
                call sfield%register(chp,'seismo comp boundary')
            endif

            call sfield%acquire

            if(mpiworld%is_master) call suformat_write('draw_'//shot%sindex,shot%dsyn,shot%nt,shot%nrcv,o_dt=shot%dt)
            
            update_wavelet=setup%get_str('UPDATE_WAVELET',o_default='per shot')

            if(update_wavelet/='no') call shot%update_wavelet !call gradient_matchfilter_data
            
            !write synthetic data
            call suformat_write('dsyn_'//shot%sindex,shot%dsyn,shot%nt,shot%nrcv,o_dt=shot%dt,o_sindex=shot%sindex)
                        
            !objective function
            call fobj%compute_dnorms
            call fobj%compute_adjsource
        
            !adjoint source
            if(update_wavelet/='no') call shot%update_adjsource

            call propagator%init_field(rfield,'rfield','rcv')

            call rfield%add_RHS(ois_adjoint=.true.)
            
            call alloc(cb%kernel,cb%mz,cb%mx,cb%my,propagator%ngrad)

            !adjoint modeling
            if(.not.cb%is_registered(chp,'kernel')) then
                call propagator%adjoint(rfield,o_sf=sfield,o_grad=cb%kernel)
                call cb%register(chp,'kernel')
            endif

            !put cb%kernel into global gradient
            call cb%project_back(fobj%gradient,cb%kernel,propagator%ngrad)
            
        enddo
        
        call hud('        END LOOP OVER SHOTS        ')
        
        !collect global objective function value
        call mpi_allreduce(MPI_IN_PLACE, fobj%dnorms, fobj%n_dnorms, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        call fobj%print_dnorms

        !collect global gradient
        call mpi_allreduce(MPI_IN_PLACE, fobj%gradient, m%n*propagator%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)

        call mpiworld%barrier

    end subroutine


    subroutine gradient_modeling_approximate(fobj)
        type(t_fobjective) :: fobj
    end subroutine

end
