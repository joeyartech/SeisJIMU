subroutine modeling_gradient_ip
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

    character(:),allocatable :: s_update_wavelet
    type(t_weighter) :: wei_wl

    type(t_field) :: fld_u, fld_a
    type(t_correlate) :: a_star_u

    type(t_suformat) :: onlydir
    
    !character(:),allocatable :: s_job
    
    !RWI misfit
    fobj%misfit=0.
    fobj%reflection=0.
    fobj%diving=0.

    call alloc(correlate_gradient,m%nz,m%nx,m%ny,ppg%ngrad)
    
    call hud('------------------------')
    call hud('     Build Ip model     ')
    call hud('------------------------')

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


        call hud('----  Solving A(m₀)u₀=s  ----')
        call ppg%init_field(fld_u, name='fld_u0');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire; call shot%write('Ru0_',shot%dsyn)

        if(setup%get_str('JOB')=='forward') cycle

        call ppg%init_field(fld_a,name='fld_a0',ois_adjoint=.true.)

        s_update_wavelet=setup%get_str('UPDATE_WAVELET')
        if(s_update_wavelet/='') call wei_wl%update(o_suffix='_4WAVELET')

        call hud('----  Computing obj func & dadj  ----')

            ! call sepa%update
            call wei%update!('_4IMAGING')
            call alloc(shot%dadj,shot%nt,shot%nrcv)

            select case (setup%get_str('DATA_NORM','DNORM',o_default='L2'))

            case ('L2'); call hud('0.5|| W(f*u - d)||² => adjsrc = f★ W W(f*u - d)')

if(s_update_wavelet/='') then
call hud('update wavelet')
call shot%update_wavelet(wei_wl%weight) !call matchfilter_apply_to_data(shot%dsyn)
call shot%write('updated_Ru0_',shot%dsyn)
call suformat_write('updated_wavelet_'//shot%sindex,shot%wavelet,shot%nt,1,shot%dt)
endif

                call wei%update
                call alloc(shot%dadj,shot%nt,shot%nrcv)

                fobj%misfit = fobj%misfit &                                                                   
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-shot%dsyn, shot%dt)  

                call kernel_L2sq(shot%dadj)
          
if(s_update_wavelet/='') then
call hud('update adjoint source')
call shot%update_adjsource
endif


            case('L2_nodir'); call hud('0.5|| W(u-u0 - d)||² => adjsrc = W W(u-u0 - d)')
                if(is_first_in) call shot%write('onlydir_',shot%dsyn)
                call onlydir%read(dir_out//'onlydir_'//shot%sindex//'.su',shot%sindex)
                
                call wei%update
                call alloc(shot%dadj,shot%nt,shot%nrcv)

                fobj%misfit = fobj%misfit &
                    + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, shot%dobs-(shot%dsyn-onlydir%trs), shot%dt)

                call kernel_L2sq(shot%dadj)


            case default
                call error('No DNORM specified!')

            end select



        call hud('remute') !0 should be still 0 in the adjoint source, esp below first arrivals
        shot%dadj=shot%dadj*wei%weight / (wei%weight+r_epsilon)

        call fld_a%ignite(o_wavelet=shot%dadj)
        call shot%write('dadj_',shot%dadj)


        call hud('----  Solving A(m₀)ᴴa₀ = Rʳ(d-Ru₀) & xcorrelate  ----')

        call ppg%init_correlate(a_star_u,'a0_star_u0')
        call ppg%adjoint(fld_a,fld_u,a_star_u)

if(allocated(a_star_u%gepsr)) call sysio_write('g0epsr_'//shot%sindex,a_star_u%gepsr,cb%mz*cb%mx)
if(allocated(a_star_u%gmur )) call sysio_write('g0mur_' //shot%sindex,a_star_u%gmur ,cb%mz*cb%mx)
! if(allocated(a_star_u%gsgma)) call sysio_write('g0sgma_'//shot%sindex,a_star_u%gsgma,cb%mz*cb%mx)

        call hud('----  Assemble  ----')
        call ppg%assemble(a_star_u)

        call hud('---------------------------------')

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')


    !allreduce misfit values
    call mpi_allreduce(mpi_in_place, [fobj%misfit], 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked misfit '//num2str(fobj%misfit))

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
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_u0, fld_u, fld_a0, fld_a
    type(t_correlate) :: a0_star_u0, a_star_u

    real,dimension(:,:),allocatable :: tmp_dsyn, gmwindow
    real,dimension(:,:),allocatable :: dadj_refl

    ! character(:),allocatable :: s_job
    ! real,dimension(:,:),allocatable :: Wdres
    !real,dimension(:,:,:),allocatable :: Ddt2
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff
    
    !RWI misfit
    fobj%misfit=0.
    fobj%reflection=0.
    fobj%diving=0.

    call alloc(correlate_gradient,m%nz,m%nx,m%ny,ppg%ngrad)    

    call hud('-------------------------')
    call hud('     Update Vp model     ')
    call hud('-------------------------')

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
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire;  call shot%write('Ru_',shot%dsyn)

        call mpiworld%barrier

        call hud('----  Computing reflection obj func & dadj  ----')
        ! call sepa%update
        call wei%update
        call alloc(shot%dadj,shot%nt,shot%nrcv)

        select case (setup%get_str('DATA_NORM','DNORM',o_default='L2_grad'))
        case ('L2_grad'); call hud('0.5|| W(Du - d)||² => adjsrc = DᵀW²(Du - d)')

            call compute_L2_grad()

        case ('L2_slope'); call hud('0.5|| W(pDu - pDd)||² => adjsrc = Dᵀ[∂p/∂Du]W²Δp')

            call compute_L2_slope()

        case default
            call error('No DNORM specified!')

        end select
                    
        call shot%write('dadj_refl_',shot%dadj)

        call hud('----  Solving A(m)ᴴa = Rʳ(d-Ru) & xcorrelate  ----')
        call ppg%init_field(fld_a ,name='fld_a' ,ois_adjoint=.true.); call fld_a%ignite
        call ppg%init_correlate(a_star_u,'a_star_u') !a★u
        call ppg%adjoint(fld_a,fld_u,a_star_u)


        call hud('----  Solving A(m₀)u₀=s  ----')
        call cb%project(ois_background=.true.)
        call ppg%init
        call ppg%init_field(fld_u0,name='fld_u0');    call fld_u0%ignite
        call ppg%forward(fld_u0)
        call fld_u0%acquire(o_seismo=shot%dsyn_aux);  call shot%write('Ru0_',shot%dsyn_aux)

        call mpiworld%barrier

        call ppg%init_field(fld_a,name='fld_a0',ois_adjoint=.true.)

        call hud('----  Computing diving obj func & dadj  ----')

        select case (setup%get_str('DATA_NORM','DNORM',o_default='L2'))
        case ('L2_grad'); !call hud('0.5|| W(Du - d)||² => adjsrc = DᵀW²(Du - d)')

            fobj%diving =0.
            shot%dadj = -shot%dadj

        case ('L2_slope'); !call hud('0.5|| W(pDu - pDd)||² => adjsrc = Dᵀ[∂p/∂Du]W²Δp')

            fobj%diving =0.
            shot%dadj = -shot%dadj

        case default
            call error('No DNORM specified!')

        end select

        call shot%write('dadj_div-refl_',shot%dadj)
        
        call hud('----  Solving A(m₀)ᴴa₀ = Rᵈ(d-Ru)-Rʳ(d-Ru) & xcorrelate  ----')
        call ppg%init_field(fld_a0,name='fld_a0',ois_adjoint=.true.); call fld_a0%ignite
        call ppg%init_correlate(a0_star_u0,'a0_star_u0') !a₀★u₀
        call ppg%adjoint(fld_a0,fld_u0,a0_star_u0)
        

call sysio_write('gepsr_'//shot%sindex,   a_star_u%gepsr,cb%mz*cb%mx)
call sysio_write('g0epsr_'//shot%sindex,a0_star_u0%gepsr,cb%mz*cb%mx)
call sysio_write('gmur_' //shot%sindex,   a_star_u%gmur ,cb%mz*cb%mx)
call sysio_write('g0mur_' //shot%sindex,a0_star_u0%gmur ,cb%mz*cb%mx)
! call sysio_write('gsgma_'//shot%sindex, a_star_u%gsgma  ,cb%mz*cb%mx)

        a_star_u%gepsr = a_star_u%gepsr + a0_star_u0%gepsr
        a_star_u%gmur  = a_star_u%gmur  + a0_star_u0%gmur
        a_star_u%gsgma = a_star_u%gsgma + a0_star_u0%gsgma

        call hud('----  Assemble  ----')
        call ppg%assemble(a_star_u)

        !produce final misfit
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
        ! call a0_star_u0%write
    endif

    !allreduce energy, gradient
    ! call mpi_allreduce(mpi_in_place, correlate_energy  , m%n          , mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, correlate_gradient, m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale by shotlist
    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)

    if(mpiworld%is_master) call sysio_write('correlate_gradient',correlate_gradient,m%n*ppg%ngrad)

    call mpiworld%barrier

    is_first_in=.false.



contains

    subroutine compute_L2_grad()

        character(:),allocatable :: s_m1, s_p1, fname_data_prefix

        type(t_suformat) ::    Ru_m1,    Ru_p1
        type(t_suformat) ::     d_m1,     d_p1
        type(t_suformat) :: W2res_m1, W2res_p1

        real,dimension(:,:),allocatable :: DRu, Dd, W2res
        
        !L2_grad: '0.5|| W(Du - d)||² => adjsrc = DᵀW²(Du - d)')
        !     [u2-u1]   [-1 1      ][u1]        [-1 -1         ][u1]   [-u1-u2]
        !     |u3-u1|   |-1 0 1    ||u2|        | 1  0 -1      ||u2|   | u1-u3|
        !Du = |u4-u2| = |  -1 0 1  ||u3|, Dᵀu = |    1  0 -1   ||u3| = | u2-u4|
        !     |u5-u3|   |    -1 0 1||u4|        |       1  0 -1||u4|   | u3-u5|
        !     [u5-u4]   [      -1 1][u5]        [          1  1][u5]   [ u4+u5]

        !residuals
        i_m1 = shot%index-1; if(i_m1==0           ) i_m1=1
        i_p1 = shot%index+1; if(i_p1==shls%nshot+1) i_p1=shls%nshot

        s_m1 = num2str(i_m1,'(i0.4)')
        s_p1 = num2str(i_p1,'(i0.4)')

        fname_data_prefix=setup%get_str('FILE_DATA_PREFIX')

        call hud('Will read Shot#'//s_m1//' & '//s_p1,mpiworld%iproc)

        call Ru_m1%read(dir_out//'Ru_Shot'//s_m1//'.su',shot%sindex)
        call Ru_p1%read(dir_out//'Ru_Shot'//s_p1//'.su',shot%sindex)

        call d_m1%read(fname_data_prefix//s_m1//'.su')
        call d_p1%read(fname_data_prefix//s_p1//'.su')

        DRu = Ru_p1%trs -Ru_m1%trs
        Dd  =  d_p1%trs - d_m1%trs

        fobj%reflection = fobj%reflection &
            + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, Dd-DRu, shot%dt)

        call alloc(W2res,shot%nt,shot%nrcv)
        call kernel_L2sq(W2res)
        call shot%write('W2res_',W2res)

        call mpiworld%barrier


        !adjoint sources
        call W2res_m1%read(dir_out//'W2res_Shot'//s_m1//'.su')
        call W2res_p1%read(dir_out//'W2res_Shot'//s_p1//'.su')

        if     (shot%index==1) then
            shot%dadj =-W2res_m1%trs - W2res_p1%trs
        else if(shot%index==shls%nshot) then
            shot%dadj = W2res_m1%trs + W2res_p1%trs
        else
            shot%dadj = W2res_m1%trs - W2res_p1%trs
        endif

    end subroutine
    

    subroutine compute_L2_slope()
    use m_median2d
        real,dimension(:,:),allocatable :: Ru, d, DRu, Dd, pu, pd, W2res, pua, dadj

        call alloc(Ru,shot%nt,shls%nshot)
        call alloc( d,shot%nt,shls%nshot)
    
        call MPI_Gather(shot%dsyn, shot%nt, MPI_REAL, Ru, shot%nt, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(shot%dobs, shot%nt, MPI_REAL,  d, shot%nt, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

        if(mpiworld%is_master) then
            DRu = grad(Ru,2)
            Dd  = grad(d, 2)

            dtr = setup%get_real('DTR',o_mandatory=1)
            
            pu = median2d(slope(DRu,shot%dt,dtr),20) !20 is half width
            pd = median2d(slope(Dd ,shot%dt,dtr),20)

            !L2_slope: '0.5|| W(pDu - pDd)||² => adjsrc = Dᵀ[∂p/∂Du]W²Δp')
            fobj%reflection = fobj%reflection &
                + L2sq(0.5, shot%nrcv*shot%nt, wei%weight, pd-pu, shot%dt)

            call alloc(W2res,shot%nt,shls%nshot)
            call kernel_L2sq(W2res)
            call sysio_write('W2res',W2res,size(W2res))

            pua= part_p_part_u(W2res,DRu,shot%dt,dtr)
            call sysio_write('pua',pua,size(pua))

            !
            !     [u2-u1]   [-1 1      ][u1]        [-1 -1         ][u1]   [-u1-u2]
            !     |u3-u1|   |-1 0 1    ||u2|        | 1  0 -1      ||u2|   | u1-u3|
            !Du = |u4-u2| = |  -1 0 1  ||u3|, Dᵀu = |    1  0 -1   ||u3| = | u2-u4|
            !     |u5-u3|   |    -1 0 1||u4|        |       1  0 -1||u4|   | u3-u5|
            !     [u5-u4]   [      -1 1][u5]        [          1  1][u5]   [ u4+u5]

            call alloc(dadj,shot%nt,shls%nshot)
                dadj(:,1)  = -pua(:,1)    -pua(:,2)
            do ix=2,shls%nshot-1
                dadj(:,ix) =  pua(:,ix-1) -pua(:,ix+1)
            enddo
                dadj(:,shls%nshot) = pua(:,shls%nshot-1) +pua(:,shls%nshot)

        endif

        if(.not.mpiworld%is_master) allocate(dadj(1,1)) !this is nonsense..

        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call MPI_Scatter(dadj, shot%nt, MPI_REAL, shot%dadj, shot%nt, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

        deallocate(dadj)

    end subroutine


    function grad(u,iaxis) result(g)
        real,dimension(:,:)             :: u
        real,dimension(:,:),allocatable :: g
        
        n1=size(u,1)
        n2=size(u,2)
        call alloc(g,n1,n2)

        if(iaxis==1) then
                g(1,:) = u(2,:)   -u(1,:)
            do i=2,n2-1
                g(i,:) = u(i+1,:) -u(i-1,:)
            enddo
                g(n2,:)= u(n2,:)  -u(n2-1,:)
        endif

        if(iaxis==2) then
                g(:,1) = u(:,2)   -u(:,1)
            do i=2,n2-1
                g(:,i) = u(:,i+1) -u(:,i-1)
            enddo
                g(:,n2)= u(:,n2)  -u(:,n2-1)
        endif

    end function

    function slope(u,d1,d2) result(p)
        real,dimension(:,:)             :: u
        real,dimension(:,:),allocatable :: p

        n1=size(u,1)
        n2=size(u,2)
        call alloc(p,n1,n2)

        p = -grad(u,2)/grad(u,1)/d2*d1

    end function

    function part_p_part_u(Dp,u,d1,d2) result(pua)
        real,dimension(:,:)             :: Dp, u
        real,dimension(:,:),allocatable :: pua

        real,dimension(:,:),allocatable :: ut,ux, Dp_ut,ux_ut

        ! ∂p/∂u = ∂ₓ(Δp/∂ₜu) -∂ₜ(Δp*∂ₓu/∂ₜu/∂ₜu) 
        ut = grad(u,1)/d1
        water = 1e-5*maxval(abs(ut))
        Dp_ut = ut*Dp / (ut*ut + water)
        ux = grad(u,2)/d2
        ux_ut = ut*ux / (ut*ut + water)

        n1=size(u,1)
        n2=size(u,2)
        call alloc(pua,n1,n2)
        pua = grad( Dp_ut ,2)/d2 - grad( Dp_ut*ux_ut ,1)/d1

    end function


end subroutine
