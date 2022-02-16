program main
use m_System
use m_Modeling
use m_Kernel
use m_Optimization

    type(t_checkpoint) :: chp_shls, chp_qp
    type(t_querypoint),target :: qp0
    character(:),allocatable :: job

    !mpiworld lives in t_mpienv
    call mpiworld%init(name='MPIWorld')

    call hud('======================================'//s_NL// &
             '       WELCOME TO SeisJIMU WPI        '//s_NL// &
             '======================================')

    call setup%init
    call sysio_init
    
    if(.not. setup%exist) then
        if(mpiworld%is_master) then
            call hud('No input setup file given. Print manual..')
            call print_manual
            call ppg%print_info
        endif
        call mpiworld%final
        stop
    endif

    !print propagator info
    call ppg%print_info

    ! !estimate required memory
    ! call m%estim_RAM
    ! call cb%estim_RAM
    ! call sfield%estim_RAM
    ! call rfield%estim_RAM
    
    !checkpoint
    call checkpoint_init

    !model
    call m%init
    call m%read
    call ppg%check_model

    !shotlist
    call shls%read_from_data
    call shls%build
    ! call chp_shls%init('FWI_shotlist_gradient',oif_fuse=.true.)
    ! if(.not.shls%is_registered(chp_shls,'sampled_shots')) then
        call shls%sample
    !     call shls%register(chp_shls,'sampled_shots')
    ! endif
    call shls%assign

    !if preconditioner needs energy terms
    if(index(setup%get_str('PRECONDITIONING','PRECO'),'energy')>0) then
        ppg%if_compute_engy=.true.
    endif

    !parametrizer
    call param%init

    !initial (model) parameters as querypoint
    call qp0%init('qp0')

    !objective function and gradient
    call fobj%init
    ! call chp_qp%init('FWI_querypoint_gradient')
    ! if(.not.qp0%is_registered(chp_qp)) then
        call fobj%eval(qp0,oif_update_m=.false.)
    !     call qp0%register(chp_qp)
    ! endif

    if(.not. qp0%is_fitting_data) then
        call hud('Negate the sign of qp0 to ensure fitting the data')
        call qp0%set_negative
    endif

    call sysio_write('qp0%g',qp0%g,size(qp0%g))
    call sysio_write('qp0%pg',qp0%pg,size(qp0%pg))

    !scale problem by linesearcher
    call ls%init
    call ls%scale(qp0)

    call hud('qp0%f, ||g||_L2 = '//num2str(qp0%f)//', '//num2str(norm2(qp0%g)))

    !if just estimate the wavelet or compute the gradient then this is it.
    if(setup%get_str('JOB')=='gradient') then
        call mpiworld%final
        stop
    endif

    !optimization
    call optimizer_init(qp0)
    call optimizer_loop
    
    call mpiworld%final
    
end

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

    ! if(setup%get_str('JOB')=='skip_imaging') then
    !     return
    ! endif
    
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
        
        tmp = L2norm_sq(size(wei%weight),wei%weight,shot%dsyn-shot%dobs,shot%dt,0.5/m%ref_kpa*m%cell_volume/shot%dt)                     ![Pa]^2    [s]   [1/Pa m^3/s]=[Nm]
        
        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call nabla_L2norm_sq(shot%dadj)
        shot%dadj=-shot%dadj
        call shot%write('Imag_dadj_',shot%dadj)

        !adjoint modeling on a
        !A^H a = R^H*Delta_d
        !C_field=0.5*(u_star_a)^2
        !grad = u \star a
        call ppg%init_field(fld_a,name='Imag_fld_a',ois_adjoint=.true.); call fld_a%ignite
        call ppg%adjoint(fld_a,fld_u,oif_compute_imag=.true.)
        call hud('---------------------------------')

        call cb%project_back

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    !collect
    call mpi_allreduce(mpi_in_place, m%image,  size(m%image), mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)

    call mpiworld%barrier

    if(mpiworld%is_master) then
        call sysio_write('unmask_I',m%image(:,:,:,1),size(m%image(:,:,:,1)))

        !apply freeze zone
        where(m%is_freeze_zone) m%image(:,:,:,1)=0.
        call sysio_write('I',m%image(:,:,:,1),size(m%image(:,:,:,1)))
    endif

    if(setup%get_str('JOB')=='imaging') then
        call mpiworld%final
        stop
    endif

end subroutine


subroutine modeling_gradient(is_fitting_data)
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

    character(:),allocatable :: corrs
    type(t_field) :: fld_u, fld_du, fld_Adj_du, fld_a, fld_da
    ! real,dimension(:,:),allocatable :: Wdres
    ! real,dimension(:,:,:),allocatable :: WImag, 
    real,dimension(:),allocatable :: tmp
    real,dimension(:,:,:),allocatable :: Imag_as_adjsrc
    real,dimension(:,:),allocatable :: tmpdsyn
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff


    if(setup%get_file('IMAGE')/='') then
        call alloc(m%image,m%nz,m%nx,m%ny,1)
        call sysio_read(setup%get_file('IMAGE'),m%image,size(m%image))
        call sysio_write('loaded_I',m%image,size(m%image))    
        
        !get shot%dt
        call shot%init(shls%yield(1))
        call shot%read_from_data
        call shot%set_var_time
    else
        call hud('Will do imaging')
        call modeling_imaging
    endif

    !compute fobjective
    call iwei%update

    !L1 & L2 norms of image
    !fobj%dnorms(1) = integral(size(iwei%weight),abs(iwei%weight),abs(m%image),1, m%cell_volume,1./m%ref_kpa/shot%dt)   ![Pa^2 s] [m^3] [1/Pa 1/s]=[Nm]
    fobj%dnorms(1) = L1norm   (size(iwei%weight),iwei%weight,m%image,m%cell_volume,1./m%ref_kpa/shot%dt)   !                                [Pa^2 s] [m^3]        [1/Pa 1/s]=[Nm]

    !fobj%dnorms(2) = integral(size(iwei%weight),iwei%weight**2,m%image,2,        m%cell_volume,0.5/m%ref_kpa**3/shot%dt**2)  ![Pa^2 s]^2 [m^3] [Pa^-3 s^-2]=[Nm]
    fobj%dnorms(2) = L2norm_sq(size(iwei%weight),iwei%weight,m%image,m%cell_volume,0.5/m%ref_kpa**3/shot%dt**2)                            ![Pa^2 s]^2 [m^3]      [Pa^-3 s^-2]=[Nm]
    
    deallocate(m%image)

    !Use L2 norm for adjsrc
    call alloc(Imag_as_adjsrc,m%nz,m%nx,m%ny)
    !Imag_as_adjsrc =reshape( -nabla_integral(), [m%nz,m%nx,m%ny])
    call nabla_L2norm_sq(Imag_as_adjsrc)
    Imag_as_adjsrc = -Imag_as_adjsrc

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

        call hud('-------- FWI SK (u star a) --------')
        !forward modeling on u
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite       
        call ppg%forward(fld_u);                    call fld_u%acquire

        !dnorm
        !C_data=Ru-d
        call shot%write('Ru_',shot%dsyn)

        call wei%update!('_4FWI')

        !compute FWI data misfit and adjsrc
        fobj%FWI_misfit = fobj%FWI_misfit &
            + L2norm_sq(size(wei%weight),wei%weight,shot%dsyn-shot%dobs,shot%dt,0.5/m%ref_kpa*m%cell_volume/shot%dt)                     ![Pa]^2    [s]   [1/Pa m^3/s]=[Nm]

        call alloc(shot%dadj,shot%nt,shot%nrcv)
        call nabla_L2norm_sq(shot%dadj); shot%dadj=-shot%dadj
        
        call shot%write('FWI_dadj_',shot%dadj)

        !adjoint modeling on a
        !A^H a = R^H*Delta_d
        !C_field=0.5*(u_star_a)^2
        !grad = u \star a
        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)

        call fld_a%ignite

        call ppg%adjoint(fld_a,fld_u,oif_compute_grad=.true.)

        cb%corr(:,:,:,1)=cb%grad(:,:,:,2) !=gkpa, propto gvp under Vp-Rho
        call hud('---------------------------------')


        if(index(corrs,'RE')>0) then

            call hud('--- extd right-side Rabbit Ears (du star a) ---')
            call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
            call ppg%init_field(fld_du,name='fld_du')
            call ppg%forward_scattering(fld_du,fld_u,Imag_as_adjsrc); call fld_du%acquire
            call shot%write('Rdu_',shot%dsyn)
            !
            !A delta_u = W_I^2*Iu
            !grad = delta_u \star a
            call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.) !adjsrc for a should not be initialized
            call ppg%adjoint_du_star_a(fld_a,fld_du,fld_u,Imag_as_adjsrc,oif_compute_grad=.true.) !,oif_compute_imag=.true.)

            cb%corr(:,:,:,2)=cb%grad(:,:,:,2) !=gkpa, propto gvp under Vp-Rho
            call hud('---------------------------------')
    ! pause

            call hud('--- extd left-side Rabbit Ears (u star da) ---')
            !re-model u
            call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
            call ppg%forward(fld_u)
            !
            !A^H delta_a = W_I^2*Ia
            !grad3 = u \star delta_a
            call ppg%init_field(fld_a, name='fld_a', ois_adjoint=.true.) !adjsrc for a should not be initialized
            call ppg%init_field(fld_da,name='fld_da',ois_adjoint=.true.)
            call ppg%adjoint_scattering_u_star_da(fld_da,fld_a,Imag_as_adjsrc,fld_u,oif_compute_grad=.true.)

            cb%corr(:,:,:,3)=cb%grad(:,:,:,2) !=gkpa, propto gvp under Vp-Rho
            call hud('---------------------------------')

        endif

        if(index(corrs,'2ndMI')>0) then

            call hud('-------- du_star_da --- (just for curiosity)')
            call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
            call ppg%init_field(fld_du,name='fld_du')
            call ppg%forward_scattering(fld_du,fld_u,Imag_as_adjsrc)

            call ppg%init_field(fld_a, name='fld_a', ois_adjoint=.true.)
            call ppg%init_field(fld_da,name='fld_da',ois_adjoint=.true.)
            call ppg%adjoint_scattering_du_star_da(fld_da,fld_a,fld_du,fld_u,Imag_as_adjsrc,oif_compute_grad=.true.)

            cb%corr(:,:,:,4)=cb%grad(:,:,:,2) !=gkpa, propto gvp under Vp-Rho
            call hud('---------------------------------')

        endif


        if(index(corrs,'DR')>0) then

            call hud('--- demig-remig (u_star_Adj(du)) ---')
            call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
            call ppg%init_field(fld_du,name='fld_du')
            call ppg%forward_scattering(fld_du,fld_u,Imag_as_adjsrc); call fld_du%acquire !shot%dsyn = R*fld_du
            !adjoint modeling on new a with adjsource=R*delta_u
            !A^H a = R^H*R*delta_u
            !grad4 = u \star a

            !!manual weighting and resampling
            shot%dsyn=shot%dsyn*wei%weight
            if(mpiworld%is_master) call shot%write('Wei_Rdu_',shot%dsyn)
            !!cheating: shot%dsyn=shot%dobs*wei%weight
            call alloc(tmpdsyn,ppg%nt,shot%nrcv)
            call resampler(shot%dsyn,tmpdsyn, &
                        ntr=shot%nrcv,&
                        din=shot%dt,nin=shot%nt,&
                        dout=ppg%dt,nout=ppg%nt)

            call ppg%init_field(fld_Adj_du,name='fld_Adj_du',ois_adjoint=.true.)
            call fld_Adj_du%ignite(o_wavelet=tmpdsyn)
            call ppg%adjoint(fld_Adj_du,fld_u,oif_compute_grad=.true.)

            cb%corr(:,:,:,5)=cb%grad(:,:,:,2) !=gkpa, propto gvp under Vp-Rho
            call hud('---------------------------------')

        endif
        
        call cb%project_back

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')


    !collect global objective function value
    call mpi_allreduce(mpi_in_place, fobj%dnorms, fobj%n_dnorms, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call fobj%print_dnorms('Stacked but not yet linesearch-scaled','')

    call mpi_allreduce(mpi_in_place, fobj%FWI_misfit, 1, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call hud('Stacked FWI_misfit '//num2str(fobj%FWI_misfit))

    !collect global correlations
    call mpi_allreduce(mpi_in_place, m%correlate,  m%n*5, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)


    if(mpiworld%is_master) then
            call sysio_write('u_star_a',      m%correlate(:,:,:,1),m%n)
        if(index(corrs,'RE')>0) then
            call sysio_write('du_star_a',     m%correlate(:,:,:,2),m%n)
            call sysio_write('u_star_da',     m%correlate(:,:,:,3),m%n)
        endif
        if(index(corrs,'2ndMI')>0) then
            call sysio_write('du_star_da',    m%correlate(:,:,:,4),m%n)
        endif
        if(index(corrs,'DR')>0) then
            call sysio_write('u_star_Adj(du)',m%correlate(:,:,:,5),m%n)
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

    deallocate(m%correlate)
    
    ! if(either(oif_gradient,.true.,present(oif_gradient))) then
        !collect global gradient
        !call mpi_allreduce(mpi_in_place, m%gradient,  m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        ! if(ppg%if_compute_engy) then
        !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        ! endif
    ! endif

    call mpiworld%barrier

end subroutine


! subroutine local_write(fname,grad,fullgrad,mm)
! use mpi
! use m_System

!     character(*) :: fname
!     real,dimension(mm) :: grad,fullgrad

    
!     call sysio_write(fname,grad,mm,o_mode='stack')
!     fullgrad=fullgrad+grad
    
! end subroutine

subroutine print_manual

    write(*,'(a)') ""
    write(*,'(a)') "----------------------------"
    write(*,'(a)') "To launch the program, do:"
    write(*,'(a)') "----------------------------"
    write(*,'(a)') ""
    write(*,'(a)') "bash $ export OMP_NUM_THREADS=1 #(use source in csh like shell)"
    write(*,'(a)') "bash $ mpirun -np $np FWI setup.in"
    write(*,'(a)') ""
    write(*,'(a)') "----------------------------"
    write(*,'(a)') "Mandatory items in setup.in:"
    write(*,'(a)') "----------------------------"
    write(*,'(a)') ""
    write(*,'(a)') "FILE_MODELS             'model'            #Prefix of input model name"
    write(*,'(a)') "MODEL_DIMENSION         '201  201  201'    #nz nx ny"
    write(*,'(a)') "MODEL_SPACING           ' 10   10   10'    #dz dx dy"
    write(*,'(a)') "MODEL_ORIGIN            '  0    0    0'    #oz ox oy"
    write(*,'(a)') "                        #all models should have same dimension, spacing and origin"
    write(*,'(a)') ""
    write(*,'(a)') "NBOUNDARYLAYER          60                 #Number of convolutional Perfect-Match Layers added to edges of computebox"
    write(*,'(a)') ""
    write(*,'(a)') "PEAK_FREQUENCY          7                  #Peak frequency of source wavelet"
    write(*,'(a)') ""
    write(*,'(a)') "FILE_DATA               'obs'              #Prefix of input obs data name"
    write(*,'(a)') "                                           #SeisJIMU will look for SU data file named as obs????.su,"
    write(*,'(a)') "                                           #where ???? represents 4 digits"
    write(*,'(a)') "                                           #Provide SHOTNO if only invert a subset of data (see below)"
    write(*,'(a)') ""
    write(*,'(a)') "---------------------------"
    write(*,'(a)') "Optional items in setup.in:"
    write(*,'(a)') "---------------------------"
    write(*,'(a)') ""
    write(*,'(a)') "APERTURE                '-99999 99999 -99999 99999'    #Define size of computebox"
    write(*,'(a)') "                                                       #='0 0 0 0' computebox bounded by farthest source & receiver"
    write(*,'(a)') "                                                       # (i.e. leftmost & rightmost device defines fx & lx, similar for fy & ly)"
    write(*,'(a)') "                                                       # (smallest possible computebox otherwise some devices will be truncated)"
    write(*,'(a)') "                                                       #='-99999 99999 -99999 99999' use whole model as computebox"
    write(*,'(a)') "                                                       # (largest possible computebox otherwise it will be larger than the model)"
    write(*,'(a)') "                                                       #Well-chosen values between above two extremes can give a good balance"
    write(*,'(a)') "                                                       # between edge-effect and computational cost"
    write(*,'(a)') "IF_ISOTROPIC            F                   #If T, omit model_eps model_del"
    write(*,'(a)') "IF_TOPO_FROM_VS         F                   #If T, decide topography from model_vs when model_vs exists and model_topo doesn't exist"
    write(*,'(a)') "IF_FREESURFACE          T                   #Free surface condition (T) or absorbing surface condition (F)"
    write(*,'(a)') ""
    write(*,'(a)') "SHOTNO                  'fshot:dshot:lshot' #Will invert only shots with index from fshot to lshot, with increments dshot"
    write(*,'(a)') "                                            #If not given, SeisJIMU will look for all data with sequential index (start from 0001)"
    write(*,'(a)') "                                            #the prefix of the filename has been defined in FILE_DATA"
    write(*,'(a)') ""
    write(*,'(a)') "FILE_WAVELET            'fricker'           #If not given, SeisJIMU will use the wavelet generated from WAVELET_TYPE"
    write(*,'(a)') "WAVELET_TYPE            'sinexp'            #Damped sine function: sin(2pi*fpeak*t)*exp(-3.33*fpeak*t)"
    write(*,'(a)') "                                            #s.t. the wavelet has ~2 peaks before dying out, followed by lowpass filtering."
    write(*,'(a)') "                                            #If not given, use Ricker wavelet:"
    write(*,'(a)') "                                            #x(t) = -(pi*fpeak*(t-t0))^2, w(t) = (1+2*x)*exp(x)"
    write(*,'(a)') "                                            #t0 = 1/fpeak by default"
    write(*,'(a)') "RICKER_DELAYTIME        '1/PEAK_FREQUENCY'  #Change t0 in Ricker wavelet"
    write(*,'(a)') "SCALE_WAVELET           'no'                #If scale source wavelet when read in"
    write(*,'(a)') "                                            #='no'       no scale"
    write(*,'(a)') "                                            #='by dxdt'  multiply by dt/dx/dy/dz s.t. the amplitude is independent on discretization"
    write(*,'(a)') "                                            #=number     user-given scaler on wavelet"
    write(*,'(a)') "UPDATE_WAVELET          'per shot'          #If update source wavelet after each forward modeling (using matching filter)"
    write(*,'(a)') "                                            #='per shot' update each shot's wavelet independently"
    write(*,'(a)') "                                            #='stack'    stack each shot's wavelet and use the same stacked wavelet for all shots"
    write(*,'(a)') "                                            #='no'       no update"
    write(*,'(a)') ""
    write(*,'(a)') "IF_HICKS                T    #Use Hicks method to interpolate source & receiver positions when not colocate with grid points."
    write(*,'(a)') "                             #If F, sources & receivers will be positioned to nearest grid points."
    write(*,'(a)') ""
    write(*,'(a)') "IF_SNAPSHOT             F                   #If T, save field snapshot to disk"
    write(*,'(a)') "SNAPSHOT_DELTA_IT       50                  #Save snapshot every SNAPSHOT_DELTA_IT time step"
    write(*,'(a)') ""
    write(*,'(a)') "FILE_WEIGHT_POLYGON     ''                  #ASCII file giving polygon-like weighting applied to data residuals in time-offset domain"
    write(*,'(a)') "                                            #partition the gather by polygons with different weights based on interpretation of phases"
    write(*,'(a)') "                                            #(e.g. if you want to separate reflections from diving waves) (modified from SU program sumute)"
    write(*,'(a)') "                                            #modified from SU program sumute. If not given, weight operator is unchanged."
    write(*,'(a)') "                                            #Example of file content (e.g. muting reflections below direct waves):"
    write(*,'(a)') "                                            #  #define a polygon that separates reflections from direct waves"
    write(*,'(a)') "                                            #  0    1000    #offsets of vertices (suppose max offset=1000m)"
    write(*,'(a)') "                                            #  0.2  1.7     #time instants of vertices (bottom of direct waves)"
    write(*,'(a)') "                                            #  0    0.1     #weight=0 below this polygon (0.1s transition)"
    write(*,'(a)') "                                            #  #more polygons can be added in the same way as above"
    write(*,'(a)') "                                            #  END          #end of file"
    write(*,'(a)') "FILE_WEIGHT_TABLE       ''                  #ASCII file giving table-like weighting applied to data residual in time-offset domain"
    write(*,'(a)') "                                            #If not given, weight operator is unchanged."
    write(*,'(a)') "                                            #Example of file content (e.g. increasing weight with absolute offset and time):"
    write(*,'(a)') "                                            #      100  500  1000  5000    #abs offset values"
    write(*,'(a)') "                                            #  #-----------------------"
    write(*,'(a)') "                                            #  0 |   1    5    10    50    #time instant (0s) and weights at 100m, 500m, 1000m & 5000m"
    write(*,'(a)') "                                            #  4 |  11   15    20    60    #time instant (0s) and weights at same positions"
    write(*,'(a)') "                                            #  #more time instants and associated weights can be added in the same way as above"
    write(*,'(a)') "                                            #  #weight at un-specified positions will be interpolated from this table"
    write(*,'(a)') "                                            #  END          #end of file"
    write(*,'(a)') ""
    write(*,'(a)') "IF_SMOOTHING            T                   #If smooth gradient by Laplacian operator"
    write(*,'(a)') "SMOOTHING_ADDMIRROR     dz                  #Smoothing could cause leakage into the mask zone (e.g. sea column)."
    write(*,'(a)') "                                            #If this is significant, increase this value to a multiple of dz (spacing in z dir)"
    write(*,'(a)') "SMOOTH_GRADIENT_WAVELENGTH_FRACTION   '1 1 1'         #Smoothing operator length versus wavelength in z,x,y dirs"
    write(*,'(a)') "SMOOTH_GRADIENT_PRESERVE_MAGNITUDE    'nopreserve'    #If restore raw gradient magnitude after smoothing"
    write(*,'(a)') "                                                      #='L1norm' restore magnitude defined with L1 norm"
    write(*,'(a)') "                                                      #='L2norm' restore magnitude defined with L2 norm"
    write(*,'(a)') "                                                      #='scale1' restore point-wise magnitude"
    write(*,'(a)') ""
    write(*,'(a)') "ZPOWER                  1                   #Depth preconditioner to gradient: multiply each column by z^ZPOWER"
    write(*,'(a)') ""
    write(*,'(a)') "EMPIRICAL_LAW           'gardner'           #If specified, use interparameter empirical law when converting gradient (parameterization)"
    write(*,'(a)') "ACTIVE_PARAMETER        'vp1500:3400'       #Specify which parameter to be inferred, followed by range of allowable values"
    write(*,'(a)') ""
    write(*,'(a)') "JOB                     'optimization'      #='optimization' continues to FWI iterations"
    write(*,'(a)') "                                            #otherwise stop after the first gradient is computed"
    write(*,'(a)') ""
    write(*,'(a)') "MIN_UPDATE              1e-8                #Minimal parameter update allowed, otherwise stop FWI iteration"
    write(*,'(a)') "MAX_ITERATE             30                  #Maximum number of allowed FWI iterations"
    write(*,'(a)') "MAX_MODELING            60                  #Maximum number of allowed forward modeling"
    write(*,'(a)') "IF_REINITIALIZE_ALPHA   T/F                 #If reset alpha=1 at each new line search procedure (Default=T for l-BFGS, Default=F for NLCG)"
    write(*,'(a)') "NPAIRS                  5                   #For l-BFGS optimzer, number of saved gradients in history"
    write(*,'(a)') "                                            #(older than NPAIRS gradients will be flushed out)"
    write(*,'(a)') ""
    write(*,'(a)') "---------------------------"
    write(*,'(a)') "Notes:"
    write(*,'(a)') "---------------------------"
    write(*,'(a)') ""
    write(*,'(a)') "Acquisition geometry, source & receiver positions, time step (nt) and time interval (dt) are read directly from obs input data,"
    write(*,'(a)') "which have to be in SU format (other formats or read method can be developed in future)."
    write(*,'(a)') "You can use SU program suaddhead, or add_suheader in ./Tools, to convert data into SU format"
    write(*,'(a)') ""

end subroutine
