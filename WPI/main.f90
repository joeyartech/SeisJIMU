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
             '       WELCOME TO SeisJIMU WFI        '//s_NL// &
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
    call chp_shls%init('FWI_shotlist_gradient',oif_fuse=.true.)
    if(.not.shls%is_registered(chp_shls,'sampled_shots')) then
        call shls%sample
        call shls%register(chp_shls,'sampled_shots')
    endif
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
    call chp_qp%init('FWI_querypoint_gradient')
    if(.not.qp0%is_registered(chp_qp)) then
        call fobj%eval(qp0,oif_update_m=.false.)
        call qp0%register(chp_qp)
    endif

    call sysio_write('qp0%g',qp0%g,size(qp0%g))
    call sysio_write('qp0%pg',qp0%pg,size(qp0%pg))

    !scale problem by linesearcher
    call ls%init
    call ls%scale(qp0)

    call hud('qp0%f, ||g||_L2 = '//num2str(qp0%f)//', '//num2str(norm2(qp0%g)))

    !if just estimate the wavelet or compute the gradient then this is it.
    ! if(setup%get_str('JOB')=='gradient') then
        call mpiworld%final
        stop
    ! endif

    ! !optimization
    ! call optimizer_init(qp0)
    ! call optimizer_loop
    
    ! call mpiworld%final
    
end

subroutine modeling_imaging
use mpi
use m_System
use m_Modeling
use m_weighter
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse

    type(t_field) :: fld_u, fld_a, fld_du, fld_da
    real,dimension(:,:),allocatable :: Wdres
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff

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
        call hud('build images')
        call hud('---------------------------------')
        !forward modeling on u
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite       
        call ppg%forward(fld_u);                    call fld_u%acquire
        call shot%write('Ru_imag_',shot%dsyn)

call wei%update('IMAGING'); !wei%weight(1:200,:)=0.
call alloc(Wdres,shot%nt,shot%nrcv)
Wdres = (shot%dsyn-shot%dobs)*wei%weight
shot%dadj = -wei%weight*Wdres*m%cell_volume/shot%dt
shot%dadj = shot%dadj/m%ref_kpa
call shot%write('dadj_',shot%dadj)

        !adjoint modeling on a
        !A^H a = R^H*Delta_d
        !C_field=0.5*(u_star_a)^2
        !grad = u \star a
        call ppg%init_field(fld_a,name='fld_a_imag');    call fld_a%ignite(ois_adjoint=.true.)
        call ppg%adjoint(fld_a,fld_u,oif_compute_imag=.true.)
        call sysio_write('I',cb%imag,size(cb%imag))
        call hud('---------------------------------')

        call cb%project_back

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    !collect global objective function value
    ! call mpi_allreduce(mpi_in_place, fobj%dnorms, fobj%n_dnorms, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    ! call fobj%print_dnorms('Unscaled stacked','')

    ! if(either(oif_gradient,.true.,present(oif_gradient))) then
        !collect global gradient
        ! call mpi_allreduce(mpi_in_place, m%gradient,  m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        ! if(ppg%if_compute_engy) then
        !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        ! endif
    ! endif

    call mpiworld%barrier

end subroutine


subroutine modeling_gradient
use mpi
use m_System
use m_Modeling
use m_weighter
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse
use m_resampler

    type(t_field) :: fld_u, fld_a, fld_du, fld_da
    real,dimension(:,:,:),allocatable :: tmpgrad
    real,dimension(:,:),allocatable :: tmpdsyn
    ! character(:),allocatable :: update_wavelet
    ! real,dimension(:,:),allocatable :: Rmu, Rdiff

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

        call hud('--- classical FWI Sensitivity Kernels ---')
        call alloc(tmpgrad,cb%nz,cb%nx,cb%ny)
        !forward modeling on u
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite       
        call ppg%forward(fld_u);                    call fld_u%acquire

        !dnorm
        !C_data=Ru-d
        call shot%write('Ru_',shot%dsyn)
        call wei%update('GRADIENT')
        call fobj%compute_dnorms

        !adjoint modeling on a
        !A^H a = R^H*Delta_d
        !C_field=0.5*(u_star_a)^2
        !grad = u \star a
        call ppg%init_field(fld_a,name='fld_a');    call fld_a%ignite(ois_adjoint=.true.)
        call ppg%adjoint(fld_a,fld_u,oif_compute_grad=.true.)
        call sysio_write('u_star_a',cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
        tmpgrad=cb%grad(:,:,:,2)
        
        call hud('--- load prior image ---')
        call alloc(cb%imag,cb%mz,cb%mx,cb%my,1)
        call sysio_read(setup%get_file('IMAGE'),cb%imag,size(cb%imag))
        call sysio_write('loaded_I',cb%imag,size(cb%imag))
        W=1.
        C_imag=0.5*sum(W*cb%imag)
        call hud('---------------------------------')
! pause
        call hud('--- extd left-side Rabbit Ears ---')
        !model u and delta_u
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
        call ppg%init_field(fld_du,name='fld_du')
        call ppg%forward_ext(fld_du,fld_u,W*W*cb%imag); call fld_du%acquire
        call shot%write('Rdu_',shot%dsyn)
        !
        !A delta_u = W^2*Iu
        !grad = delta_u \star a
        call ppg%init_field(fld_a,name='fld_a') !adjsrc for a should not be initialized
        call ppg%adjoint_ext1(fld_a,fld_du,fld_u,W*W*cb%imag,oif_compute_grad=.true.) !,oif_compute_imag=.true.)
        ! call sysio_write('image_du_star_a',cb%imag,size(cb%imag))
        call sysio_write('du_star_a',cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
        tmpgrad=tmpgrad+cb%grad(:,:,:,2)
        call hud('---------------------------------')
! pause
        call hud('--- extd right-side Rabbit Ears ---')
        !re-model u
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        !
        !A^H delta_a = W^2*Ia
        !grad3 = u \star delta_a
        call ppg%init_field(fld_a, name='fld_a') !adjsrc for a should not be initialized
        call ppg%init_field(fld_da,name='fld_da')
        call ppg%adjoint_ext2(fld_da,fld_a,W*W*cb%imag,fld_u,oif_record_adjseismo=.true.,oif_compute_grad=.true.)
        call suformat_write('da_adjseismo',fld_da%seismo,ppg%nt,1,ppg%dt)
        call sysio_write('u_star_da',cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
        tmpgrad=tmpgrad+cb%grad(:,:,:,2)
        call hud('---------------------------------')

        call hud('--- high-order term ---')
        call ppg%init_field(fld_u,name='fld_u');    call fld_u%ignite
        call ppg%init_field(fld_du,name='fld_du')
        call ppg%forward_ext(fld_du,fld_u,W*W*cb%imag); call fld_du%acquire
        !adjoint modeling on new a with adjsource=R*delta_u
        !A^H a = R^H*R*delta_u
        !grad4 = u \star a
!!manual resampling
call alloc(tmpdsyn,ppg%nt,shot%nrcv)
call resampler(shot%dsyn,tmpdsyn, &
            ntr=shot%nrcv,&
            din=shot%dt,nin=shot%nt,&
            dout=ppg%dt,nout=ppg%nt)

        call ppg%init_field(fld_a,name='fld_a');    call fld_a%ignite(ois_adjoint=.true.,o_wavelet=tmpdsyn)
        call ppg%adjoint_ext1(fld_a,fld_du,fld_u,W*W*cb%imag,oif_compute_grad=.true.)
        call sysio_write('u_star_adj(du)',cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))
        tmpgrad=tmpgrad+cb%grad(:,:,:,2)
        call hud('---------------------------------')

        cb%grad(:,:,:,2)=tmpgrad
        call sysio_write('grad',cb%grad(:,:,:,2),size(cb%grad(:,:,:,2)))

        call cb%project_back

        fobj%dnorms=fobj%dnorms+C_imag

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    !collect global objective function value
    call mpi_allreduce(mpi_in_place, fobj%dnorms, fobj%n_dnorms, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call fobj%print_dnorms('Unscaled stacked','')

    ! if(either(oif_gradient,.true.,present(oif_gradient))) then
        !collect global gradient
        call mpi_allreduce(mpi_in_place, m%gradient,  m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        ! if(ppg%if_compute_engy) then
        !     call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        ! endif
    ! endif

    call mpiworld%barrier

end subroutine

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
