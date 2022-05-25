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
             '       WELCOME TO SeisJIMU FWI        '//s_NL// &
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
    ! call fld_u%estim_RAM
    ! call fld_a%estim_RAM
    
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

    if(index(param%info,'pseudotime')>0) then
        call sysio_write('m%gradient',m%gradient,size(m%gradient))
    endif

    !scale problem by linesearcher
    call ls%init
    call ls%scale(qp0)
    
    call hud('qp0%f, ║g║₁ = '//num2str(qp0%f)//', '//num2str(sum(abs(qp0%g))))

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

subroutine modeling_gradient(is_fitting_data)!(oif_gradient)
use mpi
use m_System
use m_Modeling
use m_weighter
use m_fobjective
use m_matchfilter
use m_smoother_laplacian_sparse

    !logical,optional :: oif_gradient
    logical :: is_fitting_data

    logical,save :: is_first_in=.true.
    type(t_field) :: fld_u, fld_a
    character(:),allocatable :: update_wavelet

    is_fitting_data=.true.

    fobj%dnorms=0.
    fobj%xnorms=0.
    
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
        
        !forward modeling
        call ppg%forward(fld_u);   call fld_u%acquire

        if(mpiworld%is_master) call shot%write('draw_',shot%dsyn)

        update_wavelet=setup%get_str('UPDATE_WAVELET',o_default='per shot')

        if(update_wavelet/='no') call shot%update_wavelet !call gradient_matchfilter_data
        
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

        !objective function and adjoint source
        call fobj%stack_dnorms
        shot%dadj=-shot%dadj
        
        if(mpiworld%is_master) call fobj%print_dnorms('Shotloop-stacked','upto Shot#'//shot%sindex)
        ! if(either(oif_gradient,.true.,present(oif_gradient))) then
        
            !adjoint source
            if(update_wavelet/='no') call shot%update_adjsource
            
            call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)
            
            call fld_a%ignite
            
            !adjoint modeling
            call ppg%adjoint(fld_a,fld_u,oif_compute_grad=.true.)
            
            call cb%project_back
            
        ! endif

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    !collect global objective function value
    call mpi_allreduce(mpi_in_place, fobj%dnorms, fobj%n_dnorms, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    !scale by shotlist
    call shls%scale(fobj%n_dnorms,o_from_sampled=fobj%dnorms)

    call fobj%print_dnorms('Shotloop-stacked, shotlist-scaled, but yet linesearch-scaled','')

    ! if(either(oif_gradient,.true.,present(oif_gradient))) then
        !collect global gradient
        call mpi_allreduce(mpi_in_place, m%gradient,  m%n*ppg%ngrad, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        !scale
        call shls%scale(m%n*ppg%ngrad,o_from_sampled=m%gradient)

        if(ppg%if_compute_engy) then
            call mpi_allreduce(mpi_in_place, m%energy  ,  m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
        endif
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
