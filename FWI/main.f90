program main
use m_gradient
use m_preconditioner
use m_parameterization
use m_optimizer

    character(:),allocatable :: job
    
    call init_mpiworld

    call hud('========================================')
    call hud('     WELCOME TO LEGO MODELING CODE      ')
    call hud('========================================')
    
    call init_setup(istat)
    
    if(istat==0) then !print manual
        call fwi_print_manual
        call mpiworld_finalize
        stop
    endif
    
    
    !read initial model
    call init_model
    call field_print_info
    
    call gradient_modeling(if_gradient=.true.)
    

    open(12,file='gradient',action='write',access='stream')
    write(12) gradient
    close(12)
        
    
!     open(12,file='precond_gradient',action='write',access='direct',recl=4*m%n*2)
!     write(12,rec=1) precond_gradient
!     close(12)
    
    ! job=get_setup_char('JOB',default='optimization')
    ! !if just estimate the wavelet or compute the gradient then this is it.
    ! if(job/='optimization') then
    !     call mpiworld_finalize
    !     stop
    ! endif
    
    
    !initialize parameterization
    call init_parameterization
    
    !initialize optimizer
    call init_optimizer(m%n*npar)
    
    open(12,file='initial%pg',action='write',access='direct',recl=4*m%n*npar)
    write(12,rec=1) current%pg
    close(12)

    job=get_setup_char('JOB',default='optimization')
    !if just estimate the wavelet or compute the gradient then this is it.
    if(job/='optimization') then
        call mpiworld_finalize
        stop
    endif
    
    call optimizer
    
    
    call mpiworld_finalize
    
end


subroutine fwi_print_manual
use m_mpienv
use m_field

    if(mpiworld%is_master) then
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
        write(*,'(a)') "                                           #LEGO will look for SU data file named as obs????.su,"
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
        write(*,'(a)') "                                            #If not given, LEGO will look for all data with sequential index (start from 0001)"
        write(*,'(a)') "                                            #the prefix of the filename has been defined in FILE_DATA"
        write(*,'(a)') ""
        write(*,'(a)') "FILE_WAVELET            'fricker'           #If not given, LEGO will use the wavelet generated from WAVELET_TYPE"
        write(*,'(a)') "WAVELET_TYPE            'sinexp'            #Damped sine function: sin(2pi*fpeak*t)*exp(-3.33*fpeak*t)"
        write(*,'(a)') "                                            #s.t. the wavelet has ~2 peaks before dying out, followed by lowpass filtering."
        write(*,'(a)') "                                            #If not given, use Ricker wavelet:"
        write(*,'(a)') "                                            #x(t) = -(pi*fpeak*(t-t0))^2, w(t) = (1+2*x)*exp(x)"
        write(*,'(a)') "                                            #t0 = 1/fpeak by default"
        write(*,'(a)') "RICKER_DELAYTIME        '1/PEAK_FREQUENCY'  #Change t0 in Ricker wavelet"
        write(*,'(a)') "UPDATE_WAVELET          'per shot'          #If update source wavelet after each forward modeling (using matching filter)"
        write(*,'(a)') "                                            #='per shot' update each shot's wavelet independently"
        write(*,'(a)') "                                            #='stack'    stack each shot's wavelet and use the same stacked wavelet for all shots"
        write(*,'(a)') "                                            #='no'       no update"
        write(*,'(a)') ""
        write(*,'(a)') "IF_HICKS                T    #Use Hicks method to interpolate source & receiver positions when not colocate with grid points."
        write(*,'(a)') ""
        write(*,'(a)') "                             #If F, sources & receivers will be positioned to nearest grid points."
        write(*,'(a)') "IF_SNAPSHOT             F                   #If T, save field snapshot to disk"
        write(*,'(a)') "SNAPSHOT_DELTA_IT       50                  #Save snapshot every SNAPSHOT_DELTA_IT time step"
        write(*,'(a)') ""
        write(*,'(a)') "FILE_WEIGHT             'fweight'           #ASCII file specifying weighting operator applied to data residuals in time-offset domain;"
        write(*,'(a)') "                                            #if not given, use weight=1"
        write(*,'(a)') "                        #Two types of weighting can be specified:"
        write(*,'(a)') "                        #1) partition the gather by polygonal chains with different weights based on interpretation of phases"
        write(*,'(a)') "                        # (e.g. if you want to separate reflections from diving waves) (modified from SU program sumute)"
        write(*,'(a)') "                        #2) arbitrary weighting by anchor points (particularly time or offset)"
        write(*,'(a)') "                        #Example of fweight that mutes reflections"
        write(*,'(a)') "                        #1          #number of polygonal chains"
        write(*,'(a)') "                        #2  0  50   #1st chain has 2 vertices (ie. endpoints), weight=1 above and =0 below this chain (50-point transition)"
        write(*,'(a)') "                        #0    1000  #offset of vertices (suppose the gather's max offset=1000m)"
        write(*,'(a)') "                        #0.2  1.7   #time instant of vertices (suppose the direct wave has a period of 0.2s "
        write(*,'(a)') "                        #           # and uses 1.5s to travel from 0 to 1000m offset)"
        write(*,'(a)') "                        #more polygonal chains can be added; each is defined by 3 lines like above"
        write(*,'(a)') "                        #Example of fweight that has linearly increasing weight with offset:"
        write(*,'(a)') "                        #0               #no polygonal chains"
        write(*,'(a)') "                        #'weight_table'  #binary-format file providing the weights, e.g. 0 0.2 0.4 0.6 0.8 1"
        write(*,'(a)') "                        #'1 1. 6 200.'   #number of samples (1) and spacing (1s) in time dir, and"
        write(*,'(a)') "                        #                #number of samples (6) and spacing (200m) in space dir, given by 'weight_table'"
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

    endif

    if(mpiworld%is_master) then
        write(*,'(a)') "---------------------------"
        write(*,'(a)') "Forward operator info:"
        write(*,'(a)') "---------------------------"
        call field_print_info
        write(*,'(a)') ""
    endif
end subroutine
