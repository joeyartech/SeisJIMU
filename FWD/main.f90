program main
use m_System
use m_Modeling

    type(t_field) :: field

    !character(:),allocatable :: job  

    !mpiworld lives in t_mpienv
    call mpiworld%init(name='MPIWorld')

    call hud('======================================'//s_NL// &
             '   WELCOME TO SeisJIMU FWD MODELING   '//s_NL// &
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

    !estimate required memory
    ! call m%estim_RAM
    ! call cb%estim_RAM
    ! call field%estim_RAM
    ! call mpiworld%fin

    !checkpoint
    call checkpoint_init

    !model
    call m%init
    call m%read
    call ppg%check_model

    !shotlist
    call shls%read_from_setup
    call shls%build(o_batchsize=shls%nshot)
    call shls%sample
    call shls%assign

    call hud('===== START LOOP OVER SHOTS =====')

    do i=1,shls%nshots_per_processor

        call shot%init(shls%yield(i))
        call shot%read_from_setup
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        call hud('Modeling Shot# '//shot%sindex)

        call cb%init(ppg%nbndlayer)
        call cb%project

        call ppg%check_discretization
        call ppg%init
        call ppg%init_field(field,name='field',origin='src')
        call ppg%init_abslayer

        call field%ignite

        !forward modeling
        call ppg%forward(field)

        call field%acquire

        !write synthetic data
        call shot%write('dsyn_',shot%dsyn)

    enddo
    
    call hud('        END LOOP OVER SHOTS        ')
    
    call mpiworld%final
    
    ! stop
    
end


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
    write(*,'(a)') "ACQUI_TYPE              'spread'           #Type of acquisition geometry"
    write(*,'(a)') "                                           #2D:"
    write(*,'(a)') "                                           #='spread' split-spread geometry, i.e. sources share all receivers, typical in land acquisition"
    write(*,'(a)') "                                           #='streamer' limited offset range for each shot, typical in marine acquisition"
    write(*,'(a)') "                                           #='irregularOBN' irregular source & receiver positions, to be read from ascii files specified in SOURCE_LINE & RECEIVER_LINE"
    write(*,'(a)') "                                           #3D:"
    write(*,'(a)') "                                           #='spread3D' split-spread geometry, i.e. sources share all receivers, typical in land acquisition"
    write(*,'(a)') "                        #fz fx  fy    lz lx   ly    n"
    write(*,'(a)') "SOURCE_LINE             '10 100 1000  10 1900 1000  20'"
    write(*,'(a)') "RECEIVER_LINE           ' 5   0 1000   5  100 1000  11'"
    write(*,'(a)') "                        #Source & receiver positions:"
    write(*,'(a)') "                        #ACQUI_TYPE=='spread': Needing 7 numbers in SOURCE_LINE:"
    write(*,'(a)') "                        #                  first 3 numbers give (z,x,y) coordinate of first source,"
    write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of last source,"
    write(*,'(a)') "                        #                  followed by 1 integer number gives total number of sources."
    write(*,'(a)') "                        #                  Position of intermediate sources will be interpolated."
    write(*,'(a)') "                        #                  Same for RECEIVER_LINE"
    write(*,'(a)') "                        #ACQUI_TYPE=='streamer': Needing 7 numbers in SOURCE_LINE:"
    write(*,'(a)') "                        #                  first 3 numbers give (z,x,y) coordinate of first source,"
    write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of last source,"
    write(*,'(a)') "                        #                  followed by 1 integer number gives total number of sources."
    write(*,'(a)') "                        #                  Position of intermediate sources will be interpolated."
    write(*,'(a)') "                        #                  Needing 7 numbers in RECEIVER_LINE:"
    write(*,'(a)') "                        #                  first 3 numbers give (z,x,y) coordinate of nearest receiver (smallest offset),"
    write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of farthest receiver (largest offset),"
    write(*,'(a)') "                        #                  followed by 1 integer number gives number of receiver for each shot."
    write(*,'(a)') "                        #                  Position of intermediate receivers will be interpolated."
    write(*,'(a)') "                        #ACQUI_TYPE=='irregularOBN': Needing file name in SOURCE_LINE, which contains a Nx3 matrix (columns separated by white spaces),"
    write(*,'(a)') "                        #                  where N=total number of sources, and 3 columns specify (z,x,y) coordinate of sources"
    write(*,'(a)') "                        #                  Same for RECEIVER_LINE"
    write(*,'(a)') "                        #ACQUI_TYPE=='spread3D': We need 4 anchor points (P1-4) for a quadrilateral geometry of sources."
    write(*,'(a)') "                        #                  Thus needing 14 numbers in SOURCE_LINE:"
    write(*,'(a)') "                        #                  first 3 numbers give (z,x,y) coordinate of P1,"
    write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of P2,"
    write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of P3,"
    write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of P4,"
    write(*,'(a)') "                        #                  followed by 1 integer gives total number of points in inline direction,"
    write(*,'(a)') "                        #                  and 1 integer gives total number of points in crossline direction."
    write(*,'(a)') "                        #                  An intermediate source will be interpolated from the two endpoints of the corresponding line where the source sits,"
    write(*,'(a)') "                        #                  where the endpoints are also interpolated by distances to P1-P3 and P2-P4, respectively."
    write(*,'(a)') "                        #                  Same for RECEIVER_LINE"
    write(*,'(a)') "#Other acqui types or source & receiver lines will be considered and developed in future."
    write(*,'(a)') ""
    write(*,'(a)') "SOURCE_COMPONENT        1"
    write(*,'(a)') "RECEIVER_COMPONENT      1"
    write(*,'(a)') "                        #source & receiver component: 1=Isotropic Pressure, 2=Particle velocity in x (vx), 3=in y (vy), 4=in z (vz)"
    write(*,'(a)') "#Multi-component data will be considered and developed in future."
    write(*,'(a)') ""
    write(*,'(a)') "PEAK_FREQUENCY          7                  #Peak frequency of source wavelet"
    write(*,'(a)') ""
    write(*,'(a)') "TIME_STEP               500                #Total number of time step (nt)"
    write(*,'(a)') "TIME_INTERVAL           0.006              #Time step interval (dt)"
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
    write(*,'(a)') ""
    write(*,'(a)') "IF_ISOTROPIC            F                   #If T, omit model_eps model_del"
    write(*,'(a)') "IF_TOPO_FROM_VS         F                   #If T, decide topography from model_vs when model_vs exists and model_topo doesn't exist"
    write(*,'(a)') "IF_FREESURFACE          T                   #Free surface condition (T) or absorbing surface condition (F)"
    write(*,'(a)') ""
    write(*,'(a)') "FILE_WAVELET            'fricker_dt6000'    #If not given, SeisJIMU will use wavelet as specified in WAVELET_TYPE"
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
    write(*,'(a)') ""
    write(*,'(a)') "IF_HICKS                T    #Use Hicks method to interpolate source & receiver positions when not colocate with grid points;"
    write(*,'(a)') "                             #if F, sources & receivers will be positioned to nearest grid points."
    write(*,'(a)') ""
    write(*,'(a)') "IF_SNAPSHOT             F    #If T, save field snapshot to disk"
    write(*,'(a)') "SNAPSHOT_DELTA_IT       50   #Save snapshot every SNAPSHOT_DELTA_IT time step"
    write(*,'(a)') ""
    write(*,'(a)') "---------------------------"
    write(*,'(a)') "Notes:"
    write(*,'(a)') "---------------------------"
    write(*,'(a)') ""
    write(*,'(a)') "Modeled seismic data are in binary format (other format can be developed in future)."
    write(*,'(a)') ""

end subroutine
