program main
use m_model
use m_gen_acquisition
use m_geometry
use m_freqlist
use m_computebox
use m_field

    call init_mpiworld

    call hud('========================================')
    call hud('     WELCOME TO LEGO MODELING CODE      ')
    call hud('========================================')
    
    call init_setup(istat)
    
    if(istat==0) then !print manual
        ! call fwd_print_manual
        call mpiworld_finalize
        stop
    endif
    
    
if(mpiworld%is_master) then
    !read initial model
    call init_model
    
    !generate acquisition
    call gen_acquisition
    !Ideally, computebox can be a portion of the whole model, like in time-domain FWI
    !to reduce the size of RAM required for LU decomposition.
    !in this case, an inner loop over shots is needed inside frequency loop,
    !and computebox, S/R position 1D arrays, matrix analysis etc. 
    !should be re-done inside the inner loop.
    !However, I'm not sure how frequently freq-domain FWI will be employed,
    !and how large the problem can be.
    !To keep things simple, this inner loop is not implemented for now
    !and thus the matrix shape, S/R arrays etc. stay same (no need to rebuild)

    !loop over shots to extract acqui geom for RHS
    call build_geom_acqui
endif

    !initialize mumps
    call init_field_mumps

    !read frequency list
    if(mpiworld%is_master) call build_freqlist
    
    call hud('      START LOOP OVER FREQUENCIES         ')

    do i=1,nfreq

        if(mpiworld%is_master) then
            write(*,*) 'Modeling freq# ',i,freq(i),'Hz'
        endif
        if(mpiworld%is_master) call build_computebox(freq(i))

        !initialize field & extmodel
        call init_field_extmodel(freq(i))

        call field_matrix_factorize

        call field_RHS_substitute

        call field_write(i)

    enddo
    
    call hud('        END LOOP OVER FREQUENCIES         ')
    
    call field_mumps_finalize
    call mpiworld_finalize
    
end


! subroutine fwd_print_manual
! use m_mpienv
! use m_field

!     if(mpiworld%is_master) then
!         write(*,'(a)') ""
!         write(*,'(a)') "----------------------------"
!         write(*,'(a)') "To launch the program, do:"
!         write(*,'(a)') "----------------------------"
!         write(*,'(a)') ""
!         write(*,'(a)') "bash $ export OMP_NUM_THREADS=1 #(use source in csh like shell)"
!         write(*,'(a)') "bash $ mpirun -np $np FWI setup.in"
!         write(*,'(a)') ""
!         write(*,'(a)') "----------------------------"
!         write(*,'(a)') "Mandatory items in setup.in:"
!         write(*,'(a)') "----------------------------"
!         write(*,'(a)') ""
!         write(*,'(a)') "FILE_MODELS             'model'            #Prefix of input model name"
!         write(*,'(a)') "MODEL_DIMENSION         '201  201  201'    #nz nx ny"
!         write(*,'(a)') "MODEL_SPACING           ' 10   10   10'    #dz dx dy"
!         write(*,'(a)') "MODEL_ORIGIN            '  0    0    0'    #oz ox oy"
!         write(*,'(a)') "                        #all models should have same dimension, spacing and origin"
!         write(*,'(a)') ""
!         write(*,'(a)') "NBOUNDARYLAYER          60                 #Number of convolutional Perfect-Match Layers added to edges of computebox"
!         write(*,'(a)') ""
!         write(*,'(a)') "ACQUI_TYPE              'spread'           #Type of acquisition geometry"
!         write(*,'(a)') "                                           #2D:"
!         write(*,'(a)') "                                           #='spread' split-spread geometry, i.e. sources share all receivers, typical in land acquisition"
!         write(*,'(a)') "                                           #='streamer' limited offset range for each shot, typical in marine acquisition"
!         write(*,'(a)') "                                           #='irregularOBN' irregular source & receiver positions, to be read from ascii files specified in SOURCE_LINE & RECEIVER_LINE"
!         write(*,'(a)') "                                           #3D:"
!         write(*,'(a)') "                                           #='spread3D' split-spread geometry, i.e. sources share all receivers, typical in land acquisition"
!         write(*,'(a)') "                        #fz fx  fy    lz lx   ly    n"
!         write(*,'(a)') "SOURCE_LINE             '10 100 1000  10 1900 1000  20'"
!         write(*,'(a)') "RECEIVER_LINE           ' 5   0 1000   5  100 1000  11'"
!         write(*,'(a)') "                        #Source & receiver positions:"
!         write(*,'(a)') "                        #ACQUI_TYPE=='spread': Needing 7 numbers in SOURCE_LINE:"
!         write(*,'(a)') "                        #                  first 3 numbers give (z,x,y) coordinate of first source,"
!         write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of last source,"
!         write(*,'(a)') "                        #                  followed by 1 integer number gives total number of sources."
!         write(*,'(a)') "                        #                  Position of intermediate sources will be interpolated."
!         write(*,'(a)') "                        #                  Same for RECEIVER_LINE"
!         write(*,'(a)') "                        #ACQUI_TYPE=='streamer': Needing 7 numbers in SOURCE_LINE:"
!         write(*,'(a)') "                        #                  first 3 numbers give (z,x,y) coordinate of first source,"
!         write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of last source,"
!         write(*,'(a)') "                        #                  followed by 1 integer number gives total number of sources."
!         write(*,'(a)') "                        #                  Position of intermediate sources will be interpolated."
!         write(*,'(a)') "                        #                  Needing 7 numbers in RECEIVER_LINE:"
!         write(*,'(a)') "                        #                  first 3 numbers give (z,x,y) coordinate of nearest receiver (smallest offset),"
!         write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of farthest receiver (largest offset),"
!         write(*,'(a)') "                        #                  followed by 1 integer number gives number of receiver for each shot."
!         write(*,'(a)') "                        #                  Position of intermediate receivers will be interpolated."
!         write(*,'(a)') "                        #ACQUI_TYPE=='irregularOBN': Needing file name in SOURCE_LINE, which contains a Nx3 matrix (columns separated by white spaces),"
!         write(*,'(a)') "                        #                  where N=total number of sources, and 3 columns specify (z,x,y) coordinate of sources"
!         write(*,'(a)') "                        #                  Same for RECEIVER_LINE"
!         write(*,'(a)') "                        #ACQUI_TYPE=='spread3D': We need 4 anchor points (P1-4) for a quadrilateral geometry of sources."
!         write(*,'(a)') "                        #                  Thus needing 14 numbers in SOURCE_LINE:"
!         write(*,'(a)') "                        #                  first 3 numbers give (z,x,y) coordinate of P1,"
!         write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of P2,"
!         write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of P3,"
!         write(*,'(a)') "                        #                  then  3 numbers give (z,x,y) coordinate of P4,"
!         write(*,'(a)') "                        #                  followed by 1 integer gives total number of points in inline direction,"
!         write(*,'(a)') "                        #                  and 1 integer gives total number of points in crossline direction."
!         write(*,'(a)') "                        #                  An intermediate source will be interpolated from the two endpoints of the corresponding line where the source sits,"
!         write(*,'(a)') "                        #                  where the endpoints are also interpolated by distances to P1-P3 and P2-P4, respectively."
!         write(*,'(a)') "                        #                  Same for RECEIVER_LINE"
!         write(*,'(a)') "#Other acqui types or source & receiver lines will be considered and developed in future."
!         write(*,'(a)') ""
!         write(*,'(a)') "SOURCE_COMPONENT        1"
!         write(*,'(a)') "RECEIVER_COMPONENT      1"
!         write(*,'(a)') "                        #source & receiver component: 1=Isotropic Pressure, 2=Particle velocity in x (vx), 3=in y (vy), 4=in z (vz)"
!         write(*,'(a)') "#Multi-component data will be considered and developed in future."
!         write(*,'(a)') ""
!         write(*,'(a)') "PEAK_FREQUENCY          7                  #Peak frequency of source wavelet"
!         write(*,'(a)') ""
!         write(*,'(a)') "TIME_STEP               500                #Total number of time step (nt)"
!         write(*,'(a)') "TIME_INTERVAL           0.006              #Time step interval (dt)"
!         write(*,'(a)') ""
!         write(*,'(a)') "---------------------------"
!         write(*,'(a)') "Optional items in setup.in:"
!         write(*,'(a)') "---------------------------"
!         write(*,'(a)') ""
!         write(*,'(a)') "APERTURE                '-99999 99999 -99999 99999'    #Define size of computebox"
!         write(*,'(a)') "                                                       #='0 0 0 0' computebox bounded by farthest source & receiver"
!         write(*,'(a)') "                                                       # (i.e. leftmost & rightmost device defines fx & lx, similar for fy & ly)"
!         write(*,'(a)') "                                                       # (smallest possible computebox otherwise some devices will be truncated)"
!         write(*,'(a)') "                                                       #='-99999 99999 -99999 99999' use whole model as computebox"
!         write(*,'(a)') "                                                       # (largest possible computebox otherwise it will be larger than the model)"
!         write(*,'(a)') "                                                       #Well-chosen values between above two extremes can give a good balance"
!         write(*,'(a)') "                                                       # between edge-effect and computational cost"
!         write(*,'(a)') ""
!         write(*,'(a)') "IF_ISOTROPIC            F                   #If T, omit model_eps model_del"
!         write(*,'(a)') "IF_TOPO_FROM_VS         F                   #If T, decide topography from model_vs when model_vs exists and model_topo doesn't exist"
!         write(*,'(a)') "IF_FREESURFACE          T                   #Free surface condition (T) or absorbing surface condition (F)"
!         write(*,'(a)') ""
!         write(*,'(a)') "FILE_WAVELET            'fricker_dt6000'    #If not given, LEGO will use wavelet as specified in WAVELET_TYPE"
!         write(*,'(a)') "WAVELET_TYPE            'sinexp'            #Damped sine function: sin(2pi*fpeak*t)*exp(-3.33*fpeak*t)"
!         write(*,'(a)') "                                            #s.t. the wavelet has ~2 peaks before dying out, followed by lowpass filtering."
!         write(*,'(a)') "                                            #If not given, use Ricker wavelet:"
!         write(*,'(a)') "                                            #x(t) = -(pi*fpeak*(t-t0))^2, w(t) = (1+2*x)*exp(x)"
!         write(*,'(a)') "                                            #t0 = 1/fpeak by default"
!         write(*,'(a)') "RICKER_DELAYTIME        '1/PEAK_FREQUENCY'  #Change t0 in Ricker wavelet"
!         write(*,'(a)') "SCALE_WAVELET           'no'                #If scale source wavelet when read in"
!         write(*,'(a)') "                                            #='no'       no scale"
!         write(*,'(a)') "                                            #='by dxdt'  multiply by dt/dx/dy/dz s.t. the amplitude is independent on discretization"
!         write(*,'(a)') "                                            #=number     user-given scaler on wavelet"
!         write(*,'(a)') ""
!         write(*,'(a)') "IF_HICKS                T    #Use Hicks method to interpolate source & receiver positions when not colocate with grid points;"
!         write(*,'(a)') "                             #if F, sources & receivers will be positioned to nearest grid points."
!         write(*,'(a)') ""
!         write(*,'(a)') "IF_SNAPSHOT             F    #If T, save field snapshot to disk"
!         write(*,'(a)') "SNAPSHOT_DELTA_IT       50   #Save snapshot every SNAPSHOT_DELTA_IT time step"
!         write(*,'(a)') ""
!         write(*,'(a)') "---------------------------"
!         write(*,'(a)') "Notes:"
!         write(*,'(a)') "---------------------------"
!         write(*,'(a)') ""
!         write(*,'(a)') "Modeled seismic data are in binary format (other format can be developed in future)."
!         write(*,'(a)') ""

!     endif

!     if(mpiworld%is_master) then
!         write(*,'(a)') "---------------------------"
!         write(*,'(a)') "Forward operator info:"
!         write(*,'(a)') "---------------------------"
!         call field_print_info
!         write(*,'(a)') ""
!     endif

! end subroutine
