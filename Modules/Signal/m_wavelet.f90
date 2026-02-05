module m_wavelet
use m_System
use m_math, only: r_pi
use m_butterworth
    
    contains

    function wavelet_sinexp(nt,dt,fpeak) result(wavelet)
        ! Input parameters
        integer, intent(in) :: nt      ! Number of time samples
        real, intent(in) :: dt         ! Time step (sampling interval)
        real, intent(in) :: fpeak      ! Peak frequency of the wavelet

        ! Output
        real,dimension(:),allocatable :: wavelet

        ! Local variables
        integer :: it                  ! Loop counter for time samples
        real :: t                      ! Current time
        real :: a                      ! Exponential decay factor

        ! Validate input parameters to avoid floating-point exceptions
        if (nt <= 0) call hud('ERROR: wavelet_sinexp received invalid nt='//num2str(nt))
        if (dt <= 0. .or. isnan(dt)) call hud('ERROR: wavelet_sinexp received invalid dt='//num2str(dt))
        if (fpeak <= 0. .or. isnan(fpeak)) call hud('ERROR: wavelet_sinexp received invalid fpeak='//num2str(fpeak))

        a=-3.3333333*fpeak

        call alloc(wavelet,nt)

        do it=1,nt
            t=(it-1)*dt

            ! Calculate exponential decay, clamping to avoid underflow/denormal
            ! When a*t is very negative, exp(a*t) approaches 0
            if (a*t < -50.) then
                wavelet(it) = 0.
            else
                wavelet(it)=sin(2.*r_pi*fpeak*t)*exp(a*t)
            endif
        enddo

        !butterworth filtering to mitigate spectrum high-end tail
        call butterworth(wavelet,nt,1,dt, &
                        ois_zerophase=.false.,oif_locut=.false., &
                        o_fpasshi=fpeak,o_fstophi=2.*fpeak,o_astophi=0.1)

    end function
    
    function wavelet_ricker(nt,dt,fpeak) result(wavelet)
        ! Input parameters
        integer, intent(in) :: nt      ! Number of time samples
        real, intent(in) :: dt         ! Time step (sampling interval)
        real, intent(in) :: fpeak      ! Peak frequency of the wavelet

        ! Output
        real,dimension(:),allocatable :: wavelet

        ! Local variables
        integer :: it                  ! Loop counter for time samples
        real :: t                      ! Current time relative to delay time
        real :: t0                     ! Delay time (time shift for wavelet)
        real :: x                      ! Temporary variable for calculation

        ! Validate input parameters to avoid floating-point exceptions
        if (nt <= 0) call hud('ERROR: wavelet_ricker received invalid nt='//num2str(nt))
        if (dt <= 0. .or. isnan(dt)) call hud('ERROR: wavelet_ricker received invalid dt='//num2str(dt))
        if (fpeak <= 0. .or. isnan(fpeak)) call hud('ERROR: wavelet_ricker received invalid fpeak='//num2str(fpeak))

        t0=setup%get_real('RICKER_DELAYTIME','T0',o_default=num2str(1./fpeak))

        if (fpeak*2.5 > 1./dt) then
            call warn('Ricker wavelet peak frequency too high (fpeak*2.5 > 1/dt). Reduce it.')
        endif

        call alloc(wavelet,nt)

        do it=1,nt

            t=(it-1)*dt-t0

            x=r_pi*fpeak*t
            x=-x*x
            wavelet(it)=(1.+2.*x)*exp(x)

        enddo

    end function

end