module m_gen_wavelet
use m_arrayop
use m_sysio
use m_butterworth

    private pi

    real,parameter :: pi=3.1415927
    
    contains

    function gen_wavelet_sinexp(nt,dt,fpeak) result(wl)
        integer :: nt
        real :: dt, fpeak
        real,dimension(:),allocatable :: wl
                
        a=-3.3333333*fpeak
        
        call alloc(wl,nt)
        
        do it=1,nt
            t=(it-1)*dt
            
            wl(it)=sin(2.*pi*fpeak*t)*exp(a*t)
        enddo
        
        !butterworth filtering to mitigate spectrum high-end tail
        call butterworth(nt,dt,wl,&
        o_zerophase=.false.,o_locut=.false.,&
        o_fpasshi=fpeak,o_fstophi=2.*fpeak,&
                        o_astophi=0.1)
        
    end function
    
    function gen_wavelet_ricker(nt,dt,fpeak) result(wl)
        integer :: nt
        real :: dt, fpeak
        real,dimension(:),allocatable :: wl
        
        real t0
        
        t0=get_setup_real('RICKER_DELAYTIME',default=1./fpeak)
        
        if (fpeak*2.5 > 1./dt) then
            call hud('Ricker wavelet peak frequency too high (fpeak*2.5 > 1/dt). Reduce it.')
        endif

        call alloc(wl,nt)

        do it=1,nt

            t=(it-1)*dt-t0

            x=-pi*pi*fpeak*fpeak*t*t
            wl(it)=(1.+2.*x)*exp(x)

        enddo

    end function

end