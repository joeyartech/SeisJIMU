module m_wavelets
use m_global
use m_string
use m_setup
use m_butterworth
    
    contains

    function sinexp(nt,dt,fpeak) result(wavelet)
        real,dimension(:),allocatable :: wavelet

        a=-3.3333333*fpeak
        
        call alloc(wavelet,nt)
        
        do it=1,nt
            t=(it-1)*dt
            
            wavelet(it)=sin(2.*r_pi*fpeak*t)*exp(a*t)
        enddo
        
        !butterworth filtering to mitigate spectrum high-end tail
        call butterworth(1,nt,dt,wavelet,&
        o_zerophase=.false.,o_locut=.false.,&
        o_fpasshi=fpeak,o_fstophi=2.*fpeak,&
                        o_astophi=0.1)
        
    end function
    
    function ricker(nt,dt,fpeak) result(wavelet)
        real,dimension(:),allocatable :: wavelet

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