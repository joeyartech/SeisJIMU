program main
use m_hilbert

    integer,parameter :: nt=500
    real :: dt=0.004, fpeak=7
    
    real,dimension(nt,1) :: wl, wlconj, wlenv, wlph
    
    wl(:,1)=ricker(nt,dt,fpeak)
    
    call hilbert_transform(wl,wlconj,nt,1)
    call hilbert_envelope(wl,wlenv,nt,1)
    call hilbert_phase(wl,wlph,nt,1)
    
    open(12,file='wavelets',access='direct',recl=4*nt)
    write(12,rec=1)wl
    write(12,rec=2)wlconj
    write(12,rec=3)wlenv
    write(12,rec=4)wlph
    close(12)
    
    contains
    
    function ricker(nt,dt,fpeak) result(wavelet)
        real,dimension(:),allocatable :: wavelet

        t0=1/fpeak
        
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
    
end program
