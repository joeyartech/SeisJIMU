program main
use m_System
use m_math
use m_hilbert
use m_hilbert_nofft

    integer,parameter :: nt=500
    real :: dt=0.004, fpeak=7
    
    real,dimension(:,:),allocatable :: wl, wlconj, wlenv, wlph
    call alloc(wl,nt,1)
    call alloc(wlconj,nt,1)
    call alloc(wlenv,nt,1)
    call alloc(wlph,nt,1)
    
    wl(:,1)=ricker(nt,dt,fpeak)
    
    call cpu_time(tic)
    do i=1,500
        call hilbert_transform(wl,wlconj,nt,1)
    enddo
    call cpu_time(toc)
    tt0=toc-tic
    print*,'hilbert_fft elapsed time (500x):',tt0
    ! call hilbert_envelope(wl,wlenv,nt,1)
    ! call hilbert_phase(wl,wlph,nt,1)
    
    open(12,file='wavelets',access='direct',recl=4*nt)
    write(12,rec=1)wl
    write(12,rec=2)wlconj
    ! write(12,rec=3)wlenv
    ! write(12,rec=4)wlph
    !close(12)
    
    call hilbert_nofft('naive',wl,wlconj,nt,1)

    call cpu_time(tic)
    do i=1,500
        call hilbert_nofft('naive',wl,wlconj,nt,1) !incorrect result, why?
    enddo
    call cpu_time(toc)
    print*,'hilbert_nofft naive elapsed time (500x):',(toc-tic), (toc-tic)/tt0
    !open(12,file='wavelets',access='direct',recl=4*nt)
    write(12,rec=3)wlconj


    call cpu_time(tic)
    do i=1,500
        call hilbert_nofft('generic',wl,wlconj,nt,1)
    enddo
    call cpu_time(toc)
    print*,'hilbert_nofft generic elapsed time (500x):',(toc-tic), (toc-tic)/tt0
    !open(12,file='wavelets',access='direct',recl=4*nt)
    write(12,rec=4)wlconj
    !call hilbert_nofft(wlconj,wl,nt,1)

    
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
