module m_wavelets

    
    subroutine source_wavelet_sinexp

        a=-3.3333333*shot%src%fpeak
        
        call alloc(shot%src%wavelet,shot%src%nt)
        
        do it=1,shot%src%nt
            t=(it-1)*shot%src%dt
            
            shot%src%wavelet(it)=sin(2.*r_pi*shot%src%fpeak*t)*exp(a*t)
        enddo
        
        !butterworth filtering to mitigate spectrum high-end tail
        call butterworth(1,shot%src%nt, shot%src%dt, shot%src_wavelet,&
        o_zerophase=.false.,o_locut=.false.,&
        o_fpasshi=fpeak,o_fstophi=2.*fpeak,&
                        o_astophi=0.1)
        
    end subroutine
    
    subroutine source_wavelet_ricker

        t0=setup_get_real('RICKER_DELAYTIME','T0',default=1./shot%src%fpeak)
        
        if (shot%src%fpeak*2.5 > 1./shot%src%dt) then
            call hud('Ricker wavelet peak frequency too high (fpeak*2.5 > 1/dt). Reduce it.')
        endif

        call alloc(shot%src%wavelet,shot%src%nt)

        do it=1,shot%src%nt

            t=(it-1)*shot%src%dt-t0

            x=r_pi*shot%src%fpeak*t
            x=-x*x
            shot%src%wl(it)=(1.+2.*x)*exp(x)

        enddo

    end subroutine