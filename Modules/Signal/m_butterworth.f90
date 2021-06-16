module m_butterworth
use m_either

    private
    public :: butterworth

    contains

! A Fortran version of subfilt.c in Seismic Unix
! dir: SeisUnix/src/su/main/filters
!
! !! Copyright (c) Colorado School of Mines, 2011.!!
! !! All rights reserved.                       !!
! 
! !! SUBFILT: $Revision: 1.22 $ ; $Date: 2012/11/28 22:13:13 $    !!
! 
! !!********************** self documentation *********************!!
! "                                 ",
! " SUBFILT - apply Butterworth bandpass filter             ",
! "                                 ",
! " subfilt <stdin >stdout [optional parameters]            ",
! "                                     ",
! " Required parameters:                        ",
! "     if dt is not set in header, then dt is mandatory    ",
! "                                     ",
! " Optional parameters: (nyquist calculated internally)        ",
! "     zerophase=1        =0 for minimum phase filter     ",
! "     locut=1            =0 for no low cut filter     ",
! "     hicut=1            =0 for no high cut filter     ",
! "     fstoplo=0.10*(nyq)    freq(Hz) in low cut stop band    ",
! "     astoplo=0.05        upper bound on amp at fstoplo     ",
! "     fpasslo=0.15*(nyq)    freq(Hz) in low cut pass band    ",
! "     apasslo=0.95        lower bound on amp at fpasslo     ",
! "     fpasshi=0.40*(nyq)    freq(Hz) in high cut pass band    ",
! "     apasshi=0.95        lower bound on amp at fpasshi     ",
! "     fstophi=0.55*(nyq)    freq(Hz) in high cut stop band    ",
! "     astophi=0.05        upper bound on amp at fstophi     ",
! "     verbose=0        =1 for filter design info     ",
! "     dt = (from header)    time sampling interval (sec)    ",
! "                                     ",
! " ... or  set filter by defining  poles and 3db cutoff frequencies",
! "    npoleselo=calculated     number of poles of the lo pass band",
! "    npolesehi=calculated     number of poles of the lo pass band",
! "    f3dblo=calculated    frequency of 3db cutoff frequency",
! "    f3dbhi=calculated    frequency of 3db cutoff frequency",
! "                                     ",
! " Notes:                                ",
! " Butterworth filters were originally of interest because they  ",
! " can be implemented in hardware form through the combination of",
! " inductors, capacitors, and an amplifier. Such a filter can be ",
! " constructed in such a way as to have very small oscillations    ",
! " in the flat portion of the bandpass---a desireable attribute.    ",
! " Because the filters are composed of LC circuits, the impulse  ",
! " response is an ordinary differential equation, which translates",
! " into a polynomial in the transform domain. The filter is expressed",
! " as the division by this polynomial. Hence the poles of the filter",
! " are of interest.                            ",
! "                                     ",
! " The user may define low pass, high pass, and band pass filters",
! " that are either minimum phase or are zero phase.  The default    ",
! " is to let the program calculate the optimal number of poles in",
! " low and high cut bands.                     ",
! "                                     ",
! " Alternately the user may manually define the filter by the 3db",
! " frequency and by the number of poles in the low and or high    ",
! " cut region.                             ",
! "                                     ",
! " The advantage of using the alternate method is that the user  ",
! " can control the smoothness of the filter. Greater smoothness  ",
! " through a larger pole number results in a more bell shaped    ",
! " amplitude spectrum.                        ",
! "                                     ",
! " For simple zero phase filtering with sin squared tapering use ",
! " \"sufilter\".                                ",
! NULL};
! 
! !! Credits:
!  *    CWP: Dave Hale c. 1993 for bf.c subs and test drivers
!  *    CWP: Jack K. Cohen for su wrapper c. 1993
!  *      SEAM Project: Bruce Verwest 2009 added explicit pole option
!  *                    in a program called "subfiltpole"
!  *      CWP: John Stockwell (2012) combined Bruce Verwests changes
!  *           into the original subfilt.
!  *
!  * Caveat: zerophase will not do good if trace has a spike near
!  *       the end.  One could make a try at getting the "effective"
!  *       length of the causal filter, but padding the traces seems
!  *       painful in an already expensive algorithm.
!  *
!  *
!  * Theory:
!  * The 
!  *
!  * Trace header fields accessed: ns, dt, trid
!  !!
! !!*************** end self doc **********************************!!

    subroutine butterworth(ntr,nt,dt,trace,&
        ois_zerophase,oif_locut,oif_hicut,&
        o_fstoplo,o_fpasslo,o_fpasshi,o_fstophi,&
        o_astoplo,o_apasslo,o_apasshi,o_astophi,&
        o_npoleslo,o_npoleshi,&
        o_f3dblo,  o_f3dbhi)

        integer ntr !number of traces
        integer nt  !number of time samples
        real,dimension(ntr,nt) :: trace
        real dt     !sample spacing
        
        logical,optional :: ois_zerophase !flag for zero phase filtering
        logical,optional :: oif_locut     !flag for low cut filtering
        logical,optional :: oif_hicut     !flag for high cut filtering
        logical zerophase, locut, hicut
        
        real,optional :: o_fstoplo        !left lower corner frequency
        real,optional :: o_fpasslo        !left upper corner frequency
        real,optional :: o_fpasshi        !right lower corner frequency
        real,optional :: o_fstophi        !right upper corner frequency
        real,optional :: o_astoplo        !amp at fstoplo
        real,optional :: o_apasslo        !amp at fpasslo
        real,optional :: o_apasshi        !amp at fpasshi
        real,optional :: o_astophi        !amp at fstophi
        real fstoplo,fpasslo,fpasshi,fstophi,astoplo,apasslo,apasshi,astophi
        
        integer,optional :: o_npoleslo    !poles in low cut filter
        integer,optional :: o_npoleshi    !poles in high cut filter
        real,optional :: o_f3dblo    !3 db point of low cut filter
        real,optional :: o_f3dbhi    !3 db point of high cut filter
        integer npoleslo,npoleshi
        real f3dblo,f3dbhi
        
        real nyq        !nyquist frequency
        integer verbose        !design info flag


        nyq = 0.5/dt

        zerophase=either(ois_zerophase,.true.,present(ois_zerophase))

        locut=either(oif_locut,.true.,present(oif_locut))
        hicut=either(oif_hicut,.true.,present(oif_hicut))

        !Get design frequencies and normalize to [0, 0.5] for bfdesign
        fstoplo=either(o_fstoplo,.10*nyq,present(o_fstoplo)) *dt
        fpasslo=either(o_fpasslo,.15*nyq,present(o_fpasslo)) *dt
        fpasshi=either(o_fpasshi,.40*nyq,present(o_fpasshi)) *dt
        fstophi=either(o_fstophi,.55*nyq,present(o_fstophi)) *dt
        ! if (locut) {
        !         if (fstoplo <= 0.0)      err("fstoplo must be positive")
        !         if (fstoplo > fpasslo)  err("fstoplo must be < fpasslo")
        ! }
        ! if (hicut) {
        !         if (fpasshi > fstophi)  err("fpasshi must be < fstophi")
        !         if (fstophi > nyq)  err("fstophi must be < nyquist (%f)", nyq)
        ! }

        !Get design amplitudes and adapt in case of zerophase
        astoplo=either(o_astoplo,.05,present(o_astoplo))
        apasslo=either(o_apasslo,.95,present(o_apasslo))
        apasshi=either(o_apasshi,.95,present(o_apasshi))
        astophi=either(o_astophi,.05,present(o_astophi))        
    !     if (astoplo > apasslo || apasshi < astophi)
    !             err("Bad amplitude parameters")
        
        !Adapt user frequencies if zerophase selected
        if (zerophase) then
            astoplo = sqrt(astoplo)
            apasslo = sqrt(apasslo)
            astophi = sqrt(astophi)
            apasshi = sqrt(apasshi)
        endif

        if (present(o_npoleslo)) then
            npoleslo = o_npoleslo
            f3dblo = either(o_f3dblo,.15*nyq,present(o_f3dblo)) *dt

        else !Use bdesign to make lo cut filters
            if (locut) call bfdesign(fpasslo,apasslo,fstoplo,astoplo,npoleslo,f3dblo)

        endif

        if (present(o_npoleshi)) then
            npolehi = o_npoleshi
            f3dbhi = either(o_f3dbhi,.40*nyq,present(o_f3dbhi)) *dt

        else !Use bdesign to make hi cut filters
            if (hicut) call bfdesign(fpasshi,apasshi,fstophi,astophi,npoleshi,f3dbhi)

        endif
        
        ! !Give verbose info if requested
        ! if (verbose && locut) {
        !         if (zerophase) {
        !                 warn("low-cut filter: npoles = %d, 3db point = %f(Hz)",
        !                         2*npoleslo, f3dblo/dt)
        !         } else {
        !                 warn("low-cut filter: npoles = %d, 3db point = %f(Hz)",
        !                         npoleslo, f3dblo/dt)
        !         }
        ! }
        ! if (verbose && hicut) {
        !         if (zerophase) {
        !                 warn("high-cut filter: npoles = %d, 3db point = %f(Hz)",
        !                         2*npoleshi, f3dbhi/dt)
        !         } else {
        !                 warn("high-cut filter: npoles = %d, 3db point = %f(Hz)",
        !                         npoleshi, f3dbhi/dt)
    
        !         }
        ! }

        !low-cut (high pass) filter
        if (locut) then
            do itr=1,ntr
                call bfhighpass(npoleslo,f3dblo,nt,trace(itr,:))
                if (zerophase) then
                    do i=1,nt/2 !reverse trace in place
                        tmp = trace(itr,i)
                        trace(itr,i) = trace(itr,nt-i)
                        trace(itr,nt-i) = tmp
                    enddo
                    call  bfhighpass(npoleslo,f3dblo,nt,trace(itr,:))
                    do i=0,nt/2-1 !flip trace back
                        tmp = trace(itr,i)
                        trace(itr,i) = trace(itr,nt-i)
                        trace(itr,nt-i) = tmp
                    enddo
                endif
            enddo
        endif

        !high-cut (low pass) filter
        if (hicut) then
            do itr=1,ntr
                call bflowpass(npoleshi,f3dbhi,nt,trace(itr,:))
                if (zerophase) then
                    do i=1,nt/2 !reverse trace
                        tmp = trace(itr,i)
                        trace(itr,i) = trace(itr,nt-i)
                        trace(itr,nt-i) = tmp
                    enddo
                    call bflowpass(npoleshi,f3dbhi,nt,trace(itr,:))
                    do i=1,nt/2 !flip trace back
                        tmp = trace(itr,i)
                        trace(itr,i) = trace(itr,nt-i)
                        trace(itr,nt-i) = tmp
                    enddo
                endif
            enddo
        endif

    end subroutine



! A Fortran version of butterworth.c in Seismic Unix
! dir: SeisUnix/src/cwp/lib/
!
! !! Copyright (c) Colorado School of Mines, 2011.!!
! !! All rights reserved.                       !!
! 
! !!********************** self documentation *********************!!
! !!****************************************************************************
! BUTTERWORTH - Functions to design and apply Butterworth filters:
! 
! bfdesign    design a Butterworth filter
! bfhighpass    apply a high-pass Butterworth filter 
! bflowpass    apply a low-pass Butterworth filter 
! 
! ******************************************************************************
! Function Prototypes:
! void bfhighpass (int npoles, float f3db, int n, float p[], float q[]);
! void bflowpass (int npoles, float f3db, int n, float p[], float q[]);
! void bfdesign (float fpass, float apass, float fstop, float astop,
!     int *npoles, float *f3db);
! 
! ******************************************************************************
! bfdesign:
! Input:
! fpass        frequency in pass band at which amplitude is >= apass
! apass        amplitude in pass band corresponding to frequency fpass
! fstop         frequency in stop band at which amplitude is <= astop
! astop        amplitude in stop band corresponding to frequency fstop
! 
! Output:
! npoles        number of poles
! f3db        frequency at which amplitude is sqrt(0.5) (-3 db)
! 
! bfhighpass and bflowpass:
! Input:
! npoles        number of poles (and zeros); npoles>=0 is required
! f3db        3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
! n        length of p and q
! p        array[n] to be filtered
! 
! Output:
! q        filtered array[n] (may be equivalent to p)
! 
! ******************************************************************************
! Notes:
! (1) Nyquist frequency equals 0.5
! 
! (2) The following conditions must be true:
!     (0.0<fpass && fpass<0.5) &&
!     (0.0<fstop && fstop<0.5) &&
!     (fpass!=fstop) &&
!     (0.0<astop && astop<apass && apass<1.0)
! 
! (3) if (fpass<fstop)
! 
! bfdesign:
! Butterworth filter:  compute number of poles and -3 db frequency
! for a low-pass or high-pass filter, given a frequency response
! constrained at two frequencies.
! 
! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! ****************************************************************************!!
! !!*************** end self doc *******************************!!
! 
! #include "cwp.h"
! 
! void
! bfdesign (float fpass, float apass, float fstop, float astop,
!     int *npoles, float *f3db)
! !!****************************************************************************
! Butterworth filter:  compute number of poles and -3 db frequency
! for a low-pass or high-pass filter, given a frequency response
! constrained at two frequencies.
! ******************************************************************************
! Input:
! fpass        frequency in pass band at which amplitude is >= apass
! apass        amplitude in pass band corresponding to frequency fpass
! fstop         frequency in stop band at which amplitude is <= astop
! astop        amplitude in stop band corresponding to frequency fstop
! 
! Output:
! npoles        number of poles
! f3db        frequency at which amplitude is sqrt(0.5) (-3 db)
! ******************************************************************************
! Notes:
! (1) Nyquist frequency equals 0.5
! 
! (2) The following conditions must be true:
!     (0.0<fpass && fpass<0.5) &&
!     (0.0<fstop && fstop<0.5) &&
!     (fpass!=fstop) &&
!     (0.0<astop && astop<apass && apass<1.0)
! 
! (3) if (fpass<fstop)
!         a low-pass filter is assumed
!     else
!         a high-pass filter is assumed
! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! ****************************************************************************!!
    subroutine bfdesign (fpass, apass, fstop, astop, npoles, f3db)
    
        !warp frequencies according to bilinear transform
        wpass = 2.0*tan(r_pi*fpass)
        wstop = 2.0*tan(r_pi*fstop)

        !if lowpass filter, then
        if (fstop>fpass) then
                fnpoles = log((1.0/(apass*apass)-1.0)/(1.0/(astop*astop)-1.0)) &
                        / log(wpass*wpass/wstop/wstop)
                w3db = wpass/((1.0/(apass*apass)-1.0)**(0.5/fnpoles))

        !else, if highpass filter, then
        else
                fnpoles = log((1.0/(apass*apass)-1.0)/(1.0/(astop*astop)-1.0)) &
                        / log(wstop*wstop/wpass/wpass)
                w3db = wpass*((1.0/(apass*apass)-1.0)**(0.5/fnpoles))
        endif

        !determine integer number of poles
        npoles = 1+int(fnpoles)

        !determine (unwarped) -3 db frequency
        f3db = atan(0.5*w3db)/r_pi
        
    end subroutine

! void
! bfhighpass (int npoles, float f3db, int n, float p[], float q[])
! !!****************************************************************************
! Butterworth filter:  high-pass
! ******************************************************************************
! Input:
! npoles        number of poles (and zeros); npoles>=0 is required
! f3db        3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
! n        length of p and q
! p        array[n] to be filtered
! 
! Output:
! q        filtered array[n] (may be equivalent to p)
! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! ****************************************************************************!!
    subroutine bfhighpass (npoles, f3db, n, trace)
        real,dimension(n) :: trace,tmp
        
        tmp=trace
        
        r = 2.0*tan(r_pi*abs(f3db))
        if (mod(npoles,2)/=0) then
            scale = r+2.0
            a = 2.0/scale
            b1 = (r-2.0)/scale
            pj = 0.0
            qjm1 = 0.0
            
            do j=1,n
                    pjm1 = pj
                    pj = tmp(j)
                    trace(j) = a*(pj-pjm1)-b1*qjm1
                    qjm1 = trace(j)
            enddo
        endif

        do jpair=0,npoles/2-1
            theta = r_pi*(2.*jpair+1.)/(2.*npoles)
            scale = 4.0+4.0*r*sin(theta)+r*r
            a = 4.0/scale
            b1 = (2.0*r*r-8.0)/scale
            b2 = (4.0-4.0*r*sin(theta)+r*r)/scale
            
            pjm1 = 0.
            pj = 0.
            qjm2 = 0.
            qjm1 = 0.
            do j=1,n
                pjm2 = pjm1
                pjm1 = pj
                pj = trace(j)
                trace(j) = a*(pj-2.*pjm1+pjm2)-b1*qjm1-b2*qjm2
                qjm2 = qjm1
                qjm1 = trace(j)
            enddo
        enddo
        
    end subroutine

! void
! bflowpass (int npoles, float f3db, int n, float p[], float q[])
! !!****************************************************************************
! Butterworth filter:  low-pass
! ******************************************************************************
! Input:
! npoles        number of poles (and zeros); npoles>=0 is required
! f3db        3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
! n        length of p and q
! p        array[n] to be filtered
! 
! Output:
! q        filtered array[n] (may be equivalent to p)
! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! ****************************************************************************!!
    subroutine bflowpass (npoles, f3db, n, trace)
        real,dimension(n) :: trace,tmp
        
        tmp=trace
        
        r = 2.0*tan(r_pi*abs(f3db))
        if (mod(npoles,2)/=0) then
            scale = r+2.0
            a = r/scale
            b1 = (r-2.0)/scale
            pj = 0.0
            qjm1 = 0.0
            
            do j=1,n
                pjm1 = pj
                pj = tmp(j)
                trace(j) = a*(pj+pjm1)-b1*qjm1
                qjm1 = trace(j)
            enddo
        endif
        
        do jpair=0,npoles/2-1
            theta = r_pi*(2.*jpair+1.)/(2.*npoles)
            scale = 4.+4.*r*sin(theta)+r*r
            a = r*r/scale
            b1 = (2.*r*r-8.)/scale
            b2 = (4.-4.*r*sin(theta)+r*r)/scale
            
            pjm1 = 0.
            pj = 0.
            qjm2 = 0.
            qjm1 = 0.
            
            do j=1,n
                pjm2 = pjm1
                pjm1 = pj
                pj = trace(j)
                trace(j) = a*(pj+2.*pjm1+pjm2)-b1*qjm1-b2*qjm2
                qjm2 = qjm1
                qjm1 = trace(j)
            enddo
            
        enddo
        
    end subroutine

end
