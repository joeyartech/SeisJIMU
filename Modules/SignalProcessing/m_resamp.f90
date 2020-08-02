module m_resamp

    private
    public :: resamp
    
    integer,parameter :: ltable=8, ntable=513
    real,dimension(0:ltable-1,0:ntable-1) :: table !sinc interpolation coeff

    contains

! A Fortran version of suresamp.c in Seismic Unix
! dir: SeisUnix/src/su/main/stretching_moveout_resamp
! reduced functionality
!
! !! Copyright (c) Colorado School of Mines, 2011.!!
! !! All rights reserved.                         !!
!
! !! SURESAMP: $Revision: 1.16 $ ; $Date: 2011/11/16 23:21:55 $      !!
!
! /*********************** self documentation **********************/
! "                                                                   ",
! " SURESAMP - Resample in time                                       ",
! "                                                                   ",
! " suresamp <stdin >stdout  [optional parameters]                    ",
! "                                                                   ",
! " Example 1: (assume original data had dt=.004 nt=256)              ",
! "    sufilter <data f=40,50 amps=1.,0. |                            ",
! "    suresamp nt=128 dt=.008 | ...                                  ",
! "                                                                   ",
! " Note the typical anti-alias filtering before sub-sampling!        ",
! "                                                                   ",
! " Example 2: (assume original data had dt=.004 nt=256)              ",
! "    suresamp <data nt=512 dt=.002 | ...                            ",
! "                                                                   ",
! " Example 3: (assume original data had d1=.1524 nt=8192)            ",
! "    sufilter <data f=0,1,3,3.28 amps=1,1,1,0 |                     ",
! "    suresamp <data nt=4096 dt=.3048 | ...                          ",
! "                                                                   ",
! " Example 4: (assume original data had d1=.5 nt=4096)               ",
! "    suresamp <data nt=8192 dt=.25 | ...                            ",
! "                                                                   ",
!
! !! Credits:
!  *    CWP: Dave (resamp algorithm), Jack (SU adaptation)
!  *
!  * Algorithm:
!  *    Resampling is done via 8-coefficient sinc-interpolation.
!  *    See "$CWPROOT/src/cwp/lib/intsinc8.c" for technical details.
!  *
!  */
! /************************ end self doc ***************************/
    subroutine resamp(ntr,o_fin, din, nin, datain,&
                         o_fout,dout,nout,dataout)
        real,dimension(nin,ntr)  :: datain
        real,dimension(nout,ntr) :: dataout
        real,optional :: o_fin, o_fout
        
        real :: fin, fout
        real,dimension(:),allocatable :: t
        
        fin=0.; fout=0.
        
        if(present(o_fin )) fin =o_fin
        if(present(o_fout)) fout=o_fout
        
        if(allocated(t)) deallocate(t)
        allocate(t(nout))
        t=fout+[(i,i=0,nout-1)]*dout
        
        !build sinc interpolation table
        call build_table
        
        do i=1,n2
            !interpolate using tabulated coefficients
            call interp(fin,din,nin,   datain(:,i),&
                                nout,t,dataout(:,i))
        enddo
        
    end subroutine


! A Fortran version of butterworth.c in Seismic Unix
! dir: SeisUnix/src/cwp/lib/
!
!! Copyright (c) Colorado School of Mines, 2011.!!
!! All rights reserved.                         !!
! 
! /*********************** self documentation **********************/
! /*****************************************************************************
! INTSINC8 - Functions to interpolate uniformly-sampled data via 8-coeff. sinc
! 		approximations:
! 
! ints8c	interpolation of a uniformly-sampled complex function y(x) via an
! 	         8-coefficient sinc approximation.
! ints8r	Interpolation of a uniformly-sampled real function y(x) via a
! 		table of 8-coefficient sinc approximations
! 
! ******************************************************************************
! Function Prototypes:
! void ints8c (int nxin, float dxin, float fxin, complex yin[], 
! 	complex yinl, complex yinr, int nxout, float xout[], complex yout[]);
! void ints8r (int nxin, float dxin, float fxin, float yin[], 
! 	float yinl, float yinr, int nxout, float xout[], float yout[]);
! 
! ******************************************************************************
! Input:
! nxin		number of x values at which y(x) is input
! dxin		x sampling interval for input y(x)
! fxin		x value of first sample input
! yin		array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
! yinl		value used to extrapolate yin values to left of yin[0]
! yinr		value used to extrapolate yin values to right of yin[nxin-1]
! nxout		number of x values a which y(x) is output
! xout		array[nxout] of x values at which y(x) is output
! 
! Output:
! yout		array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.
! 
! ******************************************************************************
! Notes:
! Because extrapolation of the input function y(x) is defined by the
! left and right values yinl and yinr, the xout values are not restricted
! to lie within the range of sample locations defined by nxin, dxin, and
! fxin.
! 
! The maximum error for frequiencies less than 0.6 nyquist is less than
! one percent.
! 
! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! *****************************************************************************/
! /**************** end self doc ********************************/
! 
! #include "cwp.h"
! 
! /* these are used by both ints8c and ints8r */
! #define LTABLE 8
! #define NTABLE 513
! 
! void ints8r (int nxin, float dxin, float fxin, float yin[], 
! 	float yinl, float yinr, int nxout, float xout[], float yout[])
! /*****************************************************************************
! Interpolation of a uniformly-sampled real function y(x) via a
! table of 8-coefficient sinc approximations; maximum error for frequiencies
! less than 0.6 nyquist is less than one percent.
! ******************************************************************************
! Input:
! nxin		number of x values at which y(x) is input
! dxin		x sampling interval for input y(x)
! fxin		x value of first sample input
! yin		array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
! yinl		value used to extrapolate yin values to left of yin[0]
! yinr		value used to extrapolate yin values to right of yin[nxin-1]
! nxout		number of x values a which y(x) is output
! xout		array[nxout] of x values at which y(x) is output
! 
! Output:
! yout		array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.
! ******************************************************************************
! Notes:
! Because extrapolation of the input function y(x) is defined by the
! left and right values yinl and yinr, the xout values are not restricted
! to lie within the range of sample locations defined by nxin, dxin, and
! fxin.
! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! *****************************************************************************/
    subroutine build_table
    
        do jtable = 1, ntable-2
            table(:,jtable) = mksinc(real(jtable)/real(ntable-1))
        enddo
        
        table(:,0) = 0.
        table(:,ntable-1) = 0.
        
        table(ltable/2-1,0) = 1.
        table(ltable/2,ntable-1) = 1.

    end subroutine


! A Fortran version of inttable.c in Seismic Unix
! dir: SeisUnix/src/cwp/lib
!
!! Copyright (c) Colorado School of Mines, 2011.!!
!! All rights reserved.                         !!
!
! !!********************** self documentation *********************!!
! !!****************************************************************************
! INTTABLE8 -  Interpolation of a uniformly-sampled complex function y(x)
! 		via a table of 8-coefficient interpolators
! 
! intt8c	interpolation of a uniformly-sampled complex function y(x)
! 		via a table of 8-coefficient interpolators
! intt8r	interpolation of a uniformly-sampled real function y(x) via a
! 		table of 8-coefficient interpolators
! 
! ******************************************************************************
! Function Prototype:
! void intt8c (int ntable, float table[][8],
! 	int nxin, float dxin, float fxin, complex yin[], 
! 	complex yinl, complex yinr, int nxout, float xout[], complex yout[]);
! void intt8r (int ntable, float table[][8],
! 	int nxin, float dxin, float fxin, float yin[], 
! 	float yinl, float yinr, int nxout, float xout[], float yout[]);
! 
! ******************************************************************************
! Input:
! ntable		number of tabulated interpolation operators; ntable>=2
! table		array of tabulated 8-point interpolation operators
! nxin		number of x values at which y(x) is input
! dxin		x sampling interval for input y(x)
! fxin		x value of first sample input
! yin		array of input y(x) values:  yin[0] = y(fxin), etc.
! yinl		value used to extrapolate yin values to left of yin[0]
! yinr		value used to extrapolate yin values to right of yin[nxin-1]
! nxout		number of x values a which y(x) is output
! xout		array of x values at which y(x) is output
! 
! Output:
! yout		array of output y(x) values:  yout[0] = y(xout[0]), etc.
! 
! ******************************************************************************
! NOTES:
! ntable must not be less than 2.
! 
! The table of interpolation operators must be as follows:
! 
! Let d be the distance, expressed as a fraction of dxin, from a particular
! xout value to the sampled location xin just to the left of xout.  Then,
! for d = 0.0,
! 
! table[0][0:7] = 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
! 
! are the weights applied to the 8 input samples nearest xout.
! Likewise, for d = 1.0,
! 
! table[ntable-1][0:7] = 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
! 
! are the weights applied to the 8 input samples nearest xout.  In general,
! for d = (float)itable/(float)(ntable-1), table[itable][0:7] are the
! weights applied to the 8 input samples nearest xout.  If the actual sample
! distance d does not exactly equal one of the values for which interpolators
! are tabulated, then the interpolator corresponding to the nearest value of
! d is used.
! 
! Because extrapolation of the input function y(x) is defined by the left
! and right values yinl and yinr, the xout values are not restricted to lie
! within the range of sample locations defined by nxin, dxin, and fxin.
! 
! ******************************************************************************
! AUTHOR:  Dave Hale, Colorado School of Mines, 06/02/89
! *****************************************************************************/
! /**************** end self doc ********************************/
! 
! void intt8r (int ntable, float table[][8],
! 	int nxin, float dxin, float fxin, float yin[], float yinl, float yinr,
! 	int nxout, float xout[], float yout[])
    subroutine interp(fxin,dxin,nxin, yin, &
                                nxout,xout,yout)
        real,dimension(0:nxin-1)  :: yin
        real,dimension(0:nxout-1) :: xout, yout
        
        real,dimension(0:ltable-1) :: pyin, ptable
        integer,parameter :: ioutb=-3-8
        real,parameter :: yinl=0., yinr=0.

        !compute constants
        xoutf = fxin
        xouts = 1.0/dxin
        xoutb = real(ltable)-xoutf*xouts
        fntablem1 = real(ntable-1)
        nxinm8 = nxin-ltable
        yin0 = yin(0)
        table00 = table(0,0)

        !loop over output samples
        do ixout=0,nxout-1

            !determine pointers into table and yin
            xoutn = xoutb + xout(ixout)*xouts
            ixoutn = int(xoutn)
            kyin = ioutb+ixoutn
            pyin = yin0+kyin
            frac = xoutn-real(ixoutn)
            if (frac>=0.) then
                ktable = frac*fntablem1+0.5
            else
                ktable = (frac+1.0)*fntablem1-0.5
            endif
            ptable = table00+ktable*8;
            
            !if totally within input array, use fast method
            if (kyin>=0 .and. kyin<=nxinm8) then
                yout(ixout) = &
                    pyin(0)*ptable(0)+ &
                    pyin(1)*ptable(1)+ &
                    pyin(2)*ptable(2)+ &
                    pyin(3)*ptable(3)+ &
                    pyin(4)*ptable(4)+ &
                    pyin(5)*ptable(5)+ &
                    pyin(6)*ptable(6)+ &
                    pyin(7)*ptable(7)
            
            else ! handle end effects with care
                    ! sum over 8 tabulated coefficients
                    sum=0.
                    do itable=0,ltable-1
                            if (kyin<0) then
                                yini = yinl
                            elseif (kyin>=nxin) then
                                yini = yinr
                            else
                                yini = yin(kyin)
                            endif
                            sum = sum + yini*ptable(itable)
                        kyin = kyin+1
                    enddo
                    yout(ixout) = sum
            
            endif
            
	enddo

    end subroutine


! A Fortran version of mksinc.c in Seismic Unix
! dir: SeisUnix/src/cwp/lib
!
! !! Copyright (c) Colorado School of Mines, 2011.!!
! !! All rights reserved.                         !!
! 
! !!********************** self documentation *********************!!
! !!****************************************************************************
! MKSINC - Compute least-squares optimal sinc interpolation coefficients.
! 
! mksinc		Compute least-squares optimal sinc interpolation coefficients.
! 
! ******************************************************************************
! Function Prototype:
! void mksinc (float d, int lsinc, float sinc[]);
! 
! ******************************************************************************
! Input:
! d		fractional distance to interpolation point; 0.0<=d<=1.0
! lsinc		length of sinc approximation; lsinc%2==0 and lsinc<=20
! 
! Output:
! sinc		array[lsinc] containing interpolation coefficients
! 
! ******************************************************************************
! Notes:
! The coefficients are a least-squares-best approximation to the ideal
! sinc function for frequencies from zero up to a computed maximum
! frequency.  For a given interpolator length, lsinc, mksinc computes
! the maximum frequency, fmax (expressed as a fraction of the nyquist
! frequency), using the following empirically derived relation (from
! a Western Geophysical Technical Memorandum by Ken Larner):
! 
! 	fmax = min(0.066+0.265*log(lsinc),1.0)
! 
! Note that fmax increases as lsinc increases, up to a maximum of 1.0.
! Use the coefficients to interpolate a uniformly-sampled function y(i) 
! as follows:
! 
!             lsinc-1
!     y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
!               j=0
! 
! Interpolation error is greatest for d=0.5, but for frequencies less
! than fmax, the error should be less than 1.0 percent.
! 
! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! *****************************************************************************/
! /**************** end self doc ********************************/
! 
! #include "cwp.h"
! 
! void mksinc (float d, int lsinc, float sinc[])
! /*****************************************************************************
! Compute least-squares optimal sinc interpolation coefficients.
! ******************************************************************************
! Input:
! d		fractional distance to interpolation point; 0.0<=d<=1.0
! lsinc		length of sinc approximation; lsinc%2==0 and lsinc<=20
! 
! Output:
! sinc		array[lsinc] containing interpolation coefficients
! ******************************************************************************
! Notes:
! The coefficients are a least-squares-best approximation to the ideal
! sinc function for frequencies from zero up to a computed maximum
! frequency.  For a given interpolator length, lsinc, mksinc computes
! the maximum frequency, fmax (expressed as a fraction of the nyquist
! frequency), using the following empirically derived relation (from
! a Western Geophysical Technical Memorandum by Ken Larner):
! 
! 	fmax = min(0.066+0.265*log(lsinc),1.0)
! 
! Note that fmax increases as lsinc increases, up to a maximum of 1.0.
! Use the coefficients to interpolate a uniformly-sampled function y(i) 
! as follows:
! 
!             lsinc-1
!     y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
!               j=0
! 
! Interpolation error is greatest for d=0.5, but for frequencies less
! than fmax, the error should be less than 1.0 percent.
! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! *****************************************************************************/
    function mksinc(d) result(s)
        real,dimension(0:ltable-1) :: s,a,c,work
        
        !compute auto-correlation and cross-correlation arrays
        fmax = 0.066+0.265*log(real(ltable))
        if(fmax>=1.) fmax=1.

        a = r_sinc(fmax*[(j,j=0,ltable-1)])
        c = r_sinc(fmax*(lsinc/2-[(j,j=0,ltable-1)]-1+d))

        !solve symmetric Toeplitz system for the sinc approximation
        call stoep(ltable,a,c,s,work)

    end function

end