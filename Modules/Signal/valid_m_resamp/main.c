/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SURESAMP: $Revision: 1.16 $ ; $Date: 2011/11/16 23:21:55 $        */

#include "su.h"
#include "segy.h"
#include "header.h"

#define NTABLE 513
#define LTABLE 8

void test_ints8r (int nxin, float dxin, float fxin, float yin[], 
    float yinl, float yinr, int nxout, float xout[], float yout[]);
void test_mksinc (float d, int lsinc, float sinc[]);
void test_intt8r (int ntable, float table[][8],
    int nxin, float dxin, float fxin, float yin[], float yinl, float yinr,
    int nxout, float xout[], float yout[]);

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                   ",
" SURESAMP - Resample in time                                       ",
"                                                                   ",
" suresamp <stdin >stdout  [optional parameters]                    ",
"                                                                   ",
" Required parameters:                                              ",
"     none                                                          ",
"                                                                   ",
" Optional Parameters:                                              ",
"    nt=tr.ns    number of time samples on output                   ",
"    dt=         time sampling interval on output                   ",
"                default is:                                        ",
"                tr.dt/10^6     seismic data                        ",
"                tr.d1          non-seismic data                    ",
"    tmin=       time of first sample in output                     ",
"                default is:                                        ",
"                tr.delrt/10^3  seismic data                        ",
"                tr.f1          non-seismic data                    ",
"    rf=         resampling factor;                                 ",
"                if defined, set nt=nt_in*rf and dt=dt_in/rf        ",
"    verbose=0   =1 for advisory messages                           ",
"                                                                   ",
"                                                                   ",
" Example 1: (assume original data had dt=.004 nt=256)              ",
"    sufilter <data f=40,50 amps=1.,0. |                            ",
"    suresamp nt=128 dt=.008 | ...                                  ",
" Using the resampling factor rf, this example translates to:       ",
"    sufilter <data f=40,50 amps=1.,0. | suresamp rf=0.5 | ...      ",
"                                                                   ",
" Note the typical anti-alias filtering before sub-sampling!        ",
"                                                                   ",
" Example 2: (assume original data had dt=.004 nt=256)              ",
"    suresamp <data nt=512 dt=.002 | ...                            ",
" or use:                                                           ",
"    suresamp <data rf=2 | ...                                      ",
"                                                                   ",
" Example 3: (assume original data had d1=.1524 nt=8192)            ",
"    sufilter <data f=0,1,3,3.28 amps=1,1,1,0 |                     ",
"    suresamp <data nt=4096 dt=.3048 | ...                          ",
"                                                                   ",
" Example 4: (assume original data had d1=.5 nt=4096)               ",
"    suresamp <data nt=8192 dt=.25 | ...                            ",
"                                                                   ",
NULL};

/* Credits:
 *    CWP: Dave (resamp algorithm), Jack (SU adaptation)
 *    CENPET: Werner M. Heigl - modified for well log support
 *    RISSC: Nils Maercklin 2006 - minor fixes, added rf option
 *
 * Algorithm:
 *    Resampling is done via 8-coefficient sinc-interpolation.
 *    See "$CWPROOT/src/cwp/lib/intsinc8.c" for technical details.
 *
 * Trace header fields accessed:  ns, dt, delrt, d1, f1, trid
 * Trace header fields modified:  ns, dt, delrt (only when set tmin)
 *                                d1, f1 (only when set tmin)
 */
/************************ end self doc ***************************/


segy intrace, outtrace;

int
main(int argc, char **argv)
{
    int nt;            /* number of samples on output trace */
    int nt_in;         /* ... on input trace */
    float dt;          /* sample rate on output trace */
    int idt;           /* ... as integer */
    float dt_in;       /* ... on input trace */
    float tmin;        /* first time sample on output trace */
    float tmin_in;     /* ... on input trace */
    float *t;          /* array of output times */
    int tmin_is_set=0; /* flag for user-defined tmin */
    float rf;          /* resampling factor */
                       /* (rf>1 means down- and rf<1 up-sampling) */ 
    int verbose;       /* if 1(yes) display advisory messages */
    cwp_Bool seismic;  /* flag: is this seismic data? */
 
    
    /* Hook up getpar */
    initargs(argc, argv);
    requestdoc(1);

    /* Get verbose parameter */
    if (!getparint("verbose", &verbose)) verbose = 0;

    /* Get information from first trace */
    if (!gettr(&intrace)) err("can't get first trace");
    nt_in = intrace.ns;
    if (!nt_in)  err("ns not set in header");

    /* check for seismic or well log data */
    seismic = ISSEISMIC(intrace.trid);        
    if (seismic) {
        if (verbose)
            warn("input is seismic data, trid=%d",intrace.trid);
        dt_in   = ((double) intrace.dt)/1000000.0;
        tmin_in = ((double) intrace.delrt)/1000.0;  
    }
    else {
        if (verbose)
            warn("input is not seismic data, trid=%d",intrace.trid);
        dt_in   = intrace.d1;
        tmin_in = intrace.f1;
    }
    
    /* check input times */
    if (!dt_in) getparfloat("dt_in", &dt_in);
    if (!dt_in) err("dt or d1 not set in header or not getparred");
    if (!tmin_in && verbose) warn("delrt or f1 not set in header");

    /* Get parameters */
    if (!getparfloat("rf", &rf)) rf=0.0;
    if (rf<0.0) err("factor rf=%g must be positive", rf);

    if (rf) {
        if (!getparint("nt", &nt)) nt = NINT( ((float)nt_in)*rf);
        if (!getparfloat("dt", &dt)) dt = dt_in/rf;
    }
    else {
        if (!getparint("nt", &nt)) nt = nt_in;
        if (!getparfloat("dt", &dt)) dt = dt_in;
    }
    if (getparfloat("tmin", &tmin)) tmin_is_set = 1;
    checkpars();

    /* Validate user input nt and dt */
    CHECK_NT("nt",nt);
    idt = NINT(dt * 1000000.0);
        
    /* Allocate vector of output times */
    t = ealloc1float(nt);

    /* Loop over traces */    
    do {
        if (!tmin_is_set)    tmin = tmin_in;
            
        /* Compute output times */
        { register int itime;
          register float tvalue;
          for (itime=0,tvalue=tmin; itime<nt; itime++,tvalue+=dt)
              t[itime] = tvalue;
        }
    
        /* copy and adjust header */
        memcpy(&outtrace, &intrace, HDRBYTES);
        outtrace.ns = nt;
        if (seismic)
            outtrace.dt = idt;
        else
            outtrace.d1 = dt;

        if (tmin_is_set) {
            if (seismic)
                outtrace.delrt = NINT(tmin * 1000.0);
            else
                outtrace.f1 = tmin;
        }
        
        /* sinc interpolate new data */
        test_ints8r(nt_in, dt_in, tmin_in, intrace.data, 
                0.0, 0.0, nt, t, outtrace.data);
        
        puttr(&outtrace);
    } while (gettr(&intrace));


    return(CWP_Exit());
}


void test_ints8r (int nxin, float dxin, float fxin, float yin[], 
    float yinl, float yinr, int nxout, float xout[], float yout[])
/*****************************************************************************
Interpolation of a uniformly-sampled real function y(x) via a
table of 8-coefficient sinc approximations; maximum error for frequiencies
less than 0.6 nyquist is less than one percent.
******************************************************************************
Input:
nxin        number of x values at which y(x) is input
dxin        x sampling interval for input y(x)
fxin        x value of first sample input
yin     array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
yinl        value used to extrapolate yin values to left of yin[0]
yinr        value used to extrapolate yin values to right of yin[nxin-1]
nxout       number of x values a which y(x) is output
xout        array[nxout] of x values at which y(x) is output

Output:
yout        array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.
******************************************************************************
Notes:
Because extrapolation of the input function y(x) is defined by the
left and right values yinl and yinr, the xout values are not restricted
to lie within the range of sample locations defined by nxin, dxin, and
fxin.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/

{
    static float table[NTABLE][LTABLE];
    static int tabled=0;
    int jtable;
    float frac;

    /* tabulate sinc interpolation coefficients if not already tabulated */
    if (!tabled) {
        for (jtable=1; jtable<NTABLE-1; jtable++) {
            frac = (float)jtable/(float)(NTABLE-1);
            test_mksinc(frac,LTABLE,&table[jtable][0]);
        }

        for (jtable=0; jtable<LTABLE; jtable++) {
            table[0][jtable] = 0.0;
            table[NTABLE-1][jtable] = 0.0;
        }

        table[0][LTABLE/2-1] = 1.0;
        table[NTABLE-1][LTABLE/2] = 1.0;
        tabled = 1;
    }

    /* interpolate using tabulated coefficients */
    test_intt8r(NTABLE,table,nxin,dxin,fxin,yin,yinl,yinr,nxout,xout,yout);

}

void test_mksinc (float d, int lsinc, float sinc[])
{
        int j;
        double s[20],a[20],c[20],work[20],fmax;

        /* compute auto-correlation and cross-correlation arrays */
        fmax = 0.066+0.265*log((double)lsinc);
        fmax = (fmax<1.0)?fmax:1.0;
        for (j=0; j<lsinc; j++) {
                a[j] = dsinc(fmax*j);
                c[j] = dsinc(fmax*(lsinc/2-j-1+d));
        }

        /* solve symmetric Toeplitz system for the sinc approximation */
        stoepd(lsinc,a,c,s,work);
        for (j=0; j<lsinc; j++)
                sinc[j] = s[j];
}


void test_intt8r (int ntable, float table[][8],
    int nxin, float dxin, float fxin, float yin[], float yinl, float yinr,
    int nxout, float xout[], float yout[])
{
    int ioutb,nxinm8,ixout,ixoutn,kyin,ktable,itable;
    float xoutb,xoutf,xouts,xoutn,frac,fntablem1,yini,sum,
        *yin0,*table00,*pyin,*ptable;

    /* compute constants */
    ioutb = -3-8;
    xoutf = fxin;
    xouts = 1.0/dxin;
    xoutb = 8.0-xoutf*xouts;
    fntablem1 = (float)(ntable-1);
    nxinm8 = nxin-8;
    yin0 = &yin[0];
    table00 = &table[0][0];

    /* loop over output samples */
    for (ixout=0; ixout<nxout; ixout++) {

        /* determine pointers into table and yin */
        xoutn = xoutb+xout[ixout]*xouts;
        ixoutn = (int)xoutn;
        kyin = ioutb+ixoutn;
        pyin = yin0+kyin;
        frac = xoutn-(float)ixoutn;
        ktable = frac>=0.0?frac*fntablem1+0.5:(frac+1.0)*fntablem1-0.5;
        ptable = table00+ktable*8;
        
        /* if totally within input array, use fast method */
        if (kyin>=0 && kyin<=nxinm8) {
            yout[ixout] = 
                pyin[0]*ptable[0]+
                pyin[1]*ptable[1]+
                pyin[2]*ptable[2]+
                pyin[3]*ptable[3]+
                pyin[4]*ptable[4]+
                pyin[5]*ptable[5]+
                pyin[6]*ptable[6]+
                pyin[7]*ptable[7];
        
        /* else handle end effects with care */
        } else {
    
            /* sum over 8 tabulated coefficients */
            for (itable=0,sum=0.0; itable<8; itable++,kyin++) {
                if (kyin<0)
                    yini = yinl;
                else if (kyin>=nxin)
                    yini = yinr;
                else
                    yini = yin[kyin];
                sum += yini*(*ptable++);
            }
            yout[ixout] = sum;
        }
    }
}