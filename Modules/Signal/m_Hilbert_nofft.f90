module m_hilbert_nofft
use m_System
use m_math

    private
    public :: hilbert_nofft

    integer,parameter :: LHHALF=50   !/* half-length of Hilbert transform filter*/
    integer,parameter :: LH=2*LHHALF+1   !/* filter length must be odd */

    real,dimension(0:LH-1) :: h

    logical :: if_madeh=.false.

    contains

! A Fortran version of hilbert.c in Seismic Unix
! dir: SeisUnix/src/cwp/lib/

! /* Copyright (c) Colorado School of Mines, 2011.*/
! /* All rights reserved.                       */

! /*********************** self documentation **********************/
! /*****************************************************************************
! HILBERT - Compute Hilbert transform y of x

! hilbert     compute the Hilbert transform

! ******************************************************************************
! Function Prototype:
! void hilbert (int n, float x[], float y[]);

! ******************************************************************************
! Input:
! n       length of x and y
! x       array[n] to be Hilbert transformed

! Output:
! y       array[n] containing Hilbert transform of x

! ******************************************************************************
! Notes:
! The Hilbert transform is computed by convolving x with a
! windowed (approximate) version of the ideal Hilbert transformer.

! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! *****************************************************************************/
! /**************** end self doc ********************************/
    subroutine hilbert_nofft(algo,din,dout,nt,ntr)
        character(*) :: algo
        real,dimension(nt,ntr) :: din
        real,dimension(nt,ntr) :: dout
    
        !/* if not made h, make Hilbert transform filter; use Hamming window */
        if (.not.if_madeh) then
! print*,'if_madeh',if_madeh
            h(LHHALF) = 0.0
            do i=1,LHHALF
                taper = 0.54+0.46*cos(r_pi*i/LHHALF)
                h(LHHALF+i) = taper*(mod(i,2)*2.0/(r_pi*i))
                h(LHHALF-i) = -h(LHHALF+i)
            enddo
            if_madeh = .true.
        endif

        !/* convolve Hilbert transform with input array */
        select case(algo)
        case ('naive')
        do itr=1,ntr
            call convolve_cwp_naive(LH,-LHHALF,h, nt,0,din(:,itr), nt,0,dout(:,itr))
        enddo

        case ('generic')
        do itr=1,ntr
            !call convolve_cwp_short(LH,-LHHALF,h, nt,0,din(:,itr), nt,0,dout(:,itr))
            call convolve_cwp_generic(LH,-LHHALF,h, nt,0,din(:,itr), nt,0,dout(:,itr))
        enddo

        end select

    end subroutine


! A Fortran version of convolution.c in Seismic Unix
! dir: SeisUnix/src/cwp/lib/

! /* Copyright (c) Colorado School of Mines, 2011.*/
! /* All rights reserved.                       */

! /*********************** self documentation **********************/
! /*****************************************************************************
! CONVOLUTION - Compute z = x convolved with y

! convolve_cwp    compute the convolution of two input vector arrays

! ******************************************************************************
! Input:
! lx      length of x array
! ifx     sample index of first x
! x       array[lx] to be convolved with y
! ly      length of y array
! ify     sample index of first y
! y       array[ly] with which x is to be convolved
! lz      length of z array
! ifz     sample index of first z

! Output:
! z       array[lz] containing x convolved with y

! ******************************************************************************
! Function Prototype:
! void convolve_cwp (int lx, int ifx, float *x, int ly, int ify, float *y,
!     int lz, int ifz, float *z);

! ******************************************************************************
! Notes:
! The operation z = x convolved with y is defined to be
!            ifx+lx-1
!     z[i] =   sum    x[j]*y[i-j]  ;  i = ifz,...,ifz+lz-1
!             j=ifx
! The x samples are contained in x[0], x[1], ..., x[lx-1]; likewise for
! the y and z samples.  The sample indices of the first x, y, and z values
! determine the location of the origin for each array.  For example, if
! z is to be a weighted average of the nearest 5 samples of y, one might
! use 
!     ...
!     x[0] = x[1] = x[2] = x[3] = x[4] = 1.0/5.0;
!     conv(5,-2,x,lx,0,y,ly,0,z);
!     ...
! In this example, the filter x is symmetric, with index of first sample = -2.

! This function is optimized for architectures that can simultaneously perform
! a multiply, add, and one load from memory; e.g., the IBM RISC System/6000.
! Because, for each value of i, it accumulates the convolution sum z[i] in a
! scalar, this function is not likely to be optimal for vector architectures.

! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 11/23/91
! *****************************************************************************/
! /**************** end self doc ********************************/

    subroutine convolve_cwp_naive(lx,ifx,x, ly,ify,y, lz,ifz,z)
        real,dimension(ifx:ifx+lx-1),intent(in)  :: x
        real,dimension(ify:ify+ly-1),intent(in)  :: y
        real,dimension(ifz:ifz+lz-1),intent(out) :: z

        ilx=ifx+lx-1
        ily=ify+ly-1
        ilz=ifz+lz-1
        
        ! do i=ifz,ilz
        !     jlow = i-ily;  if (jlow<ifx) jlow = ifx
        !     jhigh = i-ify;  if (jhigh>ilx) jhigh = ilx
        !     sum=0.
        !     do j=jlow,jhigh
        !         sum=sum+x(j)*y(i-j)
        !     enddo
        !     z(i)=sum
        ! enddo

        do i=ifz,ilz
            jlow = max(i-ily,ifx)
            jhigh = min(i-ify,ilx)
            z(i)=sum(x(jlow:jhigh)*y(i-jlow:i-jhigh:-1))
        enddo
        
    end subroutine

    ! !a vectorized version dedicated to summation
    ! !however, this is not speeding up
    ! !and is bugged when looping over traces (do itr=1,ntr)
    ! !(result of 1st trace is correct, but becomes zeros starting from 2nd trace)
    ! subroutine convolve_cwp_generic(lx,ifx,x, ly,ify,y, lz,ifz,z)
    !     real,dimension(ifx:ifx+lx-1),intent(in)  :: x
    !     real,dimension(0:ly-1),intent(in)  :: y
    !     real,dimension(0:lz-1),intent(out) :: z

    !     ilx=ifx+lx-1
    !     ily=ify+ly-1
    !     ilz=ifz+lz-1
        
    !     ! /* if x is longer than y, swap x and y */
    !     ! if (lx>ly) {
    !     !     i = ifx;  ifx = ify;  ify = i;
    !     !     i = ilx;  ilx = ily;  ily = i;
    !     !     i = lx;  lx = ly;  ly = i;
    !     !     t = x;  x = y;  y = t;
    !     ! }
        
    !     !/* OFF LEFT:  i < ify+ifx */
        
    !     !/* zero output for all i */
    !     ilow = ifz;
    !     ihigh = min(ify+ifx-1,ilz)
    !     z(ilow:ihigh)=0.

    !     !/* ROLLING ON:  ify+ifx <= i < ify+ilx */
        
    !     !/* if necessary, do one i so that number of j in overlap is odd */
    !     if (i<ify+ilx .and. i<=ilz) then
    !         jlow = ifx;
    !         jhigh = i-ify;
    !         if (mod(jhigh-jlow,2)==1) then
    !             z(i) = sum(x(jlow:jhigh)*y(i-jlow:i-jhigh:-1))

    !             i=i+1
    !         endif
    !     endif
        
    !     !/* loop over pairs of i and j */
    !     ilow = i;
    !     ihigh = min(ilx+ify-1,ilz)
    !     jlow=ifx
    !     jhigh = ilow-ify;

    !     do i=ilow,ihigh-1,2
    !         sa = 0.; sb = 0.
    !         xb = x(jhigh+1)
    !         yb = 0.
    !         !do j=jhigh,jlow,-2 this might be a bug in the SU code; otherwise x(j-1) in below will be over-bounded
    !         do j=jhigh,jlow+1,-2
    !             sa = sa+ xb*yb;
    !             ya = y(i-j);
    !             sb = sb+ xb*ya;
    !             xa = x(j);
    !             sa = sa+ xa*ya;
    !             yb = y(i+1-j);
    !             sb = sb+ xa*yb;
    !             xb = x(j-1);
    !         enddo
    !         z(i) = sa;
    !         z(i+1) = sb;

    !         jhigh=jhigh+2
    !     enddo
        
    !     !/* if number of i is odd */
    !     if (i==ihigh) then
    !         jlow = ifx
    !         jhigh = i-ify
    !         z(i) = sum(x(jlow:jhigh)*y(i-jlow:i-jhigh:-1))

    !         i=i+1
    !     endif
        
    !     !/* MIDDLE:  ify+ilx <= i <= ily+ifx */
        
    !     !/* determine limits for i and j */
    !     ilow = i;
    !     ihigh = min(ily+ifx,ilz)
    !     jlow = ifx;
    !     jhigh = ilx;
                
    !     if (mod(jhigh-jlow,2)==1) then !/* if number of j is even, do j in pairs with no leftover */
    !         do i=ilow,ihigh-1,2
    !             sa = 0.; sb = 0.
    !             yb = y(i+1-jlow)
    !             xa = x(jlow)
    !             do j=jlow,jhigh-1,2
    !                 sb = sb+ xa*yb;
    !                 ya = y(i-j);
    !                 sa = sa+ xa*ya;
    !                 xb = x(j+1);
    !                 sb = sb+ xb*ya;
    !                 yb = y(i-1-j);
    !                 sa = sa+ xb*yb;
    !                 xa = x(j+2);
    !             enddo
    !             z(i) = sa
    !             z(i+1) = sb
    !         enddo
                
    !     else !/* else, number of j is odd, so do j in pairs with leftover */
    !         do i=ilow,ihigh-1,2
    !             sa = 0.; sb = 0.
    !             yb = y(i+1-jlow)
    !             xa = x(jlow)
    !             do j=jlow,jhigh-1,2
    !                 sb = sb+ xa*yb;
    !                 ya = y(i-j);
    !                 sa = sa+ xa*ya;
    !                 xb = x(j+1);
    !                 sb = sb+ xb*ya;
    !                 yb = y(i-1-j);
    !                 sa = sa+ xb*yb;
    !                 xa = x(j+2);
    !             enddo
    !             z(i) = sa+x(jhigh)*y(i-jhigh)
    !             z(i+1) = sb+x(jhigh)*y(i+1-jhigh)
    !         enddo

    !     endif
        
    !     !/* if number of i is odd */
    !     if (i==ihigh) then
    !         z(i) = sum(x(jlow:jhigh)*y(i-jlow:i-jhigh:-1))

    !         i=i+1
    !     endif

    !     ! /* ROLLING OFF:  ily+ifx < i <= ily+ilx */

    !     ! /* if necessary, do one i so that number of j in overlap is even */
    !     if (i<=ily+ilx .and. i<=ilz) then
    !         jlow = i-ily;
    !         jhigh = ilx;
    !         if (mod((jhigh-jlow),2)==0) then
    !             z(i) = sum(x(jlow:jhigh)*y(i-jlow:i-jhigh:-1))
    !             i=i+1
    !         endif
    !     endif

    !     ! /* number of j is now even, so loop over both i and j in pairs */
    !     ilow = i;
    !     ihigh = min(ily+ilx,ilz)
    !     jlow = ilow-ily;
    !     jhigh = ilx-2; !/* Dave's new patch */
    !     do i=ilow,ihigh-1,2
    !         sa = 0.; sb = 0.;
    !         xa = x(jlow);
    !         yb = 0.;
    !         do j=jlow,jhigh-1,2
    !             sb = sb+ xa*yb;
    !             ya = y(i-j);
    !             sa = sa+ xa*ya;
    !             xb = x(j+1);
    !             sb = sb+ xb*ya;
    !             yb = y(i-1-j);
    !             sa = sa+ xb*yb;
    !             xa = x(j+2);
    !         enddo
    !         sb = sb+ xa*yb;
    !         ya = y(i-j);
    !         sa = sa+ xa*ya;
    !         xb = x(j+1);
    !         sb = sb+ xb*ya;
    !         yb = y(i-1-j);
    !         sa = sa+ xb*yb;
    !         z(i) = sa;
    !         z(i+1) = sb;

    !         jlow=jlow+2
    !     enddo

    !     ! /* if number of i is odd */
    !     if (i==ihigh) then
    !         jlow = i-ily;
    !         jhigh = ilx;
    !         z(i) = sum(x(jlow:jhigh)*y(i-jlow:i-jhigh:-1))

    !         i=i+1
    !     endif
        
    !     ! /* OFF RIGHT:  ily+ilx < i */
        
    !     ! /* zero output for all i */
    !     ilow = i;
    !     ihigh = ilz;

    !     z(ilow:ihigh)=0.
        
    ! end subroutine


    !non-vectorized version
    subroutine convolve_cwp_generic(lx,ifx,x, ly,ify,y, lz,ifz,z)
        ! real,dimension(ifx:ifx+lx-1),intent(in)  :: x
        ! real,dimension(ify:ify+ly-1),intent(in)  :: y
        ! real,dimension(ifz:ifz+lz-1),intent(out) :: z
        real,dimension(ifx:ifx+lx-1),intent(in)  :: x
        real,dimension(0:ly-1),intent(in)  :: y
        real,dimension(0:lz-1),intent(out) :: z

        ilx=ifx+lx-1
        ily=ify+ly-1
        ilz=ifz+lz-1
        
        ! /* if x is longer than y, swap x and y */
        ! if (lx>ly) {
        !     i = ifx;  ifx = ify;  ify = i;
        !     i = ilx;  ilx = ily;  ily = i;
        !     i = lx;  lx = ly;  ly = i;
        !     t = x;  x = y;  y = t;
        ! }
                     
        !/* OFF LEFT:  i < ify+ifx */
        
        !/* zero output for all i */
        ilow = ifz;
        ihigh = ify+ifx-1;  if (ihigh>ilz) ihigh = ilz;
        do i=ilow,ihigh; z(i)=0; enddo
        
        !/* ROLLING ON:  ify+ifx <= i < ify+ilx */
        
        !/* if necessary, do one i so that number of j in overlap is odd */
        if (i<ify+ilx .and. i<=ilz) then
            jlow = ifx;
            jhigh = i-ify;
            if (mod(jhigh-jlow,2)==1) then
                sa = 0.0;
                do j=jlow,jhigh; sa = sa+ x(j)*y(i-j); enddo
                z(i) = sa;
                i=i+1
            endif
        endif

        !/* loop over pairs of i and j */
        ilow = i;
        ihigh = ilx+ify-1;  if (ihigh>ilz) ihigh = ilz;
        jlow=ifx
        jhigh = ilow-ify;

        do i=ilow,ihigh-1,2
            sa = 0.0; sb = 0.0;
            xb = x(jhigh+1);
            yb = 0.0;
            !do j=jhigh,jlow,-2 this might be a bug in the SU code; otherwise x(j-1) in below will be over-bounded
            do j=jhigh,jlow+1,-2 
                sa = sa+ xb*yb;
                ya = y(i-j);
                sb = sb+ xb*ya;
                xa = x(j);
                sa = sa+ xa*ya;
                yb = y(i+1-j);
                sb = sb+ xa*yb;
                xb = x(j-1);
            enddo
            z(i) = sa;
            z(i+1) = sb;

            jhigh=jhigh+2
        enddo
        
        !/* if number of i is odd */
        if (i==ihigh) then
            jlow = ifx;
            jhigh = i-ify;
            sa = 0.0;
            do j=jlow,jhigh
                sa = sa+ x(j)*y(i-j);
            enddo                
            z(i) = sa;

            i=i+1
        endif
        
        !/* MIDDLE:  ify+ilx <= i <= ily+ifx */
        
        !/* determine limits for i and j */
        ilow = i;
        ihigh = ily+ifx;  if (ihigh>ilz) ihigh = ilz;
        jlow = ifx;
        jhigh = ilx;
                
        if (mod(jhigh-jlow,2)==1) then !/* if number of j is even, do j in pairs with no leftover */
            do i=ilow,ihigh-1,2
                sa = 0.; sb = 0.0;
                yb = y(i+1-jlow);
                xa = x(jlow);
                do j=jlow,jhigh-1,2
                    sb = sb+ xa*yb;
                    ya = y(i-j);
                    sa = sa+ xa*ya;
                    xb = x(j+1);
                    sb = sb+ xb*ya;
                    yb = y(i-1-j);
                    sa = sa+ xb*yb;
                    xa = x(j+2);
                enddo
                z(i) = sa;
                z(i+1) = sb;
            enddo
                
        else !/* else, number of j is odd, so do j in pairs with leftover */
            do i=ilow,ihigh-1,2
                sa = 0; sb = 0.0;
                yb = y(i+1-jlow);
                xa = x(jlow);
                do j=jlow,jhigh-1,2
                    sb = sb+ xa*yb;
                    ya = y(i-j);
                    sa = sa+ xa*ya;
                    xb = x(j+1);
                    sb = sb+ xb*ya;
                    yb = y(i-1-j);
                    sa = sa+ xb*yb;
                    xa = x(j+2);
                enddo
                z(i) = sa+x(jhigh)*y(i-jhigh);
                z(i+1) = sb+x(jhigh)*y(i+1-jhigh);
            enddo
        endif
        
        !/* if number of i is odd */
        if (i==ihigh) then
            sa = 0.0;
            do j=jlow,jhigh; sa = sa+ x(j)*y(i-j); enddo
            z(i) = sa

            i=i+1
        endif

        !/* ROLLING OFF:  ily+ifx < i <= ily+ilx */
        
        !/* if necessary, do one i so that number of j in overlap is even */
        if (i<=ily+ilx .and. i<=ilz) then
            jlow = i-ily;
            jhigh = ilx;
            if (mod((jhigh-jlow),2)==0) then
                sa = 0.0;
                do j=jlow,jhigh
                    sa = sa+ x(j)*y(i-j);
                enddo
                z(i) = sa

                i=i+1
            endif
        endif

        ! /* number of j is now even, so loop over both i and j in pairs */
        ilow = i;
        ihigh = ily+ilx;  if (ihigh>ilz) ihigh = ilz;
        jlow = ilow-ily;
        jhigh = ilx-2; !/* Dave's new patch */
            do i=ilow,ihigh-1,2
                sa = 0.; sb = 0.;
                xa = x(jlow);
                yb = 0.0;
                do j=jlow,jhigh-1,2
                    sb = sb +xa*yb;
                    ya = y(i-j);
                    sa = sa+ xa*ya;
                    xb = x(j+1);
                    sb = sb+ xb*ya;
                    yb = y(i-1-j);
                    sa = sa+ xb*yb;
                    xa = x(j+2);
                enddo
                sb = sb+ xa*yb;
                ya = y(i-j);
                sa = sa+ xa*ya;
                xb = x(j+1);
                sb = sb+ xb*ya;
                yb = y(i-1-j);
                sa = sa+ xb*yb;
                z(i) = sa;
                z(i+1) = sb;

                jlow=jlow+2
            enddo

        !/* if number of i is odd */
        if (i==ihigh) then
            jlow = i-ily;
            jhigh = ilx;
            sa = 0.0;
            do j=jlow,jhigh
                sa = sa+ x(j)*y(i-j);
            enddo
            z(i) = sa

            i=i+1
        endif

        !/* OFF RIGHT:  ily+ilx < i */
        
        !/* zero output for all i */
        ilow = i;
        ihigh = ilz;

        !do i=ilow,ihigh; z(i)=0.; enddo
        z(ilow:ihigh)=0.
        
    end subroutine


    !/* function optimized for short x */
    subroutine convolve_cwp_short(lx,ifx,x, ly,ify,y, lz,ifz,z)
        real,dimension(ifx:ifx+lx-1),intent(in)  :: x
        real,dimension(0:ly-1),intent(in)  :: y
        real,dimension(0:lz-1),intent(out) :: z

        ilx=ifx+lx-1
        ily=ify+ly-1
        ilz=ifz+lz-1
            
        !/* OFF LEFT:  i < ifx+ify */
        ilow = ifz
        ihigh = ify+ifx-1;  if (ihigh>ilz) ihigh = ilz
        z(ilow:ihigh)=0.
        
        !/* ROLLING ON:  ify+ifx <= i < ify+ilx */
        ilow = ify+ifx;  if (ilow<ifz) ilow = ifz
        ihigh = ify+ilx-1;  if (ihigh>ilz) ihigh = ilz
        jlow = ifx
        jhigh = ilow-ify

        do i=ilow,ihigh
            sum=0.
            do j=jlow,jhigh
                sum=sum+x(j)*y(i-j)
            enddo
            z(i)=sum
            jhigh=jhigh+1
        enddo

        !/* MIDDLE:  ify+ilx <= i <= ily+ifx */
        ilow  = ify+ilx;  if (ilow<ifz) ilow = ifz
        ihigh = ily+ifx;  if (ihigh>ilz) ihigh = ilz
        ! if (lx==1) {
        !     x0 = x[ifx];
        !     for (i=ilow; i<=ihigh-1; i+=2) {
        !         ya = y[i+1-ifx];  z1 = x0*ya;
        !         yb = y[i-ifx];  z0 = x0*yb;
        !         z[i+1] = z1;
        !         z[i] = z0;
        !     } !...
        ! } else if (lx==30) {
            do i=ilow,ihigh-1,2
                ya = y(i+1-ifx );  z1 =     x(ifx   )*ya;
                yb = y(i-ifx   );  z0 =     x(ifx   )*yb;  z1 = z1+ x(ifx+1 )*yb;
                ya = y(i-ifx-1 );  z0 = z0+ x(ifx+1 )*ya;  z1 = z1+ x(ifx+2 )*ya;
                yb = y(i-ifx-2 );  z0 = z0+ x(ifx+2 )*yb;  z1 = z1+ x(ifx+3 )*yb;
                ya = y(i-ifx-3 );  z0 = z0+ x(ifx+3 )*ya;  z1 = z1+ x(ifx+4 )*ya;
                yb = y(i-ifx-4 );  z0 = z0+ x(ifx+4 )*yb;  z1 = z1+ x(ifx+5 )*yb;
                ya = y(i-ifx-5 );  z0 = z0+ x(ifx+5 )*ya;  z1 = z1+ x(ifx+6 )*ya;
                yb = y(i-ifx-6 );  z0 = z0+ x(ifx+6 )*yb;  z1 = z1+ x(ifx+7 )*yb;
                ya = y(i-ifx-7 );  z0 = z0+ x(ifx+7 )*ya;  z1 = z1+ x(ifx+8 )*ya;
                yb = y(i-ifx-8 );  z0 = z0+ x(ifx+8 )*yb;  z1 = z1+ x(ifx+9 )*yb;
                ya = y(i-ifx-9 );  z0 = z0+ x(ifx+9 )*ya;  z1 = z1+ x(ifx+10)*ya;
                yb = y(i-ifx-10);  z0 = z0+ x(ifx+10)*yb;  z1 = z1+ x(ifx+11)*yb;
                ya = y(i-ifx-11);  z0 = z0+ x(ifx+11)*ya;  z1 = z1+ x(ifx+12)*ya;
                yb = y(i-ifx-12);  z0 = z0+ x(ifx+12)*yb;  z1 = z1+ x(ifx+13)*yb;
                ya = y(i-ifx-13);  z0 = z0+ x(ifx+13)*ya;  z1 = z1+ x(ifx+14)*ya;
                yb = y(i-ifx-14);  z0 = z0+ x(ifx+14)*yb;  z1 = z1+ x(ifx+15)*yb;
                ya = y(i-ifx-15);  z0 = z0+ x(ifx+15)*ya;  z1 = z1+ x(ifx+16)*ya;
                yb = y(i-ifx-16);  z0 = z0+ x(ifx+16)*yb;  z1 = z1+ x(ifx+17)*yb;
                ya = y(i-ifx-17);  z0 = z0+ x(ifx+17)*ya;  z1 = z1+ x(ifx+18)*ya;
                yb = y(i-ifx-18);  z0 = z0+ x(ifx+18)*yb;  z1 = z1+ x(ifx+19)*yb;
                ya = y(i-ifx-19);  z0 = z0+ x(ifx+19)*ya;  z1 = z1+ x(ifx+20)*ya;
                yb = y(i-ifx-20);  z0 = z0+ x(ifx+20)*yb;  z1 = z1+ x(ifx+21)*yb;
                ya = y(i-ifx-21);  z0 = z0+ x(ifx+21)*ya;  z1 = z1+ x(ifx+22)*ya;
                yb = y(i-ifx-22);  z0 = z0+ x(ifx+22)*yb;  z1 = z1+ x(ifx+23)*yb;
                ya = y(i-ifx-23);  z0 = z0+ x(ifx+23)*ya;  z1 = z1+ x(ifx+24)*ya;
                yb = y(i-ifx-24);  z0 = z0+ x(ifx+24)*yb;  z1 = z1+ x(ifx+25)*yb;
                ya = y(i-ifx-25);  z0 = z0+ x(ifx+25)*ya;  z1 = z1+ x(ifx+26)*ya;
                yb = y(i-ifx-26);  z0 = z0+ x(ifx+26)*yb;  z1 = z1+ x(ifx+27)*yb;
                ya = y(i-ifx-27);  z0 = z0+ x(ifx+27)*ya;  z1 = z1+ x(ifx+28)*ya;
                yb = y(i-ifx-28);  z0 = z0+ x(ifx+28)*yb;  z1 = z1+ x(ifx+29)*yb;
                ya = y(i-ifx-29);  z0 = z0+ x(ifx+29)*ya;
                z(i+1) = z1;
                z(i) = z0;
            enddo
        !}
        if (ihigh>=ilow .and. mod((ihigh-ilow),2)==0) then
            ilow = ihigh
            jlow = ifx
            jhigh = ilx

            do i=ilow,ihigh
                sum=0.
                do j=jlow,jhigh
                    sum=sum+x(j)*y(i-j)
                enddo
                z(i)=sum
            enddo
            
        endif
        
        !/* ROLLING OFF:  ily+ifx < i <= ily+ilx */
        ilow = ily+ifx+1;  if (ilow<ifz) ilow = ifz
        ihigh = ily+ilx;  if (ihigh>ilz) ihigh = ilz
        jlow = ilow-ily
        jhigh = ilx

        do i=ilow,ihigh
            sum=0.
            do j=jlow,jhigh
                sum=sum+x(j)*y(i-j)
            enddo
            z(i)=sum
            jlow=jlow+1
        enddo
        
        !/* OFF RIGHT:  ily+ilx < i */
        ilow = ily+ilx+1;  if (ilow<ifz) ilow = ifz
        ihigh = ilz
        do i=ilow,ihigh
            z(i)=0.
        enddo
        
    end subroutine

end module m_hilbert_nofft