! A Fortran version of stoep.c in Seismic Unix
! dir: SeisUnix/src/cwp/lib
!
! !! Copyright (c) Colorado School of Mines, 2011.!!
! !! All rights reserved.                         !!
! 
! !!********************** self documentation *********************!!
! !!****************************************************************************
! STOEP - Functions to solve a symmetric Toeplitz linear system of equations
!    Rf=g for f
! 
! stoepd        solve a symmetric Toeplitz system - doubles
! stoepf        solve a symmetric Toeplitz system - floats
! 
! ******************************************************************************
! Function Prototypes:
! void stoepd (int n, double r[], double g[], double f[], double a[]);
! void stoepf (int n, float r[], float g[], float f[], float a[]);
! 
! ******************************************************************************
! Input:
! n     dimension of system
! r     array[n] of top row of Toeplitz matrix
! g     array[n] of right-hand-side column vector
! 
! Output:
! f     array[n] of solution (left-hand-side) column vector
! a     array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
! 
! ******************************************************************************
! Notes:
! These routines do NOT solve the case when the main diagonal is zero, it
! just silently returns.
! 
! The left column of the Toeplitz matrix is assumed to be equal to the top
! row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.
! 
! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! *****************************************************************************/
! /**************** end self doc ********************************/
! 
! #include "cwp.h"
! 
! void stoepf (int n, float r[], float g[], float f[], float a[])
! /*****************************************************************************
! Solve a symmetric Toeplitz linear system of equations Rf=g for f
! (float version) 
! ******************************************************************************
! Input:
! n		dimension of system
! r		array[n] of top row of Toeplitz matrix
! g		array[n] of right-hand-side column vector
! 
! Output:
! f		array[n] of solution (left-hand-side) column vector
! a		array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
! ******************************************************************************
! Notes:
! This routine does NOT solve the case when the main diagonal is zero, it
! just silently returns.
! 
! The left column of the Toeplitz matrix is assumed to be equal to the top
! row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.
! ******************************************************************************
! Author:  Dave Hale, Colorado School of Mines, 06/02/89
! *****************************************************************************/
    subroutine stoep(n,r,g,f,a)
        real,dimension(0:n-1) :: r,g,f,a

        if (r(0) == 0.) return

        a(0) = 1.
        v = r(0)
        f(0) = g(0)/r(0)

        do j=1,n-1
            !solve Ra=v as in Claerbout, FGDP, p. 57
            a(j) = 0.
            f(j) = 0.
            
            e=0.
            do i=0,j-1
                e = e + a(i)*r(j-i)
            enddo
            
            c = e/v
            v = v - c*e
            
            do i=0,j/2
                bot = a(j-i) - c*a(i);
                a(i) = a(i) - c*a(j-i);
                a(j-i) = bot;
            enddo

            !use a and v above to get f[i], i = 0,1,2,...,j
            w=0.
            do i=0, j-1
                w = w + f(i)*r(j-i)
            enddo
            
            c = (w-g(j))/v
            
            do i=0,j
                f(i) = f(i) - c*a(j-i)
            enddo
            
        enddo
    
    end subroutine
