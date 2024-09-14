!Fornberg, 1988, Generation of Finite-Difference Formulas on Arbitrary Spaced Grids.

!grid points: 0,1,2,3,4
!!1st-order derivative
!negative grid points (-1,-2,-3) should negate the coefs
fd_rg_O2 = 0, 0.5
fd_rg_O4 = 0, 2./3., -1./12.
fd_rg_O6 = 0, 3./4., -3./20., 1./60.
fd_rg_08 = 0, 4./5., -1./5., 4./105., -1/280

!!2nd-order derivative
!negative grid points (-1,-2,-3) should share the same coefs
fd_rg_O2 = -2., 1.
fd_rg_O4 = -5./2., 4./3., -1./12.
fd_rg_O6 = -49./18, 3./2., -3./20., 1./90.
fd_rg_08 = -205/72, 8./5., -1./5., 8./315., -1/560


!half grid points: 0.5,1.5,2.5,3.5
!!1st-order derivative
!negative half grid points (-0.5,-1.5,-2.5,-3.5) should negate the coefs
fd_rg_O2 = 1.
fd_rg_O4 = 9./8., -1./24.
fd_rg_O6 = 75./64., -25./384., 3./640.
fd_rg_08 = 1225./1024., -245./3072., 49./5120., -5./7168.

!!2nd-order derivative
!negative half grid points (-0.5,-1.5,-2.5,-3.5) share the same coefs
fd_rg_O2 = -1./2., 1./2.
fd_rg_O4 = -17./24., 13./16., -5./48.
fd_rg_O6 = -1891./2304., 1299./1280., -499./2304., 259./11520.
