#!/bin/bash

make clean

#fortran side
gfortran -ffree-line-length-none   ../../External/sinc.f90  ../../External/stoep.f90 ../m_resamp.f90  main.f90
./a.out
xwigb < out n1=200 d1=0.001 f1=0.02 &

#su side
#delrt=100  -> start timing at 0.1s
#tmin=0.02  -> start timing at 0.02s
#dt_in=0.002   -> input trace dt=0.002s
make
suaddhead < in ns=100 | sushw key=delrt a=100 | ./main dt_in=0.002 nt=200 dt=0.001  tmin=0.02 | suxwigb
