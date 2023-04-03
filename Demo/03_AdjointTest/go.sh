#!/bin/bash

if [ "$1" == "compile" ]; then
    make cleanall ; (cd ../../ ; make fwd ) ; make
fi
# 
makevel nz=101 nx=101 v000=1500 > rho0

makevel nz=50 nx=101 v000=1500 > r1
makevel nz=1  nx=101 v000=3000 > r2
cat r1 r2 r1 > tmp && transp < tmp n1=101 > rho

cat rho0 rho > simple

rm r1 r2 tmp rho0 rho
# 

# makevel nx=201 nz=201 v000=2000 > model

../../exe/AdjointTest  setup_simple.in

# suximage <  u.su  legend=1  title='u' &
# suximage <  v.su  legend=1  title='v' &
# suximage < Lu.su  legend=1  title='Lu' &
# suximage < Ladj_v.su  legend=1  title='Ladj_v' &
# sumax < Lu.su

