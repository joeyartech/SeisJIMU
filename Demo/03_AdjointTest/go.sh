#!/bin/bash

if [ "$1" == "compile" ]; then
    make cleanall ; (cd ../../ ; make fwd ) ; make
fi
# 
# makevel nz=50 nx=101 v000=1500 > c1
# makevel nz=1  nx=101 v000=2000 > c2
# cat c1 c2 c1 > tmp && transp < tmp n1=101 > model
# rm c1 c2 tmp
# 

#makevel nx=201 nz=201 v000=2000 > simple

../../exe/AdjointTest  setup_simple.in

# suximage <  u.su  legend=1  title='u' &
# suximage <  v.su  legend=1  title='v' &
# suximage < Lu.su  legend=1  title='Lu' &
# suximage < Ladj_v.su  legend=1  title='Ladj_v' &
# sumax < Lu.su

