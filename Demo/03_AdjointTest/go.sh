#!/bin/bash

if [ "$1" == "compile" ]; then
    make cleanall ; (cd ../../ ; make fwd ) ; make
fi
# 
makevel nz=50 nx=101 v000=1500 > c1
makevel nz=1  nx=101 v000=2000 > c2
cat c1 c2 c1 > tmp && transp < tmp n1=101 > simple
rm c1 c2 tmp
# 

# makevel nx=201 nz=201 v000=2000 > model


#echo "IS_Q_ATTENUATION    F"


../../exe/AdjointTest  setup_simple.in

suxwigb < results/reS.su      legend=1  title='reS'     &
suxwigb < results/imS.su      legend=1  title='imS'     &
suxwigb < results/reLadjR.su  legend=1  title='reLadjR' &
suxwigb < results/imLadjR.su  legend=1  title='imLadjR' &

