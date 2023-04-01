#!/bin/bash

# makevel nz=100 nx=201 v000=1500 > c1
# makevel nz=1  nx=201 v000=2000 > c2
# cat c1 c2 c1 > tmp && transp < tmp n1=201 > model
# rm c1 c2 tmp

makevel nz=201 nx=201 v000=1000 > model

../../exe/fwd_AC_FDSG_O4  setup.in > out

suwind < results/dsyn_Shot0001.su key=trid min=11 max=11 | suximage perc=99 legend=1 title=p &
suwind < results/dsyn_Shot0001.su key=trid min=12 max=12 | suximage perc=99 legend=1 title=vz &
suwind < results/dsyn_Shot0001.su key=trid min=14 max=14 | suximage perc=99 legend=1 title=vx &
