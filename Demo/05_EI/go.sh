#!/bin/bash

# makevel nz=50 nx=101 v000=1500 > c1
# makevel nz=1  nx=101 v000=2000 > c2
# cat c1 c2 c1 > tmp && transp < tmp n1=101 > model
# rm c1 c2 tmp

makevel nz=101 nx=101 v000=1600 > model

##data=wave
# ../../exe/PFEI  setup.in > out

#data=envelope
suenv < ../01_ForwardModeling/results/dsyn_Shot0001.su | sushw key=trid a=11 > denv_Shot0001.su
cp setup.in setup2
echo 'FILE_DATA_PREFIX  denv_Shot' >> setup2
echo 'DIR_OUT           results_env' >> setup2
../../exe/PFEI  setup2

