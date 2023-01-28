#!/bin/bash

n=101

# FWD ##
makevel nz=49 nx=$n v000=1500 > vp1
makevel nz=3  nx=$n v000=1800 > vp2
makevel nz=49 nx=$n v000=1000 > rho1
makevel nz=3  nx=$n v000=2000 > rho2
cat vp1  vp2  vp1  > tmp1 && transp < tmp1 n1=$n > tmpvp
cat rho1 rho2 rho1 > tmp1 && transp < tmp1 n1=$n > tmprho
cat tmpvp tmprho > model
rm vp1 vp2 rho1 rho2 tmp*

../../exe/FWD setup.in > out_FWD

rm -r results_fwd
mv results  results_fwd


## FWI ##
makevel nz=$n nx=$n v000=1750 > vp
makevel nz=$n nx=$n v000=1500 > rho
cat vp rho > model; rm vp rho
makevel nz=1 nx=$n v000=100 > topo
#makevel nz=1 nx=$n v000=400 > topo

rm -r results
../../exe/GradientTest  setup.in > out_FWI

echo '            alpha    pert%f    curr%f    (pert%f-curr%f)/alpha    curr%g_dot_d    if_1st_cond'
grep '1st cond' out_FWI
