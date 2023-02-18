#!/bin/bash

n=101

# # FWD ##
# makevel nz=51 nx=$n v000=1800 > vp1
# makevel nz=50 nx=$n v000=1500 > vp2
# cat vp1  vp2  > tmp1 && transp < tmp1 n1=$n > model
# rm vp1 vp2 rho1 rho2 tmp*
#
# ../../exe/FWD setup.in > out_fwd
#
# rm -r results_fwd
# mv results  results_fwd


## PFEI gradient ##
makevel nz=51 nx=$n v000=1800 > vp1
makevel nz=50 nx=$n v000=1500 > vp2
cat vp1  vp2 > tmp1 && transp < tmp1 n1=$n > model

# makevel nz=$n nx=$n v000=1750 > model
#
# makevel nz=49 nx=$n v000=1750 > vp1
# makevel nz=52 nx=$n v000=1500 > vp2
# cat vp1  vp2 > tmp1 && transp < tmp1 n1=$n > model

makevel nz=1 nx=$n v000=100 > topo
#makevel nz=1 nx=$n v000=400 > topo

rm vp1 vp2 rho1 rho2 tmp*

rm -r results_grad
../../exe/PFEI setup.in #> out_grad
mv results results_grad

suximage < results_grad/RE0_Shot0001.su clip=1e-5 &
suximage < results_grad/RdE_Shot0001.su clip=1e-5 &
susum  results_grad/RE0_Shot0001.su  results_grad/RdE_Shot0001.su | suximage clip=1e-5 &

# ## GradientTest ##
# rm -r results
# ../../exe/GradientTest  setup.in > out
#
# echo '            alpha    pert%f    curr%f    (pert%f-curr%f)/alpha    curr%g_dot_d    if_1st_cond'
# grep '1st cond' out


