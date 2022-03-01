#!/bin/bash

# # FWD ##
# makevel nz=49 nx=101 v000=1500 > c1
# makevel nz=3  nx=101 v000=2000 > c2
# cat c1 c2 c1 > tmp && transp < tmp n1=101 > model
# rm c1 c2 tmp
# 
# ../../exe/FWD setup.in
# 
# rm -r results_fwd
# mv results  results_fwd

# FWI ##
makevel nz=101 nx=101 v000=1800 > model
makevel nz=1 nx=101 v000=400 > topo

rm -r results
../../exe/gradienttest_wpi  setup.in > out

echo '            alpha    pert%f    curr%f    (pert%f-curr%f)/alpha    curr%g_dot_d    if_1st_cond'
grep '1st cond' out
