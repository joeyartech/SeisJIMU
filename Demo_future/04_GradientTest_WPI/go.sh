#!/bin/bash

# # FWD ##
# makevel nz=49 nx=101 v000=1500 > c1
# makevel nz=3  nx=101 v000=2000 > c2
# cat c1 c2 c1 > tmp && transp < tmp n1=101 > model
# rm c1 c2 tmp
# 
# ../../exe/FWD setup.in > out
# 
# rm -r results_fwd
# mv results  results_fwd

# FWI ##
makevel nz=101 nx=201 v000=1800 > model
makevel nz=1 nx=201 v000=400 > topo

rm -r results
../../exe/gradienttest_wpi  setup.in > out

echo '            alpha    pert%f    curr%f    (pert%f-curr%f)/alpha    curr%g_dot_d    if_1st_cond'
grep '1st cond' out

# ximage < ../04_GradientTest/results/qp0%g n1=101 clip=0.0005 legend=1 title=FWI &
# farith < results/du_star_a in2=results/u_star_da op=add > results/RE
# ximage < results/RE n1=101 clip=6e9 legend=1 title=RE &
# ximage < results/u_star_Adj*  n1=101 perc=99 legend=1 title=DR
