#!/bin/bash

n=201
v0=1800

# FWD ##
if [ $n -eq 101 ]; then
   makevel nz=49 nx=$n v000=1500 > c1
   makevel nz=3  nx=$n v000=2000 > c2
else
   makevel nz=98 nx=$n v000=1500 > c1
   makevel nz=5  nx=$n v000=2000 > c2
fi
cat c1 c2 c1 > tmp && transp < tmp n1=$n > model
rm c1 c2 tmp

../../exe/FWD setup.in > out

rm -r results_fwd
mv results  results_fwd


## WPI ##
makevel nz=$n nx=$n v000=$v0 > model
makevel nz=1 nx=$n v000=100 > topo
#makevel nz=1 nx=$n v000=400 > topo

rm -r results
../../exe/gradienttest_wpi  setup.in > out

echo '            alpha    pert%f    curr%f    (pert%f-curr%f)/alpha    curr%g_dot_d    ratio'
grep '1st cond' out

# farith < results/du_star_a in2=results/u_star_da op=add > results/RE
# ximage < results/RE n1=$n perc=99  legend=1 title=RE &
# ximage < results/-u_star_Adj*  n1=$n perc=99 legend=1 title=DR
