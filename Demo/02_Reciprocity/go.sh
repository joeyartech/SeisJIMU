#!/bin/bash

# makevel nz=100 nx=201 v000=1500 > c1
# makevel nz=1  nx=201 v000=2000 > c2
# cat c1 c2 c1 > tmp && transp < tmp n1=201 > model
# rm c1 c2 tmp

#makevel nz=201 nx=201 v000=1000 > model

cp setup0.in  setup.in
echo "NR                      1" >> setup.in

echo "FS                      '26 505 500'" >> setup.in
echo "FR                      '46 10  0'" >> setup.in
echo "SCOMP                   vz " >> setup.in
echo "RCOMP                   vx " >> setup.in
../../exe/fwd_AC_FDSG_O4  setup.in
mv dsyn_0001.su s-r.su

echo "FR                      '26 505 500'" >> setup.in
echo "FS                      '46 10  0'" >> setup.in
echo "RCOMP                   vz " >> setup.in
echo "SCOMP                   vx " >> setup.in
../../exe/fwd_AC_FDSG_O4  setup.in
mv dsyn_0001.su r-s.su

cat s-r.su r-s.su | suxgraph
suop < s-r.su  op=neg > tmp.su ; cat tmp.su  r-s.su  | suxgraph
