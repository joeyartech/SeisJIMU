#!/bin/bash

if [ "$1" == "compile" ]; then
    make cleanall ; (cd ../../ ; make fwd ) ; make
    exit
fi
# 
makevel nz=100 nx=201 v000=1000 > c1
makevel nz=101 nx=201 v000=1800 > c2
cat c1 c2 > tmp && transp < tmp n1=201 > simple
rm c1 c2 tmp
# 

# makevel nx=201 nz=201 v000=2000 > model

cp setup_default setup

# echo 'pz ez'
# echo 'SCOMP pz' >> setup
# echo 'RCOMP ez' >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out


echo 'pz ex'
echo 'SCOMP pz' >> setup
echo 'RCOMP ex' >> setup
../../exe/AdjointTest  setup > out

tail -10 out


echo 'pz es'
echo 'SCOMP pz' >> setup
echo 'RCOMP es' >> setup
../../exe/AdjointTest  setup > out

tail -10 out


echo 'px ez'
echo 'SCOMP px' >> setup
echo 'RCOMP ez' >> setup
../../exe/AdjointTest  setup > out

tail -10 out


echo 'ez pz'
echo 'SCOMP ez' >> setup
echo 'RCOMP pz' >> setup
../../exe/AdjointTest  setup > out

tail -10 out


echo 'ex pz'
echo 'SCOMP ex' >> setup
echo 'RCOMP pz' >> setup
../../exe/AdjointTest  setup > out

tail -10 out


echo 'es pz'
echo 'SCOMP es' >> setup
echo 'RCOMP pz' >> setup
../../exe/AdjointTest  setup > out

tail -10 out



# suximage <  u.su  legend=1  title='u' &
# suximage <  v.su  legend=1  title='v' &
# suximage < Lu.su  legend=1  title='Lu' &
# suximage < Ladj_v.su  legend=1  title='Ladj_v' &
# sumax < Lu.su

