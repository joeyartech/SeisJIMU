#!/bin/bash

# make cleanall; make
# 
makevel nz=100 nx=201 v000=1000 > c1
makevel nz=101 nx=201 v000=1800 > c2
cat c1 c2 > tmp && transp < tmp n1=201 > simple
rm c1 c2 tmp
# 

# makevel nx=201 nz=201 v000=2000 > model

cp setup_default setup

# echo 'p p'
# echo 'SCOMP p' >> setup
# echo 'RCOMP p' >> setup
# #echo "FS_METHOD    zero_stress" >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out
#
# ########################################
#
# echo 'vz p'
# echo 'SCOMP vz' >> setup
# echo 'RCOMP p ' >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out
#
# echo 'vx p' #not successful
# echo 'SCOMP vx' >> setup
# echo 'RCOMP p ' >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out
#
# echo 'p vz' #not successful
# echo 'SCOMP p ' >> setup
# echo 'RCOMP vz' >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out
#
# echo 'p vx' #not successful
# echo 'SCOMP p ' >> setup
# echo 'RCOMP vx' >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out
#
#######################################

echo 'vz vz'
echo 'SCOMP vz' >> setup
echo 'RCOMP vz' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

# echo 'vx vx'
# echo 'SCOMP vx' >> setup
# echo 'RCOMP vx' >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out

echo 'vz vx' #not successful
echo 'SCOMP vz' >> setup
echo 'RCOMP vx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'vx vz' #not successful
echo 'SCOMP vx' >> setup
echo 'RCOMP vz' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

#######################################
#
# echo 'szx szx'
# echo 'SCOMP szx' >> setup
# echo 'RCOMP szx'  >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out
#
# echo 'szx vz'
# echo 'SCOMP szx' >> setup
# echo 'RCOMP vz'  >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out
#
# echo 'szx vx'
# echo 'SCOMP szx' >> setup
# echo 'RCOMP vx'  >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out
#
# echo 'vz szx'
# echo 'SCOMP vz'  >> setup
# echo 'RCOMP szx' >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out
#
# echo 'vx szx'
# echo 'SCOMP vx'  >> setup
# echo 'RCOMP szx' >> setup
# ../../exe/AdjointTest  setup > out
#
# tail -10 out

# suximage <  u.su  legend=1  title='u' &
# suximage <  v.su  legend=1  title='v' &
# suximage < Lu.su  legend=1  title='Lu' &
# suximage < Ladj_v.su  legend=1  title='Ladj_v' &
# sumax < Lu.su

