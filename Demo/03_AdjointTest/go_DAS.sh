#!/bin/bash

# make cleanall; make
#
# makevel nz=100 nx=201 v000=1000 > c1
# makevel nz=101 nx=201 v000=1800 > c2
# cat c1 c2 > tmp && transp < tmp n1=201 > simple
# rm c1 c2 tmp
#
#
# makevel nx=201 nz=201 v000=2000 > model

cp setup_default setup

echo 'pz pz'
echo 'SCOMP pz' >> setup
echo 'RCOMP pz' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'vx vx'
echo 'SCOMP vx' >> setup
echo 'RCOMP vx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'pz vx' #not very successful
echo 'SCOMP pz' >> setup
echo 'RCOMP vx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

#####################################

echo 'p pz' #not successful successful when close to free surface
echo 'SCOMP p ' >> setup
echo 'RCOMP pz' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'p vx' #not successful successful when close to free surface
echo 'SCOMP p ' >> setup
echo 'RCOMP vx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

###################################

echo 'szz szz'
echo 'SCOMP szz' >> setup
echo 'RCOMP szz'  >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'sxx sxx'
echo 'SCOMP sxx' >> setup
echo 'RCOMP sxx'  >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'szz sxx' #fail when on the free surface
echo 'SCOMP szz' >> setup
echo 'RCOMP sxx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'p szz' #not successful successful when close to free surface
echo 'SCOMP p' >> setup
echo 'RCOMP sxx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'p sxx' #not successful successful when close to free surface
echo 'SCOMP p' >> setup
echo 'RCOMP sxx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'pz szz'
echo 'SCOMP pz'  >> setup
echo 'RCOMP szz' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'vx szz'
echo 'SCOMP vx'  >> setup
echo 'RCOMP szz' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'pz szz' #not so successful when close to free surface
echo 'SCOMP vx'  >> setup
echo 'RCOMP szz' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'vx sxx' #not so successful when close to free surface
echo 'SCOMP vx'  >> setup
echo 'RCOMP sxx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

######################################

echo 'szx szx'
echo 'SCOMP szx' >> setup
echo 'RCOMP szx'  >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'pz szx'
echo 'SCOMP pz'  >> setup
echo 'RCOMP szx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'vx szx'
echo 'SCOMP vx'  >> setup
echo 'RCOMP szx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'szz szx' #fail when close to free surface
echo 'SCOMP szz'  >> setup
echo 'RCOMP szx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'sxx szx' #not successful successful when close to free surface
echo 'SCOMP sxx'  >> setup
echo 'RCOMP szx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

echo 'p szx' #not successful successful when close to free surfacegit
echo 'SCOMP p' >> setup
echo 'RCOMP szx' >> setup
../../exe/AdjointTest  setup > out

tail -10 out

######################################
#
# suximage <  u.su  legend=1  title='u' &
# suximage <  v.su  legend=1  title='v' &
# suximage < Lu.su  legend=1  title='Lu' &
# suximage < Ladj_v.su  legend=1  title='Ladj_v' &
# sumax < Lu.su

