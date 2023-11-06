#!/bin/bash

makevel nz=100 nx=201 v000=1000 > c1
makevel nz=101 nx=201 v000=1800 > c2
cat c1 c2 > tmp && transp < tmp n1=201 > model
#makevel nz=1  nx=201 v000=1800 > c2
# cat c1 c2 c1 > tmp && transp < tmp n1=201 > model
rm c1 c2 tmp

#makevel nz=201 nx=201 v000=1000 > model

../../exe/fwd_PSV_FDSG_O4_  setup_PSV &> out_PSV
../../exe/fwd_DAS_FDSG_O4_  setup_DAS &> out_DAS


suwind < results_PSV/Ru_Shot0001.su key=trid min=12 max=12 | suximage perc=99 legend=1 title=PSV_vz &
suwind < results_DAS/Ru_Shot0001.su key=trid min=32 max=32 | suximage perc=99 legend=1 title=DAS_pz &

suwind < results_PSV/Ru_Shot0001.su key=trid min=14 max=14 | suximage perc=99 legend=1 title=PSV_vx &
suwind < results_DAS/Ru_Shot0001.su key=trid min=34 max=34 | suximage perc=99 legend=1 title=DAS_px &
