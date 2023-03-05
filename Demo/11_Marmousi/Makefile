SHELL := /bin/bash

NZ=2801
NX=13601
DZ=1.25
DX=1.25
TOPO=450

nz=281
nx=1361
h=12.5
sig:=$(shell echo "scale=4; 0.48/( 0.5/$(h) )" | bc -l )
#sig:=$(shell echo "scale=4; 0.19/( 0.5/$(h) )" | bc -l )
fpeak=6.4
dt=0.0015
nt=4600

#How good initial model do you want?
siginit=200

.SUFFIXES: .bin .smth .$(h)m .init .mig
.bin.smth:
	echo $(sig)
	@echo $^ $@ > par
	@echo $(NZ) $(NX) 1 >> par
	@echo $(DZ) $(DX) 1 >> par
	@echo $(sig) $(sig) $(sig) >> par
	@echo 2 >> par
	@echo $(TOPO) >> par
	./smoothgauss2dfast_zhou < par

.smth.$(h)m:
	unisam2 < $^ > $@ \
	  nx1=$(NZ)  dx1=$(DZ)    nx2=$(NX)  dx2=$(DX) \
	  n1=$(nz)    d1=$(h)      n2=$(nx)   d2=$(h)

.$(h)m.init :
	echo $(siginit)
	@echo $^ $@ > par
	@echo $(nz) $(nx) 1 >> par
	@echo $(h) $(h) 1 >> par
	@echo $(siginit) $(siginit) $(siginit) >> par
	@echo 0 >> par
	@echo $(TOPO) >> par
	./smoothgauss2dfast_zhou < par


#How many MPI processors do you order?
np=1
#How many OpenMP threads do you order?
nthreads=8


###########
# Prepare #
###########

code : 
	(cd ../../ ; git pull ; git checkout teaching ; make prepare)
	(cd ../../FWD ; make cleanall ; make)
	(cd ../../FWI ; make cleanall ; make)
	gfortran smoothgauss2dfast_zhou.f90 -o smoothgauss2dfast_zhou

clone :
	#wget https://s3.amazonaws.com/open.source.geoscience/open_data/elastic-marmousi/elastic-marmousi-model.tar.gz #https://wiki.seg.org/wiki/AGL_Elastic_Marmousi
	tar -xvzf elastic-marmousi-model.tar.gz
	(cd ./elastic-marmousi-model/model/ ; for fi in `ls *.tar.gz`; do tar -xvzf $$fi; done)
	segyread tape=./elastic-marmousi-model/model/MODEL_P-WAVE_VELOCITY_1.25m.segy | sustrip > vp.bin
	segyread tape=./elastic-marmousi-model/model/MODEL_DENSITY_1.25m.segy | sustrip | farith op=scale scale=1000 > rho.bin  #convert to [kg/m3]

smooth : vp.smth rho.smth
	echo "Smooth Done!"

resamp : vp.$(h)m rho.$(h)m
	echo 'Resamp Done!'

init : vp.init
	echo 'Init Done!'

topo :
	makevel nx=$(nx) nz=1 v000=$(TOPO) > $@

model : vp.$(h)m rho.$(h)m topo vp.init rho.init
	cat vp.$(h)m rho.$(h)m > true
	cat vp.init            > init
	echo 'Model Done!'
	ximage < true n1=$(nz) d1=$(h) d2=$(h) title=true wclip=1028 bclip=4700 cmap=hsv6 legend=1 &
	ximage < init n1=$(nz) d1=$(h) d2=$(h) title=init wclip=1028 bclip=4700 cmap=hsv6 legend=1 &


#################
# Observed data #
#################

fwd :
	cp setup.in setup_fwd.in
	echo "DIR_OUT         './results_fwd' " >> setup_fwd.in
	echo "NS              1 " >> setup_fwd.in
	-rm -r ./results_fwd
	(export OMP_NUM_THREADS=$(nthreads); ulimit -s unlimited; mpirun -np $(np) ../../exe/FWD setup_fwd.in > out_fwd)
	suximage < results_fwd/dsyn_Shot0001.su title='Observed data' perc=99 legend=1 &


################
# FWI Gradient #
################

ishot=1
#ishot=1:10:80
#ishot=1:20:144

grad :
	cp setup.in setup_grad.in
	echo "DIR_OUT            './results_grad' " >> setup_grad.in
	echo "FILE_MODEL         init" >> setup_grad.in
	echo "MODEL_ATTRIBUTES   vp" >> setup_grad.in
	echo "ISHOT              '$(ishot)' " >> setup_grad.in
	-rm -r ./results_grad
	(export OMP_NUM_THREADS=$(nthreads); ulimit -s unlimited; mpirun -np $(np) ../../exe/FWI setup_grad.in > out_grad)
	suximage < results_grad/dsyn_Shot0001.su title='Calculated data' perc=99 legend=1 &
	suximage < results_grad/dadj_Shot0001.su title='Data residual' perc=99 legend=1 &
	ximage < results_grad/qp0%g n1=$(nz) d1=$(dz) d2=$(dx) title='initial gradient' legend=1 &
	ximage < results_grad/qp0%pg n1=$(nz) d1=$(dz) d2=$(dx) title='initial precond gradient' legend=1 &


#################
# FWI Iteration #
#################

ishot=1
#ishot=1:10:80
#ishot=1:20:144

fwi :
	cp setup.in setup_fwi.in
	echo "DIR_OUT            './results_fwi' " >> setup_fwi.in
	echo "FILE_MODEL         init" >> setup_fwi.in
	echo "MODEL_ATTRIBUTES   vp" >> setup_fwi.in
	echo "ISHOT              '$(ishot)' " >> setup_fwi.in
	echo "SMOOTHING          Laplacian  " >> setup_fwi.in
	echo "JOB                optimization " >> setup_fwi.in
	-rm -r ./results_fwi
	(export OMP_NUM_THREADS=$(nthreads); ulimit -s unlimited; mpirun -np $(np) ../../exe/FWI setup_fwi.in > out_fwi)
	ximage < init n1=$(nz) d1=$(dz) d2=$(dx) title='initial model' wclip=1028 bclip=4700 legend=1 &
	ximage < results_fwi/model_final n1=$(nz) d1=$(dz) d2=$(dx) title='FWI model' wclip=1028 bclip=4700 legend=1 &


##############################################################################################

clean :
	-rm -r elastic-marmousi-model
	-rm binary header par
	-rm *.smth *.$(h)m *.init
	-rm setup_*.in
