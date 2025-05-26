SHELL := /bin/bash

NZ=187
NX=801
DZ=25
DX=25

nz=187
nx=801
h=25

sig=9.6

vmin=1364
vmax=5990
fpeak=4
dt=0.0025
nt=6800

#How good initial model do you want?
siginit=200

.SUFFIXES: .bin .smth .$(h)m .init .mig
.bin.smth:
	#echo $(sig)
	#@echo $^ $@ > par
	#@echo $(NZ) $(NX) 1 >> par
	#@echo $(DZ) $(DX) 1 >> par
	#@echo $(sig) $(sig) $(sig) >> par
	#@echo 2 >> par
	#@echo 0 >> par
	#./smoothgauss2dfast_zhou < par
	cp $^ $@

.smth.$(h)m:
	unisam2 < $^ > tmp \
	  nx1=$(NZ)  dx1=$(DZ)    nx2=$(NX)  dx2=$(DX) \
	  n1=$(nz)    d1=$(h)      n2=$(nx)   d2=$(h)
	#
	@echo tmp $@ > par
	@echo $(nz) $(nx) >> par
	@echo $(vmax) >> par
	./set_vmax  < par
	-rm tmp

.$(h)m.init :
	echo $(siginit)
	@echo $^ $@ > par
	@echo $(nz) $(nx) 1 >> par
	@echo $(h) $(h) 1 >> par
	@echo $(siginit) $(siginit) $(siginit) >> par
	@echo 0 >> par
	@echo 0 >> par
	./smoothgauss2dfast_zhou < par


code : 
	#(cd ../../ ; git pull ; git checkout teaching ; make prepare)
	#(cd ../../FWD ; make cleanall ; make)
	#(cd ../../FWI ; make cleanall ; make)
	gfortran smoothgauss2dfast_zhou.f90 -o smoothgauss2dfast_zhou

clone :
	cp ~/DATA/Overthrust\ \&\ Salt/overthrust2D vp.bin
	farith < vp.bin op=scale scale=0.577  > vs.bin #assuming Poisson solid

smooth : vp.smth vs.smth
	echo "Smooth Done!"

resamp : vp.$(h)m vs.$(h)m
	echo 'Resamp Done!'

init : vp.init vs.init
	makevel v000=1 nz=$(nz) nx=$(nx) > one
	cat vp.init  vs.init  one > init
	echo 'Init Done!'

model : vp.$(h)m vs.$(h)m
	makevel v000=1 nz=$(nz) nx=$(nx) > one
	cat vp.$(h)m vs.$(h)m one  > true
	echo 'Model Done!'
	#ximage < true n1=$(nz) d1=$(h) d2=$(h) title=true wclip=1028 bclip=4700 cmap=hsv6 legend=1 &
	#ximage < init n1=$(nz) d1=$(h) d2=$(h) title=init wclip=1028 bclip=4700 cmap=hsv6 legend=1 &


##############################################################################################

clean :
	-rm -r elastic-marmousi-model
	-rm binary header par
	-rm *.smth *.$(h)m *.init
	-rm setup_*.in
