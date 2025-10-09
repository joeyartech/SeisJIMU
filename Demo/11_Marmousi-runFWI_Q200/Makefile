SHELL := /bin/bash

NZ=2396
NX=13601
DZ=1.25
DX=1.25

nz=143
nx=810
h=21
sig:=$(shell echo "scale=4; 0.48/( 0.5/$(h) )" | bc -l )
#sig:=$(shell echo "scale=4; 0.19/( 0.5/$(h) )" | bc -l )
fpeak=4
dt=0.00175

#How good initial model do you want?
siginit=200

.SUFFIXES: .bin .smth .$(h)m .init .mig
.bin.smth:
	echo $(sig)
	@echo $^ $@ > par
	@echo $(NZ) $(NX)   >> par
	@echo $(DZ) $(DX)   >> par
	@echo $(sig) $(sig) >> par
	@echo 2 >> par
	@echo 0 >> par
	./smoothgauss2dfast_zhou < par

.smth.$(h)m:
	unisam2 < $^ > $@ \
	  nx1=$(NZ)  dx1=$(DZ)    nx2=$(NX)  dx2=$(DX) \
	  n1=$(nz)    d1=$(h)      n2=$(nx)   d2=$(h)

.$(h)m.init :
	echo $(siginit)
	@echo $^ $@ > par
	@echo $(nz) $(nx) >> par
	@echo $(h)  $(h)  >> par
	@echo $(siginit) $(siginit) >> par
	@echo 0 >> par
	@echo 0 >> par
	./smoothgauss2dfast_zhou < par

true : vpt.$(h)m rhot.$(h)m qpt.$(h)m

init : vpt.init rhot.init qpt.init
