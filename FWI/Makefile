mod_dir=../mod/
include ../make.inc
include ../setup.inc


obj = \
../Modules/System/*.o \
../Modules/Etc/*.o \
../Modules/Signal/*.o \
../Modules/Modeling/*.o \
../Modules/Kernel/*.o \
../Modules/Optimization/*.o


all : mod exe

mod :
	#System
	(cd ../Modules/System; make)
	#Etc
	(cd ../Modules/Etc; make)
	#Signal
	(cd ../Modules/Signal; make)
	#Modeling
	(cd ../Modules/Modeling; make)
	#Kernel
	(cd ../Modules/Kernel; make)
	#Optimization
	(cd ../Modules/Optimization; make)

exe : main.o
	mpif90 $(FLAGF90) $(MOD) $(obj) $^ -o ../exe/fwi_$(WaveEq)_$(Solver)_$(Order)_$(LineS)_$(Optim)
	(cd ../exe; ln -sf fwi_$(WaveEq)_$(Solver)_$(Order)_$(LineS)_$(Optim) FWI )

clean :
	-rm main.o

cleanall :
	#System
	(cd ../Modules/System; make clean)
	#Etc
	(cd ../Modules/Etc; make clean)
	#Signal
	(cd ../Modules/Signal; make clean)
	#Modeling
	(cd ../Modules/Modeling; make clean)
	#Kernel
	(cd ../Modules/Kernel; make clean)
	#Optimization
	(cd ../Modules/Optimization; make clean)
	#FWI
	make clean