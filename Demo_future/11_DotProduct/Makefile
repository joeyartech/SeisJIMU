mod_dir=../../mod/
include ../../make.inc
include ../../setup.inc


obj = \
../../Modules/System/*.o \
../../Modules/Etc/*.o \
../../Modules/Signal/*.o \
../../Modules/Modeling/*.o


all : exe

exe : main.o
	mpif90 $(FLAGF90) $(MOD) $(obj) $^ -o ../../exe/dotproduct_$(WaveEq)_$(Solver)_$(Order)
	(cd ../../exe; ln -sf dotproduct_$(WaveEq)_$(Solver)_$(Order) DotProduct )

clean :
	-rm main.o

cleanall :
	#System
	(cd ../../Modules/System; make clean)
	#Etc
	(cd ../../Modules/Etc; make clean)
	#Signal
	(cd ../../Modules/Signal; make clean)
	#Modeling
	(cd ../../Modules/Modeling; make clean)
	#
	make clean
