mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc


obj = \
../../Modules/System/*.o \
../../Modules/Etc/*.o \
../../Modules/Signal/*.o \
../../Modules/Modeling/*.o


all : exe

exe : main.o
	mpif90 $(FLAGF90) $(MOD) $(obj) $^ -o ../../exe/adjointtest_$(WaveEq)_$(Solver)_$(Order)
	(cd ../../exe; ln -sf adjointtest_$(WaveEq)_$(Solver)_$(Order) AdjointTest )

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
