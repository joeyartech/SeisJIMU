mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc


obj = \
../../Modules/System/*.o \
../../Modules/Etc/*.o \
../../Modules/Signal/*.o


all : mod exe

mod :
	#System
	(cd ../../Modules/System; make)
	#Etc
	(cd ../../Modules/Etc; make)
	#Signal
	(cd ../../Modules/Signal; make)

exe : main.o
	mpif90 $(FLAGF90) $(MOD) $(obj) $^ -o ../../exe/time2depth

clean :
	-rm main.o

cleanall :
	#System
	(cd ../Modules/System; make clean)
	#Etc
	(cd ../Modules/Etc; make clean)
	#Signal
	(cd ../Modules/Signal; make clean)