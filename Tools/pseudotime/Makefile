mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc


obj = \
../../Modules/System/*.o \
../../Modules/Signal/m_pseudotime.o


all : mod exe

mod :
	(cd ../../Modules/System; make)
	(cd ../../Modules/Signal; make m_pseudotime.o)

exe : main.o
	mpif90 $(FLAGF90) $(MOD) $(obj) $^ -o ../../exe/pseudotime

clean :
	-rm main.o

cleanall :
	(cd ../../Modules/System; make clean)
	(cd ../../Modules/Signal; rm m_pseudotime.o)
	rm ../../mod/m_pseudotime.mod
	make clean