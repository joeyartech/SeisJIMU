mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc


obj = \
../../Modules/System/*.o \
../../Modules/Etc/*.o \
../../Modules/Signal/m_Hilbert.o \
../../Modules/Signal/m_Hilbert_nofft.o


all : mod exe

mod :
	(cd ../../Modules/System; make)
	(cd ../../Modules/Etc; make)
	(cd ../../Modules/Signal; make m_Hilbert.o m_Hilbert_nofft.o)

#exe : main.o
exe : test2.o
	mpif90 $(FLAGF90) $(MOD) $(obj) $^ -o ../../exe/envelope

clean :
	-rm test2.o main.o

cleanall :
	(cd ../../Modules/System; make clean)
	(cd ../../Modules/Etc; make clean)
	(cd ../../Modules/Signal; rm m_Hilbert.o)
	rm ../../mod/m_Hilbert.mod ../../mod/m_Hilbert_nofft.mod
	make clean
