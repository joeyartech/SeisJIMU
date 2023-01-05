mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc


obj = \
../../Modules/System/*.o \
../../Modules/Etc/*.o \
../../Modules/Signal/m_pseudotime.o

all : mod exe

mod :
	#System
	(cd ../../Modules/System; make)
	#Etc
	(cd ../../Modules/Etc; make)
	#Signal
	(cd ../../Modules/Signal; make m_pseudotime.o)


exe : main.o
	$(F90) $(FLAGF90) $(MOD) $(obj) $^ -o ./test #../../exe/gradienttest_$(WaveEq)_$(Solver)_$(Order)_$(Param)
	#(cd ../../exe; ln -sf gradienttest_$(WaveEq)_$(Solver)_$(Order)_$(Param) GradientTest )

clean :
	-rm main.o

cleanall :
	#System
	(cd ../../Modules/System; make clean)
	#Etc
	(cd ../../Modules/Etc; make clean)
	#Signal
	(cd ../../Modules/Signal; make clean)
	#
	(cd ../../mod; rm *.mod)
	#
	make clean
