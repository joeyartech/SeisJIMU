mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc
#make sure the correct module is used
LineS=enumerate


obj = \
../../Modules/System/*.o \
../../Modules/Etc/*.o \
../../Modules/Signal/*.o \
../../Modules/Modeling/m_model.o \
../../Modules/Modeling/m_shotlist_$(ShotDec).o \
../../Modules/Modeling/m_shot.o \
../../Modules/Modeling/m_computebox.o \
../../Modules/Modeling/m_correlate.o \
../../Modules/Modeling/m_field.o \
../../Modules/Modeling/m_cpml.o \
../../Modules/Modeling/m_propagator_$(WaveEq)_$(Solver)_$(Order).o \
../../Modules/Modeling/m_Modeling.o \
../../Modules/Kernel/m_parametrizer_$(Param).o \
../../Modules/Kernel/m_querypoint.o \
../../Modules/Kernel/m_separator.o \
../../Modules/Kernel/m_weighter.o \
../../Modules/Kernel/m_preconditioner.o \
../../Modules/Kernel/m_fobjective.o \
../../Modules/Kernel/m_Kernel.o \
../../Modules/Optimization/m_linesearcher_$(LineS).o


all : mod exe

mod :
	#System
	(cd ../../Modules/System; make)
	#Etc
	(cd ../../Modules/Etc; make)
	#Signal
	(cd ../../Modules/Signal; make)
	#Modeling
	(cd ../../Modules/Modeling; make)
	#Kernel
	(cd ../../Modules/Kernel; make App=FWI)
	#Linesearch
	(cd ../../Modules/Optimization; make m_linesearcher_$(LineS).o)


exe : ../../$(App)/sub_modeling_gradient.o main.o sub_optimizer_init_loop.o
	$(F90) $(FLAGF90) $(MOD) $(obj) $^ -o ../../exe/gradienttest_$(WaveEq)_$(Solver)_$(Order)_$(Param)
	(cd ../../exe; ln -sf gradienttest_$(WaveEq)_$(Solver)_$(Order)_$(Param) GradientTest )

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
	#Kernel
	(cd ../../Modules/Kernel; make clean)
	#Optimization
	(cd ../../Modules/Optimization; rm *.o)
	#
	(cd ../../mod; rm *.mod)
	#
	make clean
