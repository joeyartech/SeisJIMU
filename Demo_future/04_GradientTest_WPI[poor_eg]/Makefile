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
../../Modules/Modeling/m_field.o \
../../Modules/Modeling/m_cpml_$(Solver).o \
../../Modules/Modeling/m_propagator_$(WaveEq)_$(Solver)_$(Order).o \
../../Modules/Modeling/m_Modeling.f90 \
../../Modules/Kernel/m_parametrizer_$(Param).o \
../../Modules/Kernel/m_querypoint.o \
../../Modules/Kernel/m_weighter.o \
../../Modules/Kernel/m_image_weighter.o \
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
	(cd ../../Modules/Kernel; make App=WPI)
	#Linesearch
	(cd ../../Modules/Optimization; make m_linesearcher_$(LineS).o)


exe : main.o ../04_GradientTest/sub_optimizer_init_loop.o #../../WPI/sub_modeling_gradient.o
	$(F90) $(FLAGF90) $(MOD) $(obj) $^ -o ../../exe/gradienttest_wpi_$(WaveEq)_$(Solver)_$(Order)
	(cd ../../exe; ln -sf gradienttest_wpi_$(WaveEq)_$(Solver)_$(Order) GradientTest_WPI )

clean :
	-rm main.o
	-rm ../04_GradientTest/sub_optimizer_init_loop.o
	-rm ../../WPI/sub_modeling_gradient.o

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
	(cd ../../Modules/Optimization; make clean)
	rm ../../Modules/Optimization/*.o
	#
	rm ../../mod/*.mod
	#
	make clean
