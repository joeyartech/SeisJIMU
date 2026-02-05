mod_dir=../mod/
include ../compiler.inc
include ../modules.inc


APP = RWI
app = $(shell echo $(APP) | tr A-Z a-z)

obj = \
../Modules/System/*.o \
../Modules/Etc/*.o \
../Modules/Signal/*.o \
../Modules/Modeling/m_model.o \
../Modules/Modeling/m_shotlist_$(ShotDec).o \
../Modules/Modeling/m_shot.o \
../Modules/Modeling/m_computebox.o \
../Modules/Modeling/m_correlate.o \
../Modules/Modeling/m_field.o \
../Modules/Modeling/m_cpml.o \
../Modules/Modeling/m_propagator_$(WaveEq)_$(Solver)_$(Order).o \
../Modules/Modeling/m_Modeling.o \
../Modules/Kernel/m_empirical.o \
../Modules/Kernel/m_parametrizer_$(Param).o \
../Modules/Kernel/m_querypoint.o \
../Modules/Kernel/m_separator.o \
../Modules/Kernel/m_weighter.o \
../Modules/Kernel/m_preconditioner.o \
../Modules/Kernel/m_fobjective.o \
../Modules/Kernel/m_Kernel.o \
../Modules/Optimization/m_linesearcher_$(LineS).o \
../Modules/Optimization/m_optimizer_common.o \
../Modules/Optimization/m_optimizer_$(Optim).o \
../Modules/Optimization/m_Optimization.o


all : mod exe

mod :
	#System
	+(cd ../Modules/System; $(MAKE))
	#Etc
	+(cd ../Modules/Etc; $(MAKE))
	#Signal
	+(cd ../Modules/Signal; $(MAKE))
	#Modeling
	+(cd ../Modules/Modeling; $(MAKE))
	#Kernel
	+(cd ../Modules/Kernel; $(MAKE) App=$(APP))
	#Optimization
	+(cd ../Modules/Optimization; $(MAKE) App=$(APP))

exe : sub_modeling_gradient.o main.o
	$(F90) $(FLAGF90) $(MOD) $(obj) $^ -o ../exe/$(app)_$(WaveEq)_$(Solver)_$(Order)_$(Param)_$(LineS)_$(Optim)_$(TAG)
	(cd ../exe; ln -sf $(app)_$(WaveEq)_$(Solver)_$(Order)_$(Param)_$(LineS)_$(Optim)_$(TAG) $(APP) )

clean :
	-rm main.o sub_modeling_gradient.o

cleanall :
	#System
	+(cd ../Modules/System; $(MAKE) clean)
	#Etc
	+(cd ../Modules/Etc; $(MAKE) clean)
	#Signal
	+(cd ../Modules/Signal; $(MAKE) clean)
	#Modeling
	+(cd ../Modules/Modeling; $(MAKE) clean)
	#Kernel
	+(cd ../Modules/Kernel; $(MAKE) clean)
	#Optimization
	+(cd ../Modules/Optimization; $(MAKE) clean)
	#$(APP)
	$(MAKE) clean
