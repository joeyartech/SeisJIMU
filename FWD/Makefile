mod_dir=../mod/
include ../compiler.inc
include ../modules.inc


APP = FWD
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
../Modules/Modeling/m_Modeling.o


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

exe : main.o
	$(F90) $(FLAGF90) $(MOD) $(obj) $^ -o ../exe/$(app)_$(WaveEq)_$(Solver)_$(Order)_$(TAG)
	(cd ../exe; ln -sf $(app)_$(WaveEq)_$(Solver)_$(Order)_$(TAG) $(APP) )

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
	#$(APP)
	make clean
