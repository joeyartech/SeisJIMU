mod_dir=../mod/
include ../compiler.inc
include ../modules.inc


obj = \
../Modules/System/*.o \
../Modules/Etc/*.o \
../Modules/Signal/*.o \
../Modules/Modeling/m_model.o \
../Modules/Modeling/m_shotlist_$(ShotDec).o \
../Modules/Modeling/m_shot.o \
../Modules/Modeling/m_computebox.o \
../Modules/Modeling/m_field.o \
../Modules/Modeling/m_cpml_$(Solver).o \
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
	mpif90 $(FLAGF90) $(MOD) $(obj) $^ -o ../exe/fwd_$(WaveEq)_$(Solver)_$(Order)
	(cd ../exe; ln -sf fwd_$(WaveEq)_$(Solver)_$(Order) FWD )

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
	#FWD
	make clean
