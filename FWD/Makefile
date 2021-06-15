mod_dir=../mod/
include ../make.inc
include ../setup.inc


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
	mpif90 $(FLAGF90) $^ $(MOD) -o ../exe/$(app)_$(WaveEq)_$(Solver)_$(Order)
	(cd ../exe; ln -sf $(app)_$(WaveEq)_$(Solver)_$(Order) $(APP) )

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