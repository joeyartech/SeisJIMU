mod_dir=../mod/
include ../make.inc


WaveEq=AC
# WaveEq=AC_VTI
# WaveEq=PSV
Solver=FDSG
Order=O4

app=fwd
APP=FWD

mod :
	#System
	(cd ../Modules/System; make)
	#Etc
	(cd ../Modules/Etc; make)
	#Signal
	(cd ../Modules/Signal; make)
	#Modeling
	(cd ../Modules/Modeling; make WaveEq=$(WaveEq) Solver=$(Solver) Order=$(Order) )

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