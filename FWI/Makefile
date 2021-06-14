mod_dir=../mod/
include ../make.inc


WaveEq=AC
# WaveEq=AC_VTI
# WaveEq=PSV
Solver=FDSG
Order=O4
LineS=Wolfe
Optim=NLCG

app=fwi
APP=FWI

mod :
	#System
	(cd ../Modules/System; make)
	#Etc
	(cd ../Modules/Etc; make)
	#Signal
	(cd ../Modules/Signal; make)
	#Modeling
	(cd ../Modules/Modeling; make WaveEq=$(WaveEq) Solver=$(Solver) Order=$(Order) )
	#Kernel
	(cd ../Modules/Kernel; make Param=$(Param) APP=$(APP))
	#Optimization
	(cd ../Modules/Optimization; make LineS=$(LineS) Optim=$(Optim))

exe : main.o
	mpif90 $(FLAGF90) $^ $(MOD) -o ../exe/$(app)_$(WaveEq)_$(Solver)_$(Order)_$(LineS)_$(Optim)
	(cd ../exe; ln -sf $(app)_$(WaveEq)_$(Solver)_$(Order)_$(LineS)_$(Optim) $(APP) )

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