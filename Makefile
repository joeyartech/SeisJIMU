
all : dir system external signal modeling kernel optimization

dir :
	mkdir mod exe

system : 
	(cd Modules/System; make )

external :
	(cd Modules/External; make)

signal :
	(cd Modules/Signal; make)

modeling :
	(cd Modules/Modeling; make)

kernel :
	(cd Modules/Kernel; make)

optimization :
	(cd Modules/Optimization; make)

# fwd : $(OBJ_FWD) FWD/main.o $(DIR)mod exe
# 	mpif90 $(FLAGF90) $(OBJ_FWD) FWD/main.o $(MOD)  -o exe/fwd_$(WaveEq)_$(Solver)
# 	(cd exe; ln -sf fwd_$(WaveEq)_$(Solver) FWD)
# 
# fwi : $(OBJ_FWI) FWI/main.o $(DIR)mod exe
# 	mpif90 $(FLAGF90) $(OBJ_FWI) FWI/main.o $(MOD)  -o exe/fwi_$(WaveEq)_$(Solver)_$(Norm)_$(Param)_$(Preco)_$(LineS)_$(Optim)
# 	(cd exe; ln -sf fwi_$(WaveEq)_$(Solver)_$(Norm)_$(Param)_$(Preco)_$(LineS)_$(Optim) FWI)


clean :
	(cd Modules/System; make clean)
	(cd Modules/External; make clean)
	(cd Modules/Signal; make clean)
	(cd Modules/Modeling; make clean)
	(cd Modules/Kernel; make clean)
	(cd Modules/Optimization; make clean)
	-rm FWD/*.o FWI/*.o
