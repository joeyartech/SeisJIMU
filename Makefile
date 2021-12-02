Compiler=gfortran

prepare :
	-mkdir mod exe
	ln -sf compiler.inc_$(Compiler) compiler.inc

mod : system etc signal modeling kernel optimization

system : 
	(cd Modules/System; make )

etc :
	(cd Modules/Etc; make)

signal :
	(cd Modules/Signal; make)

modeling :
	(cd Modules/Modeling; make)

kernel :
	(cd Modules/Kernel; make)

optimization :
	(cd Modules/Optimization; make)

exe : fwd fwi wpi

fwd :
	(cd FWD; make)
	@printf "\n"

fwi :
	(cd FWI; make)
	@printf "\n"

wpi :
	(cd WPI; make)
	@printf "\n"


clean :
	-rm FWD/*.o FWI/*.o WPI/*.o

cleanmod :
	(cd Modules/System; make clean)
	(cd Modules/Etc; make clean)
	(cd Modules/Signal; make clean)
	(cd Modules/Modeling; make clean)
	(cd Modules/Kernel; make clean)
	(cd Modules/Optimization; make clean)
	(cd mod; rm -r *)

cleanall : cleanmod clean
