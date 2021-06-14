
all : prepare system etc signal modeling kernel optimization

compiler=gfortran

prepare :
	mkdir mod exe
	ln -sf make.inc_$(compiler) make.inc

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

fwd :
	(cd FWD; make)


cleanall :
	(cd Modules/System; make clean)
	(cd Modules/Etc; make clean)
	(cd Modules/Signal; make clean)
	(cd Modules/Modeling; make clean)
	(cd Modules/Kernel; make clean)
	(cd Modules/Optimization; make clean)
	-rm FWD/*.o FWI/*.o
