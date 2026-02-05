Compiler=ifort

prepare :
	@mkdir -p mod exe
	ln -sf compiler.inc_$(Compiler) compiler.inc

mod : system etc signal modeling kernel optimization

system : 
	+(cd Modules/System; $(MAKE) )

etc :
	+(cd Modules/Etc; $(MAKE))

signal :
	+(cd Modules/Signal; $(MAKE))

modeling :
	+(cd Modules/Modeling; $(MAKE))

kernel :
	+(cd Modules/Kernel; $(MAKE))

optimization :
	+(cd Modules/Optimization; $(MAKE))

exe : fwd fwi rwi

fwd :
	+(cd FWD; $(MAKE))
	@printf "\n"

fwi :
	+(cd FWI; $(MAKE))
	@printf "\n"

rwi :
	+(cd RWI; $(MAKE))

rtm :
	+(cd RTM; $(MAKE))
	@printf "\n"


clean :
	-rm FWD/*.o FWI/*.o RWI/*.o RTM/*.o

cleanmod :
	(cd Modules/System; $(MAKE) clean)
	(cd Modules/Etc; $(MAKE) clean)
	(cd Modules/Signal; $(MAKE) clean)
	(cd Modules/Modeling; $(MAKE) clean)
	(cd Modules/Kernel; $(MAKE) clean)
	(cd Modules/Optimization; $(MAKE) clean)
	@if [ -d mod ]; then rm -f mod/*.mod mod/*.smod; fi

cleanall : cleanmod clean
