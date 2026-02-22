# Compiler selection: ifort, gfortran, or nvfortran (GPU)
# Usage: make prepare Compiler=nvfortran  (then: make cleanall && make mod && make exe)
Compiler=ifort

include modules.inc

#===============================================================================
# Core build targets (existing workflow)
#===============================================================================

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

exe : fwd fwi rwi rtm

fwd :
	+(cd FWD; $(MAKE))
	@printf "\n"

fwi :
	+(cd FWI; $(MAKE))
	@printf "\n"

rwi :
	+(cd RWI; $(MAKE))
	@printf "\n"

rtm :
	+(cd RTM; $(MAKE))
	@printf "\n"

#===============================================================================
# Convenience targets: full CPU or GPU build in one command
#===============================================================================

cpu :
	@echo "=== Building ALL apps with CPU (ifort) ==="
	$(MAKE) prepare Compiler=ifort
	$(MAKE) cleanall
	$(MAKE) mod
	$(MAKE) exe

gpu :
	@echo "=== Building ALL apps with GPU (nvfortran) ==="
	$(MAKE) prepare Compiler=nvfortran
	$(MAKE) cleanall
	$(MAKE) mod
	$(MAKE) exe

#===============================================================================
# Per-app CPU/GPU targets
#===============================================================================

fwd-cpu :
	$(MAKE) prepare Compiler=ifort
	$(MAKE) cleanall
	$(MAKE) mod
	+(cd FWD; $(MAKE))

fwd-gpu :
	$(MAKE) prepare Compiler=nvfortran
	$(MAKE) cleanall
	$(MAKE) mod
	+(cd FWD; $(MAKE))

fwi-cpu :
	$(MAKE) prepare Compiler=ifort
	$(MAKE) cleanall
	$(MAKE) mod
	+(cd FWI; $(MAKE))

fwi-gpu :
	$(MAKE) prepare Compiler=nvfortran
	$(MAKE) cleanall
	$(MAKE) mod
	+(cd FWI; $(MAKE))

rwi-cpu :
	$(MAKE) prepare Compiler=ifort
	$(MAKE) cleanall
	$(MAKE) mod
	+(cd RWI; $(MAKE))

rwi-gpu :
	$(MAKE) prepare Compiler=nvfortran
	$(MAKE) cleanall
	$(MAKE) mod
	+(cd RWI; $(MAKE))

rtm-cpu :
	$(MAKE) prepare Compiler=ifort
	$(MAKE) cleanall
	$(MAKE) mod
	+(cd RTM; $(MAKE))

rtm-gpu :
	$(MAKE) prepare Compiler=nvfortran
	$(MAKE) cleanall
	$(MAKE) mod
	+(cd RTM; $(MAKE))

#===============================================================================
# Info target: print current build configuration
#===============================================================================

info :
	@echo "=========================================="
	@echo " SeisJIMU Build Configuration"
	@echo "=========================================="
	@echo " Compiler setting : $(Compiler)"
	@if [ -L compiler.inc ]; then \
		echo " compiler.inc     : $$(readlink compiler.inc)"; \
	elif [ -f compiler.inc ]; then \
		echo " compiler.inc     : (regular file, not a symlink)"; \
	else \
		echo " compiler.inc     : NOT SET (run 'make prepare' first)"; \
	fi
	@echo "------------------------------------------"
	@echo " WaveEq           : $(WaveEq)"
	@echo " Solver            : $(Solver)"
	@echo " Order             : $(Order)"
	@echo " ShotDec           : $(ShotDec)"
	@echo " Param             : $(Param)"
	@echo " Optim             : $(Optim)"
	@echo " LineS             : $(LineS)"
	@echo "------------------------------------------"
	@echo " Executables in exe/:"
	@if [ -d exe ]; then ls -1 exe/ 2>/dev/null || echo "  (empty)"; else echo "  (exe/ not created yet)"; fi
	@echo "=========================================="

#===============================================================================
# Cleaning
#===============================================================================

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
