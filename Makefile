OMP=openmp
##Note: ifort does not implement openmp in debug mode

#OPTRPT=-qopt-report=5  -qopt-report-file=stderr -qopt-report-phase=vec

#GIT_COMMIT:=$(shell git rev-parse HEAD)
GIT_COMMIT:=$(shell git rev-parse --short HEAD)
GIT_BRANCH:=$(shell git rev-parse --abbrev-ref HEAD)

# #for Intel
# FLAGF90= -xHost -O3 -Ofast -ipo -parallel -q$(OMP) -fno-alias -traceback -assume byterecl   $(OPTRPT)
# # FLAGF90= -xCORE-AVX512 -O3 -Ofast -ipo -parallel -q$(OMP) -fno-alias -traceback -assume byterecl   $(OPTRPT)
# # FLAGF90= -g -debug -check all -check noarg_temp_created -W1 -WB -q$(OMP) -traceback -assume byterecl
# FLAGF77= $(FLAGF90)
# MOD= -module $(DIR)mod
# FPP=-fpp -DINTEL -D"commit='$(GIT_COMMIT)'" -D"branch='$(GIT_BRANCH)'"

#for GNU
FLAGF77= -Ofast -f$(OMP)                         -fbacktrace  #-flto
FLAGF90= -Ofast -f$(OMP) -ffree-line-length-none -fbacktrace  -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans
# FLAGF77= -g -Wall -f$(OMP) -fcheck=all                         -fbacktrace
# FLAGF90= -g -Wall -f$(OMP) -fcheck=all -ffree-line-length-none -fbacktrace
MOD=-J $(DIR)mod
FPP=-cpp -DGNU -D"commit='$(GIT_COMMIT)'" -D"branch='$(GIT_BRANCH)'"


DIR=./

dir :
	mkdir $(DIR)mod exe


.SUFFIX : %.f %.f90 %.o


%.o : %.f90
	mpif90 $(FPP) $(FLAGF90) $(MOD) -c $^ -o $@

%.o : %.f
	mpif77        $(FLAGF77)        -c $^ -o $@



system = \
$(DIR)Modules/System/m_global.f90 \
$(DIR)Modules/System/m_string.f90 \
$(DIR)Modules/System/m_mpienv.f90 \
$(DIR)Modules/System/m_message.f90 \
$(DIR)Modules/System/m_setup.f90 \
$(DIR)Modules/System/m_arrayop.f90 \
$(DIR)Modules/System/m_suformat.f90

externf77 = \
$(DIR)Modules/External/sgtsv.f
externf90 = \
$(DIR)Modules/External/bessi0.f90 \
$(DIR)Modules/External/singleton.f90 \
$(DIR)Modules/External/sinc.f90 \
$(DIR)Modules/External/stoep.f90

OBJ0=$(system:.f90=.o) $(externf77:.f=.o) $(externf90:.f90=.o)


#######
# FWD #
#######

WaveEq=AC
# WaveEq=AC_VTI
# WaveEq=PSV
Solver=FDSG
Domain=TimeDomain

signal_fwd = \
$(DIR)Modules/SignalProcessing/m_hicks.f90 \
$(DIR)Modules/SignalProcessing/m_butterworth.f90 \
$(DIR)Modules/SignalProcessing/m_resamp.f90

modeling = \
$(DIR)Modules/Modeling/m_model.f90 \
$(DIR)Modules/Modeling/m_gen_acqui.f90 \
$(DIR)Modules/Modeling/$(Domain)/m_shotlist.f90 \
$(DIR)Modules/Modeling/$(Domain)/m_shot.f90 \
$(DIR)Modules/Modeling/$(Domain)/m_computebox.f90 \
$(DIR)Modules/Modeling/$(Domain)/$(Solver)/m_field_$(WaveEq).f90 \
$(DIR)Modules/Modeling/$(Domain)/$(Solver)/m_boundarystore.f90 \
$(DIR)Modules/Modeling/$(Domain)/$(Solver)/m_propagator.f90

OBJ_FWD=$(OBJ0) $(signal_fwd:.f90=.o) $(modeling:.f90=.o)


#######
# FWI #
#######

Norm=L2
Param=velocities-density
#Param=velocities-impedance
#Param=slowness-density
Preco=zpower
LineS=Wolfe
#Optim=NLCG
Optim=LBFGS

signal_fwi = \
$(DIR)Modules/SignalProcessing/m_weighter_polygon.f90 \
$(DIR)Modules/SignalProcessing/m_weighter_table.f90 \
$(DIR)Modules/SignalProcessing/m_matchfilter.f90 \
$(DIR)Modules/SignalProcessing/m_smoother_laplacian_sparse.f90

gradient = \
$(DIR)Modules/Gradient/m_objectivefunc_$(Norm).f90 \
$(DIR)Modules/Gradient/m_gradient.f90 \
$(DIR)Modules/Gradient/m_parameterization_$(Param).f90 \

optimization = \
$(DIR)Modules/Optimization/m_preconditioner_$(Preco).f90 \
$(DIR)Modules/Optimization/m_linesearcher_$(LineS).f90 \
$(DIR)Modules/Optimization/m_optimizer_$(Optim).f90 \

OBJ_FWI=$(OBJ_FWD) $(signal_fwi:.f90=.o) $(gradient:.f90=.o) $(optimization:.f90=.o)


###############################################################################################

fwd : $(OBJ_FWD) FWD/main.o $(DIR)mod exe
	mpif90 $(FLAGF90) $(OBJ_FWD) FWD/main.o $(MOD)  -o exe/fwd_$(WaveEq)_$(Solver)
	(cd exe; ln -sf fwd_$(WaveEq)_$(Solver) FWD)

fwi : $(OBJ_FWI) FWI/main.o $(DIR)mod exe
	mpif90 $(FLAGF90) $(OBJ_FWI) FWI/main.o $(MOD)  -o exe/fwi_$(WaveEq)_$(Solver)_$(Norm)_$(Param)_$(Preco)_$(LineS)_$(Optim)
	(cd exe; ln -sf fwi_$(WaveEq)_$(Solver)_$(Norm)_$(Param)_$(Preco)_$(LineS)_$(Optim) FWI)


clean :
	-rm FWD/*.o FWI/*.o
	-rm $(OBJ_FWI)
	-rm $(DIR)mod/*

cleanall :
	make clean
	-rm exe/*
