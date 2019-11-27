# #for Intel
# FLAGF90= -xHost -O3 -Ofast -ipo -parallel -qopenmp -traceback -assume byterecl
# #FLAGF90= -g -debug -W1 -WB -traceback -assume byterecl
# FLAGF77= $(FLAGF90)
# MOD= -module $(DIR)mod
#for GNU
FLAGF77= -Ofast -fopenmp -fbacktrace #-flto
FLAGF90= -Ofast -fopenmp -ffree-line-length-none -fbacktrace #-flto
# FLAGF77= -g -Wall -fcheck=all                         -fopenmp -fbacktrace
# FLAGF90= -g -Wall -fcheck=all -ffree-line-length-none -fopenmp -fbacktrace
MOD=-J $(DIR)mod


DIR=../../

dir :
	mkdir $(DIR)mod


.SUFFIX : %.f %.f90 %.o

%.o : %.f90
	mpif90 $(FLAGF90) $(MOD) -c $^ -o $@

%.o : %.f
	mpif77 $(FLAGF77)        -c $^ -o $@



system = \
$(DIR)Modules/System/m_mpienv.f90 \
$(DIR)Modules/System/m_message.f90 \
$(DIR)Modules/System/m_sysio.f90 \
$(DIR)Modules/System/m_arrayop.f90 \
$(DIR)Modules/System/m_suformat.f90

externf77 = \
$(DIR)Modules/External/F77/sgtsv.f
externf90 = \
$(DIR)Modules/External/F90/singleton.f90

signal = \
$(DIR)Modules/SignalProcessing/m_hicks.f90 \
$(DIR)Modules/SignalProcessing/m_weighting.f90 \
$(DIR)Modules/SignalProcessing/m_butterworth.f90 \
$(DIR)Modules/SignalProcessing/m_matchfilter.f90 \
$(DIR)Modules/SignalProcessing/m_laplacian_smoothing_sparse.f90

modeling = \
$(DIR)Modules/Modeling/m_model.f90 \
$(DIR)Modules/Modeling/m_gen_acquisition.f90 \
$(DIR)Modules/Modeling/m_gen_wavelet.f90 \
$(DIR)Modules/Modeling/TimeDomain/m_shotlist.f90 \
$(DIR)Modules/Modeling/TimeDomain/m_shot.f90 \
$(DIR)Modules/Modeling/TimeDomain/m_computebox.f90


OBJ=$(system:.f90=.o) $(externf77:.f=.o) $(externf90:.f90=.o) $(signal:.f90=.o) $(modeling:.f90=.o)



add_suheader : $(OBJ) main.o $(DIR)mod
	mpif90 $(MOD) $(FLAGF90) $(OBJ) main.o   -o $(DIR)/exe/ADD_SUHEADER

clean :
	-rm $(OBJ)
	-rm $(DIR)mod/*
