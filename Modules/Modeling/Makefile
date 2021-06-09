include ../../make.inc

WaveEq=AC
# WaveEq=AC_VTI
# WaveEq=PSV
Solver=FDSG_O4

src = \
m_model.f90 \
m_shotlist.f90 \
m_shot.f90 \
m_computebox.f90 \
m_field.f90 \
m_boundarystore.f90 \
m_propagator_$(WaveEq)_$(Solver).f90

all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=../../mod/%.mod)
