include ../../make.inc

Param=velocities-density
#Param=velocities-impedance
#Param=slowness-density

src = \
m_weighter.f90 \
m_fobjective.f90 \
m_nabla_FWI.f90 \
m_parameterization_$(Param).f90
# m_separator.f90 \
# m_nabla_RWI.f90 \

all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=../../mod/%.mod) ../../mod/m_parameterization.mod
