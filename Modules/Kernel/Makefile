include ../../make.inc

Param=velocities-density
#Param=velocities-impedance
#Param=slowness-density

Approach=FWI

src = \
m_fobjective.f90 \
m_gradient_modeling_$(Approach).f90 \
m_parameterization_$(Param).f90 \
m_weighter_polygon.f90 \
m_weighter_table.f90

all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=../../mod/%.mod)
