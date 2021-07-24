mod_dir=../../mod/
include ../../make.inc
include ../../setup.inc


ifeq '$(APP)' 'RTM'
	src = 
endif

src = \
	m_weighter.f90 \
	m_parametrizer_$(Param).f90 \
	m_fobjective.f90 \
	m_nabla.f90 \
	m_preconditioner.f90 \
	m_Kernel.f90

ifeq '$(APP)' 'FWI'
endif

ifeq '$(APP)' 'RWI'
	src += m_separator.f90
endif

all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=$(mod_dir)%.mod) $(mod_dir)m_parametrizer.mod $(mod_dir)m_nabla.mod
