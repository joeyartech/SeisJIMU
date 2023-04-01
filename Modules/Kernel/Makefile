mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc


ifeq '$(App)' 'RTM'
	src = 
endif

src = \
	m_parametrizer_$(Param).f90 \
	m_querypoint.f90 \
	m_weighter.f90 \
	m_preconditioner.f90 \
	m_fobjective.f90 \
	m_Kernel.f90

ifeq '$(App)' 'FWI'
endif

ifeq '$(App)' 'RWI'
	src += m_separator.f90
endif

ifeq '$(App)' 'WPI'
	src += m_image_weighter.f90
endif


all : $(src:.f90=.o)
	@echo $(src)
	@printf "\n"

clean :
	-rm *.o
	@printf "\n"
	-rm $(src:%.f90=$(mod_dir)%.mod) $(mod_dir)m_parametrizer.mod $(mod_dir)m_kernel.mod $(mod_dir)m_separator.mod $(mod_dir)m_image_weighter.mod
	@printf "\n"
