mod_dir=../../mod/
include ../../make.inc
include ../../setup.inc


src = \
m_string.f90 \
m_either.f90 \
m_message.f90 \
m_arrayop.f90 \
m_setup.f90 \
m_format_su.f90

all : $(src:.f90=.o) m_mpienv.o

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=$(mod_dir)%.mod)
