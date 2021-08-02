mod_dir=../../mod/
include ../../make.inc
include ../../setup.inc


src = \
m_either.f90 \
m_string.f90 \
m_mpienv.f90 \
m_message.f90 \
m_arrayop.f90 \
m_setup.f90 \
m_sysio.f90 \
m_checkpoint.f90 \
m_suformat.f90 \
m_System.f90

all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=$(mod_dir)%.mod) $(mod_dir)m_system.mod
