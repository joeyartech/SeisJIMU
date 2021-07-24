mod_dir=../../mod/
include ../../make.inc
include ../../setup.inc


src = \
m_linesearcher_$(LineS).f90 \
m_optimizer_$(Optim).f90 \
m_Optimization.f90

all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=$(mod_dir)%.mod)  $(mod_dir)m_linesearcher.mod $(mod_dir)m_optimizer.mod
