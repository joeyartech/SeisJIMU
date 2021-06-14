mod_dir=../../mod/
include ../../make.inc


LineS=Wolfe
Optim=NLCG
#Optim=LBFGS

src = \
m_preconditioner.f90 \
m_linesearcher_$(LineS).f90 \
m_optimizer_$(Optim).f90 \

all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=$(mod_dir)%.mod)  $(mod_dir)m_linesearcher.mod $(mod_dir)m_optimizer.mod
