mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc


src = \
m_linesearcher_$(LineS).f90 \
m_optimizer_common.f90 \
m_optimizer_$(Optim).f90 \
m_Optimization.f90


all : $(src:.f90=.o)
	@printf "\n"

clean :
	-rm $(src:.f90=.o)
	@printf "\n"
	-rm $(src:%.f90=$(mod_dir)%.mod)  $(mod_dir)m_linesearcher.mod $(mod_dir)m_optimizer.mod $(mod_dir)m_optimization.mod
	@printf "\n"
