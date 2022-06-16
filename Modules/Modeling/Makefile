mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc


src = \
m_model.f90 \
m_shotlist_$(ShotDec).f90 \
m_shot.f90 \
m_computebox.f90 \
m_field.f90 \
m_cpml_$(Solver).f90 \
m_propagator_$(WaveEq)_$(Solver)_$(Order).f90 \
m_Modeling.f90


all : $(src:.f90=.o)
	@printf "\n"

clean :
	-rm *.o
	@printf "\n"
	-rm $(src:%.f90=$(mod_dir)%.mod) $(mod_dir)m_shotlist.mod $(mod_dir)m_cpml.mod $(mod_dir)m_propagator.mod $(mod_dir)m_modeling.mod
	@printf "\n"
