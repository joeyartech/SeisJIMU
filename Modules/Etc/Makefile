mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc


srcf77 = \
sgtsv.f

srcf90 = \
m_math.f90 \
findelement.f90 \
singleton.f90


all : $(srcf77:.f=.o) $(srcf90:.f90=.o)
	@printf "\n"

clean :
	-rm $(srcf77:.f=.o) $(srcf90:.f90=.o)
	@printf "\n"
	-rm $(srcf90:%.f90=$(mod_dir)%.mod)
	@printf "\n"
