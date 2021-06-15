mod_dir=../../mod/
include ../../make.inc
include ../../setup.inc


srcf77 = \
sgtsv.f

srcf90 = \
m_math.f90 \
singleton.f90

all : $(srcf77:.f=.o) $(srcf90:.f90=.o)

clean :
	-rm $(srcf77:.f=.o) $(srcf90:.f90=.o)
	-rm $(srcf90:%.f90=$(mod_dir)%.mod)
