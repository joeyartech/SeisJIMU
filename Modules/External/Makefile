include ../../make.inc


srcf77 = \
sgtsv.f

srcf90 = \
bessi0.f90 \
singleton.f90 \
sinc.f90 \
stoep.f90

all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=../../mod/%.mod)
