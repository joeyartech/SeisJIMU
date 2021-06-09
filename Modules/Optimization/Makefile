include ../../make.inc

Preco=zpower
#Preco=scaling
LineS=Wolfe
#Optim=NLCG
Optim=LBFGS

src = \
m_preconditioner_$(Preco).f90 \
m_linesearcher_$(LineS).f90 \
m_optimizer_$(Optim).f90 \

all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=../../mod/%.mod)
