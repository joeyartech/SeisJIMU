include ../../make.inc


src = \
m_hicks.f90 \
m_butterworth.f90 \
m_resampler.f90 \
m_wavelets.f90 \
m_matchfilter.f90 \
m_smoother_laplacian_sparse.f90


all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=../../mod/%.mod)
