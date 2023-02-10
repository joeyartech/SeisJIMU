mod_dir=../../mod/
include ../../compiler.inc
include ../../modules.inc


src = \
m_hicks.f90 \
m_butterworth.f90 \
m_resampler.f90 \
m_wavelet.f90 \
m_matchfilter.f90 \
m_Hilbert.f90 \
m_Lpnorm.f90 \
m_pseudotime.f90 \
m_smoother_Laplacian_sparse.f90


all : $(src:.f90=.o)
	@printf "\n"

clean :
	-rm $(src:.f90=.o)
	@printf "\n"
	-rm $(src:%.f90=$(mod_dir)%.mod)  $(mod_dir)m_lpnorm.mod $(mod_dir)m_smoother_laplacian_sparse.mod
	@printf "\n"
