include ../../make.inc


src = \
# $(DIR)Modules/SignalProcessing/m_hicks.f90 \
# $(DIR)Modules/SignalProcessing/m_butterworth.f90 \
# $(DIR)Modules/SignalProcessing/m_resamp.f90 \
# $(DIR)Modules/SignalProcessing/m_weighter_polygon.f90 \
# $(DIR)Modules/SignalProcessing/m_weighter_table.f90 \
# $(DIR)Modules/SignalProcessing/m_matchfilter.f90 \
# $(DIR)Modules/SignalProcessing/m_smoother_laplacian_sparse.f90

all : $(src:.f90=.o)

clean :
	-rm $(src:.f90=.o) 
	-rm $(src:%.f90=../../mod/%.mod)
