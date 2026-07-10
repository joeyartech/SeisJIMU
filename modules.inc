#FPP += -D"checkpoint='F'"


############
# Modeling #
############
WaveEq=TE
#WaveEq=SH
#WaveEq=PSV
#WaveEq=DAS
#WaveEq=DASAC
Solver=FDSG
Order=O4
ShotDec=batching

##########
# Kernel #
##########
#Param=epsr-mur-sgma
#Param=celerity-permeability
Param=celerity-impedance
#Param=celerity-impedance_pseudotime
#Param=velocities-density
#Param=velocities-density_pseudotime
#Param=velocities-impedance
#Param=slowness-density
App=FWI

################
# Optimization #
################
LineS=Wolfe
#LineS=enumerate
#Optim=SD
Optim=NLCG
#Optim=lBFGS
