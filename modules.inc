#FPP += -D"checkpoint='F'"


############
# Modeling #
############
WaveEq=PSV
WaveEq=PSV
Solver=FDSG
Order=O4
ShotDec=batching

##########
# Kernel #
##########
Param=velocities-density
#Param=velocities-impedance
#Param=velocities-impedance_pseudotime
App=FWI

################
# Optimization #
################
LineS=Wolfe
#LineS=enumerate
Optim=lBFGS
#Optim=NLCG
Optim=lBFGS
