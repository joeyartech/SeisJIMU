#FPP += -D"checkpoint='F'"


############
# Modeling #
############
#WaveEq=AC
WaveEq=PSV
Solver=FDSG
Order=O4
ShotDec=batching

##########
# Kernel #
##########
#Param=velocities-density
#Param=velocities-impedance
Param=vel-imp_Abdullah
#Param=velocities-impedance_pseudotime
App=FWI

################
# Optimization #
################
LineS=Wolfe
#LineS=enumerate
#Optim=SD
#Optim=NLCG
Optim=lBFGS
