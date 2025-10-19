#FPP += -D"checkpoint='F'"


############
# Modeling #
############
#WaveEq=PSV
#WaveEq=VSS
#WaveEq=VSR

#WaveEq=PSV2nd
WaveEq=VSS2nd

Solver=FDSG
Order=O4
ShotDec=batching

##########
# Kernel #
##########
#Param=vp-rho
Param=velocities-density
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
#Optim=NLCG
Optim=lBFGS
