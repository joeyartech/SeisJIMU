#FPP += -D"checkpoint='F'"


############
# Modeling #
############
WaveEq=2ndOrd_AC
#WaveEq=AC
#WaveEq=AC_VTI
# WaveEq=PSV
Solver=FDSG
#Order=O4
Order=O2
ShotDec=batching

##########
# Kernel #
##########
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
Optim=SD
#Optim=NLCG
#Optim=lBFGS
