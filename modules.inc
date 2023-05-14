#FPP += -D"checkpoint='F'"


############
# Modeling #
############
WaveEq=BornEnv_reconstruct
#WaveEq=BornEnv_monoD
#WaveEq=AC2nd
#WaveEq=AC
#WaveEq=AC_VTI
#WaveEq=PSV
Solver=FD
#Solver=FDSG
#Order=O4
Order=O2
ShotDec=batching

##########
# Kernel #
##########
Param=vp-rho-tilD
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
Optim=SD
#Optim=NLCG
#Optim=lBFGS
