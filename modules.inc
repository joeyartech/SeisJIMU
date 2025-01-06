#FPP += -D"checkpoint='F'"


############
# Modeling #
############
#WaveEq=ExactEnv
#WaveEq=BornEnv_reconstruct
#WaveEq=BornEnv_monoD
#WaveEq=AC
#WaveEq=AC_VTI
#WaveEq=PSV

WaveEq=AC2nd
WaveEq=constrho_AC2nd_noPML

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
App=PFEI

################
# Optimization #
################
LineS=Wolfe
#LineS=enumerate
#Optim=SD
Optim=NLCG
#Optim=lBFGS
