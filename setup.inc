#FPP += -D"checkpoint='F'"


############
# Modeling #
############
WaveEq=AC
# WaveEq=AC_VTI
# WaveEq=PSV
Solver=FDSG
Order=O4
ShotDec=batching

##########
# Kernel #
##########
Param=velocities-density
#Param=velocities-impedance
#Param=slowness-density
app=fwi
APP=FWI

################
# Optimization #
################
LineS=Wolfe
Optim=NLCG
#Optim=LBFGS