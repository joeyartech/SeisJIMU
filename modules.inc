#FPP += -D"checkpoint='F'"


#===============================================================================
# Modeling
#===============================================================================
# WaveEq — Wave equation (PDE) for forward/adjoint modeling
#
#   GPU-accelerated (nvfortran + OpenACC):
#     AC   — 2D/3D Acoustic isotropic
#     PSV  — 2D Elastic P-SV system
#
#   CPU-only (ifort / gfortran):
#     AC_VTI — 2D/3D Acoustic VTI
#     DAS    — Distributed Acoustic Sensing
#     AC2nd  — 2D/3D Acoustic 2nd-order formulation
#
# Uncomment exactly ONE WaveEq line below:
#===============================================================================
WaveEq=AC
#WaveEq=PSV
#WaveEq=AC_VTI
#WaveEq=DAS
#WaveEq=AC2nd
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
