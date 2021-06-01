# SeisJIMU


	███████ ███████ ██ ███████      ██ ██ ███    ███ ██    ██ 
	██      ██      ██ ██           ██ ██ ████  ████ ██    ██ 
	███████ █████   ██ ███████      ██ ██ ██ ████ ██ ██    ██ 
	     ██ ██      ██      ██ ██   ██ ██ ██  ██  ██ ██    ██ 
	███████ ███████ ██ ███████  █████  ██ ██      ██  ██████
	                                          
SeisJIMU is a module based seismic modeling and inversion (FWI) code,
licensed under GNU GENERAL PUBLIC LICENSE Version 3 (see LICENSE).

~Version 1
 - Authors: 
	 - Wei Zhou  (joeywzhou1986@gmail.com)
	 - David Lumley (david.lumley@utdallas.edu)
 - Language: Fortran
 - Libraries: OpenMP, MPI
 - Compilers: gfortran (Ver>=5), ifort
<!---  - Tested on CPU: Intel(R) Xeon(R) Gold 6140 CPU @ 2.30GHz --->

# Quick start
## Prerequisite
- Linux OS
- Makefile
- Fortran compiler (e.g. ifort, gfortran (Ver>=5))
- OpenMP and MPI
<!-- - [Seismic Unix](https://github.com/JohnWStockwellJr/SeisUnix) for data IO and plotting -->

## Compilation
1. Edit Makefile  
    a. Choose compiler and associated compilation flags. By default we use ifort (Intel Fortran compiler). FLAGF90 and FLAGF77 are default flags to compile Fortran90 (and later) and Fortran77 (and earlier) codes, respectively.  
    b. Choose modules:
    - WaveEq: which wave equation (PDE) to be solved for forward and adjoint modeling
        - AC:     2D/3D ACoustic isotropic; 
        - AC_VTI: 2D/3D ACoustic VTI, 
        - PSV:    2D elastic P-Sv system
    - Solver: which solver to solve the PDE (so far only FDSG: Finite-Difference Staggered Grid is implemented)
    - Domain: modeling in which domain (so far only TimeDomain is implemented)
    - Norm:   which norm to be used for inversion (so far only L2 norm is implemented)
    - Param:  which medium parameterization to be used for inversion
        - velocities-density
        - velocities-impedance
        - slowness-density
    - Preco:  gradient preconditioner
    - LineS:  line search method inside optimization loop
    - Optim:  optimizer
        - LBFGS:  limited memory BFGS quasi-Newton method
        - NLCG:   nonlinear Conjugate Gradient

3. Compilation  
    a. create directory:
``` $ make dir ```  
    b. compile:
    - For FWD (forward modeling):
``` $ make fwd ```
    - For FWI (full waveform inversion):
``` $ make fwi ```

    c. The executable binary can be found in ./exe. The name of the executable is a concentration of fwd or fwi and selected modules, e.g., `fwi_AC_FDSG_L2_velocities-density_zpower_Wolfe_LBFGS`  
This allows existence of multiple executables compiled with different modules.  
The newest compiled executable has a nickname simply as FWD or FWI.

3. Cleaning
- ``` $ make clean ```
    will clean all temporary files (e.g. .mod, .o). This allows a total recompilation of executables.
- ``` $ make cleanall ```
    will clean all temporary files as well as executables (not needed in most cases).

## Running programs
A handy manual is available by just running the executable in the terminal:
``` $ ./exe/FWD ```
``` $ ./exe/FWi ```

1. Create/Edit `setup.in`:
  Check the manual to see which items are mandatory or optional in `setup.in`. Modify it for your problem.
2. Decide how many (logical) CPUs and OpenMP threads you need:
  - CPU number = how many shots are processed in parallel (no need to be a factor of total number of shots)
  - Thread number = how many threads are solving the PDE in parallel (a simple version of domain decomposition)
  - CPU number * Thread number = Total number of allocated physical cores = Number of allocated nodes * Number cores per node
3. Launch program and environment in the terminal:
  a. set unlimited stack memory size (necessary for large problems):
    ``` $ ulimit -s unlimited ```
  b. specify number of OpenMP threads:
    ``` $ export OMP_NUM_THREADS=$ThreadNumber ``` (use source in csh like shell)
  c. run program with MPI wrapper, followed by setup.in file:
    ``` $ mpirun -np $CPUnumber ./exe/FWD setup.in ```
    ``` $ mpirun -np $CPUnumber ./exe/FWI setup.in ```


License
----
GNU GENERAL PUBLIC LICENSE Version 3 (GPLv3.0)
See LICENSE


References
----
If you publish results using this code, please acknowledge and reference our paper:
- [Wei Zhou and David Lumley (2021), Central-difference time-lapse 4D seismic full waveform inversion, _Geophysics_.](https://library.seg.org/doi/10.1190/geo2019-0834.1)
- [Wei Zhou and David Lumley (2021), Non-repeatability Effects on Time-Lapse 4D Seismic Full Waveform Inversion for Ocean-Bottom Node Data, _Geophysics_.](https://library.seg.org/doi/10.1190/geo2020-0577.1)
