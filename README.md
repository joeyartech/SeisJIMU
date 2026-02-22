# SeisJIMU

SeisJIMU is a modular seismic modeling and inversion code written in Fortran.
It supports forward modeling (FWD), full waveform inversion (FWI), reverse time
migration (RTM), and reflection waveform inversion (RWI). GPU acceleration is
available for selected wave equations via OpenACC (NVIDIA HPC SDK).

- **Authors**: Wei Zhou, David Lumley
- **Language**: Fortran 95 / 2003
- **Parallelism**: MPI + OpenMP (CPU), MPI + OpenACC (GPU)
- **License**: GNU General Public License v3 (GPLv3) -- see `LICENSE`

## Prerequisites

### CPU builds (ifort or gfortran)

- Linux OS
- Intel Fortran (`ifort` / `ifx`) **or** GNU Fortran (`gfortran >= 5`)
- MPI library (e.g. Intel MPI, OpenMPI, MPICH)
- OpenMP (included with both compilers)

### GPU builds (nvfortran)

- Everything above, **plus**:
- [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) (`nvfortran`, `mpifort`)
- NVIDIA GPU with compute capability >= 7.0 (V100, A100, H100, ...)

## Quick Build

One-liner commands that run the full pipeline (`prepare` -> `cleanall` -> `mod` -> `exe`):

```bash
# CPU build (ifort)
make cpu

# GPU build (nvfortran + OpenACC)
make gpu
```

Per-app shortcuts are also available:

```bash
make fwd-cpu    make fwd-gpu
make fwi-cpu    make fwi-gpu
make rtm-cpu    make rtm-gpu
make rwi-cpu    make rwi-gpu
```

## Detailed Build

### Step 1 -- Select the compiler

```bash
make prepare Compiler=ifort       # Intel Fortran (CPU)
make prepare Compiler=gfortran    # GNU Fortran (CPU)
make prepare Compiler=nvfortran   # NVIDIA HPC SDK (GPU)
```

This creates the `mod/` and `exe/` directories and symlinks
`compiler.inc` -> `compiler.inc_<Compiler>`.

### Step 2 -- Configure compiler.inc (optional)

Each `compiler.inc_*` file contains a `BUILD_TYPE` variable. Uncomment the one
you need:

| Compiler    | BUILD_TYPE options                                          |
|-------------|-------------------------------------------------------------|
| ifort       | `release`, `release_avx2` (default), `release_intel`, `debug`, `profile` |
| gfortran    | `release` (default), `release_native`, `debug`, `profile`  |
| nvfortran   | `release_gpu` (default, managed memory), `release_gpu_explicit`, `debug_gpu` |

### Step 3 -- Configure modules.inc

Edit `modules.inc` to select the wave equation and other modules. Uncomment
exactly one `WaveEq` line:

```makefile
WaveEq=AC        # Acoustic isotropic          (CPU + GPU)
#WaveEq=PSV      # Elastic P-SV                (CPU + GPU)
#WaveEq=AC_VTI   # Acoustic VTI                (CPU only)
#WaveEq=DAS      # Distributed Acoustic Sensing (CPU only)
#WaveEq=AC2nd    # Acoustic 2nd-order           (CPU only)
```

### Step 4 -- Build

```bash
make cleanall    # clean previous build artifacts
make mod         # compile modules
make exe         # compile all applications (FWD, FWI, RWI, RTM)
```

Or build a single application:

```bash
make fwd         # forward modeling only
make fwi         # full waveform inversion only
make rtm         # reverse time migration only
make rwi         # reflection waveform inversion only
```

### Inspect current configuration

```bash
make info
```

Prints the active compiler, wave equation, solver, and lists executables in
`exe/`.

## GPU Support Matrix

| WaveEq   | Description                       | CPU | GPU (OpenACC) |
|----------|-----------------------------------|:---:|:-------------:|
| `AC`     | 2D/3D Acoustic isotropic          | Yes | Yes           |
| `PSV`    | 2D Elastic P-SV                   | Yes | Yes           |
| `AC_VTI` | 2D/3D Acoustic VTI               | Yes | No            |
| `DAS`    | Distributed Acoustic Sensing      | Yes | No            |
| `AC2nd`  | 2D/3D Acoustic 2nd-order          | Yes | No            |

## Executable Naming Convention

Executables are placed in `exe/` with the pattern:

```
<app>_<WaveEq>_<Solver>_<Order>_<TAG>
```

For example: `fwi_AC_FDSG_O4_release_avx2`

A symlink with the short name (e.g. `FWI`) always points to the most recently
built executable.

## Running on a Cluster

### 1. Get an interactive node

```bash
# CPU node
srun --pty -N 1 -n 16 -p cpu_partition --time=01:00:00 bash

# GPU node
srun --pty -N 1 -n 4 --gres=gpu:1 -p gpu_partition --time=01:00:00 bash
```

### 2. Load environment modules (example)

```bash
# CPU (Intel)
module load intel-oneapi-compilers intel-oneapi-mpi

# GPU (NVIDIA)
module load nvhpc openmpi
```

### 3. Set runtime environment

```bash
ulimit -s unlimited
export OMP_NUM_THREADS=4          # adjust to your allocation
```

### 4. Run

```bash
# CPU -- MPI ranks x OMP threads should equal allocated cores
mpirun -np 4 ./exe/FWI setup.in

# GPU -- typically 1 MPI rank per GPU
mpirun -np 1 ./exe/FWI setup.in
```

A quick manual for any executable is printed by running it without arguments:

```bash
./exe/FWD
./exe/FWI
```

## Benchmark Notes

- **CPU (ifort, AVX2)**: Best throughput when MPI ranks x OMP threads fills all
  available cores. Use `release_avx2` on AMD EPYC or Intel Haswell+.
- **GPU (nvfortran, managed memory)**: Start with `release_gpu` (managed
  memory) for correctness, then switch to `release_gpu_explicit` for best
  performance once validated. One MPI rank per GPU is the typical setup.
- For large 3D problems the GPU build can provide significant speedups over
  CPU, especially for the acoustic (`AC`) propagator.

## Cleaning

```bash
make clean       # remove .o files from FWD/ FWI/ RWI/ RTM/
make cleanmod    # remove compiled modules (.mod, .smod) and module .o files
make cleanall    # both of the above
```

## References

If you publish results using this code, please cite:

- Wei Zhou and David Lumley (2021), [Central-difference time-lapse 4D seismic full waveform inversion](https://library.seg.org/doi/10.1190/geo2019-0834.1), *Geophysics*
- Wei Zhou and David Lumley (2021), [Non-repeatability Effects on Time-Lapse 4D Seismic Full Waveform Inversion for Ocean-Bottom Node Data](https://library.seg.org/doi/10.1190/geo2020-0577.1), *Geophysics*

## License

GNU General Public License v3.0 -- see `LICENSE`.
