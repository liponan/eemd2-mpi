eemd2-mpi
=========

MPI version of two-dimensional Ensemble Empirical Mode Decomposition
by Po-Nan Li in Institute of Physics, Academia Sinica, Taiwan

# Compilation example
Require GSL (GNU Scientific Library) and MPI
Intel compiler:

```
mpicxx mpi_eemd2d.2.2.cpp -lgsl -lgslcblas -lm -O3 -o mpi-eemd2
```

# Running command example

```
mpirun -np 24 ./mpi-eemd2 lena.txt 3 100 0.1
```
## Arguments
- First argument: input file in CSV style with first number indicating the number of dimension (i.e. 3), followed by three integers indicating the size in three dimensions. Remaining values are the elements of the input data in 1-D array style.
- Second argument: numbers of modes to decompose. (Default: 3)
- Third argument: Number of ensembles. (Default: 1)
- Fourth argument: Sigma of the white noise. If 0 (default), the EEMD will degenerate to EMD.
- Fifth argument: (not shown) file name suffix for the output data. Can be any integers in 0~9999. By default a random number will be assigned. 

# File usage
- **mpi_eemd2d.2.2.cpp**: main function
- **eemd.cpp**: 1-D EEMD (ensemble empirical mode decomposition) 
- **emd_core.cpp**: 1-D EMD (empirical mode decomposition)
- **find_extrema.cpp**: subfunction for emd_core.cpp
- **spline_gsl.cpp**: subfunction for emd_core.cpp. Powered by GSL.
- **print2bin.cpp**: subfunction for the main function.
- **README.md**: this file.

## Scripts for post-EEMD data processing
- **bin2m.m**: convert the binaray data to Matlab array.

