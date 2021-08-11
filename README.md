# Parallel Ant Colony Optimization Using High Performance Computing with CUDA and MPI (2021)

# Running the code
## Prerequisites (MPI)
1. MPI library
2. g++
## Prerequisites (CUDA)
1. Nvidia CUDA Toolkit

## Ussage (MPI)
### Compile
`
mpicc -lstdc++ -lm -lpthread parallel_acop.cpp
`
### Run
`
mpirun -n 4 ./a.out rl1889
`
## Usage (CUDA)
`
./aco-cuda
`

### If you find this code useful, please consider citing:

  ```
  @article{ 
  author  = {Jorge Banda-Almida, Israel Pineda},
  title   = {Parallel Ant Colony Optimization Using High Performance Computing with CUDA and MPI},
  year    = {2021},
  }
