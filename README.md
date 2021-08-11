# Parallel Ant Colony Optimization Using High Performance Computing with CUDA and MPI (2021)

Metaheuristics approximate solutions of NP-hard combinatorial optimization problems (COP) that require significant computing resources due to the exponential growth of the search space. A parallel algorithm provides an efficient approach to solve large COPs and leverage parallel high-performance computing (HPC) environments. The Ant Colony Optimization (ACO) is a population-based metaheuristic with an inherently parallel nature with outstanding time and performance results. 

## Running the code (MPI)
- Prerequisites 
  1. MPI library
  2. g++

- Compile <br>
`
mpicc -lstdc++ -lm -lpthread parallel_acop.cpp
`
- Run <br>
`
mpirun -n 4 ./a.out rl1889
`
## Running the code (CUDA)
- Prerequisites 
  1. Nvidia CUDA Toolkit
- Usage (CUDA) <br>
`
./aco-cuda
`

## Data
The instances folder contain several TSP benchmark instances from the TSPLIB benchmark library:
- rat575
- rat783
- pr1002
- rl1889
- fl3795)


### If you find this code useful, please consider citing:

  ```
  @article{ 
  author  = {Jorge Banda-Almida, Israel Pineda},
  title   = {Parallel Ant Colony Optimization Using High Performance Computing with CUDA and MPI},
  year    = {2021},
  }
