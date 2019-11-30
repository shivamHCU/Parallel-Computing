## Parallel-Computing
This repository contains assignment for Parallel Computing using MPI, OpenMP and CUDA.
----- Assignmet for Parallel Computing using MPI, OpenMP and CUDA -----

There are total 7 programs in this repo based on
This is just a informal discription.

## [1] MPI 
1.1) MPI_bsort (Parallel Bucket Sorting Algorithm and Implement using MPI (Static Strategy).
```
mpicc -o mpi_bsort MPI_bsort.c 
mpirun -np <no of processors> mpi_bsort <array_size> <max_limit>
```
1.2) MPI_search (Parallel Searching Algorithm and Implement using MPI)
```
mpicc -o mpi_search MPI_search.c 
mpirun -np <no of processors> mpi_search
``` 
1.3) MPI_intigration (Parallel Algorithm for computing the numerical integral of an arbitrary function)
```
mpicc -o mpi_intigration MPI_intigration.c 
mpirun -np <no of processors> mpi_intigration <start> <end> <intervals>
``` 
	
## [2] openMP
2.1) OMP_prefix_sum (Parallel Prefix Summation Algorithm and Implement using OpenMP)
```
gcc -o omp_prefix_sum -fopenmp OMP_prefix_sum.c
./omp_prefix_sum
```
2.2) OMP_prime (Parallel prime number generation Algorithm and Implement using OpenMP)
```
gcc -o omp_prime -fopenmp OMP_prime.c
./omp_prime
```
2.3) OMP_monte_carlo_pi (Parallel pi calculating Algorithm and Implement using OpenMP)
```
gcc -o omp_monte_carlo_pi -fopenmp OMP_monte_carlo_pi.c
./omp_monte_carlo_pi
```
## [3] CUDA

3.1) CUDA_matrix_mul (Parallel Matrix Multiplication Algorithm and Implement using CUDA)
```
nvcc -o cude_matrix_mul CUDA_matrix_mul.cu
./cude_matrix_mul
```
