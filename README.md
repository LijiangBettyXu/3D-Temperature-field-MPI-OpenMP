# Performance Comparison of 3D Temperature Field Calculation Using MPI and OpenMP

## Introduction:

The calculation of 3D temperature fields is an essential task in various scientific and engineering applications, such as fluid dynamics, heat transfer, and climate modeling. Solving these problems often requires significant computational resources due to the large-scale nature of the systems and the necessity of capturing fine spatial and temporal details. To tackle these challenges, parallel computing techniques such as Message Passing Interface (MPI) and Open Multi-Processing (OpenMP) have been employed to enable efficient utilization of multi-core processors and distributed memory systems.

In this project, the performance of a 3D temperature field calculation have been investigated using both MPI and OpenMP. The objective is to compare the performance, scalability, and ease of implementation of these two popular parallel programming models. I have developed a numerical model that solves the heat equation on a discretized 3D domain using a finite difference method. The initial temperature field is generated with a localized heat source, and the heat diffusion process is simulated over time. The model incorporates Neumann boundary conditions to maintain a no-flux condition at the boundaries of the domain.

### Parameters:<br>
Side Length (L): 5 <br>
Discretization (N): 21, 41, 81, 161, 321 <br>
grid size(∆x = ∆y = ∆z): L /(N - 1) <br>

### Initial condition:<br>
Sphere at center with radius = 0.1<br>
Sphere high temperature = 100<br>

### Boundary condition: <br>
No-flux Neuman boundary condition at the boundary of the cube.<br>

Heat Equation and the differential format using Forward Euler method:<br>

<img width="781" alt="image" src="https://user-images.githubusercontent.com/122394634/233174175-b3efcce0-ea73-4250-9f44-023ad39db3b9.png">

### Impelementation: <br>

The OpenMP file: TempAvg.cpp <br>
g++ -std=c++11 -fopenmp TempAvg.cpp -o avg <br>
export OMP_NUM_THREADS=1 <br>
./avg ${N} ${timesteps} <br>

The MPI file: TempMPI.cpp <br>
mpic++ -std=c++11 TempMPI.cpp -o mpi <br>
mpirun -n ${num_process} ./mpi ${N} ${timesteps} <br>

For the MPI implementation, we have employed a domain decomposition strategy to distribute the workload among different processes, each running on separate processing elements. Since only temperature of previous step is used in the calculation of the temperature in the next time (Forward Euler Method), there is no need to do MPI's point-to-point communication. In contrast, the OpenMP implementation leverages shared memory parallelism, with parallel for loops and work-sharing constructs to distribute the workload among available threads. 

We have conducted a series of experiments to evaluate the performance of both implementations for different problem sizes and varying numbers of processing elements or threads. The reauslt is averaged for multiple trails for stability. 




Results are analyzed in terms of execution time, speed-up, and efficiency, taking into consideration Amdahl's law and other factors that affect parallel performance. The project also discusses the trade-offs between MPI and OpenMP in terms of code complexity, maintainability, and portability.

This comprehensive comparison between MPI and OpenMP in the context of a 3D temperature field calculation provides valuable insights for practitioners and researchers seeking to optimize their parallel computing strategies for similar scientific and engineering problems.





Regenerate response
