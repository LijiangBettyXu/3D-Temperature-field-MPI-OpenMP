# Performance Comparison of 3D Temperature Field Calculation Using MPI and OpenMP

## Introduction:

The calculation of 3D temperature fields is an essential task in various scientific and engineering applications, such as fluid dynamics, heat transfer, and climate modeling. Solving these problems often requires significant computational resources due to the large-scale nature of the systems and the necessity of capturing fine spatial and temporal details. To tackle these challenges, parallel computing techniques such as Message Passing Interface (MPI) and Open Multi-Processing (OpenMP) have been employed to enable efficient utilization of multi-core processors and distributed memory systems.

In this project, the performance of a 3D temperature field calculation have been investigated using both MPI and OpenMP. The objective is to compare the performance, scalability, and ease of implementation of these two popular parallel programming models. I have developed a numerical model that solves the heat equation on a discretized 3D domain using a finite difference method. The initial temperature field is generated with a localized heat source, and the heat diffusion process is simulated over time. The model incorporates Neumann boundary conditions to maintain a no-flux condition at the boundaries of the domain.

Results are analyzed in terms of execution time, speed-up, and efficiency, taking into consideration Amdahl's law and other factors that affect parallel performance. The project also discusses the trade-offs between MPI and OpenMP in terms of code complexity, maintainability, and portability. This comprehensive comparison between MPI and OpenMP in the context of a 3D temperature field calculation provides valuable insights for practitioners and researchers seeking to optimize their parallel computing strategies for similar scientific and engineering problems.

### Parameters:<br>
Side Length (L): 5m <br>
Discretization (N): 21, 41, 81, 161, 321 <br>
Grid size(∆x = ∆y = ∆z): L /(N - 1) <br>
Thermal Conductivity coefficient (𝛼): 0.001 W/mK <br>

### Initial condition:<br>
Sphere at center with radius = 0.1m<br>
Sphere high temperature = 100 Celcius<br>

### Boundary condition: <br>
No-flux Neuman boundary condition at the boundary of the cube.<br>

Heat Equation and the differential format using Forward Euler method:<br>

<img width="980" alt="image" src="https://user-images.githubusercontent.com/122394634/233512286-d13ef1fe-785d-4268-87bd-2e0d8bced2d3.png">

### Impelementation: <br>

The OpenMP file: <br>
TempAvg.cpp <br>
g++ -std=c++11 -fopenmp TempAvg.cpp -o avg <br>
export OMP_NUM_THREADS=1 <br>
./avg ${N} ${timesteps} <br>
The submit file: submitjobOpenmp.sb <br>

The MPI file: TempMPI.cpp <br>
mpic++ -std=c++11 TempMPI.cpp -o mpi <br>
mpirun -n ${num_process} ./mpi ${N} ${timesteps} <br>
The submit file: submitjobmpi.sb <br>

For the MPI implementation, we have employed a domain decomposition strategy to distribute the workload among different processes, each running on separate processing elements. Since only temperature of previous step is used in the calculation of the temperature in the next time (Forward Euler Method), there is no need to do MPI's point-to-point communication. In contrast, the OpenMP implementation leverages shared memory parallelism, with parallel for loops and work-sharing constructs to distribute the workload among available threads. 

We have conducted a series of experiments to evaluate the performance of both implementations for different problem sizes and varying numbers of processing elements or threads. The reauslt is averaged for multiple trails for stability. 

## Results: <br>

<img width="735" alt="image" src="https://user-images.githubusercontent.com/122394634/233512371-ba20512d-7e4f-4c78-9257-b6e6a43f20ea.png">

### For the time used in OpenMP: 

As I increase the number of parallel threads from 1, 2, 4, 8, 16, 32, 64, to 128, the execution time decreases. This decrease is more pronounced as the discretization number N increases. For N = 21, the decrease is relatively small until 64 threads, after which the time suddenly increases. For N = 41, the decrease is more noticeable until 64 threads, then there is a slightly increases afterward. For N = 81, 161, and 321, the time decrease becomes more significant, and there is almost no increase after using 64 threads. 

1) I think the speedup I observe when using more threads is more significant for larger N values because the workload is larger, leading to better utilization of the available CPU resources. The cache locality and hardware limitations play a more significant role in limiting the speedup for smaller N values. 
<pre>
Workload distribution: When the discretization number N is small, the workload distributed among the threads is also small. As you increase the number of threads, each thread has less work to do, which may not be enough to keep all threads busy, resulting in a less significant decrease in time. 

Cache locality: The overhead of synchronizing and managing the threads might outweigh the benefits of parallelism when N is small. For larger N values, the cache locality decreases, and the benefits of parallelism become more apparent.

Amdahl's Law: As N becomes larger, the effect of the sequential portion of the code (F) on the total runtime diminishes, allowing the parallel portion of the code (1 - F) to dominate the runtime. As a result, the speed-up becomes more obvious with larger discretization numbers N.

S(P) = 1 / (F + (1 - F) / P)
where F is the fraction of the program that is inherently sequential (i.e., cannot be parallelized) and (1 - F) is the fraction of the program that can be parallelized.
</pre>
2) The sudden time increase when using more than 64 threads when N is small
I think this phenomenon might due to hardware limitations (Maybe HPCC can only provide 64 threads when I requesting 128). Every CPU has a maximum number of physical and logical cores. When the number of threads exceeds the number of available cores, the operating system starts sharing CPU resources among threads, leading to context switching overhead. This could explain the sudden increase in time for smaller N values when using more than 64 threads.

### For the time used in MPI:

In the case of MPI code, I also tested with thread counts from 1, 2, 4, 8, 16, 32, 64, to 128. The execution time consistently decreases for N = 21 and 41. However, for N = 81, 161, and 321, the time slightly increases after using 32 threads, which differs from the OpenMP scenario. Additionally, the execution time for MPI is always shorter than that of OpenMP for N = 21 and 41, and the time difference continues to grow as the number of threads increases. For N = 81, the time difference becomes much smaller, and the execution times converge at 128 threads. The discrepancy further narrows for N = 161 and N = 321. In fact, there seems to be virtually no difference when N = 321. Interestingly, for N = 161, the MPI code's execution time even surpasses the OpenMP code's time at around 50 threads.

MPI uses multiple processes with separate memory spaces for parallelism, while OpenMP uses threads within a single process with shared memory space for parallelism. The nature of these divisions and their impact on performance may differ. For larger problem sizes (e.g., N = 81, 161, and 321), the discrepancy between MPI and OpenMP execution times becomes smaller. This may indicate that the workload is distributed more evenly among the threads in OpenMP, allowing for better utilization of the available resources.As N increases, the MPI implementation might encounter challenges that slow down its performance compared to the OpenMP implementation.

The reasons might be:<br>
Load imbalance: As N increases, the load imbalance between different processes in the MPI implementation might become more pronounced. This could lead to some processes taking much longer to complete their tasks, thereby slowing down the overall execution time.

Better parallelism in OpenMP: The OpenMP implementation may be better suited to handle the increased problem size as it can better utilize the available resources, leading to improved performance.

Increased overhead in MPI: As the problem size grows, the overhead associated with managing processes in MPI might increase, leading to slower execution times. In contrast, the overhead in OpenMP might not be as significant, allowing for faster execution times as N increases.

Speed up: <br>
With T1 the execution time on a single processor and Tp the time on p processors, we define the speedup as Sp = T1/Tp.

Efficiency: <br>
Ep = Sp/p. (0 < Ep ≤ 1). To measure how far we are from the ideal speedup.

<img width="1376" alt="image" src="https://user-images.githubusercontent.com/122394634/233512872-2cd45971-b3f2-41dd-9ce8-3e53a1e59d07.png">

<img width="1383" alt="image" src="https://user-images.githubusercontent.com/122394634/233512571-d8b92787-d5cd-4b7a-aa31-1efca9f4101a.png">



