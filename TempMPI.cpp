// mpic++ -std=c++11 TempMPI.cpp -o mpi
// mpirun -n ${num_process} ./mpi ${N} ${timesteps}

#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>

using tensor = std::vector<std::vector<std::vector<double>>>;

// Update Temperature field using Forward Euler method
double UpdateTemp(tensor &T, int N, int timesteps, double alpha, double dt, double dx, int rank, int num_process) {
    tensor T_new(N, std::vector<std::vector<double>>(N, std::vector<double>(N, 0)));

    // Divide the rows amongst the number of processors
    int local_rows = N / num_process;
    int remainder = N % num_process;
    int start_row = rank * local_rows + std::min(rank, remainder);
    int end_row = start_row + local_rows + (rank < remainder);

    double start_time = MPI_Wtime();

    // Main loop
    for (int t = 0; t < timesteps; ++t) {


        // Calculate the Laplacian using MPI
        for (int i = std::max(start_row, 1); i < std::min(end_row, N - 1); ++i) {
            for (int j = 1; j < N - 1; ++j) {
                for (int k = 1; k < N - 1; ++k) {
                    double laplacian =
                        (T[i + 1][j][k] - 2 * T[i][j][k] + T[i - 1][j][k]) / (dx * dx) +
                        (T[i][j + 1][k] - 2 * T[i][j][k] + T[i][j - 1][k]) / (dx * dx) +
                        (T[i][j][k + 1] - 2 * T[i][j][k] + T[i][j][k - 1]) / (dx * dx);

                    // Update the temperature using the Laplacian
                    T_new[i][j][k] = T[i][j][k] + alpha * dt * laplacian;
               }
            }
        }
    
        // Enforce no flux Neumann boundary conditions
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                T_new[0][j][k] = T_new[1][j][k];
                T_new[N - 1][j][k] = T_new[N - 2][j][k];
            }
        }

        for (int i = 1; i < N - 1; ++i) {
            for (int k = 0; k < N; ++k) {
                T_new[i][0][k] = T_new[i][1][k];
                T_new[i][N - 1][k] = T_new[i][N - 2][k];
            }
        }

        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                T_new[i][j][0] = T_new[i][j][1];
                T_new[i][j][N - 1] = T_new[i][j][N - 2];
            }
        }

        // Update the temperature field
        T = T_new;
    }

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    return elapsed_time;
}

// Test whether two temperature fields are the same or not
// return 0 when they are the same
double diff(const tensor &M1, const tensor &M2){
    size_t N = M1.size();
    // size_t W = (L > 0) ? M1[0].size() : 0;
    // size_t H = (L > 0 && W > 0) ? M1[0][0].size() : 0;

    double w2 = 0;
    for (int i = 0; i < N; ++i) 
        for (int j = 0; j < N; ++j) 
            for (int k = 0; k < N; ++k) 
                w2 += pow(M1[i][j][k] - M2[i][j][k], 2);
    return w2;
}

int main(int argc, char *argv[]) {

    int result = MPI_Init(&argc,&argv);
    if (result != MPI_SUCCESS) return 1; 

    int N = (argc > 1) ? atol(argv[1]) : 21;   // Discritization
    int timesteps = (argc > 2) ? atoi(argv[2]) : 100; // Iteration number

    int rank, num_process;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);

    // Parameters
    double L = 5.0;
    double dx = L / (N - 1);
    double alpha = 1e-3;
    double dt = dx * dx / (6 * alpha);
    double heat_source_radius = 0.1;
    double heat_source_temperature = 100;

    // Initialize temperature field
    tensor T(N, std::vector<std::vector<double>>(N, std::vector<double>(N, 0)));

    // Add Initial condition
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                double dist_from_center = std::sqrt(std::pow(i * dx - L / 2, 2) + std::pow(j * dx - L / 2, 2) + std::pow(k * dx - L / 2, 2));
                if (dist_from_center <= heat_source_radius){
                    T[i][j][k] = heat_source_temperature;
                }
            }
        }
    }

    tensor T1(T), T2(T);

    // Call UpdateTemp three times and compute the average elapsed time
    double elapsed_time1 = UpdateTemp(T, N, timesteps, alpha, dt, dx, rank, num_process);

    double elapsed_time2 = UpdateTemp(T1, N, timesteps, alpha, dt, dx, rank, num_process);
    
    double elapsed_time3 = UpdateTemp(T2, N, timesteps, alpha, dt, dx, rank, num_process);

    // Calculate the average elapsed time across all processes
    double local_average_elapsed_time = (elapsed_time1 + elapsed_time2 + elapsed_time3) / 3;
    double global_average_elapsed_time;
    MPI_Allreduce(&local_average_elapsed_time, &global_average_elapsed_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    global_average_elapsed_time /= num_process;

    if (rank == 0) {
        std::cout << "Average elapsed time: " << global_average_elapsed_time << " seconds" << std::endl;
        std::cout << "Difference: " << diff(T, T1) << " " << diff(T, T2) << std::endl;
    }

    MPI_Finalize();
    return 0;
}


