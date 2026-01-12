#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "libraries/matrix_reading.c"
#include "libraries/SpMV.c"
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, size, processes, num_iterations;
    MPI_Status status;
    char filename[256] = "";

    int *row_ptr;
    double *vals, *vector, *results;
    srand(42);
    //srand(time(NULL));

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    processes = size-1; // Rank 0 is not included


    if (size < 2) {
        if (rank == 0) {
            fprintf(stderr, "Error: run with at least 2 processes.\n");
            fflush(stderr);
        }
        MPI_Finalize();
        exit(1);
    }

    /* Check the right amount of argument and open the file */
    if (argc != 3) { // PROBABLY NEEDS TO BE CHANGED
        if (rank == 0) {
            fprintf(stderr, "Intended usage: %s [martix-market-filename] [number-of-threads]\n", argv[0]);
            fflush(stderr);
        }
        MPI_Finalize();
        exit(1);
    }

    num_iterations = atoi(argv[2]); // Number of times to repeat the sending process for averaging
    double t_start, t_end;
    double *computation_time = malloc(num_iterations * sizeof(double));
    double *communication_time = malloc(num_iterations * sizeof(double));
    double *not_par_computation_time = malloc(num_iterations * sizeof(double));
    for (int i = 0; i < num_iterations; i++) {
        computation_time[i] = 0.0;
        communication_time[i] = 0.0;
        not_par_computation_time[i] = 0.0;
    }

    for (int iter = 0; iter < num_iterations; iter++) {
        if (rank == 0) {
            int M; // Number of rows
            int N; // Number of columns
            int nz; // Total number of non-zero entries
                    
            strncpy(filename, argv[1], 256); // Copy the filename to a local variable

            /* Initial checks on the matrix */
            printf("Process %d is checking the matrix: %s\n", rank, filename);
            fflush(stdout);
            if (!check_matrix(filename, &M, &N, &nz)) {
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            /* Give the filename to other processes */
            t_start = MPI_Wtime();
            printf("Process %d is broadcasting the filename to other processes.\n", rank);
            fflush(stdout);
            MPI_Bcast(&filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
            t_end = MPI_Wtime();
            communication_time[iter] += (t_end - t_start);
            

            /* Compute rows range for each process */
            int rows_per_process = M / processes;
            int remaining_rows = M % processes;

            int *rows_distribution = (int *) malloc(size * sizeof(int));
            rows_distribution[0] = 0; // Rank 0 does not process rows

            for (int i = 0; i < processes; i++) {
                rows_distribution[i+1] = rows_distribution[i] + rows_per_process; // Similar to CSR format
                if (i < remaining_rows) {
                    rows_distribution[i+1]++;
                }
            }
            
            // print row distribution for debugging
            /*for (int i = 0; i < processes; i++) {
                printf("Process %d: rows %d to %d\n", i+1, rows_distribution[i], rows_distribution[i+1]-1);
            }*/
            

            /* Send the rows distribution to all processes */
            t_start = MPI_Wtime();
            printf("Process %d is sending rows distribution to other processes.\n", rank);
            fflush(stdout);
            for (int i = 0; i < processes; i++) {
                MPI_Send(&rows_distribution[i], 1, MPI_INT, i+1, 0, MPI_COMM_WORLD); // Start row
                MPI_Send(&rows_distribution[i+1], 1, MPI_INT, i+1, 0, MPI_COMM_WORLD); // End row
            }
            t_end = MPI_Wtime();
            communication_time[iter] += (t_end - t_start);
            
            
            /* Create vector of size M */
            vector = (double *) malloc(M * sizeof(double));
            for (int i = 0; i < M; i++) {
                vector[i] = (rand() % 9) +1; // Initialize all elements to 1.0
            }


            /* Send the right part of the vector to all processes */
            t_start = MPI_Wtime();
            printf("Process %d is sending parts of the vector to other processes.\n", rank);
            fflush(stdout);
            for (int i = 0; i < processes; i++) {
                int start_row = rows_distribution[i];
                int local_M = rows_distribution[i+1] - start_row;
                MPI_Send(&vector[start_row], local_M, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD);
            }
            t_end = MPI_Wtime();
            communication_time[iter] += (t_end - t_start);


            MPI_Barrier(MPI_COMM_WORLD);
            double local_start = MPI_Wtime();
            printf("Computation started.\n");
            fflush(stdout);

            /* Allocate memory for results */
            results = (double *) malloc(M * sizeof(double));
            double *local_results = (double *) malloc(M * sizeof(double)); // Max size needed for rank 0
            

            /* Allocate memory */
            row_ptr = (int *) malloc((M+1) * sizeof(int));
            vals = (double *) malloc(nz * sizeof(double));

            /* Read the matrix into CSR format */
            if (!matrix_to_csr_total(filename, &row_ptr, &vals)) {
                fprintf(stderr, "Process 0 failed reading the whole matrix: %s\n", filename);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            /* Printf matrix rows and values */
            /*printf("Process 0 CSR Row pointer:\n");
            fflush(stdout);
            for (int i=0; i<nz; i++) {
                printf("Val %d: %f\n", i, vals[i]);
            }
            fflush(stdout);*/


            /* Compute the SpMV result */
            SpMV_csr(M, row_ptr, vals, vector, local_results);
            double local_end = MPI_Wtime();
            not_par_computation_time[iter] = local_end - local_start;

            /* Receive back results from all processes */
            t_start = MPI_Wtime();
            for (int i = 0; i < processes; i++) {
                int start_row = rows_distribution[i];
                int local_M = rows_distribution[i+1] - start_row;
                double *temp_buffer = (double *) malloc(local_M * sizeof(double));
                MPI_Recv(temp_buffer, local_M, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD, &status);
                for (int j = 0; j < local_M; j++) {
                    results[start_row + j] = temp_buffer[j];
                }
            }
            t_end = MPI_Wtime();
            communication_time[iter] += (t_end - t_start);

            /* Print result vector */
            /*printf("Process %d results:\n", rank);
            for (int i = 0; i < M; i++) {
                printf("Rank: %d - Local results[%d]: %f\n", rank, i, local_results[i]);
            }
            fflush(stdout);*/

            /* Print result vector */
            /*printf("Process %d results:\n", rank);
            for (int i = 0; i < M; i++) {
                printf("Rank: %d - Results[%d]: %f\n", rank, i, results[i]);
            }
            fflush(stdout);*/


            // Barrier to ensure all processes finish before checking results
            MPI_Barrier(MPI_COMM_WORLD);

            /* Verify correctness */
            //printf("Process %d is checking results correctness:\n", rank);
            //fflush(stdout);
            if (check_results(local_results, results, M)) {
                printf("\tResults are correct for MPI parallelization.\n");
            } else {
                printf("\tResults are NOT correct for MPI parallelization.\n");
            }

            /* Receive computation time from processes */
            for (int i = 0; i < processes; i++) {
                double proc_comp_time;
                MPI_Recv(&proc_comp_time, 1, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD, &status);
                computation_time[iter] += proc_comp_time;
            }

            free(local_results);
            free(rows_distribution);

        } else {        
            /* Receive the filename from rank 0 */
            MPI_Bcast(&filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
            

            /* Receive the rows distribution from rank 0 */
            int start_row, end_row, local_M;
            MPI_Recv(&start_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&end_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            local_M = end_row - start_row;
            /*printf("Process %d received rows %d to %d.\n", rank, start_row, end_row-1);
            printf("Process %d local_M: %d\n", rank, local_M);
            fflush(stdout);*/

            /* Receive the part of the vector from rank 0 */
            vector = (double *) malloc(local_M * sizeof(double));
            MPI_Recv(vector, local_M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

            /* Print received vector */
            //printf("Process %d received vector:\n", rank);
            /*for (int i = 0; i < local_M; i++) {
                printf("Rank: %d - Vector[%d]: %f\n", rank, i+start_row, vector[i]);
            }
            fflush(stdout);*/

            // Wait for all processes to be ready, then start timing
            MPI_Barrier(MPI_COMM_WORLD);
            t_start = MPI_Wtime();

            /* Receive result vector to fill */
            results = (double *) malloc(local_M * sizeof(double));
            //MPI_Recv(results, local_M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

            int local_nz;
            if (!matrix_to_csr_partial(filename, start_row, end_row, &local_nz, &row_ptr, &vals)) {
                fprintf(stderr, "Process %d failed reading its part of the matrix: %s\n", rank, filename);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            //printf("Process %d read its part of the matrix with %d non-zero elements.\n", rank, local_nz);
            //fflush(stdout);
            
            /* Printf matrix rows and values */
            //printf("Process %d CSR Row pointer:\n", rank);
            /*for (int i=0; i<local_nz; i++) {
                printf("Rank: %d - Val %d: %f\n", rank, i, vals[i]);
            }
            fflush(stdout);*/

            /* Compute the SpMV result */
            printf("Process %d is computing its SpMV part.\n", rank);
            fflush(stdout);
            SpMV_csr(local_M, row_ptr, vals, vector, results);

            /* Print result vector */
            //printf("Process %d results:\n", rank);
            /*for (int i = 0; i < local_M; i++) {
                printf("Rank: %d - Results[%d]: %f\n", rank, i, results[i]);
            }
            fflush(stdout);*/

            t_end = MPI_Wtime();
            double local_comp_time = t_end - t_start;

            /* Send back results to rank 0 */
            MPI_Send(results, local_M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

            /* Send computation time to rank 0 */
            MPI_Send(&local_comp_time, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

            // Barrier to ensure all processes finish before checking results
            MPI_Barrier(MPI_COMM_WORLD);
        }

        free(row_ptr);
        free(vals);
        free(vector);
        free(results);
        // Barrier to synchronize before next iteration
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* Print overall times and speedup */
    if (rank == 0) {
        double comp_time = 0.0;
        double comm_time = 0.0;
        double total_time = 0.0;
        double not_par_comp_time = 0.0;
        for (int iter = 0; iter < num_iterations; iter++) {
            comp_time += computation_time[iter];
            comm_time += communication_time[iter];
            total_time += (computation_time[iter] + communication_time[iter]);
            not_par_comp_time += not_par_computation_time[iter];
        }
        comp_time /= num_iterations; // Average over iterations
        comm_time /= num_iterations; 
        total_time /= num_iterations; 
        not_par_comp_time /= num_iterations; 
        double avg_comp_time = comp_time / processes; // Average per process
        double avg_comm_time = comm_time / processes;
        double avg_total_time = total_time / processes;

        double speedup = not_par_comp_time / avg_comp_time;

        printf("=-=\n");
        printf("Overall computation time across processes: %f seconds.\n", comp_time);
        printf("Overall communication time across processes: %f seconds.\n", comm_time);
        printf("Overall total time across processes: %f seconds.\n", total_time);
        printf("=-=\n");
        printf("Average computation time across processes: %f seconds.\n", avg_comp_time);
        printf("Average communication time across processes: %f seconds.\n", avg_comm_time);
        printf("Average total time across processes: %f seconds.\n", avg_total_time);
        printf("=-=\n");
        printf("Unparallelized computation time: %f seconds.\n", not_par_comp_time);
        printf("Speedup achieved: %f\n", speedup);
        fflush(stdout);
    }
    

    free(computation_time);
    free(communication_time);
    free(not_par_computation_time);

    MPI_Finalize();
    return 0;
}
