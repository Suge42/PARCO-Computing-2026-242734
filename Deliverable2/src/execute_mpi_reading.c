#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "libraries/SpMV.h"
#include "libraries/data_management.h"
#include "libraries/matrix_reading.h"
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, size, processes, num_iterations;
    MPI_Status status;
    char filename[256] = "";
    char result_filename[256] = "";

    int *row_ptr = NULL; // Initialize pointers to avoid problems with free()
    double *vals = NULL, *vector = NULL, *results = NULL;
    int M; // Number of rows
    int N; // Number of columns
    int nz; // Total number of non-zero entries
    //srand(42); // For debugging purposes
    srand(time(NULL));

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
    if (argc != 4) {
        if (rank == 0) {
            fprintf(stderr, "Intended usage: %s [martix-market-filename] [number-of-threads]\n", argv[0]);
            fflush(stderr);
        }
        MPI_Finalize();
        exit(1);
    }

    snprintf(result_filename, sizeof(result_filename), "%s", argv[3]); // File to store results

    num_iterations = atoi(argv[2]); // Number of times to repeat the sending process for averaging
    
    /* Allocating memory */
    double t_start, t_end;
    double *computation_time = malloc(num_iterations * sizeof(double));
    double *communication_time = malloc(num_iterations * sizeof(double));
    double *not_par_computation_time = malloc(num_iterations * sizeof(double));
    if (!computation_time || !communication_time || !not_par_computation_time) {
        fprintf(stderr, "Process %d failed to allocate memory for timing arrays\n", rank);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i = 0; i < num_iterations; i++) {
        computation_time[i] = 0.0;
        communication_time[i] = 0.0;
        not_par_computation_time[i] = 0.0;
    }
    

    for (int iter = 0; iter < num_iterations; iter++) {
        if (rank == 0) {
                    
            snprintf(filename, sizeof(filename), "%s", argv[1]); // Copy the filename to a local variable

            /* Initial checks on the matrix */
            printf("Iteration: %d - Process %d is checking the matrix: %s\n", iter+1, rank, filename);
            fflush(stdout);
            if (!check_matrix_file(filename, &M, &N, &nz)) {
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            /* Give the filename to other processes */
            t_start = MPI_Wtime();
            printf("Iteration: %d - Process %d is broadcasting the filename to other processes.\n", iter+1, rank);
            fflush(stdout);
            MPI_Bcast(&filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
            t_end = MPI_Wtime();
            communication_time[iter] += (t_end - t_start);
            

            /* Compute rows range for each process */
            int rows_per_process = M / processes;
            int remaining_rows = M % processes;

            int *rows_distribution = (int *) malloc(size * sizeof(int));
            if (!rows_distribution) {
                fprintf(stderr, "Iteration: %d - Process %d failed to allocate memory for rows distribution\n", iter+1, rank);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
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
            printf("Iteration: %d - Process %d is sending rows distribution to other processes.\n", iter+1, rank);
            fflush(stdout);
            for (int i = 0; i < processes; i++) {
                MPI_Send(&rows_distribution[i], 1, MPI_INT, i+1, 0, MPI_COMM_WORLD); // Start row
                MPI_Send(&rows_distribution[i+1], 1, MPI_INT, i+1, 0, MPI_COMM_WORLD); // End row
            }
            t_end = MPI_Wtime();
            communication_time[iter] += (t_end - t_start);
            
            
            /* Create vector of size M */
            vector = (double *) malloc(M * sizeof(double));
            if (!vector) {
                fprintf(stderr, "Iteration: %d - Process %d failed to allocate memory for random vector\n", iter+1, rank);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            for (int i = 0; i < M; i++) {
                vector[i] = (rand() % 9) + 1; // Initialize all elements to 1.0
            }


            /* Send the right part of the vector to all processes */
            t_start = MPI_Wtime();
            printf("Iteration: %d - Process %d is sending parts of the vector to other processes.\n", iter+1, rank);
            fflush(stdout);
            for (int i = 0; i < processes; i++) {
                int start_row = rows_distribution[i];
                int local_M = rows_distribution[i+1] - start_row;
                MPI_Send(&vector[start_row], local_M, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD);
            }
            t_end = MPI_Wtime();
            communication_time[iter] += (t_end - t_start);


            MPI_Barrier(MPI_COMM_WORLD);
            printf("Iteration: %d - Computation started.\n", iter+1);
            fflush(stdout);

            /* Allocate memory for results */
            results = (double *) malloc(M * sizeof(double));
            double *local_results = (double *) malloc(M * sizeof(double)); // Max size needed for rank 0
            if (!local_results || !results) {
                fprintf(stderr, "Iteration: %d - Process %d failed to allocate memory for results\n", iter+1, rank);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            

            /* Read the matrix into CSR format */
            if (!read_matrix_to_csr_total(filename, &row_ptr, &vals)) {
                fprintf(stderr, "Process 0 failed reading the whole matrix: %s\n", filename);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            /* Print matrix rows and values */
            /*printf("Process 0 CSR Row pointer:\n");
            fflush(stdout);
            for (int i=0; i<nz; i++) {
                printf("Val %d: %f\n", i, vals[i]);
            }
            fflush(stdout);*/


            /* Compute the SpMV result */
            double local_start = MPI_Wtime();
            SpMV_csr(M, row_ptr, vals, vector, local_results);
            double local_end = MPI_Wtime();
            not_par_computation_time[iter] = local_end - local_start;


            /* Receive back results from all processes */
            t_start = MPI_Wtime();
            
            int max_M = find_max_M(rows_distribution, processes);
            double *temp_buffer = (double *) malloc(max_M * sizeof(double));
            if (!temp_buffer) {
                fprintf(stderr, "Iteration: %d - Process %d failed to allocate memory for temp buffer while receiving back results\n", iter+1, rank);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            for (int i = 0; i < processes; i++) {
                int start_row = rows_distribution[i];
                int local_M = rows_distribution[i+1] - start_row;
                if (local_M > max_M) {
                    fprintf(stderr, "Iteration: %d - Process %d found local_M %d larger than max_M %d from process %d\n", iter+1, rank, local_M, max_M, i+1);
                    fflush(stderr);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                MPI_Recv(temp_buffer, local_M, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD, &status);
                for (int j = 0; j < local_M; j++) {
                    results[start_row + j] = temp_buffer[j];
                }
            }
            
            if (temp_buffer) {
                free(temp_buffer); // cleanup
                temp_buffer = NULL;
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
                printf("\tIteration: %d - Results are correct for MPI parallelization.\n", iter+1);
            } else {
                printf("\tIteration: %d - Results are NOT correct for MPI parallelization.\n", iter+1);
            }

            /* Receive computation time from processes */
            for (int i = 0; i < processes; i++) {
                double proc_comp_time;
                MPI_Recv(&proc_comp_time, 1, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD, &status);
                computation_time[iter] += proc_comp_time;
            }

            if (local_results) {
                free(local_results);
                local_results = NULL;
            }
            if (rows_distribution) {
                free(rows_distribution);
                rows_distribution = NULL;
            }

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
            if (!vector) {
                fprintf(stderr, "Process %d failed to allocate memory for vector\n", rank);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            MPI_Recv(vector, local_M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

            /* Print received vector */
            //printf("Process %d received vector:\n", rank);
            /*for (int i = 0; i < local_M; i++) {
                printf("Rank: %d - Vector[%d]: %f\n", rank, i+start_row, vector[i]);
            }
            fflush(stdout);*/

            // Wait for all processes to be ready, then start timing
            MPI_Barrier(MPI_COMM_WORLD);

            /* Receive result vector to fill */
            results = (double *) malloc(local_M * sizeof(double));
            if (!results) {
                fprintf(stderr, "Process %d failed to allocate memory for results vector\n", rank);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            

            if (!read_matrix_to_csr_partial(filename, start_row, end_row, &nz, &row_ptr, &vals)) {
                fprintf(stderr, "Process %d failed reading its part of the matrix: %s\n", rank, filename);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            //printf("Process %d read its part of the matrix with %d non-zero elements.\n", rank, local_nz);
            //fflush(stdout);
            
            /* Printf matrix rows and values */
            //printf("Process %d CSR Row pointer:\n", rank);
            /*for (int i=0; i<nz; i++) {
                printf("Rank: %d - Val %d: %f\n", rank, i, vals[i]);
            }
            fflush(stdout);*/

            /* Compute the SpMV result */
            t_start = MPI_Wtime();
            //printf("Process %d is computing its SpMV part.\n", rank);
            //fflush(stdout);
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

        // Barrier to ensure all processes finished using heap memory before freeing
        MPI_Barrier(MPI_COMM_WORLD);

        if (row_ptr) {
            free(row_ptr);
            row_ptr = NULL;
        }
        if (vals) {
            free(vals);
            vals = NULL;
        }
        if (results) {
            free(results);
            results = NULL;
        }
        if (vector) {
            free(vector);
            vector = NULL;
        }
        // Barrier to synchronize before next iteration
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* Print overall times and speedup */
    if (rank == 0) {
        double comp_time = 0.0;
        double comm_time = 0.0;
        double not_par_comp_time = 0.0;

        remove_outlier(num_iterations, computation_time, &comp_time);
        remove_outlier(num_iterations, communication_time, &comm_time);
        remove_outlier(num_iterations, not_par_computation_time, &not_par_comp_time);
        
        double total_time = comp_time + comm_time;
        
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
        printf("=-=\n\n");
        fflush(stdout);

        /* Write results to file */
        FILE *f;
        if ((f = fopen(result_filename, "a")) == NULL) {
            fprintf(stderr, "Could not open file: %s\n", result_filename);
            fflush(stderr);
            return false;
        }
        fprintf(f, "#Matrix: %s - Row: %d - Columns: %d - Working_processes: %d\n", filename, M, N, processes);
        fprintf(f, "avg_comp_time: %f\n", avg_comp_time);
        fprintf(f, "avg_comm_time: %f\n", avg_comm_time);
        fprintf(f, "avg_total_time: %f\n", avg_total_time);
        fprintf(f, "not_par_comp_time: %f\n", not_par_comp_time);
        fprintf(f, "speedup: %f\n", speedup);
        fflush(f);
        fclose(f);
    }
    

    
    free(computation_time);
    free(communication_time);
    free(not_par_computation_time);
    

    MPI_Finalize();
    return 0;
}
