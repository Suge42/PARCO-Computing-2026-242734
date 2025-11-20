#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "libraries/mmio.h"
#include "libraries/mmio.c"
#include "libraries/bubblesort.c"
#include <omp.h>

#define REPETITIONS 10

void fill_vector(int *vec, int size) {
    srand(time(NULL));
    //srand(42); // Fixed seed for testing
    for (int i = 0; i < size; i++) {
        vec[i] = (rand() % 9) +1; // Random integers between 1 and 9
    }
}

void seq_molt(int *row_ptr, double *values, int *vec, double *result, int M, int nz) {
    // Go through each row (row_ptr has M+1 elements)
    for (int i = 0; i < M; i++) {
        // For each row, go through its non-zero elements
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) {
            //printf("Row %d, accessing value index %d\n", i, j);
            result[i] += values[j] * vec[i];
        }
    }
}

void csr_par_molt(int *row_ptr, double *values, int *vec, double *result, int M, int nz) {
    // This time we parallelize the inner loop since is the most time consuming
    for (int i = 0; i < M; i++) {
        #pragma omp parallel for reduction(+:result[i])
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) {
            //printf("Row %d, accessing value index %d\n", i, j);
            result[i] += values[j] * vec[i];
        }
    }
}

void csr2_par_molt(int *s_row, int num_sr, int *row_ptr, double *values, int *vec, double *result) {
    #pragma omp parallel for
    for(int i = 0; i < num_sr; i++) {
        int sr_start = s_row[i];
        int sr_end = s_row[i+1];

        for (int j = sr_start; j < sr_end; j++) {
            int row_start = row_ptr[j];
            int row_end = row_ptr[j+1];

            for(int k = row_start; k < row_end; k++) {
                result[j] += values[k] * vec[j];
            }
        }
    }
}

void csr3_par_molt(int *ss_row, int num_ssr, int *s_row, int *row_ptr, double *values, int *vec, double *result) {
    #pragma omp parallel for
    for(int i = 0; i < num_ssr; i++) {
        int ssr_start = ss_row[i];
        int ssr_end = ss_row[i+1];

        for (int j = ssr_start; j < ssr_end; j++) {
            int sr_start = s_row[j];
            int sr_end = s_row[j+1];

            for (int k = sr_start; k < sr_end; k++) {
                int row_start = row_ptr[k];
                int row_end = row_ptr[k+1];

                for(int l = row_start; l < row_end; l++) {
                    result[k] += values[l] * vec[k];
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M; // Number of rows
    int N; // Number of columns
    int nz; // Total number of non-zero entries
    int i, *I, *J;
    double *vals;
    int *ordered_rows;
    int *ordered_colums;
    double *ordered_val;

    // Check the right amount of argument and open the file
    if (argc != 3) {
		fprintf(stderr, "Intended usage: %s [martix-market-filename] [number-of-threads]\n", argv[0]);
		exit(1);
	} else { 
        if ((f = fopen(argv[1], "r")) == NULL) {
            printf("Could not open file: %s\n", argv[1]);
            exit(1);
        }
    }

    // Set number of threads
    omp_set_num_threads(atoi(argv[2]));

    
    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)){
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */
    I = (int *) malloc(nz * sizeof(int)); // Rows pointer
    J = (int *) malloc(nz * sizeof(int)); // Columns pointer
    vals = (double *) malloc(nz * sizeof(double)); // Values pointer


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    for (i=0; i<nz; i++) {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &vals[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin)
        fclose(f);

    // Print the matrix
    /*mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    for (i=0; i<nz; i++) {
        fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);
        if (I[i] > 5)
            return 0;
    }*/
    // rows - colums - value
    // COO ordered by columns


    // Sort by row indices
    bubbleSort(I, J, vals, nz);
    
    //Conversion from COO to CSR
    int index = 0;
    int current_row = 0;
    int *row_ptr = (int *) malloc((M+1) * sizeof(int));
    row_ptr[0] = 0;
    row_ptr[1] = 0; // This is enough to initialie the array;

    while (index < nz && current_row < M) {
        if (I[index] == current_row) {
            // If I have an element in the current row, increment the row pointer
            row_ptr[current_row+1]++;

            // I only increase the index if I have found an element in the new row;
            // Otherwise, I keep the same index to check with the next row value;
            index++;
        } else {
            // If the element is not in the current row, I move to the next
            // and copy the starting value from the previous one;
            current_row++;
            row_ptr[current_row+1] = row_ptr[current_row];
        }
    }

    // Print matrix in CSR format
    /*printf("CSR Row pointer:\n");
    for (i=0; i<=M; i++) {
        printf("Row %d: %d\n", i, a_row[i]);
    }*/


    // Creation of sr and ssr vectors for CSR-K

    // Evalueate size of sr vector
    // we choose SRS = 96 as suggested in the paper (each group of row will have at max 96 matrix elements)
    int sr_size = 96;
    int progress = 0;
    int num_sr = 1; // First element is always 0
    for (i = 0; i < M; i++) {
        // We count the elements of each rows to decide how to group them
        progress += row_ptr[i+1] - row_ptr[i];
        if (progress > sr_size) {
            num_sr++;
            //printf("progress exceeded at row %d with progress %d\n", i, progress);
            progress = row_ptr[i+1] - row_ptr[i]; // We start the new group with the current row
        }
    }
    //printf("SR vector length: %d\n", sr_length);

    int *s_row = (int *) malloc((num_sr+1) * sizeof(int));
    s_row[0] = 0;
    progress = 0;
    int j = 1;
    for (i = 0; i < M; i++) {
        progress += row_ptr[i+1] - row_ptr[i];
        //count++;
        if (progress > sr_size) {
            s_row[j] = i;
            j++;
            progress = row_ptr[i+1] - row_ptr[i];
        }
    }
    s_row[j] = M; // Last element is always M

    // Print sr vector
    /*printf("SR vector:\n");
    for (i = 0; i <= sr_length; i++) {
        printf("SR[%d]: %d\n", i, s_row[i]);
    }*/



    int ssr_size = 32; // Each sub-group will have 32 elements from row_ptr vector
    progress = 0;
    int num_ssr = 1; // First element is always 0
    for (i = 0; i < num_sr; i++) {
        progress += s_row[i+1] - s_row[i];
        if (progress > ssr_size) {
            num_ssr++;
            progress = s_row[i+1] - s_row[i];
        }
    }

    int *ss_row = (int *) malloc((num_ssr+1) * sizeof(int));
    ss_row[0] = 0;
    progress = 0;
    j = 1;
    for (i = 0; i < num_sr; i++) {
        progress += s_row[i+1] - s_row[i];
        //count++;
        if (progress > ssr_size) {
            ss_row[j] = i;
            j++;
            progress = s_row[i+1] - s_row[i];
        }
    }
    ss_row[j] = num_sr;

    /*printf("SSR vector:\n");
    for (i = 0; i <= num_ssr; i++) {
        printf("SSR[%d]: %d\n", i, ss_row[i]);
    }*/

    double start, end;
    double seq_cpu_time_used, csr_cpu_time_used, csr2_cpu_time_used, csr3_cpu_time_used;

    // Collect the speedup values to plot the graph later
    double *csr_speedup_values = (double *) malloc(REPETITIONS * sizeof(double));
    double *csr2_speedup_values = (double *) malloc(REPETITIONS * sizeof(double));
    double *csr3_speedup_values = (double *) malloc(REPETITIONS * sizeof(double));

    // NORMAL PARALLELIZATION TESTING
    for (int r = 0; r < REPETITIONS; r++) {

        // Create random vector
        int *rand_vec = (int *) malloc((M) * sizeof(int));
        fill_vector(rand_vec, M);
        // Print vector
        /*for (i=0; i<M; i++) {
            printf("Vec[%d]: %d\n", i, vec[i]);
        }*/
                
        double *seq_result = (double *) malloc((M) * sizeof(double));
        double *csr_par_result = (double *) malloc((M) * sizeof(double));
        double *csr2_par_result = (double *) malloc((M) * sizeof(double));
        double *csr3_par_result = (double *) malloc((M) * sizeof(double));
        for (i = 0; i < M; i++) {
            seq_result[i] = 0.0;
            csr_par_result[i] = 0.0;
            csr2_par_result[i] = 0.0;
        }

        start = omp_get_wtime() * 1000.0;
        // Sequential SpMV
        seq_molt(row_ptr, vals, rand_vec, seq_result, M, nz);
        end = omp_get_wtime() * 1000.0;
        seq_cpu_time_used = end - start;

        start = omp_get_wtime() * 1000.0;
        // Parallel SpMV
        csr_par_molt(row_ptr, vals, rand_vec, csr_par_result, M, nz);
        end = omp_get_wtime() * 1000.0;
        //printf("Start %f\n", start);
        //printf("End %f\n", end);
        csr_cpu_time_used = end - start;


        start = omp_get_wtime() * 1000.0;
        // CSR-2 Parallel SpMV
        csr2_par_molt(s_row, num_sr, row_ptr, vals, rand_vec, csr2_par_result);
        end = omp_get_wtime() * 1000.0;
        csr2_cpu_time_used = end - start;


        start = omp_get_wtime() * 1000.0;
        // CSR-3 Parallel SpMV
        csr3_par_molt(ss_row, num_ssr, s_row, row_ptr, vals, rand_vec, csr3_par_result);
        end = omp_get_wtime() * 1000.0;
        csr3_cpu_time_used = end - start;
        

        printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
        printf("Sequential execution time for execution %d: %f milliseconds\n", r+1, seq_cpu_time_used);
        printf("Parallel execution time for execution %d: %f milliseconds\n", r+1, csr_cpu_time_used);
        printf("CSR-2 Parallel execution time for execution %d: %f milliseconds\n", r+1, csr2_cpu_time_used);
        printf("CSR-3 Parallel execution time for execution %d: %f milliseconds\n", r+1, csr3_cpu_time_used);
        csr_speedup_values[r] = seq_cpu_time_used / csr_cpu_time_used * 100.0;
        csr2_speedup_values[r] = seq_cpu_time_used / csr2_cpu_time_used * 100.0;
        csr3_speedup_values[r] = seq_cpu_time_used / csr3_cpu_time_used * 100.0;
        printf("Speedup par for execution %d : %.2f%%\n", r+1, csr_speedup_values[r]);
        printf("Speedup CSR-2 for execution %d : %.2f%%\n", r+1, csr2_speedup_values[r]);
        printf("Speedup CSR-3 for execution %d : %.2f%%\n", r+1, csr3_speedup_values[r]);

        // Check results
        bool correct = true;
        for (i = 0; i<M; i++) {
            if ((seq_result[i] - csr_par_result[i]) > 0.00001 || (csr_par_result[i] - seq_result[i]) > 0.00001) {
                correct = false;
            }
        }
        if (correct) {
            printf("Results are correct for normal parallelization.\n");
        } else {
            printf("Results are NOT correct for normal parallelization.\n");
        }


        correct = true;
        for (i = 0; i<M; i++) {
            if ((seq_result[i] - csr2_par_result[i]) > 0.00001 || (csr2_par_result[i] - seq_result[i]) > 0.00001) {
                correct = false;
            }
        }
        if (correct) {
            printf("Results are correct for CSR-2 parallelization.\n");
        } else {
            printf("Results are NOT correct for CSR-2 parallelization.\n");
        }


        correct = true;
        for (i = 0; i<M; i++) {
            if ((seq_result[i] - csr3_par_result[i]) > 0.00001 || (csr3_par_result[i] - seq_result[i]) > 0.00001) {
                correct = false;
            }
        }
        if (correct) {
            printf("Results are correct for CSR-3 parallelization.\n");
        } else {
            printf("Results are NOT correct for CSR-3 parallelization.\n");
        }

        printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
        printf("\n");

        free(seq_result);
        free(csr_par_result);

        /*for (i = 0; i < M; i++) {
            printf("Sequential result[%d]: %f\n", i, sequential_result[i]);
        }
        for (i = 0; i < M; i++) {
            printf("Parallel result[%d]: %f\n", i, parallel_result[i]);
        }*/

        fflush(stdout);
    }

    //plot_graph(speedup_values, "plots/speedup_regular_parallelization.png");
    //plot_graph(speedup_values, "plot.png");







	return 0;
}
