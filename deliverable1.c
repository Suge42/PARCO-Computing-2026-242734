/* 
*   Matrix Market I/O example program
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   and copies it to stdout.  This porgram does nothing useful, but
*   illustrates common usage of the Matrix Matrix I/O routines.
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*   Usage:  a.out [filename] > output
*
*       
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing 
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mmio.h"
#include "bubblesort.c"
#include <omp.h>

#define REPETITIONS 1

void fill_vector(int *vec, int size) {
    srand(time(NULL));
    //srand(42); // Fixed seed for testing
    for (int i = 0; i < size; i++) {
        vec[i] = (rand() % 9) +1; // Random integers between 1 and 9
    }
}

void sequential_moltiplication(int *row_ptr, double *values, int *vec, double *result, int M, int nz) {
    // Go through each row (row_ptr has M+1 elements)
    for (int i = 0; i < M; i++) {
        // For each row, go through its non-zero elements
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) {
            //printf("Row %d, accessing value index %d\n", i, j);
            result[i] += values[j] * vec[i];
        }
    }
}

void parallel_moltiplication(int *row_ptr, double *values, int *vec, double *result, int M, int nz) {
    // This time we parallelize the inner loop since is the most time consuming
    for (int i = 0; i < M; i++) {
        #pragma omp parallel for reduction(+:result[i])
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) {
            //printf("Row %d, accessing value index %d\n", i, j);
            result[i] += values[j] * vec[i];
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
    double *val;
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
    val = (double *) malloc(nz * sizeof(double)); // Values pointer


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    for (i=0; i<nz; i++) {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) fclose(f);

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
    bubbleSort(I, J, val, nz);
    
    //Conversion from COO to CSR
    int index = 0;
    int current_row = 0;
    int *a_row = (int *) malloc((M+1) * sizeof(int));
    a_row[0] = 0;
    a_row[1] = 0; // This is enough to initialie the array;

    while (index < nz && current_row < M) {
        if (I[index] == current_row) {
            // If I have an element in the current row, increment the row pointer
            a_row[current_row+1]++;

            // I only increase the index if I have found an element in the new row;
            // Otherwise, I keep the same index to check with the next row value;
            index++;
        } else {
            // If the element is not in the current row, I move to the next
            // and copy the starting value from the previous one;
            current_row++;
            a_row[current_row+1] = a_row[current_row];
        }
    }

    // Print matrix in CSR format
    /*printf("CSR Row pointer:\n");
    for (i=0; i<=M; i++) {
        printf("Row %d: %d\n", i, a_row[i]);
    }*/


    // Create random vector
    int *vec = (int *) malloc((M) * sizeof(int));
    fill_vector(vec, M);
    // Print vector
    /*for (i=0; i<M; i++) {
        printf("Vec[%d]: %d\n", i, vec[i]);
    }*/

    for (int r = 0; r < REPETITIONS; r++) {
        // Matrix-vector multiplication
        double *sequential_result = (double *) malloc((M) * sizeof(double));
        double *parallel_result = (double *) malloc((M) * sizeof(double));
        for (i = 0; i < M; i++) {
            sequential_result[i] = 0.0;
            parallel_result[i] = 0.0;
        }

        clock_t start, end;
        double seq_cpu_time_used, par_cpu_time_used;

        start = omp_get_wtime();
        // Sequential SpMV
        sequential_moltiplication(a_row, val, vec, sequential_result, M, nz);
        end = omp_get_wtime();
        seq_cpu_time_used = end - start;

        start = omp_get_wtime();
        // Parallel SpMV
        parallel_moltiplication(a_row, val, vec, parallel_result, M, nz);
        end = omp_get_wtime();
        par_cpu_time_used = end - start;

        printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
        printf("Sequential execution time for execution %d: %f seconds\n", r, seq_cpu_time_used);
        printf("Parallel execution time for execution %d: %f seconds\n", r, par_cpu_time_used);
        double speedup = seq_cpu_time_used / par_cpu_time_used * 100.0;
        printf("Speedup for execution %d: %.2f\n", r, speedup);

        // Check results
        int match = 0;
        int mismatch = 0;
        for (i = 0; i < M; i++) {
            if (sequential_result[i] != parallel_result[i]) {
                mismatch++;
            } else {
                match++;
            }
        }
        double precision = (double)match / (double)(match + mismatch) * 100.0;
        printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
        printf("Precision of parallel result for execution %d: %.2f%% (%d matches, %d mismatches)\n", r, precision, match, mismatch);
        printf("The precision error is due to the non-associativity of floating point addition.\n");
        //printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
        printf("\n");

        free(sequential_result);
        free(parallel_result);

        /*for (i = 0; i < M; i++) {
            printf("Sequential result[%d]: %f\n", i, sequential_result[i]);
        }
        for (i = 0; i < M; i++) {
            printf("Parallel result[%d]: %f\n", i, parallel_result[i]);
        }*/

        fflush(stdout);
    }







	return 0;
}
