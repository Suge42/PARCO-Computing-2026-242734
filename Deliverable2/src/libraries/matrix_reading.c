#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix_reading.h"
#include "mmio.h"
#include "bubblesort.h"

bool check_matrix_file(char *filename, int *M, int *N, int *nz) {
    FILE *f;
    MM_typecode matcode;
    int ret_code;

    if ((f = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Could not open file: %s\n", filename);
        fflush(stderr);
        return false;
    }

    if (mm_read_banner(f, &matcode) != 0) {
        fprintf(stderr, "Could not process Matrix Market banner.\n");
        fflush(stderr);
        return false;
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)){
        fprintf(stderr, "Sorry, this application does not support ");
        fprintf(stderr, "Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        fflush(stderr);
        return false;
    }

    if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) !=0) {
        fprintf(stderr, "Error reading matrix size.\n");
        fflush(stderr);
        return false;
    }

    if (f != stdin) {
        fclose(f);
    }

    return true;
}


bool read_matrix_to_csr_total(char *filename, int **row_ptr, double **vals) {
    FILE *f;
    MM_typecode matcode;
    int ret_code;
    int M; // Number of rows
    int N; // Number of columns
    int nz; // Total number of non-zero entries
    int *I = NULL, *J = NULL;

    /* Simpler checks repeat, to ensure the file is correct */
    if ((f = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Could not open file: %s\n", filename);
        fflush(stderr);
        return false;
    }

    if (mm_read_banner(f, &matcode) != 0) {
        fprintf(stderr, "Could not process Matrix Market banner.\n");
        fflush(stderr);
        return false;
    }

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)){
        fprintf(stderr, "Sorry, this application does not support complex matrices.\n");
        fflush(stderr);
        return false;
    }


    /* find out size of sparse matrix .... */
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0) {
        fprintf(stderr, "Error reading matrix size.\n");
        fflush(stderr);
        return false;
    }

    /* reseve memory for matrices */
    I = (int *) malloc(nz * sizeof(int)); // Rows pointer
    J = (int *) malloc(nz * sizeof(int)); // Columns pointer
    *vals = (double *) malloc(nz * sizeof(double)); // Values pointer
    if (!I || !J || !(*vals)) {
        fprintf(stderr, "Failed to allocate memory for matrix data.\n");
        fflush(stderr);
        return false;
    }
    

    /* Reading the actual matrix data */
    for (int i=0; i<nz; i++) {
        fscanf(f, "%d %d %lf\n", &I[i], &J[i], &(*vals)[i]);
        I[i]--;  // adjust from 1-based to 0-based
        J[i]--;
    }

    if (f != stdin) {
        fclose(f);
    }

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
    bubbleSort(I, J, *vals, nz);
    free(J); // We don't need J anymore, we only work on rows and values
    
    //Conversion from COO to CSR
    int index = 0;
    *row_ptr = (int *) malloc((M+1) * sizeof(int));
    if (!(*row_ptr)) {
        fprintf(stderr, "Allocation failed. Needed ~%zu MB for COO format alone.\n",
                ((M+1 * sizeof(int))) / (1024 * 1024));
        fflush(stderr);
        return false;
    }
    (*row_ptr)[0] = 0; // This is enough to initialie the array;

    for (int i = 0; i < M; i++) {
        // Each row ends where it starts, unless we find elements
        (*row_ptr)[i+1] = (*row_ptr)[i]; 

        while (index < nz && I[index] == i) {
            (*row_ptr)[i+1]++;
            index++;
        }
    }

    // Print matrix in CSR format
    /*printf("CSR Row pointer:\n");
    for (i=0; i<=M; i++) {
        printf("Row %d: %d\n", i, a_row[i]);
    }*/

    free(I);

    return true;
}


bool read_matrix_to_csr_partial(char *filename, int start_row, int end_row, int *local_nz, int **row_ptr, double **vals) {
    // Similar implementation as matrix_to_csr_total but only for rows in [start_row, end_row)

    FILE *f;
    int M; // Number of rows
    int N; // Number of columns
    int nz; // Total number of non-zero entries
    int *local_I = NULL, *local_J = NULL;
    int local_M = end_row - start_row;

    /* Simpler checks repeat, to ensure the file is correct */
    if ((f = fopen(filename, "r")) == NULL) {
        return false;
    }

    /* Skip the header lines */
    char line[256];
    while (fgets(line, sizeof(line), f) != NULL) {
        // Just skip to the end of the header
        if (line[0] == '%') {
            continue;
        }
        // Found the first data line
        if (sscanf(line, "%d %d %d", &M, &N, &nz) == 3) {
            break;
        }
    }

    // Save file position after header
    long data_start_pos = ftell(f);

    /* Count elements in the specified row range */
    *local_nz = 0;
    
    for (int i=0; i<nz; i++) {
        int row_tmp, col_tmp;
        double value;
        fscanf(f, "%d %d %lf\n", &row_tmp, &col_tmp, &value);
        row_tmp--;  // adjust from 1-based to 0-based
        if (row_tmp >= start_row && row_tmp < end_row) {
            (*local_nz)++;
        }
        // Since the file is in column-major order, we can't stop early
    }


    /* reseve memory for matrices */
    local_I = (int *) malloc((*local_nz) * sizeof(int)); // Rows pointer
    local_J = (int *) malloc((*local_nz) * sizeof(int)); // Columns pointer
    *vals = (double *) malloc((*local_nz) * sizeof(double)); // Values pointer
    if (!local_I || !local_J || !(*vals)) {
        fprintf(stderr, "Failed to allocate memory for local matrix data.\n");
        fflush(stderr);
        return false;
    }


    /* Reading the actual matrix data */
    fseek(f, data_start_pos, SEEK_SET); // Reset file position to start reading data
    int index = 0;
    while (index < *local_nz && !feof(f)) {
        int row_tmp, col_tmp;
        double value;
        fscanf(f, "%d %d %lf\n", &row_tmp, &col_tmp, &value);
        row_tmp--;  // adjust from 1-based to 0-based
        col_tmp--;
        // Save only if in the specified row range
        if (row_tmp >= start_row && row_tmp < end_row) {
            local_I[index] = row_tmp - start_row; // Adjust row index for local CSR
            local_J[index] = col_tmp;
            (*vals)[index] = value;
            index++;
        }
        // Since the file is in column-major order, we can't stop early
    }

    if (f != stdin) {
        fclose(f);
    }

    // Sort by row indices
    bubbleSort(local_I, local_J, *vals, *local_nz);
    free(local_J); // We don't need J anymore, we only work on rows and values


    //Conversion from COO to CSR
    index = 0;
    *row_ptr = (int *) malloc((local_M+1) * sizeof(int));
    if (!(*row_ptr)) {
        fprintf(stderr, "Allocation failed. Needed ~%zu MB for COO format alone.\n",
                ((local_M+1 * sizeof(int))) / (1024 * 1024));
        fflush(stderr);
        return false;
    }
    (*row_ptr)[0] = 0; // This is enough to initialie the array;

    for (int i = 0; i < local_M; i++) {
        // Each row ends where it starts, unless we find elements
        (*row_ptr)[i+1] = (*row_ptr)[i]; 

        while (index < *local_nz && local_I[index] == i) {
            (*row_ptr)[i+1]++;
            index++;
        }
    }

    free(local_I);

    return true;
}
