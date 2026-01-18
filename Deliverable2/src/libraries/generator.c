#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "bubblesort.c"

int generate_matrix(int rows, int cols, int percent_nonzero, int **I, int **J, double **vals, int *nz) {
    if (rows <= 0 || cols <= 0 || percent_nonzero < 0 || percent_nonzero > 100) {
        fprintf(stderr, "Invalid parameters for matrix generation.\n");
        fflush(stderr);
        return false;
    }
    

    int total_size = rows * cols;
    double *matrix = malloc(total_size * sizeof(double));
    if (!matrix) {
        fprintf(stderr, "Allocation failed. Needed ~%zu MB for matrix alone.\n",
                total_size * sizeof(double) / (1024 * 1024));
        fflush(stderr);
        return false;
    }
    

    // Fill matrix
    int local_nz = 0;
    for (int i = 0; i < total_size; ++i) {
        if ((rand() % 100) < percent_nonzero) {
            double sign = 1.0;
            if ((rand() % 100) < 50) {
                sign = -1.0;
            }
            // Random value between -10 and 10
            matrix[i] = sign * ((rand() % 1000) / 100.0);
            local_nz++;
        } else {
            matrix[i] = 0.0;
        }
    }

    *nz = local_nz;
    /* Convert matrix to COO format */
    *I = (int *) malloc((local_nz) * sizeof(int));
    *J = (int *) malloc((local_nz) * sizeof(int));
    *vals = (double *) malloc((local_nz) * sizeof(double));
    if (!(*I) || !(*J) || !(*vals)) {
        fprintf(stderr, "Allocation failed. Needed ~%zu MB for COO format alone.\n",
                ((local_nz * sizeof(double))+(local_nz * sizeof(int) * 2)) / (1024 * 1024));
        fflush(stderr);
        free(matrix);
        return false;
    }
    

    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double value = matrix[(i * cols) + j];
            if (value != 0.0) {
                (*I)[index] = i;
                (*J)[index] = j;
                (*vals)[index] = value;
                index++;
            }
        }
    }
    
    free(matrix);
    return true;
}

int coo_to_csr(int nz, int start_row, int M, int *I, int *J, double *vals, int **row_ptr) {
    // Sort by row indices
    bubbleSort(I, J, vals, nz);
    
    // Conversion from COO to CSR
    int index = 0;
    int current_row = 0;
    (*row_ptr) = (int *) malloc((M+1) * sizeof(int));
    if (!(*row_ptr)) {
        fprintf(stderr, "Allocation failed. Needed ~%zu MB for COO format alone.\n",
                ((M+1 * sizeof(int))) / (1024 * 1024));
        fflush(stderr);
        return false;
    }
    (*row_ptr)[0] = 0; 
    (*row_ptr)[1] = 0; // This is enough to initialie the array;

    while (current_row < M) {
        if (index < nz) {
            if (I[index] == current_row + start_row) {
                // If I have an element in the current row, increment the row pointer
                (*row_ptr)[current_row+1]++;

                // I only increase the index if I have found an element in the new row;
                // Otherwise, I keep the same index to check with the next row value;
                index++;
            } else {
                // If the element is not in the current row, I move to the next
                // and copy the starting value from the previous one;
                current_row++;
                (*row_ptr)[current_row+1] = (*row_ptr)[current_row];
            }
        } else {
            // If I finish the elements because the last rows are empty, I still
            // Need to fill them
            current_row++;
            (*row_ptr)[current_row+1] = (*row_ptr)[current_row];
        }
    }

    return true;
}