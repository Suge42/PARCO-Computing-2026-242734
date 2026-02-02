#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "generator.h"
#include "bubblesort.h"

bool generate_matrix(int rows, int cols, int percent_nonzero, int **I, int **J, double **vals, int *nz) {
    int total = rows * cols;
    uint8_t *mask = malloc(total); // Keep track of filled positions
    if (!mask) {
        fprintf(stderr, "Mask allocation failed\n");
        return false;
    }
    int local_nz = 0;

    /* Count non-zero elements (no need to have the real matrix in memory) */
    for (int k = 0; k < total; ++k) {
        if ((rand() % 100) < percent_nonzero) {
            mask[k] = 1; // 1 if non-zero 
            local_nz++; // Increment count if non-zero
        } else {
            mask[k] = 0; // 0 otherwise
        }
        
    }

    *nz = local_nz;

    *I = malloc(local_nz * sizeof(int));
    *J = malloc(local_nz * sizeof(int));
    *vals = malloc(local_nz * sizeof(double));
    if (!*I || !*J || !*vals) {
        fprintf(stderr, "COO allocation failed\n");
        return false;
    }

    /* Generate COO representation */
    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int k = i * cols + j;
            if (mask[k]) {
                // Decide the value sign randomly
                double sign = ((rand() & 1) ? -1.0 : 1.0);
                (*I)[index] = i;
                (*J)[index] = j;
                // Obtain values between -10.0 and 10.0
                (*vals)[index] = sign * ((rand() % 1000) / 100.0);
                index++;
            }
        }
    }

    return true;
}


bool coo_to_csr(int nz, int start_row, int M, int *I, int *J, double *vals, int **row_ptr) {
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
