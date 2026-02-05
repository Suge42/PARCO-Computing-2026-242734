#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include "SpMV.h"

void SpMV_csr(int M, int *row_ptr, int *col_idx, double *vals, double *vector, double *result) {
    for (int i = 0; i < M; i++) { // Loop over local rows
        double sum = 0.0;
        // Go through each non-zero element in the row
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) {
            sum += vals[j] * vector[col_idx[j]];
        }
        result[i] = sum;
    }
}

bool check_results(double *result_1, double *result_2, int M) {
    double epsilon = 1e-6; // Tolerance for floating-point comparison
    for (int i = 0; i < M; i++) {
        if (fabs(result_1[i] - result_2[i]) > epsilon) {
            printf("Mismatch at index %d: result_1=%f, result_2=%f\n", i, result_1[i], result_2[i]);
            return false;
        }
    }
    return true;
}