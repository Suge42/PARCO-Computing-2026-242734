#ifndef SPMV_H
#define SPMV_H

#include <stdbool.h>

void SpMV_csr(int M, int *row_ptr, int *col_idx, double *vals, double *vector, double *result);
bool check_results(double *result_1, double *result_2, int M);

#endif