#ifndef SPMV_H
#define SPMV_H

#include <stdbool.h>

bool generate_matrix(int rows, int cols, int percent_nonzero, int **I, int **J, double **vals, int *nz);
bool coo_to_csr(int nz, int start_row, int M, int *I, int *J, double *vals, int **row_ptr);

#endif