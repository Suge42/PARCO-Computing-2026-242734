#ifndef MATRIX_READING_H
#define MATRIX_READING_H

bool check_matrix_file(char *filename, int *M, int *N, int *nz);
bool read_matrix_to_csr_total(char *filename, int **row_ptr, double **vals);
bool read_matrix_to_csr_partial(char *filename, int start_row, int end_row, int *local_nz, int **row_ptr, double **vals);

#endif