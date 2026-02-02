#ifndef SPMV_H
#define SPMV_H

double compute_avg(int num_iterations, double *time);
int find_outlier(int num_iterations, double *time, double avg);
void remove_outlier(int num_iterations, double *time, double *avg_time);
int find_max_M(int *row_dist, int processes);

#endif