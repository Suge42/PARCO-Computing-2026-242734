#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "data_management.h"

double compute_avg(int num_iterations, double *time) {
    double avg = 0.0;

    for (int iter = 0; iter < num_iterations; iter++) {
        avg += time[iter];
    }
    avg = avg / num_iterations;

    return avg;
}

int find_outlier(int num_iterations, double *time, double avg) {
    int outlier_index = 0;
    double max_deviation = fabs(time[0] - avg); // By default is the first

    for (int iter = 1; iter < num_iterations; iter++) {
        double deviation = fabs(time[iter] - avg);
        if (deviation > max_deviation) { // Choose the one with the maximum deviation
            max_deviation = deviation;
            outlier_index = iter;
        }
    }

    return outlier_index;
}

void remove_outlier(int num_iterations, double *time, double *avg_time) {
    if (num_iterations <= 2) {
        // Not enough data to remove outlier
        *avg_time = compute_avg(num_iterations, time);
        return;
    }
    double avg = compute_avg(num_iterations, time);

    /* Find the index of the outlier */
    int outlier_index = find_outlier(num_iterations, time, avg);

    for (int iter = 0; iter < num_iterations; iter++) {
        if (iter != outlier_index) {
            *avg_time += time[iter];
        }
    }
    *avg_time = *avg_time / (num_iterations-1); // Average over iterations
}

int find_max_M(int *row_dist, int processes) {
    int max_M = row_dist[1] - row_dist[0]; // Initial assumption

    for (int i = 1; i < processes; i++) {
        int local_M = row_dist[i+1] - row_dist[i];
        if (local_M > max_M) {
            max_M = local_M;
        }
    }

    return max_M;
}