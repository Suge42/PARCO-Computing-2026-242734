#include <stdbool.h>
#include <stdio.h>

void swap(int *xp, int *yp){
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void swap_double(double *xp, double *yp){
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

// An optimized version of Bubble Sort
void bubbleSort(int *arr1, int *arr2, double *arr3, int n){
    int i, j;
    bool swapped;
    for (i = 0; i < n - 1; i++) {
        swapped = false;
        for (j = 0; j < n - i - 1; j++) {
            if (arr1[j] > arr1[j + 1]) {
                swap(&arr1[j], &arr1[j + 1]);
                swap(&arr2[j], &arr2[j + 1]);
                swap_double(&arr3[j], &arr3[j + 1]);
                swapped = true;
            }
        }

        // If no two elements were swapped by inner loop,
        // then break
        if (swapped == false)
            break;
    }
}
