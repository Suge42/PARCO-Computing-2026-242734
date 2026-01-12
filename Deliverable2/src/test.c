#include <stdio.h>
#include <stdlib.h>
#include "libraries/mmio.h"
#include "libraries/mmio.c"

int main(int argc, char *argv[]) {
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M; // Number of rows
    int N; // Number of columns
    int nz; // Total number of non-zero entries
    int *I, *J;
    double *vals;

    // Check the right amount of argument and open the file
    if (argc != 2) {
        fprintf(stderr, "Intended usage: %s [matrix-market-filename]\n", argv[0]);
        return 1;
    } else { 
        if ((f = fopen(argv[1], "r")) == NULL) {
            printf("Could not open file: %s\n", argv[1]);
            return 1;
        }
    }

    // printing f value
    long offset = ftell(f);
    printf("File pointer before reading banner: %ld\n", offset);
    


    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        return 1;
    }

    offset = ftell(f);
    printf("File pointer after reading banner: %ld\n", offset);

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0) {
        exit(1);
    }

    offset = ftell(f);
    printf("File pointer after reading matrix size: %ld\n", offset);


    I = (int *) malloc(10 * sizeof(int)); // Rows pointer
        J = (int *) malloc(10 * sizeof(int)); // Columns pointer
        vals = (double *) malloc(10 * sizeof(double)); // Values pointer
    

    
    long old_offset = offset;
    for (int i=0; i<10; i++) {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &vals[i]);
        offset = ftell(f);
        long distance = offset - old_offset;
        old_offset = offset;
        printf("File pointer after reading %d entry: %ld - Distance = %ld\n", i, offset, distance);
    }


    printf("Matrix Market banner read successfully.\n");
    fclose(f);
    return 0;
}