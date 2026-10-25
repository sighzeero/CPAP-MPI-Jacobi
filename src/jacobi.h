#ifndef JACOBI_H
#define JACOBI_H
#include <stdbool.h>
#include "matrix.h"
#define precision 0.00001
#ifndef CONV_CHECK_INTERVAL
#define CONV_CHECK_INTERVAL 10
#endif
#ifndef MAX_ITERATIONS
#define MAX_ITERATIONS 10000
#endif
typedef struct {
    double *x;
    int k;
    double e;
    bool converged;
} jacobi_result;
jacobi_result *jacobi_parallel_mpi(matrix_local *m, bool verbose, int rank, int count);
#endif
