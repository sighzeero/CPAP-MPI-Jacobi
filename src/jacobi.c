#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jacobi.h"
#include "timer.h"
#include "results.h"

int main(int argc, char **argv) {
    int rank, count;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count);
    if (argc != 5 || strcmp(argv[2], "M") != 0 ||
        (strcmp(argv[4], "0") != 0 && strcmp(argv[4], "1") != 0)) {
        if (rank == 0) fprintf(stderr, "Uso: %s matrice.mtx M processi verbose(0|1)\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    char *end;
    long declared = strtol(argv[3], &end, 10);
    if (*end || declared != count) {
        if (rank == 0) fprintf(stderr, "Il numero dichiarato deve coincidere con mpirun -np (%d).\n", count);
        MPI_Finalize();
        return 1;
    }
    int verbose = atoi(argv[4]);
    matrix_local *m = matrix_load_distributed(argv[1], rank, count);
    /* Il caricamento e' escluso dai tempi, setup halo e Gatherv inclusi. */
    MPI_Barrier(MPI_COMM_WORLD);
    timer *t = rank == 0 ? start_timer() : NULL;
    double begin = MPI_Wtime();
    jacobi_result *res = jacobi_parallel_mpi(m, verbose, rank, count);
    double elapsed = MPI_Wtime() - begin, maximum;
    MPI_Reduce(&elapsed, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    int success = 0;
    if (rank == 0) {
        stop_timer(t, false);
        /* Tempi di parete: massimo sui rank, in millisecondi. */
        t->timeval_diff = maximum * 1000.0;
        t->timespec_diff = maximum * 1000.0;
        success = res->converged;
        printf("processes=%d iterations=%d relative_update=%.12g converged=%d elapsed_ms=%.9f\n",
               count, res->k, res->e, success, maximum * 1000.0);
        if (verbose) {
            printf("Results:");
            for (unsigned int i = 0; i < m->size; ++i) printf(" %.17g", res->x[i]);
            printf("\n");
        }
        if (success) write_results(t, argv[1], count, 'M', m->size);
        else fprintf(stderr, "Nessuna convergenza: misura esclusa dai risultati.\n");
        free(res->x);
        free(res);
        free(t);
    }
    MPI_Bcast(&success, 1, MPI_INT, 0, MPI_COMM_WORLD);
    matrix_local_destroy(m);
    MPI_Finalize();
    return success ? 0 : 2;
}
