#ifndef MATRIX_MPI_H
#define MATRIX_MPI_H
/* CSR locale: diagonale separata e invertita, colonne con indici globali.
 * row_ptr e' ribasato a zero; intervallo di righe [row_start,row_end].
 * proc_row_start contiene i confini globali di tutti i rank (P+1 elementi).
 * La matrice viene letta in parallelo per intervalli di byte e redistribuita
 * per righe con Alltoallv. Tutti i rank devono accedere allo stesso file.
 */
typedef struct {
    unsigned int size;
    int row_start, row_end;
    int *row_ptr, *col_idx;
    double *values, *inv_diag, *b;
    int *proc_row_start;
    int process_count;
} matrix_local;
matrix_local *matrix_load_distributed(char *filename, int rank, int count);
void matrix_local_destroy(matrix_local *m);
#endif
