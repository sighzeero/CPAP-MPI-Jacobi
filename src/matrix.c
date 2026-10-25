#include <mpi.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include "matrix.h"

/* Lettura parallela Matrix Market, partizione contigua bilanciata per
 * nonzeri fuori diagonale, redistribuzione triplette e costruzione CSR. */

    // ricerca binaria del proprietario
    static int find_owner(int global_idx, int *proc_row_start, int process_count)
{
    int lo = 0, hi = process_count;
    while (lo < hi)
    {
        int mid = (lo + hi) / 2;
        if (proc_row_start[mid + 1] <= global_idx)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

// buffer dinamico per le triplette lette localmente (row, col, val)
typedef struct
{
    int *row;
    int *col;
    double *val;
    long count;
    long capacity;
} triplet_buf;

void matrix_local_destroy(matrix_local *m)
{
    free(m->row_ptr);
    free(m->col_idx);
    free(m->values);
    free(m->inv_diag);
    free(m->b);
    free(m->proc_row_start);
    free(m);
}

static void triplet_buf_init(triplet_buf *t, long initial_capacity)
{
    if (initial_capacity < 16)
        initial_capacity = 16;
    t->row = malloc(initial_capacity * sizeof(int));
    t->col = malloc(initial_capacity * sizeof(int));
    t->val = malloc(initial_capacity * sizeof(double));
    t->count = 0;
    t->capacity = initial_capacity;
}

static void triplet_buf_push(triplet_buf *t, int row, int col, double val)
{
    if (t->count == t->capacity)
    {
        t->capacity *= 2;
        t->row = realloc(t->row, t->capacity * sizeof(int));
        t->col = realloc(t->col, t->capacity * sizeof(int));
        t->val = realloc(t->val, t->capacity * sizeof(double));
    }
    t->row[t->count] = row;
    t->col[t->count] = col;
    t->val[t->count] = val;
    t->count++;
}

static void triplet_buf_free(triplet_buf *t)
{
    free(t->row);
    free(t->col);
    free(t->val);
}

/*
 * Parsa l'header MatrixMarket. Ogni rank la chiama in autonomia,
 * senza alcuna comunicazione:
 * Ritorna anche data_offset, il byte a cui inizia la prima riga di
 * dati (serve a tutti i rank per calcolare la propria fetta).
 */
static void parse_mm_header(FILE *file, long *size, long *count,
                            bool *is_symmetric, long *data_offset)
{
    char line[1024];
    long file_pos;
    *is_symmetric = false;

    do
    {
        file_pos = ftell(file);
        if (fgets(line, sizeof(line), file) == NULL)
        {
            puts("\nUnexpected end of file while skipping header");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (line[0] == '%')
        {
            if (strncmp(line, "%%MatrixMarket", 14) == 0)
            {
                if (strstr(line, "symmetric") != NULL ||
                    strstr(line, "hermitian") != NULL)
                {
                    *is_symmetric = true;
                }
            }
        }
    } while (line[0] == '%');

    fseek(file, file_pos, SEEK_SET);
    long size_cols;
    if (fscanf(file, "%li %li %li", size, &size_cols, count) != 3 ||
        *size <= 0 || *size > INT_MAX || *count < 0) {
        fprintf(stderr, "Dimensioni Matrix Market non valide\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (*size != size_cols)
    {
        puts("\nMatrix is not square");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // consuma il resto della riga delle dimensioni per arrivare
    // esattamente all'inizio della prima riga di dati
    int c;
    while ((c = fgetc(file)) != '\n' && c != EOF)
    {
    }
    *data_offset = ftell(file);
}

/*
 * Legge e parsa SOLO la fetta di byte [chunk_start, chunk_end) del
 * file, allineandosi alle righe. Riempie tb con le triplette lette
 * (mirror-ate se is_symmetric), local_row_count[] (dimensione `size`,
 * gia' allocato dal chiamante) con il numero di entry FUORI DIAGONALE
 * che finiranno in ciascuna riga globale.
 */
static void parse_chunk(const char *filename, long chunk_start, long chunk_end,
                        long data_offset, bool is_symmetric,
                        triplet_buf *tb, int *local_row_count)
{

    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        puts("\nCannot open file");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    fseek(file, chunk_start, SEEK_SET);

    // se non partiamo esattamente sull'inizio dei dati (rank 0),
    // la nostra fetta quasi certamente comincia a meta' di una riga
    // altrui: la scartiamo e ripartiamo dalla riga successiva.
    if (chunk_start > data_offset)
    {
        fseek(file, chunk_start - 1, SEEK_SET);
        int previous = fgetc(file);
        if (previous != '\n') {
            int c;
            while ((c = fgetc(file)) != '\n' && c != EOF) {}
        }
    }

    char line[1024];
    long pos;
    while ((pos = ftell(file)) < chunk_end)
    {
        if (fgets(line, sizeof(line), file) == NULL)
            break; // EOF

        long row, col;
        double val;
        if (sscanf(line, "%li %li %lf", &row, &col, &val) != 3)
            continue;
        row--;
        col--;

        if (row == col)
        {
            // entry diagonale: viaggia comunque come tripletta, la
            // smisteremo dopo l'Alltoallv guardando row==col
            triplet_buf_push(tb, (int)row, (int)col, val);
        }
        else
        {
            triplet_buf_push(tb, (int)row, (int)col, val);
            local_row_count[row]++;

            if (is_symmetric)
            {
                triplet_buf_push(tb, (int)col, (int)row, val);
                local_row_count[col]++;

            }
        }
    }

    fclose(file);
}

matrix_local *matrix_load_distributed(char *filename, int process_num, int process_count)
{

    //--------------------------------------------------------------
    // 1) Ogni rank parsa l'header DA SOLO: nessuna comunicazione.
    //--------------------------------------------------------------
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        puts("\nCannot open file");
        puts(filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    long size, count, data_offset;
    bool is_symmetric;
    parse_mm_header(file, &size, &count, &is_symmetric, &data_offset);
    fseek(file, 0, SEEK_END);
    long filesize = ftell(file);
    fclose(file);

    //--------------------------------------------------------------
    // 2) Divisione del corpo dati in fette di byte ~uguali.
    //--------------------------------------------------------------
    long data_size = filesize - data_offset;
    long chunk_start = data_offset + (long)(((double)data_size * process_num) / process_count);
    long chunk_end = (process_num == process_count - 1)
                         ? filesize
                         : data_offset + (long)(((double)data_size * (process_num + 1)) / process_count);

    //--------------------------------------------------------------
    // 3) Parsing locale della sola fetta di competenza.
    //--------------------------------------------------------------
    triplet_buf tb;
    long estimated = (count / process_count + 16) * (is_symmetric ? 2 : 1);
    triplet_buf_init(&tb, estimated);

    int *local_row_count = calloc(size, sizeof(int));
    parse_chunk(filename, chunk_start, chunk_end, data_offset, is_symmetric,
                &tb, local_row_count);

    //--------------------------------------------------------------
    // 4) UN SOLO Allreduce per ottenere lo stesso
    //    row_ptr globale su ogni rank -> partizione bilanciata
    //    ricalcolata identica su tutti, senza comunicazione extra.
    //--------------------------------------------------------------
    int *global_row_count = malloc(size * sizeof(int));
    MPI_Allreduce(local_row_count, global_row_count, (int)size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    free(local_row_count);

    int *global_row_ptr = malloc((size + 1) * sizeof(int));
    global_row_ptr[0] = 0;
    for (long i = 0; i < size; i++)
    {
        global_row_ptr[i + 1] = global_row_ptr[i] + global_row_count[i];
    }
    free(global_row_count);

    int *proc_row_start = malloc((process_count + 1) * sizeof(int));
    double weight_tot = global_row_ptr[size];
    proc_row_start[0] = 0;
    for (int t = 1; t < process_count; t++)
    {
        double target = (weight_tot * t) / process_count;
        int lo = proc_row_start[t - 1], hi = (int)size;
        while (lo < hi)
        {
            int mid = (lo + hi) / 2;
            if (global_row_ptr[mid] < target)
                lo = mid + 1;
            else
                hi = mid;
        }
        proc_row_start[t] = lo;
    }
    proc_row_start[process_count] = (int)size;

    int row_start = proc_row_start[process_num];
    int row_end = proc_row_start[process_num + 1] - 1;
    int local_rows = row_end - row_start + 1;

    //--------------------------------------------------------------
    // 5) Redistribuzione delle triplette al legittimo proprietario
    //    con un solo Alltoallv.
    //--------------------------------------------------------------
    int *send_counts = calloc(process_count, sizeof(int));
    int *owner_of = malloc((tb.count > 0 ? tb.count : 1) * sizeof(int));
    for (long i = 0; i < tb.count; i++)
    {
        int owner = find_owner(tb.row[i], proc_row_start, process_count);
        owner_of[i] = owner;
        send_counts[owner]++;
    }

    int *send_displs = malloc(process_count * sizeof(int));
    send_displs[0] = 0;
    for (int p = 1; p < process_count; p++)
        send_displs[p] = send_displs[p - 1] + send_counts[p - 1];

    long n_send = tb.count;
    int *send_row = malloc((n_send > 0 ? n_send : 1) * sizeof(int));
    int *send_col = malloc((n_send > 0 ? n_send : 1) * sizeof(int));
    double *send_val = malloc((n_send > 0 ? n_send : 1) * sizeof(double));

    int *cursor = malloc(process_count * sizeof(int));
    memcpy(cursor, send_displs, process_count * sizeof(int));
    for (long i = 0; i < tb.count; i++)
    {
        int p = owner_of[i];
        int pos = cursor[p]++;
        send_row[pos] = tb.row[i];
        send_col[pos] = tb.col[i];
        send_val[pos] = tb.val[i];
    }
    free(cursor);
    free(owner_of);
    triplet_buf_free(&tb);

    int *recv_counts = malloc(process_count * sizeof(int));
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

    int *recv_displs = malloc(process_count * sizeof(int));
    recv_displs[0] = 0;
    for (int p = 1; p < process_count; p++)
        recv_displs[p] = recv_displs[p - 1] + recv_counts[p - 1];
    int n_recv = recv_displs[process_count - 1] + recv_counts[process_count - 1];

    int *recv_row = malloc((n_recv > 0 ? n_recv : 1) * sizeof(int));
    int *recv_col = malloc((n_recv > 0 ? n_recv : 1) * sizeof(int));
    double *recv_val = malloc((n_recv > 0 ? n_recv : 1) * sizeof(double));

    MPI_Alltoallv(send_row, send_counts, send_displs, MPI_INT,
                  recv_row, recv_counts, recv_displs, MPI_INT, MPI_COMM_WORLD);
    MPI_Alltoallv(send_col, send_counts, send_displs, MPI_INT,
                  recv_col, recv_counts, recv_displs, MPI_INT, MPI_COMM_WORLD);
    MPI_Alltoallv(send_val, send_counts, send_displs, MPI_DOUBLE,
                  recv_val, recv_counts, recv_displs, MPI_DOUBLE, MPI_COMM_WORLD);

    free(send_row);
    free(send_col);
    free(send_val);
    free(send_counts);
    free(send_displs);
    free(recv_counts);
    free(recv_displs);

    //--------------------------------------------------------------
    // 6) Costruzione della CSR locale dalle triplette ricevute.
    //--------------------------------------------------------------
    int *local_row_ptr = malloc((local_rows + 1) * sizeof(int));
    int base = global_row_ptr[row_start];
    for (int i = 0; i <= local_rows; i++)
    {
        local_row_ptr[i] = global_row_ptr[row_start + i] - base;
    }
    int local_nnz = local_row_ptr[local_rows];
    free(global_row_ptr);

    int *local_col_idx = malloc((local_nnz > 0 ? local_nnz : 1) * sizeof(int));
    double *local_values = malloc((local_nnz > 0 ? local_nnz : 1) * sizeof(double));
    double *local_inv_diag = calloc(local_rows > 0 ? local_rows : 1, sizeof(double));

    int *fill_cursor = malloc((local_rows > 0 ? local_rows : 1) * sizeof(int));
    memcpy(fill_cursor, local_row_ptr, local_rows * sizeof(int));

    for (int i = 0; i < n_recv; i++)
    {
        int grow = recv_row[i];
        int lrow = grow - row_start;
        if (recv_row[i] == recv_col[i])
        {
            local_inv_diag[lrow] = recv_val[i]; // raw, invertito sotto
        }
        else
        {
            int pos = fill_cursor[lrow]++;
            local_col_idx[pos] = recv_col[i]; // resta indice GLOBALE di colonna
            local_values[pos] = recv_val[i];
        }
    }
    free(fill_cursor);
    free(recv_row);
    free(recv_col);
    free(recv_val);

    for (int i = 0; i < local_rows; i++)
    {
        if (local_inv_diag[i] == 0.0)
        {
            printf("[rank %d] Errore: diagonale mancante o nulla sulla riga globale %d\n",
                   process_num, row_start + i);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            local_inv_diag[i] = 1.0 / local_inv_diag[i];
        }
    }

    //--------------------------------------------------------------
    // Vettore b: se esiste il file _b.mtx, rank 0 lo legge e lo
    // distribuisce con Scatterv, usando la partizione
    // che ogni rank ha già' calcolato identicamente al passo 4.
    // Se il file non esiste, b = A*1 si calcola in LOCALE (ogni
    // rank ha gia' la propria CSR): niente comunicazione affatto.
    //--------------------------------------------------------------
    size_t len = strlen(filename);
    if (len < 4 || strcmp(filename + len - 4, ".mtx") != 0) {
        fprintf(stderr, "Input richiesto: file .mtx\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    char *filename_b = malloc(len + 3);
    memcpy(filename_b, filename, len - 4);
    strcpy(filename_b + len - 4, "_b.mtx");

    double *local_b = malloc((local_rows > 0 ? local_rows : 1) * sizeof(double));

    if (access(filename_b, R_OK) != -1)
    {
        double *global_b = NULL;
        if (process_num == 0)
        {
            FILE *fb = fopen(filename_b, "r");
            if (!fb) MPI_Abort(MPI_COMM_WORLD, 1);
            long bcount, bcols;
            if (fscanf(fb, "%li %li", &bcount, &bcols) != 2 || bcount != size || bcols != 1) {
                fprintf(stderr, "File b: attesi N 1 e N valori reali\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            global_b = malloc(size * sizeof(double));
            for (long i = 0; i < bcount; i++)
            {
                if (fscanf(fb, "%lf", &global_b[i]) != 1)
                    MPI_Abort(MPI_COMM_WORLD, 1);
            }
            fclose(fb);
        }
        int *row_counts = malloc(process_count * sizeof(int));
        int *row_disps = malloc(process_count * sizeof(int));
        for (int p = 0; p < process_count; p++)
        {
            row_counts[p] = proc_row_start[p + 1] - proc_row_start[p];
            row_disps[p] = proc_row_start[p];
        }
        MPI_Scatterv(global_b, row_counts, row_disps, MPI_DOUBLE,
                     local_b, local_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        free(row_counts);
        free(row_disps);
        if (process_num == 0)
            free(global_b);
    }
    else
    {
        // b = A*1: diagonale + somma degli elementi fuori diagonale
        // della riga, calcolato interamente in locale.
        for (int i = 0; i < local_rows; i++)
        {
            double diag_val = (local_inv_diag[i] != 0.0) ? (1.0 / local_inv_diag[i]) : 0.0;
            double sum = diag_val;
            for (int idx = local_row_ptr[i]; idx < local_row_ptr[i + 1]; idx++)
            {
                sum += local_values[idx];
            }
            local_b[i] = sum;
        }
    }

    free(filename_b);
    matrix_local *ml = malloc(sizeof(matrix_local));
    ml->size = size;
    ml->row_start = row_start;
    ml->row_end = row_end;
    ml->row_ptr = local_row_ptr;
    ml->col_idx = local_col_idx;
    ml->values = local_values;
    ml->inv_diag = local_inv_diag;
    ml->b = local_b;
    ml->proc_row_start = proc_row_start;
    ml->process_count = process_count;

    return ml;
}
