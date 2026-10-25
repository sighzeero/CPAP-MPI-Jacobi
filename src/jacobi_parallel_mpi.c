#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "matrix.h"
#include "jacobi.h"
/*
 * Jacobi distribuito su una matrice CSR locale.
 *
 * Il piano halo viene costruito una volta analizzando le colonne esterne
 * al blocco di righe del processo. Gli indici sono ordinati, deduplicati
 * e raggruppati per proprietario. MPI_Alltoall scambia i conteggi delle
 * richieste; MPI_Alltoallv comunica gli indici dei valori da fornire.
 * Il piano identifica anche le righe che dipendono da componenti remote.
 *
 * Ogni iterazione avvia lo scambio dei valori con MPI_Ialltoallv e calcola
 * le righe interne. Dopo MPI_Wait, i valori ricevuti vengono inseriti nel
 * vettore di lettura e si calcolano le righe di confine. Due buffer si
 * scambiano di ruolo al termine di ciascuna iterazione.
 *
 * A intervalli definiti da CONV_CHECK_INTERVAL, MPI_Allreduce combina
 * le somme dei quadrati accumulate nell'iterazione precedente per il
 * controllo dell'errore relativo. MPI_Gatherv raccoglie infine la
 * soluzione completa sul processo 0.
 */


typedef struct
{
        int local_rows;

        // valori da RICEVERE (indici GLOBALI esterni di cui ho bisogno)
        int n_recv;
        int *recv_indices; // ordinati, raggruppati per processo proprietario
        int *recv_counts;  // recv_counts[p] = quanti valori ricevo da p
        int *recv_displs;

        // valori da INVIARE (righe MIE richieste dagli altri)
        int n_send;
        int *send_indices_local; // offset locali (0-based) delle mie righe da impacchettare
        int *send_counts;                // send_counts[p] = quanti valori invio a p
        int *send_displs;

        // righe locali il cui calcolo dipende da almeno un valore esterno
        bool *row_is_boundary;
} halo_plan;

static int cmp_int(const void *a, const void *b)
{
        int ia = *(const int *)a, ib = *(const int *)b;
        return (ia > ib) - (ia < ib);
}

// Ricerca binaria del processo proprietario di un indice globale.
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

static halo_plan *build_halo_plan(matrix_local *m, int process_num, int process_count)
{
        (void)process_num;

        int row_start = m->row_start;
        int row_end = m->row_end;
        int local_rows = row_end - row_start + 1;
        int local_nnz = m->row_ptr[local_rows];


        // 1) Raccolgo gli indici esterni e marco le righe di confine.
        int *ext_tmp = malloc((local_nnz > 0 ? local_nnz : 1) * sizeof(int));
        int ext_count = 0;
        bool *is_boundary = calloc(local_rows, sizeof(bool));

        for (int i = 0; i < local_rows; i++)
        {
                int s = m->row_ptr[i], e = m->row_ptr[i + 1];
                for (int idx = s; idx < e; idx++)
                {
                        int col = m->col_idx[idx];
                        if (col < row_start || col > row_end)
                        {
                                ext_tmp[ext_count++] = col;
                                is_boundary[i] = true;
                        }
                }
        }

        //--------------------------------------------------------------
        // 2) ordino e deduplico -> recv_indices
        //--------------------------------------------------------------
        qsort(ext_tmp, ext_count, sizeof(int), cmp_int);
        int n_recv = 0;
        for (int i = 0; i < ext_count; i++)
        {
                if (i == 0 || ext_tmp[i] != ext_tmp[i - 1])
                {
                        ext_tmp[n_recv++] = ext_tmp[i];
                }
        }
        int *recv_indices = realloc(ext_tmp, (n_recv > 0 ? n_recv : 1) * sizeof(int));

        //--------------------------------------------------------------
        // 3) raggruppo per processo proprietario (l'array e' ordinato,
        //    quindi i proprietari escono automaticamente in blocchi
        //    contigui e crescenti, essendo proc_row_start monotona)
        //--------------------------------------------------------------
        int *recv_counts = calloc(process_count, sizeof(int));
        for (int i = 0; i < n_recv; i++)
        {
                int owner = find_owner(recv_indices[i], m->proc_row_start, process_count);
                recv_counts[owner]++;
        }
        int *recv_displs = malloc(process_count * sizeof(int));
        recv_displs[0] = 0;
        for (int p = 1; p < process_count; p++)
        {
                recv_displs[p] = recv_displs[p - 1] + recv_counts[p - 1];
        }

        //--------------------------------------------------------------
        // 4) scambio COUNT: dico a ciascun processo quanti valori voglio
        //    da lui (recv_counts), scopro quanti ognuno vuole da me
        //    (send_counts) -- un solo Alltoall, una tantum
        //--------------------------------------------------------------
        int *send_counts = malloc(process_count * sizeof(int));
        MPI_Alltoall(recv_counts, 1, MPI_INT, send_counts, 1, MPI_INT, MPI_COMM_WORLD);

        int *send_displs = malloc(process_count * sizeof(int));
        send_displs[0] = 0;
        for (int p = 1; p < process_count; p++)
        {
                send_displs[p] = send_displs[p - 1] + send_counts[p - 1];
        }
        int n_send = (process_count > 0) ? send_displs[process_count - 1] + send_counts[process_count - 1] : 0;

        //--------------------------------------------------------------
        // 5) scambio INDICI: invio la mia lista di richieste (recv_indices,
        //    gia' raggruppata per proprietario) e ricevo, simmetricamente,
        //    la lista esatta delle MIE righe che ciascun altro processo
        //    vuole da me (send_indices_global)
        //--------------------------------------------------------------
        int *send_indices_global = malloc((n_send > 0 ? n_send : 1) * sizeof(int));
        MPI_Alltoallv(recv_indices, recv_counts, recv_displs, MPI_INT,
                                  send_indices_global, send_counts, send_displs, MPI_INT,
                                  MPI_COMM_WORLD);

        // converto in offset locali (sono per forza righe mie, in [row_start,row_end])
        int *send_indices_local = malloc((n_send > 0 ? n_send : 1) * sizeof(int));
        for (int i = 0; i < n_send; i++)
        {
                send_indices_local[i] = send_indices_global[i] - row_start;
        }
        free(send_indices_global);

        halo_plan *hp = malloc(sizeof(halo_plan));
        hp->local_rows = local_rows;
        hp->n_recv = n_recv;
        hp->recv_indices = recv_indices;
        hp->recv_counts = recv_counts;
        hp->recv_displs = recv_displs;
        hp->n_send = n_send;
        hp->send_indices_local = send_indices_local;
        hp->send_counts = send_counts;
        hp->send_displs = send_displs;
        hp->row_is_boundary = is_boundary;

        return hp;
}

static void halo_plan_destroy(halo_plan *hp)
{
        free(hp->recv_indices);
        free(hp->recv_counts);
        free(hp->recv_displs);
        free(hp->send_indices_local);
        free(hp->send_counts);
        free(hp->send_displs);
        free(hp->row_is_boundary);
        free(hp);
}

jacobi_result *jacobi_parallel_mpi(matrix_local *m, bool verbose, int process_num, int process_count)
{
        int size = (int)m->size;
        int row_start = m->row_start;
        int row_end = m->row_end;
        int local_rows = row_end - row_start + 1;

        int *row_ptr = m->row_ptr;
        int *col_idx = m->col_idx;
        double *values = m->values;
        double *inv_diag = m->inv_diag;
        double *b = m->b;
        int *proc_row_start = m->proc_row_start;

        int i, k = 0, terminar = 0;
        double norma = INFINITY;
        bool converged = false;

        halo_plan *hp = build_halo_plan(m, process_num, process_count);

        if (verbose)
        {
                printf("[rank %d] halo: %d valori da ricevere, %d valori da inviare (righe locali = %d)\n",
                           process_num, hp->n_recv, hp->n_send, local_rows);
        }

        double *xa = malloc(size * sizeof(double));
        double *xb = malloc(size * sizeof(double));
        for (i = 0; i < size; i++)
                xa[i] = 0.0;

        double *x0 = xa;
        double *x = xb;

        double *send_buf = malloc((hp->n_send > 0 ? hp->n_send : 1) * sizeof(double));
        double *recv_buf = malloc((hp->n_recv > 0 ? hp->n_recv : 1) * sizeof(double));

        double local_n1 = 0.0;
        double local_n2 = 0.0;

        double local_vars[2];
        double global_vars[2];

        MPI_Request halo_request;

        while (k < MAX_ITERATIONS && !terminar)
        {

                //----------------------------------------------------------
                // 1. PACK ALONE
                //----------------------------------------------------------
                for (i = 0; i < hp->n_send; i++)
                {
                        send_buf[i] = x0[row_start + hp->send_indices_local[i]];
                }

                //----------------------------------------------------------
                // 2. AVVIO COMUNICAZIONE ALONE
                //----------------------------------------------------------
                MPI_Ialltoallv(send_buf,
                                           hp->send_counts,
                                           hp->send_displs,
                                           MPI_DOUBLE,
                                           recv_buf,
                                           hp->recv_counts,
                                           hp->recv_displs,
                                           MPI_DOUBLE,
                                           MPI_COMM_WORLD,
                                           &halo_request);

                //----------------------------------------------------------
                // 3. CONTROLLO NORMA ITERAZIONE PRECEDENTE
                //----------------------------------------------------------
                if (k > 0 && (k % CONV_CHECK_INTERVAL == 0))
                {
                        MPI_Allreduce(local_vars,
                                                  global_vars,
                                                  2,
                                                  MPI_DOUBLE,
                                                  MPI_SUM,
                                                  MPI_COMM_WORLD);

                        norma = global_vars[1] > 0.0 ? sqrt(global_vars[0] / global_vars[1])
                                : (global_vars[0] == 0.0 ? 0.0 : INFINITY);

                        if (verbose && process_num == 0)
                        {
                                printf("\nk = %d, norma = %.20f\n",
                                           k, norma);
                        }

                        if (norma <= precision || !isfinite(norma))
                        {
                                terminar = 1;
                                converged = isfinite(norma) && norma <= precision;

                                MPI_Wait(&halo_request, MPI_STATUS_IGNORE);

                                break;
                        }
                }

                //----------------------------------------------------------
                // 4. CALCOLO RIGHE INTERNE
                //----------------------------------------------------------
                local_n1 = 0.0;
                local_n2 = 0.0;

                for (int local_i = 0; local_i < local_rows; local_i++)
                {

                        if (hp->row_is_boundary[local_i])
                                continue;

                        int global_i = row_start + local_i;

                        double soma = 0.0;

                        for (int idx = row_ptr[local_i];
                                 idx < row_ptr[local_i + 1];
                                 idx++)
                        {
                                soma += values[idx] * x0[col_idx[idx]];
                        }

                        x[global_i] =
                                (b[local_i] - soma) * inv_diag[local_i];

                        double diff = x[global_i] - x0[global_i];

                        local_n1 += diff * diff;
                        local_n2 += x[global_i] * x[global_i];
                }

                //----------------------------------------------------------
                // 5. ATTENDO ALONE
                //----------------------------------------------------------
                MPI_Wait(&halo_request,
                                 MPI_STATUS_IGNORE);

                //----------------------------------------------------------
                // 6. AGGIORNO VALORI RICEVUTI
                //----------------------------------------------------------
                for (i = 0; i < hp->n_recv; i++)
                {
                        x0[hp->recv_indices[i]] = recv_buf[i];
                }

                //----------------------------------------------------------
                // 7. CALCOLO RIGHE DI BORDO
                //----------------------------------------------------------
                for (int local_i = 0; local_i < local_rows; local_i++)
                {

                        if (!hp->row_is_boundary[local_i])
                                continue;

                        int global_i = row_start + local_i;

                        double soma = 0.0;

                        for (int idx = row_ptr[local_i];
                                 idx < row_ptr[local_i + 1];
                                 idx++)
                        {
                                soma += values[idx] * x0[col_idx[idx]];
                        }

                        x[global_i] =
                                (b[local_i] - soma) * inv_diag[local_i];

                        double diff = x[global_i] - x0[global_i];

                        local_n1 += diff * diff;
                        local_n2 += x[global_i] * x[global_i];
                }

                //----------------------------------------------------------
                // 8. SALVO PER ALLREDUCE ITERAZIONE SUCCESSIVA
                //----------------------------------------------------------
                local_vars[0] = local_n1;
                local_vars[1] = local_n2;

                //----------------------------------------------------------
                // 9. SWAP
                //----------------------------------------------------------
                double *tmp = x0;
                x0 = x;
                x = tmp;

                k++;
        }


        /* Al limite, controlla anche l’ultimo aggiornamento calcolato. */
        if (!terminar) {
                MPI_Allreduce(local_vars, global_vars, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                norma = global_vars[1] > 0.0 ? sqrt(global_vars[0] / global_vars[1])
                        : (global_vars[0] == 0.0 ? 0.0 : INFINITY);
                converged = isfinite(norma) && norma <= precision;
        }
        /* Dopo lo swap x0 contiene sempre l’iterata piu’ recente. */
        double *x_final = x0;

        //--------------------------------------------------------------
        // MPI_Gatherv assembla la soluzione completa sul processo 0.
        //--------------------------------------------------------------
        int *row_counts = NULL, *row_displs = NULL;
        double *gathered = NULL;

        if (process_num == 0)
        {
                row_counts = malloc(process_count * sizeof(int));
                row_displs = malloc(process_count * sizeof(int));
                for (int p = 0; p < process_count; p++)
                {
                        row_counts[p] = proc_row_start[p + 1] - proc_row_start[p];
                        row_displs[p] = proc_row_start[p];
                }
                gathered = malloc(size * sizeof(double));
        }

        MPI_Gatherv(&x_final[row_start], local_rows, MPI_DOUBLE,
                                gathered, row_counts, row_displs, MPI_DOUBLE,
                                0, MPI_COMM_WORLD);

        jacobi_result *res = NULL;
        if (process_num == 0)
        {
                res = malloc(sizeof(jacobi_result));
                res->x = gathered;
                res->e = norma;
                res->k = k;
                res->converged = converged;
                free(row_counts);
                free(row_displs);
        }

        free(xa);
        free(xb);
        free(send_buf);
        free(recv_buf);
        halo_plan_destroy(hp);

        return res;
}
