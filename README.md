# Jacobi distribuito con MPI

Il progetto deriva da [utluiz/distributed-jacobi](https://github.com/utluiz/distributed-jacobi). L'obiettivo è migliorare le prestazioni del codice per poi studiarne la scalabilità.

## Funzionamento e miglioramenti

Il metodo risolve un sistema lineare `Ax = b` aggiornando iterativamente il vettore soluzione. Ogni processo calcola un blocco di righe usando i valori della soluzione dell'iterazione precedente.

Le modifiche principali riguardano la rappresentazione della matrice e l'organizzazione del calcolo:

- **Matrice sparsa in formato CSR.** Valori, indici di colonna e puntatori di riga sono memorizzati in array contigui, al posto delle liste separate per riga. La diagonale è conservata a parte e ne viene precalcolato il reciproco, evitando divisioni nel ciclo principale.
- **Caricamento distribuito.** Ogni processo legge una porzione del file Matrix Market. Le righe vengono poi distribuite in blocchi contigui bilanciati per numero di elementi non nulli fuori diagonale; ciascun processo conserva soltanto la propria parte della matrice.
- **Scambio dei soli valori necessari.** Un piano halo, costruito una sola volta, individua le componenti remote richieste dalle righe locali. Al posto di ricostruire l'intero vettore con `MPI_Allgatherv` a ogni iterazione, vengono scambiati soltanto questi valori.
- **Sovrapposizione tra comunicazione e calcolo.** Lo scambio viene avviato con `MPI_Ialltoallv`. Nel frattempo si calcolano le righe interne, che dipendono solo da valori locali; dopo il completamento della comunicazione si calcolano le righe di confine.
- **Doppio buffer.** Due vettori si scambiano di ruolo tramite puntatori, evitando la copia completa della soluzione a ogni passo. I buffer hanno dimensione globale, ma si aggiornano soltanto le componenti locali e quelle necessarie allo scambio halo.
- **Meno sincronizzazioni.** Il controllo globale di convergenza viene eseguito ogni 5 iterazioni, con `MPI_Allreduce` sui contributi dell'iterazione precedente. La soluzione completa viene raccolta sul processo 0 con un solo `MPI_Gatherv` finale.

Il criterio di arresto è l'errore relativo tra due iterate successive, con soglia `1e-5` e limite di 10000 iterazioni.

Il codice è suddiviso in `jacobi.c` per l'avvio MPI, `matrix.c` per il caricamento e la distribuzione della matrice, `jacobi_parallel_mpi.c` per il solver, `timer.c` e `results.c` per tempi e risultati.

## Compilazione

Servono un compilatore C, MPI 3 o superiore e make. Python 2.7 o 3 è usato dai test e dallo script di riepilogo. Tutti i comandi seguenti vanno eseguiti dalla radice del progetto.

```bash
make
```

In alternativa, compilazione diretta:

```bash
mpicc -std=gnu99 -O3 -o ./src/jacobi_mpi_bin \
  ./src/jacobi.c ./src/jacobi_parallel_mpi.c ./src/matrix.c \
  ./src/timer.c ./src/results.c -lm
mkdir -p results
```

## Test

Per verificare soluzione e residuo su matrici di test con 1, 2 e 4 processi:

```bash
make test
```

Per eseguire manualmente il piccolo esempio incluso e stampare la soluzione:

```bash
mpirun -np 2 ./src/jacobi_mpi_bin ./big/test10.mtx M 2 1
```

Gli argomenti sono `matrice.mtx M numero_processi verbose`: il numero di processi deve coincidere con `-np`; `verbose` vale `1` per stampare i dettagli, `0` per l'esecuzione degli esperimenti.

## Esecuzione sul cluster

Configurare `machinefile.txt` con i nomi dei nodi. Il comando seguente avvia 64 processi, con 16 processi per nodo e binding al core:

```bash
mpirun --mca btl self,vader,openib \
  -machinefile ./machinefile.txt \
  -np 64 \
  --map-by ppr:16:node \
  --bind-to core \
  ./src/jacobi_mpi_bin ./big/ecology2.mtx M 64 0
```

Il progetto e i file di input devono essere accessibili agli stessi percorsi da tutti i nodi.

## Esperimenti con script

Lo script esegue le configurazioni con **1, 2, 4, 8, 16, 32, 50 e 64 processi**, con 5 ripetizioni per configurazione. Usa gli hostfile in `machinefiles/`, da adattare ai nodi del cluster.

```bash
bash scripts/run_experiments_MPI.sh
```

La matrice predefinita è `ecology2`: lo script la cerca automaticamente in `big/` e nelle sue sottocartelle. Accetta anche il solo nome (ad esempio `pwtk`) oppure un percorso esplicito. Se trova più copie con lo stesso nome richiede il percorso esatto. Il programma viene compilato automaticamente se necessario. Per scegliere un'altra matrice o cambiare il numero di ripetizioni:

```bash
REPEATS=3 bash scripts/run_experiments_MPI.sh ./big/pwtk.mtx
```

Ogni campagna crea una directory in `results/` con log e misure separati per esecuzione. Il file `summary.csv` riporta i tempi mediani, minimo e massimo, lo speedup `T(1)/T(p)` e l'efficienza `speedup/p`.

I tempi misurati escludono il caricamento della matrice e includono la costruzione del piano halo, il calcolo e la raccolta finale. Il tempo utilizzato è il massimo del tempo trascorso sui processi MPI.

## Matrici di input

Sono incluse `test10`, `ecology2` e `pwtk` in formato Matrix Market. Se manca un file dei termini noti `nome_b.mtx` accanto alla matrice, viene generato `b = A·1`, così la soluzione attesa è il vettore di uno. Un file dei termini noti esplicito deve contenere `N 1` seguito dagli N valori.

La licenza del progetto di origine è riportata in `LICENSE`.

Gli script sono compatibili con Python 2.7.5. Il comando predefinito è `python`; per usarne un altro: `make test PYTHON=python3` oppure `PYTHON=python3 bash scripts/run_experiments_MPI.sh`. Per passare opzioni MPI ai test sul server:

```bash
MPI_TEST_FLAGS="--mca btl self,vader,openib --bind-to none" make test
```

Per controllare quale matrice verrà usata senza compilare o avviare gli esperimenti:

```bash
bash scripts/run_experiments_MPI.sh --dry-run
bash scripts/run_experiments_MPI.sh --dry-run pwtk
```

Il raggiungimento del limite di iterazioni senza convergenza non interrompe la campagna: il tempo viene salvato e nel log resta `converged=0`. Il riepilogo riporta anche il numero di iterazioni e quante esecuzioni sono convergenti. Per confrontare i tempi a parità di lavoro controllare che il numero di iterazioni coincida; piccole differenze numeriche tra riduzioni MPI possono influenzare il punto di arresto. NaN/Inf e gli errori di esecuzione continuano a interrompere il test.

Configurazioni della campagna:

| Processi MPI | Nodi | Distribuzione dei processi |
|---|---|---|
| 1, 2, 4, 8, 16 | 1 | Tutti sul primo nodo |
| 32 | 2 | 16 per nodo |
| 50 | 3 | 17, 17, 16 |
| 64 | 4 | 16 per nodo |

Il punto a 50 processi su 3 nodi replica il numero di nodi della tabella, ma supera i 48 core fisici disponibili: usa oversubscription e consente la condivisione dei core (`--oversubscribe --bind-to core:overload-allowed`). Non equivale a 50 core fisici distinti; verificare il binding effettivo nel log. Gli altri punti usano binding al core senza oversubscription.
