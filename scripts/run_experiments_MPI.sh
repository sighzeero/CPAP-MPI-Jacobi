#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=1
    shift
fi
[[ "$#" -le 1 ]] || { echo "Uso: $0 [--dry-run] [matrice.mtx|nome]" >&2; exit 1; }
REQUEST="${1:-ecology2.mtx}"
if [[ -f "$REQUEST" ]]; then
    MATRIX="$REQUEST"
elif [[ -f "$ROOT/$REQUEST" ]]; then
    MATRIX="$ROOT/$REQUEST"
else
    # Cerca per nome in big/, anche quando le matrici sono in sottocartelle.
    NAME="$(basename "$REQUEST")"
    [[ "$NAME" == *.mtx ]] || NAME="${NAME}.mtx"
    MATCHES=()
    while IFS= read -r -d '' FILE; do
        MATCHES+=("$FILE")
    done < <(find "$ROOT/big" -type f -name "$NAME" -print0)
    if [[ "${#MATCHES[@]}" -eq 0 ]]; then
        echo "Matrice non trovata: $REQUEST (ricerca in $ROOT/big)." >&2
        exit 1
    elif [[ "${#MATCHES[@]}" -gt 1 ]]; then
        echo "Più matrici con nome $NAME: specificare il percorso esatto." >&2
        printf '  %s\n' "${MATCHES[@]}" >&2
        exit 1
    fi
    MATRIX="${MATCHES[0]}"
fi
MATRIX="$(cd "$(dirname "$MATRIX")" && pwd)/$(basename "$MATRIX")"
echo "Matrice: $MATRIX"
REPEATS="${REPEATS:-5}"
HOSTFILE_DIR="${HOSTFILE_DIR:-$ROOT/machinefiles}"
CAMPAIGN="${CAMPAIGN:-$ROOT/results/$(date '+%Y%m%d_%H%M%S')}"
BIN="$ROOT/src/jacobi_mpi_bin"
# make ricompila soltanto se necessario.
if [[ "$DRY_RUN" -eq 0 ]]; then
    make -C "$ROOT"
    [[ -x "$BIN" ]] || { echo "Eseguibile non disponibile: $BIN" >&2; exit 1; }
fi
[[ "$REPEATS" =~ ^[1-9][0-9]*$ ]] || { echo 'REPEATS deve essere positivo.' >&2; exit 1; }
if [[ "$DRY_RUN" -eq 0 ]]; then
mkdir -p "$(dirname "$CAMPAIGN")"
mkdir "$CAMPAIGN"
{
    date -u
    uname -a
    mpirun --version
    mpicc --showme:version
    shasum -a 256 "$MATRIX" "$BIN" "$ROOT"/src/* "$ROOT/Makefile"
} > "$CAMPAIGN/environment.txt"
printf 'processes,repeat,nodes,hostfile,directory\n' > "$CAMPAIGN/manifest.csv"
fi
for NP in 1 2 4 8 16 32 50 64; do
    OVERSUBSCRIBE=0
    BIND=core
    case "$NP" in
        1|2|4|8|16) NODES=1; MAP=slot ;;
        32) NODES=2; MAP=ppr:16:node ;;
        50)
            # Replica il punto della tabella: 50 rank su 3 nodi (48 core).
            # Mapping per nodo: 17, 17, 16 rank; due core sono condivisi.
            NODES=3; MAP=node
            OVERSUBSCRIBE=1
            BIND=core:overload-allowed
            ;;
        64) NODES=4; MAP=ppr:16:node ;;
    esac
    MF="$HOSTFILE_DIR/mf_${NODES}nodes.txt"
    [[ -f "$MF" ]] || { echo "Hostfile assente: $MF" >&2; exit 1; }
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "MPI: $NP processi, $NODES nodi, map=$MAP, bind=$BIND, oversubscribe=$OVERSUBSCRIBE"
        continue
    fi
    for ((REP=1; REP<=REPEATS; REP++)); do
        RUN="$CAMPAIGN/p${NP}_r${REP}"
        mkdir -p "$RUN/results"
        cp "$MF" "$RUN/hostfile.txt"
        CMD=(mpirun --mca btl self,vader,openib -np "$NP" -machinefile "$RUN/hostfile.txt"
             --map-by "$MAP" --bind-to "$BIND" --report-bindings)
        if [[ "$OVERSUBSCRIBE" -eq 1 ]]; then
            CMD+=(--oversubscribe)
        fi
        CMD+=("$BIN" "$MATRIX" M "$NP" 0)
        printf '%q ' "${CMD[@]}" > "$RUN/command.txt"
        printf '\n' >> "$RUN/command.txt"
        printf '%s,%s,%s,%s,%s\n' "$NP" "$REP" "$NODES" "$MF" "p${NP}_r${REP}" >> "$CAMPAIGN/manifest.csv"
        echo "MPI: $NP processi, $NODES nodi, ripetizione $REP/$REPEATS"
        (cd "$RUN"; "${CMD[@]}" 2>&1 | tee run.log)
    done
done
if [[ "$DRY_RUN" -eq 0 ]]; then
    "${PYTHON:-python}" "$ROOT/scripts/summarize.py" "$CAMPAIGN"
fi
