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
if [[ "$DRY_RUN" -eq 1 ]]; then
    exit 0
fi
REPEATS="${REPEATS:-5}"
HOSTFILE_DIR="${HOSTFILE_DIR:-$ROOT/machinefiles}"
CAMPAIGN="${CAMPAIGN:-$ROOT/results/$(date '+%Y%m%d_%H%M%S')}"
BIN="$ROOT/src/jacobi_mpi_bin"
# make ricompila soltanto se necessario.
make -C "$ROOT"
[[ -x "$BIN" ]] || { echo "Eseguibile non disponibile: $BIN" >&2; exit 1; }
[[ "$REPEATS" =~ ^[1-9][0-9]*$ ]] || { echo 'REPEATS deve essere positivo.' >&2; exit 1; }
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
for NP in 1 2 4 8 16 32 64; do
    case "$NP" in
        1) NODES=1; MAP=slot ;;
        2) NODES=2; MAP=node ;;
        4) NODES=4; MAP=node ;;
        *) NODES=4; MAP="ppr:$((NP / 4)):node" ;;
    esac
    MF="$HOSTFILE_DIR/mf_${NODES}nodes.txt"
    [[ -f "$MF" ]] || { echo "Hostfile assente: $MF" >&2; exit 1; }
    for ((REP=1; REP<=REPEATS; REP++)); do
        RUN="$CAMPAIGN/p${NP}_r${REP}"
        mkdir -p "$RUN/results"
        cp "$MF" "$RUN/hostfile.txt"
        CMD=(mpirun --mca btl self,vader,openib -np "$NP" -machinefile "$RUN/hostfile.txt"
             --map-by "$MAP" --bind-to core --report-bindings
             "$BIN" "$MATRIX" M "$NP" 0)
        printf '%q ' "${CMD[@]}" > "$RUN/command.txt"
        printf '\n' >> "$RUN/command.txt"
        printf '%s,%s,%s,%s,%s\n' "$NP" "$REP" "$NODES" "$MF" "p${NP}_r${REP}" >> "$CAMPAIGN/manifest.csv"
        echo "MPI: $NP processi, ripetizione $REP/$REPEATS"
        (cd "$RUN"; "${CMD[@]}" 2>&1 | tee run.log)
    done
done
"${PYTHON:-python}" "$ROOT/scripts/summarize.py" "$CAMPAIGN"
