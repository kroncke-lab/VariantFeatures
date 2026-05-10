#!/bin/bash
# Per-gene full annotation pipeline: enumerate -> AM -> ANNOVAR -> VEP -> REST drains.
# Run several in parallel by invoking this script multiple times against different genes.
#
# Usage: scripts/full_gene_pipeline.sh GENE [DB_PATH] [--skip-local] [--skip-rest]
set -euo pipefail
GENE="${1:?usage: full_gene_pipeline.sh GENE [DB_PATH] [--skip-local] [--skip-rest]}"
GENE_LOWER="$(printf '%s' "$GENE" | tr '[:upper:]' '[:lower:]')"
DB="${2:-/tmp/${GENE_LOWER}_jobs.db}"
SKIP_LOCAL=0
SKIP_REST=0
for arg in "${@:3}"; do
    case "$arg" in
        --skip-local) SKIP_LOCAL=1 ;;
        --skip-rest) SKIP_REST=1 ;;
    esac
done

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON="$PROJECT_ROOT/.venv/bin/python"
log() { printf '[%s %s] %s\n' "$(date +%H:%M:%S)" "$GENE" "$*"; }

cd "$PROJECT_ROOT"

count_variants() {
    $PYTHON -c "import sqlite3, os.path; print(sqlite3.connect('$1').execute('SELECT COUNT(*) FROM variants').fetchone()[0]) if os.path.exists('$1') else print(0)" 2>/dev/null || echo 0
}

if [ "$(count_variants "$DB")" = "0" ]; then
    log "enumerate"
    $PYTHON -m variantfeatures enumerate --gene "$GENE" --db "$DB"
else
    log "skip enumerate (DB already has variants)"
fi

# Always ensure all working sources are queued (idempotent: existing rows are left alone).
log "queue all sources"
for source in alphamissense annovar vep clingen_ar myvariant; do
    $PYTHON -m variantfeatures queue --source "$source" --gene "$GENE" --db "$DB" 2>&1 \
        | grep -v urllib3 | grep -v NotOpenSSL | grep -E "Queued|already" || true
done

if [ "$SKIP_LOCAL" = "0" ]; then
    log "alphamissense (single file pass)"
    $PYTHON -m variantfeatures alphamissense-run --db "$DB" 2>&1 | grep -v urllib3 | grep -v NotOpenSSL | tail -2

    log "annovar (refGene + clinvar in one batch)"
    ANNOVAR_HOME="$PROJECT_ROOT/annovar" ANNOVAR_DB="$PROJECT_ROOT/annovar/humandb" \
        $PYTHON -m variantfeatures annovar-run --db "$DB" \
            --build hg38 --protocols refGeneWithVer,clinvar_20240611 --operations g,f 2>&1 \
            | grep -v urllib3 | grep -v NotOpenSSL | tail -2

    log "vep (cache + plugins, one batch)"
    $PYTHON -m variantfeatures vep-run --db "$DB" 2>&1 | grep -v urllib3 | grep -v NotOpenSSL | tail -2
else
    log "skip local annotators"
fi

if [ "$SKIP_REST" = "0" ]; then
    log "clingen_ar (~0.2s/var)"
    $PYTHON -m variantfeatures annotate-pending --source clingen_ar --db "$DB" --rate-limit 0.2 2>&1 \
        | grep -v urllib3 | grep -v NotOpenSSL | tail -2

    log "myvariant (~0.1s/var)"
    $PYTHON -m variantfeatures annotate-pending --source myvariant --db "$DB" --rate-limit 0.1 2>&1 \
        | grep -v urllib3 | grep -v NotOpenSSL | tail -2
else
    log "skip rest drains"
fi

# Gene-level + NMD-rule annotators (cheap, no flags to skip).
log "gene-constraint (gnomAD pLI / LOEUF / mis_z)"
$PYTHON -m variantfeatures gene-constraint -g "$GENE" --db "$DB" 2>&1 \
    | grep -v urllib3 | grep -v NotOpenSSL | tail -2

log "nmd-rule (last-exon + 50-nt rule for stop_gained)"
$PYTHON -m variantfeatures nmd-rule -g "$GENE" --db "$DB" 2>&1 \
    | grep -v urllib3 | grep -v NotOpenSSL | tail -5

log "done. summary:"
$PYTHON -m variantfeatures jobs --db "$DB" 2>&1 | tail -20 | grep -v urllib3 | grep -v NotOpenSSL
