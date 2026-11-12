#!/usr/bin/env bash
# Prueba técnica pySCENIC 0.12.1 sobre el mini-dataset oficial de SCENICprotocol.
# Uso: bash scripts/run_test_oficial.sh {grn|ctx|aucell}
set -euo pipefail
cd "$(dirname "$0")/.."
IMG=docker.io/aertslab/pyscenic:0.12.1
RUN="podman run --rm -v $PWD/test_oficial:/in:ro,Z -v $PWD/test_oficial/v2:/db:ro,Z -v $PWD/results/test_oficial:/out:Z $IMG"
case "$1" in
  grn)    $RUN pyscenic grn /in/expr_mat_tiny.loom /in/test_TFs_tiny.txt \
            -o /out/adjacencies.tsv --num_workers 4 --seed 777 ;;
  ctx)    $RUN pyscenic ctx /out/adjacencies.tsv /db/genome-ranking.genes_vs_motifs.rankings.feather \
            --annotations_fname /in/motifs.tbl --expression_mtx_fname /in/expr_mat_tiny.loom \
            --mode custom_multiprocessing --num_workers 4 --min_genes 1 \
            --output /out/regulons.csv ${CTX_EXTRA:-} ;;
  aucell) $RUN pyscenic aucell /in/expr_mat_tiny.loom /out/regulons.csv \
            -o /out/auc_mtx.csv --num_workers 4 --seed 777 ;;
esac
