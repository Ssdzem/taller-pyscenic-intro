#!/usr/bin/env bash
# Pipeline pySCENIC 0.12.1 sobre la matriz demo PBMC10k (900 células × 10,120 genes).
# Uso: bash scripts/run_pbmc_pipeline.sh {grn|ctx|aucell}
set -euo pipefail
cd "$(dirname "$0")/.."
IMG=docker.io/aertslab/pyscenic:0.12.1
W=${WORKERS:-20}
RUN="podman run --rm -v $PWD/data:/data:ro,Z -v $PWD/resources:/res:ro,Z -v $PWD/results:/out:Z $IMG"
DB1=/res/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
DB2=/res/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
case "$1" in
  grn)    $RUN pyscenic grn /data/expression_demo.tsv /res/allTFs_hg38.txt \
            -o /out/adjacencies.tsv --num_workers $W --seed 2026 ;;
  ctx)    $RUN pyscenic ctx /out/adjacencies.tsv $DB1 $DB2 \
            --annotations_fname /res/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
            --expression_mtx_fname /data/expression_demo.tsv \
            --mode custom_multiprocessing --num_workers $W \
            --output /out/regulons.csv ;;
  aucell) $RUN pyscenic aucell /data/expression_demo.tsv /out/regulons.csv \
            -o /out/auc_mtx.csv --num_workers $W --seed 2026 ;;
esac
