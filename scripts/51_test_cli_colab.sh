#!/usr/bin/env bash
# Paso 51 (prueba de viabilidad): ¿caben grn/ctx/aucell del CLI de pySCENIC en un Colab gratuito?
# Simula Colab: python:3.12-slim limitado a 2 CPU y 12 GB. pySCENIC se instala AISLADO con --target
# (como en la sección 4.2 del notebook) para no tocar el numpy/pandas del sistema.
#   FIX=pin   -> instala con "numpy<2" "pandas<2.3" (pySCENIC 0.12.1 usa np.object, eliminado en numpy 2)
#   FIX=patch -> numpy 2 + "pandas<3" (pandas 3 cambió groupby.apply) y corrige con sed las 4 líneas np.object/np.float
# Uso: bash scripts/51_test_cli_colab.sh {r1|full} {pin|patch}
MODE=${1:-r1}; FIX=${2:-pin}
cd "$(dirname "$0")/.."
podman run --rm --cpus=2 --memory=12g -v $PWD:/w:ro,Z docker.io/library/python:3.12-slim bash -c '
  MODE='"$MODE"'; FIX='"$FIX"'
  t() { local s=$(date +%s); "$@"; local rc=$?; echo "### exit=$rc segundos=$(( $(date +%s) - s )) :: $*" >&3; }
  exec 3>&1
  mkdir -p /tmp/t && cd /tmp/t
  pip install -q --root-user-action=ignore pandas >/dev/null 2>&1   # el "Colab" ya trae pandas
  P=/tmp/t/pkgs
  if [ "$FIX" = pin ]; then
    t pip install -q --root-user-action=ignore --target $P "pyscenic==0.12.1" "setuptools<81" "numpy<2" "pandas<2.3" > pip.log 2>&1
  else
    t pip install -q --root-user-action=ignore --target $P "pyscenic==0.12.1" "setuptools<81" "pandas<3" > pip.log 2>&1
    sed -i "s/np\.object,/object,/; s/dtype=np\.float)/dtype=float)/" $P/pyscenic/transform.py $P/pyscenic/rss.py
  fi
  grep -iE "error|conflict" pip.log | head -5
  export PYTHONPATH=$P PATH=$P/bin:$PATH
  python -c "import pyscenic, numpy, pandas, dask; print(\"pyscenic\", pyscenic.__version__, \"numpy\", numpy.__version__, \"pandas\", pandas.__version__, \"dask\", dask.__version__)"
  echo "cpu.max=$(cat /sys/fs/cgroup/cpu.max) mem.max=$(cat /sys/fs/cgroup/memory.max)  FIX=$FIX MODE=$MODE"
  t python -I /w/scripts/51_prep_reduced.py /w/data/expression_demo.tsv.gz /w/resources/allTFs_hg38.txt
  E=expr_$MODE.tsv
  t timeout 600 pyscenic grn $E tfs.txt -o adj_dask.tsv --num_workers 2 --seed 2026 > grn_dask.log 2>&1
  grep -E "Error|error" grn_dask.log | tail -3; echo "filas adj_dask: $(wc -l < adj_dask.tsv 2>/dev/null)"
  t timeout 3000 arboreto_with_multiprocessing.py $E tfs.txt --method grnboost2 --output adj.tsv --num_workers 2 --seed 2026 > grn_mp.log 2>&1
  grep -E "Error|error" grn_mp.log | tail -3; echo "filas adj: $(wc -l < adj.tsv 2>/dev/null)"
  R=/w/resources; TBL=$R/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
  DB500=$R/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
  DB10k=$R/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
  nreg() { python -c "from pyscenic.cli.utils import load_signatures as l; print(len(l(\"$1\")))" 2>&1 | tail -1; }
  t pyscenic ctx adj.tsv $DB500 --annotations_fname $TBL --expression_mtx_fname $E --mode custom_multiprocessing --num_workers 2 --output reg_1db.csv > ctx1.log 2>&1
  grep -E "Error|error" ctx1.log | tail -3; echo "regulones 1 base: $(nreg reg_1db.csv)  mem.peak=$(cat /sys/fs/cgroup/memory.peak)"
  t pyscenic ctx adj.tsv $DB500 $DB10k --annotations_fname $TBL --expression_mtx_fname $E --mode custom_multiprocessing --num_workers 2 --output reg_2db.csv > ctx2.log 2>&1
  grep -E "Error|error" ctx2.log | tail -3; echo "regulones 2 bases: $(nreg reg_2db.csv)  mem.peak=$(cat /sys/fs/cgroup/memory.peak)"
  t pyscenic aucell $E reg_2db.csv -o auc.csv --num_workers 2 --seed 2026 > auc.log 2>&1
  grep -E "Error|error" auc.log | tail -3; python -c "import pandas as pd; print(\"auc\", pd.read_csv(\"auc.csv\", index_col=0).shape)" 2>&1 | tail -1
  echo "### FIN mem.peak=$(cat /sys/fs/cgroup/memory.peak)"
'
