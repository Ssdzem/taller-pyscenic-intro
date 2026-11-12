#!/usr/bin/env bash
# Paso 35 (ensayo local): pySCENIC 0.12.1 con pip en un Python "tipo Colab" + AUCell real sobre la matriz demo.
# Uso: bash scripts/35_test_pip_pyscenic.sh 3.12
PY=${1:-3.12}
cd "$(dirname "$0")/.."
podman run --rm -v $PWD:/w:ro,Z -w /w docker.io/library/python:${PY}-slim sh -c '
  python --version
  start=$(date +%s)
  pip install --no-cache-dir -q --root-user-action=ignore "pyscenic==0.12.1" "setuptools<81" > /tmp/pip.log 2>&1; rc=$?
  echo "pip exit=$rc  segundos=$(( $(date +%s) - start ))"; grep -i -E "error|conflict" /tmp/pip.log | head
  python -I scripts/35_test_pip_aucell.py 2>&1 | grep -v -i warning | tail -15
'
