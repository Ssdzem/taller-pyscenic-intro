#!/usr/bin/env bash
# Paso 52: "Ejecutar todo" de un notebook en un entorno tipo Colab gratuito (2 CPU, 12 GB, Python 3.12).
# Datos del taller desde la copia local (TALLER_LOCAL) y recursos cisTarget desde resources/ (TALLER_RECURSOS),
# para no descargar 412 MB en cada prueba. pySCENIC sí se instala desde PyPI, como en Colab.
# Uso: bash scripts/52_test_notebook_2cpu.sh notebooks/X.ipynb salida.ipynb [cpus] [núcleos, p. ej. 16-17]
NB=$1; OUT=$2; CPUS=${3:-2}; CPUSET=${4:-0-$((CPUS-1))}
cd "$(dirname "$0")/.."
podman run --rm --cpus=$CPUS --cpuset-cpus=$CPUSET --memory=12g -v $PWD:/w:Z -e TALLER_LOCAL=/w -e TALLER_RECURSOS=/w/resources \
  -e TALLER_DATOS_10X=/w/data/original -e NB_TIMEOUT=${NB_TIMEOUT:-900} localhost/taller-colabtest:py312 \
  bash -c "mkdir -p /tmp/run && cd /tmp/run && python /w/scripts/test_notebook.py /w/$NB /w/$OUT"
