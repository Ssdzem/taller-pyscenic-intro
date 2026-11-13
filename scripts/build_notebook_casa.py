"""Genera notebooks/Pipeline_completo_casa.ipynb: el pipeline completo de pySCENIC en Colab, sin límite de tiempo.

Uso: python3 scripts/build_notebook_casa.py
Dos matrices a elegir:
  - "demo":     las 900 células × 10,120 genes del taller (data/expression_demo.tsv.gz). ~30 min en Colab gratuito.
  - "completa": todo el dataset 10k PBMC de 10x tras el mismo control de calidad (~10,000 células). Varias horas.
"""
import json
import textwrap
from pathlib import Path

BASE_URL = "https://raw.githubusercontent.com/Ssdzem/taller-pyscenic-intro/main"
COLAB_URL = BASE_URL.replace("https://raw.githubusercontent.com/", "https://colab.research.google.com/github/").replace(
    "/main", "/blob/main") + "/notebooks/Pipeline_completo_casa.ipynb"
OUT = Path(__file__).resolve().parents[1] / "notebooks" / "Pipeline_completo_casa.ipynb"

cells = []


def md(src):
    cells.append({"cell_type": "markdown", "metadata": {}, "source": textwrap.dedent(src).strip("\n")})


def code(src):
    cells.append({"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
                  "source": textwrap.dedent(src).strip("\n")})


# =============================================================================
md(rf"""
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({COLAB_URL})

# 🏠 pySCENIC completo, para casa

En el taller interpretaste los resultados del servidor y corriste versiones reducidas de cada comando. Aquí corres **el pipeline completo** (`grn` → `ctx` → `aucell`), con todos los genes y las dos bases de ranking, igual que el instructor.

**Elige la matriz** en la primera celda:

| Matriz | Células × genes | Tiempo aprox. en Colab gratuito (2 núcleos) | Para qué sirve |
|---|---|---|---|
| `"demo"` | 900 × 10,120 (las del taller) | ~30 min | Reproducir los resultados del servidor y compararlos |
| `"completa"` | 10,280 células × 10,165 genes de 10x tras el control de calidad | ~2–3 h (casi todo es `grn`) | Ver el pipeline a escala real |

**Antes de empezar:**

- Usa **Entorno de ejecución → Ejecutar todo** y deja la pestaña **abierta**. Colab gratuito puede desconectarse si la pestaña queda inactiva mucho tiempo, o tras ~12 h.
- Si se desconecta, vuelve a ejecutar todo: cada paso **se salta** si su resultado ya existe en la sesión. Para no perder nada, activa `GUARDAR_EN_DRIVE`.
- Necesitas ~1.5 GB de disco y ~6 GB de RAM; Colab gratuito tiene suficiente. Si tienes más núcleos (Colab Pro o tu propia computadora), el notebook usa tantos procesos como quepan en la memoria.
- Las secciones 🖥️ del notebook del taller explican cada comando y sus opciones.
""")

code(r'''
# ✏️ Configuración
MATRIZ = "demo"            # "demo" (900 células, ~30 min) o "completa" (10,280 células, ~2–3 h)
GUARDAR_EN_DRIVE = False   # True: copia los resultados a tu Google Drive al final de cada paso
CARPETA_DRIVE = "/content/drive/MyDrive/pyscenic_casa"   # dónde guardarlos
RUTA_DRIVE_RECURSOS = ""   # si tienes las bases de ranking en tu Drive, su carpeta (si no, se descargan)
''')

# =============================================================================
md(r"""
## 0 · Preparación
""")

code(rf'''
import hashlib, os, re, shutil, subprocess, sys, tarfile, time, urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

assert MATRIZ in ("demo", "completa"), 'MATRIZ debe ser "demo" o "completa"'
BASE_URL = "{BASE_URL}"
AERTS = "https://resources.aertslab.org/cistarget"
DIEZX = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3"
NUCLEOS = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()   # núcleos asignados a esta sesión

def memoria_gb():
    """RAM disponible para esta sesión (respeta el límite del contenedor si lo hay)."""
    try:
        limite = open("/sys/fs/cgroup/memory.max").read().strip()
        if limite != "max":
            return int(limite) / 1e9
    except OSError:
        pass
    return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9
RAM_GB = memoria_gb()
os.makedirs("recursos", exist_ok=True)
os.makedirs("resultados", exist_ok=True)

# archivo -> (URL, SHA1 esperado o None, variable de entorno para pruebas del instructor)
RECURSOS = {{
    "allTFs_hg38.txt": (f"{{AERTS}}/tf_lists/allTFs_hg38.txt", "b58256aab395504ff932d563528355343e234dcd", "TALLER_RECURSOS"),
    "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl": (f"{{AERTS}}/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "8cba044ef45b389f4b85237931531906b3046182", "TALLER_RECURSOS"),
    "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather": (
        f"{{AERTS}}/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "60fe44626b5752d31bac611af6e813b1f53fa9be", "TALLER_RECURSOS"),
    "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather": (
        f"{{AERTS}}/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "d9d4e927b80b2b35344aec0a9fc0a82ae9d76759", "TALLER_RECURSOS"),
    "expression_demo.tsv.gz": (f"{{BASE_URL}}/data/expression_demo.tsv.gz", None, "TALLER_LOCAL"),
    "metadata_demo.csv": (f"{{BASE_URL}}/data/metadata_demo.csv", None, "TALLER_LOCAL"),
    "regulons_demo.csv": (f"{{BASE_URL}}/results/regulons_demo.csv", None, "TALLER_LOCAL"),
    "auc_demo.csv": (f"{{BASE_URL}}/results/auc_demo.csv", None, "TALLER_LOCAL"),
    "pbmc_10k_v3_filtered_feature_bc_matrix.h5": (f"{{DIEZX}}/pbmc_10k_v3_filtered_feature_bc_matrix.h5",
        "ed9aa3361cb5bed4032e894b4d2f9311481d7919", "TALLER_DATOS_10X"),
    "pbmc_10k_v3_analysis.tar.gz": (f"{{DIEZX}}/pbmc_10k_v3_analysis.tar.gz",
        "d36f5ee97b25b61c9f7b144b7bce2a0d60fcb5eb", "TALLER_DATOS_10X"),
}}
# en la copia local del repositorio (solo pruebas) los archivos están en subcarpetas
SUBCARPETA_LOCAL = {{"expression_demo.tsv.gz": "data", "metadata_demo.csv": "data",
                    "regulons_demo.csv": "results", "auc_demo.csv": "results"}}

def sha1(ruta):
    h = hashlib.sha1()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()

def conseguir(archivo):
    """Descarga (o copia de Drive) un archivo a recursos/, una sola vez, y comprueba su SHA1."""
    url, esperado, variable = RECURSOS[archivo]
    destino = os.path.join("recursos", archivo)
    if os.path.exists(destino) and (esperado is None or sha1(destino) == esperado):
        return destino
    t0 = time.time()
    local = os.environ.get(variable)                          # solo para pruebas del instructor
    if local:
        shutil.copy(os.path.join(local, SUBCARPETA_LOCAL.get(archivo, ""), archivo), destino)
    elif RUTA_DRIVE_RECURSOS and archivo.startswith(("hg38_", "motifs-", "allTFs")):
        if not os.path.exists("/content/drive"):
            from google.colab import drive
            drive.mount("/content/drive")
        shutil.copy(os.path.join(RUTA_DRIVE_RECURSOS, archivo), destino)
    else:
        for intento in range(3):
            try:
                urllib.request.urlretrieve(url, destino)
                break
            except Exception as e:
                if intento == 2:
                    raise RuntimeError(f"No se pudo descargar {{archivo}}: {{e}}")
                time.sleep(5)
    if esperado is not None and sha1(destino) != esperado:
        os.remove(destino)
        raise RuntimeError(f"{{archivo}} llegó incompleto o dañado (SHA1 distinto). Vuelve a ejecutar la celda.")
    print(f"  ✅ {{archivo}} ({{os.path.getsize(destino) / 1e6:.0f}} MB, {{time.time() - t0:.0f}} s)")
    return destino

def guardar(*archivos):
    """Copia resultados a Google Drive si GUARDAR_EN_DRIVE = True."""
    if not GUARDAR_EN_DRIVE:
        return
    if not os.path.exists("/content/drive"):
        from google.colab import drive
        drive.mount("/content/drive")
    os.makedirs(CARPETA_DRIVE, exist_ok=True)
    for a in archivos:
        copia = os.path.join(CARPETA_DRIVE, f"{{MATRIZ}}_{{os.path.basename(a)}}")
        shutil.copy(a, copia)
        print(f"  💾 {{a}} → {{copia}}")

if GUARDAR_EN_DRIVE:   # retomar: si ya hay resultados en Drive de una sesión anterior, los recupera
    guardar()
    for a in ["adjacencies.tsv", "regulons.csv", "auc_mtx.csv"]:
        previo = os.path.join(CARPETA_DRIVE, f"{{MATRIZ}}_{{a}}")
        if os.path.exists(previo) and not os.path.exists(os.path.join("resultados", a)):
            shutil.copy(previo, os.path.join("resultados", a))
            print(f"  ↩️  recuperado de Drive: {{a}}")

print(f"Matriz: {{MATRIZ}} | núcleos: {{NUCLEOS}} | RAM: {{RAM_GB:.0f}} GB | Python {{sys.version.split()[0]}}")
''')

md(r"""
### Instalar pySCENIC (≈ 30 s)

Se aplican los mismos arreglos que en el taller: `setuptools<81`, `pandas<3` y `np.object`/`np.float` → `object`/`float`.
""")

code(r'''
PKGS = os.path.abspath("pyscenic_pkgs")
if not os.path.exists(os.path.join(PKGS, "pyscenic")):
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "--target", PKGS,
                    "pyscenic==0.12.1", "setuptools<81", "pandas<3"], check=True, timeout=900)
for archivo, viejo, nuevo in [("pyscenic/transform.py", "np.object,", "object,"),
                              ("pyscenic/rss.py", "dtype=np.float)", "dtype=float)")]:
    ruta = os.path.join(PKGS, archivo)
    with open(ruta) as f:
        texto = f.read()
    with open(ruta, "w") as f:
        f.write(texto.replace(viejo, nuevo))
os.environ["PYTHONPATH"] = PKGS
os.environ["PYTHONWARNINGS"] = "ignore"   # oculta avisos internos de librerías viejas
os.environ["PATH"] = os.path.join(PKGS, "bin") + os.pathsep + os.environ["PATH"]
subprocess.run(["pyscenic", "--help"], check=True, timeout=300, capture_output=True)
print("✅ pySCENIC 0.12.1 listo")
''')

# =============================================================================
md(r"""
## 1 · La matriz de expresión

- **`demo`:** se descarga la matriz del taller.
- **`completa`:** se descargan los datos originales de 10x (37 MB) y se aplica **el mismo control de calidad** que usó el instructor:
  - se suman los símbolos génicos duplicados;
  - células con ≥ 200 y < 4,000 genes detectados y < 15 % de lecturas mitocondriales;
  - genes con ≥ 3 × 1 % del número de células en conteos totales y detectados en ≥ 1 % de las células.

En `completa`, las poblaciones se asignan **solo con marcadores** sobre los clusters de 10x: 1, 4 y 12 → monocitos CD14+; 7 y 8 → células B; 9 → NK; el resto, "otras". El taller además verificó cada célula con CellTypist; aquí no.
""")

code(r'''
MAPA_MARCADORES = {1: "Mono_CD14", 4: "Mono_CD14", 12: "Mono_CD14", 7: "B", 8: "B", 9: "NK"}

if os.path.exists("expr.tsv") and os.path.exists("poblaciones.csv"):
    print("La matriz ya estaba preparada en esta sesión.")
elif MATRIZ == "demo":
    expr = pd.read_csv(conseguir("expression_demo.tsv.gz"), sep="\t", index_col="cell_id")
    meta = pd.read_csv(conseguir("metadata_demo.csv"), index_col="cell_id")
    expr.to_csv("expr.tsv", sep="\t")
    meta[["population"]].to_csv("poblaciones.csv")
else:
    try:
        import h5py, scipy.sparse as sp
    except ImportError:
        sin_pkgs = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}   # instalar en el Python del notebook
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", "h5py", "scipy"], check=True, env=sin_pkgs)
        import h5py, scipy.sparse as sp
    with h5py.File(conseguir("pbmc_10k_v3_filtered_feature_bc_matrix.h5"), "r") as f:
        m = f["matrix"]
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:]).T.tocsr()
        celulas = m["barcodes"][:].astype(str)
        simbolos = m["features"]["name"][:].astype(str)
    print("Original:", X.shape)
    unicos, inv = np.unique(simbolos, return_inverse=True)            # suma símbolos duplicados
    X = (X @ sp.csr_matrix((np.ones(len(simbolos)), (np.arange(len(simbolos)), inv)),
                           shape=(len(simbolos), len(unicos)))).tocsr()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    total = np.asarray(X.sum(axis=1)).ravel()
    mito = np.asarray(X[:, np.char.startswith(unicos.astype(str), "MT-")].sum(axis=1)).ravel()
    ok = (n_genes >= 200) & (n_genes < 4000) & (100 * mito / total < 15)
    X, celulas = X[ok], celulas[ok]
    gen_ok = np.asarray((X > 0).sum(axis=0)).ravel() >= 3
    X, unicos = X[:, gen_ok], unicos[gen_ok]
    print("Tras QC de células:", X.shape)
    n = X.shape[0]
    gen_ok = (np.asarray(X.sum(axis=0)).ravel() >= 3 * 0.01 * n) & (np.asarray((X > 0).sum(axis=0)).ravel() >= 0.01 * n)
    X, unicos = X[:, gen_ok], unicos[gen_ok]
    print("Tras filtrar genes:", X.shape)
    expr = pd.DataFrame(X.toarray().astype(np.int32), index=pd.Index(celulas, name="cell_id"), columns=unicos)
    expr.to_csv("expr.tsv", sep="\t")
    with tarfile.open(conseguir("pbmc_10k_v3_analysis.tar.gz")) as t:
        clusters = pd.read_csv(t.extractfile("analysis/clustering/graphclust/clusters.csv"), index_col=0)["Cluster"]
    pd.DataFrame({"population": clusters.reindex(expr.index).map(MAPA_MARCADORES).fillna("otras")},
                 index=expr.index).to_csv("poblaciones.csv")

poblaciones = pd.read_csv("poblaciones.csv", index_col=0)["population"]
print(f"✅ expr.tsv: {len(poblaciones):,} células |", poblaciones.value_counts().to_dict())
''')

# =============================================================================
md(r"""
## 2 · Recursos de aertslab (≈ 720 MB)

La lista de TF, la tabla motivo → TF y las **dos** bases de ranking. Cada archivo se comprueba con su SHA1.

**Si la descarga falla**, puede ser que el servidor de aertslab esté caído (a veces pasa). Tienes dos salidas:

- esperar y volver a ejecutar la celda;
- copiar los archivos a una carpeta de tu Drive y escribir su ruta en `RUTA_DRIVE_RECURSOS`.
""")

code(r'''
for archivo in ["allTFs_hg38.txt", "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
                "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]:
    conseguir(archivo)
print("✅ Recursos listos")
''')

# =============================================================================
md(r"""
## 3 · GRNBoost2 (`grn`)

Es el paso largo. Con 2 núcleos, `demo` tarda ~16 min y `completa` ~1.5–3 h (en una prueba con 5 procesos tardó 81 min). Usamos `arboreto_with_multiprocessing.py`, equivalente a `pyscenic grn` con los mismos argumentos, porque `pyscenic grn` falla con las versiones actuales de `dask`.
""")

code(r'''
%%time
n_genes = len(open("expr.tsv").readline().split("\t")) - 1
gb_por_proceso = 2.5 * len(poblaciones) * n_genes * 8 / 1e9        # cada proceso guarda su copia de la matriz
PROCESOS = max(1, min(NUCLEOS, int((RAM_GB - 2) // gb_por_proceso)))
print(f"GRNBoost2 con {PROCESOS} procesos ({len(poblaciones):,} células × {n_genes:,} genes, RAM {RAM_GB:.0f} GB)")
if not os.path.exists("resultados/adjacencies.tsv"):
    !arboreto_with_multiprocessing.py expr.tsv recursos/allTFs_hg38.txt --method grnboost2 --output resultados/adjacencies.tsv --num_workers {PROCESOS} --seed 2026
    guardar("resultados/adjacencies.tsv") if os.path.exists("resultados/adjacencies.tsv") else None
adj = pd.read_csv("resultados/adjacencies.tsv", sep="\t")
print(f"✅ {len(adj):,} asociaciones TF–gen de {adj.TF.nunique()} TF")
''')

# =============================================================================
md(r"""
## 4 · cisTarget (`ctx`)

Con las dos bases de ranking, igual que el servidor (~8 min con 2 núcleos).
""")

code(r'''
%%time
PROCESOS = max(1, min(NUCLEOS, int((RAM_GB - 1) // 3)))   # cada proceso carga las bases de ranking (~3 GB)
print(f"cisTarget con {PROCESOS} procesos")
if not os.path.exists("resultados/regulons.csv"):
    !pyscenic ctx resultados/adjacencies.tsv recursos/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather recursos/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname recursos/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname expr.tsv --mode custom_multiprocessing --num_workers {PROCESOS} --output resultados/regulons.csv
    guardar("resultados/regulons.csv") if os.path.exists("resultados/regulons.csv") else None

import ast
def leer_regulones_ctx(ruta):
    """Lee la salida de `pyscenic ctx` y devuelve {regulón: genes}, uniendo los motivos de cada TF."""
    tabla = pd.read_csv(ruta, header=[0, 1], index_col=[0, 1], skiprows=[2])
    regulones = {}
    for (tf, _motivo), dianas in tabla[("Enrichment", "TargetGenes")].items():
        dianas = re.sub(r"np\.float\d*\(([^)]*)\)", r"\1", dianas)   # numpy 2 escribe np.float64(0.5)
        regulones.setdefault(f"{tf}(+)", set()).update(gen for gen, _peso in ast.literal_eval(dianas))
    return regulones

regulones = leer_regulones_ctx("resultados/regulons.csv")
tamanos = pd.Series({r: len(g) for r, g in regulones.items()}).sort_values(ascending=False)
print(f"✅ {len(regulones)} regulones | tamaño mediano {tamanos.median():.0f} genes | más grande: "
      f"{tamanos.index[0]} ({tamanos.iloc[0]} genes)")
''')

# =============================================================================
md(r"""
## 5 · AUCell (`aucell`)
""")

code(r'''
%%time
if not os.path.exists("resultados/auc_mtx.csv"):
    !pyscenic aucell expr.tsv resultados/regulons.csv -o resultados/auc_mtx.csv --num_workers {NUCLEOS} --seed 2026
    guardar("resultados/auc_mtx.csv") if os.path.exists("resultados/auc_mtx.csv") else None
auc = pd.read_csv("resultados/auc_mtx.csv", index_col=0)
print(f"✅ Matriz AUCell: {auc.shape[0]:,} células × {auc.shape[1]} regulones")
''')

# =============================================================================
md(r"""
## 6 · Resultados

### Regulones más específicos de cada población

Diferencia de mediana de AUC frente a la siguiente población (la misma medida que la sección 5 del taller).
""")

code(r'''
ORDEN = ["Mono_CD14", "B", "NK"]
datos = auc.join(poblaciones)
med = datos[datos.population.isin(ORDEN)].groupby("population").median().T[ORDEN]
especificidad = pd.DataFrame({p: med[p] - med.drop(columns=p).max(axis=1) for p in ORDEN})
top = pd.DataFrame({p: especificidad[p].nlargest(8).index for p in ORDEN})
display(top)

elegidos = list(dict.fromkeys(top.to_numpy().T.ravel()))
z = med.loc[elegidos].sub(med.loc[elegidos].mean(axis=1), axis=0).div(med.loc[elegidos].std(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(5, 0.35 * len(elegidos) + 1.5))
sns.heatmap(z, cmap="RdBu_r", center=0, ax=ax, cbar_kws={"label": "z-score por regulón"})
ax.set(xlabel="", ylabel="", title=f"Mediana de AUC por población ({MATRIZ})")
plt.show()
''')

md(r"""
### Tu corrida frente a la del servidor

El servidor obtuvo **272 regulones** con la matriz `demo`, usando 20 núcleos y la imagen oficial de pySCENIC. GRNBoost2 tiene un componente aleatorio y las versiones de las librerías difieren, así que **no esperes una coincidencia exacta**. Lo importante es si los regulones con sentido biológico se repiten.
""")

code(r'''
servidor = pd.read_csv(conseguir("regulons_demo.csv")).groupby("regulon").target.apply(set).to_dict()
comunes = sorted(set(regulones) & set(servidor))
print(f"Tu corrida: {len(regulones)} regulones | servidor: {len(servidor)} | en común: {len(comunes)}")
print("Solo en tu corrida:", len(set(regulones) - set(servidor)), "| solo en el servidor:", len(set(servidor) - set(regulones)))

filas = []
for r in ["MAFB(+)", "CEBPB(+)", "SPI1(+)", "EBF1(+)", "SPIB(+)", "POU2AF1(+)", "TBX21(+)", "EOMES(+)", "FOS(+)"]:
    a, b = regulones.get(r, set()), servidor.get(r, set())
    filas.append({"regulón": r, "genes (tú)": len(a), "genes (servidor)": len(b),
                  "Jaccard": round(len(a & b) / len(a | b), 2) if a | b else np.nan})
display(pd.DataFrame(filas).set_index("regulón"))

if MATRIZ == "demo":   # mismas células: comparamos la actividad célula a célula
    auc_serv = pd.read_csv(conseguir("auc_demo.csv"), index_col="cell_id")
    r_celula = pd.Series({r: np.corrcoef(auc.loc[auc_serv.index, r], auc_serv[r])[0, 1]
                          for r in comunes if auc[r].std() > 0 and auc_serv[r].std() > 0})
    print(f"Correlación célula a célula del AUC (regulones en común): mediana {r_celula.median():.2f}; "
          f"{(r_celula > 0.8).mean():.0%} de los regulones con r > 0.8")
''')

md(r"""
### Para pensar

1. ¿Qué regulones aparecen en **ambas** corridas? ¿Coinciden con los que discutimos en el taller (MAFB, EBF1, TBX21…)?
2. Con `completa`, ¿cambian los regulones más específicos al usar ~10,000 células en lugar de 900? ¿Por qué podría pasar? *(Más células = más poder estadístico, pero también otras poblaciones —células T, plaquetas…— que aportan coexpresión.)*
3. ¿Qué harías para decidir qué regulones son **robustos**? *(Pista: repetir `grn` con otras semillas y quedarte con los regulones recurrentes.)*

**Créditos:** datos *10k PBMCs from a Healthy Donor (v3 chemistry)*, 10x Genomics (CC BY 4.0). pySCENIC: Aibar *et al.* 2017; Van de Sande *et al.* 2020. Recursos cisTarget: <https://resources.aertslab.org/cistarget/>.
""")

nb = {
    "cells": cells,
    "metadata": {
        "colab": {"provenance": [], "toc_visible": True},
        "kernelspec": {"display_name": "Python 3", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 0,
}
for c in nb["cells"]:
    c["source"] = c["source"].splitlines(keepends=True)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print(f"Escrito {OUT} ({len(cells)} celdas)")
