"""Genera notebooks/Taller_pySCENIC_2h.ipynb (pasos 28–36).

Uso: python3 scripts/build_notebook.py
El notebook es la fuente que verán los alumnos; este script es la fuente que editamos nosotros.
"""
import json
import textwrap
from pathlib import Path

BASE_URL = "https://raw.githubusercontent.com/Ssdzem/taller-pyscenic-intro/main"  # confirmar en el paso 37/39
OUT = Path(__file__).resolve().parents[1] / "notebooks" / "Taller_pySCENIC_2h.ipynb"

cells = []


def md(src):
    cells.append({"cell_type": "markdown", "metadata": {}, "source": textwrap.dedent(src).strip("\n")})


def code(src):
    cells.append({"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
                  "source": textwrap.dedent(src).strip("\n")})


# =============================================================================
COLAB_CASA = BASE_URL.replace("https://raw.githubusercontent.com/", "https://colab.research.google.com/github/").replace(
    "/main", "/blob/main") + "/notebooks/Pipeline_completo_casa.ipynb"
COLAB_URL = BASE_URL.replace("https://raw.githubusercontent.com/", "https://colab.research.google.com/github/").replace(
    "/main", "/blob/main") + "/notebooks/Taller_pySCENIC_2h.ipynb"
md(rf"""
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({COLAB_URL})

# 🧬 Taller pySCENIC: ¿qué factores de transcripción distinguen a monocitos, células B y células NK?

**Duración:** 2 horas · **Requisitos:** cuenta de Google + navegador · **No necesitas instalar nada.**

👉 **Cómo usar este notebook:** ejecuta las celdas **en orden** con el botón ▶ (o `Shift + Enter`). Las celdas marcadas con ✏️ tienen una variable que puedes cambiar.

En cada ejercicio encontrarás dos tipos de bloque nuevos:

- 🖥️ **El comando:** el comando real de pySCENIC y una tabla que explica cada opción (*flag*), como en la documentación oficial.
- ▶ **Córrelo tú:** ejecutas ese comando aquí mismo, en una versión reducida que cabe en Colab.

---

## La pregunta biológica

> **¿Qué regulones muestran distinta actividad entre monocitos CD14+, células B y células NK de sangre humana?**

## Conceptos en un minuto

- **Factor de transcripción (TF):** proteína que se une al ADN cerca de ciertos genes y activa o reprime su expresión.
- **Red de regulación génica (GRN):** mapa de *qué TF regulan a qué genes*.
- **Regulón:** un TF **+** el conjunto de genes diana que probablemente regula de forma directa.

## pySCENIC en tres pasos

```
   Matriz de expresión (células × genes)
               │
   ① GRNBoost2 ── ¿qué genes se coexpresan con cada TF?            →  asociaciones TF–gen
               │
   ② cisTarget ── ¿esos genes tienen el motivo de unión del TF?     →  regulones (podados)
               │
   ③ AUCell ───── ¿qué tan activo está cada regulón en cada célula?  →  matriz de actividad (AUC)
```

> ⚠️ **Importante:** coexpresión + presencia de un motivo = **hipótesis de regulación**, **no prueba causal**.
> Para demostrar que un TF regula un gen hacen falta experimentos (ChIP-seq, CRISPR, perturbaciones…).

## Los datos

- **Origen:** *10k PBMCs from a Healthy Donor (v3)*, 10x Genomics (dataset público). Un solo donante sano.
- **Ya calculado por el instructor** con pySCENIC 0.12.1 (en un servidor con 20 núcleos, ~7 min): **900 células** (300 por población) × **10,120 genes**. Los ejercicios interpretan **esos resultados**.
- **Además, tú correrás los tres comandos** en una versión reducida (Colab gratuito tiene solo 2 núcleos) y compararás tu resultado con el del servidor.
- **Etiquetas de población:** el dataset no trae tipos celulares; las asignamos con **genes marcadores canónicos** y las verificamos con **CellTypist**. Solo usamos células donde ambos métodos coinciden.

## Mapa del notebook

| Sección | Qué harás | Minutos aprox. |
|---|---|---|
| 0 · Preparación | Descargar datos, instalar pySCENIC y ver las poblaciones | 15–25 |
| 1 · GRNBoost2 | Correr `grn` (≈ 2 min) y explorar asociaciones TF–gen | 25–45 |
| 2 · cisTarget | Correr `ctx` (≈ 5 min) y ver cómo se podan los regulones | 45–65 |
| 3 · AUCell | Correr `aucell` (≈ 10 s) y comparar actividad entre poblaciones | 65–90 |
| 4 · *(Opcional)* AUCell por dentro | Programar AUCell desde cero | si da tiempo |
| 5 · Interpretación en parejas | Proponer una hipótesis | 90–110 |
| 6 · Limitaciones y recursos | Qué se puede y qué no se puede concluir | 110–120 |

💡 **Si una celda falla, no te preocupes:** cada sección vuelve a cargar lo que necesita. Avísale al instructor y continúa con la siguiente.
""")

# =============================================================================
md(r"""
## 0 · Preparación

Esta celda importa las librerías y define una función `cargar()` que **descarga automáticamente** los archivos del taller desde GitHub (≈ 10 MB). No tienes que subir nada.
""")

code(rf'''
# ▶ Ejecuta esta celda (tarda unos segundos).
import os, shutil, sys, time, urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

BASE_URL = "{BASE_URL}"
CARPETA = "datos_taller"

# nombre corto -> (ruta en el repositorio, opciones de lectura)
ARCHIVOS = {{
    "metadatos":   ("data/metadata_demo.csv",              dict(index_col="cell_id")),
    "expresion":   ("data/expression_demo.tsv.gz",         dict(sep="\t", index_col="cell_id")),
    "adyacencias": ("results/adjacencies_demo.tsv",        dict(sep="\t")),
    "resumen_tf":  ("results/adjacencies_tf_summary.csv",  dict(index_col="TF")),
    "regulones":   ("results/regulons_demo.csv",           dict()),
    "motivos":     ("results/motifs_summary.csv",          dict()),
    "auc":         ("results/auc_demo.csv",                dict(index_col="cell_id")),
}}

def descargar(ruta):
    """Descarga un archivo del repositorio (una sola vez) y devuelve su ruta local."""
    destino = os.path.join(CARPETA, ruta)
    if os.path.exists(destino) and os.path.getsize(destino) > 0:
        return destino
    os.makedirs(os.path.dirname(destino), exist_ok=True)
    local = os.environ.get("TALLER_LOCAL")          # solo para pruebas del instructor
    if local and os.path.exists(os.path.join(local, ruta)):
        shutil.copy(os.path.join(local, ruta), destino)
        return destino
    for intento in range(3):
        try:
            urllib.request.urlretrieve(f"{{BASE_URL}}/{{ruta}}", destino)
            return destino
        except Exception as e:
            ultimo_error = e
            time.sleep(2)
    raise RuntimeError(f"No se pudo descargar {{ruta}} ({{ultimo_error}}). "
                       "Revisa tu conexión o pide al instructor el ZIP de respaldo.")

def cargar(nombre):
    """Devuelve una tabla (DataFrame) del taller, p. ej. cargar('metadatos')."""
    ruta, opciones = ARCHIVOS[nombre]
    return pd.read_csv(descargar(ruta), **opciones)

# Estilo de gráficos y colores fijos por población
sns.set_theme(style="white", context="notebook")
ORDEN = ["Mono_CD14", "B", "NK"]
NOMBRE = {{"Mono_CD14": "Monocitos CD14+", "B": "Células B", "NK": "Células NK"}}
COLOR = {{"Mono_CD14": "#d95f02", "B": "#1b9e77", "NK": "#7570b3"}}

for nombre, (ruta, _) in ARCHIVOS.items():
    p = descargar(ruta)
    print(f"✅ {{nombre:<12}} {{ruta:<38}} {{os.path.getsize(p)/1e6:5.1f}} MB")
print("\nPython", sys.version.split()[0], "| pandas", pd.__version__, "| numpy", np.__version__)
''')

md(r"""
### Las tres poblaciones

Cada fila de `metadatos` es **una célula**. Fíjate en la columna `population` y en las coordenadas `TSNE1`, `TSNE2` (una proyección 2D calculada por 10x Genomics donde células parecidas quedan cerca).
""")

code(r'''
meta = cargar("metadatos")
display(meta.head())
print("Células por población:")
print(meta["population"].value_counts())

fig, ax = plt.subplots(figsize=(5, 4.5))
for p in ORDEN:
    d = meta[meta.population == p]
    ax.scatter(d.TSNE1, d.TSNE2, s=6, c=COLOR[p], label=NOMBRE[p])
ax.set(xticks=[], yticks=[], xlabel="t-SNE 1", ylabel="t-SNE 2", title="900 células del taller")
ax.legend(markerscale=3, frameon=False)
plt.show()
''')

code(r'''
# ✅ Comprobación de la sección 0
assert meta.shape[0] == 900, "Deberían ser 900 células"
assert meta.index.is_unique, "Los IDs de célula deben ser únicos"
assert set(meta.population) == set(ORDEN)
print("✅ Sección 0 lista:", meta.shape[0], "células,", meta.population.nunique(), "poblaciones.")
''')

md(r"""
### Instalar pySCENIC (≈ 30 s)

En los ejercicios vas a **correr los comandos reales de pySCENIC**. Esta celda lo instala en una carpeta aparte (`pyscenic_pkgs`) para no alterar el resto de Colab. Ejecútala y sigue leyendo mientras termina.

> 🧩 **Detalle técnico:** pySCENIC 0.12.1 es de 2022 y choca con librerías más nuevas que trae Colab. La celda aplica tres arreglos:
> - instala `setuptools<81` (sin él falla `import ctxcore`);
> - instala `pandas<3` (pandas 3 rompe `ctx`);
> - cambia `np.object`/`np.float` por `object`/`float` en dos archivos, porque numpy 2 eliminó esos nombres.
>
> Si la instalación falla, no pasa nada: todos los ejercicios funcionan con los resultados del servidor.
""")

code(r'''
# ▶ Instala pySCENIC 0.12.1 (≈ 30 s)
import re, subprocess
PKGS = os.path.abspath("pyscenic_pkgs")

def instalar_pyscenic():
    if not os.path.exists(os.path.join(PKGS, "pyscenic")):
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", "--target", PKGS,
                        "pyscenic==0.12.1", "setuptools<81", "pandas<3"],
                       check=True, timeout=600, capture_output=True, text=True)
    for archivo, viejo, nuevo in [("pyscenic/transform.py", "np.object,", "object,"),
                                  ("pyscenic/rss.py", "dtype=np.float)", "dtype=float)")]:
        ruta = os.path.join(PKGS, archivo)
        with open(ruta) as f:
            texto = f.read()
        with open(ruta, "w") as f:
            f.write(texto.replace(viejo, nuevo))
    # Los comandos "!" de las celdas siguientes usarán esta instalación
    os.environ["PYTHONPATH"] = PKGS
    os.environ["PYTHONWARNINGS"] = "ignore"   # oculta avisos internos de librerías viejas
    os.environ["PATH"] = os.path.join(PKGS, "bin") + os.pathsep + os.environ["PATH"]
    subprocess.run(["pyscenic", "--help"], check=True, timeout=300, capture_output=True, text=True)

t0 = time.time()
try:
    instalar_pyscenic()
    PYSCENIC_LISTO = True
    print(f"✅ pySCENIC 0.12.1 listo en {time.time() - t0:.0f} s. Comandos: pyscenic, arboreto_with_multiprocessing.py")
except Exception as e:
    PYSCENIC_LISTO = False
    detalle = re.sub(r"\x1b\[[0-9;]*m", "", str(getattr(e, "stderr", None) or e)).strip().splitlines()
    print("⚠️ No se pudo instalar pySCENIC. No pasa nada: usarás los resultados del servidor.")
    print("   Detalle:", detalle[-1] if detalle else type(e).__name__)
NUCLEOS = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()   # núcleos asignados a esta sesión
print("Núcleos de CPU disponibles:", NUCLEOS)
''')

# =============================================================================
md(r"""
---
## 1 · Ejercicio 1 — GRNBoost2: asociaciones TF–gen

**Qué hizo GRNBoost2:** para **cada gen**, entrenó un modelo (árboles de decisión) que intenta **predecir la expresión de ese gen a partir de la expresión de los ~1,000 TF** de la lista. La **importancia** (`importance`) dice cuánto ayudó cada TF a predecir ese gen.

- Importancia alta ⇒ el TF y el gen están **fuertemente asociados** (coexpresados) en estas células.
- No dice **dirección** (¿quién regula a quién?) ni **causalidad**. Es una medida **relativa**, no un p-valor.

La tabla completa del servidor tiene **2.3 millones** de asociaciones (75 MB): demasiado para descargarla en clase. Para explorarla incluimos **todas las asociaciones, sin filtrar, de 13 TF** elegidos porque cuentan la historia:

- reguladores clásicos de cada población: SPI1, CEBPB, CEBPA e IRF8 (monocitos); PAX5, EBF1, SPIB y POU2AF1 (células B); TBX21, EOMES, RUNX3 y PRDM1 (NK);
- un caso trampa: RPS4X.

Los **272 regulones** y la **matriz AUCell** sí van completos.

### 🖥️ El comando: `pyscenic grn`

Así lo corrió el instructor en el servidor:

```bash
pyscenic grn expression_demo.tsv allTFs_hg38.txt -o adjacencies.tsv --num_workers 20 --seed 2026
```

| Parte | Qué es | Valor usado | Por qué |
|---|---|---|---|
| `expression_demo.tsv` | 1.er argumento: matriz de expresión. CSV/TSV con **filas = células y columnas = genes**, o un archivo `.loom` | 900 células × 10,120 genes, conteos crudos | pySCENIC usa expresión sin escalar por gen |
| `allTFs_hg38.txt` | 2.º argumento: lista de TF candidatos, uno por línea | 1,892 TF humanos de [aertslab](https://resources.aertslab.org/cistarget/tf_lists/); 1,018 están en la matriz | Solo los genes de esta lista pueden actuar como "reguladores" |
| `-o`, `--output` | Archivo de salida: tabla `TF`, `target`, `importance` | `adjacencies.tsv` | — |
| `--num_workers` | Núcleos de CPU que trabajan en paralelo | 20 en el servidor; 2 en Colab | Más núcleos = más rápido |
| `--seed` | Semilla aleatoria de los árboles de decisión | 2026 | Para poder repetir el resultado; sin semilla, cada corrida da algo distinto |
| `-m`, `--method` | Algoritmo: `grnboost2` (por defecto) o `genie3` | `grnboost2` | Mucho más rápido que GENIE3, con resultados parecidos |
| `-t`, `--transpose` | Úsalo si tu tabla tiene **genes en filas** | no | La nuestra ya tiene células en filas |
| `--sparse` | Carga la matriz como dispersa para ahorrar memoria | no | Nuestra matriz es pequeña |

> 🧩 **Detalle técnico:** con versiones recientes de la librería `dask`, `pyscenic grn` falla con *"Must supply at least one delayed object"*. pySCENIC incluye el script `arboreto_with_multiprocessing.py`, que hace **lo mismo con los mismos argumentos** sin usar `dask`. Es el que usaremos en Colab.

### ▶ Córrelo tú (≈ 2 min)

Colab gratuito tiene **2 núcleos**, frente a los 20 del servidor. Con la matriz completa GRNBoost2 tardaría ~16 min, así que usamos una **versión reducida**: las mismas 900 células, pero solo los **2,000 genes más variables** más los **1,018 TF** (2,829 genes en total).

La primera celda prepara los archivos y la segunda corre el comando. **Mientras corre, repasen la tabla de arriba.**
""")

code(r'''
# ▶ Prepara los archivos de entrada del comando
expr = cargar("expresion")
tfs = cargar("resumen_tf").index.tolist()     # los 1,018 TF de allTFs_hg38.txt presentes en la matriz
lognorm = np.log1p(expr.div(expr.sum(axis=1), axis=0) * 1e4)
variables = set(lognorm.var().nlargest(2000).index)
genes = [g for g in expr.columns if g in variables or g in set(tfs)]
expr[genes].to_csv("expr_reducida.tsv", sep="\t")
pd.Series(tfs).to_csv("tfs.txt", index=False, header=False)
print(f"expr_reducida.tsv: {len(expr)} células × {len(genes)} genes  |  tfs.txt: {len(tfs)} TF")
''')

code(r'''
%%time
# ▶ GRNBoost2 de verdad. Equivale a:  pyscenic grn expr_reducida.tsv tfs.txt -o mi_adyacencias.tsv --num_workers 2 --seed 2026
!arboreto_with_multiprocessing.py expr_reducida.tsv tfs.txt --method grnboost2 --output mi_adyacencias.tsv --num_workers {NUCLEOS} --seed 2026
''')

code(r'''
# ¿Qué obtuviste? Tu tabla frente a la del servidor
if os.path.exists("mi_adyacencias.tsv") and os.path.getsize("mi_adyacencias.tsv") > 0:
    mia = pd.read_csv("mi_adyacencias.tsv", sep="\t")
    adj = cargar("adyacencias")
    print(f"Tu tabla: {len(mia):,} asociaciones de {mia.TF.nunique()} TF "
          f"(servidor: 2,339,328 asociaciones con 10,120 genes)")
    display(mia.head())
    tf = "TBX21"
    tuyo = mia[mia.TF == tf].nlargest(15, "importance").target.tolist()
    servidor = adj[adj.TF == tf].nlargest(15, "importance").target.tolist()
    print(f"Top 15 de {tf}, tu corrida: {', '.join(tuyo)}")
    print(f"Top 15 de {tf}, servidor:   {', '.join(servidor)}")
    print(f"En común: {len(set(tuyo) & set(servidor))} de 15")
else:
    print("No se generó mi_adyacencias.tsv (¿falló la instalación?). No pasa nada: seguimos con la tabla del servidor.")
''')

md(r"""
### Explorar las asociaciones del servidor

Ahora sí, a interpretar los resultados completos del servidor.
""")

code(r'''
adj = cargar("adyacencias")        # TF, target (gen diana), importance (peso)
resumen_tf = cargar("resumen_tf")  # nº de dianas candidatas por TF en la tabla COMPLETA
print("Filas en la tabla del taller:", len(adj))
print("TF incluidos:", sorted(adj.TF.unique()))
display(adj.head())
''')

code(r'''
# ✏️ Cambia el TF si quieres explorar otro de la lista de arriba
TF = "TBX21"

dianas = adj[adj.TF == TF].sort_values("importance", ascending=False)
print(f"{TF}: {len(dianas)} genes diana candidatos")
display(dianas.head(15))

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.hist(dianas.importance, bins=60, color="grey")
ax.set_yscale("log")
ax.set(xlabel="importancia (peso de la arista)", ylabel="nº de genes (escala log)",
       title=f"Distribución de pesos de {TF}: pocos genes fuertes, muchos débiles")
plt.show()
''')

md(r"""
**Preguntas — Ejercicio 1**

1. ¿Cuántos genes diana **candidatos** tiene TBX21? ¿Te parece realista que un TF regule *directamente* tantos genes?
2. Mira el top 15: ¿reconoces genes típicos de **células NK** (pista: `NKG7`, `GNLY`, `PRF1`, `GZMA`)?
3. ¿Qué representa el **peso** de la arista? ¿Significa que TBX21 *activa* esos genes?
4. Compara **tu** top 15 de TBX21 (de la celda ▶) con el del servidor. Usamos la misma semilla: ¿por qué no coinciden del todo? *(Pista: ¿qué genes había en tu matriz?)*

✍️ *Tus respuestas:*
""")

md(r"""
### Un caso curioso: RPS4X

La lista oficial de TF de pySCENIC incluye a **RPS4X** (una proteína ribosomal). Veamos sus asociaciones y qué TF tienen más dianas en la tabla completa.
""")

code(r'''
display(adj[adj.TF == "RPS4X"].sort_values("importance", ascending=False).head(10))
print("Asociación más fuerte de TODO el taller:")
display(adj.sort_values("importance", ascending=False).head(3))
print("TF con más dianas candidatas en la tabla completa (2.3 M filas):")
display(resumen_tf.head(10))
''')

md(r"""
**Pregunta:** RPS4X tiene las asociaciones más fuertes de toda la red, casi todas con otros genes ribosomales (`RPL…`, `RPS…`). ¿Crees que es un regulador? ¿Qué otra explicación hay? *(Lo comprobaremos en el Ejercicio 2.)*

✍️ *Tu respuesta:*
""")

code(r'''
# ✅ Comprobación del Ejercicio 1
assert {"TF", "target", "importance"} <= set(adj.columns)
assert len(adj[adj.TF == "TBX21"]) > 0
print(f"✅ Ejercicio 1 listo: {len(adj):,} asociaciones de {adj.TF.nunique()} TF;",
      f"TBX21 tiene {len(adj[adj.TF == 'TBX21']):,} candidatos.")
''')

# =============================================================================
md(r"""
---
## 2 · Ejercicio 2 — cisTarget: de asociaciones a regulones

GRNBoost2 encuentra genes **coexpresados** con un TF, pero muchos pueden ser **efectos indirectos**. cisTarget añade una segunda evidencia, independiente de la expresión:

1. Agrupa las dianas de cada TF en **módulos** (p. ej. "el top 50", "las que tienen peso > percentil 90").
2. Para cada módulo pregunta: **¿las regiones reguladoras de estos genes están enriquecidas en el motivo de unión del TF?** (se usan bases precalculadas que ordenan ~27,000 genes según su puntuación para cada motivo, en ventanas de ±10 kb y −500/+100 pb alrededor del inicio de transcripción).
3. Si el motivo está enriquecido (**NES ≥ 3**), conserva **solo los genes que tienen el motivo** → **regulón**.

Primero lo corres tú con **tus** asociaciones; después comparamos el *antes* y el *después* con los resultados del servidor.

### 🖥️ El comando: `pyscenic ctx`

```bash
pyscenic ctx adjacencies.tsv \
    hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname expression_demo.tsv \
    --mode custom_multiprocessing --num_workers 20 --output regulons.csv
```

*(La `\` al final de una línea solo indica que el comando continúa en la siguiente.)*

| Parte | Qué es | Valor usado | Por qué |
|---|---|---|---|
| `adjacencies.tsv` | 1.er argumento: las asociaciones que produjo `grn` | — | De aquí salen los módulos de cada TF |
| `*.rankings.feather` | Una o varias **bases de ranking**: para cada motivo, los ~27,000 genes humanos ordenados según qué tan fuerte aparece el motivo cerca de su inicio de transcripción | ±10 kb y −500/+100 pb, genoma hg38 (~311 MB cada una) | Dos ventanas: una región amplia y el promotor cercano |
| `--annotations_fname` | Tabla **motivo → TF** | `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` (99 MB) | Dice qué TF puede unirse a cada motivo |
| `--expression_mtx_fname` | La misma matriz de expresión de `grn` | `expression_demo.tsv` | Para calcular correlaciones y armar los módulos |
| `--mode` | Cómo repartir el trabajo entre núcleos | `custom_multiprocessing` | El más sencillo; no necesita `dask` |
| `--num_workers` | Núcleos de CPU | 20 en el servidor; 2 en Colab | — |
| `-o`, `--output` | Tabla de motivos enriquecidos con sus genes diana | `regulons.csv` | Cada TF con algún motivo enriquecido se convierte en un regulón |
| `--nes_threshold` | NES mínimo para aceptar un motivo como enriquecido | 3.0 (por defecto) | Más alto = menos regulones, pero más confiables |
| `--auc_threshold` | Fracción del ranking usada para calcular el enriquecimiento | 0.05 (por defecto) | — |
| `--rank_threshold` | Hasta qué posición del ranking se buscan dianas con el motivo | 5000 (por defecto) | — |
| `--thresholds`, `--top_n_targets`, `--top_n_regulators` | Cómo se arman los módulos: dianas por encima de los percentiles 75 y 90 de peso; las 50 mejores dianas de cada TF; los 5, 10 y 50 mejores TF de cada gen | valores por defecto | Varias formas de agrupar dan más oportunidades de encontrar el motivo |
| `--min_genes` | Tamaño mínimo de un módulo | 20 (por defecto) | Los módulos pequeños son poco fiables |
| `-a`, `--all_modules` | Incluir también regulones represores `(-)` | no (por defecto) | Por eso todos nuestros regulones son `(+)`, activadores |
| `--mask_dropouts` | Ignorar las células con cero al calcular correlaciones | no (por defecto) | Comportamiento por defecto desde pySCENIC 0.9.16 |

### ▶ Córrelo tú (≈ 1 min de descarga + 4 min de cálculo)

Para ahorrar tiempo usamos **una sola base** (−500/+100 pb) con **tus** asociaciones de la sección 1. Primero hay que conseguir la base (≈ 313 MB) y la tabla de motivos (≈ 99 MB). Hay dos maneras:

- **Opción 1 (por defecto): descarga directa** desde los servidores de aertslab. La descarga la hace Colab en la nube de Google, así que **no pasa por el wifi del salón**.
- **Opción 2: Google Drive**, si el instructor compartió una carpeta con los dos archivos (o si la descarga directa falla):
  1. Abre el enlace de la carpeta, clic derecho → *Organizar* → *Agregar acceso directo* → *Mi unidad*.
  2. Escribe su ruta en `RUTA_DRIVE`, por ejemplo `"/content/drive/MyDrive/taller_pyscenic_recursos"`.
  3. Ejecuta la celda. Colab te pedirá permiso para ver tu Drive.
""")

code(r'''
# ▶ Consigue la base de ranking y la tabla de motivos (≈ 412 MB)
import hashlib
RUTA_DRIVE = ""   # ✏️ Opción 2: ruta a la carpeta compartida en tu Drive (déjala vacía para descargar)

AERTS = "https://resources.aertslab.org/cistarget"
RECURSOS = {   # archivo -> (carpeta en aertslab, SHA1 para comprobar que llegó completo)
    "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather":
        (f"{AERTS}/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based",
         "d9d4e927b80b2b35344aec0a9fc0a82ae9d76759"),
    "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl":
        (f"{AERTS}/motif2tf", "8cba044ef45b389f4b85237931531906b3046182"),
}

def sha1(ruta):
    h = hashlib.sha1()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()

def conseguir(archivo):
    carpeta, esperado = RECURSOS[archivo]
    destino = os.path.join("recursos", archivo)
    os.makedirs("recursos", exist_ok=True)
    if os.path.exists(destino) and sha1(destino) == esperado:
        return destino
    if RUTA_DRIVE:
        if not os.path.exists("/content/drive"):
            from google.colab import drive
            drive.mount("/content/drive")
        shutil.copy(os.path.join(RUTA_DRIVE, archivo), destino)
    elif os.environ.get("TALLER_RECURSOS"):            # solo para pruebas del instructor
        shutil.copy(os.path.join(os.environ["TALLER_RECURSOS"], archivo), destino)
    else:
        urllib.request.urlretrieve(f"{carpeta}/{archivo}", destino)
    if sha1(destino) != esperado:
        os.remove(destino)
        raise RuntimeError(f"{archivo} llegó incompleto o dañado; vuelve a ejecutar la celda")
    return destino

try:
    for archivo in RECURSOS:
        t0 = time.time()
        ruta = conseguir(archivo)
        print(f"✅ {archivo}  ({os.path.getsize(ruta) / 1e6:.0f} MB, {time.time() - t0:.0f} s)")
    RECURSOS_LISTOS = True
except Exception as e:
    RECURSOS_LISTOS = False
    print(f"⚠️ No se pudieron conseguir los recursos: {e}")
    print("   Prueba la Opción 2 (Drive). Si tampoco funciona, no pasa nada: sigue con los regulones del servidor.")
''')

code(r'''
%%time
# ▶ cisTarget de verdad (≈ 4 min con 2 núcleos)
if RECURSOS_LISTOS and os.path.exists("mi_adyacencias.tsv"):
    !pyscenic ctx mi_adyacencias.tsv recursos/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname recursos/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname expr_reducida.tsv --mode custom_multiprocessing --num_workers {NUCLEOS} --output mi_regulones.csv
else:
    print("Falta mi_adyacencias.tsv (sección 1) o la base de ranking. Seguimos con los regulones del servidor.")
''')

code(r'''
# ¿Cuántos regulones obtuviste? Tus regulones frente a los del servidor
import ast

def leer_regulones_ctx(ruta):
    """Lee la salida de `pyscenic ctx` y devuelve {regulón: genes}, uniendo los motivos de cada TF."""
    tabla = pd.read_csv(ruta, header=[0, 1], index_col=[0, 1], skiprows=[2])
    regulones = {}
    for (tf, _motivo), dianas in tabla[("Enrichment", "TargetGenes")].items():
        dianas = re.sub(r"np\.float\d*\(([^)]*)\)", r"\1", dianas)   # numpy 2 escribe np.float64(0.5)
        regulones.setdefault(f"{tf}(+)", set()).update(gen for gen, _peso in ast.literal_eval(dianas))
    return regulones

if os.path.exists("mi_regulones.csv") and os.path.getsize("mi_regulones.csv") > 0:
    mios = leer_regulones_ctx("mi_regulones.csv")
    servidor = cargar("regulones").groupby("regulon").target.apply(set).to_dict()
    print(f"Tus regulones: {len(mios)}  |  servidor: {len(servidor)}  |  en común: {len(set(mios) & set(servidor))}")
    display(pd.DataFrame([{"regulón": r, "genes (tu corrida)": len(mios.get(r, set())),
                           "genes (servidor)": len(servidor[r]),
                           "genes en común": len(mios.get(r, set()) & servidor[r])}
                          for r in ["MAFB(+)", "EBF1(+)", "TBX21(+)", "FOS(+)", "SPI1(+)"]]).set_index("regulón"))
else:
    print("No hay mi_regulones.csv. No pasa nada: seguimos con los regulones del servidor.")
''')

md(r"""
### Los regulones del servidor
""")

code(r'''
reg = cargar("regulones")      # regulon, TF, target, weight
motivos = cargar("motivos")    # TF, MotifID, NES, AUC, n_targets, database
tamanos = reg.groupby("regulon").size().sort_values(ascending=False)
print(f"Regulones que sobrevivieron a cisTarget: {len(tamanos)}")
print(f"Tamaño (nº de genes): mediana {tamanos.median():.0f}, mínimo {tamanos.min()}, máximo {tamanos.max()}")
display(reg.head())

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.hist(np.log10(tamanos), bins=30, color="grey")
ax.set(xlabel="log10(nº de genes en el regulón)", ylabel="nº de regulones", title="Tamaño de los regulones")
plt.show()
''')

code(r'''
# ✏️ Compara candidatos (GRNBoost2) vs. regulón (cisTarget) para un TF
TF = "TBX21"

candidatos = set(adj.loc[adj.TF == TF, "target"])
regulon = set(reg.loc[reg.TF == TF, "target"])
print(f"{TF}: {len(candidatos):,} candidatos en GRNBoost2 → {len(regulon):,} genes en el regulón "
      f"({100 * len(regulon) / max(len(candidatos), 1):.1f}% conservados)")
print(f"¿Todos los genes del regulón venían de los candidatos? {regulon <= candidatos}")

display(reg[reg.TF == TF].head(10))
print(f"Motivos enriquecidos que justifican el regulón {TF}(+):")
display(motivos[motivos.TF == TF].sort_values("NES", ascending=False).head(5))
''')

code(r'''
# ✏️ Dibuja la red de un regulón (TF en el centro, dianas alrededor)
REGULON = "MAFB(+)"
N_MAX = 40                        # para regulones grandes, solo las N dianas de mayor peso

r = reg[reg.regulon == REGULON].sort_values("weight", ascending=False)
tf = r.TF.iloc[0]
dianas_red = [g for g in r.target if g != tf][:N_MAX]
G = nx.DiGraph()
G.add_edges_from((tf, g, {"weight": w}) for g, w in zip(r.target, r.weight) if g in dianas_red)
angulos = np.linspace(0, 2 * np.pi, len(dianas_red), endpoint=False)
pos = {tf: (0, 0), **{g: (np.cos(a), np.sin(a)) for g, a in zip(dianas_red, angulos)}}
pesos = np.array([G[u][v]["weight"] for u, v in G.edges()])

fig, ax = plt.subplots(figsize=(7, 7))
nx.draw_networkx_edges(G, pos, ax=ax, width=0.5 + 3 * (pesos - pesos.min()) / (np.ptp(pesos) or 1),
                       edge_color="#aaaaaa", arrowsize=8)
nx.draw_networkx_nodes(G, pos, nodelist=[tf], node_color="#d95f02", node_size=1500, ax=ax)
nx.draw_networkx_nodes(G, pos, nodelist=dianas_red, node_color="#fde0c5", node_size=450, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
ax.set_title(f"{REGULON}: {len(r)} genes (se muestran {len(dianas_red)}); grosor = importancia GRNBoost2")
ax.axis("off")
plt.show()
''')

code(r'''
# ¿Quién sobrevivió y quién no? (los 13 TF del taller)
tf_taller = adj.TF.unique()
tabla = pd.DataFrame({
    "candidatos_GRNBoost2": adj.groupby("TF").size(),
    "genes_en_regulon": reg.groupby("TF").size().reindex(tf_taller),
}).fillna(0).astype(int).sort_values("candidatos_GRNBoost2", ascending=False)
tabla["¿tiene regulón?"] = np.where(tabla.genes_en_regulon > 0, "sí", "NO")
display(tabla)
print("Candidatos de PAX5 en GRNBoost2 (top 10):",
      adj[adj.TF == "PAX5"].nlargest(10, "importance").target.tolist())
''')

md(r"""
**Preguntas — Ejercicio 2**

1. ¿Qué fracción de los candidatos de TBX21 sobrevivió? ¿Qué criterio usó cisTarget para descartar el resto?
2. **RPS4X** tenía las asociaciones más fuertes de toda la red. ¿Por qué no tiene regulón? ¿Qué nos dice sobre usar solo coexpresión?
3. **PAX5** es un TF esencial de células B y sus candidatos son genes B "de libro" (`IGHM`, `CD79B`, `CD40`…), pero **no tiene regulón**. ¿Significa que PAX5 no regula esos genes? *(Pistas: calidad/cobertura de los motivos, umbral de NES, ventanas genómicas usadas, elementos reguladores lejanos.)*
4. **Tu corrida** dio un número de regulones distinto de los **272** del servidor. ¿Por qué? Piensa en:
   - cuántos genes usaste (2,829 frente a 10,120);
   - cuántas bases de ranking (1 frente a 2);
   - que GRNBoost2 tiene un componente aleatorio.

   ¿Qué implica para confiar en un regulón concreto? *(Por eso SCENIC recomienda repetir el análisis varias veces y quedarse con los regulones que aparecen siempre.)*

✍️ *Tus respuestas:*
""")

code(r'''
# ✅ Comprobación del Ejercicio 2
assert {"regulon", "TF", "target", "weight"} <= set(reg.columns)
assert "MAFB(+)" in set(reg.regulon) and "RPS4X(+)" not in set(reg.regulon)
print(f"✅ Ejercicio 2 listo: {reg.regulon.nunique()} regulones; MAFB(+) tiene "
      f"{(reg.regulon == 'MAFB(+)').sum()} genes; RPS4X y PAX5 no tienen regulón.")
''')

# =============================================================================
md(r"""
---
## 3 · Ejercicio 3 — AUCell: actividad de regulones por célula

**Qué hace AUCell:** en **cada célula**, ordena todos los genes de más a menos expresado y pregunta: *¿los genes del regulón aparecen concentrados arriba del ranking?*

- El resultado es un **AUC entre 0 y 1** por célula y regulón: alto ⇒ regulón **activo** en esa célula.
- Como usa **rankings dentro de cada célula**, no depende de cómo se normalice la expresión.
- Compara **el mismo regulón entre células/poblaciones**; los valores absolutos de regulones distintos no son directamente comparables.

> 🔑 **Regla de oro:** une tablas **por ID de célula**, nunca suponiendo que las filas están en el mismo orden.

### 🖥️ El comando: `pyscenic aucell`

```bash
pyscenic aucell expression_demo.tsv regulons.csv -o auc_mtx.csv --num_workers 20 --seed 2026
```

| Parte | Qué es | Valor usado | Por qué |
|---|---|---|---|
| `expression_demo.tsv` | 1.er argumento: matriz de expresión (células × genes) | la completa, 10,120 genes | El ranking de cada célula usa **todos** los genes |
| `regulons.csv` | 2.º argumento: los regulones. Acepta la salida de `ctx` (`.csv`) o archivos `.gmt`, `.yaml` y `.dat` | los 272 regulones | — |
| `-o`, `--output` | Matriz de actividad (células × regulones), en `.csv` o `.loom` | `auc_mtx.csv` | — |
| `--auc_threshold` | Fracción superior del ranking que se mira | 0.05 (por defecto): el top 5 %, ~505 genes | Solo importa si las dianas están **arriba** del ranking |
| `--seed` | Semilla para desempatar genes con la misma expresión (muy común: muchos ceros) | 2026 | Sin semilla, los valores cambian un poco en cada corrida |
| `-w`, `--weights` | Usar los pesos de los genes del regulón | no (por defecto) | Todos los genes del regulón cuentan igual |
| `--num_workers` | Núcleos de CPU | 20 en el servidor; 2 en Colab | — |
| `-t`, `--transpose` | Si la matriz tiene genes en filas | no | — |

### ▶ Córrelo tú (≈ 10 s)

Aquí usamos la **matriz completa** y los **272 regulones del servidor**, guardados en formato `.gmt`: una línea por regulón con su nombre, una descripción y sus genes. Con la misma semilla deberías obtener **exactamente** la matriz del servidor.
""")

code(r'''
# ▶ Prepara las entradas: matriz completa (sin comprimir) y regulones del servidor en formato .gmt
expr = cargar("expresion")
expr.to_csv("expression_demo.tsv", sep="\t")
reg = cargar("regulones")
with open("regulones_servidor.gmt", "w") as f:
    for nombre, d in reg.groupby("regulon"):
        f.write("\t".join([nombre, "regulon pySCENIC"] + d.target.tolist()) + "\n")
print(f"expression_demo.tsv: {expr.shape[0]} células × {expr.shape[1]} genes  |  "
      f"regulones_servidor.gmt: {reg.regulon.nunique()} regulones")
''')

code(r'''
%%time
# ▶ AUCell de verdad
!pyscenic aucell expression_demo.tsv regulones_servidor.gmt -o mi_auc.csv --num_workers {NUCLEOS} --seed 2026
''')

code(r'''
# ¿Obtuviste lo mismo que el servidor?
if os.path.exists("mi_auc.csv") and os.path.getsize("mi_auc.csv") > 0:
    mi_auc = pd.read_csv("mi_auc.csv", index_col=0)
    auc_servidor = cargar("auc")
    dif = (mi_auc.loc[auc_servidor.index, auc_servidor.columns] - auc_servidor).abs().to_numpy().max()
    print(f"Tu matriz AUCell: {mi_auc.shape[0]} células × {mi_auc.shape[1]} regulones")
    print(f"Diferencia máxima con el servidor: {dif:.1e}",
          "✅ idéntico" if dif < 1e-9 else "⚠️ distinto (avisa al instructor)")
else:
    print("No se generó mi_auc.csv. No pasa nada: seguimos con la matriz del servidor.")
''')

md(r"""
**(Opcional, si corriste cisTarget en la sección 2)** Calcula la actividad de **tus propios** regulones y compárala con la de los del servidor:
""")

code(r'''
# ▶ (Opcional) AUCell con TUS regulones
if os.path.exists("mi_regulones.csv") and os.path.getsize("mi_regulones.csv") > 0:
    !pyscenic aucell expression_demo.tsv mi_regulones.csv -o mi_auc_propio.csv --num_workers {NUCLEOS} --seed 2026
if os.path.exists("mi_auc_propio.csv"):
    propio = pd.read_csv("mi_auc_propio.csv", index_col=0)
    servidor = cargar("auc")
    pob = cargar("metadatos")["population"]
    comparar = [r for r in ["MAFB(+)", "EBF1(+)", "TBX21(+)", "FOS(+)", "SPI1(+)"] if r in propio.columns]
    tabla = pd.concat({"tu corrida": propio[comparar].groupby(pob).median().loc[ORDEN].T,
                       "servidor": servidor[comparar].groupby(pob).median().loc[ORDEN].T}, axis=1)
    print("Mediana de AUC por población (solo regulones que aparecen en ambas corridas):")
    display(tabla.round(3))
else:
    print("No hay regulones propios (no corriste cisTarget). No pasa nada.")
''')

md(r"""
### Actividad de los regulones del servidor
""")

code(r'''
auc = cargar("auc")            # filas = células, columnas = regulones
meta = cargar("metadatos")
print("Matriz AUCell:", auc.shape, "(células × regulones)")

datos = meta[["population"]].join(auc, how="inner")   # unión POR ID de célula
assert len(datos) == len(meta) == len(auc)
display(datos.iloc[:5, :6])
''')

code(r'''
# ✏️ Regulones a comparar (puedes añadir o quitar; usa el nombre exacto con "(+)")
REGULONES = ["MAFB(+)", "EBF1(+)", "TBX21(+)", "FOS(+)", "SPI1(+)"]

mediana = datos.groupby("population")[REGULONES].median().loc[ORDEN].T
z = mediana.sub(mediana.mean(axis=1), axis=0).div(mediana.std(axis=1), axis=0)   # escala por regulón

fig, ax = plt.subplots(figsize=(6, 0.6 * len(REGULONES) + 1.2))
sns.heatmap(z, annot=mediana.round(3), fmt="", cmap="RdBu_r", center=0, ax=ax,
            xticklabels=[NOMBRE[p] for p in ORDEN], cbar_kws={"label": "z-score por regulón"})
ax.set(xlabel="", ylabel="", title="Mediana de AUC por población (número = mediana real)")
plt.show()
''')

code(r'''
# ✏️ Distribución por célula para dos regulones
PAR = ["MAFB(+)", "TBX21(+)"]

fig, axs = plt.subplots(1, len(PAR), figsize=(5 * len(PAR), 3.8))
for ax, rg in zip(np.atleast_1d(axs), PAR):
    sns.boxplot(data=datos, x="population", y=rg, order=ORDEN, hue="population", palette=COLOR,
                legend=False, ax=ax, fliersize=1)
    ax.set(xlabel="", ylabel="AUC", title=rg, xticks=range(3), xticklabels=[NOMBRE[p] for p in ORDEN])
plt.tight_layout()
plt.show()
''')

md(r"""
### ¿Actividad del regulón = expresión del TF?

Ahora comparamos, para los mismos TF, su **expresión** (normalizada: log de conteos por 10,000) con la **actividad** de su regulón.
""")

code(r'''
expr = cargar("expresion")                    # conteos crudos, filas = células, columnas = genes
TFS = [rg.replace("(+)", "") for rg in REGULONES]
expr_norm = np.log1p(expr[TFS].div(meta.loc[expr.index, "n_counts"], axis=0) * 1e4)

fig, axs = plt.subplots(2, len(TFS), figsize=(3.2 * len(TFS), 6), sharex=True)
for j, tf in enumerate(TFS):
    for i, (tabla, etiqueta) in enumerate([(expr_norm[tf], f"expresión {tf}"), (datos[f"{tf}(+)"], f"AUC {tf}(+)")]):
        ax = axs[i, j]
        sns.boxplot(x=meta.loc[tabla.index, "population"], y=tabla, order=ORDEN, hue=meta.loc[tabla.index, "population"],
                    palette=COLOR, legend=False, ax=ax, fliersize=1)
        ax.set(xlabel="", ylabel="", title=etiqueta, xticks=range(3), xticklabels=["Mono", "B", "NK"])
plt.tight_layout()
plt.show()
''')

md(r"""
**Preguntas — Ejercicio 3**

1. ¿En qué población es más activo cada regulón? ¿Coincide con lo que sabes de **MAFB** (mieloides), **EBF1** (linaje B) y **TBX21/T-bet** (NK y T citotóxicas)?
2. Compara **FOS**: el TF se expresa en las tres poblaciones, pero ¿y su regulón? ¿Actividad de regulón y expresión del TF son lo mismo? ¿Por qué podrían diferir?
3. **SPI1** es un "regulador maestro" mieloide, pero su regulón tiene **>3,000 genes**. ¿Cómo afecta el tamaño del regulón a la diferencia de AUC entre poblaciones?
4. *(Para pensar)* FOS, JUNB y FOSB son **genes de respuesta temprana** que pueden subir por el estrés de procesar la muestra. ¿Cómo cambiaría eso tu interpretación del regulón FOS(+)?

✍️ *Tus respuestas:*
""")

code(r'''
# ✅ Comprobación del Ejercicio 3
assert auc.shape == (900, auc.shape[1]) and auc.index.is_unique
assert set(auc.index) == set(meta.index), "AUCell y metadatos deben tener las mismas células"
assert auc.values.min() >= 0 and auc.values.max() <= 1
print(f"✅ Ejercicio 3 listo: AUCell {auc.shape[0]} células × {auc.shape[1]} regulones, unido por ID.")
''')

# =============================================================================
md(r"""
---
## 4 · *(Opcional)* AUCell por dentro

En la sección 3 corriste `pyscenic aucell` como una caja negra. Aquí lo **programas desde cero** en unas pocas líneas de `pandas`/`numpy` y compruebas que da **exactamente** lo mismo.

### 4.1 · AUCell desde cero

Para una célula y un regulón:

1. Ordena los genes de la célula de **más a menos expresado** (rango 0 = el más expresado; los empates se desempatan al azar con una semilla fija).
2. Mira solo el **top 5 %** del ranking (aquí, los primeros 505 genes).
3. **Curva de recuperación:** al ir bajando por el ranking, ¿cuántos genes del regulón has encontrado? Si aparecen pronto, la curva sube rápido.
4. **AUC** = área bajo esa curva ÷ área máxima posible.
""")

code(r'''
expr = cargar("expresion")
reg = cargar("regulones")
meta = cargar("metadatos")
regulones = reg.groupby("regulon")["target"].apply(list).to_dict()

UMBRAL_AUC = 0.05     # fracción del ranking que se mira (valor por defecto de pySCENIC)
SEMILLA = 2026        # misma semilla que usó el instructor

# Paso 1: ranking por célula (0 = gen más expresado); se barajan columnas antes para desempatar al azar
rankings = (expr.sample(frac=1.0, axis=1, random_state=SEMILLA)
                .rank(axis=1, ascending=False, method="first") - 1)[expr.columns]
corte = int(round(UMBRAL_AUC * expr.shape[1])) - 1   # igual que pySCENIC
print("Rankings:", rankings.shape, "| se mira el top", corte, "de", expr.shape[1], "genes")

# Curva de recuperación de MAFB(+) en un monocito y en una célula B
REGULON = "MAFB(+)"
genes = [g for g in regulones[REGULON] if g in expr.columns]
fig, ax = plt.subplots(figsize=(6, 4))
for p in ["Mono_CD14", "B"]:
    celula = meta.index[meta.population == p][0]
    r = np.sort(rankings.loc[celula, genes].values)
    x = np.arange(corte + 1)
    y = np.searchsorted(r, x, side="right")          # nº de genes del regulón encontrados hasta el rango x
    auc_manual = np.clip(corte - r, 0, None).sum() / ((corte + 1) * len(genes))
    ax.step(x, y, where="post", color=COLOR[p], label=f"{NOMBRE[p]}  (AUC = {auc_manual:.3f})")
    ax.fill_between(x, y, step="post", color=COLOR[p], alpha=0.15)
ax.set(xlabel="posición en el ranking de la célula", ylabel=f"genes de {REGULON} recuperados",
       title=f"Curva de recuperación de {REGULON} ({len(genes)} genes)")
ax.legend(frameon=False)
plt.show()
''')

code(r'''
def aucell_desde_cero(expr, regulones, umbral=0.05, semilla=2026):
    """AUCell en pocas líneas (réplica de pySCENIC 0.12.1 con sus opciones por defecto)."""
    rnk = (expr.sample(frac=1.0, axis=1, random_state=semilla)
               .rank(axis=1, ascending=False, method="first") - 1)[expr.columns].to_numpy()
    corte = int(round(umbral * expr.shape[1])) - 1
    columna = {g: i for i, g in enumerate(expr.columns)}
    resultado = {}
    for nombre, genes in regulones.items():
        idx = [columna[g] for g in genes if g in columna]
        if len(idx) < 0.8 * len(genes):                       # pySCENIC exige ≥80% de los genes presentes
            resultado[nombre] = np.zeros(len(expr)); continue
        # cada gen del regulón dentro del top suma (corte − rango); se divide por el máximo posible
        resultado[nombre] = np.clip(corte - rnk[:, idx], 0, None).sum(axis=1) / ((corte + 1) * len(idx))
    return pd.DataFrame(resultado, index=expr.index)

t0 = time.time()
auc_cero = aucell_desde_cero(expr, regulones)
print(f"AUCell desde cero: {auc_cero.shape} en {time.time() - t0:.1f} s")
auc_pre = cargar("auc")
dif = (auc_cero[auc_pre.columns].loc[auc_pre.index] - auc_pre).abs().to_numpy().max()
print(f"Diferencia máxima con el resultado del servidor: {dif:.1e}",
      "✅ (idéntico)" if dif < 1e-9 else "⚠️ (distinto; avisa al instructor)")
if os.path.exists("mi_auc.csv"):
    mi_auc = pd.read_csv("mi_auc.csv", index_col=0)
    dif2 = (auc_cero[mi_auc.columns].loc[mi_auc.index] - mi_auc).abs().to_numpy().max()
    print(f"Diferencia máxima con tu `pyscenic aucell` de la sección 3: {dif2:.1e}")
''')

# =============================================================================
md(r"""
---
## 5 · Interpretación en parejas (20 min)

Elijan **un regulón** (pueden usar la lista de abajo como inspiración), generen su "ficha" y completen la plantilla.
""")

code(r'''
auc = cargar("auc"); meta = cargar("metadatos"); reg = cargar("regulones")
motivos = cargar("motivos"); expr = cargar("expresion")
datos = meta[["population"]].join(auc, how="inner")
med = datos.groupby("population").median().T[ORDEN]
especificidad = pd.DataFrame({p: med[p] - med.drop(columns=p).max(axis=1) for p in ORDEN})
print("Regulones más específicos de cada población (diferencia de mediana de AUC vs. la siguiente población):")
display(pd.DataFrame({NOMBRE[p]: especificidad[p].nlargest(8).index for p in ORDEN}))
''')

code(r'''
# ✏️ Escribe aquí el regulón que eligieron
MI_REGULON = "EBF1(+)"

def ficha(regulon):
    tf = regulon.split("(")[0]
    r = reg[reg.regulon == regulon].sort_values("weight", ascending=False)
    print(f"■ {regulon}: {len(r)} genes diana. Top 15: {', '.join(r.target.head(15))}")
    m = motivos[motivos.TF == tf].sort_values("NES", ascending=False).head(1)
    if len(m):
        print(f"■ Motivo con mayor enriquecimiento: {m.MotifID.iloc[0]} (NES = {m.NES.iloc[0]:.2f})")
    tabla = pd.DataFrame({"mediana AUC": datos.groupby("population")[regulon].median()})
    if tf in expr.columns:
        tabla[f"expresión media {tf} (log CP10k)"] = np.log1p(
            expr[tf] / meta.loc[expr.index, "n_counts"] * 1e4).groupby(meta.population).mean()
    display(tabla.loc[ORDEN].rename(index=NOMBRE).round(3))

ficha(MI_REGULON)
''')

md(r"""
### Plantilla (complétenla aquí mismo con doble clic)

| | |
|---|---|
| **Regulón elegido** | |
| **Población con mayor actividad** (mediana AUC) | |
| **¿El TF se expresa en esa misma población?** | |
| **Hipótesis biológica** (1–2 frases) | |
| **Evidencia observada** (medida directamente en los datos) | |
| **Evidencia inferida** (deducida por el algoritmo) | |
| **Una limitación metodológica** | |
| **¿Qué experimento pondría a prueba la hipótesis?** | |
""")

# =============================================================================
md(r"""
---
## 6 · Limitaciones y recursos

**Qué SÍ podemos decir:** que ciertos regulones (TF + dianas con motivo) muestran **actividad diferencial** entre poblaciones, de forma coherente con la biología conocida.

**Qué NO podemos decir sin más evidencia:**

- **Causalidad:** coexpresión + motivo ≠ regulación demostrada.
- **Ausencia:** que un TF no tenga regulón (p. ej. PAX5) **no** significa que no regule nada; las bases de motivos y las ventanas genómicas son incompletas.
- **Generalización:** **un solo donante**, **900 células** submuestreadas y etiquetas de población **inferidas** por nosotros.
- **Artefactos técnicos:** genes de estrés/respuesta temprana (FOS, JUN), número de genes detectados por célula y *dropout* pueden influir.
- **Variabilidad:** GRNBoost2 es estocástico y el resultado depende de los genes y las bases que uses. Lo viste en la sección 2: tu corrida no dio los mismos regulones que el servidor. Por eso SCENIC recomienda repetir el análisis varias veces.

### Resumen: el pipeline completo en tres comandos

Esto es lo que corrió el instructor en el servidor (20 núcleos). En tu computadora o servidor usarías los mismos comandos:

```bash
# 1) GRNBoost2   (~2 min con 20 núcleos)
pyscenic grn expression_demo.tsv allTFs_hg38.txt -o adjacencies.tsv --num_workers 20 --seed 2026
# 2) cisTarget   (~5 min; dos bases de ranking hg38 v10, ~620 MB)
pyscenic ctx adjacencies.tsv hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname expression_demo.tsv --mode custom_multiprocessing --num_workers 20 --output regulons.csv
# 3) AUCell      (~5 s)
pyscenic aucell expression_demo.tsv regulons.csv -o auc_mtx.csv --num_workers 20 --seed 2026
```

- **Dónde se corrió:** en la imagen oficial de contenedor `aertslab/pyscenic:0.12.1`, que trae las versiones de librerías con las que pySCENIC funciona sin arreglos.
- **En Colab o con un Python reciente** hacen falta los arreglos de la sección 0, y `arboreto_with_multiprocessing.py` en lugar de `pyscenic grn`.

### 🏠 Para casa: el pipeline completo, sin límite de tiempo

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](COLAB_CASA)

El notebook **`Pipeline_completo_casa.ipynb`** corre los tres pasos completos en Colab. Puedes elegir entre dos matrices:

- **demo:** las mismas 900 células y 10,120 genes de hoy. Tarda ~30 min y debería reproducir los resultados del servidor.
- **completa:** el dataset entero de 10x, 10,280 células tras el control de calidad. Tarda ~2–3 h.

### Para seguir aprendiendo

- Documentación de pySCENIC: <https://pyscenic.readthedocs.io/>
- Protocolo SCENIC (Van de Sande *et al.*, 2020, *Nature Protocols*) y notebooks: <https://github.com/aertslab/SCENICprotocol>
- Recursos (listas de TF, motivos, bases de ranking): <https://resources.aertslab.org/cistarget/>
- Datos originales: <https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3>
""".replace("COLAB_CASA", COLAB_CASA))

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
for c in nb["cells"]:   # formato estándar: lista de líneas
    c["source"] = c["source"].splitlines(keepends=True)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print(f"Escrito {OUT} ({len(cells)} celdas)")
