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
COLAB_URL = BASE_URL.replace("https://raw.githubusercontent.com/", "https://colab.research.google.com/github/").replace(
    "/main", "/blob/main") + "/notebooks/Taller_pySCENIC_2h.ipynb"
md(rf"""
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({COLAB_URL})

# 🧬 Taller pySCENIC: ¿qué factores de transcripción distinguen a monocitos, células B y células NK?

**Duración:** 2 horas · **Requisitos:** cuenta de Google + navegador · **No necesitas instalar nada.**

👉 **Cómo usar este notebook:** ejecuta las celdas **en orden** con el botón ▶ (o `Shift + Enter`). Las celdas marcadas con ✏️ tienen una variable que puedes cambiar.

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
- **Ya calculado por el instructor** con pySCENIC 0.12.1 (en un servidor, ~7 min): **900 células** (300 por población) × **10,120 genes**.
- **Etiquetas de población:** el dataset no trae tipos celulares; las asignamos con **genes marcadores canónicos** y las verificamos con **CellTypist**. Solo usamos células donde ambos métodos coinciden.

## Mapa del notebook

| Sección | Qué harás | Minutos aprox. |
|---|---|---|
| 0 · Preparación | Descargar datos y ver las poblaciones | 15–25 |
| 1 · GRNBoost2 | Explorar asociaciones TF–gen | 25–45 |
| 2 · cisTarget | Ver cómo se podan los regulones | 45–65 |
| 3 · AUCell | Comparar actividad entre poblaciones | 65–90 |
| 4 · *(Opcional)* AUCell en vivo | Calcular AUCell tú mismo | si da tiempo |
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

# =============================================================================
md(r"""
---
## 1 · Ejercicio 1 — GRNBoost2: asociaciones TF–gen

**Qué hizo GRNBoost2:** para **cada gen**, entrenó un modelo (árboles de decisión) que intenta **predecir la expresión de ese gen a partir de la expresión de los ~1,000 TF** de la lista. La **importancia** (`importance`) dice cuánto ayudó cada TF a predecir ese gen.

- Importancia alta ⇒ el TF y el gen están **fuertemente asociados** (coexpresados) en estas células.
- No dice **dirección** (¿quién regula a quién?) ni **causalidad**. Es una medida **relativa**, no un p-valor.

La tabla completa tiene **2.3 millones** de asociaciones; para el taller incluimos **todas las asociaciones de 13 TF** seleccionados.
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

El cálculo ya está hecho; aquí comparamos el *antes* y el *después*.
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
## 4 · *(Opcional)* AUCell en vivo

Hasta ahora usamos resultados precalculados. Aquí **calculas AUCell tú mismo** con la misma matriz y los mismos regulones, y comprobamos que obtienes **lo mismo**.

- **4.1 AUCell desde cero** (solo `pandas`/`numpy`, ~1 s): siempre funciona.
- **4.2 pySCENIC real** (instala el paquete oficial, ~1 min): puede fallar según la versión de Colab; si falla, no pasa nada.
- **4.3 Comparación** de lo que hayas podido calcular contra el resultado precalculado.

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
print(f"Diferencia máxima con el resultado precalculado de pySCENIC: {dif:.1e}",
      "✅ (idéntico)" if dif < 1e-9 else "⚠️ (distinto; avisa al instructor)")
''')

md(r"""
### 4.2 · *(Opcional)* pySCENIC real

La siguiente celda **escribe un pequeño script** que usa el paquete oficial `pyscenic`. La celda de después lo **instala** (en una carpeta aparte, para no afectar al resto del notebook) y lo **ejecuta**. Tarda ~1 minuto. Si falla, simplemente continúa: ya tienes el resultado de 4.1.
""")

code(r'''
%%writefile aucell_pyscenic.py
# AUCell con el paquete oficial pySCENIC (mismas opciones que el CLI: sin pesos, umbral 0.05, semilla 2026)
import sys, pandas as pd
from ctxcore.genesig import Regulon
from pyscenic.aucell import aucell

carpeta = sys.argv[1]
expr = pd.read_csv(f"{carpeta}/data/expression_demo.tsv.gz", sep="\t", index_col="cell_id")
reg = pd.read_csv(f"{carpeta}/results/regulons_demo.csv")
firmas = [Regulon(name=n, gene2weight={g: 1.0 for g in d.target}, transcription_factor=d.TF.iloc[0],
                  gene2occurrence=[]) for n, d in reg.groupby("regulon")]
auc = aucell(expr, firmas, auc_threshold=0.05, noweights=True, seed=2026, num_workers=2)
auc.to_csv("auc_pyscenic_vivo.csv")
print("pySCENIC AUCell OK:", auc.shape)
''')

code(r'''
import subprocess
auc_vivo = None
try:
    if not os.path.exists("pyscenic_pkgs/pyscenic"):
        print("Instalando pySCENIC 0.12.1 en ./pyscenic_pkgs (≈1 min)…")
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", "--target", "pyscenic_pkgs",
                        "pyscenic==0.12.1", "setuptools<81"], check=True, timeout=420,
                       capture_output=True, text=True)
    t0 = time.time()
    salida = subprocess.run([sys.executable, "aucell_pyscenic.py", CARPETA], check=True, timeout=300,
                            capture_output=True, text=True,
                            env={**os.environ, "PYTHONPATH": os.path.abspath("pyscenic_pkgs")})
    print(salida.stdout.strip(), f"({time.time() - t0:.0f} s)")
    auc_vivo = pd.read_csv("auc_pyscenic_vivo.csv", index_col=0)
except subprocess.TimeoutExpired:
    print("⏱️ Tardó demasiado. No pasa nada: continúa con el resultado de 4.1.")
except subprocess.CalledProcessError as e:
    print("⚠️ pySCENIC no funcionó en este Colab. No pasa nada: continúa con el resultado de 4.1.")
    import re
    lineas = [l for l in re.sub(r"\x1b\[[0-9;]*m", "", e.stderr or "").splitlines() if l.strip()]
    print("   Últimas líneas del error:\n   " + "\n   ".join(lineas[-2:]))
''')

md(r"""
### 4.3 · Comparación
""")

code(r'''
auc_pre = cargar("auc")
fuentes = {"precalculado (servidor)": auc_pre}
if "auc_cero" in globals():
    fuentes["desde cero (4.1)"] = auc_cero
if globals().get("auc_vivo") is not None:
    fuentes["pySCENIC en vivo (4.2)"] = auc_vivo

for nombre, tabla in fuentes.items():
    d = (tabla.loc[auc_pre.index, auc_pre.columns] - auc_pre).abs().to_numpy().max()
    print(f"{nombre:<26} diferencia máxima vs precalculado = {d:.1e}")

RG = "MAFB(+)"
otros = [n for n in fuentes if n != "precalculado (servidor)"]
if otros:
    fig, axs = plt.subplots(1, len(otros), figsize=(4.2 * len(otros), 4), squeeze=False)
    for ax, n in zip(axs[0], otros):
        ax.scatter(auc_pre[RG], fuentes[n].loc[auc_pre.index, RG], s=5,
                   c=[COLOR[p] for p in meta.loc[auc_pre.index, "population"]])
        lim = [0, auc_pre[RG].max() * 1.05]
        ax.plot(lim, lim, "k--", lw=0.8)
        ax.set(xlabel=f"{RG} precalculado", ylabel=f"{RG} {n}", title="Cada punto = una célula")
    plt.tight_layout()
    plt.show()
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
- **Variabilidad:** GRNBoost2 es estocástico; los regulones pueden cambiar algo entre ejecuciones (por eso SCENIC recomienda repetir varias veces).

### Cómo se calcularon estos datos (no lo ejecutes aquí)

```bash
# 1) GRNBoost2                     (~2 min con 20 núcleos)
pyscenic grn expression_demo.tsv allTFs_hg38.txt -o adjacencies.tsv --num_workers 20 --seed 2026
# 2) cisTarget                     (~5 min; bases de ranking hg38 v10, ~620 MB)
pyscenic ctx adjacencies.tsv hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname expression_demo.tsv --mode custom_multiprocessing --output regulons.csv
# 3) AUCell                        (~5 s)
pyscenic aucell expression_demo.tsv regulons.csv -o auc_mtx.csv --seed 2026
```

### Para seguir aprendiendo

- Documentación de pySCENIC: <https://pyscenic.readthedocs.io/>
- Protocolo SCENIC (Van de Sande *et al.*, 2020, *Nature Protocols*) y notebooks: <https://github.com/aertslab/SCENICprotocol>
- Recursos (listas de TF, motivos, bases de ranking): <https://resources.aertslab.org/cistarget/>
- Datos originales: <https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3>
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
for c in nb["cells"]:   # formato estándar: lista de líneas
    c["source"] = c["source"].splitlines(keepends=True)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print(f"Escrito {OUT} ({len(cells)} celdas)")
