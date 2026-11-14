# STATUS — Taller pySCENIC

Máquina de trabajo ("Stargazer"): `cachyos-x8664`, usuario `emilio-work`, carpeta `~/repos/taller-pyscenic`.

## Decisiones

| Fecha | Decisión | Motivo |
|---|---|---|
| 2026-11-11 | Contenedor oficial con **Podman rootless** (no Docker) | Docker exigía meter a `emilio-work` al grupo `docker` (≈ root). Podman sin root no da privilegios extra. |
| 2026-11-11 | Notebook = **Opción A garantizada + Opción B opcional** | A (resultados precalculados) siempre funciona; B (AUCell en vivo) en 3 niveles: pySCENIC real → AUCell reimplementado en numpy (validado contra pySCENIC) → CSV precalculado. |

## Entorno

- Imagen: `docker.io/aertslab/pyscenic:0.12.1`, digest `sha256:0c06b8b0a00117e1a4b61303e9ad775bf53c9e410d4c344151c15bf0d143a288`, 1.01 GB.
- Dentro: Python 3.10.6, pyscenic 0.12.1, ctxcore 0.2.0, arboreto 0.1.6, numpy 1.23.5, pandas 1.5.1, dask/distributed 2022.11.1, pyarrow 10.0.0, loompy 3.0.7 (`logs/06_versions.log`).
- Podman 6.2.0 rootless.

## Punto de control A — ✅ COMPLETADO (2026-11-11)

| Paso | Estado | Evidencia |
|---|---|---|
| 1. Docs oficiales | ✅ | Subcomandos: `grn`, `add_cor`, `ctx`, `aucell` (`logs/06_help.log`) |
| 2. SCENICprotocol quick start | ✅ | Su ejemplo usa Nextflow + `pyscenic:0.9.19`, no 0.12.1 (ver Problemas) |
| 3. Carpetas | ✅ | `data/ resources/ results/ notebooks/ figures/ scripts/ logs/` |
| 4. Motor de contenedores | ✅ | Podman 6.2.0 rootless |
| 5. Pull de la imagen | ✅ | `logs/05_pull.log` |
| 6. `pyscenic -h` | ✅ | `logs/06_help.log` |
| 7. Mini-dataset oficial | ✅ | `test_oficial/`, sha256 en `logs/07_download_test.log` |
| 8. Ensayo grn → ctx → aucell | ✅ | `scripts/run_test_oficial.sh`, logs `logs/20_*`, `22_*`, `24_*` |
| 9. Tres salidas legibles | ✅ | `results/test_oficial/` (abajo) |

### Salidas de la prueba técnica (mini-dataset: 100 células × 200 genes, 1 TF = BRF1)

- (a) `adjacencies.tsv`: 199 filas `TF, target, importance` (importance 25–77). ~4 s.
- (b) `regulons.csv`: 8 filas motivo-enriquecido (2 motivos únicos, varios módulos) → **1 regulón `BRF1(+)` con 57 genes**. ~5 s.
- (c) `auc_mtx.csv`: 100 células × 1 regulón, AUC 0–0.039. ~4 s.
- Coherencia comprobada: genes del regulón = unión de las dianas de los motivos ⊂ adyacencias; células de AUCell = células del loom (mismo orden, 0 duplicados, 0 NaN).
- ⚠️ Es una prueba técnica; **no tiene interpretación biológica**.

### Formatos aprendidos (útiles para el notebook)

- `adjacencies.tsv`: TSV con encabezado `TF  target  importance`.
- `regulons.csv`: CSV con **encabezado de 3 filas** (multi-índice `TF, MotifID` × columnas `Enrichment/AUC, NES, …, TargetGenes`). `TargetGenes` es una lista Python serializada como texto `[('GEN', peso), …]`. Un regulón = unión por TF (y por signo, `(+)`) de las dianas de todos sus motivos. Para clase conviene exportar una versión "larga" simple (`regulon, gene`).
- `auc_mtx.csv`: filas = células (`Cell`), columnas = regulones `TF(+)`.

## Problemas encontrados y resueltos

1. **Docker sin permiso** (`/var/run/docker.sock` en `root:root 0600`) → se cambió a Podman rootless (lo instaló admin).
2. **`ctx` falló: base Feather v1 no soportada.** `genome-ranking.feather` del ejemplo oficial está en formato v1; ctxcore 0.2.0 solo lee v2:
   `ValueError: ... Feather v1 format, which is not supported anymore. Convert them with "convert_cistarget_databases_v1_to_v2.py"`.
   Solución: script oficial `aertslab/create_cisTarget_databases` (commit `304d5dc`, revisado: solo lee/escribe Feather) →
   `test_oficial/v2/genome-ranking.genes_vs_motifs.rankings.feather` (`logs/08_convert_v1_v2.log`).
   **Lección para el paso 18:** bajar directamente bases **Feather v2** (`*.genes_vs_motifs.rankings.feather`, carpeta `refseq_r80/mc_v10_clust` o similar) con su tabla de motivos v10 correspondiente.
3. Advertencias `FutureWarning` de pandas en `ctx` (groupby/apply): inofensivas.

## Encargo 2 — datos PBMC10k (pasos 10–13 ✅, 14 pendiente de ti)

| Paso | Estado | Evidencia |
|---|---|---|
| 10. Notebooks oficiales PBMC10k | ✅ | `PBMC10k_SCENIC-protocol-CLI.ipynb`, `PBMC10k_downstream-analysis.ipynb` (SCENICprotocol, master) |
| 11. Fuente de datos | ✅ | 10x Genomics *10k PBMCs from a Healthy Donor (v3 chemistry)*, Cell Ranger 3.0.0, `filtered_feature_bc_matrix`. Página: https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3 |
| 12. Descarga | ✅ | `data/original/` (sha256 en `logs/12_download_pbmc10k.log`): `pbmc_10k_v3_filtered_feature_bc_matrix.h5` (37 MB), `pbmc_10k_v3_analysis.tar.gz` (33 MB), `pbmc_10k_v3_metrics_summary.csv` |
| 13. Etiquetas reales | ✅ | `scripts/13_inspect_clusters.py`, `logs/13_clusters_markers.log` |

**Matriz:** 11,769 células × 33,538 genes, conteos UMI enteros, símbolos HGNC (24 símbolos duplicados → resolver en el paso 15).

**Qué etiquetas existen de verdad (paso 13):**
- 10x publica solo **clusters sin nombre** (`graphclust`: 17 clusters; `kmeans_2…10`) y **t-SNE** + PCA. **No hay UMAP ni tipos celulares.**
- El notebook CLI oficial solo calcula Louvain sin nombres. El downstream usa `Celltype_Garnett` (CD4 T, CD8 T, T, NK, B, Monocytes, Dendritic cells, CD34+, Unknown) leído de un `PBMC10k.loom` integrado que **no está publicado** (ruta interna del clúster de Lovaina). No es recuperable.
- Cualquier nombre de tipo celular será por tanto **anotación nuestra por marcadores canónicos**, documentada como tal.

Evidencia de marcadores (media log1p CP10k) → lectura propuesta:

| graphclust | Marcadores | Lectura | n |
|---|---|---|---|
| 1, 4, 12 | CD14, LYZ | Monocitos CD14+ | 3,226 |
| 13 | FCGR3A, LYZ | Monocitos FCGR3A+ | 363 |
| 7, 8 | MS4A1, CD79A | Células B | 1,549 |
| 9 | GNLY, NKG7, FCGR3A | NK | 605 |
| 3 | CD3E, CD8A, NKG7 | T CD8/citotóxicas | 1,112 |
| 2, 5, 6, 11 | CD3E, IL7R (CCR7) | T CD4 | 3,423 |
| 15 | FCER1A, CLEC10A, LILRA4 | DC (cDC+pDC) | 290 |
| 17 | PPBP | Plaquetas | 185 |
| 10, 14, 16 | mixtos/bajos | Ambiguos (16 ≈ dobletes T+mono) | 1,016 |

### Pasos 14–16 ✅

| Paso | Estado | Evidencia |
|---|---|---|
| 14. Contraste (TÚ, 2026-11-11) | ✅ | **Monocitos CD14+ vs Células B vs NK.** Pregunta: *¿qué regulones muestran distinta actividad entre monocitos CD14+, células B y células NK?* Anotación: marcadores canónicos + verificación CellTypist, solo células concordantes. |
| 15. Matriz reducida | ✅ | `scripts/15_prepare_demo.py`, `logs/15_prepare_demo.log`, `logs/15_prepare_demo_report.json` |
| 16. Coherencia | ✅ | `scripts/16_check_demo.py`, `logs/16_check_demo.log` → TODO OK |

**Procedimiento del paso 15** (todo *antes* de inferir la red):
1. 11,769 × 33,538 → suma de 24 símbolos duplicados (ninguno es TF clave).
2. QC (umbrales del notebook oficial): ≥200 y <4000 genes/célula, <15% mito; genes en ≥3 células → **10,280 células × 19,911 genes** (1,489 células removidas).
3. Marcadores: clusters graphclust 1/4/12 → `Mono_CD14`, 7/8 → `B`, 9 → `NK`.
4. CellTypist 1.6.3, modelo `Immune_All_Low.pkl` v2 (2022-07-16), por célula, sin majority voting; imagen `localhost/taller-celltypist:1.6.3` (`scripts/containers/celltypist/Containerfile`).
5. Concordancia: Mono_CD14 2,864/2,963 (96.7%; discordantes sobre todo "Non-classical monocytes"), B 1,446/1,463 (98.8%), NK 537/589 (91.2%; discordantes sobre todo "Tem/Temra cytotoxic T", bien excluidas).
6. Submuestreo estratificado 300/población, `random_state=2026`.
7. Filtro de genes sobre el subconjunto (umbrales oficiales): ≥27 conteos totales y detectado en ≥9 células → **10,120 genes**.

**Salidas:**
- `data/expression_demo.tsv` — 900 células × 10,120 genes, conteos UMI crudos enteros, filas = células (18.4 MB; comprimir para GitHub).
- `data/metadata_demo.csv` — `cell_id, population, cluster_10x_graphclust, celltypist_label, celltypist_conf, n_genes, n_counts, pct_mito, TSNE1, TSNE2` (75 KB).
- Coordenadas: **t-SNE de 10x** (calculado por Cell Ranger sobre las 11,769 células; verificadas, no recalculadas). No hay UMAP publicado → la figura 3 será t-SNE (o UMAP calculado por nosotros, documentado).

**Sanity check biológico (log1p CP10k medio):** SPI1 1.48 / CEBPB 1.64 en Mono; PAX5 0.62 / EBF1 0.38 / SPIB 0.84 en B; TBX21 1.04 / EOMES 0.44 en NK; ≈0 en las demás poblaciones.

⚠️ Nota para interpretación: los monocitos detectan más genes por célula (mediana 2,819) que NK (2,023) y B (1,595). AUCell es por ranking dentro de cada célula, pero conviene mencionarlo como limitación.

## Encargo 2 — pipeline real (pasos 17–27)

| Paso | Estado | Evidencia |
|---|---|---|
| 17. Lista TF | ✅ | `resources/allTFs_hg38.txt` (1,892 TF) — https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt |
| 18. Bases + motivos (aprobado 2026-11-11) | ✅ | ver tabla abajo; SHA1 = oficiales; sha256 en `logs/18_resources_sha256.log` |
| 19. Compatibilidad | ✅ | `scripts/19_check_compat.py`, `logs/19_compat.log` |
| 20. `grn` | ✅ | `results/adjacencies.tsv`, `logs/20_grn.log` (1 min 50 s, 20 workers, `--seed 2026`) |
| 21. Revisión adyacencias | ✅ | 2,339,328 filas `TF, target, importance`; 1,018 TF; 10,120 dianas |
| 22. `ctx` | ✅ | `results/regulons.csv`, `logs/22_ctx.log` (4 min 30 s; 2 bases v10, opciones por defecto) |
| 23. Regulones sobrevivientes | ✅ | **272 regulones** (todos activadores `(+)`), tamaño mediana 41 genes (3–3,106); 2,372 filas motivo-enriquecido |
| 24. `aucell` | ✅ | `results/auc_mtx.csv` 900 células × 272 regulones (5 s, `--seed 2026`); células = matriz, mismo orden; 0 NaN |
| 25. Unión por ID + resumen | ✅ | `scripts/25_postprocess.py`, `logs/25_postprocess_report.json`, `logs/25_candidates.log` |
| 26. Figuras (AMBOS) | ✅ · ⏳ revisión | `scripts/26_figures.py` → `figures/red_regulon`, `heatmap`, `tsne` (.png + .pdf) |
| 27. Exportables ligeros | ✅ | ver lista abajo |

**Resultado de cisTarget para los TF didácticos (nº genes en el regulón):** SPI1 3,106 · CEBPB 1,543 · CEBPA 177 · IRF8 243 · SPIB 413 · TBX21 1,014 · EOMES 607 · RUNX3 501 · PRDM1 554 · POU2AF1 20 · EBF1 14 · **PAX5 — (no sobrevivió)** · **RPS4X — (podado, era coexpresión ribosomal)**.

**Regulones más específicos (RSS):** B: BACH2, EBF1, POU2AF1, TFEB, FOXO1 · Mono_CD14: NFE2, TFEC, FOS, MAFB, KLF10 · NK: STAT4, GFI1, HIC1, HOXB4, EOMES.

**Regulones propuestos para clase (pendiente de tu validación):**

| Regulón | Genes | AUC mediana B / Mono / NK | Expr. TF B / Mono / NK | Uso didáctico |
|---|---|---|---|---|
| MAFB(+) | 37 | 0.113 / **0.311** / 0.109 | 0.00 / 0.84 / 0.02 | Monocitos; red completa dibujable (CSF1R, CD68, CD36, FCGR2A…) |
| EBF1(+) | 14 | **0.049** / 0.000 / 0.000 | 0.38 / 0.01 / 0.01 | Linaje B |
| TBX21(+) | 1,014 | 0.028 / 0.024 / **0.075** | 0.01 / 0.00 / 1.04 | NK (T-bet) |
| FOS(+) | 40 | 0.192 / **0.694** / 0.217 | 2.12 / 4.33 / 2.03 | Expresión ≠ actividad; posible artefacto de disociación (genes de respuesta temprana) |
| SPI1(+) (opc.) | 3,106 | 0.044 / 0.060 / 0.042 | 0.42 / 1.48 / 0.04 | Regulón enorme → diferencia de AUC diluida |

**Exportables (paso 27):**
- `results/adjacencies_demo.tsv` (0.9 MB): **todas** las aristas GRNBoost2 de 13 TF didácticos (SPI1, CEBPB, CEBPA, IRF8, PAX5, EBF1, SPIB, POU2AF1, TBX21, EOMES, RUNX3, PRDM1, RPS4X) = 45,880 filas. Filtro solo por TF; importancias sin modificar (4 decimales).
- `results/adjacencies_tf_summary.csv`: nº de dianas candidatas por TF en la tabla completa (1,018 TF).
- `results/regulons_demo.csv` (1.8 MB): formato largo `regulon, TF, target, weight`, los 272 regulones completos.
- `results/motifs_summary.csv`: por TF/motivo → NES, AUC, nº de dianas, base usada.
- `results/auc_demo.csv` (3.3 MB): AUCell completo, valores sin alterar (dif. máx. vs `auc_mtx.csv` = 3e-21, redondeo de texto), mismo orden de células que la matriz.
- `results/auc_median_by_population.csv`, `results/rss_by_population.csv`, `results/regulon_sizes.csv`.
- `data/expression_demo.tsv` → ~2 MB en `.gz` (para AUCell en vivo en Colab).
- **Se quedan solo en Stargazer:** `results/adjacencies.tsv` (75 MB), `results/regulons.csv` (5.9 MB, formato CLI), `resources/*.feather`, `data/original/`.

## Punto de control C — ✅ (2026-11-11)

Asociaciones, regulones, AUCell, metadatos y 3 figuras coherentes, todos del mismo análisis de 900 células.

**Recursos cisTarget (hg38, RefSeq r80, motivos v10 clust, Feather v2)** — base `https://resources.aertslab.org/cistarget/`:
- `motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` (99 MB; 27,072 motivos; 1,605 TF anotados)
- `databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` (311 MB; SHA1 `60fe4462…`)
- `databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` (313 MB; SHA1 `d9d4e927…`)

**Compatibilidad (paso 19):** 94.2% / 94.0% de los genes de la matriz están en cada base; 1,018 TF de la lista están en la matriz (1,015 en las bases); los 12 TF clave (SPI1, CEBPB, CEBPA, IRF8, PAX5, EBF1, SPIB, POU2AF1, TBX21, EOMES, RUNX3, PRDM1) están en lista TF, matriz y tabla de motivos.

**Adyacencias (paso 21), dianas top por importancia:**
- PAX5 → IGHM, CD79B, AFF3, BCL11A, CD40, VPREB3 · EBF1 → FCRL5, BLK, MS4A1, TNFRSF13C · SPIB → MS4A1, CD79B, BANK1, CD79A
- TBX21 → NKG7, CD247, PRF1, GNLY, GZMA, KLRD1 · EOMES → CD7, STAT4, CD247, GZMA
- SPI1 → FOS, CEBPB, FTH1, PSAP, AIF1 · CEBPB → FTH1, NEAT1, CTSS, AIF1
- 🎓 Ejemplo didáctico: las aristas más fuertes de toda la tabla son **RPS4X → genes ribosomales** (RPS4X está en la lista oficial de TF). Es coexpresión, no regulación → cisTarget debería podarlo.

## Encargo 3 — notebook (pasos 28–36)

| Paso | Estado | Evidencia |
|---|---|---|
| 28. Notebook único | ✅ | `notebooks/Taller_pySCENIC_2h.ipynb` (42 celdas, sin outputs), generado por `scripts/build_notebook.py`. Solo pandas, numpy, matplotlib, seaborn, networkx. |
| 29. Contexto (AMBOS) | ✅ · ⏳ tu revisión | Celda 1: pregunta, TF/GRN/regulón, diagrama GRNBoost2→cisTarget→AUCell, aviso explícito "no prueba causalidad", origen de datos y de etiquetas. |
| 30. Carga automática | ✅ | `cargar()` descarga de `BASE_URL` (provisional: `https://raw.githubusercontent.com/Ssdzem/taller-pyscenic-intro/main`), 3 reintentos, caché local; `TALLER_LOCAL` solo para pruebas. 7 archivos, ≈8.4 MB. |
| 31. Ejercicio 1 (GRNBoost2) | ✅ | TBX21: 3,299 candidatos (top: NKG7, CD247, PRF1, GNLY, GZMA); histograma de pesos; caso RPS4X. |
| 32. Ejercicio 2 (cisTarget) | ✅ | TBX21 3,299 → 1,014 (30.7%); motivos/NES; red MAFB(+); tabla 13 TF (RPS4X y PAX5 sin regulón). |
| 33. Ejercicio 3 (AUCell) | ✅ | Unión por ID; heatmap de medianas; boxplots; expresión del TF vs AUC (FOS, SPI1). |
| 34. Interpretación (AMBOS) | ✅ · ⏳ tu revisión | Lista de regulones específicos + `ficha(regulón)` + plantilla (hipótesis, evidencia observada vs inferida, limitación, experimento). |
| 35. AUCell en vivo (AMBOS) | ✅ local · ⏳ Colab real | Ver abajo. |
| 36. Comprobaciones | ✅ | Celda ✅ al final de cada sección; cada sección recarga sus datos con `cargar()`. |

**AUCell en vivo (sección 4, opcional), 3 niveles:**
1. **Desde cero (numpy/pandas):** réplica exacta de pySCENIC 0.12.1 (ranking con desempate aleatorio por semilla, corte = round(0.05·n) − 1 = 505, sin pesos, regla ≥80% genes). **Idéntico** al precalculado: máx |dif| = 1e-16 con pandas 1.5/numpy 1.23 y pandas 2.3/numpy 2.4 (`scripts/35_validate_numpy_aucell.py`, `logs/35_validate_numpy_aucell.log`). 0.3 s. Incluye curva de recuperación MAFB(+) en un monocito (AUC 0.340) vs una célula B (0.123).
2. **pySCENIC real:** `pip install --target pyscenic_pkgs pyscenic==0.12.1 "setuptools<81"` + ejecución en subproceso aislado (timeouts 420/300 s). **Arreglo necesario:** sin `setuptools<81` falla `import ctxcore` (`No module named 'pkg_resources'`). Probado en Python 3.11, 3.12 y 3.13 (contenedores limpios, numpy 2.4–2.5, pandas 3.1): instala en ~20 s, AUCell 1–4 s, **idéntico** (1e-16) (`scripts/35_test_pip_pyscenic.sh`, `logs/35_pip_py3*.log`).
3. **Precalculado:** si 1 y 2 fallan, la comparación usa lo que haya; los ejercicios 1–3 nunca dependen de la sección 4.

**Ejecución completa "Run all" (contenedor tipo Colab: Python 3.12, pandas 3.1.1, numpy 2.5.4, matplotlib 3.11, seaborn 0.13.2, networkx 3.7):**
- Con red: **28 s, 0 errores**, 13 figuras (`logs/36_notebook_run_local.log`; pySCENIC en vivo 25 s).
- Sin red (pip falla a propósito): **11 s, 0 errores**; mensaje amable y continúa (`logs/36_notebook_run_offline.log`).
- Imagen de prueba: `scripts/containers/colabtest/Containerfile` → `localhost/taller-colabtest:py312`; ejecutor `scripts/test_notebook.py`.
- ⚠️ Limitación: un contenedor limpio no es Colab (Colab trae cientos de paquetes preinstalados). La prueba en Colab real es tu paso 35/42, una vez publicado el repo.

## Punto de control D — ✅ local (2026-11-12)

El notebook corre de principio a fin con "Run all" en un Python limpio tipo Colab, sin Conda, sin rutas locales y sin secretos. Falta confirmarlo en Colab real tras publicar.

## Encargo 4 — publicación (pasos 37–44)

| Paso | Estado | Evidencia |
|---|---|---|
| 37. Repo público (autorizado 2026-11-12) | ✅ | https://github.com/Ssdzem/taller-pyscenic-intro (PUBLIC, rama `main`) |
| 38. `.gitignore` + licencias | ✅ | Lista blanca (todo ignorado salvo lo nombrado) + `.claude/ .codex/ .omx/`. Licencia datos 10x: **CC BY 4.0** (vía buscador; la página de 10x bloquea clientes automáticos con 403/429 → **confírmala tú en el navegador**). `LICENSE`: código MIT, material didáctico y datos derivados CC BY 4.0 (titular: Ssdzem). |
| 39. Subida + enlaces | ✅ | 36 archivos, 9.9 MB. Las 8 URLs raw responden 200 sin autenticación y su sha256 = local. Notebook bajado de GitHub + datos de GitHub (sin copia local): **30 s, 0 errores** (`logs/39_notebook_run_github.log`). |
| 40. README de una pantalla | ✅ | Público objetivo, objetivo, duración, requisitos, botón *Open in Colab*, contenido, respaldo, créditos y licencia. Botón también en la celda 1 del notebook. |
| 41. Ver repo sin sesión | ⏳ TÚ | |
| 42. Colab desde cero, Run all | ⏳ TÚ | |
| 43. Cronometrar (AMBOS) | ✅ cómputo · ⏳ humano | Tiempo de cómputo por sección: todas < 2 s salvo la 4.2 (pySCENIC en vivo ≈ 25 s). El tiempo real lo marca la lectura/discusión: cronométralo con una persona novata. |
| 44. Respaldo | ✅ | Release v1.1 (actual): `taller-pyscenic-respaldo.zip` (5.6 MB, 51 archivos, sha256 `5058da85…`), https://github.com/Ssdzem/taller-pyscenic-intro/releases/latest/download/taller-pyscenic-respaldo.zip (v1.0 sigue publicado con el notebook anterior); figuras PNG+PDF en `figures/`. |

Commits: `f2294f1` (inicial), `f88cbed` (enlace de respaldo). Local: `dist/` (ZIP, no versionado).

## Encargo 5 — docencia (pasos 45–50)

| Paso | Estado | Evidencia |
|---|---|---|
| 45. Diapositivas (AMBOS) | ✅ borrador · ⏳ tu revisión | 6 diapositivas con notas de orador (portada+objetivo, TF/GRN/regulón, tres pasos con cifras reales, abrir Colab + reloj, ejercicio en parejas observado vs inferido, limitaciones + recursos). Artifact privado: https://claude.ai/artifact/BUsry7HyPH8F3zh1gai9vr (descargable PDF/PPTX). Fuentes: `instructor/slides/project/`. Sin imágenes externas; colores por población = los del notebook. |
| 46. Respuestas esperadas (AMBOS) | ✅ borrador · ⏳ tu revisión | `instructor/SOLUCIONES.md`: todas las preguntas del notebook + 3 plantillas resueltas (EBF1, MAFB, FOS), cada afirmación etiquetada 🔎 observado / 🧮 inferido / 📚 conocimiento previo. Todas las cifras de `scripts/46_answers.py` → `logs/46_answers.log` (contenedor `taller-colabtest:py312`, exit 0). |
| Guion | ✅ borrador | `instructor/GUION.md`: checklist previo (41–43, 47, 48), mensaje para alumnos, minuto a minuto con puntos clave, tabla de contingencias y orden de recorte. |
| 47–49 | ⏳ TÚ | Enviar enlace, simulacro, impartir. |
| 50. Compartir respuestas (AMBOS) | ⏳ después del taller | `instructor/` está ignorado por git (lista blanca) → no se publica hasta que lo autorices. `scripts/46_answers.py` sin commit (pendiente de autorización). |

Hallazgos nuevos verificados (paso 46, para la discusión):
- RPS4X: **0 motivos anotados** en `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` → cisTarget no puede darle regulón por construcción.
- PAX5: 125 motivos anotados (45 directos) → la falta de regulón **no** es por falta de motivo; ningún módulo llegó a NES ≥ 3.
- Dropout: EBF1 = 0 en 206/300 células B, pero 173 de ellas tienen AUC EBF1(+) > 0; las 80 NK sin TBX21 tienen todas AUC TBX21(+) > 0.
- FOS(+) contiene FOSB, JUNB, JUND, DUSP1, KLF6, ZFP36 (respuesta temprana) y genes de identidad monocítica (S100A8, LYZ, CTSS, MNDA, VCAN, SPI1).
- SPI1(+) (3,106 genes): rango de medianas 0.018 (puesto 65/272) pero cuartiles sin solapamiento.
- Error menor: `Series.corr(method="spearman")` requiere scipy (no está en la imagen de prueba) → reemplazado por correlación de rangos.

## Siguiente acción

- **TÚ (29, 34):** revisión del contexto y del ejercicio de parejas (en curso). Si cambias algo, edito `scripts/build_notebook.py`, regenero y hago push.
- **TÚ (38):** confirmar en el navegador la línea de licencia CC BY 4.0 en la página del dataset de 10x.
- **TÚ (41):** abrir https://github.com/Ssdzem/taller-pyscenic-intro en ventana privada.
- **TÚ (35/42):** botón *Open in Colab* → *Entorno de ejecución → Ejecutar todo* en sesión limpia; anotar si la sección 4.2 (pySCENIC en vivo) funciona y cuánto tarda; repetir 2–3 veces antes de decidir si se hace en vivo.
- **TÚ (45, 46):** revisar diapositivas, `instructor/GUION.md` y `instructor/SOLUCIONES.md`.
- **TÚ (47–49):** enviar enlace, simulacro de 10 min, impartir.
- **50 (después del taller, con tu autorización):** copiar `SOLUCIONES.md` al repo, commit + push (junto con `scripts/46_answers.py`).

## Paso 51 (nuevo, a petición tuya) — ¿se puede correr el CLI de pySCENIC en Colab? (2026-11-12)

Prueba: `scripts/51_test_cli_colab.sh {r1|full} {pin|patch}` en python:3.12-slim limitado a **2 CPU / 12 GB** (como Colab gratuito); pySCENIC instalado aislado con `pip --target`. Matrices: `scripts/51_prep_reduced.py` (full = 900×10,120; r1 = 2,000 genes más variables + los 1,018 TF = 900×2,829). Logs: `logs/51_cli_*.log`.

Errores encontrados (y diagnóstico):
1. `pyscenic` CLI con numpy 2: `AttributeError: module 'numpy' has no attribute 'object'` (pyscenic/transform.py:42-44 y rss.py:21 usan `np.object`/`np.float`). La sección 4.2 no lo sufría porque importa solo `pyscenic.aucell`.
2. Arreglo "pin" (`numpy<2`): **descartado** — pyarrow reciente exige numpy ≥ 2 (`ImportError: pyarrow requires NumPy 2.0 or newer`).
3. Arreglo "patch" con pandas 3: `ctx` falla `KeyError: 'target'` (pandas 3 cambió `groupby().apply`, excluye la columna de agrupación; pyscenic/utils.py:244).
4. `pyscenic grn` (dask 2026.8): `TypeError: Must supply at least one delayed object` (problema conocido de arboreto+dask moderno) → usar `arboreto_with_multiprocessing.py` (incluido en pySCENIC, mismas opciones).
5. Error mío: la marca de tiempo se redirigía al log del paso; corregido con fd 3.

**Receta que funciona:** `pip install --target pkgs "pyscenic==0.12.1" "setuptools<81" "pandas<3"` + `sed` de 4 líneas (np.object→object, np.float→float) + `arboreto_with_multiprocessing.py` en lugar de `pyscenic grn`. numpy 2.5.4, pandas 2.3.3.

| Paso (2 CPU) | r1 (2,829 genes) | full (10,120 genes) |
|---|---|---|
| pip install | 24 s | 26 s |
| grn (arboreto_with_multiprocessing) | **121 s** | **982 s (16 min)** |
| ctx, 1 base (500 bp) | 229 s → 148 regulones | 241 s → 191 regulones |
| ctx, 2 bases | 360 s → 216 regulones | 455 s → 255 regulones (servidor: 272) |
| aucell | 7 s | 6 s |
| RAM máx. | 5.7 GB | 5.5 GB |

No medido: descarga de las bases (311 MB c/u) desde resources.aertslab.org dentro de Colab, ni el efecto de 20–30 alumnos descargando a la vez.

## Paso 52 — CLI en vivo en el notebook + notebook para casa (2026-11-12/13, autorizado: "corre todo")

Decisión del usuario: documentar cada comando (🖥️ tabla de flags) y que los alumnos corran `grn`, `ctx` y `aucell` en vivo (▶); notebook opcional `Pipeline_completo_casa.ipynb` con matriz `demo` o `completa` (dataset entero de 10x).

Cambios:
- `scripts/build_notebook.py` → notebook del taller (53 celdas): instalación de pySCENIC en sección 0; en cada sección 🖥️ el comando + flags (de `logs/51_cli_help.txt`, ayuda oficial) y ▶ córrelo tú (matriz reducida 900×2,829, 1 base de ranking); explicación de por qué 13 TF; pregunta 4 en Ej. 1 y Ej. 2; sección 4 pasa a "AUCell por dentro" (se quita la vieja 4.2/4.3); sección 6 con resumen de comandos y enlace al notebook de casa.
- `scripts/build_notebook_casa.py` → `notebooks/Pipeline_completo_casa.ipynb` (21 celdas): QC idéntico a `15_prepare_demo.py` para `completa`; recursos con SHA1; reanudación y copia a Drive opcionales; comparación con el servidor.
- Plan B para la base de ranking: carpeta compartida de Google Drive (`RUTA_DRIVE`). **Motivo:** el 2026-11-12 ~21:50 resources.aertslab.org rechazaba conexiones (`curl: (7) Failed to connect … port 443`; DNS → 134.58.50.9).
- README, `docs/METODOS.md` (sección 6b), `instructor/GUION.md` y notas de las diapositivas 3 y 4 actualizados.

Errores encontrados en las pruebas (y arreglo):
1. `os.cpu_count()` devuelve los núcleos del host (32) aunque el contenedor tenga 2 → 32 procesos: `grn` de `completa` muerto por memoria (`Killed`) y `ctx` > 15 min. Arreglo: `len(os.sched_getaffinity(0))`. En Colab ambos dan 2, pero así es robusto. Prueba con `--cpuset-cpus` (`scripts/52_test_notebook_2cpu.sh`).
2. Con numpy 2, `ctx` escribe los pesos como `np.float64(…)` → `ast.literal_eval` falla en mi lector de regulones. Arreglo: regex antes de `literal_eval`. (El `aucell` de pySCENIC sí lo lee.)
3. `pip install h5py` desde el notebook veía `PYTHONPATH=pyscenic_pkgs` (que ya trae h5py) y no lo instalaba en el kernel. Arreglo: pip sin `PYTHONPATH`. (En Colab h5py viene preinstalado.)
4. El ejecutor de pruebas cortaba celdas a 900 s → `NB_TIMEOUT` configurable.
5. Avisos `pkg_resources is deprecated` en cada comando → `PYTHONWARNINGS=ignore` para los subprocesos.

Resultados de la prueba del notebook del taller (2 CPU, 12 GB, núcleos fijos; `logs/52_taller_2cpu.log`):
- Instalación 29 s; `grn` 118 s; `ctx` 218 s; `aucell` 9 s; cómputo total ≈ 7 min.
- `aucell` real con los regulones del servidor = `auc_demo.csv` (diferencia 1e-16).
- Top 15 de TBX21: 12/15 en común con el servidor.
- En la corrida reducida no aparece MAFB(+) → ejemplo para la pregunta de robustez.
- Matriz `completa` (notebook de casa): 10,280 células × 10,165 genes (mismo nº de células que el QC original), 210 MB.
- **Prueba del notebook para casa:**
  - `demo` (2 CPU): 0 errores, 1,505 s de cómputo (`grn` 943 s, `ctx` 517 s) → 255 regulones, 212 en común con el servidor; MAFB(+) no aparece y PAX5(+) sí.
  - `completa` (6 CPU, 5 procesos, 12 GB): 0 errores, 5,268 s (`grn` 81 min, `ctx` 6 min, `aucell` 15 s) → 10,280 células, 237 regulones, 168 en común; MAFB(+) reaparece (75 genes).
  - CPU total de `grn completa` ≈ 12,200 s (40 % de sistema) → en Colab con 2 núcleos, ~1.5–3 h.
- Error 6: con 16 procesos, `grn completa` llenó los 12 GB (cada proceso copia la matriz) y avanzaba con 1.5 núcleos de 16. Arreglo: el notebook para casa elige los procesos según la RAM disponible (`memoria_gb()`, respeta `memory.max` del cgroup).
- Robustez (para la discusión): EBF1(+) y TBX21(+) encabezan B y NK en las 4 corridas; MAFB(+) aparece en 2 de 4. Diapositiva de cierre, guion y `SOLUCIONES.md` ajustados.
- Push y Release v1.1 (autorizado 2026-11-13: "ok adelante"):
  - `git push origin main` → `f88cbed..6bacc0c` (commits `2ab1876` y `6bacc0c`; este último cambia el enlace del README a `releases/latest/download/...` para que siempre apunte al último respaldo).
  - `git archive --format=zip --prefix=taller-pyscenic-respaldo/ -o dist/v1.1/taller-pyscenic-respaldo.zip HEAD` → 5.6 MB, 51 archivos, sha256 `5058da85b81f8b84f9ae46f605130e845f8175060534787e76549a3d60fa36f6`; sin `instructor/` ni `46_answers.py`.
  - `gh release create v1.1 ... --latest` → https://github.com/Ssdzem/taller-pyscenic-intro/releases/tag/v1.1
  - Verificación anónima (curl, sin credenciales): raw de los 2 notebooks, matriz y regulones → 200; enlaces de Colab de ambos notebooks → 200; ZIP de `releases/latest` con el mismo sha256; notebook remoto idéntico al local; el badge 🏠 apunta a `Pipeline_completo_casa.ipynb`.
  - Pendiente: probar "Ejecutar todo" en un Colab real (pasos 35/42).

## Paso 53 — aertslab caído: copia de los recursos en el Drive del instructor (2026-11-13)

- Síntoma (Colab real, celda 21 de la sección 2): `urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] ... Hostname mismatch, certificate is not valid for 'resources.aertslab.org'`.
- Diagnóstico (desde Stargazer):
  - `openssl s_client -connect resources.aertslab.org:443 -servername resources.aertslab.org` → certificado `CN=aertslab.org` (SAN: aertslab.org, www.aertslab.org), emitido el 2026-11-13 08:27 UTC. Le falta `resources`.
  - `curl -k` (ignorando el certificado) a la tabla de motivos y a la base 500bp → HTTP 404. El servidor sirve el sitio principal en lugar de los recursos: falla de configuración del lado de aertslab.
  - El espejo oficial `resources-mirror.aertslab.org` resuelve a la misma IP (134.58.50.9) y falla igual (curl rc=60).
  - Reportes abiertos el mismo día: aertslab/pySCENIC#644, #645; antes #641/#642 y scenicplus#643 (caída del 2026-10-26, que se recuperó ese mismo día).
  - Búsqueda de copias alternativas: API de Zenodo (0 resultados para los nombres de archivo), API de Hugging Face (datasets/modelos: nada) y búsqueda web: no hay espejos de terceros. Tampoco encontré una licencia explícita de redistribución → se descartó publicarlos en un Release de GitHub.
- Decisión del usuario: copia en su Drive, compartida "cualquier persona con el enlace"; el notebook baja con `gdown` sin montar Drive.
- Cambios:
  - `scripts/recursos_comun.py` (nuevo): `bajar()` prueba Drive del instructor (gdown) → aertslab → espejo, con SHA1 en cada fuente (si una fuente entrega un archivo equivocado se pasa a la siguiente), timeout de 30 s y 2 intentos por URL.
  - `scripts/recursos_drive.json` (nuevo): enlaces de Drive por archivo; **vacíos hasta que el usuario comparta la carpeta**.
  - `build_notebook.py` (sección 2) y `build_notebook_casa.py` (sección 2) usan esas funciones; la celda muestra `desde: <fuente>`. Se conserva `RUTA_DRIVE` / `RUTA_DRIVE_RECURSOS` (Drive montado) como último recurso.
  - `docs/METODOS.md` (recursos) e `instructor/GUION.md` (lista previa y contingencia) actualizados.
- Pruebas (`localhost/taller-colabtest:py312`):
  - aertslab caído y sin enlace de Drive → falla en 11 s con un mensaje que lista ambas fuentes.
  - Una fuente entrega un archivo equivocado → la descarta por SHA1 y usa la siguiente.
  - Enlace de Drive inexistente → error capturado y pasa a la URL siguiente. gdown no venía en la imagen y se instaló solo (6.5.1); Colab ya lo trae.
- Carpeta del usuario: `rclone link first:taller_pyscenic/resources` → https://drive.google.com/open?id=1a8HVL1UjsJKWQ8SHsKL0ySKSB7o-yX4P. IDs obtenidos con `rclone lsjson --hash`; el SHA1 de Drive coincide con el esperado en los 4 archivos. Los enlaces quedaron en `recursos_drive.json` y se regeneraron ambos notebooks.
- Descarga anónima real (contenedor, sin credenciales): 500bp 313 MB en 13 s, tabla 99 MB en 4 s, 10kbp 311 MB en 11 s, allTFs en 2 s, todas `desde: Drive del instructor`. La celda real del notebook de clase → `RECURSOS_LISTOS = True`.
- Notebook de clase completo (2 CPU, con copia local de recursos): 0 errores, 1,480 s. Fue más lento de lo normal (395 s) porque otro proceso del usuario (stellarscope, 16 procesos) ocupaba la máquina.
- Push (autorizado 2026-11-13, "yes go ahead"): `6bacc0c..91857a8`. Comprobado: el notebook de `main` en raw.githubusercontent trae el ID de Drive de la base 500bp. No se creó Release v1.2.

## Paso 54 — Teoría intercalada en diapositivas y guion (2026-11-13)

- Pedido del usuario: primero qué es una GRN, para qué sirve y por qué pySCENIC; después GRNBoost2, cisTarget y AUCell por dentro, y la biología de los ejemplos, yendo y viniendo entre el notebook y la presentación.
- Diapositivas (artefacto v11, https://claude.ai/artifact/BUsry7HyPH8F3zh1gai9vr): de 6 a 13.
  - Nuevas: `grn`, `para-que`, `por-que`, `grnboost`, `cistarget`, `poda`, `aucell` y `protagonistas`.
  - `conceptos` se eliminó: su contenido pasó a `grn`.
  - Ajustadas: notas de `cover` (ahora sí corren comandos), `pasos`, reloj de `colab`, minutos de `parejas` (18 min) y de `limites`.
  - CSS: `.grid2`, `.steplist`, `.figure`, `table.data`, etc.
  - Los diagramas (curvas de recuperación, árbol) están marcados como "esquema".
- Números nuevos, todos calculados de `results/` y `data/` (scripts en línea, no guardados):
  - Ejemplo de GRNBoost2: importancias de NKG7 ← TBX21 66.1, RORA 49.2, MYBL1 37.4; CD79A ← ANXA1 45.1, CEBPD 33.2, SPIB 27.6.
  - Spearman sobre log(CP10k): +0.73, +0.67, +0.60; −0.73, −0.65, +0.65. ANXA1 se detecta en 98 % de los monocitos y en 13 % de las B.
  - Bases de ranking (leídas con la imagen oficial): 5,876 motivos × 27,015 genes (500bp) / 27,090 genes (10kbp).
  - Candidatos → regulón: SPI1 5,829 → 3,106; TBX21 3,299 → 1,014; EBF1 1,625 → 14 (1 motivo, NES 3.37); MAFB 4,712 → 37; PAX5 2,158 → —; ANXA1 6,532 → — (0 filas en la tabla de anotación); RPS4X 6,978 → —.
  - EBF1(+) tiene AUC > 0 en 259/300 células B y en 134/600 del resto. El gen EBF1 se detecta en 31 % de las B.
  - Aristas del diagrama de red verificadas en `regulons_demo.csv`: SPI1 → LYZ y FCN1; EBF1 → CD24 y MZB1; TBX21 → NKG7 y GZMB.
- Corrección: la mediana de EBF1(+) en B es 0.04948 → **0.049**. Decía 0.050 por doble redondeo; corregido en `SOLUCIONES.md` (2 sitios) y en `GUION.md`.
- `GUION.md`:
  - Minuto a minuto reescrito, con tabla de 17 bloques: 0–15 teoría; en los pasos 1–2 se lanza el comando y se explica mientras corre; en el paso 3 se explica antes; 84–90 protagonistas; 90–108 parejas; 108–120 cierre.
  - Nueva fila de ensayo en la lista previa; orden de recorte actualizado; enlace de respaldo cambiado a `releases/latest`.
- Notebook (local, **sin push**): mapa de tiempos alineado al nuevo reloj y nota "🔁 mientras corren grn y ctx…"; parejas 18 min.

## Paso 55 — GRNBoost2 y cisTarget paso a paso, con datos reales (2026-11-13)

- Pedido del usuario: las slides de GRNBoost2 y cisTarget no se entendían ("hay unos cuates que no saben ni qué es un motivo"). GRNBoost2: árbol de regresión → boosting → importancia. cisTarget: motivo → bases → ranking → AUC y NES → poda. Una idea por slide.
- Cálculos (todo en `docker.io/aertslab/pyscenic:0.12.1`, scripts y JSON en `instructor/slides/datos/`, no versionado):
  - `grnboost_demo.py`, `grnboost_extra.py`: reentrenan el modelo de NKG7 y CD79A con `arboreto.core.fit_model` (SGBM_KWARGS, ventana 25, semilla 2026) sobre `expression_demo.tsv` y los 1,018 TF.
    - Reproduce exactamente `adjacencies.tsv`: NKG7 ← TBX21 66.09, RORA 49.18, MYBL1 37.38; CD79A ← ANXA1 45.11, CEBPD 33.16, SPIB 27.56.
    - NKG7: 482 árboles (suma de importancias = 482); TBX21 sale en 105 árboles y aporta el 13.7 %. Cada árbol: profundidad 3, 101 TF, 90 % de células, learning rate 0.01.
    - Árbol de una pregunta (DecisionTreeRegressor, profundidad 1): TBX21 > 0.5 → 225 células (220 NK, 5 B), NKG7 = 32.6; resto 675 (300 Mono, 295 B, 80 NK), NKG7 = 3.8. ECM 277 → 121.
    - Célula NK ACAGGGACAACTGCTA-1 (38 UMI): predicción 10.98 → 23.85 (100 árboles) → 32.92 (482). Monocito GGGACAAGTGAGATCG-1 (0): 10.98 → 0.20.
  - `cistarget_demo.py`: módulos con `modules_from_adjacencies` (como el CLI), base del promotor, `ctxcore.recovery`.
    - SPI1, módulo weight>90.0 % (2,694 genes, 2,600 en la base), motivo metacluster_167.5 (PU.1, anotación directa): AUC 0.05721, NES 4.3066, primero de 5,876; 7 motivos con NES ≥ 3; promedio de AUC 0.03202, desviación 0.00585. Corte del leading edge en el rango 4,495 → 740 genes. Coincide con la fila de `regulons.csv`.
    - 280 candidatas en los primeros 1,351 lugares (esperadas al azar 130). Lugares en el ranking de PU.1: CSF1R 100, LYZ 762, S100A8 1,320, CD79A 21,861. CSF1R, LYZ y S100A8 están en el leading edge y en SPI1(+).
    - Motivos de PU.1 en la tabla: AGAGGAAGTG_PU.1 (HOMER) y SPI1_full_AAAAAGCGGAAGTN (Taipale 2013). 390 motivos distintos anotados a SPI1.
- **Hallazgo (`check_le.py`):** en pySCENIC 0.12.1 las columnas TargetGenes y RankAtMax de `regulons.csv` pueden corresponder a otro motivo del mismo módulo.
  - Causa probable (transform.py, `module2features_auc1st_impl`): `annotated_features` se reordena con `sort_values` y la máscara booleana resultante se aplica a `rccs` en el orden original.
  - EBF1, módulo weight>90 %: la fila de tfdimers__MD00215 (RankAtMax 479, 14 genes) es exactamente el leading edge de jaspar__MA1509.1 (IRF6). El de MD00215 sería 801 → 18 genes.
  - SPI1 top5perTarget: 10 de 10 filas desplazadas. SPI1 weight>90 %: 3 de 4 desplazadas; la de metacluster_167.5 coincide.
  - AUC y NES sí son correctos. Afecta qué dianas entran a cada regulón (siempre dentro de las candidatas del módulo).
  - No se corrigió nada del pipeline; queda anotado en el GUION y en las notas de la slide 17. Pendiente: buscar si ya está reportado en GitHub (aertslab/pySCENIC).
- Diapositivas (artefacto v15): de 13 a 21.
  - Nuevas, generadas con `instructor/slides/datos/gen_slides_teoria.py`: `arbol`, `boosting`, `importancia`, `motivo`, `bases`, `ranking`, `curva`, `nes` y `podacurva`. `grnboost` rehecha (tabla y ANXA1). `cistarget` eliminada.
  - Revisadas con chrome-headless-shell (Playwright) a 1920×1080: ninguna rebasa el margen inferior (máx. 949 px de 952).
  - Minutos ajustados en `colab` (reloj), `poda`, `aucell`, `protagonistas` y `parejas` (17 min).
- Nuevo reloj: 0–15 teoría · 15–22 Colab · 22–44 paso 1 (slides 7–10 de 24 a 33) · 44–68 paso 2 (11–12 de 44 a 48, `ctx` a los 48, 13–17 de 50 a 60) · 68–86 paso 3 · 86–91 protagonistas · 91–108 parejas · 108–120 cierre.
- Actualizados: `GUION.md` (tabla, pasos 1 y 2, recortes) y `scripts/build_notebook.py` con el notebook de clase regenerado (mapa de tiempos; parejas 17 min). Local, **sin push**.

## Paso 56 · La red explicada en 4 diapositivas (2026-11-14)

- Origen: comentario de la autora en la diapositiva 2 del artifact ("Ni siquiera comentas qué es una red de regulación génica").
- La antigua diapositiva 2 se divide en cuatro (deck: 21 → 24 diapositivas, artifact versión 17):
  1. `grn.html` · Un gen tiene un interruptor: activador / represor; TF, motivo, diana.
  2. `red.html` · ¿Qué es una GRN?: nodos, aristas (dirección y signo), regulón. Subred real: SPI1 ⇄ CEBPB → LYZ, FCN1; SPIB → EBF1 → MZB1, CD24; TBX21 → EOMES → NKG7, GZMB.
  3. `circuitos.html` · Autorregulación, hacia adelante, retroalimentación ("motivos de red", aclarado frente al motivo de ADN).
  4. `estados.html` · Tipos celulares como estados estables (paisaje de Waddington, esquema); C/EBPα convierte células B en macrófagos (Xie et al. 2004).
- Datos (de `results/regulons_long.tsv`, sin recalcular nada): 272 regulones, 8,929 genes, 46,339 aristas, todas (+); 1,579 entre TF; SPI1 3,106 dianas (mediana 41); mediana de 5 TF por gen diana. Hacia adelante reales: SPI1→CEBPB→KLF4 y SPI1→KLF4; SPIB→POU2AF1→EBF1 y SPIB→EBF1; TBX21→EOMES→RUNX3 y TBX21→RUNX3. Mutua: SPI1→CEBPB (38.2) y CEBPB→SPI1 (21.5).
- Ojo: 153 de 272 regulones contienen a su propio TF, pero es un artefacto del método: `pyscenic/utils.py` (`add_tf`, línea 396) añade el TF a su módulo con peso 1.0, y solo sobrevive si su promotor tiene el motivo. Anotado en las notas de `circuitos.html` y en el GUION.
- Reloj: teoría 0–16 (2–4, 4–6, 6–8, 8–10 la red; 10–12 para qué; 12–14 por qué; 14–16 pasos); Colab 16–22. Actualizados: `colab.html`, `cover.html`, `para-que.html`, `por-que.html`, `pasos.html`, `limites.html`, `deck.json`, `GUION.md` (numeración +3) y `scripts/build_notebook.py` (16–22), con el notebook regenerado (53 celdas). Local, **sin push**.
- Render 1920×1080: grn 894, red 875, circuitos 935, estados 952 (límite 952).

## Paso 57 · El árbol de GRNBoost2 en tres diapositivas (2026-11-14)

- Origen: comentario de la autora en la diapositiva del árbol (pedía explicar cómo se calcula la predicción, cuántos niveles y cuántos árboles por gen, y mostrar las células en una recta con el corte del TF).
- Nuevo script `instructor/slides/datos/arbol_recta.py` (contenedor `aertslab/pyscenic:0.12.1`, matriz `data/expression_demo.tsv` de 900 × 10,120 y 1,018 TF) → `arbol_recta.json`, log en `arbol_recta.log`. Calcula:
  - TBX21 por célula: en 0 hay 300 Mono, 295 B y 80 NK; en 1–6 hay 220 NK y 5 B.
  - Error (MSE) de cada corte de TBX21: ≥1 → 121, ≥2 → 191, ≥3 → 244, ≥4 → 266, ≥5 → 271, ≥6 → 273. Sin cortar: 277.
  - Mejor corte de cada uno de los 1,018 TF: TBX21 121, RORA 145, KLRG1 162, MYBL1 168, MEF2C 174, STAT4 186, …, SPI1 197, EBF1 264.
  - Árbol real de profundidad 3 (todas las células y todos los TF): TBX21 → KLRG1 / UBB → RORA, MEF2C, RPL6, SNRPB2. Ocho hojas, de 0.7 (592 células) a 58.6 (8 NK). MSE 52.
  - Árboles por gen, mismo `fit_model` que grn (semilla 2026, parada temprana a 25): NKG7 482 (igual que antes). En 60 genes al azar: mediana 64, mínimo 25, máximo 481. Marcadores: LYZ 390, MS4A1 416, GZMB 368, CD79A 374, S100A8 372.
- `gen_slides_teoria.py`:
  - El bloque del árbol se reemplaza por tres diapositivas: `arbol` (recta de TBX21 con las 900 células, el corte y el promedio de cada lado), `corte` (todos los cortes, todos los TF) y `arbol3` (árbol de 3 niveles).
  - Boosting corregido: `max_features` de sklearn se sortea en cada pregunta (split), no una vez por árbol. Se añaden los árboles por gen.
  - Los rótulos "x de 4" pasan a "x de 6".
- Deck: 24 → 26 diapositivas (artifact versión 18). GRNBoost2 ocupa 24–35 y el Ejercicio 1 35–44. GUION renumerado (+2 desde la 11). Las notas de `grn`, `circuitos` y `red` apuntan a las diapositivas 16 y 26.
- Render 1920×1080: arbol 894, corte 935, arbol3 941, boosting 921, importancia 948, grnboost 854.

## Paso 58 · Boosting explicado con residuos (2026-11-14)
- Comentario en la diapositiva 13: el paso 2 ahora dice que cada árbol predice el residuo (real − predicción) y la predicción avanza 1 % (learning rate 0.01). Recuadro: "promedio + 482 pequeñas correcciones del residuo".
- Notas: residuo = gradiente del error cuadrático; ejemplo real de la célula NK (38 UMI, inicio 11.0, residuo 27.0).
- Cambio en `instructor/slides/datos/gen_slides_teoria.py`; render 921 px (≤ 952). Publicado como versión 20.

## Paso 59 · Importancia como suma de la baja del error (2026-11-14)
- Comentario en la diapositiva 14: la importancia se enseña como la suma, en todos los árboles, de la baja del error que aporta cada TF. Las barras se reemplazan por una tabla árbol por árbol de TBX21.
- Nuevo `instructor/slides/datos/importancia_tabla.py` (contenedor aertslab/pyscenic:0.12.1, sklearn 1.1.3, semilla 2026) → `importancia_tabla.json` y `importancia_tabla.log`. Usa `tree_.compute_feature_importances(normalize=False)` en cada árbol.
- Resultado NKG7: TBX21 en 105 de 482 árboles; suma de bajas 1,802.34; baja promedio de un árbol 27.269; 1,802.34 ÷ 27.269 = 66.09 = importancia de adjacencies.tsv. Las 8 primeras se reproducen con diferencia máxima 7e-15 (arboreto = fracción de la baja total × número de árboles).
- GUION (diapositiva 14) actualizado. Render 945 px (≤ 952). Publicado como versión 21.

## Paso 60 · Corte de la diapositiva 10 rotulado (2026-11-14)
- Comentario: la línea de corte no parecía estar en 1. Es correcto que vaya entre 0 y 1 (conteos enteros; sklearn guarda el umbral como 0.5, `arbol_recta.py` lo reporta como floor(umbral)+1 = 1). Se rotula "corte: ¿TBX21 ≥ 1?", con "0 | ≥ 1" y "← no | sí →", y el paso 3 lo explica. Render 901 px. Versión 22.
- Segundo comentario (Paso 60b): eje X uniforme (1 UMI = 128 px, marcas 0–6); cada célula se dibuja entre v y v + 0.9, así la línea del corte queda exactamente en la marca del 1. Render 940 px. Versión 23.
