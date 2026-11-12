# Métodos y procedencia de los datos

Todo lo que aparece en el notebook proviene de **un único análisis** de pySCENIC sobre la misma matriz de 900 células. Este documento describe cómo se generó y cómo reproducirlo.

## 1. Datos de partida

- **Dataset:** *10k PBMCs from a Healthy Donor (v3 chemistry)*, 10x Genomics, Cell Ranger 3.0.0 (publicado 2018-11-19).
  Página: <https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0>
- **Licencia:** Creative Commons Attribution 4.0 International (**CC BY 4.0**). Los archivos de `data/` y `results/` son **obras derivadas** (QC, submuestreo, filtrado de genes y resultados de pySCENIC); ver sección 3.
- Archivos usados (no incluidos en este repositorio):
  - `pbmc_10k_v3_filtered_feature_bc_matrix.h5` — sha256 `ebc5dedc938830e20f8e1aafb893b0a9e0bf88584f3ef2d00b232dd277a188af`
  - `pbmc_10k_v3_analysis.tar.gz` (clusters *graphclust*, t-SNE) — sha256 `78c5565a92ab8e5e781f10b8ba438da455ac453f7b93ad154e0b0737b77dc907`
  - Descarga: `https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/<archivo>`
- **Caso de referencia:** notebooks `PBMC10k_SCENIC-protocol-CLI.ipynb` y `PBMC10k_downstream-analysis.ipynb` de <https://github.com/aertslab/SCENICprotocol>.

## 2. Anotación de poblaciones

El dataset **no incluye tipos celulares** (solo clusters numerados y t-SNE). Las etiquetas `Mono_CD14`, `B` y `NK` son **anotación nuestra**:

1. **Marcadores canónicos** sobre los clusters *graphclust* de 10x: 1/4/12 → monocitos CD14+ (CD14, LYZ); 7/8 → células B (MS4A1, CD79A); 9 → NK (GNLY, NKG7, FCGR3A).
2. **Verificación independiente** célula a célula con CellTypist 1.6.3, modelo `Immune_All_Low.pkl` v2 (2022-07-16), sin *majority voting*.
3. Se conservan **solo las células donde ambos métodos coinciden**: Mono_CD14 96.7 %, B 98.8 %, NK 91.2 % de concordancia (las NK discordantes eran sobre todo "Tem/Temra cytotoxic T").

## 3. Matriz demo (`data/`)

Script: `scripts/15_prepare_demo.py` (imagen `scripts/containers/celltypist/Containerfile`).

1. 11,769 células × 33,538 genes; se suman los conteos de 24 símbolos génicos duplicados.
2. **QC** (umbrales del notebook oficial): ≥200 y <4,000 genes por célula, <15 % mitocondrial; genes detectados en ≥3 células → 10,280 células.
3. Anotación (sección 2).
4. **Submuestreo estratificado:** 300 células por población, `random_state=2026`.
5. **Filtro de genes** en el subconjunto (umbrales oficiales): ≥27 conteos totales (3 × 1 % de 900) y detectado en ≥9 células (1 %) → **10,120 genes**.
6. Salidas:
   - `expression_demo.tsv.gz`: 900 células (filas) × 10,120 genes (columnas), **conteos UMI crudos** enteros.
   - `metadata_demo.csv`: `cell_id, population, cluster_10x_graphclust, celltypist_label, celltypist_conf, n_genes, n_counts, pct_mito, TSNE1, TSNE2`. Las coordenadas t-SNE son las de Cell Ranger (calculadas sobre todas las células).

Verificación (`scripts/16_check_demo.py`, con el lector del CLI de pySCENIC): IDs únicos e idénticos y en el mismo orden en matriz y metadatos; sin NaN; ninguna célula ni gen vacío.

## 4. pySCENIC

- **Software:** imagen oficial `docker.io/aertslab/pyscenic:0.12.1` (digest `sha256:0c06b8b0a00117e1a4b61303e9ad775bf53c9e410d4c344151c15bf0d143a288`; Python 3.10.6, numpy 1.23.5, pandas 1.5.1, dask 2022.11.1), ejecutada con Podman.
- **Recursos** (<https://resources.aertslab.org/cistarget/>; hg38, RefSeq r80, motivos v10 *clust*, Feather v2; SHA1 verificados):
  - `tf_lists/allTFs_hg38.txt` (1,892 TF; 1,018 presentes en la matriz)
  - `motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`
  - `databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
  - `databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- **Comandos** (`scripts/run_pbmc_pipeline.sh`; opciones no mostradas = valores por defecto):

```bash
pyscenic grn expression_demo.tsv allTFs_hg38.txt -o adjacencies.tsv --num_workers 20 --seed 2026
pyscenic ctx adjacencies.tsv <rankings_10kbp>.feather <rankings_500bp>.feather \
    --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname expression_demo.tsv --mode custom_multiprocessing --num_workers 20 \
    --output regulons.csv
pyscenic aucell expression_demo.tsv regulons.csv -o auc_mtx.csv --num_workers 20 --seed 2026
```

- **Resultados:** 2,339,328 asociaciones TF–gen (1,018 TF); **272 regulones** (todos activadores, mediana 41 genes); AUCell 900 × 272.

## 5. Archivos de `results/`

| Archivo | Contenido | Filtrado |
|---|---|---|
| `adjacencies_demo.tsv` | Todas las aristas GRNBoost2 de 13 TF didácticos (SPI1, CEBPB, CEBPA, IRF8, PAX5, EBF1, SPIB, POU2AF1, TBX21, EOMES, RUNX3, PRDM1, RPS4X) | Solo por TF; importancias sin modificar (4 decimales) |
| `adjacencies_tf_summary.csv` | Nº de dianas candidatas por TF en la tabla completa | — |
| `regulons_demo.csv` | Los 272 regulones en formato largo: `regulon, TF, target, weight` | Ninguno (reformateo de `regulons.csv`) |
| `motifs_summary.csv` | Por TF y motivo enriquecido: NES, AUC, nº de dianas, base | Sin las listas de genes |
| `auc_demo.csv` | Matriz AUCell completa (células × regulones) | **Ninguno; valores sin alterar** |
| `auc_median_by_population.csv`, `rss_by_population.csv`, `regulon_sizes.csv` | Resúmenes por población (mediana de AUC, *regulon specificity score*) y tamaños | — |

La tabla completa de adyacencias (75 MB), el `regulons.csv` en formato CLI y las bases cisTarget no se incluyen por tamaño; se regeneran con los scripts.

## 6. AUCell en el notebook

- **Reimplementación en numpy/pandas** (sección 4.1 del notebook): réplica de `pyscenic.aucell` 0.12.1 con las opciones del CLI (sin pesos, `auc_threshold = 0.05`, corte de ranking `round(0.05·n) − 1`, desempate aleatorio con semilla 2026, regla de ≥80 % de genes presentes). Diferencia máxima con `auc_demo.csv`: 1e-16 (`scripts/35_validate_numpy_aucell.py`).
- **pySCENIC por pip** (sección 4.2): `pyscenic==0.12.1` + `setuptools<81` (sin este último, `ctxcore` falla al importar `pkg_resources`). Probado en Python 3.11–3.13 (`scripts/35_test_pip_pyscenic.sh`); resultado idéntico.

## 7. Reproducir

```bash
bash scripts/run_test_oficial.sh {grn|ctx|aucell}      # prueba técnica con el mini-dataset oficial
python scripts/15_prepare_demo.py                     # en la imagen CellTypist (requiere data/original/)
bash scripts/run_pbmc_pipeline.sh grn && bash scripts/run_pbmc_pipeline.sh ctx && bash scripts/run_pbmc_pipeline.sh aucell
python scripts/25_postprocess.py                      # en la imagen pySCENIC
python scripts/26_figures.py                          # en la imagen CellTypist
python3 scripts/build_notebook.py                     # regenera el notebook
```

Nota: el ejemplo oficial de SCENICprotocol trae una base de ranking en formato Feather v1, que pySCENIC 0.12.1 ya no lee; se convirtió con `convert_cistarget_databases_v1_to_v2.py` de <https://github.com/aertslab/create_cisTarget_databases>.

## 8. Referencias

- Aibar S. *et al.* (2017). SCENIC: single-cell regulatory network inference and clustering. *Nature Methods* 14, 1083–1086.
- Van de Sande B. *et al.* (2020). A scalable SCENIC workflow for single-cell gene regulatory network analysis. *Nature Protocols* 15, 2247–2276.
- Domínguez Conde C. *et al.* (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science* 376, eabl5197 (CellTypist).
