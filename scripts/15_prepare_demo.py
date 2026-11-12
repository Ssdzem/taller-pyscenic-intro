"""Pasos 15–16: matriz demo reproducible para pySCENIC (PBMC10k, 10x Genomics v3).

Contraste (paso 14, elegido por el instructor): Monocitos CD14+ vs Células B vs NK.

Flujo (todo antes de inferir la red; pySCENIC se corre después SOBRE ESTE MISMO subconjunto):
  1. Carga conteos UMI crudos; suma conteos de símbolos génicos duplicados.
  2. QC con los umbrales del notebook oficial PBMC10k_SCENIC-protocol-CLI:
       células: >=200 genes, <4000 genes, <15% mitocondrial; genes: detectados en >=3 células.
  3. Etiqueta por marcadores canónicos = mapa de clusters graphclust de 10x (logs/13_clusters_markers.log).
  4. Verificación independiente con CellTypist (Immune_All_Low, predicción por célula, sin majority voting).
  5. Conserva solo células donde ambos métodos coinciden.
  6. Submuestreo estratificado: N_PER_GROUP células por población, semilla fija.
  7. Filtro de genes en el subconjunto con los umbrales del notebook oficial:
       conteo total >= 3 * 0.01 * n_células  y  detectado en >= 0.01 * n_células.
  8. Exporta data/expression_demo.tsv (células × genes, conteos crudos) y data/metadata_demo.csv.

Uso (dentro de localhost/taller-celltypist:1.6.3):
  python scripts/15_prepare_demo.py
"""
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

SEED = 2026
N_PER_GROUP = 300
H5 = "data/original/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
CLUSTERS = "data/original/analysis/clustering/graphclust/clusters.csv"
TSNE = "data/original/analysis/tsne/2_components/projection.csv"
OUT_EXPR = "data/expression_demo.tsv"
OUT_META = "data/metadata_demo.csv"
OUT_REPORT = "logs/15_prepare_demo_report.json"

# Paso 3: anotación por marcadores canónicos (evidencia en logs/13_clusters_markers.log).
MARKER_MAP = {1: "Mono_CD14", 4: "Mono_CD14", 12: "Mono_CD14", 7: "B", 8: "B", 9: "NK"}


def celltypist_broad(label: str) -> str:
    """Colapsa etiquetas finas de Immune_All_Low a las tres poblaciones del contraste."""
    if label == "Classical monocytes":
        return "Mono_CD14"
    if "B cells" in label and "plasma" not in label.lower():
        return "B"
    if "NK" in label and "NKT" not in label and "T cells" not in label:
        return "NK"
    return "other"


def load_10x_h5(path):
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:]).T.tocsr()
        barcodes = m["barcodes"][:].astype(str)
        symbols = m["features"]["name"][:].astype(str)
    return X.astype(np.float32), barcodes, symbols


def main():
    rep = {"seed": SEED, "n_per_group": N_PER_GROUP}

    # 1. Carga y colapso de duplicados.
    X, barcodes, symbols = load_10x_h5(H5)
    rep["raw_shape"] = list(X.shape)
    uniq, inv = np.unique(symbols, return_inverse=True)
    agg = sp.csr_matrix((np.ones(len(symbols), np.float32), (np.arange(len(symbols)), inv)), shape=(len(symbols), len(uniq)))
    X = (X @ agg).tocsr()
    rep["duplicated_symbols_summed"] = int(len(symbols) - len(uniq))
    adata = sc.AnnData(X=X, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=uniq))
    assert np.allclose(adata.X.data, np.round(adata.X.data)), "se esperaban conteos enteros"

    # 2. QC (notebook oficial).
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    keep = (adata.obs.n_genes_by_counts >= 200) & (adata.obs.n_genes_by_counts < 4000) & (adata.obs.pct_counts_mt < 15)
    rep["qc_cells_removed"] = int((~keep).sum())
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    rep["post_qc_shape"] = list(adata.shape)

    # 3. Marcadores (clusters de 10x).
    clu = pd.read_csv(CLUSTERS, index_col=0)["Cluster"]
    assert set(adata.obs_names) <= set(clu.index)
    adata.obs["graphclust"] = clu.reindex(adata.obs_names).values
    adata.obs["population_markers"] = adata.obs["graphclust"].map(MARKER_MAP).fillna("other")

    # 4. CellTypist (requiere log1p de CP10k sobre todos los genes).
    import celltypist
    from celltypist import models
    models.download_models(model="Immune_All_Low.pkl", force_update=False)
    model = models.Model.load("Immune_All_Low.pkl")
    rep["celltypist_model"] = {"name": "Immune_All_Low.pkl", "version": model.description.get("version", "?"),
                               "date": model.description.get("date", "?"), "celltypist": celltypist.__version__}
    norm = adata.copy()
    sc.pp.normalize_total(norm, target_sum=1e4)
    sc.pp.log1p(norm)
    pred = celltypist.annotate(norm, model=model, majority_voting=False)
    lab = pred.predicted_labels["predicted_labels"].astype(str)
    adata.obs["celltypist_label"] = lab.reindex(adata.obs_names).values
    adata.obs["celltypist_conf"] = pred.probability_matrix.max(axis=1).reindex(adata.obs_names).values
    adata.obs["population_celltypist"] = adata.obs["celltypist_label"].map(celltypist_broad)

    # Tabla de concordancia (solo clusters del contraste).
    sub = adata.obs[adata.obs.population_markers != "other"]
    rep["concordance_by_marker_population"] = {
        p: sub.loc[sub.population_markers == p, "celltypist_label"].value_counts().head(8).to_dict()
        for p in ["Mono_CD14", "B", "NK"]}

    # 5. Concordancia.
    agree = adata.obs.population_markers.eq(adata.obs.population_celltypist) & adata.obs.population_markers.ne("other")
    rep["agree_counts"] = adata.obs.loc[agree, "population_markers"].value_counts().to_dict()
    rep["marker_counts_post_qc"] = adata.obs.population_markers.value_counts().to_dict()

    # 6. Submuestreo estratificado reproducible.
    chosen = (adata.obs[agree].groupby("population_markers", group_keys=False)
              .apply(lambda d: d.sample(n=N_PER_GROUP, random_state=SEED)).index)
    demo = adata[chosen].copy()

    # 7. Filtro de genes (umbrales del notebook oficial, sobre el subconjunto).
    n = demo.n_obs
    counts = np.asarray(demo.X.sum(0)).ravel()
    ncells = np.asarray((demo.X > 0).sum(0)).ravel()
    gkeep = (counts >= 3 * 0.01 * n) & (ncells >= 0.01 * n)
    demo = demo[:, gkeep].copy()
    rep["gene_filter"] = {"min_counts": 3 * 0.01 * n, "min_cells": 0.01 * n, "genes_kept": int(gkeep.sum())}
    rep["demo_shape"] = list(demo.shape)

    # 8. Exportación.
    tsne = pd.read_csv(TSNE, index_col=0)
    meta = demo.obs[["population_markers", "graphclust", "celltypist_label", "celltypist_conf",
                     "n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
    meta.columns = ["population", "cluster_10x_graphclust", "celltypist_label", "celltypist_conf",
                    "n_genes", "n_counts", "pct_mito"]
    meta["TSNE1"] = tsne.loc[meta.index, "TSNE-1"].values
    meta["TSNE2"] = tsne.loc[meta.index, "TSNE-2"].values
    meta.index.name = "cell_id"
    meta = meta.sort_values(["population"]).round({"celltypist_conf": 4, "pct_mito": 3, "TSNE1": 4, "TSNE2": 4})
    expr = pd.DataFrame(demo[meta.index].X.toarray().astype(np.int32), index=meta.index, columns=demo.var_names)
    expr.index.name = "cell_id"
    expr.to_csv(OUT_EXPR, sep="\t")
    meta.to_csv(OUT_META)
    rep["files"] = {OUT_EXPR: os.path.getsize(OUT_EXPR), OUT_META: os.path.getsize(OUT_META)}

    json.dump(rep, open(OUT_REPORT, "w"), indent=2, ensure_ascii=False, default=str)
    print(json.dumps(rep, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    sys.exit(main())
