"""Pasos 23, 25 y 27: revisar regulones, unir AUCell con metadatos POR ID, resumir y exportar versiones ligeras.

Corre en la imagen oficial pySCENIC 0.12.1 (usa sus lectores):
  podman run --rm -v $PWD:/w:Z -w /w docker.io/aertslab/pyscenic:0.12.1 python scripts/25_postprocess.py

Exporta (results/):
  regulons_long.tsv           regulon, TF, target, weight (TODOS los regulones; formato simple para clase)
  motifs_summary.csv          por TF/motivo: NES, AUC, nº de dianas (resumen de regulons.csv, sin listas)
  regulon_sizes.csv           regulon, n_targets
  auc_median_by_population.csv
  rss_by_population.csv       regulon specificity score (pySCENIC) por población
  adjacencies_demo.tsv        TODAS las aristas GRNBoost2 de los TF didácticos (filtro documentado abajo)
  adjacencies_tf_summary.csv  nº de dianas candidatas por TF en adjacencies.tsv completo
  regulons_demo.csv           = regulons_long.tsv (nombre del roadmap)
  auc_demo.csv                = auc_mtx.csv sin alterar valores (copia con índice 'cell_id')
"""
import json

import numpy as np
import pandas as pd
from pyscenic.cli.utils import load_signatures
from pyscenic.rss import regulon_specificity_scores

R = "results"
TEACH_TFS = ["SPI1", "CEBPB", "CEBPA", "IRF8", "PAX5", "EBF1", "SPIB", "POU2AF1",
             "TBX21", "EOMES", "RUNX3", "PRDM1", "RPS4X"]
rep = {}

# --- Paso 23: regulones -------------------------------------------------------
regs = load_signatures(f"{R}/regulons.csv")
long = pd.DataFrame([(r.name, r.transcription_factor, g, w)
                     for r in regs for g, w in r.gene2weight.items()],
                    columns=["regulon", "TF", "target", "weight"]).sort_values(["regulon", "weight"], ascending=[True, False])
long.to_csv(f"{R}/regulons_long.tsv", sep="\t", index=False)
long.to_csv(f"{R}/regulons_demo.csv", index=False)
sizes = long.groupby("regulon").size().rename("n_targets").sort_values(ascending=False)
sizes.to_csv(f"{R}/regulon_sizes.csv")
rep["n_regulons"] = len(regs)
rep["regulon_size"] = sizes.describe().round(1).to_dict()
rep["n_activating"] = int(sum(n.endswith("(+)") for n in sizes.index))
rep["n_repressing"] = int(sum(n.endswith("(-)") for n in sizes.index))
rep["teach_tfs_with_regulon"] = {tf: int(sizes.get(f"{tf}(+)", 0)) for tf in TEACH_TFS}

mot = pd.read_csv(f"{R}/regulons.csv", header=[0, 1], index_col=[0, 1])
mot.columns = mot.columns.droplevel(0)
ms = mot[["NES", "AUC", "Context"]].copy()
ms["n_targets"] = mot["TargetGenes"].map(lambda s: s.count("('"))
ms["database"] = ms["Context"].str.extract(r"(hg38_[^']+?_v10_clust)")[0]
ms = ms.drop(columns="Context").reset_index().sort_values(["TF", "NES"], ascending=[True, False])
ms.to_csv(f"{R}/motifs_summary.csv", index=False)
rep["n_enriched_motif_rows"] = len(ms)

# --- Paso 21/27: adyacencias ligeras -----------------------------------------
adj = pd.read_csv(f"{R}/adjacencies.tsv", sep="\t")
adj.groupby("TF").size().rename("n_candidate_targets").sort_values(ascending=False).to_csv(f"{R}/adjacencies_tf_summary.csv")
adj_demo = adj[adj.TF.isin(TEACH_TFS)].sort_values(["TF", "importance"], ascending=[True, False])
adj_demo.to_csv(f"{R}/adjacencies_demo.tsv", sep="\t", index=False, float_format="%.4f")
rep["adjacencies_demo_rows"] = len(adj_demo)

# --- Paso 24/25: AUCell + metadatos por ID -----------------------------------
auc = pd.read_csv(f"{R}/auc_mtx.csv", index_col=0)
meta = pd.read_csv("data/metadata_demo.csv", index_col="cell_id")
expr_ids = pd.read_csv("data/expression_demo.tsv", sep="\t", usecols=[0]).iloc[:, 0]
rep["auc_shape"] = list(auc.shape)
rep["auc_cells_eq_expression_cells_in_order"] = bool(list(auc.index) == list(expr_ids))
rep["auc_cols_eq_regulons"] = bool(set(auc.columns) == set(sizes.index))
rep["auc_nan"] = int(auc.isna().sum().sum())
assert set(auc.index) == set(meta.index) and not auc.index.duplicated().any()
joined = meta[["population"]].join(auc, how="inner")          # unión POR ID, nunca por posición
assert len(joined) == len(meta)
med = joined.groupby("population").median().T
med.to_csv(f"{R}/auc_median_by_population.csv", float_format="%.5f")
rss = regulon_specificity_scores(auc.loc[meta.index], meta["population"]).T
rss.to_csv(f"{R}/rss_by_population.csv", float_format="%.5f")
rep["top5_rss_by_population"] = {p: rss[p].sort_values(ascending=False).head(5).round(3).to_dict() for p in rss.columns}
# Diferencia de mediana: población vs. máximo de las otras dos (para elegir regulones visualmente claros).
delta = pd.DataFrame({p: med[p] - med.drop(columns=p).max(axis=1) for p in med.columns})
rep["top5_median_delta_by_population"] = {p: delta[p].sort_values(ascending=False).head(5).round(4).to_dict() for p in delta.columns}

out = auc.loc[expr_ids.values].copy()
out.index.name = "cell_id"
out.to_csv(f"{R}/auc_demo.csv")                                 # mismos valores, mismo orden que la matriz
chk = pd.read_csv(f"{R}/auc_demo.csv", index_col=0)
rep["auc_demo_max_abs_diff_vs_auc_mtx"] = float(np.abs(chk.values - auc.loc[chk.index, chk.columns].values).max())

json.dump(rep, open("logs/25_postprocess_report.json", "w"), indent=2, ensure_ascii=False, default=str)
print(json.dumps(rep, indent=2, ensure_ascii=False, default=str))
