"""Paso 16: verificación independiente de data/expression_demo.tsv vs data/metadata_demo.csv,
usando el mismo lector que el CLI de pySCENIC 0.12.1."""
import numpy as np, pandas as pd
from pyscenic.cli.utils import load_exp_matrix

ex = load_exp_matrix("data/expression_demo.tsv")          # lector oficial del CLI (células × genes)
meta = pd.read_csv("data/metadata_demo.csv", index_col="cell_id")
ok = lambda c, msg: print(("✅ " if c else "❌ ") + msg) or c
res = [
  ok(ex.shape == (900, 10120), f"forma matriz {ex.shape}"),
  ok(not ex.index.duplicated().any() and not meta.index.duplicated().any(), "IDs celulares sin duplicados"),
  ok(list(ex.index) == list(meta.index), "mismas células y mismo orden en matriz y metadatos"),
  ok(not ex.columns.duplicated().any(), "símbolos génicos únicos"),
  ok(not ex.isna().any().any() and not meta.drop(columns=[]).isna().any().any(), "sin NaN en matriz ni metadatos"),
  ok(np.allclose(ex.values, np.round(ex.values)) and (ex.values >= 0).all(), "conteos crudos enteros ≥0 (no escalados)"),
  ok((ex.sum(1) > 0).all() and (ex.sum(0) > 0).all(), "ninguna célula ni gen vacío"),
  ok(meta.population.value_counts().to_dict() == {"B": 300, "Mono_CD14": 300, "NK": 300}, f"grupos {meta.population.value_counts().to_dict()}"),
]
key = ["SPI1","CEBPB","CEBPA","IRF8","PAX5","EBF1","POU2AF1","SPIB","TBX21","EOMES","RUNX3","PRDM1","CD14","LYZ","MS4A1","CD79A","NKG7","GNLY"]
miss = [g for g in key if g not in ex.columns]
res.append(ok(not miss, f"genes clave presentes (faltan: {miss})"))
cpm = np.log1p(ex.div(ex.sum(1), axis=0) * 1e4)
print("\nMedia log1p(CP10k) de genes clave por población (sanity check biológico):")
print(cpm[[g for g in key if g in ex.columns]].groupby(meta.population).mean().round(2).T.to_string())
print("\nmediana genes/célula:", meta.groupby("population").n_genes.median().to_dict())
print("\nRESULTADO:", "TODO OK" if all(res) else "HAY FALLOS")
