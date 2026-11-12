"""Paso 13: ¿qué etiquetas existen? Resume clusters graphclust de 10x con marcadores canónicos de PBMC.
No asigna nombres: solo reporta evidencia (media log-normalizada y % de células que expresan)."""
import sys, h5py, numpy as np, pandas as pd, scipy.sparse as sp

h5, clu = sys.argv[1], sys.argv[2]
with h5py.File(h5, "r") as f:
    m = f["matrix"]
    X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:]).T.tocsr()  # células × genes
    bc = m["barcodes"][:].astype(str); genes = m["features"]["name"][:].astype(str)
    ftype = m["features"]["feature_type"][:].astype(str)
print("matriz:", X.shape, "| tipos de feature:", dict(zip(*np.unique(ftype, return_counts=True))),
      "| símbolos duplicados:", pd.Index(genes).duplicated().sum(), "| ¿enteros?:", np.allclose(X.data, np.round(X.data)))
c = pd.read_csv(clu, index_col=0)["Cluster"]
print("barcodes en clusters == matriz:", set(c.index) == set(bc), "| n clusters:", c.nunique())
lib = np.asarray(X.sum(1)).ravel()
L = sp.diags(1e4 / lib) @ X; L.data = np.log1p(L.data)
markers = ["CD3E","CD4","IL7R","CCR7","CD8A","GNLY","NKG7","MS4A1","CD79A","CD14","LYZ","FCGR3A","FCER1A","CLEC10A","LILRA4","PPBP","CD34","HBB"]
gi = {g: np.where(genes == g)[0][0] for g in markers}
clv = c.reindex(bc).values
rows = []
for k in sorted(c.unique()):
    idx = clv == k
    sub = L[idx][:, [gi[g] for g in markers]].toarray()
    rows.append(pd.Series(np.round(sub.mean(0), 2), index=markers, name=f"{k} (n={idx.sum()})"))
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print("\nMedia log1p(CP10k) por cluster graphclust:\n", pd.DataFrame(rows).to_string())
