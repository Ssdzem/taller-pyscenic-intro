"""Valida el AUCell reimplementado en numpy/pandas (el del notebook, nivel 2) contra results/auc_demo.csv de pySCENIC 0.12.1."""
import time, numpy as np, pandas as pd

def aucell_numpy(expr, regulons, auc_threshold=0.05, seed=2026):
    """Réplica de pySCENIC 0.12.1 (CLI por defecto: sin pesos)."""
    # 1) ranking por célula: mayor expresión = rango 0; empates rotos al azar (semilla), como pySCENIC
    rnk = expr.sample(frac=1.0, axis=1, random_state=seed).rank(axis=1, ascending=False, method="first") - 1
    rnk = rnk[expr.columns].to_numpy()
    cutoff = int(round(auc_threshold * expr.shape[1])) - 1   # igual que ctxcore.derive_rank_cutoff
    col = {g: i for i, g in enumerate(expr.columns)}
    out = {}
    for name, genes in regulons.items():
        idx = [col[g] for g in genes if g in col]
        if len(idx) < 0.8 * len(genes):
            out[name] = np.zeros(len(expr)); continue
        r = rnk[:, idx]
        # 2) cada gen del regulón en el top `cutoff` suma (cutoff - rango); 3) normalizar por el máximo posible
        out[name] = np.clip(cutoff - r, 0, None).sum(axis=1) / ((cutoff + 1) * len(idx))
    return pd.DataFrame(out, index=expr.index)

t = time.time()
expr = pd.read_csv("data/expression_demo.tsv.gz", sep="\t", index_col="cell_id")
regs = pd.read_csv("results/regulons_demo.csv").groupby("regulon")["target"].apply(list).to_dict()
ref = pd.read_csv("results/auc_demo.csv", index_col="cell_id")
t1 = time.time()
mine = aucell_numpy(expr, regs)[ref.columns].loc[ref.index]
t2 = time.time()
d = np.abs(mine.values - ref.values)
corr = np.mean([np.corrcoef(mine[c], ref[c])[0, 1] for c in ref.columns if ref[c].std() > 0])
print(f"pandas {pd.__version__} numpy {np.__version__} | carga {t1-t:.1f}s | AUCell numpy {t2-t1:.1f}s")
print(f"max |dif| = {d.max():.3e} | media |dif| = {d.mean():.3e} | correlación media por regulón = {corr:.6f}")
print("IDENTICO (<=1e-12):", d.max() <= 1e-12)
