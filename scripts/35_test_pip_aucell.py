"""Nivel 1 del notebook: AUCell con pySCENIC instalado por pip, comparado con results/auc_demo.csv."""
import time, numpy as np, pandas as pd
t = time.time()
from ctxcore.genesig import Regulon
from pyscenic.aucell import aucell
expr = pd.read_csv("data/expression_demo.tsv.gz", sep="\t", index_col="cell_id")
regs = pd.read_csv("results/regulons_demo.csv")
sigs = [Regulon(name=n, gene2weight={g: 1.0 for g in d.target}, transcription_factor=d.TF.iloc[0], gene2occurrence=[])
        for n, d in regs.groupby("regulon")]
for nw in (1, 2):
    t1 = time.time()
    auc = aucell(expr, sigs, auc_threshold=0.05, noweights=True, seed=2026, num_workers=nw)
    ref = pd.read_csv("results/auc_demo.csv", index_col="cell_id")
    d = np.abs(auc.loc[ref.index, ref.columns].values - ref.values).max()
    print(f"num_workers={nw}: {time.time()-t1:.1f}s | forma {auc.shape} | max|dif| vs precalculado = {d:.2e}")
print(f"numpy {np.__version__} pandas {pd.__version__} | total {time.time()-t:.1f}s")
