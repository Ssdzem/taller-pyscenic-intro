"""Paso 51 (prueba de viabilidad): matrices para correr el CLI de pySCENIC en un entorno tipo Colab.
Escribe en el directorio actual: expr_full.tsv (900 x 10,120), expr_r1.tsv (2,000 genes más variables
+ todos los TF presentes) y tfs.txt (TF de allTFs_hg38 presentes en la matriz)."""
import sys
import numpy as np
import pandas as pd

expr_gz, tf_list = sys.argv[1], sys.argv[2]
expr = pd.read_csv(expr_gz, sep="\t", index_col=0)
tfs = [t.strip() for t in open(tf_list) if t.strip() in expr.columns]
lognorm = np.log1p(expr.div(expr.sum(axis=1), axis=0) * 1e4)
hvg = lognorm.var().sort_values(ascending=False).index[:2000]
genes = [g for g in expr.columns if g in set(hvg) | set(tfs)]
expr.to_csv("expr_full.tsv", sep="\t")
expr[genes].to_csv("expr_r1.tsv", sep="\t")
open("tfs.txt", "w").write("\n".join(tfs) + "\n")
print(f"full {expr.shape}; r1 {expr[genes].shape}; TF presentes {len(tfs)}")
