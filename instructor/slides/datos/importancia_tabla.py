# Importancia de GRNBoost2 como suma de la baja del error, árbol por árbol (NKG7, semilla 2026).
import json, numpy as np, pandas as pd, sklearn, arboreto
from arboreto.core import to_tf_matrix, fit_model, to_feature_importances, SGBM_KWARGS, clean
ex = pd.read_csv('/w/data/expression_demo.tsv', sep='\t', index_col=0)
tfs = [l.strip() for l in open('/w/resources/allTFs_hg38.txt') if l.strip()]
tfm, tfn = to_tf_matrix(ex.to_numpy(), list(ex.columns), tfs)
y = ex['NKG7'].to_numpy(); X, names = clean(tfm, tfn, 'NKG7')
reg = fit_model('GBM', SGBM_KWARGS, X, y, 25, 2026)
imp = to_feature_importances('GBM', SGBM_KWARGS, reg)
# baja del error (MSE ponderado por células, sin normalizar) que cada TF aporta en cada árbol
D = np.array([e[0].tree_.compute_feature_importances(normalize=False) for e in reg.estimators_])  # árboles x TF
n = D.shape[0]; tot_tree = D.sum(1); mean_tot = float(tot_tree.mean())
k = names.index('TBX21')
rows = [(i + 1, round(float(D[i, k]), 3)) for i in range(n) if D[i, k] > 0]
order = np.argsort(-imp)[:8]
top = [{'tf': names[j], 'n_trees': int((D[:, j] > 0).sum()), 'sum_drop': round(float(D[:, j].sum()), 2),
        'imp': round(float(imp[j]), 2), 'check': round(float(D[:, j].sum() / mean_tot), 2)} for j in order]
out = {'versions': {'sklearn': sklearn.__version__, 'arboreto': getattr(arboreto, '__version__', '?')},
       'n_trees': n, 'mean_drop_per_tree': round(mean_tot, 3), 'total_drop_all_tf': round(float(tot_tree.sum()), 2),
       'tbx21_rows': rows, 'tbx21_n': len(rows), 'tbx21_sum': round(float(D[:, k].sum()), 2),
       'tbx21_imp': round(float(imp[k]), 2), 'top8': top,
       'max_abs_diff_check': float(np.max(np.abs(D.sum(0) / mean_tot - imp)))}
json.dump(out, open('/s/importancia_tabla.json', 'w'), indent=1, ensure_ascii=False)
print(json.dumps({k2: v for k2, v in out.items() if k2 != 'tbx21_rows'}, ensure_ascii=False, indent=1)); print(rows[:12], rows[-3:])
