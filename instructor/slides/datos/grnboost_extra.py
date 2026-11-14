import json, numpy as np, pandas as pd
from arboreto.core import to_tf_matrix, fit_model, to_feature_importances, SGBM_KWARGS, clean
ex = pd.read_csv('/w/data/expression_demo.tsv', sep='\t', index_col=0)
meta = pd.read_csv('/w/data/metadata_demo.csv', index_col=0).loc[ex.index]
tfs = [l.strip() for l in open('/w/resources/allTFs_hg38.txt') if l.strip()]
tfm, tfn = to_tf_matrix(ex.to_numpy(), list(ex.columns), tfs)
y = ex['NKG7'].to_numpy(); X, names = clean(tfm, tfn, 'NKG7')
reg = fit_model('GBM', SGBM_KWARGS, X, y, 25, 2026)
imp = to_feature_importances('GBM', SGBM_KWARGS, reg)
k = names.index('TBX21')
used = sum(1 for e in reg.estimators_ if k in set(e[0].tree_.feature))
staged = np.array([p for p in reg.staged_predict(X)])  # trees x cells
nk = np.where((meta['population'].values == 'NK') & (y >= 30) & (y <= 40))[0][0]
mono = np.where((meta['population'].values == 'Mono_CD14') & (y == 0))[0][0]
pts = [0, 1, 10, 50, 100, 200, 300, 400, len(reg.estimators_)]
def traj(c): return {p: round(float(reg.init_.predict(X[[c]])[0] if p == 0 else staged[p - 1, c]), 2) for p in pts}
mse = [float(np.mean((y - staged[i]) ** 2)) for i in range(len(staged))]
out = {'sum_importances': float(imp.sum()), 'n_trees': len(reg.estimators_), 'tbx21_imp': float(imp[k]), 'tbx21_frac': float(imp[k] / imp.sum()),
       'trees_using_tbx21': used, 'nk_cell': {'id': ex.index[nk], 'real': int(y[nk]), 'traj': traj(nk)}, 'mono_cell': {'id': ex.index[mono], 'real': int(y[mono]), 'traj': traj(mono)},
       'mean_y': float(y.mean()), 'mse_curve_every10': [round(m, 1) for m in mse[::10]] + [round(mse[-1], 1)],
       'oob_improvement_last30': [round(float(v), 4) for v in reg.oob_improvement_[-30:]],
       'max_features_n': int(reg.max_features_), 'subsample': reg.subsample, 'lr': reg.learning_rate, 'max_depth': reg.max_depth}
json.dump(out, open('/s/grnboost_extra.json', 'w'), indent=1, ensure_ascii=False); print(json.dumps(out, ensure_ascii=False)[:3000])
