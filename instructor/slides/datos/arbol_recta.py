"""Datos para las diapositivas del árbol (recta de TBX21, cortes, árbol de 3 niveles, árboles por gen)."""
import json, time, numpy as np, pandas as pd
from arboreto.core import to_tf_matrix, fit_model, SGBM_KWARGS, clean
from sklearn.tree import DecisionTreeRegressor
ex = pd.read_csv('/w/data/expression_demo.tsv', sep='\t', index_col=0)
meta = pd.read_csv('/w/data/metadata_demo.csv', index_col=0).loc[ex.index]
pop = meta['population'].to_numpy()
tfs = [l.strip() for l in open('/w/resources/allTFs_hg38.txt') if l.strip()]
genes = list(ex.columns)
tfm, tfn = to_tf_matrix(ex.to_numpy(), genes, tfs)
y = ex['NKG7'].to_numpy().astype(float)
X, names = clean(tfm, tfn, 'NKG7')
out = {'n_genes': len(genes), 'n_tf': len(names)}
tb = ex['TBX21'].to_numpy()
out['tbx21_cells'] = [[int(a), p, int(b)] for a, p, b in zip(tb, pop, y)]
out['tbx21_counts'] = pd.crosstab(tb, pop).to_dict(orient='index')
mse0 = float(np.mean((y - y.mean()) ** 2)); out['mse0'] = mse0
def mse_split(g):
    p = np.where(g, y[g].mean(), y[~g].mean()); return float(np.mean((y - p) ** 2))
out['tbx21_cuts'] = [{'cut': c, 'n_yes': int((tb >= c).sum()), 'mean_yes': float(y[tb >= c].mean()), 'mean_no': float(y[tb < c].mean()), 'mse': mse_split(tb >= c)} for c in range(1, int(tb.max()) + 1) if 0 < (tb >= c).sum() < len(tb)]
best = []
for j, n in enumerate(names):
    x = X[:, j]
    if x.max() == x.min(): continue
    s = DecisionTreeRegressor(max_depth=1).fit(x[:, None], y).tree_
    thr = float(s.threshold[0]); g = x > thr
    best.append((n, int(np.floor(thr)) + 1, mse_split(g)))
best.sort(key=lambda r: r[2])
out['best_cut_per_tf_top10'] = best[:10]
out['best_cut_selected'] = [r for r in best if r[0] in ('SPI1', 'EBF1', 'PAX5', 'CEBPB', 'RORA', 'EOMES')]
out['n_tf_tested'] = len(best)
# árbol de 3 niveles con todos los TF y todas las células
t3 = DecisionTreeRegressor(max_depth=3, random_state=0).fit(X, y)
tr = t3.tree_; leaf = t3.apply(X)
nodes = []
for i in range(tr.node_count):
    m = leaf == i if tr.children_left[i] == -1 else None
    nodes.append({'id': i, 'left': int(tr.children_left[i]), 'right': int(tr.children_right[i]),
                  'tf': names[tr.feature[i]] if tr.feature[i] >= 0 else None, 'cut': int(np.floor(tr.threshold[i])) + 1 if tr.feature[i] >= 0 else None,
                  'n': int(tr.n_node_samples[i]), 'mean': float(tr.value[i][0][0]),
                  'pops': pd.Series(pop[m]).value_counts().to_dict() if m is not None else None})
out['tree3'] = nodes; out['tree3_mse'] = float(np.mean((y - t3.predict(X)) ** 2))
# árboles por gen: muestra aleatoria de genes, mismo modelo que grn
rng = np.random.default_rng(2026)
sample = list(rng.choice(genes, 60, replace=False)) + ['NKG7', 'CD79A', 'LYZ', 'GZMB', 'MS4A1', 'S100A8']
nt = {}; t0 = time.time()
for g in sample:
    Xg, _ = clean(tfm, tfn, g)
    reg = fit_model('GBM', SGBM_KWARGS, Xg, ex[g].to_numpy(), 25, 2026)
    nt[g] = len(reg.estimators_)
out['trees_per_gene'] = nt; out['secs'] = time.time() - t0
v = np.array([nt[g] for g in sample[:60]])
out['trees_random60'] = {'min': int(v.min()), 'median': float(np.median(v)), 'max': int(v.max()), 'p10': float(np.percentile(v, 10)), 'p90': float(np.percentile(v, 90))}
json.dump(out, open('/s/arbol_recta.json', 'w'), indent=1, ensure_ascii=False)
o = dict(out); o.pop('tbx21_cells'); print(json.dumps(o, ensure_ascii=False, indent=0)[:6000])
