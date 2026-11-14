import json, numpy as np, pandas as pd
from arboreto.core import to_tf_matrix, fit_model, to_feature_importances, SGBM_KWARGS, clean
from sklearn.tree import DecisionTreeRegressor
ex = pd.read_csv('/w/data/expression_demo.tsv', sep='\t', index_col=0)
meta = pd.read_csv('/w/data/metadata_demo.csv', index_col=0).loc[ex.index]
tfs = [l.strip() for l in open('/w/resources/allTFs_hg38.txt') if l.strip()]
genes = list(ex.columns)
M = ex.to_numpy()
tfm, tfn = to_tf_matrix(M, genes, tfs)
out = {'n_tf': len(tfn)}
for target in ['NKG7', 'CD79A']:
    y = ex[target].to_numpy()
    X, names = clean(tfm, tfn, target)
    reg = fit_model('GBM', SGBM_KWARGS, X, y, 25, 2026)
    imp = to_feature_importances('GBM', SGBM_KWARGS, reg)
    top = sorted(zip(names, imp), key=lambda t: -t[1])[:8]
    staged = [float(np.mean((y - p) ** 2)) for p in reg.staged_predict(X)]
    roots = [names[est[0].tree_.feature[0]] for est in reg.estimators_[:5]]
    d = {'n_trees': len(reg.estimators_), 'top': [(n, round(float(v), 2)) for n, v in top],
         'mse_0': float(np.mean((y - y.mean()) ** 2)), 'mse_by_tree': {k: staged[k - 1] for k in [1, 10, 50, 100, 200, 400, len(staged)] if k <= len(staged)},
         'roots_first5': roots, 'mean_by_pop': ex[target].groupby(meta['population']).mean().round(2).to_dict(),
         'first_tree_depth': int(reg.estimators_[0][0].get_depth())}
    # árbol único (tocón) con todos los TF
    st = DecisionTreeRegressor(max_depth=1, random_state=0).fit(X, y)
    t = st.tree_
    f = names[t.feature[0]]; thr = float(t.threshold[0])
    left = X[:, t.feature[0]] <= thr
    d['stump'] = {'tf': f, 'thr': thr, 'n_left': int(left.sum()), 'n_right': int((~left).sum()),
                  'mean_left': float(y[left].mean()), 'mean_right': float(y[~left].mean()),
                  'mse_after': float(np.mean((y - st.predict(X)) ** 2)),
                  'pop_right': meta['population'][~left].value_counts().to_dict(),
                  'pop_left': meta['population'][left].value_counts().to_dict()}
    # tocón fijo: TBX21 > 0
    if target == 'NKG7':
        g = ex['TBX21'].to_numpy() > 0
        d['tbx21_split'] = {'n_yes': int(g.sum()), 'n_no': int((~g).sum()), 'mean_yes': float(y[g].mean()), 'mean_no': float(y[~g].mean()),
                            'pop_yes': meta['population'][g].value_counts().to_dict()}
        d2 = DecisionTreeRegressor(max_depth=2, random_state=0).fit(X, y).tree_
        d['depth2'] = [(names[d2.feature[i]] if d2.feature[i] >= 0 else None, float(d2.threshold[i]), int(d2.n_node_samples[i]), float(d2.value[i][0][0])) for i in range(d2.node_count)]
    out[target] = d
adj = pd.read_csv('/w/results/adjacencies.tsv', sep='\t')
for target in ['NKG7', 'CD79A']:
    out[target]['server_top'] = adj[adj.target == target].nlargest(8, 'importance')[['TF', 'importance']].round(2).values.tolist()
json.dump(out, open('/s/grnboost_demo.json', 'w'), indent=1, ensure_ascii=False)
print(json.dumps(out, indent=1, ensure_ascii=False))
