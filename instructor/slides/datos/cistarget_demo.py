import json, numpy as np, pandas as pd, pyarrow.feather as pf
from pyscenic.utils import modules_from_adjacencies
from ctxcore.rnkdb import FeatherRankingDatabase
from ctxcore.recovery import aucs as calc_aucs, recovery, leading_edge, derive_rank_cutoff
R = '/w/resources/'
DB500 = R + 'hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
ex = pd.read_csv('/w/data/expression_demo.tsv', sep='\t', index_col=0)
adj = pd.read_csv('/w/results/adjacencies.tsv', sep='\t')
mods = modules_from_adjacencies(adj, ex)
print('modulos', len(mods))
db = FeatherRankingDatabase(fname=DB500, name='hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings')
out = {'total_genes': db.total_genes, 'rank_cutoff': int(derive_rank_cutoff(0.05, db.total_genes))}
full = None
def analiza(tf, motif, ctxkey):
    global full
    m = [x for x in mods if x.transcription_factor == tf and ctxkey in x.context and 'activating' in x.context]
    assert len(m) == 1, [x.context for x in m]
    m = m[0]
    df = db.load(m)
    genes = df.columns.values; w = np.ones(len(genes))
    a = calc_aucs(df, db.total_genes, w, 0.05); nes = (a - a.mean()) / a.std()
    i = list(df.index).index(motif)
    rccs, _ = recovery(df, db.total_genes, w, 5000, 0.05, no_auc=True)
    avg = rccs.mean(axis=0); avg2 = avg + 2 * rccs.std(axis=0)
    le, ram = leading_edge(rccs[i], avg2, df.values[i], genes, None)
    ranks = df.values[i]
    d = {'module_genes': len(m.genes), 'in_db': len(genes), 'auc': float(a[i]), 'nes': float(nes[i]), 'auc_mean': float(a.mean()), 'auc_sd': float(a.std()),
         'n_motifs': len(a), 'n_motifs_nes3': int((nes >= 3).sum()), 'motif_rank_by_nes': int((nes > nes[i]).sum()) + 1,
         'rank_at_max': int(ram), 'leading_edge_n': len(le), 'leading_edge': sorted(le, key=lambda t: t[1]),
         'n_in_top_cutoff': int((ranks <= out['rank_cutoff']).sum()), 'expected_in_top_cutoff': len(genes) * (out['rank_cutoff'] + 1) / db.total_genes,
         'curve': [int(x) for x in rccs[i][::25]], 'avg': [round(float(x), 2) for x in avg[::25]], 'avg2': [round(float(x), 2) for x in avg2[::25]],
         'auc_hist': np.histogram(a, bins=40)[0].tolist(), 'auc_edges': np.histogram(a, bins=40)[1].round(5).tolist(),
         'module_gene_ranks_sorted': sorted(int(r) for r in ranks)}
    if full is None:
        full = pf.read_table(DB500).to_pandas().set_index('motifs')
    row = full.loc[motif].sort_values()
    d['motif_top20'] = [(g, int(r), g in set(genes)) for g, r in row.head(20).items()]
    d['known'] = {g: int(row[g]) for g in ['LYZ', 'FCN1', 'S100A8', 'CD14', 'CSF1R', 'ITGAM', 'TYROBP', 'FCER1G', 'CD79A', 'NKG7', 'GAPDH', 'ACTB'] if g in row.index}
    return d
for tf, motif, ck in [('SPI1', 'metacluster_167.5', 'top5perTarget'), ('SPI1', 'metacluster_167.5', 'weight>90.0%'), ('EBF1', 'tfdimers__MD00215', 'weight>90.0%')]:
    k = f'{tf}|{motif}|{ck}'; out[k] = analiza(tf, motif, ck)
    v = out[k]; print(k, {x: v[x] for x in ['module_genes', 'in_db', 'auc', 'nes', 'rank_at_max', 'leading_edge_n', 'n_in_top_cutoff', 'expected_in_top_cutoff', 'n_motifs', 'n_motifs_nes3', 'motif_rank_by_nes']})
    print('  top20', v['motif_top20']); print('  known', v['known']); print('  LE', [g for g, _ in v['leading_edge']][:40])
json.dump(out, open('/s/cistarget_demo.json', 'w'), ensure_ascii=False, default=lambda o: o.item() if hasattr(o,'item') else str(o))
