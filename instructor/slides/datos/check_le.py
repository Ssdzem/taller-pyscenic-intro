import ast, numpy as np, pandas as pd
from pyscenic.utils import modules_from_adjacencies, load_motifs
from ctxcore.rnkdb import FeatherRankingDatabase
from ctxcore.recovery import aucs as calc_aucs, recovery, leading_edge
R = '/w/resources/'
DB500 = R + 'hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
ex = pd.read_csv('/w/data/expression_demo.tsv', sep='\t', index_col=0)
adj = pd.read_csv('/w/results/adjacencies.tsv', sep='\t')
adj = adj  # completo: los umbrales dependen de todo
mods = modules_from_adjacencies(adj, ex)
db = FeatherRankingDatabase(fname=DB500, name='x')
mot = load_motifs('/w/results/regulons.csv')
for tf, ck in [('EBF1', 'weight>90.0%'), ('SPI1', 'top5perTarget'), ('SPI1', 'weight>90.0%')]:
    m = [x for x in mods if x.transcription_factor == tf and ck in x.context and 'activating' in x.context][0]
    df = db.load(m); genes = df.columns.values; w = np.ones(len(genes))
    a = calc_aucs(df, db.total_genes, w, 0.05); nes = (a - a.mean()) / a.std()
    rccs, _ = recovery(df, db.total_genes, w, 5000, 0.05, no_auc=True)
    avg2 = rccs.mean(axis=0) + 2 * rccs.std(axis=0)
    sel = mot.loc[tf]
    sel = sel[[ck in c and any('500bp' in s for s in c) for c in sel[('Enrichment', 'Context')]]]
    srv = {mid: (int(r[('Enrichment', 'RankAtMax')]), len(r[('Enrichment', 'TargetGenes')]), set(g for g, _ in r[('Enrichment', 'TargetGenes')])) for mid, r in sel.iterrows()}
    print(tf, ck, 'servidor:', {k: v[:2] for k, v in srv.items()})
    for i in np.where(nes >= 3)[0]:
        le, ram = leading_edge(rccs[i], avg2, df.values[i], genes, None)
        hit = [k for k, v in srv.items() if v[0] == ram and v[2] == set(g for g, _ in le)]
        print('   ', df.index[i], round(float(nes[i]), 3), 'rank_at_max', int(ram), 'LE', len(le), '<- coincide con fila del servidor:' if hit else '', hit)
