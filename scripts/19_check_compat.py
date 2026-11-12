"""Paso 19: compatibilidad de especie/nomenclatura entre matriz demo, lista de TF, bases ranking y tabla de motivos."""
import pandas as pd
from ctxcore.rnkdb import opendb
genes = set(pd.read_csv("data/expression_demo.tsv", sep="\t", nrows=0, index_col=0).columns)
tfs = set(open("resources/allTFs_hg38.txt").read().split())
mot = pd.read_csv("resources/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", sep="\t", usecols=["#motif_id", "gene_name"])
print(f"matriz: {len(genes)} genes | lista TF: {len(tfs)} | TF en matriz: {len(tfs & genes)}")
print(f"tabla motivos: {len(mot)} filas, {mot['#motif_id'].nunique()} motivos, {mot.gene_name.nunique()} TFs anotados; TF de la lista con motivo: {len(tfs & set(mot.gene_name))}")
for f in ["hg38_10kbp_up_10kbp_down_full_tx_v10_clust", "hg38_500bp_up_100bp_down_full_tx_v10_clust"]:
    db = opendb(f"resources/{f}.genes_vs_motifs.rankings.feather", name=f)
    g = set(db.genes)
    print(f"{f}: {len(g)} genes en DB, {db.total_genes if hasattr(db,'total_genes') else ''} | genes matriz en DB: {len(genes & g)} ({100*len(genes & g)/len(genes):.1f}%) | TF-en-matriz en DB: {len(tfs & genes & g)}")
key = ["SPI1","CEBPB","CEBPA","IRF8","PAX5","EBF1","SPIB","POU2AF1","TBX21","EOMES","RUNX3","PRDM1"]
print("TF clave -> en lista TF / en matriz / con motivo:", {k: (k in tfs, k in genes, k in set(mot.gene_name)) for k in key})
