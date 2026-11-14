"""Paso 46: calcula, con los mismos archivos que usa el notebook, todos los números que
aparecen en las respuestas esperadas (SOLUCIONES.md). Nada de las soluciones se escribe
a mano: cada cifra sale de aquí.

Uso (cualquier Python con pandas/numpy):  python scripts/46_answers.py > logs/46_answers.log
"""
from pathlib import Path
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
pd.set_option("display.width", 200, "display.max_columns", 20)
ORDEN = ["Mono_CD14", "B", "NK"]

meta = pd.read_csv(R / "data/metadata_demo.csv", index_col="cell_id")
expr = pd.read_csv(R / "data/expression_demo.tsv.gz", sep="\t", index_col="cell_id")
adj = pd.read_csv(R / "results/adjacencies_demo.tsv", sep="\t")
res = pd.read_csv(R / "results/adjacencies_tf_summary.csv", index_col="TF")
reg = pd.read_csv(R / "results/regulons_demo.csv")
mot = pd.read_csv(R / "results/motifs_summary.csv")
auc = pd.read_csv(R / "results/auc_demo.csv", index_col="cell_id")
datos = meta[["population"]].join(auc, how="inner")
assert len(datos) == 900


def h(t):
    print(f"\n{'=' * 8} {t} {'=' * 8}")


h("Sección 0")
print(meta.population.value_counts().to_dict())

h("Ej1 · TBX21")
d = adj[adj.TF == "TBX21"].sort_values("importance", ascending=False)
print("candidatos TBX21:", len(d))
print("top15:", d.target.head(15).tolist())
for g in ["NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "KLRD1", "FGFBP2", "CST7", "CTSW"]:
    pos = np.where(d.target.values == g)[0]
    print(f"  rango de {g}: {pos[0] + 1 if len(pos) else '—'}")
q = d.importance.quantile([0.5, 0.9, 0.99]).round(3).to_dict()
print("cuantiles importancia:", q, "max", round(d.importance.max(), 2))
print("dianas con importancia > 1:", int((d.importance > 1).sum()))

h("Ej1 · RPS4X")
r = adj[adj.TF == "RPS4X"].sort_values("importance", ascending=False)
print("candidatos RPS4X:", len(r))
print(r.head(10).to_string(index=False))
top100 = r.target.head(100)
print("fracción ribosomal (RPL/RPS) en top 100:", top100.str.match(r"^RP[LS]\d").mean())
print("top3 global:\n", adj.sort_values("importance", ascending=False).head(3).to_string(index=False))
print("TF con más candidatos (tabla completa):\n", res.head(10).to_string())

h("Ej2 · regulones")
tam = reg.groupby("regulon").size().sort_values(ascending=False)
print("n regulones:", len(tam), "mediana", tam.median(), "min", tam.min(), "max", tam.max(), tam.idxmax())
print("todos (+):", reg.regulon.str.endswith("(+)").all())
cand = set(adj.loc[adj.TF == "TBX21", "target"]); rg = set(reg.loc[reg.TF == "TBX21", "target"])
print(f"TBX21: {len(cand)} -> {len(rg)} ({100 * len(rg) / len(cand):.1f}%); subconjunto: {rg <= cand}")
print("TBX21 regulón top 10:", reg[reg.TF == "TBX21"].sort_values("weight", ascending=False).target.head(10).tolist())
print(mot[mot.TF == "TBX21"].sort_values("NES", ascending=False).head(5).to_string(index=False))
m = reg[reg.regulon == "MAFB(+)"].sort_values("weight", ascending=False)
print("MAFB(+):", len(m), "genes; top 15:", m.target.head(15).tolist(), "; autorregulación:", "MAFB" in set(m.target))
tf_taller = adj.TF.unique()
tabla = pd.DataFrame({"candidatos": adj.groupby("TF").size(),
                      "regulon": reg.groupby("TF").size().reindex(tf_taller)}).fillna(0).astype(int)
tabla["pct"] = (100 * tabla.regulon / tabla.candidatos).round(1)
print(tabla.sort_values("candidatos", ascending=False).to_string())
print("PAX5 top10:", adj[adj.TF == "PAX5"].nlargest(10, "importance").target.tolist())
print("motivos RPS4X/PAX5 en tabla de motivos enriquecidos:",
      len(mot[mot.TF == "RPS4X"]), len(mot[mot.TF == "PAX5"]))
print("regulones con TF del linaje B en la tabla:", sorted(set(reg.TF) & {"PAX5", "EBF1", "SPIB", "POU2AF1", "IRF8", "BCL11A", "TCF4"}))

h("Ej3 · medianas AUC")
REG5 = ["MAFB(+)", "EBF1(+)", "TBX21(+)", "FOS(+)", "SPI1(+)"]
med = datos.groupby("population")[REG5].median().loc[ORDEN].T
print(med.round(4).to_string())
iqr = datos.groupby("population")[REG5].quantile([0.25, 0.75]).unstack().round(4)
print("cuartiles:\n", iqr.T.to_string())
for rgn in REG5:
    print(f"  tamaño {rgn}: {tam[rgn]}")
expr_norm = np.log1p(expr.div(meta.loc[expr.index, "n_counts"], axis=0) * 1e4)
tfs = [x[:-3] for x in REG5]
print("expresión media log CP10k:\n", expr_norm[tfs].groupby(meta.population).mean().loc[ORDEN].T.round(3).to_string())
print("% células con conteo>0:\n",
      (expr[tfs] > 0).groupby(meta.population).mean().loc[ORDEN].T.mul(100).round(1).to_string())
ieg = [g for g in ["FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "IER2", "DUSP1", "HSPA1A", "HSPA1B"] if g in expr.columns]
print("genes de respuesta temprana, media log CP10k:\n",
      expr_norm[ieg].groupby(meta.population).mean().loc[ORDEN].T.round(2).to_string())
f = reg[reg.regulon == "FOS(+)"].sort_values("weight", ascending=False)
print("FOS(+) top 15:", f.target.head(15).tolist())
print("FOS(+) contiene IEG:", sorted(set(f.target) & set(ieg + ["IER3", "NR4A1", "NR4A2", "KLF6", "ZFP36", "ATF3"])))
# rango por población (max-min de medianas) vs tamaño
allmed = datos.groupby("population").median().T[ORDEN]
rango = (allmed.max(axis=1) - allmed.min(axis=1))
print("rango de medianas (max-min) por regulón; ranking dentro de 272:")
rk = rango.rank(ascending=False)
for rgn in REG5:
    print(f"  {rgn}: rango {rango[rgn]:.4f}, puesto {int(rk[rgn])}/272")
print("corr Spearman log10(tamaño) vs rango:",
      round(tam.rank().corr(rango.reindex(tam.index).rank()), 3))   # Spearman sin scipy
# efecto en tamaño: Cohen-like = diferencia de medianas / sd global
for rgn in REG5:
    s = datos[rgn]
    print(f"  {rgn}: rango/sd = {rango[rgn] / s.std():.2f}")

h("Ej3 · RSS")
rss = pd.read_csv(R / "results/rss_by_population.csv", index_col=0)
print(rss.loc[REG5].round(3).to_string())
for p in ORDEN:
    print(p, "top5 RSS:", rss[p].nlargest(5).round(3).to_dict())

h("Ej5 · especificidad (lo que ve el alumno)")
esp = pd.DataFrame({p: allmed[p] - allmed.drop(columns=p).max(axis=1) for p in ORDEN})
print(pd.DataFrame({p: esp[p].nlargest(8).index for p in ORDEN}).to_string())
print(pd.DataFrame({p: esp[p].nlargest(8).round(3).values for p in ORDEN}).to_string())
cands = sorted({r for p in ORDEN for r in esp[p].nlargest(8).index})
print("Fichas de los candidatos de la lista:")
for rgn in cands + ["EBF1(+)"]:
    tf = rgn.split("(")[0]
    rr = reg[reg.regulon == rgn].sort_values("weight", ascending=False)
    mm = mot[mot.TF == tf].sort_values("NES", ascending=False).head(1)
    e = expr_norm[tf].groupby(meta.population).mean().loc[ORDEN].round(2).tolist() if tf in expr.columns else None
    print(f"  {rgn:<11} n={len(rr):<5} medianas={allmed.loc[rgn].round(3).tolist()} expr={e} "
          f"motivo={mm.MotifID.iloc[0] if len(mm) else '—'} NES={mm.NES.iloc[0] if len(mm) else float('nan'):.2f} "
          f"top8={rr.target.head(8).tolist()}")

h("Sección 4.1 · recuperación")
print("(ver logs/35_validate_numpy_aucell.log para el ejemplo MAFB monocito vs B)")
