"""Paso 26: tres figuras de referencia a partir de results/ y data/ (todos del mismo análisis de 900 células).

  figures/red_regulon.png   red TF→dianas de MAFB(+) (todas sus dianas; grosor = importancia GRNBoost2)
  figures/heatmap.png       AUCell por célula (z-score por regulón), top-5 RSS por población + regulones didácticos
  figures/tsne.png          t-SNE de 10x coloreado por población y por actividad de MAFB(+), EBF1(+), TBX21(+)

Corre en localhost/taller-celltypist:1.6.3 (matplotlib/seaborn/networkx recientes, parecido a Colab):
  podman run --rm -v $PWD:/w:Z -w /w localhost/taller-celltypist:1.6.3 python -I scripts/26_figures.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

POP_ORDER = ["Mono_CD14", "B", "NK"]
POP_LABEL = {"Mono_CD14": "Monocitos CD14+", "B": "Células B", "NK": "Células NK"}
POP_COLOR = {"Mono_CD14": "#d95f02", "B": "#1b9e77", "NK": "#7570b3"}
SHOW = ["MAFB(+)", "EBF1(+)", "TBX21(+)"]
NET_REGULON = "MAFB(+)"

meta = pd.read_csv("data/metadata_demo.csv", index_col="cell_id")
auc = pd.read_csv("results/auc_demo.csv", index_col="cell_id")
regs = pd.read_csv("results/regulons_demo.csv")
rss = pd.read_csv("results/rss_by_population.csv", index_col=0)
assert set(auc.index) == set(meta.index)
df = meta.join(auc, how="inner")                      # unión por ID celular

def save(fig, name):
    fig.savefig(f"figures/{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"figures/{name}.pdf", bbox_inches="tight")
    plt.close(fig)

# 1) Red de un regulón ---------------------------------------------------------
r = regs[regs.regulon == NET_REGULON].sort_values("weight", ascending=False)
tf = r.TF.iloc[0]
self_loop = tf in set(r.target)                       # el TF puede estar entre sus propias dianas
targets = [t for t in r.target if t != tf]
G = nx.DiGraph()
for _, row in r[r.target != tf].iterrows():
    G.add_edge(tf, row.target, weight=row.weight)
ang = np.linspace(0, 2 * np.pi, len(targets), endpoint=False)
pos = {tf: np.array([0.0, 0.0]), **{t: np.array([np.cos(a), np.sin(a)]) for t, a in zip(targets, ang)}}
w = np.array([G[u][v]["weight"] for u, v in G.edges()])
fig, ax = plt.subplots(figsize=(8, 8))
nx.draw_networkx_edges(G, pos, ax=ax, width=0.5 + 3.5 * (w - w.min()) / (np.ptp(w) or 1), edge_color="#999999",
                       arrows=True, arrowsize=8, alpha=0.8)
nx.draw_networkx_nodes(G, pos, nodelist=[tf], node_color=POP_COLOR["Mono_CD14"], node_size=1600, ax=ax)
nx.draw_networkx_nodes(G, pos, nodelist=targets, node_color="#fde0c5", node_size=500, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
loop_txt = f" (incluye a {tf} mismo: autorregulación)" if self_loop else ""
ax.set_title(f"Regulón {NET_REGULON}: {len(r)} genes tras cisTarget{loop_txt}\n"
             "grosor de arista = importancia GRNBoost2 (asociación, no causalidad)", fontsize=11)
ax.axis("off")
save(fig, "red_regulon")

# 2) Heatmap ------------------------------------------------------------------
top = []
for p in POP_ORDER:
    top += [x for x in rss[p].sort_values(ascending=False).index[:5] if x not in top]
top += [x for x in SHOW + ["FOS(+)"] if x not in top]
cells = df.sort_values("population", key=lambda s: s.map({p: i for i, p in enumerate(POP_ORDER)})).index
z = df.loc[cells, top]
z = ((z - z.mean()) / z.std(ddof=0)).clip(-2.5, 2.5).T
fig = plt.figure(figsize=(11, 6.5))
gs = fig.add_gridspec(2, 2, height_ratios=[0.4, 10], width_ratios=[40, 1], hspace=0.03, wspace=0.03)
axc, axh, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
axc.imshow([[list(POP_COLOR).index(p) for p in df.loc[cells, "population"]]], aspect="auto",
           cmap=matplotlib.colors.ListedColormap([POP_COLOR[p] for p in POP_COLOR]))
axc.set_xticks([]); axc.set_yticks([])
for i, p in enumerate(POP_ORDER):
    axc.text(150 + 300 * i, 0, POP_LABEL[p], ha="center", va="center", color="white", fontsize=9, weight="bold")
sns.heatmap(z, ax=axh, cbar_ax=axb, cmap="RdBu_r", center=0, xticklabels=False, yticklabels=True,
            cbar_kws={"label": "AUCell (z-score por regulón)"})
axh.set_xlabel("900 células (300 por población)"); axh.set_ylabel("")
axc.set_title("Actividad de regulones (AUCell): top-5 por especificidad (RSS) en cada población + regulones didácticos", fontsize=10)
save(fig, "heatmap")

# 3) t-SNE --------------------------------------------------------------------
fig, axs = plt.subplots(1, 4, figsize=(18, 4.5))
for p in POP_ORDER:
    d = df[df.population == p]
    axs[0].scatter(d.TSNE1, d.TSNE2, s=5, c=POP_COLOR[p], label=POP_LABEL[p])
axs[0].legend(markerscale=3, fontsize=8, frameon=False); axs[0].set_title("Población")
for ax, reg in zip(axs[1:], SHOW):
    o = df.sort_values(reg)
    sc = ax.scatter(o.TSNE1, o.TSNE2, s=5, c=o[reg], cmap="viridis")
    fig.colorbar(sc, ax=ax, shrink=0.8, label="AUC"); ax.set_title(f"Actividad {reg}")
for ax in axs:
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
fig.suptitle("t-SNE de Cell Ranger (10x, calculado sobre todas las células) — se muestran las 900 células del análisis", fontsize=10)
save(fig, "tsne")
print("figuras:", ", ".join(f"figures/{n}.png/.pdf" for n in ["red_regulon", "heatmap", "tsne"]))
