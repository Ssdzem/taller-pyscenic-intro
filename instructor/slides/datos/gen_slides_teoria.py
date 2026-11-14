"""Genera las diapositivas de teoría de GRNBoost2 y cisTarget con datos reales.

Datos: grnboost_demo.json, grnboost_extra.json y cistarget_demo.json (misma carpeta),
producidos por grnboost_demo.py, grnboost_extra.py y cistarget_demo.py dentro de
docker.io/aertslab/pyscenic:0.12.1 sobre data/expression_demo.tsv y results/.
Uso: python3 instructor/slides/datos/gen_slides_teoria.py
"""
import json
from pathlib import Path

AQUI = Path(__file__).resolve().parent
SLIDES = AQUI.parent / "project" / "slides"
G = json.loads((AQUI / "grnboost_demo.json").read_text())
GX = json.loads((AQUI / "grnboost_extra.json").read_text())
C = json.loads((AQUI / "cistarget_demo.json").read_text())
SP = C["SPI1|metacluster_167.5|weight>90.0%"]
EB = C["EBF1|tfdimers__MD00215|weight>90.0%"]

INK, SOFT, LINE, MONO, B, NK = "#1E2433", "#4A5263", "#CFC7B8", "#B34A00", "#147A5C", "#5E58A0"
PAPER, ACC, TEAL, LILA = "#F7F5F0", "#F4A36B", "#7FCDB2", "#B9B4E6"
FONT = 'font-family="IBM Plex Sans, Arial, sans-serif"'


def fmt(n):
    return f"{n:,}"


def svg(label, w, vw, vh, body, vy=0):
    h = round(w * vh / vw)
    return (f'<svg aria-label="{label}" width="{w}" height="{h}" viewBox="0 {vy} {vw} {vh}" '
            f'xmlns="http://www.w3.org/2000/svg" {FONT}>\n{body}\n</svg>')


def t(x, y, s, size=28, fill=INK, anchor="start", weight=None, extra=""):
    wt = f' font-weight="{weight}"' if weight else ""
    return f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}" fill="{fill}"{wt}{extra}>{s}</text>'


def slide(sid, eyebrow, title, body, notes, cls=""):
    c = f' class="{cls}"' if cls else ""
    return (f'<section id="{sid}"{c} data-transition="fade">\n<div class="head">\n'
            f'<p class="eyebrow">{eyebrow}</p>\n<h2>{title}</h2>\n</div>\n{body}\n<aside>{notes}</aside>\n</section>\n')


def steps(items, style="flex:1;"):
    lis = "\n".join(f"<li>{i}</li>" for i in items)
    return f'<ol class="steplist" style="{style}">\n{lis}\n</ol>'


def figure(s, width, i=0):
    return f'<div class="figure" style="--i:{i}; width:{width}px;">\n{s}\n</div>'


def callout(s, i=2):
    return f'<p class="callout" style="--i:{i}">{s}</p>'


out = {}

# ---------------------------------------------------------------- 1. Árbol: un corte en la recta de TBX21
import random
A = json.loads((AQUI / "arbol_recta.json").read_text())
POPC = {"Mono_CD14": MONO, "B": B, "NK": NK}
cut1 = A["tbx21_cuts"][0]
mse0, mse1, mse3 = A["mse0"], cut1["mse"], A["tree3_mse"]
n_no, n_yes = len(A["tbx21_cells"]) - cut1["n_yes"], cut1["n_yes"]


UX, X0A = 128, 96  # eje X uniforme: 1 UMI = 128 px; el 0 está en x = 96


def xcol(v):
    return X0A + v * UX


def ynk(n):
    return 380 - n * 4.6


rnd = random.Random(2026)
groups = {k: [] for k in ("B", "Mono_CD14", "NK")}  # orden de dibujo: NK encima
for tb, p, nk7 in A["tbx21_cells"]:
    # cada célula se dibuja a la derecha de su valor (entre v y v + 0.9): así el corte cae justo en la marca del 1
    groups[p].append(f'<circle cx="{xcol(tb) + rnd.uniform(8, 0.9 * UX):.0f}" cy="{ynk(nk7):.0f}" r="4"/>')
# un <g> por población: el SVG queda por debajo del límite de 53 KB del formato
dots = [f'<g fill="{POPC[k]}" fill-opacity="0.7">' + "".join(v) + "</g>" for k, v in groups.items()]
cnt = {int(k): sum(v.values()) for k, v in A["tbx21_counts"].items()}
parts = [
    f'<line x1="70" y1="386" x2="990" y2="386" stroke="{INK}" stroke-width="3"/>',
    f'<line x1="70" y1="40" x2="70" y2="386" stroke="{INK}" stroke-width="3"/>',
    t(20, 213, "NKG7 (UMI)", 26, SOFT, "middle", None, ' transform="rotate(-90 20 213)"'),
    *[t(62, ynk(v) + 9, str(v), 24, SOFT, "end") for v in (20, 40, 60)],
    *[f'<line x1="64" y1="{ynk(v)}" x2="70" y2="{ynk(v)}" stroke="{INK}" stroke-width="2"/>' for v in (20, 40, 60)],
    *dots,
    *[f'<line x1="{xcol(v)}" y1="386" x2="{xcol(v)}" y2="396" stroke="{INK}" stroke-width="2"/>' for v in range(8)],
    f'<line x1="{xcol(1)}" y1="40" x2="{xcol(1)}" y2="396" stroke="{INK}" stroke-width="4" stroke-dasharray="12 8"/>',
    t(xcol(1) - 10, 72, "← no", 24, SOFT, "end"), t(xcol(1) + 10, 72, "sí →", 24, SOFT, "start"),
    f'<line x1="74" y1="{ynk(cut1["mean_no"]):.1f}" x2="{xcol(1) - 4}" y2="{ynk(cut1["mean_no"]):.1f}" stroke="{INK}" stroke-width="6"/>',
    f'<line x1="{xcol(1) + 4}" y1="{ynk(cut1["mean_yes"]):.1f}" x2="990" y2="{ynk(cut1["mean_yes"]):.1f}" stroke="{INK}" stroke-width="6"/>',
    f'<rect x="770" y="{ynk(cut1["mean_yes"]) - 50:.1f}" width="220" height="40" rx="8" fill="{PAPER}"/>',
    t(980, ynk(cut1["mean_yes"]) - 20, f"promedio {cut1['mean_yes']:.1f}", 28, INK, "end", 700),
    t(xcol(1) - 40, 26, "corte en 1: ¿TBX21 ≥ 1?", 26, INK, "start", 600),
    *[t(xcol(v), 416, str(v), 26, INK, "middle", 600) for v in range(7)],
    t(990, 446, "TBX21 de cada célula (UMI) →", 26, SOFT, "end"),
    t(74, 446, f"{cnt[0]} células en 0", 24, SOFT, "start"),
    f'<rect x="74" y="466" width="{xcol(1) - 78}" height="76" rx="14" fill="{LINE}"/>',
    t((74 + xcol(1) - 4) / 2, 498, f"no · {n_no}", 26, INK, "middle", 600), t((74 + xcol(1) - 4) / 2, 530, f"predice {cut1['mean_no']:.1f}", 24, INK, "middle"),
    f'<rect x="{xcol(1) + 4}" y="466" width="{986 - xcol(1)}" height="76" rx="14" fill="{NK}"/>',
    t((xcol(1) + 990) / 2, 498, f"sí · {n_yes} células", 26, PAPER, "middle", 600), t((xcol(1) + 990) / 2, 530, f"predice {cut1['mean_yes']:.1f}", 26, PAPER, "middle"),
    t(800, 70, "● Mono", 26, MONO, "start", 600), t(800, 104, "● B", 26, B, "start", 600), t(800, 138, "● NK", 26, NK, "start", 600),
]
fig = svg(f"Cada célula es un punto: a lo largo de la recta, cuánto TBX21 tiene; hacia arriba, cuánto NKG7. El corte entre 0 y 1 UMI de TBX21 deja {n_no} células del lado no, con promedio de NKG7 {cut1['mean_no']:.1f}, y {n_yes} del lado sí, con promedio {cut1['mean_yes']:.1f}.", 1000, 1000, 550, "\n".join(parts))
out["arbol"] = slide(
    "arbol", "Paso ① por dentro · 1 de 6 · mientras corre grn", "Una pregunta = un corte en la recta",
    '<div class="row" style="gap:48px; align-items:center;">\n' + figure(fig, 1048, 0) + "\n" + steps([
        "Meta: adivinar el <b>NKG7</b> de cada célula usando solo los TF.",
        "Cada punto es una célula: a la derecha, más TBX21; arriba, más NKG7.",
        "La pregunta <b>¿TBX21 ≥ 1?</b> es una línea en el 1: a la izquierda las células con 0, a la derecha las que tienen 1 o más.",
        f"Predicción de cada lado = el <b>promedio</b> de NKG7 de sus células: {cut1['mean_no']:.1f} y {cut1['mean_yes']:.1f}.",
        f"Error = promedio de (real − predicción)². Sin cortar, {mse0:.0f}; con el corte, {mse1:.0f}.",
    ]) + "\n</div>",
    f"""Minutos 24–26. Ya lanzaron la preparación y grn (~2 min). Primera de seis diapositivas sobre GRNBoost2. Todo es real: nuestras 900 células, conteos crudos, igual que grn.
La gráfica: el eje horizontal es la recta de TBX21 (UMI de cada célula, de 0 a 6), con escala uniforme. Cada célula se dibuja a la derecha de su valor (las de 1 UMI entre la marca del 1 y la del 2), así la línea del corte cae justo en el 1. La altura de cada punto es su NKG7. Colores = población.
Distribución real de TBX21: en 0 hay 300 monocitos, 295 B y 80 NK; en 1 o más solo hay NK (220) y 5 B.
El corte: la pregunta "¿TBX21 ≥ 1 UMI?" es una línea vertical en la recta. Como los conteos son enteros, cualquier línea entre 0 y 1 separa igual (sklearn guarda el umbral como 0.5); "≥ 1" es la forma fácil de decirlo. Todo lo de la izquierda cae en "no", lo de la derecha en "sí".
La predicción: cada lado predice el promedio de NKG7 de sus células. Lado no: {cut1['mean_no']:.2f}; lado sí: {cut1['mean_yes']:.2f} (las dos líneas horizontales). A una célula nueva solo se le pregunta de qué lado cae: con TBX21 = 1 o con TBX21 = 6 la predicción es la misma, {cut1['mean_yes']:.1f}.
El error: para cada célula, (real − predicción)², promediado. Sin cortar (todas valen el promedio general, {GX['mean_y']:.1f}): {mse0:.0f}. Con el corte: {mse1:.0f} (−{100 * (1 - mse1 / mse0):.0f} %).
Del lado "no" quedan 80 NK con TBX21 = 0 y NKG7 alto (puntos morados arriba a la izquierda): dropout, el gen está pero no se capturó. Por eso el árbol seguirá preguntando (diapositiva 12).""")

# ---------------------------------------------------------------- 1b. ¿Dónde cortar?
cuts = A["tbx21_cuts"]
best = A["best_cut_per_tf_top10"][:6] + [r for r in A["best_cut_selected"] if r[0] in ("SPI1", "EBF1")]
SC = 1.6
parts = [t(0, 30, "TBX21: cada corte posible", 30, INK, "start", 600), t(0, 64, "error después de cortar · menos es mejor", 26, SOFT)]
yy = 90
for c in cuts:
    w = c["mse"] * SC
    parts += [t(130, yy + 30, f"≥ {c['cut']}", 28, INK, "end", 600 if c["cut"] == 1 else None),
              f'<rect x="146" y="{yy + 4}" width="{w:.1f}" height="36" rx="6" fill="{NK if c["cut"] == 1 else LINE}"/>',
              t(146 + w - 12, yy + 31, f"{c['mse']:.0f}", 26, PAPER if c["cut"] == 1 else INK, "end", 600)]
    yy += 50
parts += [f'<line x1="{146 + mse0 * SC:.1f}" y1="84" x2="{146 + mse0 * SC:.1f}" y2="{yy}" stroke="{INK}" stroke-width="3" stroke-dasharray="8 6"/>',
          t(146 + mse0 * SC, yy + 30, f"sin cortar: {mse0:.0f}", 26, INK, "middle")]
X2 = 860
parts += [t(X2, 30, f"El mejor corte de cada TF ({fmt(A['n_tf_tested'])} probados)", 30, INK, "start", 600), t(X2, 64, "los 6 mejores y dos de otras poblaciones", 26, SOFT)]
yy = 90
for tf, c, m in best:
    if tf == "SPI1":
        parts.append(t(X2 + 120, yy + 20, "⋮", 26, SOFT, "middle")); yy += 22
    w = m * SC
    parts += [t(X2 + 200, yy + 28, f"{tf} ≥ {c}", 26, INK, "end", 600 if tf == "TBX21" else None),
              f'<rect x="{X2 + 216}" y="{yy + 4}" width="{w:.1f}" height="32" rx="6" fill="{NK if tf == "TBX21" else LINE}"/>',
              t(X2 + 216 + w - 12, yy + 29, f"{m:.0f}", 26, PAPER if tf == "TBX21" else INK, "end", 600)]
    yy += 40
parts.append(f'<line x1="{X2 + 216 + mse0 * SC:.1f}" y1="84" x2="{X2 + 216 + mse0 * SC:.1f}" y2="{yy}" stroke="{INK}" stroke-width="3" stroke-dasharray="8 6"/>')
yy = max(yy, 90 + 6 * 50 + 40)
fig = svg("Barras reales del error después de cortar. Izquierda: los seis cortes posibles de TBX21; el corte en 1 deja el error más bajo, 121, frente a 277 sin cortar. Derecha: el mejor corte de cada TF; TBX21 gana con 121, luego RORA 145, KLRG1 162, MYBL1 168, MEF2C 174, STAT4 186; SPI1 197 y EBF1 264.", 1616, 1600, yy + 6, "\n".join(parts))
out["corte"] = slide(
    "corte", "Paso ① por dentro · 2 de 6", "¿Dónde cortar? Se prueba todo",
    figure(fig, 1664, 0) + "\n" + callout(f"Gana el corte que deja el error más bajo: <b>TBX21 ≥ 1</b> ({mse0:.0f} → {mse1:.0f}). Esa es la primera pregunta del árbol.", 1),
    f"""Minutos 26–27. Cómo se escoge la pregunta: fuerza bruta.
Izquierda: para TBX21 hay seis cortes posibles (≥ 1, ≥ 2, … ≥ 6). Para cada uno se calculan los dos promedios y el error. El mejor es ≥ 1 (error {mse1:.0f}); ≥ 2 deja del lado "no" a las 102 NK con 1 UMI y el error sube a {cuts[1]['mse']:.0f}.
Derecha: se hace lo mismo con los {fmt(A['n_tf_tested'])} TF y se compara el mejor corte de cada uno. TBX21 gana; luego RORA ({best[1][2]:.0f}), KLRG1, MYBL1, MEF2C, STAT4: todos de NK. SPI1 ({[r for r in best if r[0] == 'SPI1'][0][2]:.0f}) también ayuda, porque separa a los monocitos (sin NKG7). EBF1 casi no baja el error ({[r for r in best if r[0] == 'EBF1'][0][2]:.0f}): está en pocas células B.
La línea punteada es el error sin cortar ({mse0:.0f}).
En GRNBoost2 hay un detalle: en cada pregunta solo se prueban 101 TF al azar, no los 1,018 (lo vemos en boosting).""")

# ---------------------------------------------------------------- 1c. Árbol de 3 niveles
T = {n["id"]: n for n in A["tree3"]}
POPN = {"Mono_CD14": "Mono", "B": "B", "NK": "NK"}


def comp(p):
    s = sorted(p.items(), key=lambda kv: -kv[1])
    if len(s) == 1 or s[0][1] / sum(p.values()) >= 0.93:
        return f"casi todas {POPN[s[0][0]]}" if len(s) > 1 else f"todas {POPN[s[0][0]]}"
    if sum(p.values()) > 100:
        return " y ".join(POPN[k] for k, v in s[:2])
    return " · ".join(f"{v} {POPN[k]}" for k, v in s[:2])


leaves_x = [100 + i * 200 for i in range(8)]
pos, lvl_x = {}, {}


def place(i, lo, hi, depth):
    nd = T[i]
    pos[i] = ((lo + hi) / 2, depth)
    if nd["left"] != -1:
        mid = (lo + hi) / 2
        place(nd["left"], lo, mid, depth + 1)
        place(nd["right"], mid, hi, depth + 1)


place(0, 0, 1600, 0)
YL = [0, 112, 224, 336]
parts = []
for i, nd in T.items():
    if nd["left"] == -1:
        continue
    x, d = pos[i]
    for ch, lab in ((nd["left"], "no"), (nd["right"], "sí")):
        cx, cd = pos[ch]
        parts.append(f'<line x1="{x}" y1="{YL[d] + 56}" x2="{cx}" y2="{YL[cd]}" stroke="{SOFT}" stroke-width="3"/>')
        parts.append(t((x + cx) / 2 + (-20 if lab == "no" else 20), (YL[d] + 56 + YL[cd]) / 2 + 8, lab, 24, SOFT, "end" if lab == "no" else "start", 600))
for i, nd in T.items():
    x, d = pos[i]
    if nd["left"] != -1:
        w = 300 if d < 2 else 280
        parts.append(f'<rect x="{x - w / 2}" y="{YL[d]}" width="{w}" height="56" rx="12" fill="{INK}"/>')
        parts.append(t(x, YL[d] + 38, f"¿{nd['tf']} ≥ {nd['cut']}?", 28, PAPER, "middle", 600))
    else:
        p = nd["pops"]
        top = max(p, key=p.get)
        fill = POPC[top] if p[top] / nd["n"] >= 0.9 else LINE
        tc = PAPER if fill != LINE else INK
        parts.append(f'<rect x="{x - 92}" y="{YL[d]}" width="184" height="104" rx="16" fill="{fill}"/>')
        parts.append(t(x, YL[d] + 30, f"{nd['n']} células", 24, tc, "middle"))
        parts.append(t(x, YL[d] + 66, f"{nd['mean']:.1f}", 34, tc, "middle", 700))
        parts.append(t(x, YL[d] + 94, comp(p), 24, tc, "middle"))
fig = svg("Árbol de regresión real de tres niveles para NKG7. Primera pregunta: TBX21 ≥ 1. Del lado no: KLRG1 ≥ 1 y luego RORA ≥ 2 o MEF2C ≥ 1. Del lado sí: UBB ≥ 8 y luego RPL6 ≥ 28 o SNRPB2 ≥ 3. Ocho grupos finales, cada uno con su promedio de NKG7, de 0.7 (592 monocitos y B) a 58.6 (8 NK).", 1616, 1600, 444, "\n".join(parts))
out["arbol3"] = slide(
    "arbol3", "Paso ① por dentro · 3 de 6", "Un árbol real: 3 niveles de preguntas, 8 grupos",
    figure(fig, 1664, 0) + "\n" + callout(f"Error: {mse0:.0f} → {mse1:.0f} (1 nivel) → <b>{mse3:.0f} (3 niveles)</b>. GRNBoost2 hace cientos así por gen.", 1),
    f"""Minutos 27–29. Dentro de cada lado se vuelve a hacer exactamente lo mismo: probar todos los cortes de todos los TF y quedarse con el mejor, pero solo con las células de ese lado. Cada nivel duplica los grupos: 3 niveles = hasta 8 grupos. La predicción de una célula es el promedio de NKG7 del grupo donde termina.
Árbol real (todas las células y todos los TF, profundidad 3):
- Lado no de TBX21 (675): ¿KLRG1 ≥ 1? Si no, ¿RORA ≥ 2? → 592 células (297 Mono, 283 B, 12 NK) con 0.7, o 14 (13 NK) con 30.0. Si sí, ¿MEF2C ≥ 1? → 56 (54 NK) con 29.8, o 13 (9 B, 3 Mono, 1 NK) con 2.5. Aquí rescata a las NK con dropout de TBX21.
- Lado sí (225): ¿UBB ≥ 8? y luego RPL6 o SNRPB2: grupos de NK con más o menos NKG7 (31.3, 38.3, 58.6) y uno de 4 células (3 B) con 8.5.
Ojo: UBB, RPL6 y SNRPB2 están en la lista de 1,018 "TF" que usa pySCENIC, aunque no son TF clásicos. Las preguntas de abajo afinan grupos pequeños; un árbol grande aprende detalles de pocas células. Por eso GRNBoost2 no usa un árbol grande, sino muchos árboles de 3 niveles que se corrigen (siguiente diapositiva).
Error: {mse0:.0f} → {mse1:.0f} → {mse3:.0f}.""")

# ---------------------------------------------------------------- 2. Boosting
TR = A["trees_random60"]
TRN25 = sum(1 for g, v in list(A["trees_per_gene"].items())[:60] if v == 25)
nk, mo = GX["nk_cell"]["traj"], GX["mono_cell"]["traj"]
pts = sorted(int(k) for k in nk)
X0, X1, Y0, PX = 120, 820, 420, 8.5
xs = lambda k: X0 + k * (X1 - X0) / GX["n_trees"]
ys = lambda v: Y0 - v * PX
ax = [f'<line x1="{X0}" y1="{Y0}" x2="{X1 + 10}" y2="{Y0}" stroke="{INK}" stroke-width="3"/>',
      f'<line x1="{X0}" y1="{Y0}" x2="{X0}" y2="60" stroke="{INK}" stroke-width="3"/>']
for v in (0, 10, 20, 30, 40):
    ax.append(t(X0 - 14, ys(v) + 9, v, 26, SOFT, "end"))
    if v:
        ax.append(f'<line x1="{X0}" y1="{ys(v)}" x2="{X1}" y2="{ys(v)}" stroke="{LINE}" stroke-width="1.5"/>')
for k in (0, 100, 200, 300, 400, GX["n_trees"]):
    ax.append(t(xs(k), Y0 + 36, k, 26, SOFT, "middle"))
ax.append(t((X0 + X1) / 2, 494, "número de árboles →", 28, SOFT, "middle"))
ax.append(t(40, 240, "NKG7 predicho (UMI)", 28, SOFT, "middle", extra=' transform="rotate(-90 40 240)"'))
real_nk = GX["nk_cell"]["real"]
ax.append(f'<line x1="{X0}" y1="{ys(real_nk)}" x2="{X1}" y2="{ys(real_nk)}" stroke="{NK}" stroke-width="3" stroke-dasharray="10 8"/>')
ax.append(t(X0 + 12, ys(real_nk) - 12, f"valor real de la célula NK: {real_nk}", 28, NK, "start", 600))
for traj, col in ((nk, NK), (mo, MONO)):
    d = " ".join(f"{'M' if i == 0 else 'L'}{xs(k):.1f},{ys(traj[str(k)]):.1f}" for i, k in enumerate(pts))
    ax.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="5"/>')
    ax += [f'<circle cx="{xs(k):.1f}" cy="{ys(traj[str(k)]):.1f}" r="6" fill="{col}"/>' for k in pts]
ax.append(t(xs(GX["n_trees"]), ys(nk[str(GX['n_trees'])]) + 40, f"célula NK: {nk[str(GX['n_trees'])]:.1f}", 28, NK, "end", 600))
ax.append(t(xs(250), ys(mo["200"]) - 22, f"monocito (real: 0): {mo[str(GX['n_trees'])]:.1f}", 28, MONO, "start", 600))
fig = svg(f"Predicción real de NKG7 para dos células conforme se suman árboles: una célula NK con {real_nk} UMI sube de {nk['0']:.1f} a {nk[str(GX['n_trees'])]:.1f}; un monocito con 0 baja de {mo['0']:.1f} a {mo[str(GX['n_trees'])]:.1f}.", 800, 860, 460, "\n".join(ax), vy=50)
out["boosting"] = slide(
    "boosting", "Paso ① por dentro · 4 de 6", "Boosting: muchos árboles pequeños que se corrigen",
    '<div class="row" style="gap:56px; align-items:center;">\n' + figure(fig, 848, 0) + "\n" + steps([
        f"Se empieza con el promedio: {nk['0']:.1f} para todas las células.",
        "Desde el primer árbol, cada árbol predice el <b>residuo</b> (real − predicción) y la predicción avanza solo 1 % hacia él.",
        f"Cada árbol tiene 3 niveles y ve el 90 % de las células; en cada pregunta solo escoge entre {GX['max_features_n']} TF al azar.",
        f"Se detiene cuando deja de mejorar. NKG7: {GX['n_trees']} árboles; en 60 genes al azar, de {TR['min']} a {TR['max']}.",
    ]) + "\n</div>\n" + callout(f"Predicción final = promedio + <b>{GX['n_trees']} pequeñas correcciones del residuo</b>."),
    f"""Minutos 29–31. Datos reales: el mismo modelo que corrió grn para NKG7 (semilla 2026), seguido célula por célula.
Boosting = sumar muchos árboles débiles. Se parte del promedio. Desde ahí, cada árbol se entrena para predecir el residuo de cada célula, lo que todavía falta (real − predicción actual), y su respuesta se suma multiplicada por 0.01 (learning rate). Con error cuadrático, el residuo es justo la dirección en que más baja el error (el gradiente, de ahí "gradient boosting").
Ejemplo: la célula NK tiene {real_nk} UMI y parte de {nk['0']:.1f}, residuo {real_nk - nk['0']:.1f}. Cada árbol cierra solo un poquito de ese hueco, el residuo se achica y la línea sube despacio.
Célula NK con {real_nk} UMI: {nk['0']:.1f} → {nk['10']:.1f} (10 árboles) → {nk['100']:.1f} (100) → {nk[str(GX['n_trees'])]:.1f} ({GX['n_trees']}). Monocito con 0: {mo['0']:.1f} → {mo['100']:.1f} → {mo[str(GX['n_trees'])]:.1f}.
Error cuadrático medio de todas las células: {G['NKG7']['mse_0']:.0f} al inicio, {G['NKG7']['mse_by_tree']['100']:.0f} con 100 árboles, {G['NKG7']['mse_by_tree'][str(GX['n_trees'])]:.0f} al final.
Azar a propósito: cada árbol se entrena con el 90 % de las células, y en cada pregunta solo puede escoger entre {GX['max_features_n']} de los {fmt(G['n_tf'] - 1)} TF, sorteados de nuevo en cada pregunta. Así distintos TF tienen oportunidad; los primeros cinco árboles empezaron preguntando por {', '.join(G['NKG7']['roots_first5'])}.
Parada temprana: con el 10 % de células que cada árbol no vio se mide si sigue mejorando; si en 25 árboles seguidos no mejora, se detiene. Por eso cada gen termina con un número distinto de árboles (máximo 5,000). CD79A: {G['CD79A']['n_trees']}. En 60 genes al azar de la matriz: mediana {TR['median']:.0f}, de {TR['min']} a {TR['max']}; {TRN25} de 60 se quedan en 25, el mínimo: el primer árbol ya no mejora nada, o sea que ningún TF predice ese gen. Marcadores: LYZ {A['trees_per_gene']['LYZ']}, MS4A1 {A['trees_per_gene']['MS4A1']}, GZMB {A['trees_per_gene']['GZMB']}.""")

# ---------------------------------------------------------------- 3. Importancia
IT = json.loads((AQUI / "importancia_tabla.json").read_text())
TOP = {r["tf"]: r for r in IT["top8"]}
show = IT["tbx21_rows"][:4]
rows, y = [], 44
rows.append(t(20, y, "Árbol", 28, SOFT, "start", 600))
rows.append(t(700, y, "Baja del error con TBX21", 28, SOFT, "end", 600))
rows.append(f'<line x1="20" y1="{y + 16}" x2="700" y2="{y + 16}" stroke="{INK}" stroke-width="2"/>')
y += 60
for i, (tree, d) in enumerate(show):
    if i % 2 == 0:
        rows.append(f'<rect x="10" y="{y - 31}" width="700" height="42" rx="6" fill="{LINE}" fill-opacity="0.35"/>')
    rows.append(t(20, y, fmt(tree), 28, INK))
    rows.append(t(700, y, f"{d:.1f}", 28, INK, "end"))
    y += 42
rows.append(t(20, y, f"… {IT['tbx21_n'] - len(show) - 1} árboles más …", 28, SOFT))
y += 42
last_tree, last_d = IT["tbx21_rows"][-1]
rows.append(t(20, y, fmt(last_tree), 28, INK))
rows.append(t(700, y, f"{last_d:.1f}", 28, INK, "end"))
rows.append(f'<line x1="20" y1="{y + 18}" x2="700" y2="{y + 18}" stroke="{INK}" stroke-width="2"/>')
y += 58
rows.append(f'<rect x="10" y="{y - 36}" width="700" height="52" rx="8" fill="{NK}"/>')
rows.append(t(20 + 10, y, f"Suma ({IT['tbx21_n']} árboles)", 30, PAPER, "start", 600))
rows.append(t(690, y, fmt(round(IT["tbx21_sum"])), 30, PAPER, "end", 600))
y += 60
rows.append(t(20, y, f"Importancia = {fmt(round(IT['tbx21_sum']))} ÷ {IT['mean_drop_per_tree']:.1f}", 30, INK, "start", 600))
rows.append(t(700, y, f"{IT['tbx21_imp']:.1f}", 30, NK, "end", 600))
H = y + 20
fig = svg(f"Tabla real: en el modelo de NKG7, TBX21 se usa en {IT['tbx21_n']} de {IT['n_trees']} árboles. Bajas del error en los primeros: " + ", ".join(f"árbol {a}: {b:.1f}" for a, b in show) + f"; el último, árbol {last_tree}: {last_d:.1f}. Suma: {fmt(round(IT['tbx21_sum']))}. Dividida entre lo que baja un árbol promedio, {IT['mean_drop_per_tree']:.1f}, da la importancia {IT['tbx21_imp']:.1f}.", 720, 720, H, "\n".join(rows))
R, M, K = TOP["RORA"], TOP["MYBL1"], TOP["SPI1"]
out["importancia"] = slide(
    "importancia", "Paso ① por dentro · 5 de 6", "Importancia = suma de lo que baja el error",
    '<div class="row" style="gap:56px; align-items:center;">\n' + figure(fig, 768, 0) + "\n" + steps([
        "Cada pregunta de un árbol baja el error. Esa baja se le apunta al TF de la pregunta.",
        f"Se suman las bajas de cada TF en todos los árboles. TBX21: <b>{fmt(round(IT['tbx21_sum']))}</b>.",
        f"Para tener números chicos se divide entre lo que baja un árbol promedio: TBX21 = <b>{IT['tbx21_imp']:.1f}</b>.",
        f"Igual con cada TF: RORA {fmt(round(R['sum_drop']))} → {R['imp']:.1f}; MYBL1 {fmt(round(M['sum_drop']))} → {M['imp']:.1f}.",
    ]) + "\n</div>\n" + callout("Importancia alta = <b>ayuda mucho a predecir</b>. Aún no sabemos si lo activa, lo reprime o solo coincide."),
    f"""Minutos 31–33. Así nace el número de la tabla de adyacencias. La importancia es la suma, en todos los árboles, de cuánto baja el error cada vez que se pregunta por ese TF.
La "baja del error" de una pregunta: error cuadrático del grupo antes de partirlo menos el de los dos lados después, ponderado por la fracción de células que llegan a esa pregunta (es la importancia por impureza de sklearn; cada árbol mide su error sobre el residuo que está aprendiendo). Unidades: UMI².
Tabla real (modelo de NKG7, semilla 2026, instructor/slides/datos/importancia_tabla.log): TBX21 aparece en {IT['tbx21_n']} de {IT['n_trees']} árboles. Las bajas grandes vienen al principio (árbol 4: {dict(IT['tbx21_rows'])[4]:.1f}), cuando el residuo todavía es grande; al final son chicas (árbol {last_tree}: {last_d:.1f}). Suma: {IT['tbx21_sum']:.2f}.
El último paso es solo escala: arboreto (el paquete detrás de GRNBoost2) divide entre la baja total de todos los TF ({fmt(round(IT['total_drop_all_tf']))}) y multiplica por el número de árboles ({IT['n_trees']}). Es lo mismo que dividir entre lo que baja un árbol promedio ({IT['mean_drop_per_tree']:.2f}): {IT['tbx21_sum']:.1f} ÷ {IT['mean_drop_per_tree']:.2f} = {IT['tbx21_imp']:.2f}. Por eso las importancias de un gen suman {GX['sum_importances']:.0f}, su número de árboles. El orden de los TF no cambia.
Los 8 primeros de NKG7 (suma → importancia): """ + "; ".join(f"{r['tf']} {r['sum_drop']:.0f} → {r['imp']:.1f}" for r in IT["top8"]) + f""".
Comprobación: reentrenamos el modelo de NKG7 con la semilla del servidor (2026) y salieron exactamente los mismos valores de adjacencies.tsv (TBX21 66.09, RORA 49.18, MYBL1 37.38…); la suma ÷ {IT['mean_drop_per_tree']:.2f} reproduce las 8 importancias (diferencia máxima {IT['max_abs_diff_check']:.0e}). La tabla no es una caja negra.
Ojo con SPI1 en el lugar 7 ({K['sum_drop']:.0f} → {K['imp']:.1f}): SPI1 casi no está en las NK; su AUSENCIA ayuda a predecir NKG7. Es el puente a la siguiente slide.
Se repite para los 10,120 genes: 2,339,328 filas TF – gen – importancia. Consecuencia: los genes con muchos árboles (fáciles de predecir y con mucha variación) tienen importancias más grandes. Por eso pySCENIC usa percentiles y rankings en vez de un corte fijo.""")

# ---------------------------------------------------------------- 4. Tabla (grnboost)
out["grnboost"] = slide(
    "grnboost", "Paso ① por dentro · 6 de 6 · la tabla de adyacencias", "El resultado: ¿quién predice a quién?",
    """<div class="row" style="gap:56px; align-items:flex-start;">
<div style="flex:1;">
<table class="data">
<tr><th>Gen</th><th>TF que más lo predice</th><th class="n">Import.</th><th class="n">ρ</th></tr>
<tr><td>NKG7</td><td class="tf">TBX21</td><td class="n">66.1</td><td class="n">+0.73</td></tr>
<tr><td>NKG7</td><td class="tf">RORA</td><td class="n">49.2</td><td class="n">+0.67</td></tr>
<tr><td>NKG7</td><td class="tf">MYBL1</td><td class="n">37.4</td><td class="n">+0.60</td></tr>
<tr class="hl"><td>CD79A</td><td>ANXA1</td><td class="n">45.1</td><td class="n neg">−0.73</td></tr>
<tr class="hl"><td>CD79A</td><td>CEBPD</td><td class="n">33.2</td><td class="n neg">−0.65</td></tr>
<tr><td>CD79A</td><td class="tf">SPIB</td><td class="n">27.6</td><td class="n">+0.65</td></tr>
</table>
</div>
<div class="card" style="--i:1; flex:0 0 560px;">
<p class="tag">¿Por qué ANXA1?</p>
<p class="txt">ANXA1 está en el 98 % de los monocitos y solo en el 13 % de las células B.</p>
<p class="txt">Si falta ANXA1, probablemente es una B y hay CD79A. Predice por <b>ausencia</b>.</p>
<p class="txt">ρ = correlación, calculada aparte: GRNBoost2 no la da.</p>
</div>
</div>
""" + callout("La importancia <b>no tiene signo</b>: un TF puede predecir a un gen porque falta donde el gen está."),
    """Minutos 33–35. La tabla es real (servidor, 900 células). ρ = correlación de Spearman sobre log(CP10k), calculada aparte; GRNBoost2 no la reporta.
- NKG7 (gen de citotoxicidad NK): sus mejores predictores son TBX21, RORA y MYBL1, todos con correlación positiva. Tiene sentido.
- CD79A (gen de células B): el mejor predictor es ANXA1, con correlación NEGATIVA. ANXA1 está en 98 % de los monocitos y solo en 13 % de las B: su AUSENCIA predice CD79A. Lo mismo CEBPD (un TF de monocitos).
- En el árbol real de una sola pregunta para CD79A, la pregunta ganadora es "¿ANXA1 ≥ 2 UMI?": del lado "no" quedan 294 de las 300 B.
Moraleja: importancia ≠ activación ≠ causalidad. cisTarget separa dianas con correlación positiva (regulones "+") y exige el motivo. ANXA1 no tiene ningún motivo anotado, así que nunca tendrá regulón.
Luego regresan al notebook: comparan su tabla con la del servidor (pregunta del top 15 de TBX21).""")

# ---------------------------------------------------------------- 5. Motivo
seq = "TTCCAAAGAGGAAGTGACTCGT"
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
m0, m1 = 6, 16  # AGAGGAAGTG
sx, dx = 52, 37
parts = [f'<ellipse cx="{sx + dx * (m0 + m1 - 1) / 2:.0f}" cy="70" rx="190" ry="48" fill="{MONO}"/>',
         t(sx + dx * (m0 + m1 - 1) / 2, 82, "PU.1 (proteína de SPI1)", 30, PAPER, "middle", 600),
         f'<rect x="{sx + dx * m0 - 20}" y="128" width="{dx * (m1 - m0)}" height="112" rx="10" fill="{ACC}" opacity="0.45"/>',
         f'<line x1="20" y1="184" x2="880" y2="184" stroke="{LINE}" stroke-width="2"/>']
for i, ch in enumerate(seq):
    x = sx + dx * i
    w = 700 if m0 <= i < m1 else None
    parts.append(t(x, 168, ch, 34, INK, "middle", w, ' font-family="IBM Plex Mono, Courier New, monospace"'))
    parts.append(t(x, 224, comp[ch], 34, SOFT, "middle", None, ' font-family="IBM Plex Mono, Courier New, monospace"'))
parts.append(t(20, 290, "El mismo motivo en dos colecciones:", 28, SOFT, "start", 600))


def consenso(y, etiqueta, s, core):
    out_ = [t(20, y, etiqueta, 28, INK)]
    x = 470
    i = s.index(core)
    for j, ch in enumerate(s):
        col = MONO if i <= j < i + len(core) else INK
        out_.append(t(x + j * 30, y, ch, 32, col, "middle", 700 if col == MONO else None, ' font-family="IBM Plex Mono, Courier New, monospace"'))
    return out_


parts += consenso(345, "HOMER (ChIP-seq)", "AGAGGAAGTG", "GGAA")
parts += consenso(400, "Taipale 2013 (in vitro)", "AAAAAGCGGAAGTN", "GGAA")
parts.append(t(20, 452, "esquema · la secuencia de arriba es ilustrativa; los motivos son reales", 26, SOFT))
fig = svg("Esquema: la proteína PU.1 se pega a una secuencia corta de ADN, AGAGGAAGTG. Abajo, el mismo motivo de PU.1 en dos colecciones reales, AGAGGAAGTG y AAAAAGCGGAAGTN, con el núcleo GGAA resaltado.", 860, 900, 470, "\n".join(parts))
out["motivo"] = slide(
    "motivo", "Paso ② por dentro · 1 de 6 · antes de correr ctx", "¿Qué es un motivo?",
    '<div class="row" style="gap:56px; align-items:center;">\n' + figure(fig, 908, 0) + "\n" + steps([
        "Un TF no se pega en cualquier parte del ADN: reconoce una <b>palabra corta</b>, de unas 6 a 20 letras. Esa palabra es su motivo.",
        "No es exacta: algunas posiciones aceptan varias letras (N = cualquiera).",
        "PU.1 reconoce …<b>GGAA</b>…; ese núcleo lo comparte toda su familia de TF (ETS).",
        "Un TF puede tener muchos motivos: a SPI1 se le asignan 390 en nuestra tabla.",
    ]) + "\n</div>\n" + callout("Si el motivo de un TF está cerca de un gen, ese TF <b>puede</b> pegarse ahí y regularlo."),
    """Minutos 44–46. Antes, todos lanzaron la celda de recursos (baja ~400 MB desde el Drive en ~20 s). Todavía NO lanzan ctx.
Arranca desde cero: el ADN son letras (A, C, G, T). Un factor de transcripción es una proteína con una parte que encaja en una secuencia corta. Esa secuencia preferida es su motivo.
No es una palabra exacta: el TF tolera cambios en algunas posiciones. Por eso los motivos se guardan como tablas de probabilidad por posición; para hablar los escribimos como un consenso (N = cualquier letra, R = A o G).
Los dos consensos son reales y vienen de nuestra tabla de anotación: AGAGGAAGTG_PU.1 (HOMER, ChIP-seq en macrófagos) y SPI1_full_AAAAAGCGGAAGTN (Taipale et al. 2013, in vitro). Distintas técnicas dan versiones un poco distintas; el núcleo GGAA es el sello de la familia ETS.
La secuencia de 22 letras de arriba es inventada para ilustrar.
"390 motivos": son motivos (o grupos de motivos parecidos) que la tabla asigna a SPI1, ya sea directamente o por parecido con otro motivo de SPI1. Varios TF de la misma familia comparten motivos: cisTarget no siempre puede distinguir PU.1 de SPIB.
Remate: tener el motivo cerca no garantiza que el TF se pegue (depende de la cromatina abierta, de otros TF...). Es evidencia, no prueba.""")

# ---------------------------------------------------------------- 6. Bases
k = SP["known"]
ke = EB["known"]
parts = [f'<rect x="120" y="14" width="700" height="18" rx="9" fill="{LILA}"/>',
         t(470, 66, "ventana amplia: 10 kb a cada lado del inicio", 28, NK, "middle", 600),
         f'<path d="M470,128 L470,104 L520,104" fill="none" stroke="{INK}" stroke-width="3"/>',
         f'<path d="M514,97 L526,104 L514,111 z" fill="{INK}"/>',
         t(535, 113, "inicio", 26, INK),
         f'<line x1="20" y1="150" x2="880" y2="150" stroke="{INK}" stroke-width="4"/>',
         f'<rect x="470" y="128" width="360" height="44" rx="8" fill="{TEAL}"/>',
         t(650, 159, "gen", 28, INK, "middle", 600),
         f'<rect x="420" y="182" width="66" height="18" rx="9" fill="{MONO}"/>',
         t(453, 234, "promotor: 500 pb antes y 100 pb después del inicio", 28, MONO, "middle", 600),
         f'<line x1="20" y1="262" x2="880" y2="262" stroke="{LINE}" stroke-width="2" stroke-dasharray="6 6"/>',
         t(20, 300, "La base: lugar de cada gen en el ranking de cada motivo", 28, SOFT, "start", 600)]
cols = [("CSF1R", 470), ("LYZ", 620), ("CD79A", 770)]
parts += [t(x, 344, g, 28, INK, "middle", 600) for g, x in cols]
parts.append(f'<line x1="20" y1="358" x2="880" y2="358" stroke="{LINE}" stroke-width="2"/>')
parts.append(t(20, 396, "motivo de PU.1", 28, INK))
parts += [t(x, 396, fmt(k[g] + 1), 28, INK, "middle", None, ' font-family="IBM Plex Mono, Courier New, monospace"') for g, x in cols]
parts.append(t(20, 438, "un motivo de EBF1", 28, INK))
parts += [t(x, 438, fmt(ke[g] + 1), 28, INK, "middle", None, ' font-family="IBM Plex Mono, Courier New, monospace"') for g, x in cols]
parts.append(t(20, 476, "… 5,876 motivos × 27,015 genes · dibujo de arriba no a escala", 26, SOFT))
fig = svg(f"Arriba, un gen con dos ventanas de búsqueda: el promotor, de 500 pares de bases antes a 100 después del inicio, y una ventana amplia de 10 kb a cada lado. Abajo, posiciones reales en la base: con el motivo de PU.1, CSF1R está en el lugar {k['CSF1R'] + 1}, LYZ en el {fmt(k['LYZ'] + 1)} y CD79A en el {fmt(k['CD79A'] + 1)}.", 840, 900, 486, "\n".join(parts))
out["bases"] = slide(
    "bases", "Paso ② por dentro · 2 de 6", "Las bases: el genoma ya escaneado",
    '<div class="row" style="gap:56px; align-items:center;">\n' + figure(fig, 888, 0) + "\n" + steps([
        "Alguien ya buscó cada motivo cerca de cada gen humano. Nosotros solo bajamos el resultado (313 MB).",
        "¿Qué tan cerca? Dos bases: el <b>promotor</b> o una <b>ventana amplia</b> de 10 kb a cada lado.",
        "Por motivo, los 27,015 genes quedan ordenados del que mejor tiene el motivo al que peor: un <b>ranking</b>.",
        "Otra tabla dice de qué TF es cada motivo (1,605 TF en total).",
    ]) + "\n</div>\n" + callout("En clase usamos solo la base del <b>promotor</b>; la de 10 kb queda para casa."),
    f"""Minutos 46–48. Esta es la base que acaban de bajar. Al terminar la slide, lanzan ctx (48–50).
Quién la hizo: el laboratorio de Stein Aerts. Para cada motivo calificaron qué tan bien (y cuántas veces) aparece cerca de cada gen humano, con el programa Cluster-Buster, y convirtieron la calificación en un ranking. Hacerlo en clase tomaría días.
Dos ventanas: el promotor (archivo hg38_500bp_up_100bp_down…) y una ventana amplia de 10 kb a cada lado (hg38_10kbp_up_10kbp_down…), que atrapa también regiones reguladoras cercanas. Cada base pesa ~310 MB.
La tabla de abajo son posiciones reales del archivo del promotor (1 = el gen que mejor tiene el motivo):
- Motivo de PU.1: CSF1R en el lugar {k['CSF1R'] + 1} (CSF1R es una diana clásica de PU.1: Zhang et al. 1994, Mol Cell Biol), LYZ en el {fmt(k['LYZ'] + 1)}, CD79A (gen de células B) en el {fmt(k['CD79A'] + 1)}.
- Un motivo que la tabla asigna a EBF1 (tfdimers__MD00215): CSF1R {ke['CSF1R'] + 1}, LYZ {fmt(ke['LYZ'] + 1)}, CD79A {fmt(ke['CD79A'] + 1)}.
Mismo gen, distinta posición según el motivo: cada motivo tiene su propio ranking.
Tabla de anotación (archivo .tbl): traduce motivo → TF, directamente (experimento con ese TF), por parecido con otro motivo, o por ortología con ratón.""")

# ---------------------------------------------------------------- 7. Ranking
r = SP["module_gene_ranks_sorted"]
cut = C["rank_cutoff"] + 1
bins = [0] * 20
for x in r:
    bins[min(x // cut, 19)] += 1
exp = SP["in_db"] * cut / C["total_genes"]
X0, BW, Y0 = 90, 33, 400
parts = [t(20, 36, f"Dónde caen las {fmt(SP['in_db'])} candidatas", 32, INK, "start", 600)]
for i, c in enumerate(bins):
    h = c * 1.0
    col = MONO if i == 0 else LINE
    parts.append(f'<rect x="{X0 + i * BW + 3}" y="{Y0 - h:.0f}" width="{BW - 6}" height="{h:.0f}" fill="{col}"/>')
parts.append(f'<line x1="{X0}" y1="{Y0}" x2="{X0 + 20 * BW}" y2="{Y0}" stroke="{INK}" stroke-width="3"/>')
parts.append(f'<line x1="{X0}" y1="{Y0 - exp:.0f}" x2="{X0 + 20 * BW}" y2="{Y0 - exp:.0f}" stroke="{INK}" stroke-width="3" stroke-dasharray="10 8"/>')
parts.append(t(X0 + 20 * BW, Y0 - exp - 12, f"al azar: {exp:.0f}", 30, INK, "end", 600))
parts.append(t(X0 + BW / 2, Y0 - bins[0] - 12, bins[0], 32, MONO, "middle", 700))
parts.append(t(X0, Y0 + 38, "1", 28, SOFT, "start"))
parts.append(t(X0 + 20 * BW, Y0 + 38, "27,015", 28, SOFT, "end"))
parts.append(t(X0 + 10 * BW, Y0 + 80, f"lugar en el ranking de PU.1 (bloques de {fmt(cut)})", 28, SOFT, "middle"))
fig = svg(f"Histograma real: de las {fmt(SP['in_db'])} candidatas de SPI1, {bins[0]} caen en el primer bloque de {fmt(cut)} lugares del ranking del motivo de PU.1, contra {exp:.0f} esperadas al azar; los bloques siguientes van bajando.", 660, 760, 492, "\n".join(parts), vy=8)
tabla = f"""<div style="display:flex; flex-direction:column; gap:20px; width:640px;">
<p style="font-size:28px; color:var(--soft);">Candidatas = las {fmt(SP['module_genes'])} posibles dianas de SPI1 que dio GRNBoost2.</p>
<table class="data">
<tr><th class="n">Lugar</th><th>Gen</th><th>¿Candidata?</th></tr>
<tr><td class="n">1</td><td>PCSK6-AS1</td><td>no</td></tr>
<tr class="hl"><td class="n">7</td><td>DENND4B</td><td>sí</td></tr>
<tr class="hl"><td class="n">{k['CSF1R'] + 1}</td><td>CSF1R</td><td>sí</td></tr>
<tr class="hl"><td class="n">{k['LYZ'] + 1}</td><td>LYZ</td><td>sí</td></tr>
<tr class="hl"><td class="n">{fmt(k['S100A8'] + 1)}</td><td>S100A8</td><td>sí</td></tr>
<tr><td class="n">{fmt(k['CD79A'] + 1)}</td><td>CD79A</td><td>no</td></tr>
</table>
</div>"""
out["ranking"] = slide(
    "ranking", "Paso ② por dentro · 3 de 6 · mientras corre ctx", "El ranking: ¿las candidatas de SPI1 tienen el motivo?",
    '<div class="row" style="gap:64px; align-items:center;">\n' + tabla + "\n" + figure(fig, 708, 0) + "\n</div>\n"
    + callout(f"En el primer 5 % del ranking ({fmt(cut)} genes) caen <b>{bins[0]}</b> candidatas. Al azar esperaríamos {exp:.0f}."),
    f"""Minutos 50–52. Ya lanzaron ctx (~5 min). Todo es real: base del promotor, motivo de PU.1 (grupo de motivos metacluster_167.5, anotado directamente a SPI1).
Primero, las candidatas: de la tabla de GRNBoost2, pySCENIC arma para cada TF seis grupos ("módulos"). Este es uno: los pares SPI1 → gen con importancia en el 10 % más alto de toda la tabla y correlación positiva. Son {fmt(SP['module_genes'])} genes; {fmt(SP['in_db'])} están en la base. Los otros cinco grupos: importancia en el 25 % más alto, las 50 mejores dianas del TF, y los genes donde el TF está entre sus 5, 10 o 50 mejores predictores.
La tabla de la izquierda es el ranking del motivo de PU.1 (los 27,015 genes). Los primeros lugares no son candidatas (el motivo está en muchos genes que SPI1 no regula aquí). Pero CSF1R (lugar {k['CSF1R'] + 1}), LYZ ({k['LYZ'] + 1}) y S100A8 ({fmt(k['S100A8'] + 1)}), genes de monocito, son candidatas y están arriba. CD79A, gen de células B, está casi al final ({fmt(k['CD79A'] + 1)}) y no es candidata (su correlación con SPI1 es negativa).
Histograma: las candidatas se amontonan al principio del ranking: {bins[0]} en el primer bloque contra {exp:.0f} esperadas si estuvieran repartidas al azar. Más del doble.
Pregunta al grupo: si SPI1 de verdad regula a estas candidatas, ¿qué esperarían ver? (que tengan su motivo → que queden arriba).""")

# ---------------------------------------------------------------- 8 y 9. Curva y poda
cur, avg, avg2 = SP["curve"], SP["avg"], SP["avg2"]
X0, X1, Y0, YS = 110, 820, 420, 0.45
xr = lambda rk: X0 + rk * (X1 - X0) / 5000
yr = lambda v: Y0 - v * YS


def ejes():
    a = [f'<line x1="{X0}" y1="{Y0}" x2="{X1 + 10}" y2="{Y0}" stroke="{INK}" stroke-width="3"/>',
         f'<line x1="{X0}" y1="{Y0}" x2="{X0}" y2="50" stroke="{INK}" stroke-width="3"/>']
    for v in (0, 200, 400, 600, 800):
        a.append(t(X0 - 14, yr(v) + 9, v, 26, SOFT, "end"))
    for rk in (0, 1000, 2000, 3000, 4000, 5000):
        a.append(t(xr(rk), Y0 + 36, fmt(rk), 26, SOFT, "middle"))
    a.append(t((X0 + X1) / 2, 494, "posición en el ranking del motivo →", 28, SOFT, "middle"))
    a.append(t(36, 240, "candidatas halladas", 28, SOFT, "middle", extra=' transform="rotate(-90 36 240)"'))
    return a


def linea(vals, **kw):
    d = " ".join(f"{'M' if i == 0 else 'L'}{xr(i * 25):.1f},{yr(v):.1f}" for i, v in enumerate(vals))
    attrs = " ".join(f'{k_.replace("_", "-")}="{v_}"' for k_, v_ in kw.items())
    return f'<path d="{d}" fill="none" {attrs}/>'


ic = C["rank_cutoff"] // 25 + 1
area = f"M{xr(0)},{Y0} " + " ".join(f"L{xr(i * 25):.1f},{yr(cur[i]):.1f}" for i in range(ic)) + f" L{xr((ic - 1) * 25):.1f},{Y0} Z"
parts = ejes() + [
    f'<path d="{area}" fill="{B}" opacity="0.22"/>',
    linea(avg, stroke=SOFT, stroke_width=4, stroke_dasharray="10 8"),
    linea(cur, stroke=B, stroke_width=5),
    f'<line x1="{xr(C["rank_cutoff"]):.1f}" y1="{Y0}" x2="{xr(C["rank_cutoff"]):.1f}" y2="70" stroke="{MONO}" stroke-width="3" stroke-dasharray="6 6"/>',
    t(xr(C["rank_cutoff"]) + 12, 92, f"5 % = {fmt(cut)}", 28, MONO, "start", 600),
    t(xr(1500), 150, "motivo de PU.1 (verde)", 30, B, "start", 700),
    t(xr(1750), Y0 - 26, "promedio de los 5,876 motivos", 28, SOFT, "start", 600),
    t(xr(1000), Y0 - 14, "AUC", 30, B, "middle", 700),
]
fig = svg(f"Curva de recuperación real de las candidatas de SPI1 en el ranking del motivo de PU.1, comparada con el promedio de los 5,876 motivos. En la posición {fmt(cut)} la curva de PU.1 lleva {cur[ic - 1]} candidatas y el promedio {avg[ic - 1]:.0f}. El área sombreada hasta esa posición es el AUC.", 800, 860, 466, "\n".join(parts), vy=44)
out["curva"] = slide(
    "curva", "Paso ② por dentro · 4 de 6", "La curva de recuperación y su AUC",
    '<div class="row" style="gap:56px; align-items:center;">\n' + figure(fig, 848, 0) + "\n" + steps([
        "Se baja por el ranking del motivo, gen por gen, y se cuentan las candidatas que van apareciendo.",
        "Si las candidatas tienen el motivo, aparecen pronto: la curva sube rápido.",
        f"<b>AUC</b> = el área bajo la curva en el primer 5 % ({fmt(cut)} genes). Aquí: {SP['auc']:.3f}.",
        "Se repite con los 5,876 motivos de la base: cada uno da su propia curva y su AUC.",
    ]) + "\n</div>\n" + callout(f"En la posición {fmt(cut)}, PU.1 lleva <b>{cur[ic - 1]}</b> candidatas; el motivo promedio, {avg[ic - 1]:.0f}."),
    f"""Minutos 52–54. Curva real, mismas candidatas y mismo motivo que en la slide anterior.
Recuperación = cuántas candidatas "recupero" conforme bajo por el ranking. El eje horizontal llega a 5,000 posiciones (el máximo que usa cisTarget); el área solo se mide en el primer 5 % ({fmt(cut)} posiciones; opción auc_threshold 0.05).
AUC: el área sombreada dividida entre el área máxima posible (todas las candidatas en los primeros lugares). Por eso los AUC son números pequeños: {SP['auc']:.4f} para PU.1, {SP['auc_mean']:.4f} el promedio de todos los motivos.
La línea punteada es el promedio de los 5,876 motivos con las mismas candidatas: lo que pasa con un motivo cualquiera. Que el promedio también suba es normal: algunos genes salen arriba en muchos rankings.
Es la misma lógica que AUCell más adelante: allá el ranking será de los genes de una célula.""")

hist, edges = SP["auc_hist"], SP["auc_edges"]
X0, X1, Y0, YS = 100, 830, 400, 0.82
xa = lambda a: X0 + (a - 0.015) * (X1 - X0) / 0.045
parts = [f'<line x1="{X0}" y1="{Y0}" x2="{X1}" y2="{Y0}" stroke="{INK}" stroke-width="3"/>']
for i, c in enumerate(hist):
    x0, x1 = xa(edges[i]), xa(edges[i + 1])
    parts.append(f'<rect x="{x0 + 1:.1f}" y="{Y0 - c * YS:.1f}" width="{x1 - x0 - 2:.1f}" height="{c * YS:.1f}" fill="{LINE}"/>')
for a in (0.02, 0.03, 0.04, 0.05, 0.06):
    parts.append(t(xa(a), Y0 + 36, f"{a:.2f}", 26, SOFT, "middle"))
parts.append(t((X0 + X1) / 2, Y0 + 80, "AUC de cada motivo con las mismas candidatas →", 28, SOFT, "middle"))
mu, sd = SP["auc_mean"], SP["auc_sd"]
for a, lab, col, yy in ((mu, f"promedio {mu:.3f}", INK, 40), (mu + 3 * sd, "NES = 3", MONO, 40)):
    parts.append(f'<line x1="{xa(a):.1f}" y1="{Y0}" x2="{xa(a):.1f}" y2="{yy + 14}" stroke="{col}" stroke-width="3" stroke-dasharray="8 6"/>')
    parts.append(t(xa(a), yy, lab, 28, col, "middle", 600))
xp = xa(SP["auc"])
parts.append(f'<circle cx="{xp:.1f}" cy="{Y0 - 10}" r="10" fill="{B}"/>')
parts.append(f'<line x1="{xp:.1f}" y1="{Y0 - 26}" x2="{xp:.1f}" y2="150" stroke="{B}" stroke-width="3"/>')
parts.append(t(xp, 104, "PU.1", 30, B, "middle", 700))
parts.append(t(xp, 138, f"NES {SP['nes']:.1f}", 28, B, "middle", 600))
fig = svg(f"Histograma real de los AUC de los 5,876 motivos para las candidatas de SPI1. Promedio {mu:.3f}; el umbral NES 3 está en {mu + 3 * sd:.3f}; el motivo de PU.1 tiene el AUC más alto, {SP['auc']:.3f}, con NES {SP['nes']:.1f}.", 800, 860, 490, "\n".join(parts), vy=10)
out["nes"] = slide(
    "nes", "Paso ② por dentro · 5 de 6", "NES: ¿este motivo destaca frente a los demás?",
    '<div class="row" style="gap:56px; align-items:center;">\n' + figure(fig, 848, 0) + "\n" + steps([
        f"¿{SP['auc']:.3f} es mucho? Depende de lo que logran los otros motivos con las mismas candidatas.",
        f"Promedio de los 5,876 AUC: {mu:.3f}. Desviación estándar: {sd:.4f}.",
        f"NES = (AUC − promedio) / desviación = ({SP['auc']:.3f} − {mu:.3f}) / {sd:.4f} ≈ <b>{SP['nes']:.1f}</b>.",
        f"Pasa si NES ≥ 3. Aquí pasan solo {SP['n_motifs_nes3']} de 5,876, y PU.1 es el primero.",
    ]) + "\n</div>\n" + callout("NES alto = las candidatas tienen ese motivo <b>mucho más</b> de lo esperado por azar."),
    f"""Minutos 54–56. Histograma real: un AUC por motivo, las mismas {fmt(SP['in_db'])} candidatas de SPI1.
NES (normalized enrichment score) = cuántas desviaciones estándar está un motivo por encima del promedio de todos. Valores exactos: AUC {SP['auc']:.5f}, promedio {mu:.5f}, desviación {sd:.5f} → NES {SP['nes']:.3f} (igual al valor de regulons.csv del servidor).
Umbral 3 (opción nes_threshold): de 5,876 motivos, {SP['n_motifs_nes3']} lo pasan, y el de PU.1 queda en primer lugar. Es decir: de todas las "palabras" posibles, la más sobrerrepresentada cerca de las candidatas de SPI1 es justo la palabra de SPI1. Esa coincidencia es la evidencia.
Si el motivo con NES alto fuera de otro TF, no cuenta para SPI1 (se revisa en la tabla de anotación, siguiente slide).
Contraste para después (slide de la poda): el motivo de EBF1 pasa apenas, con NES 3.37.""")

ram = SP["rank_at_max"]
X0, X1, Y0, YS = 110, 820, 420, 0.45
ia = round(ram / 25)
fill_between = (f"M{xr(0)},{yr(cur[0])} " + " ".join(f"L{xr(i * 25):.1f},{yr(cur[i]):.1f}" for i in range(ia + 1))
                + " " + " ".join(f"L{xr(i * 25):.1f},{yr(max(avg2[i], 0)):.1f}" for i in range(ia, -1, -1)) + " Z")
parts = ejes() + [
    f'<path d="{fill_between}" fill="{B}" opacity="0.15"/>',
    linea(avg2, stroke=SOFT, stroke_width=4, stroke_dasharray="10 8"),
    linea(cur, stroke=B, stroke_width=5),
    f'<line x1="{xr(ram):.1f}" y1="{Y0}" x2="{xr(ram):.1f}" y2="60" stroke="{MONO}" stroke-width="4"/>',
    t(X0 + 20, 92, f"corte (línea roja): lugar {fmt(ram + 1)}", 28, MONO, "start", 700),
    t(X0 + 20, 128, f"{SP['leading_edge_n']} candidatas antes del corte", 28, MONO, "start", 600),
    t(xr(2000), yr(cur[80]) - 30, "motivo de PU.1", 30, B, "end", 700),
    t(xr(2600), yr(avg2[104]) + 44, "lo esperado", 28, SOFT, "start", 600),
]
fig = svg(f"Curva real de PU.1 contra lo esperado (promedio más dos desviaciones). La mayor separación está en el lugar {fmt(ram + 1)}; las {SP['leading_edge_n']} candidatas que aparecen antes se quedan.", 800, 860, 466, "\n".join(parts), vy=44)
out["podacurva"] = slide(
    "podacurva", "Paso ② por dentro · 6 de 6", "La poda: quedarse con las que tienen el motivo",
    '<div class="row" style="gap:56px; align-items:center;">\n' + figure(fig, 848, 0) + "\n" + steps([
        "¿El motivo es del mismo TF? La tabla lo confirma: este motivo está anotado directamente a SPI1.",
        "Se busca dónde la curva <b>más se separa</b> de lo esperado (promedio + 2 desviaciones).",
        f"Se quedan solo las candidatas antes de ese punto: <b>{SP['leading_edge_n']} de {fmt(SP['in_db'])}</b>, entre ellas CSF1R, LYZ y S100A8.",
        "Se repite con cada motivo y cada grupo de candidatas; el regulón junta todo. SPI1(+): 3,106 genes.",
    ]) + "\n</div>\n" + callout("Regulón = candidatas de GRNBoost2 que <b>además</b> tienen el motivo de su TF."),
    f"""Minutos 56–58. Curva real (la misma de antes) contra la línea de "lo esperado": el promedio de los 5,876 motivos más dos desviaciones estándar en cada posición.
El punto de corte es donde la curva de PU.1 más se aleja de esa línea: lugar {fmt(ram + 1)}. Las candidatas que aparecen antes de ese punto ("leading edge") son las que se quedan: {SP['leading_edge_n']} de {fmt(SP['in_db'])}. Coincide con la fila del servidor (RankAtMax {ram}, {SP['leading_edge_n']} genes).
Requisito previo: que el motivo esté anotado al TF del módulo (directamente, por parecido o por ortología). Un motivo enriquecido pero de otro TF no aporta dianas a SPI1.
El regulón final SPI1(+) es la unión de las dianas que sobrevivieron con todos sus motivos enriquecidos, en sus seis grupos de candidatas y en las dos bases: 3,106 genes, entre ellos LYZ, FCN1, CSF1R y S100A8.
Lo que se tira: genes que GRNBoost2 asociaba a SPI1 pero que no tienen su motivo cerca. Pueden ser regulación indirecta (SPI1 enciende otro TF que enciende al gen) o simple coincidencia.""")

for sid, html in out.items():
    (SLIDES / f"{sid}.html").write_text(html)
    print("escrita", sid)
