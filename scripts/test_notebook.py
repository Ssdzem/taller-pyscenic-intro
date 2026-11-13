"""Ejecuta el notebook de principio a fin como "Run all" (paso 36/42, ensayo local).
Uso (en localhost/taller-colabtest:py312):  python test_notebook.py NOTEBOOK SALIDA.ipynb
Lee los datos desde TALLER_LOCAL (copia local del repo) si está definido; si no, desde GitHub."""
import os, sys, time, json, nbformat
from nbclient import NotebookClient

nb = nbformat.read(sys.argv[1], as_version=4)
client = NotebookClient(nb, timeout=int(os.environ.get("NB_TIMEOUT", 900)), kernel_name="python3", allow_errors=True)
t0 = time.time()
with client.setup_kernel():
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        t = time.time()
        client.execute_cell(cell, i)
        errs = [o for o in cell.outputs if o.output_type == "error"]
        imgs = sum("image/png" in o.get("data", {}) for o in cell.outputs)
        first = cell.source.splitlines()[0][:60]
        print(f"[{i:2d}] {time.time()-t:5.1f}s  {'❌ ' + errs[0].ename + ': ' + errs[0].evalue[:120] if errs else '✅'}"
              f"{f'  ({imgs} fig)' if imgs else ''}  | {first}", flush=True)
nbformat.write(nb, sys.argv[2])
n_err = sum(o.output_type == "error" for c in nb.cells if c.cell_type == "code" for o in c.outputs)
print(f"\nTOTAL {time.time()-t0:.0f}s | celdas con error: {n_err}")
