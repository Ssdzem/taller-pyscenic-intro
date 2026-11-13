"""Fuentes de los recursos de aertslab, compartidas por build_notebook.py y build_notebook_casa.py.

Desde el 2026-11-13 resources.aertslab.org (y su espejo, que vive en el mismo servidor) responde con
un certificado de otro dominio y HTTP 404. Por eso los notebooks prueban primero una copia en el
Drive del instructor, compartida "cualquier persona con el enlace", y después aertslab.

Para cambiar los enlaces edita scripts/recursos_drive.json y regenera ambos notebooks.
"""
import json
from pathlib import Path

_JSON = Path(__file__).with_name("recursos_drive.json")
ENLACES_DRIVE = {k: v for k, v in json.loads(_JSON.read_text()).items() if not k.startswith("_")}

AERTS = "https://resources.aertslab.org/cistarget"
ESPEJO = "https://resources-mirror.aertslab.org/cistarget"

# Código que se inserta tal cual en los notebooks (no es un f-string: las llaves son literales).
FUNCIONES_DESCARGA = r'''
def id_drive(enlace):
    """Acepta un enlace de Drive (…/file/d/ID/view, …?id=ID) o el ID solo."""
    m = re.search(r"(?:/d/|[?&]id=)([\w-]{20,})", enlace)
    return m.group(1) if m else enlace.strip()

def sha1(ruta):
    h = hashlib.sha1()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()

def bajar(destino, urls, enlace_drive="", esperado=None, intentos=2):
    """Prueba cada fuente en orden (Drive y luego las URLs) hasta que una dé el archivo correcto.
    Devuelve el nombre de la fuente que funcionó."""
    errores = []
    def correcto():
        if os.path.exists(destino) and (esperado is None or sha1(destino) == esperado):
            return True
        if os.path.exists(destino):
            os.remove(destino)
            errores.append("llegó dañado (SHA1 distinto)")
        return False
    if enlace_drive:
        try:
            try:
                import gdown
            except ImportError:   # Colab ya lo trae; esto es por si acaso
                sin_pp = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
                subprocess.run([sys.executable, "-m", "pip", "install", "-q", "gdown"], check=True, env=sin_pp)
                import gdown
            gdown.download(id=id_drive(enlace_drive), output=destino, quiet=True)
            if correcto():
                return "Drive del instructor"
        except Exception as e:
            errores.append(f"Drive: {str(e).strip().splitlines()[0][:150]}")
    for url in urls:
        servidor = url.split("/")[2]
        for intento in range(intentos):
            try:
                with urllib.request.urlopen(url, timeout=30) as r, open(destino, "wb") as f:
                    shutil.copyfileobj(r, f, 1 << 20)
                if correcto():
                    return servidor
            except Exception as e:
                errores.append(f"{servidor}: {str(e)[:150]}")
            if intento + 1 < intentos:
                time.sleep(5)
    raise RuntimeError("ninguna fuente funcionó:\n     - " + "\n     - ".join(dict.fromkeys(errores)))
'''.strip("\n")
