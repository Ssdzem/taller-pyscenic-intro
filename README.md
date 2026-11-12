# Taller pySCENIC: introducción en 2 horas

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Ssdzem/taller-pyscenic-intro/blob/main/notebooks/Taller_pySCENIC_2h.ipynb)

- **Para quién:** estudiantes de biología o biomedicina **sin experiencia** en programación, Conda ni terminal.
- **Objetivo:** reconocer los tres pasos de pySCENIC (**GRNBoost2 → cisTarget → AUCell**), explorar asociaciones TF–gen y regulones reales, graficar su actividad y discutir **qué se puede y qué no se puede concluir**.
- **Pregunta biológica:** ¿qué regulones muestran distinta actividad entre **monocitos CD14+, células B y células NK** de sangre humana?
- **Duración:** 2 horas.
- **Requisitos:** cuenta de Google, navegador actualizado e internet. **No hay que instalar nada.**

## Cómo empezar

1. Haz clic en **Open in Colab** (arriba).
2. *(Opcional, para guardar tus respuestas)* **Archivo → Guardar una copia en Drive**.
3. Ejecuta las celdas **en orden** con ▶ o `Shift + Enter`. Los datos se descargan solos.

## Qué contiene

| Carpeta | Contenido |
|---|---|
| `notebooks/` | `Taller_pySCENIC_2h.ipynb`: el único notebook del taller |
| `data/` | Matriz de expresión demo (900 células × 10,120 genes, conteos crudos) y metadatos |
| `results/` | Resultados reales de pySCENIC 0.12.1: asociaciones TF–gen, 272 regulones, matriz AUCell y resúmenes |
| `figures/` | Figuras de referencia (red de un regulón, heatmap, t-SNE) en PNG y PDF, por si Colab falla |
| `scripts/` | Código para reproducir todo el análisis |
| `docs/METODOS.md` | Procedencia de los datos, parámetros, versiones y verificaciones |

Los resultados se calcularon previamente en un servidor; en clase **no** se ejecuta el pipeline completo (aunque la sección opcional 4 calcula AUCell en vivo).

**¿Falla Colab o internet?** Descarga el respaldo completo (notebook, datos, resultados y figuras): [`taller-pyscenic-respaldo.zip`](https://github.com/Ssdzem/taller-pyscenic-intro/releases/download/v1.0/taller-pyscenic-respaldo.zip) (5.6 MB).

## Datos y créditos

- Datos: *10k PBMCs from a Healthy Donor (v3 chemistry)*, **10x Genomics**, bajo licencia [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). **Modificados:** control de calidad, anotación de poblaciones, submuestreo a 900 células, filtrado de genes y análisis con pySCENIC (detalles en [`docs/METODOS.md`](docs/METODOS.md)).
- Software: pySCENIC (Aibar *et al.* 2017, *Nat Methods*; Van de Sande *et al.* 2020, *Nat Protoc*) y CellTypist (Domínguez Conde *et al.* 2022, *Science*).
- Licencia de este material: código bajo **MIT**; texto, figuras y datos derivados bajo **CC BY 4.0** (ver [`LICENSE`](LICENSE)).
