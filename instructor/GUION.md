# Guion del instructor — Taller pySCENIC (120 min)

Material:

- **Diapositivas:** 26, con notas de orador (en las notas están las fuentes y los números exactos). Están en Claude, son privadas y se descargan como PDF o PPTX: <https://claude.ai/artifact/BUsry7HyPH8F3zh1gai9vr>. Fuentes locales en `instructor/slides/`.
- **Notebook:** <https://colab.research.google.com/github/Ssdzem/taller-pyscenic-intro/blob/main/notebooks/Taller_pySCENIC_2h.ipynb>
- **Respuestas:** `instructor/SOLUCIONES.md`. Se publica después del taller (paso 50).
- **Respaldo:** <https://github.com/Ssdzem/taller-pyscenic-intro/releases/latest/download/taller-pyscenic-respaldo.zip>, más las figuras PNG y PDF de `figures/`.

---

## Antes del día

| Cuándo | Qué | Paso |
|---|---|---|
| 1 semana antes | Abrir el repo en ventana privada; *Open in Colab* → *Ejecutar todo* en sesión limpia, 2–3 veces. Anotar si la 4.2 (pySCENIC en vivo) funciona y cuánto tarda. | 41, 42 |
| 1 semana antes | Cronometrar el notebook con una persona novata. | 43 |
| 3–5 días antes | Enviar el mensaje a los alumnos (abajo). | 47 |
| 1–2 días antes | Simulacro de los primeros 10 minutos con otra cuenta de Google: enlace → Colab → Preparación → *Instalar pySCENIC* (✅ en ~30 s) → ver la tabla de poblaciones. | 48 |
| 1–2 días antes | **Copia de los recursos en tu Drive (obligatoria desde el 2026-11-13):** sube a una carpeta de tu Google Drive los 4 archivos de `~/repos/taller-pyscenic/resources/` en Stargazer: `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` (313 MB), `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` (99 MB), `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` (311 MB, solo para casa) y `allTFs_hg38.txt` (12 KB, solo para casa). Comparte la carpeta como "cualquier persona con el enlace (lector)". Claude pone los enlaces en `scripts/recursos_drive.json`, regenera los notebooks y hace push; el notebook prueba primero tu Drive, luego aertslab. | — |
| 1–2 días antes | Ensaya los dos bloques de teoría: lanzar `grn` y explicar las diapositivas 10–15 (~11 min); explicar la 16–17 (motivo y bases), lanzar `ctx` y explicar las 18–22 (~10 min). Si `ctx` termina antes, termina la explicación de todos modos. | — |
| El día antes | Corre la celda de recursos de la sección 2 en un Colab limpio: debe decir `desde: Drive del instructor` en los dos archivos. | — |
| El día | Llevar el ZIP de respaldo en una USB y las figuras abiertas en otra pestaña. | — |

**Mensaje para los alumnos (paso 47).** Solo el enlace; nada de Conda:

> Hola. Para el taller de pySCENIC solo necesitas una **cuenta de Google** y un **navegador actualizado** (Chrome, Firefox o Edge). **No hay que instalar nada.**
> El día del taller abriremos este enlace: https://github.com/Ssdzem/taller-pyscenic-intro
> Si quieres, ábrelo antes y comprueba que puedes iniciar sesión en Google Colab.

---

## Minuto a minuto

**La idea:** 16 minutos de teoría al inicio (la red explicada en 4 diapositivas: interruptor, red, circuitos, estados estables) y, después, ir y venir entre el notebook y las diapositivas. GRNBoost2 y cisTarget se explican **una idea por diapositiva** (6 diapositivas cada uno), pensando en quien no sabe qué es un motivo. En los pasos 1 y 2 los alumnos **lanzan el comando y, mientras corre**, tú explicas el paso por dentro. En el paso 3 explicas primero, porque `aucell` tarda segundos.

**Todas las gráficas de las diapositivas 10–14 y 17–21 son datos reales** (y las aristas de la 3 y la 4 salen de nuestros regulones) de nuestras 900 células, no dibujos: la recta de TBX21, los cortes y el árbol de NKG7, el boosting de dos células, la importancia (reproduce la del servidor: TBX21 → NKG7 = 66.09) y la curva, el NES y la poda de SPI1 con el motivo de PU.1 (reproducen la fila del servidor: NES 4.31, corte en 4,495, 740 genes). Scripts y datos en `instructor/slides/datos/`.

| Min | Dónde | Diapositiva | Qué pasa |
|---|---|---|---|
| 0–2 | Proyector | 1 · Portada | La pregunta: un genoma, tres tipos celulares |
| 2–4 | Proyector | 2 · Un gen tiene un interruptor | TF, motivo, diana; activador y represor |
| 4–6 | Proyector | 3 · ¿Qué es una GRN? | Nodos, aristas con dirección y signo, regulón; red real |
| 6–8 | Proyector | 4 · Circuitos | Autorregulación, hacia adelante, retroalimentación |
| 8–10 | Proyector | 5 · Estados estables | Cada tipo celular es un valle; C/EBPα convierte B en macrófago |
| 10–12 | Proyector | 6 · ¿Para qué sirve? | Identidad, Yamanaka, GATA2, hipótesis |
| 12–14 | Proyector | 7 · ¿Por qué pySCENIC? | Frente a expresión diferencial y coexpresión |
| 14–16 | Proyector | 8 · Tres pasos | El mapa del taller |
| 16–22 | Notebook | 9 · Colab (proyectada) | Sección 0 + instalar pySCENIC |
| 22–24 | Notebook | — | Lanzan preparación + `grn` |
| 24–35 | Proyector | 10 · Un corte en la recta, 11 · ¿Dónde cortar?, 12 · Árbol de 3 niveles, 13 · Boosting, 14 · Importancia, 15 · La tabla | Mientras corre `grn` |
| 35–44 | Notebook | — | Ejercicio 1 |
| 44–48 | Proyector | 16 · ¿Qué es un motivo?, 17 · Las bases | Antes lanzan la celda de recursos (~20 s) |
| 48–50 | Notebook | — | Lanzan `ctx` |
| 50–60 | Proyector | 18 · Ranking, 19 · Curva y AUC, 20 · NES, 21 · Poda, 22 · ¿Quién sobrevive? | Mientras corre `ctx` |
| 60–68 | Notebook | — | Ejercicio 2 |
| 68–73 | Proyector | 23 · AUCell por dentro | Antes de `aucell` |
| 73–86 | Notebook | — | `aucell` + Ejercicio 3 |
| 86–91 | Proyector | 24 · Los protagonistas | Biología de los TF |
| 91–108 | Notebook | 25 · Parejas | Sección 5 |
| 108–120 | Proyector | 26 · Límites | Sección 6 y cierre |

### 0–16 · Teoría y método (diapositivas 1–8)

- **1 · Portada (2 min):** una persona, un genoma, tres tipos celulares: ¿qué los hace distintos? Hoy no se instala nada en su computadora, pero **sí corren los tres comandos** en Colab.
- **2–5 · Qué es una red de regulación, de lo concreto a lo abstracto (2 min cada una):**
  - **2 · Un gen tiene un interruptor:** el TF reconoce el motivo antes del gen; si es activador el gen se transcribe, si es represor queda apagado. Ejemplo nuestro: TBX21 → NKG7.
  - **3 · ¿Qué es una GRN?:** el mapa de quién regula a quién. Nodos = genes; aristas = flechas con dirección y signo; regulón = un TF y sus flechas de salida. El dibujo son aristas reales de nuestros regulones, incluidas flechas entre TF (SPI1 ⇄ CEBPB, SPIB → EBF1, TBX21 → EOMES). Nuestra red: 8,929 genes y 46,339 flechas, todas (+), porque pySCENIC solo reporta activaciones. Cierra con *"es un mapa propuesto, no medido"*.
  - **4 · Circuitos que se repiten:** autorregulación (memoria; PU.1 activa su propio gen), hacia adelante (filtro; SPI1 → CEBPB → KLF4 y SPI1 → KLF4, hay uno real por población) y retroalimentación (candado si se activan, elección si se reprimen; SPI1 ⇄ CEBPB). Aclara que "motivo de red" no es el motivo de ADN.
  - **5 · Estados estables:** los circuitos dejan pocas combinaciones estables de TF, los valles del paisaje de Waddington; cada tipo celular es una. Con C/EBPα forzado una célula B se vuelve macrófago (Xie et al. 2004). Por eso esperamos regulones distintos en cada población.
  - Ojo con la autorregulación en nuestros datos: 153 de 272 regulones contienen a su propio TF, pero pySCENIC lo mete a propósito (peso 1.0). No es evidencia de coexpresión; está en las notas de la 4.
- **6 · ¿Para qué sirve? (2 min):** una frase por tarjeta; la identidad ya la vieron con los valles.
  - Identidad: mismo ADN, distinto programa.
  - Reguladores maestros: 4 TF reprograman una célula (Takahashi y Yamanaka, 2006).
  - Enfermedad: la deficiencia de GATA2 deja muy pocos monocitos, células B y NK, justo nuestras poblaciones.
  - Hipótesis: qué TF perturbar.
- **7 · ¿Por qué pySCENIC? (2 min):** la expresión diferencial dice *qué* cambia, no *quién* lo controla; la coexpresión da miles de pares indirectos y sin signo. pySCENIC pide dos evidencias (coexpresión + motivo) y mide la actividad en cada célula. Dato del recuadro: el gen EBF1 se detecta en 31 % de las células B, pero su regulón está activo en 86 % (259 de 300).
- **Pregunta al grupo:** "Si dos genes suben y bajan juntos, ¿uno regula al otro?" Recoge 2–3 respuestas y no las corrijas: se contesta con ANXA1 (min 34, diapositiva 15) y RPS4X (min 58, diapositiva 22).
- **8 · Tres pasos (2 min):** el mapa (2,339,328 → 272 → 900 × 272). Avisa que cada paso tendrá su diapositiva "por dentro".

### 16–22 · Abrir Colab (diapositiva 9 proyectada)

- **Meta:** todos ven los 7 ✅ de Preparación, el t-SNE de las tres poblaciones y **"✅ pySCENIC 0.12.1 listo"** (~30 s).
- **Circula por la sala. Problemas típicos:**
  - No iniciaron sesión.
  - Aviso "no lo creó Google": elegir *Ejecutar de todos modos*.
  - Ejecutaron una celda de abajo primero: *Entorno de ejecución → Ejecutar anteriores*.
- **Si alguien no puede:** que trabaje en pareja desde ya.

### 22–44 · Paso 1: GRNBoost2

- **22–24, notebook:** que ejecuten la celda ▶ de preparación y la de `arboreto_with_multiprocessing.py` (~2 min). Señala de pasada el bloque 🖥️ de flags (`--seed`, `--num_workers`) y por qué no `pyscenic grn` (versiones de `dask`).
- **24–35, diapositivas 10–15, mientras corre.** Un solo ejemplo en las cuatro primeras: predecir **NKG7** (gen de citotoxicidad NK).
  - **10 · Un corte en la recta (2 min):** cada punto es una célula: en la recta, su TBX21 (0 a 6 UMI); hacia arriba, su NKG7. La pregunta "¿TBX21 ≥ 1?" es un corte: 675 células del lado no y 225 del lado sí. **La predicción de cada lado es el promedio de NKG7 de sus células:** 3.8 y 32.6. Error = promedio de (real − predicción)²: 277 sin cortar, 121 con el corte.
  - **11 · ¿Dónde cortar? (1 min):** fuerza bruta. Los 6 cortes de TBX21 (≥ 1 gana con 121; ≥ 2 ya da 191) y el mejor corte de cada uno de los 1,018 TF (TBX21 121, RORA 145, KLRG1 162…; SPI1 197, EBF1 264). Gana el que deja menos error.
  - **12 · Árbol de 3 niveles (2 min):** dentro de cada lado se repite lo mismo. Árbol real: TBX21 → KLRG1 / UBB → RORA, MEF2C, RPL6, SNRPB2 → 8 grupos, de 0.7 (592 Mono y B) a 58.6 (8 NK). Error 277 → 121 → 52. Rescata a las NK con dropout de TBX21 (KLRG1, MEF2C). UBB y RPL6 están en la lista de "TF" aunque no lo son.
  - **13 · Boosting (2 min):** se empieza con el promedio y cada árbol (de 3 niveles, como el de la 12) corrige 1 % de lo que falta. Una célula NK real (38 UMI): 11.0 → 23.9 con 100 árboles → 32.9 con 482. Un monocito (0): 11.0 → 0.2. Cada árbol ve el 90 % de las células y, en cada pregunta, solo escoge entre 101 TF al azar; para cuando deja de mejorar: 482 árboles para NKG7, y en 60 genes al azar de 25 (ningún TF lo predice) a 481, mediana 64.
  - **14 · Importancia (2 min):** la baja del error de cada pregunta se le apunta al TF y se suma en todos los árboles. Tabla real: TBX21 sale en 105 de 482 árboles y sus bajas suman 1,802; ÷ 27.3 (lo que baja un árbol promedio, solo escala) = 66.1. Reproducido con la semilla del servidor: idéntico a la tabla.
  - **15 · La tabla (3 min):** NKG7 lo predicen TBX21, RORA y MYBL1 (ρ positiva). CD79A lo predice **ANXA1 con ρ = −0.73**: ANXA1 está en 98 % de los monocitos y solo en 13 % de las B. **Contesta la pregunta del minuto 12:** la importancia no tiene signo ni dirección.
- **35–44, notebook:**
  - Comparan su top 15 de TBX21 con el del servidor (pregunta 4).
  - Exploran la tabla del servidor y responden las preguntas 1–3, unos 7 minutos.
  - **Puesta en común (4 min):** 3,299 candidatos de TBX21, muchos débiles (solo 871 con importancia > 1); el top 15 es el programa citotóxico NK (PRF1, GNLY, GZMA/B, NKG7). **RPS4X:** las aristas más fuertes de la red son de una proteína ribosomal. "¿Regulador?" Déjalo abierto.

### 44–68 · Paso 2: cisTarget

- **44, notebook (30 s):** que lancen la celda de recursos (~20 s desde tu Drive; si falla, aertslab). Baja justo las bases que vas a explicar.
- **44–48, diapositivas 16–17, antes de `ctx`:**
  - **16 · ¿Qué es un motivo? (2 min):** desde cero. Un TF reconoce una palabra corta del ADN (6–20 letras), no exacta. PU.1 (proteína de SPI1): AGAGGAAGTG (HOMER) y AAAAAGCGGAAGTN (Taipale 2013), mismo núcleo GGAA (familia ETS). A SPI1 se le asignan 390 motivos en la tabla. Tener el motivo cerca = *puede* regularlo.
  - **17 · Las bases (2 min):** el genoma ya escaneado por el laboratorio Aerts. Dos ventanas: promotor (−500/+100 pb) y ±10 kb. Por motivo, los 27,015 genes ordenados = un ranking. Ejemplo real: con el motivo de PU.1, CSF1R está en el lugar 100 y CD79A en el 21,861. Otra tabla traduce motivo → TF.
- **48–50, notebook:** lanzan `pyscenic ctx` (~4 min).
- **50–58, diapositivas 18–21, mientras corre.** Un solo ejemplo: las candidatas de SPI1 y el motivo de PU.1, base del promotor.
  - **18 · El ranking (2 min):** candidatas = 2,694 posibles dianas de SPI1 según GRNBoost2 (pares en el 10 % de importancias más altas, con correlación positiva; aquí cae ANXA1 → CD79A). En el primer 5 % del ranking (1,351 genes) caen 280; al azar, 130. CSF1R (lugar 100), LYZ (762) y S100A8 (1,320) son candidatas.
  - **19 · Curva y AUC (2 min):** bajar por el ranking contando candidatas. En el lugar 1,351, PU.1 lleva 280 y el motivo promedio 164. AUC = área en ese primer 5 %: 0.057.
  - **20 · NES (2 min):** promedio de los 5,876 AUC 0.032, desviación 0.0059 → NES 4.3. Solo 7 motivos pasan NES ≥ 3, y PU.1 es el primero: la palabra más sobrerrepresentada es la del propio SPI1.
  - **21 · Poda (2 min):** el motivo está anotado directamente a SPI1; corte donde la curva más se separa de promedio + 2 desviaciones (lugar 4,496) → se quedan 740 de 2,600, entre ellas CSF1R, LYZ y S100A8. El regulón junta todos los motivos y grupos: SPI1(+) 3,106. **Regulón = candidatas que además tienen el motivo de su TF.**
- **58–60, diapositiva 22 (¿quién sobrevive?), con datos del servidor:**
  - SPI1 3,106 dianas; TBX21 1,014; EBF1 solo 14 (un motivo, NES 3.37).
  - MAFB: 12 de 13 motivos de la base de ±10 kb → probablemente no les aparezca.
  - PAX5: 45 motivos propios y ningún NES ≥ 3. **Ausencia de regulón ≠ ausencia de regulación.**
  - ANXA1 y RPS4X: sin motivo anotado → sin regulón. **Contesta lo de RPS4X.**
  - Dos bases (promotor y ±10 kb); en clase solo una y 2,829 genes → ~150 regulones.
- **60–68, notebook:** comparan su número de regulones con los 272 (pregunta 4), y ejecutan la comparación de TBX21, la red de MAFB(+) y la tabla de 13 TF.
  - **Puesta en común:** TBX21 conserva el 30.7 % (1,014 de 3,299). Robustez, con honestidad: en las pruebas, MAFB(+) apareció en 2 de 4 corridas (servidor y dataset entero) y PAX5(+) en una; EBF1, TBX21, FOS y SPI1 conservan su población (celda opcional de la sección 3).
- **Si vas tarde:** no corren `ctx` (lo muestras tú) y la diapositiva 22 se vuelve un minuto.
- **Solo para ti (por si preguntan por EBF1):** al reproducir cisTarget encontramos que pySCENIC 0.12.1 puede asignar a un motivo las dianas de otro motivo del mismo grupo de candidatas. Las 16 dianas de EBF1(+) salen de la curva de un motivo de IRF6, no de la del motivo de EBF1 (que daría 18). Los NES sí son correctos. Detalle en `STATUS.md`, paso 55.

### 68–86 · Paso 3: AUCell

- **68–73, diapositiva 23 (AUCell por dentro), ANTES de correr:**
  - Ranking por célula: solo importa el orden dentro de la célula; no hay que normalizar.
  - Recorrer el 5 % de arriba (506 genes) contando dianas del regulón: la misma curva de recuperación, ahora por célula.
  - AUC de 0 a 1. Dato real: EBF1(+) con AUC > 0 en 259 de 300 células B y en 134 de 600 de las otras.
- **73–86, notebook:**
  - **▶ `aucell`** (~10 s): con los regulones del servidor debe salir "✅ idéntico". Si corrieron `ctx`, la celda opcional calcula la actividad de *sus* regulones. Recuerda unir por ID de célula.
  - **Heatmap y boxplots, con las medianas:** MAFB(+) en monocitos 0.311; EBF1(+) en B 0.049 y 0 en las demás; TBX21(+) en NK 0.075.
  - **Expresión frente a actividad:** FOS se expresa en todas las poblaciones, pero FOS(+) es activo sobre todo en monocitos (0.694). EBF1 no se detecta en 206 de 300 células B y aun así su regulón está activo.
  - **SPI1(+):** con 3,106 genes, la diferencia absoluta es mínima (0.060 frente a 0.044). El z-score del heatmap la exagera.
  - **Pregunta 4, respuesta temprana:** abre la discusión de artefactos técnicos.
  - **Sección 4 (opcional, AUCell desde cero):** para quien vaya adelantado; es la diapositiva 23 en código.

### 86–91 · Los protagonistas (diapositiva 24)

- Una columna por población; lee los TF y su función en una frase. Las fuentes están en las notas.
  - **Monocitos:** SPI1 (PU.1) maestro mieloide; CEBPB/CEBPA; KLF4; MAFB.
  - **B:** EBF1 enciende el programa y PAX5 lo fija; POU2AF1 coactiva genes de anticuerpos; SPIB.
  - **NK:** T-bet (TBX21) y EOMES, maduración y citotoxicidad; RUNX3; PRDM1.
- Cierra: "ahora elijan un regulón sabiendo qué esperar, y vean si los datos coinciden".
- Recuerda que FOS(+) sale "de monocitos" pero es sospechoso de artefacto.

### 91–108 · Parejas (diapositiva 25)

- Elegir un regulón de la tabla de específicos (o EBF1), ejecutar `ficha()` y llenar la plantilla.
- **Circula. La corrección típica:** la actividad AUC y la etiqueta de población **no** son evidencia observada.
- **Minutos 103–108:** 2–3 parejas leen su hipótesis.
- **Si eligen regulones pequeños o raros** (STAT4 con 7 genes, HOXB4 en NK, ZNF701): son buenos para hablar de fragilidad.

### 108–120 · Limitaciones y cierre (diapositiva 26)

- Qué sí y qué no. Un ejemplo de hoy para cada "no": RPS4X, PAX5, un donante, FOS(+) y la variabilidad.
- Cierra volviendo a la diapositiva 3: regulón = hipótesis.
- Recursos para seguir (sección 6 del notebook): el resumen de los tres comandos y el **notebook para casa** (`Pipeline_completo_casa.ipynb`, ~30 min con la matriz demo o ~2–3 h con la completa). Avisa que compartirás las respuestas sugeridas.

---

## Si algo sale mal

| Falla | Qué hacer |
|---|---|
| Colab no abre o no hay cuenta de Google | Trabajar en pareja. Mostrar en el proyector tu propia ejecución. |
| La red de la sala es lenta | La descarga es de ~10 MB y se hace una sola vez. Si falla: ZIP de respaldo por USB → subirlo a Colab (icono de carpeta) o abrir las figuras. |
| Se cae internet del todo | Seguir con las figuras PDF/PNG de `figures/` y la discusión biológica. Las diapositivas descargadas en PDF funcionan sin red. |
| Falla *Instalar pySCENIC* | Las celdas ▶ muestran un aviso y todo sigue con los resultados del servidor. Muestra tu propia ejecución en el proyector. |
| No bajan los recursos (Drive dice "demasiadas descargas" y aertslab está caído) | Que agreguen tu carpeta a su Drive (clic derecho → *Organizar* → *Agregar acceso directo*) y usen la Opción 2 (`RUTA_DRIVE = "/content/drive/MyDrive/<carpeta>"`). Si tampoco, `ctx` solo en tu proyector o para casa. |
| Vas tarde | Recorta, en este orden: la sección 4 completa, la diapositiva 22 (¿quién sobrevive?) a un minuto, la 14 (importancia) a un minuto, el `ctx` en vivo (tú lo muestras), la red de MAFB(+), la diapositiva 6 (¿para qué sirve?) a un minuto, la 4 (circuitos) a un minuto y la pregunta 4 del Ejercicio 3. Las 10–13 y 18–21 no se recortan: ocurren mientras corren los comandos. La 16 (motivo) tampoco: sin ella no se entiende cisTarget. **Nunca** recortes la sección 5 ni el cierre. |
