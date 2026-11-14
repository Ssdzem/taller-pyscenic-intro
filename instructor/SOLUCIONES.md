# Respuestas sugeridas — Taller pySCENIC (2 h)

> **Para el instructor.** No compartir antes del taller: los alumnos verían las respuestas. Paso 50: publicarlo en el repositorio **después** de la sesión.
>
> Todas las cifras salen de los mismos archivos que usa el notebook y se recalculan con `scripts/46_answers.py` (salida en `logs/46_answers.log`). Son las cifras que verá cualquier alumno que ejecute el notebook sin cambiar nada.

**Cómo leer las etiquetas de evidencia**

- 🔎 **Observado:** medido directamente en los datos (conteos UMI, cuántas células expresan un gen).
- 🧮 **Inferido:** lo propone un algoritmo o una decisión nuestra:
  - asociación TF–gen (GRNBoost2);
  - pertenencia al regulón (cisTarget);
  - etiqueta de población (marcadores + CellTypist).
- La **actividad AUC** es mixta: se calcula sobre la expresión medida, pero con un conjunto de genes inferido. Cuéntala como 🧮.
- 📚 **Conocimiento previo:** literatura, no estos datos.

---

## Respuesta corta a la pregunta del taller

> *¿Qué regulones muestran distinta actividad entre monocitos CD14+, células B y células NK?*

Regulones más específicos por población, según la tabla de la sección 5: diferencia de mediana de AUC respecto a la siguiente población.

| Monocitos CD14+ | Células B | Células NK |
|---|---|---|
| FOS(+) ⚠️ · **MAFB(+)** · TFEC(+) · JUNB(+) ⚠️ · FOSB(+) ⚠️ · NFE2(+) · CEBPD(+) · ZBTB43(+) | TFEB(+) · **POU2AF1(+)** · ZNF432(+) · TCF4(+) · BACH2(+) · MYC(+) · **EBF1(+)** · SPIB(+) | **STAT4(+)** · ZNF701(+) · PRDM1(+) · NFATC2(+) · HOXB4(+) · RUNX3(+) · EOMES(+) · GFI1(+) |

- **Coherentes con la biología conocida 📚:**
  - Monocitos: MAFB, CEBPD, TFEC.
  - Células B: EBF1, POU2AF1, SPIB, TCF4, BACH2.
  - Células NK: TBX21, EOMES, PRDM1, RUNX3, STAT4.
- ⚠️ **Robustez entre corridas** (comprobado al repetir el pipeline; ver Ejercicio 2, pregunta 4):
  - **MAFB(+) aparece solo en 2 de 4 corridas:** el servidor y la del dataset entero (10,280 células). No aparece en la reducida de clase ni en la demo con `arboreto_with_multiprocessing.py`.
  - EBF1(+), TBX21(+), SPI1(+) y FOS(+) aparecen en las cuatro, aunque POU2AF1(+) falta en la de 10,280 células.
  - Presenta a MAFB como el caso más claro **de esta corrida**, no como un hallazgo robusto.
- ⚠️ **Posible artefacto:** FOS, JUNB y FOSB son genes de respuesta temprana; ver Ejercicio 3, pregunta 4.
- **Regulones muy pequeños** (≤ 9 genes): STAT4(+) 7, TFEB(+) 7, ZNF432(+) 5, ZNF701(+) 8, MYC(+) 9, ZBTB43(+) 9. Son frágiles y muy sensibles a la variabilidad de GRNBoost2.
- **Poco conocidos en estos linajes:** HOXB4(+) en NK y ZNF701(+) en NK. Son buenos ejemplos de *no sobreinterpretar*.

---

## Sección 0 · Preparación

- 300 células por población (Mono_CD14, B, NK); 900 en total. 🔎 (conteo de filas)
- La población de cada célula es 🧮: la asignamos con marcadores canónicos sobre los clusters de 10x y la confirmamos con CellTypist. Solo se usan células donde ambos métodos coinciden: 96.7 % de los monocitos, 98.8 % de las B y 91.2 % de las NK.

---

## Ejercicio 1 · GRNBoost2

**1. ¿Cuántos genes diana candidatos tiene TBX21? ¿Es realista que regule directamente tantos?**

- **3,299 candidatos.** No es realista como regulación directa:
  - GRNBoost2 lista **todo** gen al que el TF ayudó a predecir algo, por poco que sea.
  - La distribución está muy sesgada: mediana de importancia 0.31, solo **871** dianas con importancia > 1, percentil 99 = 26.5 y máximo 66.1 (NKG7). Lo ve el histograma en escala log.
- 🧮 Las asociaciones son inferidas. Que un TF real tenga cientos o miles de dianas directas no es imposible 📚, pero esta lista mezcla efectos directos, indirectos y de composición celular.

**2. ¿Reconoces genes típicos de células NK en el top 15?**

- Sí. Top 15:
  - NKG7, CD247, PRF1, GNLY, GZMA, CST7, CTSW, ARL4C, HOPX, KLRD1, GZMM, IFITM2, APMAP, ABHD17A, GZMB.
- Entre ellos hay perforina (PRF1, puesto 3), granulisina (GNLY, 4), granzimas (GZMA 5, GZMB 15, GZMM), NKG7 (1), KLRD1 (CD94, 10) y CD247 (CD3ζ, también presente en NK). Es el programa citotóxico 📚.
- 🧮 Asociación inferida. 🔎 Que esos genes se expresan sobre todo en las NK de este dataset sí es observable.

**3. ¿Qué representa el peso? ¿Significa que TBX21 activa esos genes?**

- Es la **importancia** del TF en un modelo de árboles que predice la expresión del gen diana. Es una medida **relativa** de asociación.
- Lo que el peso **no** es:
  - **No tiene signo:** no distingue activar de reprimir.
  - **No tiene dirección:** el gen podría regular al TF, o ambos responder a un tercero.
  - **No es un p-valor.**
- Trampa clave: con tres poblaciones muy distintas, TBX21 funciona como un **marcador de célula NK**. Cualquier gen específico de NK se asocia con él, lo regule o no (confusión por tipo celular).

**4. Tu top 15 de TBX21 frente al del servidor: ¿por qué no coinciden del todo, si la semilla es la misma?**

- En la prueba en un entorno tipo Colab (2 núcleos) coincidieron **12 de 15**. El grupo puede ver cifras ligeramente distintas.
  - Tu corrida: NKG7, CD247, FCGR3A, ABHD17A, CTSW, GNLY, HOPX, PRF1, APOBEC3G, APMAP, CST7, GZMM, CMC1, KLRD1, GZMA.
  - Servidor: NKG7, CD247, PRF1, GNLY, GZMA, CST7, CTSW, ARL4C, HOPX, KLRD1, GZMM, IFITM2, APMAP, ABHD17A, GZMB.
- **Motivo principal:** la matriz reducida solo tiene **2,829 genes** (los 2,000 más variables + los TF).
  - Genes como ARL4C, IFITM2 o GZMB pueden no estar entre los "más variables", o competir con otros genes.
  - La importancia es **relativa**: depende de qué otros TF y genes hay en el modelo.
- **Motivos secundarios:** distinta implementación (`arboreto_with_multiprocessing.py` frente a `pyscenic grn` con dask) y versiones de librerías.
- **Mensaje:** el núcleo del programa citotóxico (NKG7, PRF1, GNLY, GZMA, KLRD1…) aparece en ambas. Es una señal **robusta**; los puestos exactos no lo son.

**RPS4X: ¿es un regulador? ¿Qué otra explicación hay?**

- Tiene **6,978** candidatos (el máximo de la tabla completa) y las tres asociaciones más fuertes de toda la red: RPL13 97.4, RPL41 95.3 y RPS3 93.7.
- El **70 %** de su top 100 son proteínas ribosomales (RPL…/RPS…).
- Es una proteína ribosomal, no un TF clásico 📚. Aparece porque la lista oficial de pySCENIC (`allTFs_hg38.txt`, 1,892 genes) la incluye.
- **Explicación más probable:** los genes ribosomales se coexpresan entre sí porque comparten programa transcripcional y responden al tamaño o actividad celular y a la profundidad de secuenciación. Es coexpresión **sin** regulación por RPS4X.
- Se resuelve en el Ejercicio 2.

---

## Ejercicio 2 · cisTarget

Contexto: **272 regulones**, todos activadores "(+)". Tamaño: mediana 41 genes, mínimo 3, máximo 3,106 (SPI1). 🧮

| TF | Candidatos GRNBoost2 | Genes en el regulón | % conservado |
|---|---:|---:|---:|
| RPS4X | 6,978 | **0** | 0 |
| SPI1 | 5,829 | 3,106 | 53.3 |
| CEBPB | 5,503 | 1,543 | 28.0 |
| CEBPA | 3,914 | 177 | 4.5 |
| RUNX3 | 3,721 | 501 | 13.5 |
| IRF8 | 3,423 | 243 | 7.1 |
| TBX21 | 3,299 | 1,014 | 30.7 |
| SPIB | 2,847 | 413 | 14.5 |
| PRDM1 | 2,386 | 554 | 23.2 |
| EOMES | 2,193 | 607 | 27.7 |
| PAX5 | 2,158 | **0** | 0 |
| POU2AF1 | 2,004 | 20 | 1.0 |
| EBF1 | 1,625 | 14 | 0.9 |

**1. ¿Qué fracción de los candidatos de TBX21 sobrevivió? ¿Qué criterio usó cisTarget?**

- **1,014 de 3,299 = 30.7 %.** Todos los genes del regulón vienen de los candidatos (es un subconjunto).
- **Criterio:**
  1. Agrupa las dianas en módulos.
  2. Para cada módulo, comprueba si las regiones reguladoras de esos genes (±10 kb o −500/+100 pb del inicio de transcripción) están enriquecidas en un motivo asociado al TF.
  3. Si el enriquecimiento es suficiente (**NES ≥ 3**), conserva solo los genes que puntúan alto para ese motivo.
- Motivo de TBX21 con mayor enriquecimiento: `metacluster_191.3`, NES 4.59, base −500/+100 pb.
- 🧮 El enriquecimiento de motivo es evidencia **independiente de la expresión** (secuencia de ADN), pero sigue siendo una predicción. No mide unión real.

**2. ¿Por qué RPS4X no tiene regulón? ¿Qué nos dice sobre usar solo coexpresión?**

- En la tabla de anotación motivo→TF usada (`motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`) hay **0 motivos anotados para RPS4X** 🔎 (dato del recurso). Sin motivo, cisTarget **no puede** construirle un regulón, por muy fuertes que sean sus asociaciones.
- **Moraleja:** la asociación más fuerte de toda la red era un falso regulador. La coexpresión sola no distingue regulación de programas compartidos.

**3. PAX5 no tiene regulón. ¿Significa que no regula esos genes?**

- **No.** PAX5 es esencial para la identidad B 📚 y sus candidatos son genes B de libro: IGHM, CD79B, AFF3, BCL11A, FAM129C, CBX7, IGHD, CD40, CXCR4, VPREB3.
- **Lo que sí sabemos 🔎:** la tabla de anotación tiene **125 motivos asociados a PAX5** (45 anotados directamente). El problema **no** es que falte su motivo.
- **Lo que inferimos:** ningún módulo de PAX5 alcanzó NES ≥ 3 en estas dos bases. Explicaciones posibles, ninguna demostrada aquí:
  - PAX5 se une a potenciadores **distales**, fuera de las ventanas de ±10 kb.
  - Su motivo es largo o degenerado y puntúa peor en los rankings.
  - GRNBoost2 reparte la importancia entre TF B muy correlacionados (EBF1, SPIB, POU2AF1, TCF4…). Los módulos de PAX5 quedan "mezclados" y diluyen el enriquecimiento.
  - El umbral NES es arbitrario y el tamaño es pequeño (300 células B).
- Otros regulones B sí sobrevivieron: EBF1(+) 14 genes, POU2AF1(+) 20, TCF4(+) 34, BACH2(+) 14, SPIB(+) 413.
- **Moraleja:** **ausencia de regulón ≠ ausencia de regulación.**
- **Dato extra para cerrar la discusión:** al repetir el pipeline completo con otra implementación de GRNBoost2 (notebook para casa), **PAX5(+) sí aparece como regulón**. Su ausencia en el servidor fue, en parte, cuestión de azar y umbrales.

**4. ¿Por qué tu número de regulones no coincide con los 272 del servidor?**

Prueba en un entorno tipo Colab: matriz reducida y 1 base (−500/+100 pb) → **148 regulones**, de los cuales **119** también están entre los 272 del servidor. 🧮

| Regulón | Genes (corrida reducida) | Genes (servidor) | En común |
|---|---:|---:|---:|
| MAFB(+) | **0, no aparece** | 37 | 0 |
| EBF1(+) | 39 | 14 | 6 |
| TBX21(+) | 171 | 1,014 | 155 |
| FOS(+) | 30 | 40 | 10 |
| SPI1(+) | 1,073 | 3,106 | 994 |

- **Causas:**
  - **Menos genes** (2,829 frente a 10,120): menos dianas candidatas y módulos más pequeños.
  - **Una base** en lugar de dos: los motivos que solo se detectan en la ventana de ±10 kb se pierden. En el servidor, **12 de los 13 motivos** enriquecidos de MAFB(+) vinieron de la base de ±10 kb (`results/motifs_summary.csv`), que en clase no se usa. Con solo la de −500/+100 pb, ningún módulo de MAFB llegó a NES ≥ 3. Sugiere que MAFB actúa más bien sobre regiones alejadas del promotor (inferido).
  - **El componente aleatorio de GRNBoost2** y la implementación distinta.
- **Qué implica:**
  - La **pertenencia exacta** de genes a un regulón es frágil.
  - **El patrón por población** es mucho más estable. Con sus propios regulones (celda opcional de la sección 3), EBF1(+) sigue siendo más activo en B, TBX21(+) en NK, y FOS(+) y SPI1(+) en monocitos.

    | Mediana AUC (Mono / B / NK) | Corrida reducida | Servidor |
    |---|---|---|
    | EBF1(+) | 0.054 / **0.200** / 0.061 | 0.000 / **0.049** / 0.000 |
    | TBX21(+) | 0.034 / 0.037 / **0.168** | 0.024 / 0.028 / **0.075** |
    | FOS(+) | **0.300** / 0.132 / 0.112 | **0.694** / 0.192 / 0.217 |
    | SPI1(+) | **0.111** / 0.070 / 0.067 | **0.060** / 0.044 / 0.042 |

- **Una tercera corrida** lo confirma: el notebook para casa con la matriz **demo completa** (10,120 genes), **las dos bases** y la misma semilla. Solo cambia la implementación de GRNBoost2 (`arboreto_with_multiprocessing.py` en lugar de `pyscenic grn` con dask).
  - Resultado: **255 regulones**, 212 en común con el servidor.
  - **MAFB(+) tampoco aparece.**
  - **PAX5(+) sí aparece**, aunque en el servidor no existía.
  - Correlación célula a célula del AUC en los regulones comunes: mediana 0.45.
  - Jaccard de genes frente al servidor: SPI1(+) 0.78, CEBPB(+) 0.52, EOMES(+) 0.49, TBX21(+) 0.33, SPIB(+) 0.19, EBF1(+) 0.12, FOS(+) 0.12.
- **Cuarta corrida:** el notebook para casa con el **dataset entero** (10,280 células × 10,165 genes, 2 bases).
  - Resultado: **237 regulones**, 168 en común con el servidor.
  - **MAFB(+) reaparece**, con 75 genes, y es el 4.º más específico de monocitos. Ahí también están FOS, CEBPA, NFE2 y CEBPD.
  - EBF1(+) encabeza las células B y TBX21(+) las NK.
  - POU2AF1(+) desaparece.
- **Por eso** SCENIC recomienda repetir GRNBoost2 varias veces y quedarse con los regulones recurrentes. Y por eso no hay que sacar conclusiones de un regulón que solo aparece en una corrida.

---

## Ejercicio 3 · AUCell

Mediana de AUC por población, con los regulones por defecto del notebook. 🧮

| Regulón | Genes | Mono CD14+ | B | NK | Población más activa |
|---|---:|---:|---:|---:|---|
| MAFB(+) | 37 | **0.311** | 0.113 | 0.110 | Monocitos |
| EBF1(+) | 14 | 0.000 | **0.049** | 0.000 | B |
| TBX21(+) | 1,014 | 0.024 | 0.028 | **0.075** | NK |
| FOS(+) | 40 | **0.694** | 0.192 | 0.217 | Monocitos |
| SPI1(+) | 3,106 | **0.060** | 0.044 | 0.042 | Monocitos |

Expresión del TF: media de log(CP10k), con el % de células donde se detecta entre paréntesis. 🔎

| TF | Mono CD14+ | B | NK |
|---|---:|---:|---:|
| MAFB | 0.84 (65 %) | 0.00 (0.3 %) | 0.02 (2 %) |
| EBF1 | 0.01 (1 %) | 0.38 (31 %) | 0.01 (1 %) |
| TBX21 | 0.00 (0 %) | 0.01 (2 %) | 1.04 (73 %) |
| FOS | 4.33 (100 %) | 2.12 (95 %) | 2.03 (90 %) |
| SPI1 | 1.48 (93 %) | 0.42 (37 %) | 0.04 (3 %) |

**1. ¿En qué población es más activo cada regulón? ¿Coincide con lo que sabes?**

- Sí, coincide con el conocimiento previo 📚: MAFB(+) en monocitos (mieloide), EBF1(+) en células B (linaje B) y TBX21(+) en NK (T-bet, citotóxicas).
- El TF también se expresa sobre todo en esa población 🔎. Las dos evidencias apuntan al mismo lado.
- **Matiz interesante:** el regulón se detecta **aunque el TF no se detecte** (*dropout*):
  - EBF1 no se detecta en **206 de 300** células B. Aun así, **173** de ellas tienen AUC > 0 para EBF1(+); mediana 0.047, frente a 0.055 en las B que sí expresan EBF1.
  - Lo mismo pasa con TBX21: las 80 NK sin conteos de TBX21 tienen todas AUC > 0, mediana 0.070.
  - AUCell mira las **dianas**, no solo al TF, y por eso es más robusto a la escasez de conteos.

**2. FOS se expresa en las tres poblaciones; ¿y su regulón? ¿Actividad de regulón = expresión del TF?**

- **No son lo mismo.** FOS se detecta en el 90–100 % de las células de las tres poblaciones 🔎. En cambio, FOS(+) es mucho más activo en monocitos: mediana 0.694 frente a 0.192 y 0.217 🧮. Es el **2.º regulón** con mayor rango entre poblaciones de los 272.
- **Por qué pueden diferir:**
  - La actividad depende de que se expresen las **dianas**, no solo el TF.
  - Hay regulación postraduccional del TF (fosforilación, localización, socios como JUN).
  - Hay cofactores específicos de cada tipo celular.
- **Matiz importante:** el regulón FOS(+) contiene genes de identidad monocítica (S100A8, LYZ, CTSS, MNDA, VCAN, SPI1).
  - En un dataset de tres poblaciones, la coexpresión está dominada por las diferencias *entre* tipos celulares.
  - Parte de la "actividad" de FOS(+) puede ser simplemente "ser monocito".

**3. SPI1 tiene un regulón de más de 3,000 genes. ¿Cómo afecta el tamaño a la diferencia de AUC?**

- **SPI1(+)** tiene 3,106 genes; es el regulón más grande.
- La diferencia **absoluta** de AUC entre poblaciones es mínima: 0.060 frente a 0.044 y 0.042, rango 0.018, puesto 65 de 272. Eso pasa aunque SPI1 se exprese claramente más en monocitos 🔎: 1.48 frente a 0.42 y 0.04.
- **Motivo:** con miles de genes, el regulón incluye muchos genes poco específicos. Su AUC tiende a un valor intermedio parecido en todas las células.
- **Pero la diferencia es consistente:** los cuartiles no se solapan (Mono 0.059–0.061; B 0.042–0.045; NK 0.041–0.043).
- **Para la clase:**
  - Compara AUC **dentro del mismo regulón**, nunca entre regulones de distinto tamaño.
  - El heatmap usa z-score por regulón y hace que SPI1(+) parezca tan "rojo" como MAFB(+). Hay que mirar también los números reales anotados.
  - En la puntuación de especificidad (RSS, en `results/rss_by_population.csv`), SPI1(+) en monocitos vale 0.483, frente a 0.588 de MAFB(+).

**4. FOS, JUNB y FOSB son genes de respuesta temprana. ¿Cómo cambia la interpretación de FOS(+)?**

- El regulón FOS(+) contiene varios genes de respuesta temprana: **FOSB, JUNB, JUND, DUSP1, KLF6, ZFP36** (y el propio FOS en la lista). 🧮
- Estos genes se expresan más en monocitos 🔎 (media de log CP10k, Mono / B / NK):

  | Gen | Mono | B | NK |
  |---|---:|---:|---:|
  | FOS | 4.33 | 2.12 | 2.03 |
  | DUSP1 | 3.57 | 1.84 | 1.54 |
  | FOSB | 2.65 | 1.71 | 1.51 |
  | JUNB | 2.69 | 1.94 | 1.66 |

- **Interpretación prudente:**
  - Parte de la señal FOS(+) puede reflejar el **estrés de procesar la muestra** (disociación, tiempo, temperatura) más que un programa regulador fisiológico 📚.
  - Los monocitos podrían ser más sensibles a ese estrés, o expresar estos genes de forma basal. Con estos datos **no podemos distinguirlo**: haría falta una muestra fijada al instante o un control de disociación.
  - Lo mismo aplica a JUNB(+) y FOSB(+) en la lista de específicos de monocitos.

---

## Sección 3 · `pyscenic aucell` en vivo

- Con los 272 regulones del servidor (`.gmt`), la matriz completa y `--seed 2026`, el resultado es **idéntico** al del servidor: diferencia máxima 1e-16.
- **Mensaje:** AUCell es determinista con semilla. Lo variable está en los pasos anteriores (GRNBoost2 y cisTarget).

## Sección 4 (opcional) · AUCell por dentro

- La reimplementación en numpy reproduce **exactamente** el CSV precalculado: diferencia máxima 1e-16 (`scripts/35_validate_numpy_aucell.py`).
- Curva de recuperación de MAFB(+): AUC **0.340** en un monocito frente a **0.123** en una célula B.
- También coincide exactamente (diferencia 0) con el `pyscenic aucell` que corrieron en la sección 3.

---

## Sección 5 · Plantillas resueltas (ejemplos)

### EBF1(+), el ejemplo que trae el notebook

| | |
|---|---|
| **Regulón elegido** | EBF1(+): 14 genes (CD24, MZB1, TCF4, STRBP, ARHGAP25, FAM129C, DENND4A, SMAGP…) |
| **Población con mayor actividad** | Células B: mediana 0.049; 0.000 en monocitos y NK |
| **¿El TF se expresa ahí?** | Sí, pero poco: log CP10k 0.38, detectado en el 31 % de las B (≈ 0 en las otras) |
| **Hipótesis** | EBF1 mantiene un programa de identidad B que incluye a CD24, MZB1 y TCF4. |
| **Evidencia observada 🔎** | EBF1 se expresa casi solo en células B; sus dianas (CD24, MZB1…) también se expresan en B |
| **Evidencia inferida 🧮** | Asociación EBF1–diana (GRNBoost2); motivo `tfdimers__MD00215` enriquecido (NES 3.37); actividad AUC; la etiqueta "B" |
| **Limitación** | Regulón muy pequeño (14 de 1,625 candidatos, 0.9 %); NES cerca del umbral; un solo donante |
| **Experimento** | CRISPR o knockdown de EBF1 en células B y medir CD24/MZB1; ChIP-seq o CUT&RUN de EBF1 para ver unión directa |

### MAFB(+), el caso más claro de monocitos

| | |
|---|---|
| **Regulón elegido** | MAFB(+): 37 genes (NCF2, FTH1, PSAP, DUSP1, FTL, TNFSF10, ALDH2, FGL2, CD36…); incluye al propio MAFB (autorregulación) |
| **Población con mayor actividad** | Monocitos: 0.311, frente a 0.113 (B) y 0.110 (NK); 5.º regulón con mayor rango de los 272 |
| **¿El TF se expresa ahí?** | Sí: log CP10k 0.84, detectado en el 65 % de los monocitos (≤ 2 % en B y NK) |
| **Hipótesis** | MAFB activa genes de función mieloide (NADPH oxidasa NCF2, ferritinas FTH1/FTL, receptor CD36) en monocitos. |
| **Evidencia observada 🔎** | Expresión de MAFB y de sus dianas, concentrada en monocitos |
| **Evidencia inferida 🧮** | Asociaciones; motivo `metacluster_22.36` (NES 4.67); AUC; la etiqueta "Mono_CD14" |
| **Limitación** | Varias dianas (FTL, FTH1, PSAP) son genes muy abundantes en cualquier monocito: pueden ser marcadores de tipo celular más que dianas directas |
| **Experimento** | Silenciar MAFB en monocitos primarios o en macrófagos derivados y medir NCF2/CD36; ChIP-seq de MAFB |

### FOS(+), un ejemplo de qué **no** concluir sin más

| | |
|---|---|
| **Población con mayor actividad** | Monocitos: 0.694, frente a 0.192 y 0.217 |
| **¿El TF se expresa ahí?** | Sí, pero también en las otras dos poblaciones (90–100 % de las células) |
| **Hipótesis prudente** | La diferencia puede deberse a la respuesta al estrés de la muestra y a la identidad monocítica, más que a una regulación por FOS propia de los monocitos. |
| **Limitación** | Genes de respuesta temprana en el regulón; confusión con el tipo celular |
| **Experimento** | Comparar con una muestra fijada sin disociación prolongada, o con un protocolo de disociación en frío |

**Qué buscar al corregir la plantilla:**

1. La hipótesis está formulada como **hipótesis** ("podría", "propone").
2. **No** ponen la actividad AUC ni la etiqueta de población como evidencia "observada".
3. La limitación es concreta: tamaño del regulón, NES, un donante, artefacto…
4. El experimento pone a prueba la **causalidad**: perturbación o unión directa.

---

## Sección 6 · Limitaciones (lo que deberían poder decir al final)

- **Sí:** ciertos regulones muestran actividad diferencial entre poblaciones, coherente con la biología conocida.
- **No sin más evidencia:**
  - **Causalidad.** Caso RPS4X.
  - **Ausencia de regulación.** Caso PAX5.
  - **Generalización.** Un donante, 900 células, etiquetas inferidas.
  - **Artefactos técnicos.** Caso FOS(+).
  - **Estabilidad.** GRNBoost2 es estocástico; SCENIC recomienda varias ejecuciones y quedarse con los regulones recurrentes.
