# Pauta de informe – Módulo de RNA-seq

## Objetivo

Elaborar un informe científico que integre los resultados obtenidos durante los prácticos de RNA-seq, desde el control de calidad de las secuencias hasta el análisis de expresión diferencial y su interpretación biológica.

El informe debe presentar de manera sintética **el racional de los análisis, las principales decisiones metodológicas, los resultados obtenidos y su interpretación**. No se espera responder individualmente todas las preguntas de los prácticos ni describir cada comando ejecutado.

**Extensión máxima:** 10 páginas, excluyendo portada, referencias y anexos. Los scripts, tablas extensas y figuras complementarias deberán incluirse como anexos.

## 1. Introducción y objetivos
*Extensión sugerida: 1 página.*

Presentar brevemente el contexto biológico del experimento, el diseño experimental y la pregunta que se busca responder mediante RNA-seq.

Incluir las características de las muestras, condiciones comparadas y objetivos del análisis.

## 2. Estrategia de análisis
*Extensión sugerida: 1 página.*

Presentar un esquema del flujo bioinformático desarrollado, indicando las herramientas utilizadas en cada etapa y su propósito.

Describir brevemente los criterios metodológicos más relevantes, incluyendo los parámetros seleccionados, el genoma de referencia y las estrategias de cuantificación, normalización y análisis estadístico.

Los comandos completos deben incorporarse en anexos.

## 3. Resultados y discusión
*Extensión sugerida: 5–6 páginas.*

Organizar los resultados en las siguientes secciones. Se recomienda utilizar figuras y tablas que permitan comparar las seis muestras y destacar los hallazgos principales.

### 3.1. Control de calidad y preprocesamiento
- Resumir la calidad inicial de los datos, identificando los principales problemas detectados.
- Explicar y justificar la estrategia de preprocesamiento seleccionada.
- Comparar los resultados antes y después del procesamiento, incluyendo calidad, retención de reads y cambios relevantes.
- Discutir si el procesamiento fue adecuado y qué problemas persistieron.

### 3.2. Alineamiento al genoma de referencia
- Describir brevemente la estrategia de alineamiento y su racional.
- Presentar las principales estadísticas de alineamiento de las seis muestras.
- Evaluar la consistencia entre réplicas, posibles muestras atípicas y la relación con los resultados del control de calidad.

### 3.3. Cuantificación y expresión diferencial
- Explicar brevemente cómo se obtuvieron los conteos génicos, los criterios de asignación de fragmentos y el fundamento de la normalización.
- Resumir los resultados de cuantificación y la proporción de fragmentos asignados a genes.
- Presentar el número de genes diferencialmente expresados, distinguiendo aquellos con mayor expresión en WT y mutantes.
- Incorporar e interpretar el **PCA y el MA plot**, considerando el comportamiento de las réplicas y las diferencias entre condiciones.

### 3.4. Interpretación biológica y enriquecimiento funcional
- Seleccionar **tres genes con mayor expresión en mutantes y tres con mayor expresión en WT**, utilizando preferentemente los criterios `padj < 0,05` y `|log2FoldChange| ≥ 1`.
- Presentar una tabla con sus identificadores, nombres, `baseMean`, `log2FoldChange`, `padj` y funciones biológicas.
- Investigar y discutir su posible relación con el contexto biológico del experimento, utilizando literatura científica.
- Realizar un **análisis de enriquecimiento funcional con DAVID o g:Profiler**, utilizando los genes diferencialmente expresados.
- Incorporar las tablas de enriquecimiento obtenidas, indicando las categorías funcionales, vías o procesos biológicos identificados y su significancia estadística. Las tablas completas pueden quedar en anexos, destacando los resultados principales en el informe.
- Discutir qué procesos biológicos podrían estar asociados a las diferencias transcriptómicas observadas y cómo se relacionan con la función de *pax6b*.

Se debe indicar qué conjunto de genes se utilizó para el enriquecimiento, cuál fue el universo de referencia y qué criterio estadístico se aplicó. Evitar interpretar el enriquecimiento como demostración de activación o inhibición de una vía sin evidencia adicional.

## 4. Conclusiones
*Extensión sugerida: 0,5–1 página.*

Integrar los principales hallazgos metodológicos y biológicos del análisis.

Discutir brevemente las limitaciones de los resultados y plantear las conclusiones más relevantes respecto de las diferencias transcriptómicas entre WT y mutantes *pax6b-/-*.

## 5. Referencias

Incluir las referencias científicas utilizadas para contextualizar el experimento, justificar aspectos metodológicos relevantes e interpretar los resultados biológicos.

## 6. Anexos

Los anexos no se consideran dentro del límite de páginas.

Deberán incluir:

1. **Scripts y comandos:** archivos utilizados en cada etapa del análisis (Bash y R), organizados e identificados para permitir su reproducción.
2. **Tablas de resultados:** matrices de conteos crudos y normalizados, resultados completos de expresión diferencial y tablas de enriquecimiento funcional de DAVID o g:Profiler.
3. **Figuras complementarias:** resultados adicionales de FastQC, MultiQC, Cutadapt, HISAT2 u otros análisis que respalden las interpretaciones del informe.

Los scripts deben contener las instrucciones y parámetros realmente utilizados. No es necesario incluir archivos FASTQ, SAM o BAM.

---

## Consideraciones generales

- **Priorizar la interpretación sobre la descripción:** no basta con mostrar una figura; se debe explicar qué representa y qué se concluye de ella.
- **Integrar las etapas:** relacionar los resultados del control de calidad y preprocesamiento con el alineamiento, la cuantificación y la expresión diferencial.
- **Justificar las decisiones:** explicar brevemente por qué se utilizaron determinados parámetros y criterios.
- **Sintetizar la información:** seleccionar las figuras y tablas más relevantes y evitar repetir resultados.
- **Mantener la reproducibilidad:** todos los análisis deben poder reconstruirse a partir de los scripts y archivos entregados.

**El informe será evaluado principalmente por la capacidad de seleccionar, interpretar e integrar los resultados, así como por la claridad de la argumentación científica y la reproducibilidad de los análisis.**

## Rúbrica de evaluación del informe de RNA-seq

El informe será evaluado considerando la comprensión del flujo bioinformático, la justificación de las decisiones metodológicas, la interpretación de los resultados y su integración biológica.

| Criterio de evaluación | Puntaje máximo |
|---|---:|
| **1. Introducción y objetivos** | **5** |
| Contextualiza el problema biológico, describe correctamente el diseño experimental y plantea objetivos claros. | 5 |
| **2. Estrategia de análisis** | **15** |
| Presenta un flujo de análisis coherente, identificando las herramientas utilizadas y el propósito de cada etapa. | 5 |
| Explica el racional de los análisis y justifica los principales parámetros y decisiones metodológicas. | 10 |
| **3. Resultados y discusión** | **60** |
| **3.1. Control de calidad y preprocesamiento:** interpreta las métricas de calidad, identifica problemas y justifica la estrategia de procesamiento. Evalúa críticamente los cambios obtenidos. | 15 |
| **3.2. Alineamiento:** interpreta las estadísticas de alineamiento, compara las muestras y relaciona los resultados con el control de calidad previo. | 10 |
| **3.3. Cuantificación y expresión diferencial:** comprende los criterios de cuantificación, normalización y análisis estadístico. Interpreta correctamente los genes diferencialmente expresados, el PCA y el MA plot. | 20 |
| **3.4. Interpretación biológica:** analiza genes seleccionados e interpreta las tablas de enriquecimiento funcional de DAVID o g:Profiler. Relaciona los resultados con el contexto biológico utilizando literatura científica. | 15 |
| **4. Conclusiones** | **10** |
| Integra los principales resultados metodológicos y biológicos, presenta conclusiones sustentadas en evidencia e identifica limitaciones del estudio. | 10 |
| **5. Referencias y reproducibilidad** | **5** |
| Utiliza referencias pertinentes y entrega scripts, comandos y archivos complementarios completos y organizados que permitan reproducir los análisis. | 5 |
| **6. Presentación científica** | **5** |
| Presenta una redacción clara y precisa, figuras y tablas correctamente identificadas, organización coherente y adecuada capacidad de síntesis. | 5 |
| **TOTAL** | **100 puntos** |

### Niveles de desempeño

Para cada criterio, el puntaje se asignará de acuerdo con los siguientes niveles:

| Nivel | Porcentaje del puntaje | Descripción |
|---|---:|---|
| **Excelente** | 90–100% | Demuestra comprensión sólida, interpreta correctamente los resultados, justifica sus decisiones y establece relaciones críticas entre los análisis. |
| **Bueno** | 75–89% | Presenta resultados e interpretaciones correctas, con algunas debilidades menores en la justificación o integración. |
| **Suficiente** | 60–74% | Presenta los resultados principales, pero predominan las descripciones sobre el análisis crítico o existen explicaciones incompletas. |
| **Insuficiente** | 0–59% | Presenta omisiones importantes, errores de interpretación o evidencia limitada de comprensión de los análisis. |

### Consideraciones para la evaluación

- Se privilegiará la **comprensión y el razonamiento bioinformático** por sobre la simple ejecución de herramientas.
- Los resultados deben estar acompañados de una interpretación y una justificación del análisis realizado.
- Se valorará especialmente la capacidad de relacionar los resultados entre las diferentes etapas del flujo de RNA-seq.
- La ausencia de resultados biológicos esperados, como una separación clara en el PCA o categorías funcionales significativamente enriquecidas, **no será penalizada si los resultados son correctamente analizados y discutidos**.
- La inclusión de numerosas figuras, tablas o resultados sin interpretación no constituye evidencia suficiente de aprendizaje.
- La reproducibilidad del análisis debe quedar respaldada mediante los scripts y comandos entregados en los anexos.

**El criterio central será demostrar que el estudiante comprende qué hizo, por qué lo hizo, qué información obtuvo y qué conclusiones puede (o no puede) extraer de sus resultados.**
