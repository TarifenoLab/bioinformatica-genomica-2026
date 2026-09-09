# Práctico 2. Preprocesamiento de datos RNA-seq

## Introducción

En el práctico anterior examinamos los archivos FASTQ sin procesar de seis muestras de células endocrinas pancreáticas —pancreatic endocrine cells, PECs— aisladas mediante FACS desde embriones de pez cebra de 27 horas post-fertilización. Tres muestras corresponden a organismos wild-type y tres a mutantes `pax6b-/-`.

Mediante SeqKit, FastQ Screen, FastQC y MultiQC evaluamos la estructura, composición, posible origen biológico y calidad de los reads. A partir de esos resultados formulamos un diagnóstico y propusimos posibles intervenciones.

En este práctico utilizaremos ese diagnóstico para diseñar y evaluar una estrategia de preprocesamiento. Trabajaremos con [Cutadapt](https://cutadapt.readthedocs.io/en/stable/), una herramienta que permite identificar y remover adaptadores, recortar bases desde los extremos y filtrar reads según distintos criterios.

El preprocesamiento no consiste en aplicar automáticamente todas las modificaciones disponibles ni en conseguir que todos los módulos de FastQC alcancen un estado `PASS`. Cada intervención debe responder a un problema observado y debe conservar la mayor cantidad posible de información biológicamente útil.

> Los parámetros de preprocesamiento no estarán resueltos en las instrucciones. Deberás seleccionarlos y justificarlos utilizando los resultados obtenidos en el Práctico 1, la estructura de las bibliotecas y la documentación de Cutadapt.

Los archivos procesados que resulten de esta actividad serán utilizados posteriormente para el alineamiento al genoma de referencia de *Danio rerio*.

---

## Objetivos de aprendizaje

Al finalizar este práctico, podrás:

* Diferenciar entre eliminación de adaptadores, recorte y filtrado de reads.
* Relacionar problemas observados mediante FastQC y MultiQC con posibles intervenciones de preprocesamiento.
* Reconocer los principales tipos de adaptadores utilizados por Cutadapt.
* Seleccionar y justificar parámetros de preprocesamiento.
* Procesar simultáneamente R1 y R2 sin perder la sincronización de los pares.
* Interpretar el reporte generado por Cutadapt.
* Evaluar una estrategia mediante un procesamiento piloto.
* Comparar cuantitativamente los datos crudos y procesados.
* Reconocer los riesgos asociados al sobreprocesamiento.
* Documentar un flujo de trabajo reproducible.

---

# 1. Punto de partida

## 1.1. Diseño experimental

Continuaremos trabajando con las seis muestras utilizadas en el Práctico 1:

| Muestra     | Genotipo   | Réplica | R1                      | R2                      |
| ----------- | ---------- | ------: | ----------------------- | ----------------------- |
| `PEC_WT_1`  | Wild-type  |       1 | `PEC_WT_1_R1.fastq.gz`  | `PEC_WT_1_R2.fastq.gz`  |
| `PEC_WT_2`  | Wild-type  |       2 | `PEC_WT_2_R1.fastq.gz`  | `PEC_WT_2_R2.fastq.gz`  |
| `PEC_WT_3`  | Wild-type  |       3 | `PEC_WT_3_R1.fastq.gz`  | `PEC_WT_3_R2.fastq.gz`  |
| `PEC_MUT_1` | `pax6b-/-` |       1 | `PEC_MUT_1_R1.fastq.gz` | `PEC_MUT_1_R2.fastq.gz` |
| `PEC_MUT_2` | `pax6b-/-` |       2 | `PEC_MUT_2_R1.fastq.gz` | `PEC_MUT_2_R2.fastq.gz` |
| `PEC_MUT_3` | `pax6b-/-` |       3 | `PEC_MUT_3_R1.fastq.gz` | `PEC_MUT_3_R2.fastq.gz` |

Cada fila representa una muestra biológica. Los archivos R1 y R2 no son réplicas independientes: corresponden a los dos extremos secuenciados de los mismos fragmentos de biblioteca.

## 1.2. Recuperación del diagnóstico

Revisa el reporte `Practico_01_MultiQC.html` y la tabla de diagnóstico que completaste al finalizar el Práctico 1.

Antes de seleccionar cualquier parámetro, resume la evidencia disponible:

| Problema potencial             | Evidencia observada | Muestras afectadas | R1, R2 o ambos | ¿Requiere intervención? | Justificación |
| ------------------------------ | ------------------- | ------------------ | -------------- | ----------------------- | ------------- |
| Presencia de adaptadores       |                     |                    |                |                         |               |
| Calidad reducida en 3′         |                     |                    |                |                         |               |
| Calidad reducida en 5′         |                     |                    |                |                         |               |
| Bases `N`                      |                     |                    |                |                         |               |
| Reads de baja calidad global   |                     |                    |                |                         |               |
| Sesgo de composición inicial   |                     |                    |                |                         |               |
| Secuencias sobrerrepresentadas |                     |                    |                |                         |               |
| Otro                           |                     |                    |                |                         |               |

### Preguntas iniciales

1. ¿Qué problemas fueron comunes a todas las muestras?
2. ¿Qué diferencias sistemáticas observaste entre R1 y R2?
3. ¿Qué resultados requieren una intervención y cuáles corresponden a características esperables de la biblioteca?
4. ¿Qué información adicional necesitas antes de seleccionar los parámetros?
5. ¿Qué consecuencias podría tener modificar los reads sin evidencia suficiente?

> Un estado `WARN` o `FAIL` de FastQC no constituye por sí solo una indicación para modificar los reads. La interpretación debe considerar el protocolo de preparación de bibliotecas, la plataforma de secuenciación y el objetivo del análisis.

---

# 2. Preparación del espacio de trabajo

## 2.1. Conexión y activación del ambiente

Conéctate a la estación de cálculo mediante SSH:

```bash
ssh usuario@152.74.15.46
```

Reemplaza `usuario` por tu nombre de usuario.

Carga Conda y activa el ambiente `RNAseq`:

```bash
source /usr/local/miniconda3/etc/profile.d/conda.sh
conda activate RNAseq
```

Dirígete al directorio del módulo:

```bash
cd ~/Bioinformatica_Omicas/RNAseq
```

Comprueba tu ubicación y verifica que los datos crudos continúan disponibles:

```bash
pwd
ls -lh data/raw/
```

## 2.2. Creación de directorios

Crea los directorios para los datos procesados y los nuevos resultados:

```bash
mkdir -p data/processed
mkdir -p results/cutadapt
mkdir -p results/fastqc_processed
mkdir -p results/multiqc_processed
```

Los archivos originales permanecerán en `data/raw/` y no deben ser modificados. Cutadapt escribirá nuevos archivos dentro de `data/processed/`.

## 2.3. Comprobación de Cutadapt

Consulta la ayuda y registra la versión instalada:

```bash
cutadapt --help
cutadapt --version
```

La versión utilizada deberá indicarse posteriormente en la sección de materiales y métodos del informe.

---

# Actividad 1. Principios del preprocesamiento

Antes de construir el comando, distingue las operaciones que puede realizar una herramienta de preprocesamiento.

| Operación                       | Pregunta que responde                                                    |
| ------------------------------- | ------------------------------------------------------------------------ |
| Eliminación de adaptadores      | ¿Existe una secuencia técnica dentro del read que debe removerse?        |
| Recorte por calidad             | ¿Existen bases de baja calidad en uno o ambos extremos?                  |
| Eliminación fija de bases       | ¿Existe una razón técnica para eliminar siempre determinadas posiciones? |
| Recorte de bases `N` terminales | ¿Existen bases ambiguas en los extremos?                                 |
| Filtrado por contenido de `N`   | ¿Cuántas bases ambiguas puede contener un read antes de ser descartado?  |
| Filtrado por longitud           | ¿Cuál es la longitud mínima útil después del recorte?                    |
| Filtrado por errores esperados  | ¿Cuál es el error acumulado tolerable dentro de un read?                 |

## 1.1. Modificación y filtrado no son equivalentes

Una opción de **modificación** cambia la secuencia o la longitud de un read. Por ejemplo, puede remover un adaptador o recortar bases de baja calidad desde un extremo.

Una opción de **filtrado** determina si el read procesado será conservado o descartado. Por ejemplo, después del recorte puede eliminar los reads que quedaron por debajo de una longitud mínima.

### Ejercicios

1. Explica la diferencia entre remover un adaptador y eliminar un número fijo de bases.
2. ¿Por qué el recorte por calidad suele aplicarse en los extremos y no en cualquier posición interna?
3. ¿Qué diferencia existe entre remover bases `N` terminales y descartar un read que contiene bases `N`?
4. ¿Por qué el filtrado por longitud debe evaluarse después de las operaciones de recorte?
5. ¿Qué información puede perderse si el umbral de calidad es demasiado exigente?
6. ¿Cómo podría afectar al alineamiento conservar reads excesivamente cortos?
7. ¿Por qué conseguir más bases Q30 no garantiza por sí solo un mejor conjunto de datos?

---

# Actividad 2. Reconocimiento de adaptadores con Cutadapt

Consulta la [guía de Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html), especialmente las secciones **Adapter types**, **Adapter search parameters**, **Quality trimming**, **Filtering reads** y **Trimming paired-end reads**.

## 2.1. Tipos de adaptadores

Cutadapt utiliza distintas opciones según la posición esperada de una secuencia:

| Tipo de secuencia                       | R1   | R2   |
| --------------------------------------- | ---- | ---- |
| Adaptador 3′                            | `-a` | `-A` |
| Adaptador 5′                            | `-g` | `-G` |
| Adaptador en cualquiera de los extremos | `-b` | `-B` |

La posición esperada puede especificarse con mayor precisión mediante adaptadores anclados u otras formas descritas en la documentación.

### Ejercicios

Explica con tus propias palabras:

1. Adaptador regular 3′.
2. Adaptador regular 5′.
3. Adaptador anclado 3′.
4. Adaptador anclado 5′.
5. Adaptador no interno.
6. Adaptador enlazado o *linked adapter*.
7. Diferencia entre las opciones `-a`, `-g` y `-b`.
8. Diferencia entre las opciones en minúscula y sus equivalentes en mayúscula.

## 2.2. Tolerancia de error y solapamiento

Cutadapt permite reconocer adaptadores aun cuando la coincidencia no sea perfecta. Investiga las opciones asociadas con:

* Tasa máxima de error.
* Número de errores permitidos.
* Longitud mínima de solapamiento.
* Inclusión o exclusión de inserciones y deleciones.

### Ejercicios

1. ¿Cómo calcula Cutadapt la tasa de error de una coincidencia?
2. ¿Por qué el número de errores permitidos depende de la longitud de la coincidencia?
3. ¿Qué riesgo existe al aceptar un solapamiento demasiado corto?
4. ¿Qué riesgo existe al exigir un solapamiento demasiado largo?
5. Para una coincidencia de 10 nucleótidos con un error, calcula la tasa de error.
6. Para una coincidencia de 8 nucleótidos con un error, calcula la tasa de error.
7. ¿Ambas coincidencias serían aceptadas utilizando la tasa máxima predeterminada? Fundamenta.

## 2.3. Wildcards y repeticiones

Cutadapt permite representar bases ambiguas mediante códigos IUPAC y escribir repeticiones de manera abreviada.

Describe la secuencia:

```text
YAGTTA{8}TTCGA
```

En tu respuesta indica:

1. Qué nucleótidos puede representar `Y`.
2. Qué significa `{8}`.
3. Cuál es la longitud total de la secuencia representada.
4. En qué situación sería útil utilizar una base degenerada al definir un adaptador o primer.

---

# Actividad 3. Diseño de la estrategia de preprocesamiento

## 3.1. Identificación de las secuencias técnicas

Utiliza los módulos **Overrepresented sequences** y **Adapter Content** de FastQC, el reporte MultiQC y la información del protocolo de preparación de bibliotecas para investigar:

1. Qué secuencias técnicas fueron detectadas.
2. Si aparecen en R1, R2 o ambos.
3. Si corresponden efectivamente a adaptadores, primers u otra secuencia.
4. En qué extremo del read se encuentran.
5. Si aparecen completas o parcialmente.
6. Qué secuencia debería proporcionarse a Cutadapt.

> No asumas que toda secuencia sobrerrepresentada es un adaptador. Puede corresponder a un transcrito biológicamente abundante, una secuencia repetitiva u otro componente de la biblioteca.

Registra las fuentes utilizadas para identificar las secuencias técnicas.

## 3.2. Selección de operaciones y parámetros

Completa la siguiente tabla antes de ejecutar Cutadapt:

| Problema observado            | Evidencia | Intervención propuesta | Opción de Cutadapt | Parámetro seleccionado | Resultado esperado | Riesgo |
| ----------------------------- | --------- | ---------------------- | ------------------ | ---------------------- | ------------------ | ------ |
| Adaptadores                   |           |                        |                    |                        |                    |        |
| Calidad 3′ de R1              |           |                        |                    |                        |                    |        |
| Calidad 3′ de R2              |           |                        |                    |                        |                    |        |
| Calidad 5′ de R1              |           |                        |                    |                        |                    |        |
| Calidad 5′ de R2              |           |                        |                    |                        |                    |        |
| Bases `N`                     |           |                        |                    |                        |                    |        |
| Longitud posterior al recorte |           |                        |                    |                        |                    |        |
| Otro                          |           |                        |                    |                        |                    |        |

### Preguntas para orientar la selección

1. ¿Qué adaptador o adaptadores utilizarás?
2. ¿Qué tipo de adaptador representa cada secuencia?
3. ¿Qué solapamiento mínimo exigirás?
4. ¿Mantendrás la tolerancia de error predeterminada o la modificarás?
5. ¿Aplicarás recorte por calidad?
6. ¿R1 y R2 requieren los mismos umbrales?
7. ¿Existe evidencia para recortar el extremo 5′?
8. ¿Eliminarás un número fijo de bases? ¿Qué evidencia lo justifica?
9. ¿Cómo tratarás las bases `N`?
10. ¿Cuál será la longitud mínima aceptada?
11. ¿Cómo afectaría esa longitud a la capacidad de alineamiento?
12. ¿Qué pérdida de reads o nucleótidos considerarías excesiva?

> No elimines bases iniciales únicamente para corregir un sesgo de composición observado por FastQC. Antes debes determinar si el patrón es compatible con el protocolo de preparación de bibliotecas y si representa realmente un problema para el análisis posterior.

## 3.3. Coherencia entre muestras

Las seis bibliotecas pertenecen al mismo experimento y fueron generadas utilizando el mismo protocolo. Por esta razón, una vez validada la estrategia piloto, deberá aplicarse el mismo conjunto de parámetros a todas las muestras.

Cualquier procesamiento diferente entre muestras debe estar respaldado por evidencia y discutirse como una posible fuente de sesgo técnico.

---

# Actividad 4. Procesamiento piloto de una muestra

Antes de procesar el conjunto completo, evalúa tu estrategia utilizando el par R1–R2 de `PEC_WT_1`.

## 4.1. Sintaxis paired-end

La estructura general de un comando paired-end es:

```bash
cutadapt OPCIONES_JUSTIFICADAS \
    -o SALIDA_R1.fastq.gz \
    -p SALIDA_R2.fastq.gz \
    ENTRADA_R1.fastq.gz \
    ENTRADA_R2.fastq.gz
```

Este esquema no es un comando ejecutable. Debes reemplazar los nombres genéricos y agregar únicamente las opciones justificadas en la Actividad 3.

En modo paired-end:

* `-o` indica el archivo de salida para R1.
* `-p` indica el archivo de salida para R2.
* Ambos archivos de entrada deben pertenecer a la misma muestra.
* R1 y R2 deben procesarse simultáneamente.
* Cutadapt comprueba la compatibilidad de los identificadores de ambos mates.
* Si un criterio de filtrado determina que un par debe eliminarse, Cutadapt mantiene sincronizadas las salidas.

Investiga el funcionamiento de `--pair-filter` y determina qué comportamiento utilizarás.

## 4.2. Construcción del comando

Construye el comando para `PEC_WT_1` utilizando:

* Como entradas: los FASTQ ubicados en `data/raw/`.
* Como salidas: nuevos FASTQ dentro de `data/processed/`.
* Un nombre que identifique claramente la muestra, el mate y que los reads fueron procesados.
* Las opciones y valores definidos en la Actividad 3.
* Un reporte de Cutadapt guardado dentro de `results/cutadapt/`.

Antes de ejecutarlo, comprueba:

* Que R1 y R2 pertenecen a `PEC_WT_1`.
* Que las salidas no se escribirán dentro de `data/raw/`.
* Que los archivos de salida R1 y R2 tienen nombres distintos.
* Que las opciones dirigidas a R1 y R2 utilizan las formas correctas.
* Que cada parámetro posee una justificación.

> Guarda el comando exacto utilizado. El análisis debe poder reproducirse posteriormente sin depender del historial de la terminal.

## 4.3. Interpretación del reporte de Cutadapt

Después de la ejecución, examina el reporte y registra:

| Métrica                             | R1 | R2 | Interpretación |
| ----------------------------------- | -: | -: | -------------- |
| Pares procesados                    |    |    |                |
| Reads con adaptadores               |    |    |                |
| Bases removidas por adaptadores     |    |    |                |
| Bases removidas por calidad         |    |    |                |
| Pares descartados por longitud      |    |    |                |
| Pares descartados por otro criterio |    |    |                |
| Pares conservados                   |    |    |                |
| Porcentaje de pares conservados     |    |    |                |

Dependiendo de las opciones utilizadas, algunas categorías podrían no aparecer en el reporte.

### Ejercicios

1. ¿Qué proporción de reads contenía adaptadores?
2. ¿La frecuencia fue similar entre R1 y R2?
3. ¿Cuántas bases fueron removidas por calidad?
4. ¿Qué criterio produjo la mayor pérdida de pares?
5. ¿La pérdida de información fue mayor o menor que la esperada?
6. ¿La distribución de las longitudes de coincidencia apoya que la secuencia removida corresponde a un adaptador real?
7. ¿Detectas evidencia de coincidencias inespecíficas?
8. ¿Qué resultado del reporte te haría reconsiderar un parámetro?

---

# Actividad 5. Evaluación del procesamiento piloto

## 5.1. Caracterización con SeqKit

Utiliza `seqkit stats -a` para caracterizar los archivos procesados de `PEC_WT_1` y guarda el resultado dentro de `results/seqkit/`.

Compara:

* Número de reads.
* Número total de nucleótidos.
* Longitud mínima.
* Longitud promedio.
* Longitud máxima.
* Porcentaje de GC.
* Porcentaje de bases Q20.
* Porcentaje de bases Q30.

## 5.2. Nuevo control de calidad

Ejecuta FastQC sobre R1 y R2 procesados de `PEC_WT_1` y guarda los reportes dentro de:

```text
results/fastqc_processed/
```

Utiliza la estructura general:

```bash
fastqc [OPCIONES] ARCHIVO_R1_PROCESADO ARCHIVO_R2_PROCESADO
```

Selecciona las opciones necesarias revisando `fastqc --help` y utilizando como referencia el comando empleado en el Práctico 1.

## 5.3. Comparación antes y después

Completa la tabla:

| Métrica o módulo             | Datos crudos | Datos procesados | Cambio observado | ¿Representa una mejora? |
| ---------------------------- | ------------ | ---------------- | ---------------- | ----------------------- |
| Número de reads              |              |                  |                  |                         |
| Nucleótidos totales          |              |                  |                  |                         |
| Longitud de los reads        |              |                  |                  |                         |
| %Q20                         |              |                  |                  |                         |
| %Q30                         |              |                  |                  |                         |
| Per base sequence quality    |              |                  |                  |                         |
| Per base sequence content    |              |                  |                  |                         |
| Per base N content           |              |                  |                  |                         |
| Sequence Length Distribution |              |                  |                  |                         |
| Overrepresented sequences    |              |                  |                  |                         |
| Adapter Content              |              |                  |                  |                         |

### Evaluación del piloto

1. ¿Disminuyó la señal correspondiente a adaptadores?
2. ¿Mejoró la calidad de los extremos?
3. ¿Cambió la longitud promedio de los reads?
4. ¿Qué proporción de reads y nucleótidos se perdió?
5. ¿El procesamiento eliminó información que parecía de buena calidad?
6. ¿Apareció algún problema nuevo en FastQC?
7. ¿Existe algún módulo que no cambió? ¿Era esperable que cambiara?
8. ¿Mantendrías todos los parámetros utilizados?
9. ¿Modificarías algún parámetro antes de procesar las demás muestras?
10. ¿La estrategia produce un equilibrio razonable entre calidad y retención de información?

> El objetivo no es maximizar el número de módulos `PASS`. Una estrategia adecuada debe resolver los problemas relevantes sin provocar una pérdida desproporcionada de reads, longitud o complejidad.

---

# Actividad 6. Definición de la estrategia final

A partir del procesamiento piloto, registra la estrategia definitiva:

| Componente                | Decisión final | Evidencia | Parámetro definitivo |
| ------------------------- | -------------- | --------- | -------------------- |
| Adaptador de R1           |                |           |                      |
| Adaptador de R2           |                |           |                      |
| Tolerancia de error       |                |           |                      |
| Solapamiento mínimo       |                |           |                      |
| Calidad 5′ de R1          |                |           |                      |
| Calidad 3′ de R1          |                |           |                      |
| Calidad 5′ de R2          |                |           |                      |
| Calidad 3′ de R2          |                |           |                      |
| Tratamiento de bases `N`  |                |           |                      |
| Longitud mínima           |                |           |                      |
| Comportamiento paired-end |                |           |                      |
| Otro                      |                |           |                      |

Discute la propuesta con el curso antes de continuar. Una vez acordada y validada la estrategia, utiliza los mismos parámetros para las seis muestras.

---

# Actividad 7. Procesamiento de las seis muestras

Procesa cada par de archivos de manera independiente:

| Orden | Muestra     | Entrada R1              | Entrada R2              |
| ----: | ----------- | ----------------------- | ----------------------- |
|     1 | `PEC_WT_1`  | `PEC_WT_1_R1.fastq.gz`  | `PEC_WT_1_R2.fastq.gz`  |
|     2 | `PEC_WT_2`  | `PEC_WT_2_R1.fastq.gz`  | `PEC_WT_2_R2.fastq.gz`  |
|     3 | `PEC_WT_3`  | `PEC_WT_3_R1.fastq.gz`  | `PEC_WT_3_R2.fastq.gz`  |
|     4 | `PEC_MUT_1` | `PEC_MUT_1_R1.fastq.gz` | `PEC_MUT_1_R2.fastq.gz` |
|     5 | `PEC_MUT_2` | `PEC_MUT_2_R1.fastq.gz` | `PEC_MUT_2_R2.fastq.gz` |
|     6 | `PEC_MUT_3` | `PEC_MUT_3_R1.fastq.gz` | `PEC_MUT_3_R2.fastq.gz` |

Para cada muestra:

1. Comprueba los nombres de entrada.
2. Utiliza la estrategia definitiva.
3. Especifica una salida para R1 y otra para R2.
4. Guarda un reporte independiente de Cutadapt.
5. Comprueba que la ejecución finalizó sin errores.
6. Registra el comando utilizado.

No sobrescribas los archivos crudos ni reutilices el mismo nombre de salida para muestras diferentes.

## 7.1. Registro de procesamiento

Completa la siguiente tabla:

| Muestra     | Pares iniciales | Pares conservados | Pares descartados | Retención (%) | Reads R1 finales | Reads R2 finales | Observaciones |
| ----------- | --------------: | ----------------: | ----------------: | ------------: | ---------------: | ---------------: | ------------- |
| `PEC_WT_1`  |                 |                   |                   |               |                  |                  |               |
| `PEC_WT_2`  |                 |                   |                   |               |                  |                  |               |
| `PEC_WT_3`  |                 |                   |                   |               |                  |                  |               |
| `PEC_MUT_1` |                 |                   |                   |               |                  |                  |               |
| `PEC_MUT_2` |                 |                   |                   |               |                  |                  |               |
| `PEC_MUT_3` |                 |                   |                   |               |                  |                  |               |

### Ejercicios

1. ¿Todas las muestras conservaron una proporción similar de pares?
2. ¿Alguna muestra perdió considerablemente más reads?
3. ¿La pérdida se relaciona con adaptadores, calidad, longitud u otro criterio?
4. ¿Las diferencias observadas ya eran evidentes en el control de calidad inicial?
5. ¿Existe alguna muestra que deba examinarse nuevamente antes del alineamiento?

---

# Actividad 8. Verificación de los archivos procesados

## 8.1. Número y nombres de archivos

Cuenta los FASTQ procesados:

```bash
find data/processed -name "*.fastq.gz" | wc -l
```

Cuenta R1 y R2 por separado, adaptando los patrones de búsqueda a la nomenclatura que utilizaste:

```bash
find data/processed -name "*_R1*.fastq.gz" | wc -l
find data/processed -name "*_R2*.fastq.gz" | wc -l
```

### Ejercicios

1. ¿Cuántos archivos procesados esperabas obtener?
2. ¿Cada muestra posee una salida R1 y una salida R2?
3. ¿Los nombres permiten distinguir inequívocamente los datos crudos de los procesados?
4. ¿Detectas archivos incompletos, vacíos o con nombres inconsistentes?

## 8.2. Sincronización de R1 y R2

Ejecuta `seqkit stats -a` sobre todos los archivos procesados:

```bash
seqkit stats -a data/processed/*.fastq.gz \
    > results/seqkit/seqkit_stats_processed.txt
```

Para cada muestra comprueba que R1 y R2 contienen el mismo número de reads.

### Ejercicios

1. ¿Todos los pares permanecieron sincronizados?
2. ¿Por qué el número de reads de R1 debe coincidir con el de R2?
3. ¿Qué error en el procesamiento podría producir archivos desincronizados?
4. ¿Por qué no se recomienda procesar R1 y R2 por separado cuando se aplican criterios de filtrado?

---

# Actividad 9. Control de calidad posterior al preprocesamiento

## 9.1. FastQC

Ejecuta FastQC sobre los doce archivos procesados y guarda los resultados en:

```text
results/fastqc_processed/
```

Construye el comando utilizando como referencia el Práctico 1 y la ayuda de FastQC.

## 9.2. MultiQC

Genera un reporte MultiQC que integre:

* Los reportes FastQC de los datos crudos.
* Los reportes FastQC de los datos procesados.
* Los reportes de Cutadapt.

Guarda el nuevo reporte dentro de:

```text
results/multiqc_processed/
```

Utiliza el nombre:

```text
Practico_02_MultiQC.html
```

Revisa la ayuda de MultiQC y construye el comando necesario:

```bash
multiqc --help
```

> Los nombres de los archivos procesados deben permitir distinguirlos de los archivos crudos dentro del reporte integrado.

---

# Actividad 10. Comparación global antes y después

Utiliza SeqKit, los reportes de Cutadapt y MultiQC para completar la tabla:

| Muestra     | Pares crudos | Pares procesados | Retención (%) | Longitud cruda | Longitud procesada | %Q30 crudo | %Q30 procesado | Adaptadores antes/después | Observaciones |
| ----------- | -----------: | ---------------: | ------------: | -------------- | ------------------ | ---------: | -------------: | ------------------------- | ------------- |
| `PEC_WT_1`  |              |                  |               |                |                    |            |                |                           |               |
| `PEC_WT_2`  |              |                  |               |                |                    |            |                |                           |               |
| `PEC_WT_3`  |              |                  |               |                |                    |            |                |                           |               |
| `PEC_MUT_1` |              |                  |               |                |                    |            |                |                           |               |
| `PEC_MUT_2` |              |                  |               |                |                    |            |                |                           |               |
| `PEC_MUT_3` |              |                  |               |                |                    |            |                |                           |               |

## Preguntas de análisis

1. ¿Qué problemas fueron corregidos o reducidos?
2. ¿Qué módulos de FastQC cambiaron de manera importante?
3. ¿Qué módulos permanecieron sin cambios?
4. ¿Los cambios fueron similares en R1 y R2?
5. ¿Las seis muestras respondieron de manera comparable al procesamiento?
6. ¿Alguna muestra se convirtió en un outlier después del procesamiento?
7. ¿Qué porcentaje de pares y nucleótidos se perdió?
8. ¿La pérdida fue proporcional al problema observado inicialmente?
9. ¿El aumento de Q20 o Q30 se debe en parte a que se eliminaron bases?
10. ¿Cómo cambió la distribución de longitudes?
11. ¿La reducción de adaptadores justifica la pérdida de información?
12. ¿Existe evidencia de sobreprocesamiento?
13. ¿Qué problemas no pueden resolverse mediante Cutadapt?
14. ¿Consideras que los archivos procesados son apropiados para el alineamiento? Fundamenta.

> Una comparación antes/después debe considerar simultáneamente la calidad y la cantidad de información retenida. Reportar solamente que aumentó el porcentaje de bases Q30 es insuficiente.

---

# Actividad 11. Comparación de herramientas

Además de Cutadapt, existen otras herramientas para el preprocesamiento de datos de secuenciación. Revisa la documentación oficial de:

* [fastp](https://github.com/OpenGene/fastp)
* [Trimmomatic](https://github.com/usadellab/Trimmomatic)
* [BBDuk](https://bbmap.org/tools/bbduk)

Completa la tabla:

| Característica                      | Cutadapt | fastp | Trimmomatic | BBDuk |
| ----------------------------------- | -------- | ----- | ----------- | ----- |
| Eliminación de adaptadores          |          |       |             |       |
| Detección automática de adaptadores |          |       |             |       |
| Recorte por calidad                 |          |       |             |       |
| Filtrado por longitud               |          |       |             |       |
| Tratamiento de bases `N`            |          |       |             |       |
| Manejo paired-end                   |          |       |             |       |
| Filtrado por k-mers                 |          |       |             |       |
| Reporte HTML o JSON                 |          |       |             |       |
| Integración con MultiQC             |          |       |             |       |
| Principal fortaleza                 |          |       |             |       |
| Principal limitación                |          |       |             |       |

### Preguntas

1. ¿Qué herramienta ofrece el mayor control explícito sobre la búsqueda de adaptadores?
2. ¿Cuál automatiza una mayor proporción del análisis?
3. ¿Qué ventajas y riesgos tiene la detección automática de adaptadores?
4. ¿Qué herramienta seleccionarías para procesar este experimento y por qué?
5. ¿Utilizar herramientas diferentes con parámetros predeterminados debería producir resultados equivalentes? Fundamenta.

---

# Actividad 12. Documentación del análisis

Redacta un párrafo para la sección de materiales y métodos que incluya:

* Herramienta y versión.
* Tipo de datos procesados.
* Forma de procesamiento paired-end.
* Adaptadores evaluados y removidos.
* Criterios de recorte por calidad.
* Criterios de filtrado.
* Longitud mínima conservada.
* Estrategia utilizada para evaluar el resultado.

No escribas solamente una lista de comandos. La metodología debe permitir comprender qué se hizo, por qué se hizo y con qué parámetros.

Redacta además una breve sección de resultados que explique:

* Cuánto material fue conservado.
* Qué problemas se redujeron.
* Qué diferencias persistieron.
* Qué limitaciones posee la estrategia.
* Por qué los archivos finales fueron seleccionados —o no— para el alineamiento.

---

# Transferencia de resultados

Los comandos `scp` deben ejecutarse desde una terminal de tu computador personal, no desde la sesión SSH abierta en la estación.

Sal de la estación:

```bash
exit
```

Transfiere los reportes necesarios para su revisión:

```bash
scp -r usuario@152.74.15.46:~/Bioinformatica_Omicas/RNAseq/results/cutadapt .
scp -r usuario@152.74.15.46:~/Bioinformatica_Omicas/RNAseq/results/fastqc_processed .
scp -r usuario@152.74.15.46:~/Bioinformatica_Omicas/RNAseq/results/multiqc_processed .
```

Reemplaza `usuario` por tu nombre de usuario.

No es necesario transferir los FASTQ completos a tu computador personal para abrir los reportes HTML.

---

# Presentación y discusión de resultados

Los resultados serán discutidos en clase utilizando el mismo PPT colaborativo del módulo de RNA-seq.

Cada presentación deberá incluir:

* Diagnóstico que motivó el procesamiento.
* Intervenciones seleccionadas.
* Parámetros y justificación.
* Resultado del procesamiento piloto.
* Porcentaje de pares y nucleótidos conservados.
* Comparación antes/después.
* Problemas corregidos y problemas persistentes.
* Riesgos o limitaciones de la estrategia.
* Decisión respecto del alineamiento.

Las actividades desarrolladas desde el Práctico 1 formarán parte de un único informe asociado al módulo de RNA-seq. No se debe entregar un informe independiente para este práctico.

---

# Resultados mínimos esperados

Al finalizar debes contar con una estructura equivalente a:

```text
RNAseq/
├── data/
│   ├── raw/
│   │   ├── PEC_WT_1_R1.fastq.gz
│   │   ├── PEC_WT_1_R2.fastq.gz
│   │   ├── ...
│   │   ├── PEC_MUT_3_R1.fastq.gz
│   │   └── PEC_MUT_3_R2.fastq.gz
│   └── processed/
│       ├── PEC_WT_1_R1.processed.fastq.gz
│       ├── PEC_WT_1_R2.processed.fastq.gz
│       ├── ...
│       ├── PEC_MUT_3_R1.processed.fastq.gz
│       └── PEC_MUT_3_R2.processed.fastq.gz
└── results/
    ├── seqkit/
    │   ├── seqkit_stats_raw.txt
    │   └── seqkit_stats_processed.txt
    ├── fastqc/
    │   └── ...
    ├── multiqc/
    │   └── Practico_01_MultiQC.html
    ├── cutadapt/
    │   ├── PEC_WT_1.cutadapt.log
    │   ├── ...
    │   └── PEC_MUT_3.cutadapt.log
    ├── fastqc_processed/
    │   └── ...
    └── multiqc_processed/
        └── Practico_02_MultiQC.html
```

Los nombres exactos dependerán de la nomenclatura definida durante el práctico, pero deben ser consistentes e inequívocos.

Antes de finalizar, comprueba que:

* Utilizaste los resultados del Práctico 1 para seleccionar los parámetros.
* Justificaste cada intervención.
* Procesaste R1 y R2 simultáneamente.
* Validaste primero la estrategia con una muestra piloto.
* Aplicaste la estrategia definitiva de manera coherente a las seis muestras.
* No modificaste ni sobrescribiste los datos crudos.
* Conservaste una salida R1 y una salida R2 por muestra.
* Verificaste la sincronización de los pares.
* Registraste la versión de Cutadapt.
* Guardaste los comandos y reportes.
* Repetiste FastQC sobre los archivos procesados.
* Integraste los resultados mediante MultiQC.
* Comparaste calidad y cantidad de información antes y después.
* Fundamentaste si los datos pueden continuar al alineamiento.
