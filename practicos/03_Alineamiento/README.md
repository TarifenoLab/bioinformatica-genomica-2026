# Práctico 3. Alineamiento de datos RNA-seq

## Introducción

En los prácticos anteriores evaluamos la calidad de seis bibliotecas de RNA-seq paired-end correspondientes a células endocrinas pancreáticas aisladas mediante FACS desde embriones de pez cebra de 27 horas post-fertilización.

Trabajamos con tres réplicas biológicas *wild-type* y tres réplicas mutantes `pax6b-/-`.

Después del control de calidad inicial, diseñamos una estrategia de preprocesamiento, la evaluamos utilizando una muestra piloto y posteriormente aplicamos la estrategia definitiva a las seis bibliotecas.

En este práctico utilizaremos los archivos FASTQ procesados para realizar el **alineamiento de las secuencias contra el genoma de referencia de *Danio rerio***.

Para ello utilizaremos **HISAT2**, un alineador capaz de trabajar con datos de RNA-seq y reconocer alineamientos que atraviesan uniones de *splicing*.

Posteriormente utilizaremos **SAMtools** para manipular los archivos generados durante el alineamiento y obtener archivos BAM ordenados e indexados que puedan utilizarse en los análisis posteriores.

El objetivo no es solamente obtener archivos alineados. Durante este práctico deberás comprender:

- Qué significa alinear un read contra un genoma de referencia.
- Qué características particulares posee el alineamiento de RNA-seq.
- Qué información contienen los índices utilizados por HISAT2.
- Cómo incorporar información conocida sobre uniones de *splicing*.
- Cómo interpretar las estadísticas de alineamiento.
- Qué información contiene un archivo SAM/BAM.
- Cómo transformar los alineamientos en archivos adecuados para análisis posteriores.
- Cómo evaluar si las seis bibliotecas presentan resultados de alineamiento consistentes.

---

## Objetivos de aprendizaje

Al finalizar este práctico podrás:

- Diferenciar entre genoma de referencia, anotación génica e índice de alineamiento.
- Explicar por qué el alineamiento de RNA-seq requiere considerar el *splicing*.
- Reconocer los archivos necesarios para realizar un alineamiento con HISAT2.
- Consultar la documentación de HISAT2 y seleccionar los parámetros apropiados.
- Construir un comando de alineamiento para datos paired-end.
- Interpretar el resumen de alineamiento generado por HISAT2.
- Reconocer la estructura y los principales campos de un archivo SAM.
- Convertir archivos SAM a BAM utilizando SAMtools.
- Ordenar e indexar archivos BAM.
- Evaluar estadísticas básicas de los alineamientos.
- Comparar el desempeño del alineamiento entre muestras.
- Identificar posibles muestras problemáticas antes de continuar con la cuantificación.

---

# 1. Punto de partida

## 1.1. Diseño experimental

Continuaremos trabajando con las seis muestras utilizadas en los prácticos anteriores:

| Muestra | Genotipo | Réplica |
|---|---|---:|
| `PEC_WT_1` | Wild-type | 1 |
| `PEC_WT_2` | Wild-type | 2 |
| `PEC_WT_3` | Wild-type | 3 |
| `PEC_MUT_1` | `pax6b-/-` | 1 |
| `PEC_MUT_2` | `pax6b-/-` | 2 |
| `PEC_MUT_3` | `pax6b-/-` | 3 |

Cada muestra posee dos archivos FASTQ porque las bibliotecas fueron secuenciadas en modalidad **paired-end**.

Los archivos R1 y R2 corresponden a los dos extremos secuenciados de los mismos fragmentos de biblioteca y **no representan réplicas independientes**.

---

## 1.2. Datos de entrada

Utilizaremos los archivos obtenidos al finalizar el Práctico 2.

Estos archivos deberían encontrarse dentro de:

```bash
~/Bioinformatica_Omicas/RNAseq/data/processed/
```

Comprueba su contenido:

```bash
cd ~/Bioinformatica_Omicas/RNAseq

ls -lh data/processed/
```

Deberías encontrar dos archivos por muestra, correspondientes a R1 y R2.

Antes de continuar verifica:

- Que existen seis muestras.
- Que cada muestra posee R1 y R2.
- Que los archivos corresponden a los FASTQ **procesados**.
- Que R1 y R2 poseen el mismo número de reads.
- Que los nombres permiten identificar inequívocamente la muestra y el mate.

### Preguntas iniciales

1. ¿Por qué utilizaremos los FASTQ procesados y no los archivos crudos?
2. ¿Qué consecuencias tendría intercambiar accidentalmente R1 y R2 de muestras diferentes?
3. ¿Por qué es importante que R1 y R2 permanezcan sincronizados?
4. ¿Qué información obtenida en los prácticos anteriores podría ayudarnos posteriormente a interpretar una tasa de alineamiento inesperadamente baja?

---

# 2. Genoma de referencia y anotación

En este práctico alinearemos las secuencias contra el genoma de referencia de:

**Danio rerio – GRCz11**

Los archivos necesarios ya se encuentran disponibles en la estación de cálculo.

No deberás descargar ni construir nuevamente estos archivos.

Las anotaciones se encuentran en:

```text
/home/Materiales_Curso/Anotaciones/Danio_rerio_114/
```

Explora el directorio:

```bash
ls -lh /home/Materiales_Curso/Anotaciones/Danio_rerio_114/
```

Identifica los archivos presentes antes de continuar.

---

## Actividad 1. Reconocimiento de los archivos de referencia

Dentro del directorio busca los archivos correspondientes a:

- Genoma de referencia en formato FASTA.
- Anotación génica en formato GTF.
- Índices de HISAT2.
- Archivo de sitios conocidos de *splicing*.

Para cada uno registra:

| Archivo | Formato | ¿Qué información contiene? | ¿Para qué se utiliza? |
|---|---|---|---|
| Genoma | | | |
| Anotación | | | |
| Índices HISAT2 | | | |
| Sitios de *splicing* | | | |

### Preguntas

1. ¿Cuál es la diferencia entre el archivo FASTA del genoma y el archivo GTF?
2. ¿Qué información contiene un archivo GTF?
3. ¿Qué representan las coordenadas presentes en un GTF?
4. ¿Qué tipos de elementos pueden encontrarse anotados?
5. ¿Por qué el archivo GTF no reemplaza al FASTA del genoma?
6. ¿Qué relación existe entre ambos archivos?

---

# 3. ¿Por qué el alineamiento de RNA-seq es especial?

En una biblioteca de RNA-seq, los reads derivan principalmente de moléculas de RNA maduras.

En organismos eucariontes, muchos transcritos han sufrido *splicing*. Como consecuencia, un read puede contener secuencias provenientes de dos exones que se encuentran separados por un intrón en el genoma.

Conceptualmente podemos tener:

```text
RNA:

EXÓN 1 | EXÓN 2
--------READ--------

Genoma:

EXÓN 1 -------- INTRÓN -------- EXÓN 2
```

Por esta razón, el alineamiento de RNA-seq requiere herramientas capaces de identificar **alineamientos divididos o spliced alignments**.

HISAT2 puede identificar este tipo de alineamientos y además utilizar información conocida sobre sitios de *splicing*.

---

## Actividad 2. Alineamiento y splicing

Responde:

1. ¿Por qué un alineador diseñado únicamente para secuencias genómicas podría tener dificultades con datos de RNA-seq?
2. ¿Qué ocurre cuando un read atraviesa una unión exón–exón?
3. ¿Qué significa que un alineador sea *splice-aware*?
4. ¿Qué ventaja podría proporcionar entregar a HISAT2 sitios de *splicing* conocidos?
5. ¿Significa esto que HISAT2 solamente puede encontrar uniones previamente anotadas? Investiga y fundamenta tu respuesta.

---

# 4. Índices de HISAT2

HISAT2 no busca directamente cada read recorriendo desde cero el archivo FASTA completo.

Para realizar búsquedas eficientes utiliza un conjunto de archivos de **índice**, construidos previamente a partir del genoma de referencia.

Estos archivos ya fueron generados para este práctico.

Busca los archivos asociados al índice dentro de:

```bash
/home/Materiales_Curso/Anotaciones/Danio_rerio_114/
```

Puedes ayudarte utilizando:

```bash
ls -lh /home/Materiales_Curso/Anotaciones/Danio_rerio_114/*.ht2
```

o, dependiendo del tipo de índice disponible:

```bash
ls -lh /home/Materiales_Curso/Anotaciones/Danio_rerio_114/*.ht2l
```

---

## Actividad 3. Reconocimiento del índice

Examina los nombres de los archivos.

### Preguntas

1. ¿Cuántos archivos forman el índice?
2. ¿Qué extensión poseen?
3. ¿Qué parte del nombre comparten?
4. ¿Qué significa el **basename** de un índice HISAT2?
5. ¿Se entrega a HISAT2 el nombre de uno de los archivos `.ht2` o el basename común?
6. ¿Qué ocurriría si se entrega como `-x` una ruta o basename que no corresponde a un índice existente?

---

## Información: ¿cómo se construyeron estos índices?

> **NO EJECUTES LOS SIGUIENTES COMANDOS.**
>
> Los índices ya están construidos y disponibles en la estación.
>
> Esta sección tiene como objetivo mostrar cómo se prepara una referencia desde cero.

La estructura general utilizada por HISAT2 es:

```bash
hisat2-build [opciones] <reference_in> <ht2_base>
```

donde:

- `<reference_in>` corresponde al genoma de referencia.
- `<ht2_base>` corresponde al nombre base que compartirán los archivos del índice.

Consulta la documentación:

```bash
hisat2-build --help
```

### Preguntas

1. ¿Por qué construir un índice puede tomar considerablemente más tiempo que simplemente abrir el FASTA?
2. ¿Por qué el índice se construye una sola vez y luego puede reutilizarse para múltiples muestras?
3. ¿Qué problemas podrían aparecer si se utiliza un índice construido a partir de una versión diferente del genoma que la utilizada para la anotación GTF?

---

# 5. Sitios conocidos de splicing

HISAT2 puede recibir un archivo que contenga sitios de *splicing* conocidos.

Este archivo puede generarse a partir de una anotación GTF utilizando una herramienta incluida con HISAT2.

La estructura general es:

```bash
hisat2_extract_splice_sites.py genes.gtf > splicesites.txt
```

En nuestra estación este archivo **ya fue generado**.

Busca dentro del directorio de anotaciones el archivo:

```text
Danio_rerio.GRCz11.114.ss
```

Examina sus primeras líneas:

```bash
head /home/Materiales_Curso/Anotaciones/Danio_rerio_114/Danio_rerio.GRCz11.114.ss
```

---

## Actividad 4. Sitios de splicing

1. ¿Qué información contiene cada línea?
2. ¿De dónde se obtuvo esta información?
3. ¿Qué relación existe entre este archivo y el GTF?
4. ¿Por qué entregar esta información puede facilitar el alineamiento de RNA-seq?
5. ¿Qué diferencia conceptual existe entre el **índice del genoma** y el **archivo de sitios de splicing**?

---

# 6. Preparación del espacio de trabajo

Regresa al directorio del análisis:

```bash
cd ~/Bioinformatica_Omicas/RNAseq
```

Crea los directorios que utilizaremos:

```bash
mkdir -p results/hisat2
mkdir -p results/sam
mkdir -p results/bam
mkdir -p results/samtools
```

Comprueba la estructura:

```bash
find results -maxdepth 2 -type d
```

Los archivos tendrán funciones diferentes:

```text
results/
├── hisat2/      → reportes generados durante el alineamiento
├── sam/         → archivos SAM
├── bam/         → archivos BAM ordenados e indexados
└── samtools/    → estadísticas generadas con SAMtools
```

---

# 7. Ambiente de HISAT2

HISAT2 se encuentra instalado en un ambiente Conda independiente.

Carga Conda:

```bash
source /usr/local/miniconda3/etc/profile.d/conda.sh
```

Activa el ambiente:

```bash
conda activate Hisat2
```

Comprueba que HISAT2 está disponible:

```bash
hisat2 --version
```

Consulta su documentación:

```bash
hisat2 --help
```

Registra la versión utilizada.

---

# 8. Construcción del comando de alineamiento

Antes de ejecutar HISAT2 deberás comprender los parámetros que utilizarás.

Considera el siguiente esquema:

```bash
hisat2 \
    -x INDEX \
    -1 READ1 \
    -2 READ2 \
    -S OUTPUT.sam \
    --known-splicesite-infile SPLICE_SITES \
    [OTRAS_OPCIONES]
```

> Este esquema **no es un comando ejecutable**.
>
> Debes identificar las rutas, nombres de archivos y opciones apropiadas para nuestro experimento.

---

## Actividad 5. Parámetros de HISAT2

Consulta:

```bash
hisat2 --help
```

y la documentación oficial de HISAT2.

Investiga qué función cumplen:

```text
-x
-1
-2
-S
-q
-p
--known-splicesite-infile
--rna-strandness
--dta
```
Nota: para informarte sobre direccionalidad de la librerias revisa el siguiente link: [Xiaofei Carl Zang](https://x-zang.github.io/blog/check-strandness/)

Completa:

| Parámetro | Función | ¿Lo utilizaremos? | Justificación |
|---|---|---|---|
| `-x` | | | |
| `-1` | | | |
| `-2` | | | |
| `-S` | | | |
| `-q` | | | |
| `-p` | | | |
| `--known-splicesite-infile` | | | |
| `--rna-strandness` | | | |
| `--dta` | | | |

### Preguntas

1. ¿Qué información espera HISAT2 después de `-x`?
2. ¿Qué diferencia existe entre `-1` y `-2`?
3. ¿Qué formato de salida genera `-S`?
4. ¿Es necesario indicar `-q` cuando nuestros archivos son FASTQ? ¿Qué ocurre si no se especifica?
5. ¿Qué controla `-p`?
6. ¿Qué información necesita conocerse antes de utilizar `--rna-strandness`?
7. ¿Podemos asumir que una biblioteca es stranded solamente porque estamos trabajando con RNA-seq?
8. ¿Para qué tipo de análisis posterior está diseñada la opción `--dta`?
9. ¿Todos los parámetros disponibles deben utilizarse? Fundamenta.

---

# 9. Alineamiento piloto

Antes de procesar las seis muestras realizaremos un alineamiento piloto utilizando:

```text
PEC_WT_1
```

Los archivos de entrada deben corresponder a los FASTQ procesados obtenidos en el Práctico 2.

---

## Actividad 6. Construcción del comando piloto

Construye el comando de HISAT2 para `PEC_WT_1`.

Deberás definir correctamente:

- Ruta al índice.
- Basename del índice.
- FASTQ procesado correspondiente a R1.
- FASTQ procesado correspondiente a R2.
- Archivo de sitios conocidos de *splicing*.
- Archivo SAM de salida.
- Número de threads.
- Parámetros adicionales que estén justificados.

El archivo SAM deberá guardarse dentro de:

```text
results/sam/
```

El reporte de HISAT2 deberá guardarse dentro de:

```text
results/hisat2/
```

### Antes de ejecutar

Comprueba:

- Que R1 y R2 pertenecen a `PEC_WT_1`.
- Que estás utilizando los archivos procesados.
- Que el basename del índice existe.
- Que el archivo de sitios de *splicing* existe.
- Que la salida se escribirá dentro de `results/sam/`.
- Que no sobrescribirás ningún FASTQ.
- Que el reporte de HISAT2 quedará almacenado.

---

# 10. STDOUT y STDERR

Antes de ejecutar el alineamiento es importante comprender dónde escribe cada programa sus resultados.

En Linux existen diferentes flujos de salida.

Investiga la diferencia entre:

```text
>
2>
&>
```

### Preguntas

1. ¿Qué significa `stdout`?
2. ¿Qué significa `stderr`?
3. ¿Cuál utiliza HISAT2 para imprimir su resumen de alineamiento?
4. ¿Qué diferencia existe entre:

```bash
> reporte.txt
```

y:

```bash
2> reporte.txt
```

5. ¿Qué hace:

```bash
&> reporte.txt
```

6. ¿Cuál utilizarías para guardar el reporte de HISAT2 y por qué?

Una vez decidido, ejecuta el alineamiento piloto y guarda el reporte.

---

# 11. Interpretación del reporte de HISAT2

Cuando HISAT2 finaliza entrega un resumen del alineamiento.

Abre el reporte generado:

```bash
cat results/hisat2/ARCHIVO_REPORTE
```

---

## Actividad 7. Interpretación del alineamiento piloto

Registra los resultados principales.

| Métrica | Resultado | Interpretación |
|---|---:|---|
| Total de pares | | |
| Concordant 0 times | | |
| Concordant exactly 1 time | | |
| Concordant >1 times | | |
| Discordant alignments | | |
| Mates alineados individualmente | | |
| Overall alignment rate | | |

### Preguntas

1. ¿Qué significa que un par se alinee concordantemente?
2. ¿Qué significa `concordantly exactly 1 time`?
3. ¿Qué significa `concordantly >1 times`?
4. ¿Qué características del genoma podrían producir múltiples alineamientos?
5. ¿Qué significa que un par no se alinee concordantemente?
6. ¿Puede uno de los mates alinearse aunque el par completo no lo haga concordantemente?
7. ¿Qué representa `overall alignment rate`?
8. ¿Por qué una tasa de alineamiento alta no garantiza por sí sola que el experimento sea de buena calidad?
9. ¿Qué causas biológicas o técnicas podrían producir una tasa de alineamiento baja?
10. ¿Cómo podrían ayudarte los resultados de FastQC, FastQ Screen y Cutadapt a interpretar este resultado?

---

# 12. Exploración del archivo SAM

HISAT2 genera un archivo **SAM —Sequence Alignment/Map—**.

Antes de transformarlo, examina su contenido.

```bash
head results/sam/ARCHIVO.sam
```

Las líneas que comienzan con:

```text
@
```

corresponden al **header**.

Puedes observar solamente registros de alineamiento utilizando:

```bash
grep -v "^@" results/sam/ARCHIVO.sam | head
```

---

## Actividad 8. Estructura de un archivo SAM

Investiga los campos obligatorios del formato SAM:

| Campo | Significado |
|---|---|
| QNAME | |
| FLAG | |
| RNAME | |
| POS | |
| MAPQ | |
| CIGAR | |
| RNEXT | |
| PNEXT | |
| TLEN | |
| SEQ | |
| QUAL | |

### Preguntas

1. ¿Qué representa `RNAME`?
2. ¿Qué representa `POS`?
3. ¿Qué información entrega `MAPQ`?
4. ¿MAPQ corresponde a la calidad de secuenciación de las bases? Explica.
5. ¿Qué representa el campo CIGAR?
6. ¿Cómo puede CIGAR representar un alineamiento que atraviesa un intrón?
7. ¿Por qué el campo FLAG no debe interpretarse simplemente como un número decimal independiente?
8. ¿Qué información relacionada con paired-end puede codificarse mediante FLAG?
9. ¿Qué diferencia existe entre los campos obligatorios y los tags opcionales del formato SAM?

---

# 13. De SAM a BAM

Los archivos SAM contienen la información de alineamiento en formato de texto.

Para los análisis posteriores utilizaremos archivos **BAM**, que almacenan esta información en una representación binaria más compacta.

La manipulación de estos archivos se realizará utilizando **SAMtools**.

HISAT2 y SAMtools se encuentran instalados en ambientes Conda diferentes.

Primero desactiva el ambiente actual:

```bash
conda deactivate
```

Activa:

```bash
conda activate samtools
```

Comprueba la instalación:

```bash
samtools --version
```

Consulta:

```bash
samtools --help
```

Registra la versión utilizada.

---

# 14. Conversión, ordenamiento e indexación

Investiga los siguientes comandos:

```text
samtools view
samtools sort
samtools index
```

---

## Actividad 9. SAMtools

Completa:

| Comando | Función | Entrada | Salida |
|---|---|---|---|
| `samtools view` | | | |
| `samtools sort` | | | |
| `samtools index` | | | |

### Preguntas

1. ¿Qué diferencia existe entre SAM y BAM?
2. ¿Qué ventaja posee BAM respecto de SAM?
3. ¿Qué significa ordenar un BAM por coordenadas?
4. ¿Por qué muchos programas requieren BAM ordenados?
5. ¿Qué información contiene un archivo `.bai`?
6. ¿Por qué el índice `.bai` permite acceder rápidamente a regiones específicas del genoma?

---

# 15. Procesamiento del alineamiento piloto

## 15.1. Conversión SAM → BAM

Consulta:

```bash
samtools view --help
```

Construye un comando que transforme:

```text
results/sam/PEC_WT_1.sam
```

en un archivo BAM.

Guarda el resultado dentro de:

```text
results/bam/
```

No elimines todavía el SAM.

---

## 15.2. Ordenamiento

Consulta:

```bash
samtools sort --help
```

Ordena el BAM por coordenadas.

Utiliza una nomenclatura que permita reconocer claramente que se trata de un archivo ordenado.

Por ejemplo:

```text
PEC_WT_1.sorted.bam
```

---

## 15.3. Indexación

Consulta:

```bash
samtools index --help
```

Genera el índice correspondiente al BAM ordenado.

Después del proceso deberías disponer al menos de:

```text
PEC_WT_1.sorted.bam
PEC_WT_1.sorted.bam.bai
```

Comprueba:

```bash
ls -lh results/bam/
```

---

# 16. Evaluación del BAM piloto

SAMtools posee distintas herramientas para obtener estadísticas sobre los alineamientos.

Investiga:

```bash
samtools flagstat
samtools stats
samtools idxstats
```

---

## Actividad 10. Estadísticas del alineamiento

Ejecuta `samtools flagstat` sobre el BAM ordenado de `PEC_WT_1`.

Guarda el resultado dentro de:

```text
results/samtools/
```

Utiliza una nomenclatura informativa, por ejemplo:

```text
PEC_WT_1.flagstat.txt
```

Examina el reporte.

### Preguntas

1. ¿Cuántos reads aparecen en el archivo?
2. ¿Cuántos están mapeados?
3. ¿Qué porcentaje está mapeado?
4. ¿Cuántos reads están correctamente pareados?
5. ¿Qué significa `properly paired`?
6. ¿Cuántos reads corresponden a read1?
7. ¿Cuántos corresponden a read2?
8. ¿Por qué esperarías números similares para read1 y read2?
9. ¿Qué diferencias existen entre la información proporcionada por HISAT2 y `samtools flagstat`?

---

# 17. Alineamiento de las seis muestras

Una vez validado el procedimiento utilizando `PEC_WT_1`, aplica la misma estrategia a las seis muestras.

Recuerda que **cada muestra debe alinearse independientemente**.

| Orden | Muestra |
|---:|---|
| 1 | `PEC_WT_1` |
| 2 | `PEC_WT_2` |
| 3 | `PEC_WT_3` |
| 4 | `PEC_MUT_1` |
| 5 | `PEC_MUT_2` |
| 6 | `PEC_MUT_3` |

Para cada muestra deberás obtener:

1. Reporte de HISAT2.
2. Archivo SAM.
3. BAM.
4. BAM ordenado.
5. Índice del BAM.
6. Reporte `flagstat`.

---

## Actividad 11. Registro de los alineamientos

Completa:

| Muestra | Pares procesados | Concordant 0 veces | Concordant 1 vez | Concordant >1 vez | Overall alignment rate |
|---|---:|---:|---:|---:|---:|
| `PEC_WT_1` | | | | | |
| `PEC_WT_2` | | | | | |
| `PEC_WT_3` | | | | | |
| `PEC_MUT_1` | | | | | |
| `PEC_MUT_2` | | | | | |
| `PEC_MUT_3` | | | | | |

### Preguntas

1. ¿Las seis muestras presentan tasas de alineamiento similares?
2. ¿Existe alguna muestra que se comporte como un posible outlier?
3. ¿Las diferencias entre réplicas son mayores o menores que las diferencias entre genotipos?
4. ¿Alguna muestra presenta una proporción particularmente alta de reads no alineados?
5. ¿Alguna muestra presenta una proporción particularmente alta de alineamientos múltiples?
6. ¿Qué información adicional revisarías antes de concluir que una muestra posee un problema?

---

# 18. Comparación HISAT2 vs SAMtools

Completa una segunda tabla utilizando `flagstat`:

| Muestra | Reads totales | Reads mapeados | Mapeados (%) | Properly paired | Properly paired (%) |
|---|---:|---:|---:|---:|---:|
| `PEC_WT_1` | | | | | |
| `PEC_WT_2` | | | | | |
| `PEC_WT_3` | | | | | |
| `PEC_MUT_1` | | | | | |
| `PEC_MUT_2` | | | | | |
| `PEC_MUT_3` | | | | | |

### Preguntas

1. ¿Qué mide HISAT2 cuando informa `overall alignment rate`?
2. ¿Qué mide SAMtools cuando informa reads mapeados?
3. ¿Por qué estas métricas no deben interpretarse automáticamente como exactamente la misma cantidad?
4. ¿Qué información adicional entrega `flagstat` sobre la estructura paired-end de los datos?

---

# 19. Verificación de archivos

Comprueba los BAM:

```bash
ls -lh results/bam/
```

Utiliza además:

```bash
samtools quickcheck -v results/bam/*.sorted.bam
```

Investiga qué evalúa `samtools quickcheck`.

### Preguntas

1. ¿Qué significa que `samtools quickcheck` no produzca ninguna salida?
2. ¿Garantiza esto que todos los alineamientos sean biológicamente correctos?
3. ¿Qué tipo de problemas puede detectar?
4. ¿Qué tipo de problemas no puede detectar?

---

# 20. Integración de resultados con MultiQC

Los reportes generados por HISAT2 pueden integrarse mediante MultiQC para facilitar la comparación entre muestras.

Recuerda que MultiQC se encuentra disponible en el ambiente:

```text
RNAseq
```

Cambia nuevamente de ambiente:

```bash
conda deactivate
conda activate RNAseq
```

Consulta:

```bash
multiqc --help
```

Genera un reporte que incorpore los resultados de alineamiento.

Guarda el reporte dentro de un nuevo directorio:

```bash
mkdir -p results/multiqc_alignment
```

El reporte deberá llamarse:

```text
Practico_03_MultiQC.html
```

---

## Actividad 12. Evaluación global del alineamiento

Utiliza los reportes individuales y MultiQC para comparar las seis muestras.

### Preguntas

1. ¿Las seis bibliotecas muestran un comportamiento comparable?
2. ¿Existe alguna muestra que se separe claramente del resto?
3. ¿La proporción de reads alineados de forma única es similar?
4. ¿Existen diferencias importantes en alineamientos múltiples?
5. ¿Los resultados son coherentes con la calidad observada antes y después del preprocesamiento?
6. ¿Algún problema observado en FastQC o FastQ Screen parece reflejarse en el alineamiento?
7. ¿Consideras que todas las muestras pueden continuar a la etapa de cuantificación? Fundamenta utilizando evidencia.

---

# 21. Discusión: ¿qué significa un buen alineamiento?

Una tasa de alineamiento alta es deseable, pero **no constituye por sí sola una demostración de calidad**.

La interpretación debe considerar conjuntamente:

- Calidad inicial de los reads.
- Preprocesamiento realizado.
- Organismo y genoma de referencia.
- Versión del ensamblaje.
- Protocolo de preparación de bibliotecas.
- Estructura paired-end.
- Proporción de alineamientos únicos.
- Proporción de alineamientos múltiples.
- Proporción de pares concordantes.
- Posible contaminación.
- Complejidad de la biblioteca.

### Preguntas de discusión

1. ¿Una tasa de alineamiento de 100% sería necesariamente ideal?
2. ¿Por qué algunos reads pueden alinearse en múltiples regiones?
3. ¿Cómo afectan las regiones repetitivas al alineamiento?
4. ¿Qué efecto podría tener utilizar el genoma de una especie incorrecta?
5. ¿Qué efecto podría tener utilizar una versión diferente del ensamblaje?
6. ¿Qué ocurriría si el GTF y el genoma corresponden a versiones incompatibles?
7. ¿Cómo podría una contaminación detectada mediante FastQ Screen afectar el alineamiento?
8. ¿Por qué no deberíamos eliminar automáticamente una muestra solamente porque posee una tasa de alineamiento menor que las demás?

---

# 22. Documentación del análisis

Redacta un párrafo para la sección de **Materiales y Métodos** que incluya:

- Herramienta de alineamiento.
- Versión de HISAT2.
- Genoma de referencia utilizado.
- Versión del ensamblaje.
- Tipo de datos.
- Estrategia paired-end.
- Utilización de sitios conocidos de *splicing*.
- Parámetros relevantes utilizados.
- Conversión a BAM.
- Ordenamiento.
- Indexación.
- Versión de SAMtools.
- Métodos utilizados para evaluar el alineamiento.

No escribas solamente una lista de comandos.

La metodología debe permitir que otra persona comprenda y reproduzca el análisis.

---

Redacta además una breve sección de **Resultados** que describa:

- Rango de tasas de alineamiento observado.
- Proporción de alineamientos únicos y múltiples.
- Comportamiento de los pares.
- Consistencia entre réplicas.
- Posibles muestras atípicas.
- Relación entre los resultados de alineamiento y el QC previo.
- Decisión respecto de qué muestras continuarán a cuantificación.

---

# 23. Transferencia de resultados

No es necesario descargar los archivos SAM o BAM al computador personal para realizar el informe.

Los reportes de texto y el reporte MultiQC son suficientes para la discusión de los resultados.

Los comandos `scp` deben ejecutarse desde una terminal de tu computador personal y no desde la sesión SSH abierta en la estación.

Primero sal de la estación:

```bash
exit
```

Luego puedes transferir, por ejemplo:

```bash
scp -r usuario@152.74.15.46:~/Bioinformatica_Omicas/RNAseq/results/hisat2 .
```

```bash
scp -r usuario@152.74.15.46:~/Bioinformatica_Omicas/RNAseq/results/samtools .
```

```bash
scp -r usuario@152.74.15.46:~/Bioinformatica_Omicas/RNAseq/results/multiqc_alignment .
```

Reemplaza `usuario` por tu nombre de usuario.

---

# 24. Presentación y discusión de resultados

Los resultados serán discutidos en clase utilizando el mismo PPT colaborativo del módulo de RNA-seq.

Cada presentación deberá incluir:

- Estrategia de alineamiento.
- Genoma de referencia utilizado.
- Parámetros relevantes de HISAT2.
- Resultado del alineamiento piloto.
- Comparación de las seis muestras.
- Porcentaje de alineamiento.
- Alineamientos únicos y múltiples.
- Comportamiento paired-end.
- Resultados de `samtools flagstat`.
- Posibles muestras atípicas.
- Relación con el QC de los prácticos anteriores.
- Limitaciones del análisis.
- Decisión respecto de la cuantificación.

Las actividades desarrolladas desde el Práctico 1 formarán parte de un **único informe asociado al módulo de RNA-seq**.

No se debe entregar un informe independiente para este práctico.

---

# Resultados mínimos esperados

Al finalizar el práctico deberías contar con una estructura equivalente a:

```text
RNAseq/
├── data/
│   ├── raw/
│   │   └── ...
│   └── processed/
│       ├── PEC_WT_1_R1.processed.fastq.gz
│       ├── PEC_WT_1_R2.processed.fastq.gz
│       ├── ...
│       ├── PEC_MUT_3_R1.processed.fastq.gz
│       └── PEC_MUT_3_R2.processed.fastq.gz
│
└── results/
    ├── hisat2/
    │   ├── PEC_WT_1.hisat2.log
    │   ├── PEC_WT_2.hisat2.log
    │   ├── ...
    │   └── PEC_MUT_3.hisat2.log
    │
    ├── sam/
    │   ├── PEC_WT_1.sam
    │   ├── ...
    │   └── PEC_MUT_3.sam
    │
    ├── bam/
    │   ├── PEC_WT_1.sorted.bam
    │   ├── PEC_WT_1.sorted.bam.bai
    │   ├── ...
    │   ├── PEC_MUT_3.sorted.bam
    │   └── PEC_MUT_3.sorted.bam.bai
    │
    ├── samtools/
    │   ├── PEC_WT_1.flagstat.txt
    │   ├── ...
    │   └── PEC_MUT_3.flagstat.txt
    │
    └── multiqc_alignment/
        └── Practico_03_MultiQC.html
```

Los nombres exactos dependerán de la nomenclatura utilizada durante el práctico, pero deben ser **consistentes e inequívocos**.

---

# Antes de finalizar

Comprueba que:

- [ ] Utilizaste los FASTQ procesados obtenidos en el Práctico 2.
- [ ] Identificaste correctamente el genoma de referencia.
- [ ] Identificaste el GTF correspondiente.
- [ ] Identificaste el basename del índice HISAT2.
- [ ] Identificaste el archivo de sitios conocidos de *splicing*.
- [ ] Registraste la versión de HISAT2.
- [ ] Comprendiste los parámetros utilizados en el alineamiento.
- [ ] Alineaste R1 y R2 conjuntamente.
- [ ] Guardaste un reporte de HISAT2 para cada muestra.
- [ ] Evaluaste primero una muestra piloto.
- [ ] Interpretaste el resumen de alineamiento.
- [ ] Examinaste la estructura de un archivo SAM.
- [ ] Registraste la versión de SAMtools.
- [ ] Convertiste los alineamientos a BAM.
- [ ] Ordenaste los BAM por coordenadas.
- [ ] Indexaste los BAM.
- [ ] Generaste estadísticas con `samtools flagstat`.
- [ ] Comparaste las seis muestras.
- [ ] Verificaste la integridad básica de los BAM.
- [ ] Integraste los resultados mediante MultiQC.
- [ ] Evaluaste los resultados junto con el QC de los prácticos anteriores.
- [ ] Fundamentaste qué muestras pueden continuar a cuantificación.
- [ ] Conservaste los comandos utilizados para garantizar la reproducibilidad.
