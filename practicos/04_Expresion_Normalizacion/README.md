# Práctico 4. Estimación de expresión génica y análisis de expresión diferencial

## Introducción

En este práctico trabajaremos con los archivos BAM obtenidos en el práctico de alineamiento. Utilizaremos **HTSeq-count** para estimar la expresión a nivel de gen y **DESeq2** para identificar genes diferencialmente expresados entre las muestras wild-type y mutantes *pax6b-/-*.

## Actividad 1. Estimación de la expresión génica

HTSeq-count cruza los alineamientos con la anotación génica para obtener conteos por gen. En datos paired-end, cada par se cuenta como una unidad correspondiente a un fragmento.

Consulta el [manual de HTSeq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html). En este práctico utilizaremos el modo `union`, contando alineamientos sobre exones agrupados por `gene_id`.

### Preparación de los archivos

Desde el directorio de trabajo, crea las carpetas de salida:

```bash
cd ~/Bioinformatica_Omicas/RNAseq
mkdir -p results/htseq_results/deseq2
```

Los BAM del práctico anterior están ordenados por coordenadas. HTSeq-count admite ese orden, pero aquí utilizaremos archivos ordenados por nombre para mantener juntos los alineamientos de cada par.

Genera una nueva copia ordenada por nombre:

```bash
samtools sort -n \
    -o results/bam/PEC_WT_1.name_sorted.bam \
    results/bam/PEC_WT_1.sorted.bam
```

Conserva el BAM original ordenado por coordenadas. Repite este procedimiento para las seis muestras.

### Conteo por gen

Activa el ambiente de HTSeq-count:

```bash
conda activate htseq-count
```

Identifica la anotación GTF correspondiente al genoma utilizado durante el alineamiento:

```bash
ls -lh /home/Materiales_Curso/Anotaciones/Danio_rerio_114/*.gtf*
```

Ejecuta el siguiente comando reemplazando `ANOTACION.gtf` por el nombre del archivo correspondiente:

```bash
htseq-count \
    --format=bam \
    --order=name \
    --stranded=no \
    --type=exon \
    --idattr=gene_id \
    --mode=union \
    results/bam/PEC_WT_1.name_sorted.bam \
    /home/Materiales_Curso/Anotaciones/Danio_rerio_114/ANOTACION.gtf \
    > results/htseq/PEC_WT_1.count
```

Repite el conteo para las seis muestras, manteniendo los mismos parámetros y utilizando nombres de salida que permitan identificar cada muestra.

### Ejercicio 1. Parámetros y resultados del conteo

1. Explica las diferencias entre los modos `union`, `intersection-strict` e `intersection-nonempty`.
2. Justifica el uso de `--type=exon`, `--idattr=gene_id` y `--mode=union`.
3. ¿Qué indica `--stranded=no`? ¿Por qué la modalidad paired-end no determina por sí sola la direccionalidad de una biblioteca?
4. ¿Qué ocurre con los fragmentos que solapan más de un gen? ¿Y con los que se encuentran completamente dentro de intrones?
5. Examina un archivo `.count` y explica qué representan sus dos columnas.
6. Compara entre muestras las categorías `__no_feature`, `__ambiguous` y `__alignment_not_unique`. ¿Observas alguna muestra con resultados diferentes de las demás?

## Actividad 2. Preparación de los datos para DESeq2

Crea un archivo llamado `metadata.txt` en el directorio de trabajo. Debe contener las siguientes columnas separadas por **tabulaciones**:

```text
Sample	File	Condition
PEC_WT_1	PEC_WT_1.count	WT
PEC_WT_2	PEC_WT_2.count	WT
PEC_WT_3	PEC_WT_3.count	WT
PEC_MUT_1	PEC_MUT_1.count	MUT
PEC_MUT_2	PEC_MUT_2.count	MUT
PEC_MUT_3	PEC_MUT_3.count	MUT
```

Inicia R en el ambiente donde se encuentre instalado DESeq2. Carga el paquete e importa los conteos:

```r
library(DESeq2)

setwd("~/Bioinformatica_Omicas/RNAseq")

metadata <- read.delim("metadata.txt")

metadata$Condition <- factor(
    metadata$Condition,
    levels = c("WT", "MUT")
)

dds <- DESeqDataSetFromHTSeqCount(
    sampleTable = metadata,
    directory = "results/htseq",
    design = ~ Condition
)
```

DESeq2 debe recibir **conteos crudos**, ya que realiza la normalización durante el análisis. Consulta la [documentación de DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

Exporta la matriz de conteos:

```r
write.csv(
    counts(dds),
    "results/deseq2/conteos_crudos.csv"
)
```

### Ejercicio 2. Matriz de expresión

1. Comprueba que la tabla contiene una fila por gen y una columna por muestra.
2. ¿Por qué las categorías especiales de HTSeq que comienzan con `__` no deben tratarse como genes?
3. ¿Por qué no debemos entregar conteos normalizados como entrada a DESeq2?
4. Adjunta la tabla de conteos crudos como anexo al informe.

## Actividad 3. Análisis de expresión diferencial

Utilizaremos wild-type como condición de referencia. Ejecuta el análisis y obtiene los resultados para la comparación **mutante respecto de wild-type**:

```r
dds <- DESeq(dds)

res <- results(
    dds,
    contrast = c("Condition", "MUT", "WT"),
    alpha = 0.05
)

summary(res)
```

En esta comparación:

- Un `log2FoldChange` positivo indica mayor expresión en mutantes.
- Un `log2FoldChange` negativo indica mayor expresión en wild-type.
- Consideraremos diferencialmente expresados los genes con `padj < 0.05`.

Ordena y guarda los resultados:

```r
resOrdered <- res[order(res$padj), ]

write.csv(
    as.data.frame(resOrdered),
    "results/deseq2/expresion_diferencial.csv"
)

sig <- subset(
    resOrdered,
    !is.na(padj) & padj < 0.05
)

write.csv(
    as.data.frame(sig),
    "results/deseq2/genes_significativos.csv"
)
```

Consulta los factores de normalización y exporta los conteos normalizados:

```r
sizeFactors(dds)

write.csv(
    counts(dds, normalized = TRUE),
    "results/deseq2/conteos_normalizados.csv"
)
```

### Ejercicio 3. Interpretación de los resultados

1. Explica qué representan `baseMean`, `log2FoldChange`, `pvalue` y `padj`.
2. ¿Cuántos genes presentan expresión diferencial con `padj < 0.05`?
3. ¿Cuántos tienen mayor expresión en mutantes y cuántos en wild-type?
4. ¿Por qué utilizamos el valor de *p* ajustado?
5. ¿Qué representan los factores de normalización? Compara los conteos crudos y normalizados de un mismo gen entre muestras.

## Actividad 4. Visualización e interpretación biológica

### Gráfico MA

Genera un gráfico que relacione la abundancia media de cada gen con su cambio de expresión:

```r
pdf("results/deseq2/MAplot.pdf")

plotMA(
    res,
    alpha = 0.05,
    main = "Mutante vs wild-type",
    ylim = c(-5, 5)
)

dev.off()
```

### Análisis de componentes principales

Utiliza una transformación estabilizadora de la varianza para explorar la relación entre las muestras:

```r
vsd <- varianceStabilizingTransformation(
    dds,
    blind = TRUE
)

pdf("results/deseq2/PCA.pdf")

print(plotPCA(vsd, intgroup = "Condition"))

dev.off()
```

Esta transformación se utiliza para explorar y visualizar los datos. El análisis de expresión diferencial se realiza a partir de los conteos crudos.

### Ejercicio 4. Visualización y selección de genes

1. Describe los ejes del gráfico MA. ¿Qué representan los puntos destacados?
2. Examina el PCA. ¿Las réplicas de cada condición se agrupan? ¿Existe separación entre genotipos?
3. Si alguna muestra se separa de sus réplicas, revisa sus resultados de control de calidad y alineamiento. Discute posibles explicaciones sin excluirla únicamente por su posición en el PCA.
4. Selecciona tres genes con mayor expresión en mutantes y tres con mayor expresión en wild-type, utilizando los criterios `padj < 0.05` y `|log2FoldChange| ≥ 1`. Si no hay suficientes genes que cumplan ambos criterios, indícalo.
5. Para cada gen seleccionado, registra su identificador, nombre si está disponible, `baseMean`, `log2FoldChange` y `padj`. Revisa también sus conteos normalizados en las seis muestras.
6. Investiga sus funciones utilizando literatura científica y discute su posible relación con el contexto del experimento.

**Recuerda:** un cambio de expresión grande no implica necesariamente una expresión abundante. Considera ambas características al interpretar los genes seleccionados.

## Presentación de resultados

Integra los resultados en el informe del módulo de RNA-seq y en el PPT colaborativo. No se requiere un informe independiente para este práctico.

Incluye:

- Los parámetros utilizados para el conteo y su justificación.
- La comparación de los resultados de HTSeq entre muestras.
- Las tablas de conteos crudos y normalizados como anexos.
- El número de genes diferencialmente expresados en cada dirección.
- Los gráficos MA y PCA con su interpretación.
- La tabla de genes seleccionados y su discusión biológica, con referencias.

Conserva los comandos y el script de R utilizados para garantizar la reproducibilidad del análisis.
