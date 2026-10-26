# Práctico N°5: Ensamble de novo y evaluación de transcriptoma

## Introducción

En este práctico realizaremos el **ensamble de novo de un transcriptoma a partir de datos de RNA-seq paired-end**, utilizando **Trinity**.

A diferencia de los análisis realizados previamente, aquí no utilizaremos un genoma de referencia para reconstruir los transcritos. El ensamble será generado directamente a partir de la información contenida en los reads.

Trabajaremos mezclando las muestras de WT y MUT con las que hemos trabajado hasta ahora. 

Posteriormente evaluaremos la calidad del ensamble utilizando distintas aproximaciones. Es importante considerar que **no existe una única métrica capaz de definir si un transcriptoma está bien ensamblado**, por lo que los resultados deberán interpretarse de manera integrada.

### Objetivos

Al finalizar este práctico deberás ser capaz de:

- Comprender la lógica general de un ensamble de novo de transcriptoma.
- Construir un ensamble utilizando Trinity.
- Reconocer la estructura de los identificadores generados por Trinity.
- Evaluar estadísticas básicas del ensamble.
- Determinar qué proporción de los reads originales está representada en el transcriptoma.
- Calcular e interpretar N50 y ExN50.
- Evaluar la completitud del transcriptoma utilizando BUSCO.
- Integrar distintas métricas para emitir una conclusión sobre la calidad del ensamble.

---

# Actividad 1. Ensamble de novo con Trinity

Utilizaremos **Trinity**, programa diseñado para reconstruir transcriptomas a partir de datos de RNA-seq.

Trinity implementa tres componentes principales: **Inchworm, Chrysalis y Butterfly**, que participan secuencialmente en la reconstrucción y resolución de los transcritos.

Puedes consultar la documentación oficial:

[Trinity RNA-seq de novo assembly](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

## 1.1 Preparación de los datos

Trabajaremos con el set de 6 muestras con las que hemos trabajado hasta ahora. Utiliza los archivos **paired-end preprocesados**.

En esta etapa las muestras serán utilizadas conjuntamente para generar **un único ensamble de referencia**.

Por lo tanto, deberás entregar a Trinity:

- los seis archivos correspondientes al read R1 mediante `--left`
- los seis archivos correspondientes al read R2 mediante `--right`

Los nombres de los archivos deben separarse mediante comas.

> **Importante:** las librerías utilizadas en este curso son **unstranded**.

Para un ensamble de novo, la calidad del preprocesamiento es particularmente importante. Errores, adaptadores residuales o secuencias contaminantes pueden incorporarse directamente al transcriptoma ensamblado.

---

## 1.2 Construcción del comando

Trinity se encuentra instalado en el ambiente Conda:

```bash
trinity
```

Activa el ambiente y consulta la ayuda del programa:

```bash
conda activate trinity

Trinity --help
```

A partir de la documentación, construye un script que permita ejecutar el ensamble.

La estructura general será:

```bash
Trinity \
    --seqType <TIPO_SECUENCIA> \
    --left <R1_REP1,R1_REP2,R1_REP3> \
    --right <R2_REP1,R2_REP2,R2_REP3> \
    --CPU <NUMERO_CPU> \
    --max_memory <MEMORIA> \
    --output <DIRECTORIO_SALIDA>
```

### Restricción de recursos

Debido a que el servidor es compartido:

```text
--max_memory no debe superar 10G
```

Selecciona además un número razonable de CPU.

No copies simplemente la estructura anterior: revisa qué significa cada parámetro y reemplaza los campos indicados.

---

## 1.3 Ejecución

Como en los prácticos anteriores:

1. Genera un script Bash.
2. Hazlo ejecutable.
3. Ejecuta el análisis utilizando `nohup`.
4. Guarda la salida estándar y los mensajes del programa en un archivo `.log`.

Por ejemplo, la estructura general de ejecución será:

```bash
nohup ./ensamble_trinity.sh > ensamble_trinity.log 2>&1 &
```

Puedes monitorear la ejecución utilizando:

```bash
tail -f ensamble_trinity.log
```

El principal resultado del análisis será el archivo:

```text
Trinity.fasta
```

---

## 1.4 Explorar el ensamble

Observa las primeras secuencias:

```bash
head Trinity.fasta
```

Examina los identificadores generados por Trinity.

Por ejemplo:

```text
TRINITY_DN100_c0_g1_i1
```

Investiga qué representan:

```text
DN
c
g
i
```

y distingue conceptualmente entre:

- transcriptos/isoformas Trinity
- “genes” Trinity

> Los grupos denominados *genes* por Trinity corresponden a agrupaciones computacionales de transcritos relacionados y no deben interpretarse directamente como genes biológicos anotados.

### Preguntas

1. ¿Por qué utilizamos conjuntamente todas las muestras para construir el transcriptoma?
2. ¿Qué diferencia conceptual existe entre ensamblar RNA-seq con y sin genoma de referencia?
3. ¿Qué podrían representar varias isoformas asociadas a un mismo “gene” de Trinity?

---

# Actividad 2. Estadísticas generales del ensamble

Trinity incluye herramientas para caracterizar el transcriptoma generado.

Utiliza:

```text
TrinityStats.pl
```

Localiza el programa dentro de la instalación de Trinity y ejecútalo sobre:

```text
Trinity.fasta
```

El análisis entrega, entre otros:

- número de “genes” Trinity
- número de transcritos
- porcentaje GC
- longitud media y mediana
- número total de bases ensambladas
- N50

Además, Trinity calcula estadísticas considerando:

1. todos los transcritos ensamblados
2. solamente la isoforma más larga de cada “gene”

### Preguntas

1. ¿Cuántos “genes” Trinity fueron reconstruidos?
2. ¿Cuántos transcritos/isoformas fueron reconstruidos?
3. ¿Cuál es la relación aproximada entre transcritos y “genes”?
4. ¿Cuál es el N50 considerando todos los transcritos?
5. ¿Cuál es el N50 considerando solamente la isoforma más larga por “gene”?
6. ¿Por qué estos dos valores pueden ser diferentes?

> **Importante:** un N50 mayor no significa necesariamente que el transcriptoma sea biológicamente mejor. En transcriptómica, una gran cantidad de isoformas o transcritos largos puede aumentar artificialmente esta métrica.

---

# Actividad 3. Representación de los reads en el ensamble

Una forma de evaluar un transcriptoma es determinar qué proporción de los reads originales puede volver a alinearse contra las secuencias reconstruidas.

Conceptualmente realizaremos:

```text
reads originales
      ↓
Trinity.fasta
      ↓
alineamiento
      ↓
porcentaje de reads representados
```

Revisa la sección correspondiente de:

[Transcriptome Assembly Quality Assessment](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment)

Identifica la herramienta recomendada por Trinity para alinear los reads nuevamente contra el transcriptoma.

Realiza el análisis utilizando los reads de las tres réplicas empleadas para construir el ensamble.

Registra:

- porcentaje de alineamiento global
- porcentaje concordante de los reads paired-end, cuando corresponda

### Preguntas

1. ¿Qué información entrega esta métrica sobre el ensamble?
2. ¿Qué podría significar que una fracción importante de los reads no pueda alinearse nuevamente?
3. ¿Un porcentaje alto de alineamiento garantiza que el transcriptoma sea completo? Justifica.

---

# Actividad 4. ExN50

El **N50 tradicional** considera únicamente las longitudes de los contigs.

Sin embargo, un transcriptoma puede contener una gran cantidad de transcritos de muy baja abundancia o pobremente ensamblados. Por esta razón Trinity propone utilizar el **ExN50**, que incorpora información de abundancia.

Para calcularlo necesitamos primero estimar la expresión de los transcritos.

## 4.1 Estimación de abundancia

Utilizaremos **RSEM**, mediante las herramientas disponibles en Trinity.

A diferencia del ensamble, la cuantificación debe realizarse **por separado para cada réplica biológica**.

La estructura general del análisis será:

```text
Rep1 → RSEM
Rep2 → RSEM
Rep3 → RSEM
          ↓
matriz de abundancia
          ↓
ExN50
```

Consulta:

[Trinity Transcript Quantification](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification)

Identifica y utiliza:

```text
align_and_estimate_abundance.pl
```

con:

```text
--est_method RSEM
```

Realiza la cuantificación independientemente para cada réplica.

Luego identifica la herramienta de Trinity que permite combinar los resultados de las tres muestras en una matriz de expresión.

---

## 4.2 Cálculo de ExN50

Utilizando la matriz de expresión y `Trinity.fasta`, calcula el perfil ExN50.

Identifica dentro de Trinity el script:

```text
contig_ExN50_statistic.pl
```

Realiza el cálculo a nivel de **gene**, dado que reduce parte del efecto producido por isoformas de muy baja expresión.

Genera además el gráfico:

```text
Ex vs ExN50
```

utilizando la herramienta proporcionada por Trinity.

### Preguntas

1. ¿Cuál es el E90N50 de tu ensamble?
2. ¿En qué valor de Ex se alcanza aproximadamente el máximo ExN50?
3. ¿Por qué el E90N50 puede resultar más informativo que el N50 convencional?
4. ¿Qué ocurre con el ExN50 al acercarse a E100 y por qué?

---

# Actividad 5. Evaluación de completitud con BUSCO

Hasta ahora hemos evaluado principalmente características estructurales del ensamble.

Sin embargo, también necesitamos responder:

> **¿Nuestro transcriptoma contiene los genes que biológicamente esperaríamos encontrar?**

Para ello utilizaremos **BUSCO** (*Benchmarking Universal Single-Copy Orthologs*).

BUSCO busca genes ortólogos altamente conservados que se espera estén presentes en organismos pertenecientes a un determinado linaje evolutivo.

Consulta la documentación:

[BUSCO User Guide](https://busco.ezlab.org/busco_userguide)

## 5.1 Selección del modo de análisis

Determina:

- cuál debe ser el `--mode` para analizar un transcriptoma
- qué significa `--lineage_dataset`
- qué linaje corresponde utilizar para los datos empleados en este curso

Construye el comando correspondiente utilizando:

```text
Trinity.fasta
```

como entrada.

```
conda activate busco

busco \
  -i Trinity.fasta \
  -o BUSCO_metazoa \
  -m transcriptome \
  -l /opt/busco_datasets/lineages/metazoa_odb12.2 \
  --offline \
  -c 8
```


Registra las siguientes categorías:

- Complete BUSCOs
- Complete and single-copy BUSCOs
- Complete and duplicated BUSCOs
- Fragmented BUSCOs
- Missing BUSCOs

### Preguntas

1. ¿Qué porcentaje de BUSCOs completos presenta tu transcriptoma?
2. ¿Qué porcentaje aparece duplicado?
3. ¿Qué podrían representar los BUSCOs fragmentados?
4. ¿Qué podrían representar los BUSCOs ausentes?
5. ¿Por qué en un transcriptoma pueden aparecer BUSCOs duplicados sin que necesariamente exista una duplicación genómica?

---

# Actividad 6. Evaluación integrada del transcriptoma

Resume tus principales resultados en una tabla.

| Métrica | Resultado | ¿Qué evalúa? | Interpretación |
|---|---:|---|---|
| N° “genes” Trinity | | Complejidad del ensamble | |
| N° transcritos | | Complejidad/isoformas | |
| N50 | | Continuidad | |
| E90N50 | | Continuidad ponderada por expresión | |
| Reads representados | | Soporte de los datos originales | |
| BUSCO completos | | Completitud biológica | |
| BUSCO fragmentados | | Reconstrucción parcial | |
| BUSCO ausentes | | Potencial incompletitud | |

Finalmente responde:

> **¿Consideras que el transcriptoma generado presenta una calidad suficiente para ser utilizado en análisis posteriores?**

Tu conclusión debe integrar las diferentes métricas. **No evalúes la calidad utilizando únicamente N50 o BUSCO.**

---

# Material que debes conservar

Guarda los siguientes archivos, ya que podrán ser utilizados posteriormente:

```text
Trinity.fasta
```

Resultados de cuantificación RSEM de cada réplica.

Matriz de abundancia generada a partir de RSEM.

Resultados de BUSCO.

Resultados de ExN50.

Scripts utilizados.

Archivos `.log` de los análisis.

> La reproducibilidad del análisis es parte del trabajo bioinformático. Los scripts deben permitir reconstruir claramente cada etapa realizada.
