# Práctico 1. Formato FASTQ y control de calidad de datos RNA-seq

## Introducción

En este práctico comenzaremos el análisis de un experimento de RNA-seq trabajando directamente con los archivos de secuenciación sin procesar.

Los datos provienen del proyecto [[GSE149081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149081)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149081), disponible en el repositorio Gene Expression Omnibus (GEO). El estudio investigó las similitudes transcriptómicas entre las células endocrinas pancreáticas y las células enteroendocrinas del pez cebra, así como las redes regulatorias controladas por el factor de transcripción **Pax6b**.

Durante el curso utilizaremos seis muestras correspondientes a células endocrinas pancreáticas —pancreatic endocrine cells, PECs— aisladas mediante FACS desde embriones de pez cebra de 27 horas post-fertilización:

* Tres réplicas biológicas wild-type.
* Tres réplicas biológicas mutantes `pax6b-/-`.

Las bibliotecas fueron preparadas a partir de RNA total utilizando el protocolo Smart-seq2 y Nextera XT. La secuenciación se realizó en una plataforma Illumina HiSeq 2000, generando aproximadamente 60 millones de reads paired-end de 100 nucleótidos por muestra.

A diferencia de versiones anteriores de este práctico, trabajaremos con los archivos completos. Estos mismos datos se utilizarán en las siguientes etapas del módulo:

1. Control de calidad.
2. Preprocesamiento.
3. Alineamiento al genoma de referencia.
4. Cuantificación de la expresión génica.
5. Análisis de expresión diferencial.
6. Anotación y enriquecimiento funcional.

---

## Objetivos de aprendizaje

Al finalizar este práctico, podrás:

* Reconocer la estructura y los componentes de un archivo FASTQ.
* Comprender la relación entre un fragmento secuenciado y sus reads R1 y R2.
* Interpretar los identificadores y puntajes de calidad generados por una plataforma Illumina.
* Explorar archivos FASTQ comprimidos utilizando Unix y `SeqKit`.
* Describir el diseño experimental y verificar la organización de los datos.
* Evaluar el posible origen biológico y la presencia de contaminación mediante FastQ Screen.
* Evaluar la calidad de los reads mediante FastQC.
* Integrar los resultados de múltiples muestras mediante MultiQC.
* Identificar problemas potenciales que deberán abordarse durante el preprocesamiento.

> En este práctico no modificaremos ni filtraremos las secuencias. Primero realizaremos un diagnóstico de los datos sin procesar. Las decisiones y procedimientos de preprocesamiento se desarrollarán en el siguiente práctico.

---

# 1. Diseño experimental

Las seis muestras analizadas corresponden al siguiente diseño:

| Muestra     | Tipo celular                 |   Genotipo | Réplica | Archivos                                          |
| ----------- | ---------------------------- | ---------: | ------: | ------------------------------------------------- |
| `PEC_WT_1`  | Célula endocrina pancreática |  Wild-type |       1 | `PEC_WT_1_R1.fastq.gz` y `PEC_WT_1_R2.fastq.gz`   |
| `PEC_WT_2`  | Célula endocrina pancreática |  Wild-type |       2 | `PEC_WT_2_R1.fastq.gz` y `PEC_WT_2_R2.fastq.gz`   |
| `PEC_WT_3`  | Célula endocrina pancreática |  Wild-type |       3 | `PEC_WT_3_R1.fastq.gz` y `PEC_WT_3_R2.fastq.gz`   |
| `PEC_MUT_1` | Célula endocrina pancreática | `pax6b-/-` |       1 | `PEC_MUT_1_R1.fastq.gz` y `PEC_MUT_1_R2.fastq.gz` |
| `PEC_MUT_2` | Célula endocrina pancreática | `pax6b-/-` |       2 | `PEC_MUT_2_R1.fastq.gz` y `PEC_MUT_2_R2.fastq.gz` |
| `PEC_MUT_3` | Célula endocrina pancreática | `pax6b-/-` |       3 | `PEC_MUT_3_R1.fastq.gz` y `PEC_MUT_3_R2.fastq.gz` |

En una secuenciación paired-end, cada fragmento de la biblioteca se secuencia desde ambos extremos:

* El archivo R1 contiene la primera lectura de cada fragmento.
* El archivo R2 contiene la segunda lectura del mismo fragmento.

Por lo tanto, cada muestra biológica está representada por dos archivos FASTQ que deben mantenerse sincronizados durante todo el análisis.

## Preguntas iniciales

1. ¿Cuál es la unidad de replicación biológica de este experimento?
2. ¿Cuántas condiciones experimentales se compararán?
3. ¿Cuántas réplicas biológicas existen por condición?
4. ¿Por qué cada muestra está representada por dos archivos?
5. ¿Cuántos archivos FASTQ esperas encontrar en total?
6. ¿Los archivos R1 y R2 representan réplicas experimentales independientes? Fundamenta tu respuesta.

---

# 2. Preparación del espacio de trabajo

## 2.1. Conexión a la estación de cálculo

Conéctate mediante SSH utilizando las credenciales enviadas por correo:

```bash
ssh usuario@152.74.15.46
```

Reemplaza `usuario` por tu nombre de usuario.

## 2.2. Creación de directorios

Crea un directorio para el módulo de RNA-seq:

```bash
mkdir -p ~/Bioinformatica_Omicas/RNAseq/{data/raw,results}
cd ~/Bioinformatica_Omicas/RNAseq
```

Comprueba tu ubicación:

```bash
pwd
```

Examina la estructura creada:

```bash
ls -R
```

La estructura inicial debería ser:

```text
RNAseq/
├── data/
│   └── raw/
└── results/
```

Los datos sin procesar se mantendrán siempre dentro de `data/raw`. Los resultados generados por las distintas herramientas deberán almacenarse dentro de `results`.

## 2.3. Acceso a los datos completos

Los archivos se encuentran en:

```text
/home/Materiales_Curso/Data_RNAseq
```

Examina los archivos disponibles:

```bash
ls -lh /home/Materiales_Curso/Data_RNAseq/
```

Debido al tamaño de los datos completos, utilizaremos enlaces simbólicos en lugar de generar una copia adicional de cada archivo:

```bash
ln -s /home/Materiales_Curso/Data_RNAseq/*.fastq.gz data/raw/
```

Comprueba que los enlaces fueron creados:

```bash
ls -lh data/raw/
```

Los enlaces simbólicos permiten utilizar los archivos compartidos como datos de entrada sin duplicar su contenido dentro de cada cuenta.

> Nunca modifiques, muevas ni elimines los archivos originales almacenados en `/home/Materiales_Curso/Data_RNAseq`. Todos los resultados y archivos procesados deberán guardarse dentro de tu directorio personal.

## 2.4. Verificación del número de archivos

Cuenta los archivos disponibles:

```bash
find data/raw -name "*.fastq.gz" | wc -l
```

Cuenta por separado los archivos R1 y R2:

```bash
find data/raw -name "*_R1.fastq.gz" | wc -l
find data/raw -name "*_R2.fastq.gz" | wc -l
```

### Ejercicios

1. ¿Cuántos archivos FASTQ encontraste?
2. ¿El número de archivos coincide con el diseño experimental?
3. ¿Cada muestra posee un archivo R1 y un archivo R2?
4. ¿Detectas inconsistencias en los nombres?
5. ¿Por qué es importante utilizar una nomenclatura uniforme desde el comienzo del análisis?

---

# 3. Activación del ambiente bioinformático

Las herramientas utilizadas se encuentran instaladas dentro del ambiente Conda denominado `RNAseq`.

Carga Conda mediante:

```bash
source /usr/local/miniconda3/etc/profile.d/conda.sh
```

Activa el ambiente:

```bash
conda activate RNAseq
```

Comprueba las herramientas y registra sus versiones:

```bash
seqkit version
fastqc --version
fastq_screen --version
multiqc --version
```

Las versiones utilizadas deben indicarse posteriormente en la sección de materiales y métodos del informe.

Para salir del ambiente utiliza:

```bash
conda deactivate
```

---

# Actividad 1. Estructura del formato FASTQ

Los datos se encuentran en archivos con extensión `.fastq.gz`.

La extensión `.fastq` indica que el archivo utiliza el formato FASTQ, mientras que `.gz` indica que se encuentra comprimido mediante `gzip`.

Las herramientas bioinformáticas que utilizaremos pueden leer directamente archivos comprimidos, por lo que no es necesario descomprimirlos.

## 1.1. Exploración de un archivo comprimido

Selecciona uno de los archivos y examínalo mediante:

```bash
zless data/raw/PEC_WT_1_R1.fastq.gz
```

Utiliza las flechas para desplazarte y presiona `q` para salir.

También puedes visualizar los primeros dos registros mediante:

```bash
zcat data/raw/PEC_WT_1_R1.fastq.gz | head -n 8
```

> No utilices `more` o `less` directamente sobre archivos comprimidos. Para archivos `.gz` debes utilizar `zmore`, `zless` o `zcat`.

## 1.2. Estructura de un registro FASTQ

Examina cuidadosamente los primeros registros e identifica la función de cada línea.

### Ejercicios

1. ¿Cuántas líneas conforman un registro FASTQ?

2. Describe la información contenida en cada línea:

   * Línea 1.
   * Línea 2.
   * Línea 3.
   * Línea 4.

3. ¿Con qué carácter comienza el identificador de un registro FASTQ?

4. ¿Cuál es la función de la línea que comienza con `+`?

5. ¿Qué caracteres pueden aparecer en la secuencia de nucleótidos?

6. ¿Qué significa la presencia de una base `N`?

7. Compara la longitud de la secuencia con la longitud de la línea de calidad. ¿Qué relación debería existir entre ellas?

8. ¿Por qué cada nucleótido necesita un carácter de calidad asociado?

---

# Actividad 2. Identificadores Illumina y datos paired-end

## 2.1. Interpretación de los encabezados

Examina el encabezado del primer registro de R1 y R2 de la misma muestra:

```bash
zcat data/raw/PEC_WT_1_R1.fastq.gz | head -n 1
zcat data/raw/PEC_WT_1_R2.fastq.gz | head -n 1
```

Dependiendo de la versión del formato Illumina, el encabezado puede contener información como:

* Instrumento de secuenciación.
* Número de corrida.
* Identificador de la flow cell.
* Lane.
* Tile.
* Coordenadas `x` e `y` del clúster.
* Número del read.
* Estado del filtro de calidad.
* Índice o barcode de la muestra.

### Ejercicios

1. Identifica todos los campos que puedas reconocer en el encabezado.

2. Indica el tile y las coordenadas del clúster del primer fragmento.

3. Compara el identificador del primer read de R1 con el primero de R2.

4. ¿Qué campos se mantienen y cuáles cambian?

5. ¿Qué información permite determinar que ambos reads provienen del mismo fragmento?

6. Repite la comparación utilizando el último registro de ambos archivos.

Para examinar el último registro puedes investigar el uso de:

```bash
seqkit range
```

## 2.2. Diferencia entre read, fragmento y par

Explica con tus propias palabras los siguientes conceptos:

* Read.
* Fragmento de biblioteca.
* Read 1.
* Read 2.
* Par de reads.
* Réplica biológica.

> R1 y R2 no constituyen réplicas experimentales. Son dos observaciones del mismo fragmento de la biblioteca.

---

# Actividad 3. Caracterización general con [SeqKit](https://bioinf.shenwei.me/seqkit/)

[`SeqKit`](https://bioinf.shenwei.me/seqkit/) es una herramienta mantenida actualmente que permite explorar y manipular archivos FASTA y FASTQ.

Consulta la ayuda:

```bash
seqkit --help
seqkit stats --help
```

Crea un directorio para los resultados:

```bash
mkdir -p results/seqkit
```

Ejecuta el resumen estadístico sobre todos los archivos:

```bash
seqkit stats -a data/raw/*.fastq.gz
```

Guarda también el resultado en un archivo:

```bash
seqkit stats -a data/raw/*.fastq.gz \
    > results/seqkit/seqkit_stats_raw.txt
```

## Ejercicios

1. Para cada archivo registra:

   * Número de reads.
   * Número total de nucleótidos.
   * Longitud mínima.
   * Longitud promedio.
   * Longitud máxima.
   * Porcentaje de GC.
   * Porcentaje de bases Q20.
   * Porcentaje de bases Q30.

2. ¿Todos los reads tienen la misma longitud?

3. ¿R1 y R2 contienen el mismo número de reads dentro de cada muestra?

4. ¿Existen muestras con un número de reads considerablemente diferente?

5. ¿La proporción de GC es similar entre las réplicas?

6. ¿Se observan diferencias sistemáticas entre R1 y R2?

7. ¿Qué consecuencias tendría que R1 y R2 tuvieran distinto número de registros?

8. Considerando que cada registro FASTQ ocupa cuatro líneas, estima cuántas líneas tendría uno de estos archivos si estuviera descomprimido.

9. ¿Por qué `seqkit stats` es una estrategia más robusta para contar secuencias que contar líneas y dividir por cuatro?

---

# Actividad 4. Evaluación del posible origen de las secuencias

FastQ Screen permite comparar una muestra de reads contra un panel de genomas de referencia. Este análisis se utiliza para evaluar si la composición de una biblioteca es compatible con su origen esperado y para detectar posibles fuentes de contaminación.

Consulta la [[documentación de FastQ Screen](https://stevenwingett.github.io/FastQ-Screen/)](https://stevenwingett.github.io/FastQ-Screen/) y la ayuda del programa:

```bash
fastq_screen --help
```

El archivo de configuración y los índices de referencia ya se encuentran preparados en la estación de cálculo.

## 4.1. Ejecución

Crea el directorio de resultados:

```bash
mkdir -p results/fastq_screen
```

Ejecuta FastQ Screen sobre todos los archivos:

```bash
fastq_screen \
    --outdir results/fastq_screen \
    data/raw/*.fastq.gz
```

Este análisis puede tardar algunos minutos.

Mientras se ejecuta, investiga:

* Qué alineador utiliza FastQ Screen.
* Qué genomas contiene el archivo de configuración.
* Cuántos reads analiza por defecto.
* Qué función cumple la opción `--subset`.

> Aunque estamos trabajando con los archivos completos, FastQ Screen normalmente evalúa un subconjunto de reads para reducir el tiempo de cómputo. Es importante distinguir entre proporcionar un archivo completo como entrada y analizar todos sus registros.

## 4.2. Interpretación

Para cada genoma de referencia, FastQ Screen clasifica los reads según si:

* No mapearon.
* Mapearon una sola vez.
* Mapearon en múltiples posiciones.
* Mapearon también contra otros genomas del panel.

### Ejercicios

1. ¿Cuál es el organismo al que mapea la mayor proporción de reads?

2. ¿El resultado es consistente con la metadata del experimento?

3. ¿R1 y R2 presentan resultados similares?

4. ¿Existen reads compatibles con otros organismos?

5. ¿La proporción observada sugiere una contaminación relevante?

6. ¿Las seis muestras presentan un perfil comparable?

7. ¿Existe alguna muestra que se comporte como un outlier?

8. ¿Qué representan los reads que no mapearon contra ninguno de los genomas incluidos?

9. ¿Un read que mapea contra dos especies necesariamente corresponde a contaminación?

10. ¿Cómo podrían influir en el resultado:

    * Las regiones evolutivamente conservadas.
    * Las secuencias repetitivas.
    * La completitud de los genomas de referencia.
    * La distancia evolutiva entre las especies.
    * El hecho de que se trate de datos de RNA-seq?

> FastQ Screen evalúa compatibilidad con un panel limitado de referencias. Un alto porcentaje de mapeo apoya un posible origen biológico, pero no demuestra por sí solo la identidad taxonómica de la muestra.

---

# Actividad 5. Evaluación de calidad con [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) genera un reporte gráfico con diferentes métricas de calidad de los datos de secuenciación.

Consulta la ayuda:

```bash
fastqc --help
```

También puedes revisar la [[documentación de los módulos de FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).

## 5.1. Ejecución

Crea un directorio para los resultados:

```bash
mkdir -p results/fastqc
```

Ejecuta FastQC sobre los doce archivos:

```bash
fastqc \
    --threads 2 \
    --outdir results/fastqc \
    data/raw/*.fastq.gz
```

Por cada archivo de entrada se generarán:

* Un reporte en formato HTML.
* Un archivo comprimido con los resultados y datos utilizados en el reporte.

## 5.2. Interpretación de los módulos

Revisa los siguientes módulos:

1. **Basic Statistics**
2. **Per base sequence quality**
3. **Per sequence quality scores**
4. **Per tile sequence quality**
5. **Per base sequence content**
6. **Per sequence GC content**
7. **Per base N content**
8. **Sequence Length Distribution**
9. **Sequence Duplication Levels**
10. **Overrepresented sequences**
11. **Adapter Content**

Para cada módulo, determina:

* Qué propiedad de los datos evalúa.
* Qué representa el eje X.
* Qué representa el eje Y.
* Qué patrón se considera esperable.
* Qué patrón presentan tus datos.
* Si existe una diferencia entre R1 y R2.
* Si existe una diferencia entre muestras o condiciones.

> Los estados `PASS`, `WARN` y `FAIL` no deben interpretarse automáticamente como “muestra buena” o “muestra mala”. FastQC aplica umbrales generales que no consideran completamente el tipo de biblioteca ni el protocolo experimental. Por ejemplo, algunas bibliotecas RNA-seq pueden presentar sesgos de composición en los primeros nucleótidos sin que esto implique necesariamente un error técnico.

## 5.3. Ejercicios de interpretación

1. ¿Cuál es la longitud de los reads?

2. ¿Qué codificación se utilizó para los puntajes de calidad?

3. ¿La calidad cambia a medida que aumenta la posición dentro del read?

4. ¿R2 presenta una calidad menor que R1?

5. ¿Existen posiciones con un aumento importante de bases `N`?

6. ¿La composición de A, C, G y T es uniforme a lo largo de los reads?

7. ¿Existe un sesgo en los primeros nucleótidos?

8. ¿Qué procesos técnicos o biológicos podrían explicar ese sesgo?

9. ¿La distribución de GC es unimodal?

10. ¿Alguna muestra muestra una distribución de GC diferente del resto?

11. ¿Se detectan secuencias sobrerrepresentadas?

12. ¿Las secuencias sobrerrepresentadas corresponden a adaptadores, secuencias biológicas abundantes u otro tipo de secuencia?

13. ¿Existe evidencia de contaminación por adaptadores?

14. ¿El nivel de duplicación es similar entre las muestras?

15. ¿Una duplicación elevada constituye necesariamente un problema en RNA-seq de células obtenidas mediante FACS y amplificadas con Smart-seq2? Fundamenta tu respuesta.

---

# Actividad 6. Integración de resultados con [MultiQC](https://seqera.io/multiqc/)

Cuando se analizan muchas muestras, revisar los reportes individualmente dificulta detectar patrones globales y muestras atípicas.

[`MultiQC`](https://seqera.io/multiqc/) permite integrar los resultados generados por diferentes herramientas en un único reporte interactivo.

FastQ Screen es compatible con MultiQC y su documentación recomienda utilizarlo cuando se analizan múltiples muestras.

## 6.1. Ejecución

Crea el directorio de salida:

```bash
mkdir -p results/multiqc
```

Genera el reporte integrado:

```bash
multiqc results/fastqc results/fastq_screen \
    --outdir results/multiqc \
    --filename Practico_01_MultiQC.html
```

## 6.2. Comparación global

Utiliza el reporte MultiQC para completar una tabla como la siguiente:

| Muestra     | Reads R1 | Reads R2 | Longitud | %GC | %Q30 | Calidad R1/R2 | Adaptadores | Secuencias sobrerrepresentadas | FastQ Screen | Observaciones |
| ----------- | -------: | -------: | -------: | --: | ---: | ------------- | ----------- | ------------------------------ | ------------ | ------------- |
| `PEC_WT_1`  |          |          |          |     |      |               |             |                                |              |               |
| `PEC_WT_2`  |          |          |          |     |      |               |             |                                |              |               |
| `PEC_WT_3`  |          |          |          |     |      |               |             |                                |              |               |
| `PEC_MUT_1` |          |          |          |     |      |               |             |                                |              |               |
| `PEC_MUT_2` |          |          |          |     |      |               |             |                                |              |               |
| `PEC_MUT_3` |          |          |          |     |      |               |             |                                |              |               |

### Ejercicios

1. ¿Las réplicas de una misma condición presentan perfiles de calidad comparables?

2. ¿Existe alguna muestra que se comporte como un outlier?

3. ¿Las diferencias principales se asocian con:

   * La muestra biológica.
   * El genotipo.
   * El mate R1/R2.
   * La posición dentro del read.
   * La corrida o tile de secuenciación.

4. ¿Los problemas detectados afectan a todas las muestras o solamente a algunas?

5. ¿Sería metodológicamente apropiado aplicar exactamente el mismo preprocesamiento a todos los archivos? Fundamenta tu respuesta.

---

# Actividad 7. Diagnóstico previo al preprocesamiento

En esta etapa no debes modificar los reads. El objetivo es formular un diagnóstico basado en evidencia que será utilizado en el siguiente práctico.

Completa la siguiente tabla:

| Problema potencial                 | Evidencia observada | Muestras afectadas | R1, R2 o ambos | Posible intervención | Resultado esperado | Riesgo de la intervención |
| ---------------------------------- | ------------------- | ------------------ | -------------- | -------------------- | ------------------ | ------------------------- |
| Calidad reducida al final del read |                     |                    |                |                      |                    |                           |
| Presencia de adaptadores           |                     |                    |                |                      |                    |                           |
| Bases `N`                          |                     |                    |                |                      |                    |                           |
| Reads de baja calidad global       |                     |                    |                |                      |                    |                           |
| Sesgo de composición inicial       |                     |                    |                |                      |                    |                           |
| Secuencias sobrerrepresentadas     |                     |                    |                |                      |                    |                           |
| Contaminación biológica            |                     |                    |                |                      |                    |                           |

Para cada posible intervención, indica:

1. Qué evidencia justifica aplicarla.
2. A qué muestras debería aplicarse.
3. Qué parámetro de calidad esperas mejorar.
4. Qué información podría perderse.
5. Cómo comprobarías posteriormente si la intervención fue beneficiosa.

> Estas propuestas constituyen hipótesis de trabajo. En el siguiente práctico se seleccionarán herramientas y parámetros de preprocesamiento, y se evaluará si las modificaciones mejoran efectivamente los datos.

---

# Transferencia de resultados

Los comandos `scp` deben ejecutarse desde una terminal de tu computador personal, no desde la sesión SSH abierta en la estación.

Sal de la estación:

```bash
exit
```

Transfiere los resultados:

```bash
scp -r usuario@152.74.15.46:~/Bioinformatica_Omicas/RNAseq/results .
```

Reemplaza `usuario` por tu nombre de usuario.

Abre los reportes `.html` utilizando un navegador web.

---

# Presentación y discusión de resultados

Los resultados serán discutidos en clase utilizando el siguiente [[PPT colaborativo](https://docs.google.com/presentation/d/1MHL3N8V4S7amSoAlNC3wo0_cMoCgJkYEf6PqwatC4h4/edit?usp=sharing)](https://docs.google.com/presentation/d/1MHL3N8V4S7amSoAlNC3wo0_cMoCgJkYEf6PqwatC4h4/edit?usp=sharing).

Cada resultado presentado debe incluir:

* Pregunta u objetivo del análisis.
* Herramienta utilizada y versión.
* Resultado principal.
* Interpretación.
* Limitaciones.
* Implicancias para el preprocesamiento.

Las actividades desarrolladas desde este práctico formarán parte de un único informe asociado al módulo de RNA-seq. No se debe entregar un informe independiente para este práctico.

---

# Resultados mínimos esperados

Al finalizar debes contar con:

```text
RNAseq/
├── data/
│   └── raw/
│       ├── PEC_WT_1_R1.fastq.gz
│       ├── PEC_WT_1_R2.fastq.gz
│       ├── ...
│       ├── PEC_MUT_3_R1.fastq.gz
│       └── PEC_MUT_3_R2.fastq.gz
└── results/
    ├── seqkit/
    │   └── seqkit_stats_raw.txt
    ├── fastq_screen/
    │   └── ...
    ├── fastqc/
    │   └── ...
    └── multiqc/
        └── Practico_01_MultiQC.html
```

Antes de finalizar, comprueba que:

* Trabajaste con las seis muestras completas.
* Cada muestra posee R1 y R2.
* No modificaste los datos originales.
* Registraste las versiones de las herramientas.
* Guardaste los resultados en los directorios correspondientes.
* Integraste los resultados mediante MultiQC.
* Elaboraste un diagnóstico previo al preprocesamiento.
