# Práctico N°6: Enriquecimiento funcional de datos de RNA-seq

## INTRODUCCIÓN

En los prácticos anteriores identificamos genes diferencialmente expresados entre las condiciones experimentales. Sin embargo, una lista de genes diferencialmente expresados por sí sola no necesariamente permite comprender qué procesos biológicos están siendo afectados.

El **análisis de enriquecimiento funcional** busca determinar si determinados procesos biológicos, funciones moleculares, componentes celulares o vías de señalización se encuentran representados con mayor frecuencia de la esperada dentro de una lista de genes.

Una estrategia ampliamente utilizada es el **Over-Representation Analysis (ORA)**. En este análisis se compara nuestra lista de genes de interés —por ejemplo, genes diferencialmente expresados— contra una población de referencia o **background**.

En este práctico utilizaremos dos herramientas web para explorar funcionalmente nuestros resultados de RNA-seq:

- **DAVID**
- **g:Profiler**

Ambas herramientas permiten analizar categorías funcionales como **Gene Ontology (GO)** y diferentes bases de datos de vías biológicas.

---

## OBJETIVOS

Al finalizar este práctico deberías ser capaz de:

1. Generar listas de genes adecuadas para análisis de enriquecimiento funcional.
2. Interpretar enriquecimientos de Gene Ontology.
3. Identificar vías biológicas asociadas a los genes diferencialmente expresados.
4. Comprender la importancia de definir correctamente el **background** del análisis.
5. Comparar resultados obtenidos utilizando DAVID y g:Profiler.
6. Integrar los resultados del enriquecimiento funcional con los resultados del análisis de expresión diferencial.

---

# Actividad 1. Preparación de las listas de genes

A partir de los resultados obtenidos en el análisis de expresión diferencial, genera al menos dos listas:

- **Genes significativamente sobreexpresados (UP)**
- **Genes significativamente subexpresados (DOWN)**

Utiliza los criterios de significancia definidos durante los prácticos anteriores.

Cada lista debe contener únicamente los identificadores de los genes, idealmente **Ensembl Gene IDs**.

Por ejemplo:

```text
ENSMUSG00000000001
ENSMUSG00000000028
ENSMUSG00000000037
ENSMUSG00000000049
```

> **Importante:** los genes UP y DOWN deben analizarse inicialmente **por separado**.

### ¿Por qué?

Un mismo proceso biológico puede contener genes cuya expresión aumenta y otros cuya expresión disminuye. Si juntamos ambas listas desde el inicio, podemos perder información respecto de la dirección de la respuesta transcriptómica.

Antes de continuar, registra:

| Lista | Número de genes |
|---|---:|
| Genes UP | |
| Genes DOWN | |
| Total de genes expresados | |

---

# Actividad 2. Enriquecimiento funcional utilizando DAVID

DAVID (**Database for Annotation, Visualization and Integrated Discovery**) permite identificar anotaciones funcionales sobrerrepresentadas dentro de una lista de genes.

Ingresa a la página de **DAVID Bioinformatics Resources** y selecciona la herramienta de análisis funcional.

## 2.1 Cargar una lista de genes

Comenzaremos analizando los **genes DOWN**.

En el menú para cargar genes:

1. Ingresa la lista de identificadores Ensembl correspondiente a los genes DOWN.
2. Selecciona el organismo correspondiente a nuestro experimento.
3. Selecciona **ENSEMBL_GENE_ID** como tipo de identificador, si es solicitado.
4. Define la lista como **Gene List**.
5. Envía la lista.

Revisa cuántos genes fueron reconocidos correctamente por DAVID.

> **Pregunta:** ¿DAVID reconoció todos los genes de tu lista?  
> Si algunos genes no fueron reconocidos, ¿qué razones podrían explicar esta diferencia?

---

## 2.2 Functional Annotation

Selecciona la herramienta **Functional Annotation**.

Explora el panel **Functional Annotation Summary**.

DAVID organiza las anotaciones en diferentes categorías. Para este práctico nos concentraremos principalmente en:

### Gene Ontology

Explora:

- **Biological Process (BP)**
- **Molecular Function (MF)**
- **Cellular Component (CC)**

Selecciona una categoría y abre el **Functional Annotation Chart**.

Observa las columnas disponibles.

Presta especial atención a:

- número de genes asociados al término;
- porcentaje de genes;
- significancia estadística;
- corrección por comparaciones múltiples;
- genes responsables del enriquecimiento.

Selecciona algunos de los términos más significativos y determina **qué genes de tu lista están generando el enriquecimiento**.

### Preguntas

- ¿Cuál es el proceso biológico significativamente enriquecido que presenta mayor evidencia estadística?
- ¿Cuántos genes de tu lista participan en ese proceso?
- ¿Son los mismos genes responsables de varios términos enriquecidos?
- ¿Existen términos GO muy similares o redundantes?

> **Importante:** encontrar muchos términos GO significativos no significa necesariamente que existan muchos procesos biológicos independientes. La estructura jerárquica de Gene Ontology genera frecuentemente términos relacionados o parcialmente redundantes.

---

# Actividad 3. Exploración de pathways

Dentro de **Functional Annotation**, explora las categorías asociadas a vías biológicas.

Por ejemplo:

- KEGG
- Reactome
- WikiPathways

Las bases disponibles pueden depender del organismo analizado.

Identifica las vías significativamente enriquecidas.

Selecciona al menos **tres pathways** que consideres relevantes y registra:

| Pathway | Genes observados | Significancia ajustada | Interpretación |
|---|---:|---:|---|
| | | | |
| | | | |
| | | | |

Cuando esté disponible, abre el pathway y observa dónde se localizan los genes diferencialmente expresados dentro de la vía.

### Pregunta

> ¿Los pathways identificados son coherentes con los resultados de Gene Ontology y con la biología del experimento?

No basta con que una vía sea estadísticamente significativa. Debes evaluar si existe una relación razonable entre:

**genes diferencialmente expresados → funciones enriquecidas → vías biológicas → fenómeno experimental**

---

# Actividad 4. ¿Qué estamos utilizando como background?

Hasta ahora hemos realizado el análisis utilizando el **background definido por defecto por DAVID**.

Pero consideremos nuestro experimento de RNA-seq.

No todos los genes del genoma tienen necesariamente la misma posibilidad de aparecer como diferencialmente expresados. Por ejemplo, genes que no se expresan en nuestras muestras difícilmente podrán ser detectados como DEGs.

Por esta razón, para un experimento transcriptómico puede ser más apropiado utilizar como universo estadístico los **genes que fueron efectivamente evaluados en el análisis diferencial**.

Es decir:

```text
Genes diferencialmente expresados
              ↓
        Query gene list

Genes evaluados/expresados
              ↓
       Statistical background
```

Utiliza ahora como background la lista de genes expresados proporcionada para el práctico.

Carga esta lista en DAVID como **Background** o **Population** y repite el análisis de los genes DOWN.

Compara los resultados obtenidos utilizando:

1. background por defecto;
2. background experimental.

### Preguntas

- ¿Cambió el número de términos significativos?
- ¿Cambió el orden de los términos?
- ¿Aparecieron o desaparecieron procesos biológicos?
- ¿Por qué el background puede modificar el resultado estadístico del enriquecimiento?

---

# Actividad 5. Repetir el análisis con los genes UP

Realiza ahora el mismo procedimiento utilizando los **genes UP**.

Analiza:

- GO Biological Process;
- GO Molecular Function;
- GO Cellular Component;
- pathways.

Utiliza como background los genes expresados en el experimento.

### Pregunta

> ¿Los genes UP y DOWN representan procesos biológicos diferentes?

Identifica procesos que parezcan:

- activados o aumentados;
- disminuidos;
- compartidos entre ambas listas.

Recuerda que el análisis de enriquecimiento **no demuestra directamente que una vía esté activada o inhibida**. El resultado indica que genes asociados a una determinada función se encuentran sobrerrepresentados en nuestra lista.

---

# Actividad 6. Enriquecimiento funcional utilizando g:Profiler

Ahora repetiremos parte del análisis utilizando **g:Profiler**.

g:Profiler permite analizar múltiples fuentes de anotación funcional y además facilita la comparación entre diferentes listas de genes.

Ingresa a **g:Profiler** y selecciona **g:GOSt – Functional profiling**.

---

## 6.1 Análisis inicial

Carga primero la lista de genes DOWN.

Selecciona el organismo correspondiente al experimento.

Ejecuta el análisis.

g:Profiler mostrará un gráfico general de enriquecimiento y una tabla con los términos significativamente enriquecidos.

Explora las siguientes fuentes:

- **GO:BP** — Biological Process
- **GO:MF** — Molecular Function
- **GO:CC** — Cellular Component
- **KEGG**
- **Reactome**
- **WikiPathways**

Dependiendo del organismo, algunas bases pueden no estar disponibles.

---

## 6.2 Interpretar la tabla de resultados

Observa las columnas disponibles.

Entre ellas encontrarás información relacionada con:

- término funcional;
- fuente de anotación;
- significancia estadística;
- tamaño del término;
- número de genes de tu lista asociados al término;
- genes responsables del enriquecimiento.

Selecciona **3–5 términos** que consideres biológicamente relevantes.

No selecciones términos únicamente porque presentan el menor valor de *p*.

Considera también:

- especificidad del término;
- número de genes involucrados;
- redundancia con otros términos;
- relación con el fenómeno experimental.

---

# Actividad 7. Incorporar un background experimental en g:Profiler

Abre las opciones avanzadas de g:Profiler.

Busca la opción asociada al:

**Statistical domain scope**

Selecciona un **Custom background**.

Ingresa como background la lista de genes expresados/evaluados en nuestro experimento.

Repite el análisis.

### Compara:

| Resultado | Background por defecto | Background experimental |
|---|---|---|
| Número de términos significativos | | |
| Principal GO:BP | | |
| Principal pathway | | |
| Cambios relevantes | | |

### Pregunta

> ¿Qué análisis consideras estadísticamente más apropiado para interpretar este experimento de RNA-seq? Justifica.

---

# Actividad 8. Comparar genes UP y DOWN simultáneamente

Una característica útil de g:Profiler es que permite analizar más de una lista de genes.

Carga simultáneamente:

```text
>UP
gene1
gene2
gene3
...

>DOWN
gene4
gene5
gene6
...
```

Ejecuta el análisis utilizando el background experimental.

Utiliza la visualización para identificar términos:

- exclusivos de genes UP;
- exclusivos de genes DOWN;
- enriquecidos en ambas listas.

### Pregunta

> ¿Existen funciones biológicas compartidas entre UP y DOWN? Si las hay, ¿significa necesariamente que el proceso está siendo regulado en la misma dirección?

---

# Actividad 9. Comparación DAVID vs g:Profiler

Compara los resultados obtenidos utilizando ambas herramientas.

No es esperable que sean idénticos.

Completa la siguiente tabla:

| Característica | DAVID | g:Profiler |
|---|---|---|
| GO Biological Process | | |
| Pathway principal | | |
| Número de términos significativos | | |
| Genes asociados a los términos principales | | |
| Facilidad de interpretación | | |

### Discusión

- ¿Qué resultados fueron consistentes entre las dos herramientas?
- ¿Qué diferencias observaste?
- ¿Por qué dos herramientas de enriquecimiento funcional podrían producir resultados diferentes usando la misma lista de genes?
- ¿Qué aspectos deberías reportar en un artículo científico para permitir que el análisis sea reproducible?

---

# Actividad 10. Integración biológica

Finalmente, integra los resultados obtenidos en los prácticos de RNA-seq.

Considera conjuntamente:

```text
Control de calidad
        ↓
Alineamiento / cuantificación
        ↓
Expresión diferencial
        ↓
Genes UP / DOWN
        ↓
Enriquecimiento funcional
        ↓
Interpretación biológica
```

Selecciona los resultados que consideres más importantes y construye una interpretación global del experimento.

Tu interpretación debería intentar responder:

**¿Qué procesos celulares o moleculares parecen modificarse entre las condiciones experimentales y qué genes sustentan esta conclusión?**

---

# Resultados que debes guardar

Guarda los resultados necesarios para incluir posteriormente en tu informe de RNA-seq.

Como mínimo:

1. Tabla de enriquecimiento funcional para genes **UP**.
2. Tabla de enriquecimiento funcional para genes **DOWN**.
3. Resultados utilizando el **background experimental**.
4. Una visualización representativa de DAVID o g:Profiler.
5. Una tabla resumida con los principales procesos o pathways.
6. Los genes que sustentan los enriquecimientos que consideres biológicamente más relevantes.

No es necesario incluir en el informe todas las categorías significativas obtenidas.

La selección de resultados debe ser **justificada biológicamente**.

---

## Para pensar

El enriquecimiento funcional es una herramienta de apoyo para interpretar una lista de genes y no constituye por sí solo una demostración de mecanismo.

Antes de concluir que un determinado proceso está alterado, pregúntate:

> **¿Qué genes están generando este enriquecimiento?**

> **¿Cuál fue el universo de genes utilizado como background?**

> **¿El término está sustentado por muchos genes o solamente por unos pocos?**

> **¿Existen términos redundantes describiendo esencialmente el mismo fenómeno?**

> **¿La dirección de cambio de los genes es consistente con la interpretación que estoy proponiendo?**

> **¿Los resultados son coherentes con el contexto biológico del experimento?**
