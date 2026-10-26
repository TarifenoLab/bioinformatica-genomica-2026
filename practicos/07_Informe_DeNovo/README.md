# Informe: Ensamble de novo y evaluación del transcriptoma

## Extensión

**Máximo 5 páginas**, sin considerar:

- portada
- referencias
- anexos

Los scripts utilizados deberán incluirse como **anexo** y no dentro del cuerpo principal del informe.

El informe debe privilegiar la presentación e interpretación de resultados por sobre la descripción detallada de comandos.

---

# 1. Introducción y objetivo

**Máximo ½ página.**

Explicar brevemente:

- qué es un ensamble de novo de transcriptoma
- en qué situación se utiliza
- qué desafíos presenta su evaluación

Finalizar indicando el objetivo del análisis realizado.

---

# 2. Estrategia de análisis

Presentar un esquema simple del flujo de trabajo realizado.

Por ejemplo:

```text
RNA-seq preprocesado
        ↓
      Trinity
        ↓
   Transcriptoma
        ↓
 ┌──────┼────────┬─────────┐
 ↓      ↓        ↓         ↓
N50   Reads    ExN50     BUSCO
      mapping
 └──────┴────────┴─────────┘
        ↓
Evaluación integrada
```

Explicar brevemente la **racional** de cada análisis:

- ¿qué evalúa?
- ¿por qué es necesario?

No describir comandos ni parámetros técnicos en esta sección.

---

# 3. Resultados

Esta corresponde a la **sección principal del informe**.

## 3.1 Ensamble

Reportar al menos:

| Métrica | Resultado |
|---|---:|
| “Genes” Trinity | |
| Transcritos/isoformas | |
| Bases ensambladas | |
| Longitud media | |
| Longitud mediana | |
| N50 | |
| N50, isoforma más larga/gene | |

Explicar brevemente qué indica la relación entre número de transcritos y número de “genes” Trinity.

---

## 3.2 Representación de reads

Presentar el porcentaje de reads que pueden alinearse nuevamente contra el transcriptoma.

Interpretar el resultado:

- ¿el transcriptoma representa adecuadamente los reads utilizados para construirlo?
- ¿qué limitaciones tiene esta métrica?

---

## 3.3 N50 y ExN50

Presentar:

- N50
- E90N50
- gráfico Ex vs ExN50

Explicar las diferencias conceptuales entre ambas métricas.

Interpretar la forma general de la curva ExN50 y discutir si el resultado es consistente con un transcriptoma adecuadamente reconstruido.

---

## 3.4 BUSCO

Presentar gráficamente o mediante una tabla:

- Complete
- Single-copy
- Duplicated
- Fragmented
- Missing

Interpretar qué indican estos resultados sobre la **completitud biológica** del transcriptoma.

No limitarse a describir porcentajes.

---

# 4. Integración y conclusión

Esta sección debe responder:

> **¿Qué tan bueno es el transcriptoma generado y qué evidencia respalda esta conclusión?**

Integrar al menos:

- continuidad del ensamble
- representación de los reads
- expresión/ExN50
- completitud BUSCO

Discutir posibles inconsistencias entre métricas.

Por ejemplo:

- N50 alto pero BUSCO bajo
- mapping alto pero BUSCO incompleto
- BUSCO alto pero abundantes transcritos cortos
- elevada proporción de BUSCO duplicados

La conclusión debe basarse en el **conjunto de resultados**, no en una sola métrica.

---

# 5. Anexos

## Anexo 1. Scripts

Incluir los scripts utilizados para:

- Trinity
- alineamiento de reads contra el transcriptoma
- RSEM
- construcción de matriz de expresión
- ExN50
- BUSCO

Los scripts deben estar ordenados y claramente identificados.

No es necesario incluir archivos `.log` completos, pero pueden incorporarse fragmentos relevantes si permiten demostrar o explicar algún resultado.

---

# Evaluación sugerida

| Sección | Puntaje |
|---|---:|
| Introducción y objetivo | 10 |
| Racional y estrategia de análisis | 15 |
| Resultados del ensamble | 15 |
| Representación de reads | 10 |
| N50 y ExN50 | 15 |
| BUSCO | 15 |
| Integración y conclusión | 15 |
| Reproducibilidad, presentación y anexos | 5 |
| **Total** | **100** |

### Criterio central de evaluación

Se evaluará principalmente la capacidad de **comprender qué información entrega cada análisis, interpretar correctamente los resultados e integrar distintas métricas para evaluar la calidad del transcriptoma**.

La mera ejecución correcta de los programas o presentación de sus outputs, sin interpretación, no será suficiente.
