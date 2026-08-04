# Grafo de la red de datos de la Universidad de Cuenca

Grafo reconstruido a partir de los 34 diagramas del documento
**"Diagramas de red final.docx"** (informe técnico de topologías y conexiones de red).
Cada imagen del documento (weathermaps de monitoreo + diagrama MPLS) fue transcrita
a una lista de nodos y aristas, fusionada en un único grafo no dirigido y
posteriormente **saneada** (ver *Correcciones aplicadas*).

## Archivos

| Archivo | Contenido |
|---|---|
| `red_ucuenca.graphml` | Grafo completo con atributos — importable en Gephi, yEd, Cytoscape, NetworkX, igraph, Graphs.jl (vía GraphIO.jl) |
| `red_ucuenca_nodes.csv` | 177 nodos: `id, label, campus, capa, diagrams` |
| `red_ucuenca_edges.csv` | 209 aristas: `source, target, label, trafico_mbps, capacidad_mbps, rol, diagrams` |
| `red_ucuenca_correcciones.py` | Script que produjo la versión saneada a partir de la transcripción original |

El grafo es **simple y no dirigido**: 177 nodos, 209 aristas, una sola componente conexa.

## Atributos

### Nodos

- **`id`** / **`label`**: nombre del equipo tal como aparece en el diagrama. Una
  sola convención: mayúsculas y guion medio, sin espacios. Los identificadores
  son **idénticos** en los CSV y en el GraphML.
- **`campus`**: ubicación física, valor único. `Campus Central` (75),
  `Campus Paraiso` (44), `Campus Balzay` (35), `Campus Yanuncay` (13),
  `Campus Hospitalidad` (5), `Sede Museo` (2), `Sede Centro Historico` (2),
  `Nube MPLS` (1).
- **`capa`**: posición en la jerarquía descrita por el informe. `acceso` (132),
  `agregacion` (27), `core` (5), `wan` (11), `interconexion` (2, los firewalls
  perimetrales).
- **`diagrams`**: número(s) del diagrama de origen; el diagrama *N* corresponde a
  la imagen *N* del documento.

### Aristas

- **`label`**: etiqueta original del enlace (puerto o *lag* y lectura del monitor).
- **`trafico_mbps`**: tráfico **medido** en el instante de la captura, en Mbps.
  Definido en 170 aristas.
- **`capacidad_mbps`**: capacidad **nominal** declarada explícitamente en la
  etiqueta del enlace (10 Gbps, 2×10 Gbps, 1 Gbps). Definida en 28 aristas, todas
  del diagrama 34. En el resto de enlaces la capacidad no está documentada en la
  fuente y debe estimarse.
- **`rol`**: función del enlace.
  - `principal` (157): enlace activo, miembro primario.
  - `wan` (32): enlace de la nube MPLS y de la interconexión entre campus.
  - `respaldo` (12): miembro secundario de un par redundante, sin tráfico
    (etiquetas `lag NNN 0 bps`).
  - `secundario` (5): miembro secundario de un par redundante que sí cursa tráfico.
  - `inferido` (3): ver nota 3 más abajo.
- **`diagrams`**: diagrama(s) de procedencia, separados por `;`.

## Correcciones aplicadas

La transcripción original contenía inconsistencias heredadas del documento
fuente. `red_ucuenca_correcciones.py` las resuelve de forma reproducible
(61 correcciones); ejecutarlo con `--dry-run` lista una por una.

| Id | Corrección |
|---|---|
| C1 | **Identificadores unificados.** Seis nodos usaban espacios (`INTERNET MPLS`, `ROUTER CAMPUS *`) y el GraphML usaba guion bajo mientras los CSV usaban guion medio. Ahora hay una sola convención, idéntica en los tres archivos. |
| C2 | **Aristas duplicadas colapsadas.** Seis pares de aristas «paralelas» no eran enlaces redundantes sino **el mismo enlace físico transcrito desde dos diagramas distintos**, con dos lecturas de tráfico diferentes (p. ej. `BAL-CENTEC-D2 – VLIR-0A-A107`: 13,75 Mbps en el diagrama 17 y 2,23 Mbps en el 24). Se conserva la medición máxima y ambos diagramas en `diagrams`. El grafo pasa de 215 a 209 aristas y deja de ser multigrafo. |
| C3 | **Tráfico separado de capacidad.** El antiguo `bandwidth_mbps` mezclaba dos magnitudes: en el diagrama 34 la etiqueta declara la capacidad del enlace («Trafico MPLS 10 Gbps»), mientras que en el resto es una lectura del weathermap. Se separó en `trafico_mbps` y `capacidad_mbps`. |
| C4 | **`style` sustituido por `rol`.** El atributo de dibujo era ambiguo: `dashed` significaba a la vez «respaldo sin tráfico» y «enlace lógico del MPLS», y 23 aristas `dashed` cursaban tráfico. `rol` tiene semántica explícita y determinada por los datos. |
| C5 | **Etiquetas `*UNKNOWN*`** del monitor sustituidas por campo vacío (2 aristas). |
| C6 | **`campus` con un valor único.** Cinco nodos frontera declaraban dos campus; ahora indican su ubicación física. Las sedes **Museo** y **Centro Histórico**, que no aparecen en la taxonomía de campus del informe y solo existen en el diagrama 34, tienen su propio valor. |
| C7 | **Nuevo atributo `capa`**, derivado de la convención de nombres (`-C<n>` core, `-D<n>` agregación, `-A<n>` acceso) y del rol de los equipos WAN y de los firewalls perimetrales. |

## Decisiones de reconstrucción (transcripción original)

1. Los nombres de edificio del diagrama general de cada campus se unificaron con el
   switch de distribución de su diagrama de detalle (p. ej. "CC-ARQUITECTURA" →
   `CC-ARQUITECTURA-D107`, "AULARIO 1" → `BAL-AUL1-D4`, "A4 ENFERMERIA" →
   `CP-ENFERMERIA-D1`, "C12" del diagrama MPLS → `DT-0A-C12`).
2. "Huayna Cápac" (diagrama 27) y el router Campus Huayna Cápac (diagrama 34) se
   trataron como el mismo nodo.
3. Tres enlaces con `rol=inferido` (diagrama `0`) conectan Yanuncay, Hospitalidad y
   el Consultorio Jurídico a la nube MPLS: el informe indica esa interconexión, pero
   los diagramas no dibujan explícitamente el enlace router↔switch de núcleo.

## Notas sobre el documento fuente

Detectadas al sanear, **no** corregibles en los datos:

- El `.docx` rotula dos veces "Diagrama 21" y su texto termina en el 33, mientras
  que las imágenes llegan a 34. La columna `diagrams` usa el índice real de
  imagen; a partir del segundo "Diagrama 21" el rótulo del texto va desfasado en
  uno.
- La descripción del "Diagrama 23" (CEA y VLIR) repite el texto del diagrama 20
  ("switch BAL-CENTEC-D2 … edificio Centro tecnológico").
- El informe afirma redundancia física completa *core*↔agregación en **Balzay y
  Paraíso**. Los datos confirman Balzay (seis equipos con doble enlace a
  `DT-0A-C12` y `DT-0A-C13`) pero **no** Paraíso: hay un solo switch de core
  (`CPAR-C10`), y el doble enlace de sus siete equipos va al mismo core, lo que
  es agregación de puertos y no redundancia de núcleo.
- El informe describe el Campus Central con enlaces simples. Eso vale para
  agregación↔acceso; a nivel *core*↔agregación los trece switches `CC-*` sí están
  doblemente conectados a `DATCC-2A-C2` y `DATCC-2A-C3`.

## Uso rápido

```python
import networkx as nx
G = nx.read_graphml("red_ucuenca.graphml")   # 177 nodos, 209 aristas, 1 componente
```

```julia
using Graphs, GraphIO, GraphIO.GraphML
G = loadgraph("red_ucuenca.graphml", "UCuenca", GraphMLFormat())
```

Para trabajar con los atributos en Julia conviene partir de los CSV: ver
`codigo_base/cargar_red.jl`.

> `red_ucuenca_grafo.png` y `red_ucuenca_noc.html` son renderizaciones de la
> transcripción original y no reflejan las correcciones C1–C7.
