# Código base — Proyecto integrador de Redes Complejas

Punto de partida para el proyecto descrito en
[`../proyecto_red_ucuenca.pdf`](../proyecto_red_ucuenca.pdf).

Ambos scripts hacen exactamente lo mismo en los dos lenguajes del módulo:
cargan la red de datos de la Universidad de Cuenca, **verifican la carga**
contra los valores del Anexo A del enunciado e imprimen un resumen estructural.
Elijan uno de los dos (o usen los dos y contrasten) y construyan las cinco
fases del proyecto encima.

## Ejecución

```bash
# Python
python3 -m pip install -r requirements.txt
python3 cargar_red.py                  # desde los CSV
python3 cargar_red.py --fuente graphml # desde el GraphML

# Julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. cargar_red.jl
```

Ambos terminan con código de salida `0` si las trece comprobaciones del
Anexo A pasan, y `1` si alguna falla. Los datos se buscan en la raíz del
repositorio; para usar otra ubicación, exporten `RED_UCUENCA_DIR` (o pasen
`--dir` en Python).

## Qué hace y qué no hace

Hace:

- Lee `red_ucuenca_nodes.csv` + `red_ucuenca_edges.csv` o `red_ucuenca.graphml`.
  Los identificadores son idénticos en las dos fuentes, así que ambas producen
  el mismo grafo.
- Devuelve un grafo **simple y no dirigido**: 177 nodos, 209 aristas, una sola
  componente conexa.
- Conserva los atributos `campus`, `capa`, `label`, `diagrams` en los nodos y
  `trafico_mbps`, `capacidad_mbps`, `rol`, `label`, `diagrams` en las aristas.

No hace, a propósito:

- **No pondera** las aristas. Elegir el modelo de peso (saltos, latencia,
  carga) es parte de P5.
- **No estima capacidades.** `capacidad_mbps` viene dado solo para los 28
  enlaces del diagrama MPLS; estimar los 181 restantes es el primer punto de
  P6. Ojo: `trafico_mbps` es el tráfico *medido* en el instante de la captura,
  no una capacidad — son columnas distintas justamente para que no se
  confundan.
- **No calcula** centralidades, comunidades, caminos mínimos, flujos ni
  percolación. Eso es el proyecto.

## Uso como biblioteca

```python
from cargar_red import cargar_red, verificar, resumen
G = cargar_red()             # networkx.Graph
verificar(G)
```

```julia
include("cargar_red.jl")
red = cargar_red()           # red.g, red.ids, red.idx,
verificar(red)               # red.campus, red.capa, red.nodos, red.aristas
```

## Material relacionado del repositorio

| Ruta | Contenido |
|---|---|
| `../red_ucuenca_README.md` | Diccionario de datos completo y correcciones aplicadas al conjunto de datos |
| `../codigo_referencia/kmeans/` | k-means con ejemplos (P4) |
| `../codigo_referencia/ford-fulkerson/` | Ford-Fulkerson con BFS/DFS, corte mínimo, animaciones (P6) |
| `../codigo_referencia/edmonds-karp/` | Edmonds-Karp con onda BFS (P6) |
| `../codigo_referencia/actividad_flujo_maximo.pdf` | Actividad de laboratorio de flujo máximo, base directa de P6 |
| `intro/codes/` (en GitHub) | Carga y visualización de grafos, distribución de grado, Louvain en Julia |

Reutilizar y adaptar ese código es correcto y esperado; **citar de dónde
proviene cada pieza es obligatorio**.
