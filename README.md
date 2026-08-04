# Proyecto integrador — Redes Complejas (1217)

Maestría en Ciencias de la Ingeniería Eléctrica · IV cohorte
Universidad de Cuenca · Dr. Fabián Astudillo-Salinas

**Caso de estudio:** la red de datos de la Universidad de Cuenca —
177 equipos y 209 enlaces reconstruidos a partir de los 34 diagramas del
informe técnico de topologías y conexiones de red.

> ### 📅 Entrega: 26 de agosto de 2026
> Informe PDF + repositorio de código + diapositivas. **25 puntos.**
> Grupos de 2 o 3 integrantes.

---

## Empiecen por aquí

1. Lean **`proyecto_red_ucuenca.pdf`** — el enunciado completo: 11 problemas
   agrupados en 5 fases, rúbrica detallada y estructura sugerida del informe.
2. Ejecuten el código base para comprobar que su entorno funciona:

```bash
cd codigo_base

# Python
python3 -m pip install -r requirements.txt
python3 cargar_red.py

# Julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. cargar_red.jl
```

Debe imprimir **«Todas las comprobaciones pasaron»** y terminar con código de
salida `0`. Si no, revisen su instalación antes de seguir.

3. Abran **`red_ucuenca_noc.html`** en el navegador para explorar la red de
   forma interactiva, y lean **`red_ucuenca_README.md`** (el diccionario de
   datos) antes de escribir la primera línea de análisis.

---

## Contenido del paquete

### Enunciado

| Archivo | Contenido |
|---|---|
| `proyecto_red_ucuenca.pdf` | **Enunciado del proyecto.** Problemas P1–P11, rúbrica, entregables y anexos |

### Conjunto de datos

| Archivo | Contenido |
|---|---|
| `red_ucuenca_nodes.csv` | 177 nodos: `id`, `label`, `campus`, `capa`, `diagrams` |
| `red_ucuenca_edges.csv` | 209 aristas: `source`, `target`, `label`, `trafico_mbps`, `capacidad_mbps`, `rol`, `diagrams` |
| `red_ucuenca.graphml` | El mismo grafo con atributos — para Gephi, yEd, Cytoscape, NetworkX, igraph |
| `red_ucuenca_README.md` | **Diccionario de datos.** Qué significa cada atributo y qué correcciones se aplicaron |
| `red_ucuenca_grafo.png` | Visualización estática de referencia |
| `red_ucuenca_noc.html` | Vista interactiva tipo NOC: mapa, filtros por sede y por rol, tabla de enlaces |

El grafo es **simple, no dirigido y conexo**. Los identificadores son idénticos
en los CSV y en el GraphML, así que ambas fuentes producen el mismo grafo.

### Código base

| Archivo | Contenido |
|---|---|
| `codigo_base/cargar_red.py` | Carga + verificación + resumen, en Python con NetworkX |
| `codigo_base/cargar_red.jl` | Lo mismo en Julia con Graphs.jl |
| `codigo_base/requirements.txt` | Dependencias de Python |
| `codigo_base/Project.toml` | Dependencias de Julia |
| `codigo_base/README.md` | Qué hace y qué **no** hace el código base |

Los dos scripts son equivalentes: elijan uno (o usen ambos y contrasten) y
construyan las cinco fases encima. **Cargan y verifican, pero no ponderan las
aristas, no estiman capacidades y no calculan ninguna métrica**: eso es el
proyecto.

### Código de referencia

Material del módulo que el enunciado les pide reutilizar. Citar de dónde
proviene cada pieza es obligatorio.

| Ruta | Se usa en |
|---|---|
| `codigo_referencia/kmeans/` | P4 — partición de la red |
| `codigo_referencia/ford-fulkerson/` | P6 — flujo máximo, con animaciones paso a paso |
| `codigo_referencia/edmonds-karp/` | P6 — flujo máximo con onda BFS |
| `codigo_referencia/actividad_flujo_maximo.pdf` | P6 — la actividad de laboratorio en la que se apoya |

### Procedencia de los datos

`red_ucuenca_correcciones.py` documenta las 61 correcciones aplicadas a la
transcripción original de los diagramas. Se incluye para que puedan verificar
cómo se construyó el conjunto de datos; **no hace falta ejecutarlo** — los
archivos que reciben ya están corregidos. (Ejecutado sobre ellos avisa de eso y
no modifica nada.)

---

## Notas

- **No hay que auditar los datos.** El conjunto que reciben es consistente:
  identificadores unificados, sin aristas duplicadas, tráfico y capacidad en
  columnas separadas. Trabajen sobre él directamente.
- **`trafico_mbps` no es capacidad.** Es la lectura del monitor en el instante
  de la captura. La capacidad nominal solo está documentada en los 28 enlaces
  del diagrama MPLS (`capacidad_mbps`); estimar la de los demás es parte de P6.
- **Todo resultado debe ser reproducible.** El docente ejecutará su
  repositorio; un número del informe que no aparezca al ejecutarlo se califica
  sobre cero.
- El material completo del módulo está en
  <https://github.com/fabianastudillo/ComplexNetworks>.
