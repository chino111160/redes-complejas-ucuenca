# Análisis del proyecto integrador — red de datos UCuenca

Código de las cinco fases (P1–P11) del enunciado (`../proyecto_red_ucuenca.pdf`),
en Julia. Reutiliza `../codigo_base/cargar_red.jl` (carga y verificación) y
`../codigo_referencia/` (Ford-Fulkerson, Edmonds-Karp, k-means) tal como pide
el enunciado — ver la cabecera de cada script para el detalle de qué se
reutiliza de dónde.

## Ejecución

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # una sola vez

julia --project=. scripts/fase1_p1_p2_caracterizacion.jl
julia --project=. scripts/fase2_p3_p4_recorrido_particion.jl
julia --project=. scripts/fase3_p5_p6_p7_optimizacion.jl
julia --project=. scripts/fase4_p8_p9_p10_percolacion.jl
julia --project=. scripts/fase5_p11_rediseno.jl
```

**El orden importa**: cada script empieza verificando la carga contra el
Anexo A y termina con código de salida `0` si todo corrió bien. La Fase 4
(P10) lee `resultados/p6_cortes_minimos_aristas.csv` (lo escribe la Fase 3) y
la Fase 5 (P11) lee `resultados/p10_top10_ficha.csv` (lo escribe la Fase 4) —
por eso no pueden ejecutarse fuera de orden ni de forma aislada.

## Estructura

| Ruta | Contenido |
|---|---|
| `src/ModeloRed.jl` | Utilidades compartidas: modelos de peso (P5), estimación de capacidad (P6.1), puentes/puntos de articulación por Tarjan (P1, reusado en P8), modelo de configuración (P2), NMI/ARI (P4.3) |
| `scripts/fase1_p1_p2_caracterizacion.jl` | P1 (medidas fundamentales) + P2 (modelos nulos, visualización) |
| `scripts/fase2_p3_p4_recorrido_particion.jl` | P3 (BFS/DFS desde cero) + P4 (Louvain desde cero, k-means espectral, NMI/ARI) |
| `scripts/fase3_p5_p6_p7_optimizacion.jl` | P5 (Dijkstra/Floyd-Warshall desde cero) + P6 (Ford-Fulkerson/Edmonds-Karp reutilizados, flujo de costo mínimo) + P7 (p-mediana/p-centro, heurística + JuMP/HiGHS) |
| `scripts/fase4_p8_p9_p10_percolacion.jl` | P8 (percolación de nodos/enlaces) + P9 (cascadas Motter-Lai, SIR) + P10 (ranking de puntos críticos) |
| `scripts/fase5_p11_rediseno.jl` | P11 (propuesta de rediseño, antes/después, comparación con alternativas) |
| `resultados/*.csv` | Tablas que produce cada script — insumo directo del informe |
| `figuras/*.png` | Figuras que produce cada script |

## Decisiones que vale la pena tener presentes al leer el código

- Las funciones "desde cero" que pide el enunciado explícitamente (BFS, DFS,
  Dijkstra, Floyd-Warshall, Louvain, k-means) están implementadas a mano, con
  verificación cruzada contra `Graphs.jl` donde el enunciado lo permite.
  Centralidades, modelos aleatorios (Erdős–Rényi, Barabási–Albert) y
  algoritmos de grafos no mencionados explícitamente usan `Graphs.jl`.
- `Ford-Fulkerson` y `Edmonds-Karp` son los de `codigo_referencia/`,
  incluidos textualmente vía `include()`; como ambos archivos definen un
  `struct RedFlujo` con el mismo nombre, cada uno se aísla en su propio
  módulo (`module FF ... end` / `module EK ... end`) para poder usar los dos
  a la vez sin colisión — ver `fase3_p5_p6_p7_optimizacion.jl`.
- La capacidad de los 181 enlaces no documentados es una ESTIMACIÓN basada en
  reglas (capa, rol, el estándar de 10 Gbps del informe); está documentada
  arista por arista en `resultados/p6_capacidad_estimada.csv`.
