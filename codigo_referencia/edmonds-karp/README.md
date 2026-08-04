# Edmonds-Karp — Flujo máximo con BFS (con animaciones)

Implementación educativa del algoritmo de Edmonds-Karp en Julia. A
diferencia del proyecto hermano [`../ford-fulkerson`](../ford-fulkerson),
aquí el foco está en lo que distingue a Edmonds-Karp:

- **Animación de la onda BFS**: se ve crecer la exploración nivel por
  nivel sobre la red residual, con cada nodo anotado con su distancia
  `d` desde la fuente.
- **Camino más corto**: el camino aumentante siempre tiene longitud
  mínima; la tabla de iteraciones muestra que las longitudes nunca
  decrecen (lema central de la cota O(V·E²)).
- **Red zigzag**: ejemplo donde una mala elección de caminos (DFS por el
  arco trampa de capacidad 1) requeriría hasta 2000 iteraciones, mientras
  BFS termina en 2.

## Archivos

| Archivo | Descripción |
|---|---|
| `edmonds_karp.jl` | Algoritmo con registro de niveles BFS, corte mínimo y funciones de dibujo/animación |
| `ejemplo1.jl` | Red clásica de CLRS (flujo máximo = 23) y red zigzag: tabla de iteraciones, GIFs y corte mínimo |
| `notebook_pluto.jl` | Notebook Pluto interactivo con selector de red y slider de pasos |

## Uso rápido

```bash
cd optimization/edmonds-karp
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # solo la primera vez
julia --project=. ejemplo1.jl
```

Genera `edmonds_karp_clrs.gif`, `edmonds_karp_zigzag.gif` y
`flujo_maximo_final.png`.

## Modo interactivo (para clase)

```julia
include("ejemplo1.jl")
edmonds_karp_interactivo(red, s, t)          # avanza con [Enter]
edmonds_karp_interactivo(red_zigzag, 1, 4)
```

## Notebook Pluto

```julia
using Pluto            # ]add Pluto si no está instalado
Pluto.run(notebook=joinpath(pwd(), "notebook_pluto.jl"))
```

## Convenciones de la visualización

- Panel izquierdo: red de flujo, etiquetas `flujo/capacidad` (azul = con
  flujo, naranja = camino aumentante, púrpura = corte mínimo).
- Panel derecho durante la búsqueda: **onda BFS** — nodos coloreados por
  nivel con su distancia `d`, árbol BFS en verde azulado (etiquetado con
  capacidades residuales), nodos no alcanzados en gris.
- Panel derecho tras aumentar: red residual (gris = capacidad restante,
  rojo punteado = arcos de retroceso).
- Último fotograma: nodos de `S` en dorado y aristas del corte mínimo en
  púrpura — teorema max-flow min-cut.
