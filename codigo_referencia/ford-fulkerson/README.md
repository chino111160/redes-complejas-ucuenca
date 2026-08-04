# Ford-Fulkerson — Flujo máximo (con animaciones)

Implementación educativa del algoritmo de Ford-Fulkerson en Julia, con
visualización paso a paso de la red de flujo y la red residual, animación
en GIF, modo interactivo en el REPL y notebook Pluto con sliders.

## Archivos

| Archivo | Descripción |
|---|---|
| `ford_fulkerson.jl` | Algoritmo (BFS/Edmonds-Karp y DFS clásico), corte mínimo y funciones de dibujo/animación |
| `ejemplo1.jl` | Red clásica de CLRS (flujo máximo = 23): traza en consola, GIFs y corte mínimo |
| `notebook_pluto.jl` | Notebook Pluto interactivo con selector de método y slider de iteraciones |

## Uso rápido

```bash
cd algoritmos/ford-fulkerson
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # solo la primera vez
julia --project=. ejemplo1.jl
```

Esto imprime la traza de iteraciones y genera `ford_fulkerson_bfs.gif`,
`ford_fulkerson_dfs.gif` y `flujo_maximo_final.png`.

## Modo interactivo (para clase)

En el REPL o VSCode, desde esta carpeta:

```julia
include("ejemplo1.jl")
ford_fulkerson_interactivo(red, s, t)              # avanza con [Enter]
ford_fulkerson_interactivo(red, s, t; metodo=:dfs) # versión clásica DFS
```

## Notebook Pluto

```julia
using Pluto            # ]add Pluto si no está instalado
Pluto.run(notebook=joinpath(pwd(), "notebook_pluto.jl"))
```

Permite cambiar el método (BFS/DFS) y recorrer las iteraciones con un
slider, mostrando la red de flujo y la red residual lado a lado.

## Convenciones de la visualización

- Etiquetas de arcos: `flujo/capacidad`. Azul = con flujo, gris = sin flujo.
- Naranja: camino aumentante de la iteración actual (punteado si usa un
  arco residual inverso, es decir, cancela flujo).
- Red residual: gris sólido = capacidad restante; rojo punteado = arco de
  retroceso (permite deshacer flujo).
- Último fotograma: nodos de `S` en dorado y aristas del corte mínimo en
  púrpura — ilustración del teorema max-flow min-cut.
