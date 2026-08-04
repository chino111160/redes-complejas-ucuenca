### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 91204fbe-970f-4512-86a8-ab170d392c2a
using Plots, PlutoUI, Printf

# ╔═╡ 8997ef38-2f09-49dd-9ed0-77d42db2be01
md"""
# Algoritmo de Ford-Fulkerson — Flujo máximo

**Redes Complejas — Universidad de Cuenca**

Este notebook interactivo permite explorar el algoritmo de Ford-Fulkerson
paso a paso sobre la red clásica de Cormen et al. (CLRS):

1. Elige el método de búsqueda de caminos aumentantes (BFS o DFS).
2. Usa el *slider* para avanzar por las iteraciones del algoritmo.
3. Observa a la izquierda la red de flujo (etiquetas `flujo/capacidad`)
   y a la derecha la **red residual** (los arcos punteados rojos permiten
   *cancelar* flujo ya enviado).
"""

# ╔═╡ 303d959e-dc3f-477d-bb8f-3bd36c90259f
# Cargamos el algoritmo y las funciones de dibujo desde el archivo local
include(joinpath(@__DIR__, "ford_fulkerson.jl"));

# ╔═╡ 1cf6ee5e-ffbe-48c5-863f-42c1ded5e121
begin
    # Red clásica de CLRS: nodos 1=s, 2=v₁, 3=v₂, 4=v₃, 5=v₄, 6=t
    nombres = ["s", "v₁", "v₂", "v₃", "v₄", "t"]
    pos = [(0.0, 1.0), (1.0, 2.0), (1.0, 0.0), (2.2, 2.0), (2.2, 0.0), (3.2, 1.0)]
    C = zeros(Int, 6, 6)
    C[1, 2] = 16; C[1, 3] = 13
    C[2, 4] = 12; C[3, 2] = 4; C[3, 5] = 14
    C[4, 3] = 9;  C[4, 6] = 20
    C[5, 4] = 7;  C[5, 6] = 4
    red = RedFlujo(C, nombres, pos)
    s, t = 1, 6
end;

# ╔═╡ c8c9631a-8069-4109-b628-3313eacf0654
md"""
**Método de búsqueda del camino aumentante:**
$(@bind metodo Select([:bfs => "BFS (Edmonds-Karp)", :dfs => "DFS (Ford-Fulkerson clásico)"]))
"""

# ╔═╡ c288e1d9-564d-4c9b-a5db-b4506adf2fe7
begin
    flujo_max, F_final, historia =
        ford_fulkerson(red, s, t; metodo=metodo, verbose=false)
    frames = _fotogramas(red, s, t, historia)
end;

# ╔═╡ 41cfbca0-3b77-4f4d-b603-721661e72304
md"""
**Paso de la ejecución:**
$(@bind k Slider(1:length(frames); default=1, show_value=true))
*(mueve el slider para avanzar; el último paso muestra el corte mínimo)*
"""

# ╔═╡ b8655b09-b4fa-44d5-b776-285a7548ad19
md"### $(frames[k].titulo)"

# ╔═╡ 9843c077-0055-4417-8885-77d8dc3464fc
dibujar_fotograma(red, frames[k]; s=s, t=t)

# ╔═╡ f214f762-fe6e-477d-9b2b-e9c1add28afe
md"""
!!! info "Teorema max-flow min-cut"
    El flujo máximo encontrado es **$(flujo_max)**, alcanzado en
    **$(length(historia))** iteraciones con el método **$(metodo)**.

    Al terminar, los nodos alcanzables desde ``s`` en la red residual
    forman el conjunto ``S`` (dorado en el último fotograma). Las aristas
    que cruzan de ``S`` a ``V \\setminus S`` forman el **corte mínimo**, y la
    suma de sus capacidades es exactamente igual al flujo máximo:
    ``|f^*| = c(S, V \\setminus S)``.
"""

# ╔═╡ 931c4cd3-70ef-4ad0-81b2-3901aa0085b4
md"""
### Preguntas para discusión en clase

- ¿Por qué son necesarios los arcos *de retroceso* (punteados rojos) en la
  red residual? Busca una iteración donde el camino aumentante cancele
  flujo ya enviado.
- Compara BFS y DFS: ¿cambia el número de iteraciones? ¿Cambia el flujo
  máximo final? ¿Por qué?
- ¿Qué pasaría con la versión DFS si las capacidades fueran números
  irracionales?
- Modifica la matriz `C` (celda de definición de la red) y observa cómo
  cambian el flujo máximo y el corte mínimo.
"""

# ╔═╡ Cell order:
# ╟─8997ef38-2f09-49dd-9ed0-77d42db2be01
# ╠═91204fbe-970f-4512-86a8-ab170d392c2a
# ╠═303d959e-dc3f-477d-bb8f-3bd36c90259f
# ╠═1cf6ee5e-ffbe-48c5-863f-42c1ded5e121
# ╟─c8c9631a-8069-4109-b628-3313eacf0654
# ╠═c288e1d9-564d-4c9b-a5db-b4506adf2fe7
# ╟─41cfbca0-3b77-4f4d-b603-721661e72304
# ╟─b8655b09-b4fa-44d5-b776-285a7548ad19
# ╠═9843c077-0055-4417-8885-77d8dc3464fc
# ╟─f214f762-fe6e-477d-9b2b-e9c1add28afe
# ╟─931c4cd3-70ef-4ad0-81b2-3901aa0085b4
