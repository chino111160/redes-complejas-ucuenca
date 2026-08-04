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

# ╔═╡ d2d3c5af-576f-4f19-8528-3694fbd80e5a
using Plots, PlutoUI, Printf

# ╔═╡ 5ba34839-9c66-4d6f-acfb-217e149e5c1b
md"""
# Algoritmo de Edmonds-Karp — Flujo máximo con BFS

**Redes Complejas — Universidad de Cuenca**

Edmonds-Karp es Ford-Fulkerson con una regla adicional: el camino
aumentante se busca **siempre con BFS**, es decir, siempre se aumenta por
el camino con *menos arcos*. Esto garantiza terminación y la cota
``O(V \\cdot E^2)``.

Usa el *slider* para recorrer la ejecución: verás la **onda BFS** crecer
nivel por nivel (derecha), el camino más corto resaltado en naranja, y la
red de flujo actualizándose (izquierda, etiquetas `flujo/capacidad`).
"""

# ╔═╡ 5cbae456-19b4-4875-97fe-3832f827222f
# Cargamos el algoritmo y las funciones de dibujo desde el archivo local
include(joinpath(@__DIR__, "edmonds_karp.jl"));

# ╔═╡ 36510c3f-bb55-4ee9-9048-824fb4780b37
begin
    # Red 1: clásica de CLRS (nodos 1=s, 2=v₁, 3=v₂, 4=v₃, 5=v₄, 6=t)
    C1 = zeros(Int, 6, 6)
    C1[1, 2] = 16; C1[1, 3] = 13
    C1[2, 4] = 12; C1[3, 2] = 4; C1[3, 5] = 14
    C1[4, 3] = 9;  C1[4, 6] = 20
    C1[5, 4] = 7;  C1[5, 6] = 4
    red_clrs = RedFlujo(C1, ["s", "v₁", "v₂", "v₃", "v₄", "t"],
                        [(0.0, 1.0), (1.0, 2.0), (1.0, 0.0),
                         (2.2, 2.0), (2.2, 0.0), (3.2, 1.0)])

    # Red 2: "zigzag" — con DFS podría tomar hasta 2000 iteraciones
    M = 1000
    C2 = zeros(Int, 4, 4)
    C2[1, 2] = M; C2[1, 3] = M; C2[2, 3] = 1; C2[2, 4] = M; C2[3, 4] = M
    red_zigzag = RedFlujo(C2, ["s", "u", "v", "t"],
                          [(0.0, 1.0), (1.2, 2.0), (1.2, 0.0), (2.4, 1.0)])
end;

# ╔═╡ 79fcf205-18ac-48b2-9b73-73b5861dfc4f
md"""
**Red de ejemplo:**
$(@bind eleccion Select(["clrs" => "Red clásica (CLRS)", "zigzag" => "Red zigzag (arco trampa u→v)"]))
"""

# ╔═╡ b7e99ee7-0eed-4026-b1d2-d18ef1cc4570
begin
    red = eleccion == "clrs" ? red_clrs : red_zigzag
    s, t = 1, size(red.C, 1)
    flujo_max, F_final, historia = edmonds_karp(red, s, t; verbose=false)
    frames = _fotogramas(red, s, t, historia)
end;

# ╔═╡ b08160b9-48df-4e23-b0f3-26cd8f1c5e7a
md"""
**Paso de la ejecución:**
$(@bind k Slider(1:length(frames); default=1, show_value=true))
*(mueve el slider; el último paso muestra el corte mínimo)*
"""

# ╔═╡ cf66aa33-5844-4297-afd1-8d43fbbe1efb
md"### $(frames[k].titulo)"

# ╔═╡ 8c22915a-cbdc-4fea-b166-f179b6eaf797
dibujar_fotograma(red, frames[k]; s=s, t=t)

# ╔═╡ d5691a96-6a1f-4897-a8f5-bd95bdf1c68b
md"""
!!! info "Resultado"
    Flujo máximo = **$(flujo_max)**, en **$(length(historia))** iteraciones.
    Longitudes de los caminos aumentantes:
    **$(join([length(p.camino)-1 for p in historia], ", "))** — obsérvese
    que nunca decrecen (lema central de la demostración de complejidad).
"""

# ╔═╡ e3c475de-e3d5-4f9a-891c-0b7771c16b38
md"""
!!! warning "¿Por qué BFS y no cualquier camino?"
    En la **red zigzag**, un Ford-Fulkerson que eligiera mal los caminos
    (alternando por el arco trampa ``u \\to v`` de capacidad 1) enviaría
    solo 1 unidad de flujo por iteración: hasta ``2M = 2000`` iteraciones.
    Edmonds-Karp ignora ese arco porque el camino que lo usa tiene
    longitud 3 > 2, y termina en **2 iteraciones**. El número de
    iteraciones queda acotado por ``O(V \\cdot E)`` *independientemente de
    las capacidades*.
"""

# ╔═╡ 25aaa697-53f8-45a0-ab14-cea9f3289a85
md"""
### Preguntas para discusión en clase

- ¿Por qué la distancia BFS de ``s`` a cada nodo nunca disminuye entre
  iteraciones? ¿Cómo implica esto la cota ``O(V \\cdot E^2)``?
- En la onda BFS, ¿qué significa que un nodo quede gris (no alcanzado)?
  Relación con el corte mínimo del último fotograma.
- Cambia ``M`` en la red zigzag: ¿cambia el número de iteraciones de
  Edmonds-Karp? ¿Y el de un DFS mal elegido?
- Modifica la matriz `C1` y observa cómo cambian los niveles BFS y el
  corte mínimo.
"""

# ╔═╡ Cell order:
# ╟─5ba34839-9c66-4d6f-acfb-217e149e5c1b
# ╠═d2d3c5af-576f-4f19-8528-3694fbd80e5a
# ╠═5cbae456-19b4-4875-97fe-3832f827222f
# ╠═36510c3f-bb55-4ee9-9048-824fb4780b37
# ╟─79fcf205-18ac-48b2-9b73-73b5861dfc4f
# ╠═b7e99ee7-0eed-4026-b1d2-d18ef1cc4570
# ╟─b08160b9-48df-4e23-b0f3-26cd8f1c5e7a
# ╟─cf66aa33-5844-4297-afd1-8d43fbbe1efb
# ╠═8c22915a-cbdc-4fea-b166-f179b6eaf797
# ╟─d5691a96-6a1f-4897-a8f5-bd95bdf1c68b
# ╟─e3c475de-e3d5-4f9a-891c-0b7771c16b38
# ╟─25aaa697-53f8-45a0-ab14-cea9f3289a85
