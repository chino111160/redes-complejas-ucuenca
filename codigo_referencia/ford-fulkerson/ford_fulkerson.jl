# ============================================================
# Algoritmo de Ford-Fulkerson (flujo máximo) en Julia
# con visualización paso a paso, animaciones y corte mínimo
# ============================================================
#
# Este archivo contiene:
#   1. Estructuras de datos (RedFlujo, Paso)
#   2. Búsqueda de caminos aumentantes (BFS y DFS)
#   3. El algoritmo de Ford-Fulkerson
#   4. Cálculo del corte mínimo (teorema max-flow min-cut)
#   5. Funciones de dibujo (red de flujo y red residual)
#   6. Animación (GIF) y modo interactivo (paso a paso)
#
# Uso típico (ver ejemplo1.jl):
#   include("ford_fulkerson.jl")
#   flujo, F, historia = ford_fulkerson(red, s, t)
#   animar_ford_fulkerson(red, s, t; archivo="ff.gif")
#   ford_fulkerson_interactivo(red, s, t)   # [Enter] para avanzar

using Plots
using Printf

# ------------------------------------------------------------
# 1. Estructuras de datos
# ------------------------------------------------------------

"""
    RedFlujo(C, nombres, pos)

Red de flujo dirigida con `n` nodos.

- `C[u,v]`: capacidad del arco u→v (0 si el arco no existe).
- `nombres`: etiquetas de los nodos (para graficar).
- `pos`: posiciones fijas (x, y) de cada nodo en el plano.
"""
struct RedFlujo
    C::Matrix{Int}
    nombres::Vector{String}
    pos::Vector{Tuple{Float64,Float64}}
end

"""
    Paso

Registro de una iteración del algoritmo: el camino aumentante
encontrado, el cuello de botella Δ, una copia de la matriz de
flujo *después* de aumentar, y el flujo total acumulado.

La matriz de flujo `F` es antisimétrica (F[u,v] = -F[v,u]), la
convención clásica que hace que la capacidad residual sea
simplemente r(u,v) = C[u,v] - F[u,v] para *cualquier* par de nodos.
"""
struct Paso
    camino::Vector{Int}
    Δ::Int
    F::Matrix{Int}
    flujo_total::Int
end

# ------------------------------------------------------------
# 2. Búsqueda de caminos aumentantes en la red residual
# ------------------------------------------------------------

"Reconstruye el camino s→t a partir del vector de padres."
function _reconstruir(padre::Vector{Int}, s::Int, t::Int)
    camino = [t]
    while camino[1] != s
        pushfirst!(camino, padre[camino[1]])
    end
    return camino
end

"""
    buscar_camino_bfs(C, F, s, t) -> Vector{Int}

Busca un camino aumentante s→t en la red residual usando BFS
(variante de Edmonds-Karp: siempre encuentra el camino aumentante
con MENOS arcos, lo que garantiza complejidad O(V·E²)).
Devuelve el camino como lista de nodos, o `Int[]` si no existe.
"""
function buscar_camino_bfs(C::Matrix{Int}, F::Matrix{Int}, s::Int, t::Int)
    n = size(C, 1)
    padre = zeros(Int, n)
    padre[s] = s
    cola = [s]
    while !isempty(cola)
        u = popfirst!(cola)
        for v in 1:n
            # r(u,v) = C - F > 0 ⇒ el arco residual u→v existe
            if padre[v] == 0 && C[u, v] - F[u, v] > 0
                padre[v] = u
                v == t && return _reconstruir(padre, s, t)
                push!(cola, v)
            end
        end
    end
    return Int[]
end

"""
    buscar_camino_dfs(C, F, s, t) -> Vector{Int}

Versión "clásica" de Ford-Fulkerson: busca cualquier camino
aumentante con DFS. Funciona, pero puede necesitar muchas más
iteraciones que BFS (con capacidades irracionales incluso podría
no terminar — buen tema de discusión en clase).
"""
function buscar_camino_dfs(C::Matrix{Int}, F::Matrix{Int}, s::Int, t::Int)
    n = size(C, 1)
    padre = zeros(Int, n)
    padre[s] = s
    pila = [s]
    while !isempty(pila)
        u = pop!(pila)
        for v in n:-1:1   # orden inverso para explorar primero los índices bajos
            if padre[v] == 0 && C[u, v] - F[u, v] > 0
                padre[v] = u
                v == t && return _reconstruir(padre, s, t)
                push!(pila, v)
            end
        end
    end
    return Int[]
end

# ------------------------------------------------------------
# 3. Algoritmo de Ford-Fulkerson
# ------------------------------------------------------------

"""
    ford_fulkerson(red, s, t; metodo=:bfs, verbose=true)
        -> (flujo_max, F, historia)

Calcula el flujo máximo de `s` a `t` en la red `red`.

Idea del algoritmo:
1. Empezar con flujo cero.
2. Mientras exista un camino aumentante s→t en la red residual:
   a. Encontrar el cuello de botella Δ = mínima capacidad residual
      a lo largo del camino.
   b. Aumentar el flujo en Δ a lo largo del camino (los arcos
      inversos "cancelan" flujo ya enviado).
3. Al no existir más caminos, el flujo es máximo (teorema
   max-flow min-cut).

- `metodo = :bfs` usa Edmonds-Karp; `:dfs` usa la versión clásica.
- `historia` es un vector de `Paso` con cada iteración, que las
  funciones de animación usan para reconstruir la ejecución.
"""
function ford_fulkerson(red::RedFlujo, s::Int, t::Int;
                        metodo::Symbol=:bfs, verbose::Bool=true)
    C = red.C
    n = size(C, 1)
    F = zeros(Int, n, n)
    historia = Paso[]
    flujo_total = 0
    buscar = metodo == :bfs ? buscar_camino_bfs : buscar_camino_dfs

    while true
        camino = buscar(C, F, s, t)
        isempty(camino) && break   # no hay más caminos aumentantes

        # Cuello de botella: mínima capacidad residual del camino
        Δ = minimum(C[camino[i], camino[i+1]] - F[camino[i], camino[i+1]]
                    for i in 1:length(camino)-1)

        # Aumentar el flujo (F antisimétrica: el arco inverso cancela)
        for i in 1:length(camino)-1
            u, v = camino[i], camino[i+1]
            F[u, v] += Δ
            F[v, u] -= Δ
        end
        flujo_total += Δ

        push!(historia, Paso(camino, Δ, copy(F), flujo_total))
        if verbose
            @printf("Iteración %d: camino %s,  Δ = %d,  flujo total = %d\n",
                    length(historia), join(red.nombres[camino], " → "),
                    Δ, flujo_total)
        end
    end

    verbose && @printf("\nFlujo máximo: %d (en %d iteraciones)\n",
                       flujo_total, length(historia))
    return flujo_total, F, historia
end

# ------------------------------------------------------------
# 4. Corte mínimo (teorema max-flow min-cut)
# ------------------------------------------------------------

"""
    corte_minimo(C, F, s) -> (S, aristas_corte)

Con el flujo máximo `F`, el conjunto `S` de nodos alcanzables desde
`s` en la red residual define el corte mínimo (S, V∖S). Las aristas
del corte son los arcos originales que van de S hacia V∖S; la suma
de sus capacidades es exactamente el flujo máximo.
"""
function corte_minimo(C::Matrix{Int}, F::Matrix{Int}, s::Int)
    n = size(C, 1)
    visitado = falses(n)
    visitado[s] = true
    cola = [s]
    while !isempty(cola)
        u = popfirst!(cola)
        for v in 1:n
            if !visitado[v] && C[u, v] - F[u, v] > 0
                visitado[v] = true
                push!(cola, v)
            end
        end
    end
    S = findall(visitado)
    aristas = [(u, v) for u in S, v in findall(.!visitado) if C[u, v] > 0]
    return S, vec(aristas)
end

# ------------------------------------------------------------
# 5. Funciones de dibujo
# ------------------------------------------------------------

"""
Dibuja una flecha de `p1` a `p2`, acortada para no tapar los nodos,
con desplazamiento perpendicular `offset` (para arcos antiparalelos)
y una etiqueta opcional junto al punto medio.
"""
function _flecha!(plt, p1, p2; color=:gray55, lw=1.5, estilo=:solid,
                  etiqueta="", lab_color=:gray30, offset=0.0, radio=0.17)
    dx, dy = p2[1] - p1[1], p2[2] - p1[2]
    L = hypot(dx, dy)
    ux, uy = dx / L, dy / L          # dirección unitaria
    nx, ny = -uy, ux                 # perpendicular (izquierda)
    ax, ay = p1[1] + ux * radio + nx * offset, p1[2] + uy * radio + ny * offset
    bx, by = p2[1] - ux * radio + nx * offset, p2[2] - uy * radio + ny * offset
    plot!(plt, [ax, bx], [ay, by];
          color=color, lw=lw, linestyle=estilo, arrow=true)
    if etiqueta != ""
        mx = (ax + bx) / 2 + nx * 0.14
        my = (ay + by) / 2 + ny * 0.14
        annotate!(plt, mx, my, text(etiqueta, 9, lab_color, :center))
    end
    return plt
end

"Dibuja los nodos con sus nombres. Los nodos de `S` se pintan dorados."
function _nodos!(plt, red::RedFlujo, s::Int, t::Int, S::Vector{Int})
    n = length(red.nombres)
    for i in 1:n
        color = i in S       ? :gold :
                i == s       ? :palegreen :
                i == t       ? :lightsalmon : :lightblue
        scatter!(plt, [red.pos[i][1]], [red.pos[i][2]];
                 markersize=17, color=color, markerstrokecolor=:black,
                 markerstrokewidth=1.5)
        annotate!(plt, red.pos[i][1], red.pos[i][2],
                  text(red.nombres[i], 10, :black, :center))
    end
    return plt
end

"Crea el lienzo base con límites calculados a partir de las posiciones."
function _lienzo(red::RedFlujo, titulo::String)
    xs = [p[1] for p in red.pos]; ys = [p[2] for p in red.pos]
    return plot(; legend=false, axis=false, grid=false, ticks=false,
                aspect_ratio=:equal, title=titulo, titlefontsize=10,
                xlims=(minimum(xs) - 0.45, maximum(xs) + 0.45),
                ylims=(minimum(ys) - 0.45, maximum(ys) + 0.45))
end

"""
    dibujar_red(red, F; camino, titulo, S, s, t)

Dibuja la red de flujo con etiquetas "flujo/capacidad".
- Arcos con flujo positivo en azul; sin flujo en gris.
- El `camino` aumentante se resalta en naranja (si el camino usa un
  arco residual inverso, se dibuja punteado).
- Si `S` no está vacío, se resaltan los nodos de S y las aristas del
  corte mínimo en púrpura.
"""
function dibujar_red(red::RedFlujo, F::Matrix{Int};
                     camino::Vector{Int}=Int[], titulo::String="",
                     S::Vector{Int}=Int[], s::Int=0, t::Int=0)
    plt = _lienzo(red, titulo)
    n = size(red.C, 1)
    arcos_camino = Set(zip(camino[1:max(end - 1, 0)], camino[2:end]))

    for u in 1:n, v in 1:n
        red.C[u, v] > 0 || continue
        off = red.C[v, u] > 0 ? 0.07 : 0.0   # separar arcos antiparalelos
        f = max(F[u, v], 0)
        en_camino = (u, v) in arcos_camino
        en_corte  = !isempty(S) && (u in S) && !(v in S)
        color = en_camino ? :orangered :
                en_corte  ? :purple :
                f > 0     ? :steelblue : :gray60
        lw = (en_camino || en_corte) ? 4 : (f > 0 ? 2.5 : 1.5)
        _flecha!(plt, red.pos[u], red.pos[v];
                 color=color, lw=lw, offset=off,
                 etiqueta="$f/$(red.C[u, v])",
                 lab_color=f > 0 ? :steelblue : :gray45)
    end

    # Arcos residuales inversos usados por el camino (cancelan flujo)
    for (u, v) in arcos_camino
        if red.C[u, v] == 0
            _flecha!(plt, red.pos[u], red.pos[v];
                     color=:orangered, lw=3, estilo=:dash, offset=0.09)
        end
    end

    return _nodos!(plt, red, s, t, S)
end

"""
    dibujar_residual(red, F; titulo, s, t)

Dibuja la red residual: cada arco con capacidad residual r > 0.
Los arcos "de avance" (parte de la capacidad original sin usar) van
en gris sólido; los arcos "de retroceso" (que permiten cancelar
flujo ya enviado) van punteados en rojo.
"""
function dibujar_residual(red::RedFlujo, F::Matrix{Int};
                          titulo::String="Red residual", s::Int=0, t::Int=0)
    plt = _lienzo(red, titulo)
    n = size(red.C, 1)
    for u in 1:n, v in 1:n
        r = red.C[u, v] - F[u, v]
        r > 0 || continue
        # ¿existe también el arco residual contrario? → separar
        off = (red.C[v, u] - F[v, u] > 0) ? 0.07 : 0.0
        if red.C[u, v] > 0
            _flecha!(plt, red.pos[u], red.pos[v];
                     color=:gray50, lw=1.8, offset=off,
                     etiqueta="$r", lab_color=:gray35)
        else
            _flecha!(plt, red.pos[u], red.pos[v];
                     color=:indianred, lw=1.8, estilo=:dash, offset=off,
                     etiqueta="$r", lab_color=:indianred)
        end
    end
    return _nodos!(plt, red, s, t, Int[])
end

# ------------------------------------------------------------
# 6. Animación y modo interactivo
# ------------------------------------------------------------

"""
Genera la secuencia de "fotogramas" lógicos de una ejecución:
estado inicial → (camino resaltado, flujo actualizado) por cada
iteración → corte mínimo final.
"""
function _fotogramas(red::RedFlujo, s::Int, t::Int, historia::Vector{Paso})
    n = size(red.C, 1)
    frames = NamedTuple[]
    F_prev = zeros(Int, n, n)
    push!(frames, (F=F_prev, camino=Int[], S=Int[],
                   titulo="Red inicial — flujo = 0"))
    for (i, p) in enumerate(historia)
        ruta = join(red.nombres[p.camino], " → ")
        push!(frames, (F=F_prev, camino=p.camino, S=Int[],
                       titulo="Iteración $i: camino $ruta   (Δ = $(p.Δ))"))
        push!(frames, (F=p.F, camino=Int[], S=Int[],
                       titulo="Iteración $i: flujo total = $(p.flujo_total)"))
        F_prev = p.F
    end
    S, _ = corte_minimo(red.C, F_prev, s)
    flujo = isempty(historia) ? 0 : historia[end].flujo_total
    push!(frames, (F=F_prev, camino=Int[], S=S,
                   titulo="Flujo máximo = $flujo = capacidad del corte mínimo (S en dorado)"))
    return frames
end

"Dibuja un fotograma: red de flujo y, opcionalmente, la red residual al lado."
function dibujar_fotograma(red::RedFlujo, fr; s::Int, t::Int, residual::Bool=true)
    p1 = dibujar_red(red, fr.F; camino=fr.camino, titulo=fr.titulo,
                     S=fr.S, s=s, t=t)
    residual || return plot(p1, size=(680, 500))
    p2 = dibujar_residual(red, fr.F; s=s, t=t)
    return plot(p1, p2; layout=(1, 2), size=(1250, 500))
end

"""
    animar_ford_fulkerson(red, s, t; metodo=:bfs, archivo="ford_fulkerson.gif",
                          fps=0.6, residual=true)

Ejecuta el algoritmo y guarda un GIF animado: cada iteración muestra
primero el camino aumentante resaltado y luego el flujo actualizado,
junto con la red residual. El último fotograma muestra el corte mínimo.
"""
function animar_ford_fulkerson(red::RedFlujo, s::Int, t::Int;
                               metodo::Symbol=:bfs,
                               archivo::String="ford_fulkerson.gif",
                               fps::Real=0.6, residual::Bool=true,
                               verbose::Bool=true)
    _, _, historia = ford_fulkerson(red, s, t; metodo=metodo, verbose=verbose)
    frames = _fotogramas(red, s, t, historia)
    anim = @animate for fr in frames
        dibujar_fotograma(red, fr; s=s, t=t, residual=residual)
    end
    return gif(anim, archivo; fps=fps)
end

"""
    ford_fulkerson_interactivo(red, s, t; metodo=:bfs, residual=true)

Modo interactivo para clase: muestra la ejecución fotograma a
fotograma en la ventana de gráficos; presiona [Enter] en la consola
para avanzar al siguiente paso.
"""
function ford_fulkerson_interactivo(red::RedFlujo, s::Int, t::Int;
                                    metodo::Symbol=:bfs, residual::Bool=true)
    _, _, historia = ford_fulkerson(red, s, t; metodo=metodo, verbose=false)
    frames = _fotogramas(red, s, t, historia)
    for (k, fr) in enumerate(frames)
        display(dibujar_fotograma(red, fr; s=s, t=t, residual=residual))
        println("[$k/$(length(frames))] ", fr.titulo)
        if k < length(frames)
            print("    [Enter] para continuar... ")
            readline()
        end
    end
    return nothing
end
