# ============================================================
# Algoritmo de Edmonds-Karp (flujo máximo) en Julia
# con animación de la onda BFS, camino más corto y corte mínimo
# ============================================================
#
# Edmonds-Karp es la especialización de Ford-Fulkerson en la que el
# camino aumentante SIEMPRE se busca con BFS, es decir, siempre se
# aumenta por el camino con MENOS arcos. Esto garantiza:
#
#   1. Terminación incluso con capacidades irracionales.
#   2. Complejidad O(V·E²): la distancia BFS de s a cada nodo nunca
#      disminuye entre iteraciones, y cada arco puede "saturarse"
#      a lo sumo O(V) veces.
#
# Este archivo contiene:
#   1. Estructuras de datos (RedFlujo, PasoEK)
#   2. BFS con registro de niveles y árbol de exploración
#   3. El algoritmo de Edmonds-Karp (con tabla de iteraciones)
#   4. Corte mínimo (teorema max-flow min-cut)
#   5. Funciones de dibujo (red de flujo, onda BFS, red residual)
#   6. Animación (GIF) y modo interactivo (paso a paso)
#
# Uso típico (ver ejemplo1.jl):
#   include("edmonds_karp.jl")
#   flujo, F, historia = edmonds_karp(red, s, t)
#   animar_edmonds_karp(red, s, t; archivo="ek.gif")
#   edmonds_karp_interactivo(red, s, t)   # [Enter] para avanzar

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
    PasoEK

Registro de una iteración de Edmonds-Karp. Además del camino, el
cuello de botella Δ y el flujo, guarda la información de la BFS que
encontró el camino (para poder animar la "onda" de exploración):

- `nivel[v]`: distancia BFS desde s en la red residual (-1 si v no
  es alcanzable).
- `arbol`: arcos (u,v) del árbol BFS, en el orden en que se
  descubrieron los nodos.

La matriz de flujo `F` es antisimétrica (F[u,v] = -F[v,u]), de modo
que la capacidad residual es r(u,v) = C[u,v] - F[u,v] para cualquier
par de nodos.
"""
struct PasoEK
    camino::Vector{Int}
    Δ::Int
    F::Matrix{Int}
    flujo_total::Int
    nivel::Vector{Int}
    arbol::Vector{Tuple{Int,Int}}
end

# ------------------------------------------------------------
# 2. BFS con niveles (el corazón de Edmonds-Karp)
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
    bfs_niveles(C, F, s, t) -> (camino, nivel, arbol)

BFS sobre la red residual desde `s`, registrando la distancia
(`nivel`) de cada nodo y el árbol de exploración. Como BFS visita
los nodos por capas, el camino devuelto es un camino aumentante de
longitud MÍNIMA (en número de arcos) — la propiedad que define a
Edmonds-Karp.

Nota: para animar la onda completa, esta versión no se detiene al
llegar a `t`; explora todos los nodos alcanzables (una implementación
de producción se detendría en `t`).
"""
function bfs_niveles(C::Matrix{Int}, F::Matrix{Int}, s::Int, t::Int)
    n = size(C, 1)
    nivel = fill(-1, n)
    padre = zeros(Int, n)
    nivel[s] = 0
    padre[s] = s
    cola = [s]
    arbol = Tuple{Int,Int}[]
    while !isempty(cola)
        u = popfirst!(cola)
        for v in 1:n
            # r(u,v) = C - F > 0 ⇒ el arco residual u→v existe
            if nivel[v] == -1 && C[u, v] - F[u, v] > 0
                nivel[v] = nivel[u] + 1
                padre[v] = u
                push!(arbol, (u, v))
                push!(cola, v)
            end
        end
    end
    camino = nivel[t] == -1 ? Int[] : _reconstruir(padre, s, t)
    return camino, nivel, arbol
end

# ------------------------------------------------------------
# 3. Algoritmo de Edmonds-Karp
# ------------------------------------------------------------

"""
    edmonds_karp(red, s, t; verbose=true) -> (flujo_max, F, historia)

Calcula el flujo máximo de `s` a `t` con Edmonds-Karp:

1. Empezar con flujo cero.
2. Mientras BFS encuentre un camino aumentante s→t en la red
   residual (el más corto en número de arcos):
   a. Δ = mínima capacidad residual a lo largo del camino.
   b. Aumentar el flujo en Δ (los arcos inversos cancelan flujo).
3. Al no existir más caminos, el flujo es máximo.

Con `verbose=true` imprime una tabla con la longitud de cada camino
aumentante: obsérvese que las longitudes NUNCA decrecen — este es el
lema central de la demostración de la cota O(V·E²).
"""
function edmonds_karp(red::RedFlujo, s::Int, t::Int; verbose::Bool=true)
    C = red.C
    n = size(C, 1)
    F = zeros(Int, n, n)
    historia = PasoEK[]
    flujo_total = 0

    verbose && @printf("%-5s %-9s %-6s %-7s %s\n",
                       "Iter", "Longitud", "Δ", "Flujo", "Camino")
    while true
        camino, nivel, arbol = bfs_niveles(C, F, s, t)
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

        push!(historia, PasoEK(camino, Δ, copy(F), flujo_total, nivel, arbol))
        verbose && @printf("%-5d %-9d %-6d %-7d %s\n",
                           length(historia), length(camino) - 1, Δ, flujo_total,
                           join(red.nombres[camino], " → "))
    end

    if verbose
        longitudes = [length(p.camino) - 1 for p in historia]
        @printf("\nFlujo máximo: %d (en %d iteraciones)\n",
                flujo_total, length(historia))
        println("Longitudes de los caminos: ", join(longitudes, ", "),
                "  → no decrecientes ✓ (lema de Edmonds-Karp)")
    end
    return flujo_total, F, historia
end

# ------------------------------------------------------------
# 4. Corte mínimo (teorema max-flow min-cut)
# ------------------------------------------------------------

"""
    corte_minimo(C, F, s) -> (S, aristas_corte)

Con el flujo máximo `F`, el conjunto `S` de nodos alcanzables desde
`s` en la red residual define el corte mínimo (S, V∖S); la suma de
las capacidades de sus aristas es exactamente el flujo máximo.
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

# Paleta de colores para los niveles BFS (nivel 0 = fuente)
const PALETA_NIVELES = [:gold, :palegreen, :skyblue, :plum,
                        :lightsalmon, :khaki, :lightpink]

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

"Crea el lienzo base con límites calculados a partir de las posiciones."
function _lienzo(red::RedFlujo, titulo::String)
    xs = [p[1] for p in red.pos]; ys = [p[2] for p in red.pos]
    return plot(; legend=false, axis=false, grid=false, ticks=false,
                aspect_ratio=:equal, title=titulo, titlefontsize=10,
                xlims=(minimum(xs) - 0.45, maximum(xs) + 0.45),
                ylims=(minimum(ys) - 0.55, maximum(ys) + 0.45))
end

"Dibuja los nodos con sus nombres. Los nodos de `S` se pintan dorados."
function _nodos!(plt, red::RedFlujo, s::Int, t::Int, S::Vector{Int})
    for i in 1:length(red.nombres)
        color = i in S ? :gold :
                i == s ? :palegreen :
                i == t ? :lightsalmon : :lightblue
        scatter!(plt, [red.pos[i][1]], [red.pos[i][2]];
                 markersize=17, color=color, markerstrokecolor=:black,
                 markerstrokewidth=1.5)
        annotate!(plt, red.pos[i][1], red.pos[i][2],
                  text(red.nombres[i], 10, :black, :center))
    end
    return plt
end

"""
    dibujar_red(red, F; camino, titulo, S, s, t)

Dibuja la red de flujo con etiquetas "flujo/capacidad".
Arcos con flujo en azul, sin flujo en gris; el `camino` aumentante
en naranja (punteado si usa un arco residual inverso); si `S` no
está vacío, resalta el corte mínimo en púrpura.
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
    dibujar_bfs(red, F; nivel, arbol, hasta, camino, titulo, s, t)

Dibuja la "onda" BFS sobre la red residual:
- Arcos residuales de fondo en gris claro.
- Arcos del árbol BFS descubiertos hasta el nivel `hasta` en verde
  azulado, etiquetados con su capacidad residual.
- Nodos coloreados por su nivel BFS (d = distancia desde s), con la
  distancia anotada debajo; los nodos aún no alcanzados en gris.
- Si se pasa `camino`, se resalta en naranja (el camino más corto).
"""
function dibujar_bfs(red::RedFlujo, F::Matrix{Int};
                     nivel::Vector{Int}, arbol::Vector{Tuple{Int,Int}},
                     hasta::Int, camino::Vector{Int}=Int[],
                     titulo::String="Onda BFS en la red residual",
                     s::Int=0, t::Int=0)
    plt = _lienzo(red, titulo)
    n = size(red.C, 1)
    _off(u, v) = (red.C[v, u] - F[v, u] > 0) ? 0.07 : 0.0

    # Fondo: todos los arcos residuales
    for u in 1:n, v in 1:n
        if red.C[u, v] - F[u, v] > 0
            _flecha!(plt, red.pos[u], red.pos[v];
                     color=:gray85, lw=1.2, offset=_off(u, v))
        end
    end

    # Árbol BFS hasta el nivel `hasta`
    for (u, v) in arbol
        nivel[v] <= hasta || continue
        r = red.C[u, v] - F[u, v]
        _flecha!(plt, red.pos[u], red.pos[v];
                 color=:teal, lw=2.5, offset=_off(u, v),
                 etiqueta="$r", lab_color=:teal)
    end

    # Camino más corto resaltado
    for i in 1:length(camino)-1
        u, v = camino[i], camino[i+1]
        _flecha!(plt, red.pos[u], red.pos[v];
                 color=:orangered, lw=4, offset=_off(u, v))
    end

    # Nodos coloreados por nivel
    for i in 1:n
        alcanzado = nivel[i] >= 0 && nivel[i] <= hasta
        color = alcanzado ? PALETA_NIVELES[mod1(nivel[i] + 1, length(PALETA_NIVELES))] :
                            :gray92
        scatter!(plt, [red.pos[i][1]], [red.pos[i][2]];
                 markersize=17, color=color, markerstrokecolor=:black,
                 markerstrokewidth=1.5)
        annotate!(plt, red.pos[i][1], red.pos[i][2],
                  text(red.nombres[i], 10, :black, :center))
        alcanzado && annotate!(plt, red.pos[i][1], red.pos[i][2] - 0.3,
                               text("d=$(nivel[i])", 8, :teal, :center))
    end
    return plt
end

"""
    dibujar_residual(red, F; titulo, s, t)

Dibuja la red residual: gris sólido = capacidad restante de los
arcos originales; rojo punteado = arcos de retroceso (permiten
cancelar flujo ya enviado).
"""
function dibujar_residual(red::RedFlujo, F::Matrix{Int};
                          titulo::String="Red residual", s::Int=0, t::Int=0)
    plt = _lienzo(red, titulo)
    n = size(red.C, 1)
    for u in 1:n, v in 1:n
        r = red.C[u, v] - F[u, v]
        r > 0 || continue
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
Genera la secuencia de fotogramas lógicos de una ejecución:
por cada iteración, la onda BFS creciendo nivel a nivel, el camino
más corto encontrado y el flujo actualizado; al final, el corte
mínimo. Cada fotograma es un NamedTuple con un campo `tipo`
(:bfs, :camino, :flujo o :corte).
"""
function _fotogramas(red::RedFlujo, s::Int, t::Int, historia::Vector{PasoEK})
    n = size(red.C, 1)
    frames = NamedTuple[]
    F_prev = zeros(Int, n, n)
    sin_nivel = fill(-1, n)
    push!(frames, (tipo=:flujo, F=F_prev, camino=Int[], S=Int[],
                   nivel=sin_nivel, arbol=Tuple{Int,Int}[], hasta=-1,
                   titulo="Red inicial — flujo = 0"))
    for (i, p) in enumerate(historia)
        d_t = p.nivel[t]
        for k in 1:d_t   # la onda BFS crece capa por capa
            push!(frames, (tipo=:bfs, F=F_prev, camino=Int[], S=Int[],
                           nivel=p.nivel, arbol=p.arbol, hasta=k,
                           titulo="Iteración $i — BFS explora el nivel $k"))
        end
        ruta = join(red.nombres[p.camino], " → ")
        push!(frames, (tipo=:camino, F=F_prev, camino=p.camino, S=Int[],
                       nivel=p.nivel, arbol=p.arbol, hasta=d_t,
                       titulo="Iteración $i — camino más corto (long. $d_t): $ruta, Δ = $(p.Δ)"))
        push!(frames, (tipo=:flujo, F=p.F, camino=Int[], S=Int[],
                       nivel=sin_nivel, arbol=Tuple{Int,Int}[], hasta=-1,
                       titulo="Iteración $i — flujo total = $(p.flujo_total)"))
        F_prev = p.F
    end
    S, _ = corte_minimo(red.C, F_prev, s)
    flujo = isempty(historia) ? 0 : historia[end].flujo_total
    push!(frames, (tipo=:corte, F=F_prev, camino=Int[], S=S,
                   nivel=sin_nivel, arbol=Tuple{Int,Int}[], hasta=-1,
                   titulo="Flujo máximo = $flujo = capacidad del corte mínimo (S en dorado)"))
    return frames
end

"""
Dibuja un fotograma: a la izquierda la red de flujo; a la derecha,
la onda BFS (fotogramas :bfs y :camino) o la red residual
(fotogramas :flujo y :corte).
"""
function dibujar_fotograma(red::RedFlujo, fr; s::Int, t::Int)
    izq = dibujar_red(red, fr.F;
                      camino=fr.tipo == :camino ? fr.camino : Int[],
                      titulo=fr.titulo, S=fr.S, s=s, t=t)
    der = if fr.tipo in (:bfs, :camino)
        dibujar_bfs(red, fr.F; nivel=fr.nivel, arbol=fr.arbol, hasta=fr.hasta,
                    camino=fr.tipo == :camino ? fr.camino : Int[], s=s, t=t)
    else
        dibujar_residual(red, fr.F; s=s, t=t)
    end
    return plot(izq, der; layout=(1, 2), size=(1250, 500))
end

"""
    animar_edmonds_karp(red, s, t; archivo="edmonds_karp.gif", fps=0.8)

Ejecuta el algoritmo y guarda un GIF animado: cada iteración muestra
la onda BFS creciendo capa por capa, el camino más corto resaltado y
el flujo actualizado; el último fotograma muestra el corte mínimo.
"""
function animar_edmonds_karp(red::RedFlujo, s::Int, t::Int;
                             archivo::String="edmonds_karp.gif",
                             fps::Real=0.8, verbose::Bool=true)
    _, _, historia = edmonds_karp(red, s, t; verbose=verbose)
    frames = _fotogramas(red, s, t, historia)
    anim = @animate for fr in frames
        dibujar_fotograma(red, fr; s=s, t=t)
    end
    return gif(anim, archivo; fps=fps)
end

"""
    edmonds_karp_interactivo(red, s, t)

Modo interactivo para clase: muestra la ejecución fotograma a
fotograma (incluida la onda BFS) en la ventana de gráficos; presiona
[Enter] en la consola para avanzar al siguiente paso.
"""
function edmonds_karp_interactivo(red::RedFlujo, s::Int, t::Int)
    _, _, historia = edmonds_karp(red, s, t; verbose=false)
    frames = _fotogramas(red, s, t, historia)
    for (k, fr) in enumerate(frames)
        display(dibujar_fotograma(red, fr; s=s, t=t))
        println("[$k/$(length(frames))] ", fr.titulo)
        if k < length(frames)
            print("    [Enter] para continuar... ")
            readline()
        end
    end
    return nothing
end
