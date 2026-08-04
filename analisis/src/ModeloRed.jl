#=============================================================
  ModeloRed.jl -- Utilidades compartidas entre las cinco fases
  del proyecto integrador (Redes Complejas 1217, UCuenca)

  Se incluye desde cada script de analisis/scripts/ con:
      include(joinpath(@__DIR__, "..", "src", "ModeloRed.jl"))

  Reúne lo que varias fases necesitan en común, para no duplicarlo:
    - Modelos de peso de P5 (w_saltos, w_latencia, w_carga)
    - Estimación de capacidad de P6.1 (compartida con P5, que la
      necesita para w_latencia y w_carga)
    - Puntos de articulación y puentes (Tarjan, DFS con low-link),
      usados en P1 y reutilizados en P8 (ataque dirigido a puentes)
    - Modelo de configuración con secuencia de grados exacta (P2):
      Graphs.jl no garantiza un grafo simple, así que se implementa
      la reconexión de "stubs" a mano (Newman, 2010, cap. 13)
    - NMI y ARI (P4.3), fórmulas de Newman (2010) y Hubert & Arabie
      (1985) respectivamente
    - Pequeñas utilidades de E/S (tablas top-k, guardar CSV creando
      el directorio si hace falta)
=============================================================#

using Graphs
using DataFrames
using CSV
using Random
using Statistics

# ---------------------------------------------------------------
# Modelos de peso -- Problema P5
# ---------------------------------------------------------------

"""
    w_saltos(u, v) -> Float64

Modelo de peso "número de saltos": todas las aristas valen 1. Con este
peso, Dijkstra/Floyd-Warshall reproducen BFS (camino con menos enlaces).
"""
w_saltos(u, v) = 1.0

"""
    w_latencia(c; α, β) -> Float64

Modelo de peso "latencia estimada": w = α + β / c, con `c` la capacidad
del enlace (Mbps). A mayor capacidad, menor latencia de transmisión
aportada por ese enlace; α representa una latencia fija (propagación +
procesamiento) común a todos los enlaces. Los valores de α y β deben
fijarse y justificarse en el informe (P5, nota "Antes de aplicar...").
"""
w_latencia(c::Real; α::Real, β::Real) = α + β / c

"""
    w_carga(b, c) -> Float64

Modelo de peso "carga relativa": w = b / c, la fracción de la capacidad
que ya está en uso (tráfico medido `b` sobre capacidad `c`). Es siempre
no negativo si `b, c ≥ 0`; puede superar 1 si el tráfico medido excede la
capacidad estimada (indicio de que la estimación de capacidad fue baja).
"""
w_carga(b::Real, c::Real) = b / c

# ---------------------------------------------------------------
# Estimación de capacidad -- Problema P6.1 (compartida con P5)
# ---------------------------------------------------------------

"""
    estimar_capacidad(red) -> (cap::Vector{Float64}, supuesto::Vector{String})

Completa `c(u,v)` para las 209 aristas de la red. Para las 28 aristas del
diagrama MPLS se usa `capacidad_mbps` tal cual viene documentada. Para las
181 restantes se estima a partir de la capa de sus extremos y del rol del
enlace, con las siguientes reglas (documentadas también en el vector
`supuesto`, una por arista, para citarlas en el informe):

  1. Documentada  -> se usa `capacidad_mbps` sin modificar.
  2. Troncal core-agregación -> 10 Gbps, tal como declara el informe
     técnico (capítulo 2.1: "los enlaces troncales entre capas operan a
     10 Gbps").
  3. `rol` ∈ {wan, inferido} (enlaces hacia la nube MPLS o inferidos como
     su interconexión) -> 10 Gbps, por analogía con los enlaces MPLS que
     sí están documentados a esa velocidad.
  4. Agregación-acceso (o core-acceso) -> 1 Gbps, uplink Ethernet
     estándar de un switch de acceso/sub-laboratorio; no está declarado
     por el informe, así que es la estimación más discutible del
     conjunto y se señala como tal.
  5. Cualquier otro caso (p. ej. firewalls de interconexión) -> 1 Gbps
     por defecto, mismo criterio que 4.

En los casos 4 y 5 se anota además si el enlace es `respaldo` o
`secundario` (miembro de un par redundante), porque eso no cambia la
capacidad nominal del enlace individual pero sí es relevante para P1/P6
al leer la tabla de supuestos.
"""
function estimar_capacidad(red)
    n = nrow(red.aristas)
    cap = Vector{Float64}(undef, n)
    supuesto = Vector{String}(undef, n)
    for (i, r) in enumerate(eachrow(red.aristas))
        if !ismissing(r.capacidad_mbps)
            cap[i] = Float64(r.capacidad_mbps)
            supuesto[i] = "documentada (diagrama MPLS)"
            continue
        end
        cu, cv = red.capa[r.u], red.capa[r.v]
        capas = Set((cu, cv))
        if capas == Set(("core", "agregacion"))
            cap[i] = 10_000.0
            supuesto[i] = "troncal core-agregacion, 10 Gbps (informe, cap. 2.1)"
        elseif r.rol in ("wan", "inferido")
            cap[i] = 10_000.0
            supuesto[i] = "wan/inferido hacia la nube MPLS, asumido 10 Gbps por analogía con los enlaces MPLS documentados"
        elseif capas == Set(("agregacion", "acceso")) || capas == Set(("core", "acceso"))
            cap[i] = 1_000.0
            supuesto[i] = "agregacion-acceso, asumido 1 Gbps (uplink Ethernet estandar, no documentado por el informe)"
        else
            cap[i] = 1_000.0
            supuesto[i] = "capa no troncal, asumido 1 Gbps por defecto"
        end
        if r.rol in ("respaldo", "secundario")
            supuesto[i] = supuesto[i] * " [miembro de par redundante]"
        end
    end
    return cap, supuesto
end

# ---------------------------------------------------------------
# Puntos de articulación y puentes -- Tarjan (low-link DFS)
# ---------------------------------------------------------------

"""
    puentes_y_articulacion(g::SimpleGraph) -> (articulacion::Vector{Int}, puentes::Vector{Tuple{Int,Int}})

Algoritmo de Tarjan (DFS con tiempos de descubrimiento `disc` y
"low-link" `low`) para hallar puntos de articulación y puentes en un
grafo simple no dirigido en O(V+E). Un nodo `u` (no raíz) es de
articulación si tiene un hijo `v` con `low[v] >= disc[u]`; la raíz lo es
si tiene más de un hijo en el árbol DFS. Una arista `(u,v)` (con `v`
hijo de `u`) es puente si `low[v] > disc[u]` (no hay atajo hacia `u` o
más arriba desde el subárbol de `v`).

Se implementa a mano porque el enunciado pide contabilizar puentes y
puntos de articulación por campus y por capa (P1.5) y reutilizarlos como
blanco de ataque dirigido en P8; tenerlos como función propia evita
depender de que `Graphs.jl` exponga ambos en la misma versión.
"""
function puentes_y_articulacion(g::SimpleGraph)
    n = nv(g)
    visitado = falses(n)
    disc = zeros(Int, n)
    low = zeros(Int, n)
    padre = zeros(Int, n)
    es_articulacion = falses(n)
    puentes = Tuple{Int,Int}[]
    tiempo = Ref(0)

    function dfs(raiz::Int)
        pila = [(raiz, 1)]  # (nodo, índice del próximo vecino a explorar), DFS iterativo
        vecinos_cache = Dict{Int,Vector{Int}}()
        hijos_raiz = 0
        visitado[raiz] = true
        tiempo[] += 1
        disc[raiz] = low[raiz] = tiempo[]
        while !isempty(pila)
            u, k = pila[end]
            vs = get!(vecinos_cache, u) do
                collect(neighbors(g, u))
            end
            if k <= length(vs)
                pila[end] = (u, k + 1)
                v = vs[k]
                if !visitado[v]
                    visitado[v] = true
                    padre[v] = u
                    tiempo[] += 1
                    disc[v] = low[v] = tiempo[]
                    u == raiz && (hijos_raiz += 1)
                    push!(pila, (v, 1))
                elseif v != padre[u]
                    low[u] = min(low[u], disc[v])
                end
            else
                pop!(pila)
                if !isempty(pila)
                    p = pila[end][1]
                    low[p] = min(low[p], low[u])
                    if low[u] >= disc[p] && p != raiz
                        es_articulacion[p] = true
                    end
                    if low[u] > disc[p]
                        push!(puentes, (p, u))
                    end
                end
            end
        end
        hijos_raiz > 1 && (es_articulacion[raiz] = true)
    end

    for s in 1:n
        visitado[s] || dfs(s)
    end
    return findall(es_articulacion), puentes
end

# ---------------------------------------------------------------
# Modelo de configuración con secuencia de grados exacta -- P2
# ---------------------------------------------------------------

"""
    modelo_configuracion(secuencia::Vector{Int}) -> SimpleGraph

Construye una realización del modelo de configuración (Newman, 2010,
cap. 13) para la secuencia de grados dada, mediante emparejamiento
aleatorio de "medios enlaces" (*stub matching*). Se descartan los pares
que producirían un self-loop o una arista repetida (para mantener el
grafo simple, como exige el resto del análisis); los pares descartados
se reintentan en unas pocas rondas adicionales con el resto de stubs. Es
normal terminar con unas pocas aristas menos que `sum(secuencia)/2`; el
script que llama a esta función reporta esa diferencia.
"""
function modelo_configuracion(secuencia::Vector{Int})
    n = length(secuencia)
    g = SimpleGraph(n)
    stubs = Int[]
    for (i, d) in enumerate(secuencia)
        append!(stubs, fill(i, d))
    end
    shuffle!(stubs)
    restantes = Int[]
    for i in 1:2:length(stubs)-1
        u, v = stubs[i], stubs[i+1]
        if u == v || has_edge(g, u, v)
            push!(restantes, u, v)
        else
            add_edge!(g, u, v)
        end
    end
    for _ in 1:5
        length(restantes) < 2 && break
        shuffle!(restantes)
        nuevos = Int[]
        for i in 1:2:length(restantes)-1
            u, v = restantes[i], restantes[i+1]
            if u == v || has_edge(g, u, v)
                push!(nuevos, u, v)
            else
                add_edge!(g, u, v)
            end
        end
        restantes = nuevos
    end
    return g
end

# ---------------------------------------------------------------
# NMI y ARI -- Problema P4.3
# ---------------------------------------------------------------

"Tabla de contingencia (conteo) entre dos particiones dadas como vectores de enteros de igual longitud."
function tabla_contingencia(a::Vector{Int}, b::Vector{Int})
    ua, ub = sort(unique(a)), sort(unique(b))
    ia = Dict(v => i for (i, v) in enumerate(ua))
    ib = Dict(v => i for (i, v) in enumerate(ub))
    T = zeros(Int, length(ua), length(ub))
    for (x, y) in zip(a, b)
        T[ia[x], ib[y]] += 1
    end
    return T, ua, ub
end

"""
    nmi(a, b) -> Float64

Información mutua normalizada entre dos particiones (Newman, 2010, ec.
11.13): `I(a,b) / ((H(a)+H(b))/2)`, en `[0, 1]`. `1` indica particiones
idénticas salvo relabeling; `0`, independencia estadística.
"""
function nmi(a::Vector{Int}, b::Vector{Int})
    T, _, _ = tabla_contingencia(a, b)
    n = sum(T)
    Pa = vec(sum(T, dims=2)) ./ n
    Pb = vec(sum(T, dims=1)) ./ n
    Ha = -sum(p * log(p) for p in Pa if p > 0)
    Hb = -sum(p * log(p) for p in Pb if p > 0)
    I = 0.0
    for i in axes(T, 1), j in axes(T, 2)
        pij = T[i, j] / n
        pij > 0 && (I += pij * log(pij / (Pa[i] * Pb[j])))
    end
    denom = (Ha + Hb) / 2
    return denom == 0 ? 1.0 : I / denom
end

"""
    ari(a, b) -> Float64

Índice de Rand ajustado (Hubert & Arabie, 1985): corrige el índice de
Rand por el acuerdo esperado al azar dado el tamaño de los grupos. Vale
`1` para particiones idénticas y ≈`0` para particiones tan parecidas
como dos particiones aleatorias independientes; puede ser negativo.
"""
function ari(a::Vector{Int}, b::Vector{Int})
    T, _, _ = tabla_contingencia(a, b)
    n = sum(T)
    comb2(x) = x * (x - 1) / 2
    sum_ij = sum(comb2, T)
    ai = vec(sum(T, dims=2))
    bj = vec(sum(T, dims=1))
    sum_a, sum_b = sum(comb2, ai), sum(comb2, bj)
    total = comb2(n)
    esperado = total == 0 ? 0.0 : sum_a * sum_b / total
    maximo = (sum_a + sum_b) / 2
    denom = maximo - esperado
    return denom == 0 ? 1.0 : (sum_ij - esperado) / denom
end

# ---------------------------------------------------------------
# Utilidades de E/S
# ---------------------------------------------------------------

"Escribe `df` en `ruta` como CSV, creando el directorio contenedor si hace falta."
function guardar_csv(df::DataFrame, ruta::AbstractString)
    mkpath(dirname(ruta))
    CSV.write(ruta, df)
    return ruta
end

"""
    tabla_topk(ids, valores; k=10, nombre_valor="valor") -> DataFrame

Tabla de las `k` mayores entradas de `valores`, con su identificador y
su puesto (1 = mayor). Se usa para los top-10 de centralidades, P1.3.
"""
function tabla_topk(ids::Vector{String}, valores::Vector{<:Real}; k::Int=10, nombre_valor::String="valor")
    orden = sortperm(valores; rev=true)[1:min(k, length(valores))]
    DataFrame(puesto=1:length(orden), id=ids[orden], Symbol(nombre_valor) => valores[orden])
end
