#=============================================================
  fase3_p5_p6_p7_optimizacion.jl -- Fase 3: Optimización en Redes
  Proyecto integrador de Redes Complejas (1217) -- Universidad de Cuenca

  Resuelve P5 (caminos más cortos), P6 (flujo máximo y corte mínimo) y
  P7 (localización de instalaciones) del enunciado
  (proyecto_red_ucuenca.pdf, cap. 5).

  Uso:
      julia --project=analisis analisis/scripts/fase3_p5_p6_p7_optimizacion.jl
=============================================================#

using Graphs
using DataFrames
using CSV
using Statistics
using StatsBase
using Printf
using Random
using LinearAlgebra
using JuMP
using HiGHS

const DIR_ESTE_SCRIPT = @__DIR__
const DIR_ANALISIS = dirname(DIR_ESTE_SCRIPT)
const DIR_RAIZ = dirname(DIR_ANALISIS)

include(joinpath(DIR_RAIZ, "codigo_base", "cargar_red.jl"))
include(joinpath(DIR_ANALISIS, "src", "ModeloRed.jl"))

const DIR_RESULTADOS = joinpath(DIR_ANALISIS, "resultados")
const DIR_FIGURAS = joinpath(DIR_ANALISIS, "figuras")
mkpath(DIR_RESULTADOS); mkpath(DIR_FIGURAS)

# Ford-Fulkerson y Edmonds-Karp de codigo_referencia/ definen AMBOS un
# `struct RedFlujo` y funciones auxiliares con el mismo nombre (son scripts
# educativos independientes, no módulos entre sí). Se aíslan cada uno en su
# propio módulo para poder incluir los dos archivos sin colisión de
# nombres, sin modificarles una sola línea.
module FF
    include(joinpath(@__DIR__, "..", "..", "codigo_referencia", "ford-fulkerson", "ford_fulkerson.jl"))
end
module EK
    include(joinpath(@__DIR__, "..", "..", "codigo_referencia", "edmonds-karp", "edmonds_karp.jl"))
end

println("="^70)
println("FASE 3 -- Optimización en Redes (P5-P7)")
println("="^70)

red = cargar_red()
verificar(red) || error("La carga no reproduce el Anexo A: revise el conjunto de datos.")
g = red.g
n, m = nv(g), ne(g)

# =================================================================
# P6.1 -- Estimación de capacidad (se necesita ya para P5)
# =================================================================
println("\n--- P6.1 Estimación de capacidad c(u,v) ---")
capacidad_estimada, supuesto_capacidad = estimar_capacidad(red)
trafico_raw = red.aristas.trafico_mbps
trafico_estimado = [ismissing(x) ? 0.0 : Float64(x) for x in trafico_raw]
n_sin_trafico = count(ismissing, trafico_raw)
@printf("Capacidad estimada para %d aristas (%d documentadas + %d estimadas).\n",
        m, count(!ismissing, red.aristas.capacidad_mbps), m - count(!ismissing, red.aristas.capacidad_mbps))
@printf("Tráfico medido ausente en %d aristas: se asume 0 Mbps para w_carga (enlaces sin lectura ", n_sin_trafico)
println("del monitor en el instante de la captura -- ver supuesto documentado por arista).")

df_capacidad = DataFrame(
    origen=red.aristas.source, destino=red.aristas.target, rol=red.aristas.rol,
    capacidad_mbps=capacidad_estimada, supuesto=supuesto_capacidad, trafico_mbps=trafico_estimado,
)
guardar_csv(df_capacidad, joinpath(DIR_RESULTADOS, "p6_capacidad_estimada.csv"))
println("Resumen de supuestos de capacidad:")
for (s, k) in sort(collect(countmap([replace(s, r" \[.*\]" => "") for s in supuesto_capacidad])); by=x -> -x[2])
    @printf("  %-90s %3d\n", s, k)
end

# =================================================================
# P5.1 -- Dijkstra (cola de prioridad propia) y Floyd-Warshall
# =================================================================
println("\n--- P5.1 Dijkstra y Floyd-Warshall ---")

"Montículo binario mínimo minimalista: `push!`/`pop!` en O(log n), suficiente para Dijkstra."
struct MiniHeap
    datos::Vector{Tuple{Float64,Int}}
end
MiniHeap() = MiniHeap(Tuple{Float64,Int}[])
function heap_push!(h::MiniHeap, prioridad::Float64, item::Int)
    push!(h.datos, (prioridad, item))
    i = length(h.datos)
    while i > 1
        p = i ÷ 2
        h.datos[p][1] <= h.datos[i][1] && break
        h.datos[p], h.datos[i] = h.datos[i], h.datos[p]
        i = p
    end
end
function heap_pop!(h::MiniHeap)
    tope = h.datos[1]
    h.datos[1] = h.datos[end]
    pop!(h.datos)
    i, nn = 1, length(h.datos)
    while true
        l, r, menor = 2i, 2i + 1, i
        l <= nn && h.datos[l][1] < h.datos[menor][1] && (menor = l)
        r <= nn && h.datos[r][1] < h.datos[menor][1] && (menor = r)
        menor == i && break
        h.datos[i], h.datos[menor] = h.datos[menor], h.datos[i]
        i = menor
    end
    return tope
end
heap_vacio(h::MiniHeap) = isempty(h.datos)

"""
    dijkstra_desde_cero(g, W, origen) -> (dist, padre)

Dijkstra con montículo binario propio (`MiniHeap`), O((V+E) log V). `W` es
la matriz de pesos densa (n×n, `Inf` donde no hay arista). Requiere pesos
no negativos -- ver discusión P5.5: ninguno de los tres modelos de esta
red puede producir un peso negativo, así que la precondición siempre se
cumple.
"""
function dijkstra_desde_cero(g::SimpleGraph, W::Matrix{Float64}, origen::Int)
    n = nv(g)
    dist = fill(Inf, n)
    padre = zeros(Int, n)
    visitado = falses(n)
    dist[origen] = 0.0
    h = MiniHeap()
    heap_push!(h, 0.0, origen)
    while !heap_vacio(h)
        d, u = heap_pop!(h)
        visitado[u] && continue
        visitado[u] = true
        for v in neighbors(g, u)
            nd = d + W[u, v]
            if nd < dist[v]
                dist[v] = nd
                padre[v] = u
                heap_push!(h, nd, v)
            end
        end
    end
    return dist, padre
end

"Floyd-Warshall clásico, O(V³), sobre la matriz de pesos densa `W` (con 0 en la diagonal e `Inf` donde no hay arista)."
function floyd_warshall_desde_cero(W::Matrix{Float64})
    D = copy(W)
    n = size(D, 1)
    for k in 1:n, i in 1:n, j in 1:n
        nuevo = D[i, k] + D[k, j]
        nuevo < D[i, j] && (D[i, j] = nuevo)
    end
    return D
end

"Matriz de pesos densa n×n para uno de los tres modelos de P5."
function matriz_pesos(red, modelo::Symbol; α::Float64=0.05, β::Float64=12.0)
    n = nv(red.g)
    W = fill(Inf, n, n)
    for i in 1:n
        W[i, i] = 0.0
    end
    for (k, (u, v)) in enumerate(zip(red.aristas.u, red.aristas.v))
        w = modelo == :saltos ? w_saltos(u, v) :
            modelo == :latencia ? w_latencia(capacidad_estimada[k]; α=α, β=β) :
            w_carga(trafico_estimado[k], capacidad_estimada[k])
        W[u, v] = w
        W[v, u] = w
    end
    return W
end

# α (ms): latencia fija de propagación+conmutación por salto, campus-escala
#   (decenas de km a lo sumo, más el retardo de un switch/router L2-L3).
# β (Mbit): tamaño de una trama Ethernet típica (~1500 B ≈ 12 kbit ⇒ β=12
#   da el tiempo de serialización en ms cuando c está en Mbps).
const α_LATENCIA, β_LATENCIA = 0.05, 12.0
W_saltos = matriz_pesos(red, :saltos)
W_latencia = matriz_pesos(red, :latencia; α=α_LATENCIA, β=β_LATENCIA)
W_carga = matriz_pesos(red, :carga)

D_saltos = floyd_warshall_desde_cero(W_saltos)
D_latencia = floyd_warshall_desde_cero(W_latencia)
D_carga = floyd_warshall_desde_cero(W_carga)

println("Verificando Dijkstra vs. Floyd-Warshall en 20 pares al azar, para los 3 modelos...")
Random.seed!(1217)
pares_prueba = [(rand(1:n), rand(1:n)) for _ in 1:20]
for (nombre, W, D) in [("saltos", W_saltos, D_saltos), ("latencia", W_latencia, D_latencia), ("carga", W_carga, D_carga)]
    ok = true
    for (i, j) in pares_prueba
        dist_i, _ = dijkstra_desde_cero(g, W, i)
        ok &= isapprox(dist_i[j], D[i, j]; atol=1e-9)
    end
    @printf("  Modelo %-10s: %s\n", nombre, ok ? "Dijkstra y Floyd-Warshall COINCIDEN en los 20 pares" : "DISCREPANCIA -- revisar")
end

# =================================================================
# P5.2 -- Tiempos de ejecución: Dijkstra vs. Floyd-Warshall
# =================================================================
println("\n--- P5.2 Comparación empírica de tiempos ---")
idx_core_central = red.idx["DATCC-2A-C2"]
orden_bfs, _, _ = let
    visitado = falses(n); distancia = fill(-1, n); orden = Int[]
    cola = Int[idx_core_central]; visitado[idx_core_central] = true; distancia[idx_core_central] = 0
    while !isempty(cola)
        u = popfirst!(cola); push!(orden, u)
        for v in neighbors(g, u)
            if !visitado[v]
                visitado[v] = true; distancia[v] = distancia[u] + 1; push!(cola, v)
            end
        end
    end
    orden, distancia, nothing
end

tamaños = [20, 40, 60, 80, 100, 130, n]
filas_tiempo = NamedTuple[]
for tam in tamaños
    verts = orden_bfs[1:tam]
    sg, _ = induced_subgraph(g, verts)
    Wsg = fill(Inf, tam, tam)
    for i in 1:tam; Wsg[i, i] = 0.0; end
    for e in edges(sg)
        Wsg[src(e), dst(e)] = 1.0; Wsg[dst(e), src(e)] = 1.0
    end
    t_fw = @elapsed floyd_warshall_desde_cero(Wsg)
    t_dij_una = @elapsed dijkstra_desde_cero(sg, Wsg, 1)
    crossover_k = t_dij_una > 0 ? t_fw / t_dij_una : NaN
    push!(filas_tiempo, (n=tam, t_floyd_warshall_s=t_fw, t_dijkstra_1fuente_s=t_dij_una,
                          consultas_equilibrio=crossover_k))
    @printf("  n=%4d   FW=%.5fs (todas las parejas)   Dijkstra(1 fuente)=%.6fs   equilibrio≈%.0f fuentes\n",
            tam, t_fw, t_dij_una, crossover_k)
end
df_tiempos = DataFrame(filas_tiempo)
guardar_csv(df_tiempos, joinpath(DIR_RESULTADOS, "p5_tiempos_dijkstra_fw.csv"))
println("""
Lectura: Floyd-Warshall calcula TODAS las distancias en una sola corrida
(O(V³)); Dijkstra da todas las distancias DESDE UNA fuente (O((V+E) log V)).
"consultas_equilibrio" ≈ cuántas fuentes de Dijkstra cuestan lo mismo que
UNA corrida completa de Floyd-Warshall: por debajo de ese número de fuentes
que realmente interesan, conviene Dijkstra repetido; con más, conviene
Floyd-Warshall una sola vez (esto es justo lo que predice la teoría: para
V grande, V·(V+E)logV crece más lento que V³ en la práctica salvo que casi
todas las V fuentes hagan falta).
""")

# =================================================================
# P5.3-P5.4 -- Ranking por cercanía bajo cada modelo, par más distante
# =================================================================
println("--- P5.3-P5.4 Cercanía por modelo y par de acceso más distante ---")
function top10_cercania(D::Matrix{Float64})
    n = size(D, 1)
    cercania = [1.0 / sum(D[i, :]) for i in 1:n]  # cercanía = 1/suma de distancias (Freeman)
    orden = sortperm(cercania; rev=true)[1:10]
    return red.ids[orden], round.(cercania[orden]; digits=6)
end
id_s, val_s = top10_cercania(D_saltos)
id_l, val_l = top10_cercania(D_latencia)
id_c, val_c = top10_cercania(D_carga)
tabla_cercania_modelos = DataFrame(puesto=1:10, id_saltos=id_s, cercania_saltos=val_s,
                                    id_latencia=id_l, cercania_latencia=val_l,
                                    id_carga=id_c, cercania_carga=val_c)
guardar_csv(tabla_cercania_modelos, joinpath(DIR_RESULTADOS, "p5_top10_cercania_por_modelo.csv"))
println(tabla_cercania_modelos)
println("¿Cambia el ranking según el peso? Coincidencia top-10 saltos∩latencia: ",
        length(intersect(Set(id_s), Set(id_l))), "/10   saltos∩carga: ",
        length(intersect(Set(id_s), Set(id_c))), "/10")
n_pesos_carga_cero = count(k -> trafico_estimado[k] == 0.0, 1:m)
println("Nota sobre w_carga: ", n_pesos_carga_cero, " de $m aristas tienen tráfico medido/asumido en 0 Mbps,",
        "\ny por lo tanto peso EXACTO 0 en este modelo. Eso hace que muchos nodos de acceso conectados por",
        "\nenlaces sin tráfico registrado queden a distancia 0 de su switch y compartan exactamente la misma",
        "\ncercanía agregada (de ahí los empates en la columna cercania_carga): es una propiedad real del",
        "\nmodelo -mide congestión, no topología- y no un error de cálculo.")

idx_acceso = findall(==("acceso"), red.capa)
function par_mas_distante(D::Matrix{Float64})
    mejor = (0.0, 0, 0)
    for i in idx_acceso, j in idx_acceso
        i < j && D[i, j] > mejor[1] && (mejor = (D[i, j], i, j))
    end
    return mejor
end
function reconstruir_camino(padre::Vector{Int}, origen::Int, destino::Int)
    camino = [destino]
    while camino[1] != origen
        camino[1] == 0 && return Int[]
        pushfirst!(camino, padre[camino[1]])
    end
    return camino
end
for (nombre, W, D) in [("saltos", W_saltos, D_saltos), ("latencia", W_latencia, D_latencia), ("carga", W_carga, D_carga)]
    dmax, i, j = par_mas_distante(D)
    _, padre_i = dijkstra_desde_cero(g, W, i)
    camino = reconstruir_camino(padre_i, i, j)
    @printf("\nModelo %s: par de acceso más distante = %s <-> %s  (distancia=%.3f)\n",
            nombre, red.ids[i], red.ids[j], dmax)
    println("  Ruta: ", join(red.ids[camino], " -> "))
end

# =================================================================
# P5.5 -- Pesos negativos y elección de protocolo
# =================================================================
println("\n--- P5.5 Pesos negativos y protocolos reales ---")
println("""
Ninguno de los tres modelos puede producir un peso negativo aquí:
w_saltos=1 siempre; w_latencia=α+β/c con α,β>0 y c>0 (capacidad estimada
siempre positiva); w_carga=b/c con b,c≥0. Dijkstra es válido en los tres
casos sin necesidad de Bellman-Ford.

OSPF/IS-IS usan en la práctica un costo ESTÁTICO (típicamente ∝1/ancho de
banda nominal, muy parecido a w_latencia con tráfico fuera de la fórmula),
NO un costo que dependa del tráfico instantáneo como w_carga. Si el peso
dependiera del tráfico en tiempo real, cada cambio de carga dispararía un
recálculo de rutas (LSA/actualización de estado de enlace), y el tráfico
que huye de un enlace cargado podría sobrecargar la ruta alternativa,
haciendo que ESA se vuelva "cara" y el tráfico regrese -- oscilación de
rutas (route flapping) e inestabilidad, un problema documentado desde los
primeros diseños de enrutamiento adaptativo de ARPANET.
""")

# =================================================================
# P6.2-P6.4 -- Flujo máximo y corte mínimo por campus
# =================================================================
println("="^70); println("P6 -- Flujo máximo y corte mínimo"); println("="^70)

idx_sink = red.idx["INTERNET-MPLS"]
cap_int = round.(Int, capacidad_estimada)
const CAP_SUPERFUENTE = 200_000  # mucho mayor que cualquier capacidad real (10 Gbps=10 000): nunca es el cuello de botella

"Construye la red de flujo (177+1 nodos: +1 super-fuente) para un campus dado, en el struct RedFlujo de FF/EK."
function red_flujo_campus(campus_objetivo::String)
    N = n + 1
    C = zeros(Int, N, N)
    for (k, (u, v)) in enumerate(zip(red.aristas.u, red.aristas.v))
        C[u, v] = cap_int[k]
        C[v, u] = cap_int[k]
    end
    idx_acceso_campus = findall(i -> red.capa[i] == "acceso" && red.campus[i] == campus_objetivo, 1:n)
    for i in idx_acceso_campus
        C[N, i] = CAP_SUPERFUENTE
    end
    nombres = vcat(red.ids, ["SUPERFUENTE"])
    pos = [(0.0, 0.0) for _ in 1:N]  # no se usan las funciones de dibujo, solo las de cómputo
    # FF.RedFlujo y EK.RedFlujo son tipos DISTINTOS (cada módulo aísla su propia
    # definición del struct, ver comentario más arriba): se construye una
    # instancia de cada uno con los mismos datos, en vez de compartir un objeto.
    return FF.RedFlujo(C, nombres, pos), EK.RedFlujo(C, nombres, pos), N, idx_acceso_campus
end

campus_con_acceso = [c for c in sort(unique(red.campus))
                      if c != "Nube MPLS" && count(i -> red.capa[i] == "acceso" && red.campus[i] == c, 1:n) > 0]
println("Campus con nodos de acceso propios (se les calcula flujo máximo hacia INTERNET-MPLS): ",
        join(campus_con_acceso, ", "))

filas_flujo = NamedTuple[]
filas_cortes = NamedTuple[]
for campus in campus_con_acceso
    redflujo_ff, redflujo_ek, idx_fuente, idx_acceso_campus = red_flujo_campus(campus)
    println("\n>>> Campus: ", campus, " (", length(idx_acceso_campus), " nodos de acceso)")

    flujo_ff, F_ff, hist_ff = FF.ford_fulkerson(redflujo_ff, idx_fuente, idx_sink; metodo=:dfs, verbose=false)
    flujo_ek, F_ek, hist_ek = EK.edmonds_karp(redflujo_ek, idx_fuente, idx_sink; verbose=false)
    @assert flujo_ff == flujo_ek "Ford-Fulkerson y Edmonds-Karp deberían dar el mismo flujo máximo"

    S, aristas_corte_todas = FF.corte_minimo(redflujo_ff.C, F_ff, idx_fuente)
    capacidad_corte = sum(redflujo_ff.C[u, v] for (u, v) in aristas_corte_todas)
    @assert capacidad_corte == flujo_ff "Max-flow != capacidad del corte mínimo -- violación del teorema"
    # El arco sintético super-fuente->acceso tiene capacidad artificialmente
    # enorme (CAP_SUPERFUENTE) para no ser nunca el cuello de botella real;
    # se filtra explícitamente en vez de asumir que aparece (o no) en el corte.
    aristas_corte = [(u, v) for (u, v) in aristas_corte_todas if u != idx_fuente]

    longs_ff = [length(p.camino) - 1 for p in hist_ff]
    longs_ek = [length(p.camino) - 1 for p in hist_ek]
    # Propiedad de Edmonds-Karp (BFS): la longitud del camino aumentante NUNCA
    # decrece de una iteración a la siguiente -- es la base de su cota
    # O(V*E). Se verifica explícitamente en vez de solo asumirla.
    @assert issorted(longs_ek) "Edmonds-Karp debería producir longitudes de camino no decrecientes -- violación de la propiedad que garantiza su cota O(V·E)"
    # Sin `init` en minimum/maximum: hist_ff/hist_ek nunca están vacíos aquí
    # (flujo_ff > 0 en los siete campus implica al menos una iteración), así
    # que un `init=0` sería incorrecto (reportaría 0 como longitud mínima
    # aunque ningún camino real tenga longitud 0) y se prefiere que un caso
    # verdaderamente vacío falle con un error claro en vez de mentir con un 0.
    @printf("  Flujo máximo: %d Mbps   FF(DFS): %d iteraciones (long. %d-%d)   EK(BFS): %d iteraciones (long. %d-%d)\n",
            flujo_ff, length(hist_ff), minimum(longs_ff), maximum(longs_ff),
            length(hist_ek), minimum(longs_ek), maximum(longs_ek))
    @printf("  Longitudes de los caminos aumentantes -- FF(DFS): %s\n", join(longs_ff, ", "))
    @printf("  Longitudes de los caminos aumentantes -- EK(BFS): %s (no decreciente ✓)\n", join(longs_ek, ", "))
    @printf("  Corte mínimo verificado a mano: %d aristas físicas, capacidad total = %d Mbps (= flujo máximo ✓)\n",
            length(aristas_corte), capacidad_corte)

    for (u, v) in aristas_corte
        @printf("    corte: %-28s -- %-28s  cap=%5d Mbps\n", red.ids[u], red.ids[v], redflujo_ff.C[u, v])
        push!(filas_cortes, (campus=campus, origen=red.ids[u], destino=red.ids[v],
                              capacidad_mbps=redflujo_ff.C[u, v]))
    end

    push!(filas_flujo, (campus=campus, nodos_acceso=length(idx_acceso_campus),
                         flujo_maximo_mbps=flujo_ff, iteraciones_ff=length(hist_ff),
                         iteraciones_ek=length(hist_ek), aristas_corte=length(aristas_corte),
                         capacidad_corte_mbps=capacidad_corte,
                         longitud_ff_min=minimum(longs_ff), longitud_ff_max=maximum(longs_ff),
                         longitud_ek_min=minimum(longs_ek), longitud_ek_max=maximum(longs_ek),
                         longitudes_ff=join(longs_ff, ";"), longitudes_ek=join(longs_ek, ";")))
end
df_flujo = DataFrame(filas_flujo)
guardar_csv(df_flujo, joinpath(DIR_RESULTADOS, "p6_flujo_maximo_por_campus.csv"))
guardar_csv(DataFrame(filas_cortes), joinpath(DIR_RESULTADOS, "p6_cortes_minimos_aristas.csv"))
println("\nResumen de flujo máximo por campus:")
println(df_flujo)

_, puentes_p1 = puentes_y_articulacion(g)
set_puentes = Set([(min(u, v), max(u, v)) for (u, v) in puentes_p1])
println("\nInterpretación P6.4: los cortes mínimos anteriores están compuestos por enlaces reales;",
        "\ncuando ese enlace también aparece en la lista de puentes de P1 (p6_flujo... vs p1_puentes.csv),",
        "\nsignifica que la falla de ESE ÚNICO enlace ya reduce a cero el flujo de salida del campus --",
        "\nel cuello de botella físico y el punto de fragilidad topológica coinciden.")

# =================================================================
# P6.5 -- Flujo de costo mínimo
# =================================================================
println("\n--- P6.5 Flujo de costo mínimo ---")
"""
    mincostflow_ssp(C, costo, s, t, objetivo) -> (flujo, costo_total, F)

Flujo de costo mínimo por caminos aumentantes sucesivos más baratos
(Bellman-Ford sobre la red residual, que admite arcos residuales de costo
negativo -- por eso no se usa Dijkstra aquí). Se detiene al alcanzar
`objetivo` unidades de flujo o al agotarse los caminos aumentantes.
Complejidad O(objetivo · V · E) en el peor caso (cada aumento manda al
menos 1 unidad); aceptable para las demandas pequeñas de este problema.
"""
function mincostflow_ssp(C::Matrix{Int}, costo::Matrix{Int}, s::Int, t::Int, objetivo::Int)
    Nn = size(C, 1)
    F = zeros(Int, Nn, Nn)
    flujo_total = 0
    costo_total = 0
    while flujo_total < objetivo
        dist = fill(Inf, Nn); dist[s] = 0.0
        padre = zeros(Int, Nn)
        for _ in 1:Nn-1
            cambiado = false
            for u in 1:Nn, v in 1:Nn
                (C[u, v] - F[u, v]) > 0 || continue
                c_uv = C[u, v] > 0 ? costo[u, v] : -costo[v, u]
                if dist[u] + c_uv < dist[v]
                    dist[v] = dist[u] + c_uv; padre[v] = u; cambiado = true
                end
            end
            cambiado || break
        end
        dist[t] == Inf && break
        camino = [t]
        while camino[1] != s
            pushfirst!(camino, padre[camino[1]])
        end
        Δ = minimum(C[camino[i], camino[i+1]] - F[camino[i], camino[i+1]] for i in 1:length(camino)-1)
        Δ = min(Δ, objetivo - flujo_total)
        for i in 1:length(camino)-1
            u, v = camino[i], camino[i+1]
            c_uv = C[u, v] > 0 ? costo[u, v] : -costo[v, u]
            F[u, v] += Δ; F[v, u] -= Δ
            costo_total += Δ * c_uv
        end
        flujo_total += Δ
    end
    return flujo_total, costo_total, F
end

# Demanda fija desde dos campus simultáneamente hacia Internet/MPLS, con
# costo por salto (1 por defecto, 5 si el enlace es de respaldo -- para
# penalizar el uso de rutas de respaldo salvo que sean estrictamente
# necesarias). Red combinada: una super-fuente para cada campus, cada una
# con una arista de capacidad EXACTA igual a su demanda (para forzar que
# se sirva esa demanda, no más).
campus_demanda = campus_con_acceso[1:min(2, length(campus_con_acceso))]
N2 = n + length(campus_demanda)
C2 = zeros(Int, N2, N2)
costo2 = ones(Int, N2, N2)
for (k, (u, v)) in enumerate(zip(red.aristas.u, red.aristas.v))
    C2[u, v] = cap_int[k]; C2[v, u] = cap_int[k]
    c = red.aristas.rol[k] == "respaldo" ? 5 : 1
    costo2[u, v] = c; costo2[v, u] = c
end
demandas = Int[]
for (i, campus) in enumerate(campus_demanda)
    idx_fs = n + i
    idx_acceso_campus = findall(j -> red.capa[j] == "acceso" && red.campus[j] == campus, 1:n)
    flujo_campus = df_flujo[df_flujo.campus .== campus, :flujo_maximo_mbps][1]
    demanda_campus = max(1, flujo_campus ÷ 2)  # la mitad del flujo máximo individual, para que la red combinada sí alcance a servirla
    push!(demandas, demanda_campus)
    for j in idx_acceso_campus
        C2[idx_fs, j] = CAP_SUPERFUENTE
    end
end
# Súper-súper-fuente con arcos de capacidad EXACTA = demanda de cada campus
N3 = N2 + 1
C3 = zeros(Int, N3, N3); C3[1:N2, 1:N2] = C2
costo3 = ones(Int, N3, N3); costo3[1:N2, 1:N2] = costo2
for (i, d) in enumerate(demandas)
    C3[N3, n+i] = d
end
objetivo_total = sum(demandas)
flujo_mcf, costo_mcf, _ = mincostflow_ssp(C3, costo3, N3, idx_sink, objetivo_total)
@printf("Flujo de costo mínimo: demanda fija de %s desde %s\n",
        join(["$(campus_demanda[i])=$(demandas[i])" for i in 1:length(demandas)], " y "), "hacia INTERNET-MPLS")
@printf("  Flujo servido: %d de %d Mbps solicitados   Costo total: %d (unidades de costo·Mbps)\n",
        flujo_mcf, objetivo_total, costo_mcf)
flujo_max_puro_combinado = sum(df_flujo[in.(df_flujo.campus, Ref(campus_demanda)), :flujo_maximo_mbps])
@printf("  Comparación: la suma de los flujos máximos individuales de esos campus es %d Mbps;\n", flujo_max_puro_combinado)
println("  el flujo de costo mínimo sirve la demanda solicitada (menor) al menor costo posible,",
        "\n  evitando enlaces de respaldo salvo que no quede alternativa -- un objetivo distinto",
        "\n  al de flujo máximo puro, que no le importa el costo, solo maximizar el caudal.")

guardar_csv(DataFrame(campus=campus_demanda, demanda_mbps=demandas), joinpath(DIR_RESULTADOS, "p6_demanda_costo_minimo.csv"))

# =================================================================
# P7 -- Localización de instalaciones (p-mediana, p-centro)
# =================================================================
println("\n" * "="^70); println("P7 -- Localización de instalaciones"); println("="^70)
println("Distancia usada: número de saltos (D_saltos de P5) -- el criterio estándar de p-mediana/p-centro.")

function greedy_p_mediana(D::Matrix{Float64}, p::Int)
    nn = size(D, 1)
    elegidos = Int[]
    for _ in 1:p
        mejor_nodo, mejor_costo = 0, Inf
        for cand in setdiff(1:nn, elegidos)
            conjunto = vcat(elegidos, cand)
            costo = sum(minimum(D[i, conjunto]) for i in 1:nn)
            if costo < mejor_costo
                mejor_costo, mejor_nodo = costo, cand
            end
        end
        push!(elegidos, mejor_nodo)
    end
    return elegidos
end
function greedy_p_centro(D::Matrix{Float64}, p::Int)
    nn = size(D, 1)
    elegidos = Int[]
    for _ in 1:p
        mejor_nodo, mejor_costo = 0, Inf
        for cand in setdiff(1:nn, elegidos)
            conjunto = vcat(elegidos, cand)
            costo = maximum(minimum(D[i, conjunto]) for i in 1:nn)
            if costo < mejor_costo
                mejor_costo, mejor_nodo = costo, cand
            end
        end
        push!(elegidos, mejor_nodo)
    end
    return elegidos
end

function ilp_p_mediana(D::Matrix{Float64}, p::Int)
    nn = size(D, 1)
    modelo = Model(HiGHS.Optimizer); set_silent(modelo)
    @variable(modelo, y[1:nn], Bin)
    @variable(modelo, x[1:nn, 1:nn], Bin)
    @constraint(modelo, sum(y) == p)
    @constraint(modelo, [i = 1:nn], sum(x[i, :]) == 1)
    @constraint(modelo, [i = 1:nn, j = 1:nn], x[i, j] <= y[j])
    @objective(modelo, Min, sum(D[i, j] * x[i, j] for i in 1:nn, j in 1:nn))
    optimize!(modelo)
    return findall(k -> value(y[k]) > 0.5, 1:nn), objective_value(modelo)
end
function ilp_p_centro(D::Matrix{Float64}, p::Int)
    nn = size(D, 1)
    modelo = Model(HiGHS.Optimizer); set_silent(modelo)
    @variable(modelo, y[1:nn], Bin)
    @variable(modelo, x[1:nn, 1:nn], Bin)
    @variable(modelo, z >= 0)
    @constraint(modelo, sum(y) == p)
    @constraint(modelo, [i = 1:nn], sum(x[i, :]) == 1)
    @constraint(modelo, [i = 1:nn, j = 1:nn], x[i, j] <= y[j])
    @constraint(modelo, [i = 1:nn], sum(D[i, j] * x[i, j] for j in 1:nn) <= z)
    @objective(modelo, Min, z)
    optimize!(modelo)
    return findall(k -> value(y[k]) > 0.5, 1:nn), objective_value(modelo)
end

println("""
Formulación p-mediana: min Σᵢ Σⱼ d(i,j)·xᵢⱼ  s.a. Σⱼ yⱼ=p, xᵢⱼ≤yⱼ, Σⱼxᵢⱼ=1, x,y∈{0,1}
Formulación p-centro:  min z            s.a. Σⱼ d(i,j)·xᵢⱼ ≤ z ∀i, (mismas restricciones que arriba)
""")

filas_p7 = NamedTuple[]
for p in [1, 2, 3, 5]
    g_med = greedy_p_mediana(D_saltos, p)
    g_cen = greedy_p_centro(D_saltos, p)
    i_med, obj_med = ilp_p_mediana(D_saltos, p)
    i_cen, obj_cen = ilp_p_centro(D_saltos, p)
    costo_greedy_med = sum(minimum(D_saltos[i, g_med]) for i in 1:n)
    costo_greedy_cen = maximum(minimum(D_saltos[i, g_cen]) for i in 1:n)
    @printf("\np=%d\n", p)
    @printf("  p-mediana  greedy: %s (costo=%.1f)   ILP: %s (costo=%.1f)\n",
            join(red.ids[g_med], ", "), costo_greedy_med, join(red.ids[i_med], ", "), obj_med)
    @printf("  p-centro   greedy: %s (radio=%.1f)   ILP: %s (radio=%.1f)\n",
            join(red.ids[g_cen], ", "), costo_greedy_cen, join(red.ids[i_cen], ", "), obj_cen)
    push!(filas_p7, (p=p, mediana_greedy=join(red.ids[g_med], ";"), mediana_ilp=join(red.ids[i_med], ";"),
                      costo_mediana_ilp=obj_med, centro_greedy=join(red.ids[g_cen], ";"),
                      centro_ilp=join(red.ids[i_cen], ";"), radio_centro_ilp=obj_cen))
end
df_p7 = DataFrame(filas_p7)
guardar_csv(df_p7, joinpath(DIR_RESULTADOS, "p7_localizacion.csv"))

cent_intermediacion = betweenness_centrality(g)
cent_cercania = closeness_centrality(g)
top5_intermediacion = red.ids[sortperm(cent_intermediacion; rev=true)[1:5]]
top5_cercania = red.ids[sortperm(cent_cercania; rev=true)[1:5]]
println("\nP7.3 -- Comparación con P1: top-5 por intermediación: ", join(top5_intermediacion, ", "))
println("                              top-5 por cercanía:       ", join(top5_cercania, ", "))
println("""
Un colector de telemetría no necesita estar en el camino de MUCHOS pares
(alta intermediación) ni "en el medio de todo" en distancia agregada (alta
cercanía): necesita estar DONDE HAY EQUIPOS QUE MONITOREAR sin quedar lejos
de ninguno. La p-mediana minimiza distancia PROMEDIO (bueno para telemetría
agregada, tolera algunos equipos algo lejos); el p-centro minimiza la
distancia PEOR CASO (garantiza un techo de latencia de sondeo para TODOS los
equipos). Los nodos de mayor grado/intermediación suelen ser el CORE, que
está muy cerca de pocos equipos de acceso pero lejos de la mayoría de los
equipos hoja -- por eso raramente coincide con la solución óptima de P7.

P7.4 -- Restricciones que el modelo omite: espacio en rack y alimentación
eléctrica redundante en el punto elegido, puertos de administración
disponibles en el equipo huésped, seguridad física del cuarto de
telecomunicaciones, y el costo/licenciamiento del colector (p.ej. sondas
NetFlow por puerto vs. una sonda SNMP centralizada). Incorporarlas
requeriría, como mínimo, restringir el conjunto de candidatos a equipos con
espacio de rack certificado y sumar un costo fijo de licencia por sitio a
la función objetivo (una variante de "localización con costos fijos").
""")

println("\n" * "="^70)
println("FASE 3 completa. Resultados en ", DIR_RESULTADOS)
println("="^70)
