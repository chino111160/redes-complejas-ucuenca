#=============================================================
  fase2_p3_p4_recorrido_particion.jl -- Fase 2: Recorrido y Partición
  Proyecto integrador de Redes Complejas (1217) -- Universidad de Cuenca

  Resuelve P3 (BFS/DFS) y P4 (comunidades y modularidad) del enunciado
  (proyecto_red_ucuenca.pdf, cap. 4).

  Uso:
      julia --project=analisis analisis/scripts/fase2_p3_p4_recorrido_particion.jl
=============================================================#

using Graphs
using DataFrames
using CSV
using Statistics
using StatsBase
using Printf
using Random
using LinearAlgebra

const DIR_ESTE_SCRIPT = @__DIR__
const DIR_ANALISIS = dirname(DIR_ESTE_SCRIPT)
const DIR_RAIZ = dirname(DIR_ANALISIS)

include(joinpath(DIR_RAIZ, "codigo_base", "cargar_red.jl"))
include(joinpath(DIR_ANALISIS, "src", "ModeloRed.jl"))

const DIR_RESULTADOS = joinpath(DIR_ANALISIS, "resultados")
const DIR_FIGURAS = joinpath(DIR_ANALISIS, "figuras")
mkpath(DIR_RESULTADOS); mkpath(DIR_FIGURAS)

println("="^70)
println("FASE 2 -- Recorrido y Partición (P3-P4)")
println("="^70)

red = cargar_red()
verificar(red) || error("La carga no reproduce el Anexo A: revise el conjunto de datos.")
g = red.g
n, m = nv(g), ne(g)

# =================================================================
# P3.1 -- BFS y DFS desde cero
# =================================================================
println("\n--- P3.1 BFS y DFS desde cero ---")

"""
    bfs_desde_cero(g, origen) -> (orden, distancia, padre)

BFS con una cola FIFO (`Vector{Int}` + `popfirst!`) sobre la lista de
adyacencia de `g`. Complejidad O(V+E): cada nodo se encola una vez y cada
arista se examina como máximo dos veces (una desde cada extremo).
`distancia[v] = -1` si `v` no es alcanzable desde `origen`.
"""
function bfs_desde_cero(g::SimpleGraph, origen::Int)
    n = nv(g)
    visitado = falses(n)
    distancia = fill(-1, n)
    padre = zeros(Int, n)
    orden = Int[]
    cola = Int[origen]
    visitado[origen] = true
    distancia[origen] = 0
    while !isempty(cola)
        u = popfirst!(cola)
        push!(orden, u)
        for v in neighbors(g, u)
            if !visitado[v]
                visitado[v] = true
                distancia[v] = distancia[u] + 1
                padre[v] = u
                push!(cola, v)
            end
        end
    end
    return orden, distancia, padre
end

"""
    dfs_desde_cero(g, origen) -> (orden, padre)

DFS iterativo con una pila explícita de `(nodo, siguiente_vecino_a_probar)`
(evita recursión para no toparse con el límite de pila de Julia en grafos
grandes). Misma complejidad O(V+E) que BFS; la diferencia es la estructura
de datos (pila LIFO en vez de cola FIFO), que determina el orden de
visita y, con ello, el árbol de recorrido resultante.
"""
function dfs_desde_cero(g::SimpleGraph, origen::Int)
    n = nv(g)
    visitado = falses(n)
    padre = zeros(Int, n)
    orden = Int[]
    pila = [(origen, 1)]
    visitado[origen] = true
    vecinos_cache = Dict{Int,Vector{Int}}()
    while !isempty(pila)
        u, k = pila[end]
        vs = get!(() -> collect(neighbors(g, u)), vecinos_cache, u)
        if k == 1
            push!(orden, u)
        end
        if k <= length(vs)
            pila[end] = (u, k + 1)
            v = vs[k]
            if !visitado[v]
                visitado[v] = true
                padre[v] = u
                push!(pila, (v, 1))
            end
        else
            pop!(pila)
        end
    end
    return orden, padre
end

# Verificación cruzada contra Graphs.jl (permitido solo para verificar, no
# para reemplazar la implementación propia).
_, dist_propia, _ = bfs_desde_cero(g, 1)
dist_biblioteca = gdistances(g, 1)
@assert dist_propia == dist_biblioteca "BFS propio no coincide con Graphs.jl"
orden_dfs_propia, _ = dfs_desde_cero(g, 1)
@assert Set(orden_dfs_propia) == Set(1:n) "DFS propio no visita todos los nodos"
println("BFS y DFS propios verificados contra Graphs.jl (mismas distancias / mismo conjunto visitado).")

# =================================================================
# P3.2-P3.3 -- Perfil de profundidad desde el core y desde el nodo MPLS
# =================================================================
println("\n--- P3.2-P3.3 Perfiles de profundidad (BFS) ---")

idx_core_central = red.idx["DATCC-2A-C2"]  # uno de los dos switches de core redundantes del Campus Central
idx_core_central_b = red.idx["DATCC-2A-C3"]
idx_mpls = red.idx["INTERNET-MPLS"]

function perfil_profundidad(distancia::Vector{Int})
    válidas = filter(>=(0), distancia)
    conteo = countmap(válidas)
    DataFrame(distancia=sort(collect(keys(conteo))),
              nodos=[conteo[d] for d in sort(collect(keys(conteo)))])
end

_, dist_core, _ = bfs_desde_cero(g, idx_core_central)
_, dist_core_b, _ = bfs_desde_cero(g, idx_core_central_b)
_, dist_mpls, _ = bfs_desde_cero(g, idx_mpls)

perfil_core = perfil_profundidad(dist_core)
perfil_mpls = perfil_profundidad(dist_mpls)
guardar_csv(perfil_core, joinpath(DIR_RESULTADOS, "p3_perfil_core_central.csv"))
guardar_csv(perfil_mpls, joinpath(DIR_RESULTADOS, "p3_perfil_mpls.csv"))

println("Perfil de profundidad desde DATCC-2A-C2 (core, Campus Central):")
println(perfil_core)
println("(desde el otro switch de core redundante, DATCC-2A-C3, ",
        dist_core == dist_core_b ? "el perfil es IDÉNTICO" :
        "el perfil difiere en $(count(dist_core .!= dist_core_b)) nodos", ")")
println("\nPerfil de profundidad desde INTERNET-MPLS (nube MPLS):")
println(perfil_mpls)

capa_por_distancia_core = combine(groupby(
    DataFrame(distancia=dist_core, capa=red.capa)[dist_core .>= 0, :], :distancia),
    :capa => (x -> join(sort(collect(countmap(x)); by=p->-p[2]) .|> p -> "$(p[1])=$(p[2])", ", ")) => :capas)
println("\nCapa jerárquica por nivel de distancia desde el core (Campus Central):")
println(capa_por_distancia_core)
println("Interpretación P3.2: si la jerarquía core(0)->agregación(~1)->acceso(~2) declarada por el",
        "\ninforme se cumpliera al pie de la letra, casi todos los nodos deberían quedar a distancia",
        "\n<=2. Distancias mayores indican acceso indirecto (p.ej. un equipo de acceso colgado de otro",
        "\nequipo de acceso) o campus remotos que solo llegan al core vía la nube MPLS (más saltos).")

campus_más_lejos_core = combine(groupby(
    DataFrame(campus=red.campus, distancia=dist_core)[dist_core .>= 0, :], :campus),
    :distancia => maximum => :distancia_maxima, :distancia => mean => :distancia_media)
println("\nDistancia máxima/media desde el core de Campus Central, por campus:")
println(sort(campus_más_lejos_core, :distancia_maxima; rev=true))

campus_más_lejos_mpls = combine(groupby(
    DataFrame(campus=red.campus, distancia=dist_mpls)[dist_mpls .>= 0, :], :campus),
    :distancia => maximum => :distancia_maxima, :distancia => mean => :distancia_media)
println("\nDistancia máxima/media desde la nube MPLS, por campus (P3.3 -- ¿qué campus queda más 'lejos'?):")
println(sort(campus_más_lejos_mpls, :distancia_maxima; rev=true))

# =================================================================
# P3.4 -- Ciclos del grafo (DFS)
# =================================================================
println("\n--- P3.4 Ciclos del grafo ---")
# Con un árbol de expansión DFS de n-1 aristas, cada una de las m-(n-1)
# aristas restantes ("de retroceso", no usadas por el árbol) cierra
# exactamente un ciclo fundamental. El número de aristas de retroceso es
# el número ciclomático m-n+componentes = 209-177+1 = 33 y NO depende del
# árbol elegido; sí depende del árbol cuáles aristas concretas se marcan
# como "de retroceso" (las demás combinaciones de árbol muestran otro
# representante del mismo ciclo).
_, padre_dfs = dfs_desde_cero(g, idx_core_central)
es_arista_arbol(u, v) = padre_dfs[v] == u || padre_dfs[u] == v
aristas_ciclo = [(src(e), dst(e)) for e in edges(g) if !es_arista_arbol(src(e), dst(e))]
num_ciclomático = m - n + length(connected_components(g))
@assert length(aristas_ciclo) == num_ciclomático "Conteo de ciclos inconsistente"
@printf("Número ciclomático (m - n + componentes): %d\n", num_ciclomático)

df_ciclos = DataFrame(
    origen = [red.ids[u] for (u, v) in aristas_ciclo],
    destino = [red.ids[v] for (u, v) in aristas_ciclo],
    campus_origen = [red.campus[u] for (u, v) in aristas_ciclo],
    campus_destino = [red.campus[v] for (u, v) in aristas_ciclo],
)
guardar_csv(df_ciclos, joinpath(DIR_RESULTADOS, "p3_aristas_no_arbol_ciclos.csv"))
println("Aristas 'de retroceso' (generadoras de ciclo) por campus:")
for (c, k) in sort(collect(countmap(vcat(df_ciclos.campus_origen, df_ciclos.campus_destino))); by=x -> -x[2])
    @printf("  %-24s %3d\n", c, k)
end
articulacion, puentes = puentes_y_articulacion(g)
println("De las $num_ciclomático aristas que NO son de árbol, ninguna puede ser puente (por definición: ",
        "un puente no está en ningún ciclo). De las ", n - 1, " aristas SÍ de árbol, ",
        length(puentes), " resultan ser puentes -- las otras ", (n - 1) - length(puentes),
        " son de árbol pero quedan 'cubiertas' por algún ciclo que las hace prescindibles.")
println("Ausencia de ciclos en una zona = ausencia de camino alternativo = zona sin redundancia,",
        "\nexactamente donde P1 ya localizó los puentes: son el mismo fenómeno visto desde dos algoritmos distintos.")

# =================================================================
# P3.5 -- BFS vs DFS para una tarea física
# =================================================================
println("\n--- P3.5 ¿Qué recorrido modela mejor una visita física? ---")
println("""
Un técnico que camina por los edificios no puede "teletransportarse" entre
ramas del árbol como sí lo hace BFS (que visita todos los nodos a distancia
1, después todos los de distancia 2, saltando de rama en rama y obligando a
volver físicamente al core entre visitas). DFS, en cambio, avanza por una
rama hasta el final antes de retroceder -- igual que caminar un pasillo o
un edificio completo antes de pasar al siguiente. El ORDEN que produce DFS
es directamente una ruta física razonable (con retrocesos = caminar de
vuelta el mismo tramo); BFS produce un orden por niveles que, llevado a
recorrido físico real, multiplicaría los trayectos de ida y vuelta al core.
""")

# =================================================================
# P4.1 -- Louvain (implementación propia, Blondel et al. 2008)
# =================================================================
println("\n--- P4.1 Louvain ---")

"""
Representación de grafo ponderado para Louvain: aristas `(u<v, peso)`
entre nodos distintos más un vector de "autolazos" (peso interno de cada
super-nodo, acumulado al agregar comunidades en fases sucesivas). El
autolazo se guarda como el peso interno real (sin duplicar); el grado
ponderado de un nodo suma sus aristas incidentes MÁS dos veces su
autolazo, como es estándar en la formulación de modularidad con self-loops.
"""
struct GrafoLouvain
    n::Int
    aristas::Vector{Tuple{Int,Int,Float64}}
    autolazo::Vector{Float64}
end
GrafoLouvain(g::SimpleGraph) = GrafoLouvain(
    nv(g), [(min(src(e), dst(e)), max(src(e), dst(e)), 1.0) for e in edges(g)],
    zeros(Float64, nv(g)))

function grado_ponderado(gl::GrafoLouvain)
    k = 2 .* copy(gl.autolazo)
    for (u, v, w) in gl.aristas
        k[u] += w; k[v] += w
    end
    return k
end
peso_total(gl::GrafoLouvain) = sum(w for (_, _, w) in gl.aristas; init=0.0) + sum(gl.autolazo)
function lista_adyacencia(gl::GrafoLouvain)
    adj = [Tuple{Int,Float64}[] for _ in 1:gl.n]
    for (u, v, w) in gl.aristas
        push!(adj[u], (v, w)); push!(adj[v], (u, w))
    end
    return adj
end

"""
Fase de movimiento local: cada nodo se reasigna, de forma voraz, a la
comunidad vecina que maximiza la ganancia de modularidad `ΔQ` (fórmula
simplificada de Blondel et al. 2008, ec. 2), hasta que una pasada
completa no mueve a nadie.
"""
function louvain_local(gl::GrafoLouvain; rng=Random.default_rng())
    n = gl.n
    k = grado_ponderado(gl)
    m = peso_total(gl)
    adj = lista_adyacencia(gl)
    comunidad = collect(1:n)
    Σtot = copy(k)
    mejora_global = true
    while mejora_global
        mejora_global = false
        for i in shuffle(rng, 1:n)
            c0 = comunidad[i]
            Σtot[c0] -= k[i]
            pesos = Dict{Int,Float64}()
            for (v, w) in adj[i]
                v == i && continue
                pesos[comunidad[v]] = get(pesos, comunidad[v], 0.0) + w
            end
            # Maximizar ΔQ ∝ peso_hacia_c - Σtot[c]*k[i]/m (factor 1/(2m) común omitido)
            mejor_c, mejor_score = c0, get(pesos, c0, 0.0) - Σtot[c0] * k[i] / m
            for (c, win) in pesos
                c == c0 && continue
                score = win - Σtot[c] * k[i] / m
                if score > mejor_score + 1e-10
                    mejor_score, mejor_c = score, c
                end
            end
            Σtot[mejor_c] += k[i]
            if mejor_c != c0
                comunidad[i] = mejor_c
                mejora_global = true
            end
        end
    end
    return comunidad
end

"Agrega cada comunidad en un super-nodo: aristas entre comunidades se suman; aristas internas pasan a autolazo."
function louvain_agregar(gl::GrafoLouvain, comunidad::Vector{Int})
    etiquetas = sort(unique(comunidad))
    idx = Dict(c => i for (i, c) in enumerate(etiquetas))
    nc = length(etiquetas)
    peso_entre = Dict{Tuple{Int,Int},Float64}()
    autolazo2 = zeros(Float64, nc)
    for (u, v, w) in gl.aristas
        cu, cv = idx[comunidad[u]], idx[comunidad[v]]
        if cu == cv
            autolazo2[cu] += w
        else
            clave = cu < cv ? (cu, cv) : (cv, cu)
            peso_entre[clave] = get(peso_entre, clave, 0.0) + w
        end
    end
    for i in 1:gl.n
        autolazo2[idx[comunidad[i]]] += gl.autolazo[i]
    end
    aristas2 = [(k[1], k[2], w) for (k, w) in peso_entre]
    return GrafoLouvain(nc, aristas2, autolazo2), idx
end

"""
    louvain(g; rng, max_niveles) -> Vector{Int}

Louvain completo: alterna movimiento local y agregación (Blondel et al.
2008) hasta que agregar ya no reduce el número de super-nodos. Devuelve
la comunidad final de cada nodo ORIGINAL de `g` (no de los super-nodos).
"""
function louvain(g::SimpleGraph; rng=Random.default_rng(), max_niveles::Int=50)
    n0 = nv(g)
    gl = GrafoLouvain(g)
    pertenece = collect(1:n0)
    for _ in 1:max_niveles
        com = louvain_local(gl; rng=rng)
        for i in 1:n0
            pertenece[i] = com[pertenece[i]]
        end
        gl2, idx = louvain_agregar(gl, com)
        for i in 1:n0
            pertenece[i] = idx[pertenece[i]]
        end
        gl2.n == gl.n && break
        gl = gl2
    end
    return pertenece
end

modularidad(g::SimpleGraph, comunidad::Vector{Int}) = let
    m2 = 2 * ne(g); grados_ = degree(g); Q = 0.0
    for c in unique(comunidad)
        miembros = findall(==(c), comunidad)
        e_in = sum(count(==(c), comunidad[neighbors(g, u)]) for u in miembros; init=0)
        sum_deg = sum(grados_[miembros])
        Q += e_in / m2 - (sum_deg / m2)^2
    end
    Q
end

Random.seed!(1217)
comunidad_louvain = louvain(g)
Q_louvain = modularidad(g, comunidad_louvain)
@printf("Louvain: %d comunidades, Q = %.4f\n", length(unique(comunidad_louvain)), Q_louvain)

println("\nEstabilidad con 5 semillas distintas:")
resultados_semillas = NamedTuple[]
for semilla in [1, 2, 3, 42, 2024]
    Random.seed!(semilla)
    com_i = louvain(g)
    push!(resultados_semillas, (semilla=semilla, comunidades=length(unique(com_i)),
                                 Q=round(modularidad(g, com_i); digits=4),
                                 nmi_vs_primera=round(nmi(comunidad_louvain, com_i); digits=4)))
end
df_estabilidad = DataFrame(resultados_semillas)
guardar_csv(df_estabilidad, joinpath(DIR_RESULTADOS, "p4_louvain_estabilidad.csv"))
println(df_estabilidad)

# =================================================================
# P4.2-P4.3 -- Comparación con la partición por campus (NMI, ARI, confusión)
# =================================================================
println("\n--- P4.2-P4.3 Louvain vs. partición por campus ---")
campus_como_entero = Dict(c => i for (i, c) in enumerate(sort(unique(red.campus))))
etiqueta_campus = [campus_como_entero[c] for c in red.campus]

valor_nmi = nmi(comunidad_louvain, etiqueta_campus)
valor_ari = ari(comunidad_louvain, etiqueta_campus)
@printf("NMI(Louvain, campus) = %.4f    ARI(Louvain, campus) = %.4f\n", valor_nmi, valor_ari)

T, com_ids, campus_ids = tabla_contingencia(comunidad_louvain, etiqueta_campus)
campus_por_id = Dict(v => k for (k, v) in campus_como_entero)
df_confusion = DataFrame(T, Symbol.(campus_por_id[c] for c in campus_ids))
insertcols!(df_confusion, 1, :comunidad_louvain => com_ids)
guardar_csv(df_confusion, joinpath(DIR_RESULTADOS, "p4_matriz_confusion_comunidad_campus.csv"))
println("\nMatriz de confusión comunidad de Louvain × campus:")
println(df_confusion)

# Nodos "discrepantes": su comunidad de Louvain está dominada por un campus distinto al propio
dominante_por_comunidad = Dict{Int,Int}()
for c in unique(comunidad_louvain)
    miembros = findall(==(c), comunidad_louvain)
    dominante_por_comunidad[c] = mode(etiqueta_campus[miembros])
end
idx_discrepantes = [i for i in 1:n if dominante_por_comunidad[comunidad_louvain[i]] != etiqueta_campus[i]]
df_discrepantes = DataFrame(
    id = red.ids[idx_discrepantes],
    campus_propio = red.campus[idx_discrepantes],
    comunidad = comunidad_louvain[idx_discrepantes],
    campus_dominante_de_su_comunidad = [campus_por_id[dominante_por_comunidad[comunidad_louvain[i]]] for i in idx_discrepantes],
)
guardar_csv(df_discrepantes, joinpath(DIR_RESULTADOS, "p4_nodos_discrepantes.csv"))
@printf("\nNodos cuya comunidad de Louvain está dominada por un campus distinto al propio: %d de %d\n",
        length(idx_discrepantes), n)
println("En ingeniería de red esto señala equipos cuyo tráfico/estructura los liga MÁS a otro campus",
        "\nque al suyo propio (típicamente enlaces WAN/MPLS o equipos frontera) -- ver tabla completa en",
        "\np4_nodos_discrepantes.csv.")

# =================================================================
# P4.4 -- k-means sobre el espectro del laplaciano
# =================================================================
println("\n--- P4.4 k-means sobre el espectro del laplaciano ---")
# k-means desde cero: funciones copiadas de
# codigo_referencia/kmeans/ejemplo1.jl (líneas 1-149; se omite la demo
# final de ese archivo, que genera sus propios datos sintéticos).
dist2(x::AbstractVector, μ::AbstractVector)::Float64 = dot(x .- μ, x .- μ)
function init_plusplus(X::Matrix, K::Int)
    nfil = size(X, 1)
    μ = [X[rand(1:nfil), :]]
    for _ in 2:K
        D = [minimum(dist2(X[i, :], c) for c in μ) for i in 1:nfil]
        probs = D ./ sum(D)
        j = findfirst(cumsum(probs) .≥ rand())
        push!(μ, X[j, :])
    end
    return μ
end
function assign_clusters(X::Matrix, μ::Vector)::Vector{Int}
    nfil = size(X, 1)
    labels = Vector{Int}(undef, nfil)
    for i in 1:nfil
        labels[i] = argmin([dist2(X[i, :], c) for c in μ])
    end
    return labels
end
function update_centroids(X::Matrix, labels::Vector{Int}, K::Int)::Vector
    μ_new = Vector(undef, K)
    for k in 1:K
        mask = findall(labels .== k)
        μ_new[k] = isempty(mask) ? X[rand(1:size(X, 1)), :] : vec(mean(X[mask, :], dims=1))
    end
    return μ_new
end
function my_kmeans(X::Matrix, K::Int; max_iter=300, tol=1e-6, seed=42)
    Random.seed!(seed)
    μ = init_plusplus(X, K)
    labels = assign_clusters(X, μ)
    for _ in 1:max_iter
        μ_old = deepcopy(μ)
        labels = assign_clusters(X, μ)
        μ = update_centroids(X, labels, K)
        all(norm(μ[k] .- μ_old[k]) < tol for k in 1:K) && break
    end
    return labels
end

K = length(unique(comunidad_louvain))
L = Matrix(laplacian_matrix(g))
autoval, autovec = eigen(Symmetric(Float64.(L)))
# Se omite el primer autovector (autovalor ≈0, constante: no separa nada,
# es el "modo trivial" del laplaciano de un grafo conexo). Se usan los
# siguientes K componentes como coordenadas espectrales de cada nodo
# (relajación del corte espectral normalizado, Newman 2010, cap. 11).
X_espectral = autovec[:, 2:K+1]
etiquetas_kmeans = my_kmeans(X_espectral, K)

valor_nmi_km = nmi(etiquetas_kmeans, comunidad_louvain)
valor_ari_km = ari(etiquetas_kmeans, comunidad_louvain)
@printf("k-means espectral (K=%d) vs. Louvain: NMI=%.4f  ARI=%.4f\n", K, valor_nmi_km, valor_ari_km)
println("""
El espacio espectral usa distancia EUCLÍDEA entre coordenadas derivadas del
laplaciano; Louvain optimiza directamente la modularidad sobre la
TOPOLOGÍA (conteo de aristas), sin necesitar ninguna noción de distancia
entre nodos. Cuando el grafo es muy disperso y jerárquico como este (grado
medio 2.36), los primeros K autovectores capturan sobre todo las ramas más
grandes del árbol de expansión, no necesariamente los mismos bloques que
Louvain encuentra por densidad local de aristas: coincidencias parciales
(NMI/ARI intermedios) son el resultado esperable, no un error.
""")

# =================================================================
# P4.5 -- Límite de resolución de la modularidad
# =================================================================
println("--- P4.5 Límite de resolución ---")
tam_comunidades = combine(groupby(DataFrame(c=comunidad_louvain), :c), nrow => :tamaño)
@printf("Tamaño de comunidades de Louvain: mínimo=%d, máximo=%d, mediana=%.0f\n",
        minimum(tam_comunidades.tamaño), maximum(tam_comunidades.tamaño), median(tam_comunidades.tamaño))
println("""
La modularidad tiene un límite de resolución conocido (Fortunato & Barthélemy,
2007): con m=$(m) aristas, dos comunidades solo se separan si la suma de sus
grados internos supera aproximadamente √(2m) ≈ $(round(sqrt(2m); digits=1)).
Bloques pequeños y débilmente conectados a la red (p.ej. un sub-laboratorio
con 3-4 equipos de acceso, todos colgando de un mismo switch) pueden
fusionarse con un bloque vecino más grande AUNQUE un administrador de red
-que sí distingue subredes/VLANs- los consideraría dominios separados. La
tabla de tamaños de arriba (y la lista completa en p4_louvain_estabilidad.csv)
es el punto de partida para revisar si eso ocurrió aquí.
""")

println("\n" * "="^70)
println("FASE 2 completa. Resultados en ", DIR_RESULTADOS)
println("="^70)
