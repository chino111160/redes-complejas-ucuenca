#=============================================================
  fase4_p8_p9_p10_percolacion.jl -- Fase 4: Percolación y Robustez
  Proyecto integrador de Redes Complejas (1217) -- Universidad de Cuenca

  Resuelve P8 (percolación), P9 (cascadas y epidemias) y P10 (diagnóstico
  de puntos críticos) del enunciado (proyecto_red_ucuenca.pdf, cap. 6).

  Uso:
      julia --project=analisis analisis/scripts/fase4_p8_p9_p10_percolacion.jl
=============================================================#

using Graphs
using DataFrames
using CSV
using Statistics
using StatsBase
using Printf
using Random
using Plots

const DIR_ESTE_SCRIPT = @__DIR__
const DIR_ANALISIS = dirname(DIR_ESTE_SCRIPT)
const DIR_RAIZ = dirname(DIR_ANALISIS)

include(joinpath(DIR_RAIZ, "codigo_base", "cargar_red.jl"))
include(joinpath(DIR_ANALISIS, "src", "ModeloRed.jl"))

const DIR_RESULTADOS = joinpath(DIR_ANALISIS, "resultados")
const DIR_FIGURAS = joinpath(DIR_ANALISIS, "figuras")
mkpath(DIR_RESULTADOS); mkpath(DIR_FIGURAS)

println("="^70)
println("FASE 4 -- Percolación y Robustez (P8-P10)")
println("="^70)

red = cargar_red()
verificar(red) || error("La carga no reproduce el Anexo A: revise el conjunto de datos.")
g = red.g
n, m = nv(g), ne(g)

# =================================================================
# P8.1-P8.2 -- Percolación de nodos, 4 estrategias
# =================================================================
println("\n--- P8.1-P8.2 Percolación de nodos ---")

"Tamaño relativo de la componente gigante entre los nodos aún vivos (según `vivo`, alineado 1:n)."
function tam_componente_gigante(g::SimpleGraph, vivo::BitVector, n::Int)
    idx = findall(vivo)
    isempty(idx) && return 0.0
    sg, _ = induced_subgraph(g, idx)
    comps = connected_components(sg)
    return maximum(length.(comps)) / n
end

"S(f) para una secuencia FIJA de eliminación (nodos en el orden en que se eliminan)."
function percolacion_nodos(g::SimpleGraph, orden::Vector{Int})
    n = nv(g)
    vivo = trues(n)
    S = zeros(Float64, n + 1)
    S[1] = 1.0
    for k in 1:n
        vivo[orden[k]] = false
        S[k+1] = tam_componente_gigante(g, vivo, n)
    end
    return S
end

"S(f) recalculando el criterio (p.ej. intermediación) sobre la red YA REDUCIDA antes de cada eliminación."
function percolacion_nodos_recalculada(g::SimpleGraph, criterio::Function)
    n = nv(g)
    vivo = trues(n)
    S = zeros(Float64, n + 1)
    S[1] = 1.0
    for k in 1:n
        idx_vivos = findall(vivo)
        sg, mapa = induced_subgraph(g, idx_vivos)
        valores = criterio(sg)
        peor_local = argmax(valores)
        vivo[mapa[peor_local]] = false
        S[k+1] = tam_componente_gigante(g, vivo, n)
    end
    return S
end

Random.seed!(1217)
const N_REAL_PERCOLACION = 100
S_aleatorio_muestras = zeros(N_REAL_PERCOLACION, n + 1)
for r in 1:N_REAL_PERCOLACION
    S_aleatorio_muestras[r, :] = percolacion_nodos(g, shuffle(1:n))
end
S_aleatorio_media = vec(mean(S_aleatorio_muestras; dims=1))
S_aleatorio_std = vec(std(S_aleatorio_muestras; dims=1))

grados_p8 = degree(g)
orden_grado = sortperm(grados_p8; rev=true)
S_grado = percolacion_nodos(g, orden_grado)

intermediacion_p8 = betweenness_centrality(g)
orden_intermediacion = sortperm(intermediacion_p8; rev=true)
S_intermediacion_estatica = percolacion_nodos(g, orden_intermediacion)

println("Calculando percolación con intermediación RECALCULADA (recomputa en cada paso; puede tardar unos segundos)...")
S_intermediacion_recalculada = percolacion_nodos_recalculada(g, sg -> betweenness_centrality(sg))

f_vals = (0:n) ./ n
function estimar_fc(S::Vector{Float64}, f::AbstractVector; umbral::Float64=0.5)
    k = findfirst(<(umbral), S)
    return isnothing(k) ? NaN : f[k]
end
fc_aleatorio = estimar_fc(S_aleatorio_media, f_vals)
fc_grado = estimar_fc(S_grado, f_vals)
fc_intermediacion_estatica = estimar_fc(S_intermediacion_estatica, f_vals)
fc_intermediacion_recalculada = estimar_fc(S_intermediacion_recalculada, f_vals)
@printf("f_c (S cae por debajo de 0.5): aleatorio=%.3f  grado=%.3f  intermediación(estática)=%.3f  intermediación(recalculada)=%.3f\n",
        fc_aleatorio, fc_grado, fc_intermediacion_estatica, fc_intermediacion_recalculada)

df_percolacion_nodos = DataFrame(
    f=f_vals, S_aleatorio_media=S_aleatorio_media, S_aleatorio_std=S_aleatorio_std,
    S_grado=S_grado, S_intermediacion_estatica=S_intermediacion_estatica,
    S_intermediacion_recalculada=S_intermediacion_recalculada,
)
guardar_csv(df_percolacion_nodos, joinpath(DIR_RESULTADOS, "p8_percolacion_nodos.csv"))

plt_perc_nodos = plot(f_vals, S_aleatorio_media; ribbon=S_aleatorio_std, label="Aleatorio (media±std, 100 real.)",
                       xlabel="fracción de nodos eliminados (f)", ylabel="S(f) -- componente gigante (relativa)",
                       title="Percolación de nodos -- red UCuenca", legend=:topright, lw=2,
                       left_margin=6Plots.mm, top_margin=3Plots.mm)
plot!(plt_perc_nodos, f_vals, S_grado; label="Ataque por grado", lw=2)
plot!(plt_perc_nodos, f_vals, S_intermediacion_estatica; label="Ataque por intermediación (estática)", lw=2)
plot!(plt_perc_nodos, f_vals, S_intermediacion_recalculada; label="Ataque por intermediación (recalculada)", lw=2)
savefig(plt_perc_nodos, joinpath(DIR_FIGURAS, "p8_percolacion_nodos.png"))

println("""
Paradoja de la robustez: si f_c(aleatorio) >> f_c(ataques dirigidos), la red
tolera mucho mejor fallos al azar que ataques deliberados -- el patrón
típico de redes con grado heterogéneo (unos pocos switches de core/agregación
concentran muchísimas conexiones). Consecuencia operativa: un plan de
MANTENIMIENTO puede tratar los fallos como si fueran aproximadamente
aleatorios (tolerancia razonable); un plan de RESPUESTA A INCIDENTES DE
SEGURIDAD no puede asumir eso -- un atacante que apunta a los nodos de
mayor grado/intermediación colapsa la red con muchas menos bajas que el
peor caso aleatorio, así que esos nodos necesitan protección desproporcionada
(redundancia, monitoreo, control de acceso físico) respecto a su número.
""")

# =================================================================
# P8.3 -- Percolación de enlaces, incluyendo ataque a los puentes de P1
# =================================================================
println("--- P8.3 Percolación de enlaces ---")
_, puentes_p1 = puentes_y_articulacion(g)

"S(f) eliminando aristas en el orden dado (lista de (u,v)); no relabела vértices, solo quita aristas."
function percolacion_enlaces(g::SimpleGraph, orden_aristas::Vector{Tuple{Int,Int}})
    n = nv(g)
    h = SimpleGraph(n)
    for e in edges(g)
        add_edge!(h, src(e), dst(e))
    end
    mm = length(orden_aristas)
    S = zeros(Float64, mm + 1)
    S[1] = maximum(length.(connected_components(h))) / n
    for k in 1:mm
        rem_edge!(h, orden_aristas[k][1], orden_aristas[k][2])
        S[k+1] = maximum(length.(connected_components(h))) / n
    end
    return S
end

lista_aristas = [(src(e), dst(e)) for e in edges(g)]
Random.seed!(2024)
S_enlaces_aleatorio_muestras = zeros(N_REAL_PERCOLACION, m + 1)
for r in 1:N_REAL_PERCOLACION
    S_enlaces_aleatorio_muestras[r, :] = percolacion_enlaces(g, shuffle(lista_aristas))
end
S_enlaces_aleatorio_media = vec(mean(S_enlaces_aleatorio_muestras; dims=1))

# Ataque dirigido: primero TODOS los puentes de P1 (en orden aleatorio entre
# ellos, porque quitar un puente no cambia que los demás lo sigan siendo),
# después el resto de las aristas ordenadas por intermediación de arista
# descendente (estática).
function intermediacion_aristas(g::SimpleGraph)
    n = nv(g)
    ib = Dict{Tuple{Int,Int},Float64}()
    for i in 1:n
        d, padre = let
            dist = fill(-1, n); padre_ = zeros(Int, n); dist[i] = 0
            cola = [i]
            while !isempty(cola)
                u = popfirst!(cola)
                for v in neighbors(g, u)
                    if dist[v] == -1
                        dist[v] = dist[u] + 1; padre_[v] = u; push!(cola, v)
                    end
                end
            end
            dist, padre_
        end
        for v in 1:n
            v == i && continue
            u = v
            while padre[u] != 0
                clave = padre[u] < u ? (padre[u], u) : (u, padre[u])
                ib[clave] = get(ib, clave, 0.0) + 1.0
                u = padre[u]
            end
        end
    end
    return ib
end
ib = intermediacion_aristas(g)
puentes_set = Set([(min(u, v), max(u, v)) for (u, v) in puentes_p1])
resto_ordenado = sort([e for e in lista_aristas if !((min(e[1], e[2]), max(e[1], e[2])) in puentes_set)];
                       by=e -> -get(ib, (min(e[1], e[2]), max(e[1], e[2])), 0.0))
Random.seed!(7)
orden_puentes_primero = vcat(shuffle(collect(puentes_set)), resto_ordenado)
S_enlaces_puentes = percolacion_enlaces(g, orden_puentes_primero)

f_vals_enlaces = (0:m) ./ m
df_percolacion_enlaces = DataFrame(f=f_vals_enlaces, S_aleatorio_media=S_enlaces_aleatorio_media,
                                    S_ataque_puentes_primero=S_enlaces_puentes)
guardar_csv(df_percolacion_enlaces, joinpath(DIR_RESULTADOS, "p8_percolacion_enlaces.csv"))

plt_perc_enlaces = plot(f_vals_enlaces, S_enlaces_aleatorio_media; label="Aleatorio (media, 100 real.)",
                         xlabel="fracción de enlaces eliminados (f)", ylabel="S(f)",
                         title="Percolación de enlaces -- red UCuenca", lw=2)
plot!(plt_perc_enlaces, f_vals_enlaces, S_enlaces_puentes; label="Puentes primero, luego intermediación de arista", lw=2)
savefig(plt_perc_enlaces, joinpath(DIR_FIGURAS, "p8_percolacion_enlaces.png"))
@printf("f_c de enlaces: aleatorio=%.3f   puentes-primero=%.3f (con %d puentes de %d aristas, %.1f%% del total)\n",
        estimar_fc(S_enlaces_aleatorio_media, f_vals_enlaces), estimar_fc(S_enlaces_puentes, f_vals_enlaces),
        length(puentes_p1), m, 100 * length(puentes_p1) / m)

# =================================================================
# P8.4 -- Eficiencia global E(f)
# =================================================================
println("\n--- P8.4 Eficiencia global E(f) ---")
"Eficiencia global E = (1/(n(n-1))) Σ_{i≠j} 1/d(i,j), sobre los nodos vivos (d=Inf entre componentes -> aporta 0)."
function eficiencia_global(g::SimpleGraph, vivo::BitVector, n_total::Int)
    idx = findall(vivo)
    length(idx) < 2 && return 0.0
    sg, _ = induced_subgraph(g, idx)
    total = 0.0
    for u in 1:nv(sg)
        d = gdistances(sg, u)
        for v in 1:nv(sg)
            v == u && continue
            d[v] > 0 && (total += 1.0 / d[v])
        end
    end
    return total / (n_total * (n_total - 1))
end

function curva_eficiencia(g::SimpleGraph, orden::Vector{Int})
    n = nv(g)
    vivo = trues(n)
    E = zeros(Float64, n + 1)
    E[1] = eficiencia_global(g, vivo, n)
    for k in 1:n
        vivo[orden[k]] = false
        E[k+1] = eficiencia_global(g, vivo, n)
    end
    return E
end

# CORRECCIÓN (auditoría de reproducibilidad): igual que en P8.5, un solo
# shuffle(1:n) es una única muestra aleatoria con varianza alta para esta
# métrica; se promedia sobre N_EFICIENCIA réplicas, igual que S_aleatorio_media.
const N_EFICIENCIA = 20
Random.seed!(2020)
E_aleatorio_muestras = zeros(N_EFICIENCIA, n + 1)
for r in 1:N_EFICIENCIA
    E_aleatorio_muestras[r, :] = curva_eficiencia(g, shuffle(1:n))
end
E_aleatorio = vec(mean(E_aleatorio_muestras; dims=1))
E_intermediacion = curva_eficiencia(g, orden_intermediacion)
df_eficiencia = DataFrame(f=f_vals, E_aleatorio=E_aleatorio, E_ataque_intermediacion=E_intermediacion)
guardar_csv(df_eficiencia, joinpath(DIR_RESULTADOS, "p8_eficiencia_global.csv"))

plt_eficiencia = plot(f_vals, E_aleatorio ./ E_aleatorio[1]; label="Aleatorio", lw=2,
                       xlabel="fracción de nodos eliminados (f)", ylabel="E(f)/E(0)",
                       title="Eficiencia global relativa vs. percolación de nodos")
plot!(plt_eficiencia, f_vals, E_intermediacion ./ E_intermediacion[1]; label="Ataque por intermediación", lw=2)
plot!(plt_eficiencia, f_vals, S_aleatorio_media; label="S(f) aleatorio (referencia)", lw=1, linestyle=:dash, color=:gray)
plot!(plt_eficiencia, f_vals, S_intermediacion_estatica; label="S(f) ataque (referencia)", lw=1, linestyle=:dash, color=:black)
savefig(plt_eficiencia, joinpath(DIR_FIGURAS, "p8_eficiencia_vs_percolacion.png"))
println("""
E(f) puede caer mucho antes de que S(f) (componente gigante) colapse porque
la eficiencia pondera el INVERSO de la distancia: un nodo que sigue
"técnicamente conectado" pero ahora a 8 saltos en vez de a 2 aporta casi
nada a E aunque siga contando en la componente gigante. Para el diagnóstico
operativo esto importa: la red puede seguir "arriba" (conectada) mucho
después de que el SERVICIO percibido (latencia efectiva) ya se haya
degradado severamente.
""")

# =================================================================
# P8.5 -- Comparación con los modelos nulos de P2
# =================================================================
# CORRECCIÓN (auditoría de reproducibilidad): la versión anterior generaba
# un solo Erdős-Rényi y un solo modelo de configuración (una única muestra
# aleatoria), con varianza demasiado alta para reportar f_c como un punto
# fijo -- inconsistente con el resto de P8, que sí promedia (100 réplicas
# para el fallo aleatorio de nodos) y con P2 (100 réplicas de cada modelo
# nulo). Aquí se promedia sobre N_NULOS_P8 réplicas de cada modelo, cada
# una con su propio orden aleatorio de eliminación, igual que S_aleatorio.
println("--- P8.5 Comparación con modelos nulos ---")
grados_reales = collect(grados_p8)
const N_NULOS_P8 = 50
Random.seed!(8500)
S_er_muestras = zeros(N_NULOS_P8, n + 1)
S_conf_muestras = zeros(N_NULOS_P8, n + 1)
fc_er_muestras = zeros(N_NULOS_P8)
fc_conf_muestras = zeros(N_NULOS_P8)
for r in 1:N_NULOS_P8
    g_er_r = erdos_renyi(n, m)
    g_conf_r = modelo_configuracion(grados_reales)
    S_er_r = percolacion_nodos(g_er_r, shuffle(1:n))
    S_conf_r = percolacion_nodos(g_conf_r, shuffle(1:n))
    S_er_muestras[r, :] = S_er_r
    S_conf_muestras[r, :] = S_conf_r
    fc_er_muestras[r] = estimar_fc(S_er_r, f_vals)
    fc_conf_muestras[r] = estimar_fc(S_conf_r, f_vals)
end
S_er_media = vec(mean(S_er_muestras; dims=1))
S_conf_media = vec(mean(S_conf_muestras; dims=1))
fc_er_media, fc_er_std = mean(fc_er_muestras), std(fc_er_muestras)
fc_conf_media, fc_conf_std = mean(fc_conf_muestras), std(fc_conf_muestras)
df_percolacion_nulos = DataFrame(f=f_vals, S_real=S_aleatorio_media, S_ER=S_er_media, S_configuracion=S_conf_media)
guardar_csv(df_percolacion_nulos, joinpath(DIR_RESULTADOS, "p8_percolacion_vs_nulos.csv"))
@printf("f_c bajo fallo aleatorio (promedio de %d réplicas por modelo): red real=%.3f   Erdős-Rényi=%.3f±%.3f   configuración=%.3f±%.3f\n",
        N_NULOS_P8, fc_aleatorio, fc_er_media, fc_er_std, fc_conf_media, fc_conf_std)
if fc_aleatorio < fc_conf_media - fc_conf_std
    println("f_c(real) cae POR DEBAJO de f_c(configuración) incluso preservando exactamente la secuencia",
            "\nde grados: la robustez ante fallos aleatorios NO se explica solo por la heterogeneidad de",
            "\ngrados -- el cableado específico de UCuenca concentra fragilidad de forma que el recableado",
            "\naleatorio (con los mismos grados) tiende a evitar. Consistente con el hallazgo de P6",
            "\n(Paraíso depende de un único enlace) y con la asimetría de redundancia ya señalada en P1.")
else
    println("Si f_c(real) ≈ f_c(configuración) « f_c(ER), la robustez de UCuenca ante fallos aleatorios",
            "\nse explica por su secuencia de grados heterogénea (unos pocos hubs, muchas hojas) y NO por",
            "\nalgún otro patrón estructural adicional -- coherente con lo ya visto en P1/P2.")
end

# =================================================================
# P9.1-P9.2 -- Cascadas de fallos, margen crítico τ_c
# =================================================================
println("\n" * "="^70); println("P9 -- Propagación de fallos y epidemias"); println("="^70)
println("\n--- P9.1-P9.2 Cascadas de fallos ---")

L0 = betweenness_centrality(g; normalize=false)  # carga de referencia = intermediación SIN normalizar (conteo de caminos más cortos)

"""
    simular_cascada(g, nodo_inicial, τ) -> Vector{Int}

Modelo de Motter & Lai (2002): capacidad fija Cᵢ=(1+τ)Lᵢ (Lᵢ=carga INICIAL,
intermediación no normalizada). Al eliminar un nodo, la carga se
redistribuye por los caminos más cortos de la red reducida (se recalcula la
intermediación ahí); cualquier nodo cuya carga recalculada supere su
capacidad FIJA también falla, y la cascada continúa hasta que nadie más
supera su capacidad.
"""
function simular_cascada(g::SimpleGraph, nodo_inicial::Int, τ::Float64, capacidad_nodo::Vector{Float64})
    n = nv(g)
    vivo = trues(n)
    vivo[nodo_inicial] = false
    caidos = Set([nodo_inicial])
    cambiado = true
    while cambiado
        cambiado = false
        idx_vivos = findall(vivo)
        length(idx_vivos) < 2 && break
        sg, mapa = induced_subgraph(g, idx_vivos)
        carga_actual = betweenness_centrality(sg; normalize=false)
        for (i_local, i_orig) in enumerate(mapa)
            if carga_actual[i_local] > capacidad_nodo[i_orig]
                vivo[i_orig] = false
                push!(caidos, i_orig)
                cambiado = true
            end
        end
    end
    return caidos
end

# Candidatos a disparador: los de mayor carga inicial -- son los únicos con
# capacidad de desencadenar una cascada relevante (un nodo hoja, al caer,
# no redistribuye carga significativa hacia nadie). Restringir a estos
# candidatos mantiene la búsqueda de τ_c computacionalmente tratable.
const N_CANDIDATOS_CASCADA = 15
candidatos = sortperm(L0; rev=true)[1:N_CANDIDATOS_CASCADA]
# Grilla de τ ampliada respecto de la versión original (que partía en 0,05):
# se agregan 0,01/0,02/0,03 para acotar más finamente el umbral justo por
# encima del piso del modelo. τ < 0 no se explora porque no es físicamente
# interpretable aquí: C_i=(1+τ)L_i con τ<0 implicaría que un nodo arranca
# YA por encima de su propia capacidad (falla en t=0 sin ninguna perturbación),
# lo cual no es un "margen de tolerancia" sino un estado inicial inválido.
# τ=0 (capacidad exactamente igual a la carga inicial, el caso más frágil
# posible dentro del modelo) ya estaba cubierto y sigue siendo el piso real.
grilla_τ = [0.0, 0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0]

# Criterio de daño alternativo: además de la fracción de NODOS caídos
# (criterio original), se mide la fracción de CAPACIDAD (Mbps) de la red
# que queda incidente a algún nodo caído -- un criterio operativo distinto,
# consistente con que el resto del informe (P6) mide todo en Mbps. Los dos
# criterios pueden discrepar: pocos nodos de alta capacidad caídos bastan
# para dañar mucha capacidad, o muchos nodos de acceso de baja capacidad
# pueden caer sin comprometer mucho ancho de banda total.
capacidad_arista_p9, _ = estimar_capacidad(red)
capacidad_total_red = sum(capacidad_arista_p9)
aristas_por_nodo_p9 = [Int[] for _ in 1:n]
for (k, (u, v)) in enumerate(zip(red.aristas.u, red.aristas.v))
    push!(aristas_por_nodo_p9[u], k)
    push!(aristas_por_nodo_p9[v], k)
end
"Fracción de la capacidad total de la red (Mbps) incidente a algún nodo del conjunto `caidos`."
function daño_capacidad(caidos)
    idxs_arista = Set{Int}()
    for i in caidos, k in aristas_por_nodo_p9[i]
        push!(idxs_arista, k)
    end
    return sum(capacidad_arista_p9[k] for k in idxs_arista; init=0.0) / capacidad_total_red
end

filas_cascada = NamedTuple[]
for τ in grilla_τ
    capacidad_nodo = (1 + τ) .* L0
    for cand in candidatos
        caidos = simular_cascada(g, cand, τ, capacidad_nodo)
        push!(filas_cascada, (τ=τ, nodo=red.ids[cand], caidos=length(caidos),
                               fraccion_afectada=length(caidos) / n,
                               fraccion_capacidad_afectada=daño_capacidad(caidos)))
    end
end
df_cascada = DataFrame(filas_cascada)
guardar_csv(df_cascada, joinpath(DIR_RESULTADOS, "p9_cascadas_barrido_tau.csv"))

const UMBRAL_CASCADA_GRAVE = 0.2
τ_c_por_nodo = combine(groupby(df_cascada, :nodo)) do sub
    peligrosos_nodos = sub[sub.fraccion_afectada .> UMBRAL_CASCADA_GRAVE, :τ]
    peligrosos_cap = sub[sub.fraccion_capacidad_afectada .> UMBRAL_CASCADA_GRAVE, :τ]
    (τ_c=isempty(peligrosos_nodos) ? missing : maximum(peligrosos_nodos),
     max_fraccion_observada=maximum(sub.fraccion_afectada),
     τ_c_capacidad=isempty(peligrosos_cap) ? missing : maximum(peligrosos_cap),
     max_fraccion_capacidad_observada=maximum(sub.fraccion_capacidad_afectada))
end
sort!(τ_c_por_nodo, :max_fraccion_observada; rev=true)
guardar_csv(τ_c_por_nodo, joinpath(DIR_RESULTADOS, "p9_tau_critico_por_nodo.csv"))
existe_disparador = any(!ismissing, τ_c_por_nodo.τ_c)
existe_disparador_capacidad = any(!ismissing, τ_c_por_nodo.τ_c_capacidad)
# τ_c tiene sentido ("el margen por debajo del cual UN nodo dispara una
# cascada grave") solo si algún candidato lo logra en la grilla probada; si
# ninguno lo hace ni siquiera a margen cero, no hay un τ_c positivo que
# reportar y hay que decirlo explícitamente en vez de imprimir 0.0 (que se
# leería, al revés, como "la red es frágil incluso sin margen"). Se repite
# el mismo razonamiento para el criterio alternativo de capacidad.
τ_c_red = existe_disparador ? maximum(skipmissing(τ_c_por_nodo.τ_c)) : 0.0
τ_c_red_capacidad = existe_disparador_capacidad ? maximum(skipmissing(τ_c_por_nodo.τ_c_capacidad)) : 0.0
if existe_disparador
    @printf("Margen crítico τ_c de la red (criterio nodos, peor nodo candidato): %.2f\n", τ_c_red)
else
    peor = first(τ_c_por_nodo)
    @printf("Criterio NODOS: ningún nodo candidato desencadena una cascada > %.0f%% de la red, ni siquiera con margen τ=0.\n",
            100 * UMBRAL_CASCADA_GRAVE)
    @printf("  Peor caso observado en la grilla: %s alcanza %.1f%% de la red afectada (máximo entre los τ probados).\n",
            peor.nodo, 100 * peor.max_fraccion_observada)
end
if existe_disparador_capacidad
    @printf("Margen crítico τ_c de la red (criterio capacidad Mbps, peor nodo candidato): %.2f\n", τ_c_red_capacidad)
else
    peor_cap = sort(τ_c_por_nodo, :max_fraccion_capacidad_observada; rev=true) |> first
    @printf("Criterio CAPACIDAD: ningún nodo candidato desencadena una cascada que comprometa > %.0f%% de la capacidad (Mbps) de la red, ni siquiera con margen τ=0.\n",
            100 * UMBRAL_CASCADA_GRAVE)
    @printf("  Peor caso observado en la grilla: %s compromete %.1f%% de la capacidad de la red (máximo entre los τ probados).\n",
            peor_cap.nodo, 100 * peor_cap.max_fraccion_capacidad_observada)
end
if !existe_disparador && !existe_disparador_capacidad
    println("Los DOS criterios de daño (fracción de nodos, fracción de capacidad Mbps) coinciden: ninguno",
            "\nencuentra un τ_c positivo en la grilla ampliada. Esto refuerza -- por una vía adicional e",
            "\nindependiente -- que la conclusión no es un artefacto de cómo se mide el daño.")
    println("Esto NO contradice la fragilidad hallada en P8: ahí la métrica es CONECTIVIDAD (aristas que",
            "\ndesconectan partes de la red); acá es SOBRECARGA (redistribución de camino más corto). Esta red",
            "\ntiene poca redundancia, pero las rutas alternativas tras una falla son mayormente CORTAS y NO",
            "\nconcentran suficiente tráfico redistribuido en un tercer nodo como para tumbarlo en cadena.")
end
println("Nodos disparadores más peligrosos (mayor daño máximo observado en la grilla de τ, criterio nodos):")
println(first(τ_c_por_nodo, 5))

# =================================================================
# P9.3 -- SIR y umbral epidémico
# =================================================================
println("\n--- P9.3 SIR y umbral epidémico ---")
"""
    simular_sir(g, semilla, β, γ; inmunes) -> tamaño final del brote

SIR discreto de tiempo síncrono: cada infectado contagia a cada vecino
susceptible con probabilidad β por paso, y se recupera con probabilidad γ
por paso. Los nodos en `inmunes` empiezan en R (no pueden contagiarse ni
contagiar).
"""
function simular_sir(g::SimpleGraph, semilla::Int, β::Float64, γ::Float64;
                      inmunes::Set{Int}=Set{Int}(), max_pasos::Int=300)
    n = nv(g)
    estado = fill(:S, n)
    for i in inmunes
        estado[i] = :R
    end
    estado[semilla] == :R && return 0
    estado[semilla] = :I
    for _ in 1:max_pasos
        infectados = findall(==(:I), estado)
        isempty(infectados) && break
        nuevas = Int[]
        for i in infectados, v in neighbors(g, i)
            estado[v] == :S && rand() < β && push!(nuevas, v)
        end
        for i in infectados
            rand() < γ && (estado[i] = :R)
        end
        for v in unique(nuevas)
            estado[v] == :S && (estado[v] = :I)
        end
    end
    return count(x -> x == :R || x == :I, estado)
end

grados_g = degree(g)
k_medio = mean(grados_g)
k2_medio = mean(grados_g .^ 2)
λ_c = k_medio / k2_medio
@printf("Predicción de campo medio: <k>=%.3f  <k²>=%.3f  λ_c=<k>/<k²>=%.4f\n", k_medio, k2_medio, λ_c)

const γ_SIR = 0.3
grilla_λ = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0] .* λ_c
const N_REAL_SIR = 60
Random.seed!(99)
filas_sir = NamedTuple[]
for λ in grilla_λ
    β = λ * γ_SIR
    tamaños = Float64[]
    for _ in 1:N_REAL_SIR
        semilla = rand(1:n)
        push!(tamaños, simular_sir(g, semilla, β, γ_SIR) / n)
    end
    push!(filas_sir, (λ=λ, λ_sobre_λc=λ / λ_c, tamaño_medio=mean(tamaños), tamaño_std=std(tamaños)))
end
df_sir = DataFrame(filas_sir)
guardar_csv(df_sir, joinpath(DIR_RESULTADOS, "p9_sir_umbral.csv"))
println(df_sir)

plt_sir = plot(df_sir.λ_sobre_λc, df_sir.tamaño_medio; yerror=df_sir.tamaño_std, marker=:circle,
               xlabel="λ/λ_c", ylabel="tamaño final del brote (fracción de la red)",
               title="SIR -- umbral epidémico", label="simulación (media±std, $N_REAL_SIR real.)")
vline!(plt_sir, [1.0]; label="λ_c (campo medio)", linestyle=:dash, color=:red)
savefig(plt_sir, joinpath(DIR_FIGURAS, "p9_sir_umbral.png"))

# =================================================================
# P9.4 -- Inmunización: aleatoria vs. por centralidad
# =================================================================
println("\n--- P9.4 Estrategias de inmunización ---")
const M_INMUNIZADOS = max(5, round(Int, 0.06 * n))  # ~6% de la red ("presupuesto" de parches)
# λ=2λ_c (usado en el barrido de P9.3) queda demasiado cerca del umbral: ahí
# el tamaño del brote es dominado por ruido de qué tan lejos cae la semilla
# de la componente epidémica, no por la inmunización -- con pocas
# realizaciones eso puede incluso mostrar una reducción "negativa" por puro
# azar. Para que la comparación sea legible se usa un λ bien supercrítico y
# más realizaciones (menos ruido, efecto de la inmunización dominante).
λ_prueba = 8.0 * λ_c
β_prueba = λ_prueba * γ_SIR
const N_REAL_INMUNIZACION = 200

orden_centralidad_inmunizacion = sortperm(intermediacion_p8; rev=true)[1:M_INMUNIZADOS]
inmunes_centralidad = Set(orden_centralidad_inmunizacion)

Random.seed!(555)
tam_sin_inmunizar = [simular_sir(g, rand(1:n), β_prueba, γ_SIR) / n for _ in 1:N_REAL_INMUNIZACION]
tam_inmune_aleatorio = Float64[]
for _ in 1:N_REAL_INMUNIZACION
    inmunes_al = Set(sample(1:n, M_INMUNIZADOS; replace=false))
    semilla = rand(setdiff(1:n, inmunes_al))
    push!(tam_inmune_aleatorio, simular_sir(g, semilla, β_prueba, γ_SIR; inmunes=inmunes_al) / n)
end
tam_inmune_centralidad = Float64[]
for _ in 1:N_REAL_INMUNIZACION
    semilla = rand(setdiff(1:n, inmunes_centralidad))
    push!(tam_inmune_centralidad, simular_sir(g, semilla, β_prueba, γ_SIR; inmunes=inmunes_centralidad) / n)
end

df_inmunizacion = DataFrame(
    estrategia=["sin inmunizar", "aleatoria ($M_INMUNIZADOS nodos)", "por centralidad ($M_INMUNIZADOS nodos)"],
    tamaño_medio_brote=[mean(tam_sin_inmunizar), mean(tam_inmune_aleatorio), mean(tam_inmune_centralidad)],
    tamaño_std=[std(tam_sin_inmunizar), std(tam_inmune_aleatorio), std(tam_inmune_centralidad)],
)
guardar_csv(df_inmunizacion, joinpath(DIR_RESULTADOS, "p9_inmunizacion.csv"))
println(df_inmunizacion)
@printf("Reducción del tamaño medio del brote: aleatoria=%.1f%%   por centralidad=%.1f%%\n",
        100 * (1 - mean(tam_inmune_aleatorio) / mean(tam_sin_inmunizar)),
        100 * (1 - mean(tam_inmune_centralidad) / mean(tam_sin_inmunizar)))

# =================================================================
# P9.5 -- Analogía con sistemas de energía eléctrica
# =================================================================
println("""
--- P9.5 Analogía con redes de transmisión eléctrica ---
En una red de transmisión, la "carga" de una línea es el flujo de potencia
que transporta (determinado por las leyes de Kirchhoff, no por caminos más
cortos, pero igualmente concentrado en unas pocas líneas troncales); la
"capacidad" es el límite térmico/de estabilidad de esa línea; la "cascada"
es exactamente el mecanismo de un apagón en cascada (p.ej. el apagón del
noreste de EE. UU.-Canadá de 2003): al desconectarse una línea sobrecargada,
su flujo se redistribuye instantáneamente por las líneas restantes según
las mismas leyes físicas, que pueden entonces sobrecargarse y disparar la
siguiente desconexión. Motter & Lai (2002, "Cascade-based attacks on
complex networks", Phys. Rev. E) formalizan precisamente esta analogía y es
la referencia directa del modelo de cascada usado arriba.
""")

# =================================================================
# P10 -- Ranking de puntos críticos
# =================================================================
println("="^70); println("P10 -- Diagnóstico de puntos críticos"); println("="^70)

articulacion_p10, _ = puentes_y_articulacion(g)
es_articulacion = falses(n); es_articulacion[articulacion_p10] .= true

# Participación en cortes mínimos: cuántos de los cortes mínimos de P6
# (uno por campus) incluyen a este nodo como uno de sus dos extremos.
ruta_cortes = joinpath(DIR_RESULTADOS, "p6_cortes_minimos_aristas.csv")
participacion_corte = zeros(Int, n)
if isfile(ruta_cortes)
    df_cortes = CSV.read(ruta_cortes, DataFrame)
    for r in eachrow(df_cortes)
        haskey(red.idx, r.origen) && (participacion_corte[red.idx[r.origen]] += 1)
        haskey(red.idx, r.destino) && (participacion_corte[red.idx[r.destino]] += 1)
    end
else
    println("Aviso: no se encontró p6_cortes_minimos_aristas.csv (ejecute fase3 primero) -- participación en cortes queda en 0.")
end

# Daño en cascada: fracción de la red que cae si ESTE nodo falla, a un τ de
# referencia fijo (se usa el τ_c estimado en P9.2, la frontera de seguridad
# ya diagnosticada para la red).
τ_referencia = max(τ_c_red, 0.1)
capacidad_referencia = (1 + τ_referencia) .* L0
println("Calculando daño en cascada de cada nodo a τ=$(round(τ_referencia; digits=2)) (puede tardar unos segundos)...")
daño_cascada = zeros(Float64, n)
for i in 1:n
    daño_cascada[i] = length(simular_cascada(g, i, τ_referencia, capacidad_referencia)) / n
end

normalizar(v) = maximum(v) > 0 ? v ./ maximum(v) : v
indice = (normalizar(intermediacion_p8) .+ Float64.(es_articulacion) .+
          normalizar(Float64.(participacion_corte)) .+ normalizar(daño_cascada)) ./ 4

df_ranking = DataFrame(
    id=red.ids, campus=red.campus, capa=red.capa,
    intermediacion_norm=round.(normalizar(intermediacion_p8); digits=3),
    es_articulacion=es_articulacion,
    participacion_corte=participacion_corte,
    daño_cascada=round.(daño_cascada; digits=3),
    indice_compuesto=round.(indice; digits=4),
)
sort!(df_ranking, :indice_compuesto; rev=true)
guardar_csv(df_ranking, joinpath(DIR_RESULTADOS, "p10_ranking_puntos_criticos.csv"))

println("""
Índice compuesto = promedio simple (pesos iguales, punto de partida
razonable a falta de una prioridad de negocio declarada por la institución)
de 4 componentes normalizados a [0,1]: intermediación, condición de punto
de articulación (0/1), participación en los cortes mínimos de P6 y daño
causado en la cascada de P9. Top-10:
""")
top10_ranking = first(df_ranking, 10)
println(top10_ranking)
guardar_csv(top10_ranking, joinpath(DIR_RESULTADOS, "p10_top10_ficha.csv"))

for r in eachrow(top10_ranking)
    @printf("\n%-28s  campus=%-20s capa=%-12s  índice=%.3f\n", r.id, r.campus, r.capa, r.indice_compuesto)
    @printf("  intermediación_norm=%.2f  articulación=%s  participa en %d corte(s)  daño_cascada=%.2f\n",
            r.intermediacion_norm, r.es_articulacion ? "sí" : "no", r.participacion_corte, r.daño_cascada)
end

println("\n" * "="^70)
println("FASE 4 completa. Resultados en ", DIR_RESULTADOS, " y figuras en ", DIR_FIGURAS)
println("="^70)
