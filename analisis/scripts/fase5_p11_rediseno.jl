#=============================================================
  fase5_p11_rediseno.jl -- Fase 5: Propuesta de Rediseño (P11)
  Proyecto integrador de Redes Complejas (1217) -- Universidad de Cuenca

  Parte del diagnóstico real de las Fases 1-4 (lee
  resultados/p10_ranking_puntos_criticos.csv, p1_puentes.csv,
  p6_flujo_maximo_por_campus.csv) para proponer hasta 5 intervenciones,
  cuantificar su efecto sobre la red modificada y compararlas con dos
  alternativas ingenuas.

  Uso:
      julia --project=analisis analisis/scripts/fase5_p11_rediseno.jl
=============================================================#

using Graphs
using DataFrames
using CSV
using Statistics
using Printf
using Random
using Plots

const DIR_ESTE_SCRIPT = @__DIR__
const DIR_ANALISIS = dirname(DIR_ESTE_SCRIPT)
const DIR_RAIZ = dirname(DIR_ANALISIS)

include(joinpath(DIR_RAIZ, "codigo_base", "cargar_red.jl"))
include(joinpath(DIR_ANALISIS, "src", "ModeloRed.jl"))
module FF
    include(joinpath(@__DIR__, "..", "..", "codigo_referencia", "ford-fulkerson", "ford_fulkerson.jl"))
end

const DIR_RESULTADOS = joinpath(DIR_ANALISIS, "resultados")
const DIR_FIGURAS = joinpath(DIR_ANALISIS, "figuras")
mkpath(DIR_RESULTADOS); mkpath(DIR_FIGURAS)

println("="^70)
println("FASE 5 -- Propuesta de Rediseño (P11)")
println("="^70)

red = cargar_red()
verificar(red) || error("La carga no reproduce el Anexo A: revise el conjunto de datos.")
g0 = red.g
n = nv(g0)
cap0, _ = estimar_capacidad(red)

ruta_top10 = joinpath(DIR_RESULTADOS, "p10_top10_ficha.csv")
isfile(ruta_top10) || error("Ejecute fase4 primero: falta $ruta_top10")
top10 = CSV.read(ruta_top10, DataFrame)
println("\nDiagnóstico de partida (P10, top-10 completo en $ruta_top10):")
println(top10)

# =================================================================
# Métricas de referencia (idénticas a P1/P6/P8, reevaluadas aquí para
# poder recalcularlas sobre red base y red modificada con el mismo código)
# =================================================================
"Distancia media y diámetro por BFS exhaustivo (igual método que fase1)."
function distancia_media_y_diametro(g::SimpleGraph)
    nn = nv(g); suma = 0; npares = 0; dmax = 0
    for u in 1:nn
        d = gdistances(g, u)
        for v in (u+1):nn
            suma += d[v]; npares += 1; dmax = max(dmax, d[v])
        end
    end
    return suma / npares, dmax
end
function eficiencia_global(g::SimpleGraph)
    nn = nv(g); total = 0.0
    for u in 1:nn
        d = gdistances(g, u)
        for v in 1:nn
            v == u && continue
            d[v] > 0 && (total += 1.0 / d[v])
        end
    end
    return total / (nn * (nn - 1))
end
function auc_percolacion_grado(g::SimpleGraph)
    nn = nv(g)
    orden = sortperm(degree(g); rev=true)
    vivo = trues(nn)
    S = zeros(Float64, nn + 1); S[1] = 1.0
    for k in 1:nn
        vivo[orden[k]] = false
        idx = findall(vivo)
        S[k+1] = isempty(idx) ? 0.0 : maximum(length.(connected_components(induced_subgraph(g, idx)[1]))) / nn
    end
    return sum(S) / length(S)  # área bajo la curva -- resumen escalar de robustez ante ataque por grado
end

"""
Flujo máximo campus->INTERNET-MPLS con super-fuente, reutilizando
FF.ford_fulkerson (metodo BFS = Edmonds-Karp). `cap_dict` da la capacidad
por par de nodos (sin importar el orden `(u,v)`/`(v,u)`) -- se usa un
diccionario en vez de un vector posicional porque el orden en que
`edges(g)` itera las aristas NO coincide, en general, con el orden de
filas de `red.aristas` (verificado: son el mismo conjunto pero en orden
distinto); emparejar por POSICIÓN escondería una capacidad mal asignada.
"""
function flujo_maximo_campus(g::SimpleGraph, cap_dict::Dict{Tuple{Int,Int},Float64}, campus_objetivo::String, idx_sink::Int)
    nn = nv(g)
    N = nn + 1
    C = zeros(Int, N, N)
    for e in edges(g)
        u, v = src(e), dst(e)
        c = get(cap_dict, (u, v), get(cap_dict, (v, u), 0.0))
        C[u, v] = round(Int, c); C[v, u] = round(Int, c)
    end
    idx_acceso = findall(i -> red.capa[i] == "acceso" && red.campus[i] == campus_objetivo, 1:nn)
    for i in idx_acceso
        C[N, i] = 200_000
    end
    rf = FF.RedFlujo(C, vcat(red.ids, ["SUPERFUENTE"]), [(0.0, 0.0) for _ in 1:N])
    flujo, _, _ = FF.ford_fulkerson(rf, N, idx_sink; metodo=:bfs, verbose=false)
    return flujo
end

"Recalcula el capítulo completo de métricas (P11.2) sobre una red (grafo + diccionario de capacidades (u,v)=>Mbps) dada."
function evaluar_red(g::SimpleGraph, cap_dict::Dict{Tuple{Int,Int},Float64}, idx_sink::Int)
    articulacion, puentes = puentes_y_articulacion(g)
    dist_media, _ = distancia_media_y_diametro(g)
    efic = eficiencia_global(g)
    auc = auc_percolacion_grado(g)
    flujo_paraiso = flujo_maximo_campus(g, cap_dict, "Campus Paraiso", idx_sink)
    flujo_hospitalidad = flujo_maximo_campus(g, cap_dict, "Campus Hospitalidad", idx_sink)
    flujo_yanuncay = flujo_maximo_campus(g, cap_dict, "Campus Yanuncay", idx_sink)
    return (puentes=length(puentes), articulacion=length(articulacion), distancia_media=dist_media,
            eficiencia_global=efic, robustez_auc_percolacion=auc,
            flujo_paraiso_mbps=flujo_paraiso, flujo_hospitalidad_mbps=flujo_hospitalidad,
            flujo_yanuncay_mbps=flujo_yanuncay)
end

idx_sink = red.idx["INTERNET-MPLS"]
# Diccionario (u,v)=>capacidad para la red base, construido directamente
# desde red.aristas.u/.v (la fuente de verdad, alineada por fila con cap0)
# -- no desde edges(g), cuyo orden de iteración es distinto.
cap_dict0 = Dict((min(u, v), max(u, v)) => c for (u, v, c) in zip(red.aristas.u, red.aristas.v, cap0))
println("\nEvaluando la red BASE (puede tardar unos segundos por los flujos máximos)...")
metricas_base = evaluar_red(g0, cap_dict0, idx_sink)
println(metricas_base)

# =================================================================
# P11.1 -- Intervenciones propuestas (derivadas del diagnóstico real)
# =================================================================
println("\n--- P11.1 Intervenciones propuestas ---")
println("""
Del top-10 de P10 y de p6_flujo_maximo_por_campus.csv surgen dos patrones:
(a) Campus Paraíso tiene el flujo máximo MÁS BAJO de la institución (1000 Mbps,
    35 nodos de acceso) porque su único core (CPAR-C10) sale por un ÚNICO
    enlace de 1 Gbps hacia ROUTER-CAMPUS-HUAYNA-CAPAC -- el informe afirma
    redundancia en Balzay/Paraíso, pero P1/P6 ya mostraron que Paraíso en
    realidad no la tiene.
(b) Hospitalidad y el Consultorio Jurídico (HOS-0A-D05, CCJ-CJURIDICO-D4)
    dependen cada uno de un ÚNICO enlace rol=inferido directo a
    INTERNET-MPLS -- ambos son puntos de articulación en P1.
(c) Yanuncay depende también de un único enlace (AGRPRI-1A-D10 --
    ROUTER-CAMPUS-YANUNCAY) para sus 11 equipos de acceso.
(d) BAL-AUL2-D1 protege 10 equipos de acceso tras un único switch de
    agregación (aunque SU propio uplink al core ya es redundante) -- el
    equipo de mayor tráfico detrás de él (AUL2-0A-A121, 81.71 Mbps) queda
    sin ruta alterna si BAL-AUL2-D1 falla.
""")

struct Intervencion
    tipo::Symbol  # :nueva o :duplicar
    origen::String
    destino::String
    capacidad_mbps::Float64
    justificacion::String
end

intervenciones = [
    Intervencion(:duplicar, "CPAR-C10", "ROUTER-CAMPUS-HUAYNA-CAPAC", 10_000.0,
        "Único enlace de salida de todo Campus Paraíso (35 nodos de acceso); hoy 1 Gbps, el " *
        "cuello de botella de flujo más severo de la institución (P6). Se sube al estándar " *
        "de 10 Gbps que ya usan el resto de troncales WAN documentados."),
    Intervencion(:nueva, "HOS-0A-D05", "FORTIGATE-1800F-CENTRAL", 1_000.0,
        "Hospitalidad depende hoy de un único enlace rol=inferido directo a INTERNET-MPLS " *
        "(punto de articulación en P1). Ruta alterna vía el firewall perimetral de Central."),
    Intervencion(:nueva, "CCJ-CJURIDICO-D4", "DATCC-2A-C2", 1_000.0,
        "Mismo patrón que Hospitalidad: el Consultorio Jurídico solo llega a la red por un " *
        "enlace inferido a INTERNET-MPLS. Ruta alterna vía un core de Campus Central."),
    Intervencion(:nueva, "ROUTER-CAMPUS-YANUNCAY", "DATCC-2A-C3", 10_000.0,
        "Yanuncay (11 nodos de acceso) depende de un único enlace (AGRPRI-1A-D10 -- " *
        "ROUTER-CAMPUS-YANUNCAY) para salir de campus; se añade una segunda ruta troncal."),
    Intervencion(:nueva, "AUL2-0A-A121", "BAL-CENTEC-D2", 1_000.0,
        "Mitigación puntual: el equipo de mayor tráfico tras BAL-AUL2-D1 (81.71 Mbps, el " *
        "switch que más nodos de acceso protege en el top-10 de P10) obtiene una ruta " *
        "alterna hacia otro switch de agregación del mismo campus. No resuelve el patrón " *
        "general de acceso de único homing (132 nodos en toda la red, fuera de presupuesto)."),
]
for (i, iv) in enumerate(intervenciones)
    @printf("%d. [%s] %s -- %s (%.0f Mbps)\n   %s\n", i, iv.tipo, iv.origen, iv.destino, iv.capacidad_mbps, iv.justificacion)
end
guardar_csv(DataFrame(intervenciones), joinpath(DIR_RESULTADOS, "p11_intervenciones_propuestas.csv"))

"""
    aplicar_intervenciones(g, cap_dict, intervenciones) -> (g2, cap_dict2)

`:nueva` agrega una arista real (cambia topología: afecta puentes,
articulación, distancia, eficiencia, percolación y flujo). `:duplicar`
NO cambia la topología del grafo simple -- una segunda fibra entre el
MISMO PAR de equipos no crea una ruta alternativa a nivel de nodos, así
que deliberadamente NO se cuenta como una arista nueva -- solo reemplaza
la capacidad de la arista existente (afecta flujo, no puentes/articulación
de ESE enlace). Esta distinción se discute explícitamente en el informe.
Parte de `cap_dict` (no de un vector posicional) por la misma razón que
`flujo_maximo_campus`: evitar cualquier dependencia del orden de
iteración de `edges(g)`.
"""
function aplicar_intervenciones(g::SimpleGraph, cap_dict::Dict{Tuple{Int,Int},Float64}, intervenciones::Vector{Intervencion})
    g2 = SimpleGraph(nv(g))
    for e in edges(g)
        add_edge!(g2, src(e), dst(e))
    end
    cap_dict2 = copy(cap_dict)
    for iv in intervenciones
        u, v = red.idx[iv.origen], red.idx[iv.destino]
        clave = u < v ? (u, v) : (v, u)
        iv.tipo == :nueva && add_edge!(g2, u, v)
        cap_dict2[clave] = iv.capacidad_mbps
    end
    return g2, cap_dict2
end

g_prop, cap_dict_prop = aplicar_intervenciones(g0, cap_dict0, intervenciones)
println("\nEvaluando la red MODIFICADA...")
metricas_prop = evaluar_red(g_prop, cap_dict_prop, idx_sink)
println(metricas_prop)

# =================================================================
# P11.2 -- Tabla antes/después/variación
# =================================================================
println("\n--- P11.2 Antes / después / variación ---")
claves = collect(keys(metricas_base))
df_comparacion = DataFrame(
    métrica=String.(claves),
    antes=[getfield(metricas_base, k) for k in claves],
    despues=[getfield(metricas_prop, k) for k in claves],
)
df_comparacion.variacion = df_comparacion.despues .- df_comparacion.antes
df_comparacion.variacion_pct = round.(100 .* df_comparacion.variacion ./ abs.(df_comparacion.antes); digits=1)
guardar_csv(df_comparacion, joinpath(DIR_RESULTADOS, "p11_antes_despues.csv"))
println(df_comparacion)

f_vals = (0:n) ./ n
function curva_percolacion_grado(g::SimpleGraph)
    nn = nv(g)
    orden = sortperm(degree(g); rev=true)
    vivo = trues(nn)
    S = zeros(Float64, nn + 1); S[1] = 1.0
    for k in 1:nn
        vivo[orden[k]] = false
        idx = findall(vivo)
        S[k+1] = isempty(idx) ? 0.0 : maximum(length.(connected_components(induced_subgraph(g, idx)[1]))) / nn
    end
    return S
end
S_base = curva_percolacion_grado(g0)
S_prop = curva_percolacion_grado(g_prop)
plt_comparacion = plot(f_vals, S_base; label="Red base", lw=2,
                        xlabel="fracción de nodos eliminados (f, ataque por grado)",
                        ylabel="S(f)", title="Percolación bajo ataque dirigido (P11): antes vs. después",
                        titlefontsize=11, top_margin=3Plots.mm)
plot!(plt_comparacion, f_vals, S_prop; label="Red con las 5 intervenciones", lw=2)
savefig(plt_comparacion, joinpath(DIR_FIGURAS, "p11_percolacion_antes_despues.png"))

# =================================================================
# P11.3 -- Comparación con alternativas ingenuas
# =================================================================
println("\n--- P11.3 Comparación con alternativas ---")
"5 aristas entre los pares de mayor grado combinado que aún no están conectados (criterio ingenuo)."
function alternativa_mayor_grado(g::SimpleGraph, k::Int=5)
    nn = nv(g)
    grados_ = degree(g)
    candidatos = sortperm(grados_; rev=true)[1:min(30, nn)]
    pares = Tuple{Int,Int}[]
    for i in candidatos, j in candidatos
        i < j && !has_edge(g, i, j) && push!(pares, (i, j))
    end
    sort!(pares; by=p -> -(grados_[p[1]] + grados_[p[2]]))
    return pares[1:k]
end
"5 aristas al azar entre pares no conectados (línea base nula)."
function alternativa_aleatoria(g::SimpleGraph, k::Int=5; semilla::Int=0)
    Random.seed!(semilla)
    nn = nv(g)
    elegidas = Tuple{Int,Int}[]
    while length(elegidas) < k
        u, v = rand(1:nn), rand(1:nn)
        u == v && continue
        e = u < v ? (u, v) : (v, u)
        (has_edge(g, u, v) || e in elegidas) && continue
        push!(elegidas, e)
    end
    return elegidas
end

function evaluar_alternativa(g::SimpleGraph, cap_dict::Dict{Tuple{Int,Int},Float64}, pares::Vector{Tuple{Int,Int}})
    g2 = SimpleGraph(nv(g))
    for e in edges(g); add_edge!(g2, src(e), dst(e)); end
    for (u, v) in pares; add_edge!(g2, u, v); end
    cap_dict2 = copy(cap_dict)
    for (u, v) in pares; cap_dict2[(u, v)] = 1_000.0; end
    return evaluar_red(g2, cap_dict2, idx_sink)
end

pares_grado = alternativa_mayor_grado(g0)
pares_aleatorios = alternativa_aleatoria(g0)
println("Alternativa A (ingenua, mayor grado): ", [(red.ids[u], red.ids[v]) for (u, v) in pares_grado])
println("Alternativa B (aleatoria): ", [(red.ids[u], red.ids[v]) for (u, v) in pares_aleatorios])

metricas_alt_grado = evaluar_alternativa(g0, cap_dict0, pares_grado)
metricas_alt_aleatoria = evaluar_alternativa(g0, cap_dict0, pares_aleatorios)

df_alternativas = DataFrame(
    propuesta=["Red base (sin cambios)", "P11 (propuesta justificada)", "Alternativa: mayor grado", "Alternativa: aleatoria"],
    puentes=[metricas_base.puentes, metricas_prop.puentes, metricas_alt_grado.puentes, metricas_alt_aleatoria.puentes],
    articulacion=[metricas_base.articulacion, metricas_prop.articulacion, metricas_alt_grado.articulacion, metricas_alt_aleatoria.articulacion],
    distancia_media=round.([metricas_base.distancia_media, metricas_prop.distancia_media, metricas_alt_grado.distancia_media, metricas_alt_aleatoria.distancia_media]; digits=3),
    eficiencia_global=round.([metricas_base.eficiencia_global, metricas_prop.eficiencia_global, metricas_alt_grado.eficiencia_global, metricas_alt_aleatoria.eficiencia_global]; digits=4),
    robustez_auc=round.([metricas_base.robustez_auc_percolacion, metricas_prop.robustez_auc_percolacion, metricas_alt_grado.robustez_auc_percolacion, metricas_alt_aleatoria.robustez_auc_percolacion]; digits=4),
    flujo_paraiso_mbps=[metricas_base.flujo_paraiso_mbps, metricas_prop.flujo_paraiso_mbps, metricas_alt_grado.flujo_paraiso_mbps, metricas_alt_aleatoria.flujo_paraiso_mbps],
    flujo_hospitalidad_mbps=[metricas_base.flujo_hospitalidad_mbps, metricas_prop.flujo_hospitalidad_mbps, metricas_alt_grado.flujo_hospitalidad_mbps, metricas_alt_aleatoria.flujo_hospitalidad_mbps],
)
guardar_csv(df_alternativas, joinpath(DIR_RESULTADOS, "p11_comparacion_alternativas.csv"))
println(df_alternativas)
println("""
La propuesta justificada ataca los CUELLOS DE BOTELLA REALES identificados en
P1/P6/P10 (Paraíso, Hospitalidad, Consultorio Jurídico, Yanuncay, un tramo de
acceso de Balzay). "Unir los dos nodos de mayor grado" no repara ningún
punto de articulación real -- casi siempre conecta dos switches de CORE que
YA tienen múltiples caminos entre sí, así que apenas mueve puentes/
articulación/distancia/eficiencia; el flujo de Paraíso e Hospitalidad, en
particular, NO mejora con ninguna alternativa ingenua porque ninguna de
las dos toca esos campus.
""")

# =================================================================
# P11.4-P11.5 -- Factibilidad y limitaciones
# =================================================================
println("""
--- P11.4 Factibilidad ---
Las 5 intervenciones son, en términos de obra civil, de complejidad muy
distinta: (1) y (4) son upgrades de un enlace WAN YA EXISTENTE entre
equipos que ya están físicamente enlazados (solo cambio de óptica/SFP y
posible contrato de mayor ancho de banda con el proveedor MPLS) -- bajo
costo, alta factibilidad. (2) y (3) requieren tender fibra nueva entre
Hospitalidad/Consultorio Jurídico y Campus Central -- si la distancia
física excede lo razonable para fibra propia, en la práctica esto se
resolvería con un segundo circuito del proveedor, no con un enlace físico
directo (limitación del modelo: no conocemos la distancia real ni el
trazado disponible). (5) es la de menor costo (cableado interno del mismo
campus) pero también la de menor alcance (un solo equipo).

--- P11.5 Limitaciones del estudio ---
- Las capacidades no documentadas (181 de 209 enlaces) son ESTIMADAS por
  reglas heurísticas (P6.1); la magnitud del beneficio de cada intervención
  depende de esas estimaciones, no de mediciones reales.
- El modelo no incorpora energía, espacio físico en rack, ni licenciamiento
  -- ver P7.4. Una intervención "barata" en el grafo puede no serlo en la
  práctica.
- La red se trata como estática: no hay series de tiempo de tráfico, así
  que no se puede distinguir un cuello de botella permanente de uno que
  solo aparece en horas pico.
- Los resultados de percolación / cascada / SIR son promedios de muchas
  realizaciones; una falla real concreta puede comportarse distinto a la
  media reportada.
- No puede concluirse de este análisis cuál sería el ROI económico de cada
  intervención, solo su impacto estructural sobre las métricas de red.
""")

println("="^70)
println("FASE 5 completa. Resultados en ", DIR_RESULTADOS, " y figuras en ", DIR_FIGURAS)
println("="^70)
