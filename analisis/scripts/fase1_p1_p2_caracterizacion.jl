#=============================================================
  fase1_p1_p2_caracterizacion.jl -- Fase 1: Modelado y Caracterización
  Proyecto integrador de Redes Complejas (1217) -- Universidad de Cuenca

  Resuelve P1 (medidas fundamentales) y P2 (modelos nulos y
  visualización) del enunciado (proyecto_red_ucuenca.pdf, cap. 3).

  Uso:
      julia --project=analisis analisis/scripts/fase1_p1_p2_caracterizacion.jl

  Reutiliza cargar_red()/verificar() de codigo_base/cargar_red.jl (no se
  reescribe la carga) y las utilidades de analisis/src/ModeloRed.jl.
=============================================================#

using Graphs
using DataFrames
using CSV
using Statistics
using StatsBase
using Printf
using Random
using Plots
using NetworkLayout

const DIR_ESTE_SCRIPT = @__DIR__
const DIR_ANALISIS = dirname(DIR_ESTE_SCRIPT)
const DIR_RAIZ = dirname(DIR_ANALISIS)

include(joinpath(DIR_RAIZ, "codigo_base", "cargar_red.jl"))
include(joinpath(DIR_ANALISIS, "src", "ModeloRed.jl"))

const DIR_RESULTADOS = joinpath(DIR_ANALISIS, "resultados")
const DIR_FIGURAS = joinpath(DIR_ANALISIS, "figuras")
mkpath(DIR_RESULTADOS)
mkpath(DIR_FIGURAS)

println("="^70)
println("FASE 1 -- Modelado y Caracterización (P1-P2)")
println("="^70)

red = cargar_red()
verificar(red) || error("La carga no reproduce el Anexo A: revise el conjunto de datos.")
g = red.g
n, m = nv(g), ne(g)

# =================================================================
# P1.1 -- Medidas fundamentales
# =================================================================
println("\n--- P1.1 Medidas fundamentales ---")
componentes = connected_components(g)
densidad = 2m / (n * (n - 1))
@printf("Nodos: %d   Aristas: %d   Densidad: %.4f\n", n, m, densidad)
@printf("Componentes conexas: %d (tamaño de la mayor: %d)\n",
        length(componentes), maximum(length.(componentes)))

# =================================================================
# P1.2 -- Distribución de grado
# =================================================================
println("\n--- P1.2 Distribución de grado ---")
grados = degree(g)
@printf("Grado medio: %.2f   máximo: %d   mínimo: %d\n",
        mean(grados), maximum(grados), minimum(grados))

conteo_grados = countmap(grados)
df_grados = sort(DataFrame(grado=collect(keys(conteo_grados)),
                            frecuencia=collect(values(conteo_grados))), :grado)
guardar_csv(df_grados, joinpath(DIR_RESULTADOS, "p1_distribucion_grado.csv"))

plt_hist = histogram(grados; bins=maximum(grados), legend=false,
                      xlabel="grado k", ylabel="número de nodos",
                      title="Distribución de grado -- red UCuenca")
plt_loglog = scatter(df_grados.grado, df_grados.frecuencia;
                      xscale=:log10, yscale=:log10, legend=false,
                      xlabel="grado k (log)", ylabel="frecuencia (log)",
                      title="Distribución de grado, escala log-log",
                      markersize=5)
savefig(plot(plt_hist, plt_loglog; layout=(1, 2), size=(1100, 450)),
        joinpath(DIR_FIGURAS, "p1_distribucion_grado.png"))

# Ajuste de ley de potencias discreta por máxima verosimilitud (Clauset,
# Shalizi & Newman, 2009, ec. 3), solo como referencia exploratoria: con
# n = 177 nodos la muestra es demasiado pequeña para afirmar "libre de
# escala" sin una prueba formal (bondad de ajuste Kolmogorov-Smirnov +
# comparación por razón de verosimilitud contra alternativas como la
# log-normal, ver Clauset et al. 2009, sec. 3-4). Aquí solo se reporta el
# exponente estimado para la cola k >= kmin, como insumo de la discusión.
function ajuste_powerlaw_discreto(datos::Vector{<:Integer}, kmin::Integer)
    cola = filter(>=(kmin), datos)
    isempty(cola) && return NaN, 0
    α = 1 + length(cola) / sum(log.(cola ./ (kmin - 0.5)))
    return α, length(cola)
end
kmin_cola = round(Int, quantile(grados, 0.75))
α_hat, n_cola = ajuste_powerlaw_discreto(grados, kmin_cola)
@printf("Ajuste exploratorio de cola de potencias: kmin=%d, α≈%.2f (n_cola=%d de %d)\n",
        kmin_cola, α_hat, n_cola, n)
println("Con n=177 esta estimación NO alcanza para afirmar 'libre de escala'; ",
        "haría falta la prueba de bondad de ajuste (KS) y la comparación por ",
        "razón de verosimilitud de Clauset-Shalizi-Newman (2009) frente a ",
        "alternativas (log-normal, exponencial).")

# =================================================================
# P1.3 -- Centralidades
# =================================================================
println("\n--- P1.3 Centralidades ---")
cent_grado = Float64.(grados)
cent_cercania = closeness_centrality(g)
cent_intermediacion = betweenness_centrality(g)
cent_vecprop = eigenvector_centrality(g)

tabla_cent = DataFrame(
    puesto = 1:10,
    id_grado = red.ids[sortperm(cent_grado; rev=true)[1:10]],
    grado = round.(sort(cent_grado; rev=true)[1:10]; digits=2),
    id_cercania = red.ids[sortperm(cent_cercania; rev=true)[1:10]],
    cercania = round.(sort(cent_cercania; rev=true)[1:10]; digits=4),
    id_intermediacion = red.ids[sortperm(cent_intermediacion; rev=true)[1:10]],
    intermediacion = round.(sort(cent_intermediacion; rev=true)[1:10]; digits=4),
    id_vector_propio = red.ids[sortperm(cent_vecprop; rev=true)[1:10]],
    vector_propio = round.(sort(cent_vecprop; rev=true)[1:10]; digits=4),
)
guardar_csv(tabla_cent, joinpath(DIR_RESULTADOS, "p1_centralidades_top10.csv"))
println("Top-3 por grado: ", join(tabla_cent.id_grado[1:3], ", "))
println("Top-3 por intermediación: ", join(tabla_cent.id_intermediacion[1:3], ", "))
coinciden = length(intersect(Set(tabla_cent.id_grado), Set(tabla_cent.id_intermediacion)))
println("Nodos que coinciden en el top-10 de grado E intermediación: ", coinciden, "/10")

# =================================================================
# P1.4 -- Clustering, diámetro, distancia media, asortatividad
# =================================================================
println("\n--- P1.4 Clustering, diámetro, distancia media, asortatividad ---")
clustering_medio = global_clustering_coefficient(g)
diam = diameter(g)

function distancia_media_y_diametro(g::SimpleGraph)
    n = nv(g)
    suma = 0
    npares = 0
    dmax = 0
    for u in 1:n
        d = gdistances(g, u)
        for v in (u+1):n
            suma += d[v]
            npares += 1
            dmax = max(dmax, d[v])
        end
    end
    return suma / npares, dmax
end
dist_media, diam_verificado = distancia_media_y_diametro(g)
@assert diam_verificado == diam "Diámetro inconsistente entre Graphs.jl y el cálculo directo"

function asortatividad_grado(g::SimpleGraph)
    # Newman (2010), "Networks: An Introduction", ec. 8.28: coeficiente de
    # correlación de Pearson entre los grados de los dos extremos de cada
    # arista (forma simétrica, cada arista contribuye una vez).
    grados_locales = degree(g)
    M = ne(g)
    S1 = S2 = S3 = 0.0
    for e in edges(g)
        j, k = grados_locales[src(e)], grados_locales[dst(e)]
        S1 += j * k
        S2 += 0.5 * (j + k)
        S3 += 0.5 * (j^2 + k^2)
    end
    S1 /= M; S2 /= M; S3 /= M
    return (S1 - S2^2) / (S3 - S2^2)
end
asort = asortatividad_grado(g)

@printf("Clustering medio: %.4f\n", clustering_medio)
@printf("Diámetro: %d   Distancia media: %.3f\n", diam, dist_media)
@printf("Asortatividad por grado: %.4f\n", asort)

# =================================================================
# P1.5 -- Puntos de articulación y puentes
# =================================================================
println("\n--- P1.5 Puntos de articulación y puentes ---")
articulacion, puentes = puentes_y_articulacion(g)
@printf("Puntos de articulación: %d   Puentes: %d\n", length(articulacion), length(puentes))

df_articulacion = DataFrame(
    id = red.ids[articulacion],
    campus = red.campus[articulacion],
    capa = red.capa[articulacion],
    grado = grados[articulacion],
)
guardar_csv(df_articulacion, joinpath(DIR_RESULTADOS, "p1_puntos_articulacion.csv"))

df_puentes = DataFrame(
    origen = [red.ids[u] for (u, v) in puentes],
    destino = [red.ids[v] for (u, v) in puentes],
    campus_origen = [red.campus[u] for (u, v) in puentes],
    campus_destino = [red.campus[v] for (u, v) in puentes],
    capa_origen = [red.capa[u] for (u, v) in puentes],
    capa_destino = [red.capa[v] for (u, v) in puentes],
)
guardar_csv(df_puentes, joinpath(DIR_RESULTADOS, "p1_puentes.csv"))

println("Puntos de articulación por campus:")
for (c, k) in sort(collect(countmap(red.campus[articulacion])); by=x -> -x[2])
    @printf("  %-24s %3d\n", c, k)
end
println("Puentes por campus (contando el campus de cada extremo):")
campus_puentes = vcat(df_puentes.campus_origen, df_puentes.campus_destino)
for (c, k) in sort(collect(countmap(campus_puentes)); by=x -> -x[2])
    @printf("  %-24s %3d\n", c, k)
end

# =================================================================
# P1.6 -- Contraste con el informe técnico (redundancia core-agregación)
# =================================================================
println("\n--- P1.6 Contraste con el informe técnico ---")
# El informe (cap. 2.1) afirma redundancia física completa core-agregación
# en Balzay y Paraíso, y enlaces simples en Campus Central. Se verifica
# contando, por campus, cuántos nodos de agregación tienen >=2 aristas
# hacia nodos de capa `core` (redundancia real) frente a los que tienen
# exactamente 1 (enlace simple).
idx_agregacion = findall(==("agregacion"), red.capa)
filas = NamedTuple[]
for i in idx_agregacion
    vecinos_core = count(v -> red.capa[v] == "core", neighbors(g, i))
    push!(filas, (id=red.ids[i], campus=red.campus[i], enlaces_a_core=vecinos_core))
end
df_redundancia = DataFrame(filas)
guardar_csv(df_redundancia, joinpath(DIR_RESULTADOS, "p1_redundancia_core_agregacion.csv"))

println("Nodos de agregación con >=2 enlaces a core (redundancia real), por campus:")
resumen_redundancia = combine(groupby(df_redundancia, :campus),
    :enlaces_a_core => (x -> count(>=(2), x)) => :con_redundancia,
    :enlaces_a_core => (x -> count(==(1), x)) => :enlace_simple,
    nrow => :total)
println(resumen_redundancia)
guardar_csv(resumen_redundancia, joinpath(DIR_RESULTADOS, "p1_redundancia_resumen.csv"))
println("Balzay y Paraíso deberían mostrar mayoría con >=2 enlaces si el informe es exacto; ",
        "el diccionario de datos (red_ucuenca_README.md) ya adelanta que Paraíso en realidad ",
        "usa agregación de puertos hacia un único core (CPAR-C10), no redundancia de núcleo: ",
        "la tabla de arriba debe confirmarlo con las aristas reales.")

# Tabla resumen global de P1, para el informe
tabla_global = DataFrame(
    métrica = ["Nodos", "Aristas", "Densidad", "Componentes conexas",
               "Grado medio", "Grado máximo", "Grado mínimo",
               "Clustering medio", "Diámetro", "Distancia media",
               "Asortatividad por grado", "Puntos de articulación", "Puentes"],
    valor = [n, m, round(densidad; digits=4), length(componentes),
             round(mean(grados); digits=2), maximum(grados), minimum(grados),
             round(clustering_medio; digits=4), diam, round(dist_media; digits=3),
             round(asort; digits=4), length(articulacion), length(puentes)],
)
guardar_csv(tabla_global, joinpath(DIR_RESULTADOS, "p1_medidas_globales.csv"))

# =================================================================
# P2 -- Modelos nulos y visualización
# =================================================================
println("\n" * "="^70)
println("P2 -- Modelos nulos y visualización")
println("="^70)

const N_REALIZACIONES = 100

function métricas_red(h::SimpleGraph)
    comps = connected_components(h)
    hg = induced_subgraph(h, comps[argmax(length.(comps))])[1]  # componente gigante, para que diámetro/distancia estén definidos
    dm, dg = distancia_media_y_diametro(hg)
    return (clustering=global_clustering_coefficient(h), dist_media=dm,
            diametro=dg, asortatividad=asortatividad_grado(h))
end

function resumir(nombre::String, muestras::Vector)
    cl = [s.clustering for s in muestras]
    dm = [s.dist_media for s in muestras]
    dg = [s.diametro for s in muestras]
    at = [s.asortatividad for s in muestras]
    (modelo=nombre,
     clustering_media=mean(cl), clustering_std=std(cl),
     dist_media_media=mean(dm), dist_media_std=std(dm),
     diametro_media=mean(dg), diametro_std=std(dg),
     asortatividad_media=mean(at), asortatividad_std=std(at))
end

println("\nGenerando $N_REALIZACIONES realizaciones de Erdős–Rényi G(n,m)...")
Random.seed!(1217)
muestras_er = [métricas_red(erdos_renyi(n, m)) for _ in 1:N_REALIZACIONES]

println("Generando $N_REALIZACIONES realizaciones del modelo de configuración...")
muestras_conf = Vector{NamedTuple}(undef, N_REALIZACIONES)
aristas_conf_logradas = zeros(Int, N_REALIZACIONES)
for i in 1:N_REALIZACIONES
    h = modelo_configuracion(collect(grados))
    aristas_conf_logradas[i] = ne(h)
    muestras_conf[i] = métricas_red(h)
end
@printf("Aristas logradas por el modelo de configuración: %.1f de %d en promedio (%.1f%%)\n",
        mean(aristas_conf_logradas), m, 100 * mean(aristas_conf_logradas) / m)

resumen_real = (modelo="UCuenca (real)",
                clustering_media=clustering_medio, clustering_std=0.0,
                dist_media_media=dist_media, dist_media_std=0.0,
                diametro_media=diam, diametro_std=0.0,
                asortatividad_media=asort, asortatividad_std=0.0)

tabla_nulos = DataFrame([resumen_real, resumir("Erdős–Rényi G(n,m)", muestras_er),
                          resumir("Configuración (grados exactos)", muestras_conf)])
guardar_csv(tabla_nulos, joinpath(DIR_RESULTADOS, "p2_modelos_nulos.csv"))
println(tabla_nulos)

println("\nInterpretación P2.1: si clustering_media y asortatividad de UCuenca caen lejos")
println("(varias desviaciones estándar) de ER, esa propiedad NO se explica solo por n y m;")
println("si además el modelo de configuración -que sí fija la secuencia de grados- se acerca")
println("más al valor real que ER, la propiedad es explicada en gran parte por la heterogeneidad")
println("de grados (p.ej. los switches de core con grado muy alto) y no por más estructura.")

# --- Barabási-Albert, n y m comparables ---
k_ba = max(1, round(Int, m / n))  # cada nodo nuevo se conecta a k_ba existentes -> m ≈ k_ba*(n-k_ba)
g_ba = barabasi_albert(n, k_ba)
met_ba = métricas_red(g_ba)
@printf("\nBarabási-Albert (n=%d, k=%d, %d aristas): clustering=%.4f, dist_media=%.3f, diametro=%d, asort=%.4f\n",
        n, k_ba, ne(g_ba), met_ba.clustering, met_ba.dist_media, met_ba.diametro, met_ba.asortatividad)
println("Una red de infraestructura física como UCuenca NO crece por conexión preferencial libre:")
println("cada switch se conecta a un número acotado de puertos y a la jerarquía core/agregación/acceso")
println("que dicta el diseño, no a 'quien ya tiene más enlaces'. Barabási-Albert sirve de referencia")
println("de cuán distinto sería el grado más extremo bajo ese mecanismo de crecimiento.")

guardar_csv(DataFrame([merge((modelo="Barabási-Albert (1 realización)",), met_ba)]),
            joinpath(DIR_RESULTADOS, "p2_barabasi_albert.csv"))

# --- Visualizaciones ---
println("\nGenerando visualizaciones...")
campus_unicos = sort(unique(red.campus))
paleta = palette(:tab10, length(campus_unicos))
color_por_campus = Dict(c => paleta[i] for (i, c) in enumerate(campus_unicos))
colores_nodos = [color_por_campus[c] for c in red.campus]

# Layout de resortes (Fruchterman-Reingold vía NetworkLayout.jl): apropiado
# para un grafo disperso jerárquico como este, separa visualmente los
# distintos campus sin coordenadas geográficas. Con 177 nodos, etiquetar
# cada uno sería ilegible: se usa color/tamaño y una leyenda aparte en vez
# de texto por nodo. Se dibuja a mano con NetworkLayout + Plots.scatter (en
# vez de GraphRecipes.graphplot, cuyo recipe en la versión instalada falla
# al recibir un tamaño de nodo distinto por nodo -- ver discusión en el
# informe/README; el resto del código de referencia del módulo ya dibuja
# redes a mano con Plots de la misma manera, ver
# codigo_referencia/ford-fulkerson/ford_fulkerson.jl).
Random.seed!(1217)
pos = NetworkLayout.spring(g; iterations=200)
xs = [p[1] for p in pos]
ys = [p[2] for p in pos]

function lienzo_red(; titulo="")
    plt = plot(; legend=false, axis=false, grid=false, ticks=false,
               title=titulo, size=(1000, 700))
    for e in edges(g)
        u, v = src(e), dst(e)
        plot!(plt, [xs[u], xs[v]], [ys[u], ys[v]]; color=:gray70, alpha=0.5, lw=0.6)
    end
    return plt
end

plt_campus = lienzo_red(; titulo="Red UCuenca coloreada por campus")
scatter!(plt_campus, xs, ys; markersize=4, color=colores_nodos, markerstrokewidth=0.3)
leyenda_campus = plot(; legend=:left, framestyle=:none, grid=false, axis=false)
for c in campus_unicos  # una serie por campus: Plots.jl solo agrega una entrada de leyenda por serie, no por color de punto
    scatter!(leyenda_campus, [NaN], [NaN]; color=color_por_campus[c], label=c, markersize=6)
end
savefig(plot(plt_campus, leyenda_campus; layout=@layout([a{0.82w} b]), size=(1200, 700)),
        joinpath(DIR_FIGURAS, "p2_red_por_campus.png"))

tam_nodos = 2 .+ 18 .* (cent_intermediacion ./ max(maximum(cent_intermediacion), eps()))
plt_intermediacion = lienzo_red(; titulo="Red UCuenca -- tamaño de nodo ∝ intermediación")
scatter!(plt_intermediacion, xs, ys; markersize=tam_nodos, color=:steelblue, markerstrokewidth=0.3)
savefig(plt_intermediacion, joinpath(DIR_FIGURAS, "p2_red_por_intermediacion.png"))
println("Figuras guardadas en ", DIR_FIGURAS)

println("\n" * "="^70)
println("FASE 1 completa. Resultados en ", DIR_RESULTADOS, " y figuras en ", DIR_FIGURAS)
println("="^70)
