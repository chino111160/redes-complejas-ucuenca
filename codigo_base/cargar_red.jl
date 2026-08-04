#=============================================================
  cargar_red.jl -- Código base del proyecto integrador
  Módulo de Redes Complejas (1217) -- Universidad de Cuenca

  Carga la red de datos de la Universidad de Cuenca desde los CSV,
  verifica el pipeline contra los valores del Anexo A del enunciado e
  imprime un resumen estructural.

  Uso:
      julia --project=. cargar_red.jl
      julia --project=. -e 'include("cargar_red.jl"); red = cargar_red()'
=============================================================#

using CSV
using DataFrames
using Graphs
using Printf
using Statistics

const DIR_SCRIPT = @__DIR__
const TESTIGO = "red_ucuenca_nodes.csv"

"""
    buscar_dir_datos() -> String

Busca hacia arriba el directorio que contiene el conjunto de datos. Así el
código funciona igual dentro del repositorio (`proyecto/codigo_base/` con los
datos en la raíz) que en el paquete entregado a los estudiantes
(`codigo_base/` con los datos un nivel arriba), sin editar nada.
"""
function buscar_dir_datos()
    dir = DIR_SCRIPT
    while true
        isfile(joinpath(dir, TESTIGO)) && return dir
        padre = dirname(dir)
        padre == dir && break
        dir = padre
    end
    dirname(DIR_SCRIPT)  # sin datos a la vista: falla luego con un error claro
end

"""
    ruta_datos(nombre) -> String

Ruta absoluta a un archivo del conjunto de datos. El directorio puede fijarse
con la variable de entorno `RED_UCUENCA_DIR`; si no, se busca hacia arriba
desde la ubicación de este script.
"""
function ruta_datos(nombre::AbstractString)
    base = get(ENV, "RED_UCUENCA_DIR", buscar_dir_datos())
    joinpath(base, nombre)
end

"""
Estructura que agrupa el grafo y sus atributos.

Campos:
  - `g`       : `SimpleGraph` no dirigido, 177 nodos y 209 aristas.
  - `ids`     : identificadores, indexados por vértice.
  - `idx`     : diccionario identificador -> índice de vértice.
  - `campus`  : ubicación física de cada vértice.
  - `capa`    : capa jerárquica de cada vértice (core / agregacion / acceso /
                wan / interconexion).
  - `nodos`   : `DataFrame` de nodos.
  - `aristas` : `DataFrame` de aristas, con las columnas `u` y `v` añadidas
                (índices de vértice).
"""
struct RedUCuenca
    g::SimpleGraph{Int}
    ids::Vector{String}
    idx::Dict{String,Int}
    campus::Vector{String}
    capa::Vector{String}
    nodos::DataFrame
    aristas::DataFrame
end

"""
    cargar_red() -> RedUCuenca

Construye el grafo a partir de `red_ucuenca_nodes.csv` y
`red_ucuenca_edges.csv`.

El grafo es **simple y no dirigido**. Los identificadores son idénticos a los
del GraphML, de modo que ambas fuentes producen el mismo grafo.

DECISIÓN DE MODELADO: el grafo se entrega SIN pesos. Elegir el modelo de peso
(saltos, latencia, carga) es parte del Problema P5; estimar la capacidad de
los enlaces cuya capacidad nominal no está documentada es parte del P6.
"""
function cargar_red()
    nodos = CSV.read(ruta_datos("red_ucuenca_nodes.csv"), DataFrame;
                     types=Dict(:id => String, :label => String,
                                :campus => String, :capa => String,
                                :diagrams => String))
    # Los tipos se fijan explícitamente: si se dejaran inferir, una columna
    # cuyos valores sean todos enteros (capacidad_mbps) se leería como Int64 y
    # la otra como Float64, con semánticas de división distintas.
    # Las celdas vacías quedan como `missing`, la convención de Julia.
    aristas = CSV.read(ruta_datos("red_ucuenca_edges.csv"), DataFrame;
                       types=Dict(:source => String, :target => String,
                                  :label => String,
                                  :trafico_mbps => Float64,
                                  :capacidad_mbps => Float64,
                                  :rol => String,
                                  :diagrams => String))

    ids = String.(nodos.id)
    @assert length(unique(ids)) == length(ids) "Identificadores de nodo duplicados"
    idx = Dict(id => i for (i, id) in enumerate(ids))

    g = SimpleGraph(length(ids))
    u = Vector{Int}(undef, nrow(aristas))
    v = Vector{Int}(undef, nrow(aristas))
    for (k, r) in enumerate(eachrow(aristas))
        haskey(idx, r.source) || error("Arista $k: nodo origen desconocido: $(r.source)")
        haskey(idx, r.target) || error("Arista $k: nodo destino desconocido: $(r.target)")
        u[k], v[k] = idx[r.source], idx[r.target]
        add_edge!(g, u[k], v[k])
    end
    aristas.u = u
    aristas.v = v

    RedUCuenca(g, ids, idx, String.(nodos.campus), String.(nodos.capa),
               nodos, aristas)
end

"""
    trafico(red), capacidad(red)

Vectores alineados con `red.aristas`.

  - `trafico`   : tráfico MEDIDO en el instante de la captura, en Mbps
                  (definido en 170 de las 209 aristas).
  - `capacidad` : capacidad NOMINAL declarada en la etiqueta del enlace
                  (definida en 28 aristas, todas del diagrama MPLS).

Son magnitudes distintas: no las mezcle. En el Problema P6 deberá estimar la
capacidad de los enlaces donde `capacidad` no está definida y justificar el
criterio empleado.
"""
trafico(red::RedUCuenca)   = [ismissing(x) ? missing : Float64(x) for x in red.aristas.trafico_mbps]
capacidad(red::RedUCuenca) = [ismissing(x) ? missing : Float64(x) for x in red.aristas.capacidad_mbps]

"Conteo por valor para un vector de cadenas (evita una dependencia extra)."
function countmap_str(v::AbstractVector{<:AbstractString})
    d = Dict{String,Int}()
    for x in v
        d[x] = get(d, x, 0) + 1
    end
    d
end

"""
    verificar(red) -> Bool

Comprueba la carga contra los valores de referencia del Anexo A del enunciado.
Devuelve `true` si todas las comprobaciones pasan.
"""
function verificar(red::RedUCuenca)
    roles = countmap_str(red.aristas.rol)
    capas = countmap_str(red.capa)
    comprobaciones = [
        ("Nodos",                       nv(red.g),                              177),
        ("Aristas",                     ne(red.g),                              209),
        ("Componentes conexas",         length(connected_components(red.g)),      1),
        ("Nodos en capa acceso",        get(capas, "acceso", 0),                132),
        ("Nodos en capa agregacion",    get(capas, "agregacion", 0),             27),
        ("Nodos en capa core",          get(capas, "core", 0),                    5),
        ("Aristas rol=principal",       get(roles, "principal", 0),             157),
        ("Aristas rol=wan",             get(roles, "wan", 0),                    32),
        ("Aristas rol=respaldo",        get(roles, "respaldo", 0),               12),
        ("Aristas rol=secundario",      get(roles, "secundario", 0),              5),
        ("Aristas rol=inferido",        get(roles, "inferido", 0),                3),
        ("Aristas con trafico_mbps",    count(!ismissing, trafico(red)),        170),
        ("Aristas con capacidad_mbps",  count(!ismissing, capacidad(red)),       28),
    ]

    println("\n", "="^62)
    println("VERIFICACIÓN DEL PIPELINE (Anexo A del enunciado)")
    println("="^62)
    todo_ok = true
    for (nombre, obtenido, esperado) in comprobaciones
        ok = obtenido == esperado
        todo_ok &= ok
        @printf("  %-30s %8d  esperado %6d  %s\n",
                nombre, obtenido, esperado, ok ? "OK" : "FALLA")
    end
    println("="^62)
    println(todo_ok ? "  Todas las comprobaciones pasaron." :
                      "  HAY COMPROBACIONES FALLIDAS: revise la carga antes de continuar.")
    println("="^62)
    todo_ok
end

"""
    resumen(red)

Imprime un resumen estructural mínimo. Es solo un punto de partida: el
Problema P1 pide bastante más que esto.
"""
function resumen(red::RedUCuenca)
    g = red.g
    n, m = nv(g), ne(g)
    grados = degree(g)

    println("\n", "="^62)
    println("RESUMEN ESTRUCTURAL (punto de partida del Problema P1)")
    println("="^62)
    @printf("  Nodos                        %d\n", n)
    @printf("  Aristas                      %d\n", m)
    @printf("  Densidad                     %.4f\n", 2m / (n * (n - 1)))
    @printf("  Grado medio                  %.2f\n", mean(grados))
    @printf("  Grado máximo / mínimo        %d / %d\n", maximum(grados), minimum(grados))
    @printf("  Nodos de grado 1             %d\n", count(==(1), grados))
    @printf("  Componentes conexas          %d\n", length(connected_components(g)))

    println("\n  Nodos por campus:")
    for (c, k) in sort(collect(countmap_str(red.campus)), by = x -> -x[2])
        @printf("    %-42s %4d\n", c, k)
    end

    println("\n  Top-10 por grado:")
    for i in sortperm(grados; rev=true)[1:min(10, n)]
        @printf("    %-36s %-12s %4d\n", red.ids[i], red.capa[i], grados[i])
    end
    println("="^62)

    println("""

    SIGUIENTES PASOS (no implementados aquí a propósito):
      P1  centralidades, clustering, diámetro, asortatividad,
          puntos de articulación y puentes
      P2  modelos nulos (Erdos-Renyi, configuración, Barabasi-Albert)
      P3  BFS y DFS implementados desde cero
      P4  Louvain y k-means, comparados con la partición por campus
      P5  Dijkstra y Floyd-Warshall con los tres modelos de peso
      P6  capacidad de los enlaces no documentados + Ford-Fulkerson /
          Edmonds-Karp (reutilizar optimization/ford-fulkerson y
          optimization/edmonds-karp)
      P7  p-mediana y p-centro
      P8  percolación de nodos y de enlaces
      P9  cascadas de fallos y SIR
      P10 ranking de puntos críticos
      P11 propuesta de rediseño con a lo sumo 5 enlaces nuevos
    """)
end

# --- Punto de entrada -----------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    red = cargar_red()
    ok = verificar(red)
    resumen(red)
    exit(ok ? 0 : 1)
end
