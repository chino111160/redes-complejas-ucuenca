# ============================================================
# K-Means sobre el dataset Palmer Penguins
# ============================================================
#
# Palmer Penguins es un dataset recopilado por la Dra. Kristen Gorman
# en la estación Palmer, Antártida (LTER). Contiene mediciones de
# 344 pingüinos de 3 especies:
#
#   - Adelie      (152 individuos)  — pico corto y ancho
#   - Chinstrap    (68 individuos)  — pico largo y delgado
#   - Gentoo      (124 individuos)  — pico intermedio, mayor masa corporal
#
# Variables numéricas utilizadas como features para K-means:
#   1. bill_length_mm    — largo del pico (mm)
#   2. bill_depth_mm     — profundidad/alto del pico (mm)
#   3. flipper_length_mm — largo de la aleta (mm)
#   4. body_mass_g       — masa corporal (g)
#
# Las escalas son muy distintas (mm vs. g), por lo que se aplica
# normalización Z-score antes de ejecutar K-means.
#
# Objetivo: verificar si K-means (sin usar las etiquetas de especie)
# puede recuperar los 3 grupos naturales solo a partir de las
# mediciones morfológicas.
#
# Referencia:
#   Gorman KB, Williams TD, Fraser WR (2014). Ecological Sexual
#   Dimorphism and Environmental Variability within a Community of
#   Antarctic Penguins (Genus Pygoscelis). PLoS ONE 9(3): e90081.
# ============================================================

using LinearAlgebra
using Statistics
using Random
using DataFrames
using PalmerPenguins
using Plots

# ------------------------------------------------------------
# 1. Funciones del algoritmo (reutilizadas de ejemplo1)
# ------------------------------------------------------------

"""
    dist2(x, μ) -> Float64

Distancia euclidiana al cuadrado: ‖x − μ‖².
"""
dist2(x::AbstractVector, μ::AbstractVector)::Float64 = dot(x .- μ, x .- μ)

"""
    init_plusplus(X, K) -> Vector{Vector{Float64}}

Inicialización K-means++: elige centroides con probabilidad
proporcional a D(x)² para garantizar buena separación inicial.
"""
function init_plusplus(X::Matrix, K::Int)
    n = size(X, 1)
    μ = [X[rand(1:n), :]]
    for _ in 2:K
        D = [minimum(dist2(X[i, :], c) for c in μ) for i in 1:n]
        probs = D ./ sum(D)
        j = findfirst(cumsum(probs) .≥ rand())
        push!(μ, X[j, :])
    end
    return μ
end

"""
    assign_clusters(X, μ) -> Vector{Int}

Asigna cada punto al centroide más cercano.
"""
function assign_clusters(X::Matrix, μ::Vector)::Vector{Int}
    n = size(X, 1)
    labels = Vector{Int}(undef, n)
    for i in 1:n
        dists = [dist2(X[i, :], c) for c in μ]
        labels[i] = argmin(dists)
    end
    return labels
end

"""
    update_centroids(X, labels, K) -> Vector{Vector{Float64}}

Recalcula centroides como la media de cada cluster.
"""
function update_centroids(X::Matrix, labels::Vector{Int}, K::Int)::Vector
    μ_new = Vector(undef, K)
    for k in 1:K
        mask = findall(labels .== k)
        if isempty(mask)
            μ_new[k] = X[rand(1:size(X, 1)), :]
        else
            μ_new[k] = vec(mean(X[mask, :], dims=1))
        end
    end
    return μ_new
end

"""
    wcss(X, μ, labels) -> Float64

Within-Cluster Sum of Squares (inercia).
"""
function wcss(X::Matrix, μ::Vector, labels::Vector{Int})::Float64
    return sum(dist2(X[i, :], μ[labels[i]]) for i in 1:size(X, 1))
end

"""
    my_kmeans(X, K; max_iter=300, tol=1e-6, seed=42) -> NamedTuple

Ejecuta K-means completo con inicialización K-means++.
Retorna `(labels, centroids, wcss_history)`.
"""
function my_kmeans(X::Matrix, K::Int; max_iter=300, tol=1e-6, seed=42)
    Random.seed!(seed)
    μ = init_plusplus(X, K)
    labels = assign_clusters(X, μ)
    J_hist = Float64[]

    for iter in 1:max_iter
        μ_old = deepcopy(μ)
        labels = assign_clusters(X, μ)
        μ = update_centroids(X, labels, K)
        J = wcss(X, μ, labels)
        push!(J_hist, J)

        if all(norm(μ[k] .- μ_old[k]) < tol for k in 1:K)
            println("Converge en iteración $iter")
            break
        end
    end
    return (labels=labels, centroids=μ, wcss_history=J_hist)
end

# ------------------------------------------------------------
# 2. Cargar y preparar datos
# ------------------------------------------------------------

df = DataFrame(PalmerPenguins.load())
dropmissing!(df)

features = [:bill_length_mm, :bill_depth_mm, :flipper_length_mm, :body_mass_g]
X_raw = Matrix{Float64}(df[:, features])

"""
    normalize(X) -> (X_norm, μ_col, σ_col)

Normalización Z-score por columna: x' = (x − μ) / σ.
Necesaria porque las features tienen escalas muy distintas
(mm vs. gramos). Sin normalizar, body_mass_g dominaría la distancia.
"""
function normalize(X::Matrix)
    μ_col = mean(X, dims=1)
    σ_col = std(X, dims=1)
    X_norm = (X .- μ_col) ./ σ_col
    return X_norm, μ_col, σ_col
end

X, μ_col, σ_col = normalize(X_raw)

species_true = df.species
species_names = unique(species_true)
println("Dataset: $(size(X, 1)) pingüinos, $(length(features)) features")
println("Especies reales: ", species_names)

# ------------------------------------------------------------
# 3. Ejecutar K-means (K=3, una por especie)
# ------------------------------------------------------------

K = 3
result = my_kmeans(X, K)
println("WCSS final: ", round(result.wcss_history[end], digits=2))

# ------------------------------------------------------------
# 4. Método del codo: evaluar K=1..8
# ------------------------------------------------------------

"""
Calcula WCSS para distintos valores de K.
El "codo" en la curva sugiere el K óptimo.
"""
K_range = 1:8
wcss_values = Float64[]
for k in K_range
    r = my_kmeans(X, k)
    push!(wcss_values, r.wcss_history[end])
end

# ------------------------------------------------------------
# 5. Visualización
# ------------------------------------------------------------

# Plot 1: Clusters sobre bill_length vs bill_depth (normalizado)
p1 = scatter(X[:, 1], X[:, 2],
    group=result.labels,
    palette=:Set1,
    markersize=5,
    alpha=0.7,
    xlabel="Bill Length (norm)",
    ylabel="Bill Depth (norm)",
    title="K-Means (K=$K)",
    legend=:topright,
    label=permutedims(["Cluster $k" for k in 1:K])
)
cx = [c[1] for c in result.centroids]
cy = [c[2] for c in result.centroids]
scatter!(p1, cx, cy,
    color=:black, markersize=10, markershape=:star5, label="Centroides"
)

# Plot 2: Especies reales (ground truth) para comparar
p2 = scatter(X[:, 1], X[:, 2],
    group=species_true,
    palette=:Dark2,
    markersize=5,
    alpha=0.7,
    xlabel="Bill Length (norm)",
    ylabel="Bill Depth (norm)",
    title="Especies reales",
    legend=:topright
)

# Plot 3: Método del codo
p3 = plot(collect(K_range), wcss_values,
    linewidth=2,
    marker=:circle,
    markersize=5,
    color=:steelblue,
    xlabel="K",
    ylabel="WCSS",
    title="Método del Codo",
    legend=false
)
vline!(p3, [3], linestyle=:dash, color=:red, label="K=3")

# Plot 4: Convergencia
p4 = plot(result.wcss_history,
    linewidth=2,
    marker=:circle,
    markersize=4,
    color=:steelblue,
    xlabel="Iteración",
    ylabel="WCSS",
    title="Convergencia (K=$K)",
    legend=false
)

fig = plot(p1, p2, p3, p4, layout=(2, 2), size=(1100, 900))
savefig(fig, joinpath(@__DIR__, "penguins_kmeans.png"))
println(joinpath(@__DIR__, "penguins_kmeans.png"), " generado")
display(fig)
