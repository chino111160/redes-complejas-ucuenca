# ============================================================
# K-Means Clustering desde cero en Julia
# ============================================================

using LinearAlgebra
using Statistics
using Random
using Plots

# ------------------------------------------------------------
# 1. Funciones de distancia
# ------------------------------------------------------------

"""
    dist2(x, μ) -> Float64

Calcula la distancia euclidiana al cuadrado entre dos vectores: ‖x − μ‖².
Se usa el cuadrado para evitar la raíz cuadrada (no afecta la comparación
de distancias y es más eficiente).
"""
dist2(x::AbstractVector, μ::AbstractVector)::Float64 = dot(x .- μ, x .- μ)

# ------------------------------------------------------------
# 2. Inicialización de centroides
# ------------------------------------------------------------

"""
    init_centroids(X, K) -> Vector{Vector{Float64}}

Inicialización aleatoria simple: selecciona K puntos del dataset
al azar como centroides iniciales. Es rápida pero puede producir
inicializaciones malas si los puntos elegidos quedan en el mismo cluster.
"""
function init_centroids(X::Matrix, K::Int)
    n = size(X, 1)
    idx = randperm(n)[1:K]
    return [X[i, :] for i in idx]
end

"""
    init_plusplus(X, K) -> Vector{Vector{Float64}}

Inicialización K-means++: el primer centroide se elige al azar.
Los siguientes se eligen con probabilidad proporcional a D(x)²,
donde D(x) es la distancia al centroide más cercano ya elegido.
Esto garantiza centroides bien separados y mejora la convergencia.
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

# ------------------------------------------------------------
# 3. Pasos del algoritmo
# ------------------------------------------------------------

"""
    assign_clusters(X, μ) -> Vector{Int}

Paso de asignación: para cada punto xᵢ, calcula la distancia a todos
los centroides y le asigna la etiqueta del centroide más cercano.
Complejidad: O(n × K × d).
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

Paso de actualización: recalcula cada centroide como la media aritmética
de todos los puntos asignados a su cluster. Si un cluster queda vacío,
se reinicializa eligiendo un punto aleatorio del dataset.
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

# ------------------------------------------------------------
# 4. Métrica de evaluación
# ------------------------------------------------------------

"""
    wcss(X, μ, labels) -> Float64

Within-Cluster Sum of Squares (inercia): suma de las distancias al
cuadrado de cada punto a su centroide asignado.
J = Σᵢ ‖xᵢ − μ_label(i)‖². Un valor menor indica clusters más compactos.
"""
function wcss(X::Matrix, μ::Vector, labels::Vector{Int})::Float64
    return sum(dist2(X[i, :], μ[labels[i]]) for i in 1:size(X, 1))
end

# ------------------------------------------------------------
# 5. Algoritmo principal
# ------------------------------------------------------------

"""
    my_kmeans(X, K; max_iter=300, tol=1e-6, seed=42) -> NamedTuple

Ejecuta el algoritmo K-means completo:
  1. Inicializa centroides con K-means++.
  2. Itera asignación → actualización hasta convergencia.
  3. Converge cuando ningún centroide se mueve más que `tol`.

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
# 6. Generación de datos (blobs) y ejecución
# ------------------------------------------------------------

Random.seed!(42)
n, d, K = 60, 2, 3

centers = [0.0 0.0;     # centro cluster 1
           5.0 5.0;     # centro cluster 2
          -5.0 5.0]     # centro cluster 3
points_per_cluster = n ÷ K
X = vcat([randn(points_per_cluster, d) .+ centers[k, :]' for k in 1:K]...)

result = my_kmeans(X, K)
println("Labels:     ", result.labels)
println("Centroids:  ", result.centroids)
println("WCSS final: ", result.wcss_history[end])

# ------------------------------------------------------------
# 7. Visualización
# ------------------------------------------------------------

# Plot 1: Clusters con centroides
p1 = scatter(X[:, 1], X[:, 2],
    group=result.labels,
    palette=:Set1,
    markersize=6,
    alpha=0.7,
    xlabel="x₁", ylabel="x₂",
    title="K-Means Clustering (K=$K)",
    legend=:topright,
    label=permutedims(["Cluster $k" for k in 1:K])
)
cx = [c[1] for c in result.centroids]
cy = [c[2] for c in result.centroids]
scatter!(p1, cx, cy,
    color=:black,
    markersize=12,
    markershape=:star5,
    label="Centroides"
)

# Plot 2: Convergencia de WCSS por iteración
p2 = plot(result.wcss_history,
    linewidth=2,
    marker=:circle,
    markersize=4,
    color=:steelblue,
    xlabel="Iteración",
    ylabel="WCSS (inercia)",
    title="Convergencia de K-Means",
    legend=false
)

# Combinar ambos plots
fig = plot(p1, p2, layout=(1, 2), size=(1000, 450))
savefig(fig, joinpath(@__DIR__, "kmeans_resultado.png"))
println(@__DIR__, "/kmeans_resultado.png generado")
display(fig)
