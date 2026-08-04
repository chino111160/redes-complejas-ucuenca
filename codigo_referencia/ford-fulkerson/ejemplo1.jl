# ============================================================
# Ejemplo 1: Flujo máximo con Ford-Fulkerson
# Red clásica de Cormen et al. (CLRS), flujo máximo = 23
# ============================================================
#
# Ejecutar desde esta carpeta:
#   julia --project=. ejemplo1.jl
#
# O en el REPL / VSCode (recomendado para clase):
#   include("ejemplo1.jl")
#   ford_fulkerson_interactivo(red, s, t)          # paso a paso con [Enter]
#   ford_fulkerson_interactivo(red, s, t; metodo=:dfs)  # versión clásica DFS

include("ford_fulkerson.jl")

# ------------------------------------------------------------
# Definición de la red
# ------------------------------------------------------------
# Nodos: 1=s, 2=v₁, 3=v₂, 4=v₃, 5=v₄, 6=t
nombres = ["s", "v₁", "v₂", "v₃", "v₄", "t"]
pos = [(0.0, 1.0),   # s
       (1.0, 2.0),   # v₁
       (1.0, 0.0),   # v₂
       (2.2, 2.0),   # v₃
       (2.2, 0.0),   # v₄
       (3.2, 1.0)]   # t

C = zeros(Int, 6, 6)
C[1, 2] = 16   # s  → v₁
C[1, 3] = 13   # s  → v₂
C[2, 4] = 12   # v₁ → v₃
C[3, 2] = 4    # v₂ → v₁
C[3, 5] = 14   # v₂ → v₄
C[4, 3] = 9    # v₃ → v₂
C[4, 6] = 20   # v₃ → t
C[5, 4] = 7    # v₄ → v₃
C[5, 6] = 4    # v₄ → t

red = RedFlujo(C, nombres, pos)
s, t = 1, 6

# ------------------------------------------------------------
# 1) Ejecución con traza en consola
# ------------------------------------------------------------
println("=== Ford-Fulkerson (Edmonds-Karp, BFS) ===")
flujo, F, historia = ford_fulkerson(red, s, t)

S, aristas_corte = corte_minimo(red.C, F, s)
println("\nCorte mínimo: S = {", join(nombres[S], ", "), "}")
println("Aristas del corte: ",
        join(["$(nombres[u])→$(nombres[v]) ($(C[u,v]))" for (u, v) in aristas_corte], ", "))
println("Capacidad del corte: ", sum(C[u, v] for (u, v) in aristas_corte),
        "  (= flujo máximo, teorema max-flow min-cut)")

# ------------------------------------------------------------
# 2) Animación → GIF (para proyectar en clase)
# ------------------------------------------------------------
println("\nGenerando animación...")
animar_ford_fulkerson(red, s, t; archivo="ford_fulkerson_bfs.gif", verbose=false)
println("GIF guardado en ford_fulkerson_bfs.gif")

# Comparación: la versión clásica con DFS puede tomar otro camino
animar_ford_fulkerson(red, s, t; metodo=:dfs,
                      archivo="ford_fulkerson_dfs.gif", verbose=false)
println("GIF guardado en ford_fulkerson_dfs.gif")

# Imagen estática del estado final con el corte mínimo (para diapositivas)
frames = _fotogramas(red, s, t, historia)
savefig(dibujar_fotograma(red, frames[end]; s=s, t=t), "flujo_maximo_final.png")
println("Imagen final guardada en flujo_maximo_final.png")

# ------------------------------------------------------------
# 3) Modo interactivo (descomentar o llamar desde el REPL)
# ------------------------------------------------------------
# ford_fulkerson_interactivo(red, s, t)
