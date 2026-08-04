# ============================================================
# Ejemplo 1: Flujo máximo con Edmonds-Karp
#   Red 1: red clásica de Cormen et al. (CLRS), flujo máximo = 23
#   Red 2: red "zigzag" donde elegir mal los caminos es catastrófico
# ============================================================
#
# Ejecutar desde esta carpeta:
#   julia --project=. ejemplo1.jl
#
# O en el REPL / VSCode (recomendado para clase):
#   include("ejemplo1.jl")
#   edmonds_karp_interactivo(red, s, t)          # paso a paso con [Enter]
#   edmonds_karp_interactivo(red_zigzag, 1, 4)

include("edmonds_karp.jl")

# ------------------------------------------------------------
# Red 1: la red clásica de CLRS
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

println("=== Edmonds-Karp: red clásica de CLRS ===\n")
flujo, F, historia = edmonds_karp(red, s, t)

S, aristas_corte = corte_minimo(red.C, F, s)
println("\nCorte mínimo: S = {", join(nombres[S], ", "), "}")
println("Aristas del corte: ",
        join(["$(nombres[u])→$(nombres[v]) ($(C[u,v]))" for (u, v) in aristas_corte], ", "))
println("Capacidad del corte: ", sum(C[u, v] for (u, v) in aristas_corte),
        "  (= flujo máximo, teorema max-flow min-cut)")

println("\nGenerando animación (incluye la onda BFS nivel por nivel)...")
animar_edmonds_karp(red, s, t; archivo="edmonds_karp_clrs.gif", verbose=false)
println("GIF guardado en edmonds_karp_clrs.gif")

# Imagen estática del estado final (para diapositivas)
frames = _fotogramas(red, s, t, historia)
savefig(dibujar_fotograma(red, frames[end]; s=s, t=t), "flujo_maximo_final.png")
println("Imagen final guardada en flujo_maximo_final.png")

# ------------------------------------------------------------
# Red 2: la red "zigzag" — por qué importa elegir el camino corto
# ------------------------------------------------------------
# Con Ford-Fulkerson "puro" (p. ej. DFS), una mala elección puede
# alternar por el arco central u→v (capacidad 1) y necesitar hasta
# 2·M iteraciones (aquí M = 1000, ¡2000 iteraciones!). Edmonds-Karp
# usa siempre el camino más corto y termina en 2 iteraciones.
M = 1000
nombres_z = ["s", "u", "v", "t"]
pos_z = [(0.0, 1.0), (1.2, 2.0), (1.2, 0.0), (2.4, 1.0)]
Cz = zeros(Int, 4, 4)
Cz[1, 2] = M   # s → u
Cz[1, 3] = M   # s → v
Cz[2, 3] = 1   # u → v  (el arco "trampa")
Cz[2, 4] = M   # u → t
Cz[3, 4] = M   # v → t
red_zigzag = RedFlujo(Cz, nombres_z, pos_z)

println("\n=== Edmonds-Karp: red zigzag (M = $M) ===\n")
flujo_z, _, historia_z = edmonds_karp(red_zigzag, 1, 4)
println("\nCon caminos mal elegidos (DFS alternando por u→v), este mismo")
println("flujo podría requerir hasta $(2M) iteraciones; BFS usó $(length(historia_z)).")

animar_edmonds_karp(red_zigzag, 1, 4; archivo="edmonds_karp_zigzag.gif", verbose=false)
println("GIF guardado en edmonds_karp_zigzag.gif")

# ------------------------------------------------------------
# Modo interactivo (descomentar o llamar desde el REPL)
# ------------------------------------------------------------
# edmonds_karp_interactivo(red, s, t)
# edmonds_karp_interactivo(red_zigzag, 1, 4)
