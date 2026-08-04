#!/usr/bin/env python3
"""Código base del proyecto integrador -- Módulo de Redes Complejas (1217).

Carga la red de datos de la Universidad de Cuenca, verifica el pipeline contra
los valores del Anexo A del enunciado e imprime un resumen estructural.

Uso:
    python3 cargar_red.py                # carga desde los CSV
    python3 cargar_red.py --fuente graphml
    python3 cargar_red.py --dir /ruta/a/los/datos

Como módulo:
    from cargar_red import cargar_red, verificar, resumen
    G = cargar_red()
"""

from __future__ import annotations

import argparse
import collections
import os
import pathlib
import sys

import networkx as nx

DIR_SCRIPT = pathlib.Path(__file__).resolve().parent
TESTIGO = "red_ucuenca_nodes.csv"


def _buscar_dir_datos() -> pathlib.Path:
    """Busca hacia arriba el directorio que contiene el conjunto de datos.

    Así el código funciona igual dentro del repositorio
    (``proyecto/codigo_base/`` con los datos en la raíz) que en el paquete
    entregado a los estudiantes (``codigo_base/`` con los datos un nivel
    arriba), sin tener que editar nada.
    """
    for directorio in (DIR_SCRIPT, *DIR_SCRIPT.parents):
        if (directorio / TESTIGO).is_file():
            return directorio
    return DIR_SCRIPT.parent  # sin datos a la vista: falla luego con un error claro


def ruta_datos(nombre: str, base: pathlib.Path | None = None) -> pathlib.Path:
    """Ruta a un archivo del conjunto de datos.

    El directorio puede sobreescribirse con el argumento ``base`` o con la
    variable de entorno ``RED_UCUENCA_DIR``; si no, se busca hacia arriba
    desde la ubicación de este script.
    """
    if base is not None:
        return base / nombre
    entorno = os.environ.get("RED_UCUENCA_DIR")
    raiz = pathlib.Path(entorno) if entorno else _buscar_dir_datos()
    return raiz / nombre


def _leer_csv(ruta: pathlib.Path) -> list[dict[str, str]]:
    import csv

    with open(ruta, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _num(s: str) -> float | None:
    s = (s or "").strip()
    return float(s) if s else None


def cargar_red(fuente: str = "csv", base: pathlib.Path | None = None) -> nx.Graph:
    """Construye la red y la devuelve como ``networkx.Graph``.

    El grafo es **simple y no dirigido**: 177 nodos y 209 aristas, una sola
    componente conexa. Los identificadores son idénticos en los CSV y en el
    GraphML, así que ambas fuentes producen exactamente el mismo grafo.

    Atributos de nodo:  ``label``, ``campus``, ``capa``, ``diagrams``.
    Atributos de arista: ``label``, ``trafico_mbps``, ``capacidad_mbps``,
    ``rol``, ``diagrams``.

    DECISIÓN DE MODELADO: se devuelve el grafo SIN pesos. Elegir el modelo de
    peso (saltos, latencia, carga) es parte del Problema P5; definir la
    capacidad de los enlaces no documentados es parte del Problema P6.

    Parameters
    ----------
    fuente : {"csv", "graphml"}
        Origen de los datos. Ambos deben producir el mismo grafo.
    """
    if fuente not in ("csv", "graphml"):
        raise ValueError(f"fuente desconocida: {fuente!r}")

    if fuente == "graphml":
        G = nx.Graph(nx.read_graphml(ruta_datos("red_ucuenca.graphml", base)))
        # GraphML omite los atributos vacíos y añade un `id` de arista propio.
        # Se normaliza para que ambas fuentes entreguen exactamente el mismo
        # diccionario de atributos: así el código escrito contra una funciona
        # sin cambios contra la otra.
        for _, _, d in G.edges(data=True):
            d.pop("id", None)
            d.setdefault("label", "")
            d.setdefault("trafico_mbps", None)
            d.setdefault("capacidad_mbps", None)
        return G

    G = nx.Graph()
    for fila in _leer_csv(ruta_datos("red_ucuenca_nodes.csv", base)):
        G.add_node(
            fila["id"],
            label=fila["label"],
            campus=fila["campus"],
            capa=fila["capa"],
            diagrams=fila["diagrams"],
        )
    for k, fila in enumerate(_leer_csv(ruta_datos("red_ucuenca_edges.csv", base))):
        u, v = fila["source"], fila["target"]
        for extremo in (u, v):
            if extremo not in G:
                raise KeyError(f"Arista {k}: nodo desconocido {extremo!r}")
        G.add_edge(
            u,
            v,
            label=fila["label"],
            trafico_mbps=_num(fila["trafico_mbps"]),
            capacidad_mbps=_num(fila["capacidad_mbps"]),
            rol=fila["rol"],
            diagrams=fila["diagrams"],
        )
    return G


def verificar(G: nx.Graph) -> bool:
    """Compara la carga con los valores del Anexo A del enunciado."""
    roles = collections.Counter(d["rol"] for _, _, d in G.edges(data=True))
    capas = collections.Counter(nx.get_node_attributes(G, "capa").values())
    con_trafico = sum(
        1 for _, _, d in G.edges(data=True) if d.get("trafico_mbps") is not None
    )
    con_capacidad = sum(
        1 for _, _, d in G.edges(data=True) if d.get("capacidad_mbps") is not None
    )

    comprobaciones = [
        ("Nodos", G.number_of_nodes(), 177),
        ("Aristas", G.number_of_edges(), 209),
        ("Componentes conexas", nx.number_connected_components(G), 1),
        ("Nodos en capa acceso", capas["acceso"], 132),
        ("Nodos en capa agregacion", capas["agregacion"], 27),
        ("Nodos en capa core", capas["core"], 5),
        ("Aristas rol=principal", roles["principal"], 157),
        ("Aristas rol=wan", roles["wan"], 32),
        ("Aristas rol=respaldo", roles["respaldo"], 12),
        ("Aristas rol=secundario", roles["secundario"], 5),
        ("Aristas rol=inferido", roles["inferido"], 3),
        ("Aristas con trafico_mbps", con_trafico, 170),
        ("Aristas con capacidad_mbps", con_capacidad, 28),
    ]

    print("\n" + "=" * 62)
    print("VERIFICACIÓN DEL PIPELINE (Anexo A del enunciado)")
    print("=" * 62)
    todo_ok = True
    for nombre, obtenido, esperado in comprobaciones:
        ok = obtenido == esperado
        todo_ok &= ok
        print(
            f"  {nombre:<30} {obtenido:>8}  esperado {esperado:>6}  "
            f"{'OK' if ok else 'FALLA'}"
        )
    print("=" * 62)
    print(
        "  Todas las comprobaciones pasaron."
        if todo_ok
        else "  HAY COMPROBACIONES FALLIDAS: revise la carga antes de continuar."
    )
    print("=" * 62)
    return todo_ok


def resumen(G: nx.Graph) -> None:
    """Resumen estructural mínimo: punto de partida del Problema P1."""
    n, m = G.number_of_nodes(), G.number_of_edges()
    grados = dict(G.degree())
    campus = collections.Counter(nx.get_node_attributes(G, "campus").values())

    print("\n" + "=" * 62)
    print("RESUMEN ESTRUCTURAL (punto de partida del Problema P1)")
    print("=" * 62)
    print(f"  Nodos                        {n}")
    print(f"  Aristas                      {m}")
    print(f"  Densidad                     {nx.density(G):.4f}")
    print(f"  Grado medio                  {2 * m / n:.2f}")
    print(f"  Grado máximo / mínimo        {max(grados.values())} / {min(grados.values())}")
    print(f"  Nodos de grado 1             {sum(1 for d in grados.values() if d == 1)}")
    print(f"  Componentes conexas          {nx.number_connected_components(G)}")

    print("\n  Nodos por campus:")
    for c, k in campus.most_common():
        print(f"    {c:<42} {k:>4}")

    print("\n  Top-10 por grado:")
    for nodo, d in sorted(grados.items(), key=lambda kv: -kv[1])[:10]:
        capa = G.nodes[nodo].get("capa", "?")
        print(f"    {nodo:<36} {capa:<12} {d:>4}")
    print("=" * 62)

    print(
        """
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
    """
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--fuente",
        choices=("csv", "graphml"),
        default="csv",
        help="origen de los datos (por defecto: csv)",
    )
    parser.add_argument(
        "--dir",
        type=pathlib.Path,
        default=None,
        help="directorio que contiene los archivos red_ucuenca*",
    )
    args = parser.parse_args(argv)

    G = cargar_red(fuente=args.fuente, base=args.dir)
    ok = verificar(G)
    resumen(G)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
