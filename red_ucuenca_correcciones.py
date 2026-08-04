#!/usr/bin/env python3
"""Corrección del conjunto de datos de la red de la Universidad de Cuenca.

Toma la transcripción original de los 34 diagramas y produce la versión
saneada que consume el proyecto integrador del módulo de Redes Complejas.
Se conserva en el repositorio para que cada corrección sea auditable y
reproducible.

Correcciones aplicadas (ver red_ucuenca_README.md):

  C1  Identificadores unificados: una sola convención con guion medio, sin
      espacios, idéntica en los CSV y en el GraphML.
  C2  Aristas duplicadas colapsadas: seis pares eran el MISMO enlace físico
      transcrito desde dos diagramas distintos, no enlaces redundantes.
      El grafo pasa a ser simple (209 aristas).
  C3  `bandwidth_mbps` se separa en `trafico_mbps` (medición del weathermap)
      y `capacidad_mbps` (capacidad nominal declarada en la etiqueta). El
      atributo original mezclaba ambas magnitudes.
  C4  `style` (atributo de dibujo, ambiguo) se sustituye por `rol`, con
      semántica explícita: principal / secundario / respaldo / wan / inferido.
  C5  Etiquetas `*UNKNOWN*` del monitor sustituidas por campo vacío.
  C6  `campus` pasa a tener un único valor (ubicación física). Las sedes
      Museo y Centro Histórico dejan de estar bajo «MPLS/Interconexion».
  C7  Nuevo atributo `capa`: core / agregacion / acceso / interconexion / wan.

Uso:
    python3 red_ucuenca_correcciones.py            # escribe los archivos
    python3 red_ucuenca_correcciones.py --dry-run  # solo informa
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import sys
import xml.etree.ElementTree as ET

RAIZ = pathlib.Path(__file__).resolve().parent
ORIG_NODES = RAIZ / "red_ucuenca_nodes.csv"
ORIG_EDGES = RAIZ / "red_ucuenca_edges.csv"

# --- C1: identificadores ---------------------------------------------------


def id_canonico(s: str) -> str:
    """Una sola convención: mayúsculas y guion medio, sin espacios."""
    return re.sub(r"\s+", "-", s.strip())


# --- C6: ubicación física única --------------------------------------------

CAMPUS_EXPLICITO = {
    "PE1-CENTRAL": "Campus Central",
    "PE2-CENTRAL": "Campus Central",
    "FORTIGATE-1800F-CENTRAL": "Campus Central",
    "PE1-BALZAY": "Campus Balzay",
    "PE2-BALZAY": "Campus Balzay",
    "FORTIGATE-1800F-BALZAY": "Campus Balzay",
    "ROUTER-L2TP-BALZAY": "Campus Balzay",
    "ROUTER-CAMPUS-YANUNCAY": "Campus Yanuncay",
    "ROUTER-CAMPUS-PARAISO": "Campus Paraiso",
    "ROUTER-CAMPUS-HUAYNA-CAPAC": "Campus Paraiso",
    "ROUTER-CAMPUS-MUSEO": "Sede Museo",
    "SW-ARUBA-MUSEO": "Sede Museo",
    "ROUTER-CAMPUS-CENTRO-HISTORICO": "Sede Centro Historico",
    "SW-ARUBA-CENTRO-HISTORICO": "Sede Centro Historico",
    "INTERNET-MPLS": "Nube MPLS",
}


def campus_unico(nid: str, campus_original: str) -> str:
    if nid in CAMPUS_EXPLICITO:
        return CAMPUS_EXPLICITO[nid]
    # Los nodos frontera declaraban dos campus; el primero es la ubicación física
    return campus_original.split(",")[0].strip()


# --- C7: capa jerárquica ---------------------------------------------------


def capa(nid: str) -> str:
    if nid == "INTERNET-MPLS" or nid.startswith(("PE1-", "PE2-", "ROUTER-")):
        return "wan"
    if nid.startswith("FORTIGATE-"):
        return "interconexion"
    if nid.startswith("SW-ARUBA-"):
        return "acceso"
    if re.search(r"-C\d+$", nid):
        return "core"
    if re.search(r"-D\d+$", nid) or nid == "CC-BLOQUE-IQ":
        return "agregacion"
    return "acceso"


# --- C3: tráfico medido frente a capacidad nominal -------------------------

# En el diagrama 34 (nube MPLS) la etiqueta declara la CAPACIDAD del enlace
# ("Trafico MPLS 10 Gbps", "L2TP 1 Gbps backup"). En el resto de diagramas la
# etiqueta es una lectura del weathermap, es decir TRÁFICO instantáneo.
RE_GBPS = re.compile(r"(?:(\d+)\s*x\s*)?(\d+(?:[.,]\d+)?)\s*Gbps", re.IGNORECASE)


def capacidad_desde_etiqueta(label: str) -> float | None:
    m = RE_GBPS.search(label or "")
    if not m:
        return None
    factor = int(m.group(1)) if m.group(1) else 1
    return factor * float(m.group(2).replace(",", ".")) * 1000


# --- C4: rol del enlace ----------------------------------------------------

PRIORIDAD_ROL = ["principal", "secundario", "wan", "respaldo", "inferido"]


def rol_enlace(style: str, diagramas: str, trafico: float | None) -> str:
    if style == "inferido":
        return "inferido"
    if "34" in [d.strip() for d in re.split(r"[;,]", diagramas)]:
        return "wan"
    if style == "dashed":
        # Miembro secundario de un par redundante; si no cursa tráfico, está
        # en espera (todas estas etiquetas son "lag NNN 0 bps").
        return "secundario" if trafico else "respaldo"
    return "principal"


# --- Proceso ---------------------------------------------------------------


def cargar_original():
    with open(ORIG_NODES, newline="", encoding="utf-8") as fh:
        nodos = list(csv.DictReader(fh))
    with open(ORIG_EDGES, newline="", encoding="utf-8") as fh:
        aristas = list(csv.DictReader(fh))

    # El script reescribe sus propias entradas: si ya se ejecutó, los archivos
    # tienen el esquema corregido y volver a procesarlos no tiene sentido.
    if aristas and "bandwidth_mbps" not in aristas[0]:
        raise SystemExit(
            "Los archivos ya están corregidos (esquema con trafico_mbps/rol).\n"
            "Para volver a generarlos, recupere la transcripción original:\n"
            "    git checkout HEAD -- red_ucuenca_nodes.csv red_ucuenca_edges.csv\n"
            "y ejecute de nuevo este script."
        )
    return nodos, aristas


def corregir(nodos, aristas):
    informe = []

    # --- nodos ---
    nodos_out = []
    for n in nodos:
        nid = id_canonico(n["id"])
        if nid != n["id"]:
            informe.append(f"C1 identificador: {n['id']!r} -> {nid!r}")
        camp = campus_unico(nid, n["campus"])
        if camp != n["campus"]:
            informe.append(f"C6 campus de {nid}: {n['campus']!r} -> {camp!r}")
        nodos_out.append(
            {
                "id": nid,
                "label": nid,
                "campus": camp,
                "capa": capa(nid),
                "diagrams": n["diagrams"],
            }
        )
    ids = {n["id"] for n in nodos_out}

    # --- aristas: normalizar y reclasificar ---
    intermedias = []
    for k, e in enumerate(aristas):
        u, v = id_canonico(e["source"]), id_canonico(e["target"])
        for extremo in (u, v):
            if extremo not in ids:
                raise KeyError(f"arista {k}: nodo desconocido {extremo!r}")

        label = (e["label"] or "").strip()
        if "*UNKNOWN*" in label:
            informe.append(f"C5 etiqueta {u}--{v}: {label!r} -> ''")
            label = ""

        bruto = (e["bandwidth_mbps"] or "").strip()
        valor = float(bruto) if bruto else None
        rol = rol_enlace(e["style"], e["diagrams"], valor)

        if rol == "wan":
            # En el diagrama MPLS el número de la etiqueta es capacidad
            capacidad = capacidad_desde_etiqueta(label) or valor
            trafico = None
            if valor is not None:
                informe.append(
                    f"C3 {u}--{v}: {valor:g} Mbps reclasificado de "
                    f"trafico a capacidad ({capacidad:g} Mbps)"
                )
        else:
            capacidad = None
            trafico = valor

        intermedias.append(
            {
                "clave": frozenset((u, v)),
                "source": u,
                "target": v,
                "label": label,
                "trafico_mbps": trafico,
                "capacidad_mbps": capacidad,
                "rol": rol,
                "diagrams": e["diagrams"],
            }
        )

    # --- C2: colapsar duplicados ---
    grupos: dict[frozenset, list[dict]] = {}
    for a in intermedias:
        grupos.setdefault(a["clave"], []).append(a)

    aristas_out = []
    vistos = set()
    for a in intermedias:  # preserva el orden original
        if a["clave"] in vistos:
            continue
        vistos.add(a["clave"])
        grupo = grupos[a["clave"]]
        if len(grupo) == 1:
            aristas_out.append(grupo[0])
            continue

        traficos = [g["trafico_mbps"] for g in grupo if g["trafico_mbps"] is not None]
        capacidades = [g["capacidad_mbps"] for g in grupo if g["capacidad_mbps"] is not None]
        diagramas = sorted(
            {d.strip() for g in grupo for d in re.split(r"[;,]", g["diagrams"])},
            key=int,
        )
        # Se conserva la medición máxima: son lecturas del mismo enlace en
        # instantes distintos, no dos enlaces sumables.
        mejor = max(grupo, key=lambda g: (g["trafico_mbps"] or -1))
        informe.append(
            f"C2 duplicado {a['source']}--{a['target']} (diagramas "
            f"{'+'.join(diagramas)}): {len(grupo)} filas -> 1, "
            f"trafico {'/'.join(f'{t:g}' for t in traficos) or '-'} -> "
            f"{max(traficos) if traficos else '-'}"
        )
        aristas_out.append(
            {
                "clave": a["clave"],
                "source": a["source"],
                "target": a["target"],
                "label": mejor["label"],
                "trafico_mbps": max(traficos) if traficos else None,
                "capacidad_mbps": max(capacidades) if capacidades else None,
                "rol": min(
                    (g["rol"] for g in grupo), key=lambda r: PRIORIDAD_ROL.index(r)
                ),
                "diagrams": ";".join(diagramas),
            }
        )

    for a in aristas_out:
        a.pop("clave")
    return nodos_out, aristas_out, informe


def escribir_csv(nodos, aristas):
    with open(RAIZ / "red_ucuenca_nodes.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "label", "campus", "capa", "diagrams"])
        w.writeheader()
        w.writerows(nodos)

    with open(RAIZ / "red_ucuenca_edges.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "source",
                "target",
                "label",
                "trafico_mbps",
                "capacidad_mbps",
                "rol",
                "diagrams",
            ],
        )
        w.writeheader()
        for a in aristas:
            w.writerow(
                {
                    **a,
                    "trafico_mbps": "" if a["trafico_mbps"] is None else f"{a['trafico_mbps']:g}",
                    "capacidad_mbps": "" if a["capacidad_mbps"] is None else f"{a['capacidad_mbps']:g}",
                }
            )


def escribir_graphml(nodos, aristas):
    NS = "http://graphml.graphdrawing.org/xmlns"
    ET.register_namespace("", NS)
    root = ET.Element(f"{{{NS}}}graphml")

    claves = [
        ("label", "node", "label", "string"),
        ("campus", "node", "campus", "string"),
        ("capa", "node", "capa", "string"),
        ("ndiag", "node", "diagrams", "string"),
        ("elabel", "edge", "label", "string"),
        ("trafico", "edge", "trafico_mbps", "double"),
        ("capacidad", "edge", "capacidad_mbps", "double"),
        ("rol", "edge", "rol", "string"),
        ("ediag", "edge", "diagrams", "string"),
    ]
    for kid, dominio, nombre, tipo in claves:
        ET.SubElement(
            root,
            f"{{{NS}}}key",
            {"id": kid, "for": dominio, "attr.name": nombre, "attr.type": tipo},
        )

    g = ET.SubElement(root, f"{{{NS}}}graph", {"id": "UCuenca", "edgedefault": "undirected"})
    for n in nodos:
        el = ET.SubElement(g, f"{{{NS}}}node", {"id": n["id"]})
        for kid, campo in (("label", "label"), ("campus", "campus"),
                           ("capa", "capa"), ("ndiag", "diagrams")):
            ET.SubElement(el, f"{{{NS}}}data", {"key": kid}).text = n[campo]

    for i, a in enumerate(aristas):
        el = ET.SubElement(
            g, f"{{{NS}}}edge", {"id": f"e{i}", "source": a["source"], "target": a["target"]}
        )
        if a["label"]:
            ET.SubElement(el, f"{{{NS}}}data", {"key": "elabel"}).text = a["label"]
        if a["trafico_mbps"] is not None:
            ET.SubElement(el, f"{{{NS}}}data", {"key": "trafico"}).text = f"{a['trafico_mbps']:g}"
        if a["capacidad_mbps"] is not None:
            ET.SubElement(el, f"{{{NS}}}data", {"key": "capacidad"}).text = f"{a['capacidad_mbps']:g}"
        ET.SubElement(el, f"{{{NS}}}data", {"key": "rol"}).text = a["rol"]
        ET.SubElement(el, f"{{{NS}}}data", {"key": "ediag"}).text = a["diagrams"]

    ET.indent(root, space="  ")
    ET.ElementTree(root).write(
        RAIZ / "red_ucuenca.graphml", encoding="utf-8", xml_declaration=True
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true", help="no escribe archivos")
    args = ap.parse_args(argv)

    nodos, aristas = cargar_original()
    nodos_out, aristas_out, informe = corregir(nodos, aristas)

    print(f"nodos:   {len(nodos)} -> {len(nodos_out)}")
    print(f"aristas: {len(aristas)} -> {len(aristas_out)}")
    print(f"\ncorrecciones aplicadas: {len(informe)}")
    for linea in informe:
        print("  ", linea)

    if not args.dry_run:
        escribir_csv(nodos_out, aristas_out)
        escribir_graphml(nodos_out, aristas_out)
        print("\nArchivos reescritos: red_ucuenca_nodes.csv, red_ucuenca_edges.csv, red_ucuenca.graphml")
    return 0


if __name__ == "__main__":
    sys.exit(main())
