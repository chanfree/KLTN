#!/usr/bin/env python3
"""
Export the TF/DEG gene interaction network as a vector PDF
suitable for editing in Adobe Illustrator.

Each node, edge, label, and legend item is a separate editable vector object.
Uses matplotlib's PDF backend with TrueType font embedding (fonttype=42)
so all text remains live and editable in Illustrator.

Usage:
    python3 scripts/export_vector_pdf.py \
        top100_network.cys \
        output/data/top100_tf_go_classification.csv \
        output/pdf/top100_network_GO_TF_network.pdf
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import zipfile
from pathlib import Path

# Add scripts directory to path so we can reuse parsing/layout logic
SCRIPTS_DIR = Path(__file__).parent.resolve()
sys.path.insert(0, str(SCRIPTS_DIR))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle

import build_go_tf_artifacts as bga

# Visual constants — must match build_go_tf_artifacts.py exactly
TF_FILL = "#F3C43C"
TF_BORDER = "#D6A51E"
DEG_FILL = "#B45A9D"
EDGE_COLOR = "#5E8FC7"
LABEL_COLOR = "#164E96"
EDGE_ALPHA = 65 / 255  # ~0.255, matching edge_rgba in render_network

# Layout scale: original plot area is 2110 px wide for 2.0 normalized units
PIXELS_PER_UNIT = 1055.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export the gene interaction network as an Illustrator-editable vector PDF."
    )
    parser.add_argument("input_cys", type=Path, help="Cytoscape session file (.cys)")
    parser.add_argument("classification_csv", type=Path, help="TF/DEG classification CSV")
    parser.add_argument("output_pdf", type=Path, help="Output vector PDF path")
    return parser.parse_args()


def find_entry_flexible(names: list[str], *needles: str) -> str:
    """Find a unique archive entry that contains all needle substrings."""
    matches = [n for n in names if all(needle in n for needle in needles)]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one entry for {needles!r}, got {matches}")
    return matches[0]


def render_vector_pdf(
    output_pdf: Path,
    edge_weights,
    positions: dict[str, tuple[float, float]],
    tf_nodes: list[str],
    deg_nodes: list[str],
    gene_by_suid: dict[str, str],
    degree,
    tf_sizes: dict[str, float],
) -> None:
    """Render the network as a fully vector PDF using matplotlib.

    Every artist (Line2D, Circle, Rectangle, Text) becomes a separate
    editable object when the PDF is opened in Adobe Illustrator.
    """
    # TrueType embedding → live, editable text in Illustrator
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

    # Figure size: 24 x 17 inches matches the original 2400x1700 px canvas at 100 dpi
    fig, ax = plt.subplots(figsize=(24, 17))
    ax.set_xlim(-1.08, 1.08)
    ax.set_ylim(-1.08, 1.08)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # --- Edges: one Line2D per unique node pair ---
    for source, target in edge_weights:
        sx, sy = positions[source]
        tx, ty = positions[target]
        ax.add_line(
            Line2D(
                [sx, tx],
                [sy, ty],
                color=EDGE_COLOR,
                alpha=EDGE_ALPHA,
                linewidth=1.0,
                zorder=1,
            )
        )

    # --- DEG nodes: small purple squares ---
    deg_half = 6.0 / PIXELS_PER_UNIT
    for node in deg_nodes:
        x, y = positions[node]
        ax.add_patch(
            Rectangle(
                (x - deg_half, y - deg_half),
                2 * deg_half,
                2 * deg_half,
                facecolor=DEG_FILL,
                edgecolor="none",
                zorder=2,
            )
        )

    # --- TF nodes: gold circles with darker border ---
    for node in tf_nodes:
        gene = gene_by_suid[node]
        x, y = positions[node]
        radius = tf_sizes[gene] * 0.72 / PIXELS_PER_UNIT
        ax.add_patch(
            Circle(
                (x, y),
                radius,
                facecolor=TF_FILL,
                edgecolor=TF_BORDER,
                linewidth=1.5,
                zorder=3,
            )
        )

    # --- TF labels: live text centred on each TF circle ---
    for node in tf_nodes:
        gene = gene_by_suid[node]
        x, y = positions[node]
        # Original pixel font-size → points (1 px @ 100 dpi = 0.72 pt)
        font_size_pt = (27 + min(13, degree[node] / 2.8)) * 0.72
        ax.text(
            x,
            y,
            gene,
            ha="center",
            va="center",
            fontsize=font_size_pt,
            color=LABEL_COLOR,
            zorder=4,
        )

    # --- Title and subtitle ---
    ax.text(
        -1.04,
        1.055,
        "Regulatory networks of top TFs",
        fontsize=28,
        color="#171717",
        fontweight="bold",
        ha="left",
        va="center",
        zorder=5,
    )
    ax.text(
        -1.035,
        1.005,
        "Gene Ontology Molecular Function: GO:0003700 and descendants",
        fontsize=14,
        color="#555555",
        ha="left",
        va="center",
        zorder=5,
    )

    # --- Legend ---
    leg_x, leg_y = 0.80, -0.85
    # Use an average TF radius for the legend swatch
    avg_tf_radius = sum(tf_sizes.values()) / len(tf_sizes) * 0.72 / PIXELS_PER_UNIT
    ax.add_patch(
        Circle(
            (leg_x, leg_y),
            avg_tf_radius,
            facecolor=TF_FILL,
            edgecolor=TF_BORDER,
            linewidth=1.5,
            zorder=5,
        )
    )
    ax.text(
        leg_x + avg_tf_radius + 0.03,
        leg_y,
        "TFs",
        fontsize=16,
        color="#171717",
        va="center",
        ha="left",
        zorder=5,
    )

    deg_leg_y = leg_y - 0.10
    ax.add_patch(
        Rectangle(
            (leg_x - deg_half, deg_leg_y - deg_half),
            2 * deg_half,
            2 * deg_half,
            facecolor=DEG_FILL,
            edgecolor="none",
            zorder=5,
        )
    )
    ax.text(
        leg_x + avg_tf_radius + 0.03,
        deg_leg_y,
        "DEGs",
        fontsize=16,
        color="#171717",
        va="center",
        ha="left",
        zorder=5,
    )

    # Save as vector PDF (no rasterisation)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        str(output_pdf),
        format="pdf",
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.2,
    )
    plt.close(fig)


def main() -> None:
    args = parse_args()
    classification = bga.load_classification(args.classification_csv)

    with tempfile.TemporaryDirectory(prefix="vector-pdf-") as temp_dir:
        extracted = Path(temp_dir)
        with zipfile.ZipFile(args.input_cys) as archive:
            archive.extractall(extracted)
        names = [
            p.relative_to(extracted).as_posix()
            for p in extracted.rglob("*")
            if p.is_file()
        ]

        node_table_path = extracted / find_entry_flexible(
            names,
            "M.+musculus+%281%29",
            "SHARED_ATTRS",
            "CyNode",
            "shared+node",
        )
        network_xgmml_path = extracted / find_entry_flexible(
            names,
            "/networks/",
            "M.+musculus+%281%29",
            ".xgmml",
        )

        _, nodes = bga.parse_node_table(node_table_path)
        graph_nodes, induced_edges, degree, edge_weights = bga.parse_network(
            network_xgmml_path, nodes
        )

        # Same deterministic layout as build_go_tf_artifacts.py (seed=42)
        positions = bga.compute_layout(graph_nodes, edge_weights)
        bga.separate_tf_positions(positions, nodes, classification)

        gene_by_suid = {suid: row.get("gene name", "") for suid, row in nodes.items()}
        tf_sizes: dict[str, float] = {}
        tf_node_set = set()
        for suid in graph_nodes:
            gene = gene_by_suid[suid]
            info = classification[gene]
            if info["node_type"] == "TF":
                degree_value = float(degree.get(suid, 0))
                tf_sizes[gene] = 42.0 + 0.9 * degree_value
                tf_node_set.add(suid)

        tf_nodes = [s for s in graph_nodes if s in tf_node_set]
        deg_nodes = [s for s in graph_nodes if s not in tf_node_set]

        render_vector_pdf(
            args.output_pdf,
            edge_weights,
            positions,
            tf_nodes,
            deg_nodes,
            gene_by_suid,
            degree,
            tf_sizes,
        )

    print(f"Vector PDF (Illustrator-editable): {args.output_pdf}")
    print(f"  Nodes: {len(graph_nodes)}  (TF: {len(tf_nodes)}, DEG: {len(deg_nodes)})")
    print(f"  Edges: {len(edge_weights)}")
    print(f"  TFs: {', '.join(sorted(tf_sizes))}")


if __name__ == "__main__":
    main()
