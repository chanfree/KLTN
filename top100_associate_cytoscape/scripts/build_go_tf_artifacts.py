#!/usr/bin/env python3

from __future__ import annotations

import argparse
import collections
import csv
import math
import shutil
import tempfile
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from PIL import Image as PILImage
from PIL import ImageDraw, ImageFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


XGMML_NS = "http://www.cs.rpi.edu/XGMML"
CY_NS = "http://www.cytoscape.org"
ET.register_namespace("", XGMML_NS)
ET.register_namespace("cy", CY_NS)

TF_FILL = "#F3C43C"
TF_BORDER = "#D6A51E"
DEG_FILL = "#B45A9D"
EDGE_COLOR = "#5E8FC7"
LABEL_COLOR = "#164E96"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_cys", type=Path)
    parser.add_argument("classification_csv", type=Path)
    parser.add_argument("--output-cys", type=Path, required=True)
    parser.add_argument("--backup-cys", type=Path, required=True)
    parser.add_argument("--preview-png", type=Path, required=True)
    parser.add_argument("--report-pdf", type=Path, required=True)
    return parser.parse_args()


def find_entry(names: list[str], needle: str) -> str:
    matches = [name for name in names if needle in name]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one archive entry for {needle!r}, got {matches}")
    return matches[0]


def read_cycsv(path: Path) -> list[list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.reader(stream))


def write_cycsv(path: Path, rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerows(rows)


def load_classification(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return {row["gene_name"]: row for row in csv.DictReader(stream)}


def parse_node_table(path: Path) -> tuple[list[list[str]], dict[str, dict[str, str]]]:
    rows = read_cycsv(path)
    header = rows[1]
    data = {}
    for row in rows[5:]:
        row += [""] * (len(header) - len(row))
        data[row[0]] = dict(zip(header, row))
    return rows, data


def parse_network(
    network_path: Path,
    nodes: dict[str, dict[str, str]],
) -> tuple[set[str], list[tuple[str, str, str]], collections.Counter[str], collections.Counter[tuple[str, str]]]:
    tree = ET.parse(network_path)
    root = tree.getroot()
    query_suids = {
        suid for suid, row in nodes.items() if row.get("node type") == "query"
    }
    induced_edges = []
    degree = collections.Counter()
    edge_weights = collections.Counter()
    for edge in root.findall(f".//{{{XGMML_NS}}}edge"):
        source = edge.get("source", "")
        target = edge.get("target", "")
        if source not in query_suids or target not in query_suids:
            continue
        induced_edges.append((source, target, edge.get("id", "")))
        degree[source] += 1
        degree[target] += 1
        edge_weights[tuple(sorted((source, target)))] += 1
    return query_suids, induced_edges, degree, edge_weights


def edge_type_counts(
    edge_table_path: Path,
    induced_edges: list[tuple[str, str, str]],
) -> collections.Counter[str]:
    rows = read_cycsv(edge_table_path)
    header = rows[1]
    type_by_suid = {}
    for row in rows[5:]:
        row += [""] * (len(header) - len(row))
        record = dict(zip(header, row))
        type_by_suid[record.get("SUID", "")] = record.get("data type", "Other") or "Other"
    return collections.Counter(type_by_suid.get(edge_id, "Other") for _, _, edge_id in induced_edges)


def update_node_table(
    path: Path,
    classification: dict[str, dict[str, str]],
    degree: collections.Counter[str],
) -> dict[str, float]:
    rows = read_cycsv(path)
    header = rows[1]
    original_len = len(header)
    new_columns = [
        "node_type",
        "tf_database_type",
        "tf_database_version",
        "tf_reference_rule",
        "tf_go_terms",
        "tf_go_term_names",
        "tf_evidence_codes",
        "tf_classification_date",
        "degree_val",
        "display_size",
        "label_text",
    ]
    for column in new_columns:
        if column in header:
            raise RuntimeError(f"Column already exists: {column}")
    rows[1].extend(new_columns)
    rows[2].extend([
        "java.lang.String",
        "java.lang.String",
        "java.lang.String",
        "java.lang.String",
        "java.lang.String",
        "java.lang.String",
        "java.lang.String",
        "java.lang.String",
        "java.lang.Double",
        "java.lang.Double",
        "java.lang.String",
    ])
    rows[3].extend(["mutable"] * len(new_columns))
    rows[4].extend([""] * len(new_columns))

    sizes = {}
    gene_index = header.index("gene name")
    source_type_index = header.index("node type")
    for row in rows[5:]:
        row += [""] * (original_len - len(row))
        suid = row[0]
        gene = row[gene_index]
        source_type = row[source_type_index]
        if source_type == "query":
            info = classification[gene]
            node_type = info["node_type"]
            degree_value = float(degree.get(suid, 0))
            display_size = 42.0 + 0.9 * degree_value if node_type == "TF" else 8.0
            values = [
                node_type,
                info["tf_database_type"],
                info["tf_database_version"],
                info["tf_reference_rule"],
                info["tf_go_terms"],
                info["tf_go_term_names"],
                info["tf_evidence_codes"],
                info["classification_date"],
                f"{degree_value:.1f}",
                f"{display_size:.1f}",
                gene if node_type == "TF" else "",
            ]
            if node_type == "TF":
                sizes[gene] = display_size
        else:
            values = ["BACKGROUND", "", "", "", "", "", "", "", "0.0", "0.0", ""]
        row.extend(values)
    write_cycsv(path, rows)
    return sizes


def update_network_name(path: Path, new_name: str) -> None:
    rows = read_cycsv(path)
    header = rows[1]
    for row in rows[5:]:
        row += [""] * (len(header) - len(row))
        for column in ("name", "shared name"):
            if column in header:
                row[header.index(column)] = new_name
    write_cycsv(path, rows)


def replace_mapping(
    visual_property: ET.Element,
    mapping_type: str,
    attribute_name: str,
    attribute_type: str,
    entries: list[tuple[str, str]] | None = None,
) -> None:
    for child in list(visual_property):
        visual_property.remove(child)
    mapping = ET.SubElement(
        visual_property,
        mapping_type,
        {"attributeName": attribute_name, "attributeType": attribute_type},
    )
    if mapping_type == "discreteMapping":
        for key, value in entries or []:
            ET.SubElement(
                mapping,
                "discreteMappingEntry",
                {"attributeValue": key, "value": value},
            )


def update_visual_style(path: Path, tf_sizes: dict[str, float]) -> str:
    tree = ET.parse(path)
    root = tree.getroot()
    style = next(s for s in root.findall("visualStyle") if s.get("name") == "M musculus (1)")
    style_name = "Top100 TF - GO MF"
    style.set("name", style_name)

    def vp(section: str, name: str) -> ET.Element:
        parent = style.find(section)
        assert parent is not None
        existing = next(
            (x for x in parent.findall("visualProperty") if x.get("name") == name),
            None,
        )
        if existing is not None:
            return existing
        return ET.SubElement(parent, "visualProperty", {"name": name})

    vp("network", "NETWORK_BACKGROUND_PAINT").set("default", "#FFFFFF")

    node_visible = vp("node", "NODE_VISIBLE")
    node_visible.set("default", "true")
    replace_mapping(
        node_visible,
        "discreteMapping",
        "node_type",
        "string",
        [("TF", "true"), ("DEG", "true"), ("BACKGROUND", "false")],
    )

    node_shape = vp("node", "NODE_SHAPE")
    node_shape.set("default", "RECTANGLE")
    replace_mapping(
        node_shape,
        "discreteMapping",
        "node_type",
        "string",
        [("TF", "ELLIPSE"), ("DEG", "RECTANGLE"), ("BACKGROUND", "RECTANGLE")],
    )

    node_fill = vp("node", "NODE_FILL_COLOR")
    node_fill.set("default", DEG_FILL)
    replace_mapping(
        node_fill,
        "discreteMapping",
        "node_type",
        "string",
        [("TF", TF_FILL), ("DEG", DEG_FILL), ("BACKGROUND", "#D5D8DC")],
    )

    node_border = vp("node", "NODE_BORDER_PAINT")
    node_border.set("default", DEG_FILL)
    replace_mapping(
        node_border,
        "discreteMapping",
        "node_type",
        "string",
        [("TF", TF_BORDER), ("DEG", DEG_FILL), ("BACKGROUND", "#D5D8DC")],
    )

    border_width = vp("node", "NODE_BORDER_WIDTH")
    border_width.set("default", "0.0")
    replace_mapping(
        border_width,
        "discreteMapping",
        "node_type",
        "string",
        [("TF", "1.5"), ("DEG", "0.0"), ("BACKGROUND", "0.0")],
    )

    node_label = vp("node", "NODE_LABEL")
    node_label.set("default", "")
    replace_mapping(node_label, "passthroughMapping", "label_text", "string")
    vp("node", "NODE_LABEL_COLOR").set("default", LABEL_COLOR)
    vp("node", "NODE_LABEL_POSITION").set("default", "C,C,c,0.00,0.00")

    label_size = vp("node", "NODE_LABEL_FONT_SIZE")
    label_size.set("default", "1")
    replace_mapping(
        label_size,
        "discreteMapping",
        "node_type",
        "string",
        [("TF", "16"), ("DEG", "1"), ("BACKGROUND", "1")],
    )

    node_size = vp("node", "NODE_SIZE")
    node_size.set("default", "8.0")
    replace_mapping(
        node_size,
        "discreteMapping",
        "gene name",
        "string",
        [(gene, f"{size:.1f}") for gene, size in sorted(tf_sizes.items())],
    )
    vp("node", "NODE_WIDTH").set("default", "8.0")
    vp("node", "NODE_HEIGHT").set("default", "8.0")
    vp("node", "NODE_TRANSPARENCY").set("default", "255")

    for name in ("EDGE_UNSELECTED_PAINT", "EDGE_STROKE_UNSELECTED_PAINT"):
        prop = vp("edge", name)
        prop.set("default", EDGE_COLOR)
        for child in list(prop):
            prop.remove(child)
    vp("edge", "EDGE_TRANSPARENCY").set("default", "105")
    width = vp("edge", "EDGE_WIDTH")
    width.set("default", "0.75")
    for child in list(width):
        width.remove(child)
    vp("edge", "EDGE_CURVED").set("default", "true")
    vp("edge", "EDGE_SOURCE_ARROW_SHAPE").set("default", "NONE")
    vp("edge", "EDGE_TARGET_ARROW_SHAPE").set("default", "NONE")

    tree.write(path, encoding="UTF-8", xml_declaration=True)
    return style_name


def compute_layout(
    graph_nodes: set[str],
    edge_weights: collections.Counter[tuple[str, str]],
) -> dict[str, tuple[float, float]]:
    connected_nodes = {node for edge in edge_weights for node in edge}
    isolated_nodes = sorted(graph_nodes - connected_nodes, key=int)
    node_order = sorted(connected_nodes, key=int)
    index = {node: i for i, node in enumerate(node_order)}
    rng = np.random.default_rng(42)
    pos = rng.normal(0.0, 0.25, size=(len(node_order), 2))
    edge_pairs = np.array(
        [(index[source], index[target]) for source, target in edge_weights],
        dtype=int,
    )
    weights = np.array(
        [math.sqrt(edge_weights[edge]) for edge in edge_weights],
        dtype=float,
    )
    k = 0.30
    temperature = 0.32
    attraction_strength = 0.03

    for iteration in range(650):
        delta = pos[:, None, :] - pos[None, :, :]
        distance = np.linalg.norm(delta, axis=2)
        np.fill_diagonal(distance, np.inf)
        repulsive = np.sum(
            delta / distance[:, :, None] * (k * k / distance)[:, :, None],
            axis=1,
        )

        displacement = repulsive
        if len(edge_pairs):
            source_index = edge_pairs[:, 0]
            target_index = edge_pairs[:, 1]
            edge_delta = pos[source_index] - pos[target_index]
            edge_distance = np.maximum(np.linalg.norm(edge_delta, axis=1), 1e-6)
            attraction = (
                edge_delta
                / edge_distance[:, None]
                * ((edge_distance * edge_distance / k) * weights * attraction_strength)[:, None]
            )
            np.add.at(displacement, source_index, -attraction)
            np.add.at(displacement, target_index, attraction)

        displacement_norm = np.maximum(np.linalg.norm(displacement, axis=1), 1e-9)
        step = displacement / displacement_norm[:, None] * np.minimum(displacement_norm, temperature)[:, None]
        pos += step
        pos -= pos.mean(axis=0)
        temperature = max(0.004, 0.32 * (1.0 - (iteration + 1) / 650.0))

    max_abs = max(float(np.max(np.abs(pos))), 1e-9)
    pos = pos / max_abs * 0.82
    result = {node: (float(pos[i, 0]), float(pos[i, 1])) for i, node in enumerate(node_order)}
    for i, node in enumerate(isolated_nodes):
        angle = -math.pi / 3.0 + i * (2.0 * math.pi / max(1, len(isolated_nodes)))
        result[node] = (0.92 * math.cos(angle), 0.92 * math.sin(angle))
    return result


def separate_tf_positions(
    positions: dict[str, tuple[float, float]],
    nodes: dict[str, dict[str, str]],
    classification: dict[str, dict[str, str]],
) -> None:
    tf_nodes = [
        suid
        for suid, row in nodes.items()
        if row.get("gene name") in classification
        and classification[row["gene name"]]["node_type"] == "TF"
        and suid in positions
    ]
    for _ in range(40):
        changed = False
        for i, source in enumerate(tf_nodes):
            for target in tf_nodes[i + 1 :]:
                sx, sy = positions[source]
                tx, ty = positions[target]
                dx, dy = sx - tx, sy - ty
                distance = math.hypot(dx, dy)
                minimum_distance = 0.15
                if distance >= minimum_distance:
                    continue
                if distance < 1e-6:
                    dx, dy, distance = 1.0, 0.0, 1.0
                shift = (minimum_distance - distance) / 2.0
                ux, uy = dx / distance, dy / distance
                positions[source] = (
                    max(-0.90, min(0.90, sx + ux * shift)),
                    max(-0.90, min(0.90, sy + uy * shift)),
                )
                positions[target] = (
                    max(-0.90, min(0.90, tx - ux * shift)),
                    max(-0.90, min(0.90, ty - uy * shift)),
                )
                changed = True
        if not changed:
            break


def update_view(path: Path, positions: dict[str, tuple[float, float]], style_name: str) -> None:
    tree = ET.parse(path)
    root = tree.getroot()
    root.set(f"{{{CY_NS}}}visualStyle", style_name)
    graphics = root.find(f"{{{XGMML_NS}}}graphics")
    if graphics is not None:
        for att in graphics.findall(f"{{{XGMML_NS}}}att"):
            name = att.get("name")
            if name == "NETWORK_CENTER_X_LOCATION" or name == "NETWORK_CENTER_Y_LOCATION":
                att.set("value", "0.0")
            elif name == "NETWORK_SCALE_FACTOR":
                att.set("value", "0.50")
            elif name == "NETWORK_WIDTH":
                att.set("value", "900.0")
            elif name == "NETWORK_HEIGHT":
                att.set("value", "700.0")

    for node in root.findall(f"{{{XGMML_NS}}}node"):
        suid = node.get(f"{{{CY_NS}}}nodeId", "")
        if suid not in positions:
            continue
        x, y = positions[suid]
        node_graphics = node.find(f"{{{XGMML_NS}}}graphics")
        if node_graphics is not None:
            node_graphics.set("x", f"{x * 520.0:.6f}")
            node_graphics.set("y", f"{y * 470.0:.6f}")
    tree.write(path, encoding="UTF-8", xml_declaration=True)


def render_network(
    output_png: Path,
    graph_nodes: set[str],
    edge_weights: collections.Counter[tuple[str, str]],
    positions: dict[str, tuple[float, float]],
    nodes: dict[str, dict[str, str]],
    classification: dict[str, dict[str, str]],
    degree: collections.Counter[str],
    tf_sizes: dict[str, float],
) -> None:
    gene_by_suid = {suid: row.get("gene name", "") for suid, row in nodes.items()}
    tf_nodes = [suid for suid in graph_nodes if classification[gene_by_suid[suid]]["node_type"] == "TF"]
    deg_nodes = [suid for suid in graph_nodes if suid not in tf_nodes]

    width, height = 2400, 1700
    image = PILImage.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    overlay = PILImage.new("RGBA", (width, height), (255, 255, 255, 0))
    edge_draw = ImageDraw.Draw(overlay)
    plot_left, plot_top, plot_right, plot_bottom = 130, 190, 2240, 1570

    def project(node: str) -> tuple[float, float]:
        x, y = positions[node]
        px = plot_left + (x + 1.0) * 0.5 * (plot_right - plot_left)
        py = plot_top + (1.0 - (y + 1.0) * 0.5) * (plot_bottom - plot_top)
        return px, py

    edge_rgba = (94, 143, 199, 65)
    for source, target in edge_weights:
        edge_draw.line([project(source), project(target)], fill=edge_rgba, width=2)
    image = PILImage.alpha_composite(image.convert("RGBA"), overlay).convert("RGB")
    draw = ImageDraw.Draw(image)

    for node in deg_nodes:
        x, y = project(node)
        half = 6
        draw.rectangle((x - half, y - half, x + half, y + half), fill=DEG_FILL)
    for node in tf_nodes:
        gene = gene_by_suid[node]
        x, y = project(node)
        radius = tf_sizes[gene] * 0.72
        draw.ellipse(
            (x - radius, y - radius, x + radius, y + radius),
            fill=TF_FILL,
            outline=TF_BORDER,
            width=3,
        )

    font_regular_path = Path("/System/Library/Fonts/Supplemental/Arial.ttf")
    font_bold_path = Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf")
    if not font_regular_path.exists():
        font_regular_path = Path("/System/Library/Fonts/Helvetica.ttc")
    if not font_bold_path.exists():
        font_bold_path = font_regular_path
    title_font = ImageFont.truetype(str(font_regular_path), 58)
    subtitle_font = ImageFont.truetype(str(font_regular_path), 25)
    legend_font = ImageFont.truetype(str(font_regular_path), 34)
    draw.text((70, 42), "Regulatory networks of top TFs", fill="#171717", font=title_font)
    draw.text(
        (72, 112),
        "Gene Ontology Molecular Function: GO:0003700 and descendants",
        fill="#555555",
        font=subtitle_font,
    )
    for node in tf_nodes:
        gene = gene_by_suid[node]
        x, y = project(node)
        font_size = int(27 + min(13, degree[node] / 2.8))
        label_font = ImageFont.truetype(str(font_regular_path), font_size)
        box = draw.textbbox((0, 0), gene, font=label_font)
        text_width = box[2] - box[0]
        text_height = box[3] - box[1]
        draw.text(
            (x - text_width / 2, y - text_height / 2 - box[1]),
            gene,
            fill=LABEL_COLOR,
            font=label_font,
        )

    legend_x, legend_y = 2040, 1390
    draw.ellipse((legend_x, legend_y, legend_x + 48, legend_y + 32), fill=TF_FILL, outline=TF_BORDER, width=2)
    draw.text((legend_x + 65, legend_y - 5), "TFs", fill="#171717", font=legend_font)
    draw.rectangle((legend_x + 7, legend_y + 72, legend_x + 39, legend_y + 96), fill=DEG_FILL)
    draw.text((legend_x + 65, legend_y + 62), "DEGs", fill="#171717", font=legend_font)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    image.save(output_png, format="PNG", optimize=True)


def make_report(
    output_pdf: Path,
    preview_png: Path,
    classification: dict[str, dict[str, str]],
    nodes: dict[str, dict[str, str]],
    degree: collections.Counter[str],
    edge_counts: collections.Counter[str],
    tf_sizes: dict[str, float],
    graph_nodes: set[str],
    induced_edge_count: int,
) -> None:
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))
    page_size = landscape(A4)
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "ChineseTitle",
        parent=styles["Title"],
        fontName="STSong-Light",
        fontSize=22,
        leading=28,
        textColor=colors.HexColor("#17202A"),
        alignment=TA_LEFT,
        spaceAfter=8,
    )
    h2_style = ParagraphStyle(
        "ChineseH2",
        parent=styles["Heading2"],
        fontName="STSong-Light",
        fontSize=15,
        leading=20,
        textColor=colors.HexColor("#164E96"),
        spaceBefore=6,
        spaceAfter=7,
    )
    body_style = ParagraphStyle(
        "ChineseBody",
        parent=styles["BodyText"],
        fontName="STSong-Light",
        fontSize=10.2,
        leading=15,
        textColor=colors.HexColor("#263238"),
        alignment=TA_LEFT,
    )
    small_style = ParagraphStyle(
        "ChineseSmall",
        parent=body_style,
        fontSize=8.2,
        leading=11,
        textColor=colors.HexColor("#52606D"),
    )
    table_header = ParagraphStyle(
        "TableHeader",
        parent=body_style,
        fontSize=8.5,
        leading=10,
        textColor=colors.white,
        alignment=TA_CENTER,
    )
    table_cell = ParagraphStyle(
        "TableCell",
        parent=body_style,
        fontSize=8.2,
        leading=10,
    )

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor("#D9E2EC"))
        canvas.line(18 * mm, 12 * mm, page_size[0] - 18 * mm, 12 * mm)
        canvas.setFont("STSong-Light", 8)
        canvas.setFillColor(colors.HexColor("#64748B"))
        canvas.drawString(18 * mm, 7.5 * mm, "top100_network.cys - GO Molecular Function TF analysis")
        canvas.drawRightString(page_size[0] - 18 * mm, 7.5 * mm, f"第 {doc.page} 页")
        canvas.restoreState()

    doc = SimpleDocTemplate(
        str(output_pdf),
        pagesize=page_size,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=14 * mm,
        bottomMargin=17 * mm,
        title="top100_network GO Molecular Function TF analysis",
        author="Codex",
    )
    story = []

    intro_title = Paragraph("基因互作网络与 TF 核实结果", title_style)
    intro_text = Paragraph(
        "节点分类数据库：Gene Ontology - Molecular Function；严格 TF 判据为 GO:0003700 "
        "DNA-binding transcription factor activity 及其所有子条目。",
        body_style,
    )
    image = Image(str(preview_png))
    max_w = page_size[0] - 36 * mm
    max_h = 125 * mm
    ratio = min(max_w / image.imageWidth, max_h / image.imageHeight)
    image.drawWidth = image.imageWidth * ratio
    image.drawHeight = image.imageHeight * ratio
    story.append(KeepTogether([intro_title, intro_text, Spacer(1, 3 * mm), image]))
    story.append(PageBreak())

    tf_rows = []
    gene_to_suid = {
        row.get("gene name", ""): suid
        for suid, row in nodes.items()
        if row.get("node type") == "query"
    }
    for gene, info in classification.items():
        if info["node_type"] != "TF":
            continue
        suid = gene_to_suid[gene]
        tf_rows.append({
            "gene": gene,
            "degree": degree[suid],
            "centrality": degree[suid] / max(1, len(graph_nodes) - 1),
            "size": tf_sizes[gene],
            "go": info["tf_go_terms"],
            "evidence": info["tf_evidence_codes"],
        })
    tf_rows.sort(key=lambda row: (-row["degree"], row["gene"]))

    story.append(Paragraph("分析总结", title_style))
    metrics = [
        ["候选基因", "可见互作边", "严格 TF", "其他 DEG"],
        [
            str(len(graph_nodes)),
            f"{induced_edge_count:,}",
            str(len(tf_rows)),
            str(len(graph_nodes) - len(tf_rows)),
        ],
    ]
    metric_table = Table(metrics, colWidths=[48 * mm] * 4, rowHeights=[10 * mm, 14 * mm])
    metric_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#164E96")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("BACKGROUND", (0, 1), (-1, 1), colors.HexColor("#F5F8FB")),
        ("TEXTCOLOR", (0, 1), (-1, 1), colors.HexColor("#17202A")),
        ("FONTNAME", (0, 0), (-1, -1), "STSong-Light"),
        ("FONTSIZE", (0, 0), (-1, 0), 9),
        ("FONTSIZE", (0, 1), (-1, 1), 16),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#D9E2EC")),
    ]))
    story.append(metric_table)
    story.append(Spacer(1, 5 * mm))
    story.append(Paragraph(
        "重新核实后，严格 TF 集合为：" + "、".join(row["gene"] for row in tf_rows) + "。"
        "此前标注中的 Nupr1 和 Eya1 未满足 GO:0003700 严格判据；Etv1、Pax9 和 Sox21 被新增为 TF。",
        body_style,
    ))
    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph("TF 节点及证据", h2_style))
    tf_table_data = [[
        Paragraph("基因", table_header),
        Paragraph("度数", table_header),
        Paragraph("度中心性", table_header),
        Paragraph("节点大小", table_header),
        Paragraph("GO 条目", table_header),
        Paragraph("证据代码", table_header),
    ]]
    for row in tf_rows:
        tf_table_data.append([
            Paragraph(row["gene"], table_cell),
            Paragraph(str(row["degree"]), table_cell),
            Paragraph(f"{row['centrality']:.3f}", table_cell),
            Paragraph(f"{row['size']:.1f}", table_cell),
            Paragraph(row["go"].replace("; ", "<br/>"), table_cell),
            Paragraph(row["evidence"].replace("; ", ", "), table_cell),
        ])
    tf_table = Table(
        tf_table_data,
        colWidths=[24 * mm, 18 * mm, 24 * mm, 24 * mm, 72 * mm, 34 * mm],
        repeatRows=1,
    )
    tf_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#164E96")),
        ("BACKGROUND", (0, 1), (-1, -1), colors.white),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7F9FB")]),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#D9E2EC")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (1, 1), (3, -1), "CENTER"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    story.append(tf_table)
    story.append(PageBreak())

    story.append(Paragraph("方法、数据库和网络结构", title_style))
    database_version = next(iter(classification.values()))["tf_database_version"]
    method_rows = [
        [Paragraph("项目", table_header), Paragraph("定义", table_header)],
        [Paragraph("物种", table_cell), Paragraph("Mus musculus（小鼠）", table_cell)],
        [Paragraph("数据库类型", table_cell), Paragraph("Gene Ontology - Molecular Function", table_cell)],
        [Paragraph("本地数据库版本", table_cell), Paragraph(database_version, table_cell)],
        [Paragraph("严格 TF 根条目", table_cell), Paragraph("GO:0003700 DNA-binding transcription factor activity", table_cell)],
        [Paragraph("纳入范围", table_cell), Paragraph("GO:0003700 及其全部后代条目；只纳入 MF 注释。", table_cell)],
        [Paragraph("排除范围", table_cell), Paragraph("仅为转录共因子、染色质调节因子或一般转录调节因子的基因，不自动判为 TF。", table_cell)],
        [Paragraph("节点大小", table_cell), Paragraph("DEG 固定为 8；TF 大小 = 42 + 0.9 x 诱导网络度数。", table_cell)],
        [Paragraph("边定义", table_cell), Paragraph("97 个候选基因在原 GeneMANIA 网络中的诱导子图，共 1,181 条证据边。", table_cell)],
    ]
    method_table = Table(method_rows, colWidths=[44 * mm, 154 * mm], repeatRows=1)
    method_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#164E96")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7F9FB")]),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#D9E2EC")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(method_table)
    story.append(Spacer(1, 3 * mm))

    story.append(Paragraph("互作证据构成", h2_style))
    interaction_data = [[Paragraph("互作类型", table_header), Paragraph("边数", table_header), Paragraph("占比", table_header)]]
    for edge_type, count in edge_counts.most_common():
        interaction_data.append([
            Paragraph(edge_type, table_cell),
            Paragraph(f"{count:,}", table_cell),
            Paragraph(f"{count / induced_edge_count:.1%}", table_cell),
        ])
    interaction_table = Table(interaction_data, colWidths=[80 * mm, 32 * mm, 32 * mm], repeatRows=1)
    interaction_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#164E96")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7F9FB")]),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#D9E2EC")),
        ("ALIGN", (1, 1), (-1, -1), "CENTER"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ]))
    story.append(interaction_table)
    story.append(Spacer(1, 3 * mm))

    story.append(Paragraph("解释与限制", h2_style))
    notes = [
        "GO 严格判据强调 DNA 结合型转录因子，因此适合制作可复现的 TF/DEG 网络，但可能漏掉尚未被 GO 完整注释的 TF。",
        "Eya1 通常作为转录共激活因子处理；Nupr1 具有应激和染色质相关调控功能，但在本数据库版本中未满足 GO:0003700 严格条件。",
        "节点大小表示网络连接度，不代表表达倍数、差异显著性或调控效应大小。",
        "GeneMANIA 边为综合证据，不能直接解释为 TF 对靶基因的方向性调控。",
    ]
    for note in notes:
        story.append(Paragraph("- " + note, body_style))
        story.append(Spacer(1, 0.5 * mm))

    story.append(Spacer(1, 0.5 * mm))
    story.append(Paragraph(
        "参考信息：Gene Ontology GO:0003700 - https://amigo.geneontology.org/amigo/term/GO:0003700；"
        "Bioconductor org.Mm.eg.db / GO.db；网络来源为 top100_network.cys 中的 GeneMANIA Mus musculus 会话数据。",
        small_style,
    ))

    doc.build(story, onFirstPage=footer, onLaterPages=footer)


def replace_thumbnail(extracted_root: Path, preview_png: Path) -> None:
    thumbnail = next(extracted_root.glob("*/session_thumbnail.png"))
    shutil.copy2(preview_png, thumbnail)


def repack(extracted_root: Path, output_cys: Path) -> None:
    output_cys.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output_cys, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(extracted_root.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(extracted_root).as_posix())


def main() -> None:
    args = parse_args()
    classification = load_classification(args.classification_csv)
    args.backup_cys.parent.mkdir(parents=True, exist_ok=True)
    if not args.backup_cys.exists():
        shutil.copy2(args.input_cys, args.backup_cys)

    with tempfile.TemporaryDirectory(prefix="top100-go-tf-") as temp_dir:
        extracted = Path(temp_dir)
        with zipfile.ZipFile(args.input_cys) as archive:
            archive.extractall(extracted)
        names = [path.relative_to(extracted).as_posix() for path in extracted.rglob("*") if path.is_file()]

        node_table = extracted / find_entry(
            names,
            "tables/87-M.+musculus+%281%29/SHARED_ATTRS-org.cytoscape.model.CyNode",
        )
        network_xgmml = extracted / find_entry(names, "/networks/87-M.+musculus+%281%29.xgmml")
        edge_table = extracted / find_entry(
            names,
            "tables/87-M.+musculus+%281%29/SHARED_ATTRS-org.cytoscape.model.CyEdge",
        )
        vizmap = extracted / find_entry(names, "/session_vizmap.xml")
        view_xgmml = extracted / find_entry(names, "/views/128-7954-M.+musculus+%281%29.xgmml")

        _, nodes = parse_node_table(node_table)
        graph_nodes, induced_edges, degree, edge_weights = parse_network(network_xgmml, nodes)
        edge_counts = edge_type_counts(edge_table, induced_edges)
        tf_sizes = update_node_table(node_table, classification, degree)

        for needle in (
            "tables/87-M.+musculus+%281%29/SHARED_ATTRS-org.cytoscape.model.CyNetwork",
            "tables/128-M.+musculus+%281%29/LOCAL_ATTRS-org.cytoscape.model.CyNetwork",
        ):
            update_network_name(extracted / find_entry(names, needle), "top100_network_GO_TF")

        style_name = update_visual_style(vizmap, tf_sizes)
        positions = compute_layout(graph_nodes, edge_weights)
        separate_tf_positions(positions, nodes, classification)
        update_view(view_xgmml, positions, style_name)
        render_network(
            args.preview_png,
            graph_nodes,
            edge_weights,
            positions,
            nodes,
            classification,
            degree,
            tf_sizes,
        )
        make_report(
            args.report_pdf,
            args.preview_png,
            classification,
            nodes,
            degree,
            edge_counts,
            tf_sizes,
            graph_nodes,
            len(induced_edges),
        )
        replace_thumbnail(extracted, args.preview_png)
        repack(extracted, args.output_cys)

    print(f"Modified session: {args.output_cys}")
    print(f"Report PDF: {args.report_pdf}")
    print(f"Preview PNG: {args.preview_png}")
    print(f"Query nodes: {len(graph_nodes)}")
    print(f"Induced edges: {len(induced_edges)}")
    print("TFs: " + ", ".join(sorted(tf_sizes)))


if __name__ == "__main__":
    main()
