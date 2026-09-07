#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export the rendered gene interaction network as a one-page PDF."
    )
    parser.add_argument("input_png", type=Path)
    parser.add_argument("output_pdf", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)

    page_width, page_height = landscape(A4)
    margin = 6.0
    available_width = page_width - 2 * margin
    available_height = page_height - 2 * margin

    image = ImageReader(str(args.input_png))
    image_width, image_height = image.getSize()
    scale = min(available_width / image_width, available_height / image_height)
    draw_width = image_width * scale
    draw_height = image_height * scale
    x = (page_width - draw_width) / 2
    y = (page_height - draw_height) / 2

    pdf = canvas.Canvas(str(args.output_pdf), pagesize=(page_width, page_height))
    pdf.setTitle("Regulatory networks of top TFs")
    pdf.setAuthor("Codex")
    pdf.setSubject("GO Molecular Function TF gene interaction network")
    pdf.setFillColorRGB(1, 1, 1)
    pdf.rect(0, 0, page_width, page_height, stroke=0, fill=1)
    pdf.drawImage(
        image,
        x,
        y,
        width=draw_width,
        height=draw_height,
        preserveAspectRatio=True,
        mask="auto",
    )
    pdf.showPage()
    pdf.save()

    print(f"Network PDF: {args.output_pdf}")


if __name__ == "__main__":
    main()
