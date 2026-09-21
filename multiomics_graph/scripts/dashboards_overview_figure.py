"""Compose the per-strain summary dashboards into two overview figures.

Figure 1: B36, MS_14384, MS_14385
Figure 2: MS_14386, MS_14387

Run from multiomics_graph/:  python scripts/dashboards_overview_figure.py
"""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "outputs"
REPORTS = BASE / "reports"

GROUPS = [
    ("dashboards_example.png", ["B36", "MS_14386"]),
    ("dashboards_overview_1.png", ["B36", "MS_14384", "MS_14385"]),
    ("dashboards_overview_2.png", ["MS_14386", "MS_14387"]),
]
COLS = 3
PAD = 18


def compose(out_name, strains):
    imgs = []
    for s in strains:
        f = OUT / s / "figures" / "summary_dashboard.png"
        if not f.exists():
            f = OUT / s / "figures" / f"{s}_summary_dashboard.png"
        imgs.append((s, Image.open(f).convert("RGB")))

    w = max(im.width for _, im in imgs)
    h = max(im.height for _, im in imgs)
    rows = (len(imgs) + COLS - 1) // COLS
    canvas_w = COLS * w + (COLS + 1) * PAD
    canvas_h = rows * h + (rows + 1) * PAD

    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")

    for i, (s, im) in enumerate(imgs):
        r, c = divmod(i, COLS)
        x = PAD + c * (w + PAD)
        y = PAD + r * (h + PAD)
        canvas.paste(im, (x, y))

    out = REPORTS / out_name
    canvas.save(out)
    print("Wrote", out, f"({canvas_w}x{canvas_h})")


def main():
    for out_name, strains in GROUPS:
        compose(out_name, strains)


if __name__ == "__main__":
    main()