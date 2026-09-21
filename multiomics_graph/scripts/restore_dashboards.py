"""Restore the original pipeline-generated dashboards from the earlier montage.

dashboards_overview.png was composed from the original outputs dashboards (before
any regeneration), so each original summary_dashboard.png can be cropped back out
of it. B36 is copied from the untouched B36_summary_dashboard.png.

Run from multiomics_graph/:  python scripts/restore_dashboards.py
"""

from pathlib import Path

from PIL import Image

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "outputs"
MONTAGE = BASE / "reports" / "dashboards_overview.png"

W, H = 2388, 1521  # montage cell geometry (PAD=18, LABEL_H=30)
OFFSETS = {
    "B36": (18, 18),
    "MS_14384": (18 + 1 * (W + 18), 18),
    "MS_14385": (18 + 2 * (W + 18), 18),
    "MS_14386": (18, 18 + 1 * (H + 30 + 18)),
    "MS_14387": (18 + 1 * (W + 18), 18 + 1 * (H + 30 + 18)),
}


def main():
    montage = Image.open(MONTAGE).convert("RGB")
    for strain, (x, y) in OFFSETS.items():
        cell = montage.crop((x, y, x + W, y + H))
        bbox = cell.getbbox()
        if bbox is None:
            print(f"  {strain}: empty cell, skipping")
            continue
        img = cell.crop(bbox)
        out = OUT / strain / "figures" / "summary_dashboard.png"
        img.save(out)
        print(f"Restored {strain} -> {out} ({img.width}x{img.height})")


if __name__ == "__main__":
    main()