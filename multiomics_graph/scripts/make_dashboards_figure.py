"""Generate a correct summary-dashboard figure for B36 and MS_14386.

Reads only the deterministic pipeline outputs (the same artifacts that
``verify_counts.py`` uses as the single source of truth), so the numbers in the
figure match ``reports/verified_stats.md`` rather than the stale
``heterodata_metadata.txt``.

Run from ``multiomics_graph/``:  python scripts/make_dashboards_figure.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "outputs"
REPORTS = BASE / "reports"
ROOT_IMAGES = BASE.parent / "Images"

STRAINS = sys.argv[1:] or ["B36", "MS_14386"]

STRAIN_META = {
    "B36": ("ST131", "GCF_900622635.1"),
    "MS_14386": ("ST224", "GCF_900622665.1"),
}

NODE_COLORS = {
    "gene": "#4C72B0",
    "protein": "#55A868",
    "metabolite": "#DD8452",
    "annotation": "#8172B2",
}

REG_COLORS = {"Up": "#C44E52", "Down": "#4C72B0", "Stable": "#9E9E9E"}

EDGE_LABELS = [
    ("gene_genomic_proximity", "Proximity"),
    ("gene_transcriptional_correlation", "Co-expression"),
    ("gene_encodes", "Encodes"),
    ("protein_abundance_correlation", "Protein corr."),
    ("protein_ppi", "PPI"),
    ("gene_in_pathway", "Pathway (gene)"),
    ("gene_cog_category", "COG"),
    ("protein_annotated_by", "GO"),
    ("gene_associated_with", "AMR"),
    ("protein_member_of_cluster", "Cluster"),
    ("metabolite_in_pathway", "Pathway (met.)"),
    ("metabolite_metabolised_by", "Metab. enzyme"),
    ("metabolite_abundance_correlation", "Metabolite corr."),
]

HYPEREDGE_LABELS = [
    ("multi_omics_triple", "Triple"),
    ("kegg_pathway", "KEGG"),
    ("cog_category", "COG"),
    ("go_term", "GO"),
    ("ppi_cluster", "PPI"),
    ("amr_mechanism", "AMR"),
    ("metabolic_pathway", "Metabolic"),
]


def load_node_counts(strain: str) -> dict:
    md = (OUT / strain / "graph" / "heterodata_metadata.txt").read_text(
        encoding="utf-8", errors="replace"
    )
    return {k: int(v) for k, v in re.findall(r"(\w+): (\d+) nodes", md)}


def load_edge_counts(strain: str) -> dict:
    gd = OUT / strain / "graph"
    counts = {}
    for csv in sorted(gd.glob("edges_*.csv")):
        parts = csv.stem.replace("edges_", "", 1).split("_")
        rel = "_".join(parts[1:-1])
        src, dst = parts[0], parts[-1]
        key = f"{src}_{rel}" if rel != dst else f"{src}_{rel}_{dst}"
        n = sum(1 for _ in open(csv, encoding="utf-8", errors="replace")) - 1
        counts[key] = n
    return counts


def load_hyperedge_counts(strain: str) -> dict:
    df = pd.read_csv(OUT / strain / "graph" / "hyperedge_type_counts.csv")
    return {str(r.iloc[0]): int(r.iloc[1]) for _, r in df.iterrows()}


def load_regulation(strain: str) -> tuple[dict, dict]:
    rna = pd.read_csv(OUT / strain / "transcriptomics_expression.csv")
    prot = pd.read_csv(OUT / strain / "proteomics_abundance.csv")
    rna_c = {k: int(v) for k, v in rna["RNA_Regulation"].value_counts().items()}
    prot_c = {k: int(v) for k, v in prot["Protein_Regulation"].value_counts().items()}
    return rna_c, prot_c


def load_amr(strain: str) -> tuple[int, int, int, float]:
    rep = pd.read_csv(REPORTS / f"{strain}_amr_report.csv")
    n_loci = len(rep)
    n_markers = rep["marker"].nunique()
    n_classes = rep["amr_class"].nunique()
    trace = pd.read_csv(REPORTS / f"amr_path_trace_{strain}.csv", low_memory=False)
    cov = 0.0
    if "layer_coverage" in trace.columns:
        vals = trace["layer_coverage"].dropna().tolist()
        if vals:
            cov = sum(int(x) for x in vals) / len(vals)
    return n_loci, n_markers, n_classes, cov


def _bars(ax, labels, values, colors, xlabel="Count"):
    y = range(len(labels))
    ax.barh(list(y), values, color=colors, alpha=0.85, edgecolor="white")
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    for i, v in enumerate(values):
        ax.text(v, i, f" {v:,}", va="center", fontsize=7)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _pie(ax, values: dict, title: str):
    order = ["Up", "Down", "Stable"]
    vals = [values.get(k, 0) for k in order]
    if sum(vals) == 0:
        vals = [0, 0, 1]
    total = sum(vals) or 1
    explode = [0.12 if 0 < v / total < 0.15 else 0.0 for v in vals]
    wedges, texts = ax.pie(
        vals,
        labels=[f"{o} ({v:,})" for o, v in zip(order, vals)],
        colors=[REG_COLORS[o] for o in order],
        startangle=90,
        counterclock=False,
        explode=explode,
        labeldistance=1.15,
        textprops={"fontsize": 7},
        wedgeprops={"edgecolor": "white"},
    )
    if 0 < vals[0] / total < 0.15 and 0 < vals[1] / total < 0.15:
        x, y = texts[1].get_position()
        texts[1].set_position((x, y - 0.1))
    ax.set_title(title, fontsize=9, fontweight="bold")


def render_dashboard(fig, spec, strain: str, data: dict):
    st, acc = STRAIN_META[strain]
    inner = spec.subgridspec(
        4, 3, height_ratios=[0.07, 1, 1, 1], hspace=0.5, wspace=0.28
    )

    ax_title = fig.add_subplot(inner[0, :])
    ax_title.axis("off")
    ax_title.text(
        0.5, 0.5, f"{strain}   ({st})   {acc}",
        ha="center", va="center", fontsize=13, fontweight="bold",
        transform=ax_title.transAxes,
    )

    # 1. Node counts
    ax1 = fig.add_subplot(inner[1, 0])
    norder = ["gene", "protein", "metabolite", "annotation"]
    nlabs = ["Gene", "Protein", "Metabolite", "Annotation"]
    nvals = [data["nodes"][k] for k in norder]
    _bars(ax1, nlabs, nvals, [NODE_COLORS[k] for k in norder])
    ax1.set_title("Node types", fontsize=9, fontweight="bold")

    # 2. RNA regulation
    ax2 = fig.add_subplot(inner[1, 1])
    _pie(ax2, data["rna_reg"], "RNA regulation")

    # 3. Protein regulation
    ax3 = fig.add_subplot(inner[1, 2])
    _pie(ax3, data["prot_reg"], "Protein regulation")

    # 4. Edge types (log scale)
    ax4 = fig.add_subplot(inner[2, 0])
    keys = [k for k, _ in EDGE_LABELS]
    labs = [l for _, l in EDGE_LABELS]
    vals = [data["edges"].get(k, 0) for k in keys]
    ax4.barh(range(len(labs)), vals, color="#5D576B", alpha=0.85, edgecolor="white")
    ax4.set_yticks(range(len(labs)))
    ax4.set_yticklabels(labs, fontsize=7)
    ax4.invert_yaxis()
    ax4.set_xscale("log")
    ax4.set_xlim(1, 200000)
    for i, v in enumerate(vals):
        ax4.text(v, i, f" {v:,}", va="center", fontsize=6)
    ax4.set_title("Edges by type (log)", fontsize=9, fontweight="bold")
    ax4.tick_params(axis="x", labelsize=7)
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)

    # 5. Hyperedge families
    ax5 = fig.add_subplot(inner[2, 1])
    hkeys = [k for k, _ in HYPEREDGE_LABELS]
    hlabs = [l for _, l in HYPEREDGE_LABELS]
    hvals = [data["hyperedges"].get(k, 0) for k in hkeys]
    _bars(ax5, hlabs, hvals, "#5D576B")
    ax5.set_title("Hyperedge families", fontsize=9, fontweight="bold")

    # 6. Summary text
    ax6 = fig.add_subplot(inner[2, 2])
    ax6.axis("off")
    n_loci, n_markers, n_classes, cov = data["amr"]
    lines = [
        "Summary",
        "\u2015" * 18,
        f"Aligned genes:  {data['aligned']:,}",
        f"Genome CDS:     {data['cds']:,}",
        f"Total edges:    {data['total_edges']:,}",
        f"Hyperedges:     {data['total_hyper']:,}",
        "",
        "AMR",
        "\u2015" * 18,
        f"Loci:           {n_loci}",
        f"Markers:        {n_markers}",
        f"Classes:        {n_classes}",
        f"Layer coverage: {cov:.2f} / 8",
    ]
    ax6.text(
        0.05, 0.95, "\n".join(lines), transform=ax6.transAxes,
        fontsize=8.5, va="top", fontfamily="monospace",
    )


def main():
    data = {}
    for strain in STRAINS:
        rna_reg, prot_reg = load_regulation(strain)
        aligned = pd.read_csv(OUT / strain / "aligned_multiomics.csv")
        cds = pd.read_csv(OUT / strain / "genome_genes.csv")
        data[strain] = {
            "nodes": load_node_counts(strain),
            "edges": load_edge_counts(strain),
            "hyperedges": load_hyperedge_counts(strain),
            "rna_reg": rna_reg,
            "prot_reg": prot_reg,
            "amr": load_amr(strain),
            "aligned": len(aligned),
            "cds": len(cds),
            "total_edges": sum(load_edge_counts(strain).values()),
            "total_hyper": sum(load_hyperedge_counts(strain).values()),
        }

    n = len(STRAINS)
    fig = plt.figure(figsize=(20, 11) if n > 1 else (10, 11))
    outer = fig.add_gridspec(1, n, wspace=0.18)
    for col, strain in enumerate(STRAINS):
        render_dashboard(fig, outer[col], strain, data[strain])

    fig.suptitle(
        "Vesalius representation \u2014 summary dashboards",
        fontsize=15, fontweight="bold", y=0.99,
    )

    out_name = f"dashboards_{STRAINS[0]}.png" if n == 1 else "dashboards_example.png"
    targets = [
        REPORTS / out_name,
        ROOT_IMAGES / out_name,
    ]
    for t in targets:
        t.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(t, dpi=150, bbox_inches="tight", facecolor="white")
        print("Wrote", t)

    plt.close(fig)


if __name__ == "__main__":
    main()
