"""Figure of representative hyperedges, one per family, from the saved B36 hypergraph.

Reads outputs/B36/graph/hyperedges.csv and hypergraph_nodes.csv, resolves readable
names for gene/protein/metabolite nodes from the saved alignment tables, and draws
each selected hyperedge as a hub-and-spoke diagram (the hyperedge hub in the centre,
member nodes arranged on a circle).

Run from multiomics_graph/:  python scripts/hyperedge_samples_figure.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "outputs" / "B36"
REPORTS = BASE / "reports"

FAMILY_ORDER = [
    "multi_omics_triple",
    "kegg_pathway",
    "cog_category",
    "go_term",
    "ppi_cluster",
    "amr_mechanism",
    "metabolic_pathway",
]
FAMILY_COLORS = {
    "multi_omics_triple": "#1f77b4",
    "kegg_pathway": "#ff7f0e",
    "cog_category": "#2ca02c",
    "go_term": "#d62728",
    "ppi_cluster": "#9467bd",
    "amr_mechanism": "#8c564b",
    "metabolic_pathway": "#17becf",
}
NODE_TYPE_COLORS = {
    "genome": "#4C72B0",
    "transcriptome": "#DD8452",
    "protein": "#55A868",
    "metabolite": "#C44E52",
    "annotation": "#8172B2",
}

MAX_SHOWN = 12


def load_name_maps():
    gene = {}
    prot = {}
    if (OUT / "genome_genes.csv").exists():
        g = pd.read_csv(OUT / "genome_genes.csv", dtype=str)
        for _, r in g.iterrows():
            gid = r.get("GeneID", "")
            name = r.get("GeneName", "") or gid
            gene[gid] = name
            pid = r.get("ProteinID", "")
            if pid:
                prot[pid] = name
    metab = {}
    if (OUT / "metabolomics_abundance.csv").exists():
        m = pd.read_csv(OUT / "metabolomics_abundance.csv", dtype=str)
        for _, r in m.iterrows():
            mid = r.get("MetaboliteID", "")
            name = r.get("MetaboliteName", "") or mid
            if mid:
                metab[mid] = name
    return gene, prot, metab


def node_type(node, node_types):
    if node in node_types:
        return node_types[node]
    if node.startswith("WP_"):
        return "protein"
    if node.endswith("_RNA"):
        return "transcriptome"
    if node.startswith("CHEBI:"):
        return "metabolite"
    if node.startswith("AMR:"):
        return "annotation"
    if node.startswith("EW"):
        return "genome"
    return "annotation"


def readable(node, node_types, gene, prot, metab):
    if node in gene:
        return gene[node]
    if node in prot:
        return prot[node]
    if node in metab:
        return metab[node]
    if node.endswith("_RNA"):
        base = node[:-4]
        return f"{readable(base, node_types, gene, prot, metab)} RNA"
    return node


def pick(hyperedges, family, preferred=None):
    sub = hyperedges[hyperedges["type"] == family]
    if sub.empty:
        return None
    if preferred:
        sel = sub[sub["name"] == preferred]
        if not sel.empty:
            return sel.iloc[0]
    sub = sub[sub["n_nodes"] >= 3].sort_values("n_nodes")
    if sub.empty:
        sub = hyperedges[hyperedges["type"] == family].sort_values("n_nodes")
    return sub.iloc[0]


def draw_panel(ax, nodes, hub_name, family, node_types, gene, prot, metab, title):
    nodes = list(dict.fromkeys(nodes))
    truncated = len(nodes) - MAX_SHOWN
    shown = nodes[:MAX_SHOWN]
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")

    n = len(shown)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False) if n else []
    for i, nd in enumerate(shown):
        a = angles[i]
        x, y = np.cos(a), np.sin(a)
        ax.plot([0, x], [0, y], color="#cccccc", lw=1, zorder=1)
        t = node_type(nd, node_types)
        color = NODE_TYPE_COLORS.get(t, "#999999")
        ax.scatter([x], [y], s=120, color=color, edgecolors="white", lw=0.6, zorder=3)
        label = str(readable(nd, node_types, gene, prot, metab))
        if label.lower() in ("nan", "none", ""):
            label = nd
        if len(label) > 16:
            label = label[:15] + "~"
        ax.text(x * 1.18, y * 1.18, label, ha="center", va="center",
                fontsize=6.5, zorder=4)

    ax.scatter([0], [0], s=320, color=FAMILY_COLORS[family], edgecolors="white",
               lw=1.2, zorder=3)
    hub_label = hub_name if len(hub_name) <= 24 else hub_name[:23] + "~"
    ax.text(0, 0, hub_label, ha="center", va="center", color="white",
            fontsize=6.5, fontweight="bold", zorder=5)
    if truncated > 0:
        ax.text(1.15, -1.35, f"+{truncated} more", ha="right", va="center",
                fontsize=6, color="#666666")
    ax.set_title(title, fontsize=9, pad=2)


def main():
    hyperedges = pd.read_csv(OUT / "graph" / "hyperedges.csv")
    node_df = pd.read_csv(OUT / "graph" / "hypergraph_nodes.csv", dtype=str)
    node_types = dict(zip(node_df["node_id"], node_df["node_type"]))
    gene, prot, metab = load_name_maps()

    preferred = {
        "multi_omics_triple": "Triple_EW036_RS26140",
        "amr_mechanism": "beta_lactam_resistance",
    }

    sel = {}
    for fam in FAMILY_ORDER:
        h = pick(hyperedges, fam, preferred.get(fam))
        if h is None:
            continue
        sel[fam] = h

    n_rows = int(np.ceil(len(sel) / 3))
    fig, axes = plt.subplots(n_rows, 3, figsize=(12, 4.0 * n_rows))
    axes = np.atleast_2d(axes)
    for k, fam in enumerate(FAMILY_ORDER):
        r, c = divmod(k, 3)
        if fam not in sel:
            axes[r, c].axis("off")
            continue
        h = sel[fam]
        draw_panel(axes[r, c], h["nodes"].split(";"), h["name"], fam, node_types,
                   gene, prot, metab, fam.replace("_", " ").title())

    for idx in range(len(sel), len(axes.ravel())):
        axes.ravel()[idx].axis("off")

    handles = [
        plt.Line2D([], [], marker="o", ls="", color=col, label=name)
        for name, col in NODE_TYPE_COLORS.items()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=len(NODE_TYPE_COLORS),
               frameon=False, fontsize=8)
    fig.suptitle("Representative hyperedges in the B36 hypergraph "
                 "(hub = hyperedge, spokes = member nodes)",
                 fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])

    out = REPORTS / "hyperedge_samples.png"
    fig.savefig(out, dpi=200)
    print("Wrote", out)
    for fam, h in sel.items():
        print(f"  {fam:20s} n={h['n_nodes']:3d}  {h['name']}")


if __name__ == "__main__":
    main()