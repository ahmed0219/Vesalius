"""Generate the comparative and discovery figures used in the results chapter.

All figures read only the verified saved artifacts and write to
``reports/``. Run from ``multiomics_graph/``:

    python scripts/results_figures.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "outputs"
REPORTS = BASE / "reports"
STRAINS = ["B36", "MS_14384", "MS_14385", "MS_14386", "MS_14387"]
AMR_STRAINS = ["B36", "MS_14384", "MS_14386", "MS_14387"]

GENE_REL = {
    "genomic_proximity": "genomic proximity",
    "encodes": "encodes",
    "in_pathway": "gene-pathway",
    "cog_category": "gene-COG",
    "associated_with": "gene-AMR",
}
PROT_REL = {
    "ppi": "PPI",
    "annotated_by": "protein-GO",
    "member_of_cluster": "cluster",
}
MET_REL = {
    "in_pathway": "metab-pathway",
    "metabolised_by": "metab-enzyme",
}


def edge_table(strain):
    gd = OUT / strain / "graph"
    out = {}
    for csv in gd.glob("edges_*.csv"):
        parts = csv.stem.replace("edges_", "", 1).split("_")
        rel = "_".join(parts[1:-1])
        src, dst = parts[0], parts[-1]
        n = sum(1 for _ in open(csv, encoding="utf-8", errors="replace")) - 1
        out[(src, rel, dst)] = n
    return out


def hyperedge_table(strain):
    f = OUT / strain / "graph" / "hyperedge_type_counts.csv"
    out = {}
    if f.exists():
        df = pd.read_csv(f)
        for _, r in df.iterrows():
            out[str(r.iloc[0])] = int(r.iloc[1])
    return out


def fig1_cross_strain_edge_comparison():
    """Cross-strain comparative bar chart of realized edge-type counts."""
    rows = []
    for s in STRAINS:
        et = edge_table(s)
        for (src, rel, dst), n in et.items():
            if rel in GENE_REL:
                label = GENE_REL[rel]
            elif rel in PROT_REL:
                label = PROT_REL[rel]
            elif rel in MET_REL:
                label = MET_REL[rel]
            else:
                continue
            rows.append({"strain": s, "edge_type": label, "count": n})
    df = pd.DataFrame(rows)
    # aggregate (a relation may appear under multiple (src,rel,dst) keys that
    # map to the same display label, e.g. gene-pathway)
    df = df.groupby(["strain", "edge_type"], as_index=False)["count"].sum()
    order = ["genomic proximity", "encodes", "gene-pathway", "gene-COG",
             "gene-AMR", "PPI", "protein-GO", "cluster", "metab-pathway",
             "metab-enzyme"]
    order = [o for o in order if o in set(df["edge_type"])]
    df["edge_type"] = pd.Categorical(df["edge_type"], categories=order, ordered=True)
    piv = df.pivot(index="edge_type", columns="strain", values="count").fillna(0)
    piv = piv.loc[order]

    fig, ax = plt.subplots(figsize=(11, 5))
    x = np.arange(len(piv))
    width = 0.15
    colors = plt.cm.tab10(np.linspace(0, 1, len(STRAINS)))
    for i, s in enumerate(STRAINS):
        ax.bar(x + (i - 2) * width, piv[s].values, width, label=s,
               color=colors[i])
    ax.set_yscale("symlog")
    ax.set_xticks(x)
    ax.set_xticklabels(piv.index, rotation=30, ha="right")
    ax.set_ylabel("edge count (log scale)")
    ax.set_title("Realized edge-type counts across the five strain graphs")
    ax.legend(ncol=5, fontsize=8)
    fig.tight_layout()
    out = REPORTS / "results_edge_type_comparison.png"
    fig.savefig(out, dpi=150)
    print("Wrote", out)


def fig2_layer_coverage():
    """Mean layer coverage of AMR determinants per strain."""
    vals, labels = [], []
    for s in AMR_STRAINS:
        f = REPORTS / f"amr_path_trace_{s}.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f, low_memory=False)
        if "layer_coverage" in df.columns:
            vals.append(df["layer_coverage"].dropna().astype(float).mean())
            labels.append(s)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(labels, vals, color="#2a9d8f")
    ax.axhline(8, color="grey", ls="--", lw=0.8)
    ax.set_ylabel("mean layer coverage / 8")
    ax.set_ylim(0, 9)
    ax.set_title("Mean AMR evidence-layer coverage per strain")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.1, f"{v:.2f}", ha="center", fontsize=9)
    fig.tight_layout()
    out = REPORTS / "results_layer_coverage.png"
    fig.savefig(out, dpi=150)
    print("Wrote", out)


def fig3_amr_support_heatmap():
    """AMR determinant x strain: multi-omics-supported vs genome-only."""
    rows = []
    for s in AMR_STRAINS:
        f = REPORTS / f"{s}_amr_report.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f)
        for _, r in df.iterrows():
            rows.append({
                "strain": s,
                "marker": str(r["marker"]),
                "support": "multi-omics" if r["in_aligned_omics"] else "genome-only",
            })
    df = pd.DataFrame(rows)
    markers = sorted(df["marker"].unique())
    mat = np.zeros((len(markers), len(AMR_STRAINS)))
    for i, m in enumerate(markers):
        for j, s in enumerate(AMR_STRAINS):
            sub = df[(df["marker"] == m) & (df["strain"] == s)]
            if len(sub):
                mat[i, j] = 1 if (sub["support"] == "multi-omics").any() else 0.5
    fig, ax = plt.subplots(figsize=(8, 8))
    cmap = matplotlib.colors.ListedColormap(["#ffffff", "#f4a261", "#2a9d8f"])
    im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(AMR_STRAINS)))
    ax.set_xticklabels(AMR_STRAINS)
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels(markers, fontsize=8)
    ax.set_title("AMR determinant multi-omics support per strain",
                 fontsize=11)
    cbar = fig.colorbar(im, ticks=[0, 0.5, 1], ax=ax)
    cbar.ax.set_yticklabels(["absent", "genome-only", "multi-omics"])
    fig.tight_layout()
    fig.subplots_adjust(top=0.94)
    out = REPORTS / "results_amr_support_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("Wrote", out)


def fig4_concordance():
    """RNA vs protein logFC concordance for multi-omics AMR determinants."""
    pts = []
    for s in AMR_STRAINS:
        f = REPORTS / f"amr_path_trace_{s}.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f, low_memory=False)
        if not {"L2_RNA_log2_derived", "L3_Protein_logFC"}.issubset(df.columns):
            continue
        sub = df[df["L3_Protein_logFC"].notna()].copy()
        for _, r in sub.iterrows():
            pts.append({
                "strain": s,
                "marker": str(r["marker"]),
                "rna": float(r["L2_RNA_log2_derived"]),
                "protein": float(r["L3_Protein_logFC"]),
            })
    if not pts:
        return
    df = pd.DataFrame(pts)
    fig, ax = plt.subplots(figsize=(7, 6))
    for s in AMR_STRAINS:
        sub = df[df["strain"] == s]
        ax.scatter(sub["rna"], sub["protein"], label=s, s=60)
        for _, r in sub.iterrows():
            ax.annotate(r["marker"], (r["rna"], r["protein"]), fontsize=7,
                        xytext=(4, 4), textcoords="offset points")
    lim = max(abs(df["rna"].min()), abs(df["rna"].max()),
              abs(df["protein"].min()), abs(df["protein"].max())) * 1.1
    ax.plot([-lim, lim], [-lim, lim], ls="--", color="grey", lw=0.8)
    ax.axhline(0, color="grey", lw=0.5)
    ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlabel("RNA log2(Sera/RPMI)")
    ax.set_ylabel("Protein logFC")
    ax.set_title("RNA-protein concordance of multi-omics AMR determinants")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = REPORTS / "results_rna_protein_concordance.png"
    fig.savefig(out, dpi=150)
    print("Wrote", out)


def fig5_candidate_network():
    """Gene candidates recurring in AMR molecular context.

    A tiered selection spanning the determinant-count distribution: the
    top-ranked gene candidates plus the amino-acid/aromatic-biosynthesis and
    central-metabolism genes discussed in the results chapter.
    """
    cand = REPORTS / "amr_path_candidates.csv"
    if not cand.exists():
        return
    cdf = pd.read_csv(cand)
    gene = cdf[cdf["entity_type"] == "gene"]
    if gene.empty:
        return

    wanted = ["ilvB", "mrcA", "mdh", "accB", "nanA", "gltD",
              "asd", "aroB", "pabA", "argD", "pckA", "gpmM", "ilvN"]
    by_name = {}
    for _, r in gene.iterrows():
        by_name.setdefault(str(r["candidate_name"]), r)
    rows, seen = [], set()
    for nm in wanted:
        if nm in by_name and nm not in seen:
            rows.append(by_name[nm])
            seen.add(nm)
    for _, r in gene.iterrows():
        nm = str(r["candidate_name"])
        if nm not in seen:
            rows.append(r)
            seen.add(nm)
        if len(rows) >= 12:
            break
    top = pd.DataFrame(rows).head(12)

    fig, ax = plt.subplots(figsize=(8, 5))
    names = [str(r["candidate_name"]) for _, r in top.iterrows()]
    n_det = top["number_of_amr_determinants"].astype(float).values
    y = np.arange(len(names))
    bars = ax.barh(y, n_det, color="#457b9d")
    for yy, v in zip(y, n_det):
        ax.text(v + 0.1, yy, str(int(v)), va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, n_det.max() + 1.5)
    ax.set_xlabel("Number of connected AMR determinants")
    ax.set_title("Gene candidates recurring in AMR molecular context")
    fig.tight_layout()
    out = REPORTS / "results_candidate_genes.png"
    fig.savefig(out, dpi=150)
    print("Wrote", out)


def fig6_hyperedge_comparison():
    """Cross-strain hyperedge family comparison (grouped bar)."""
    families = ["multi_omics_triple", "kegg_pathway", "cog_category",
                "go_term", "ppi_cluster", "amr_mechanism", "metabolic_pathway"]
    vals = {}
    for s in STRAINS:
        vals[s] = hyperedge_table(s)
    present = [f for f in families if any(vals[s].get(f, 0) for s in STRAINS)]
    fig, ax = plt.subplots(figsize=(11, 5))
    x = np.arange(len(present))
    width = 0.15
    colors = plt.cm.tab10(np.linspace(0, 1, len(STRAINS)))
    for i, s in enumerate(STRAINS):
        ax.bar(x + (i - 2) * width, [vals[s].get(f, 0) for f in present],
               width, label=s, color=colors[i])
    ax.set_yscale("symlog")
    ax.set_xticks(x)
    ax.set_xticklabels(present, rotation=30, ha="right")
    ax.set_ylabel("hyperedge count (log scale)")
    ax.set_title("Hyperedge family counts across the five strain graphs")
    ax.legend(ncol=5, fontsize=8)
    fig.tight_layout()
    out = REPORTS / "results_hyperedge_comparison.png"
    fig.savefig(out, dpi=150)
    print("Wrote", out)


def main():
    fig1_cross_strain_edge_comparison()
    fig2_layer_coverage()
    fig3_amr_support_heatmap()
    fig4_concordance()
    fig5_candidate_network()
    fig6_hyperedge_comparison()


if __name__ == "__main__":
    main()