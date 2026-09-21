"""Topology analysis: are AMR determinants structurally distinguishable in the
multi-omics graph?

Reads the saved (enriched) graph artifacts and asks whether curated AMR
determinants occupy different positions in the graph than the non-AMR
background, on a set of network features:

  - knowledge_degree:  annotation edges on the gene node
                       (in_pathway + cog_category + associated_with)
  - protein_go_degree: GO edges on the encoded protein
  - protein_ppi_degree: STRING PPI edges on the encoded protein
  - proximity_degree:  genomic-proximity edges on the gene node
  - hyperedge_participation: number of hyperedges containing the gene or its
                       encoded protein (all families)
  - rna_measured / protein_measured / in_aligned: multi-omics coverage

Analysis:
  1. Descriptive per-strain comparison (AMR vs non-AMR median).
  2. Pooled stratified permutation test (AMR labels permuted within strain)
     per feature (two-sided, 10000 draws).
  3. Genome-only vs multi-omics-supported AMR determinants (within the AMR set).

Outputs:
  - reports/amr_topology_summary.md
  - reports/amr_topology_boxplot.png   (AMR vs background, key features)
  - reports/amr_topology_hyperedge_families.png (hyperedge-family participation)

Run from ``multiomics_graph/``:  python scripts/amr_topology_analysis.py
"""

import math
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
RNG = np.random.default_rng(42)

GENE_EDGE_RELS = {
    "genomic_proximity": "proximity_degree",
    "in_pathway": "knowledge_degree",
    "cog_category": "knowledge_degree",
    "associated_with": "knowledge_degree",
}
PROTEIN_EDGE_RELS = {
    "annotated_by": "protein_go_degree",
    "ppi": "protein_ppi_degree",
}


def load_edges(strain):
    gd = OUT / strain / "graph"
    gene = {}
    protein = {}
    if gd.exists():
        for csv in gd.glob("edges_*.csv"):
            parts = csv.stem.replace("edges_", "", 1).split("_")
            rel = "_".join(parts[1:-1])
            src, dst = parts[0], parts[-1]
            df = pd.read_csv(csv)
            cols = [c for c in df.columns if c != "relation"]
            if len(cols) < 2:
                continue
            s, d = cols[0], cols[1]
            for r, f in GENE_EDGE_RELS.items():
                if rel == r and src == "gene":
                    for x in df[s].dropna():
                        gene.setdefault(str(x), {})[f] = gene.setdefault(str(x), {}).get(f, 0) + 1
                    break
            for r, f in PROTEIN_EDGE_RELS.items():
                if rel == r and src == "protein":
                    for x in df[s].dropna():
                        protein.setdefault(str(x), {})[f] = protein.setdefault(str(x), {}).get(f, 0) + 1
                    break
            if rel == "encodes" and src == "gene":
                enc = dict(zip(df[s].astype(str), df[d].astype(str)))
    # encodes map
    enc = {}
    f = gd / "edges_gene_encodes_protein.csv"
    if f.exists():
        df = pd.read_csv(f)
        cols = [c for c in df.columns if c != "relation"]
        if len(cols) >= 2:
            enc = dict(zip(df[cols[0]].astype(str), df[cols[1]].astype(str)))
    return gene, protein, enc


def load_hyperedge_membership(strain):
    """{node_id: set of hyperedge families} for gene and protein nodes.

    Hyperedge ``nodes`` may carry a type prefix (``gene_<locus>``,
    ``<locus>_RNA``, ``WP_...``) or be a bare locus tag; we index both the
    raw node id and its stripped locus-tag form so AMR loci resolve correctly
    across all families.
    """
    f = OUT / strain / "graph" / "hyperedges.csv"
    out = {}
    if not f.exists():
        return out
    df = pd.read_csv(f)
    for _, r in df.iterrows():
        fam = str(r["type"])
        for node in str(r["nodes"]).split(";"):
            node = node.strip()
            if not node:
                continue
            keys = {node}
            if node.startswith("gene_"):
                keys.add(node[len("gene_"):])
            elif node.endswith("_RNA"):
                keys.add(node[: -len("_RNA")])
            for k in keys:
                d = out.setdefault(k, set())
                d.add(fam)
    return out


def is_amr_map(strain):
    f = OUT / strain / "gene_features.csv"
    amr = {}
    if f.exists():
        df = pd.read_csv(f)
        if "is_amr_gene" in df.columns and "GeneID" in df.columns:
            for _, r in df.iterrows():
                if bool(r["is_amr_gene"]):
                    amr[str(r["GeneID"])] = True
    return amr


def aligned_set(strain):
    f = OUT / strain / "aligned_multiomics.csv"
    if not f.exists():
        return set()
    df = pd.read_csv(f, low_memory=False)
    return set(str(x) for x in df["GeneID"].dropna())


def build_features(strain):
    gene_edges, protein_edges, enc = load_edges(strain)
    he = load_hyperedge_membership(strain)
    amr = is_amr_map(strain)
    aligned = aligned_set(strain)

    rows = []
    # gene universe from genome_genes (gene nodes)
    gf = OUT / strain / "genome_genes.csv"
    universe = set(gene_edges)
    if gf.exists():
        df = pd.read_csv(gf)
        if "GeneID" in df.columns:
            universe |= set(str(x) for x in df["GeneID"].dropna())

    for g in universe:
        ge = gene_edges.get(g, {})
        prot = enc.get(g)
        pe = protein_edges.get(prot, {}) if prot else {}
        he_count = len(he.get(g, set())) + (len(he.get(prot, set())) if prot else 0)
        row = {
            "strain": strain,
            "gene": g,
            "amr": bool(amr.get(g, False)),
            "in_aligned": g in aligned,
            "knowledge_degree": ge.get("knowledge_degree", 0),
            "proximity_degree": ge.get("proximity_degree", 0),
            "protein_go_degree": pe.get("protein_go_degree", 0),
            "protein_ppi_degree": pe.get("protein_ppi_degree", 0),
            "hyperedge_participation": he_count,
        }
        rows.append(row)
    return pd.DataFrame(rows)


def pooled_perm_test(df, feature, n_perm=10000, seed=42):
    """Stratified permutation test: permute feature values within strain.

    Streaming over permutations (memory-safe): each draw permutes feature
    values independently within each strain while AMR labels stay fixed, then
    compares the pooled AMR-vs-background mean difference to the observed one.
    """
    rng = np.random.default_rng(seed)
    strains = df["strain"].values
    y = df[feature].values.astype(float)
    amr_flag = df["amr"].values.astype(bool)

    groups = [np.where(strains == s)[0] for s in np.unique(strains)]
    obs_diff = y[amr_flag].mean() - y[~amr_flag].mean()

    y_perm = y.copy()
    perm_diff = np.empty(n_perm)
    for k in range(n_perm):
        for g in groups:
            y_perm[g] = rng.permutation(y_perm[g])
        perm_diff[k] = y_perm[amr_flag].mean() - y_perm[~amr_flag].mean()
    bigger = int((np.abs(perm_diff) >= abs(obs_diff)).sum())
    p = (bigger + 1) / (n_perm + 1)
    return obs_diff, p


def _bh_fdr_correct(pvalues):
    """Benjamini-Hochberg FDR correction across features.

    ``pvalues`` is a dict {feature: p}; returns {feature: q}. Monotonicity is
    enforced in the standard direction (cumulative min propagated from the
    LARGEST p down to the smallest), so q is always >= p.
    """
    feats = list(pvalues)
    n = len(feats)
    if n == 0:
        return {}
    order = sorted(feats, key=lambda f: pvalues[f])
    raw = {f: pvalues[f] * n / (r + 1) for r, f in enumerate(order)}
    q = {}
    prev = float('inf')
    for f in reversed(order):
        v = min(raw[f], prev)
        q[f] = min(1.0, v)
        prev = v
    return q


def _cohens_d(df, feature):
    """Effect size: AMR-vs-background mean difference in pooled-SD units."""
    amr = df.loc[df["amr"], feature].values.astype(float)
    bg = df.loc[~df["amr"], feature].values.astype(float)
    na, nb = len(amr), len(bg)
    if na == 0 or nb == 0:
        return float("nan")
    va, vb = amr.var(ddof=1), bg.var(ddof=1)
    sp = math.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    if sp == 0:
        return 0.0
    return (amr.mean() - bg.mean()) / sp


def _bootstrap_ci_mean_diff(df, feature, n_boot=2000, seed=42, ci=0.95):
    """Bootstrap CI for the AMR-vs-background mean difference, stratified by
    strain within each resample (so the panel structure is preserved).

    Chunked & vectorized: resample indices are drawn per strain in chunks of
    ``CHUNK`` draws (avoiding a giant (n_boot, n_rows) allocation), and pooled
    AMR/bg means are accumulated across chunks.

    Returns (low, high) as 100*ci% percentile interval.
    """
    CHUNK = 250
    rng = np.random.default_rng(seed)
    amr_sum = np.zeros(n_boot)
    amr_n = np.zeros(n_boot)
    bg_sum = np.zeros(n_boot)
    bg_n = np.zeros(n_boot)
    for s in np.unique(df["strain"].values):
        sub = df[df["strain"] == s]
        vals = sub[feature].values.astype(float)
        amr = sub["amr"].values.astype(bool)
        n_s = len(vals)
        if n_s == 0:
            continue
        for start in range(0, n_boot, CHUNK):
            stop = min(start + CHUNK, n_boot)
            m = stop - start
            idx = rng.integers(0, n_s, size=(m, n_s))
            v = vals[idx]              # (m, n_s)
            a = amr[idx]
            amr_sum[start:stop] += (v * a).sum(axis=1)
            amr_n[start:stop] += a.sum(axis=1)
            bg_sum[start:stop] += (v * (~a)).sum(axis=1)
            bg_n[start:stop] += (~a).sum(axis=1)
    valid = (amr_n > 0) & (bg_n > 0)
    diffs = np.full(n_boot, np.nan)
    diffs[valid] = (amr_sum[valid] / amr_n[valid]) - (bg_sum[valid] / bg_n[valid])
    diffs = diffs[~np.isnan(diffs)]
    if len(diffs) == 0:
        return float("nan"), float("nan")
    lo = float(np.percentile(diffs, (1 - ci) / 2 * 100))
    hi = float(np.percentile(diffs, (1 + ci) / 2 * 100))
    return lo, hi


def main():
    all_rows = []
    for s in STRAINS:
        df = build_features(s)
        all_rows.append(df)
        # descriptive
        if df["amr"].any():
            amr_df = df[df["amr"]]
            bg_df = df[~df["amr"]]
            print(f"[{s}] AMR n={len(amr_df)}, background n={len(bg_df)}")
            for feat in ["knowledge_degree", "proximity_degree",
                         "protein_go_degree", "protein_ppi_degree",
                         "hyperedge_participation"]:
                print(f"    {feat}: AMR med={amr_df[feat].median():.1f} "
                      f"bg med={bg_df[feat].median():.1f}")
    pool = pd.concat(all_rows, ignore_index=True)

    features = ["knowledge_degree", "proximity_degree", "protein_go_degree",
                "protein_ppi_degree", "hyperedge_participation"]
    lines = ["# AMR Topology Analysis", ""]
    lines.append(
        "Stratified permutation test (AMR labels permuted within strain, "
        "10,000 draws, two-sided). AMR determinants = curated manifest loci "
        "present in the gene-node universe of each strain."
    )
    lines.append("")
    lines.append("| feature | AMR mean | background mean | diff | p (perm) | "
                 "p_adj (BH) | Cohen's d | CI95 low | CI95 high | n_AMR | n_bg |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    perm_results = {}
    effect_sizes = {}
    cis = {}
    diffs = {}
    for feat in features:
        diff, p = pooled_perm_test(pool, feat)
        perm_results[feat] = p
        diffs[feat] = diff
        effect_sizes[feat] = _cohens_d(pool, feat)
        cis[feat] = _bootstrap_ci_mean_diff(pool, feat)
    p_adj = _bh_fdr_correct(perm_results)
    amr_pool = pool[pool["amr"]]
    bg_pool = pool[~pool["amr"]]
    for feat in features:
        lo, hi = cis[feat]
        lines.append(
            f"| {feat} | {amr_pool[feat].mean():.3f} | {bg_pool[feat].mean():.3f} "
            f"| {diffs[feat]:+.3f} | {perm_results[feat]:.4f} | {p_adj[feat]:.4f} "
            f"| {effect_sizes[feat]:+.3f} | {lo:+.3f} | {hi:+.3f} "
            f"| {len(amr_pool)} | {len(bg_pool)} |"
        )
    lines.append("")
    lines.append("`p (perm)` = two-sided stratified permutation test "
                 "(AMR labels permuted within strain, 10,000 draws). "
                 "`p_adj (BH)` = Benjamini-Hochberg FDR across the 5 features. "
                 "`Cohen's d` = standardized mean difference (effect size). "
                 "`CI95` = bootstrap 95% percentile interval (2,000 resamples, "
                 "stratified by strain).")
    lines.append("")

    # genome-only vs multi-omics within AMR set
    lines.append("## Genome-only vs multi-omics-supported AMR determinants")
    lines.append("")
    lines.append("| feature | genome-only mean | multi-omics mean | n_go | n_mo |")
    lines.append("|---|---|---|---|---|")
    for s in AMR_STRAINS:
        rep = REPORTS / f"{s}_amr_report.csv"
        if not rep.exists():
            continue
        rdf = pd.read_csv(rep)
        go_loci = set(str(x) for x in rdf.loc[rdf["in_aligned_omics"] == False, "locus_tag"])
        sub = pool[(pool["strain"] == s) & (pool["amr"])].copy()
        sub["go"] = sub["gene"].isin(go_loci)
        if len(sub) and sub["go"].any() and (~sub["go"]).any():
            lines.append(f"### {s}")
            for feat in features:
                g = sub[sub["go"]][feat].mean()
                m = sub[~sub["go"]][feat].mean()
                lines.append(
                    f"| {feat} | {g:.2f} | {m:.2f} | "
                    f"{int(sub['go'].sum())} | {int((~sub['go']).sum())} |"
                )
    lines.append("")

    # hyperedge-family participation of AMR determinants
    lines.append("## Hyperedge-family participation of AMR determinants")
    lines.append("")
    fam_counts = {}
    for s in AMR_STRAINS:
        amr = is_amr_map(s)
        enc = load_edges(s)[2]
        he = load_hyperedge_membership(s)
        for g in amr:
            fams = set(he.get(g, set()))
            prot = enc.get(g)
            if prot:
                fams |= set(he.get(prot, set()))
            for f_ in fams:
                fam_counts[f_] = fam_counts.get(f_, 0) + 1
    lines.append("| family | determinant-memberships |")
    lines.append("|---|---|")
    for f_, c in sorted(fam_counts.items(), key=lambda kv: -kv[1]):
        lines.append(f"| {f_} | {c} |")
    lines.append("")

    out = REPORTS / "amr_topology_summary.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out}")

    # Figures
    make_boxplot(pool, features, perm_results, p_adj)
    make_boxplot_stacked(pool, features, perm_results, p_adj)
    make_hyperedge_figure(fam_counts)


def make_boxplot(pool, features, perm_results, p_adj):
    fig, axes = plt.subplots(1, len(features), figsize=(16, 4), sharey=False)
    for ax, feat in zip(axes, features):
        amr_vals = pool.loc[pool["amr"], feat].values
        bg_vals = pool.loc[~pool["amr"], feat].values
        bp = ax.boxplot([bg_vals, amr_vals], tick_labels=["background", "AMR"],
                        widths=0.6, patch_artist=True,
                        medianprops=dict(color="black"))
        for patch, c in zip(bp["boxes"], ["#c9d6e3", "#f4a261"]):
            patch.set_facecolor(c)
        ax.set_title(f"{feat}\np={perm_results[feat]:.4f}\nFDR={p_adj[feat]:.4f}",
                     fontsize=9)
        ax.set_yscale("symlog")
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("AMR determinants vs background: multi-omics graph topology "
                 "(BH-corrected permutation p-values)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = REPORTS / "amr_topology_boxplot.png"
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


def make_boxplot_stacked(pool, features, perm_results, p_adj):
    fig, axes = plt.subplots(len(features), 1, figsize=(7, 13), sharex=True)
    for ax, feat in zip(axes, features):
        amr_vals = pool.loc[pool["amr"], feat].values
        bg_vals = pool.loc[~pool["amr"], feat].values
        bp = ax.boxplot([bg_vals, amr_vals], tick_labels=["background", "AMR"],
                        widths=0.6, patch_artist=True, vert=False,
                        medianprops=dict(color="black"))
        for patch, c in zip(bp["boxes"], ["#c9d6e3", "#f4a261"]):
            patch.set_facecolor(c)
        ax.set_title(f"{feat}\np={perm_results[feat]:.4f}\nFDR={p_adj[feat]:.4f}",
                     fontsize=9)
        ax.set_xscale("symlog")
        ax.grid(axis="x", alpha=0.3)
    fig.suptitle("AMR determinants vs background: multi-omics graph topology "
                 "(BH-corrected permutation p-values)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = REPORTS / "amr_topology_boxplot_stacked.png"
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


def make_hyperedge_figure(fam_counts):
    if not fam_counts:
        return
    names = sorted(fam_counts, key=lambda k: -fam_counts[k])
    vals = [fam_counts[k] for k in names]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(names, vals, color="#457b9d")
    ax.set_xlabel("hyperedge family")
    ax.set_ylabel("AMR determinant-memberships (all strains)")
    ax.set_title("Higher-order context of AMR determinants")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    out = REPORTS / "amr_topology_hyperedge_families.png"
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()