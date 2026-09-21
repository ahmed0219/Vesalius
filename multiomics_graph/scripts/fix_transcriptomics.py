"""Correct the transcriptomics-derived outputs for a strain.

The pipeline originally loaded the raw-counts file and computed RNA_logFC as a
raw-count difference, which produced nonsensical regulation calls (e.g. only
~1% "Stable") and co-expression edges on raw counts. This script recomputes:

  * RNA_RPMI / RNA_Sera / RNA_logFC / RNA_Regulation from the log2 CPM matrix
  * the gene-gene co-expression edges (|Pearson r| >= 0.7 + mutual top-k) on the
    log2 CPM replicate profiles

and overwrites ``transcriptomics_expression.csv``, the RNA columns of
``aligned_multiomics.csv``, and the co-expression edge CSVs.

Run from ``multiomics_graph/``:  python scripts/fix_transcriptomics.py [strain ...]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent          # multiomics_graph/
OUT = BASE / "outputs"
STRAINS_DIR = BASE.parent / "strains"                  # repo-root strains/

CORR_THRESHOLD = 0.7
CORR_TOP_K = 10
REGULATION_THRESHOLD = 1.0


def find_cpmlog2(strain: str) -> Path:
    d = STRAINS_DIR / strain / "transcriptomic"
    matches = sorted(d.rglob("*cpmlog2*")) if d.exists() else []
    if not matches:
        raise FileNotFoundError(f"no *cpmlog2* file for {strain}")
    return matches[0]


def load_config(strain: str) -> dict:
    with open(STRAINS_DIR / strain / "config.json") as f:
        return json.load(f)


def load_cpmlog2_matrix(strain: str) -> tuple[pd.DataFrame, list[str], list[str]]:
    cfg = load_config(strain)
    rpmi = cfg["transcriptomics"]["rpmi_samples"]
    sera = cfg["transcriptomics"]["sera_samples"]

    path = find_cpmlog2(strain)
    df = pd.read_csv(path, sep="\t")
    df.columns = [str(c).strip('"') for c in df.columns]

    # Column names may be sample IDs directly or embedded in a longer name.
    def resolve(ids):
        out = []
        for sid in ids:
            hit = [c for c in df.columns if str(sid) in str(c)]
            if not hit:
                raise KeyError(f"sample {sid} not found in {path.name}")
            out.append(hit[0])
        return out

    rpmi_cols = resolve(rpmi)
    sera_cols = resolve(sera)
    for c in rpmi_cols + sera_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df, rpmi_cols, sera_cols


def correlation_edges(corr: np.ndarray, threshold: float, top_k: int):
    """Replicate HeterogeneousGraphBuilder._correlation_edges."""
    abs_c = np.abs(corr)
    n = corr.shape[0]
    mask = (abs_c >= threshold) & np.triu(np.ones((n, n), dtype=bool), k=1)
    if top_k and top_k > 0:
        order = np.argsort(-abs_c, axis=1)
        mutual = np.zeros((n, n), dtype=bool)
        for i in range(n):
            for j in order[i, :top_k]:
                if j != i:
                    mutual[i, j] = True
                    mutual[j, i] = True
        mask &= mutual
    i_idx, j_idx = np.where(mask)
    values = [float(corr[i, j]) for i, j in zip(i_idx, j_idx)]
    return i_idx, j_idx, values


def process(strain: str) -> dict:
    out = OUT / strain

    df, rpmi_cols, sera_cols = load_cpmlog2_matrix(strain)
    gene_col = "GeneID" if "GeneID" in df.columns else df.columns[0]
    name_col = "GeneName" if "GeneName" in df.columns else None

    expr = pd.DataFrame()
    expr["GeneID"] = df[gene_col]
    if name_col:
        expr["GeneName"] = df[name_col]
    else:
        expr["GeneName"] = df[gene_col]
    expr["RNA_RPMI"] = df[rpmi_cols].mean(axis=1)
    expr["RNA_Sera"] = df[sera_cols].mean(axis=1)
    expr["RNA_logFC"] = expr["RNA_Sera"] - expr["RNA_RPMI"]
    expr["RNA_Regulation"] = np.select(
        [expr["RNA_logFC"] > REGULATION_THRESHOLD,
         expr["RNA_logFC"] < -REGULATION_THRESHOLD],
        ["Up", "Down"], default="Stable",
    )
    expr.to_csv(out / "transcriptomics_expression.csv", index=False)

    # Update the RNA columns of the aligned table.
    aligned_path = out / "aligned_multiomics.csv"
    aligned = pd.read_csv(aligned_path)
    rna_map = expr.set_index("GeneID")[["RNA_RPMI", "RNA_Sera", "RNA_logFC"]]
    aligned = aligned.drop(columns=["RNA_RPMI", "RNA_Sera", "RNA_logFC"])
    aligned = aligned.merge(rna_map, left_on="GeneID", right_index=True, how="left")
    aligned.to_csv(aligned_path, index=False)

    # Recompute co-expression edges from the log2 CPM replicate matrix.
    genome = pd.read_csv(out / "genome_genes.csv")
    gene_ids = sorted(set(genome["GeneID"].astype(str)) | set(aligned["GeneID"].astype(str)))
    rep_cols = rpmi_cols + sera_cols
    rna_idx = df.set_index(gene_col)[rep_cols]
    gene_set = [g for g in gene_ids if g in rna_idx.index]
    rep_mat = rna_idx.loc[gene_set].values.astype(float)
    col_mean = np.nanmean(rep_mat, axis=0)
    inds = np.where(np.isnan(rep_mat))
    rep_mat[inds] = np.take(col_mean, inds[1])

    corr = np.corrcoef(rep_mat)
    i_idx, j_idx, _ = correlation_edges(corr, CORR_THRESHOLD, CORR_TOP_K)

    rows = pd.DataFrame({
        "gene_id": [gene_set[i] for i in i_idx],
        "gene_target": [gene_set[j] for j in j_idx],
        "relation": "transcriptional_correlation",
    })
    rows.to_csv(out / "graph" / "edges_gene_transcriptional_correlation_gene.csv",
                index=False)

    reg = expr["RNA_Regulation"].value_counts().to_dict()
    return {
        "strain": strain,
        "n_genes": len(expr),
        "regulation": {k: int(reg.get(k, 0)) for k in ("Up", "Down", "Stable")},
        "n_coexpression": len(rows),
    }


def main():
    strains = sys.argv[1:] or ["B36", "MS_14386"]
    for s in strains:
        r = process(s)
        print(f"{r['strain']}: genes={r['n_genes']} "
              f"Up={r['regulation']['Up']} Down={r['regulation']['Down']} "
              f"Stable={r['regulation']['Stable']} "
              f"coexpression={r['n_coexpression']}")


if __name__ == "__main__":
    main()
