"""
scripts/amr_report.py
Phase 3 - biological interpretation of the AMR knowledge layer.

Reads the already-generated pipeline outputs for each strain and produces:
  1. A per-strain AMR report (marker presence, locus, COG, KEGG, PPI
     neighborhood, omics availability).
  2. A cross-strain comparison matrix (marker x strain presence, shared vs
     strain-specific mechanisms) + a comparison figure.

This is read-only over `multiomics_graph/outputs/<strain>/` and the AMR
manifest. It does NOT re-run the pipeline.

Usage (from multiomics_graph/):
    python scripts/amr_report.py [--output reports]
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from amr.amr import AMRManifest

REPO = Path(__file__).resolve().parent.parent
STRAINS_OUT = REPO / 'outputs'
MANIFEST = AMRManifest()

STRAIN_ORDER = ['B36', 'MS_14384', 'MS_14386', 'MS_14387']

CLASS_SHORT = {
    'beta_lactam': 'Bla',
    'aminoglycoside': 'AAC/APH',
    'sulfonamide': 'Sul',
    'trimethoprim': 'DfrA',
    'tetracycline': 'Tet',
    'macrolide': 'Mph',
    'fluoroquinolone': 'FQ',
    'phenicol': 'FloR',
    'fosfomycin': 'FosA',
}


def _read(path: Path):
    return pd.read_csv(path, low_memory=False) if path.exists() else None


def load_genome_loci(strain: Path) -> set:
    gg = _read(strain / 'genome_genes.csv')
    if gg is None:
        gg = _read(strain / 'gene_features.csv')
    if gg is None or 'GeneID' not in gg.columns:
        return set()
    return set(gg['GeneID'])


def build_strain_report(name: str):
    sd = STRAINS_OUT / name
    if not sd.exists():
        return None

    genome_loci = load_genome_loci(sd)
    aligned = _read(sd / 'aligned_multiomics.csv')
    aligned_loci = set(aligned['GeneID']) if aligned is not None else set()

    # COG categories per locus (gene -> annotation)
    cog_of = defaultdict(set)
    ce = _read(sd / 'graph' / 'edges_gene_cog_category_annotation.csv')
    if ce is not None:
        for _, r in ce.iterrows():
            cog_of.setdefault(r['gene_id'], set()).add(r['annotation_id'])

    # KEGG pathways per locus
    kegg_of = defaultdict(set)
    ke = _read(sd / 'graph' / 'edges_gene_in_pathway_annotation.csv')
    if ke is None:
        ke = _read(sd / 'edges_gene_in_pathway_annotation.csv')
    if ke is not None and len(ke.columns) >= 2:
        cols = list(ke.columns)
        gcol = next((c for c in cols if 'gene' in c), None)
        acol = next((c for c in cols if c not in ('relation',) and c != gcol), None)
        for _, r in ke.iterrows():
            if gcol and acol:
                kegg_of.setdefault(r[gcol], set()).add(r[acol])

    # PPI neighbor count per protein
    n_ppi_of = {}
    pf = _read(sd / 'graph' / 'edges_protein_ppi_protein.csv')
    pdf = None
    if pf is not None and len(pf.columns) >= 2:
        pcols = list(pf.columns)
        pdf = pf
        ps, pd_ = pcols[0], pcols[1]

    markers = []
    for m in MANIFEST.markers(name):
        prot = m['protein_ids'][0] if m['protein_ids'] else None
        n_ppi = 0
        if prot and pdf is not None:
            n_ppi = int((pdf[ps] == prot).sum() + (pdf[pd_] == prot).sum())
        for locus in m['locus_tags']:
            markers.append({
                'strain': name,
                'marker': m['name'],
                'amr_class': m['amr_class'],
                'locus_tag': locus,
                'protein_id': ';'.join(m['protein_ids']) or '(pseudogene)',
                'COG': ','.join(sorted(cog_of.get(locus, []))) or '-',
                'n_kegg_pathways': len(kegg_of.get(locus, [])),
                'n_ppi_neighbors': n_ppi,
                'in_genome': locus in genome_loci,
                'in_aligned_omics': locus in aligned_loci,
            })

    return {
        'name': name,
        'markers': markers,
        'genome_only': [m['marker'] for m in markers if not m['in_aligned_omics']],
        'quantified': [m['marker'] for m in markers if m['in_aligned_omics']],
    }


def write_report(rec: dict, out_dir: Path):
    df = pd.DataFrame(rec['markers'])
    if df.empty:
        return
    df.to_csv(out_dir / f"{rec['name']}_amr_report.csv", index=False)


def build_comparison(records):
    all_markers = []
    for r in records:
        for m in r['markers']:
            if m['marker'] not in all_markers:
                all_markers.append(m['marker'])
    heat = pd.DataFrame(
        0, index=all_markers, columns=[r['name'] for r in records]
    )
    for r in records:
        for m in r['markers']:
            heat.loc[m['marker'], r['name']] = 1
    return heat


def plot_heatmap(heat: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(7, max(4, 0.35 * len(heat))))
    im = ax.imshow(heat.values.astype(int), cmap='Reds', aspect='auto')
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=45, ha='right')
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index, fontsize=8)
    ax.set_title('AMR determinant presence across strains', fontsize=11)
    for i in range(len(heat.index)):
        for j in range(len(heat.columns)):
            ax.text(j, i, int(heat.iloc[i, j]),
                    ha='center', va='center', color='black', fontsize=7)
    fig.colorbar(im, ax=ax, label='present')
    plt.tight_layout()
    path = out / 'amr_cross_strain_heatmap.png'
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [Saved] {path}")


def symmetries_summary(records):
    marker_present = defaultdict(set)
    class_present = defaultdict(set)
    all_sets = {r['name'] for r in records}
    for r in records:
        for m in r['markers']:
            marker_present[m['marker']].add(r['name'])
            class_present[m['amr_class']].add(r['name'])

    core_mechanisms = sorted(k for k, v in class_present.items()
                             if v == all_sets)
    shared_markers = sorted(k for k, v in marker_present.items() if len(v) >= 2)
    strain_specific = sorted(k for k, v in marker_present.items() if len(v) == 1)
    only_one = [k for k, v in marker_present.items() if len(v) == 1]
    return {
        'core_mechanisms': core_mechanisms,
        'shared_markers': shared_markers,
        'strain_specific_markers': strain_specific,
        'singletons': only_one,
    }


def main():
    ap = argparse.ArgumentParser(description='AMR knowledge-layer report')
    ap.add_argument('--output', default='reports')
    ap.add_argument('--strains', nargs='*', default=STRAIN_ORDER)
    args = ap.parse_args()

    out = REPO / args.output
    out.mkdir(parents=True, exist_ok=True)

    records = []
    for name in args.strains:
        rec = build_strain_report(name)
        if rec is None:
            print(f"  [{name}] no outputs dir, skipping")
            continue
        records.append(rec)
        write_report(rec, out)
        print(f"  [{name}] {len(rec['markers'])} markers | "
              f"{len(rec['genome_only'])} genome-only | "
              f"{len(rec['quantified'])} in multi-omics")

    if len(records) >= 2:
        heat = build_comparison(records)
        heat.to_csv(out / 'amr_cross_strain_matrix.csv')
        plot_heatmap(heat, out)

        syn = symmetries_summary(records)
        txt = (
            "AMR synergies across strains\n"
            "===========================\n"
            f"Core mechanisms (present in every strain):\n"
            f"  {syn['core_mechanisms']}\n\n"
            f"Shared markers (>=2 strains):\n"
            f"  {syn['shared_markers']}\n\n"
            f"Strain-specific markers (single strain):\n"
            f"  {syn['strain_specific_markers']}\n"
        )
        (out / 'amr_synergies.txt').write_text(txt, encoding='utf-8')
        print("  cross-strain matrix, heatmap and synergies written")

    print("\nReports directory:", out.resolve())


if __name__ == '__main__':
    main()