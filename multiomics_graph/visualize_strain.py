"""
visualize_strain.py
Re-generate plots for a completed strain run from saved outputs.

Usage:
    python visualize_strain.py B36
    python visualize_strain.py MS_14384
"""

import sys
import json
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')

# Ensure project root is in path
PROJECT_ROOT = Path(__file__).parent.absolute()
sys.path.insert(0, str(PROJECT_ROOT))
from main import load_config
from graph.hypergraph import BiologicalHypergraphBuilder
from visualization.network_plots import NetworkVisualizer
from annotation.eggnog import COG_CATEGORIES


def main():
    parser = argparse.ArgumentParser(description="Visualize a completed strain run")
    parser.add_argument('strain', type=str, help='Strain name (e.g. B36, MS_14384)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory for figures (default: outputs/<strain>/figures)')
    parser.add_argument('--show', action='store_true',
                        help='Display plots interactively instead of saving to files')
    args = parser.parse_args()

    config = load_config(args.strain)
    base_dir = Path('outputs') / args.strain

    if args.output_dir:
        fig_dir = Path(args.output_dir)
    else:
        fig_dir = base_dir / 'figures'

    viz = NetworkVisualizer(str(fig_dir))

    # 1. Load aligned table
    aligned_path = base_dir / 'aligned_multiomics.csv'
    if aligned_path.exists():
        aligned = pd.read_csv(aligned_path)
        viz.plot_regulation_concordance(aligned)
        print(f"  Regulation concordance → {fig_dir / 'regulation_concordance.png'}")
    else:
        aligned = None
        print("  [skip] aligned_multiomics.csv not found")

    # 2. Load COG annotations
    cog_path = base_dir / 'gene_features.csv'
    if cog_path.exists():
        features = pd.read_csv(cog_path)
        cog_cols = [c for c in features.columns if c.startswith('COG_')]
        if cog_cols:
            cog_counts = {}
            for col in cog_cols:
                count = features[col].sum()
                if count > 0:
                    cog_code = col.replace('COG_', '')
                    cog_counts[cog_code] = int(count)
            viz.plot_cog_distribution(cog_counts, cog_categories=COG_CATEGORIES)
            print(f"  COG distribution → {fig_dir / 'cog_distribution.png'}")

    # 3. Load hypergraph
    hyperedge_path = base_dir / 'graph' / 'hyperedges.csv'
    if hyperedge_path.exists():
        he_df = pd.read_csv(hyperedge_path)
        hyperedge_type_map = {}
        for he_type in he_df['hyperedge_type'].unique() if 'hyperedge_type' in he_df.columns else []:
            ids = he_df[he_df['hyperedge_type'] == he_type]['hyperedge_id'].unique()
            hyperedge_type_map[he_type] = list(ids)
        if hyperedge_type_map:
            he_counts = {k: len(v) for k, v in hyperedge_type_map.items()}
            viz.plot_hyperedge_type_distribution(he_counts)
            print(f"  Hyperedge types → {fig_dir / 'hyperedge_types.png'}")

        # Load incidence matrix
        inc_path = base_dir / 'graph' / 'incidence_matrix.npy'
        if inc_path.exists():
            H = np.load(inc_path)
            viz.plot_incidence_matrix(H, hyperedge_types=hyperedge_type_map)
            print(f"  Incidence matrix → {fig_dir / 'incidence_matrix.png'}")

    # 4. Summary dashboard
    genome_s = {'genes': None}
    rna_s = {'genes': None}
    protein_s = {'proteins': None}
    graph_s = {}
    hg_s = {'n_hyperedges': None}

    if aligned is not None:
        graph_s['gene_nodes'] = len(aligned)
        graph_s['protein_nodes'] = aligned.get('ProteinID', pd.Series()).nunique() if 'ProteinID' in aligned.columns else 0

    try:
        he_df = pd.read_csv(base_dir / 'graph' / 'hyperedges.csv')
        hg_s['n_hyperedges'] = len(he_df)
        if 'n_nodes' in he_df.columns:
            hg_s['mean_hyperedge_size'] = float(he_df['n_nodes'].mean())
        nodes = pd.read_csv(base_dir / 'graph' / 'hypergraph_nodes.csv')
        hg_s['n_nodes'] = int(nodes['node_id'].nunique())
        H = np.load(base_dir / 'graph' / 'incidence_matrix.npy')
        hg_s['incidence_density'] = float(np.mean(H > 0))
    except Exception:
        pass

    try:
        genome_genes = pd.read_csv(base_dir / 'genome_genes.csv')
        genome_s['genes'] = len(genome_genes)
    except Exception:
        pass
    try:
        rna = pd.read_csv(base_dir / 'transcriptomics_expression.csv')
        rna_s['genes'] = len(rna)
        if 'RNA_Regulation' in rna.columns:
            vc = rna['RNA_Regulation'].value_counts()
            rna_s['up_regulated'] = int(vc.get('Up', 0))
            rna_s['down_regulated'] = int(vc.get('Down', 0))
            rna_s['stable'] = int(vc.get('Stable', 0))
    except Exception:
        pass
    try:
        prot = pd.read_csv(base_dir / 'proteomics_abundance.csv')
        protein_s['proteins'] = len(prot)
        if 'Protein_Regulation' in prot.columns:
            vc = prot['Protein_Regulation'].value_counts()
            protein_s['up_regulated'] = int(vc.get('Up', 0))
            protein_s['down_regulated'] = int(vc.get('Down', 0))
            protein_s['stable'] = int(vc.get('Stable', 0))
    except Exception:
        pass
    try:
        with open(base_dir / 'graph' / 'hyperedge_type_counts.csv') as f:
            counts = pd.read_csv(f)
            if 'count' in counts.columns:
                hg_s['n_hyperedges'] = counts['count'].sum()
                type_col = counts.columns[0]
                hg_s['hyperedge_types'] = dict(
                    zip(counts[type_col], counts['count'])
                )
    except Exception:
        pass

    viz.plot_summary_dashboard(genome_s, rna_s, protein_s, graph_s, hg_s,
                               title=args.strain)
    print(f"  Summary dashboard → {fig_dir / 'summary_dashboard.png'}")

    print(f"\n  All figures saved to: {fig_dir.resolve()}")


if __name__ == '__main__':
    main()
