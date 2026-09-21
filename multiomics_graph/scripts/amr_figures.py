"""
scripts/amr_figures.py
Regenerate the AMR highlight figures from saved pipeline outputs.

Every known AMR determinant from the paper (AMR manifest) is kept as a gene
node even when it has no RNA/protein quantitation. Genome-only markers are
drawn with an outlined/hatched style, multi-omics-supported markers as solid
nodes; the AMR mechanism/class hyperedge nodes are shown as diamonds. The
script then automatically verifies that 100% of the manifest's AMR loci are
present in the generated figure's node set (raising on any absence).

Read-only over `multiomics_graph/outputs/<strain>/` and the AMR manifest; it
does NOT re-run the pipeline.

Usage (from multiomics_graph/):
    python scripts/amr_figures.py
    python scripts/amr_figures.py --strains B36 MS_14386
    python scripts/amr_figures.py --output outputs --visualizer-output outputs
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from amr.amr import AMRManifest
from visualization.network_plots import NetworkVisualizer

REPO = Path(__file__).resolve().parent.parent
STRAINS_OUT = REPO / 'outputs'
MANIFEST = AMRManifest()

STRAIN_ORDER = ['B36', 'MS_14384', 'MS_14386', 'MS_14387']

COL_TYPE_MAP = {
    'gene': 'gene', 'protein': 'protein', 'annotation': 'annotation',
    'metabolite': 'metabolite',
    'gene_target': 'gene', 'protein_target': 'protein',
}


def load_draw_graph(sd: Path):
    """Rebuild the pairwise draw graph exactly as main.py does for figures."""
    import networkx as nx
    node_types = {'gene': set(), 'protein': set(),
                  'annotation': set(), 'metabolite': set()}
    edge_types = {}
    graph_dir = sd / 'graph'
    edge_csvs = sorted(graph_dir.glob('edges_*.csv'))
    for csv_path in edge_csvs:
        df = pd.read_csv(csv_path, low_memory=False)
        cols = [c for c in df.columns if c != 'relation']
        if len(cols) < 2:
            continue
        src_col, dst_col = cols[0], cols[1]
        parts = csv_path.stem.replace('edges_', '', 1).split('_')
        rel = '_'.join(parts[1:-1]) if len(parts) >= 3 else parts[0]

        def get_type(c):
            for prefix, t in COL_TYPE_MAP.items():
                if str(c).startswith(prefix):
                    return t
            return 'gene'

        stype, dtype = get_type(src_col), get_type(dst_col)
        key = (stype, rel, dtype)
        edge_types[key] = []
        for _, row in df.iterrows():
            s, d = str(row[src_col]).strip(), str(row[dst_col]).strip()
            if s and d:
                edge_types[key].append((f"{stype}_{s}", f"{dtype}_{d}"))
                node_types[stype].add(f"{stype}_{s}")
                node_types[dtype].add(f"{dtype}_{d}")
    return {k: sorted(v) for k, v in node_types.items()}, edge_types


def aligned_loci(sd: Path) -> set:
    al = sd / 'aligned_multiomics.csv'
    if not al.exists():
        return set()
    df = pd.read_csv(al, low_memory=False)
    return set(str(g) for g in df['GeneID'].dropna().tolist())


def build_amr_draw(name: str, sd: Path):
    node_types, edge_types = load_draw_graph(sd)
    aligned = aligned_loci(sd)

    locus_to_amr = MANIFEST.locus_to_amr(name)
    amr_node_types = {k: set(v) for k, v in node_types.items()}
    amr_edge_types = {k: list(v) for k, v in edge_types.items()}

    amr_class_map = {}
    genome_only_loci = set()
    for locus in locus_to_amr:
        gnode = f"gene_{locus}"
        amr_class_map[gnode] = locus_to_amr[locus]['amr_class']
        amr_node_types['gene'].add(gnode)
        if locus not in aligned:
            genome_only_loci.add(gnode)

    mechanism_nodes = []
    amr_mech_key = ('gene', 'amr_mechanism', 'amr_mechanism')
    amr_edge_types.setdefault(amr_mech_key, [])
    amr_node_types.setdefault('amr_mechanism', set())
    for locus, rec in locus_to_amr.items():
        mech = f"AMR:{rec['amr_class']}"
        mechanism_nodes.append(mech)
        amr_node_types['amr_mechanism'].add(mech)
        amr_edge_types[amr_mech_key].append((f"gene_{locus}", mech))

    amr_node_type_lists = {k: sorted(v) for k, v in amr_node_types.items()}
    return {
        'node_types': amr_node_type_lists,
        'edge_types': amr_edge_types,
        'amr_loci': list(locus_to_amr),
        'amr_class_map': amr_class_map,
        'genome_only_loci': genome_only_loci,
        'mechanism_nodes': mechanism_nodes,
    }


def main():
    ap = argparse.ArgumentParser(description='Regenerate AMR highlight figures')
    ap.add_argument('--strains', nargs='*', default=STRAIN_ORDER)
    ap.add_argument('--output', default='outputs',
                    help='Base output directory (default: outputs)')
    ap.add_argument('--no-verify', action='store_true',
                    help='Skip the 100%-coverage assertion')
    args = ap.parse_args()

    out_base = REPO / args.output
    ok = True
    for name in args.strains:
        sd = out_base / name
        if not sd.exists():
            print(f"  [{name}] no outputs dir, skipping")
            continue
        n_markers = len(MANIFEST.markers(name))
        n_loci = len(MANIFEST.amr_loci(name))
        if n_markers == 0:
            print(f"  [{name}] no AMR markers in manifest (0) — nothing to plot")
            continue

        draw = build_amr_draw(name, sd)
        viz = NetworkVisualizer(output_dir=sd / 'figures')

        amr_state = viz.plot_amr_highlight(
            draw['node_types'], draw['edge_types'],
            draw['amr_loci'], draw['amr_class_map'],
            title=f"{name} — AMR Determinants in the Graph",
            filename="amr_highlight.png",
            max_nodes=0,
            genome_only_loci=draw['genome_only_loci'],
            mechanism_nodes=draw['mechanism_nodes'],
            verify=not args.no_verify,
        )
        if amr_state is None:
            continue

        n_total = len(amr_state)
        n_multi = sum(1 for s in amr_state.values() if s == 'multi_omics')
        n_genome = sum(1 for s in amr_state.values() if s == 'genome_only')
        status = "PASS" if n_total == n_loci else "FAIL"
        if n_total != n_loci:
            ok = False
        print(f"  [{name}] {status}  {n_total}/{n_loci} AMR loci drawn "
              f"({n_multi} multi-omics, {n_genome} genome-only)")

    if not ok:
        print("\nERROR: some strains did not reach 100% AMR locus coverage.")
        sys.exit(1)
    print("\nAll AMR loci from the paper present in the figures (100%).")


if __name__ == '__main__':
    main()