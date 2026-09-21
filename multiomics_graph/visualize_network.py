"""
visualize_network.py
Render the full heterogeneous graph for a completed strain run.

Usage:
    python visualize_network.py B36
    python visualize_network.py MS_14384 --max-nodes 500
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

PROJECT_ROOT = Path(__file__).parent.absolute()
sys.path.insert(0, str(PROJECT_ROOT))
from main import load_config

# Colors
COLORS = {
    'gene': '#2E86AB',
    'protein': '#A23B72',
    'annotation': '#F18F01',
    'genomic_proximity': '#73AB84',
    'encodes': '#C44536',
    'ppi': '#5D576B',
    'in_pathway': '#FF6B6B',
    'cog_category': '#4ECDC4',
}
EDGE_STYLES = {
    'genomic_proximity': 'solid',
    'encodes': 'dotted',
    'ppi': 'dashdot',
    'in_pathway': 'dashed',
    'cog_category': 'dashed',
}


def prefixed_id(node_id: str, ntype: str) -> str:
    """Create a type-prefixed node ID to avoid collisions across types."""
    return f"{ntype}_{node_id}"


def load_strain_data(strain_dir: Path):
    """
    Load all edge CSVs and build node/edge structures.

    Node IDs are prefixed with their type (e.g. gene_0, protein_0)
    to disambiguate when different node types share the same numeric IDs.
    """
    node_types = {'gene': set(), 'protein': set(), 'annotation': set()}
    edge_types = {}

    # Column prefix > node type mapping
    COL_PREFIX_MAP = {
        'gene': 'gene',
        'protein': 'protein',
        'annotation': 'annotation',
    }

    for csv_path in sorted(strain_dir.glob('edges_*.csv')):
        df = pd.read_csv(csv_path)
        cols = df.columns.tolist()

        # Determine edge type from filename: edges_<src>_<relation>_<dst>.csv
        name = csv_path.stem.replace('edges_', '', 1)
        parts = name.split('_')
        if len(parts) >= 3:
            edge_type = '_'.join(parts[1:-1])
        else:
            edge_type = name

        # Find source and target columns (first two data columns, exclude 'relation')
        data_cols = [c for c in cols if c != 'relation']
        if len(data_cols) < 2:
            continue
        src_col, dst_col = data_cols[0], data_cols[1]

        # Infer node types from column name prefixes
        def get_type(col_name):
            for prefix, ntype in COL_PREFIX_MAP.items():
                if col_name.startswith(prefix):
                    return ntype
            return 'gene'

        src_type = get_type(src_col)
        dst_type = get_type(dst_col)

        # Collect nodes with type-prefixed IDs
        for col in data_cols:
            ctype = get_type(col)
            node_types[ctype].update(
                prefixed_id(str(v).strip(), ctype)
                for v in df[col].dropna().unique()
                if pd.notna(v) and str(v).strip()
            )

        # Collect edges with prefixed IDs
        edges = []
        for _, row in df.iterrows():
            src = row[src_col]
            dst = row[dst_col]
            if pd.notna(src) and pd.notna(dst):
                s = prefixed_id(str(src).strip(), src_type)
                d = prefixed_id(str(dst).strip(), dst_type)
                if s and d:
                    edges.append((s, d))
        edge_types[edge_type] = edges

    return node_types, edge_types


def build_graph(node_types: dict, edge_types: dict, max_nodes: int = None):
    """Build a NetworkX graph from node/edge structures."""
    G = nx.Graph()

    if max_nodes:
        # Degree-based sampling: keep highest-degree nodes for connectivity
        degree_count = {}
        for rel, edges in edge_types.items():
            for s, d in edges:
                degree_count[s] = degree_count.get(s, 0) + 1
                degree_count[d] = degree_count.get(d, 0) + 1

        # Sort all nodes by degree, take top max_nodes
        sorted_nodes = sorted(degree_count.items(), key=lambda x: -x[1])
        kept = set(n for n, _ in sorted_nodes[:max_nodes])

        # Ensure at least one node per type
        for ntype, nodes in node_types.items():
            if not any(n in kept for n in nodes):
                # Add the highest-degree node of this missing type
                for n, _ in sorted_nodes:
                    if n in nodes:
                        kept.add(n)
                        break

        for n in kept:
            G.add_node(n, type=next(
                (t for t, ns in node_types.items() if n in ns), 'unknown'
            ))
        for rel, edges in edge_types.items():
            kept_edges = [(s, d) for s, d in edges if s in kept and d in kept]
            if kept_edges:
                G.add_edges_from(kept_edges, relation=rel)
    else:
        # Full graph
        for ntype, nodes in node_types.items():
            for n in nodes:
                G.add_node(n, type=ntype)
        for rel, edges in edge_types.items():
            if edges:
                G.add_edges_from(edges, relation=rel)

    return G


def plot_graph(G, node_types, edge_types, output_path: Path, title: str):
    """Render the graph at high resolution."""
    if G.number_of_nodes() == 0:
        print("  Empty graph, nothing to plot")
        return

    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    print(f"  Computing layout (may take a minute)...")

    # Layout: spring with tuned parameters for large graphs
    k_val = 1.5 / (G.number_of_nodes() ** 0.3)
    pos = nx.spring_layout(G, seed=42, k=k_val, iterations=50)

    # Map node colors by type
    node_colors = []
    for n in G.nodes():
        ntype = G.nodes[n].get('type', 'gene')
        node_colors.append(COLORS.get(ntype, '#888888'))

    # Edge colors and styles by relation
    edge_colors = []
    edge_styles = []
    for u, v, d in G.edges(data=True):
        rel = d.get('relation', 'ppi')
        edge_colors.append(COLORS.get(rel, '#888888'))
        edge_styles.append(EDGE_STYLES.get(rel, 'solid'))

    fig, ax = plt.subplots(1, 1, figsize=(32, 24))

    # Draw edges
    for rel in set(d.get('relation', '') for _, _, d in G.edges(data=True)):
        edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('relation') == rel]
        if edges:
            nx.draw_networkx_edges(
                G, pos, edgelist=edges, ax=ax,
                edge_color=COLORS.get(rel, '#888888'),
                style=EDGE_STYLES.get(rel, 'solid'),
                alpha=0.15, width=0.5,
            )

    # Draw nodes (no labels for full graph — too many)
    for ntype, color in COLORS.items():
        if ntype in node_types and node_types[ntype]:
            nodes_of_type = [n for n in G.nodes() if G.nodes[n].get('type') == ntype]
            if nodes_of_type:
                nx.draw_networkx_nodes(
                    G, pos, nodelist=nodes_of_type, ax=ax,
                    node_size=8, node_color=color, alpha=0.7,
                    edgecolors='none',
                )

    # Legend
    legend_elements = []
    for ntype, color in COLORS.items():
        if ntype in node_types and node_types[ntype]:
            legend_elements.append(
                mpatches.Patch(color=color, label=ntype, alpha=0.8)
            )
    # Add edge type legend
    for rel, color in COLORS.items():
        if rel in edge_types and edge_types[rel]:
            legend_elements.append(
                plt.Line2D([0], [0], color=color, linestyle=EDGE_STYLES.get(rel, 'solid'),
                           label=rel, linewidth=1.5)
            )

    ax.legend(handles=legend_elements, fontsize=10, loc='upper right',
              title=f'{G.number_of_nodes()} nodes, {G.number_of_edges()} edges',
              title_fontsize=11)

    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()

    # Save Hi-Res
    fig.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_path} ({output_path.stat().st_size // 1024} KB)")


def main():
    parser = argparse.ArgumentParser(description="Render full heterogeneous graph")
    parser.add_argument('strain', type=str, help='Strain name')
    parser.add_argument('--max-nodes', type=int, default=0,
                        help='Max nodes (0 = all)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output PNG path')
    parser.add_argument('--export-gexf', type=str, default=None,
                        help='Export graph to GEXF file (for Gephi)')
    parser.add_argument('--export-sif', type=str, default=None,
                        help='Export graph to SIF file (for Cytoscape)')
    args = parser.parse_args()

    strain_dir = Path('outputs') / args.strain
    if not strain_dir.exists():
        print(f"Output directory not found: {strain_dir.resolve()}")
        print("Run the pipeline first: python main.py --strain", args.strain)
        sys.exit(1)

    print(f"Loading strain data from {strain_dir}...")
    node_types, edge_types = load_strain_data(strain_dir)

    total_nodes = sum(len(n) for n in node_types.values())
    total_edges = sum(len(e) for e in edge_types.values())
    print(f"  Found {total_nodes} nodes, {total_edges} edges across {len(edge_types)} types")

    max_nodes = args.max_nodes if args.max_nodes > 0 else None
    G = build_graph(node_types, edge_types, max_nodes=max_nodes)

    label = f" (sampled {max_nodes} nodes)" if max_nodes else " (full)"
    title_fragment = f"{args.strain} — Heterogeneous Multi-Omics Graph{label}"

    # ── Export GEXF (for Gephi) ──
    if args.export_gexf:
        gexf_path = Path(args.export_gexf)
        nx.write_gexf(G, str(gexf_path))
        print(f"  Exported GEXF: {gexf_path} ({gexf_path.stat().st_size // 1024} KB)")
        print(f"  Open in Gephi: File > Open > {gexf_path}")

    # ── Export SIF (for Cytoscape) ──
    if args.export_sif:
        sif_path = Path(args.export_sif)
        with open(sif_path, 'w') as f:
            for u, v, d in G.edges(data=True):
                rel = d.get('relation', 'interacts_with')
                node1_type = G.nodes[u].get('type', 'gene')
                node2_type = G.nodes[v].get('type', 'gene')
                f.write(f"{u}\t{rel}\t{v}\n")
        print(f"  Exported SIF: {sif_path} ({sif_path.stat().st_size // 1024} KB)")
        print(f"  Open in Cytoscape: File > Import > Network > File > {sif_path}")

    output_path = Path(args.output) if args.output else strain_dir / 'figures' / 'heterogeneous_graph_full.png'

    plot_graph(G, node_types, edge_types, output_path, title_fragment)

    # Also save a sampled version for quick viewing
    if max_nodes is None and total_nodes > 1000:
        sampled_path = output_path.parent / 'heterogeneous_graph_sampled.png'
        G_sample = build_graph(node_types, edge_types, max_nodes=500)
        plot_graph(G_sample, node_types, edge_types, sampled_path,
                   f"{args.strain} — Heterogeneous Graph (sampled 500 nodes)")


if __name__ == '__main__':
    main()
