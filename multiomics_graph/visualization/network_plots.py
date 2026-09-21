"""
network_plots.py
Visualisation functions for the heterogeneous multi-omics graph.

Generates:
  - Heterogeneous graph layout (colored by node type)
  - Hypergraph incidence matrix heatmap
  - Edge type distribution bar chart
  - Node feature distributions
  - Regulation summary (RNA vs Protein concordance)
  - COG category frequency plot
"""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
from collections import Counter


# Color scheme
COLORS = {
    'gene': '#2E86AB',
    'protein': '#A23B72',
    'annotation': '#F18F01',
    'metabolite': '#1B998B',
    'genomic_proximity': '#73AB84',
    'transcriptional_correlation': '#E8A838',
    'encodes': '#C44536',
    'ppi': '#5D576B',
    'multi_omics_triple': '#1B998B',
    'kegg_pathway': '#FF6B6B',
    'cog_category': '#4ECDC4',
    'go_term': '#FFE66D',
    'rna': '#3B8EA5',
    'protein_omics': '#F08080',
    'amr_mechanism': '#E63946',
    'amr_gene': '#D90429',
}

AMR_CLASS_COLORS = {
    'beta_lactam': '#D90429',
    'aminoglycoside': '#F4A261',
    'sulfonamide': '#2A9D8F',
    'trimethoprim': '#8D99AE',
    'tetracycline': '#E9C46A',
    'macrolide': '#B5179E',
    'fluoroquinolone': '#4C956C',
    'phenicol': '#F77F00',
    'fosfomycin': '#003049',
}

EDGE_STYLES = {
    'genomic_proximity': 'solid',
    'transcriptional_correlation': 'dashed',
    'encodes': 'dotted',
    'ppi': 'dashdot',
}


class NetworkVisualizer:
    """
    Visualise the heterogeneous multi-omics graph and hypergraph.

    Parameters
    ----------
    output_dir : str
        Directory to save figures
    """

    def __init__(self, output_dir: str = 'outputs/figures'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def plot_heterogeneous_graph(
        self,
        node_types: dict,
        edge_types: dict,
        title: str = "Heterogeneous Multi-Omics Graph",
        filename: str = "heterogeneous_graph.png",
        max_nodes: int = 0,
    ):
        """
        Plot the heterogeneous multi-omics network.

        Parameters
        ----------
        node_types : dict of {type: list of node_ids}
        edge_types : dict of {(src, rel, dst): [(s, d), ...]}
        title : str
        filename : str
        max_nodes : int
            Max nodes (0 = all). When > 0, samples highest-degree nodes.
        """
        G = nx.Graph()

        # Degree-based sampling when max_nodes is set
        if max_nodes > 0:
            degree_count = {}
            for (_, _, _), edges in edge_types.items():
                for s, d in edges:
                    degree_count[s] = degree_count.get(s, 0) + 1
                    degree_count[d] = degree_count.get(d, 0) + 1
            sorted_nodes = sorted(degree_count.items(), key=lambda x: -x[1])
            kept = set(n for n, _ in sorted_nodes[:max_nodes])
            # Ensure at least one node per type
            for ntype, nodes in node_types.items():
                if not any(n in kept for n in nodes):
                    for n, _ in sorted_nodes:
                        if n in nodes:
                            kept.add(n)
                            break
        else:
            kept = set()
            for nodes in node_types.values():
                kept.update(nodes)

        node_color_map = []
        for ntype, nodes in node_types.items():
            for n in nodes:
                if n in kept:
                    G.add_node(n, type=ntype)
                    node_color_map.append(COLORS.get(ntype, '#888888'))

        for (src, rel, dst), edges in edge_types.items():
            for s, d in edges:
                if s in G and d in G:
                    G.add_edge(s, d, relation=rel)

        if G.number_of_nodes() == 0:
            print("  [Viz] Empty graph, skipping")
            return

        # Layout — scale k by node count
        k_val = 1.0 / (G.number_of_nodes() ** 0.3)
        pos = nx.spring_layout(G, seed=42, k=k_val, iterations=50)

        figsize = (20, 16) if G.number_of_nodes() > 1000 else (14, 10)
        fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Draw edges grouped by relation type
        for rel in set(
            data.get('relation', '') for _, _, data in G.edges(data=True)
        ):
            edges = [
                (u, v) for u, v, d in G.edges(data=True)
                if d.get('relation') == rel
            ]
            if edges:
                nx.draw_networkx_edges(
                    G, pos, edgelist=edges, ax=ax,
                    edge_color=COLORS.get(rel, '#888888'),
                    style=EDGE_STYLES.get(rel, 'solid'),
                    alpha=0.4, width=0.8,
                )

        # Draw nodes
        nx.draw_networkx_nodes(
            G, pos, ax=ax, node_size=40,
            node_color=node_color_map, alpha=0.8,
            edgecolors='white', linewidths=0.3,
        )

        # Legend
        legend_elements = []
        for ntype, color in COLORS.items():
            if ntype in node_types and node_types[ntype]:
                legend_elements.append(
                    mpatches.Patch(color=color, label=ntype, alpha=0.8)
                )
        ax.legend(
            handles=legend_elements, fontsize=8, loc='upper right',
            title='Node Types', title_fontsize=9,
        )

        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.axis('off')

        plt.tight_layout()
        path = self.output_dir / filename
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [Saved] {path}")

    def plot_amr_highlight(
        self,
        node_types: dict,
        edge_types: dict,
        amr_loci: list,
        amr_class_map: dict,
        title: str = "AMR Determinants in the Multi-Omics Graph",
        filename: str = "amr_highlight.png",
        max_nodes: int = 0,
        genome_only_loci: set = None,
        mechanism_nodes: list = None,
        verify: bool = True,
    ):
        """
        Plot the heterogeneous graph highlighting AMR determinant genes.

        AMR genes are enlarged and colored by their resistance mechanism
        class; all other nodes are faint greyed out so resistance loci stand
        out. Requires networkx spring layout.

        Every known AMR locus from the manifest is kept as a gene node even
        when it has no RNA/protein quantitation (a genome-only determinant is
        still a real genomic marker). Such nodes are drawn with an outlined /
        hatched style to distinguish them from multi-omics-supported AMR genes
        (solid fill). Their omics status is never invented.

        Parameters
        ----------
        node_types : dict of {type: [node_ids]}
        edge_types : dict of {(src, rel, dst): [(s, d), ...]}
        amr_loci : list
            Locus tags flagged as AMR determinants.
        amr_class_map : dict
            {locus_tag: amr_class} resolved from the manifest.
        genome_only_loci : set, optional
            Locus tags (same keys as ``amr_class_map``) present only at the
            genome level (no RNA/protein quantitation).
        mechanism_nodes : list, optional
            ``AMR:<class>`` mechanism node ids to draw as diamonds.
        verify : bool
            When True, assert every ``amr_class_map`` key is present in the
            drawn graph and return the per-node state map. Default True.

        Returns
        -------
        dict or None
            {node_id: 'multi_omics'|'genome_only'} for every AMR locus node in
            the figure; None if the graph was empty.
        """
        G = nx.Graph()
        kept = set()
        for nodes in node_types.values():
            kept.update(nodes)

        if not amr_class_map and not mechanism_nodes:
            print("  [Viz] No AMR loci, skipping amr_highlight")
            return None

        for ntype, nodes in node_types.items():
            for n in nodes:
                if n in kept:
                    G.add_node(n, type=ntype)
        for (src, rel, dst), edges in edge_types.items():
            for s, d in edges:
                if s in G and d in G:
                    G.add_edge(s, d, relation=rel)

        if G.number_of_nodes() == 0:
            print("  [Viz] Empty graph, skipping")
            return None

        # For very large graphs, restrict the drawn/layout subgraph to the
        # AMR molecular context: every known AMR determinant, its 1-hop
        # neighbourhood, and the mechanism nodes. This keeps the figure fast
        # and readable while guaranteeing every AMR locus stays in the figure.
        amr_set = set(amr_class_map)
        mech_set = set(mechanism_nodes or [])
        if len(G) > 2500 and (amr_set or mech_set):
            keep = set(amr_set)
            keep.update(mech_set)
            for n in list(amr_set):
                keep.update(G.neighbors(n))
            G = G.subgraph(keep).copy()

        k_val = 1.0 / (G.number_of_nodes() ** 0.3)
        pos = nx.spring_layout(G, seed=42, k=k_val, iterations=50)

        fig, ax = plt.subplots(1, 1, figsize=(20, 16) if len(G) > 1000 else (14, 10))

        # Faint background nodes + edges
        nx.draw_networkx_edges(
            G, pos, ax=ax, edge_color='#CCCCCC', alpha=0.15, width=0.5,
        )
        non_amr = [n for n in G.nodes if n not in amr_class_map
                   and n not in (mechanism_nodes or [])]
        nx.draw_networkx_nodes(
            G, pos, ax=ax, nodelist=non_amr,
            node_size=30, node_color='#DDDDDD', alpha=0.6,
            edgecolors='white', linewidths=0.2,
        )

        # Degree-sampled subset if requested
        if max_nodes > 0 and len(amr_class_map) > max_nodes:
            deg = dict(G.degree())
            amr_class_map = dict(sorted(
                amr_class_map.items(), key=lambda kv: -deg.get(kv[0], 0)
            )[:max_nodes])

        # Never silently drop a known AMR determinant: all locus keys were made
        # part of ``node_types['gene']`` by the caller before this point. Nodes
        # not reachable through a pairwise edge are still valid genome-only
        # determinants and are drawn (outlined) below.
        genome_only_loci = set(genome_only_loci or [])
        amr_state = {}
        by_class = {}
        for locus, cls in amr_class_map.items():
            if locus not in G:
                G.add_node(locus, type='gene')
            by_class.setdefault(cls, []).append(locus)
            amr_state[locus] = 'genome_only' if locus in genome_only_loci \
                else 'multi_omics'

        legend = []
        for cls, loci in by_class.items():
            color = AMR_CLASS_COLORS.get(cls, '#D90429')
            solid = [n for n in loci if amr_state[n] == 'multi_omics']
            genome_only = [n for n in loci if amr_state[n] == 'genome_only']
            if solid:
                nx.draw_networkx_nodes(
                    G, pos, ax=ax, nodelist=solid,
                    node_size=260, node_color=color, alpha=0.95,
                    edgecolors='black', linewidths=1.0,
                )
            if genome_only:
                # Outlined / hatched -> genome-only determinant (no invented
                # RNA/protein values; status stays explicitly genome_only).
                ax.scatter(
                    [pos[n][0] for n in genome_only],
                    [pos[n][1] for n in genome_only],
                    s=210, facecolors='none',
                    edgecolors=color,
                    linewidths=1.2, hatch='///', alpha=0.9,
                )
            legend.append(mpatches.Patch(color=color, label=cls, alpha=0.95))

        # AMR mechanism / class nodes (diamonds) — representation of the
        # curated resistance hyperedge classes (beta-lactam, aminoglycoside,
        # sulfonamide, trimethoprim, tetracycline, macrolide, ...).
        mech = mechanism_nodes or []
        if mech:
            mech_drawn = [n for n in mech if n in G]
            nx.draw_networkx_nodes(
                G, pos, ax=ax, nodelist=mech_drawn,
                node_size=340, node_color=COLORS.get('amr_mechanism', '#E63946'),
                node_shape='D', alpha=0.9, edgecolors='black', linewidths=0.8,
            )
            legend.append(mpatches.Patch(
                color=COLORS.get('amr_mechanism', '#E63946'),
                label='AMR mechanism (class)', alpha=0.9,
            ))

        # Legend: class colors + genome-only marker
        legend.append(mpatches.Patch(
            facecolor='white', edgecolor='grey', hatch='///',
            label='AMR gene (genome-only)',
        ))
        legend.append(mpatches.Patch(
            facecolor='white', edgecolor='grey',
            label='AMR gene (multi-omics)',
        ))
        ax.legend(
            handles=legend, fontsize=9, loc='upper right',
            title='AMR Class',
        )
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.axis('off')
        plt.tight_layout()
        path = self.output_dir / filename
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [Saved] {path}")

        if verify:
            missing = [n for n in amr_class_map if n not in amr_state]
            if missing:
                raise RuntimeError(
                    f"AMR loci absent from figure {filename}: {missing}"
                )
        return amr_state

    def plot_incidence_matrix(
        self,
        H: np.ndarray,
        hyperedge_types: dict = None,
        filename: str = "incidence_matrix.png",
    ):
        """
        Plot hypergraph incidence matrix as a heatmap.

        Parameters
        ----------
        H : np.ndarray (n_nodes, n_hyperedges)
        hyperedge_types : dict, optional
            {type_name: [hyperedge_indices]}
        filename : str
        """
        if H.size == 0:
            return

        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

        # Downsample if too large
        max_display = 500
        if H.shape[0] > max_display or H.shape[1] > max_display:
            H_display = H[:max_display, :max_display]
            suffix = f" ({min(max_display, H.shape[0])} nodes × {min(max_display, H.shape[1])} hyperedges)"
        else:
            H_display = H
            suffix = ""

        im = ax.imshow(H_display, aspect='auto', cmap='Blues',
                       interpolation='nearest')
        ax.set_xlabel('Hyperedges', fontsize=11)
        ax.set_ylabel('Nodes', fontsize=11)
        ax.set_title(f'Hypergraph Incidence Matrix{suffix}',
                     fontsize=12, fontweight='bold')

        cbar = plt.colorbar(im, ax=ax, shrink=0.6)
        cbar.set_label('H[i,j]', fontsize=10)

        # Mark hyperedge type boundaries if available
        if hyperedge_types:
            current = 0
            for htype, indices in hyperedge_types.items():
                if indices:
                    current += len(indices)
                    ax.axvline(x=current - 0.5, color='red',
                               linewidth=0.5, alpha=0.5)

        plt.tight_layout()
        path = self.output_dir / filename
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [Saved] {path}")

    def plot_regulation_concordance(
        self,
        aligned_data: pd.DataFrame,
        filename: str = "regulation_concordance.png",
    ):
        """
        Plot RNA vs Protein regulation concordance scatter.

        Parameters
        ----------
        aligned_data : pd.DataFrame
            Must have RNA_logFC, Protein_logFC columns
        filename : str
        """
        if 'RNA_logFC' not in aligned_data or 'Protein_logFC' not in aligned_data:
            return

        mask = (
            aligned_data['RNA_logFC'].notna()
            & aligned_data['Protein_logFC'].notna()
        )
        if mask.sum() == 0:
            return

        rna = aligned_data.loc[mask, 'RNA_logFC'].values
        prot = aligned_data.loc[mask, 'Protein_logFC'].values
        agreement = rna * prot

        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

        # Scatter
        sc = axes[0].scatter(rna, prot, c=agreement, cmap='RdBu_r',
                             alpha=0.6, s=20, edgecolors='grey', linewidth=0.2)
        axes[0].axhline(y=0, color='grey', linestyle='--', linewidth=0.8)
        axes[0].axvline(x=0, color='grey', linestyle='--', linewidth=0.8)
        axes[0].set_xlabel('RNA logFC (Serum vs RPMI)', fontsize=11)
        axes[0].set_ylabel('Protein logFC (Serum vs RPMI)', fontsize=11)
        axes[0].set_title('RNA-Protein Regulation Concordance',
                          fontsize=12, fontweight='bold')
        cbar = plt.colorbar(sc, ax=axes[0])
        cbar.set_label('Agreement Score\n(RNA_logFC × Protein_logFC)',
                       fontsize=9)

        # Distribution of agreement scores
        axes[1].hist(agreement, bins=50, color='#2E86AB', alpha=0.7,
                     edgecolor='white', linewidth=0.5)
        axes[1].axvline(x=0, color='red', linestyle='--', linewidth=1.5)
        axes[1].set_xlabel('Agreement Score', fontsize=11)
        axes[1].set_ylabel('Frequency', fontsize=11)
        axes[1].set_title('Agreement Score Distribution',
                          fontsize=12, fontweight='bold')

        concordant = (agreement > 0).sum()
        discordant = (agreement < 0).sum()
        axes[1].text(
            0.95, 0.95,
            f"Concordant: {concordant}\nDiscordant: {discordant}",
            transform=axes[1].transAxes, fontsize=9, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
        )

        plt.tight_layout()
        path = self.output_dir / filename
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [Saved] {path}")

    def plot_cog_distribution(
        self,
        cog_counts: dict,
        filename: str = "cog_distribution.png",
        cog_categories: dict = None,
    ):
        if not cog_counts:
            return

        cat = cog_categories or {}
        letters = sorted(cog_counts.keys())
        counts = [cog_counts.get(l, 0) for l in letters]
        names = [cat.get(l, f'COG {l}') for l in letters]

        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        bars = ax.bar(range(len(letters)), counts, color='#4ECDC4',
                      alpha=0.8, edgecolor='white', linewidth=0.5)

        ax.set_xticks(range(len(letters)))
        ax.set_xticklabels(letters, fontsize=10)
        ax.set_xlabel('COG Category', fontsize=11)
        ax.set_ylabel('Number of Genes', fontsize=11)
        ax.set_title('COG Functional Category Distribution',
                     fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Annotate
        for bar, count, name in zip(bars, counts, names):
            if count > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + max(counts) * 0.01,
                        str(count), ha='center', va='bottom', fontsize=7,
                        rotation=90)

        plt.tight_layout()
        path = self.output_dir / filename
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [Saved] {path}")

    def plot_hyperedge_type_distribution(
        self,
        hyperedge_type_counts: dict,
        filename: str = "hyperedge_types.png",
    ):
        """
        Plot bar chart of hyperedge types.

        Parameters
        ----------
        hyperedge_type_counts : dict
            {type_name: count}
        filename : str
        """
        if not hyperedge_type_counts:
            return

        types = list(hyperedge_type_counts.keys())
        counts = list(hyperedge_type_counts.values())
        colors = [COLORS.get(t, '#888888') for t in types]

        fig, ax = plt.subplots(1, 1, figsize=(8, 5))
        bars = ax.barh(range(len(types)), counts, color=colors,
                       alpha=0.8, edgecolor='white')

        ax.set_yticks(range(len(types)))
        ax.set_yticklabels(types, fontsize=10)
        ax.set_xlabel('Number of Hyperedges', fontsize=11)
        ax.set_title('Hyperedge Type Distribution',
                     fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        for bar, count in zip(bars, counts):
            ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                    str(count), ha='left', va='center', fontsize=9)

        plt.tight_layout()
        path = self.output_dir / filename
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [Saved] {path}")

    def plot_metabolite_regulation(
        self,
        metabolite_table: pd.DataFrame,
        filename: str = "metabolite_regulation.png",
    ):
        """
        Plot metabolite logFC distribution and top-changing metabolites.

        Parameters
        ----------
        metabolite_table : pd.DataFrame
            Must have MetaboliteName, Metabolite_logFC columns
        filename : str
        """
        if 'Metabolite_logFC' not in metabolite_table or \
                'MetaboliteName' not in metabolite_table:
            return
        t = metabolite_table.dropna(subset=['Metabolite_logFC'])
        if t.empty:
            return

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # logFC histogram
        axes[0].hist(t['Metabolite_logFC'], bins=40, color='#1B998B',
                     alpha=0.7, edgecolor='white', linewidth=0.5)
        axes[0].axvline(x=0, color='grey', linestyle='--', linewidth=1)
        axes[0].set_xlabel('Metabolite logFC (Sera vs RPMI)', fontsize=11)
        axes[0].set_ylabel('Number of Metabolites', fontsize=11)
        axes[0].set_title('Metabolite Regulation Distribution',
                          fontsize=12, fontweight='bold')
        axes[0].spines['top'].set_visible(False)
        axes[0].spines['right'].set_visible(False)

        # Top up/down
        top = t.nlargest(15, 'Metabolite_logFC')
        bottom = t.nsmallest(15, 'Metabolite_logFC')
        sel = pd.concat([top, bottom].iloc[0:0] if False else [top, bottom])
        sel = sel.sort_values('Metabolite_logFC')
        names = [str(n)[:28] for n in sel['MetaboliteName']]
        colors = ['#C44536' if v < 0 else '#2E86AB'
                  for v in sel['Metabolite_logFC']]
        axes[1].barh(range(len(sel)), sel['Metabolite_logFC'], color=colors,
                     alpha=0.8, edgecolor='white')
        axes[1].set_yticks(range(len(sel)))
        axes[1].set_yticklabels(names, fontsize=7)
        axes[1].set_xlabel('logFC', fontsize=11)
        axes[1].set_title('Top Changed Metabolites', fontsize=12,
                          fontweight='bold')
        axes[1].axvline(x=0, color='grey', linestyle='--', linewidth=1)
        axes[1].spines['top'].set_visible(False)
        axes[1].spines['right'].set_visible(False)

        plt.tight_layout()
        path = self.output_dir / filename
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [Saved] {path}")

    def plot_summary_dashboard(
        self, genome_summary: dict, rna_summary: dict,
        protein_summary: dict, graph_summary: dict,
        hypergraph_summary: dict,
        filename: str = "summary_dashboard.png",
        title: str = None,
    ):
        """
        Create a comprehensive summary dashboard figure.

        Parameters
        ----------
        genome_summary, rna_summary, protein_summary : dict
        graph_summary, hypergraph_summary : dict
        filename : str
        """
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))

        # 1. Omics layer sizes
        ax = axes[0, 0]
        labels = ['Genome', 'Transcriptome', 'Proteome']
        sizes = [
            int(genome_summary.get('genes') or 0),
            int(rna_summary.get('genes') or 0),
            int(protein_summary.get('proteins') or 0),
        ]
        colors = [COLORS['gene'], '#3B8EA5', COLORS['protein']]
        bars = ax.bar(labels, sizes, color=colors, alpha=0.8, edgecolor='white')
        for bar, s in zip(bars, sizes):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                    str(s), ha='center', va='bottom', fontsize=9)
        ax.set_title('Data Layer Sizes', fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # 2. RNA regulation
        ax = axes[0, 1]
        reg = ['Up', 'Down', 'Stable']
        rna_counts = [rna_summary.get(f'{r.lower()}_regulated', 0) for r in reg]
        if rna_summary.get('stable', 0) > 0:
            rna_counts[2] = rna_summary['stable']
        if sum(rna_counts) == 0:
            rna_counts = [0, 0, max(int(rna_summary.get('genes') or 0), 1)]
        rna_colors = ['#C44536', '#2E86AB', '#888888']
        ax.pie(rna_counts, labels=reg, autopct='%1.1f%%',
               colors=rna_colors, startangle=90)
        ax.set_title('RNA Regulation', fontweight='bold')

        # 3. Protein regulation
        ax = axes[0, 2]
        prot_counts = [
            protein_summary.get('up_regulated', 0),
            protein_summary.get('down_regulated', 0),
            protein_summary.get('stable', 0),
        ]
        if sum(prot_counts) == 0:
            prot_counts = [0, 0, max(int(protein_summary.get('proteins') or 0), 1)]
        ax.pie(prot_counts, labels=reg, autopct='%1.1f%%',
               colors=rna_colors, startangle=90)
        ax.set_title('Protein Regulation', fontweight='bold')

        # 4. Graph summary
        ax = axes[1, 0]
        graph_labels = list(graph_summary.keys())
        graph_vals = list(graph_summary.values())
        ax.barh(range(len(graph_labels)), graph_vals, color='#5D576B',
                alpha=0.8, edgecolor='white')
        ax.set_yticks(range(len(graph_labels)))
        ax.set_yticklabels(graph_labels, fontsize=9)
        ax.set_title('Heterogeneous Graph', fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # 5. Hyperedge types
        ax = axes[1, 1]
        he_types = hypergraph_summary.get('hyperedge_types', {})
        if he_types:
            he_labels = list(he_types.keys())
            he_counts = list(he_types.values())
            he_colors = [COLORS.get(t, '#888888') for t in he_labels]
            ax.barh(range(len(he_labels)), he_counts, color=he_colors,
                    alpha=0.8, edgecolor='white')
            ax.set_yticks(range(len(he_labels)))
            ax.set_yticklabels(he_labels, fontsize=8)
            ax.set_title('Hyperedge Types', fontweight='bold')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            if he_counts:
                ax.set_xlim(0, max(he_counts) * 1.14)
            for i, v in enumerate(he_counts):
                ax.text(v, i, f' {int(v)}', va='center', fontsize=7)

        # 6. Hypergraph stats
        ax = axes[1, 2]
        ax.axis('off')
        stats_text = (
            f"Hypergraph Summary\n"
            f"{'─' * 20}\n"
            f"Nodes: {hypergraph_summary.get('n_nodes', 0):,}\n"
            f"Hyperedges: {hypergraph_summary.get('n_hyperedges', 0):,}\n"
            f"Mean size: {hypergraph_summary.get('mean_hyperedge_size', 0):.1f}\n"
            f"Incidence density: {hypergraph_summary.get('incidence_density', 0):.6f}"
        )
        ax.text(0.1, 0.5, stats_text, transform=ax.transAxes,
                fontsize=11, va='center', fontfamily='monospace')

        plt.suptitle('Multi-Omics Heterogeneous Graph Architecture'
                     + (f' - {title}' if title else ''),
                     fontsize=14, fontweight='bold', y=1.01)
        plt.tight_layout()
        path = self.output_dir / filename
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [Saved] {path}")
