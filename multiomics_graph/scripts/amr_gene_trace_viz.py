"""
scripts/amr_gene_trace_viz.py
Trace a single AMR determinant across the eight omics layers and visualise
how it propagates through the multi-omics graph.

For each requested (strain, marker) the script renders one two-panel figure:

  Panel A (top)  - Layer trace ladder. The determinant gene is the hub on the
                   left; each of the 8 layers (genome, transcriptome, proteome,
                   metabolome, KEGG, COG, GO, PPI/AMR) is a box on the right.
                   Box border = layer colour, box fill = layer status
                   (observed / not_observed / unavailable), box text = the key
                   value carried by that layer (logFC, counts, ids).

  Panel B (bottom) - Propagation network. The determinant's molecular context
                   in the entity graph (genes/proteins/metabolites/annotation
                   hubs) up to ``--max-depth`` hops, coloured by layer and with
                   relationship-labelled edges. Members of very large annotation
                   hubs are capped so broad pathways (e.g. ko01130) do not flood
                   the figure; the hubs themselves are always shown.

Read-only over `multiomics_graph/outputs/<strain>/` (plus the AMR manifest and
the saved reports); it never re-runs the pipeline. Reuses the tracing logic of
`scripts/amr_path_analysis.py`.

Usage (from multiomics_graph/):
    python scripts/amr_gene_trace_viz.py
    python scripts/amr_gene_trace_viz.py --strains B36 --marker blaOXA-1
    python scripts/amr_gene_trace_viz.py --strains MS_14386 --max-depth 3
    python scripts/amr_gene_trace_viz.py --output reports
"""

import argparse
import sys
from collections import deque
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import networkx as nx
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import amr_path_analysis as ap

# Layer-specific colours for boxes / nodes / edges.
LAYER_COLORS = {
    'genome': '#2E86AB',
    'transcriptome': '#3B8EA5',
    'proteome': '#A23B72',
    'metabolome': '#1B998B',
    'kegg': '#FF6B6B',
    'cog': '#4ECDC4',
    'go': '#F2C14E',
    'ppi_amr': '#5D576B',
    'amr_mechanism': '#E63946',
    'annotation': '#F18F01',
    'hypergraph': '#8E7DBE',
}

EDGE_COLORS = {
    'encodes': '#C44536',
    'genomic_proximity': '#73AB84',
    'in_pathway': '#FF6B6B',
    'cog_category': '#4ECDC4',
    'annotated_by': '#F2C14E',
    'associated_with': '#E63946',
    'member_of_cluster': '#8E7DBE',
    'ppi': '#5D576B',
    'metabolised_by': '#2A9D8F',
    'member_of': '#BBBBBB',
}

STATUS_STYLE = {
    'observed': {'face': '#E8F5E9', 'edge': '#2E7D32', 'ls': 'solid'},
    'not_observed': {'face': '#F5F5F5', 'edge': '#9E9E9E', 'ls': 'dashed'},
    'unavailable': {'face': '#FFFFFF', 'edge': '#E53935', 'ls': 'dotted'},
}

MAX_HUB_MEMBERS = 5
HUB_SIZE_LIMIT = 60
MAX_HUBS = 12


def _trunc(s, n=72):
    s = str(s).strip()
    return s if len(s) <= n else s[:n - 1] + '\u2026'


def _fmt(v, nd=2):
    try:
        f = float(v)
        if pd.isna(f):
            return 'n/a'
        return f'{f:.{nd}f}'
    except (TypeError, ValueError):
        return str(v)


def _layer_value(rec, layer):
    """Human-readable value for each layer in the trace record."""
    if layer == 'genome':
        return f"locus {rec['locus_tag']}\n{rec['L1_genomic_neighbors']} genomic neighbors"
    if layer == 'transcriptome':
        if rec['L2_rna_present']:
            return (f"RNA log2FC {_fmt(rec.get('L2_RNA_log2_derived'))} "
                    f"({rec['L2_RNA_direction']})")
        return 'no RNA quantitation'
    if layer == 'proteome':
        if rec['L3_protein_present']:
            return (f"Protein logFC {_fmt(rec.get('L3_Protein_logFC'))} "
                    f"({rec.get('L3_Protein_Regulation') or 'n/a'})")
        return 'no protein quantitation'
    if layer == 'metabolome':
        if rec['L4_metabolites']:
            return (f"{rec['L4_metabolites']} metabolites\n"
                    f"{_trunc(rec.get('L4_metabolite_names') or '', 68)}")
        return 'no metabolite link'
    if layer == 'kegg':
        if rec['L5_kegg_pathways']:
            return (f"{rec['L5_kegg_pathways']} KEGG pathways\n"
                    f"{_trunc(rec.get('L5_kegg_names') or '', 68)}")
        return 'no KEGG pathway'
    if layer == 'cog':
        return _trunc(rec.get('L6_cog_categories') or 'no COG assignment', 72)
    if layer == 'go':
        if rec['L7_go_terms']:
            return (f"{rec['L7_go_terms']} GO terms\n"
                    f"{_trunc(rec.get('L7_go_ids') or '', 68)}")
        return 'no GO annotation'
    if layer == 'ppi_amr':
        parts = [f"{rec['L8_ppi_neighbors']} PPI neighbors"]
        if rec.get('L8_ppi_clusters'):
            parts.append(_trunc(rec['L8_ppi_clusters'], 60))
        if rec.get('L8_amr_mechanisms'):
            parts.append('; '.join(rec['L8_amr_mechanisms']))
        return '\n'.join(parts)
    return ''


def _layer_status(rec, layer):
    return str(rec.get(ap._layer_key(layer), 'unavailable'))


def _prepare_strain_data(name):
    """Load a strain exactly like amr_path_analysis.main() does."""
    data, issues = ap._load_strain(name)
    if not data:
        return None, issues
    data['in_pathway_map'] = ap._edge_map(data['in_pathway'], 'gene_id', 'annotation_id')
    data['cog_map'] = ap._edge_map(data['cog'], 'gene_id', 'annotation_id')
    data['go_map'] = ap._wp_edge_map(data['go'], 'protein_id', 'annotation_id')
    data['clusters_map'] = ap._wp_edge_map(data['clusters'], 'protein_id', 'annotation_id')
    data['amr_edges_map'] = ap._edge_map(data['amr_edges'], 'gene_id', 'annotation_id')
    data['genomic_map'] = ap._edge_map(data['genomic'], 'gene_id', 'gene_target')
    data['_entity_graph'] = ap._build_entity_graph(data)
    data['_amr_loci'] = set(ap.MANIFEST.amr_loci(name))
    data['_amr_proteins'] = {
        p for rec in ap.MANIFEST.markers(name) for p in rec['protein_ids']}
    return data, issues


# ---------------------------------------------------------------------------
# Propagation graph (Panel B)
# ---------------------------------------------------------------------------

def _node_layer(node, data):
    if node in data.get('_metabolite_ids', set()):
        return 'metabolome'
    if node in data.get('_protein_ids', set()):
        return 'proteome'
    if node in data.get('_gene_ids', set()):
        return 'genome'
    if str(node).startswith('AMR:'):
        return 'amr_mechanism'
    return 'annotation'


def _node_label(node, layer):
    if layer == 'metabolome':
        return ap.CHEBI_NAMES.get(node, node)[:18]
    if layer in ('genome', 'proteome'):
        nm = (ap.ENTITY_NAMES.get(node)
              or ap.ENTITY_NAMES.get(ap._norm_wp(node)) or node)
        return str(nm)[:18]
    return str(node)[:24]


def _build_context_graph(data, locus, protein, max_depth):
    """Bounded context subgraph: determinant -> 1-hop -> capped 2-hop members."""
    g = data.get('_entity_graph') or {}
    merged = {}
    for s in (locus, protein):
        if not s:
            continue
        for node, (d, seq, eseq) in ap._bfs(g, s, max_depth).items():
            if node not in merged or d < merged[node][0]:
                merged[node] = (d, seq, eseq)

    measured = data.get('_quantified', {})
    measured_ids = set(measured.get('rna', set())) | set(measured.get('protein', set())) \
        | set(measured.get('metabolites', set()))

    layer_of = {n: _node_layer(n, data) for n in merged}

    # Annotation hubs directly attached to the determinant (d == 1).
    # Sort by *specificity* first: hubs with few biological co-members
    # (e.g. map00311 penicillin/cephalosporin biosynthesis -> 3 metabolites)
    # carry real information, while giant hubs (ko01130, 162 members) are
    # generic and are only shown if budget remains.
    def _hub_members(hub):
        return [(m[0], m[1]) for m in g.get(hub, [])
                if layer_of.get(m[0]) not in (None, 'annotation')
                and m[0] != locus and m[0] != protein]

    hubs = sorted(
        ((n, info, _hub_members(n)) for n, info in merged.items()
         if layer_of[n] == 'annotation' and info[0] == 1),
        key=lambda kv: len(kv[2]),
    )[:MAX_HUBS]

    # Keep direct biological neighbours (d == 1) in full.
    keep = set()
    for n, (d, _seq, _eseq) in merged.items():
        if d == 1 and layer_of[n] != 'annotation':
            keep.add(n)

    # Cap members reached through a hub (d == 2): skip giant hub pathways,
    # keep only a few measured-first co-members per hub.
    for hub, _info, hub_members in hubs:
        members = [m for m, _rel in hub_members
                   if layer_of.get(m) not in (None, 'annotation')
                   and m != locus and m != protein]
        if len(members) > HUB_SIZE_LIMIT:
            continue
        ordered = sorted(
            members,
            key=lambda n: (n not in measured_ids, len(g.get(n, []))),
        )
        for m in ordered[:MAX_HUB_MEMBERS]:
            if m in merged and merged[m][0] <= 2:
                keep.add(m)

    # Assemble the networkx graph from the surviving paths.
    G = nx.Graph()
    for n, (d, seq, eseq) in merged.items():
        if n in keep or n == locus or n == protein:
            G.add_node(n, layer=layer_of.get(n, 'annotation'), depth=d)
    for hub, _info, _hub_members in hubs:
        G.add_node(hub, layer='annotation', depth=1)
    for n, (d, seq, eseq) in merged.items():
        if n not in G:
            continue
        for prev, rel in zip(seq, eseq):
            if prev in G:
                G.add_edge(prev, n, relation=rel)
    return G


# ---------------------------------------------------------------------------
# Figure rendering
# ---------------------------------------------------------------------------

def _draw_layer_ladder(ax, rec, marker, strain):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    title = (f"{marker}  |  {strain}  |  class: {rec['amr_class']}  |  "
             f"status: {rec['status']}  |  layer coverage: "
             f"{rec['layer_coverage']}/8  |  {rec['layer_coverage_pct']}")
    ax.set_title(title, fontsize=11, fontweight='bold', pad=8)

    # Determinant gene hub.
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.02, 0.36), 0.21, 0.28, boxstyle='round,pad=0.008',
        fc='#D90429', ec='black', lw=1.2))
    ax.text(0.125, 0.58, marker, ha='center', va='center',
            fontsize=11, fontweight='bold', color='white')
    ax.text(0.125, 0.47, rec['locus_tag'], ha='center', va='center',
            fontsize=8, color='white')
    ax.text(0.125, 0.40, rec['protein_id'] or 'no protein',
            ha='center', va='center', fontsize=6.5, color='#FFE0B2')

    evidence = (rec['evidence_types'] or '').replace('; ', '\n')
    ax.text(0.02, 0.05, f"Evidence:\n{evidence}", fontsize=8,
            va='bottom', ha='left', family='monospace',
            bbox=dict(boxstyle='round', fc='#FAFAFA', ec='#BBBBBB'))

    rows = ['genome', 'transcriptome', 'proteome', 'metabolome',
            'kegg', 'cog', 'go', 'ppi_amr']
    for i, layer in enumerate(rows):
        col, row = i // 4, i % 4
        x = 0.36 + col * 0.31
        y = 0.88 - row * 0.20
        status = _layer_status(rec, layer)
        style = STATUS_STYLE.get(status, STATUS_STYLE['unavailable'])
        colour = LAYER_COLORS.get(layer, '#888888')

        # Gene -> layer connector.
        ax.annotate(
            '', xy=(x + 0.005, y + 0.075), xytext=(0.235, 0.50),
            arrowprops=dict(arrowstyle='-|>', color=colour, lw=1.4,
                            shrinkA=0, shrinkB=0),
        )

        ax.add_patch(mpatches.FancyBboxPatch(
            (x, y), 0.295, 0.17, boxstyle='round,pad=0.004',
            fc=style['face'], ec=colour, lw=1.4,
            linestyle=style['ls']))

        ax.text(x + 0.005, y + 0.145, f"{layer}  [{status}]",
                fontsize=8.5, fontweight='bold', color=colour,
                va='top', family='monospace')
        ax.text(x + 0.005, y + 0.115, _layer_value(rec, layer),
                fontsize=6.5, va='top', color='#222222', family='monospace')


def _draw_propagation(ax, G, marker, locus, protein):
    if G.number_of_nodes() == 0:
        ax.axis('off')
        ax.text(0.5, 0.5, 'no graph context', ha='center', va='center')
        return

    pos = nx.spring_layout(G, seed=42, k=1.0 / (G.number_of_nodes() ** 0.4),
                           iterations=80)
    nodes = list(G.nodes(data=True))

    # Edges grouped by relation.
    rels = sorted({d['relation'] for _, _, d in G.edges(data=True)})
    for rel in rels:
        edges = [(u, v) for u, v, d in G.edges(data=True)
                 if d['relation'] == rel]
        nx.draw_networkx_edges(
            G, pos, ax=ax, edgelist=edges,
            edge_color=EDGE_COLORS.get(rel, '#BBBBBB'),
            alpha=0.55, width=1.0,
        )

    # Nodes grouped by layer.
    layers = sorted({d['layer'] for _, d in nodes})
    for layer in layers:
        nlist = [n for n, d in nodes if d['layer'] == layer]
        size = 380 if layer == 'amr_mechanism' else \
            260 if layer == 'annotation' else 180
        nx.draw_networkx_nodes(
            G, pos, ax=ax, nodelist=nlist, node_size=size,
            node_color=LAYER_COLORS.get(layer, '#888888'),
            alpha=0.9, edgecolors='white', linewidths=0.5,
        )

    # Emphasise the determinant gene.
    seed = [n for n in (locus, protein) if n in G]
    if seed:
        nx.draw_networkx_nodes(
            G, pos, ax=ax, nodelist=seed, node_size=520,
            node_color='#D90429', alpha=0.95, edgecolors='black',
            linewidths=1.4, node_shape='s',
        )

    # Labels.
    labels = {}
    for n, d in nodes:
        if n == locus or n == protein:
            labels[n] = marker if n == locus else ('protein' if n == protein else n)
        else:
            labels[n] = _node_label(n, d['layer'])
    nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=6.5,
                            font_color='#111111')

    # Edge labels only for determinant-incident propagation edges.
    edge_labels = {
        (u, v): d['relation'] for u, v, d in G.edges(data=True)
        if u in seed or v in seed
    }
    if edge_labels:
        nx.draw_networkx_edge_labels(
            G, pos, edge_labels=edge_labels, ax=ax, font_size=6,
            font_color='#333333', bbox=dict(boxstyle='round,pad=0.1',
                                            fc='white', ec='none', alpha=0.7),
        )

    # Legends.
    node_legend = [mpatches.Patch(color=LAYER_COLORS.get(l, '#888888'),
                                  label=l) for l in layers]
    edge_legend = [Line2D([], [], color=EDGE_COLORS.get(r, '#888888'),
                          label=r, alpha=0.7) for r in rels]
    lg1 = ax.legend(handles=node_legend, fontsize=7, loc='upper left',
                    title='Node layer', title_fontsize=8,
                    bbox_to_anchor=(0.0, 1.0))
    ax.add_artist(lg1)
    ax.legend(handles=edge_legend, fontsize=7, loc='lower left',
              title='Relationship', title_fontsize=8,
              bbox_to_anchor=(0.0, 0.0))

    ax.set_title(f"Propagation: {marker} molecular context "
                 f"(direct links + capped co-annotation members)",
                 fontsize=11, fontweight='bold')
    ax.axis('off')


def render_trace(rec, data, out_path, max_depth=2, dpi=150):
    marker = rec['marker']
    strain = rec['strain']
    locus = rec['locus_tag']
    protein = rec['protein_id'] or None

    fig, (axA, axB) = plt.subplots(
        2, 1, figsize=(18, 15),
        gridspec_kw={'height_ratios': [1, 1.55]})

    _draw_layer_ladder(axA, rec, marker, strain)

    G = _build_context_graph(data, locus, protein, max_depth)
    _draw_propagation(axB, G, marker, locus, protein)

    fig.suptitle(f"AMR determinant trace — {marker} ({strain})",
                 fontsize=14, fontweight='bold', y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  [Saved] {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap_parser = argparse.ArgumentParser(
        description='Trace AMR determinants across omics layers + propagation')
    ap_parser.add_argument('--strains', nargs='*', default=ap.STRAIN_ORDER)
    ap_parser.add_argument('--marker', default=None,
                           help='Marker name (default: all markers of the strain)')
    ap_parser.add_argument('--output', default='reports',
                           help='Output directory (default: reports)')
    ap_parser.add_argument('--max-depth', type=int, default=2)
    ap_parser.add_argument('--dpi', type=int, default=150)
    args = ap_parser.parse_args()

    out_dir = ap.REPO / args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    made = 0
    for name in args.strains:
        data, issues = _prepare_strain_data(name)
        if not data:
            print(f"  [{name}] no outputs dir, skipping")
            continue
        markers = ap.MANIFEST.markers(name)
        if not markers:
            print(f"  [{name}] no AMR markers, skipping")
            continue
        if args.marker:
            markers = [m for m in markers
                       if m['name'].lower() == args.marker.lower()]
            if not markers:
                print(f"  [{name}] marker '{args.marker}' not found, skipping")
                continue
        for m in markers:
            rec = ap.trace_marker(name, m, data)
            rec['hyperedges'] = ap._hyperedge_memberships(
                data, rec['locus_tag'], rec['protein_id'] or None)
            safe = rec['marker'].replace('/', '_').replace(' ', '_')
            path = out_dir / f"amr_trace_{name}_{safe}.png"
            render_trace(rec, data, path, args.max_depth, args.dpi)
            made += 1

    print(f"\n{made} figure(s) written to: {out_dir.resolve()}")


if __name__ == '__main__':
    main()