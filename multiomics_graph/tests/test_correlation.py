"""
test_correlation.py
Co-expression (transcriptional_correlation) and abundance-correlation edges.

Ensures:
  - replicate-level data produces correlation edges (mutual top-k capped)
  - condition-means-only input degrades to zero correlation edges
  - the mutual top-k cap keeps the graph sparse vs the raw threshold rule
  - protein abundance correlation uses the same replicate-column detection
"""
import numpy as np
import pandas as pd
import torch

from graph.heterogeneous_graph import HeterogeneousGraphBuilder

try:
    from torch_geometric.data import HeteroData
    _HAS_PYG = True
except Exception:
    _HAS_PYG = False


def _make_builder(top_k=2, threshold=0.5):
    b = HeterogeneousGraphBuilder(
        correlation_threshold=threshold, correlation_top_k=top_k)
    if _HAS_PYG:
        b.data = HeteroData()
    else:
        b.data = {}
    return b


def _replicate_frame(gids, n_reps=6, seed=1):
    rng = np.random.default_rng(seed)
    cols = {f'rep{i}' for i in range(n_reps)}
    data = {g: [float(v) for v in rng.normal(0, 1, n_reps)] for g in gids}
    # Force two strongly co-expressed clusters so edges are expected.
    for g in gids[:4]:
        data[g] = [float(v) + 5.0 * i for i, v in enumerate(data[g])]
    return pd.DataFrame.from_dict(data, orient='index', columns=sorted(cols))


def test_replicates_produce_correlation_edges():
    gids = [f'g{i:04d}' for i in range(12)]
    rep = _replicate_frame(gids)
    b = _make_builder(top_k=3)
    b.add_gene_nodes(gids, pd.DataFrame(index=gids), pd.DataFrame(index=gids))
    b.add_transcriptional_correlation_edges(gids, rep)
    key = ('gene', 'transcriptional_correlation', 'gene')
    if not _HAS_PYG:
        return
    ei = b.data[key].edge_index
    assert ei.size(1) > 0
    assert b.data[key].edge_attr is not None
    assert b.data[key].edge_attr.size(0) == ei.size(1)


def test_means_only_degrades_to_zero_edges():
    gids = [f'g{i:04d}' for i in range(12)]
    means = pd.DataFrame({
        'RNA_RPMI': np.zeros(len(gids)),
        'RNA_Sera': np.ones(len(gids)),
        'RNA_logFC': np.ones(len(gids)),
    }, index=gids)
    b = _make_builder()
    b.add_gene_nodes(gids, pd.DataFrame(index=gids), pd.DataFrame(index=gids))
    b.add_transcriptional_correlation_edges(gids, means)
    if not _HAS_PYG:
        return
    key = ('gene', 'transcriptional_correlation', 'gene')
    # Means-only input must not fabricate correlation edges.
    assert key not in b.data


def test_top_k_caps_edge_count():
    rng = np.random.default_rng(0)
    gids = [f'g{i:04d}' for i in range(60)]
    rep = pd.DataFrame(
        rng.normal(0, 1, (60, 8)).astype(float),
        index=gids, columns=[f'rep{i}' for i in range(8)])

    sparse = _make_builder(top_k=2)
    sparse.add_gene_nodes(gids, pd.DataFrame(index=gids),
                          pd.DataFrame(index=gids))
    sparse.add_transcriptional_correlation_edges(gids, rep)
    dense = _make_builder(top_k=0)
    dense.add_gene_nodes(gids, pd.DataFrame(index=gids),
                         pd.DataFrame(index=gids))
    dense.add_transcriptional_correlation_edges(gids, rep)
    if not _HAS_PYG:
        return
    key = ('gene', 'transcriptional_correlation', 'gene')
    n_sparse = sparse.data[key].edge_index.size(1)
    n_dense = dense.data[key].edge_index.size(1)
    assert n_sparse < n_dense
    # Mutual top-2 over 60 nodes is at most ~2*60 undirected pairs.
    assert n_sparse <= 60 * 2


def test_protein_abundance_correlation():
    pids = [f'WP_{i:09d}.1' for i in range(12)]
    rep = _replicate_frame(pids)
    b = _make_builder(top_k=3)
    b.add_protein_nodes(pids, pd.DataFrame(index=pids))
    b.add_abundance_correlation_edges(pids, rep)
    if not _HAS_PYG:
        return
    key = ('protein', 'abundance_correlation', 'protein')
    assert b.data[key].edge_index.size(1) > 0


def test_metabolite_correlation():
    mids = [f'CHEBI:{i}' for i in range(12)]
    rep = _replicate_frame(mids)
    b = _make_builder(top_k=3)
    b.add_metabolite_nodes(mids, pd.DataFrame(index=mids))
    b.add_metabolite_correlation_edges(mids, rep)
    if not _HAS_PYG:
        return
    key = ('metabolite', 'abundance_correlation', 'metabolite')
    assert b.data[key].edge_index.size(1) > 0
    assert b.data[key].edge_attr is not None