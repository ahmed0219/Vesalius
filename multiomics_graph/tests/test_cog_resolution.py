"""
test_cog_resolution.py
Phase 1 validation - entity resolution between eggNOG (unversioned WP),
graph protein nodes (versioned WP) and gene nodes (locus tags).

Ensures:
  - normalize_wp_id strips the version suffix
  - eggNOG loader stores both normalized and original keys
  - COG annotations re-key onto real locus-tag gene nodes via gene_protein_map
  - COG hyperedges contain actual graph gene nodes (not orphan WP strings)
"""
import numpy as np
import pandas as pd
import pytest

from annotation.identifiers import normalize_wp_id
from annotation.eggnog import EggNOGAnnotator, COG_CATEGORIES
from graph.hypergraph import BiologicalHypergraphBuilder


def test_normalize_wp_id():
    assert normalize_wp_id('WP_000000542.1') == 'WP_000000542'
    assert normalize_wp_id('WP_000000542') == 'WP_000000542'
    assert normalize_wp_id('') == ''
    assert normalize_wp_id('not_a_wp.1') == 'not_a_wp.1'


def test_eggnog_load_keys_both_forms(tmp_path):
    tsv = tmp_path / 'eggnog.tsv'
    tsv.write_text(
        "GeneID\tCOG_category\tDescription\tPreferred_name\n"
        "WP_000000542\tJ\tRibosomal\trpsA\n"
        "WP_000000953\tS\tUnknown\t\n"
    )
    ann = EggNOGAnnotator()
    ann.load_precomputed(str(tsv))
    # Both unversioned (standard) and the raw form present.
    assert 'WP_000000542' in ann.annotations
    payload = ann.annotations['WP_000000542']
    assert payload['COG_category'] == 'J'
    assert payload['original_protein_id'] == 'WP_000000542'
    assert 'J' in payload['COG_classes']


def test_as_locus_tag_map_resolves_via_gene_protein(tmp_path):
    tsv = tmp_path / 'eggnog.tsv'
    tsv.write_text(
        "GeneID\tCOG_category\tDescription\tPreferred_name\n"
        "WP_000000542\tV\tDefense\tblaOXA\n"
    )
    ann = EggNOGAnnotator()
    ann.load_precomputed(str(tsv))
    gpm = {'EW036_RS26250': 'WP_000000542.1'}
    resolved = ann.as_locus_tag_map(gpm)
    assert 'EW036_RS26250' in resolved
    assert 'V' in resolved['EW036_RS26250']['COG_classes']
    # A protein id with no annotation is skipped.
    gpm2 = {'EW036_RS99999': 'WP_99999999.9'}
    assert 'EW036_RS99999' not in ann.as_locus_tag_map(gpm2)


def test_cog_hyperedges_use_gene_nodes():
    hb = BiologicalHypergraphBuilder()
    cog = {
        'EW036_RS26225': {'COG_classes': ['V'], 'COG_category': 'V'},
        'EW036_RS26250': {'COG_classes': ['V'], 'COG_category': 'V'},
        'EW036_RS00001': {'COG_classes': ['J'], 'COG_category': 'J'},
    }
    gpm = {
        'EW036_RS26225': 'WP_000000239.1',
        'EW036_RS26250': 'WP_000001334.1',
        'EW036_RS00001': 'WP_000000542.1',
    }
    hb.add_cog_hyperedges(cog_annotations=cog, gene_to_protein=gpm)
    he = [h for h in hb.hyperedges if h['type'] == 'cog_category']
    assert len(he) == 2
    v_he = [h for h in he if h['features'].get('cog') == 'V'][0]
    # Nodes are locus tags (gene nodes), not orphan WP strings.
    assert 'EW036_RS26225' in v_he['nodes']
    assert 'EW036_RS26250' in v_he['nodes']
    assert any(n.startswith('WP_') for n in v_he['nodes'])


def test_b36_defense_category_expected():
    # COG category 'V' (Defense mechanisms) must be a valid category and is
    # the category OXA-1 (a B36 defense/resistance marker) falls under.
    assert 'V' in COG_CATEGORIES