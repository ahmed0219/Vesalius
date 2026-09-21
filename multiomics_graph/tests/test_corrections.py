"""
test_corrections.py
Targeted tests for the pipeline corrections:
  - 2-layer (RNA-less) multi-omics triples are masked (no NaN agreement)
  - per-instrument metabolomics normalisation
  - KEGG name-fallback prefers an exact-name match
"""
import numpy as np
import pandas as pd
import pytest

from graph.hypergraph import BiologicalHypergraphBuilder
from preprocessing.metabolomics import MetabolomicsProcessor
from annotation.kegg_compound import KEGGCompoundResolver


# ── 2-layer triples ───────────────────────────────────────────
def test_2layer_triple_agreement_masked():
    hb = BiologicalHypergraphBuilder()
    aligned = pd.DataFrame({
        'GeneID': ['g1', 'g2'],
        'ProteinID': ['p1', 'p2'],
        'RNA_logFC': [np.nan, np.nan],
        'Protein_logFC': [1.5, -0.5],
    })
    hb.add_multi_omics_triples(aligned)

    triples = [h for h in hb.hyperedges if h['type'] == 'multi_omics_triple']
    assert len(triples) == 2
    for t in triples:
        f = t['features']
        assert f['has_rna'] is False
        # Agreement must be explicitly masked, never a numeric NaN.
        assert f['agreement_score'] is None
        assert np.isfinite(f['Protein_logFC'])
        assert f['RNA_logFC'] == 0.0


def test_3layer_triple_agreement_computed():
    hb = BiologicalHypergraphBuilder()
    aligned = pd.DataFrame({
        'GeneID': ['g1'],
        'ProteinID': ['p1'],
        'RNA_logFC': [2.0],
        'Protein_logFC': [3.0],
    })
    hb.add_multi_omics_triples(aligned)
    t = hb.hyperedges[0]['features']
    assert t['has_rna'] is True
    assert t['agreement_score'] == 6.0


# ── Metabolomics per-instrument normalisation ───────────────────
def test_metabolomics_merges_cross_instrument(tmp_path):
    # Two "instruments": one on a small scale, one on a large scale.
    # Both detect Citrate; each also detects a private metabolite.
    gc = tmp_path / 'gc_ms.tsv'
    lc = tmp_path / 'lc_ms.tsv'
    gc.write_text(
        "metabolite_identification\tdatabase_identifier\tdatabase\tchemical_formula\tmass_to_charge\tretention_time\t50913\t50914\n"
        "Citrate\tCHEBI:16947\tMetabolomicsDB\tC6H8O7\t191.0\t1.0\t1\t3\n"
        "Malate\tCHEBI:6650\tMetabolomicsDB\tC4H6O5\t133.0\t1.0\t1\t3\n"
    )
    lc.write_text(
        "metabolite_identification\tdatabase_identifier\tdatabase\tchemical_formula\tmass_to_charge\tretention_time\t50913\t50914\n"
        "Citrate\tCHEBI:16947\tMetabolomicsDB\tC6H8O7\t191.0\t1.0\t100\t300\n"
        "Alanine\tCHEBI:5\tMetabolomicsDB\tC3H7NO2\t89.0\t1.0\t1000\t3000\n"
    )
    proc = MetabolomicsProcessor([str(gc), str(lc)])
    cond = {'50913': 'Sera', '50914': 'Sera'}
    table = proc.load(condition_map=cond)

    # Citrate detected by both instruments -> merged into one row.
    cit = table[table['MetaboliteName'] == 'Citrate']
    assert len(cit) == 1
    assert int(cit.iloc[0]['n_instruments']) == 2
    # Private metabolites survive with a single-instrument provenance.
    assert int(table[table['MetaboliteName'] == 'Alanine'].iloc[0]['n_instruments']) == 1
    assert int(table[table['MetaboliteName'] == 'Malate'].iloc[0]['n_instruments']) == 1
    # Resulting abundance is finite and ordered sensibly (Sera > RPMI absent:
    # all samples are Sera here, so both conditions reflect the scale).
    assert np.isfinite(cit.iloc[0]['Metabolite_Sera'])


# ── KEGG name fallback prefers exact match ──────────────────────
def test_find_by_name_prefers_exact_match(tmp_path):
    resolver = KEGGCompoundResolver(cache_dir=tmp_path)

    fake_find = (
        "C00036\tOxaloacetate\tC4H4O5\n"
        "C00042\tSuccinic acid\tC4H6O4\n"
        "C00424\tSuccinate monoamide\tC4H7NO3\n"
    )

    def fake_get(url: str) -> str:
        if '/find/compound/' in url:
            return fake_find
        # fetch_compound -> a raw KEGG record for the requested cpd
        cid = url.split('cpd:')[-1]
        return (
            f"ENTRY {cid}\n"
            f"NAME  {'Succinic acid' if cid == 'C00042' else 'Some compound'}\n"
            "FORMULA C4H6O4\n"
            "//\n"
        )

    resolver._kegg_get = fake_get

    out = resolver.find_by_name(['Succinic acid'], top_k=10)
    assert out['Succinic acid']['kegg_id'] == 'C00042'
    assert out['Succinic acid']['name'] == 'Succinic acid'


def test_find_by_name_falls_back_to_first_hit(tmp_path):
    resolver = KEGGCompoundResolver(cache_dir=tmp_path)

    def fake_get(url: str) -> str:
        if '/find/compound/' in url:
            return "C00036\tOxaloacetate\tC4H4O5\nC00123\tAmino acid X\n"
        cid = url.split('cpd:')[-1]
        return f"ENTRY {cid}\nNAME  Whatever\n//\n"

    resolver._kegg_get = fake_get
    out = resolver.find_by_name(['Ambigous name'], top_k=10)
    # No exact name match: first candidate (legacy behaviour) is used.
    assert out['Ambigous name']['kegg_id'] == 'C00036'