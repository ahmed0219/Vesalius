"""
test_amr_path_analysis.py
Validation for scripts/amr_path_analysis.py core logic:
  - layer coverage scoring over synthetic multi-omics edge data
  - RNA/protein concordance classification
  - genome-only vs multi-omics status assignment
  - connected-metabolite detection via shared pathways
  - evidence-level classification of relationships
  - explicit layer status (observed / not_observed / unavailable)
  - candidate AMR-associated element aggregation (excludes known AMR loci)
"""
import pandas as pd
from pathlib import Path

from scripts.amr_path_analysis import (
    LAYERS,
    EVIDENCE_ANNOTATION,
    EVIDENCE_DIRECT,
    EVIDENCE_GRAPH,
    EVIDENCE_INFERRED,
    _connected_metabolites,
    _edge_map,
    _hyperedge_memberships,
    _quantified_sets,
    _strain_layer_availability,
    _wp_edge_map,
    build_candidates,
    build_evidence_rows,
    trace_marker,
)
from scripts.enrich_amr_graph import _update_hyperedges


def _sample_data():
    aligned = pd.DataFrame([
        {'GeneID': 'X_00001', 'RNA_RPMI': 100.0, 'RNA_Sera': 50.0,
         'RNA_logFC': -50.0, 'ProteinID': 'WP_0001.1',
         'Protein_RPMI': 10.0, 'Protein_Sera': 12.0,
         'Protein_logFC': 2.0, 'Protein_Regulation': 'Up'},
        {'GeneID': 'X_00002', 'RNA_RPMI': 1.0, 'RNA_Sera': 1.0,
         'RNA_logFC': 0.0, 'ProteinID': 'WP_0002.1',
         'Protein_RPMI': 1.0, 'Protein_Sera': 1.0,
         'Protein_logFC': 0.0, 'Protein_Regulation': 'Stable'},
    ])
    return {
        'aligned': aligned,
        'genome_genes': pd.DataFrame({'GeneID': ['X_00001', 'X_00002']}),
        'gene_features': pd.DataFrame({'GeneID': ['X_00001', 'X_00002']}),
        'metabolomics': pd.DataFrame([
            {'MetaboliteID': 'CHEBI:1', 'MetaboliteName': 'TestM',
             'Metabolite_logFC': 1.0},
        ]),
        'in_pathway': pd.DataFrame([
            {'gene_id': 'X_00001', 'annotation_id': 'eco00010'},
        ]),
        'cog': pd.DataFrame([
            {'gene_id': 'X_00001', 'annotation_id': 'V'},
        ]),
        'go': pd.DataFrame([
            {'protein_id': 'WP_0001.1', 'annotation_id': 'GO:0008800'},
        ]),
        'ppi': pd.DataFrame([
            {'protein_id': 'WP_0001.1', 'protein_target': 'WP_9999.1'},
        ]),
        'clusters': pd.DataFrame([
            {'protein_id': 'WP_0001.1', 'annotation_id': 'Cluster_0'},
        ]),
        'amr_edges': pd.DataFrame([
            {'gene_id': 'X_00001', 'annotation_id': 'beta_lactam'},
        ]),
        'genomic': pd.DataFrame([
            {'gene_id': 'X_00001', 'gene_target': 'X_00002'},
        ]),
        'metab_in_pathway': pd.DataFrame([
            {'metabolite_id': 'CHEBI:1', 'annotation_id': 'eco00010'},
        ]),
        'metab_enzymes': pd.DataFrame([
            {'protein_id': 'WP_0001.1', 'metabolite_id': 'CHEBI:1'},
        ]),
        'hyperedges': pd.DataFrame([
            {'hyperedge_id': 0, 'type': 'multi_omics_triple',
             'name': 'Triple_X_00001', 'n_nodes': 3,
             'nodes': 'X_00001;X_00001_RNA;WP_0001.1'},
            {'hyperedge_id': 1, 'type': 'amr_mechanism',
             'name': 'beta_lactam_resistance', 'n_nodes': 2,
             'nodes': 'X_00001;AMR:beta_lactam'},
        ]),
    }


def _marker():
    return {
        'name': 'blaTest', 'amr_class': 'beta_lactam',
        'locus_tags': ['X_00001'], 'protein_ids': ['WP_0001.1'],
        'source': 'test', 'note': '',
    }


def test_trace_marker_coverage_all_layers():
    data = _sample_data()
    data['in_pathway_map'] = {'X_00001': {'eco00010'}}
    data['cog_map'] = {'X_00001': {'V'}}
    data['go_map'] = {'WP_0001.1': {'GO:0008800'}}
    data['clusters_map'] = {'WP_0001.1': {'Cluster_0'}}
    data['amr_edges_map'] = {'X_00001': {'beta_lactam'}}
    data['genomic_map'] = {'X_00001': {'X_00002'}}
    rec = trace_marker('STRAIN', _marker(), data)
    assert rec['layer_coverage'] == len(LAYERS) == 8
    assert rec['status'] == 'multi_omics'
    assert rec['L5_kegg_pathways'] == 1
    assert rec['L6_cog_categories'] == 'V (Defense mechanisms)'
    assert rec['L7_go_terms'] == 1
    assert rec['L8_ppi_neighbors'] == 1
    assert rec['L4_metabolites'] == 1


def test_concordance_and_genome_only():
    data = _sample_data()
    data['in_pathway_map'] = {}
    data['cog_map'] = {}
    data['go_map'] = {}
    data['clusters_map'] = {}
    data['amr_edges_map'] = {}
    data['genomic_map'] = {}
    marker = _marker()
    rec = trace_marker('STRAIN', marker, data)
    # RNA down (log2(Sera/RPMI) = -1) vs protein up (+2) -> discordant
    assert rec['rna_protein_concordance'] == 'discordant'

    # A marker in the genome but with no omics quantitation -> genome_only
    marker['locus_tags'] = ['X_00099']
    marker['protein_ids'] = []
    data['genome_genes'] = pd.DataFrame({'GeneID': ['X_00001', 'X_00002', 'X_00099']})
    data['gene_features'] = pd.DataFrame({'GeneID': ['X_00001', 'X_00002', 'X_00099']})
    rec2 = trace_marker('STRAIN', marker, data)
    assert rec2['status'] == 'genome_only'
    assert rec2['layer_coverage'] == 1  # genome layer only


def _prepared_data():
    data = _sample_data()
    data['in_pathway_map'] = _edge_map(data['in_pathway'], 'gene_id', 'annotation_id')
    data['cog_map'] = _edge_map(data['cog'], 'gene_id', 'annotation_id')
    data['go_map'] = _wp_edge_map(data['go'], 'protein_id', 'annotation_id')
    data['clusters_map'] = _wp_edge_map(data['clusters'], 'protein_id', 'annotation_id')
    data['amr_edges_map'] = _edge_map(data['amr_edges'], 'gene_id', 'annotation_id')
    data['genomic_map'] = _edge_map(data['genomic'], 'gene_id', 'gene_target')
    data['_layer_available'] = _strain_layer_availability(data)
    data['_quantified'] = _quantified_sets(data)
    data['_amr_loci'] = {'X_00001'}
    data['_amr_proteins'] = {'WP_0001.1'}
    return data


def test_evidence_classification():
    data = _prepared_data()
    rec = trace_marker('STRAIN', _marker(), data)
    rec['hyperedges'] = _hyperedge_memberships(data, rec['locus_tag'],
                                               rec['protein_id'] or None)
    rows = build_evidence_rows(rec, data)
    by_rel = {r['relationship']: r for r in rows}
    assert by_rel['rna_quantified']['evidence_type'] == EVIDENCE_DIRECT
    assert by_rel['protein_quantified']['evidence_type'] == EVIDENCE_DIRECT
    assert by_rel['genomic_proximity']['evidence_type'] == EVIDENCE_GRAPH
    assert by_rel['ppi']['evidence_type'] == EVIDENCE_GRAPH
    assert by_rel['in_pathway']['evidence_type'] == EVIDENCE_ANNOTATION
    assert by_rel['go_annotation']['evidence_type'] == EVIDENCE_ANNOTATION
    assert by_rel['cog_category']['evidence_type'] == EVIDENCE_ANNOTATION
    assert by_rel['shares_pathway']['evidence_type'] == EVIDENCE_INFERRED


def test_layer_status_unavailable():
    data = _prepared_data()
    data['aligned'] = None
    data['_layer_available'] = _strain_layer_availability(data)
    data['_quantified'] = _quantified_sets(data)
    rec = trace_marker('STRAIN', _marker(), data)
    # transcriptome/proteome data do not exist for this strain -> unavailable
    assert rec['L2_transcriptome_status'] == 'unavailable'
    assert rec['L3_proteome_status'] == 'unavailable'
    assert rec['n_missing_layers'] >= 2
    assert rec['status'] == 'genome_only'


def test_candidates_exclude_amr_and_aggregate():
    def mk_ev(strain, det, tgt, plen):
        return {'strain': strain, 'amr_determinant': det, 'amr_class': 'beta_lactam',
                'source_entity': 'X_00001', 'source_layer': 'genome',
                'target_entity': tgt, 'target_layer': 'genome',
                'relationship': 'genomic_proximity', 'path_length': plen,
                'node_sequence': '', 'edge_sequence': '',
                'evidence_type': EVIDENCE_GRAPH, 'evidence_source': 'test',
                'observed_status': EVIDENCE_GRAPH}

    strain_datas = {
        'S1': {'_amr_loci': {'X_00001', 'X_00004'}, '_amr_proteins': set(),
               '_quantified': {}},
        'S2': {'_amr_loci': {'X_00001'}, '_amr_proteins': set(),
               '_quantified': {}},
    }
    evidence = {
        'S1': [mk_ev('S1', 'bla1', 'X_00002', 1),   # non-AMR neighbor
               mk_ev('S1', 'bla1', 'X_00004', 1)],  # AMR locus -> excluded
        'S2': [mk_ev('S2', 'bla2', 'X_00002', 1)],  # same neighbor, 2nd strain
    }
    cand = build_candidates(strain_datas, evidence)
    assert len(cand) == 1
    row = cand.iloc[0]
    assert row['candidate'] == 'X_00002'
    assert row['number_of_amr_determinants'] == 2
    assert row['number_of_strains'] == 2
    assert row['repeated_across_strains'] == True


def test_candidates_debias_hubs_behind_specific():
    """Candidates reaching >= HUB_FRACTION of the max determinant count are
    flagged hub_candidate and ordered after specific (non-hub) candidates."""
    def mk_ev(strain, det, tgt, plen, tl='genome'):
        return {'strain': strain, 'amr_determinant': det, 'amr_class': 'beta_lactam',
                'source_entity': 'X_00001', 'source_layer': 'genome',
                'target_entity': tgt, 'target_layer': tl,
                'relationship': 'genomic_proximity', 'path_length': plen,
                'node_sequence': '', 'edge_sequence': '',
                'evidence_type': EVIDENCE_GRAPH, 'evidence_source': 'test',
                'observed_status': EVIDENCE_GRAPH}

    strain_datas = {
        'S1': {'_amr_loci': {'X_00001'}, '_amr_proteins': set(),
               '_quantified': {}},
        'S2': {'_amr_loci': {'X_00001'}, '_amr_proteins': set(),
               '_quantified': {}},
        'S3': {'_amr_loci': {'X_00001'}, '_amr_proteins': set(),
               '_quantified': {}},
    }
    # 4 distinct determinants; the hub candidate reaches all 4, the specific
    # candidate only 1.
    hub_evs = []
    for det in ('bla1', 'bla2', 'bla3', 'bla4'):
        for s in ('S1', 'S2', 'S3'):
            hub_evs.append(mk_ev(s, det, 'CHEBI:HUB', 2, tl='metabolome'))
    specific_evs = [mk_ev('S1', 'bla1', 'X_00002', 1)]
    evidence = {
        'S1': hub_evs + [mk_ev('S1', 'bla1', 'X_00002', 1)],
        'S2': hub_evs,
        'S3': hub_evs,
    }
    cand = build_candidates(strain_datas, evidence)
    assert set(cand['candidate']) == {'CHEBI:HUB', 'X_00002'}
    hub_row = cand[cand['candidate'] == 'CHEBI:HUB'].iloc[0]
    spec_row = cand[cand['candidate'] == 'X_00002'].iloc[0]
    assert hub_row['hub_candidate'] == True
    assert spec_row['hub_candidate'] == False
    assert spec_row['rank'] < hub_row['rank']
    # candidate_name resolves only for metabolites
    assert 'candidate_name' in cand.columns


def _tmp_hyperedge_df():
    return pd.DataFrame([
        {'hyperedge_id': 10, 'type': 'kegg_pathway', 'name': 'ko01100',
         'nodes': 'GENE_A;GENE_B', 'n_nodes': 2},
        {'hyperedge_id': 11, 'type': 'go_term', 'name': 'GO:0001: x',
         'nodes': 'WP_0001.1', 'n_nodes': 1},
    ])


def test_update_hyperedges_bumps_existing(tmp_path: Path):
    """A hyperedge already present must gain the new locus node (int keys).
    Regression: nodes_map was keyed by string ids while bump() looked up by
    int, silently never updating existing hyperedges."""
    he = _tmp_hyperedge_df()
    loci = ['X_00009']
    gene_to_protein = {'X_00009': 'WP_0002.1'}
    eggnog = {'WP_0002': {'cog': [], 'kegg': ['ko01100']}}
    added, updated = _update_hyperedges(he, tmp_path, loci,
                                        gene_to_protein, eggnog, {})
    assert added == 0  # ko01100 already exists -> no new hyperedge
    assert updated >= 1  # but the existing hyperedge was bumped
    # the existing kegg_pathway hyperedge now contains the new locus
    out = pd.read_csv(tmp_path / 'hyperedges.csv', dtype=str)
    row = out[out['hyperedge_id'] == '10'].iloc[0]
    assert 'X_00009' in str(row['nodes'])
