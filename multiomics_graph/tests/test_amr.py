"""
test_amr.py
Phase 2 validation - AMR knowledge layer:
  - manifest loads for all strains
  - class grouping / locus mapping is correct
  - AMR hyperedges include gene nodes + AMR_mechanism node
  - genome-only markers (no protein) are still represented
  - every manifest locus is a drawn gene node (100% coverage)
"""
from amr.amr import AMRManifest
from graph.hypergraph import BiologicalHypergraphBuilder


def _merged_amr_draw_nodes(man, strain, aligned_loci):
    """Mirror the AMR draw-graph merge shared by main.py / scripts/amr_figures.py."""
    locus_to_amr = man.locus_to_amr(strain)
    amr_class_map = {}
    genome_only = set()
    for locus in locus_to_amr:
        gnode = f"gene_{locus}"
        amr_class_map[gnode] = locus_to_amr[locus]['amr_class']
        if locus not in aligned_loci:
            genome_only.add(gnode)
    return amr_class_map, genome_only


def test_manifest_loads_all_strains():
    man = AMRManifest()
    expected = {'B36', 'MS_14384', 'MS_14386', 'MS_14387'}
    assert expected.issubset(set(man.strains()))


def test_b36_marker_count():
    man = AMRManifest()
    assert man.summary('B36')['n_markers'] == 8
    assert set(man.amr_loci('B36')) == {
        'EW036_RS26095', 'EW036_RS26140', 'EW036_RS26145',
        'EW036_RS26155', 'EW036_RS26210', 'EW036_RS26225',
        'EW036_RS26250', 'EW036_RS26255',
    }


def test_genome_only_marker_has_no_protein():
    # aadA1 in MS_14386 is a frameshifted pseudogene (protein_id null)
    man = AMRManifest()
    aad = [m for m in man.markers('MS_14386') if m['name'] == 'aadA1'][0]
    assert aad['protein_ids'] == []
    assert aad['locus_tags'] == ['EW035_RS23490']


def test_all_manifest_loci_kept_even_genome_only():
    # No manifest locus may be dropped just because it lacks RNA/protein
    # quantification: every one must become a drawn gene node.
    man = AMRManifest()
    # Empty aligned set -> every marker is a genome-only determinant and all
    # must still be represented as gene nodes.
    amr_class_map, genome_only = _merged_amr_draw_nodes(man, 'B36', set())
    assert len(amr_class_map) == 8
    assert len(genome_only) == 8


def test_genome_only_classification_matches_report_level():
    # Genes absent from aligned_multiomics are classified genome_only.
    man = AMRManifest()
    # B36: tet(A) EW036_RS26095 is genome-only (no omics quantification).
    aligned = set(man.amr_loci('B36')) - {'EW036_RS26095'}
    amr_class_map, genome_only = _merged_amr_draw_nodes(man, 'B36', aligned)
    assert 'gene_EW036_RS26095' in genome_only
    assert 'gene_EW036_RS26140' not in genome_only


def test_amr_hyperedges_include_mechanism_node():
    man = AMRManifest()
    hb = BiologicalHypergraphBuilder()
    hb.add_amr_hyperedges(
        man.class_to_loci('B36'),
        man.amr_loci('B36'),
        gene_to_protein={'EW036_RS26225': 'WP_000000239.1'},
    )
    amr = [h for h in hb.hyperedges if h['type'] == 'amr_mechanism']
    assert len(amr) == 7  # 7 distinct classes in B36

    beta = [h for h in amr if h['name'] == 'beta_lactam_resistance'][0]
    # Contains both beta-lactam genes + the mechanism node.
    assert 'EW036_RS26225' in beta['nodes']
    assert 'EW036_RS26250' in beta['nodes']
    assert 'AMR:beta_lactam' in beta['nodes']
    # Genome-only class (tetracycline) is still represented via mechanism node.
    tet = [h for h in amr if h['name'] == 'tetracycline_resistance'][0]
    assert 'EW036_RS26095' in tet['nodes']
    assert 'AMR:tetracycline' in tet['nodes']