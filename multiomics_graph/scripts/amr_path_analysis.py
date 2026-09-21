"""
scripts/amr_path_analysis.py
Multi-omics molecular-context tracing for every curated AMR determinant.

Scope
-----
This is an EXPLORATORY / mechanism-characterization framework. It does NOT
train a model, does NOT classify or predict AMR, and does NOT claim causal
relationships. For each known AMR determinant it describes:

  - which of the 8 multi-omics graph layers contain information about it,
  - what evidence supports each relationship (direct measurement, graph edge,
    annotation, inference), keeping observed vs inferred information distinct,
  - bounded (<=4 hop) paths through the biological entity graph,
  - per-strain differences and conserved context, and
  - candidate AMR-associated elements (non-AMR entities that repeatedly occur
    in the molecular context of known determinants) -- hypotheses, not
    confirmed resistance determinants.

The script reads ONLY saved pipeline outputs under `outputs/<strain>/` and
`outputs/<strain>/graph/`. It never reruns the upstream pipeline.

Evidence-level taxonomy
-----------------------
  directly_observed     measured abundance (RNA/protein/metabolite) or a
                        directly represented genomic locus
  graph_derived         graph connectivity: genomic-proximity edge, PPI edge,
                        cluster membership, hyperedge membership, bounded
                        multi-hop path
  annotation_mediated   KEGG pathway membership, COG assignment, GO annotation,
                        curated AMR mechanism
  inferred              metabolite relationship mediated by pathway/enzyme
                        information, or multi-hop context reaching a metabolite
  unavailable           the omics layer itself does not exist for the strain

Missing data are never fabricated: a missing measurement is reported as
'not_observed' and a missing layer as 'unavailable', never as a measured zero.

Outputs (written under `reports/`):
  amr_path_trace_<strain>.csv        one row per AMR determinant (all layers)
  amr_path_evidence_<strain>.csv     one row per determinant x evidence item
  amr_path_hyperedge_<strain>.csv    hyperedge memberships per determinant
  amr_path_cross_strain.csv          determinant x strain comparison matrix
  amr_path_candidates.csv            cross-strain candidate AMR-associated elements
  amr_path_report.md                 biological interpretation report

Usage (from multiomics_graph/):
    python scripts/amr_path_analysis.py
    python scripts/amr_path_analysis.py --strains B36 MS_14386
    python scripts/amr_path_analysis.py --output reports
"""

import argparse
import math
import re
import sys
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from amr.amr import AMRManifest

REPO = Path(__file__).resolve().parent.parent
STRAINS_OUT = REPO / 'outputs'
MANIFEST = AMRManifest()

STRAIN_ORDER = ['B36', 'MS_14384', 'MS_14386', 'MS_14387']

# Bounded path strategy: multi-hop context is limited to <=4 hops and each
# marker contributes at most MAX_BFS_ROWS path records (shortest paths first).
MAX_PATH_DEPTH = 4
MAX_BFS_ROWS = 400
TOP_CANDIDATES = 20
# Statistical enrichment of candidates is reported alongside the descriptive
# ranking so the "top candidates" table distinguishes hypotheses with a real
# over-representation in AMR context from chance associations. These stats
# never change the descriptive ordering; they add p-value / FDR columns.
ENRICHMENT_ENABLED = True

LAYERS = [
    'genome', 'transcriptome', 'proteome', 'metabolome',
    'kegg', 'cog', 'go', 'ppi_amr',
]

COG_DESC = {
    'J': 'Translation', 'A': 'RNA processing', 'K': 'Transcription',
    'L': 'Replication/repair', 'B': 'Chromatin', 'D': 'Cell division',
    'Y': 'Nuclear structure', 'V': 'Defense mechanisms', 'T': 'Signal transduction',
    'M': 'Cell wall/membrane', 'N': 'Cell motility', 'Z': 'Cytoskeleton',
    'W': 'Extracellular', 'U': 'Intracellular trafficking', 'O': 'Posttranslational mod.',
    'X': 'Mobilome/prophage', 'C': 'Energy production', 'G': 'Carbohydrate metabolism',
    'E': 'Amino acid metabolism', 'F': 'Nucleotide metabolism', 'H': 'Coenzyme metabolism',
    'I': 'Lipid metabolism', 'P': 'Inorganic ion transport', 'Q': 'Secondary metabolites',
    'R': 'General function prediction', 'S': 'Unknown function',
}

# Evidence-level taxonomy (see module docstring)
EVIDENCE_DIRECT = 'directly_observed'
EVIDENCE_GRAPH = 'graph_derived'
EVIDENCE_ANNOTATION = 'annotation_mediated'
EVIDENCE_INFERRED = 'inferred'
EVIDENCE_UNAVAILABLE = 'unavailable'
EVIDENCE_ORDER = [EVIDENCE_DIRECT, EVIDENCE_GRAPH, EVIDENCE_ANNOTATION,
                  EVIDENCE_INFERRED, EVIDENCE_UNAVAILABLE]

# Default evidence level attributed to a *present* layer at the marker level.
LAYER_EVIDENCE = {
    'genome': EVIDENCE_DIRECT,
    'transcriptome': EVIDENCE_DIRECT,
    'proteome': EVIDENCE_DIRECT,
    'metabolome': EVIDENCE_INFERRED,
    'kegg': EVIDENCE_ANNOTATION,
    'cog': EVIDENCE_ANNOTATION,
    'go': EVIDENCE_ANNOTATION,
    'ppi_amr': EVIDENCE_GRAPH,
}

_WP_RE = re.compile(r'^WP_\d+(\.\d+)?$')

# A candidate whose number of connected determinants reaches this fraction of
# the MAXIMUM determinant count reached by any candidate is treated as a
# generic "hub" entity (e.g. a ubiquitous pathway-shared metabolite); hubs
# carry little discriminative information and are down-ranked below specific
# candidates. Used only to order the hypothesis list, not as a scored filter.
HUB_FRACTION = 0.75


def _chebi_names() -> dict:
    """Metabolite id -> name from the saved metabolomics tables (read-only)."""
    names = {}
    for sd in STRAINS_OUT.iterdir():
        mtab = sd / 'metabolomics_abundance.csv'
        if not mtab.exists():
            continue
        try:
            df = pd.read_csv(mtab, low_memory=False)
        except Exception:
            continue
        idcol = ('MetaboliteID' if 'MetaboliteID' in df.columns
                 else df.columns[0])
        ncol = 'MetaboliteName' if 'MetaboliteName' in df.columns else None
        for _, r in df.iterrows():
            mid = str(r[idcol]).strip()
            if not mid:
                continue
            if ncol and pd.notna(r[ncol]) and str(r[ncol]).strip():
                names.setdefault(mid, str(r[ncol]).strip())
    return names


CHEBI_NAMES = _chebi_names()


def _norm_wp(pid):
    """Return the unversioned WP identifier (WP_000002542.1 -> WP_000002542).

    Some eggNOG/COG resources use the unversioned form while the graph uses the
    versioned form. Normalizing to a single canonical key prevents orphan or
    duplicate nodes when matching annotations to proteins.
    """
    if not isinstance(pid, str):
        return ''
    s = pid.strip()
    m = re.match(r'^WP_(\d+)(?:\.\d+)?$', s)
    return f'WP_{m.group(1)}' if m else s


def _entity_names() -> dict:
    """Entity id -> human-readable name for gene/protein candidates.

    Gene symbols and product descriptions come from each strain's
    `genome_genes.csv`; proteins are mapped back to their encoding gene so a
    WP accession reports the gene product name.
    """
    names = {}
    gene_of_protein = {}
    for sd in STRAINS_OUT.iterdir():
        gg = sd / 'genome_genes.csv'
        if not gg.exists():
            continue
        try:
            df = pd.read_csv(gg, low_memory=False)
        except Exception:
            continue
        gcol = 'GeneID' if 'GeneID' in df.columns else None
        if not gcol:
            continue
        for _, r in df.iterrows():
            gid = str(r[gcol]).strip()
            if not gid:
                continue
            sym = str(r.get('GeneName', '') or '').strip()
            prod = str(r.get('Product', '') or '').strip()
            label = sym if sym else prod
            if label:
                names.setdefault(gid, label)
            pid = str(r.get('ProteinID', '') or '').strip()
            if pid and label:
                names.setdefault(pid, label)
                gene_of_protein.setdefault(_norm_wp(pid), gid)
    for pid, _g in gene_of_protein.items():
        if pid not in names and _g in names:
            names[pid] = names[_g]
    return names


ENTITY_NAMES = _entity_names()


def _candidate_name(entity: str, etype: str) -> str:
    if etype == 'metabolite':
        return CHEBI_NAMES.get(entity, '')
    return ENTITY_NAMES.get(entity, '') or ENTITY_NAMES.get(_norm_wp(entity), '')


def _paths_from_cache() -> dict:
    """KEGG pathway id -> name, from the saved pipeline cache (read-only)."""
    cache = STRAINS_OUT / 'cache' / 'kegg_list_pathways_eco.json'
    if not cache.exists():
        return {}
    try:
        raw = cache.read_text(encoding='utf-8')
        return {k: (v.get('name', k) if isinstance(v, dict) else k)
                for k, v in pd.read_json(raw).to_dict().items()}
    except Exception:
        return {}


PATHWAY_NAMES = _paths_from_cache()


# --------------------------------------------------------------------------
# Low-level readers
# --------------------------------------------------------------------------

def _read(path: Path):
    return pd.read_csv(path, low_memory=False) if path.exists() else None


def _read_edges(path: Path, issues: list, req_cols=None):
    """Read an edge CSV and record data-correctness problems.

    Duplicate columns (a known historical bug in some same-type edge files)
    are detected here so corrupted tables are never silently trusted.
    """
    if not path.exists():
        return None
    df = pd.read_csv(path, low_memory=False)
    dup = sorted({c for c in df.columns if list(df.columns).count(c) > 1})
    if dup:
        issues.append(f'{path.name}: duplicate columns {dup} -> not trusted')
    if req_cols:
        missing = [c for c in req_cols if c not in df.columns]
        if missing:
            issues.append(f'{path.name}: missing required columns {missing}')
    return df


def _collect(df, col):
    if df is None or col not in df.columns:
        return set()
    return set(df[col].dropna().astype(str))


# --------------------------------------------------------------------------
# Recovered UniProt annotations (per-ID lookup, bypasses batch truncation)
# --------------------------------------------------------------------------

RECOVERY_PATH = REPO / 'reports' / 'amr_uniprot_recovery.csv'


def _load_uniprot_recovery() -> dict:
    """Load per-ID recovered UniProt annotations for AMR determinants.

    Produced by `scripts/recover_amr_uniprot.py` because the upstream
    UniProt batch search (`_search_refseq_batch`) caps results at 500 entries
    per request, silently dropping later IDs in each ~1000-ID batch. A single
    WP query is not truncated, so the canonical UniProt entry (and its GO
    terms) can be recovered. Keyed by both versioned and unversioned WP ids.
    """
    df = _read(RECOVERY_PATH)
    if df is None or df.empty or 'protein_id' not in df.columns:
        return {}
    out = {}
    for _, r in df.iterrows():
        pid = str(r['protein_id']).strip()
        go = {g.strip() for g in str(r.get('go_ids', '')).split(';') if g.strip()}
        rec = {
            'go': go,
            'uniprot_id': str(r.get('uniprot_id', '')),
            'protein_name': str(r.get('protein_name', '')),
            'function': str(r.get('function', '')),
            'recovered': bool(r.get('recovered', False)),
        }
        out[pid] = rec
        out[_norm_wp(pid)] = rec
    return out


UNIPROT_RECOVERY = _load_uniprot_recovery()


def _recovered_for(protein: str, data: dict) -> dict:
    """Recovered UniProt record for a protein id (versioned or unversioned)."""
    if not protein:
        return {}
    rec = (data.get('_uniprot_recovery') or {}).get(protein)
    if rec is None:
        rec = (data.get('_uniprot_recovery') or {}).get(_norm_wp(protein))
    return rec or {}


# --------------------------------------------------------------------------
# Layer availability / status
# --------------------------------------------------------------------------

def _strain_layer_availability(data: dict) -> dict:
    return {
        'genome': data.get('genome_genes') is not None,
        'transcriptome': (data.get('aligned') is not None
                          and 'RNA_RPMI' in data['aligned'].columns),
        'proteome': (data.get('aligned') is not None
                     and 'ProteinID' in data['aligned'].columns),
        'metabolome': data.get('metabolomics') is not None,
        'kegg': data.get('in_pathway') is not None,
        'cog': data.get('cog') is not None,
        'go': data.get('go') is not None,
        'ppi_amr': (data.get('ppi') is not None
                    or data.get('clusters') is not None
                    or data.get('amr_edges') is not None),
    }


def _layer_status(available: bool, present: bool) -> str:
    """observed / not_observed / unavailable -- never a fabricated zero."""
    if not available:
        return 'unavailable'
    return 'observed' if present else 'not_observed'


# --------------------------------------------------------------------------
# Strain data loading
# --------------------------------------------------------------------------

def _load_strain(name: str):
    """Load saved pipeline outputs for a strain.

    Returns (data, issues) where issues is a list of data-correctness warnings
    that must be surfaced in the report.
    """
    sd = STRAINS_OUT / name
    if not sd.exists():
        return {}, []
    issues = []
    gd = sd / 'graph'
    data = {
        'aligned': _read(sd / 'aligned_multiomics.csv'),
        'transcriptomics': _read(sd / 'transcriptomics_expression.csv'),
        'genome_genes': _read(sd / 'genome_genes.csv'),
        'gene_features': _read(sd / 'gene_features.csv'),
        'metabolomics': _read(sd / 'metabolomics_abundance.csv'),
        'in_pathway': _read_edges(gd / 'edges_gene_in_pathway_annotation.csv',
                                  issues, ['gene_id', 'annotation_id']),
        'cog': _read_edges(gd / 'edges_gene_cog_category_annotation.csv',
                           issues, ['gene_id', 'annotation_id']),
        'go': _read_edges(gd / 'edges_protein_annotated_by_annotation.csv',
                          issues, ['protein_id', 'annotation_id']),
        'ppi': _read_edges(gd / 'edges_protein_ppi_protein.csv',
                           issues, ['protein_id', 'protein_target']),
        'clusters': _read_edges(gd / 'edges_protein_member_of_cluster_annotation.csv',
                                issues, ['protein_id', 'annotation_id']),
        'amr_edges': _read_edges(gd / 'edges_gene_associated_with_annotation.csv',
                                 issues, ['gene_id', 'annotation_id']),
        'genomic': _read_edges(gd / 'edges_gene_genomic_proximity_gene.csv',
                               issues, ['gene_id', 'gene_target']),
        'metab_in_pathway': _read_edges(gd / 'edges_metabolite_in_pathway_annotation.csv',
                                        issues, ['metabolite_id', 'annotation_id']),
        'metab_enzymes': _read_edges(gd / 'edges_metabolite_metabolised_by_protein.csv',
                                     issues, ['protein_id', 'metabolite_id']),
        'encodes': _read_edges(gd / 'edges_gene_encodes_protein.csv',
                               issues, ['gene_id', 'protein_id']),
        'hyperedges': _read(gd / 'hyperedges.csv'),
    }

    # Recovered per-ID UniProt annotations (see _load_uniprot_recovery)
    data['_uniprot_recovery'] = UNIPROT_RECOVERY

    data['_gene_ids'] = (
        _collect(data['aligned'], 'GeneID')
        | _collect(data['genome_genes'], 'GeneID')
        | _collect(data['gene_features'], 'GeneID')
        | _collect(data['encodes'], 'gene_id')
        | _collect(data['in_pathway'], 'gene_id')
        | _collect(data['cog'], 'gene_id')
        | _collect(data['genomic'], 'gene_id')
        | _collect(data['genomic'], 'gene_target')
        | _collect(data['amr_edges'], 'gene_id'))
    data['_protein_ids'] = (
        _collect(data['aligned'], 'ProteinID')
        | _collect(data['encodes'], 'protein_id')
        | _collect(data['ppi'], 'protein_id')
        | _collect(data['ppi'], 'protein_target')
        | _collect(data['go'], 'protein_id')
        | _collect(data['clusters'], 'protein_id')
        | _collect(data['metab_enzymes'], 'protein_id'))
    if data['metabolomics'] is not None:
        idcol = ('MetaboliteID' if 'MetaboliteID' in data['metabolomics'].columns
                 else data['metabolomics'].columns[0])
        data['_metabolite_ids'] = _collect(data['metabolomics'], idcol)
    else:
        data['_metabolite_ids'] = set()
    data['_metabolite_ids'] |= _collect(data['metab_enzymes'], 'metabolite_id')
    data['_metabolite_ids'] |= _collect(data['metab_in_pathway'], 'metabolite_id')

    data['_annotation_ids'] = (
        _collect(data['in_pathway'], 'annotation_id')
        | _collect(data['cog'], 'annotation_id')
        | _collect(data['go'], 'annotation_id')
        | _collect(data['clusters'], 'annotation_id')
        | _collect(data['amr_edges'], 'annotation_id'))

    data['_quantified'] = _quantified_sets(data)
    data['_layer_available'] = _strain_layer_availability(data)
    return data, issues


# --------------------------------------------------------------------------
# Quantitation helpers
# --------------------------------------------------------------------------

def _log2_ratio(numer, denom):
    """log2(numer/denom) for raw abundances; NaN when either is missing/<=0.

    RNA_logFC in the pipeline is stored as a raw-count difference
    (RNA_Sera - RNA_RPMI). This derives the proper log2 fold ratio so it is
    comparable to the log-scale Protein_logFC.
    """
    try:
        n, d = float(numer), float(denom)
        if n <= 0 or d <= 0:
            return float('nan')
        return math.log2(n / d)
    except (TypeError, ValueError):
        return float('nan')


def _quantified_sets(data: dict) -> dict:
    """Sets of genes/proteins/metabolites that have a direct measurement."""
    rna, prot, met = set(), set(), set()
    a = data.get('aligned')
    if a is not None and 'GeneID' in a.columns:
        for col in ('RNA_RPMI', 'RNA_Sera'):
            if col in a.columns:
                rna |= set(a.loc[a[col].notna(), 'GeneID'].astype(str))
        if 'Protein_RPMI' in a.columns:
            prot |= set(a.loc[a['Protein_RPMI'].notna(), 'GeneID'].astype(str))
    if data.get('metabolomics') is not None:
        idcol = ('MetaboliteID' if 'MetaboliteID' in data['metabolomics'].columns
                 else data['metabolomics'].columns[0])
        met |= set(data['metabolomics'][idcol].dropna().astype(str))
    return {'rna': rna, 'protein': prot, 'metabolites': met}


def _rna_row(aligned, locus: str, transcriptomics=None):
    """RNA expression for a locus.

    Looks the locus up in the aligned multi-omics table first (it carries the
    expression values used across the pipeline), then falls back to the raw
    transcriptomics table so genome-only AMR determinants (no RNA measured in
    the aligned set) still get a transcriptome layer instead of a fabricated
    zero.
    """
    if aligned is not None and 'GeneID' in aligned.columns:
        hit = aligned[aligned['GeneID'].astype(str) == locus]
        if not hit.empty:
            r = hit.iloc[0]
            return {
                'RNA_RPMI': r.get('RNA_RPMI'), 'RNA_Sera': r.get('RNA_Sera'),
                'RNA_logFC': r.get('RNA_logFC'),
            }
    if transcriptomics is not None and len(transcriptomics.columns) >= 1:
        idcol = ('GeneID' if 'GeneID' in transcriptomics.columns
                 else transcriptomics.columns[0])
        if idcol:
            hit = transcriptomics[transcriptomics[idcol].astype(str) == locus]
            if not hit.empty:
                r = hit.iloc[0]
                return {
                    'RNA_RPMI': r.get('RNA_RPMI'), 'RNA_Sera': r.get('RNA_Sera'),
                    'RNA_logFC': r.get('RNA_logFC'),
                }
    return {}


def _protein_row(aligned, locus: str):
    if aligned is None or 'GeneID' not in aligned.columns:
        return {}
    hit = aligned[aligned['GeneID'].astype(str) == locus]
    if hit.empty:
        return {}
    r = hit.iloc[0]
    return {
        'ProteinID': r.get('ProteinID'),
        'Protein_RPMI': r.get('Protein_RPMI'), 'Protein_Sera': r.get('Protein_Sera'),
        'Protein_logFC': r.get('Protein_logFC'),
        'Protein_Regulation': r.get('Protein_Regulation'),
    }


def _edge_map(df, src_col, dst_col):
    out = {}
    if df is None:
        return out
    for _, r in df.iterrows():
        s, d = r[src_col], r[dst_col]
        if pd.isna(s) or pd.isna(d):
            continue
        s, d = str(s).strip(), str(d).strip()
        out.setdefault(s, set()).add(d)
    return out


def _wp_edge_map(df, src_col, dst_col):
    """Edge map keyed by both versioned and unversioned WP identifiers.

    Prevents duplicate/orphan biological nodes caused by identifier
    versioning (WP_000002542.1 vs WP_000002542).
    """
    out = {}
    if df is None:
        return out
    for _, r in df.iterrows():
        s, d = r[src_col], r[dst_col]
        if pd.isna(s) or pd.isna(d):
            continue
        s, d = str(s).strip(), str(d).strip()
        out.setdefault(s, set()).add(d)
        out.setdefault(_norm_wp(s), set()).add(d)
    return out


def _go_annotations(go_map, protein):
    """GO/cluster annotations for a protein using versioned+unversioned keys."""
    if not protein or not go_map:
        return []
    ids = set(go_map.get(protein, set()))
    ids |= set(go_map.get(_norm_wp(protein), set()))
    return sorted(ids)


# --------------------------------------------------------------------------
# Metabolite connections
# --------------------------------------------------------------------------

def _metabolite_connections(data: dict, locus: str, protein: str) -> dict:
    """Map each connected metabolite to its connection metadata.

    Connections are INFERRED (mediated by pathway/enzyme information), except
    that 'measured' records whether the metabolite itself has a direct
    abundance measurement.
    """
    conns = {}
    measured = data.get('_quantified', _quantified_sets(data))['metabolites']
    gene_paths = set(data.get('in_pathway_map', {}).get(locus, set()))
    if data.get('metab_in_pathway') is not None:
        for _, r in data['metab_in_pathway'].iterrows():
            mid, pid = r.get('metabolite_id'), r.get('annotation_id')
            if pd.isna(mid) or pd.isna(pid):
                continue
            pid = str(pid).strip()
            if pid in gene_paths:
                mid = str(mid).strip()
                conns.setdefault(mid, {
                    'via': 'pathway', 'path_length': 2, 'pathway': pid,
                    'measured': mid in measured})
    if protein and data.get('metab_enzymes') is not None:
        for _, r in data['metab_enzymes'].iterrows():
            if str(r.get('protein_id', '')) == protein and pd.notna(r.get('metabolite_id')):
                mid = str(r['metabolite_id']).strip()
                conns.setdefault(mid, {
                    'via': 'enzyme', 'path_length': 1, 'pathway': '',
                    'measured': mid in measured})
    return conns


def _connected_metabolites(data: dict, locus: str, protein: str) -> set:
    """Set of connected metabolite ids (compat helper)."""
    return set(_metabolite_connections(data, locus, protein))


def _metabolite_abundance(data, mids):
    out = {}
    if data.get('metabolomics') is None:
        return out
    mtab = data['metabolomics']
    idcol = 'MetaboliteID' if 'MetaboliteID' in mtab.columns else mtab.columns[0]
    for mid in mids:
        hit = mtab[mtab[idcol].astype(str) == mid]
        if not hit.empty:
            r = hit.iloc[0]
            out[mid] = {
                'name': r.get('MetaboliteName', ''),
                'Metabolite_logFC': r.get('Metabolite_logFC'),
            }
    return out


# --------------------------------------------------------------------------
# Bounded path context (<=4 hops over biological entities)
# --------------------------------------------------------------------------

def _build_entity_graph(data: dict) -> dict:
    """Undirected adjacency over biological entities (genes/proteins/metabolites)
    plus KEGG pathway nodes used as connector hubs.

    Direct biological relationships (co-localization, protein interaction,
    encoding, enzyme-substrate) are always traversed. KEGG pathway nodes and
    annotation hyperedges (GO/KEGG/COG terms, PPI clusters) are traversed as
    *annotation-mediated* connector hubs so a determinant can be related to
    co-annotated / co-pathway genes and proteins -- the primary way
    plasmid-borne AMR genes connect to the broader molecular state. The AMR
    mechanism node is a leaf (curated label). Traversal is bounded (<=4 hops,
    capped per determinant).
    """
    adj = {}

    def add(a, b, rel, la, lb):
        if a is None or b is None or pd.isna(a) or pd.isna(b):
            return
        adj.setdefault(str(a), []).append((str(b), rel, la, lb))
        adj.setdefault(str(b), []).append((str(a), rel, lb, la))

    if data.get('encodes') is not None:
        for _, r in data['encodes'].iterrows():
            add(r.get('gene_id'), r.get('protein_id'), 'encodes', 'genome', 'proteome')
    if data.get('genomic') is not None:
        for _, r in data['genomic'].iterrows():
            add(r.get('gene_id'), r.get('gene_target'), 'genomic_proximity', 'genome', 'genome')
    if data.get('ppi') is not None:
        for _, r in data['ppi'].iterrows():
            add(r.get('protein_id'), r.get('protein_target'), 'ppi', 'proteome', 'proteome')
    if data.get('metab_enzymes') is not None:
        for _, r in data['metab_enzymes'].iterrows():
            add(r.get('protein_id'), r.get('metabolite_id'), 'metabolised_by', 'proteome', 'metabolome')
    # KEGG pathway connector hubs (annotation-mediated context)
    if data.get('in_pathway') is not None:
        for _, r in data['in_pathway'].iterrows():
            add(r.get('gene_id'), r.get('annotation_id'), 'in_pathway', 'genome', 'kegg')
    if data.get('metab_in_pathway') is not None:
        for _, r in data['metab_in_pathway'].iterrows():
            add(r.get('metabolite_id'), r.get('annotation_id'), 'in_pathway', 'metabolome', 'kegg')
    # Hyperedge connector hubs (co-annotation / co-cluster context)
    if data.get('hyperedges') is not None:
        metab_ids = data.get('_metabolite_ids', set())
        for _, r in data['hyperedges'].iterrows():
            he_type = str(r.get('type', ''))
            he_name = str(r.get('name', ''))
            if not he_name:
                continue
            rel = f'member_of_{he_type}'
            for n in str(r.get('nodes', '')).split(';'):
                n = n.strip()
                if not n:
                    continue
                if n.startswith('WP_'):
                    add(n, he_name, rel, 'proteome', 'hypergraph')
                elif n in metab_ids:
                    add(n, he_name, rel, 'metabolome', 'hypergraph')
                else:
                    add(n, he_name, rel, 'genome', 'hypergraph')
    return adj


def _bfs(graph: dict, start: str, max_depth: int) -> dict:
    """Shortest-path BFS from start; returns node -> (depth, node_seq, edge_seq)."""
    paths = {start: (0, [start], [])}
    dq = deque([start])
    while dq:
        u = dq.popleft()
        du, seq, eseq = paths[u]
        if du >= max_depth:
            continue
        for v, rel, _la, _lb in graph.get(u, []):
            if v in paths:
                continue
            paths[v] = (du + 1, seq + [v], eseq + [rel])
            dq.append(v)
    return paths


def _bfs_context(data: dict, locus: str, protein: str, max_depth: int) -> list:
    """Multi-hop molecular context (depth 2..max_depth) for one determinant.

    Annotation nodes (pathways/GO/COG/cluster/AMR) are traversal hubs, not
    endpoints: only biological entities (genes/proteins/metabolites) are
    reported. Pathways serve as connectors; the other annotation types are
    leaf evidence.
    """
    graph = data.get('_entity_graph')
    if not graph:
        return []
    gene_ids = data.get('_gene_ids', set())
    protein_ids = data.get('_protein_ids', set())
    metab_ids = data.get('_metabolite_ids', set())
    annotation_ids = data.get('_annotation_ids', set())
    own = {x for x in (locus, protein) if x}
    merged = {}
    for s in own:
        for node, (d, seq, eseq) in _bfs(graph, s, max_depth).items():
            cur = merged.get(node)
            if cur is None or d < cur[0]:
                merged[node] = (d, seq, eseq)
    rows = []
    for node, (d, seq, eseq) in merged.items():
        if d < 2 or node in own or node in annotation_ids:
            continue
        if node in metab_ids:
            layer = 'metabolome'
        elif node in protein_ids:
            layer = 'proteome'
        elif node in gene_ids:
            layer = 'genome'
        else:
            continue
        rows.append({'endpoint': node, 'layer': layer, 'path_length': d,
                     'node_sequence': seq, 'edge_sequence': eseq})
    rows.sort(key=lambda x: (x['path_length'], x['endpoint']))
    return rows[:MAX_BFS_ROWS]


# --------------------------------------------------------------------------
# Hyperedge memberships
# --------------------------------------------------------------------------

def _hyperedge_memberships(data: dict, locus: str, protein: str) -> list:
    memberships = []
    if data.get('hyperedges') is None:
        return memberships
    for _, r in data['hyperedges'].iterrows():
        nodes = set(str(n).strip() for n in str(r.get('nodes', '')).split(';'))
        if locus in nodes or (protein and protein in nodes):
            memberships.append({
                'hyperedge_id': r.get('hyperedge_id'),
                'type': r.get('type'),
                'name': r.get('name'),
                'n_nodes': r.get('n_nodes'),
            })
    return memberships


# --------------------------------------------------------------------------
# Per-marker tracing
# --------------------------------------------------------------------------

def trace_marker(name: str, marker: dict, data: dict) -> dict:
    locus = marker['locus_tags'][0] if marker['locus_tags'] else ''
    protein = marker['protein_ids'][0] if marker['protein_ids'] else None

    available = data.get('_layer_available') or _strain_layer_availability(data)

    # Layer 1 -- Genome
    genome_loci = set(data.get('_gene_ids') or [])
    if not genome_loci:
        for key in ('genome_genes', 'gene_features'):
            if data.get(key) is not None and 'GeneID' in data[key].columns:
                genome_loci |= set(str(g) for g in data[key]['GeneID'].dropna())
    in_genome = locus in genome_loci
    genomic_neighbors = (sorted(data['genomic_map'].get(locus, set()))
                         if data.get('genomic_map') else [])

    # Layer 2 -- Transcriptome
    rna = _rna_row(data.get('aligned'), locus, data.get('transcriptomics'))
    has_rna = (available['transcriptome'] and 'RNA_RPMI' in rna
               and pd.notna(rna['RNA_RPMI']))
    rna_log2 = _log2_ratio(rna.get('RNA_Sera'), rna.get('RNA_RPMI'))
    if math.isnan(rna_log2):
        rna_log2 = None
    if rna_log2 is None:
        rna_dir = ''
    elif rna_log2 > 0:
        rna_dir = 'up'
    elif rna_log2 < 0:
        rna_dir = 'down'
    else:
        rna_dir = 'stable'

    # Layer 3 -- Proteome
    prot = _protein_row(data.get('aligned'), locus)
    if protein is None and prot.get('ProteinID') is not None:
        protein = prot['ProteinID']
    has_protein = (available['proteome'] and 'ProteinID' in prot
                   and pd.notna(prot['ProteinID']))
    if not has_protein:
        protein = None

    # Layer 4 -- Metabolome
    metab_conn = _metabolite_connections(data, locus, protein)
    mids = set(metab_conn)
    metab = _metabolite_abundance(data, mids)

    # Layer 5 -- KEGG
    kegg_ids = sorted(data['in_pathway_map'].get(locus, set()))
    kegg_names = [PATHWAY_NAMES.get(k, k) for k in kegg_ids]

    # Layer 6 -- COG
    cog_ids = sorted(data['cog_map'].get(locus, set()))
    cog_names = [f"{c} ({COG_DESC.get(c, 'Unknown')})" for c in cog_ids]

    # Layer 7 -- GO
    go_ids_graph = _go_annotations(data.get('go_map'), protein)
    recovered = _recovered_for(
        protein or (marker['protein_ids'][0] if marker['protein_ids'] else None),
        data,
    )
    # Recovered GO terms are applied only when the protein is quantified
    # (multi-omics). For genome-only determinants the recovered UniProt
    # identity is still reported, but no proteome-layer GO edges are claimed.
    recovered_go = sorted(set(recovered.get('go', []))) if protein else []
    go_ids = sorted(set(go_ids_graph) | set(recovered_go))

    # Layer 8 -- PPI / clusters / AMR mechanism
    ppi_neighbors = []
    if protein and data.get('ppi') is not None:
        pcols = [c for c in data['ppi'].columns if c != 'relation']
        if len(pcols) >= 2:
            src, dst = pcols[0], pcols[1]
            for _, r in data['ppi'].iterrows():
                if str(r[src]) == protein:
                    ppi_neighbors.append(str(r[dst]))
                elif str(r[dst]) == protein:
                    ppi_neighbors.append(str(r[src]))
    cluster_ids = _go_annotations(data.get('clusters_map'), protein)
    amr_classes = sorted(data['amr_edges_map'].get(locus, set()))

    # Hyperedge memberships (all layers)
    memberships = _hyperedge_memberships(data, locus, protein)
    he_types = sorted({m['type'] for m in memberships})

    # Bounded multi-hop molecular context
    bfs_paths = _bfs_context(data, locus, protein, MAX_PATH_DEPTH)

    # RNA/protein concordance on compatible log scales
    rna_lfc = rna_log2
    prot_lfc = prot.get('Protein_logFC')
    concordance = 'unavailable'
    if rna_lfc is not None and pd.notna(prot_lfc):
        try:
            rl, pl = float(rna_lfc), float(prot_lfc)
            if rl == 0 or pl == 0:
                concordance = 'neutral'
            elif (rl > 0) == (pl > 0):
                concordance = 'concordant'
            else:
                concordance = 'discordant'
        except (TypeError, ValueError):
            concordance = 'unavailable'

    layers_present = {
        'genome': in_genome,
        'transcriptome': has_rna,
        'proteome': has_protein,
        'metabolome': len(mids) > 0,
        'kegg': len(kegg_ids) > 0,
        'cog': len(cog_ids) > 0,
        'go': len(go_ids) > 0,
        'ppi_amr': (len(ppi_neighbors) > 0 or len(cluster_ids) > 0
                    or len(amr_classes) > 0),
    }
    layer_status = {L: _layer_status(available.get(L, False), layers_present[L])
                    for L in LAYERS}

    # Evidence-level counts (machine-readable)
    n_direct = int(has_rna) + int(has_protein) \
        + sum(1 for c in metab_conn.values() if c.get('measured'))
    n_graph = (len(genomic_neighbors) + len(ppi_neighbors)
               + len(cluster_ids) + len(memberships) + len(bfs_paths))
    n_annotation = (len(kegg_ids) + len(cog_ids) + len(go_ids)
                    + len(amr_classes))
    n_inferred = len(metab_conn)
    n_missing = sum(1 for L in LAYERS if not available.get(L, False))
    evidence_types = []
    for ev, n in [(EVIDENCE_DIRECT, n_direct), (EVIDENCE_GRAPH, n_graph),
                  (EVIDENCE_ANNOTATION, n_annotation), (EVIDENCE_INFERRED, n_inferred)]:
        if n > 0:
            evidence_types.append(ev)
    if n_missing:
        evidence_types.append(EVIDENCE_UNAVAILABLE)

    if has_rna or has_protein:
        status = 'multi_omics'
    elif in_genome:
        status = 'genome_only'
    else:
        status = 'unavailable'

    rec = {
        'strain': name,
        'marker': marker['name'],
        'amr_class': marker['amr_class'],
        'locus_tag': locus,
        'protein_id': protein or '',
        # Layer 1
        'L1_genome_present': in_genome,
        'L1_genome_status': layer_status['genome'],
        'L1_genomic_neighbors': len(genomic_neighbors),
        # Layer 2
        'L2_rna_present': has_rna,
        'L2_transcriptome_status': layer_status['transcriptome'],
        'L2_RNA_RPMI': rna.get('RNA_RPMI'),
        'L2_RNA_Sera': rna.get('RNA_Sera'),
        'L2_RNA_logFC_rawdiff': rna.get('RNA_logFC'),
        'L2_RNA_log2_derived': rna_log2,
        'L2_RNA_direction': rna_dir,
        # Layer 3
        'L3_protein_present': has_protein,
        'L3_proteome_status': layer_status['proteome'],
        'L3_Protein_RPMI': prot.get('Protein_RPMI'),
        'L3_Protein_Sera': prot.get('Protein_Sera'),
        'L3_Protein_logFC': prot.get('Protein_logFC'),
        'L3_Protein_Regulation': prot.get('Protein_Regulation', ''),
        # Layer 4
        'L4_metabolites': len(mids),
        'L4_metabolome_status': layer_status['metabolome'],
        'L4_metabolite_names': '; '.join(m['name'] for m in metab.values() if m.get('name')),
        # Layer 5
        'L5_kegg_pathways': len(kegg_ids),
        'L5_kegg_status': layer_status['kegg'],
        'L5_kegg_names': '; '.join(kegg_names),
        # Layer 6
        'L6_cog_categories': '; '.join(cog_names),
        'L6_cog_status': layer_status['cog'],
        # Layer 7
        'L7_go_terms': len(go_ids),
        'L7_go_status': layer_status['go'],
        'L7_go_ids': '; '.join(go_ids),
        'L7_go_recovered_count': len(recovered_go),
        'L7_go_recovered_ids': '; '.join(recovered_go),
        'recovered_uniprot_id': recovered.get('uniprot_id', ''),
        'recovered_protein_name': recovered.get('protein_name', ''),
        'recovered_function': recovered.get('function', ''),
        # Layer 8
        'L8_ppi_neighbors': len(ppi_neighbors),
        'L8_ppi_amr_status': layer_status['ppi_amr'],
        'L8_ppi_clusters': '; '.join(cluster_ids),
        'L8_amr_mechanisms': '; '.join(amr_classes),
        # Cross-layer
        'rna_protein_concordance': concordance,
        'hyperedge_types': '; '.join(he_types),
        'n_hyperedges': len(memberships),
        # Descriptive integration measure (NOT an importance/evidence score)
        'layer_coverage': sum(layers_present.values()),
        'multi_omics_layer_coverage': sum(layers_present.values()),
        'layer_coverage_pct': f"{100.0 * sum(layers_present.values()) / len(LAYERS):.1f}%",
        # Evidence summary
        'evidence_types': '; '.join(evidence_types),
        'n_direct_observations': n_direct,
        'n_graph_derived': n_graph,
        'n_annotation_mediated': n_annotation,
        'n_inferred': n_inferred,
        'n_missing_layers': n_missing,
        'status': status,
        # Internal (excluded from the trace CSV)
        '_genomic_neighbors': genomic_neighbors,
        '_ppi_neighbors': ppi_neighbors,
        '_kegg_ids': kegg_ids,
        '_cog_ids': cog_ids,
        '_go_ids': go_ids_graph,
        '_go_ids_merged': go_ids,
        '_recovered_go': recovered_go,
        '_cluster_ids': cluster_ids,
        '_amr_classes': amr_classes,
        '_metab_connections': metab_conn,
        '_bfs_paths': bfs_paths,
    }
    return rec


# --------------------------------------------------------------------------
# Evidence-level CSV (one row per determinant x evidence/relationship item)
# --------------------------------------------------------------------------

EVIDENCE_COLUMNS = [
    'strain', 'amr_determinant', 'amr_class',
    'source_entity', 'source_layer', 'target_entity', 'target_layer',
    'relationship', 'path_length', 'node_sequence', 'edge_sequence',
    'evidence_type', 'evidence_source', 'observed_status',
]


def build_evidence_rows(rec: dict, data: dict) -> list:
    locus = rec['locus_tag']
    protein = rec['protein_id'] or None
    rows = []

    def add(src, sl, tgt, tl, rel, plen, ev, esrc, nseq=None, eseq=None, obs=None):
        rows.append({
            'strain': rec['strain'], 'amr_determinant': rec['marker'],
            'amr_class': rec['amr_class'],
            'source_entity': src, 'source_layer': sl,
            'target_entity': tgt, 'target_layer': tl,
            'relationship': rel, 'path_length': plen,
            'node_sequence': ';'.join(nseq) if nseq else '',
            'edge_sequence': ';'.join(eseq) if eseq else '',
            'evidence_type': ev, 'evidence_source': esrc,
            'observed_status': obs or ev,
        })

    # Direct quantitative observations
    if rec.get('L2_rna_present'):
        add(locus, 'genome', f'{locus}_RNA', 'transcriptome', 'rna_quantified', 0,
            EVIDENCE_DIRECT, 'aligned_multiomics.csv (RNA_RPMI/RNA_Sera)',
            obs='observed')
    if rec.get('L3_protein_present') and protein:
        add(locus, 'genome', protein, 'proteome', 'protein_quantified', 0,
            EVIDENCE_DIRECT, 'aligned_multiomics.csv (Protein_RPMI/Protein_Sera)',
            obs='observed')

    # Graph-derived: direct relationships
    for nb in rec.get('_genomic_neighbors', []):
        add(locus, 'genome', nb, 'genome', 'genomic_proximity', 1,
            EVIDENCE_GRAPH, 'edges_gene_genomic_proximity_gene.csv',
            [locus, nb], ['genomic_proximity'])
    for nb in rec.get('_ppi_neighbors', []):
        add(protein, 'proteome', nb, 'proteome', 'ppi', 1,
            EVIDENCE_GRAPH, 'edges_protein_ppi_protein.csv',
            [protein, nb], ['ppi'])
    for c in rec.get('_cluster_ids', []):
        add(protein, 'proteome', c, 'ppi_amr', 'member_of_cluster', 1,
            EVIDENCE_GRAPH, 'edges_protein_member_of_cluster_annotation.csv',
            [protein, c], ['member_of_cluster'])
    for he in rec.get('hyperedges', []):
        src = protein or locus
        add(src, 'hypergraph', he['name'], 'hypergraph', 'hyperedge_membership', 1,
            EVIDENCE_GRAPH, 'hyperedges.csv',
            [src, he['name']], [f"member_of:{he['type']}"])

    # Annotation-mediated
    for k in rec.get('_kegg_ids', []):
        add(locus, 'genome', k, 'kegg', 'in_pathway', 1,
            EVIDENCE_ANNOTATION, 'edges_gene_in_pathway_annotation.csv (KEGG)',
            [locus, k], ['in_pathway'])
    for c in rec.get('_cog_ids', []):
        add(locus, 'genome', c, 'cog', 'cog_category', 1,
            EVIDENCE_ANNOTATION, 'edges_gene_cog_category_annotation.csv (eggNOG/COG)',
            [locus, c], ['cog_category'])
    for g in rec.get('_go_ids', []):
        add(protein, 'proteome', g, 'go', 'go_annotation', 1,
            EVIDENCE_ANNOTATION, 'edges_protein_annotated_by_annotation.csv (GO)',
            [protein, g], ['annotated_by'])
    for g in rec.get('_recovered_go', []):
        add(protein, 'proteome', g, 'go', 'go_annotation', 1,
            EVIDENCE_ANNOTATION,
            'amr_uniprot_recovery.csv (per-ID UniProt lookup; '
            'bypasses upstream batch-search truncation)',
            [protein, g], ['annotated_by'])
    for m in rec.get('_amr_classes', []):
        add(locus, 'genome', m, 'ppi_amr', 'amr_mechanism', 1,
            EVIDENCE_ANNOTATION, 'amr_manifest / edges_gene_associated_with_annotation.csv',
            [locus, m], ['associated_with'])

    # Inferred: metabolite relationships (pathway/enzyme mediated)
    for mid, conn in rec.get('_metab_connections', {}).items():
        if conn['via'] == 'pathway':
            add(locus, 'genome', mid, 'metabolome', 'shares_pathway', 2,
                EVIDENCE_INFERRED, 'pathway-mediated (gene->pathway->metabolite)',
                [locus, conn['pathway'], mid], ['in_pathway', 'in_pathway'])
        else:
            add(protein, 'proteome', mid, 'metabolome', 'metabolised_by', 1,
                EVIDENCE_INFERRED, 'edges_metabolite_metabolised_by_protein.csv',
                [protein, mid], ['metabolised_by'])
        if conn.get('measured'):
            add(locus, 'genome', mid, 'metabolome', 'metabolite_quantified',
                conn['path_length'], EVIDENCE_DIRECT, 'metabolomics_abundance.csv',
                obs='observed')

    # Graph/annotation-derived: bounded multi-hop context (depth 2..4)
    known_metabs = set(rec.get('_metab_connections', {}))
    for p in rec.get('_bfs_paths', []):
        if p['layer'] == 'metabolome' and p['endpoint'] in known_metabs:
            continue  # already reported as a direct L4 metabolite connection
        eseq = p['edge_sequence']
        if p['layer'] == 'metabolome':
            ev = EVIDENCE_INFERRED  # metabolite relationship reached via context
        elif any(e == 'member_of_ppi_cluster' for e in eseq):
            ev = EVIDENCE_GRAPH  # co-cluster is network evidence
        elif any(e.startswith('member_of_') for e in eseq) or 'in_pathway' in eseq:
            ev = EVIDENCE_ANNOTATION  # co-annotation / co-pathway association
        else:
            ev = EVIDENCE_GRAPH
        add(locus, 'genome', p['endpoint'], p['layer'], 'graph_context',
            p['path_length'], ev, f'bounded_BFS (<= {MAX_PATH_DEPTH} hops)',
            p['node_sequence'], p['edge_sequence'])
    return rows


# --------------------------------------------------------------------------
# Per-strain CSV writers
# --------------------------------------------------------------------------

def write_trace_csv(recs, out_dir: Path, name: str):
    public = [{k: v for k, v in r.items() if not k.startswith('_')} for r in recs]
    if public:
        pd.DataFrame(public).to_csv(out_dir / f'amr_path_trace_{name}.csv', index=False)


def write_evidence_csv(rows, out_dir: Path, name: str):
    if rows:
        pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).to_csv(
            out_dir / f'amr_path_evidence_{name}.csv', index=False)


def write_hyperedge_detail(name: str, records, out_dir: Path):
    rows = []
    for rec in records:
        for he in rec.get('hyperedges', []):
            rows.append({
                'strain': name,
                'marker': rec['marker'],
                'amr_class': rec['amr_class'],
                'locus_tag': rec['locus_tag'],
                'protein_id': rec['protein_id'],
                'hyperedge_id': he['hyperedge_id'],
                'hyperedge_type': he['type'],
                'hyperedge_name': he['name'],
                'hyperedge_n_nodes': he['n_nodes'],
            })
    if rows:
        pd.DataFrame(rows).to_csv(out_dir / f'amr_path_hyperedge_{name}.csv', index=False)


# --------------------------------------------------------------------------
# Cross-strain matrix
# --------------------------------------------------------------------------

def build_cross_strain(traces: dict) -> pd.DataFrame:
    records = [r for recs in traces.values() for r in recs]
    if not records:
        return pd.DataFrame()
    strains = list(traces.keys())
    all_markers = sorted({r['marker'] for r in records})
    rows = []
    for m in all_markers:
        first = next((r for r in records if r['marker'] == m), None)
        row = {'marker': m, 'amr_class': first['amr_class'] if first else ''}
        for s in strains:
            hits = [r for r in records if r['marker'] == m and r['strain'] == s]
            if not hits:
                row[f'{s}_presence'] = 'absent'
                row[f'{s}_status'] = 'absent'
                for suf in ('coverage', 'coverage_pct', 'n_direct', 'n_graph_derived',
                            'n_annotation_mediated', 'n_inferred', 'n_hyperedges',
                            'n_ppi', 'concordance'):
                    row[f'{s}_{suf}'] = ''
                continue
            r = hits[0]
            row[f'{s}_presence'] = 'present'
            row[f'{s}_status'] = r['status']
            row[f'{s}_coverage'] = r['layer_coverage']
            row[f'{s}_coverage_pct'] = r['layer_coverage_pct']
            row[f'{s}_n_direct'] = r['n_direct_observations']
            row[f'{s}_n_graph_derived'] = r['n_graph_derived']
            row[f'{s}_n_annotation_mediated'] = r['n_annotation_mediated']
            row[f'{s}_n_inferred'] = r['n_inferred']
            row[f'{s}_n_hyperedges'] = r['n_hyperedges']
            row[f'{s}_n_ppi'] = r['L8_ppi_neighbors']
            row[f'{s}_concordance'] = r['rna_protein_concordance']
        rows.append(row)
    cols = ['marker', 'amr_class']
    for s in strains:
        cols += [f'{s}_presence', f'{s}_status', f'{s}_coverage', f'{s}_coverage_pct',
                 f'{s}_n_direct', f'{s}_n_graph_derived', f'{s}_n_annotation_mediated',
                 f'{s}_n_inferred', f'{s}_n_hyperedges', f'{s}_n_ppi', f'{s}_concordance']
    return pd.DataFrame(rows, columns=cols)


# --------------------------------------------------------------------------
# Candidate AMR-associated elements
# --------------------------------------------------------------------------

def _candidate_entity_type(tl: str) -> str:
    return {'genome': 'gene', 'proteome': 'protein',
            'metabolome': 'metabolite'}.get(tl, '')


# --------------------------------------------------------------------------
# Statistical enrichment of candidate elements
# --------------------------------------------------------------------------
#
# Candidates are ranked by the descriptive criteria in build_candidates().
# Independently, each candidate gets a one-sided over-representation test for
# "does this entity connect to MORE AMR determinants than expected by chance,
# given how many entities of the same type exist and how many determinant
# connections that entity type carries overall?".
#
# Null model (per entity type): the S total (entity -> determinant) incidence
# pairs of that type are distributed across the N eligible entities of that
# type, each incidence landing on entity e with probability 1/N. The tail
# P(X >= k_e) is therefore a binomial tail, computed in log-space for numerical
# stability. Hypergeometric (draws without replacement) is a slightly more
# exact choice, but with hundreds of entities the binomial approximation is
# negligible and much easier to interpret/defend.
#
# p-values are corrected across ALL candidates (all entity types pooled) with
# the Benjamini-Hochberg (BH) false-discovery-rate procedure.


def _log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float('-inf')
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _binom_tail_pvalue(k: int, S: int, N: int) -> float:
    """One-sided P(X >= k) for X ~ Binomial(S, 1/N).

    Answers: probability that an entity connects to >= k determinants by
    chance, when the ``S`` determinant-incidences of its type are uniformly
    spread over the ``N`` candidate entities of that type (the entities that
    touch AMR context at least once).

    Uses ``scipy.stats.binom.sf`` for speed and numerical accuracy when scipy
    is available; otherwise falls back to an exact log-space summation (slower
    but dependency-free).
    """
    if k <= 0:
        return 1.0
    if N <= 1:
        return 1.0 if k <= S else 0.0
    if S <= 0:
        return 1.0
    p = 1.0 / N
    try:
        from scipy.stats import binom
        # Survival P(X >= k) = sf(k - 1).
        return float(binom.sf(k - 1, S, p))
    except Exception:
        pass
    # Fallback: exact log-space summation of the tail k..S. Could be slow for
    # very large S; kept only as a dependency-free safety net.
    logp = math.log(p)
    log1mp = math.log1p(-p)
    total = 0.0
    for i in range(k, S + 1):
        log_term = _log_choose(S, i) + i * logp + (S - i) * log1mp
        total += math.exp(log_term)
    return max(0.0, min(1.0, total))


def _bh_fdr(pvalues: list) -> list:
    """Benjamini-Hochberg FDR correction on a list of p-values.

    Returns q-values in the same order as ``pvalues``. NaNs are left as NaN.
    Monotonicity is applied in the standard direction: q-values are forced to
    be non-decreasing as p increases (the cumulative minimum propagates from
    the LARGEST p down to the smallest).
    """
    n = len(pvalues)
    q = [float('nan')] * n
    # Keep only non-NaN indices for ordering.
    valid = [(i, pvalues[i]) for i in range(n)
             if not (isinstance(pvalues[i], float) and math.isnan(pvalues[i]))]
    order = [i for i, _ in sorted(valid, key=lambda t: t[1])]
    m = len(order)
    if m == 0:
        return q
    # Raw BH q = p * n_total_pval / rank (rank over ordered valid p-values).
    # Use n (including NaN) as the denominator is the conservative choice; we
    # use the number of valid tests m to avoid deflating when NaNs are present.
    raw = {i: pvalues[i] * m / (r + 1) for r, i in enumerate(order)}
    # Enforce monotonicity: iterate from the largest p to the smallest.
    prev = float('inf')
    for i in reversed(order):
        v = min(raw[i], prev)
        q[i] = min(1.0, v)
        prev = v
    return q


def build_candidates(strain_datas: dict, evidence_by_strain: dict) -> pd.DataFrame:
    """Non-AMR biological entities recurring in AMR molecular context.

    Candidates are HYPOTHESES, ranked by transparent descriptive criteria:
      1. number of distinct AMR determinants connected (desc)
      2. number of strains (desc)
      3. minimum path length (asc)
      4. number of directly observed items (desc)
    No black-box biological score is used.

    Hub de-biasing: an entity connected to >= HUB_FRACTION of the *maximum*
    number of determinants reached by any candidate (i.e. effectively every
    determinant that can reach it, e.g. a ubiquitous pathway-shared metabolite)
    has low discriminative value, so non-hub candidates are ordered first and
    hubs are appended afterward. Hub status is reported explicitly, not
    silently filtered.
    """
    agg = {}
    for strain, evs in evidence_by_strain.items():
        data = strain_datas.get(strain, {})
        amr_loci = data.get('_amr_loci', set())
        amr_proteins = {_norm_wp(p) for p in data.get('_amr_proteins', set())}
        for ev in evs:
            tl = ev['target_layer']
            etype = _candidate_entity_type(tl)
            if not etype:
                continue  # annotations/hyperedges are not candidate elements
            tgt = ev['target_entity']
            if not tgt or tgt == ev['source_entity']:
                continue
            if tl == 'genome' and tgt in amr_loci:
                continue  # known AMR locus
            if tl == 'proteome' and _norm_wp(tgt) in amr_proteins:
                continue  # protein of a known AMR determinant
            key = (tgt, etype)
            a = agg.setdefault(key, {
                'candidate': tgt, 'entity_type': etype,
                'determinants': set(), 'strains': set(), 'layers': set(),
                'min_path': 99, 'ev_types': set(),
                'n_direct': 0, 'n_graph': 0, 'n_annotation': 0, 'n_inferred': 0,
                'measured_strains': set(),
            })
            a['determinants'].add(f"{strain}|{ev['amr_determinant']}")
            a['strains'].add(strain)
            a['layers'].add(tl)
            if ev['path_length'] is not None and isinstance(ev['path_length'], (int, float)):
                a['min_path'] = min(a['min_path'], int(ev['path_length']))
            a['ev_types'].add(ev['evidence_type'])
            a['n_direct'] += ev['evidence_type'] == EVIDENCE_DIRECT
            a['n_graph'] += ev['evidence_type'] == EVIDENCE_GRAPH
            a['n_annotation'] += ev['evidence_type'] == EVIDENCE_ANNOTATION
            a['n_inferred'] += ev['evidence_type'] == EVIDENCE_INFERRED
            if _directly_measured(data, tgt, etype):
                a['measured_strains'].add(strain)

    rows = []
    for (tgt, etype), a in agg.items():
        n_det = len(a['determinants'])
        rows.append({
            'candidate': tgt,
            'candidate_name': _candidate_name(tgt, etype),
            'entity_type': etype,
            'number_of_amr_determinants': n_det,
            'number_of_strains': len(a['strains']),
            'layers_supported': '; '.join(sorted(a['layers'])),
            'minimum_path_length': a['min_path'],
            'evidence_types': '; '.join(sorted(a['ev_types'])),
            'n_direct_observations': a['n_direct'],
            'n_graph_derived': a['n_graph'],
            'n_annotation_mediated': a['n_annotation'],
            'n_inferred': a['n_inferred'],
            'n_directly_measured_strains': len(a['measured_strains']),
            'strains': '; '.join(sorted(a['strains'])),
            'determinants': '; '.join(sorted({d.split('|', 1)[1] for d in a['determinants']})),
            'repeated_across_strains': len(a['strains']) >= 2,
        })
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    # Hub = connected to ~every determinant that can reach it (uses the max
    # candidate determinant count as the reference, not the raw total which
    # includes determinants with no reachable candidate).
    max_det = int(df['number_of_amr_determinants'].max())
    df['determinant_fraction'] = df['number_of_amr_determinants'] / max_det \
        if max_det else 0.0
    df['hub_candidate'] = df['number_of_amr_determinants'] >= HUB_FRACTION * max_det

    # ── Statistical enrichment (over-representation vs. chance) ──
    # The comparison pool is the OTHER CANDIDATE entities of the SAME type:
    # i.e., among all genes/proteins/metabolites that appear at least once in
    # AMR context, does this entity connect to unusually MANY determinants?
    #
    # Null (per entity type): the ``S`` total (entity -> determinant) incidence
    # pairs of that type are allocated uniformly across the ``N`` candidate
    # entities of that type, so an entity's determinant count is
    # Binomial(S, 1/N). ``p_value`` is the tail P(X >= k_e). This distinguishes
    # a specific, strongly-connected candidate (e.g. an efflux pump) from a
    # generic hub that merely touches everything; the existing `hub_candidate`
    # flag captures the orthogonal "low discriminative value" notion.
    if ENRICHMENT_ENABLED:
        # N = number of candidates of each type; S = total incidence per type.
        S_by_type = {
            t: int(df.loc[df['entity_type'] == t, 'number_of_amr_determinants'].sum())
            for t in ('gene', 'protein', 'metabolite')
        }
        N_by_type = {
            t: int((df['entity_type'] == t).sum())
            for t in ('gene', 'protein', 'metabolite')
        }

        pvals = []
        if len(df):
            # Vectorized binomial tails (one scipy call per entity type).
            k_arr = df['number_of_amr_determinants'].astype(int).values
            etype_arr = df['entity_type'].values
            p_arr = np.full(len(df), np.nan)
            for t in ('gene', 'protein', 'metabolite'):
                sel = etype_arr == t
                S = S_by_type.get(t, 0)
                N = N_by_type.get(t, 0)
                if not sel.any() or N <= 1 or S <= 0:
                    continue
                try:
                    from scipy.stats import binom
                    p_arr[sel] = binom.sf(k_arr[sel] - 1, S, 1.0 / N)
                except Exception:
                    # dependency-free fallback (per element)
                    for i in np.where(sel)[0]:
                        p_arr[i] = _binom_tail_pvalue(int(k_arr[i]), S, N)
            pvals = list(p_arr)
        df['p_value'] = pvals
        df['q_value'] = _bh_fdr(pvals)
        df['enrichment_significant'] = df['q_value'] <= 0.05

    df = df.sort_values(
        ['hub_candidate', 'number_of_amr_determinants', 'number_of_strains',
         'minimum_path_length', 'n_direct_observations'],
        ascending=[True, False, False, True, False]).reset_index(drop=True)
    df.insert(0, 'rank', range(1, len(df) + 1))
    return df


def _directly_measured(data: dict, entity: str, etype: str) -> bool:
    q = data.get('_quantified', {})
    if etype == 'gene':
        return entity in q.get('rna', set()) or entity in q.get('protein', set())
    if etype == 'protein':
        return entity in q.get('protein', set())
    if etype == 'metabolite':
        return entity in q.get('metabolites', set())
    return False


# --------------------------------------------------------------------------
# Data-correctness checks
# --------------------------------------------------------------------------

def _kegg_quality(data: dict, name: str):
    """Detect suspicious/inconsistent KEGG outputs without normalizing them away."""
    df = data.get('in_pathway')
    if df is None:
        return f'{name}: KEGG layer unavailable (no edges_gene_in_pathway_annotation.csv)'
    if 'annotation_id' not in df.columns:
        return f'{name}: KEGG file missing annotation_id column (inconsistent)'
    ann = df['annotation_id'].astype(str)
    n = len(df)
    if n == 0:
        return f'{name}: KEGG file present but empty (inconsistent)'
    distinct = ann.nunique()
    if distinct < 20:
        return (f'{name}: KEGG has only {distinct} distinct pathway annotations -- '
                f'suspiciously degenerate; pathway memberships should be read with caution')
    umbrella = ann[ann.str.startswith(('map0110', 'ko0110'))]
    frac = len(umbrella) / n
    if frac > 0.5:
        return (f'{name}: {frac:.0%} of KEGG rows are global umbrella maps '
                f'(map01100/ko01100); pathway specificity is low')
    return None


def _wp_versioning_note(data: dict, name: str):
    found = []
    for df, col in [(data.get('go'), 'protein_id'),
                    (data.get('clusters'), 'protein_id'),
                    (data.get('ppi'), 'protein_id'),
                    (data.get('encodes'), 'protein_id')]:
        if df is None or col not in df.columns:
            continue
        vals = df[col].dropna().astype(str)
        uv = vals[vals.str.match(r'^WP_\d+$', na=False)]
        if len(uv):
            found.append(f'{col}: {len(uv)} unversioned WP ids')
    if found:
        return (f'{name}: WP versioning -- {"; ".join(found)}; '
                f'annotations were version-normalized to the canonical node')
    return f'{name}: WP versioning -- all protein ids are versioned (WP_....N); no normalization needed'


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

def write_report(traces, cross, candidates, issues_by_strain,
                 strain_datas, out_dir: Path):
    lines = []
    lines.append('# AMR Multi-Omics Path Analysis')
    lines.append('')
    lines.append('> **Scope.** This is an **exploratory, mechanism/context '
                 'characterization** framework. It is **not an AMR classifier** '
                 'and **not a causal inference system**. Graph connectivity does '
                 'not imply causality; annotation-mediated relationships are not '
                 'equivalent to experimental measurements; candidate elements '
                 'are hypotheses, not confirmed resistance determinants.')
    lines.append('')
    lines.append('## A. Dataset / analysis summary')
    lines.append('')
    lines.append(f'- Strains analyzed: {", ".join(traces.keys())}')
    lines.append('- The analysis reads only saved pipeline outputs '
                 '(`outputs/<strain>/` and `outputs/<strain>/graph/`).')
    lines.append(f'- Bounded molecular-context paths: <= {MAX_PATH_DEPTH} hops '
                 f'(shortest first, max {MAX_BFS_ROWS} path records per determinant).')
    lines.append('- Eight layers are traced per determinant; layer coverage is a '
                 '**descriptive integration measure** and is **not** an evidence '
                 'strength or biological importance score.')
    lines.append('')

    lines.append('## B. Per-strain AMR determinant table')
    lines.append('')
    lines.append('| Strain | Determinants | Multi-omics | Genome-only | '
                 'Mean layer coverage |')
    lines.append('|---|---|---|---|---|')
    for s, recs in traces.items():
        n_multi = sum(1 for r in recs if r['status'] == 'multi_omics')
        n_genome = sum(1 for r in recs if r['status'] == 'genome_only')
        mean_cov = (sum(r['layer_coverage'] for r in recs) / len(recs)
                    if recs else 0)
        lines.append(f"| {s} | {len(recs)} | {n_multi} | {n_genome} | "
                     f"{mean_cov:.2f}/8 |")
    lines.append('')

    lines.append('## C. Layer coverage summary')
    lines.append('')
    lines.append('Layer coverage counts how many of the 8 layers contain '
                 'information about a determinant. The layers are heterogeneous '
                 'and not equivalent: a determinant with 6/8 layers is **not** '
                 'necessarily biologically more important than one with 4/8.')
    lines.append('')
    lines.append('| Strain | Determinants | Mean layer coverage |')
    lines.append('|---|---|---|')
    for s, recs in traces.items():
        mean_cov = (sum(r['layer_coverage'] for r in recs) / len(recs)
                    if recs else 0)
        lines.append(f'| {s} | {len(recs)} | {mean_cov:.2f}/8 |')
    lines.append('')

    lines.append('## D. Evidence-level summary')
    lines.append('')
    lines.append('Evidence is classified as **directly observed** (measured '
                 'abundance / represented locus), **graph-derived** (edges, '
                 'hyperedge membership, bounded paths), **annotation-mediated** '
                 '(KEGG/COG/GO/curated mechanism), or **inferred** (metabolite '
                 'relationships via pathway/enzyme information). Missing layers '
                 'are reported as **unavailable**, never fabricated as zeros.')
    lines.append('')
    lines.append('| Strain | Direct | Graph-derived | Annotation-mediated | '
                 'Inferred | Missing layers (marker-level) |')
    lines.append('|---|---|---|---|---|---|')
    for s, recs in traces.items():
        lines.append(f"| {s} | {sum(r['n_direct_observations'] for r in recs)} | "
                     f"{sum(r['n_graph_derived'] for r in recs)} | "
                     f"{sum(r['n_annotation_mediated'] for r in recs)} | "
                     f"{sum(r['n_inferred'] for r in recs)} | "
                     f"{sum(r['n_missing_layers'] for r in recs)} |")
    lines.append('')
    lines.append('**Important:** annotation-mediated evidence (KEGG membership, '
                 'GO annotation, COG assignment) indicates association in a '
                 'reference database, not observed activity in the strain under '
                 'study. STRING PPI edges are database/network evidence, not '
                 'strain-specific experimental validation.')
    lines.append('')

    lines.append('## E. RNA/protein concordance summary')
    lines.append('')
    lines.append('Concordance compares the sign of `RNA log2(Sera/RPMI)` '
                 '(derived) with the sign of `Protein_logFC`; both are log '
                 'scales. It is a descriptive observation, not proof of '
                 'post-transcriptional regulation.')
    lines.append('')
    lines.append('| Strain | Concordant | Discordant | Neutral | Unavailable |')
    lines.append('|---|---|---|---|---|')
    for s, recs in traces.items():
        counts = {v: sum(1 for r in recs
                         if r['rna_protein_concordance'] == v)
                  for v in ('concordant', 'discordant', 'neutral', 'unavailable')}
        lines.append(f"| {s} | {counts['concordant']} | {counts['discordant']} | "
                     f"{counts['neutral']} | {counts['unavailable']} |")
    lines.append('')

    lines.append('## F. Hyperedge summary')
    lines.append('')
    he_counts = {}
    for s, recs in traces.items():
        for r in recs:
            for t in r['hyperedge_types'].split('; '):
                if t:
                    he_counts[t] = he_counts.get(t, 0) + 1
    if he_counts:
        lines.append('| Hyperedge type | Determinant-memberships across strains |')
        lines.append('|---|---|')
        for t, n in sorted(he_counts.items(), key=lambda x: -x[1]):
            lines.append(f'| {t} | {n} |')
    else:
        lines.append('(no hyperedge memberships)')
    lines.append('')
    lines.append('Hyperedge membership is graph-derived context, not evidence '
                 'of functional activity.')
    lines.append('')

    lines.append('## G. Top candidate AMR-associated elements')
    lines.append('')
    lines.append('Candidates are **non-AMR** genes, proteins, or metabolites that '
                 'repeatedly occur in the molecular context of known AMR '
                 'determinants. They are ranked by transparent descriptive '
                 'criteria (number of distinct AMR determinants connected, '
                 'number of strains, minimum path length, direct observations) '
                 'and are **hypotheses for follow-up, not confirmed resistance '
                 'determinants**. Candidates connected to >= '
                 f'{int(HUB_FRACTION * 100)}% of the maximum determinant count '
                 'reached by any candidate are labelled **hub candidates** '
                 '(e.g. ubiquitous pathway-shared metabolites) and are ordered '
                  'after the specific candidates.')
    lines.append('')
    lines.append('Candidates also carry a one-sided **over-representation test**: '
                 'the probability that another entity of the same type would '
                 'connect to >= that many determinants purely by chance '
                 '(binomial tail, `p_value`, corrected across all candidates '
                 'with **Benjamini-Hochberg FDR** into `q_value`). A '
                 '`q_value <= 0.05` marks a candidate as **significantly '
                 'enriched** among the entities of its type that touch AMR '
                 'context. Significance is reported alongside the descriptive '
                 'ranking and does not change the ordering; it is distinct '
                 'from the `hub` flag, which measures non-specificity.')
    lines.append('')
    if candidates is not None and not candidates.empty:
        top = candidates.head(TOP_CANDIDATES)
        lines.append('| Rank | Candidate | Name | Type | #AMR determinants | '
                     '#Strains | Min path | Hub ? | p-value | q-value (FDR) | '
                     'Significant | Evidence types |')
        lines.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
        for _, r in top.iterrows():
            name = r.get('candidate_name', '')
            hub = 'yes' if r.get('hub_candidate') else ''
            pv = f"{r['p_value']:.3g}" if pd.notna(r.get('p_value')) else ''
            qv = f"{r['q_value']:.3g}" if pd.notna(r.get('q_value')) else ''
            sig = 'yes' if r.get('enrichment_significant') else ''
            lines.append(f"| {r['rank']} | {r['candidate']} | {name} | "
                         f"{r['entity_type']} | {r['number_of_amr_determinants']} | "
                         f"{r['number_of_strains']} | {r['minimum_path_length']} | "
                         f"{hub} | {pv} | {qv} | {sig} | {r['evidence_types']} |")
    else:
        lines.append('(no candidate elements identified)')
    lines.append('')

    lines.append('## H. Cross-strain observations')
    lines.append('')
    if cross is not None and not cross.empty:
        multi = cross[cross[[c for c in cross.columns
                              if c.endswith('_presence')]].eq('present').sum(axis=1) >= 2]
        lines.append('Determinants present in >= 2 strains (differences are '
                     'descriptive, not evidence of biological causation):')
        lines.append('')
        lines.append('| Marker | Class | ' + ' | '.join(
            f'{s} coverage' for s in traces.keys()) + ' |')
        lines.append('|' + '---|' * (len(traces.keys()) + 2))
        for _, r in multi.iterrows():
            vals = [str(r.get(f'{s}_coverage', '')) for s in traces.keys()]
            lines.append(f"| {r['marker']} | {r['amr_class']} | " + ' | '.join(vals) + ' |')
    else:
        lines.append('(no cross-strain data)')
    lines.append('')

    lines.append('## I. Missing-data and data-quality notes')
    lines.append('')
    all_issues = []
    for s in traces:
        for L in LAYERS:
            unavail = any(r[f'{_layer_key(L)}'] == 'unavailable' for r in traces[s])
            if unavail:
                all_issues.append(f'{s}: layer {L} unavailable')
    if all_issues:
        lines.append('Unavailable layers (distinct from absence):')
        lines.append('')
        for msg in sorted(set(all_issues)):
            lines.append(f'- {msg}')
        lines.append('')
    else:
        lines.append('No layer is unavailable for any strain; all 8 layers exist '
                     'for every strain in this analysis.')
        lines.append('')
    lines.append('Data-correctness checks:')
    lines.append('')
    for s in traces:
        for msg in issues_by_strain.get(s, []):
            lines.append(f'- **{s}** {msg}')
        q = _kegg_quality(strain_datas.get(s, {}), s)
        if q:
            lines.append(f'- **{s}** {q}')
        else:
            df = strain_datas.get(s, {}).get('in_pathway')
            n_distinct = df['annotation_id'].nunique() if df is not None else 0
            lines.append(f'- **{s}** KEGG consistency: healthy '
                         f'({n_distinct} distinct pathway annotations)')
        lines.append(f'- {_wp_versioning_note(strain_datas.get(s, {}), s)}')
    lines.append('')
    lines.append('Where a saved edge table has duplicate columns, it is flagged '
                 'and not silently trusted. Suspicious KEGG results are reported, '
                 'not silently normalized away.')
    lines.append('')

    # Recovered UniProt annotations (per-ID lookup, bypasses batch truncation)
    n_rec = sum(1 for s in traces for r in traces[s]
                if r.get('L7_go_recovered_count', 0) > 0)
    if n_rec:
        lines.append('Recovered UniProt annotations:')
        lines.append('')
        lines.append(f'- **{n_rec}** determinant rows have GO terms recovered by per-ID '
                     'UniProt lookup (`reports/amr_uniprot_recovery.csv`). The upstream '
                     'pipeline maps RefSeq -> UniProt in ~1000-ID batches but only keeps '
                     'the first 500 matching UniProt entries per request; later IDs in a '
                     'batch are silently dropped, which is why some multi-omics-'
                     'quantified AMR determinants previously had no GO annotation. '
                     'Re-querying each WP accession individually (same search, same '
                     'candidate resolution and parsing) recovers those annotations. The '
                     'recovered terms are real UniProt GO annotations and are labelled '
                     'annotation-mediated with their own evidence source.')
        lines.append('')
    lines.append('Recovered eggNOG KEGG pathways:')
    lines.append('')
    lines.append('- AMR determinants now carry KEGG pathway memberships merged from '
                 'the eggNOG-mapper `KEGG_Pathway` column (regenerated TSVs via '
                 '`scripts/convert_emapper_to_tsv.py`, applied to the saved graphs by '
                 '`scripts/enrich_amr_graph.py`). This closes the previous coverage gap '
                 'where `eco`-bridged UniProt cross-references skipped plasmid-borne '
                 'genes.')
    lines.append('')
    lines.append('Known reference-coverage limitations (KEGG / PPI):')
    lines.append('- **KEGG for acquired AMR determinants.** The `eco`-bridged UniProt '
                 'KEGG cross-references do not cover plasmid-borne genes, but the '
                 'eggNOG-mapper output now carries a `KEGG_Pathway` column '
                 '(see `scripts/convert_emapper_to_tsv.py`), and those memberships are '
                 'merged into the graph (`scripts/enrich_amr_graph.py`). KEGG is still '
                 'genuinely absent where eggNOG assigns no pathway (e.g. many '
                 'efflux/modification determinants with only a `-` KEGG_Pathway).')
    lines.append('- **PPI = 0 for AMR determinants.** STRING is queried against '
                 'E. coli K-12/UPEC reference proteins; the clinical-strain WP '
                 'accessions are not hosted by STRING, and the acquired resistance genes '
                 'have no K-12 orthologs. Zero STRING edges for a determinant reflects '
                 'reference coverage, not an absence of physical interactions.')
    lines.append('')

    lines.append('## J. Interpretation guidelines and limitations')
    lines.append('')
    for line in [
        '- This is exploratory **molecular-context characterization**, not AMR '
        'prediction and not causal inference.',
        '- Graph connectivity does not imply causality. A path between a '
        'determinant and an element means a *graph-derived association*, '
        'nothing more.',
        '- KEGG membership does **not** mean the pathway is active; GO '
        'annotation does **not** mean functional activity was observed; STRING '
        'PPI is **not** strain-specific experimental validation.',
        '- A metabolite sharing a pathway does **not** mean it directly '
        'interacts with the AMR gene; such links are labeled *inferred*.',
        '- A determinant with higher layer coverage is **not** automatically '
        'more important; coverage is a descriptive integration measure.',
        '- Candidate AMR-associated elements are **hypotheses**, not confirmed '
        'resistance determinants.',
        '- Discordant RNA/protein signs are descriptive observations and are '
        '**not** interpreted as proof of post-transcriptional regulation.',
        '- Missing measurements are never converted to zeros; they are reported '
        'as *not observed* (layer exists, no value) or *unavailable* (layer '
        'data do not exist).',
    ]:
        lines.append(line)
    lines.append('')

    (out_dir / 'amr_path_report.md').write_text('\n'.join(lines), encoding='utf-8')


def _layer_key(layer: str) -> str:
    return {'genome': 'L1_genome_status', 'transcriptome': 'L2_transcriptome_status',
            'proteome': 'L3_proteome_status', 'metabolome': 'L4_metabolome_status',
            'kegg': 'L5_kegg_status', 'cog': 'L6_cog_status',
            'go': 'L7_go_status', 'ppi_amr': 'L8_ppi_amr_status'}[layer]


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description='AMR multi-omics molecular-context analysis')
    ap.add_argument('--strains', nargs='*', default=STRAIN_ORDER)
    ap.add_argument('--output', default='reports')
    args = ap.parse_args()

    out_dir = REPO / args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    traces = {}
    issues_by_strain = {}
    evidence_by_strain = {}
    strain_datas = {}
    for name in args.strains:
        data, issues = _load_strain(name)
        if not data:
            print(f"  [{name}] no outputs dir, skipping")
            continue
        markers = MANIFEST.markers(name)
        if not markers:
            print(f"  [{name}] no AMR markers, skipping")
            continue
        issues_by_strain[name] = issues

        data['in_pathway_map'] = _edge_map(data['in_pathway'], 'gene_id', 'annotation_id')
        data['cog_map'] = _edge_map(data['cog'], 'gene_id', 'annotation_id')
        data['go_map'] = _wp_edge_map(data['go'], 'protein_id', 'annotation_id')
        data['clusters_map'] = _wp_edge_map(data['clusters'], 'protein_id', 'annotation_id')
        data['amr_edges_map'] = _edge_map(data['amr_edges'], 'gene_id', 'annotation_id')
        data['genomic_map'] = _edge_map(data['genomic'], 'gene_id', 'gene_target')
        data['_entity_graph'] = _build_entity_graph(data)
        data['_amr_loci'] = set(MANIFEST.amr_loci(name))
        data['_amr_proteins'] = {p for rec in markers for p in rec['protein_ids']}
        strain_datas[name] = data

        strain_traces = []
        strain_evidence = []
        for m in markers:
            rec = trace_marker(name, m, data)
            rec['hyperedges'] = _hyperedge_memberships(data, rec['locus_tag'],
                                                      rec['protein_id'] or None)
            strain_traces.append(rec)
            strain_evidence.extend(build_evidence_rows(rec, data))

        write_trace_csv(strain_traces, out_dir, name)
        write_evidence_csv(strain_evidence, out_dir, name)
        write_hyperedge_detail(name, strain_traces, out_dir)
        traces[name] = strain_traces
        evidence_by_strain[name] = strain_evidence
        mean_cov = (sum(r['layer_coverage'] for r in strain_traces)
                    / len(strain_traces))
        print(f"  [{name}] {len(strain_traces)} markers traced | "
              f"mean layer coverage {mean_cov:.1f}/8 | "
              f"{len(strain_evidence)} evidence rows")

    cross = build_cross_strain(traces)
    if not cross.empty:
        cross.to_csv(out_dir / 'amr_path_cross_strain.csv', index=False)

    candidates = build_candidates(strain_datas, evidence_by_strain)
    if candidates is not None and not candidates.empty:
        candidates.to_csv(out_dir / 'amr_path_candidates.csv', index=False)
        print(f"\nCandidates: {len(candidates)} non-AMR elements ranked")

    write_report(traces, cross, candidates, issues_by_strain,
                 strain_datas, out_dir)
    print("Reports written to:", out_dir.resolve())


if __name__ == '__main__':
    main()
