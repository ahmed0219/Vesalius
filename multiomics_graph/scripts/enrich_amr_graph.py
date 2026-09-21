"""
scripts/enrich_amr_graph.py
Surgically enrich saved per-strain graph outputs for AMR determinants.

Why
---
A full pipeline re-run (network-bound annotation) is not needed to close the
KEGG/COG/GO/encodes/proximity gaps for AMR determinants. The saved outputs in
`outputs/<strain>/graph/` are CSV edge files that the read-only analysis
script (`scripts/amr_path_analysis.py`) consumes directly. This script adds
the missing rows for the determinants so the trace shows the recovered layers:

  * encodes          -- gene -> protein edges for genome-only determinants
                        (tet(A), aadA5, blaCTX-M-15 in B36; aadA1, floR in
                        MS_14386) that have no aligned multi-omics row but DO
                        carry a RefSeq protein in `genome_genes.csv`.
  * COG              -- gene -> COG class edges from the eggNOG annotation TSV
                        (which now carries KEGG columns; see
                        scripts/convert_emapper_to_tsv.py).
  * KEGG in_pathway  -- gene -> KEGG pathway (ko/map) edges from the eggNOG
                        KEGG_Pathway column. The UniProt `eco` bridge used by
                        the pipeline drops acquired/plasmid-borne genes.
  * GO               -- protein -> GO edges from
                        `reports/amr_uniprot_recovery.csv`.
  * genomic proximity-- gene -> gene edges from the *unfiltered*
                        `genome_genomic_edges.csv` for determinants the graph
                        did not previously include.

It also updates `hyperedges.csv` / `hypergraph_nodes.csv` so the new
annotation memberships become hyperedge connector hubs (so BFS path tracing
can reach co-annotated genes via COG/GO/KEGG). Existing hyperedges are
bumped with the determinant's locus node (the `nodes_map` uses integer
`hyperedge_id` keys so already-present hyperedges are updated, not just
newly created ones).

Usage (from multiomics_graph/):
    python scripts/enrich_amr_graph.py
    python scripts/enrich_amr_graph.py --strains B36 MS_14386
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from amr.amr import AMRManifest

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / 'outputs'
STRAINS = REPO.parent / 'strains'
RECOVERY_PATH = REPO / 'reports' / 'amr_uniprot_recovery.csv'
STRAIN_ORDER = ['B36', 'MS_14384', 'MS_14386', 'MS_14387']
COG_CATEGORIES = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')


def _read(path):
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, dtype=str)
    except Exception:
        return None


def _norm_wp(pid: str) -> str:
    if not pid:
        return ''
    return pid.split('.')[0]


def _unique(df, cols):
    if df is None or df.empty:
        return pd.DataFrame(columns=cols)
    return df.drop_duplicates(subset=cols)


def _append_unique(base, rows, cols):
    new = pd.DataFrame(rows, columns=cols)
    merged = pd.concat([base, new], ignore_index=True) if base is not None \
        else new
    return merged.drop_duplicates(subset=cols, keep='first')


def load_eggnog_lookup(strain: str) -> dict:
    """Return {unversioned_wp: {'cog': [letters], 'kegg': [pathways]}}."""
    tsv = STRAINS / strain / 'eggnog_annotations.tsv'
    if not tsv.exists():
        return {}
    df = pd.read_csv(tsv, sep='\t', dtype=str)
    if 'GeneID' not in df.columns:
        return {}
    out = {}
    for _, r in df.iterrows():
        gid = str(r.get('GeneID', '')).strip()
        if not gid:
            continue
        cog = str(r.get('COG_category', ''))
        classes = [c for c in cog if c in COG_CATEGORIES]
        kegg_raw = str(r.get('KEGG_Pathway', ''))
        kegg = [p for p in (x.strip() for x in kegg_raw.split(','))
                if p and p.startswith(('ko', 'map'))]
        if classes or kegg:
            out[gid] = {'cog': classes, 'kegg': kegg}
    return out


def load_recovery() -> dict:
    """Return {protein_id: [GO ids]} from the recovery CSV (versioned keys)."""
    rec = _read(RECOVERY_PATH)
    if rec is None:
        return {}
    out = {}
    for _, r in rec.iterrows():
        pid = str(r.get('protein_id', '') or '').strip()
        go_raw = str(r.get('go_ids', '') or '')
        if not pid or not go_raw or go_raw.lower() == 'nan':
            continue
        gos = [g.strip() for g in go_raw.split(';') if g.strip()]
        if gos:
            out.setdefault(pid, []).extend(gos)
    return out


def enrich_strain(strain: str) -> dict:
    gd = OUT_DIR / strain / 'graph'
    sd = OUT_DIR / strain
    report = {'strain': strain, 'encodes': 0, 'cog': 0, 'kegg': 0,
              'go': 0, 'proximity': 0, 'hyperedges_added': 0,
              'hyperedges_updated': 0}
    if not gd.exists():
        return report

    manifest = AMRManifest()
    markers = manifest.markers(strain)
    loci = [mk['locus_tags'][0] for mk in markers if mk['locus_tags']]

    genome = _read(sd / 'genome_genes.csv')
    gene_to_protein = {}
    if genome is not None and 'GeneID' in genome.columns:
        pidcol = 'ProteinID' if 'ProteinID' in genome.columns else None
        for _, r in genome.iterrows():
            g = str(r.get('GeneID', '')).strip()
            if not g:
                continue
            p = str(r.get(pidcol, '') or '').strip() if pidcol else ''
            if p and p.lower() != 'nan':
                gene_to_protein[g] = p

    eggnog = load_eggnog_lookup(strain)
    recovery = load_recovery()

    # --- edge files ------------------------------------------------------
    enc_path = gd / 'edges_gene_encodes_protein.csv'
    cog_path = gd / 'edges_gene_cog_category_annotation.csv'
    kegg_path = gd / 'edges_gene_in_pathway_annotation.csv'
    go_path = gd / 'edges_protein_annotated_by_annotation.csv'
    prox_path = gd / 'edges_gene_genomic_proximity_gene.csv'

    enc = _read(enc_path)
    cog = _read(cog_path)
    kegg = _read(kegg_path)
    go = _read(go_path)
    prox = _read(prox_path)

    enc_rows, cog_rows, kegg_rows, go_rows, prox_rows = [], [], [], [], []

    for locus in loci:
        protein = gene_to_protein.get(locus, '')
        wp = _norm_wp(protein)

        # encodes edge
        if protein and enc is not None and not \
                ((enc['gene_id'].astype(str) == locus) &
                 (enc['protein_id'].astype(str) == protein)).any():
            enc_rows.append({'gene_id': locus, 'protein_id': protein,
                             'relation': 'encodes'})
            report['encodes'] += 1

        # COG edges (only if the graph does not already carry them)
        if cog is not None:
            existing = set(cog.loc[cog['gene_id'].astype(str) == locus,
                                   'annotation_id'].astype(str))
        else:
            existing = set()
        for cls in eggnog.get(wp, {}).get('cog', []):
            if cls not in existing:
                cog_rows.append({'gene_id': locus,
                                 'annotation_id': cls,
                                 'relation': 'cog_category'})
                report['cog'] += 1

        # KEGG in_pathway edges
        if kegg is not None:
            existing = set(kegg.loc[kegg['gene_id'].astype(str) == locus,
                                    'annotation_id'].astype(str))
        else:
            existing = set()
        for pw in eggnog.get(wp, {}).get('kegg', []):
            if pw not in existing:
                kegg_rows.append({'gene_id': locus,
                                  'annotation_id': pw,
                                  'relation': 'in_pathway'})
                report['kegg'] += 1

        # GO edges (protein -> GO from recovery)
        if protein and go is not None:
            existing = set(go.loc[go['protein_id'].astype(str) == protein,
                                  'annotation_id'].astype(str))
            for goid in recovery.get(protein, []) + recovery.get(wp, []):
                if goid not in existing:
                    go_rows.append({'protein_id': protein,
                                    'annotation_id': goid,
                                    'relation': 'annotated_by'})
                    report['go'] += 1

        # genomic proximity (unfiltered genome edges)
        if prox is not None:
            existing = set(prox.loc[prox['gene_id'].astype(str) == locus,
                                    'gene_target'].astype(str))
        else:
            existing = set()
        gge = _read(sd / 'genome_genomic_edges.csv')
        if gge is not None and 'GeneID_A' in gge.columns:
            hits = gge[gge['GeneID_A'].astype(str) == locus]
            hits = pd.concat([hits, gge[gge['GeneID_B'].astype(str) == locus]],
                             ignore_index=True)
            for _, r in hits.iterrows():
                a = str(r.get('GeneID_A', '')).strip()
                b = str(r.get('GeneID_B', '')).strip()
                tgt = b if a == locus else a
                if not tgt or tgt == locus:
                    continue
                if tgt in existing:
                    continue
                # only connect to genes that exist in the graph edge files
                if _in_graph(tgt, enc, cog, kegg, prox):
                    prox_rows.append({'gene_id': locus,
                                      'gene_target': tgt,
                                      'relation': 'genomic_proximity'})
                    report['proximity'] += 1

    # --- write back ------------------------------------------------------
    enc = _append_unique(enc, enc_rows, ['gene_id', 'protein_id', 'relation'])
    cog = _append_unique(cog, cog_rows, ['gene_id', 'annotation_id', 'relation'])
    kegg = _append_unique(kegg, kegg_rows, ['gene_id', 'annotation_id', 'relation'])
    go = _append_unique(go, go_rows, ['protein_id', 'annotation_id', 'relation'])
    prox = _append_unique(prox, prox_rows, ['gene_id', 'gene_target', 'relation'])
    enc.to_csv(enc_path, index=False)
    cog.to_csv(cog_path, index=False)
    kegg.to_csv(kegg_path, index=False)
    go.to_csv(go_path, index=False)
    prox.to_csv(prox_path, index=False)

    # --- hyperedges ------------------------------------------------------
    he_path = gd / 'hyperedges.csv'
    he = _read(he_path)
    if he is not None:
        added, updated = _update_hyperedges(
            he, gd, loci, gene_to_protein, eggnog, recovery)
        report['hyperedges_added'] = added
        report['hyperedges_updated'] = updated
    return report


def _in_graph(node: str, *dfs) -> bool:
    for df in dfs:
        if df is None:
            continue
        for col in df.columns:
            if node in df[col].astype(str).values:
                return True
    return False


def _update_hyperedges(he: pd.DataFrame, gd: Path, loci, gene_to_protein,
                       eggnog, recovery) -> tuple:
    """Add determinant memberships to kegg_pathway / cog_category / go_term
    hyperedges, creating new hyperedges when a term is absent."""
    he_path = gd / 'hyperedges.csv'
    he = he.copy()
    nodes_map = {}   # hyperedge_id -> set(nodes)
    for _, r in he.iterrows():
        nodes_map[int(r['hyperedge_id'])] = set(
            str(n).strip() for n in str(r.get('nodes', '')).split(';') if n)
    next_id = int(he['hyperedge_id'].astype(int).max()) + 1
    changed = {}     # hyperedge_id -> added nodes
    new_rows = []

    def bump(he_id, nid):
        if he_id in nodes_map:
            if nid not in nodes_map[he_id]:
                nodes_map[he_id].add(nid)
                changed.setdefault(he_id, []).append(nid)
        else:
            pass

    for locus in loci:
        protein = gene_to_protein.get(locus, '')
        wp = _norm_wp(protein)

        for pw in eggnog.get(wp, {}).get('kegg', []):
            hit = he[(he['type'] == 'kegg_pathway') & (he['name'] == pw)]
            if hit.empty:
                new_rows.append({'hyperedge_id': next_id, 'type': 'kegg_pathway',
                                 'name': pw, 'nodes': locus, 'n_nodes': 1})
                nodes_map[next_id] = {locus}
                next_id += 1
            else:
                bump(int(hit.iloc[0]['hyperedge_id']), locus)

        for cls in eggnog.get(wp, {}).get('cog', []):
            hit = he[(he['type'] == 'cog_category')
                     & (he['name'].str.startswith(f'COG_{cls}:'))]
            if hit.empty:
                name = f'COG_{cls}: {cls}'
                new_rows.append({'hyperedge_id': next_id, 'type': 'cog_category',
                                 'name': name, 'nodes': locus, 'n_nodes': 1})
                nodes_map[next_id] = {locus}
                next_id += 1
            else:
                bump(int(hit.iloc[0]['hyperedge_id']), locus)

        for goid in recovery.get(protein, []) + recovery.get(wp, []):
            if not protein:
                continue
            hit = he[(he['type'] == 'go_term')
                     & (he['name'].str.startswith(goid + ':'))]
            if hit.empty:
                name = f'{goid}: {goid}'
                node = protein or locus
                new_rows.append({'hyperedge_id': next_id, 'type': 'go_term',
                                 'name': name, 'nodes': node, 'n_nodes': 1})
                nodes_map[next_id] = {node}
                next_id += 1
            else:
                node = protein or locus
                if node:
                    bump(int(hit.iloc[0]['hyperedge_id']), node)

    rows_out = []
    for _, r in he.iterrows():
        he_id = int(r['hyperedge_id'])
        nodes = nodes_map.get(he_id, set())
        rows_out.append({'hyperedge_id': he_id, 'type': r['type'],
                         'name': r['name'],
                         'nodes': ';'.join(sorted(nodes)),
                         'n_nodes': len(nodes)})
    for row in new_rows:
        rows_out.append(row)

    pd.DataFrame(rows_out).to_csv(he_path, index=False)

    # hypergraph_nodes.csv sidecar
    hn_path = gd / 'hypergraph_nodes.csv'
    hn = _read(hn_path)
    all_nodes = set()
    if hn is not None:
        all_nodes = set(hn['node_id'].astype(str))
    for nid in nodes_map.values():
        all_nodes |= nid
    if hn is None:
        hn = pd.DataFrame({'node_id': sorted(all_nodes),
                           'node_type': 'genome'})
    else:
        existing = set(hn['node_id'].astype(str))
        add = sorted(all_nodes - existing)
        if add:
            hn = pd.concat([
                hn,
                pd.DataFrame({'node_id': add,
                              'node_type': 'genome'}),
            ], ignore_index=True)
    hn.to_csv(hn_path, index=False)

    return len(new_rows), len(changed)


def main():
    ap = argparse.ArgumentParser(description='Enrich AMR graph outputs')
    ap.add_argument('--strains', nargs='*', default=STRAIN_ORDER)
    args = ap.parse_args()

    print("Enriching saved graph outputs for AMR determinants\n")
    for name in args.strains:
        r = enrich_strain(name)
        print(f"  [{name}] encodes+{r['encodes']} cog+{r['cog']} "
              f"kegg+{r['kegg']} go+{r['go']} proximity+{r['proximity']} | "
              f"hyperedges +{r['hyperedges_added']} updated {r['hyperedges_updated']}")


if __name__ == '__main__':
    main()
