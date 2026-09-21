"""
scripts/recover_amr_uniprot.py
Recover UniProt annotations for curated AMR determinants.

Why
---
The upstream pipeline maps RefSeq -> UniProt in batches of ~1000 IDs with a
search result cap of 500 entries per request (`_search_refseq_batch` in
annotation/uniprot.py). When a batch's matching UniProt entries exceed 500,
later IDs in the batch are silently dropped from the mapping. As a result,
AMR determinants that ARE represented in the multi-omics graph (e.g. dfrA17,
sul1, blaOXA-1, mph(A), aac(6')-Ib-cr5 in B36) received no UniProt/GO
annotation at all.

This script re-runs the SAME UniProtKB search -- one WP accession per request,
so there is no truncation -- reusing the pipeline's own candidate-resolution
and entry-parsing code. It writes the recovered annotations to
`reports/amr_uniprot_recovery.csv`. The read-only analysis script
(`scripts/amr_path_analysis.py`) consumes this file to enrich the GO layer.

It does NOT modify or re-run the upstream pipeline.

Usage (from multiomics_graph/):
    python scripts/recover_amr_uniprot.py
    python scripts/recover_amr_uniprot.py --strains B36 MS_14386
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from annotation.uniprot import UniProtAnnotator
from amr.amr import AMRManifest

REPO = Path(__file__).resolve().parent.parent
MANIFEST = AMRManifest()
STRAIN_ORDER = ['B36', 'MS_14384', 'MS_14386', 'MS_14387']


def recover_for_wp(ann: UniProtAnnotator, wp: str):
    """Recover the best UniProt entry for a single WP accession (no truncation).

    Uses the pipeline's own search / candidate-resolution / parsing code.
    """
    norm = ann._normalize_wp(wp)
    if not norm:
        return None
    entries = ann._search_refseq_batch([norm])
    if not entries:
        return None
    index = ann._build_refseq_index(entries)
    best = ann._resolve_candidates(index.get(norm, []))
    if not best:
        return None
    entry = next((e for e in entries if e.get('primaryAccession') == best), None)
    if entry is None:
        return None
    parsed = ann._parse_entry(entry)
    return parsed


def go_ids_of(parsed: dict) -> list:
    out = []
    for col in ('GO_BP', 'GO_MF', 'GO_CC'):
        for g in str(parsed.get(col, '')).split(';'):
            g = g.strip()
            if g:
                out.append(g)
    return sorted(set(out))


def main():
    ap = argparse.ArgumentParser(description='Recover UniProt annotations for AMR determinants')
    ap.add_argument('--output', default='reports')
    ap.add_argument('--strains', nargs='*', default=STRAIN_ORDER)
    args = ap.parse_args()

    out_dir = REPO / args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    ann = UniProtAnnotator(organism_id='562')

    # Collect unique WP accessions across requested strains
    wanted = []
    seen = set()
    for name in args.strains:
        for mk in MANIFEST.markers(name):
            for p in mk.get('protein_ids', []):
                if not p or '(pseudogene)' in p:
                    continue
                norm = ann._normalize_wp(p)
                if norm not in seen:
                    seen.add(norm)
                    wanted.append((name, mk, p, norm))

    # Recover per unique WP
    recovered = {}
    for norm in sorted({w[3] for w in wanted}):
        parsed = recover_for_wp(ann, norm)
        if parsed:
            recovered[norm] = parsed
            n_go = len(go_ids_of(parsed))
            print(f"  {norm} -> {parsed.get('UniProtID')} | "
                  f"{parsed.get('ProteinName', '')[:40]} | GO={n_go}")
        else:
            print(f"  {norm} -> NO UniProt entry found")

    rows = []
    for name, mk, p, norm in wanted:
        parsed = recovered.get(norm)
        if parsed is None:
            rows.append({
                'strain': name, 'marker': mk['name'], 'amr_class': mk['amr_class'],
                'locus_tag': mk['locus_tags'][0] if mk['locus_tags'] else '',
                'protein_id': p, 'uniprot_id': '',
                'protein_name': '', 'function': '', 'ec_number': '',
                'go_ids': '', 'go_bp': '', 'go_mf': '', 'go_cc': '',
                'n_go': 0, 'recovered': False,
            })
            continue
        rows.append({
            'strain': name, 'marker': mk['name'], 'amr_class': mk['amr_class'],
            'locus_tag': mk['locus_tags'][0] if mk['locus_tags'] else '',
            'protein_id': p, 'uniprot_id': parsed.get('UniProtID', ''),
            'protein_name': parsed.get('ProteinName', ''),
            'function': parsed.get('Function', ''),
            'ec_number': parsed.get('EC_number', ''),
            'go_ids': '; '.join(go_ids_of(parsed)),
            'go_bp': parsed.get('GO_BP', ''), 'go_mf': parsed.get('GO_MF', ''),
            'go_cc': parsed.get('GO_CC', ''),
            'n_go': len(go_ids_of(parsed)),
            'recovered': True,
        })

    df = pd.DataFrame(rows)
    path = out_dir / 'amr_uniprot_recovery.csv'
    df.to_csv(path, index=False)
    n_ok = int(df['recovered'].sum())
    n_go = int(df['n_go'].sum())
    print(f"\nSaved {len(df)} determinant rows to {path} "
          f"({n_ok} recovered, {n_go} GO terms total)")


if __name__ == '__main__':
    main()
