"""Convert eggNOG-mapper .emapper.annotations output into the pipeline TSV format.

Produces the `eggnog_annotations.tsv` expected by
`annotation/eggnog.py::load_precomputed`:

    GeneID<TAB>COG_category<TAB>Description<TAB>Preferred_name
           <TAB>KEGG_ko<TAB>KEGG_Pathway

Matches the format of the reference files committed for B36 / MS_14384:
versioned WP accessions (WP_000002283.1) are stripped to unversioned
(WP_000002283), '-' COG categories become 'S' (Function unknown), and only
rows present in the emapper output are kept (the pipeline fills placeholders
for any unannotated gene via `annotate_by_locus_tag`).

KEGG_ko (e.g. ko:K17836) and KEGG_Pathway (e.g. ko00311,map01130) are carried
through from the emapper output. These enrich the graph's KEGG pathway layer
for acquired resistance genes that the `eco`-bridged UniProt cross-reference
route cannot cover (plasmid-borne determinants).

Usage:
    python scripts/convert_emapper_to_tsv.py <strain> [--in FILE] [--out FILE]
"""

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
STRAINS_DIR = REPO_ROOT / 'strains'


def strip_wp_version(gene_id: str) -> str:
    """WP_000002283.1 -> WP_000002283 (strip numeric version suffix)."""
    parts = gene_id.split('.')
    if len(parts) > 1 and parts[-1].isdigit():
        return '.'.join(parts[:-1])
    return gene_id


def convert(emapper_in: Path, out_path: Path) -> int:
    rows = []
    cols = None
    with emapper_in.open(encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#') and 'seed_ortholog' in stripped:
                cols = stripped.lstrip('#').split('\t')
                continue
            if stripped.startswith('#'):
                continue
            parts = stripped.split('\t')
            if cols is None:
                cols = parts
                continue
            rec = dict(zip(cols, parts))
            query = rec.get('query_name', rec.get('query', ''))
            cog = rec.get('COG_category', '-')
            desc = rec.get('Description', '')
            pref = rec.get('Preferred_name', '')
            kegg_ko = rec.get('KEGG_ko', rec.get('KEGG_KO', ''))
            kegg_path = rec.get('KEGG_Pathway', '')
            cog = 'S' if (cog == '-' or not cog) else cog
            desc = desc if desc else 'Unknown'
            rows.append((strip_wp_version(query), cog, desc, pref, kegg_ko, kegg_path))

    # Deduplicate by gene accession (last wins), sort for stable output.
    seen = {}
    for q, cog, desc, pref, kegg_ko, kegg_path in rows:
        seen[q] = (cog, desc, pref, kegg_ko, kegg_path)
    with out_path.open('w', encoding='utf-8', newline='') as fh:
        fh.write('GeneID\tCOG_category\tDescription\tPreferred_name\tKEGG_ko\tKEGG_Pathway\n')
        for q in sorted(seen):
            cog, desc, pref, kegg_ko, kegg_path = seen[q]
            fh.write(f"{q}\t{cog}\t{desc}\t{pref}\t{kegg_ko}\t{kegg_path}\n")
    return len(seen)


def main():
    ap = argparse.ArgumentParser(description='Convert emapper annotations to pipeline TSV')
    ap.add_argument('strain', help='Strain name (uses strains/<strain>/ by default)')
    ap.add_argument('--in', dest='in_file', default=None,
                    help='Path to .emapper.annotations (default: strains/<strain>/eggnog_annotations.emapper.annotations)')
    ap.add_argument('--out', dest='out_file', default=None,
                    help='Output TSV (default: strains/<strain>/eggnog_annotations.tsv)')
    args = ap.parse_args()

    strain_dir = STRAINS_DIR / args.strain
    in_file = Path(args.in_file) if args.in_file else \
        strain_dir / 'eggnog_annotations.emapper.annotations'
    out_file = Path(args.out_file) if args.out_file else \
        strain_dir / 'eggnog_annotations.tsv'

    if not in_file.exists():
        raise SystemExit(f"input not found: {in_file}")

    n = convert(in_file, out_file)
    print(f"wrote {out_file} with {n} genes")


if __name__ == '__main__':
    main()
