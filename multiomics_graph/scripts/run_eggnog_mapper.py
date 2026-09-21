"""
run_eggnog_mapper.py
Reproducibly generate `eggnog_annotations.tsv` for every strain using
eggNOG-mapper (external dependency, installed separately).

eggNOG-mapper must be available as `emapper.py` (v2, e.g. via conda):

    conda create -n eggnog python=3.8
    conda activate eggnog
    pip install eggnog-mapper==2.1.15
    conda install -c bioconda diamond
    download_eggnog_data.py --data_dir <EGGNOG_DB> --file eggnog_proteins.dmnd
    download_eggnog_data.py --data_dir <EGGNOG_DB> -y

This script mirrors the format currently hand-produced for B36/MS_14384:

    GeneID<TAB>COG_category<TAB>Description<TAB>Preferred_name

where GeneID is the WP protein accession present in each strain's
`protein.faa`. Only rows with a mapped COG category are kept; unmapped
proteins get `COG_category=S` (Function unknown) so every protein has an
entry and the pipeline can build COG hyperedges.

Run (from repo root):
    python scripts/run_eggnog_mapper.py [--data-dir <EGGNOG_DB>] [--strains MS_14385 MS_14386]
"""

import argparse
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
STRAINS_DIR = REPO_ROOT / 'strains'


def strip_wp_version(gene_id: str) -> str:
    """WP_000002283.1 -> WP_000002283 (strip numeric version suffix)."""
    parts = gene_id.split('.')
    if len(parts) > 1 and parts[-1].isdigit():
        return '.'.join(parts[:-1])
    return gene_id


def protein_faa(strain: str) -> Path:
    """Locate protein.faa for a strain (handles the capital `Genome/` dir)."""
    cfg = STRAINS_DIR / strain / 'config.json'
    if not cfg.exists():
        return Path('')
    data = json.loads(cfg.read_text(encoding='utf-8'))
    paths = data.get('paths', {})
    if paths.get('proteins_faa'):
        p = STRAINS_DIR / strain / paths['proteins_faa']
        return p if p.exists() else Path('')
    # Fall back to any NCBI protein.faa under the strain dir
    for d in sorted((STRAINS_DIR / strain).rglob('protein.faa')):
        return d
    return Path('')


def convert_annotations(emapper_out: Path, out_path: Path):
    """Convert raw eggNOG-mapper `.annotations` into the pipeline TSV.

    The emapper annotations file has leading `# ...` provenance comment
    lines and a single-hash data header line (`#query ...`). Column order
    for eggNOG-mapper v2:
        query seed_ortholog evalue score eggNOG_OGs max_annot_lvl
        COG_category Description Preferred_name GOs EC KEGG_ko ...
    This reader detects the data header by the presence of `seed_ortholog`,
    then extracts the query, COG category, description and preferred name
    regardless of exact names. WP accessions are stripped to their
    unversioned form (WP_000002283.1 -> WP_000002283) to match the reference
    files; '-' COG categories become 'S' (Function unknown).
    """
    rows = []
    cols = None
    with emapper_out.open(encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            # emapper v2.1.13 data header row is a single `#query` line
            # containing `seed_ortholog`; other `#` lines are provenance.
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
            cog = 'S' if (cog == '-' or not cog) else cog
            desc = desc if desc else 'Unknown'
            rows.append((strip_wp_version(query), cog, desc, pref))

    # Deduplicate by query accession (last wins), sort for stable output.
    seen = {}
    for q, cog, desc, pref in rows:
        seen[q] = (cog, desc, pref)
    with out_path.open('w', encoding='utf-8', newline='') as fh:
        fh.write('GeneID\tCOG_category\tDescription\tPreferred_name\n')
        for q in sorted(seen):
            cog, desc, pref = seen[q]
            fh.write(f"{q}\t{cog}\t{desc}\t{pref}\n")


def main():
    ap = argparse.ArgumentParser(description='Generate eggNOG annotations per strain')
    ap.add_argument('--data-dir', default=os.environ.get('eggNOG_DB', ''),
                    help='eggNOG data dir (where eggnog_proteins.dmnd lives). '
                         'If empty, uses the conda env default.')
    ap.add_argument('--strains', nargs='*', default=None,
                    help='Strains to process (default: all in strains/)')
    ap.add_argument('--emapper', default='emapper.py',
                    help='eggNOG-mapper executable name/path')
    ap.add_argument('--keep-tmp', action='store_true',
                    help='Do not delete the temporary output dir')
    args = ap.parse_args()

    strains = args.strains or [d.name for d in STRAINS_DIR.iterdir()
                               if (d / 'config.json').exists()]
    for strain in strains:
        fasta = protein_faa(strain)
        out = STRAINS_DIR / strain / 'eggnog_annotations.tsv'
        if not fasta.exists():
            print(f"  [{strain}] protein.faa not found, skipping")
            continue
        print(f"  [{strain}] annotating {fasta.name}")

        tmp = Path(tempfile.mkdtemp(prefix=f"eggnog_{strain}_"))
        try:
            # NOTE: emapper's `-o/--output` is the *base name* for output
            # files, not a format list. Files are written into `--output_dir`
            # as `{prefix}.emapper.annotations`. Pass both explicitly so the
            # annotations land in our temp dir regardless of the CWD.
            cmd = [
                args.emapper,
                '-i', str(fasta),
                '-o', 'out',
                '--output_dir', str(tmp),
                '--dmnd_iterate', 'no',
                '--cpu', '8',
                '--override',
            ]
            if args.data_dir:
                cmd += ['--data_dir', str(args.data_dir)]
            print('   running:', ' '.join(cmd))
            subprocess.run(cmd, check=True)
            ann = tmp / 'out.emapper.annotations'
            if not ann.exists():
                raise FileNotFoundError(f"expected {ann}")
            convert_annotations(ann, out)
            n = sum(1 for _ in out.open(encoding='utf-8')) - 1
            print(f"  [{strain}] wrote {out.name} with {n} proteins")
        finally:
            if not args.keep_tmp:
                shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    main()