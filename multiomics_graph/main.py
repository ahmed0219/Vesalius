#!/usr/bin/env python3
"""
main.py
Multi-strain heterogeneous multi-omics graph construction pipeline.

Orchestrates:
  1. Genomics   — GFF parsing → GeneID/ProteinID mapping
  2. Transcriptomics — RNA-seq log2 CPM → condition means + logFC
  3. Proteomics — SWATH-MS protein areas → condition means + logFC
  4. Integration — cross-omics alignment via central dogma
  5. Annotation — UniProt/STRING/KEGG/eggNOG retrieval
  6. Graph construction — PyG HeteroData with multiple edge types
  7. Hypergraph construction — biological module hyperedges
  8. Visualisation — summary plots and network figures
  9. Export — CSV, NPY, PT files for downstream GNN training

Usage:
    python main.py --strain B36 --genome path/to/gff --rna path/to/expression --prot path/to/protein_groups.xlsx
    python main.py --all-strains   # process all discovered strains
"""

import sys
import json
import time
import argparse
import warnings
import hashlib
import platform
import importlib.metadata
from pathlib import Path
from typing import Dict, List, Optional

import torch
import pandas as pd
import numpy as np

warnings.filterwarnings('ignore')

# Ensure project root is in path
PROJECT_ROOT = Path(__file__).parent.absolute()
sys.path.insert(0, str(PROJECT_ROOT))

# Ensure the RNG families used anywhere in the pipeline (Louvain clustering,
# correlation ties, etc.) share a single reproducible seed so full-pipeline
# reruns are deterministic. Recorded in every per-strain manifest.json.
PIPELINE_SEED = 42

from preprocessing import GenomicsProcessor, TranscriptomicsProcessor, ProteomicsProcessor
from annotation.uniprot import UniProtAnnotator
from annotation.string_api import STRINGClient
from annotation.kegg import KEGGAnnotator
from annotation.eggnog import EggNOGAnnotator, COG_CATEGORIES
from amr.amr import AMRManifest, AMR_CLASS_NAMES
from graph.heterogeneous_graph import HeterogeneousGraphBuilder
from graph.hypergraph import BiologicalHypergraphBuilder
from visualization.network_plots import NetworkVisualizer

BASE_DIR = Path(__file__).parent.parent.absolute()  # agent/


# ──────────────────────────────────────────────
# Strain configuration — per-strain folder or JSON file
# ──────────────────────────────────────────────

STRAINS_DIR = BASE_DIR / 'strains'
STRAINS_CONFIG_PATH = Path(__file__).parent / 'strains_config.json'


def list_available_strains() -> list:
    """Return list of available strain names from both layouts."""
    strains = set()

    # Layout 1: per-strain folders (strains/<name>/config.json)
    if STRAINS_DIR.exists():
        for entry in STRAINS_DIR.iterdir():
            if entry.is_dir() and (entry / 'config.json').exists() and entry.name != 'TEMPLATE':
                strains.add(entry.name)

    # Layout 2: centralized JSON (strains_config.json)
    if STRAINS_CONFIG_PATH.exists():
        with open(STRAINS_CONFIG_PATH, 'r') as f:
            strains.update(json.load(f).keys())

    return sorted(strains)


def _find_strain_dir(strain: str) -> Path | None:
    """Check if strain has a per-strain folder with config.json."""
    candidate = STRAINS_DIR / strain / 'config.json'
    return candidate if candidate.exists() else None


def _load_strain_from_folder(strain: str) -> dict:
    """Load config from strains/<name>/config.json, resolve paths relative to that directory."""
    config_path = STRAINS_DIR / strain / 'config.json'
    with open(config_path, 'r') as f:
        cfg = json.load(f)

    strain_root = config_path.parent

    # Flatten into pipeline config dict
    flat = {
        'name': cfg.get('name', strain),
        'ncbi_tax_id': cfg.get('ncbi_tax_id', 562),
        'kegg_code': cfg.get('kegg_code', 'eco'),
        'regulation_std_multiplier': cfg.get('regulation_std_multiplier', 0.0),
        'string_species_primary': cfg.get('string_species_primary', 511145),
        'string_species_fallback': cfg.get('string_species_fallback', 511145),
        'string_species': cfg.get('string_species_primary', 511145),
    }

    # ── Resolve file paths ─────
    # Helper: resolve explicit path (relative to strain dir) or discover by convention
    def _resolve(key: str, subdir: str, pattern: str) -> Path:
        """Try explicit path from config, then scan strain_dir/subdir for pattern."""
        # Check explicit paths section
        paths = cfg.get('paths', {})
        if key in paths:
            p = Path(paths[key])
            return (strain_root / p).resolve() if not p.is_absolute() else p
        # Scan subdirectory by convention
        scan_dir = strain_root / subdir
        if scan_dir.exists():
            matches = list(scan_dir.rglob(pattern)) if '*' in pattern else list(scan_dir.rglob(pattern))
            if matches:
                return matches[0]
        return Path('')

    flat['genome_gff'] = _resolve('genome_gff', 'genome', '*.gff')
    flat['genome_fna'] = _resolve('genome_fna', 'genome', '*.fna')
    flat['cds_fna'] = _resolve('cds_fna', 'genome', 'cds_from_genomic.fna')
    flat['protein_faa'] = _resolve('protein_faa', 'genome', '*.faa')
    # RNA counts: prefer files matching *counts* over other files
    rna_explicit = cfg.get('paths', {}).get('rna_counts')
    if rna_explicit:
        flat['rna_counts_file'] = _resolve('rna_counts', 'transcriptomic', '*')
    else:
        rna_dir = strain_root / 'transcriptomic'
        if rna_dir.exists():
            files = sorted(rna_dir.iterdir())
            counts = [f for f in files if 'count' in f.name.lower()]
            flat['rna_counts_file'] = counts[0] if counts else files[0] if files else Path('')
        else:
            flat['rna_counts_file'] = Path('')

    # RNA log2 CPM (preferred for logFC/regulation): discover *cpmlog2* file.
    # Raw counts are not comparable across samples, so logFC/regulation must be
    # computed from the log2 CPM matrix rather than from raw count differences.
    rna_dir = strain_root / 'transcriptomic'
    cpmlog2 = []
    if rna_dir.exists():
        cpmlog2 = sorted(rna_dir.rglob('*cpmlog2*'))
    flat['rna_file'] = cpmlog2[0] if cpmlog2 else Path('')

    # Transcriptomics sample IDs
    tx = cfg.get('transcriptomics', {})
    flat['rna_rpmi_samples'] = tx.get('rpmi_samples', [])
    flat['rna_sera_samples'] = tx.get('sera_samples', [])

    # Proteomics: use subdir or explicit paths
    px = cfg.get('proteomics', {})
    flat['prot_rpmi_ids'] = px.get('rpmi_ids', [])
    flat['prot_sera_ids'] = px.get('sera_ids', [])
    flat['prot_gene_col'] = px.get('gene_col', None)

    # Metabolomics
    meta = cfg.get('metabolomics', {})
    flat['metabolomics_enabled'] = meta.get('enabled', False)
    flat['metabolomics_dir'] = meta.get('dir', '')
    flat['metabolomics_maf_files'] = meta.get('maf_files', [])
    flat['metabolomics_sample_file'] = meta.get('sample_file', '')
    flat['metabolomics_rpmi_samples'] = meta.get('rpmi_samples', [])
    flat['metabolomics_sera_samples'] = meta.get('sera_samples', [])
    # Resolve metabolomics paths relative to the strain dir
    if flat['metabolomics_dir']:
        meta_dir = Path(flat['metabolomics_dir'])
        flat['metabolomics_dir'] = str(meta_dir if meta_dir.is_absolute() else (strain_root / meta_dir).resolve())
    if flat['metabolomics_sample_file']:
        p = Path(flat['metabolomics_sample_file'])
        flat['metabolomics_sample_file'] = str(p if p.is_absolute() else (strain_root / p).resolve())
    resolved_maf = []
    for mf in flat['metabolomics_maf_files']:
        p = Path(mf)
        resolved_maf.append(str(p if p.is_absolute() else (strain_root / p).resolve()))
    flat['metabolomics_maf_files'] = resolved_maf

    # SWATH
    sw = px.get('swath', {})
    flat['prot_swath_enabled'] = sw.get('enabled', True)
    swath_file = sw.get('file', '')
    if swath_file:
        p = Path(swath_file)
        flat['prot_file'] = p if p.is_absolute() else (strain_root / p).resolve()
    else:
        prot_dir = strain_root / 'proteomics'
        xlsx_files = list(prot_dir.glob('*.xlsx')) if prot_dir.exists() else []
        flat['prot_file'] = xlsx_files[0] if xlsx_files else Path('')
    flat['prot_sheet'] = sw.get('sheet', 'Area - proteins')
    flat['prot_protein_col'] = sw.get('protein_col', 'Protein')
    flat['prot_swath_rpmi_template'] = sw.get('rpmi_template', '')
    flat['prot_swath_sera_template'] = sw.get('sera_template', '')

    return flat


def _load_strain_from_central_config(strain: str) -> dict:
    """Load config from the centralized strains_config.json (legacy layout)."""
    with open(STRAINS_CONFIG_PATH, 'r') as f:
        all_configs = json.load(f)

    cfg = all_configs[strain]

    flat = {
        'name': cfg.get('name', strain),
        'ncbi_tax_id': cfg.get('ncbi_tax_id', 562),
        'kegg_code': cfg.get('kegg_code', 'eco'),
        'regulation_std_multiplier': cfg.get('regulation_std_multiplier', 0.0),
        'string_species_primary': cfg.get('string_species_primary', 511145),
        'string_species_fallback': cfg.get('string_species_fallback', 511145),
        'string_species': cfg.get('string_species_primary', 511145),
    }

    # Resolve file paths relative to BASE_DIR
    paths = cfg.get('paths', {})
    def _resolve(p):
        p = Path(p or '')
        return p if p.is_absolute() else BASE_DIR / p

    flat['genome_gff'] = _resolve(paths.get('genome_gff', ''))
    flat['genome_fna'] = _resolve(paths.get('genome_fna', ''))
    flat['cds_fna'] = _resolve(paths.get('cds_fna', ''))
    flat['protein_faa'] = _resolve(paths.get('protein_faa', ''))
    flat['rna_file'] = _resolve(paths.get('rna_cpmlog2', ''))
    flat['rna_counts_file'] = _resolve(paths.get('rna_counts', ''))
    flat['prot_file'] = _resolve(paths.get('prot_file', ''))

    # Transcriptomics
    tx = cfg.get('transcriptomics', {})
    flat['rna_rpmi_samples'] = tx.get('rpmi_samples', [])
    flat['rna_sera_samples'] = tx.get('sera_samples', [])

    # Proteomics
    px = cfg.get('proteomics', {})
    flat['prot_sheet'] = px.get('sheet', 'Area - proteins')
    flat['prot_protein_col'] = px.get('protein_col', 'Protein')
    flat['prot_gene_col'] = px.get('gene_col', None)
    flat['prot_rpmi_ids'] = px.get('rpmi_ids', [])
    flat['prot_sera_ids'] = px.get('sera_ids', [])

    # Metabolomics (legacy central config)
    meta = cfg.get('metabolomics', {})
    flat['metabolomics_enabled'] = meta.get('enabled', False)
    flat['metabolomics_dir'] = meta.get('dir', '')
    flat['metabolomics_maf_files'] = meta.get('maf_files', [])
    flat['metabolomics_sample_file'] = meta.get('sample_file', '')
    flat['metabolomics_rpmi_samples'] = meta.get('rpmi_samples', [])
    flat['metabolomics_sera_samples'] = meta.get('sera_samples', [])

    sw = px.get('swath', {})
    flat['prot_swath_enabled'] = sw.get('enabled', True)
    flat['prot_swath_rpmi_template'] = sw.get('rpmi_template', '')
    flat['prot_swath_sera_template'] = sw.get('sera_template', '')

    return flat


def load_config(strain: str) -> dict:
    """Load configuration for a given strain.

    Resolution order:
      1. strains/<name>/config.json  (per-strain folder — recommended)
      2. multiomics_graph/strains_config.json  (legacy centralized config)
    """
    # Try per-strain folder first
    if _find_strain_dir(strain):
        return _load_strain_from_folder(strain)

    # Fall back to centralized config
    if STRAINS_CONFIG_PATH.exists():
        with open(STRAINS_CONFIG_PATH, 'r') as f:
            all_configs = json.load(f)
        if strain in all_configs:
            return _load_strain_from_central_config(strain)

    available = list_available_strains()
    raise ValueError(
        f"Unknown strain '{strain}'. Available: {available}\n\n"
        f"To add a new strain:\n"
        f"  1. Create strains/<name>/config.json (see strains/TEMPLATE/)\n"
        f"  2. Place data in strains/<name>/genome/, transcriptomic/, proteomics/\n"
        f"  3. Run: python main.py --strain <name>"
    )


# ──────────────────────────────────────────────
# Strain processing pipeline
# ──────────────────────────────────────────────

class StrainPipeline:
    """
    Process a single bacterial strain through the full multi-omics pipeline.

    Each stage produces intermediate outputs and passes data to the next.
    """

    def __init__(self, strain: str, config: dict,
                 output_base: str = 'outputs',
                 skip_annotation: bool = False,
                 no_cache: bool = False,
                 correlation_threshold: float = 0.7,
                 correlation_top_k: int = 10):
        self.strain = strain
        self.strain_dir = STRAINS_DIR / strain
        self.config = config
        self.skip_annotation = skip_annotation
        self.no_cache = no_cache
        self.correlation_threshold = correlation_threshold
        self.correlation_top_k = correlation_top_k

        # Output directories
        self.out_dir = Path(output_base) / strain
        self.fig_dir = self.out_dir / 'figures'
        self.cache_dir = self.out_dir / 'cache'
        self.graph_dir = self.out_dir / 'graph'
        for d in [self.out_dir, self.fig_dir, self.cache_dir, self.graph_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # If --no-cache, wipe cache directories
        if self.no_cache:
            import shutil
            cache_root = Path(output_base) / 'cache'
            if cache_root.exists():
                shutil.rmtree(cache_root)
            cache_root.mkdir(parents=True, exist_ok=True)
            for p in Path(output_base).rglob('*.pt'):
                p.unlink()
            print("  [--no-cache] All cached data cleared")

        # State: populated by each stage
        self.genome = None
        self.rna = None
        self.rna_raw = None
        self.has_rna = False
        self.protein = None
        self.protein_raw = None
        self.metabolomics = None
        self.metabolite_kegg = {}
        self.metabolite_to_pathway = {}
        self.metabolite_to_enzymes = {}
        self.aligned = None
        self.annotations = {}
        self._uniprot_annotator = None
        self.ppi_edges = None
        self.pathway_membership = {}
        self.cog_annotations = {}
        self.gene_protein_map = {}
        self.hetero_data = None
        self.hypergraph = None
        self.visualizer = NetworkVisualizer(str(self.fig_dir))

        self.timing = {}

        # Reproducible RNG across the pipeline (recorded in manifest.json).
        self.seed = PIPELINE_SEED
        np.random.seed(self.seed)
        torch.manual_seed(self.seed)

    def _log(self, stage: str, msg: str):
        elapsed = time.time() - self._stage_start
        print(f"  [{self.strain}] [{stage}] {msg} ({elapsed:.1f}s)")

    def _start_stage(self, name: str):
        self._stage_start = time.time()
        if name not in self.timing:
            self.timing[name] = 0.0
        print(f"\n{'=' * 50}")
        print(f"  {self.strain}: {name}")
        print(f"{'=' * 50}")

    def _end_stage(self, name: str):
        elapsed = time.time() - self._stage_start
        self.timing[name] = elapsed

    # ── Stage 1: Genomics ────────────────────────

    def run_genomics(self):
        """Parse genome GFF → gene table with GeneID ↔ ProteinID mapping."""
        self._start_stage("Genomics")

        gff = self.config['genome_gff']
        if not gff.exists():
            raise FileNotFoundError(f"GFF not found: {gff}")

        self.genome = GenomicsProcessor(strain=self.strain)
        self.genome.parse_gff(str(gff))

        s = self.genome.summary()
        self._log("Genomics",
                  f"{s['genes']} CDS, {s['replicons']} replicons, "
                  f"{s['with_gene_symbol']} with symbols")

        # Save gene table
        self.genome.gene_table.to_csv(
            self.out_dir / 'genome_genes.csv', index=False
        )

        # Compute genomic proximity edges
        genomic_edges = self.genome.compute_genomic_edges()
        genomic_edges.to_csv(self.out_dir / 'genome_genomic_edges.csv',
                              index=False)
        self._log("Genomics",
                  f"{len(genomic_edges)} genomic proximity edges")

        self._end_stage("Genomics")

    # ── Stage 2: Transcriptomics ─────────────────

    def run_transcriptomics(self):
        """Process RNA expression → condition means + logFC."""
        self._start_stage("Transcriptomics")

        # Prefer the log2 CPM matrix; raw counts are not comparable across
        # samples, so logFC/regulation must be computed from log2 CPM.
        rna_file = self.config.get('rna_file') or self.config.get('rna_counts_file')
        if not rna_file or str(rna_file) == '.' or not rna_file.exists() or not rna_file.is_file():
            self._log("Transcriptomics", "No transcriptomics file; skipping stage")
            self.rna = None
            self.has_rna = False
            return

        self.rna = TranscriptomicsProcessor(
            strain=self.strain,
            std_multiplier=self.config.get('regulation_std_multiplier', 0.0),
        )
        self.rna.load_log2_cpm(str(rna_file))
        self.has_rna = True

        self.rna.compute_condition_means(
            rpmi_samples=self.config['rna_rpmi_samples'],
            sera_samples=self.config['rna_sera_samples'],
        )

        s = self.rna.summary()
        self._log("Transcriptomics",
                  f"{s['genes']} genes, {s['up_regulated']} up, "
                  f"{s['down_regulated']} down")

        # Save
        self.rna.expression_table.to_csv(
            self.out_dir / 'transcriptomics_expression.csv', index=False
        )

        self._end_stage("Transcriptomics")

    # ── Stage 3: Proteomics ──────────────────────

    def run_proteomics(self):
        """Process SWATH-MS proteomics data to condition means + logFC.

        Loads the SWATH-MS Excel `Area - proteins` sheet. Column naming
        templates are read from the strain config.
        """
        self._start_stage("Proteomics")

        cfg = self.config
        rpmi_ids = cfg.get('prot_rpmi_ids', [])
        sera_ids = cfg.get('prot_sera_ids', [])

        if not rpmi_ids or not sera_ids:
            self._log("Proteomics", "No proteomics sample IDs configured, skipping")
            return

        proto_df = None
        prot_file = cfg.get('prot_file')

        # ── SWATH-MS Excel ──
        if prot_file and prot_file.exists() and cfg.get('prot_swath_enabled', True):
            try:
                self._log("Proteomics", "Loading SWATH-MS Excel...")
                df = pd.read_excel(
                    str(prot_file),
                    sheet_name=cfg.get('prot_sheet', 'Area - proteins'),
                )

                protein_col = cfg.get('prot_protein_col', 'Protein')
                rpmi_template = cfg.get('prot_swath_rpmi_template', '')
                sera_template = cfg.get('prot_swath_sera_template', '')

                rpmi_cols = [rpmi_template.format(sid=sid) for sid in rpmi_ids]
                sera_cols = [sera_template.format(sid=sid) for sid in sera_ids]

                available = set(df.columns)
                rpmi_cols = [c for c in rpmi_cols if c in available]
                sera_cols = [c for c in sera_cols if c in available]

                if rpmi_cols and sera_cols:
                    proto_df = pd.DataFrame()
                    proto_df['ProteinID'] = df[protein_col].str.replace(
                        '^ref\\|', '', regex=True
                    )
                    for col in rpmi_cols + sera_cols:
                        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

                    proto_df['Protein_RPMI'] = np.log2(df[rpmi_cols].mean(axis=1) + 1)
                    proto_df['Protein_Sera'] = np.log2(df[sera_cols].mean(axis=1) + 1)
                    proto_df['Protein_logFC'] = proto_df['Protein_Sera'] - proto_df['Protein_RPMI']

                    self._log("Proteomics", f"SWATH: {len(proto_df)} protein groups")
            except Exception as e:
                self._log("Proteomics", f"SWATH failed: {e}")

        if proto_df is None or proto_df.empty:
            self._log("Proteomics", "No proteomics data loaded")
            return

        # Aggregate by ProteinID (some appear in multiple rows)
        proto_agg = proto_df.groupby('ProteinID', as_index=False).agg({
            'Protein_RPMI': 'mean',
            'Protein_Sera': 'mean',
            'Protein_logFC': 'mean',
        })

        # Regulation call — std-based or absolute threshold
        std_mult = cfg.get('regulation_std_multiplier', 0.0)
        if std_mult > 0:
            pthresh = std_mult * proto_agg['Protein_logFC'].std()
        else:
            pthresh = 1.0
        proto_agg['Protein_Regulation'] = np.select(
            [
                proto_agg['Protein_logFC'] > pthresh,
                proto_agg['Protein_logFC'] < -pthresh,
            ],
            ['Up', 'Down'],
            default='Stable',
        )

        # Store as simple namespace
        reg_counts = proto_agg['Protein_Regulation'].value_counts()
        self.protein = type('obj', (object,), {
            'protein_table': proto_agg,
        })
        # Retain the per-sample intensity matrix (ProteinID x replicates) for
        # the protein abundance-correlation edge family.
        self.protein_raw = None
        try:
            if rpmi_cols and sera_cols and 'ProteinID' in proto_df.columns:
                prot_repl = proto_df[['ProteinID']].copy()
                for _col in rpmi_cols + sera_cols:
                    if _col in df.columns:
                        prot_repl[_col] = df[_col].values
                self.protein_raw = (
                    prot_repl.groupby('ProteinID', as_index=False)
                    .first()
                    .set_index('ProteinID')
                )
        except Exception as e:
            self._log("Proteomics", f"raw replicate matrix skipped: {e}")
        self.protein._summary_cache = {
            'proteins': len(proto_agg),
            'up_regulated': int(reg_counts.get('Up', 0)),
            'down_regulated': int(reg_counts.get('Down', 0)),
            'stable': int(reg_counts.get('Stable', 0)),
        }

        proto_agg.to_csv(self.out_dir / 'proteomics_abundance.csv', index=False)
        self._log("Proteomics", f"{len(proto_agg)} proteins quantified")
        self._end_stage("Proteomics")

    # ── Stage 3.5: Metabolomics ───────────────────

    def run_metabolomics(self):
        """Process metabolite abundance (MetaboLights MAF) → condition means + logFC.

        Assigns each metabolite sample column to RPMI or Sera based on the
        strain config, computes mean abundances and log2 fold change, then
        resolves ChEBI identifiers to KEGG compounds (enzymes + pathways).
        """
        self._start_stage("Metabolomics")

        cfg = self.config
        if not cfg.get('metabolomics_enabled', False):
            self._log("Metabolomics", "Disabled for this strain, skipping")
            self._end_stage("Metabolomics")
            return

        from preprocessing.metabolomics import MetabolomicsProcessor
        from annotation.kegg_compound import KEGGCompoundResolver

        maf_files = cfg.get('metabolomics_maf_files', [])
        if not maf_files:
            meta_dir = cfg.get('metabolomics_dir', '')
            if meta_dir:
                d = Path(meta_dir)
                maf_files = [str(f) for f in sorted(d.glob('*maf.tsv'))]
            else:
                self._log("Metabolomics", "No MAF files configured, skipping")
                self._end_stage("Metabolomics")
                return

        proc = MetabolomicsProcessor(
            maf_files,
            sample_file=cfg.get('metabolomics_sample_file') or None,
        )
        try:
            table = proc.load(
                rpmi_samples=cfg.get('metabolomics_rpmi_samples', []),
                sera_samples=cfg.get('metabolomics_sera_samples', []),
                min_detection=0.0,
            )
        except ValueError as e:
            self._log("Metabolomics", f"Processing failed: {e}")
            self._end_stage("Metabolomics")
            return

        table.to_csv(self.out_dir / 'metabolomics_abundance.csv', index=False)
        self.metabolomics = type('obj', (object,), {'table': table})
        self.metabolomics._summary_cache = proc.summary()
        self.metabolomics.raw_replicates = getattr(
            proc, 'raw_replicates', pd.DataFrame()
        )

        # Resolve ChEBI -> KEGG compound -> enzymes/pathways
        chebis = [
            mid for mid in table['MetaboliteID'].tolist()
            if isinstance(mid, str) and mid.startswith('CHEBI:')
        ]
        resolver = KEGGCompoundResolver()
        resolved = resolver.resolve(chebis) if chebis else {}
        self.metabolite_kegg = resolved
        self.metabolite_to_pathway = {}
        self.metabolite_to_enzymes = {}

        # Name-based KEGG fallback for ChEBI-less / unresolved metabolites.
        # Keys returned by find_by_name are metabolite names; map them back
        # to MetroboliteIDs so the graph edges use the canonical node ids.
        unresolved = table[
            (~table['MetaboliteID'].isin(set(resolved.keys())))
        ]
        name_records = resolver.find_by_name(
            unresolved['MetaboliteName'].astype(str).str.strip().tolist()
        ) if not unresolved.empty else {}
        name_to_id = {
            str(row['MetaboliteName']).strip(): row['MetaboliteID']
            for _, row in unresolved.iterrows()
        }

        def _apply(info, target_id):
            if info.get('pathways'):
                self.metabolite_to_pathway[target_id] = list(info['pathways'].keys())
            if info.get('enzymes'):
                self.metabolite_to_enzymes[target_id] = info['enzymes']

        for che, info in resolved.items():
            _apply(info, che)
        for nm, info in name_records.items():
            tid = name_to_id.get(nm)
            if tid:
                _apply(info, tid)

        n_resolved = len(resolved) + len(name_records)
        s = proc.summary()
        self._log("Metabolomics",
                  f"{s['metabolites']} metabolites, {n_resolved} KEGG-resolved "
                  f"({len(resolved)} ChEBI, {len(name_records)} by name), "
                  f"{s['up_regulated']} up, {s['down_regulated']} down")
        self._end_stage("Metabolomics")

    # ── Stage 4: Integration ─────────────────────

    def run_integration(self):
        """Align gene, RNA, and protein data via the central dogma bridge.

        Requires the genome. If transcriptomics (RNA) is unavailable, falls
        back to a 2-layer (genome + proteomics) alignment so the graph can
        still be built from the available omics layers.
        """
        self._start_stage("Integration")

        if self.genome is None:
            self._log("Integration", "Missing genome data")
            return

        # Build the bridge: ProteinID → GeneID from genome
        prot_to_gene = {}
        for _, row in self.genome.gene_table.iterrows():
            pid = row.get('ProteinID', '')
            if pid and not pd.isna(pid):
                prot_to_gene[pid] = row['GeneID']

        # Deduplicate genome data (some GeneIDs may have multiple rows)
        genome_dedup = self.genome.gene_table.groupby('GeneID').first().reset_index()
        genome_data = genome_dedup.set_index('GeneID')
        protein_data = self.protein.protein_table.copy()

        # Map protein data to GeneID
        if 'GeneID' not in protein_data.columns:
            protein_data['GeneID'] = protein_data['ProteinID'].map(prot_to_gene)

        # Drop unmapped proteins
        protein_data = protein_data.dropna(subset=['GeneID'])
        protein_data = protein_data.set_index('GeneID')

        # Merge RNA + protein on GeneID (handle duplicate indices)
        prot_cols = ['ProteinID', 'Protein_RPMI', 'Protein_Sera', 'Protein_logFC']
        if 'Protein_Regulation' in protein_data.columns:
            prot_cols.append('Protein_Regulation')

        rna_available = self.has_rna and self.rna is not None and self.rna.raw_data is not None
        if rna_available:
            rna_data = self.rna.expression_table.set_index('GeneID')
            aligned = rna_data[['RNA_RPMI', 'RNA_Sera', 'RNA_logFC']].join(
                protein_data[prot_cols],
                how='inner',
            )
        else:
            # 2-layer fallback: align proteomics against the genome only.
            aligned = protein_data[prot_cols].copy()
            for col in ['RNA_RPMI', 'RNA_Sera', 'RNA_logFC']:
                aligned[col] = float('nan')
            self._log("Integration",
                      "No RNA data; building 2-layer (genome + proteome) alignment")

        aligned = aligned[~aligned.index.duplicated(keep='first')]

        # Add genome features
        for col in ['start', 'end', 'strand', 'CDS_length', 'replicon']:
            if col in genome_data.columns:
                aligned[col] = genome_data.loc[aligned.index, col]

        aligned = aligned.reset_index()
        aligned = aligned.rename(columns={'index': 'GeneID'})

        # Build gene → protein mapping for graph construction
        for _, row in aligned.iterrows():
            gid = row['GeneID']
            pid = row.get('ProteinID', '')
            if pid and not pd.isna(pid):
                self.gene_protein_map[gid] = pid

        self.aligned = aligned
        aligned.to_csv(self.out_dir / 'aligned_multiomics.csv', index=False)

        n_layers = "3 omics" if rna_available else "2 omics (no RNA)"
        self._log("Integration",
                  f"{len(aligned)} genes aligned across {n_layers} layers")
        self._end_stage("Integration")

    # ── Stage 5: Biological Annotation ───────────

    def run_annotation(self):
        """
        Retrieve biological annotations (can be skipped for offline use).

        UniProt:    functional descriptions, GO terms
        STRING:     protein-protein interactions
        KEGG:       pathway membership
        eggNOG:     COG functional categories
        """
        if self.skip_annotation:
            self._log("Annotation", "Skipping (offline mode)")
            # Create empty annotations for hypergraph structure
            if self.aligned is not None:
                self.cog_annotations = {
                    gid: {'COG_classes': ['S'], 'COG_category': 'S', 'Description': ''}
                    for gid in self.aligned['GeneID']
                }
                # Create empty pathway membership from COG
                for gid in self.aligned['GeneID']:
                    self.pathway_membership.setdefault('S', []).append(gid)
            return

        self._start_stage("Annotation")

        protein_ids = self.aligned['ProteinID'].dropna().unique().tolist() if self.aligned is not None else []

        # UniProt
        if protein_ids:
            uniprot = UniProtAnnotator(organism_id=str(self.config.get('ncbi_tax_id', 562)))
            uni_df = uniprot.fetch_annotations(protein_ids)
            if not uni_df.empty:
                uni_df.to_csv(self.out_dir / 'uniprot_annotations.csv', index=False)
                self.annotations['uniprot'] = uni_df
                self._uniprot_annotator = uniprot

        # STRING — use UniProt accessions with species fallback (primary UPEC strain → K-12)
        if protein_ids:
            uniprot_map = uniprot.map_refseq_to_uniprot(protein_ids)
            uniprot_ids = list(set(uniprot_map.values())) if uniprot_map else []
            if uniprot_ids:
                primary_sp = self.config.get(
                    'string_species_primary',
                    self.config.get('string_species', 511145)
                )
                fallback_sp = self.config.get(
                    'string_species_fallback',
                    self.config.get('string_species', 511145)
                )

                # Collect PPI edges from primary species then fallback
                all_ppi = []
                for sp in [primary_sp, fallback_sp]:
                    client = STRINGClient(species_id=sp)
                    ppi_part = client.fetch_interactions(uniprot_ids)
                    if not ppi_part.empty:
                        all_ppi.append(ppi_part)

                ppi = pd.concat(all_ppi, ignore_index=True) if all_ppi else pd.DataFrame()
                if not ppi.empty:
                    # Deduplicate across merged results (same pair may appear in both species)
                    ppi['_pair_key'] = ppi.apply(
                        lambda r: tuple(sorted([str(r['ProteinID_A']), str(r['ProteinID_B'])])), axis=1
                    )
                    ppi = ppi.loc[
                        ppi.groupby('_pair_key')['combined_score'].idxmax()
                    ].drop(columns=['_pair_key']).reset_index(drop=True)

                    # Map UniProt accessions back to RefSeq WP_ for graph node IDs
                    uniprot_to_wp = {v: k for k, v in uniprot_map.items()}
                    unver_to_ver = {}
                    for pid in protein_ids:
                        unver = pid.split('.')[0]
                        if unver not in unver_to_ver:
                            unver_to_ver[unver] = pid
                    if 'ProteinID_A' in ppi.columns:
                        ppi['ProteinID_A'] = ppi['ProteinID_A'].map(uniprot_to_wp).map(unver_to_ver).fillna(ppi['ProteinID_A'])
                        ppi['ProteinID_B'] = ppi['ProteinID_B'].map(uniprot_to_wp).map(unver_to_ver).fillna(ppi['ProteinID_B'])
                    ppi.to_csv(self.out_dir / 'string_ppi.csv', index=False)
                    self.ppi_edges = ppi
                    self._log("STRING", f"{len(ppi)} PPI edges (primary={primary_sp}, fallback={fallback_sp})")
            else:
                self._log("STRING", "No UniProt mappings available for STRING query")

        # KEGG — use UniProt cross-references to bridge locus tags to KEGG gene IDs
        if self.aligned is not None:
            kegg = KEGGAnnotator(organism_code=self.config.get('kegg_code', 'eco'), strain_name=self.strain)
            gene_ids = self.aligned['GeneID'].tolist()
            # Build gene → KEGG ID map from UniProt cross-references
            gene_to_kegg = {}
            if hasattr(uniprot, 'mapping_stats') and uniprot.mapping_stats and uniprot.mapping_stats.mapped > 0:
                uniprot_map = uniprot.map_refseq_to_uniprot(protein_ids)
                if self.gene_protein_map:
                    for gid, wp_list in self.gene_protein_map.items():
                        for wp in (wp_list if isinstance(wp_list, list) else [wp_list]):
                            wp_norm = wp.split('.')[0].strip()
                            up_id = uniprot_map.get(wp_norm)
                            if up_id and 'uniprot' in self.annotations:
                                uni_row = self.annotations['uniprot']
                                row = uni_row[uni_row['UniProtID'] == up_id]
                                if not row.empty and 'KEGG' in row.columns:
                                    kegg_ids = str(row.iloc[0]['KEGG'])
                                    for kid in kegg_ids.split('; '):
                                        if kid.startswith(f"{kegg.organism_code}:"):
                                            gene_to_kegg[gid] = kid.split(':')[1]
                                            break
            pw_membership = kegg.build_pathway_membership(gene_ids, gene_to_kegg_map=gene_to_kegg if gene_to_kegg else None)
            if pw_membership:
                self.pathway_membership = pw_membership
                self._log("KEGG", f"{len(pw_membership)} pathways mapped")

        # eggNOG
        if self.aligned is not None:
            eggnog = EggNOGAnnotator()
            eggnog_file = self.strain_dir / 'eggnog_annotations.tsv'
            if eggnog_file.exists():
                eggnog.load_precomputed(str(eggnog_file))
                self._log("eggNOG", f"loaded {len(eggnog.annotations)} annotations from {eggnog_file.name}")
            else:
                gene_ids = self.aligned['GeneID'].tolist()
                eggnog.annotate_by_locus_tag(gene_ids)
                self._log("eggNOG", f"{len(eggnog.annotations)} genes annotated (placeholder)")

            # Re-key annotations from WP accessions onto locus tags via the
            # gene->protein mapping, so COG edges/hyperedges attach to real
            # gene nodes instead of orphan WP nodes.
            resolved = eggnog.as_locus_tag_map(self.gene_protein_map)
            if resolved:
                self.cog_annotations = resolved
                self._log("eggNOG",
                          f"{len(resolved)} COG annotations resolved onto graph gene nodes")
            else:
                self.cog_annotations = eggnog.annotations

            # eggNOG KEGG pathway membership: acquire/plasmid-borne AMR genes
            # are invisible to the `eco`-bridged UniProt cross-reference route,
            # but the emapper output carries KEGG_Pathway for them. Merge those
            # memberships so the graph gets a KEGG layer for these determinants.
            try:
                eggnog_kegg = eggnog.kegg_pathway_membership(self.gene_protein_map)
                if eggnog_kegg:
                    if not self.pathway_membership:
                        self.pathway_membership = {}
                    n_added = 0
                    for pw, genes in eggnog_kegg.items():
                        pw_existing = self.pathway_membership.setdefault(pw, [])
                        for g in genes:
                            if g not in pw_existing:
                                pw_existing.append(g)
                                n_added += 1
                    self._log("eggNOG KEGG",
                              f"{len(eggnog_kegg)} pathways, {n_added} gene memberships added "
                              f"from eggNOG KEGG_Pathway")
            except Exception as e:
                self._log("eggNOG KEGG", f"skip: {e}")

        # AMR knowledge layer (curated manifest, genome-derived)
        try:
            self.amr_manifest = AMRManifest()
            self.amr_summary = self.amr_manifest.summary(self.strain)
            n_amr = self.amr_summary['n_markers']
            self._log("AMR", f"{n_amr} resistance markers in manifest "
                             f"({self.amr_summary['n_loci']} loci, "
                             f"{len(self.amr_summary['classes'])} classes)")
        except FileNotFoundError:
            self.amr_manifest = None
            self.amr_summary = {}
            self._log("AMR", "no manifest, skipping AMR annotation")

        self._end_stage("Annotation")

    # ── Stage 6: Heterogeneous Graph ─────────────

    def run_graph_construction(self):
        """Build PyG HeteroData heterogeneous graph."""
        self._start_stage("Graph Construction")

        if self.aligned is None or self.aligned.empty:
            self._log("Graph", "No aligned data available")
            return

        # Gene node universe = genome-annotated genes ∪ multi-omics aligned
        # genes. Genome-only AMR determinants (e.g. tet(A), aadA5,
        # blaCTX-M-15) are NOT in `aligned` (no RNA/protein measurement) but
        # must still appear as gene nodes so they carry encodes/COG/KEGG/
        # proximity edges and connect to the broader molecular state.
        gene_ids = set(self.aligned['GeneID'].tolist())
        if self.genome is not None and self.genome.gene_table is not None \
                and 'GeneID' in self.genome.gene_table.columns:
            gene_ids |= set(self.genome.gene_table['GeneID'].dropna().astype(str))
        gene_ids = sorted(gene_ids)
        protein_ids = self.aligned['ProteinID'].dropna().unique().tolist()

        # Prepare indexed dataframes
        rna_idx = self.aligned.set_index('GeneID')[
            ['RNA_RPMI', 'RNA_Sera', 'RNA_logFC']
        ]

        # For the RNA correlation edges, we need per-replicate data
        # Re-load the raw counts for replicate-level access
        rna_raw = rna_idx
        if self.rna is not None and getattr(self.rna, 'raw_data', None) is not None \
                and not self.rna.raw_data.empty:
            raw = self.rna.raw_data
            gcol = getattr(self.rna, '_gene_id_col', None) or \
                ('GeneID' if 'GeneID' in raw.columns else None)
            if gcol and gcol in raw.columns:
                rna_raw = raw.set_index(gcol)

        # For the protein abundance-correlation edges, use the per-sample
        # intensity matrix retained by run_proteomics (ProteinID x replicates).
        protein_raw = getattr(self, 'protein_raw', None)

        genome_idx = self.genome.gene_table.groupby('GeneID').first()[
            ['start', 'end', 'strand', 'CDS_length']
        ] if self.genome else pd.DataFrame(index=gene_ids)

        protein_idx = self.aligned.set_index('ProteinID')[
            ['Protein_RPMI', 'Protein_Sera', 'Protein_logFC']
        ] if protein_ids else pd.DataFrame()

        # Genomic edges (filtered to graph gene nodes)
        all_edges = self.genome.compute_genomic_edges() if self.genome else pd.DataFrame()
        aligned_gene_set = set(gene_ids)
        if not all_edges.empty:
            all_edges = all_edges[
                all_edges['GeneID_A'].isin(aligned_gene_set) &
                all_edges['GeneID_B'].isin(aligned_gene_set)
            ]

        # Build graph
        builder = HeterogeneousGraphBuilder(
            correlation_threshold=self.correlation_threshold,
            correlation_top_k=self.correlation_top_k,
        )

        # Pathway annotations as dict for graph
        gene_to_ko = {}
        if self.pathway_membership:
            for pw, genes in self.pathway_membership.items():
                for g in genes:
                    gene_to_ko.setdefault(g, []).append(pw)

        ko_names = {pw: pw.replace('path:', '') for pw in self.pathway_membership}

        # Convert cog_annotations {GeneID: {COG_classes: [...], ...}} -> {GeneID: [class, ...]}
        gene_to_cog_map = {}
        if self.cog_annotations:
            for gid, ann in self.cog_annotations.items():
                classes = ann.get('COG_classes', [])
                if classes:
                    gene_to_cog_map[gid] = classes

        # Metabolomics layer inputs
        metabolite_ids = None
        metabolite_data = None
        metabolite_replicates = None
        if getattr(self.metabolomics, 'table', None) is not None and not self.metabolomics.table.empty:
            mtable = self.metabolomics.table
            metabolite_ids = mtable['MetaboliteID'].tolist()
            metabolite_data = mtable.set_index('MetaboliteID')
            metabolite_replicates = getattr(
                self.metabolomics, 'raw_replicates', None
            )

        # pathway id -> name map for metabolite pathway edges (from KEGG resolver)
        pathway_names = {}
        for info in self.metabolite_kegg.values():
            for pid, pname in info.get('pathways', {}).items():
                pathway_names.setdefault(pid, pname)

        # enzyme (EC) -> proteins map from UniProt EC annotations.
        # UniProt uses unversioned WP ids while graph nodes are versioned,
        # so translate unversioned -> first versioned protein id.
        unver_to_ver = {}
        for p in protein_ids:
            unver = p.split('.')[0]
            if unver not in unver_to_ver:
                unver_to_ver[unver] = p
        enzyme_to_protein = {}
        if 'uniprot' in self.annotations:
            uni = self.annotations['uniprot']
            if 'EC_number' in uni.columns:
                for _, row in uni.iterrows():
                    pid = row.get('ProteinID', '')
                    ec = str(row.get('EC_number', '')).strip()
                    if not pid or not ec or ec == 'nan':
                        continue
                    pid = unver_to_ver.get(pid.split('.')[0], pid)
                    for ec_entry in ec.split(';'):
                        ec_entry = ec_entry.strip()
                        if ec_entry:
                            enzyme_to_protein.setdefault(ec_entry, []).append(pid)

        # Prepare GO annotations for HeteroData
        protein_to_go = {}
        go_info = {}
        if 'uniprot' in self.annotations:
            uni = self.annotations['uniprot']
            for _, row in uni.iterrows():
                pid = row.get('ProteinID', '')
                if not pid:
                    continue
                pid = unver_to_ver.get(pid.split('.')[0], pid)
                for col in ('GO_BP', 'GO_MF', 'GO_CC'):
                    raw = row.get(col, '')
                    if not isinstance(raw, str) or not raw.strip():
                        continue
                    for go in raw.split(';'):
                        go = go.strip()
                        if go:
                            protein_to_go.setdefault(pid, []).append(go)
            if getattr(self._uniprot_annotator, 'go_terms', None):
                for go_id, gt in self._uniprot_annotator.go_terms.items():
                    go_info[go_id] = {'name': gt.name, 'aspect': gt.aspect}

            # Merge recovered per-ID UniProt annotations for AMR determinants
            # (see scripts/recover_amr_uniprot.py). The upstream batch mapping
            # truncated >500-hit batches, dropping AMR proteins; the recovery
            # re-runs the same search one ID at a time and records GO terms in
            # reports/amr_uniprot_recovery.csv. Merge those GO terms into the
            # graph's protein->GO map so graph GO edges/hyperedges cover the
            # determinants, not just the trace layer.
            recovery_csv = Path(__file__).resolve().parent.parent \
                / 'reports' / 'amr_uniprot_recovery.csv'
            if recovery_csv.exists():
                try:
                    rec_df = pd.read_csv(recovery_csv, dtype=str)
                    rec_n = 0
                    for _, rec in rec_df.iterrows():
                        pid = str(rec.get('protein_id', '') or '').strip()
                        go_raw = str(rec.get('go_ids', '') or '')
                        if not pid or not go_raw or go_raw.lower() == 'nan':
                            continue
                        # Resolve unversioned recovery id onto a versioned
                        # graph protein node when possible.
                        pid_ver = unver_to_ver.get(pid.split('.')[0], pid)
                        for go in go_raw.split(';'):
                            go = go.strip()
                            if not go:
                                continue
                            if go not in protein_to_go.get(pid_ver, []):
                                protein_to_go.setdefault(pid_ver, []).append(go)
                            rec_n += 1
                    if rec_n:
                        self._log("GO recovery",
                                  f"{rec_n} recovered GO terms merged into protein_to_go "
                                  f"from amr_uniprot_recovery.csv")
                except Exception as e:
                    self._log("GO recovery", f"skip: {e}")

        self.protein_to_go = protein_to_go if protein_to_go else None
        self.go_info = go_info if go_info else None

        # Prepare AMR mechanisms for HeteroData
        class_to_loci = None
        class_names = None
        if getattr(self, 'amr_manifest', None) is not None:
            class_to_loci = self.amr_manifest.class_to_loci(self.strain)
            class_names = AMR_CLASS_NAMES

        # Prepare PPI clusters (Louvain communities) for HeteroData
        ppi_clusters = []
        if self.ppi_edges is not None and not self.ppi_edges.empty:
            import networkx as nx
            G_ppi = nx.Graph()
            for _, row in self.ppi_edges.iterrows():
                a, b = row.get('ProteinID_A', ''), row.get('ProteinID_B', '')
                w = row.get('combined_score', 1.0)
                if a and b:
                    G_ppi.add_edge(a, b, weight=w)
            if G_ppi.number_of_nodes() > 0:
                ppi_clusters = [list(c) for c in nx.community.louvain_communities(G_ppi, weight='weight', resolution=1.0, seed=42) if len(c) >= 3]

        self.hetero_data = builder.build(
            gene_ids=gene_ids,
            protein_ids=protein_ids,
            rna_data=rna_idx,
            genome_data=genome_idx,
            protein_data=protein_idx,
            genomic_edges=all_edges,
            gene_protein_map=self.gene_protein_map,
            ppi_edges=self.ppi_edges,
            rna_replicates=rna_raw,
            protein_replicates=protein_raw,
            gene_to_ko=gene_to_ko if gene_to_ko else None,
            ko_to_pathway=ko_names if ko_names else None,
            gene_to_cog=gene_to_cog_map if gene_to_cog_map else None,
            metabolite_ids=metabolite_ids,
            metabolite_data=metabolite_data,
            metabolite_to_pathway=self.metabolite_to_pathway or None,
            pathway_names=pathway_names or None,
            metabolite_to_enzymes=self.metabolite_to_enzymes or None,
            enzyme_to_protein=enzyme_to_protein or None,
            metabolite_replicates=metabolite_replicates,
            protein_to_go=protein_to_go if protein_to_go else None,
            go_info=go_info if go_info else None,
            class_to_loci=class_to_loci if class_to_loci else None,
            class_names=class_names if class_names else None,
            ppi_clusters=ppi_clusters if ppi_clusters else None,
        )

        builder.save(str(self.graph_dir))

        # Count edge types
        edge_types = {}
        edge_keys = self.hetero_data.edge_types if hasattr(self.hetero_data, 'edge_types') else \
                    [k for k in self.hetero_data if isinstance(k, tuple)] if isinstance(self.hetero_data, dict) else []
        for key in edge_keys:
            store = self.hetero_data[key]
            ei = store.get('edge_index', torch.zeros((2, 0))) if isinstance(store, dict) else \
                 getattr(store, 'edge_index', torch.zeros((2, 0)))
            edge_types[key[1]] = ei.size(1)

        self._log("Graph",
                  f"{len(gene_ids)} gene nodes, {len(protein_ids)} protein nodes"
                  + (f", {len(metabolite_ids)} metabolite nodes" if metabolite_ids else ""))
        self._end_stage("Graph Construction")

        # Return for summary
        return builder.summary()

    # ── Stage 7: Hypergraph ──────────────────────

    def run_hypergraph(self):
        """Build biological module hypergraph."""
        self._start_stage("Hypergraph Construction")

        if self.aligned is None or self.aligned.empty:
            self._log("Hypergraph", "No aligned data available")
            return

        hg = BiologicalHypergraphBuilder()

        # 1. Multi-omics triples (central dogma hyperedges)
        hg.add_multi_omics_triples(self.aligned)
        self._log("Hypergraph",
                  f"{len([h for h in hg.hyperedges if h['type'] == 'multi_omics_triple'])} "
                  f"multi-omics triple hyperedges")

        # 2. KEGG pathway hyperedges
        if self.pathway_membership:
            hg.add_kegg_pathway_hyperedges(self.pathway_membership)
            n_path = len([h for h in hg.hyperedges if h['type'] == 'kegg_pathway'])
            self._log("Hypergraph", f"{n_path} KEGG pathway hyperedges")

        # 3. COG category hyperedges
        if self.cog_annotations:
            hg.add_cog_hyperedges(
                self.cog_annotations,
                gene_to_protein=self.gene_protein_map,
                cog_categories=COG_CATEGORIES,
            )
            n_cog = len([h for h in hg.hyperedges if h['type'] == 'cog_category'])
            self._log("Hypergraph", f"{n_cog} COG category hyperedges")

        # 3b. AMR mechanism hyperedges (curated resistance knowledge)
        if getattr(self, 'amr_manifest', None) is not None:
            class_to_loci = self.amr_manifest.class_to_loci(self.strain)
            amr_loci = self.amr_manifest.amr_loci(self.strain)
            if class_to_loci:
                hg.add_amr_hyperedges(
                    class_to_loci,
                    amr_loci,
                    gene_to_protein=self.gene_protein_map,
                    class_names=AMR_CLASS_NAMES,
                )
                n_amr = len([h for h in hg.hyperedges if h['type'] == 'amr_mechanism'])
                self._log("Hypergraph",
                          f"{n_amr} AMR mechanism hyperedges "
                          f"({len(amr_loci)} AMR loci)")

        # 4. PPI cluster hyperedges
        if self.ppi_edges is not None and not self.ppi_edges.empty:
            protein_ids = self.aligned['ProteinID'].dropna().unique().tolist()
            hg.add_ppi_clusters(self.ppi_edges, protein_ids)
            n_ppi = len([h for h in hg.hyperedges if h['type'] == 'ppi_cluster'])
            self._log("Hypergraph", f"{n_ppi} PPI cluster hyperedges")

        # 4b. GO term hyperedges (from UniProt GO_BP / GO_MF annotations)
        # Reuse the versioned protein->GO map prepared during graph
        # construction so hyperedge node ids match the protein nodes.
        protein_to_go = getattr(self, 'protein_to_go', None) or {}
        go_info = getattr(self, 'go_info', None) or {}
        if protein_to_go:
            hg.add_go_hyperedges(protein_to_go, go_info)
            n_go = len([h for h in hg.hyperedges if h['type'] == 'go_term'])
            self._log("Hypergraph", f"{n_go} GO term hyperedges")

        # 5. Metabolite pathway hyperedges (co-join genes + metabolites)
        if self.metabolite_to_pathway:
            hg.add_metabolite_hyperedges(
                self.metabolite_to_pathway,
                pathway_genes=self.pathway_membership or None,
            )
            n_met = len([h for h in hg.hyperedges if h['type'] == 'metabolic_pathway'])
            self._log("Hypergraph", f"{n_met} metabolic pathway hyperedges")

        # Build incidence matrix
        hg.build_incidence_matrix()

        self.hypergraph = hg
        hg.save(str(self.graph_dir))

        self._end_stage("Hypergraph Construction")

    # ── Stage 8: Visualisation ───────────────────

    def run_visualization(self):
        """Generate all figures including a full network render."""
        self._start_stage("Visualization")

        viz = self.visualizer

        # 0. Full heterogeneous network graph — load from saved edge CSVs
        edge_csvs = list(self.graph_dir.glob('edges_*.csv'))
        if edge_csvs:
            try:
                import pandas as pd
                node_types = {'gene': set(), 'protein': set(), 'annotation': set(), 'metabolite': set()}
                edge_types = {}
                col_type_map = {
                    'gene': 'gene', 'protein': 'protein', 'annotation': 'annotation',
                    'metabolite': 'metabolite',
                    'gene_target': 'gene', 'protein_target': 'protein',
                }

                for csv_path in sorted(edge_csvs):
                    df = pd.read_csv(csv_path)
                    cols = [c for c in df.columns if c != 'relation']
                    if len(cols) < 2:
                        continue
                    src_col, dst_col = cols[0], cols[1]

                    # Edge type from filename: edges_<src>_<rel>_<dst>.csv
                    parts = csv_path.stem.replace('edges_', '', 1).split('_')
                    rel = '_'.join(parts[1:-1]) if len(parts) >= 3 else parts[0]

                    def get_type(c):
                        for prefix, t in col_type_map.items():
                            if c.startswith(prefix):
                                return t
                        return 'gene'

                    stype, dtype = get_type(src_col), get_type(dst_col)
                    key = (stype, rel, dtype)
                    edge_types[key] = []
                    for _, row in df.iterrows():
                        s, d = str(row[src_col]).strip(), str(row[dst_col]).strip()
                        if s and d:
                            s_id = f"{stype}_{s}"
                            d_id = f"{dtype}_{d}"
                            edge_types[key].append((s_id, d_id))
                            node_types[stype].add(s_id)
                            node_types[dtype].add(d_id)

                # Convert to lists for the viz function
                node_type_lists = {k: list(v) for k, v in node_types.items()}
                total_viz_nodes = sum(len(v) for v in node_type_lists.values())

                # Full render is only feasible up to a few thousand nodes
                # (networkx spring layout is O(V^2) per iteration); for larger
                # graphs the 500-node sampled figure below is used instead.
                if total_viz_nodes <= 2500:
                    viz.plot_heterogeneous_graph(
                        node_type_lists, edge_types,
                        title=f"{self.strain} — Heterogeneous Multi-Omics Graph",
                        filename="heterogeneous_graph.png",
                        max_nodes=0,
                    )
                else:
                    self._log(
                        "Visualization",
                        f"skipping full heterogeneous_graph.png "
                        f"({total_viz_nodes} nodes); using sampled render",
                    )
                # Also save a sampled version
                if total_viz_nodes > 1000:
                    viz.plot_heterogeneous_graph(
                        node_type_lists, edge_types,
                        title=f"{self.strain} — Heterogeneous Graph (sampled 500)",
                        filename="heterogeneous_graph_sampled.png",
                        max_nodes=500,
                    )

                # AMR highlight plot
                if getattr(self, 'amr_manifest', None) is not None \
                        and 'gene' in node_type_lists:
                    amr = self.amr_manifest
                    locus_to_amr = amr.locus_to_amr(self.strain)

                    # Skip when there are no curated determinants for this
                    # strain (e.g. MS_14385): an empty highlight would still
                    # force an O(V^2) spring layout on the full graph.
                    if not locus_to_amr:
                        self._log(
                            "Visualization",
                            "no AMR loci, skipping amr_highlight",
                        )
                    else:
    
                        # Merge the AMR knowledge layer into the draw graph.
                        # Every manifest locus stays a gene node even without
                        # RNA/protein quantitation; classes become mechanism
                        # (hyperedge) nodes linked to their genes. No omics
                        # values are invented for genome-only determinants.
                        amr_node_types = {k: set(v) for k, v in node_types.items()}
                        amr_edge_types = {
                            k: list(v) for k, v in edge_types.items()
                        }
                        aligned_loci = set()
                        if self.aligned is not None:
                            aligned_loci = set(
                                str(g) for g in
                                self.aligned['GeneID'].dropna().tolist()
                            )
                        elif (self.out_dir / 'aligned_multiomics.csv').exists():
                            import pandas as pd
                            _al = pd.read_csv(
                                self.out_dir / 'aligned_multiomics.csv',
                                low_memory=False,
                            )
                            aligned_loci = set(
                                str(g) for g in _al['GeneID'].dropna().tolist()
                            )
    
                        amr_class_map = {}
                        genome_only_loci = set()
                        for locus in locus_to_amr:
                            gnode = f"gene_{locus}"
                            amr_class_map[gnode] = locus_to_amr[locus]['amr_class']
                            amr_node_types['gene'].add(gnode)
                            if locus not in aligned_loci:
                                genome_only_loci.add(gnode)
    
                        mechanism_nodes = []
                        amr_mech_key = ('gene', 'amr_mechanism', 'amr_mechanism')
                        amr_edge_types.setdefault(amr_mech_key, [])
                        amr_node_types.setdefault('amr_mechanism', set())
                        for locus, rec in locus_to_amr.items():
                            mech = f"AMR:{rec['amr_class']}"
                            mechanism_nodes.append(mech)
                            amr_node_types['amr_mechanism'].add(mech)
                            amr_edge_types[amr_mech_key].append(
                                (f"gene_{locus}", mech)
                            )
    
                        amr_node_type_lists = {
                            k: sorted(v) for k, v in amr_node_types.items()
                        }
                        viz.plot_amr_highlight(
                            amr_node_type_lists, amr_edge_types,
                            amr_loci=list(locus_to_amr),
                            amr_class_map=amr_class_map,
                            title=f"{self.strain} — AMR Determinants in the Graph",
                            filename="amr_highlight.png",
                            max_nodes=0,
                            genome_only_loci=genome_only_loci,
                            mechanism_nodes=mechanism_nodes,
                        )
            except Exception as e:
                self._log("Visualization", f"Network plot failed: {e}")

        # 1. Regulation concordance
        if self.aligned is not None and not self.aligned.empty:
            viz.plot_regulation_concordance(self.aligned)

        # 1b. Metabolite regulation
        if getattr(self.metabolomics, 'table', None) is not None \
                and not self.metabolomics.table.empty:
            viz.plot_metabolite_regulation(self.metabolomics.table)

        # 2. COG distribution
        if self.cog_annotations:
            cog_counts = {}
            for ann in self.cog_annotations.values():
                for cog in ann.get('COG_classes', []):
                    cog_counts[cog] = cog_counts.get(cog, 0) + 1
            viz.plot_cog_distribution(cog_counts, cog_categories=COG_CATEGORIES)

        # 3. Hyperedge type distribution
        if self.hypergraph:
            he_types = {
                k: len(v) for k, v in self.hypergraph.hyperedge_type_map.items()
            }
            viz.plot_hyperedge_type_distribution(he_types)

        # 4. Incidence matrix
        if self.hypergraph and self.hypergraph.H is not None:
            viz.plot_incidence_matrix(
                self.hypergraph.H,
                hyperedge_types=self.hypergraph.hyperedge_type_map,
            )

        # 5. Summary dashboard
        genome_s = self.genome.summary() if self.genome else {}
        rna_s = self.rna.summary() if hasattr(self.rna, 'summary') else {}
        protein_s = getattr(self.protein, '_summary_cache', {'proteins': len(getattr(self.protein, 'protein_table', []))})

        if self.hetero_data is not None:
            graph_s = HeterogeneousGraphBuilder().summary()  # rough
            graph_s = {
                'gene_nodes': len(self.aligned) if self.aligned is not None else 0,
                'protein_nodes': len(self.aligned['ProteinID'].unique()) if self.aligned is not None else 0,
            }
        else:
            graph_s = {}

        hg_s = self.hypergraph.summary() if self.hypergraph else {}

        viz.plot_summary_dashboard(
            genome_s, rna_s, protein_s, graph_s, hg_s,
            title=self.strain,
        )

        self._end_stage("Visualization")

    # ── Stage 9: Export ──────────────────────────

    def run_export(self):
        """Export all data for downstream GNN training."""
        self._start_stage("Export")

        out = self.out_dir

        # 1. Multi-omics aligned table (for graph construction)
        if self.aligned is not None:
            self.aligned.to_csv(out / 'aligned_multiomics.csv', index=False)

        # 2. Gene features (with AMR knowledge annotation layer)
        if self.genome is not None:
            gene_features = self.genome.gene_table.copy()
            amr = getattr(self, 'amr_manifest', None)
            if amr is not None:
                locus_to_amr = amr.locus_to_amr(self.strain)
                gene_ids = gene_features.get(
                    'GeneID', gene_features.index
                ).tolist()
                amr_cls = [
                    [locus_to_amr[g]['amr_class']]
                    if g in locus_to_amr else []
                    for g in gene_ids
                ]
            else:
                amr_cls = []
            if amr_cls and any(amr_cls):
                gene_features['is_amr_gene'] = [bool(c) for c in amr_cls]
                gene_features['amr_class'] = [
                    ';'.join(sorted(c)) for c in amr_cls
                ]
            gene_features.to_csv(out / 'gene_features.csv', index=False)

        # 3. Metabolomics exports
        if getattr(self.metabolomics, 'table', None) is not None:
            self.metabolomics.table.to_csv(
                out / 'metabolomics_abundance.csv', index=False
            )
        if self.metabolite_to_pathway:
            rows = []
            for mid, paths in self.metabolite_to_pathway.items():
                for p in paths:
                    rows.append({'MetaboliteID': mid, 'KEGGPathway': p})
            if rows:
                pd.DataFrame(rows).to_csv(
                    out / 'metabolite_pathway_membership.csv', index=False
                )
        if self.metabolite_to_enzymes:
            rows = []
            for mid, ecs in self.metabolite_to_enzymes.items():
                for ec in ecs:
                    rows.append({'MetaboliteID': mid, 'EC_number': ec})
            if rows:
                pd.DataFrame(rows).to_csv(
                    out / 'metabolite_enzymes.csv', index=False
                )

        # 4. Edge list format for all graph relations
        if self.hetero_data is not None:
            export_keys = self.hetero_data.edge_types if hasattr(self.hetero_data, 'edge_types') else \
                          [k for k in self.hetero_data if isinstance(k, tuple)] if isinstance(self.hetero_data, dict) else []
            for key in export_keys:
                store = self.hetero_data[key]
                ei = store.get('edge_index') if isinstance(store, dict) else getattr(store, 'edge_index', None)
                if ei is not None and ei.size(1) > 0:
                    src_type, rel, dst_type = key
                    src_col = f'{src_type}_id'
                    dst_col = f'{dst_type}_id'
                    if src_col == dst_col:
                        dst_col = f'{dst_type}_target_id'
                    edge_df = pd.DataFrame({
                        src_col: ei[0].numpy(),
                        dst_col: ei[1].numpy(),
                        'relation': rel,
                    })
                    safe_name = f"edges_{src_type}_{rel}_{dst_type}.csv"
                    edge_df.to_csv(out / safe_name, index=False)

        # 5. Hypergraph incidence matrix
        if self.hypergraph and self.hypergraph.H is not None:
            np.save(out / 'hypergraph_incidence.npy', self.hypergraph.H)

            # Also save as sparse edge list
            rows, cols = np.where(self.hypergraph.H > 0)
            inc_df = pd.DataFrame({
                'node_idx': rows,
                'hyperedge_idx': cols,
            })
            inc_df.to_csv(out / 'hypergraph_incidence_edges.csv', index=False)

        # 6. PyTorch objects
        if self.hetero_data is not None:
            import torch
            torch.save(self.hetero_data, out / 'heterodata.pt')

        if self.hypergraph and self.hypergraph.H is not None:
            import torch
            laplacian = self.hypergraph.build_hypergraph_laplacian()
            for attempt in range(5):
                try:
                    torch.save(laplacian, out / 'hypergraph_laplacian.pt')
                    break
                except RuntimeError:
                    if attempt == 4:
                        self._log("Export", "Failed to save hypergraph_laplacian.pt, skipping")
                    else:
                        time.sleep(3)

        self._log("Export", f"All data saved to {out}")
        self._end_stage("Export")

    # ── Run all stages ───────────────────────────

    def _file_sha256(self, path) -> str:
        try:
            p = Path(path)
            if not p.exists() or not p.is_file():
                return 'missing'
            h = hashlib.sha256()
            with open(p, 'rb') as f:
                for chunk in iter(lambda: f.read(65536), b''):
                    h.update(chunk)
            return h.hexdigest()
        except Exception as e:
            return f'error: {e}'

    def _git_info(self):
        try:
            import subprocess
            commit = subprocess.run(
                ['git', 'rev-parse', 'HEAD'], capture_output=True, text=True,
                cwd=str(BASE_DIR), timeout=5,
            ).stdout.strip()[:12] or 'unknown'
            dirty = subprocess.run(
                ['git', 'status', '--porcelain'], capture_output=True,
                text=True, cwd=str(BASE_DIR), timeout=5,
            ).stdout.strip() != ''
            return {'commit': commit, 'dirty': dirty}
        except Exception:
            return {'commit': 'unknown', 'dirty': False}

    def write_manifest(self, run_params: Optional[dict] = None):
        """Write ``outputs/<strain>/manifest.json`` for reproducibility.

        Captures the pipeline identity (git commit + dirty flag), CLI params,
        a SHA-256 hash per input data file, an RNG seed, and installed package
        versions, so a future reader can reconstruct exactly what produced the
        outputs in this directory.
        """
        import importlib.metadata as im

        cfg = self.config
        data_files = {
            'genome_gff': str(cfg.get('genome_gff', '')),
            'protein_faa': str(cfg.get('protein_faa', '')),
            'cds_fna': str(cfg.get('cds_fna', '')),
            'rna_counts_file': str(cfg.get('rna_counts_file', '')),
            'prot_file': str(cfg.get('prot_file', '')),
        }
        input_hashes = {
            k: self._file_sha256(v) for k, v in data_files.items() if v
        }

        scalar_params = {}
        for key, val in cfg.items():
            if isinstance(val, (str, int, float, bool)) or val is None:
                scalar_params[key] = val

        packages = {}
        for name in ('numpy', 'pandas', 'matplotlib', 'networkx', 'requests',
                     'torch', 'torch_geometric', 'scipy'):
            try:
                packages[name] = im.version(name)
            except im.PackageNotFoundError:
                packages[name] = None

        manifest = {
            'pipeline': 'multiomics_graph',
            'strain': self.strain,
            'seed': self.seed,
            'git': self._git_info(),
            'timestamp': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'python_version': platform.python_version(),
            'platform': platform.platform(),
            'packages': packages,
            'parameters': {
                'correlation_threshold': self.correlation_threshold,
                'correlation_top_k': self.correlation_top_k,
                'skip_annotation': self.skip_annotation,
                'no_cache': self.no_cache,
            },
            'run_params': run_params or {},
            'config_scalars': scalar_params,
            'input_file_hashes_sha256': input_hashes,
        }
        out_path = self.out_dir / 'manifest.json'
        with open(out_path, 'w') as f:
            json.dump(manifest, f, indent=2)
        print(f"  [{self.strain}] [Manifest] reproducibility manifest -> {out_path.name}")

    def run_all(self):
        """Execute the complete multi-omics pipeline."""
        total_start = time.time()

        print(f"\n{'#' * 60}")
        print(f"# PROCESSING STRAIN: {self.strain} ({self.config.get('name', '')})")
        print(f"{'#' * 60}")

        self.run_genomics()
        self.run_transcriptomics()
        self.run_proteomics()
        self.run_metabolomics()
        self.run_integration()
        self.run_annotation()
        self.run_graph_construction()
        self.run_hypergraph()
        self.run_visualization()
        self.run_export()

        # Reproducibility: write a per-strain manifest capturing everything
        # needed to reconstruct this run (seed, params, file hashes, versions).
        self.write_manifest(run_params={
            'strain': self.strain,
            'correlation_threshold': self.correlation_threshold,
            'correlation_top_k': self.correlation_top_k,
        })

        total_time = time.time() - total_start
        self.timing['total'] = total_time

        # Print summary
        print(f"\n{'=' * 50}")
        print(f"  {self.strain} PIPELINE COMPLETE")
        print(f"{'=' * 50}")
        print(f"  Total time: {total_time:.1f}s")
        for stage, t in self.timing.items():
            if stage != 'total':
                print(f"    {stage}: {t:.1f}s")

        # Final data summary
        if self.aligned is not None:
            print(f"\n  Aligned triples: {len(self.aligned)}")
        if self.hypergraph:
            print(f"  Hypergraph: {self.hypergraph.summary().get('n_hyperedges', 0)} hyperedges")
        print(f"  Output: {self.out_dir}")

        return self


# ──────────────────────────────────────────────
# CLI Entry Point
# ──────────────────────────────────────────────


def discover_strains() -> list:
    """Scan for potential new strains in both layouts.

    Layout 1: strains/<name>/  (per-strain folder — recommended)
    Layout 2: genome/<dir>/   (legacy layout)

    Returns list of dicts with strain info.
    """
    found = []

    # Scan per-strain folders (strains/<name>/)
    if STRAINS_DIR.exists():
        for entry in sorted(STRAINS_DIR.iterdir()):
            if not entry.is_dir() or entry.name == 'TEMPLATE':
                continue
            config_path = entry / 'config.json'
            if config_path.exists():
                with open(config_path) as f:
                    cfg = json.load(f)
                # Check data accessibility via resolved config
                try:
                    resolved_cfg = _load_strain_from_folder(entry.name)
                    gff_ok = resolved_cfg.get('genome_gff') and Path(resolved_cfg['genome_gff']).exists()
                    rna_ok = resolved_cfg.get('rna_counts_file') and Path(resolved_cfg['rna_counts_file']).exists()
                    prot_ok = resolved_cfg.get('prot_file') and Path(resolved_cfg['prot_file']).exists()
                    has_data = gff_ok or rna_ok or prot_ok
                except Exception:
                    has_data = False
                found.append({
                    'layout': 'per-strain folder',
                    'name': entry.name,
                    'display_name': cfg.get('name', entry.name),
                    'config': str(config_path),
                    'configured': True,
                    'has_data': has_data,
                })
            else:
                found.append({
                    'layout': 'per-strain folder',
                    'name': entry.name,
                    'display_name': entry.name,
                    'config': 'missing',
                    'configured': False,
                    'has_data': (entry / 'genome').exists() or (entry / 'transcriptomic').exists() or (entry / 'proteomics').exists(),
                })

    # Scan legacy genome/ directory
    genome_dir = BASE_DIR / 'genome'
    if genome_dir.exists():
        for entry in sorted(genome_dir.iterdir()):
            if not entry.is_dir():
                continue
            gff_files = list(entry.rglob('*.gff'))
            if not gff_files:
                continue
            already = any(f['name'] == entry.name for f in found)
            if not already:
                found.append({
                    'layout': 'legacy (genome/ dir)',
                    'name': entry.name,
                    'display_name': entry.name,
                    'genome_gff': str(gff_files[0].relative_to(BASE_DIR)),
                    'configured': entry.name in list_available_strains(),
                })

    return found


def main():
    parser = argparse.ArgumentParser(
        description="Multi-Omics Heterogeneous Graph Architecture for Bacterial Systems Biology"
    )
    parser.add_argument('--strain', type=str, default='B36',
                        help='Strain identifier (see --list-strains)')
    parser.add_argument('--all-strains', action='store_true',
                        help='Process all configured strains')
    parser.add_argument('--output', type=str, default='outputs',
                        help='Base output directory')
    parser.add_argument('--correlation', type=float, default=0.7,
                        help='Correlation threshold for co-expression edges')
    parser.add_argument('--correlation-top-k', type=int, default=10,
                        help='Mutual top-k cap per node for co-expression edges '
                             '(0 = no cap; plain threshold rule is very dense)')
    parser.add_argument('--skip-annotation', action='store_true', default=False,
                        help='Skip online annotation API calls (UniProt, STRING, KEGG, eggNOG)')
    parser.add_argument('--no-skip-annotation', action='store_false', dest='skip_annotation',
                        help='Enable online annotation API calls (deprecated — online is now default)')
    parser.add_argument('--no-cache', action='store_true', default=False,
                        help='Bypass all cached results and re-download everything')
    parser.add_argument('--list-strains', action='store_true',
                        help='List all configured strains and exit')
    parser.add_argument('--discover', action='store_true',
                        help='Scan directories for potential new strains and print config templates')

    args = parser.parse_args()

    # ── List strains mode ──
    if args.list_strains:
        strains = list_available_strains()
        print(f"\nConfigured strains ({len(strains)}):")
        for s in strains:
            location = f"strains/{s}/config.json" if (STRAINS_DIR / s / 'config.json').exists() else str(STRAINS_CONFIG_PATH)
            # Try to get display name from the respective config source
            display = s
            if (STRAINS_DIR / s / 'config.json').exists():
                with open(STRAINS_DIR / s / 'config.json') as f:
                    display = json.load(f).get('name', s)
            elif STRAINS_CONFIG_PATH.exists():
                with open(STRAINS_CONFIG_PATH) as f:
                    cfg = json.load(f)
                    display = cfg.get(s, {}).get('name', s)
            print(f"  {s:20s}  {display}")
            print(f"  {'':20s}  [{location}]")
        print(f"\n  Use: python main.py --strain <name>")
        return

    # ── Discover mode ──
    if args.discover:
        print(f"\nScanning for strains...\n")
        found = discover_strains()
        if not found:
            print("  No strain directories found.")
            return
        print(f"  Found {len(found)} potential strain(s):\n")
        for info in found:
            if info['layout'] == 'per-strain folder':
                status = "[configured]" if info['configured'] else "[no config.json]"
                data = "has data" if info.get('has_data') else "empty"
                print(f"  {info['name']:20s}  [strain folder]  {status}, {data}")
            else:
                status = "✓ in config" if info['configured'] else "not configured"
                print(f"  {info['name']:20s}  [legacy]  {status}")
                print(f"  {'':20s}   GFF: {info.get('genome_gff', '?')}")
        print()
        print(f"  To add a new strain:")
        print(f"    1. Copy strains/TEMPLATE/ to strains/<name>/")
        print(f"    2. Place genome/ transcriptomic/ proteomics/ data inside")
        print(f"    3. Edit strains/<name>/config.json with sample IDs")
        print(f"    4. Run: python main.py --strain <name>")
        return

    # ── Process strains ──
    if args.all_strains:
        strains = list_available_strains()
    else:
        strains = [args.strain]

    for strain in strains:
        config = load_config(strain)
        pipeline = StrainPipeline(
            strain=strain,
            config=config,
            output_base=args.output,
            skip_annotation=args.skip_annotation,
            no_cache=args.no_cache,
            correlation_threshold=args.correlation,
            correlation_top_k=args.correlation_top_k,
        )
        pipeline.run_all()

    if len(strains) > 1:
        print(f"\n{'#' * 60}")
        print(f"# ALL {len(strains)} STRAINS COMPLETE")
        print(f"{'#' * 60}")


if __name__ == '__main__':
    main()
