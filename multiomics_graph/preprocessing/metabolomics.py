"""
metabolomics.py
Parse MetaboLights MAF (metabolite abundance files) into per-strain
condition summary tables compatible with the multi-omics graph.

Each MAF file follows the MetaboLights MAF format:
  - descriptive columns (database_identifier, chemical_formula, smiles,
    inchi, metabolite_identification, mass_to_charge, retention_time, ...)
  - one abundance column per sample (header e.g. `102.100.100/50913`)

Samples are assigned to RPMI vs Sera groups from the associated MTBLS
study-design file (`s_MTBLS2015.txt`), keyed by the numeric suffix of the
sample identifier.

Output table columns (per metabolite):
  - MetaboliteID   (ChEBI id, or name fallback)
  - MetaboliteName
  - Database
  - ChemicalFormula, MassToCharge, RetentionTime
  - Metabolite_RPMI, Metabolite_Sera, Metabolite_logFC
"""

import re
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional


class MetabolomicsProcessor:
    """
    Process MetaboLights MAF files into per-strain metabolite summaries.

    Parameters
    ----------
    maf_files : list of str
        Paths to one or more MAF TSV files (GC-MS, LC-MS, ...)
    sample_file : str, optional
        Path to the study `s_*.txt` sample sheet used to assign
        RPMI / Sera condition labels.
    """

    def __init__(self, maf_files: List[str],
                 sample_file: Optional[str] = None):
        self.maf_files = [Path(f) for f in maf_files]
        self.sample_file = Path(sample_file) if sample_file else None
        self.sample_condition = {}   # {sample_key(str): 'RPMI'|'Pooled sera'}
        self.raw = None
        self.table = None
        self.raw_replicates = pd.DataFrame()   # MetaboliteID x per-sample

    # ── Sample condition mapping ─────────────────────────────
    def parse_sample_file(self) -> Dict[str, str]:
        """
        Parse the MTBLS sample sheet: map sample key -> growth media.

        The sheet has one line per sample with `Source Name` in column 0
        (e.g. `102.100.100/50913`) and the growth media in the
        `Factor Value[Growth media]` column.
        """
        if self.sample_file is None:
            return {}
        if not self.sample_file.exists():
            print(f"  [Metabolomics] sample file not found: {self.sample_file}")
            return {}

        mapping = {}
        rows = self.sample_file.read_text(encoding='utf-8').splitlines()
        if not rows:
            return mapping

        header = rows[0].split('\t')
        media_idx = None
        for i, col in enumerate(header):
            if 'Factor Value[Growth media]' in col:
                media_idx = i
                break
        if media_idx is None:
            return mapping

        for line in rows[1:]:
            cells = line.split('\t')
            if len(cells) < 2:
                continue
            source = cells[0].strip()
            media = cells[media_idx].strip() if media_idx < len(cells) else ''
            if not source or not media:
                continue
            key = source.split('/')[-1]
            mapping[key] = media

        self.sample_condition = mapping
        return mapping

    # ── MAF parsing ──────────────────────────────────────────
    @staticmethod
    def _load_maf(path: Path) -> pd.DataFrame:
        """Load a MAF TSV; returns DataFrame with raw columns."""
        df = pd.read_csv(path, sep='\t', dtype={'database_identifier': str},
                         low_memory=False)
        return df

    def load_raw(self) -> pd.DataFrame:
        """Load all MAF files and return a single long DataFrame.

        Each MAF file is tagged with its originating instrument (file stem)
        so later stages can normalise across instruments before averaging.
        """
        frames = []
        for maf in self.maf_files:
            if not maf.exists():
                print(f"  [Metabolomics] MAF file not found: {maf}")
                continue
            df = self._load_maf(maf)
            df['_instrument'] = maf.stem
            frames.append(df)
        if not frames:
            raise FileNotFoundError(f"No MAF files found in {self.maf_files}")
        self.raw = pd.concat(frames, ignore_index=True)
        return self.raw

    # ── Per-instrument scale normalisation ────────────────
    @staticmethod
    def _normalise_instruments(df: pd.DataFrame,
                               abundance_cols: List[str]) -> pd.DataFrame:
        """
        Bring per-instrument abundances onto a common scale before averaging.

        Different instruments (e.g. GC-MS relative log-abundance vs LC-MS raw
        ion counts) report values on very different scales. Averaging them raw
        lets the higher-scale instrument dominate the merged estimate. This
        scales each instrument's columns by the ratio of its median positive
        abundance to the mean of all instrument medians, which preserves
        within-instrument log fold-changes while making cross-instrument
        averaging meaningful.

        Returns a copy with `_instrument` retained for downstream QC.
        """
        if '_instrument' not in df.columns:
            return df

        out = df.copy()
        instruments = out['_instrument'].unique()
        if len(instruments) < 2:
            return out

        medians = {}
        for inst in instruments:
            sub = out[out['_instrument'] == inst]
            mat = sub[abundance_cols].to_numpy(dtype=float)
            pos = mat[mat > 0]
            medians[inst] = float(np.median(pos)) if len(pos) else np.nan

        valid = [m for m in medians.values() if m and m > 0]
        if not valid:
            return out
        target = float(np.mean(valid))

        scales = {inst: target / medians[inst] for inst in instruments
                  if medians.get(inst) and medians[inst] > 0}

        # Per-instrument scale factor, broadcast to every abundance column.
        out[abundance_cols] = out[abundance_cols].astype(float).mul(
            out['_instrument'].map(scales).to_numpy(), axis=0
        )
        return out

    # ── Canonicalise each row ────────────────────────────────
    @staticmethod
    def _metabolite_id(row) -> Optional[str]:
        """Stable metabolite id: ChEBI accession when present else name."""
        dbid = str(row.get('database_identifier', '')).strip()
        m = re.match(r'^CHEBI:(\d+)$', dbid)
        if m:
            return f"CHEBI:{m.group(1)}"
        name = str(row.get('metabolite_identification', '')).strip()
        if name and name.lower() != 'nan':
            return f"M-{name}"
        return None

    # ── Main entry point ─────────────────────────────────────
    def load(self, rpmi_samples: Optional[List] = None,
             sera_samples: Optional[List] = None,
             log2: bool = True,
             min_detection: float = 0.0,
             condition_map: Optional[Dict[str, str]] = None) -> pd.DataFrame:
        """
        Build a per-metabolite summary for one strain.

        Parameters
        ----------
        rpmi_samples : list, optional
            Sample keys to average as RPMI replicates. If None, inferred
            from the condition map (values starting with 'RPMI').
        sera_samples : list, optional
            Sample keys to average as Sera replicates. If None, inferred
            from the condition map (values starting with 'Sera').
        log2 : bool
            Whether to return log2-transformed mean abundances.
        min_detection : float
            Drop metabolites whose max mean abundance is below this value.
        condition_map : dict, optional
            Override {sample_key: condition} mapping.
        """
        if condition_map is not None:
            self.sample_condition = condition_map
        if not self.sample_condition and self.sample_file:
            self.parse_sample_file()

        if self.raw is None:
            self.load_raw()

        df = self.raw.copy()
        abundance_cols = [
            c for c in df.columns
            if re.match(r'^\d+$', c) or re.match(r'^102\.\d+\.\d+/\d+$', str(c))
        ]

        # Normalise across instruments (GC-MS vs LC-MS) before averaging so a
        # single high-scale instrument does not dominate merged abundances.
        df = self._normalise_instruments(df, abundance_cols)

        # Determine condition for each abundance column
        def _cond(col: str):
            key = str(col).split('/')[-1]
            return self.sample_condition.get(key)

        rpmi_cols = [c for c in abundance_cols
                     if (str(_cond(c)).lower().startswith('rpm'))]
        sera_cols = [c for c in abundance_cols
                     if (str(_cond(c)).lower().startswith('sera'))]

        # Respect explicit sample lists (config) when provided
        if rpmi_samples:
            rpmi_cols = [c for c in abundance_cols
                         if str(c).split('/')[-1] in {str(s) for s in rpmi_samples}]
        if sera_samples:
            sera_cols = [c for c in abundance_cols
                         if str(c).split('/')[-1] in {str(s) for s in sera_samples}]

        if not rpmi_cols and not sera_cols:
            raise ValueError(
                "No sample columns could be mapped to RPMI/Sera. "
                f"Abundance columns found: {abundance_cols[:5]}..."
            )

        def _mat(cols):
            return df[cols].apply(pd.to_numeric, errors='coerce')

        rpmi_mean = _mat(rpmi_cols).fillna(0).mean(axis=1) if rpmi_cols else 0.0
        sera_mean = _mat(sera_cols).fillna(0).mean(axis=1) if sera_cols else 0.0

        if log2:
            rpmi_l = np.log2(rpmi_mean + 1)
            sera_l = np.log2(sera_mean + 1)
        else:
            rpmi_l, sera_l = rpmi_mean, sera_mean

        df['Metabolite_RPMI'] = rpmi_l
        df['Metabolite_Sera'] = sera_l
        df['Metabolite_logFC'] = sera_l - rpmi_l
        df['MetaboliteID'] = df.apply(self._metabolite_id, axis=1)
        df['MetaboliteName'] = df['metabolite_identification'].astype(str).str.strip()
        df['Database'] = df['database'].astype(str).str.strip()
        df['ChemicalFormula'] = df['chemical_formula'].astype(str).str.strip()
        df['MassToCharge'] = pd.to_numeric(df['mass_to_charge'], errors='coerce')
        df['RetentionTime'] = pd.to_numeric(df['retention_time'], errors='coerce')
        df['Instruments'] = df['_instrument'].astype(str).str.strip() \
            if '_instrument' in df.columns else ''

        keep = ['MetaboliteID', 'MetaboliteName', 'Database',
                'ChemicalFormula', 'MassToCharge', 'RetentionTime',
                'Metabolite_RPMI', 'Metabolite_Sera', 'Metabolite_logFC',
                'Instruments']
        out = df[keep].copy()
        out = out[out['MetaboliteID'].notna() & (out['MetaboliteID'] != '')]

        # Merge identical metabolite ids across instruments: average the
        # abundance values (already log2-transformed per instrument) and
        # coalesce the descriptive metadata (first non-null value wins).
        def _first_nonnull(s: pd.Series):
            vals = s.dropna()
            return vals.iloc[0] if len(vals) else np.nan

        agg = out.groupby('MetaboliteID', as_index=False).agg({
            'Metabolite_RPMI': 'mean',
            'Metabolite_Sera': 'mean',
            'Metabolite_logFC': 'mean',
            'MetaboliteName': _first_nonnull,
            'Database': lambda s: ', '.join(str(v) for v in s.dropna().unique() if str(v) not in ('', 'nan')),
            'ChemicalFormula': _first_nonnull,
            'MassToCharge': _first_nonnull,
            'RetentionTime': _first_nonnull,
            'Instruments': lambda s: ', '.join(sorted({str(v) for v in s.dropna().unique() if str(v) not in ('', 'nan')})),
        })
        # Recompute logFC from the averaged condition means
        agg['Metabolite_logFC'] = agg['Metabolite_Sera'] - agg['Metabolite_RPMI']
        agg['n_instruments'] = agg['Instruments'].str.split(', ').apply(len)

        if min_detection > 0:
            agg = agg[
                (agg['Metabolite_RPMI'] >= min_detection) |
                (agg['Metabolite_Sera'] >= min_detection)
            ]

        out = agg.sort_values('Metabolite_RPMI', ascending=False).reset_index(drop=True)
        self.table = out

        # Retain the per-sample abundance matrix (MetaboliteID x sample
        # columns) for metabolite--metabolite abundance-correlation edges.
        # Values are instrument-normalised (so GC-MS and LC-MS are on a
        # comparable scale) and averaged across duplicate/instrument rows;
        # sample columns where a metabolite was not detected are left as NaN
        # and handled by the correlation builder's mean-imputation.
        try:
            rep = df[['MetaboliteID'] + abundance_cols].copy()
            rep = rep[rep['MetaboliteID'].notna() & (rep['MetaboliteID'] != '')]
            rep[abundance_cols] = rep[abundance_cols].apply(
                pd.to_numeric, errors='coerce')
            rep = rep.groupby('MetaboliteID', as_index=False)[
                abundance_cols].mean()
            self.raw_replicates = rep.set_index('MetaboliteID')
        except Exception:
            self.raw_replicates = pd.DataFrame()
        return out

    def summary(self) -> dict:
        """Return summary statistics."""
        t = getattr(self, 'table', None)
        if t is None:
            return {'metabolites': 0}
        n = len(t)
        ups = int((t['Metabolite_logFC'] > 0.5).sum())
        down = int((t['Metabolite_logFC'] < -0.5).sum())
        multi = int(t['n_instruments'].gt(1).sum()) \
            if 'n_instruments' in t.columns else 0
        return {
            'metabolites': n,
            'up_regulated': ups,
            'down_regulated': down,
            'stable': n - ups - down,
            'detected_on_multiple_instruments': multi,
        }
