"""
transcriptomics.py
Process bacterial RNA-seq expression data (log2 CPM format).

Handles:
  - Loading gene-level log2 CPM expression tables
  - Mapping sample IDs to conditions (RPMI vs Serum)
  - Computing per-condition means across biological replicates
  - Differential expression (logFC) and regulation calls

Biological context:
  RPMI  = normal in vitro growth (baseline)
  Serum = bloodstream-like stress environment
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional


# Sample-to-condition mapping for B36 (from GEO GSE152966)
B36_SAMPLE_MAP = {
    # RPMI replicates
    '50857': 'RPMI', '50858': 'RPMI',
    '51033': 'RPMI', '51034': 'RPMI',
    '51035': 'RPMI', '51036': 'RPMI',
    # Pooled sera replicates
    '50863': 'Sera', '50864': 'Sera',
    '51037': 'Sera', '51038': 'Sera',
    '51039': 'Sera', '51040': 'Sera',
}

B36_RPMI_SAMPLES = ['50857', '50858', '51033', '51034', '51035', '51036']
B36_SERA_SAMPLES = ['50863', '50864', '51037', '51038', '51039', '51040']


class TranscriptomicsProcessor:
    """
    Process bacterial transcriptomics data.

    The input is a gene × sample matrix of log2 CPM values.
    The output summarises expression per condition and
    identifies differentially regulated genes.

    Attributes
    ----------
    strain : str
    expression_table : pd.DataFrame
        Long-format summary with RNA_RPMI, RNA_Sera, RNA_logFC
    raw_data : pd.DataFrame
        Full expression matrix
    regulation_threshold : float
        logFC threshold for calling Up/Down regulation
    """

    def __init__(self, strain: str = "unknown",
                 regulation_threshold: float = 1.0,
                 std_multiplier: float = 0.0):
        self.strain = strain
        self.expression_table = pd.DataFrame()
        self.raw_data = pd.DataFrame()
        self.regulation_threshold = regulation_threshold
        self.std_multiplier = std_multiplier

    def load_log2_cpm(self, filepath: str,
                       sample_map: Optional[Dict[str, str]] = None,
                       sheet_name: Optional[str] = None) -> pd.DataFrame:
        """
        Load gene expression data in log2 CPM format.

        The file should have genes as rows and sample IDs as columns.
        First two columns are expected to be GeneID and GeneName.

        Parameters
        ----------
        filepath : str
            Path to expression file (CSV, TSV, or gzipped TSV)
        sample_map : dict, optional
            Maps sample column names to conditions ('RPMI' or 'Sera')
        sheet_name : str, optional
            For Excel files, the sheet name

        Returns
        -------
        pd.DataFrame with raw expression data
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"Expression file not found: {filepath}")

        # Detect file type
        suffix = str(filepath.suffix).lower()
        if filepath.name.endswith('.gz'):
            suffix = ''.join(Path(filepath.stem).suffixes).lower()

        if suffix in ('.csv',):
            df = pd.read_csv(filepath)
        elif suffix in ('.tsv', '.txt'):
            df = pd.read_csv(filepath, sep='\t')
        elif suffix in ('.xlsx', '.xls'):
            df = pd.read_excel(filepath, sheet_name=sheet_name or 0)
        else:
            # Try TSV as default
            df = pd.read_csv(filepath, sep='\t')

        # Identify gene identifier columns
        id_cols = []
        for col in df.columns:
            col_lower = col.lower().replace(' ', '_')
            if col_lower in ('geneid', 'gene_id', 'gene', 'locus_tag'):
                id_cols.append(col)
                break

        name_cols = []
        for col in df.columns:
            col_lower = col.lower().replace(' ', '_')
            if col_lower in ('genename', 'gene_name', 'name', 'symbol'):
                name_cols.append(col)
                break

        self._gene_id_col = id_cols[0] if id_cols else df.columns[0]
        self._gene_name_col = name_cols[0] if name_cols else None

        self.raw_data = df.copy()
        self.sample_map = sample_map or B36_SAMPLE_MAP
        return df

    def compute_condition_means(
        self,
        rpmi_samples: Optional[List[str]] = None,
        sera_samples: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Compute mean expression per condition across biological replicates.

        Derives:
          RNA_RPMI   = mean of RPMI replicate log2 CPM values
          RNA_Sera   = mean of Serum replicate log2 CPM values
          RNA_logFC  = RNA_Sera - RNA_RPMI
          Regulation = Up / Down / Stable (based on logFC threshold)

        Parameters
        ----------
        rpmi_samples : list, optional
            Column names for RPMI condition replicates
        sera_samples : list, optional
            Column names for Serum condition replicates

        Returns
        -------
        pd.DataFrame with summarised expression data
        """
        if self.raw_data.empty:
            raise ValueError("No expression data loaded. Call load_log2_cpm() first.")

        rpmi = rpmi_samples or B36_RPMI_SAMPLES
        sera = sera_samples or B36_SERA_SAMPLES

        # Validate columns exist
        available = set(self.raw_data.columns)
        missing_rpmi = [c for c in rpmi if c not in available]
        missing_sera = [c for c in sera if c not in available]
        if missing_rpmi:
            raise ValueError(f"RPMI columns not in data: {missing_rpmi}")
        if missing_sera:
            raise ValueError(f"Sera columns not in data: {missing_sera}")

        # Ensure numeric
        for col in rpmi + sera:
            self.raw_data[col] = pd.to_numeric(
                self.raw_data[col], errors='coerce'
            ).fillna(0)

        result = pd.DataFrame()
        result['GeneID'] = self.raw_data[self._gene_id_col]
        if self._gene_name_col:
            result['GeneName'] = self.raw_data[self._gene_name_col]
        else:
            result['GeneName'] = result['GeneID']

        result['RNA_RPMI'] = self.raw_data[rpmi].mean(axis=1)
        result['RNA_Sera'] = self.raw_data[sera].mean(axis=1)
        result['RNA_logFC'] = result['RNA_Sera'] - result['RNA_RPMI']

        # Regulation call — use std-based threshold if configured
        if self.std_multiplier > 0:
            thresh = self.std_multiplier * result['RNA_logFC'].std()
        else:
            thresh = self.regulation_threshold
        result['RNA_Regulation'] = np.select(
            [
                result['RNA_logFC'] > thresh,
                result['RNA_logFC'] < -thresh,
            ],
            ['Up', 'Down'],
            default='Stable',
        )

        self.expression_table = result
        return result

    def get_regulated_genes(self, direction: str = 'Up') -> pd.DataFrame:
        """Return genes regulated in the specified direction."""
        return self.expression_table[
            self.expression_table['RNA_Regulation'] == direction
        ].copy()

    def summary(self) -> dict:
        """Return summary statistics of transcriptome data."""
        if self.expression_table.empty:
            return {'strain': self.strain, 'genes': 0}
        reg = self.expression_table['RNA_Regulation'].value_counts()
        return {
            'strain': self.strain,
            'genes': len(self.expression_table),
            'up_regulated': int(reg.get('Up', 0)),
            'down_regulated': int(reg.get('Down', 0)),
            'stable': int(reg.get('Stable', 0)),
            'mean_logFC': float(self.expression_table['RNA_logFC'].mean()),
            'max_logFC': float(self.expression_table['RNA_logFC'].max()),
            'min_logFC': float(self.expression_table['RNA_logFC'].min()),
        }
