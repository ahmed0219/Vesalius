"""
proteomics.py
Process bacterial SWATH-MS proteomics data (protein area quantification).

Handles:
  - Loading protein abundance from Excel (Area - proteins sheet)
  - Filtering bacterial proteins (remove reverse hits, contaminants)
  - Computing per-condition means across biological replicates
  - Differential abundance analysis (logFC) and regulation calls

Biological context:
  The proteomics data comes from SWATH-MS acquisition.
  Quantification is at the protein area level.
  Each condition has 6 biological replicates.
"""

from unittest import result

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional


# Sample ID mapping for B36 proteomics
B36_RPMI_SAMPLE_IDS = [50973, 50974, 50975, 50976, 50977, 50978]
B36_SERA_SAMPLE_IDS = [50979, 50980, 50981, 50982, 50983, 50984]

# Column name pattern for the Area - proteins sheet
# Format: "{sample_id} (filename.wiff (sample 1))"
B36_RPMI_COLS = [
    f"{sid} (180525_P21506_SWATH_{sid}.wiff (sample 1))"
    for sid in B36_RPMI_SAMPLE_IDS
]
B36_SERA_COLS = [
    f"{sid} (180525_P21506_SWATH_{sid}.wiff (sample 1))"
    for sid in B36_SERA_SAMPLE_IDS
]


class ProteomicsProcessor:
    """
    Process bacterial SWATH-MS proteomics data.

    Attributes
    ----------
    strain : str
    protein_table : pd.DataFrame
        Summarised protein abundance with Protein_RPMI, Protein_Sera, Protein_logFC
    raw_data : pd.DataFrame
        Full protein quantification matrix
    regulation_threshold : float
        logFC threshold for regulation calls
    """

    def __init__(self, strain: str = "unknown",
                 regulation_threshold: float = 1.0,
                 std_multiplier: float = 0.0):
        self.strain = strain
        self.protein_table = pd.DataFrame()
        self.raw_data = pd.DataFrame()
        self.regulation_threshold = regulation_threshold
        self.std_multiplier = std_multiplier
        self._bacterial_mask = None

    def load_protein_groups(self, filepath: str,
                             sheet_name: str = "Area - proteins",
                             protein_col: str = "Protein Accession",
                             gene_col: str = "Gene Symbol") -> pd.DataFrame:
        """
        Load protein quantification data from SWATH-MS output.

        Parameters
        ----------
        filepath : str
            Path to the protein groups Excel file
        sheet_name : str
            Worksheet containing protein-level area quantification
        protein_col : str
            Column containing protein identifiers (e.g., WP_*)
        gene_col : str
            Column containing gene symbols

        Returns
        -------
        pd.DataFrame with raw protein quantification
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"Proteomics file not found: {filepath}")

        df = pd.read_excel(filepath, sheet_name=sheet_name)
        self._protein_col = protein_col if protein_col in df.columns else df.columns[0]
        self._gene_col = gene_col if gene_col in df.columns else None

        self.raw_data = df
        return df

    def filter_bacterial_proteins(self, valid_prefix: str = "WP_") -> pd.DataFrame:
        """
        Retain only bacterial proteins with RefSeq WP_ accessions.

        Removes:
          - Reverse/decoy hits (REV__ prefix)
          - Contaminants (CON__ prefix)
          - Non-bacterial proteins (no WP_ accession)

        Parameters
        ----------
        valid_prefix : str
            Accession prefix for bacterial RefSeq proteins

        Returns
        -------
        pd.DataFrame with only bacterial protein entries
        """
        if self.raw_data.empty:
            raise ValueError("No data loaded. Call load_protein_groups() first.")

        df = self.raw_data.copy()

        # Parse protein accessions (may contain semicolons for isoforms)
        def _is_bacterial(acc_str):
            if pd.isna(acc_str):
                return False
            accessions = str(acc_str).split(';')
            for acc in accessions:
                acc = acc.strip()
                if acc.startswith('REV__') or acc.startswith('CON__'):
                    return False
                if acc.startswith(valid_prefix):
                    return True
            return False

        mask = df[self._protein_col].apply(_is_bacterial)
        self._bacterial_mask = mask
        df = df[mask].copy()
        self.raw_data = df
        return df

    def compute_condition_means(
        self,
        rpmi_cols: Optional[List[str]] = None,
        sera_cols: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Compute mean abundance per condition across replicates.

        Derives:
          Protein_RPMI   = mean log2(RPMI replicate areas + 1)
          Protein_Sera   = mean log2(Serum replicate areas + 1)
          Protein_logFC  = Protein_Sera - Protein_RPMI
          Regulation     = Up / Down / Stable

        Parameters
        ----------
        rpmi_cols : list, optional
            Column names for RPMI condition
        sera_cols : list, optional
            Column names for Serum condition

        Returns
        -------
        pd.DataFrame with summarised protein data
        """
        if self.raw_data.empty:
            raise ValueError("No data loaded. Call load_protein_groups() first.")

        rpmi = rpmi_cols or B36_RPMI_COLS
        sera = sera_cols or B36_SERA_COLS

        # Validate columns
        available = set(self.raw_data.columns)
        rpmi = [c for c in rpmi if c in available]
        sera = [c for c in sera if c in available]

        if not rpmi:
            # Auto-detect RPMI columns by sample IDs in column names
            rpmi = [
                c for c in self.raw_data.columns
                if any(str(sid) in c for sid in B36_RPMI_SAMPLE_IDS)
            ]
        if not sera:
            sera = [
                c for c in self.raw_data.columns
                if any(str(sid) in c for sid in B36_SERA_SAMPLE_IDS)
            ]

        if not rpmi or not sera:
            raise ValueError(
                f"Could not identify condition columns. "
                f"Available: {list(available)[:10]}..."
            )

        # log2 transform: log2(area + 1)
        for col in rpmi + sera:
            self.raw_data[col] = pd.to_numeric(
                self.raw_data[col], errors='coerce'
            ).fillna(0)

        result = pd.DataFrame()

        # Extract primary Protein Accession (first before semicolon)
        result['ProteinID'] = self.raw_data[self._protein_col].apply(
            lambda x: str(x).split(';')[0].strip()
        )
        if self._gene_col and self._gene_col in self.raw_data.columns:
            result['GeneName'] = self.raw_data[self._gene_col]
        else:
            result['GeneName'] = result['ProteinID']

        # Mean log2 area per condition
        result['Protein_RPMI'] = np.log2(
            self.raw_data[rpmi] + 1
        ).mean(axis=1)

        result['Protein_Sera'] = np.log2(
            self.raw_data[sera] + 1
        ).mean(axis=1)
        result['Protein_logFC'] = result['Protein_Sera'] - result['Protein_RPMI']

        # Regulation call — use std-based threshold if configured
        if self.std_multiplier > 0:
            thresh = self.std_multiplier * result['Protein_logFC'].std()
        else:
            thresh = self.regulation_threshold
        result['Protein_Regulation'] = np.select(
            [
                result['Protein_logFC'] > thresh,
                result['Protein_logFC'] < -thresh,
            ],
            ['Up', 'Down'],
            default='Stable',
        )

        self.protein_table = result
        return result

    def summary(self) -> dict:
        """Return summary statistics."""
        if self.protein_table.empty:
            return {'strain': self.strain, 'proteins': 0}
        reg = self.protein_table['Protein_Regulation'].value_counts()
        return {
            'strain': self.strain,
            'proteins': len(self.protein_table),
            'up_regulated': int(reg.get('Up', 0)),
            'down_regulated': int(reg.get('Down', 0)),
            'stable': int(reg.get('Stable', 0)),
            'mean_logFC': float(self.protein_table['Protein_logFC'].mean()),
        }
