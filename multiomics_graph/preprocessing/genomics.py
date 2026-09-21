"""
genomics.py
Parse bacterial genome GFF files and extract CDS-level features.

The GFF is the central backbone linking all omics layers:
  - CDS entries connect locus_tag (GeneID) ↔ protein_id (ProteinID)
  - Genomic coordinates define gene-gene neighborhood relationships
  - The gene symbol (Name/gene attribute) provides the human-readable name
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re


class GenomicsProcessor:
    """
    Process bacterial genome GFF files.

    Extracts CDS features from NCBI RefSeq GFF3 format.
    Builds the GeneID ↔ ProteinID mapping table used for multi-omics integration.

    Attributes
    ----------
    strain : str
        Strain identifier (e.g., 'B36')
    gene_table : pd.DataFrame
        Gene-level table with columns: GeneID, GeneName, ProteinID,
        start, end, strand, replicon, CDS_length
    replicon_count : int
        Number of distinct replicons (chromosomes + plasmids)
    """

    def __init__(self, strain: str = "unknown"):
        self.strain = strain
        self.gene_table = pd.DataFrame()
        self._protein_id_map = {}

    def parse_gff(self, gff_path: str) -> pd.DataFrame:
        """
        Parse a bacterial GFF file, extracting only CDS entries.

        NCBI GFF3 CDS entries have attributes:
          - ID=cds-{locus_tag}
          - Parent=gene-{locus_tag}
          - Dbxref=GenBank:{protein_id}
          - protein_id={WP_accession}
          - locus_tag={locus_tag}
          - gene={gene_symbol}  (when available)
          - product=description

        Parameters
        ----------
        gff_path : str
            Path to the .gff file

        Returns
        -------
        pd.DataFrame with CDS-level information
        """
        gff_path = Path(gff_path)
        if not gff_path.exists():
            raise FileNotFoundError(f"GFF not found: {gff_path}")

        cds_records = []

        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                if parts[2] != 'CDS':
                    continue

                attr_str = parts[8]
                attrs = dict(
                    item.split('=', 1) for item in attr_str.split(';') if '=' in item
                )

                locus_tag = attrs.get('locus_tag', '')
                protein_id = attrs.get('protein_id', '')
                gene_symbol = attrs.get('gene', '')
                product = attrs.get('product', '')
                parent = attrs.get('Parent', '')

                # Some GFF versions store protein_id in Dbxref
                if not protein_id:
                    dbxref = attrs.get('Dbxref', '')
                    match = re.search(r'GenBank:([^\s,;]+)', dbxref)
                    if match:
                        protein_id = match.group(1)

                if not locus_tag and parent:
                    locus_tag = parent.replace('gene-', '')

                record = {
                    'GeneID': locus_tag,
                    'GeneName': gene_symbol,
                    'ProteinID': protein_id,
                    'Product': product,
                    'replicon': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'CDS_length': int(parts[4]) - int(parts[3]) + 1,
                }
                cds_records.append(record)

        df = pd.DataFrame(cds_records)

        # Remove entries missing both GeneID and ProteinID
        df = df[df['GeneID'].notna() & (df['GeneID'] != '')]
        df = df.sort_values(['replicon', 'start']).reset_index(drop=True)

        self.gene_table = df
        self._build_protein_id_map()
        return df

    def _build_protein_id_map(self):
        """Build ProteinID → GeneID lookup dictionary."""
        self._protein_id_map = {}
        for _, row in self.gene_table.iterrows():
            pid = row['ProteinID']
            if pid and not pd.isna(pid):
                self._protein_id_map[pid] = row['GeneID']

    def get_gene_by_protein(self, protein_id: str) -> str:
        """Look up GeneID from ProteinID."""
        return self._protein_id_map.get(protein_id, '')

    def get_protein_by_gene(self, gene_id: str) -> str:
        """Look up ProteinID from GeneID."""
        match = self.gene_table.loc[
            self.gene_table['GeneID'] == gene_id, 'ProteinID'
        ]
        return match.iloc[0] if not match.empty else ''

    def compute_genomic_edges(self, max_gap: int = 5000) -> pd.DataFrame:
        """
        Build gene-gene edges based on genomic proximity.

        Consecutive genes on the same replicon within `max_gap` bp
        receive an edge. The weight captures physical proximity.

        Parameters
        ----------
        max_gap : int
            Maximum intergenic distance (bp) to consider adjacent.

        Returns
        -------
        pd.DataFrame with columns: GeneID_A, GeneID_B, replicon, distance, weight
        """
        if self.gene_table.empty:
            return pd.DataFrame()

        edges = []
        for replicon, group in self.gene_table.groupby('replicon'):
            group = group.sort_values('start')
            genes = group.to_dict('records')
            for i in range(len(genes) - 1):
                g1, g2 = genes[i], genes[i + 1]
                dist = g2['start'] - g1['end']
                if 0 < dist <= max_gap:
                    weight = 1.0 / max(1.0, float(dist))
                    edges.append({
                        'GeneID_A': g1['GeneID'],
                        'GeneID_B': g2['GeneID'],
                        'replicon': replicon,
                        'distance': int(dist),
                        'weight': weight,
                    })

        return pd.DataFrame(edges)

    def summary(self) -> dict:
        """Return summary statistics of the genome."""
        if self.gene_table.empty:
            return {'strain': self.strain, 'genes': 0}
        return {
            'strain': self.strain,
            'genes': len(self.gene_table),
            'replicons': self.gene_table['replicon'].nunique(),
            'with_protein_id': self.gene_table['ProteinID'].astype(bool).sum(),
            'with_gene_symbol': self.gene_table['GeneName'].astype(bool).sum(),
            'total_cds_length': int(self.gene_table['CDS_length'].sum()),
        }
