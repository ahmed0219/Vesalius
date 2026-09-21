"""
eggnog.py
Retrieve eggNOG (evolutionary genealogy of genes: Non-supervised
Orthologous Groups) annotations for bacterial genes.

eggNOG provides:
  - COG/NOG functional categories (e.g., J, K, L, ...)
  - Orthologous groups
  - Functional descriptions
  - Taxonomic scope

COG categories (used in hypergraph construction):
  J: Translation, ribosomal structure
  K: Transcription
  L: Replication, recombination, repair
  D: Cell cycle control, division
  V: Defense mechanisms
  T: Signal transduction
  M: Cell wall/membrane biogenesis
  U: Intracellular trafficking
  O: Post-translational modification
  C: Energy production
  G: Carbohydrate transport
  E: Amino acid transport
  F: Nucleotide transport
  H: Coenzyme transport
  I: Lipid transport
  P: Inorganic ion transport
  Q: Secondary metabolites
  R: General function prediction
  S: Function unknown
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from typing import Dict, List, Optional

from .identifiers import normalize_wp_id


CACHE_DIR = Path(__file__).parent.parent / "outputs" / "cache"

COG_CATEGORIES = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'O': 'Post-translational modification, protein turnover, chaperones',
    'X': 'Mobilome: prophages, transposons',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown',
}


class EggNOGAnnotator:
    """
    Assign eggNOG/COG functional categories to bacterial genes.

    Uses a precomputed mapping file or the eggNOG API.

    Attributes
    ----------
    annotations : dict
        {GeneID: {cog_category, cog_class, description, ...}}
    category_counts : dict
        {COG_category: count_of_genes}
    """

    def __init__(self, cache_dir: Optional[Path] = None):
        self.cache_dir = Path(cache_dir) if cache_dir else CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.annotations = {}
        self.category_counts = {k: 0 for k in COG_CATEGORIES}

    def load_precomputed(self, filepath: str,
                          gene_col: str = 'GeneID',
                          cog_col: str = 'COG_category') -> pd.DataFrame:
        """
        Load eggNOG annotations from a precomputed TSV/CSV file.

        Expected format (from eggNOG-mapper output):
          query, seed_ortholog, evalue, score, eggNOG_OGs,
          COG_category, Description, ...

        Parameters
        ----------
        filepath : str
            Path to eggNOG annotations file
        gene_col : str
            Column name for gene identifier
        cog_col : str
            Column name for COG category

        Returns
        -------
        pd.DataFrame with annotations
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"eggNOG file not found: {filepath}")

        if filepath.suffix in ('.csv',):
            df = pd.read_csv(filepath)
        else:
            df = pd.read_csv(filepath, sep='\t')

        required = {gene_col, cog_col}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Missing columns: {missing}")

        for _, row in df.iterrows():
            gene = row[gene_col]
            cog = str(row[cog_col])
            description = row.get('Description', '')

            # COG category is a single letter or comma-separated
            categories = [c.strip() for c in cog if c.strip() in COG_CATEGORIES]

            payload = {
                'COG_category': cog,
                'COG_classes': categories,
                'Description': description,
                'Preferred_name': row.get('Preferred_name', ''),
                'KEGG_ko': str(row.get('KEGG_ko', '')),
                'KEGG_Pathway': str(row.get('KEGG_Pathway', '')),
                'original_protein_id': gene,
            }

            # Key annotations by the normalized (unversioned) RefSeq accession
            # so they align with versioned protein nodes in the graph, while
            # preserving the original id for traceability. The versioned form
            # is also registered so both id-spaces resolve.
            norm = normalize_wp_id(gene)
            self.annotations[norm] = payload
            if gene and gene not in self.annotations:
                self.annotations[gene] = payload

            for cat in categories:
                if cat in self.category_counts:
                    self.category_counts[cat] += 1

        return df

    def annotate_by_locus_tag(self, gene_ids: List[str]) -> Dict[str, dict]:
        """
        Create placeholder annotations based on locus tag patterns.
        In a production setting, this would call the eggNOG API.

        For the proof of concept, creates empty annotations that
        serve as a template for real data integration.

        Parameters
        ----------
        gene_ids : list of str

        Returns
        -------
        dict: {GeneID: {annotation dict}}
        """
        for gene in gene_ids:
            if gene not in self.annotations:
                self.annotations[gene] = {
                    'COG_category': 'S',  # Default: Function unknown
                    'COG_classes': ['S'],
                    'Description': '',
                    'Preferred_name': '',
                }
                self.category_counts['S'] += 1

        return self.annotations

    def get_genes_by_cog(self, cog_category: str) -> List[str]:
        """Return list of genes belonging to a COG category."""
        return [
            gene for gene, ann in self.annotations.items()
            if cog_category in ann.get('COG_classes', [])
        ]

    def as_locus_tag_map(self, gene_protein_map: Dict[str, str]) -> Dict[str, dict]:
        """Re-key annotations onto locus tags via the gene->protein mapping.

        The heterogeneous graph nodes genes are locus tags, while eggNOG keys
        are WP protein accessions (versioned or not). This resolves each
        locus_tag to its protein via ``gene_protein_map``, normalizes the
        protein accession, and looks up its COG annotation — so COG edges and
        hyperedges attach to real gene nodes, not orphan WP nodes.

        Parameters
        ----------
        gene_protein_map : dict
            {locus_tag: versioned_protein_id} (as built in ``run_integration``)

        Returns
        -------
        dict
            {locus_tag: {COG_classes, COG_category, Description, ...}} for the
            loci that map onto an annotated protein.
        """
        resolved: Dict[str, dict] = {}
        for locus, protein in gene_protein_map.items():
            if not protein:
                continue
            ann = self.annotations.get(protein) \
                or self.annotations.get(normalize_wp_id(protein))
            if ann is None:
                continue
            resolved[locus] = dict(ann)
        return resolved

    def kegg_pathway_membership(self,
                                gene_protein_map: Dict[str, str]) -> Dict[str, List[str]]:
        """Build {pathway_id: [locus_tag, ...]} from eggNOG KEGG_Pathway.

        Parses the comma-separated KEGG_Pathway column (ko00311, map01130, ...)
        carried through by `scripts/convert_emapper_to_tsv.py` and re-keys it
        onto locus tags via the gene->protein mapping, exactly like
        ``as_locus_tag_map`` does for COG. This gives acquired/plasmid-borne
        AMR determinants a KEGG pathway layer that the `eco`-bridged UniProt
        cross-reference route cannot provide.

        Parameters
        ----------
        gene_protein_map : dict
            {locus_tag: versioned_protein_id} (as built in ``run_integration``)

        Returns
        -------
        dict
            {pathway_id: [locus_tag, ...]}
        """
        membership: Dict[str, List[str]] = {}
        for locus, protein in gene_protein_map.items():
            if not protein:
                continue
            ann = self.annotations.get(protein) \
                or self.annotations.get(normalize_wp_id(protein))
            if ann is None:
                continue
            raw = ann.get('KEGG_Pathway', '')
            for pw in (p.strip() for p in str(raw).split(',') if p.strip()):
                if pw == '-' or pw.lower() == 'nan':
                    continue
                membership.setdefault(pw, []).append(locus)
        return membership

    def summary(self) -> dict:
        """Return COG annotation summary."""
        non_zero = {k: v for k, v in self.category_counts.items() if v > 0}
        return {
            'annotated_genes': len(self.annotations),
            'cog_categories_found': len(non_zero),
            'top_categories': dict(
                sorted(non_zero.items(), key=lambda x: -x[1])[:5]
            ),
        }
