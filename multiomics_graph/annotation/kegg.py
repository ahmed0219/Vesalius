"""
kegg.py
Retrieve KEGG pathway annotations for bacterial genes.

KEGG (Kyoto Encyclopedia of Genes and Genomes) provides:
  - Pathway maps (e.g., eco00010: Glycolysis)
  - KO (KEGG Orthology) identifiers
  - Module definitions

This module maps bacterial genes to KEGG pathways via KO identifiers
and builds pathway membership lists for hypergraph construction.

Uses the KEGG REST API with local caching.
"""

import pandas as pd
import numpy as np
import requests
import time
import json
from pathlib import Path
from typing import List, Dict, Optional
import re


KEGG_BASE = "https://rest.kegg.jp"
KEGG_CACHE_DIR = Path(__file__).parent.parent / "outputs" / "cache"


class KEGGAnnotator:
    """
    Retrieve KEGG pathway annotations for bacterial genes.

    Maps bacterial genes to pathways via:
      Gene → KO (KEGG Orthology) → Pathway

    Attributes
    ----------
    organism_code : str
        KEGG organism code (e.g., 'eco' for E. coli K-12)
    pathway_membership : dict
        {pathway_id: [list of GeneIDs]}
    ko_to_genes : dict
        {KO_id: [list of GeneIDs]}
    pathways : dict
        {pathway_id: {name, description, class}}
    """

    def __init__(self, organism_code: str = "eco",
                 cache_dir: Optional[Path] = None,
                 strain_name: Optional[str] = None):
        """
        Parameters
        ----------
        organism_code : str
            KEGG organism code (eco = E. coli K-12 MG1655).
            For B36, try 'ecb' (E. coli) or closest relative.
        cache_dir : Path, optional
        strain_name : str, optional
            Name of the active strain, used for strain-specific caching.
        """
        self.organism_code = organism_code
        self.cache_dir = Path(cache_dir) if cache_dir else KEGG_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.strain_name = strain_name

        self.pathway_membership = {}
        self.ko_to_genes = {}
        self.pathways = {}
        self.gene_to_pathways = {}

    def _cache_path(self, key: str) -> Path:
        safe = key.replace('/', '_').replace('.', '_')
        return self.cache_dir / f"kegg_{safe}.json"

    def _load_cache(self, key: str):
        path = self._cache_path(key)
        if path.exists():
            with open(path, 'r') as f:
                return json.load(f)
        return None

    def _save_cache(self, key: str, data):
        with open(self._cache_path(key), 'w') as f:
            json.dump(data, f, indent=2)

    def _kegg_get(self, url: str, retries: int = 3) -> str:
        """Make KEGG REST API call with rate limiting and retries."""
        for attempt in range(retries):
            try:
                resp = requests.get(url, timeout=30)
                if resp.status_code == 200:
                    return resp.text
            except requests.RequestException:
                pass
            time.sleep(1 + attempt)
        return ''

    def list_pathways(self) -> Dict[str, str]:
        """
        List all KEGG pathways for the organism.

        Returns
        -------
        dict: {pathway_id: pathway_name}
        """
        cache_key = f"list_pathways_{self.organism_code}"
        cached = self._load_cache(cache_key)
        if cached:
            self.pathways = cached
            return cached

        text = self._kegg_get(f"{KEGG_BASE}/list/pathway/{self.organism_code}")
        pathways = {}
        for line in text.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    path_id = parts[0]
                    name = parts[1]
                    pathways[path_id] = {'name': name, 'pathway_id': path_id}

        self.pathways = pathways
        self._save_cache(cache_key, pathways)
        return pathways

    def get_ko_for_genes(
        self, gene_ids: List[str],
        gene_to_kegg_map: Optional[Dict[str, str]] = None,
    ) -> Dict[str, str]:
        """
        Map bacterial gene identifiers to KO (KEGG Orthology) numbers.

        Uses KEGG gene-to-KO mapping for the organism.

        Parameters
        ----------
        gene_ids : list of str
            Bacterial gene identifiers (locus tags)
        gene_to_kegg_map : dict, optional
            Maps input gene IDs to KEGG gene IDs (e.g., {'EW036_RS00010': 'b0001'}).
            Used when the input IDs differ from KEGG's native naming.

        Returns
        -------
        dict: {gene_id: KO_id}
        """
        cache_prefix = f"gene2ko_{self.strain_name}_" if getattr(self, 'strain_name', None) else "gene2ko_"
        cache_key = f"{cache_prefix}{self.organism_code}"
        cached = self._load_cache(cache_key)
        if cached:
            self.ko_to_genes = cached
            return cached

        # Get all KO assignments for the organism
        text = self._kegg_get(
            f"{KEGG_BASE}/link/ko/{self.organism_code}"
        )
        kegg_gene_to_ko = {}
        for line in text.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                gene_ko = parts[0]
                ko = parts[1].split(':')[-1] if ':' in parts[1] else parts[1]
                # Extract gene ID from the full identifier (format: eco:b0001)
                gene_id_parts = gene_ko.split(':')
                if len(gene_id_parts) >= 2:
                    gene_id = gene_id_parts[1]
                    if ko.startswith('K'):
                        kegg_gene_to_ko[gene_id] = ko

        result = {}
        for gid in gene_ids:
            # Direct match
            if gid in kegg_gene_to_ko:
                result[gid] = kegg_gene_to_ko[gid]
                continue
            # Via gene-to-kegg bridge map
            if gene_to_kegg_map and gid in gene_to_kegg_map:
                kegg_id = gene_to_kegg_map[gid]
                if kegg_id in kegg_gene_to_ko:
                    result[gid] = kegg_gene_to_ko[kegg_id]
                    continue
            # Try stripped short name
            short = gid.split('_')[-1] if '_' in gid else gid
            if short in kegg_gene_to_ko:
                result[gid] = kegg_gene_to_ko[short]

        self.ko_to_genes = result
        self._save_cache(cache_key, result)
        return result

    def map_ko_to_pathways(self, ko_ids: List[str]) -> Dict[str, List[str]]:
        """
        Map KO numbers to KEGG pathway IDs via batch API.

        Parameters
        ----------
        ko_ids : list of str
            KO identifiers (e.g., K02305)

        Returns
        -------
        dict: {KO_id: [pathway_id, ...]}
        """
        cache_key = f"ko2pathway_{self.organism_code}"
        cached = self._load_cache(cache_key)
        if cached:
            return cached

        BATCH = 50  # KOs per batch request
        ko_to_pathways = {}
        for i in range(0, len(ko_ids), BATCH):
            batch = ko_ids[i:i + BATCH]
            query = "+".join(f"ko:{ko}" for ko in batch)
            text = self._kegg_get(f"{KEGG_BASE}/link/pathway/{query}")
            for line in text.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    # Format: ko:K02305\tpath:eco00240
                    ko = parts[0].split(':')[-1] if ':' in parts[0] else parts[0]
                    path_id = parts[1].strip()
                    if path_id.startswith('path:'):
                        path_id = path_id.replace('path:', '')
                    if ko not in ko_to_pathways:
                        ko_to_pathways[ko] = []
                    if path_id not in ko_to_pathways[ko]:
                        ko_to_pathways[ko].append(path_id)

        self._save_cache(cache_key, ko_to_pathways)
        return ko_to_pathways

    def build_pathway_membership(
        self, gene_ids: List[str],
        gene_to_kegg_map: Optional[Dict[str, str]] = None,
    ) -> Dict[str, List[str]]:
        """
        Build complete pathway → genes membership.

        1. Map genes → KO
        2. Map KO → pathways
        3. Invert to pathway → genes

        Parameters
        ----------
        gene_ids : list of str
        gene_to_kegg_map : dict, optional
            Maps input gene IDs to KEGG gene IDs (e.g., {'EW036_RS00010': 'b0001'}).

        Returns
        -------
        dict: {pathway_id: [list of GeneIDs]}
        """
        cache_prefix = f"pathway_membership_{self.strain_name}_" if getattr(self, 'strain_name', None) else "pathway_membership_"
        cache_key = f"{cache_prefix}{self.organism_code}"
        cached = self._load_cache(cache_key)
        if cached:
            self.pathway_membership = cached
            return cached

        # Step 1: Genes → KO
        gene_to_ko = self.get_ko_for_genes(gene_ids, gene_to_kegg_map)
        kos = list(set(gene_to_ko.values()))

        # Step 2: KO → Pathways
        ko_to_pathways = self.map_ko_to_pathways(kos)

        # Step 3: Invert to Pathway → Genes
        pathway_genes = {}
        for gene, ko in gene_to_ko.items():
            for pathway in ko_to_pathways.get(ko, []):
                if pathway not in pathway_genes:
                    pathway_genes[pathway] = []
                if gene not in pathway_genes[pathway]:
                    pathway_genes[pathway].append(gene)

        # Also fetch pathway names
        self.list_pathways()

        # Build gene → pathway mapping
        self.gene_to_pathways = {}
        for pw, genes in pathway_genes.items():
            for g in genes:
                if g not in self.gene_to_pathways:
                    self.gene_to_pathways[g] = []
                self.gene_to_pathways[g].append(pw)

        self.pathway_membership = pathway_genes
        self._save_cache(cache_key, pathway_genes)
        return pathway_genes

    def get_pathway_name(self, pathway_id: str) -> str:
        """Return human-readable pathway name."""
        if pathway_id in self.pathways:
            return self.pathways[pathway_id].get('name', pathway_id)
        # Fall back to API
        text = self._kegg_get(f"{KEGG_BASE}/get/{pathway_id}")
        for line in text.split('\n'):
            if line.startswith('NAME'):
                return line.replace('NAME', '').strip()
        return pathway_id

    def summary(self) -> dict:
        """Return KEGG annotation summary."""
        return {
            'organism': self.organism_code,
            'pathways': len(self.pathway_membership),
            'genes_with_pathways': len(self.gene_to_pathways),
            'ko_assigned': len(self.ko_to_genes),
        }
