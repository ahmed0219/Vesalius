"""
string.py
Retrieve protein-protein interaction data from the STRING database.

Provides:
  - PPI edges with confidence scores for bacterial proteins
  - STRING protein identifiers mapped from RefSeq WP_ accessions
  - Interaction network export for graph construction

Uses the STRING REST API with local caching.
"""

import pandas as pd
import numpy as np
import requests
import time
import json
from pathlib import Path
from typing import List, Optional


STRING_BASE = "https://string-db.org/api"
STRING_JSON = f"{STRING_BASE}/json"

CACHE_DIR = Path(__file__).parent.parent / "outputs" / "cache"


class STRINGClient:
    """
    Retrieve protein-protein interactions from STRING DB.

    Maps RefSeq WP_* accessions to STRING identifiers and
    retrieves interaction networks with confidence scores.

    Attributes
    ----------
    species_id : int
        NCBI taxonomy ID for STRING (511145 = E. coli K-12 MG1655)
    confidence_threshold : float
        Minimum STRING interaction score. API uses 0-1000 scale internally
        (400 = medium confidence). The returned scores are 0-1 floats.
    interactions : pd.DataFrame
        Retrieved PPI edges
    """

    def __init__(self, species_id: int = 511145,
                 confidence_threshold: float = 400,
                 cache_dir: Optional[Path] = None):
        """
        Parameters
        ----------
        species_id : int
            STRING species identifier (511145 = E. coli K-12)
        confidence_threshold : float
            Minimum interaction score. STRING API uses 0-1000 scale
            (400 = medium confidence). Returned scores are 0-1 floats.
        cache_dir : Path, optional
        """
        self.species_id = species_id
        self.confidence_threshold = confidence_threshold
        self.cache_dir = Path(cache_dir) if cache_dir else CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.interactions = pd.DataFrame()
        self._string_id_map = {}

    def _cache_path(self, key: str) -> Path:
        safe = key.replace('/', '_').replace('.', '_')
        return self.cache_dir / f"string_{safe}.json"

    def _load_cache(self, key: str):
        path = self._cache_path(key)
        if path.exists():
            with open(path, 'r') as f:
                return json.load(f)
        return None

    def _save_cache(self, key: str, data):
        with open(self._cache_path(key), 'w') as f:
            json.dump(data, f, indent=2)

    def get_string_ids(self, protein_ids: List[str]) -> dict:
        """
        Map RefSeq WP_* proteins to STRING identifiers.

        Parameters
        ----------
        protein_ids : list of str

        Returns
        -------
        dict: {original_id: string_id}
        """
        cache_key = f"string_ids_{self.species_id}"
        cached = self._load_cache(cache_key)
        if cached:
            self._string_id_map = cached
            return cached

        mapping = {}

        # STRING accepts up to 2000 identifiers per call
        for i in range(0, len(protein_ids), 1000):
            batch = protein_ids[i:i + 1000]
            try:
                resp = requests.post(
                    f"{STRING_JSON}/get_string_ids",
                    data={
                        'identifiers': '\r\n'.join(batch),
                        'species': self.species_id,
                        'limit': 1,
                        'echo_query': 1,
                    },
                    timeout=30,
                )
                if resp.status_code != 200:
                    continue

                for entry in resp.json():
                    query = entry.get('queryItem') or entry.get('queryString', '')
                    string_id = entry.get('stringId', '')
                    if query and string_id:
                        mapping[query] = string_id

            except Exception as e:
                print(f"  [STRING] Warning: batch failed: {e}")
                continue

        self._string_id_map = mapping
        self._save_cache(cache_key, mapping)
        return mapping

    def fetch_interactions(self, protein_ids: List[str]) -> pd.DataFrame:
        """
        Fetch PPI network for the given protein identifiers.

        Parameters
        ----------
        protein_ids : list of str
            UniProt accessions (preferred — these are passed from main.py)

        Returns
        -------
        pd.DataFrame with columns:
            ProteinID_A, ProteinID_B, combined_score,
            experimental_score, database_score, ...
        """
        cache_key = f"ppi_{self.species_id}_{len(protein_ids)}_{int(self.confidence_threshold)}"
        cached = self._load_cache(cache_key)
        if cached is not None:
            self.interactions = pd.DataFrame(cached)
            return self.interactions

        # Map to STRING IDs
        str_map = self.get_string_ids(protein_ids)
        string_ids = list(str_map.values())

        if not string_ids:
            print("  [STRING] No identifiers mapped")
            return pd.DataFrame()

        # Fetch network in batches
        all_edges = []
        for i in range(0, len(string_ids), 500):
            batch = string_ids[i:i + 500]
            try:
                resp = requests.post(
                    f"{STRING_JSON}/network",
                    data={
                        'identifiers': '\r\n'.join(batch),
                        'species': self.species_id,
                        'required_score': self.confidence_threshold,
                    },
                    timeout=60,
                )
                if resp.status_code != 200:
                    continue

                for edge in resp.json():
                    sid_a = edge.get('stringId_A', '')
                    sid_b = edge.get('stringId_B', '')
                    score = float(edge.get('score', 0))

                    all_edges.append({
                        'ProteinID_A': sid_a,
                        'ProteinID_B': sid_b,
                        'combined_score': score,
                        'experimental_score': float(
                            edge.get('escore', 0)
                        ),
                        'database_score': float(
                            edge.get('dscore', 0)
                        ),
                        'coexpression_score': float(
                            edge.get('tscore', 0)
                        ),
                        'cooccurrence_score': float(
                            edge.get('pscore', 0)
                        ),
                    })

            except Exception as e:
                print(f"  [STRING] Warning: network fetch failed: {e}")
                continue

        result = pd.DataFrame(all_edges)
        if not result.empty:
            # Sort protein pair so (A,B) and (B,A) deduplicate as the same edge,
            # keeping the entry with the highest combined_score per pair
            result['_pair_key'] = result.apply(
                lambda r: tuple(sorted([str(r['ProteinID_A']), str(r['ProteinID_B'])])), axis=1
            )
            result = result.loc[
                result.groupby('_pair_key')['combined_score'].idxmax()
            ].drop(columns=['_pair_key']).reset_index(drop=True)

        # Map STRING IDs back to original query IDs (UniProt accessions)
        reverse_map = {v: k for k, v in str_map.items()}
        if not result.empty and 'ProteinID_A' in result.columns:
            result['ProteinID_A'] = result['ProteinID_A'].map(reverse_map).fillna(
                result['ProteinID_A']
            )
            result['ProteinID_B'] = result['ProteinID_B'].map(reverse_map).fillna(
                result['ProteinID_B']
            )

        self.interactions = result
        self._save_cache(cache_key, result.to_dict('records'))
        return result

    def summary(self) -> dict:
        """Return summary of PPI network."""
        if self.interactions.empty:
            return {'ppi_edges': 0}
        return {
            'ppi_edges': len(self.interactions),
            'unique_proteins': len(set(
                list(self.interactions['ProteinID_A'].unique()) +
                list(self.interactions['ProteinID_B'].unique())
            )),
            'mean_confidence': float(
                self.interactions['combined_score'].mean()
            ),
        }
