"""
kegg_compound.py
Resolve MetaboLights ChEBI metabolite identifiers to KEGG compound entries,
then enrich with the enzymes (EC numbers) and KEGG pathways that act on or
involve each compound.

Mapping strategy (ChEBI API unreachable):
  1. Download the full KEGG ChEBI cross-reference table via
     `rest.kegg.jp/conv/chebi/cpd` (cpd -> chebi), invert to chebi -> [cpd].
  2. For each KEGG compound id, fetch the compound record
     (`rest.kegg.jp/get/cpd:Cxxxxx`) and parse the `ENZYME` and `PATHWAY`
     blocks.

All network results are cached to disk so reruns are offline.
"""

import re
import json
import time
import requests
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional

KEGG_BASE = "https://rest.kegg.jp"
KEGG_CACHE_DIR = Path(__file__).parent.parent / "outputs" / "cache"


class KEGGCompoundResolver:
    """
    Map ChEBI ids to KEGG compounds and their enzymes/pathways.

    Attributes
    ----------
    chebi_to_cpd : dict
        {chebi_accession: [kegg_compound_id, ...]}
    compound_info : dict
        {kegg_compound_id: {name, formula, enzymes: [EC...],
                            pathways: [pathway_id...]}}
    """

    def __init__(self, cache_dir: Optional[Path] = None):
        self.cache_dir = Path(cache_dir) if cache_dir else KEGG_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.chebi_to_cpd = {}
        self.compound_info = {}
        self._crossref_file = self.cache_dir / "kegg_chebi_cpd.json"
        self._compound_dir = self.cache_dir / "kegg_compounds"
        self._compound_dir.mkdir(parents=True, exist_ok=True)

    # ── Cross-reference table (chebi -> cpd) ─────────────────
    def _load_crossref(self) -> Dict[str, List[str]]:
        """Load the cached chebi->cpd map, downloading it if missing."""
        if self.chebi_to_cpd:
            return self.chebi_to_cpd
        if self._crossref_file.exists():
            with open(self._crossref_file) as f:
                self.chebi_to_cpd = json.load(f)
            return self.chebi_to_cpd

        print("  [KEGG compound] Downloading ChEBI cross-reference table...")
        text = self._kegg_get(f"{KEGG_BASE}/conv/chebi/cpd")
        rev: Dict[str, List[str]] = {}
        for line in text.strip().splitlines():
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            cpd, che = parts
            che_id = che.split(':')[1]
            rev.setdefault(che_id, []).append(cpd.split(':')[1])

        self.chebi_to_cpd = rev
        with open(self._crossref_file, 'w') as f:
            json.dump(rev, f)
        print(f"  [KEGG compound] {len(rev)} ChEBI->compound mappings cached")
        return rev

    # ── Compound record fetching ─────────────────────────────
    def _kegg_get(self, url: str, retries: int = 3) -> str:
        for attempt in range(retries):
            try:
                resp = requests.get(url, timeout=30)
                if resp.status_code == 200:
                    return resp.text
            except requests.RequestException:
                pass
            time.sleep(1 + attempt)
        return ''

    def _compound_cache_path(self, cpd_id: str) -> Path:
        return self._compound_dir / f"{cpd_id}.json"

    def _parse_compound(self, text: str) -> Dict:
        """Parse a KEGG compound record into info dict."""
        info = {'name': '', 'formula': '', 'enzymes': [], 'pathways': {}}
        block = ''
        for line in text.splitlines():
            if line.startswith('//'):
                break
            if not line or line[0] == ' ':
                continue
            key = line.split()[0]
            block = key
            if key == 'NAME':
                # first name line, strip trailing ';'
                name = line.replace('NAME', '', 1).strip().rstrip(';')
                info['name'] = name
            elif key == 'FORMULA':
                info['formula'] = line.replace('FORMULA', '', 1).strip()
        # Filter enzyme block: keep valid EC numbers (N.N.N.N), drop
        # artifacts like `2.4.1.-` or stray tokens
        def _valid_ec(tok: str) -> bool:
            return bool(re.match(r'^\d+\.\d+\.\d+\.\d+$', tok))

        # enzymes and pathways are on continuation lines
        in_block = None
        for line in text.splitlines():
            if line.startswith('//'):
                break
            if line and line[0] != ' ':
                key = line.split()[0]
                in_block = key if key in ('ENZYME', 'PATHWAY') else None
                if key == 'ENZYME':
                    rest = line[6:].strip()
                    if rest:
                        info['enzymes'].extend(rest.split())
                elif key == 'PATHWAY':
                    rest = line[6:].strip()
                    if rest:
                        parts = rest.split(maxsplit=1)
                        if len(parts) == 2:
                            info['pathways'][parts[0]] = parts[1]
                continue
            if in_block == 'ENZYME':
                info['enzymes'].extend(line.strip().split())
            elif in_block == 'PATHWAY':
                parts = line.strip().split(maxsplit=1)
                if len(parts) == 2:
                    info['pathways'][parts[0]] = parts[1]
        info['enzymes'] = [e for e in info['enzymes'] if _valid_ec(e)]
        info['pathways'] = {k: v for k, v in info['pathways'].items()
                            if k.lower().startswith('map')}
        return info

    def fetch_compound(self, cpd_id: str) -> Dict:
        """Fetch (and cache) a single KEGG compound record."""
        cache = self._compound_cache_path(cpd_id)
        if cache.exists():
            with open(cache) as f:
                return json.load(f)
        text = self._kegg_get(f"{KEGG_BASE}/get/cpd:{cpd_id}")
        info = self._parse_compound(text) if text else {}
        with open(cache, 'w') as f:
            json.dump(info, f)
        time.sleep(0.3)  # be polite to KEGG
        return info

    # ── Public API ───────────────────────────────────────────
    def resolve(self, chebi_ids: List[str]) -> Dict[str, Dict]:
        """
        Resolve a list of ChEBI accessions to KEGG compound metadata.

        Returns
        -------
        dict: {chebi_id: {kegg_id, name, formula, enzymes, pathways}}
        """
        self._load_crossref()
        result = {}
        for che_id in chebi_ids:
            che = str(che_id)
            if not re.match(r'^CHEBI:\d+$', che):
                continue
            accession = che.split(':')[1]
            cpd_ids = self.chebi_to_cpd.get(accession, [])
            if not cpd_ids:
                continue
            best = cpd_ids[0]  # first compound mapping
            info = self.fetch_compound(best)
            if not info:
                continue
            result[che] = {
                'kegg_id': best,
                'name': info.get('name', ''),
                'formula': info.get('formula', ''),
                'enzymes': info.get('enzymes', []),
                'pathways': info.get('pathways', {}),
            }
        self.resolved = result
        return result

    def find_by_name(self, names: List[str],
                     top_k: int = 1) -> Dict[str, Dict]:
        """
        Resolve metabolites by their common name via KEGG `find/compound`.

        Fallback for ChEBI-less or cross-reference-missing metabolites.
        The first (most relevant) KEGG compound per name is used.

        Parameters
        ----------
        names : list of str
            Metabolite display names (also used as result keys).
        top_k : int
            Number of candidate KEGG compounds to consider per name
            (only the first is used for enrichment).

        Returns
        -------
        dict: {name: {kegg_id, name, formula, enzymes, pathways}}
        """
        result = {}
        for name in names:
            if not name or not isinstance(name, str):
                continue
            query = name.replace(' ', '%20')
            text = self._kegg_get(f"{KEGG_BASE}/find/compound/{query}")

            # Parse all candidate KEGG compounds, keyed by their exact name.
            candidates = []  # (cpd_id, exact_name)
            for line in text.strip().splitlines():
                if not line.strip():
                    continue
                parts = line.split('\t')
                if parts and parts[0].startswith('C'):
                    candidates.append((parts[0], parts[1].strip() if len(parts) > 1 else ''))

            if not candidates:
                continue

            target = name.strip().lower()
            cpd_id = None
            # 1) Prefer an exact (case-insensitive) name match.
            for cid, cname in candidates:
                if cname.lower() == target:
                    cpd_id = cid
                    break
            # 2) Otherwise accept a candidate whose name starts with the query.
            if cpd_id is None:
                for cid, cname in candidates:
                    if cname.lower().startswith(target):
                        cpd_id = cid
                        break
            # 3) Last resort: first candidate (as before).
            if cpd_id is None and candidates:
                cpd_id = candidates[0][0]

            info = self.fetch_compound(cpd_id)
            if not info:
                continue
            result[name] = {
                'kegg_id': cpd_id,
                'name': info.get('name', ''),
                'formula': info.get('formula', ''),
                'enzymes': info.get('enzymes', []),
                'pathways': info.get('pathways', {}),
            }
        return result

    def summary(self) -> dict:
        """Return resolver summary."""
        return {
            'resolved_compounds': len(getattr(self, 'resolved', {})),
            'crossref_entries': len(self.chebi_to_cpd),
            'fetched_compound_records': len(list(self._compound_dir.glob('*.json'))),
        }
