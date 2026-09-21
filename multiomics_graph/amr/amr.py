"""
amr/amr.py
Load the curated antimicrobial-resistance (AMR) manifest and expose AMR
mechanisms as a biological annotation layer on the multi-omics graph.

The manifest (`amr_manifest.json`) records known resistance determinants per
strain as curated from the paper and verified against each strain's NCBI
genomic.gff. AMR knowledge is attached to the graph as *annotation*, not as a
predictive input: markers remain represented even when they have no
RNA/protein quantification (a genome-only resistance gene is biologically
valid).
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Union

AMR_DIR = Path(__file__).parent
DEFAULT_MANIFEST = AMR_DIR / 'amr_manifest.json'

# Canonical AMR mechanism labels used for the AMR_mechanism node type.
AMR_CLASS_NAMES = {
    'beta_lactam': 'Beta-lactam resistance',
    'aminoglycoside': 'Aminoglycoside resistance',
    'sulfonamide': 'Sulfonamide resistance',
    'trimethoprim': 'Trimethoprim resistance',
    'tetracycline': 'Tetracycline resistance',
    'macrolide': 'Macrolide resistance',
    'fluoroquinolone': 'Fluoroquinolone resistance',
    'phenicol': 'Phenicol resistance',
    'fosfomycin': 'Fosfomycin resistance',
}


def _to_list(value: Union[str, List[str], None]) -> List[str]:
    """Normalize a manifest value to a list of strings."""
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v) for v in value]
    return [str(value)]


class AMRManifest:
    """Load and query the per-strain AMR manifest."""

    def __init__(self, manifest_path: Optional[Path] = None):
        self.path = Path(manifest_path) if manifest_path else DEFAULT_MANIFEST
        self.data: Dict[str, dict] = self._load()

    def _load(self) -> Dict[str, dict]:
        if not self.path.exists():
            raise FileNotFoundError(f"AMR manifest not found: {self.path}")
        return json.loads(self.path.read_text(encoding='utf-8'))

    def strains(self) -> List[str]:
        return list(self.data.keys())

    def markers(self, strain: str) -> List[dict]:
        """Return normalized marker records for a strain."""
        return [self._normalize(name, rec)
                for name, rec in self.data.get(strain, {}).items()]

    def _normalize(self, marker_name: str, rec: dict) -> dict:
        return {
            'name': rec.get('name', marker_name),
            'amr_class': rec.get('amr_class', 'unknown'),
            'locus_tags': _to_list(rec.get('locus_tag')),
            'protein_ids': _to_list(rec.get('protein_id')),
            'source': rec.get('source', ''),
            'note': rec.get('note', ''),
        }

    def locus_to_amr(self, strain: str) -> Dict[str, dict]:
        """Map each locus_tag to its AMR marker record (one-to-one)."""
        mapping: Dict[str, dict] = {}
        for rec in self.markers(strain):
            for locus in rec['locus_tags']:
                mapping[locus] = rec
        return mapping

    def class_to_loci(self, strain: str) -> Dict[str, List[str]]:
        """Group AMR locus tags by mechanism class."""
        groups: Dict[str, List[str]] = {}
        for rec in self.markers(strain):
            groups.setdefault(rec['amr_class'], [])
            for locus in rec['locus_tags']:
                groups[rec['amr_class']].append(locus)
        return groups

    def amr_loci(self, strain: str) -> List[str]:
        """All AMR locus tags for a strain (deduplicated, ordered)."""
        seen: List[str] = []
        for rec in self.markers(strain):
            for locus in rec['locus_tags']:
                if locus and locus not in seen:
                    seen.append(locus)
        return seen

    def class_name(self, amr_class: str) -> str:
        return AMR_CLASS_NAMES.get(amr_class, amr_class.replace('_', ' '))

    def summary(self, strain: str) -> dict:
        markers = self.markers(strain)
        classes = {}
        for rec in markers:
            classes[rec['amr_class']] = classes.get(rec['amr_class'], 0) + 1
        return {
            'strain': strain,
            'n_markers': len(markers),
            'n_loci': len(self.amr_loci(strain)),
            'classes': classes,
        }