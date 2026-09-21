"""
uniprot.py
Retrieve functional annotations from UniProt for bacterial proteins.

Maps RefSeq WP_ accessions -> UniProt accessions -> functional annotations.
Extracts GO terms (BP/MF/CC), EC numbers, protein names, domains, pathways.

Uses the UniProt REST API with local caching, parallel downloads,
and proper rate limiting.
"""

import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import requests

logger = logging.getLogger(__name__)

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
CACHE_DIR = Path(__file__).resolve().parent.parent / "outputs" / "cache"
DEFAULT_NCBI_TAX_ID = 562
BATCH_SIZE = 30  # IDs per search batch (max ~32 OR conditions)
MAX_RETRIES = 5


@dataclass
class GoTerm:
    id: str
    name: str
    aspect: str


@dataclass
@dataclass
class MappingStats:
    submitted: int = 0
    mapped: int = 0
    candidates_total: int = 0
    failed_batches: int = 0
    elapsed: float = 0.0

    @property
    def coverage(self) -> float:
        return self.mapped / max(self.submitted, 1)


@dataclass
class AnnotationStats:
    requested: int = 0
    mapped: int = 0
    annotated: int = 0
    swissprot: int = 0
    trembl: int = 0
    avg_go: float = 0.0
    avg_keywords: float = 0.0
    missing: int = 0
    elapsed: float = 0.0


class _CacheManager:
    """Manages persistent cache of UniProt API responses."""

    def __init__(self, cache_dir: Path):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _path(self, key: str) -> Path:
        safe = key.replace("/", "_").replace(".", "_")
        return self.cache_dir / f"uniprot_{safe}.json"

    def load(self, key: str) -> Optional[Any]:
        path = self._path(key)
        if path.exists() and path.stat().st_size > 4:
            try:
                with open(path, "r") as f:
                    return json.load(f)
            except (json.JSONDecodeError, ValueError):
                logger.warning("Corrupt cache file: %s", path)
                return None
        return None

    def save(self, key: str, data: Any):
        path = self._path(key)
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    def invalidate(self, key: str):
        path = self._path(key)
        if path.exists():
            path.unlink()

    def is_valid(self, key: str) -> bool:
        path = self._path(key)
        if not path.exists():
            return False
        if path.stat().st_size <= 4:
            return False
        try:
            data = self.load(key)
            return data is not None and len(data) > 0 if isinstance(data, dict) else bool(data)
        except (json.JSONDecodeError, ValueError):
            return False


def _request_with_retry(
    method: str, url: str, retries: int = MAX_RETRIES, **kwargs: Any
) -> requests.Response:
    """Make an HTTP request with exponential backoff retry."""
    for attempt in range(retries):
        try:
            resp = requests.request(method, url, timeout=60, **kwargs)
            if resp.status_code == 429:
                wait = min(2 ** attempt * 5, 120)
                logger.warning("Rate limited (429), waiting %ds", wait)
                time.sleep(wait)
                continue
            resp.raise_for_status()
            return resp
        except requests.RequestException as e:
            if attempt == retries - 1:
                raise
            wait = 2 ** attempt
            logger.debug("Request failed (attempt %d/%d): %s", attempt + 1, retries, e)
            time.sleep(wait)
    raise RuntimeError(f"Request failed after {retries} retries: {url}")


class _AnnotationParser:
    """Parses UniProtKB JSON responses into structured records."""

    @staticmethod
    def _safe_get(data: dict, *keys: str, default: str = "") -> str:
        for key in keys:
            if isinstance(data, dict):
                data = data.get(key, {})
            else:
                return default
        return str(data) if not isinstance(data, (dict, list)) else default

    @staticmethod
    def _get_gene_name(data: dict) -> str:
        genes = data.get("genes") or []
        if genes:
            return genes[0].get("geneName", {}).get("value", "")
        return ""

    @staticmethod
    def _get_protein_name(data: dict) -> str:
        desc = data.get("proteinDescription", {})
        rec = desc.get("recommendedName") or desc.get("submittedName") or {}
        return rec.get("fullName", {}).get("value", "")

    @staticmethod
    def _get_alternative_names(data: dict) -> str:
        alt_names: List[str] = []
        desc = data.get("proteinDescription", {})
        for alt in desc.get("alternativeNames") or []:
            val = alt.get("fullName", {}).get("value", "")
            if val:
                alt_names.append(val)
        return "; ".join(alt_names)

    @staticmethod
    def _get_comment(data: dict, comment_type: str) -> str:
        comments = data.get("comments") or []
        values: List[str] = []
        for c in comments:
            if c.get("commentType") == comment_type:
                for t in c.get("texts") or []:
                    val = t.get("value", "")
                    if val:
                        values.append(val)
        return " ".join(values)

    @staticmethod
    def _get_ec_numbers(data: dict) -> str:
        ec_set: Set[str] = set()
        for db in data.get("uniProtKBCrossReferences") or []:
            if db.get("database") == "EC":
                ec_set.add(db.get("id", ""))
        comments = data.get("comments") or []
        for c in comments:
            if c.get("commentType") in ("CATALYTIC ACTIVITY",):
                reaction = c.get("reaction") or {}
                reactions = reaction if isinstance(reaction, list) else [reaction]
                for r in reactions:
                    if isinstance(r, dict):
                        ec = r.get("ecNumber", "")
                        if ec:
                            ec_set.add(ec)
        return "; ".join(sorted(ec_set))

    @staticmethod
    def _get_subcellular_location(data: dict) -> str:
        comments = data.get("comments") or []
        locs: List[str] = []
        for c in comments:
            if c.get("commentType") == "SUBCELLULAR LOCATION":
                for loc in c.get("subcellularLocations") or []:
                    val = loc.get("location", {}).get("value", "")
                    if val:
                        locs.append(val)
        return "; ".join(locs)

    @staticmethod
    def _get_keywords(data: dict) -> List[str]:
        return [kw["name"] for kw in data.get("keywords") or []]

    @staticmethod
    def _get_go_terms(data: dict) -> Dict[str, List[GoTerm]]:
        result: Dict[str, List[GoTerm]] = {"BP": [], "MF": [], "CC": []}
        # Try field-named keys (used with fields param)
        aspect_map = {"go_p": "BP", "go_f": "MF", "go_c": "CC"}
        found = False
        for field, aspect in aspect_map.items():
            if field in data:
                found = True
                for entry in data.get(field) or []:
                    go_id = entry.get("goId", "")
                    if go_id:
                        result[aspect].append(
                            GoTerm(id=go_id, name=entry.get("name", ""), aspect=aspect)
                        )
        if found:
            return result
        # Fallback: extract GO from uniProtKBCrossReferences
        go_prefix = {"C": "CC", "F": "MF", "P": "BP"}
        for db in data.get("uniProtKBCrossReferences") or []:
            if db.get("database") == "GO":
                go_id = db.get("id", "")
                if not go_id:
                    continue
                aspect = "BP"
                for prop in db.get("properties") or []:
                    if prop.get("key") == "GoTerm":
                        val = prop.get("value", "")
                        if val and val[0] in go_prefix:
                            aspect = go_prefix[val[0]]
                        break
                result[aspect].append(
                    GoTerm(id=go_id, name=val if val else "", aspect=aspect)
                )
        return result

    @staticmethod
    def _get_cross_references(data: dict) -> Dict[str, List[str]]:
        xrefs: Dict[str, List[str]] = {}
        for db in data.get("uniProtKBCrossReferences") or []:
            db_name = db.get("database", "")
            db_id = db.get("id", "")
            if db_name and db_id:
                xrefs.setdefault(db_name, []).append(db_id)
        return xrefs

    @staticmethod
    def _get_organism(data: dict) -> Tuple[str, str]:
        org = data.get("organism", {})
        name = org.get("scientificName", "")
        tax_id = ""
        for db in org.get("lineages") or []:
            pass
        for db in org.get("dbReferences") or []:
            if db.get("type") == "NCBI Taxonomy":
                tax_id = db.get("id", "")
        return name, tax_id

    @staticmethod
    def _get_protein_length(data: dict) -> int:
        seq = data.get("sequence", {})
        return int(seq.get("length", 0)) if seq else 0

    @staticmethod
    def _get_annotation_score(data: dict) -> float:
        return float(data.get("annotationScore", 0))

    @staticmethod
    def is_reviewed(data: dict) -> bool:
        return data.get("entryType", "").startswith("UniProtKB reviewed")


class UniProtAnnotator:
    """
    Retrieve and cache UniProt annotations for bacterial RefSeq proteins.

    Maps RefSeq WP_* accessions to UniProt entries -> rich functional annotations.
    Handles multiple UniProt candidates per RefSeq, selects the best one.

    Public interface:
        annotator = UniProtAnnotator()
        df = annotator.fetch_annotations(wp_ids)
    """

    def __init__(
        self,
        organism_id: str = "562",
        cache_dir: Optional[Path] = None,
    ):
        self.organism_id = organism_id
        self.cache = _CacheManager(Path(cache_dir) if cache_dir else CACHE_DIR)
        self.go_terms: Dict[str, GoTerm] = {}
        self.mapping_stats: Optional[MappingStats] = None
        self.annotation_stats: Optional[AnnotationStats] = None

    # ------------------------------------------------------------------
    # Identifier normalization
    # ------------------------------------------------------------------

    @staticmethod
    def _normalize_wp(wp: str) -> str:
        """Strip version suffix from RefSeq accession (WP_000673464.1 -> WP_000673464)."""
        return wp.split(".")[0].strip() if wp else ""

    @staticmethod
    def _normalize_wp_list(ids: List[str]) -> List[str]:
        seen: Set[str] = set()
        result: List[str] = []
        for raw in ids:
            norm = UniProtAnnotator._normalize_wp(raw)
            if norm and norm not in seen:
                seen.add(norm)
                result.append(norm)
        return result

    # ------------------------------------------------------------------
    # Mapping: RefSeq -> UniProt (via search API)
    # ------------------------------------------------------------------

    def _search_refseq_batch(self, batch: List[str]) -> List[dict]:
        """Search UniProtKB for entries cross-referencing any of the given RefSeq IDs."""
        query = " OR ".join(f"(database:RefSeq {pid})" for pid in batch)
        try:
            resp = _request_with_retry(
                "GET",
                f"{UNIPROT_BASE}/search",
                params={"query": query, "size": 500, "format": "json"},
            )
            return resp.json().get("results", [])
        except Exception as e:
            logger.warning("Batch search failed for %d IDs: %s", len(batch), e)
            return []

    def _build_refseq_index(self, entries: List[dict]) -> Dict[str, List[Tuple[str, dict]]]:
        """Build {norm_wp: [(uniprot_id, entry)]} from search results."""
        index: Dict[str, List[Tuple[str, dict]]] = {}
        for entry in entries:
            acc = entry.get("primaryAccession", "")
            if not acc:
                continue
            for xref in entry.get("uniProtKBCrossReferences") or []:
                if xref.get("database") == "RefSeq":
                    wp = xref.get("id", "").split(".")[0].strip()
                    if wp:
                        index.setdefault(wp, []).append((acc, entry))
        return index

    def _resolve_candidates(
        self, candidates: List[Tuple[str, dict]]
    ) -> Optional[str]:
        """Select the best UniProt entry from multiple candidates."""
        if not candidates:
            return None
        if len(candidates) == 1:
            return candidates[0][0]

        def rank_key(cand: Tuple[str, dict]) -> tuple:
            uniprot_id, entry = cand
            _, tax_id = _AnnotationParser._get_organism(entry)
            organism_match = 1 if tax_id == self.organism_id else 0
            reviewed = _AnnotationParser.is_reviewed(entry)
            score = _AnnotationParser._get_annotation_score(entry)
            return (organism_match, reviewed, score, uniprot_id)

        return max(candidates, key=rank_key)[0]

    def map_refseq_to_uniprot(
        self, wp_accessions: List[str], batch_size: int = BATCH_SIZE
    ) -> Dict[str, str]:
        """
        Map RefSeq WP_* accessions to UniProt accessions via search API.

        Returns {WP_accession: UniProt_accession}.
        """
        start = time.time()
        normalized = self._normalize_wp_list(wp_accessions)

        cache_key = f"refseq2uniprot_{self.organism_id}"

        cached = self.cache.load(cache_key)
        if cached and isinstance(cached, dict) and len(cached) > 0:
            logger.info("Loaded mapping from cache: %d RefSeq -> UniProt", len(cached))
            self.mapping_stats = MappingStats(
                submitted=len(normalized), mapped=len(cached),
                elapsed=time.time() - start,
            )
            return cached

        # Search in batches
        all_entries: List[dict] = []
        failed_batches = 0
        n_batches = (len(normalized) + batch_size - 1) // batch_size

        for i in range(0, len(normalized), batch_size):
            batch = normalized[i : i + batch_size]
            batch_idx = i // batch_size + 1
            logger.info("Search batch %d/%d (%d IDs)", batch_idx, n_batches, len(batch))
            entries = self._search_refseq_batch(batch)
            if entries:
                all_entries.extend(entries)
            else:
                failed_batches += 1
            logger.info("Search progress: %d/%d (%.0f%%)",
                        i + len(batch), len(normalized),
                        (i + len(batch)) / max(len(normalized), 1) * 100)

        # Build index: norm_wp -> [(uniprot_id, entry)]
        refseq_index = self._build_refseq_index(all_entries)

        # Resolve best candidate per input ID
        final_map: Dict[str, str] = {}
        for wp_id in normalized:
            candidates = refseq_index.get(wp_id, [])
            best = self._resolve_candidates(candidates)
            if best:
                final_map[wp_id] = best

        # Cache entries for later use by fetch_annotations
        entry_cache_key = f"entries_{self.organism_id}"
        self.cache.save(entry_cache_key, all_entries)

        elapsed = time.time() - start
        self.mapping_stats = MappingStats(
            submitted=len(normalized), mapped=len(final_map),
            candidates_total=sum(len(v) for v in refseq_index.values()),
            failed_batches=failed_batches, elapsed=elapsed,
        )

        logger.info(
            "Mapping: %d/%d proteins mapped (%.1f%%), %d candidates, %d failed batches, %.1fs",
            self.mapping_stats.mapped, self.mapping_stats.submitted,
            self.mapping_stats.coverage * 100,
            self.mapping_stats.candidates_total,
            self.mapping_stats.failed_batches,
            self.mapping_stats.elapsed,
        )

        if final_map:
            self.cache.save(cache_key, final_map)
        else:
            logger.warning("Not caching empty mapping result.")

        return final_map

    # ------------------------------------------------------------------
    # Annotation parsing
    # ------------------------------------------------------------------

    def _parse_entry(self, entry: dict) -> Optional[dict]:
        """Parse a single UniProtKB entry into a flat record."""
        try:
            primary_acc = entry.get("primaryAccession", "")
            parser = _AnnotationParser()

            gene_name = parser._get_gene_name(entry)
            protein_name = parser._get_protein_name(entry)
            alt_names = parser._get_alternative_names(entry)

            function = parser._get_comment(entry, "FUNCTION")
            catalytic_activity = parser._get_comment(entry, "CATALYTIC ACTIVITY")
            ec_numbers = parser._get_ec_numbers(entry)

            subcell = parser._get_subcellular_location(entry)
            keywords = parser._get_keywords(entry)
            organism_name, org_tax_id = parser._get_organism(entry)

            go_data = parser._get_go_terms(entry)
            xrefs = parser._get_cross_references(entry)

            length = parser._get_protein_length(entry)
            score = parser._get_annotation_score(entry)
            reviewed = parser.is_reviewed(entry)

            entry_type = entry.get("entryType", "")
            protein_existence = entry.get("proteinExistence", "")

            # Store GO terms
            go_bp: List[str] = []
            go_mf: List[str] = []
            go_cc: List[str] = []
            for aspect, go_list in go_data.items():
                for gt in go_list:
                    self.go_terms[gt.id] = gt
                    if aspect == "BP":
                        go_bp.append(gt.id)
                    elif aspect == "MF":
                        go_mf.append(gt.id)
                    elif aspect == "CC":
                        go_cc.append(gt.id)

            return {
                "UniProtID": primary_acc,
                "GeneName": gene_name,
                "ProteinName": protein_name,
                "AlternativeNames": alt_names,
                "Function": function,
                "CatalyticActivity": catalytic_activity,
                "EC_number": ec_numbers,
                "ProteinLength": length,
                "AnnotationScore": score,
                "Reviewed": reviewed,
                "EntryType": entry_type,
                "ProteinExistence": protein_existence,
                "Organism": organism_name,
                "TaxonomyID": org_tax_id,
                "SubcellularLocation": subcell,
                "Keywords": "; ".join(keywords),
                "KeywordCount": len(keywords),
                "GO_BP": "; ".join(go_bp),
                "GO_MF": "; ".join(go_mf),
                "GO_CC": "; ".join(go_cc),
                "GO_BP_count": len(go_bp),
                "GO_MF_count": len(go_mf),
                "GO_CC_count": len(go_cc),
                "GO_total": len(go_bp) + len(go_mf) + len(go_cc),
                "InterPro": "; ".join(xrefs.get("InterPro", [])),
                "Pfam": "; ".join(xrefs.get("Pfam", [])),
                "KEGG": "; ".join(xrefs.get("KEGG", [])),
                "RefSeq": "; ".join(xrefs.get("RefSeq", [])),
                "PDB": "; ".join(xrefs.get("PDB", [])),
                "AlphaFold": "; ".join(xrefs.get("AlphaFoldDB", [])),
            }
        except (KeyError, TypeError, ValueError, AttributeError) as e:
            logger.warning("Failed to parse entry %s: %s", entry.get("primaryAccession", ""), e)
            return None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def fetch_annotations(self, wp_accessions: List[str]) -> pd.DataFrame:
        """
        Fetch functional annotations for a list of RefSeq WP_ accessions.

        Returns a DataFrame with one row per protein and columns for
        UniProt ID, gene names, function, GO terms, EC numbers, etc.
        """
        start = time.time()
        normalized_input = self._normalize_wp_list(wp_accessions)

        # Step 1: Map RefSeq -> UniProt (also caches entries internally)
        uniprot_map = self.map_refseq_to_uniprot(wp_accessions)

        if not uniprot_map:
            logger.warning("No RefSeq -> UniProt mappings found. Returning empty DataFrame.")
            return pd.DataFrame()

        # Step 2: Load entries from cache (saved by map_refseq_to_uniprot)
        entry_cache_key = f"entries_{self.organism_id}"
        entries = self.cache.load(entry_cache_key) or []

        if not entries:
            logger.warning("No cached entries found for annotation parsing.")
            return pd.DataFrame()

        # Step 3: Build {uniprot_id -> entry} lookup
        id_to_entry = {e.get("primaryAccession", ""): e for e in entries}
        unique_mapped = set(uniprot_map.values())

        # Step 4: Parse only the mapped entries
        parsed: List[dict] = []
        for uid in unique_mapped:
            entry = id_to_entry.get(uid)
            if entry:
                rec = self._parse_entry(entry)
                if rec:
                    parsed.append(rec)

        # Step 5: Build result DataFrame
        result = pd.DataFrame(parsed) if parsed else pd.DataFrame()

        # Step 6: Map UniProt IDs back to RefSeq WP accessions
        if not result.empty:
            wp_reverse: Dict[str, str] = {}
            for wp, up in uniprot_map.items():
                if up not in wp_reverse:
                    wp_reverse[up] = wp
            result["ProteinID"] = result["UniProtID"].map(wp_reverse)

        # Step 7: Validation
        result = self._validate(result, uniprot_map, normalized_input)

        elapsed = time.time() - start

        # Step 8: Statistics
        self.annotation_stats = self._compute_stats(result, uniprot_map, normalized_input, elapsed)
        self._log_stats(self.annotation_stats)

        return result

    def _validate(
        self,
        result: pd.DataFrame,
        uniprot_map: Dict[str, str],
        normalized_input: List[str],
    ) -> pd.DataFrame:
        """Validate and deduplicate the result DataFrame."""
        if result.empty:
            return result

        # Drop duplicate UniProtIDs (keep first)
        before = len(result)
        result = result.drop_duplicates(subset=["UniProtID"])
        if len(result) < before:
            logger.info("Removed %d duplicate UniProt entries", before - len(result))

        # Drop duplicate ProteinIDs
        if "ProteinID" in result.columns:
            before = len(result)
            result = result.drop_duplicates(subset=["ProteinID"])
            if len(result) < before:
                logger.info("Removed %d duplicate ProteinID mappings", before - len(result))

        # Check for missing UniProt IDs
        mapped_set = set(uniprot_map.values())
        annotated_set = set(result["UniProtID"].unique())
        missing = mapped_set - annotated_set
        if missing:
            logger.warning(
                "%d mapped UniProt IDs have no annotation data", len(missing)
            )

        return result

    def _compute_stats(
        self,
        result: pd.DataFrame,
        uniprot_map: Dict[str, str],
        normalized_input: List[str],
        elapsed: float,
    ) -> AnnotationStats:
        stats = AnnotationStats(elapsed=elapsed)
        stats.requested = len(normalized_input)
        stats.mapped = len(uniprot_map)
        stats.annotated = len(result)
        stats.missing = stats.mapped - stats.annotated

        if not result.empty:
            stats.swissprot = int(result["Reviewed"].sum())
            stats.trembl = stats.annotated - stats.swissprot
            if "GO_total" in result.columns:
                stats.avg_go = float(result["GO_total"].mean())
            if "KeywordCount" in result.columns:
                stats.avg_keywords = float(result["KeywordCount"].mean())

        return stats

    def _log_stats(self, stats: AnnotationStats):
        logger.info("=" * 50)
        logger.info("UniProt Annotation Summary")
        logger.info("=" * 50)
        logger.info("Proteins requested:    %d", stats.requested)
        logger.info("Proteins mapped:       %d", stats.mapped)
        logger.info("Proteins annotated:    %d", stats.annotated)
        logger.info("Swiss-Prot entries:    %d", stats.swissprot)
        logger.info("TrEMBL entries:        %d", stats.trembl)
        logger.info("Average GO terms:      %.1f", stats.avg_go)
        logger.info("Average keywords:      %.1f", stats.avg_keywords)
        logger.info("Missing annotations:   %d", stats.missing)
        logger.info("Elapsed time:          %.1fs", stats.elapsed)
        logger.info("=" * 50)

    def get_go_term_info(self, go_id: str) -> str:
        gt = self.go_terms.get(go_id)
        if gt:
            return f"{gt.name} [{gt.aspect}]"
        return f"{go_id} [unknown]"
