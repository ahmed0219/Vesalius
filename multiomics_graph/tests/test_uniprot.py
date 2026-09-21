"""
Tests for the UniProt annotation module.
"""

from pathlib import Path
from typing import List
from unittest.mock import MagicMock

import pandas as pd
import pytest

from annotation.uniprot import (
    GoTerm,
    UniProtAnnotator,
    _AnnotationParser,
    _CacheManager,
)


# ------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------

@pytest.fixture
def annotator(tmp_path: Path) -> UniProtAnnotator:
    return UniProtAnnotator(organism_id="562", cache_dir=tmp_path / "cache")


@pytest.fixture
def sample_wp_ids() -> List[str]:
    return ["WP_000673464.1", "WP_000060112.1", "WP_000072072"]


@pytest.fixture
def sample_wp_ids_normalized() -> List[str]:
    return ["WP_000673464", "WP_000060112", "WP_000072072"]


@pytest.fixture
def sample_mapping_response() -> dict:
    return {
        "results": [
            {"from": "WP_000673464", "to": "P0A989"},
            {"from": "WP_000060112", "to": "P0A9B8"},
            {"from": "WP_000072072", "to": "P0A9E3"},
        ]
    }


@pytest.fixture
def sample_uniprot_entry() -> dict:
    return {
        "primaryAccession": "P0A989",
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "proteinExistence": "Experimental evidence at protein level",
        "annotationScore": 5.0,
        "organism": {
            "scientificName": "Escherichia coli",
            "lineages": [],
            "dbReferences": [{"type": "NCBI Taxonomy", "id": "83333"}],
        },
        "sequence": {"length": 1101, "molWeight": 123456},
        "genes": [{"geneName": {"value": "thrL"}}],
        "proteinDescription": {
            "recommendedName": {
                "fullName": {"value": "Threonine synthase"}
            },
            "alternativeNames": [{"fullName": {"value": "TS"}}],
        },
        "comments": [
            {
                "commentType": "FUNCTION",
                "texts": [{"value": "Catalyzes threonine biosynthesis"}],
            },
            {
                "commentType": "CATALYTIC ACTIVITY",
                "reaction": {"ecNumber": "4.2.3.1"},
                "texts": [{"value": "Catalytic activity description"}],
            },
            {
                "commentType": "SUBCELLULAR LOCATION",
                "subcellularLocations": [
                    {"location": {"value": "Cytoplasm"}}
                ],
            },
        ],
        "keywords": [
            {"name": "Amino-acid biosynthesis"},
            {"name": "Threonine biosynthesis"},
        ],
        "go_f": [{"goId": "GO:0003824", "name": "catalytic activity"}],
        "go_p": [
            {"goId": "GO:0009088", "name": "threonine biosynthetic process"}
        ],
        "go_c": [{"goId": "GO:0005737", "name": "cytoplasm"}],
        "uniProtKBCrossReferences": [
            {"database": "EC", "id": "4.2.3.1"},
            {"database": "InterPro", "id": "IPR001000"},
            {"database": "Pfam", "id": "PF00291"},
            {"database": "KEGG", "id": "eco:b0003"},
            {"database": "RefSeq", "id": "WP_000673464"},
        ],
    }


# ------------------------------------------------------------------
# Identifier normalization
# ------------------------------------------------------------------

class TestIdentifierNormalization:
    def test_strips_version_suffix(self):
        assert UniProtAnnotator._normalize_wp("WP_000673464.1") == "WP_000673464"

    def test_handles_unversioned(self):
        assert UniProtAnnotator._normalize_wp("WP_000673464") == "WP_000673464"

    def test_handles_empty(self):
        assert UniProtAnnotator._normalize_wp("") == ""

    def test_handles_none(self):
        assert UniProtAnnotator._normalize_wp(None) == ""  # type: ignore

    def test_normalize_list_deduplicates(self):
        ids = ["WP_000673464.1", "WP_000673464", "WP_000060112.1"]
        result = UniProtAnnotator._normalize_wp_list(ids)
        assert result == ["WP_000673464", "WP_000060112"]

    def test_normalize_list_skips_empty(self):
        ids = ["WP_000673464.1", "", "WP_000060112.1"]
        result = UniProtAnnotator._normalize_wp_list(ids)
        assert result == ["WP_000673464", "WP_000060112"]


# ------------------------------------------------------------------
# Mapping candidate ranking
# ------------------------------------------------------------------

class TestMappingCandidateRanking:
    def _make_entry(self, acc: str, tax_id: str, reviewed: bool, score: float) -> dict:
        return {
            "primaryAccession": acc,
            "entryType": "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)",
            "annotationScore": score,
            "organism": {"scientificName": "", "lineages": [], "dbReferences": [{"type": "NCBI Taxonomy", "id": tax_id}]},
            "uniProtKBCrossReferences": [],
        }

    def test_organism_match_ranked_first(self, annotator: UniProtAnnotator):
        candidates = [
            ("P0A989", self._make_entry("P0A989", "562", False, 3.0)),
            ("P0A988", self._make_entry("P0A988", "83333", False, 5.0)),
        ]
        best = annotator._resolve_candidates(candidates)
        assert best == "P0A989"

    def test_reviewed_over_unreviewed_same_organism(self, annotator: UniProtAnnotator):
        candidates = [
            ("P0A989", self._make_entry("P0A989", "562", False, 5.0)),
            ("P0A988", self._make_entry("P0A988", "562", True, 3.0)),
        ]
        best = annotator._resolve_candidates(candidates)
        assert best == "P0A988"

    def test_score_as_tiebreaker(self, annotator: UniProtAnnotator):
        candidates = [
            ("P0A989", self._make_entry("P0A989", "562", True, 3.0)),
            ("P0A988", self._make_entry("P0A988", "562", True, 5.0)),
        ]
        best = annotator._resolve_candidates(candidates)
        assert best == "P0A988"


# ------------------------------------------------------------------
# Cache manager
# ------------------------------------------------------------------

class TestCacheManager:
    def test_save_and_load(self, tmp_path: Path):
        cache = _CacheManager(tmp_path)
        cache.save("test_key", {"a": 1, "b": 2})
        loaded = cache.load("test_key")
        assert loaded == {"a": 1, "b": 2}

    def test_load_missing(self, tmp_path: Path):
        cache = _CacheManager(tmp_path)
        assert cache.load("nonexistent") is None

    def test_invalidate(self, tmp_path: Path):
        cache = _CacheManager(tmp_path)
        cache.save("test_key", {"a": 1})
        assert cache.load("test_key") is not None
        cache.invalidate("test_key")
        assert cache.load("test_key") is None

    def test_is_valid_empty_dict(self, tmp_path: Path):
        cache = _CacheManager(tmp_path)
        cache.save("empty_key", {})
        assert not cache.is_valid("empty_key")

    def test_is_valid_populated(self, tmp_path: Path):
        cache = _CacheManager(tmp_path)
        cache.save("populated", {"a": 1})
        assert cache.is_valid("populated")


# ------------------------------------------------------------------
# Annotation parser
# ------------------------------------------------------------------

class TestAnnotationParser:
    def test_get_gene_name(self, sample_uniprot_entry):
        assert _AnnotationParser._get_gene_name(sample_uniprot_entry) == "thrL"

    def test_get_gene_name_missing(self):
        assert _AnnotationParser._get_gene_name({}) == ""

    def test_get_protein_name(self, sample_uniprot_entry):
        assert _AnnotationParser._get_protein_name(sample_uniprot_entry) == "Threonine synthase"

    def test_get_protein_name_missing(self):
        assert _AnnotationParser._get_protein_name({}) == ""

    def test_get_alternative_names(self, sample_uniprot_entry):
        names = _AnnotationParser._get_alternative_names(sample_uniprot_entry)
        assert "TS" in names

    def test_get_function(self, sample_uniprot_entry):
        func = _AnnotationParser._get_comment(sample_uniprot_entry, "FUNCTION")
        assert "threonine" in func

    def test_get_ec_numbers(self, sample_uniprot_entry):
        ec = _AnnotationParser._get_ec_numbers(sample_uniprot_entry)
        assert "4.2.3.1" in ec

    def test_get_subcellular_location(self, sample_uniprot_entry):
        loc = _AnnotationParser._get_subcellular_location(sample_uniprot_entry)
        assert "Cytoplasm" in loc

    def test_get_keywords(self, sample_uniprot_entry):
        kws = _AnnotationParser._get_keywords(sample_uniprot_entry)
        assert "Amino-acid biosynthesis" in kws
        assert len(kws) == 2

    def test_get_go_terms(self, sample_uniprot_entry):
        go_data = _AnnotationParser._get_go_terms(sample_uniprot_entry)
        assert "BP" in go_data
        assert "MF" in go_data
        assert "CC" in go_data
        assert len(go_data["BP"]) == 1
        assert go_data["BP"][0].id == "GO:0009088"

    def test_get_cross_references(self, sample_uniprot_entry):
        xrefs = _AnnotationParser._get_cross_references(sample_uniprot_entry)
        assert "EC" in xrefs
        assert "InterPro" in xrefs
        assert xrefs["KEGG"] == ["eco:b0003"]

    def test_get_organism(self, sample_uniprot_entry):
        name, tax_id = _AnnotationParser._get_organism(sample_uniprot_entry)
        assert "Escherichia" in name

    def test_get_protein_length(self, sample_uniprot_entry):
        assert _AnnotationParser._get_protein_length(sample_uniprot_entry) == 1101

    def test_is_reviewed(self, sample_uniprot_entry):
        assert _AnnotationParser.is_reviewed(sample_uniprot_entry) is True

    def test_is_reviewed_false(self):
        entry = {"entryType": "UniProtKB unreviewed (TrEMBL)"}
        assert _AnnotationParser.is_reviewed(entry) is False

    def test_get_annotation_score(self, sample_uniprot_entry):
        assert _AnnotationParser._get_annotation_score(sample_uniprot_entry) == 5.0

    def test_parse_entry_full(self, annotator: UniProtAnnotator, sample_uniprot_entry: dict):
        record = annotator._parse_entry(sample_uniprot_entry)
        assert record is not None
        assert record["UniProtID"] == "P0A989"
        assert record["GeneName"] == "thrL"
        assert record["ProteinName"] == "Threonine synthase"
        assert record["EC_number"] == "4.2.3.1"
        assert record["Reviewed"] is True
        assert record["GO_total"] == 3
        assert record["InterPro"] == "IPR001000"

    def test_parse_entry_empty(self, annotator: UniProtAnnotator):
        record = annotator._parse_entry({})
        # Empty dict should produce a record with empty/default values, not crash
        assert record is not None
        assert record["UniProtID"] == ""


# ------------------------------------------------------------------
# Mapping via search API (mocked)
# ------------------------------------------------------------------

class TestSearchMapping:
    def _make_refseq_entry(self, acc: str, refseq_id: str, tax_id: str = "562",
                           reviewed: bool = True, score: float = 5.0) -> dict:
        return {
            "primaryAccession": acc,
            "entryType": "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)",
            "annotationScore": score,
            "organism": {"scientificName": "", "lineages": [],
                         "dbReferences": [{"type": "NCBI Taxonomy", "id": tax_id}]},
            "sequence": {"length": 100},
            "uniProtKBCrossReferences": [{"database": "RefSeq", "id": refseq_id}],
        }

    def test_single_match(self, annotator: UniProtAnnotator):
        annotator._search_refseq_batch = MagicMock(
            return_value=[self._make_refseq_entry("P0A989", "WP_000673464")]
        )
        result = annotator.map_refseq_to_uniprot(["WP_000673464.1"])
        assert result == {"WP_000673464": "P0A989"}

    def test_multiple_candidates_organism_match_first(self, annotator: UniProtAnnotator):
        annotator._search_refseq_batch = MagicMock(
            return_value=[
                self._make_refseq_entry("P0A989", "WP_000673464", tax_id="562", reviewed=False, score=3.0),
                self._make_refseq_entry("P0A988", "WP_000673464", tax_id="83333", reviewed=False, score=5.0),
            ]
        )
        result = annotator.map_refseq_to_uniprot(["WP_000673464.1"])
        # P0A989 matches organism 562, so ranked first
        assert result.get("WP_000673464") == "P0A989"

    def test_no_mapping_returns_empty(self, annotator: UniProtAnnotator):
        annotator._search_refseq_batch = MagicMock(return_value=[])
        result = annotator.map_refseq_to_uniprot(["WP_NONEXISTENT"])
        assert len(result) == 0


# ------------------------------------------------------------------
# Full annotation pipeline (mocked)
# ------------------------------------------------------------------

class TestAnnotationPipeline:
    def test_fetch_annotations(
        self,
        annotator: UniProtAnnotator,
        sample_uniprot_entry: dict,
    ):
        # Mock search to return sample entry (with RefSeq cross-ref to WP_000673464)
        annotator._search_refseq_batch = MagicMock(
            return_value=[sample_uniprot_entry]
        )

        df = annotator.fetch_annotations(["WP_000673464.1", "WP_000060112.1"])
        assert not df.empty
        assert "UniProtID" in df.columns
        assert "ProteinID" in df.columns
        assert "GO_BP" in df.columns
        assert "GO_MF" in df.columns
        assert "GO_CC" in df.columns
        assert "InterPro" in df.columns
        assert "Reviewed" in df.columns
        assert "Keywords" in df.columns

    def test_empty_mapping_returns_empty_df(
        self,
        annotator: UniProtAnnotator,
    ):
        annotator._search_refseq_batch = MagicMock(return_value=[])

        df = annotator.fetch_annotations(["WP_NONEXISTENT"])
        assert df.empty

    def test_go_terms_populated(
        self,
        annotator: UniProtAnnotator,
        sample_uniprot_entry: dict,
    ):
        annotator._search_refseq_batch = MagicMock(
            return_value=[sample_uniprot_entry]
        )

        annotator.fetch_annotations(["WP_000673464.1"])

        assert annotator.go_terms["GO:0009088"].name == "threonine biosynthetic process"
        assert annotator.go_terms["GO:0009088"].aspect == "BP"
        assert annotator.go_terms["GO:0003824"].aspect == "MF"

    def test_annotation_stats(
        self,
        annotator: UniProtAnnotator,
        sample_uniprot_entry: dict,
    ):
        annotator._search_refseq_batch = MagicMock(
            return_value=[sample_uniprot_entry]
        )

        annotator.fetch_annotations(["WP_000673464.1"])

        assert annotator.annotation_stats is not None
        assert annotator.annotation_stats.requested == 1
        assert annotator.annotation_stats.mapped == 1
        assert annotator.annotation_stats.annotated == 1
        assert annotator.annotation_stats.swissprot == 1
        assert annotator.annotation_stats.avg_go == 3.0
        assert annotator.annotation_stats.avg_keywords == 2.0


# ------------------------------------------------------------------
# Cache behavior
# ------------------------------------------------------------------

class TestCacheBehavior:
    def test_cache_not_overwritten_with_empty(self, tmp_path: Path):
        cache = _CacheManager(tmp_path)
        cache.save("test", {"WP_000673464": "P0A989"})
        cache.save("test", {})
        loaded = cache.load("test")
        assert loaded is None
        assert not cache.is_valid("test")

    def test_cache_invalidated_on_corrupt(self, tmp_path: Path):
        path = tmp_path / "uniprot_corrupt.json"
        path.write_text("{invalid json")
        cache = _CacheManager(tmp_path)
        assert cache.load("corrupt") is None

    def test_cache_skipped_for_empty_mapping(self, tmp_path: Path):
        cache = _CacheManager(tmp_path)
        key = "refseq2uniprot_562"
        cache.save(key, {})
        assert not cache.is_valid(key)


# ------------------------------------------------------------------
# Go term info
# ------------------------------------------------------------------

class TestGoTermInfo:
    def test_get_go_term_info(self, annotator: UniProtAnnotator):
        annotator.go_terms["GO:0009088"] = GoTerm(
            id="GO:0009088",
            name="threonine biosynthetic process",
            aspect="BP",
        )
        info = annotator.get_go_term_info("GO:0009088")
        assert "threonine biosynthetic process" in info
        assert "BP" in info

    def test_get_go_term_info_unknown(self, annotator: UniProtAnnotator):
        info = annotator.get_go_term_info("GO:0000000")
        assert "unknown" in info


# ------------------------------------------------------------------
# Edge cases
# ------------------------------------------------------------------

class TestEdgeCases:
    def test_empty_accession_list(self, annotator: UniProtAnnotator):
        df = annotator.fetch_annotations([])
        assert df.empty

    def test_all_unknown_accessions(self, annotator: UniProtAnnotator):
        annotator._search_refseq_batch = MagicMock(return_value=[])
        df = annotator.fetch_annotations(["WP_NOT_EXIST_A", "WP_NOT_EXIST_B"])
        assert isinstance(df, pd.DataFrame)

    def test_duplicate_accessions_deduped(self, annotator: UniProtAnnotator):
        normalized = annotator._normalize_wp_list(
            ["WP_000673464.1", "WP_000673464", "WP_000673464.2"]
        )
        assert len(normalized) == 1
        assert normalized[0] == "WP_000673464"

    def test_parse_entry_malformed(
        self, annotator: UniProtAnnotator, caplog: pytest.LogCaptureFixture
    ):
        malformed = {"notAnEntry": True, "primaryAccession": "P00001"}
        record = annotator._parse_entry(malformed)
        assert record is not None  # Should handle gracefully


# ------------------------------------------------------------------
# Backward compatibility
# ------------------------------------------------------------------

class TestBackwardCompatibility:
    def test_public_api_unchanged(self):
        """The public interface must remain UniProtAnnotator().fetch_annotations(ids)."""
        annotator = UniProtAnnotator()
        assert hasattr(annotator, "fetch_annotations")
        assert hasattr(annotator, "map_refseq_to_uniprot")
        assert hasattr(annotator, "get_go_term_info")

    def test_go_terms_attribute(self, annotator: UniProtAnnotator):
        """self.go_terms must remain accessible for downstream code."""
        assert hasattr(annotator, "go_terms")
        assert isinstance(annotator.go_terms, dict)
