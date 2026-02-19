"""Tests for rsID functionality."""

import pytest

import ferro_hgvs


class TestRsIdParsing:
    """Tests for rsID parsing functions."""

    def test_parse_rsid_with_prefix(self) -> None:
        value = ferro_hgvs.parse_rsid_value("rs121913529")
        assert value == 121913529

    def test_parse_rsid_without_prefix(self) -> None:
        value = ferro_hgvs.parse_rsid_value("121913529")
        assert value == 121913529

    def test_parse_rsid_invalid(self) -> None:
        with pytest.raises(ValueError, match="Invalid rsID"):
            ferro_hgvs.parse_rsid_value("invalid")

    def test_format_rsid(self) -> None:
        formatted = ferro_hgvs.format_rsid_value(121913529)
        assert formatted == "rs121913529"


class TestInMemoryRsIdLookup:
    """Tests for InMemoryRsIdLookup class."""

    def test_create_with_test_data(self) -> None:
        lookup = ferro_hgvs.InMemoryRsIdLookup.with_test_data()
        assert lookup is not None

    def test_contains(self) -> None:
        lookup = ferro_hgvs.InMemoryRsIdLookup.with_test_data()
        # Test data contains rs121913529 (BRAF V600E) and rs80357906 (BRCA1)
        assert lookup.contains("rs121913529")

    def test_lookup_existing(self) -> None:
        lookup = ferro_hgvs.InMemoryRsIdLookup.with_test_data()
        results = lookup.lookup("rs121913529")
        assert isinstance(results, list)
        assert len(results) > 0

    def test_lookup_nonexistent(self) -> None:
        lookup = ferro_hgvs.InMemoryRsIdLookup.with_test_data()
        with pytest.raises(ValueError, match="Lookup failed"):
            lookup.lookup("rs999999999999")


class TestRsIdResult:
    """Tests for RsIdResult class."""

    def test_result_attributes(self) -> None:
        lookup = ferro_hgvs.InMemoryRsIdLookup.with_test_data()
        results = lookup.lookup("rs121913529")
        assert len(results) > 0
        result = results[0]
        assert hasattr(result, "rsid")
        assert hasattr(result, "contig")
        assert hasattr(result, "position")
        assert hasattr(result, "reference")
        assert hasattr(result, "alternate")

    def test_result_is_snv(self) -> None:
        lookup = ferro_hgvs.InMemoryRsIdLookup.with_test_data()
        results = lookup.lookup("rs121913529")
        assert len(results) > 0
        result = results[0]
        # BRAF V600E is a SNV
        assert result.is_snv()

    def test_result_to_hgvs(self) -> None:
        lookup = ferro_hgvs.InMemoryRsIdLookup.with_test_data()
        results = lookup.lookup("rs121913529")
        assert len(results) > 0
        result = results[0]
        hgvs = result.to_hgvs()
        assert isinstance(hgvs, str)
        assert "g." in hgvs
