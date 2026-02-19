"""Tests for SPDI functionality."""

import pytest

import ferro_hgvs


class TestSpdiParsing:
    """Tests for SPDI parsing."""

    def test_parse_spdi_snv(self) -> None:
        spdi = ferro_hgvs.parse_spdi("NC_000001.11:12344:A:G")
        assert spdi.sequence == "NC_000001.11"
        assert spdi.position == 12344
        assert spdi.deletion == "A"
        assert spdi.insertion == "G"

    def test_parse_spdi_deletion(self) -> None:
        spdi = ferro_hgvs.parse_spdi("NC_000001.11:12344:ATG:")
        assert spdi.deletion == "ATG"
        assert spdi.insertion == ""

    def test_parse_spdi_insertion(self) -> None:
        spdi = ferro_hgvs.parse_spdi("NC_000001.11:12344::ATG")
        assert spdi.deletion == ""
        assert spdi.insertion == "ATG"

    def test_parse_spdi_invalid(self) -> None:
        with pytest.raises(ValueError):
            ferro_hgvs.parse_spdi("invalid")


class TestSpdiVariant:
    """Tests for SpdiVariant class."""

    def test_is_snv(self) -> None:
        spdi = ferro_hgvs.parse_spdi("NC_000001.11:12344:A:G")
        assert spdi.is_snv()

    def test_is_deletion(self) -> None:
        spdi = ferro_hgvs.parse_spdi("NC_000001.11:12344:ATG:")
        assert spdi.is_deletion()

    def test_is_insertion(self) -> None:
        spdi = ferro_hgvs.parse_spdi("NC_000001.11:12344::ATG")
        assert spdi.is_insertion()

    def test_str_representation(self) -> None:
        spdi = ferro_hgvs.parse_spdi("NC_000001.11:12344:A:G")
        assert str(spdi) == "NC_000001.11:12344:A:G"

    def test_equality(self) -> None:
        s1 = ferro_hgvs.parse_spdi("NC_000001.11:12344:A:G")
        s2 = ferro_hgvs.parse_spdi("NC_000001.11:12344:A:G")
        assert s1 == s2

    def test_hash(self) -> None:
        s1 = ferro_hgvs.parse_spdi("NC_000001.11:12344:A:G")
        s2 = ferro_hgvs.parse_spdi("NC_000001.11:12344:A:G")
        assert hash(s1) == hash(s2)

    def test_to_dict(self) -> None:
        spdi = ferro_hgvs.parse_spdi("NC_000001.11:12344:A:G")
        d = spdi.to_dict()
        assert d["sequence"] == "NC_000001.11"
        assert d["position"] == 12344
        assert d["deletion"] == "A"
        assert d["insertion"] == "G"
