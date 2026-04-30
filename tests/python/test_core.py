"""Tests for core ferro-hgvs functionality."""

import json

import pytest

import ferro_hgvs


class TestParsing:
    """Tests for HGVS parsing."""

    def test_parse_coding_substitution(self) -> None:
        variant = ferro_hgvs.parse("NM_000088.3:c.100A>G")
        assert variant.reference == "NM_000088.3"
        assert variant.variant_type == "coding"
        assert str(variant) == "NM_000088.3:c.100A>G"

    def test_parse_genomic_substitution(self) -> None:
        variant = ferro_hgvs.parse("NC_000001.11:g.12345A>G")
        assert variant.reference == "NC_000001.11"
        assert variant.variant_type == "genomic"

    def test_parse_protein_substitution(self) -> None:
        variant = ferro_hgvs.parse("NP_000079.2:p.Glu6Val")
        assert variant.reference == "NP_000079.2"
        assert variant.variant_type == "protein"

    def test_parse_noncoding_substitution(self) -> None:
        variant = ferro_hgvs.parse("NR_046018.2:n.100A>G")
        assert variant.reference == "NR_046018.2"
        assert variant.variant_type == "non_coding"

    def test_parse_invalid_raises_error(self) -> None:
        with pytest.raises(ValueError, match="Parse error"):
            ferro_hgvs.parse("invalid")

    def test_parse_deletion(self) -> None:
        variant = ferro_hgvs.parse("NM_000088.3:c.100del")
        assert "del" in str(variant)

    def test_parse_insertion(self) -> None:
        variant = ferro_hgvs.parse("NM_000088.3:c.100_101insATG")
        assert "ins" in str(variant)

    def test_parse_duplication(self) -> None:
        variant = ferro_hgvs.parse("NM_000088.3:c.100dup")
        assert "dup" in str(variant)

    @pytest.mark.parametrize(
        ("hgvs", "selector"),
        [
            ("MYSEQ(1):c.100A>G", "1"),
            ("MY-SEQ(GENE1):c.100A>G", "GENE1"),
            ("MYREF_SEQ(1):c.100A>G", "1"),
            ("MYREF_SEQ(1):p.(Arg8Gln)", "1"),
        ],
    )
    def test_parse_accepts_gene_selector_on_non_refseq(self, hgvs: str, selector: str) -> None:
        # Issue #69: gene selectors must parse on any valid accession AND be
        # captured on the variant — not silently dropped.
        variant = ferro_hgvs.parse(hgvs)
        # to_json() wraps the variant body under its discriminator key (e.g. "Cds").
        body = next(iter(json.loads(variant.to_json()).values()))
        assert body["gene_symbol"] == selector


class TestHgvsVariant:
    """Tests for HgvsVariant class."""

    def test_variant_equality(self) -> None:
        v1 = ferro_hgvs.parse("NM_000088.3:c.100A>G")
        v2 = ferro_hgvs.parse("NM_000088.3:c.100A>G")
        assert v1 == v2

    def test_variant_hash(self) -> None:
        v1 = ferro_hgvs.parse("NM_000088.3:c.100A>G")
        v2 = ferro_hgvs.parse("NM_000088.3:c.100A>G")
        assert hash(v1) == hash(v2)
        # Can be used in sets/dicts
        variant_set = {v1, v2}
        assert len(variant_set) == 1

    def test_variant_repr(self) -> None:
        variant = ferro_hgvs.parse("NM_000088.3:c.100A>G")
        repr_str = repr(variant)
        assert "HgvsVariant" in repr_str
        assert "NM_000088.3" in repr_str

    def test_to_dict(self) -> None:
        variant = ferro_hgvs.parse("NM_000088.3:c.100A>G")
        d = variant.to_dict()
        assert isinstance(d, dict)
        assert "reference" in d
        assert d["reference"] == "NM_000088.3"
        assert d["variant_type"] == "coding"


class TestVersion:
    """Tests for version information."""

    def test_version_exists(self) -> None:
        assert hasattr(ferro_hgvs, "__version__")
        assert isinstance(ferro_hgvs.__version__, str)
        assert len(ferro_hgvs.__version__) > 0
