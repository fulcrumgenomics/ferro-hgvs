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


class TestNormalizeMergeConsecutive:
    """Smoke tests for HGVS-spec consecutive-edit merging via the Python binding."""

    def test_consecutive_subs_collapse_to_delins(self) -> None:
        # Issue #72: two adjacent SNVs in cis must collapse to one delins.
        result = ferro_hgvs.normalize("NC_000001.11:g.[1000G>A;1001A>C]")
        assert result == "NC_000001.11:g.1000_1001delinsAC"

    def test_consecutive_dels_collapse_to_ranged_del(self) -> None:
        result = ferro_hgvs.normalize("NC_000001.11:g.[1000del;1001del]")
        assert result == "NC_000001.11:g.1000_1001del"

    def test_non_adjacent_subs_stay_separate(self) -> None:
        # One unchanged nt between variants -> spec keeps them separate.
        result = ferro_hgvs.normalize("NC_000001.11:g.[100G>A;102C>T]")
        assert "100G>A" in result
        assert "102C>T" in result
        assert ";" in result


class TestNormalizeEmptyInsertDelinsToDel:
    """Issue #81 item A3: a delins with an empty inserted sequence -> del."""

    def test_single_position_empty_delins(self) -> None:
        # NM_000088.3 c.10 = G; c.10delins is semantically a deletion of G.
        result = ferro_hgvs.normalize("NM_000088.3:c.10delins")
        assert result == "NM_000088.3:c.10del"

    def test_multi_position_empty_delins_with_3prime_shift(self) -> None:
        # NM_000088.3 c.10_11 = GT; rewriting to del then applying the spec's
        # 3'-rule shifts to c.11_12 because ref[10]=G == ref[12]=G.
        result = ferro_hgvs.normalize("NM_000088.3:c.10_11delins")
        assert result == "NM_000088.3:c.11_12del"
