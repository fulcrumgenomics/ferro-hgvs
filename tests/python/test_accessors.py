"""Tests for Python variant accessor properties and methods."""

import ferro_hgvs


class TestPositionAccessors:
    """Tests for start/end/offset properties."""

    def test_genomic_substitution(self) -> None:
        v = ferro_hgvs.parse("NC_000001.11:g.12345A>G")
        assert v.start == 12345
        assert v.end == 12345
        assert v.offset is None

    def test_range_deletion(self) -> None:
        v = ferro_hgvs.parse("NC_000001.11:g.12345_12350del")
        assert v.start == 12345
        assert v.end == 12350

    def test_coding_intronic_offset(self) -> None:
        v = ferro_hgvs.parse("NM_000088.3:c.93+1G>T")
        assert v.start == 93
        assert v.offset == 1

    def test_noncoding_intronic_offset(self) -> None:
        v = ferro_hgvs.parse("NR_046018.2:n.100+5A>G")
        assert v.start == 100
        assert v.offset == 5

    def test_rna_intronic_offset(self) -> None:
        v = ferro_hgvs.parse("NM_000088.3:r.100+5a>g")
        assert v.start == 100
        assert v.offset == 5

    def test_utr3_indistinguishable_from_cds(self) -> None:
        # The `*` UTR marker is not exposed by start/end — both report the raw
        # base value. This is documented behavior; callers needing to
        # distinguish should inspect the variant string.
        utr = ferro_hgvs.parse("NM_000088.3:c.*5A>G")
        cds = ferro_hgvs.parse("NM_000088.3:c.5A>G")
        assert utr.start == cds.start == 5

    def test_protein_returns_none(self) -> None:
        v = ferro_hgvs.parse("NP_000079.2:p.Glu6Val")
        assert v.start is None
        assert v.end is None
        assert v.offset is None
        assert v.indel_length is None


class TestEditAccessors:
    """Tests for substitution_bases, is_identity, indel_length, is_frameshift."""

    def test_substitution_bases(self) -> None:
        v = ferro_hgvs.parse("NC_000001.11:g.12345A>G")
        assert v.substitution_bases == ("A", "G")

    def test_substitution_bases_none_for_deletion(self) -> None:
        v = ferro_hgvs.parse("NC_000001.11:g.12345del")
        assert v.substitution_bases is None

    def test_is_identity_true(self) -> None:
        v = ferro_hgvs.parse("NM_000088.3:c.100=")
        assert v.is_identity()

    def test_is_identity_false(self) -> None:
        v = ferro_hgvs.parse("NC_000001.11:g.12345A>G")
        assert not v.is_identity()

    def test_indel_length_substitution(self) -> None:
        assert ferro_hgvs.parse("NC_000001.11:g.12345A>G").indel_length == 0

    def test_indel_length_deletion(self) -> None:
        assert ferro_hgvs.parse("NC_000001.11:g.12345_12350del").indel_length == -6

    def test_indel_length_insertion(self) -> None:
        assert ferro_hgvs.parse("NC_000001.11:g.12345_12346insATG").indel_length == 3

    def test_indel_length_delins(self) -> None:
        # Replaces 6 bases with 4 → net -2.
        assert ferro_hgvs.parse("NC_000001.11:g.12345_12350delinsATTT").indel_length == -2

    def test_indel_length_uncertain_insertion_is_none(self) -> None:
        # ins(10_20) has unknown exact length.
        v = ferro_hgvs.parse("NC_000001.11:g.12345_12346ins(10_20)")
        assert v.indel_length is None

    def test_is_frameshift_true(self) -> None:
        assert ferro_hgvs.parse("NC_000001.11:g.12345del").is_frameshift()

    def test_is_frameshift_in_frame(self) -> None:
        assert not ferro_hgvs.parse("NC_000001.11:g.12345_12347del").is_frameshift()


class TestAlleles:
    """Tests for num_variants, variants(), and allele start/end behavior."""

    def test_simple_variant_num_variants(self) -> None:
        v = ferro_hgvs.parse("NC_000001.11:g.12345A>G")
        assert v.num_variants == 1
        assert len(v.variants()) == 1

    def test_allele_num_variants(self) -> None:
        v = ferro_hgvs.parse("NM_000088.3:c.[100A>G;200C>T]")
        assert v.num_variants == 2
        subs = v.variants()
        assert len(subs) == 2
        assert subs[0].start == 100
        assert subs[1].start == 200

    def test_multi_allele_start_is_none(self) -> None:
        # Multi-sub-variant alleles have ambiguous start; expose None.
        v = ferro_hgvs.parse("NM_000088.3:c.[100A>G;200C>T]")
        assert v.start is None
        assert v.end is None

    def test_single_allele_delegates(self) -> None:
        v = ferro_hgvs.parse("NM_000088.3:c.[100A>G]")
        assert v.start == 100
        assert v.end == 100


class TestToDict:
    """Tests for to_dict() including new accessor keys."""

    def test_substitution_keys(self) -> None:
        d = ferro_hgvs.parse("NC_000001.11:g.12345A>G").to_dict()
        assert d["start"] == 12345
        assert d["end"] == 12345
        assert d["ref_base"] == "A"
        assert d["alt_base"] == "G"
        assert d["indel_length"] == 0
        assert d["num_variants"] == 1

    def test_deletion_omits_substitution_keys(self) -> None:
        d = ferro_hgvs.parse("NC_000001.11:g.12345_12350del").to_dict()
        assert d["start"] == 12345
        assert d["end"] == 12350
        assert d["indel_length"] == -6
        assert "ref_base" not in d
        assert "alt_base" not in d

    def test_intronic_includes_offset(self) -> None:
        d = ferro_hgvs.parse("NM_000088.3:c.93+1G>T").to_dict()
        assert d["offset"] == 1

    def test_protein_omits_position_keys(self) -> None:
        d = ferro_hgvs.parse("NP_000079.2:p.Glu6Val").to_dict()
        assert "start" not in d
        assert "end" not in d
        assert "offset" not in d
        assert "indel_length" not in d
