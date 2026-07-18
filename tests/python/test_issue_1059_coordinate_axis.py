"""Issue #1059: expose an allele's coordinate axis on the top-level variant.

Before this change, the axis predicates (``is_coding``/``is_genomic``/...) returned
``False`` for *every* allele parent, so an allele's coordinate system could only be
discovered by iterating ``variants()``. This adds a first-class ``axis`` property (an
``Axis`` enum, or ``None`` when there is no single axis) and deprecates the boolean
predicates in favour of it.
"""

import warnings

import pytest

import ferro_hgvs
from ferro_hgvs import Axis


class TestAxisProperty:
    """The ``axis`` property resolves the coordinate system of any variant."""

    @pytest.mark.parametrize(
        ("desc", "expected"),
        [
            ("NM_000088.3:c.100A>G", Axis.Coding),
            ("NC_000001.11:g.12345A>G", Axis.Genomic),
            ("NR_046018.2:n.100A>G", Axis.NonCoding),
            ("NM_000088.3:r.100a>g", Axis.Rna),
            ("NP_000079.2:p.Glu6Val", Axis.Protein),
            ("NC_012920.1:m.100A>G", Axis.Mitochondrial),
        ],
    )
    def test_leaf_variants_report_their_axis(self, desc: str, expected: Axis) -> None:
        assert ferro_hgvs.parse(desc).axis == expected

    @pytest.mark.parametrize(
        ("desc", "expected"),
        [
            # The reported bug: an allele parent must report its members' shared axis.
            ("NM_000000.1:c.[6G>C;16_18del]", Axis.Coding),
            ("NC_000000.1:g.[6G>C;16_18del]", Axis.Genomic),
            ("NP_000000.1:p.[(Arg8Gln);(Ser10Gly)]", Axis.Protein),
        ],
    )
    def test_allele_parent_reports_shared_axis(self, desc: str, expected: Axis) -> None:
        variant = ferro_hgvs.parse(desc)
        assert variant.num_variants > 1  # this is genuinely an allele
        assert variant.axis == expected

    def test_no_single_axis_returns_none(self) -> None:
        # An RNA fusion spans two partner frames -> no single coordinate axis.
        fusion = ferro_hgvs.parse("NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924")
        assert fusion.axis is None

    def test_genome_ring_is_genomic(self) -> None:
        # A genome ring (`ACC:g.[seg1::seg2]`) is a single genomic accession whose
        # `::` are intra-molecule ISCN2020 break junctions, not a two-frame fusion,
        # so it resolves to Genomic (unlike the RNA fusion above).
        ring = ferro_hgvs.parse("NC_000022.11:g.pter_1000del::2000_qterdel")
        assert ring.axis == Axis.Genomic

    def test_single_and_allele_are_consistent(self) -> None:
        # The inconsistency from the issue: wrapping an edit in an allele must not
        # change the reported axis.
        single = ferro_hgvs.parse("NM_000000.1:c.6G>C")
        allele = ferro_hgvs.parse("NM_000000.1:c.[6G>C;16_18del]")
        assert single.axis == allele.axis == Axis.Coding


class TestAxisEnum:
    """The ``Axis`` enum carries its HGVS code and molecule-type groupings."""

    @pytest.mark.parametrize(
        ("axis", "code"),
        [
            (Axis.Genomic, "g"),
            (Axis.Coding, "c"),
            (Axis.NonCoding, "n"),
            (Axis.Rna, "r"),
            (Axis.Protein, "p"),
            (Axis.Mitochondrial, "m"),
            (Axis.Circular, "o"),
        ],
    )
    def test_code(self, axis: Axis, code: str) -> None:
        assert axis.code() == code
        assert str(axis) == code

    def test_dna_grouping_includes_mito_and_circular(self) -> None:
        for axis in (Axis.Genomic, Axis.Coding, Axis.NonCoding, Axis.Mitochondrial, Axis.Circular):
            assert axis.is_dna()
            assert not axis.is_rna()
            assert not axis.is_protein()

    def test_rna_and_protein_groupings(self) -> None:
        assert Axis.Rna.is_rna()
        assert not Axis.Rna.is_dna()
        assert Axis.Protein.is_protein()
        assert not Axis.Protein.is_dna()

    def test_reporter_is_dna_usecase(self) -> None:
        # The workaround in the issue collapses to one expression, and now works
        # for coding alleles (which previously returned all-False).
        def is_dna(v: "ferro_hgvs.HgvsVariant") -> bool:
            return v.axis is not None and v.axis.is_dna()

        assert is_dna(ferro_hgvs.parse("NM_000000.1:c.[6G>C;16_18del]"))
        assert is_dna(ferro_hgvs.parse("NC_012920.1:m.100A>G"))
        assert not is_dna(ferro_hgvs.parse("NP_000000.1:p.[(Arg8Gln);(Ser10Gly)]"))


class TestLegacyBooleanDeprecation:
    """The boolean axis predicates still work (now allele-aware) but are deprecated."""

    def test_allele_predicate_now_matches_single_edit(self) -> None:
        allele = ferro_hgvs.parse("NM_000000.1:c.[6G>C;16_18del]")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            assert allele.is_coding() is True
            assert allele.is_genomic() is False

    @pytest.mark.parametrize(
        "predicate",
        ["is_genomic", "is_coding", "is_noncoding", "is_protein", "is_rna", "is_mitochondrial"],
    )
    def test_predicates_emit_deprecation_warning(self, predicate: str) -> None:
        variant = ferro_hgvs.parse("NM_000088.3:c.100A>G")
        with pytest.warns(DeprecationWarning, match="axis"):
            getattr(variant, predicate)()
