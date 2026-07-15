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


class TestNormalizeIssue160RevcompInvSubspans:
    """Issue #160 (as corrected by issue #1034): revcomp inv detection.

    Per HGVS an inversion describes a *maximal contiguous run* whose alt is the
    reverse complement of the reference (`DNA/inversion.md`, `general.md:56`).
    A delins whose WHOLE contiguous run is a revcomp becomes an `inv`; a revcomp
    SUB-run of a longer contiguous change is NOT carved out — that change stays
    a single `delins` (issue #1034 regression fix).

    Uses a temporary JSON reference so the test can pin specific genomic bases
    at the positions referenced by the variant strings (the bundled
    `with_test_data` mock has no NC_000001.11 sequence).
    """

    @staticmethod
    def _make_reference_json(tmp_path, contig: str, start_1based: int, bases: str) -> str:
        """Write a JSON reference with `bases` at `start_1based..` and return its path."""
        # Filler `A` from position 1 up to start_1based, then the test bases,
        # then more filler to cover the normalize window. JSON loader accepts
        # an object with `genomic_sequences`.
        seq = ["A"] * max(2000, start_1based + len(bases) + 200)
        for i, b in enumerate(bases):
            seq[start_1based - 1 + i] = b
        payload = {
            "transcripts": [],
            "proteins": {},
            "genomic_sequences": {contig: "".join(seq)},
        }
        path = tmp_path / "ref.json"
        path.write_text(json.dumps(payload))
        return str(path)

    def test_full_span_revcomp_in_cis_merges_to_inv(self, tmp_path) -> None:
        # Cis allele over GG at 1092-1093 → full-span inv.
        ref = self._make_reference_json(tmp_path, "NC_000001.11", 1092, "GG")
        normalizer = ferro_hgvs.Normalizer(reference_json=ref)
        result = normalizer.normalize("NC_000001.11:g.[1092G>C;1093G>C]")
        assert result == "NC_000001.11:g.1092_1093inv"

    def test_sub_span_revcomp_stays_single_delins(self, tmp_path) -> None:
        # Headline issue #160 case, spec-corrected by issue #1034: cis allele
        # over TCC at 1150-1152 merges to delinsGAG. The 2-nt sub-run 1150_1151
        # (TC->GA) is a revcomp, but the whole contiguous run TCC->GAG is not
        # (revcomp(TCC)=GGA != GAG), so it must NOT split into
        # [1150_1151inv;1152C>G] — it stays a single delins.
        ref = self._make_reference_json(tmp_path, "NC_000001.11", 1150, "TCC")
        normalizer = ferro_hgvs.Normalizer(reference_json=ref)
        result = normalizer.normalize("NC_000001.11:g.[1150T>G;1151C>A;1152C>G]")
        assert result == "NC_000001.11:g.1150_1152delinsGAG"

    def test_user_typed_delins_with_revcomp_subspan_stays_delins_symmetric(self, tmp_path) -> None:
        # User-typed `g.1150_1152delinsGAG` produces the same canonical
        # form as the cis-allele input — the rule depends only on
        # (ref, position, alt), not input shape (issue #1034).
        ref = self._make_reference_json(tmp_path, "NC_000001.11", 1150, "TCC")
        normalizer = ferro_hgvs.Normalizer(reference_json=ref)
        result = normalizer.normalize("NC_000001.11:g.1150_1152delinsGAG")
        assert result == "NC_000001.11:g.1150_1152delinsGAG"

    def test_cds_full_span_revcomp_in_cis_merges_to_inv(self) -> None:
        # NM_000088.3 c.9-10 = "GG" (in the bundled mock); revcomp = "CC".
        # Cis allele merges to delinsCC then canonicalizes to inv.
        result = ferro_hgvs.normalize("NM_000088.3:c.[9G>C;10G>C]")
        assert result == "NM_000088.3:c.9_10inv"

    def test_cds_sub_span_revcomp_stays_single_delins(self) -> None:
        # NM_000088.3 c.13-15 = "CTG"; the 2-nt sub-run c.13-14 (CT->AG) is a
        # revcomp, but the whole contiguous run CTG->AGT is not
        # (revcomp(CTG)=CAG != AGT), so no sub-run is carved out (issue #1034):
        # it stays a single delins, not [13_14inv;15G>T].
        result = ferro_hgvs.normalize("NM_000088.3:c.[13C>A;14T>G;15G>T]")
        assert result == "NM_000088.3:c.13_15delinsAGT"

    def test_rna_full_span_revcomp_in_cis_merges_to_inv(self) -> None:
        # NM_000088.3 r.9-10 over "GG" → revcomp "CC" → r.9_10inv.
        # Lower-case bases per HGVS r. spec.
        result = ferro_hgvs.normalize("NM_000088.3:r.[9g>c;10g>c]")
        assert result == "NM_000088.3:r.9_10inv"

    def test_n_user_typed_delins_with_full_span_revcomp(self) -> None:
        # Coordinate-system parity for n.: the n. allele syntax `[...]` isn't
        # accepted by the parser, so cover the n. inv-split path via a
        # user-typed delins. NM_000088.3 n.9-10 = "GG" (transcript bytes);
        # revcomp = "CC" → n.9_10inv.
        result = ferro_hgvs.normalize("NM_000088.3:n.9_10delinsCC")
        assert result == "NM_000088.3:n.9_10inv"


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


class TestNormalizeDelinsSubOnlyDecompose:
    """Issue #81 item A10 (sub-only branch, tracked in #165): a delins whose
    span contains two or more independent single-base mismatches separated
    by at least one unchanged nucleotide decomposes to the individual
    `[X>A; B>Y]` form per `general.md:56` (sub > delins). The codon-frame
    exception (`general.md:35-38`) preserves an embedded
    `[Sub; Identity; Sub]` triplet whose endpoints share a codon.

    Reference: NM_000088.3 has CDS sequence
    `ATGCCCAAGGTGCTGCCCCAGATGCTGCCAGTGCTGCTGCTGCTGCTGCTGCTGCTGCTG`,
    so c.10 = G, c.11 = T, c.12 = G, c.13 = C, c.14 = T. Codon 4 spans
    c.10..c.12 and codon 5 spans c.13..c.15."""

    def test_coding_delins_codon_frame_pair_preserved(self) -> None:
        # ref c.10..c.12 = GTG, alt = ATC → [Sub(G>A); Identity(T);
        # Sub(G>C)]. Endpoints 10 and 12 share codon 4, so the spec
        # codon-frame exception keeps the form as a 3-base delins.
        result = ferro_hgvs.normalize("NM_000088.3:c.10_12delinsATC")
        assert result == "NM_000088.3:c.10_12delinsATC"

    def test_coding_cross_codon_pair_decomposes(self) -> None:
        # ref c.12..c.14 = GCT, alt = ACA → [Sub(G>A); Identity(C);
        # Sub(T>A)]. Endpoints 12 and 14 sit in codons 4 and 5
        # respectively, so the codon-frame exception does NOT apply and
        # the spec's sub > delins priority splits the pair.
        result = ferro_hgvs.normalize("NM_000088.3:c.12_14delinsACA")
        assert "12G>A" in result
        assert "14T>A" in result
        # No 3-base delins should remain.
        assert "12_14delins" not in result

    def test_coding_embedded_codon_frame_triplet_within_longer_span(
        self,
    ) -> None:
        # ref c.10..c.13 = GTGC, alt = ATCA → [Sub(G>A); Identity(T);
        # Sub(G>C); Sub(C>A)]. The (10, 12) pair is in codon 4 and
        # preserves as a 3-base delins; c.13 sits in codon 5, so the
        # spec exception "two variants together affecting one amino
        # acid" does not extend to it and it emits as a separate sub.
        result = ferro_hgvs.normalize("NM_000088.3:c.10_13delinsATCA")
        assert "10_12delinsATC" in result
        assert "13C>A" in result


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


class TestNormalizeDelinsToInversion:
    """Issue #81 item A2: delins whose insert is the revcomp of the deleted ref -> inv."""

    def test_simple_delins_to_inv(self) -> None:
        # NM_000088.3 c.7_9 = AAG; revcomp(AAG) = CTT; outer A/G not complement
        # so no shortening — result is c.7_9inv.
        result = ferro_hgvs.normalize("NM_000088.3:c.7_9delinsCTT")
        assert result == "NM_000088.3:c.7_9inv"

    def test_delins_to_inv_with_outer_shortening(self) -> None:
        # NM_000088.3 c.3_6 = GCCC; revcomp(GCCC) = GGGC; outer G/C are
        # complement so the inversion shortens to the inner CC at c.4_5.
        result = ferro_hgvs.normalize("NM_000088.3:c.3_6delinsGGGC")
        assert result == "NM_000088.3:c.4_5inv"

    def test_delins_one_base_complement_is_substitution(self) -> None:
        # NM_000088.3 c.1 = A; insert T. A 1-base complement is a substitution
        # per HGVS spec, NEVER an inversion.
        result = ferro_hgvs.normalize("NM_000088.3:c.1delinsT")
        assert result == "NM_000088.3:c.1A>T"


class TestNormalizeDelinsSharedAffixTrimming:
    """Issue #81 item A2 (extension): shared-affix trimming reduces a delins
    to its minimal HGVS form (sub / del / ins / smaller delins) per the
    sub > del > inv > dup > ins priority."""

    def test_delins_shared_suffix_becomes_substitution(self) -> None:
        # NM_000088.3 c.1_4 = ATGC; insert AAGC shares the AGC suffix and
        # the leading A, leaving T -> A at c.2.
        result = ferro_hgvs.normalize("NM_000088.3:c.1_4delinsAAGC")
        assert result == "NM_000088.3:c.2T>A"

    def test_delins_shared_affix_pure_deletion(self) -> None:
        # NM_000088.3 c.1_4 = ATGC; insert AC shares prefix A and suffix C,
        # leaving TG deleted at c.2_3.
        result = ferro_hgvs.normalize("NM_000088.3:c.1_4delinsAC")
        assert result == "NM_000088.3:c.2_3del"

    def test_delins_shared_affix_pure_insertion(self) -> None:
        # NM_000088.3 c.1_3 = ATG; insert ATCG shares prefix AT and suffix G,
        # leaving C inserted between c.2 and c.3.
        result = ferro_hgvs.normalize("NM_000088.3:c.1_3delinsATCG")
        assert result == "NM_000088.3:c.2_3insC"
