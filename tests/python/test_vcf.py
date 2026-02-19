"""Tests for VCF functionality."""

import pytest

import ferro_hgvs


class TestVcfRecord:
    """Tests for VcfRecord class."""

    def test_create_vcf_record(self) -> None:
        record = ferro_hgvs.VcfRecord("chr1", 12345, "A", "G")
        assert record.chrom == "chr1"
        assert record.pos == 12345
        assert record.reference == "A"
        assert record.alternate == "G"

    def test_create_vcf_record_with_id(self) -> None:
        record = ferro_hgvs.VcfRecord("chr1", 12345, "A", "G", id="rs12345")
        assert record.id == "rs12345"

    def test_snv_factory(self) -> None:
        record = ferro_hgvs.VcfRecord.snv("chr1", 12345, "A", "G")
        assert record.chrom == "chr1"
        assert record.pos == 12345
        assert record.reference == "A"
        assert record.alternate == "G"

    def test_snv_empty_ref_raises(self) -> None:
        with pytest.raises(ValueError, match="ref_base must be a non-empty string"):
            ferro_hgvs.VcfRecord.snv("chr1", 12345, "", "G")

    def test_snv_empty_alt_raises(self) -> None:
        with pytest.raises(ValueError, match="alt_base must be a non-empty string"):
            ferro_hgvs.VcfRecord.snv("chr1", 12345, "A", "")

    def test_alternates_property(self) -> None:
        record = ferro_hgvs.VcfRecord("chr1", 12345, "A", "G")
        alts = record.alternates
        assert isinstance(alts, list)
        assert "G" in alts

    def test_repr(self) -> None:
        record = ferro_hgvs.VcfRecord("chr1", 12345, "A", "G")
        repr_str = repr(record)
        assert "VcfRecord" in repr_str
        assert "chr1" in repr_str


class TestVcfToHgvs:
    """Tests for VCF to HGVS conversion."""

    def test_vcf_to_genomic_hgvs(self) -> None:
        record = ferro_hgvs.VcfRecord("NC_000001.11", 12345, "A", "G")
        variant = ferro_hgvs.vcf_to_genomic_hgvs(record)
        assert variant.variant_type == "genomic"
        assert "12345" in str(variant)

    def test_vcf_to_genomic_hgvs_deletion(self) -> None:
        record = ferro_hgvs.VcfRecord("NC_000001.11", 12345, "ATG", "A")
        variant = ferro_hgvs.vcf_to_genomic_hgvs(record)
        assert "del" in str(variant)

    def test_vcf_to_genomic_hgvs_insertion(self) -> None:
        record = ferro_hgvs.VcfRecord("NC_000001.11", 12345, "A", "ATG")
        variant = ferro_hgvs.vcf_to_genomic_hgvs(record)
        assert "ins" in str(variant)
