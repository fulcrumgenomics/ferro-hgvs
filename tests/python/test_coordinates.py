"""Tests for coordinate conversion functionality."""

import pytest

import ferro_hgvs


class TestPositionClasses:
    """Tests for ZeroBasedPos and OneBasedPos."""

    def test_zero_based_pos(self) -> None:
        pos = ferro_hgvs.ZeroBasedPos(99)
        assert pos.value == 99

    def test_one_based_pos(self) -> None:
        pos = ferro_hgvs.OneBasedPos(100)
        assert pos.value == 100

    def test_zero_to_one_conversion(self) -> None:
        zero = ferro_hgvs.ZeroBasedPos(99)
        one = zero.to_one_based()
        assert one.value == 100

    def test_one_to_zero_conversion(self) -> None:
        one = ferro_hgvs.OneBasedPos(100)
        zero = one.to_zero_based()
        assert zero.value == 99

    def test_roundtrip_conversion(self) -> None:
        original = ferro_hgvs.ZeroBasedPos(50)
        roundtrip = original.to_one_based().to_zero_based()
        assert roundtrip.value == original.value


class TestCoordinateFunctions:
    """Tests for coordinate conversion functions."""

    def test_hgvs_pos_to_index(self) -> None:
        # HGVS positions are 1-based, indices are 0-based
        index = ferro_hgvs.hgvs_pos_to_index(1)
        assert index == 0

    def test_index_to_hgvs_pos(self) -> None:
        pos = ferro_hgvs.index_to_hgvs_pos(0)
        assert pos == 1


class TestCoordinateMapper:
    """Tests for CoordinateMapper class."""

    def test_create_mapper(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        assert mapper is not None

    def test_has_transcript(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        # Test data should have NM_000088.3
        assert mapper.has_transcript("NM_000088.3")
        assert not mapper.has_transcript("NONEXISTENT.1")

    def test_c_to_n(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        result = mapper.c_to_n("NM_000088.3", 100)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_c_to_p(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        protein_pos = mapper.c_to_p("NM_000088.3", 100)
        assert isinstance(protein_pos, int)
        assert protein_pos > 0

    def test_n_to_c(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        result = mapper.n_to_c("NM_000088.3", 100)
        assert isinstance(result, tuple)
        assert len(result) == 3  # (cds_position, offset, is_utr3)

    def test_n_to_c_with_downstream(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        # Should accept downstream parameter
        result = mapper.n_to_c("NM_000088.3", 100, downstream=False)
        assert isinstance(result, tuple)

    def test_get_strand(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        strand = mapper.get_strand("NM_000088.3")
        assert strand in (ferro_hgvs.Strand.Plus, ferro_hgvs.Strand.Minus)

    def test_invalid_transcript_raises(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        with pytest.raises(ValueError, match="Transcript not found"):
            mapper.c_to_n("NONEXISTENT.1", 100)
