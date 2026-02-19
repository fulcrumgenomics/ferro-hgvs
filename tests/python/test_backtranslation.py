"""Tests for backtranslation functionality."""

import ferro_hgvs


class TestCodonTable:
    """Tests for CodonTable class."""

    def test_standard_codon_table(self) -> None:
        table = ferro_hgvs.CodonTable.standard()
        assert table is not None


class TestBacktranslator:
    """Tests for Backtranslator class."""

    def test_create_backtranslator(self) -> None:
        bt = ferro_hgvs.Backtranslator()
        assert bt is not None

    def test_create_with_codon_table(self) -> None:
        table = ferro_hgvs.CodonTable.standard()
        bt = ferro_hgvs.Backtranslator(table)
        assert bt is not None

    def test_standard_factory(self) -> None:
        bt = ferro_hgvs.Backtranslator.standard()
        assert bt is not None

    def test_backtranslate_substitution(self) -> None:
        bt = ferro_hgvs.Backtranslator.standard()
        changes = bt.backtranslate_substitution("Leu", "Phe")
        assert isinstance(changes, list)
        # Leu -> Phe should have possible codon changes
        assert len(changes) > 0

    def test_backtranslate_to_stop(self) -> None:
        bt = ferro_hgvs.Backtranslator.standard()
        changes = bt.backtranslate_to_stop("Glu")
        assert isinstance(changes, list)

    def test_backtranslate_stop_loss(self) -> None:
        bt = ferro_hgvs.Backtranslator.standard()
        changes = bt.backtranslate_stop_loss("Gln")
        assert isinstance(changes, list)


class TestCodonChange:
    """Tests for CodonChange class."""

    def test_codon_change_attributes(self) -> None:
        bt = ferro_hgvs.Backtranslator.standard()
        changes = bt.backtranslate_substitution("Leu", "Phe")
        assert len(changes) > 0
        change = changes[0]
        assert hasattr(change, "ref_codon")
        assert hasattr(change, "alt_codon")
        assert hasattr(change, "changed_positions")
        assert isinstance(change.ref_codon, str)
        assert isinstance(change.alt_codon, str)
        # changed_positions is bytes containing position values
        assert isinstance(change.changed_positions, bytes)

    def test_codon_change_num_changes(self) -> None:
        bt = ferro_hgvs.Backtranslator.standard()
        changes = bt.backtranslate_substitution("Leu", "Phe")
        if changes:
            change = changes[0]
            num = change.num_changes()
            assert isinstance(num, int)
            assert num >= 1

    def test_codon_change_str(self) -> None:
        bt = ferro_hgvs.Backtranslator.standard()
        changes = bt.backtranslate_substitution("Leu", "Phe")
        if changes:
            change = changes[0]
            s = str(change)
            assert isinstance(s, str)
