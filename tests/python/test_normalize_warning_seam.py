"""The module-level `normalize()` had no warning-bearing sibling.

`ferro_hgvs.normalize(...)` returns a `str`, so it can only report that
normalization *happened* — never that it *repaired* something. Some of those
repairs are lossy in a way the returned string does not record, the clearest
being `MEMBERS_COALESCED_FROM_REPORTED_FORM`: an input describing two cis
members comes back as one `delins`, and the individually-reported form is gone
with no trace in the result.

`Normalizer.normalize_with_warnings(...)` (issue #395 item 5) already gave the
class-based route a channel. The free function did not, so the shortest Python
path — the one the README and the docstrings put first — stayed silent.

These tests cover the free function `normalize_with_warnings`, and pin the
property that makes it a *disclosure* rather than a behaviour change: it returns
the same normalized string `normalize` does.

They run under `pytest tests/python/` after `maturin develop --features python`.
"""

import warnings

import ferro_hgvs

# Two adjacent protein substitutions on a transcript the bundled test data
# carries. They coalesce into a single `delins`, so the input's two-member
# provenance is unrecoverable from the output — which is exactly what W5005
# reports. Also asserted in Rust by
# `coalesced_members_diagnostic::a_protein_cis_allele_carries_the_provenance_warning`.
COALESCING_INPUT = "NP_000079.2:p.[Ala100Gly;Ala101Gly]"
COALESCED_OUTPUT = "NP_000079.2:p.Ala100_Ala101delinsGlyGly"

# A clean input on the bundled `NM_000088.3` test transcript: `c.4C>G` matches
# the reference base and stays in the CDS, so nothing is repaired.
CLEAN_INPUT = "NM_000088.3:c.4C>G"


def _normalize_with_warnings(variant: str):
    """Call the free function, ignoring the one-time reduced-capability warning.

    Both free functions use the bundled genome-less test data and emit a
    `UserWarning` saying so (#1012). That is orthogonal to what is under test
    here and fires at most once per process, so filtering it keeps the
    assertions about `NormalizationWarning`s rather than about Python warnings.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return ferro_hgvs.normalize_with_warnings(variant)


def _normalize(variant: str) -> str:
    """The quiet free function, under the same UserWarning filter."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return ferro_hgvs.normalize(variant)


class TestFreeFunctionNormalizeWithWarnings:
    def test_the_coalesce_provenance_reaches_the_caller(self) -> None:
        """W5005 is reachable from the module-level entry point.

        Asserted on the warning's identity *and* its text, because a caller
        needs both: the code to branch on, the message to act on.
        """
        result = _normalize_with_warnings(COALESCING_INPUT)

        assert str(result.result) == COALESCED_OUTPUT, (
            "precondition: this input must actually coalesce"
        )
        assert [w.code for w in result.warnings] == ["MEMBERS_COALESCED_FROM_REPORTED_FORM"]
        assert result.has_warnings() is True
        message = result.warnings[0].message
        assert "input described 2 cis members" in message, message
        assert "normalized form describes 1" in message, message
        assert "not recoverable from the normalized string" in message, message

    def test_the_string_is_identical_to_the_quiet_entry_point(self) -> None:
        """The disclosure moves nothing.

        `normalize` and `normalize_with_warnings` share
        `normalize_core_checked`, so the normalized variant is the same object
        either way. Pinned by comparison rather than by argument, over both a
        repaired input and a clean one.
        """
        for variant in (COALESCING_INPUT, CLEAN_INPUT):
            quiet = _normalize(variant)
            loud = _normalize_with_warnings(variant)
            assert quiet == str(loud.result), (
                f"the warning-bearing exit moved the string for {variant}"
            )

    def test_a_clean_input_reports_nothing(self) -> None:
        """Guards against a channel that warns on everything, which discloses
        nothing useful."""
        result = _normalize_with_warnings(CLEAN_INPUT)
        assert result.warnings == []
        assert result.has_warnings() is False
        assert str(result.result) == CLEAN_INPUT

    def test_it_is_exported_from_the_package(self) -> None:
        """The function must be reachable as `ferro_hgvs.normalize_with_warnings`
        and listed in `__all__`, or the type stubs and the runtime disagree."""
        assert "normalize_with_warnings" in ferro_hgvs.__all__
        assert "NormalizeResultWithWarnings" in ferro_hgvs.__all__
        assert callable(ferro_hgvs.normalize_with_warnings)
