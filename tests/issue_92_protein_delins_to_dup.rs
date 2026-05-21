//! Issue #92 — `p.delins` → canonical-form rewrites via affix-trim.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/92>.
//!
//! HGVS spec mandates that a `delins` whose inserted residues share a
//! prefix or suffix with the deleted residues must be re-described
//! with the trimmed range and trimmed insert (`delins.md:53-56`). After
//! affix-trim, the residual is canonicalized per Prioritization
//! (`general.md:56-57`): `sub > del > dup > ins > delins`.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// Synthetic protein `MKVAAALELE` (10 AAs): poly-A at p.4..p.6,
/// LELE motif at p.7..p.10. Matches the fixture in
/// `tests/issue_92_protein_ins_to_dup.rs`.
fn provider_with_polya_protein() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein("NP_TESTPROT.1", "MKVAAALELE");
    provider
}

/// Affix-trim collapses a delins whose inserted sequence equals the
/// deleted sequence to a silent identity. `p.Leu7_Glu8delinsLeuGlu`
/// over ref `LE` should normalize to `p.Leu7_Glu8=`.
#[test]
fn protein_delins_collapses_to_identity_when_insert_equals_delete() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Leu7_Glu8delinsLeuGlu").expect("parse p.Leu7_Glu8delinsLeuGlu");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Leu7_Glu8delinsLeuGlu");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Leu7_Glu8=") || out.ends_with(":p.="),
        "delins where insert == delete must collapse to identity; got {out:?}",
    );
    assert!(
        !out.contains("delins"),
        "delins must NOT survive the affix-trim collapse; got {out:?}",
    );
}
