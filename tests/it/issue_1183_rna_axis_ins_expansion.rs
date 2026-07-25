//! Issue #1183 — the `r.` axis must run the #333 `ins[...]` expansion.
//!
//! Expansion wrappers existed for the genome/cds/tx/mt axes but not for RNA, so
//! an `r.` input carried its `ins[...]` payload through **unexpanded and
//! unvalidated**. Two consequences, both pinned here:
//!
//!   1. A resolvable payload stayed bracketed on `r.` while `c.`/`n.` flattened
//!      it to a literal.
//!   2. An *out-of-scope* payload — which `c.` correctly refuses — was laundered
//!      through `r.`, so projecting onto `c.` emitted a description ferro's own
//!      validation then rejected.
//!
//! `NM_000088.3` in the mock has `cds_start = 1`, so `r.N` and `c.N` address the
//! same base. That makes the two axes directly comparable: the `r.` result must
//! be the `c.` result re-rendered as RNA (lowercase, `u` for `T`), and any
//! divergence is the bug. Comparing against `c.` rather than a hardcoded string
//! is deliberate — it keeps these tests honest if the shared expansion or
//! canonicalization behavior changes later.

use ferro_hgvs::error::FerroError;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

fn normalize(input: &str) -> Result<String, FerroError> {
    let normalizer = Normalizer::new(MockProvider::with_test_data());
    let variant = parse_hgvs(input)?;
    normalizer.normalize(&variant).map(|v| format!("{}", v))
}

/// `r.4_6` is `CCC` and the payload `1_3` is `ATG` — chosen so no position of the
/// replaced range matches the inserted sequence. Coincidental matches make the
/// normalizer's canonical split factor the result into an allele, which would
/// hide the inserted literal these tests are about.
const OUTER: &str = "r.4_6";

// =============================================================================
// A resolvable payload must flatten on `r.`, exactly as it does on `c.`
// =============================================================================

/// Same-reference position-range payload (`ins[A_B]`, the #333 shape).
#[test]
fn rna_axis_expands_a_same_reference_position_range_payload() {
    let out = normalize(&format!("NM_000088.3:{OUTER}delins[1_3]"))
        .expect("r.-axis same-reference payload must normalize");
    assert_eq!(
        out, "NM_000088.3:r.4_6delinsaug",
        "the payload must flatten to the literal r.1_3 == aug",
    );
    assert!(
        !out.contains('['),
        "no bracket may survive expansion; got {out}",
    );
}

/// Cross-reference payload (`ins[ACC:axis.A_B]`, the #422/#773 shape). This is
/// the exact form in the issue report.
#[test]
fn rna_axis_expands_a_cross_reference_payload() {
    let out = normalize(&format!("NM_000088.3:{OUTER}delins[NR_000123.1:n.1_4]"))
        .expect("r.-axis cross-reference payload must normalize");
    // NR_000123.1 is "ACGTACGTACGT", so n.1_4 == ACGT -> rendered acgu on r.
    // Pin the whole description like the same-reference sibling above, not just
    // the literal: `r.4_6` is `CCC`, which shares neither its first nor its last
    // base with `acgu`, so the canonical split cannot trim either endpoint and
    // the anchor must stay `4_6`. A weaker `contains` would miss an anchor shift.
    assert_eq!(
        out, "NM_000088.3:r.4_6delinsacgu",
        "the cross-reference must flatten to the literal acgu at the unshifted anchor",
    );
}

/// The `r.` result must agree with the `c.` result for the same positions
/// (`cds_start == 1`), differing only in RNA rendering. This is the property the
/// issue is really about: `r.` must not diverge from `c.`.
#[test]
fn rna_axis_expansion_agrees_with_the_c_axis() {
    for payload in ["1_3", "NR_000123.1:n.1_4"] {
        let rna = normalize(&format!("NM_000088.3:r.4_6delins[{payload}]"))
            .unwrap_or_else(|e| panic!("r. with [{payload}] must normalize, got {e:?}"));
        let cds = normalize(&format!("NM_000088.3:c.4_6delins[{payload}]"))
            .unwrap_or_else(|e| panic!("c. with [{payload}] must normalize, got {e:?}"));
        // Re-render the c. answer as RNA: lowercase, and u for T.
        let expected = cds.replace("c.", "r.");
        let split = expected.split_at(expected.find("delins").unwrap_or_else(|| {
            panic!(
                "the c. answer for [{payload}] is expected to stay a delins; got {cds} — \
                 if canonicalization legitimately changed shape, update this re-rendering"
            )
        }));
        let (prefix, edit) = split;
        let as_rna = format!("{prefix}{}", edit.to_lowercase().replace('t', "u"));
        assert_eq!(
            rna, as_rna,
            "r. must match the c. answer re-rendered as RNA (payload [{payload}])",
        );
    }
}

// =============================================================================
// An out-of-scope payload must be REFUSED on `r.`, not laundered through
// =============================================================================

/// The laundering half of the bug. `c.` refuses each of these shapes; before the
/// fix `r.` accepted them silently, so `ferro project --axis c` emitted a `c.`
/// string that ferro itself then failed to parse/normalize.
#[test]
fn rna_axis_refuses_the_payloads_the_c_axis_refuses() {
    for payload in [
        // p. axis — structurally invalid as a DNA-insertion payload.
        "NP_000079.2:p.10_15",
        // Intronic offsets and UTR markers — spec-undefined for expansion.
        "NM_000088.3:c.244-8_249",
        "NM_000088.3:c.100_*200",
        // pter/qter decorations have no concrete coordinate.
        "NC_000022.10:g.pter_100",
    ] {
        let on_cds = normalize(&format!("NM_000088.3:c.4_6delins[{payload}]"));
        assert!(
            on_cds.is_err(),
            "premise check: c. must refuse [{payload}]; got {on_cds:?}",
        );
        let on_rna = normalize(&format!("NM_000088.3:r.4_6delins[{payload}]"));
        assert!(
            on_rna.is_err(),
            "r. must refuse [{payload}] like c. does, not launder it through; got {on_rna:?}",
        );
    }
}

/// A payload naming an accession absent from the provider must surface the
/// lookup failure on `r.` too, rather than passing the variant through.
#[test]
fn rna_axis_surfaces_a_missing_payload_accession() {
    let result = normalize(&format!("NM_000088.3:{OUTER}delins[NM_999999.9:c.1_3]"));
    assert!(
        result.is_err(),
        "an unresolvable payload accession must error on r.; got {result:?}",
    );
}

// =============================================================================
// The `con` asymmetry noted in the issue's "Related" section
// =============================================================================

/// `con` is rewritten to `delins` on `r.` already, but without the expansion
/// step the rewritten payload stayed a bare `delins1_3` where `c.` produced the
/// literal. Same root cause, so it is fixed by the same wrapper — this pins it.
#[test]
fn rna_axis_con_edit_expands_to_a_literal() {
    let out = normalize(&format!("NM_000088.3:{OUTER}con1_3")).expect("r.-axis con must normalize");
    assert_eq!(
        out, "NM_000088.3:r.4_6delinsaug",
        "con must expand to the literal, not stay a bare delins1_3",
    );
    assert!(
        !out.contains("delins1_3"),
        "the unexpanded position range must not survive; got {out}",
    );
}
