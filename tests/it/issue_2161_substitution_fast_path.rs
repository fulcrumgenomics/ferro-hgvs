//! Regression coverage for the lone-substitution fast path (#2169, part of #2161).
//!
//! `Normalizer::sequence_first_pass` short-circuits a lone single-base
//! substitution on a nucleotide axis: `is_lone_single_base_substitution`
//! returns `true` and the pass returns `None` ("already canonical") without
//! running the apply + align + partition machinery. The claim under test is
//! that this is **output-preserving**, not a heuristic: a substitution changes
//! exactly one reference base, and the minimal, unique re-derivation of a
//! one-base difference is that same substitution — there is no competitor
//! spelling to prefer and no room to shift — so the full pass would hand back
//! the identical variant. A substitution is therefore a fixed point of
//! normalization: `normalize(sub) == sub`, byte-for-byte.
//!
//! Because the output is byte-identical to the input, denoted-sequence
//! equivalence between input and output is trivially satisfied — the two are
//! the same string — so the fixed-point assertion below *is* the
//! denoted-sequence check. The remaining risk the fast path could introduce is
//! **skipping validation that lives outside the sequence-first pass** (e.g. the
//! reference-base check), so a reference-mismatch case is included and asserts
//! the `RefSeqMismatch` warning still fires.
//!
//! The expected outputs here were cross-checked against the pre-fast-path
//! behaviour: with the fast-path early-return temporarily removed from
//! `sequence_first_pass`, every case below produces the identical output and
//! diagnostics (the full pass agrees), which is what makes the short-circuit a
//! proven no-op rather than a shortcut. See #2161 / #2169.
//!
//! Coverage: lone substitutions on all five nucleotide axes — `g.`, `m.`, `c.`,
//! `n.`, `r.` — through both `normalize()` and `normalize_with_diagnostics()`,
//! plus an **uncertain-edit** substitution on `g.` and `c.` (`g.(..>..)`, which
//! also takes the fast path because `Mu::inner()` is `Some` for `Uncertain`),
//! and a **reference-mismatch** substitution. Each asserts the fixed point,
//! idempotency, and reparsability.

use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, NormalizeConfig, Normalizer};

// ---------------------------------------------------------------------------
// Fixtures — non-repetitive sequences so no test position sits in a tract a
// shuffle could touch (a substitution never shifts, but this removes all doubt
// about a homopolymer interaction). Bases at every test position were computed
// from these constants, so the stated reference base in each string is correct
// and the substitution is a clean (non-mismatch) edit unless stated otherwise.
// ---------------------------------------------------------------------------

const G_ACC: &str = "NC_TEST.1";
const G_SEQ: &str = "GCATTCAGGTACCGATTGCACGTGCTAAGCATCAGTTGCTCATGACTGGCATCAGTGA";

const M_ACC: &str = "NC_012920.1";
const M_SEQ: &str = "AACCGGTTACGTTGCACATGCTAGGCATCAGTTGCACATGACTGGCATCAGTAGTCA";

const TX_ACC: &str = "NM_TEST.1";
const TX_SEQ: &str = "ATGCAGTCGATTGCACCTAGCATCAGTTGCACATGACTGGCATCAGTAAGCTTGCAA";

/// One provider carrying the genomic contig, the mitochondrial contig, and the
/// transcript, so a single normalizer serves every axis. `cds_start = 1` makes
/// `c.N` land on transcript position `N`.
fn provider() -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence(G_ACC, G_SEQ);
    p.add_genomic_sequence(M_ACC, M_SEQ);
    let len = TX_SEQ.len() as u64;
    let tx = Transcript::new(
        TX_ACC.to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        TX_SEQ.to_string(),
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    p.add_transcript(tx);
    p
}

/// Assert that `input` is a fixed point of normalization through **both** entry
/// points, is idempotent, and reparses — the property the substitution fast
/// path must preserve. `expect_clean` requires no warnings (a clean edit);
/// otherwise warnings are left to the caller to inspect.
fn assert_substitution_is_a_fixed_point(nz: &Normalizer<MockProvider>, input: &str) {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));

    // normalize(): byte-identical output.
    let out = nz
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"));
    assert_eq!(
        out.to_string(),
        input,
        "a lone substitution must be a fixed point of normalize()"
    );

    // normalize_with_diagnostics(): same output, and no warnings for a clean edit.
    let diag = nz
        .normalize_with_diagnostics(&variant)
        .unwrap_or_else(|e| panic!("normalize_with_diagnostics {input}: {e}"));
    assert_eq!(
        diag.result.to_string(),
        input,
        "normalize_with_diagnostics() must agree with normalize() on the fixed point"
    );
    assert!(
        diag.warnings.is_empty(),
        "a clean substitution must raise no warnings; got {:?}",
        diag.warnings
    );

    // Idempotency: normalizing the output again yields the same string.
    let reparsed = parse_hgvs(&out.to_string()).expect("output reparses");
    let again = nz.normalize(&reparsed).expect("re-normalize");
    assert_eq!(
        again.to_string(),
        input,
        "normalization must be idempotent on a substitution"
    );
}

// ---------------------------------------------------------------------------
// One clean lone substitution per nucleotide axis.
// ---------------------------------------------------------------------------

#[test]
fn genomic_substitution_is_a_fixed_point() {
    assert_substitution_is_a_fixed_point(&Normalizer::new(provider()), "NC_TEST.1:g.20A>G");
}

#[test]
fn mitochondrial_substitution_is_a_fixed_point() {
    assert_substitution_is_a_fixed_point(&Normalizer::new(provider()), "NC_012920.1:m.25G>A");
}

#[test]
fn coding_substitution_is_a_fixed_point() {
    assert_substitution_is_a_fixed_point(&Normalizer::new(provider()), "NM_TEST.1:c.15A>G");
}

#[test]
fn noncoding_substitution_is_a_fixed_point() {
    assert_substitution_is_a_fixed_point(&Normalizer::new(provider()), "NM_TEST.1:n.22A>G");
}

#[test]
fn rna_substitution_is_a_fixed_point() {
    assert_substitution_is_a_fixed_point(&Normalizer::new(provider()), "NM_TEST.1:r.10a>g");
}

// ---------------------------------------------------------------------------
// Uncertain-edit substitutions also take the fast path (`Mu::inner()` is `Some`
// for `Uncertain`). They must be fixed points too — the `(...)` uncertainty
// wrapper is preserved, not stripped.
// ---------------------------------------------------------------------------

#[test]
fn uncertain_genomic_substitution_is_a_fixed_point() {
    assert_substitution_is_a_fixed_point(&Normalizer::new(provider()), "NC_TEST.1:g.(33C>T)");
}

#[test]
fn uncertain_coding_substitution_is_a_fixed_point() {
    assert_substitution_is_a_fixed_point(&Normalizer::new(provider()), "NM_TEST.1:c.(28T>C)");
}

// ---------------------------------------------------------------------------
// A reference-mismatch substitution: the fast path skips the sequence-first
// re-derivation, but the reference-base check lives elsewhere and must still
// fire. Per #1052 the substitution arm keeps the stated reference base in the
// output (`corrected: false`), so the string is still a fixed point while a
// `RefSeqMismatch` warning is raised.
// ---------------------------------------------------------------------------

#[test]
fn reference_mismatch_substitution_still_validates() {
    // g.45 is `A` on the contig; state `T` to force a mismatch.
    let input = "NC_TEST.1:g.45T>G";
    let nz = Normalizer::with_config(provider(), NormalizeConfig::lenient());
    let variant = parse_hgvs(input).expect("parse");

    let out = nz.normalize(&variant).expect("normalize (lenient)");
    assert_eq!(
        out.to_string(),
        input,
        "a substitution keeps its stated reference base (corrected: false), so the string is unchanged"
    );

    let diag = nz
        .normalize_with_diagnostics(&variant)
        .expect("normalize_with_diagnostics (lenient)");
    assert_eq!(diag.result.to_string(), input);
    let refseq_mismatches = diag
        .warnings
        .iter()
        .filter(|w| matches!(w, NormalizationWarning::RefSeqMismatch { .. }))
        .count();
    assert_eq!(
        refseq_mismatches, 1,
        "the reference-base check must still fire despite the fast path; warnings={:?}",
        diag.warnings
    );
}

// ---------------------------------------------------------------------------
// A single-member non-substitution must NOT take the fast path — it is not a
// fixed point in a repeat, so this guards that the predicate is narrow (only
// substitutions short-circuit). A deletion inside a homopolymer 3'-shifts.
// ---------------------------------------------------------------------------

#[test]
fn a_deletion_in_a_repeat_still_runs_the_full_pass() {
    // Contig with a homopolymer: a del at its 5' end must shift to the 3' end,
    // proving the fast path did not swallow a non-substitution.
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_HP.1", "GCGC AAAAAAAA GCGC".replace(' ', ""));
    let nz = Normalizer::new(p);
    let variant = parse_hgvs("NC_HP.1:g.5del").expect("parse");
    let out = nz.normalize(&variant).expect("normalize");
    assert_ne!(
        out.to_string(),
        "NC_HP.1:g.5del",
        "a deletion in a homopolymer must 3'-shift — it must not be treated as a fixed point"
    );
    // Sanity: the type is still a genome variant (guards the helper wiring).
    assert!(matches!(out, HgvsVariant::Genome(_)));
}
