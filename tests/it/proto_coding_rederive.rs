//! Coding-axis rederive prototype (issue #2155 follow-up).
//!
//! Exercises the experimental `FERRO_CODING_REDERIVE` path: a CDS-body `c.`
//! allele re-derived from its resulting sequence through the confluent
//! sequence-first core, so the `c.[..]` spellings of one variant converge on one
//! description. Like `proto_genomic_rederive`, these assert **convergence +
//! idempotence**, and each convergence case also pins the form the prototype
//! produces today. The form is still what the operator re-adjudicates; the pin
//! only makes a change of form a deliberate edit rather than a silent pass by
//! some other sequence-equivalent spelling.
//!
//! Convergence alone is a self-consistency check: two spellings could converge on
//! one *wrong* form and still pass. So each case also checks
//! **sequence-preservation against an independent oracle** —
//! `cis_apply_oracle::apply_with`, which reaches the bases through `hgvs_to_spdi`
//! and an SPDI splice rather than through the normalizer — so a merge onto the
//! wrong bases fails here even though both spellings would still "converge".
//!
//! The path is **off by default**; each test enables it via the env flag. Nextest
//! runs each test in its own process, so the flag (read once into a process-wide
//! `OnceLock`) is scoped to the test that sets it before its first normalize.

use crate::common::cis_apply_oracle::apply_with;
use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
use ferro_hgvs::reference::transcript::Strand;

/// Enable the experimental coding rederive for this test process. Must run
/// before the first `normalize`, which is why every test calls it first.
///
/// The flag is process-wide and cached, so setting it is only safe when this
/// test owns its process. nextest guarantees that and says so in
/// `NEXTEST_EXECUTION_MODE`; under `cargo test` every test in the `it` binary
/// shares one process, and setting the flag here would switch the coding path on
/// for whichever tests ran afterwards. Refuse rather than contaminate them.
fn enable_coding_rederive() {
    assert_eq!(
        std::env::var("NEXTEST_EXECUTION_MODE").as_deref(),
        Ok("process-per-test"),
        "FERRO_CODING_REDERIVE is a process-wide flag: run this test under \
         `cargo nextest` (process-per-test) so enabling it cannot leak into other tests"
    );
    // `set_var` is safe on edition 2021, and the assertion above means no other
    // test shares this process to read the flag.
    std::env::set_var("FERRO_CODING_REDERIVE", "1");
}

/// Assert that two `c.` spellings of one variant converge under coding rederive,
/// and that the convergence is **sequence-preserving** — checked against an oracle
/// independent of the normalizer, so a merge onto the wrong bases fails here.
///
/// `core` is the transcript sequence (`cds_start = 1`, so `c.N` is transcript
/// position `N`); `a`/`b` are the two spellings. Returns the converged form.
fn assert_converges(core: &str, cds_end: u64, a: &str, b: &str) -> String {
    let p = SyntheticBuilder::cds(core, 1, cds_end, Strand::Plus).build();

    let na = normalize_to_string(p.clone(), a);
    let nb = normalize_to_string(p.clone(), b);
    assert_eq!(
        na, nb,
        "the two c. spellings must converge under coding rederive:\n  {a}\n  {b}"
    );

    // Independent oracle: apply each input and the converged output to the
    // transcript sequence via SPDI (NOT the normalizer, and not the core's own
    // `verify_round_trip`, which shares the applier). Without this, a merge onto
    // the wrong sequence would still pass the equality above.
    let bases_a = apply_with(&p, core, a).expect("input A denotes a sequence");
    let bases_b = apply_with(&p, core, b).expect("input B denotes a sequence");
    let bases_out = apply_with(&p, core, &na).expect("converged output denotes a sequence");
    assert_eq!(
        bases_a, bases_b,
        "fixture sanity: the two inputs must denote the same bases"
    );
    assert_eq!(
        bases_out, bases_a,
        "coding rederive must preserve the denoted bases (independent SPDI oracle)"
    );

    // Idempotent: normalizing the converged form again is a no-op.
    assert_eq!(
        normalize_to_string(p, &na),
        na,
        "coding rederive must be idempotent"
    );
    na
}

#[test]
fn coding_rederive_converges_the_2155_pair() {
    enable_coding_rederive();
    // The #2155 window: c.10_17 = CTTAGTTA, replaced by AAACAAAC. Two spellings of
    // one variant — the minimal-alignment split and the spanning delins.
    let core = "ATGCCCGGGCTTAGTTAGGCCAATTCCGGAT";
    assert_eq!(&core[9..17], "CTTAGTTA", "fixture: c.10_17 == CTTAGTTA");
    let converged = assert_converges(
        core,
        30,
        "NM_TEST.1:c.[10_12delinsAA;14_16delinsCAA;17_18insC]",
        "NM_TEST.1:c.10_17delinsAAACAAAC",
    );
    // Pinned as well: convergence plus the oracle admit any idempotent,
    // sequence-equivalent form, and this is the one the prototype is offered with.
    assert_eq!(converged, "NM_TEST.1:c.10_17delinsAAACAAAC");
}

#[test]
fn coding_rederive_merges_the_one_codon_exception() {
    enable_coding_rederive();
    // codon 4 = c.10..c.12 = CTT; the one-amino-acid exception (`general.md:35`):
    // `[10C>T;12T>A]` and `c.10_12delinsTTA` both denote CTT -> TTA.
    let core = "ATGCCCGGGCTTAGTTAGGCCAATTCCGGAT";
    assert_eq!(&core[9..12], "CTT", "fixture: c.10_12 == CTT");
    let converged = assert_converges(
        core,
        30,
        "NM_TEST.1:c.[10C>T;12T>A]",
        "NM_TEST.1:c.10_12delinsTTA",
    );
    assert_eq!(converged, "NM_TEST.1:c.10_12delinsTTA");
}

/// `general.md:44` exempts a deletion or duplication at an exon/exon junction
/// from the 3' rule, and surface B refuses to re-derive a group whose members
/// straddle a junction (#1450). The coding rederive path must not undo that:
/// `NM_001234.1` (exons at transcript 1-15 / 16-30 / 31-44, `cds_start = 5`) has
/// a G-run over c.9-c.33, and re-deriving `c.[10del;13del]` across it would carry
/// both members past their junctions and merge them into one deletion in the
/// third exon. With the flag on, the pair must still come back as surface B's
/// clamped answer.
#[test]
fn coding_rederive_does_not_cross_an_exon_junction() {
    enable_coding_rederive();
    let normalizer =
        ferro_hgvs::Normalizer::new(ferro_hgvs::reference::MockProvider::with_test_data());
    let variant = ferro_hgvs::parse_hgvs("NM_001234.1:c.[10del;13del]").expect("input parses");
    assert_eq!(
        normalizer
            .normalize(&variant)
            .expect("normalization succeeds")
            .to_string(),
        "NM_001234.1:c.[11del;26del]"
    );
}

/// The input-side junction check above cannot see a **lone** member that sits
/// wholly inside one exon and is carried out of it by the derivation:
/// `derive_block_members` shifts within the sequence window, not the exon, so on
/// `NM_001234.1` (exon 1 ends at `c.11`, G-run over c.9-c.33) a `c.11del`
/// derived to `c.33del`, two exons on. Surface B clamps that shift to the
/// enclosing exon, and with the flag on the result must still be its answer.
#[test]
fn coding_rederive_keeps_a_lone_member_in_its_exon() {
    enable_coding_rederive();
    let normalizer =
        ferro_hgvs::Normalizer::new(ferro_hgvs::reference::MockProvider::with_test_data());
    for input in ["NM_001234.1:c.10del", "NM_001234.1:c.11del"] {
        let variant = ferro_hgvs::parse_hgvs(input).expect("input parses");
        assert_eq!(
            normalizer
                .normalize(&variant)
                .expect("normalization succeeds")
                .to_string(),
            "NM_001234.1:c.11del",
            "{input} must stay clamped to exon 1"
        );
    }
}

/// A `cds_start_NF` transcript has no confirmed `c.1` (#972), so the coding
/// frame the rederive path reasons in is undefined on it, and the flag must not
/// change what such an input normalizes to. `c.[11_12insGG;12_13insGG]` is a
/// shape where the two paths disagree on a complete transcript (the coding core
/// keeps the two insertions, surface B merges them). On a `cds_start_NF`
/// transcript the sequence-first pass declines before either runs, so with the
/// flag on the input passes through unchanged, as `normalize_cds` passes a lone
/// member through and as the flag-off path does
/// (`issue_972_cds_start_nf_decline::a_cds_start_nf_allele_passes_through_the_sequence_first_pass`).
#[test]
fn coding_rederive_passes_a_cds_start_nf_allele_through() {
    use ferro_hgvs::reference::ReferenceProvider;
    enable_coding_rederive();
    let core = "ATGCCCGGGCTTAGTTAGGCCAATTCCGGAT";
    let mut provider = SyntheticBuilder::cds(core, 1, 30, Strand::Plus).build();
    let mut transcript = (*provider
        .get_transcript("NM_TEST.1")
        .expect("fixture transcript"))
    .clone();
    transcript.cds_start_incomplete = true;
    provider.add_transcript(transcript);

    let normalizer =
        ferro_hgvs::Normalizer::with_config(provider, ferro_hgvs::NormalizeConfig::lenient());
    let variant =
        ferro_hgvs::parse_hgvs("NM_TEST.1:c.[11_12insGG;12_13insGG]").expect("input parses");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient normalization succeeds");
    assert_eq!(
        result.result.to_string(),
        "NM_TEST.1:c.[11_12insGG;12_13insGG]"
    );
}

/// The rederive path templates on the first member and renders in the coding
/// frame, so it must admit an allele only when **every** member is `c.`. A
/// mixed-axis cis allele on one transcript parses, and `[c.10C>T;n.12T>A]`
/// names the same codon the one-codon exception merges in
/// `coding_rederive_merges_the_one_codon_exception` (`cds_start = 1`, so
/// `n.12` is `c.12`). With the flag on it must come back exactly as surface B
/// leaves it, members and axes unchanged, rather than as one `c.` delins.
#[test]
fn coding_rederive_leaves_a_mixed_axis_allele_to_surface_b() {
    enable_coding_rederive();
    let core = "ATGCCCGGGCTTAGTTAGGCCAATTCCGGAT";
    let p = SyntheticBuilder::cds(core, 1, 30, Strand::Plus).build();
    for input in [
        "[NM_TEST.1:c.10C>T;NM_TEST.1:n.12T>A]",
        "[NM_TEST.1:c.10_12delinsAA;NM_TEST.1:n.14_16delinsCAA;NM_TEST.1:c.17_18insC]",
        "[NM_TEST.1:c.10_12delinsAA;NM_TEST.1:r.14_16delinscaa]",
    ] {
        assert_eq!(normalize_to_string(p.clone(), input), input);
    }
}
