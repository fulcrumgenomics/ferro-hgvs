//! `normalize` must return a description `parse_hgvs` accepts.
//!
//! This is not a canonicalization question. Several open arguments about cis
//! normalization are about *which* of several valid spellings is canonical
//! (#1260, #1261, #1262), and their answers depend on policy calls about the
//! delins corpus. This one does not: whatever the right canonical form turns
//! out to be, `normalize` is documented to return a valid HGVS description, and
//! every consumer chaining a second call depends on that.
//!
//! It is also a hole in the idempotency oracle rather than a separate concern.
//! `FERRO_ASSERT_IDEMPOTENT` verifies by re-normalizing its own output, which it
//! cannot do for an output that will not parse — so an unparseable result is
//! precisely the shape that oracle cannot see. The re-parse check therefore sits
//! at the same seam (`normalize_core_checked`'s single exit) behind its own
//! flag, `FERRO_ASSERT_REPARSE`, so it is clear which invariant broke.
//!
//! The violation that motivated it (#1263) is fixed: a duplication whose span
//! collided with a sibling's bases is now re-spelled as the equivalent
//! insertion. `FERRO_ASSERT_REPARSE=1` is clean across the whole suite and is
//! set in CI alongside `FERRO_ASSERT_IDEMPOTENT=1`.
//!
//! Two exemptions in the oracle are load-bearing for that, and both are tracked
//! by #1264 rather than left implicit:
//!
//! - `0` and `?` (null and unknown allele) are legal HGVS but not standalone
//!   parseable, since `parse_hgvs` wants an accession. A limit of the entry
//!   point, not a malformed output.
//! - A description the parser accepts inbound but rejects outbound — an empty
//!   allele rendering as `[]`, or an `ins` whose non-adjacent anchors are
//!   tolerated on the way in — is a parse/display round-trip asymmetry that
//!   normalization merely passes through. The oracle fires only when the
//!   *input* re-parses and the *output* does not.

use crate::common::synthetic::{padded, SyntheticBuilder};
use ferro_hgvs::{parse_hgvs, Normalizer};

fn normalize(core: &str, input: &str) -> String {
    let normalizer = Normalizer::new(SyntheticBuilder::genomic(core).build());
    let variant = parse_hgvs(input).expect("input parses");
    format!(
        "{}",
        normalizer.normalize(&variant).expect("normalize succeeds")
    )
}

/// The reported violation (#1263), now fixed.
///
/// `262_263insA` canonicalises to `262dup`, and the sibling-crossing clamp lets
/// a deletion sit flush against that junction — correct while the member is
/// spelled `ins`, an overlap once it is spelled `dup`, because a duplication's
/// span *is* the base it copies. ferro then rejected its own output.
///
/// The fix is a re-spelling, not a repositioning: `Xdup` and
/// `X_X+1ins<ref[X]>` denote the same edit, and the insertion form claims no
/// base. Tightening the junction clamp by one base also fixes this and breaks
/// #1135, where a deletion must reach the junction for the self-cancelling
/// del+ins merge to fire — which is why the fix reads the reference instead of
/// moving a boundary.
#[test]
fn a_del_beside_a_dup_re_spells_instead_of_colliding() {
    let core = "GCAACAGCGTAAAC";
    let input = "NC_TEST.1:g.[260_261del;262_263insA;264delinsAC]";
    let normalized = normalize(core, input);

    // The re-spelling exposes a del and an ins that then cancel further, and the
    // fixed point is the re-derivation of what remains — which is why the repair
    // runs inside the normalization loop and not after it.
    //
    // Was `NC_TEST.1:g.[261del;263_264insA]` until the input-relative weight
    // bound was deleted
    // (`rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`). Same
    // block as `issue_1284_transcript_axis_collision`, reached from an unrelated
    // input: `CAG` -> `AGA`, EQUAL length 3/3, three consecutive changed
    // columns, so `DNA/delins.md:16` gives the single `delins`. The old form was
    // the input's spelling surviving a refusal, not a derived answer.
    assert_eq!(
        normalized, "NC_TEST.1:g.261_263delinsAGA",
        "expected the colliding `262dup` to be re-spelled and then re-derived"
    );
    parse_hgvs(&normalized).unwrap_or_else(|e| {
        panic!("`{input}` normalized to `{normalized}`, which ferro cannot re-parse: {e}")
    });
    assert_eq!(
        applied(core, &normalized),
        applied(core, input),
        "`{input}` -> `{normalized}` changed the sequence"
    );
}

/// Apply `descriptor` to the synthetic reference independently of the
/// normalizer, via SPDI triples, and return the **whole** resulting sequence.
///
/// Delegates to the shared oracle rather than re-deriving the splice walk: that
/// walk is what makes "declines an overlapping description" true, and it is
/// already the single copy the sibling-crossing tests rest on.
///
/// The full sequence, not a window around the core. This previously returned
/// `skip(254).take(20)`, which is narrower than the contract its caller states
/// ("changed the sequence"): the cases here change length, so a member that
/// shifted out of the window, or content pushed past its end, compared equal
/// while the sequences differed.
fn applied(core: &str, descriptor: &str) -> String {
    let provider = SyntheticBuilder::genomic(core).build();
    let reference = padded(core);
    crate::common::cis_apply_oracle::apply_with(&provider, &reference, descriptor)
        .unwrap_or_else(|| {
            panic!("`{descriptor}` has no well-defined resulting sequence (overlapping or unconvertible)")
        })
}

/// The invariant itself, over the cases that already satisfy it — so the
/// round-trip is exercised rather than merely described, and a regression that
/// broke a currently-good shape would be caught here.
#[test]
fn normalized_output_re_parses() {
    for (core, input) in [
        ("AAAAAA", "NC_TEST.1:g.[258_259insC;259_260insC]"),
        ("AAAAAA", "NC_TEST.1:g.258_260delinsACACA"),
        ("AAAAAA", "NC_TEST.1:g.[258A>C;260del]"),
        ("AAAAAA", "NC_TEST.1:g.258_260delinsCA"),
        ("CAGTATGCAGGCAA", "NC_TEST.1:g.[258_259insA;259_260insAG]"),
        ("GCAACAGCGTAAAC", "NC_TEST.1:g.[260_261del;264delinsAC]"),
        ("GCAACAGCGTAAAC", "NC_TEST.1:g.262_263insA"),
    ] {
        let normalized = normalize(core, input);
        parse_hgvs(&normalized).unwrap_or_else(|e| {
            panic!("`{input}` normalized to `{normalized}`, which ferro cannot re-parse: {e}")
        });
    }
}

/// Each of the four sibling passes must reach a fixed point on a **compound**
/// cis allele, not merely on a lone member.
///
/// The passes run inside `normalize_allele`'s loop and rewrite members in place
/// — a clamp pulls a span back, a demotion re-spells a repeat, a respell turns a
/// `dup` into an `ins`. Any of those can expose further work for the next pass,
/// so "the output is correct" and "the output is stable" are different claims
/// and only the second one rules out a rewrite that oscillates. The suite's
/// other idempotency coverage is single-member, which cannot reach these passes
/// at all: every one of them bails below two members.
///
/// `respell_colliding_duplications` is the case that motivated this — it had to
/// move inside the loop precisely because running it once at the end left a
/// non-fixed-point result, the del and ins it exposes cancelling further on the
/// next pass.
#[test]
fn every_sibling_pass_is_idempotent_on_a_compound_allele() {
    use crate::common::cis_apply_oracle::normalize as normalize_template;

    // (pass under test, reference, input)
    let cases = [
        (
            "clamp_sibling_crossing_shifts",
            "TAATATATATAATATATATT",
            "TEMPLATE:g.[3_4del;9del]",
        ),
        (
            "demote_repeats_spanning_siblings",
            "TTTTTTTTTAATATATTTTA",
            "TEMPLATE:g.[1_2del;4del]",
        ),
        (
            "clamp_sibling_crossing_junctions",
            "ACAAAAAAAACGTACGTACG",
            "TEMPLATE:g.[2_3insA;5A>G]",
        ),
    ];
    for (pass, reference, input) in cases {
        let once = normalize_template(reference, input);
        let twice = normalize_template(reference, &once);
        assert_eq!(
            twice, once,
            "{pass}: `{input}` -> `{once}` is not a fixed point"
        );
    }

    // `respell_colliding_duplications` needs the padded builder, since the
    // shape that reaches it is the #1263 del + dup collision pinned above.
    for (core, input) in [
        ("GCAACAGCGTAAAC", "NC_TEST.1:g.[260_261del;264delinsAC]"),
        (
            "GCAACAGCGTAAAC",
            "NC_TEST.1:g.[260_261del;262_263insA;264delinsAC]",
        ),
    ] {
        let once = normalize(core, input);
        let twice = normalize(core, &once);
        assert_eq!(
            twice, once,
            "respell_colliding_duplications: `{input}` -> `{once}` is not a fixed point"
        );
    }
}

/// The sibling clamp must not depend on how a substitution is *spelled*.
///
/// `9>G` (`SubstitutionNoRef`) replaces the base at 9 exactly as `9A>G` does,
/// but it was a separate `NaEdit` variant and was missing from
/// `claims_reference_bases`, so the clamp could not see it. The two spellings of
/// one variant therefore normalized differently — and the `>G` one to a
/// *malformed* result, two members both claiming position 9:
///
/// ```text
/// g.[4del;9A>G]  ->  g.8_9delinsG    clamped, then coalesced
/// g.[4del;9>G]   ->  g.[9>G;9del]    overlapping; the SPDI oracle declines it
/// ```
///
/// Both spellings are pinned here, so the contract is "same variant, same
/// answer" rather than "this one case happens to work".
///
/// Only the `9A>G` spelling gets the SPDI sequence-preservation check:
/// `hgvs_to_spdi` has no conversion for `SubstitutionNoRef`, so `apply` declines
/// the `9>G` input and the oracle cannot speak to it. The string equality is the
/// contract under test anyway — the two inputs denote one variant, so one
/// answer — and the sequence is pinned through the spelling the oracle can serve.
#[test]
fn the_clamp_does_not_depend_on_how_a_substitution_is_spelled() {
    use crate::common::cis_apply_oracle::{apply, normalize as normalize_template};

    let seq = "GGGAAAAAAGGG";
    let expected = "TEMPLATE:g.8_9delinsG";
    for input in ["TEMPLATE:g.[4del;9A>G]", "TEMPLATE:g.[4del;9>G]"] {
        assert_eq!(
            normalize_template(seq, input),
            expected,
            "spelling changed the answer for {input}"
        );
    }
    let spelled = "TEMPLATE:g.[4del;9A>G]";
    assert_eq!(
        apply(seq, expected).expect("output applies"),
        apply(seq, spelled).expect("input applies"),
        "`{spelled}` -> `{expected}` changed the sequence"
    );
}

/// `NPaddedDeletion` (`delN[k]`) claims the bases under its span, so a sibling
/// must not shift onto them — the same contract a plain `del` gets.
///
/// Asserted on the emitted string rather than through the SPDI oracle the other
/// sibling tests use: `hgvs_to_spdi` has no conversion for `delN`, so `apply`
/// declines both sides and the oracle cannot distinguish right from wrong here.
/// The string is enough — `4del` sits in the A-run `4..9` and its standalone
/// 3'-most position is `9del`, on top of the sibling. Stopping at `8del` is the
/// clamp firing; before `NPaddedDeletion` was added to `claims_reference_bases`
/// the sibling was invisible and the deletion landed on it.
#[test]
fn an_n_padded_deletion_is_a_sibling_barrier() {
    use crate::common::cis_apply_oracle::normalize as normalize_template;

    let seq = "GGGAAAAAAGGG";
    assert_eq!(
        normalize_template(seq, "TEMPLATE:g.[4del;9delN[1]]"),
        "TEMPLATE:g.[8del;9delN[1]]"
    );
}

/// The mitochondrial axis takes the same sibling passes as `g.`, and nothing
/// else in the suite exercises `CisKind::Mt` through them.
#[test]
fn the_mitochondrial_axis_clamps_like_the_genomic_one() {
    use ferro_hgvs::MockProvider;

    let seq = "GGGAAAAAAGGG";
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_012920.1", seq.to_string());
    let normalizer = Normalizer::new(provider);
    let input = "NC_012920.1:m.[4del;9A>G]";
    let normalized = normalizer
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize")
        .to_string();
    // Same shape as the genomic `a_deletion_does_not_shift_onto_a_substituted_base`:
    // clamped to `8del`, adjacent to `9A>G`, the two render as one delins.
    assert_eq!(normalized, "NC_012920.1:m.8_9delinsG");
}
