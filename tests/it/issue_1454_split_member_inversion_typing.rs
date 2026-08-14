//! #1454 — a member the canonical split *creates* must be typed when it is
//! created, not one normalization pass later.
//!
//! ```text
//! NM_TEST.1:c.10_15delinsTAATAT
//!   on main   -> NM_TEST.1:c.[10_12delinsTAA;13_15delinsTAT]
//!   pass 2    -> NM_TEST.1:c.[10_12delinsTAA;13_15inv]
//! ```
//!
//! `c.13_15` is `ATA` and `revcomp(ATA)` is `TAT`, so the second member is an
//! exact inversion. The first output was therefore not a fixed point, which
//! `FERRO_ASSERT_IDEMPOTENT=1` catches and the other two seam oracles do not —
//! the string re-parses and every coordinate is in range.
//!
//! # The mechanism is in the split, not in the derivation
//!
//! The issue attributes this to the alternation loop in
//! `merge::canonicalize_from_sequence`. That loop never runs for this input:
//! `Normalizer::sequence_first_pass` admits a lone member only on `g.`/`m.`
//! (`is_splittable_single_member`), so a lone `c.` delins declines on its first
//! line. Instrumenting the loop against the reproducer logs nothing.
//!
//! The split comes from `apply_canonical_split` instead:
//!
//! 1. `rules::decompose_delins` reports
//!    `[Sub@10; IdentityAt@11; Sub@12; Sub@13; Sub@14; Sub@15]`. It types an
//!    inversion only when a **whole** maximal contiguous mismatch run is a
//!    reverse complement, and deliberately does not carve a reverse-complement
//!    *sub*-run out of a longer contiguous change (#1034, #1040). Positions
//!    12-15 are one run, `TATA -> ATAT`, which is not one.
//! 2. `build_split_variants` then applies the codon-frame exception
//!    (`general.md:35-38`) and **consumes** `[Sub@10; IdentityAt@11; Sub@12]`
//!    into `10_12delinsTAA`.
//! 3. That leaves `Sub@13,14,15` — a span that did not exist when step 1 typed
//!    the input, and which *is* a whole reverse complement. It fell to
//!    `flush_substitution_run`, which emitted `Delins` unconditionally.
//!
//! So the re-grouping creates spans, and nothing re-typed them. `general.md:56`
//! ranks inversion above deletion-insertion and `delins.md:5` defines a delins
//! as a replacement "which is not a substitution or inversion", so `delins` over
//! a reverse-complement span is not a description the spec admits at all.
//!
//! The fix routes both of `build_split_variants`' multi-base groupings through
//! `rules::canonicalize_delins` — the same function the *next* pass reaches via
//! `normalize_na_edit`'s Delins arm. Answering the question the second pass
//! would ask, with the function it would ask it of, is what makes the two passes
//! unable to disagree.
//!
//! # Measurements
//!
//! Over 479,682 rows — every `delins` of length 3, 4 and 5 at every position of
//! three 63/64 nt cores, on the `c.` and `n.` axes, counted as
//! `normalize(normalize(x)) != normalize(x)`:
//!
//! | | on `main` (`a530cdaf`) | with the fix |
//! |---|---|---|
//! | outputs that are not fixed points | 267 | **0** |
//!
//! All 267 are the one shape: a two-member split whose trailing member is an
//! exact reverse complement, re-typed `inv` only on the following pass. They
//! span 188 distinct payload pairs and both the 2 nt and 3 nt inversion widths
//! (`c.1_5delinsATCAA -> c.[1_3delinsATC;4_5delinsAA]`, whose `4_5` is `TT`).
//! Every row of the census sits on the `c.` axis; the `n.` half contributes
//! none, for the reason `an_axis_without_a_reading_frame_is_not_re_split` pins.
//!
//! Representation stability, via `examples/dump_normalized_corpus.rs` against
//! `a530cdaf`: **52 of 78,028 rows moved (0.1%)**, all of them in the
//! `delins_hiding_an_inversion` family and all of the same shape — a trailing
//! `delins` member becoming `inv`. Every moved row is an output that `main`
//! itself rewrote on a second pass, so this converges on the string ferro
//! already produced rather than inventing a third one; it is nonetheless a moved
//! string for a consumer that normalizes once and stores the result.

//! # What #1524 changed here
//!
//! The reproducer above no longer reproduces: `c.[10_12delinsTAA;13_15inv]`
//! puts members on the consecutive nucleotides 12 and 13, which
//! `delins.md:16` forbids, so `c.10_15delinsTAATAT` is now one member and its
//! own normal form. Four expectations in this file moved with it and each says
//! so at its assertion.
//!
//! The fix this file exists for is intact and still exercised. Only the set of
//! groupings that can *create* a span narrowed: a substitution-run member is
//! exactly one maximal mismatch run, which `decompose_delins` has already
//! typed, so the codon-frame triplet is now the sole producer — and
//! `codon_triplet_over_an_ambiguous_centre_is_an_inversion`, retargeted onto a
//! separated trailing run, is its guard.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// The issue's core. Positions 10-15 are `AATATA`; `c.13_15` is `ATA`, whose
/// reverse complement is `TAT`.
const CORE: &str = "TTTTTTTTTAATATATTTTAATATAATTAAAAAAATAATTTTTATAAATATATTATTTTAAAAA";

/// Positions 10/11/12 are `A`/`N`/`C`.
///
/// `N` is its own complement, which is what makes the codon-frame branch of the
/// fix reachable at all: that branch emits `[Sub@i; IdentityAt@i+1; Sub@i+2]` as
/// one three-base member whose centre is unchanged, so `alt[1] == ref[1]`, and a
/// reverse complement needs `alt[1] == complement(ref[1])`. Only the
/// self-complementary IUPAC codes (`W`, `S`, `N`) satisfy both. The endpoints
/// are `A` and `C` rather than a complementary pair, so the span is a genuine
/// inversion instead of collapsing to an identity.
const AMBIGUOUS_CORE: &str = "GGGGGGGGGANCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

/// One transcript carrying `core`, coding or not.
///
/// The coding arm sets the CDS to the whole transcript so `c.N`, `n.N` and `r.N`
/// all name transcript position `N` — which is what makes the cross-axis
/// comparisons below controlled: the records differ only in whether a reading
/// frame exists, never in which bases a position names.
fn provider(accession: &str, core: &str, coding: bool) -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = core.to_string();
    let length = sequence.len() as u64;
    let (cds_start, cds_end) = if coding {
        (Some(1), Some(length))
    } else {
        (None, None)
    };
    provider.add_transcript(Transcript::new(
        accession.to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        sequence,
        cds_start,
        cds_end,
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

fn normalize(provider: &MockProvider, descriptor: &str) -> String {
    let variant = parse_hgvs(descriptor).expect("fixture must parse");
    Normalizer::new(provider.clone())
        .normalize(&variant)
        .expect("fixture must normalize")
        .to_string()
}

/// Normalize, and assert the result is its own normal form.
///
/// The fixed-point half is asserted rather than left to
/// `FERRO_ASSERT_IDEMPOTENT`, so these tests fail in a plain `cargo nextest run`
/// too. CI's gating `Test` job runs without the oracles, so an oracle-only guard
/// for #1454 would be green in the one job that blocks a merge.
fn normalized_fixed_point(provider: &MockProvider, descriptor: &str) -> String {
    let once = normalize(provider, descriptor);
    let twice = normalize(provider, &once);
    assert_eq!(
        twice, once,
        "`{descriptor}` normalized to a form that is not a fixed point",
    );
    once
}

/// The original reproducer, **retired as such by #1524** and kept as the row
/// that records why.
///
/// The shape it pinned — a span the codon exception *leaves behind*, typed
/// `inv` — is no longer reachable from a substitution run. `c.13` touches
/// `c.12`, so `delins.md:16` puts them in one member and there is no leftover
/// span to type. `c.10_15delinsTAATAT` is therefore its own normal form.
///
/// This is not a coverage loss for #1454's actual fix: a member the split
/// *creates* must still be typed when it is created, and
/// `codon_triplet_over_an_ambiguous_centre_is_an_inversion` below is the guard
/// for it. What changed is which grouping can create such a member. Since
/// #1524, only the codon-frame triplet can — a substitution run member is
/// exactly one maximal mismatch run, and `decompose_delins` has already typed
/// those.
///
/// # MOVED BY #1835 AND RETURNED BY #1878
///
/// The span `c.10_15` reads `AATATA` against a payload of `TAATAT` — the same
/// bases rotated one position. Read column-wise that is five mismatches; read as
/// an alignment it is one inserted `T` at the 5' end and one deleted `A` at the
/// 3', cost two. #1835 took the second reading and pinned `c.[9dup;15del]`.
///
/// `rulings[equal-length-block-column-correspondence-is-unique]` (decided) says
/// the first reading is not a competing choice but the only one: six reference
/// bases against a six-base payload have exactly one column correspondence, so
/// what changed is a fact rather than an alignment ferro selects, and the search
/// that found the rotation had nothing to decide. `delins.md:16` then types five
/// consecutive changed columns as one `delins`, which is this row's original and
/// restored answer. The name — "is now one member" — is accurate again.
#[test]
fn a_span_the_codon_exception_used_to_leave_behind_is_now_one_member() {
    let provider = provider("NM_TEST.1", CORE, true);
    assert_eq!(
        normalized_fixed_point(&provider, "NM_TEST.1:c.10_15delinsTAATAT"),
        "NM_TEST.1:c.10_15delinsTAATAT",
    );
}

/// Every spelling of the variant reaches that one form, including the two-member
/// spelling and the already-typed one.
///
/// Pinned as convergence rather than as three independent strings: a test that
/// pinned only the first would go green if a later change moved the others onto
/// a different form, which is the opposite of what this issue is about.
///
/// #1835 moved the converged form to `c.[9dup;15del]` and #1878 returned it to
/// the spanning `c.10_15delinsTAATAT` — see
/// [`a_span_the_codon_exception_used_to_leave_behind_is_now_one_member`] for the
/// argument. **The property this test asserts survived both moves untouched, and
/// it is the one that matters here**: all three spellings, including the one
/// authored with an explicit `inv` member, reach a single form. Confluence is
/// what #1454 is about; only the representative ever moved.
#[test]
fn the_other_spellings_of_the_variant_converge_on_the_same_form() {
    let provider = provider("NM_TEST.1", CORE, true);
    for input in [
        "NM_TEST.1:c.10_15delinsTAATAT",
        "NM_TEST.1:c.[10_12delinsTAA;13_15delinsTAT]",
        "NM_TEST.1:c.[10_12delinsTAA;13_15inv]",
    ] {
        assert_eq!(
            normalized_fixed_point(&provider, input),
            "NM_TEST.1:c.10_15delinsTAATAT",
            "spelling `{input}` did not converge",
        );
    }
}

/// The `r.` axis reaches the same `inv`, which is what pins the `U`/`T` fold.
///
/// `r.` is codon-frame-aware on a coding transcript (#275), so it takes the same
/// branch as `c.`. Its alt bases come from the author's literal and so are `u`,
/// while the transcript reference is `T`; `rules::is_revcomp` maps `A` to `T`
/// and never to `U`, so without the fold this member reads as an ordinary delins
/// and the `r.` axis silently keeps the defect.
///
/// # #1835 moved this to `r.[9dup;15del]`; #1878 returned it
///
/// The `r.` axis agrees with `c.` under either pinning, which is the stronger
/// half of what this row buys. What #1878 restores is the `inv` label: with the
/// unique correspondence the member is one contiguous reverse-complemented run
/// again, so this row witnesses the `U`/`T` fold the way it was written to — an
/// unfolded comparison cannot see that `uaauau` reverse-complements `AATATA`.
#[test]
fn the_rna_axis_reaches_the_same_inversion() {
    let provider = provider("NM_TEST.1", CORE, true);
    assert_eq!(
        normalized_fixed_point(&provider, "NM_TEST.1:r.10_15delinsuaauau"),
        "NM_TEST.1:r.10_15delinsuaauau",
    );
}

/// The discriminating case for *where* the fix fires, and the reason the two
/// tests above pin a reading-frame axis specifically.
///
/// Without a reading frame there is no codon-frame exception, so nothing
/// consumes part of a maximal mismatch run and no new span is ever created —
/// `decompose_delins` has already tested the only spans that reach
/// `flush_substitution_run`, and its verdict stands. The same core and the same
/// span therefore split differently and produce no `inv` at all.
///
/// This is what an over-general fix would break: a rule that typed *any*
/// multi-base member as `inv` whenever its payload happened to be a reverse
/// complement would have to change this row too, and it must not.
///
/// # #1835 ERASED THE CONTRAST THIS ROW DRAWS; #1878 RESTORES IT
///
/// Under the flip both axes reached `[9dup;15del]` — the same form the `c.` axis
/// reached — so a frame-aware and a frameless axis produced identical output over
/// this core and removing the reading-frame condition would not have reddened
/// this row. That was recorded there as a real coverage loss.
///
/// It is a loss no longer. With the unique correspondence
/// (`rulings[equal-length-block-column-correspondence-is-unique]`) the derivation
/// no longer finds a rotation to partition around, so the codon-frame exception
/// has an adjacent pair to act on again and acts on it only where a reading frame
/// is declared. The frameless axes split at `c.11`; the coding axis does not.
#[test]
fn an_axis_without_a_reading_frame_is_not_re_split() {
    let coding = provider("NM_TEST.1", CORE, true);
    let noncoding = provider("NR_TEST.1", CORE, false);
    assert_eq!(
        normalized_fixed_point(&coding, "NM_TEST.1:n.10_15delinsTAATAT"),
        "NM_TEST.1:n.[10A>T;12_15delinsATAT]",
    );
    assert_eq!(
        normalized_fixed_point(&noncoding, "NR_TEST.1:n.10_15delinsTAATAT"),
        "NR_TEST.1:n.[10A>T;12_15delinsATAT]",
    );
}

/// The #1034 / #1040 control: a reverse-complement **sub**-run of one contiguous
/// change is still not carved out.
///
/// `c.12_15` is `TATA -> CTAT`, one contiguous mismatch run that is not itself a
/// reverse complement, even though `13_15` (`ATA -> TAT`) inside it is. The fix
/// types spans the split *creates*; it must not create one.
///
/// # #1835 moved this to `c.[11_12insC;15del]`; #1878 returned it
///
/// `TATA -> CTAT` is another equal-length rotation, and the flip read it as one
/// inserted `C` plus one deleted `A`. There is nothing to read: four reference
/// bases against a four-base payload have one correspondence, all four columns
/// differ, and `delins.md:16` makes them one `delins`. The control the row exists
/// for held under both pinnings — the interior `13_15` reverse complement is
/// carved out by neither.
#[test]
fn a_reverse_complement_sub_run_of_one_contiguous_change_stays_a_delins() {
    let provider = provider("NM_TEST.1", CORE, true);
    assert_eq!(
        normalized_fixed_point(&provider, "NM_TEST.1:c.12_15delinsCTAT"),
        "NM_TEST.1:c.12_15delinsCTAT",
    );
}

/// The second call site: the codon-frame triplet itself can be an inversion.
///
/// See [`AMBIGUOUS_CORE`] for why this needs a self-complementary centre — every
/// other base makes the merged three-base member's middle column disagree with
/// its own reverse complement, which rules the span out. Reachable rather than
/// defensive, which is why the branch is typed instead of assumed.
///
/// **The input has to be longer than the triplet itself.** A lone
/// `c.10_12delinsGNT` never reaches `build_split_variants` at all: a whole-span
/// reverse complement is typed `inv` by `normalize_na_edit`'s Delins arm before
/// any split happens, so it exercises the direct path and says nothing about the
/// codon branch. Instrumenting that branch against the three-base spelling logs
/// zero hits. Extending the input puts a second mismatch run after the triplet,
/// which is what makes the split fire — and only then does the codon exception
/// consume `[Sub@10; IdentityAt@11; Sub@12]` and hand `ANC -> GNT` to
/// [`push_typed_replacement`] as a span the split created.
///
/// **And the second run has to be separated from the triplet, not touching it
/// (#1524).** The input was `c.10_15delinsGNTAAA`, whose trailing run starts at
/// `c.13` — consecutive with the triplet's `c.12`, which `delins.md:16` says is
/// one member. That spelling now yields a single `c.10_15delinsGNTAAA` and
/// exercises nothing. Leaving `c.13` unchanged puts one nucleotide between the
/// two, `general.md:34` keeps them apart, and the triplet reaches
/// [`push_typed_replacement`] exactly as before. Without this retarget #1454's
/// only live guard on that call site would have been silently deleted.
#[test]
fn codon_triplet_over_an_ambiguous_centre_is_an_inversion() {
    let provider = provider("NM_TEST.1", AMBIGUOUS_CORE, true);
    // `10_12` is `ANC` -> `GNT`, an exact reverse complement, and it is the
    // codon-frame triplet the exception merges. `14_15` is `GG` -> `AA`, which
    // is not one, so it stays a `delins` — the two outcomes of the same call
    // site in one row.
    assert_eq!(
        normalized_fixed_point(&provider, "NM_TEST.1:c.10_15delinsGNTGAA"),
        "NM_TEST.1:c.[10_12inv;14_15delinsAA]",
    );
    // The three-base spelling is the same variant reached by the direct path,
    // and must agree. Kept as the control that the two routes to `inv` do not
    // disagree, not as coverage of the branch above.
    assert_eq!(
        normalized_fixed_point(&provider, "NM_TEST.1:c.10_12delinsGNT"),
        "NM_TEST.1:c.10_12inv",
    );
}
