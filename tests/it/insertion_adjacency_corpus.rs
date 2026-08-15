//! What a cis allele containing an insertion **next to** another member does.
//!
//! # The geometry this module is about
//!
//! An insertion occupies a zero-width junction, not the two positions it names.
//! The spec says so directly — `docs/syntax.yaml:11`: "Ranges are inclusive for
//! all variant types except insertions, which for which ranges are exclusive"
//! (sic) — and every insertion example in `DNA/insertion.md` is glossed "between
//! nucleotides". The `^` spelling (`g.123^124insG`) is rejected at
//! `DNA/insertion.md:105` purely on character economy ("a character was already
//! used to indicate a range"), not because junction semantics are wrong.
//!
//! So four arrangements have to be told apart, and conflating any two of them
//! is the defect class this module guards:
//!
//! | arrangement | example | verdict |
//! |---|---|---|
//! | junction at a span's **edge** | `[301_302insG;302T>C]` | accepted, merged |
//! | junction **strictly interior** to a span that keeps its bases | `[302_304delinsGG;302_303insG]` | refused, W5002 |
//! | junction **shared** with a span edit that writes there (`dup`) | `[302_304dup;304_305insG]` | accepted, composed |
//! | junction **shared** with a second true insertion | `[302_303insA;302_303insT]` | refused, W5002 |
//!
//! A reader who models an insertion as covering `[a, a+1]` cannot separate row 1
//! from row 2, and would refuse the spec-legal shape in row 1.
//!
//! **Row 3 was itself a refusal until this change**, and it is the one the
//! `delins-adjacent-members-when-both-consume-reference` ruling reverses — see
//! [`insertion_at_a_duplications_write_junction_composes_with_it`], which pins
//! the acceptance, and the doc comment there for the clause that carries it.
//! Rows 3 and 4 are what the narrowing separates: one of each kind is ordered by
//! `duplication.md:91`'s gloss, two of the same kind are not.
//!
//! # Why these are pinned rather than left to the sweeps
//!
//! `cis_junction_crossing_shift.rs` and friends sweep this space and assert that
//! normalization preserves the sequence. That is necessary and not sufficient:
//! it cannot see a *verdict* change (an accepted shape starting to be refused
//! costs no sequence), and it cannot see a *form* change (two spellings of one
//! sequence both preserve it). Both are what regressed historically — #1286,
//! #1290, #1292, #1301, #1316 are all this family — so the verdict and the form
//! are pinned here explicitly.
//!
//! # Status of the refusals
//!
//! Rows 2 and 4 are **ferro policy, not spec compliance**, and are labelled as
//! such at each assertion. The spec defines member conflict in exactly one
//! clause, `general.md:58`, whose stated ground ("removing part of a reference
//! sequence and replacing it with part of the same sequence") describes only its
//! own `del`+`dup` example. Nothing in the corpus addresses an insertion
//! interior to a span, and nothing addresses two *true insertions* at one
//! junction — for the latter the spec instead supplies a different designated
//! spelling, `ins[A;B]` (`general.md:78`, worked at `:81`). Refusing is
//! defensible by construction and it agrees with Mutalyzer's `EOVERLAP`; it is
//! not something the spec requires.
//!
//! Rows 1 and 3, by contrast, **are** answered by the spec — row 1 by
//! `delins.md:86-89` (see
//! `recommended_form_pins::substitution_and_abutting_insertion_merge_per_delins_md_86`)
//! and row 3 by `duplication.md:90`, which publishes the shared-junction
//! `[dup;ins]` geometry as a correct description.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer};

use crate::common::cis_apply_oracle::apply_with;
use crate::common::hg38_window::{base_at, local_desc, provider, HG38_WINDOW};

/// Normalize a local-coordinate description in strict mode.
///
/// Strict, not the library default: `NormalizeConfig::default()` is **lenient**
/// (`config.rs:129`, "for backwards compatibility"), under which every
/// `should_reject_*` rung is a no-op and the refusal rows below would silently
/// pass by being accepted. The CLI's own default is strict, so this is also the
/// configuration a caller actually meets.
fn normalize_strict(body: &str) -> Result<String, String> {
    let input = local_desc(body);
    let variant: HgvsVariant = parse_hgvs(&input).map_err(|e| format!("parse: {e}"))?;
    let normalizer = Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    );
    normalizer
        .normalize(&variant)
        .map(|v| {
            v.to_string()
                .strip_prefix(&format!("{}:g.", crate::common::hg38_window::LOCAL_CONTIG))
                .expect("rendered description keeps its accession and axis")
                .to_string()
        })
        .map_err(|e| format!("{e}"))
}

/// Assert `body` normalizes to `expected` **and** that the two denote the same
/// sequence.
///
/// The second half is the load-bearing one and it does not go through the
/// normalizer: `apply_with` converts each member to its SPDI triple and splices
/// the reference itself, so "the output means what the input meant" is not being
/// checked by the component under test. Pinning the string alone would freeze
/// whatever ferro happens to emit, including a wrong answer.
#[track_caller]
fn assert_normalizes(body: &str, expected: &str) {
    let got = normalize_strict(body).unwrap_or_else(|e| panic!("{body} was refused: {e}"));
    assert_eq!(got, expected, "form drifted for {body}");

    let p = provider();
    let before = apply_with(&p, HG38_WINDOW, &local_desc(body))
        .unwrap_or_else(|| panic!("input {body} denotes no sequence"));
    let after = apply_with(&p, HG38_WINDOW, &local_desc(expected))
        .unwrap_or_else(|| panic!("output {expected} denotes no sequence"));
    assert_eq!(
        before, after,
        "normalizing {body} changed the sequence it denotes"
    );
}

/// Assert `body` is refused, and say which of the two policy rows it is.
#[track_caller]
fn assert_refused_by_policy(body: &str, because: &str) {
    match normalize_strict(body) {
        Ok(out) => panic!("{body} was accepted as {out}; expected refusal ({because})"),
        Err(e) => assert!(
            e.contains("coincident cis-allele edits") || e.contains("W5002"),
            "{body} was refused, but not as an overlap conflict: {e}"
        ),
    }
}

// ---------------------------------------------------------------------------
// Row 1 — a junction at a span's edge is not interior, and merges.
// ---------------------------------------------------------------------------

/// An insertion immediately 5' of a substitution.
#[test]
fn insertion_abutting_a_substitution_on_its_five_prime_side_merges() {
    assert_eq!(base_at(302), b'T', "case assumes the reference base at 302");
    assert_normalizes("[301_302insG;302T>C]", "302delinsGC");
}

/// An insertion immediately 3' of a substitution.
///
/// The mirror of the case above, and the one the spec rules on directly. The
/// payload order differs (`GC` against `CG`) because the insertion lands on the
/// other side of the substituted base — a merge that ignored the junction's side
/// would render these identically.
#[test]
fn insertion_abutting_a_substitution_on_its_three_prime_side_merges() {
    assert_normalizes("[302T>C;302_303insG]", "302delinsCG");
}

/// An insertion immediately 5' of a deletion.
#[test]
fn insertion_abutting_a_deletion_on_its_five_prime_side_merges() {
    assert_normalizes("[301_302insG;302_304del]", "302_304delinsG");
}

/// An insertion immediately 3' of a deletion.
///
/// Converges with the case above: both denote "replace `302_304` with `G`",
/// and the junction side is not recoverable from the result, so one form is
/// correct for both.
#[test]
fn insertion_abutting_a_deletion_on_its_three_prime_side_merges() {
    assert_normalizes("[302_304del;304_305insG]", "302_304delinsG");
}

/// An insertion immediately 5' of an inversion.
#[test]
fn insertion_abutting_an_inversion_on_its_five_prime_side_is_accepted() {
    assert_normalizes("[301_302insG;302_306inv]", "[302_303delinsGT;305_306insG]");
}

/// An insertion immediately 3' of an inversion.
///
/// The inversion trims to `303_305`: `302` is `T` and `306` is `A`, so
/// reverse-complementing `302..306` leaves both ends unchanged and the
/// unchanged flanks are not part of the change (`DNA/complex.md:50`, the 3'
/// rule as "maintaining the longest unchanged sequence").
#[test]
fn insertion_abutting_an_inversion_on_its_three_prime_side_is_accepted() {
    assert_eq!((base_at(302), base_at(306)), (b'T', b'A'));
    assert_normalizes("[302_306inv;306_307insG]", "[303_305inv;306_307insG]");
}

/// An insertion immediately 5' of a duplication.
///
/// Not the same as row 3: a `dup` writes at the junction *after* its span, so an
/// insertion before its span collides with nothing.
#[test]
fn insertion_before_a_duplications_span_is_accepted() {
    assert_normalizes("[301_302insG;302_304dup]", "301_302insGTCT");
}

/// An insertion abutting a `delins`, on each side.
#[test]
fn insertion_abutting_a_delins_merges_from_either_side() {
    assert_normalizes("[301_302insG;302_304delinsGG]", "302_304delinsGGG");
    assert_normalizes("[302_304delinsGG;304_305insG]", "302_304delinsGGG");
}

/// Two insertions at *adjacent* junctions — distinct junctions, so not row 3.
#[test]
fn two_insertions_at_adjacent_junctions_are_accepted() {
    assert_normalizes("[301_302insG;302_303insC]", "[301dup;303dup]");
}

/// Three members, an insertion abutting the middle one.
///
/// Member count is not incidental: the junction detector runs over the whole
/// pre-merge member list, so a two-member-only guard would not cover it.
#[test]
fn a_three_member_allele_with_an_abutting_insertion_is_accepted() {
    assert_normalizes(
        "[301_302insG;302_306inv;308del]",
        "[302_303delinsGT;306_308inv]",
    );
}

// ---------------------------------------------------------------------------
// Rows 2 and 3 — POLICY refusals. See the module header: the spec is silent.
// ---------------------------------------------------------------------------

/// A junction strictly interior to a `delins`, which keeps its payload.
#[test]
fn policy_insertion_interior_to_a_delins_is_refused() {
    assert_refused_by_policy(
        "[302_304delinsGG;302_303insG]",
        "the payload survives, so the insertion's position within it is undetermined",
    );
}

/// A junction strictly interior to an `inv`, which keeps its bases.
#[test]
fn policy_insertion_interior_to_an_inversion_is_refused() {
    assert_refused_by_policy(
        "[302_306inv;303_304insG]",
        "the inverted bases survive, so the interior junction stays meaningful",
    );
}

/// A junction shared with a `dup`'s write is **accepted**, and composes.
///
/// This used to be row 3 of the table above, refused alongside two insertions at
/// one junction. It is not the same shape: `DNA/duplication.md:90` publishes
/// `NC_000001.11(NM_206933.2):c.[675-542_1211-703dup;1211-703_1211-702insGTAAA]`
/// as a correct description — an insertion at the junction the duplication
/// writes into — and glosses it as a duplication "**followed by**" the
/// insertion. That gloss supplies the order the refusal said was missing.
///
/// Keeping the pair split is also what satisfies the *other* clause in play:
/// `duplication.md:18` requires a variant describable as a duplication to be
/// described as one, which is why the ledger's
/// `delins-adjacent-members-when-both-consume-reference` declined to merge such
/// a pair into one `delins`. Accepting the split pair is the only answer that
/// neither destroys the duplication nor refuses a shape the spec publishes.
///
/// What is asserted is **convergence**, not a pinned string: the composed pair
/// must reach the same form as the single-member spelling of the same edit. A
/// pinned string would pass even if the two settled apart, which is the failure
/// this whole family keeps producing.
#[test]
fn insertion_at_a_duplications_write_junction_composes_with_it() {
    assert_eq!(
        (base_at(302), base_at(303), base_at(304)),
        (b'T', b'C', b'T'),
        "case assumes g.302_304 is TCT"
    );
    let pair = normalize_strict("[302_304dup;304_305insG]")
        .expect("duplication.md:90 publishes this geometry as correct");
    // The dup copies `TCT` into the junction after 304, then the insertion adds
    // `G` — so the single-member spelling of the same edit is `insTCTG`.
    let single = normalize_strict("304_305insTCTG").expect("the single-member spelling normalizes");
    assert_eq!(
        pair, single,
        "the composed pair settled apart from its single-member spelling"
    );
}

/// …and the composition does not depend on which member was written first.
#[test]
fn a_duplication_and_its_junction_insertion_compose_in_one_order() {
    assert_eq!(
        normalize_strict("[302_304dup;304_305insG]"),
        normalize_strict("[304_305insG;302_304dup]"),
        "member order decided the composed payload"
    );
}

/// Two insertions at one junction stay refused.
///
/// The narrowing above is deliberately not a blanket. Two payloads into one slot
/// really is undetermined — `[4_5insA;4_5insT]` is `…AT…` or `…TA…` — and
/// `general.md:78` supplies a different designated spelling for that content
/// (`ins[A;B]`), so refusing the two-member form is supported by construction.
#[test]
fn policy_two_insertions_at_one_junction_are_still_refused() {
    assert_refused_by_policy(
        "[302_303insA;302_303insT]",
        "two payloads into one slot with nothing to order them",
    );
}

/// A deletion is the documented exception to row 2 (#1406): it removes every
/// base it spans, so an interior junction has nothing left to be positioned
/// against and the composition is unique.
///
/// # Why this one checks the form only
///
/// `apply_with` declines any description whose members claim the same base —
/// that conservatism is what makes it a usable oracle elsewhere. Here the
/// insertion's junction *is* interior to the deletion, so the oracle returns
/// `None` for the input and cannot be used to check sequence preservation.
///
/// That is not a gap in the oracle; it is the whole reason #1406 needed a
/// ruling rather than a measurement. A purely geometric reader cannot tell that
/// this particular overlap composes uniquely — it takes the argument that a
/// deletion leaves nothing for the junction to be positioned against. So the
/// unique composition is asserted here by hand instead: `302_304` is `TCT`, and
/// whichever order the members apply in, the result is `G` in place of those
/// three bases.
#[test]
fn insertion_interior_to_a_deletion_is_accepted_because_the_bases_are_gone() {
    assert_eq!(
        (base_at(302), base_at(303), base_at(304)),
        (b'T', b'C', b'T'),
        "case assumes g.302_304 is TCT"
    );
    let got = normalize_strict("[302_304del;302_303insG]")
        .expect("#1406 — an interior junction in a deletion composes uniquely");
    assert_eq!(got, "302_304delinsG");
}

// ---------------------------------------------------------------------------
// Order independence — the verdict and the form are properties of the allele.
// ---------------------------------------------------------------------------

/// Permuting members changes neither the result nor the verdict.
///
/// `#1261` and `#1103` are both order-leak defects in this exact family, and the
/// authored order is what a conflict-preserving path deliberately keeps, so the
/// two interact.
#[test]
fn member_order_does_not_change_the_normalized_form() {
    for (a, b) in [
        ("[301_302insG;302_303insC]", "[302_303insC;301_302insG]"),
        ("[301_302insG;302_306inv]", "[302_306inv;301_302insG]"),
        ("[301_302insG;302T>C]", "[302T>C;301_302insG]"),
    ] {
        assert_eq!(
            normalize_strict(a),
            normalize_strict(b),
            "member order leaked into the result for {a}"
        );
    }
}

/// Every accepted output re-parses.
///
/// A normalizer that emits a string its own parser rejects has produced
/// something no consumer can read back; this is the cheap suite-wide check for
/// it, and it is the shape that catches the CDS/3'UTR seam defect pinned in
/// `insertion_adjacency_defects`.
#[test]
fn every_accepted_output_reparses() {
    let bodies = [
        "[301_302insG;302T>C]",
        "[302T>C;302_303insG]",
        "[301_302insG;302_304del]",
        "[302_304del;304_305insG]",
        "[301_302insG;302_306inv]",
        "[302_306inv;306_307insG]",
        "[301_302insG;302_304dup]",
        "[301_302insG;302_304delinsGG]",
        "[302_304delinsGG;304_305insG]",
        "[301_302insG;302_303insC]",
        "[301_302insG;302_306inv;308del]",
        "[302_304del;302_303insG]",
    ];
    for body in bodies {
        let out = normalize_strict(body).unwrap_or_else(|e| panic!("{body}: {e}"));
        let rendered = local_desc(&out);
        parse_hgvs(&rendered)
            .unwrap_or_else(|e| panic!("{body} emitted {rendered}, which does not re-parse: {e}"));
    }
}
