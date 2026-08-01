//! Issue #1249: an `inv` whose complementary outer bases all cancel down to a
//! *single* remaining base must become a substitution, not identity.
//!
//! # The bug
//!
//! `rules::shorten_inversion` peels complementary outer pairs off an inversion.
//! That peeling is sound: when the first base is the complement of the last,
//! both endpoints survive the reverse complement unchanged, so neither needs
//! describing. What was unsound is the stopping condition — it treated a
//! residue of **zero or one** base as identity:
//!
//! ```text
//! reference   CCAATGGCC
//! g.3_5inv    span AAT, revcomp ATT   ->  ferro emitted g.3_5=
//! applied     CCATTGGCC                   the sequence DID change
//! ```
//!
//! A zero-base residue really is identity: every base cancelled. A **one**-base
//! residue is not — that base is replaced by its own complement, and no
//! standard base is self-complementary. `AAT` peels its `A..T` pair and leaves
//! `A` at position 4, which inverts to `T`, so the variant is `g.4A>T`.
//!
//! # Spec
//!
//! `DNA/inversion.md:5` defines an inversion as **more than one nucleotide**,
//! and `:16` gives the replacement rule outright:
//!
//! > by definition, the region inverted (`positions_inverted`) contains **more
//! > than one nucleotide**. The description `g.234inv` is therefore not
//! > allowed; a one-nucleotide inversion should be described as a
//! > [substitution](substitution.md).
//!
//! So the residue is emitted as a substitution to its complement — never as a
//! one-base `inv`, and never as `=`.
//!
//! This is the normalizer-side counterpart to #1079, which enforces the same
//! rule in the *parser*. #1079 can only reject `g.234inv`, because recovering
//! the substitution needs the reference sequence and the parser has none. The
//! normalizer has the reference, so it performs the rewrite #1079 documents.
//!
//! # Why it matters twice
//!
//! Beyond losing the variant outright — `=` asserts the reference is unchanged,
//! so any consumer filtering identities drops it silently — it is also a
//! **confluence** break of the kind #1229–#1235 track: the `delins` and
//! substitution spellings of this very variant normalize correctly, and only
//! the `inv` spelling is lost.

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder, PAD_OFFSET};

/// 1-based HGVS position of the `n`-th (1-based) core base on a genomic build.
fn at(n: u64) -> u64 {
    PAD_OFFSET + n
}

/// Normalize `input` against a genomic provider carrying `core`.
fn norm(core: &str, input: &str) -> String {
    normalize_to_string(SyntheticBuilder::genomic(core).build(), input)
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        other => other,
    }
}

fn revcomp(seq: &str) -> String {
    seq.bytes().rev().map(|b| complement(b) as char).collect()
}

/// The sequence that `g.<s>_<e>inv` produces, computed independently of the
/// normalizer: splice the reverse complement of the span back into `seq`.
/// `s`/`e` are 1-based inclusive positions into `seq`.
fn apply_inv(seq: &str, s: usize, e: usize) -> String {
    format!("{}{}{}", &seq[..s - 1], revcomp(&seq[s - 1..e]), &seq[e..])
}

// ---------------------------------------------------------------------------
// The reported case.
// ---------------------------------------------------------------------------

#[test]
fn one_base_residue_becomes_a_substitution_to_its_complement() {
    // CCAATGGCC: the 3_5 span is AAT. The outer A..T pair is complementary and
    // cancels, leaving `A` at core position 4, which inverts to `T`.
    let input = format!("NC_TEST.1:g.{}_{}inv", at(3), at(5));
    assert_eq!(
        norm("CCAATGGCC", &input),
        format!("NC_TEST.1:g.{}A>T", at(4)),
        "an inv leaving a one-base residue is a substitution, not identity"
    );
}

#[test]
fn an_explicit_inverted_sequence_takes_the_same_path() {
    // The `invAAT` spelling must not escape the rewrite.
    let input = format!("NC_TEST.1:g.{}_{}invAAT", at(3), at(5));
    assert_eq!(
        norm("CCAATGGCC", &input),
        format!("NC_TEST.1:g.{}A>T", at(4))
    );
}

// ---------------------------------------------------------------------------
// The even-length counterpart must stay identity.
// ---------------------------------------------------------------------------

#[test]
fn a_fully_cancelling_even_span_is_still_identity() {
    // TTGGCCAA: the 3_6 span is GGCC. Both outer pairs (G..C, G..C) cancel and
    // nothing remains, so this genuinely denotes no change. The fix must not
    // turn these into substitutions.
    let input = format!("NC_TEST.1:g.{}_{}inv", at(3), at(6));
    assert_eq!(
        norm("TTGGCCAA", &input),
        format!("NC_TEST.1:g.{}_{}=", at(3), at(6))
    );
}

// ---------------------------------------------------------------------------
// Confluence: every spelling of this variant must agree (#1229–#1235 family).
// ---------------------------------------------------------------------------

#[test]
fn every_spelling_of_the_variant_normalizes_alike() {
    let core = "CCAATGGCC";
    let expected = format!("NC_TEST.1:g.{}A>T", at(4));
    let spellings = [
        format!("NC_TEST.1:g.{}_{}inv", at(3), at(5)),
        format!("NC_TEST.1:g.{}_{}invAAT", at(3), at(5)),
        format!("NC_TEST.1:g.{}_{}delinsATT", at(3), at(5)),
        format!("NC_TEST.1:g.{}_{}delAATinsATT", at(3), at(5)),
        format!("NC_TEST.1:g.{}A>T", at(4)),
    ];
    for spelling in &spellings {
        assert_eq!(
            norm(core, spelling),
            expected,
            "spelling `{spelling}` must converge with the others"
        );
    }
}

// ---------------------------------------------------------------------------
// Exhaustive parity sweep — the invariant behind the whole fix.
// ---------------------------------------------------------------------------

/// No inversion may normalize to `=` unless applying it really leaves the
/// reference unchanged, and one that does leave it unchanged must still say so.
///
/// Sweeping every 4^L core for L in 3..=6 covers both parities: an even-length
/// span whose outer pairs all cancel reduces to nothing and is identity, while
/// an odd-length one always retains a centre base that changes. Before the fix
/// this failed on 16 cores at L=3 and 64 at L=5 (`4^((L+1)/2)`), i.e. on every
/// odd-length span that cancelled completely — there were no false positives to
/// weigh against, so the old branch was wrong in every case it fired.
#[test]
fn identity_is_emitted_only_when_the_sequence_is_actually_unchanged() {
    const BASES: [u8; 4] = *b"ACGT";

    for len in 3..=6usize {
        for n in 0..4usize.pow(len as u32) {
            let mut k = n;
            let core: String = (0..len)
                .map(|_| {
                    let b = BASES[k % 4] as char;
                    k /= 4;
                    b
                })
                .collect();

            let input = format!("NC_TEST.1:g.{}_{}inv", at(1), at(len as u64));
            let out = norm(&core, &input);
            let emitted_identity = out.ends_with('=');
            let unchanged = apply_inv(&core, 1, len) == core;

            assert_eq!(
                emitted_identity,
                unchanged,
                "core {core}: normalize gave `{out}`, but applying the inversion \
                 {} the sequence",
                if unchanged { "left" } else { "changed" }
            );
        }
    }
}

// ---------------------------------------------------------------------------
// The rewrite is axis-independent: it lives in the shared `NaEdit::Inversion`
// arm, so c. and n. must behave identically to g.
// ---------------------------------------------------------------------------

#[test]
fn the_rewrite_applies_on_the_cds_axis() {
    use ferro_hgvs::reference::transcript::Strand;
    // Whole transcript is the reference for c., so positions are literal.
    let provider = SyntheticBuilder::cds("CCAATGGCC", 1, 9, Strand::Plus).build();
    assert_eq!(
        normalize_to_string(provider, "NM_TEST.1:c.3_5inv"),
        "NM_TEST.1:c.4A>T"
    );
}

#[test]
fn the_rewrite_applies_on_the_noncoding_axis() {
    use ferro_hgvs::reference::transcript::Strand;
    let provider = SyntheticBuilder::noncoding("CCAATGGCC", Strand::Plus).build();
    assert_eq!(
        normalize_to_string(provider, "NR_TEST.1:n.3_5inv"),
        "NR_TEST.1:n.4A>T"
    );
}

// ---------------------------------------------------------------------------
// IUPAC ambiguity codes in the reference.
//
// The peeling test is `complement(first) == last`, and `complement` originally
// modelled only A/C/G/T/U — every other symbol mapped to itself. That made an
// ambiguity code look self-complementary, which broke the same two ways this
// issue is about:
//
//   * `R` (A/G) inverts to `Y` (C/T), but looked like identity, so an inversion
//     resolving to a lone `R` was silently discarded — #1249 again.
//   * a lone centre base cancelled against *itself* (`ref_seq[s]` and
//     `ref_seq[e - 1]` are the same byte once the span is one base wide),
//     driving the shortened end below the shortened start.
// ---------------------------------------------------------------------------

#[test]
fn an_ambiguous_one_base_residue_becomes_a_substitution() {
    // ART: A cancels against T, leaving R. revcomp(R) == Y, so this is a real
    // change and must be described, not dropped as `=`.
    assert_eq!(
        norm("ART", &format!("NC_TEST.1:g.{}_{}inv", at(1), at(3))),
        format!("NC_TEST.1:g.{}R>Y", at(2))
    );
}

#[test]
fn a_self_complementary_residue_is_identity() {
    // W (A/T), S (C/G) and N are the only self-complementary symbols: each is
    // its own reverse complement, so the span really is unchanged.
    for centre in ["W", "S", "N"] {
        let core = format!("A{centre}T");
        assert_eq!(
            norm(&core, &format!("NC_TEST.1:g.{}_{}inv", at(1), at(3))),
            format!("NC_TEST.1:g.{}_{}=", at(1), at(3)),
            "{core}: a self-complementary centre leaves the sequence unchanged"
        );
    }
}

#[test]
fn an_ambiguous_pair_cancels_only_when_it_actually_complements() {
    // revcomp("RY") == "RY" -- a genuine no-op.
    assert_eq!(
        norm("RY", &format!("NC_TEST.1:g.{}_{}inv", at(1), at(2))),
        format!("NC_TEST.1:g.{}_{}=", at(1), at(2))
    );
    // revcomp("RR") == "YY" -- nothing cancels, so the inversion stands.
    assert_eq!(
        norm("RR", &format!("NC_TEST.1:g.{}_{}inv", at(1), at(2))),
        format!("NC_TEST.1:g.{}_{}inv", at(1), at(2))
    );
}
