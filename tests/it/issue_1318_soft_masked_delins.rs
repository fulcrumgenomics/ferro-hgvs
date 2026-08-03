//! Issue #1318 — the delins path must not depend on the reference's masking.
//!
//! #1250 case-folded the *inversion-shortening* comparisons and was
//! deliberately scope-limited; the delins→inv typing path kept its raw byte
//! tests. On a soft-masked reference `deleted` arrives lowercase while
//! `inserted` comes from the author's upper-case description, so those tests
//! saw a mismatch that is not one.
//!
//! The consequence is worse than a mis-typing. Two comparisons on the same path
//! disagreed about case:
//!
//! - `is_revcomp` complements through `complement()`, which folds to
//!   upper-case, so lower-case deleted vs upper-case inserted compared **equal**;
//! - `shared_affix_lengths` compares raw bytes, so the same pair shared **no**
//!   prefix or suffix and nothing was trimmed.
//!
//! A palindromic delins therefore entered the inversion branch with an
//! untrimmed range, `shorten_inversion` (correctly case-folding since #1250)
//! collapsed it to nothing, and the branch's `expect` — whose stated invariant
//! is that trimming has already handled that case — **panicked**.
//!
//! Reference FASTAs are routinely soft-masked, and soft-masked bytes reach
//! normalization in this codebase (`issue_1235_cis_allele_confluence.rs::
//! soft_masked_reference_yields_the_same_canonical_form`), so this is reachable
//! from ordinary input.

use std::io::Write;

use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};

/// Positions 29..34 are `CACGTG` and 35..40 are `CAATTG` — both palindromes, so
/// both are their own reverse complement.
const SEQ: &str = concat!(
    "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCG",
    "TACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA"
);

fn provider_for(sequence: &str) -> JsonProvider {
    let n = sequence.len() as u64;
    let json = serde_json::json!({
        "version": "1.0",
        "transcripts": [{
            "id": "TEMPLATE-gene.1",
            "chromosome": "TEMPLATE",
            "strand": "+",
            "sequence": sequence,
            "cds_start": 1,
            "cds_end": n - (n % 3),
            "genomic_start": 1,
            "genomic_end": n,
            "protein_id": "TEMPLATE-gene.1",
            "exons": [{"number": 1, "start": 1, "end": n,
                       "genomic_start": 1, "genomic_end": n}],
        }],
        "genomic_sequences": {"TEMPLATE": sequence},
    });
    let mut f = tempfile::NamedTempFile::new().expect("tempfile");
    write!(f, "{json}").expect("write json");
    JsonProvider::from_json(f.path()).expect("load reference")
}

fn normalize_with(sequence: &str, input: &str) -> String {
    let normalizer = Normalizer::new(provider_for(sequence));
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize `{input}`: {e}"))
        .to_string()
}

/// The crash. A palindromic delins over a soft-masked span.
///
/// On an upper-case reference the insert equals the deleted bases and the
/// identity branch answers immediately; on the lower-case one the byte
/// comparison misses that, and the range reaches the inversion branch untrimmed.
#[test]
fn a_palindromic_delins_on_a_soft_masked_reference_does_not_panic() {
    // `CACGTG` at 29..34 and `CAATTG` at 35..40 are both palindromes. Each
    // insert equals the deleted bases modulo case, so the canonical form is the
    // identity — which is the answer the uppercase reference already gave, and
    // the one the soft-masked reference panicked instead of giving.
    for (input, expected) in [
        ("TEMPLATE:g.29_34delinsCACGTG", "TEMPLATE:g.29_34="),
        ("TEMPLATE:g.35_40delinsCAATTG", "TEMPLATE:g.35_40="),
    ] {
        let masked = normalize_with(&SEQ.to_lowercase(), input);
        // Pin the output, not merely parity: two equally-wrong answers would
        // satisfy a parity check, and the whole class of defect here is one
        // comparison folding case while its sibling does not.
        assert_eq!(masked, expected, "canonical form of `{input}`");
        let upper = normalize_with(SEQ, input);
        assert_eq!(
            masked, upper,
            "soft-masking must not change the canonical form of `{input}`"
        );
    }
}

/// The spread of shapes, each reaching a different one of the case-sensitive
/// comparison sites, with its canonical form pinned.
///
/// Two assertions per case, and both are needed. The pinned output is the
/// primary check — it is what a parity-only test cannot give, since masked and
/// unmasked can be equally wrong. The parity assertion is the supplementary
/// one: it says the *reason* the output is right is not "the uppercase path
/// happens to agree today".
///
/// The eight cases land on eight distinct canonical shapes — identity, `inv`,
/// three different single substitutions, a two-member decomposition, and a
/// `dup` — which is how you can tell the table exercises separate branches
/// rather than eight spellings of one.
#[test]
fn soft_masking_never_changes_the_canonical_form_of_a_delins() {
    for (input, expected) in [
        // Palindromes — the inversion branch, which for a palindrome reduces to
        // the identity.
        ("TEMPLATE:g.29_34delinsCACGTG", "TEMPLATE:g.29_34="),
        ("TEMPLATE:g.35_40delinsCAATTG", "TEMPLATE:g.35_40="),
        // A true reverse complement that is not an identity: 63..70 is
        // `TACGTACG`, which is not its own reverse complement, and the insert
        // is that reverse complement. So `is_revcomp` must fire while the
        // identity test above it must not — the discrimination a palindrome
        // cannot show, because for a palindrome the two agree.
        ("TEMPLATE:g.63_70delinsCGTACGTA", "TEMPLATE:g.63_70inv"),
        // Shared suffix, then shared prefix, each differing only in case: the
        // trimming must consume it and leave a bare substitution.
        ("TEMPLATE:g.29_34delinsCACTTG", "TEMPLATE:g.32G>T"),
        ("TEMPLATE:g.29_34delinsCTCGTG", "TEMPLATE:g.30A>T"),
        // A genuine interior identity run flanked by mismatches on both ends:
        // 29..36 is `CACGTGCA`, and the insert changes only the first and last
        // base, leaving `ACGTGC` identical in the middle. That is what the
        // decomposition loop's identity-run detector exists to find, and it is
        // why the answer is two members rather than one delins; an insert equal
        // to the whole deleted span never reaches it, because the whole-range
        // identity test answers first.
        ("TEMPLATE:g.29_36delinsGACGTGCT", "TEMPLATE:g.[29C>G;36A>T]"),
        // A duplication (`N -> 2N`), which `canonicalize_delins` classifies
        // before trimming and with its own pair of comparisons. The 3' shift
        // moves it off the input's own span, so this also pins that the dup was
        // typed as a dup rather than trimmed down to an insertion.
        ("TEMPLATE:g.29_34delinsCACGTGCACGTG", "TEMPLATE:g.31_36dup"),
        // A plain substitution reached through trimming.
        ("TEMPLATE:g.29_31delinsCAG", "TEMPLATE:g.31C>G"),
    ] {
        let masked = normalize_with(&SEQ.to_lowercase(), input);
        assert_eq!(masked, expected, "canonical form of `{input}`");
        let upper = normalize_with(SEQ, input);
        assert_eq!(
            masked, upper,
            "soft-masking must not change the canonical form of `{input}`"
        );
    }
}

/// The mirror-image mixed-case shape: an **r.** description, whose edit
/// literals are lowercase by convention, against an **uppercase** reference.
///
/// Every case above puts the lowercase bytes on the reference side and the
/// uppercase ones in the description. The RNA axis inverts that — `to_rna`
/// lowercases the edit literals and maps `t` to `u` — so a comparison that
/// folds case in only one direction, or that folds case but not the `t`/`u`
/// spelling, is invisible to the whole table above. Both orientations reach the
/// same comparison sites, so both must answer the same shapes.
#[test]
fn an_rna_description_with_lowercase_payloads_takes_the_same_shapes() {
    for (input, expected) in [
        // Identity — the palindrome.
        (
            "TEMPLATE-gene.1:r.29_34delinscacgug",
            "TEMPLATE-gene.1:r.29_34=",
        ),
        // Inversion — the non-palindromic reverse complement.
        (
            "TEMPLATE-gene.1:r.63_70delinscguacgua",
            "TEMPLATE-gene.1:r.63_70inv",
        ),
        // The interior identity run, decomposing into two substitutions.
        (
            "TEMPLATE-gene.1:r.29_36delinsgacgugcu",
            "TEMPLATE-gene.1:r.[29c>g;36a>u]",
        ),
    ] {
        let upper = normalize_with(SEQ, input);
        assert_eq!(upper, expected, "canonical form of `{input}`");
        let masked = normalize_with(&SEQ.to_lowercase(), input);
        assert_eq!(
            masked, upper,
            "soft-masking must not change the canonical form of `{input}`"
        );
    }
}

/// Normalization must be a fixed point on a soft-masked reference, not merely
/// agree with the uppercase run once.
///
/// The parity tests above compare two *first* passes. They would still hold if
/// masked and unmasked were equally non-idempotent, which is the failure this
/// class actually produces: a comparison that folds case in one place and not
/// another can type the first pass one way and the re-read output another. So
/// this asserts `norm(norm(x)) == norm(x)` against a lowercase reference with
/// uppercase edit literals — the real mixed-case shape, since HGVS descriptions
/// are written uppercase and reference FASTAs are routinely soft-masked.
///
/// It also re-parses the intermediate, so an unreadable output fails here
/// rather than surfacing as the second pass's parse error.
#[test]
fn normalization_is_idempotent_on_a_soft_masked_reference() {
    let masked_seq = SEQ.to_lowercase();
    for input in [
        "TEMPLATE:g.29_34delinsCACGTG",
        "TEMPLATE:g.35_40delinsCAATTG",
        "TEMPLATE:g.63_70delinsCGTACGTA",
        "TEMPLATE:g.29_34delinsCACTTG",
        "TEMPLATE:g.29_34delinsCTCGTG",
        "TEMPLATE:g.29_36delinsGACGTGCT",
        "TEMPLATE:g.29_34delinsCACGTGCACGTG",
        "TEMPLATE:g.29_31delinsCAG",
        // The inverse orientation: lowercase `r.` payloads (with `u` for `t`)
        // against the same soft-masked reference.
        "TEMPLATE-gene.1:r.29_34delinscacgug",
        "TEMPLATE-gene.1:r.63_70delinscguacgua",
        "TEMPLATE-gene.1:r.29_36delinsgacgugcu",
    ] {
        let once = normalize_with(&masked_seq, input);
        parse_hgvs(&once)
            .unwrap_or_else(|e| panic!("`{input}` -> `{once}` does not re-parse: {e}"));
        let twice = normalize_with(&masked_seq, &once);
        assert_eq!(
            twice, once,
            "`{input}` -> `{once}` is not a fixed point on a soft-masked reference"
        );
    }
}
