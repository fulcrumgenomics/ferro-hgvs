//! #1517 — a whole-block reverse complement is an `inv` when the description it
//! competes with is made of `delins` members, and is **not** when that
//! description is made of substitutions.
//!
//! Both directions come from one clause, `general.md:56`:
//!
//! > when a description is possible according to several types, the preferred
//! > description is: (1) substitution, (2) deletion, (3) inversion,
//! > (4) duplication, (5) insertion.
//!
//! `delins` is absent from that list. `DNA/delins.md` defines it as the residual
//! form — "not a substitution, inversion or conversion" — so it ranks below
//! every named type, inversion included.
//!
//! That single rule settles two cases the project previously reasoned about
//! separately, and which pull in **opposite** directions. They are written here
//! side by side, on one reference, because the contrast is the whole point: the
//! blocks differ only in the *type* of the members the inversion competes with.
//!
//! | filed as | block | competing members | priority | preferred |
//! |---|---|---|---|---|
//! | #1230 | `GATG` -> `CATC` | two substitutions | 1 beats 3 | the split |
//! | #179 / #1517 | `AATGCACA` -> `TGTGCATT` | two `delins` | 3 beats residual | `inv` |
//! | #1706 | `GACAAGTG` -> `CACTTGTC` | one sub **and** one `delins` | `:56` ranks neither mixture | `inv` |
//!
//! The third row is a 2026-08-12 operator ruling and is the case the two
//! candidate gates disagreed on — see
//! `a_reverse_complement_competing_with_mixed_members_is_an_inversion`, which
//! states both readings and why the wider one was taken.
//!
//! #1230 and #179 were filed as *contradictory* bug reports — #1230 says an
//! inversion with an unchanged interior is over-recognized and should be
//! separate substitutions; #179 says an inversion whose interior is
//! self-complementary is under-recognized and should be `inv`. Read through
//! `general.md:56` they are consistent, and the apparent contradiction is an
//! artefact of neither report naming the type of the alternative it was
//! comparing against.
//!
//! # What the spec does not settle
//!
//! `general.md:110` mentions `c.76_83inv` only to illustrate the `inv` keyword;
//! it is not a worked example with a stated normalized output, and the spec has
//! no worked example of an inversion whose interior coincides. Applying a
//! *type* priority to rank descriptions of differing **arity** — one member
//! against two — is an implementer's reading of `general.md:56`, not something
//! its text compels. The counter-reading is `general.md:34`: two changes
//! separated by unchanged nucleotides are two variants, described individually.
//! This file pins the project's choice, and the reasoning for it, so that a
//! future reader can see it was a choice.

use ferro_hgvs::reference::mock::JsonProvider;
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, Normalizer};
use std::io::Write;

/// Two inversion blocks separated by non-repetitive flanks, so neither can
/// 3'-shift into the other and neither sits in a tract that would let a member
/// rotate:
///
/// ```text
///  1-10  GGTTCCAAGT   flank
/// 11-14  GATG         #1230's block; reverse complement CATC, interior `AT` unchanged
/// 15-24  GGTTCCAAGT   flank
/// 25-32  AATGCACA     #1517's block; reverse complement TGTGCATT, 4 of 8 columns unchanged
/// 33-42  GGTTCCAAGT   flank
/// 43-50  GACAAGTG     the MIXED block; reverse complement CACTTGTC (see below)
/// 51-60  GGTTCCAAGT   flank
/// ```
const SEQ: &str = "GGTTCCAAGTGATGGGTTCCAAGTAATGCACAGGTTCCAAGTGACAAGTGGGTTCCAAGT";

fn provider() -> JsonProvider {
    let n = SEQ.len() as u64;
    let json = serde_json::json!({
        "version": "1.0",
        "transcripts": [{
            "id": "TEMPLATE-gene.1",
            "chromosome": "TEMPLATE",
            "strand": "+",
            "sequence": SEQ,
            "cds_start": 1,
            "cds_end": n - (n % 3),
            "genomic_start": 1,
            "genomic_end": n,
            "protein_id": "TEMPLATE-gene.1",
            "exons": [{
                "number": 1, "start": 1, "end": n,
                "genomic_start": 1, "genomic_end": n
            }]
        }],
        "genomic_sequences": { "TEMPLATE": SEQ }
    });
    let mut f = tempfile::NamedTempFile::new().expect("tempfile");
    write!(f, "{json}").expect("write json");
    JsonProvider::from_json(f.path()).expect("load reference")
}

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(provider());
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"));
    format!("{normalized}")
}

/// Assert every encoding normalizes to `canonical`, and that `canonical` is
/// itself a fixed point — the second half is what makes this a confluence check
/// rather than a table of expected strings.
fn converges_to(canonical: &str, encodings: &[&str]) {
    for input in encodings {
        assert_eq!(
            normalize_to_string(input),
            canonical,
            "`{input}` should normalize to `{canonical}`"
        );
    }
    assert_eq!(
        normalize_to_string(canonical),
        canonical,
        "`{canonical}` must be a fixed point"
    );
}

/// The competing members are **`delins`**, so inversion wins (`general.md:56`:
/// inversion is ranked, `delins` is the unranked residual).
///
/// `AATGCACA` -> `TGTGCATT` is a true reverse complement that coincides at 4 of
/// its 8 columns, leaving two 2-base changes six bases apart. Described
/// individually those are `delins` members, which is the shape #179 reported as
/// an inversion the normalizer fails to detect.
#[test]
fn a_reverse_complement_competing_with_delins_members_is_an_inversion() {
    converges_to(
        "TEMPLATE:g.25_32inv",
        &[
            "TEMPLATE:g.25_32inv",
            "TEMPLATE:g.25_32delinsTGTGCATT",
            "TEMPLATE:g.[25_26delinsTG;31_32delinsTT]",
        ],
    );
}

/// The competing members are **substitutions**, so the split wins
/// (`general.md:56` ranks substitution first, above inversion).
///
/// This is #1230's case, and it must keep its answer: the rule introduced for
/// #1517 ranks by member *type*, so it must not disturb a block whose
/// alternative outranks inversion. Without this control the new rule could be
/// written as "always prefer `inv`" and still look correct.
#[test]
fn a_reverse_complement_competing_with_substitutions_stays_split() {
    converges_to(
        "TEMPLATE:g.[11G>C;14G>C]",
        &[
            "TEMPLATE:g.11_14inv",
            "TEMPLATE:g.11_14delinsCATC",
            "TEMPLATE:g.[11G>C;14G>C]",
        ],
    );
}

/// The competing members are **mixed** — one substitution and one `delins` — and
/// the ruling is that a mixed competitor is not a substitution competitor, so
/// inversion wins.
///
/// # This is the case the two candidate gates disagree on
///
/// Operator ruling, 2026-08-12. Two readings of `general.md:56` were live at
/// once and no test discriminated them, which meant merge order would have
/// settled the question silently:
///
/// | gate | admits an `inv` when | this block |
/// |---|---|---|
/// | narrow | **no** competing piece is a substitution | refused |
/// | wide | **not every** competing piece is a substitution | **admitted** |
///
/// The wide reading is the one this project takes, and this test is that ruling.
/// The ground is that `:56` ranks *substitution* above inversion — it says
/// nothing about a description that is only partly substitutions, and a
/// competitor holding a `delins` is not an alternative `:56` can rank, which is
/// exactly `inversion-vs-two-delins-76-83`'s holding. Under the narrow reading a
/// single coincidental substitution column vetoes the typing of the whole span,
/// which makes the notation turn on a base coincidence rather than on the event.
///
/// # The block is built so that no other route can admit it
///
/// `GACAAGTG` -> `CACTTGTC` is an exact reverse complement. Its mirror pairs
/// `(2,7)` and `(3,6)` are self-complementary and its pairs `(1,8)` and `(4,5)`
/// are not, so the changed columns are `{1, 4, 5, 8}` — three runs, `[1]`,
/// `[4,5]`, `[8]`, i.e. a substitution, a two-base `delins`, a substitution.
/// Every other admission route is deliberately starved:
///
/// - separations are **2** and **2**, so the single-base route cannot fire;
/// - changed columns are **4 of 8**, exactly half, so the density route reads it
///   as not dominant;
/// - payload columns are **4 of 8** too, so the payload route reads the same.
///
/// So this block is decided by the member-type route alone, and by which form of
/// it is in force. If a future change makes another route admit this block, this
/// test stops discriminating the two gates and should be rebuilt rather than
/// re-blessed — check the geometry above before trusting a green run.
#[test]
fn a_reverse_complement_competing_with_mixed_members_is_an_inversion() {
    converges_to(
        "TEMPLATE:g.43_50inv",
        &[
            "TEMPLATE:g.43_50inv",
            "TEMPLATE:g.43_50delinsCACTTGTC",
            "TEMPLATE:g.[43G>C;46_47delinsTT;50G>C]",
        ],
    );
}

/// The reference with `[start, end]` (1-based, inclusive) replaced by its
/// reverse complement — built here, independently of anything the normalizer
/// does, so it can serve as the expected answer.
fn inverted(start: usize, end: usize) -> String {
    let block: String = SEQ.as_bytes()[start - 1..end]
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            other => *other as char,
        })
        .collect();
    format!("{}{}{}", &SEQ[..start - 1], block, &SEQ[end..])
}

/// The sequence a description denotes, obtained by projecting each member
/// through [`hgvs_to_spdi`] and applying the triples to `SEQ`.
///
/// Projected through SPDI rather than re-implemented, so this asks what the
/// description *means* rather than re-deriving it the way the normalizer
/// would. Returns `None` when members overlap and there is no resulting
/// sequence at all — reported as a failure by the caller rather than silently
/// compared equal, which is the shape the cis-allele regressions take.
fn denotes(provider: &JsonProvider, descriptor: &str) -> Option<String> {
    let members: Vec<HgvsVariant> = match parse_hgvs(descriptor).ok()? {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let mut triples = Vec::new();
    for member in &members {
        triples.push(hgvs_to_spdi(member, provider).ok()?);
    }
    // 3'-to-5', longer deletion first at a tied position, matching
    // `common::synthetic::assert_padded_preserving`'s applier.
    triples.sort_by_key(|t| std::cmp::Reverse((t.position, t.deletion.len())));
    let mut edited = SEQ.as_bytes().to_vec();
    let mut claimed = SEQ.len();
    for triple in &triples {
        let start = usize::try_from(triple.position).ok()?;
        let end = start.checked_add(triple.deletion.len())?;
        if end > claimed {
            return None; // overlapping members
        }
        edited.splice(start..end, triple.insertion.bytes());
        claimed = start;
    }
    String::from_utf8(edited).ok()
}

/// Whatever the chosen spelling, it must denote the sequence the input denotes.
///
/// The two tests above compare only strings, so a re-spelling pass that changed
/// the resulting bases would satisfy both. This asserts the denotation instead:
/// every encoding is **normalized**, and the normalizer's output is projected
/// back to a sequence and compared byte-for-byte against [`inverted`], which is
/// built here from `SEQ` and never consults the normalizer.
///
/// It is byte equality, not a substring check: a payload written at the wrong
/// offset still contains the right bases, so `contains` cannot see it.
#[test]
fn both_choices_preserve_the_sequence_they_denote() {
    let reference = provider();
    let normalizer = Normalizer::new(provider());

    let check = |encodings: &[&str], expected: String, block: &str| {
        for input in encodings {
            let output = normalizer
                .normalize(&parse_hgvs(input).expect("parse"))
                .expect("normalize")
                .to_string();
            let denoted = denotes(&reference, &output)
                .unwrap_or_else(|| panic!("`{input}` -> `{output}` has no resulting sequence"));
            assert_eq!(
                denoted, expected,
                "{block}: `{input}` -> `{output}` denotes a different sequence"
            );
        }
    };

    // The #1517 block: 25..32 inverts to TGTGCATT, and every encoding must
    // still denote that after normalization types it as one `inv`.
    check(
        &[
            "TEMPLATE:g.25_32inv",
            "TEMPLATE:g.25_32delinsTGTGCATT",
            "TEMPLATE:g.[25_26delinsTG;31_32delinsTT]",
        ],
        inverted(25, 32),
        "#1517",
    );

    // The #1230 block: 11..14 inverts to CATC, and the split that wins here
    // must denote exactly what the `inv` spelling did.
    check(
        &[
            "TEMPLATE:g.11_14inv",
            "TEMPLATE:g.11_14delinsCATC",
            "TEMPLATE:g.[11G>C;14G>C]",
        ],
        inverted(11, 14),
        "#1230",
    );

    // The mixed block: 43..50 inverts to CACTTGTC. The gate ruling decides how
    // it is SPELLED; it must not decide what it denotes, and this is the half
    // that would catch a widening that typed the wrong span as the inversion.
    check(
        &[
            "TEMPLATE:g.43_50inv",
            "TEMPLATE:g.43_50delinsCACTTGTC",
            "TEMPLATE:g.[43G>C;46_47delinsTT;50G>C]",
        ],
        inverted(43, 50),
        "mixed",
    );
}
