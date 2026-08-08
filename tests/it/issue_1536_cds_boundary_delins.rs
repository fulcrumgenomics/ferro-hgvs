//! #1536 — a lone `delins` whose span crosses a CDS boundary is never
//! re-derived from the sequence it denotes, so two spellings of one variant
//! reach two different fixed points.
//!
//! The identical block placed entirely inside the CDS converges correctly. The
//! discriminator is nothing but where the start or stop codon happens to fall.
//!
//! # Mechanism
//!
//! `join_pos` (`src/normalize/merge.rs:1085-1087`) refuses any interval whose
//! two ends are in different regions:
//!
//! ```text
//! if rs != re {
//!     return None;
//! }
//! ```
//!
//! `simple_cds_pos` classifies a position as `Cds`, `FivePrimeUtr` or
//! `ThreePrimeUtr` and returns the axis coordinate *within that region* — `*1`
//! is `1`, `-3` is `-3`, `9` is `9`. Those numbers are not on one line, so
//! `edit_span_union` and the `w_lo`/`w_hi` window arithmetic cannot be done in
//! axis coordinates across a boundary. `join_pos` therefore declines,
//! `collect_canonical_edits` returns `None`, and the whole sequence-first pass
//! is skipped for the variant.
//!
//! Unlike the exon-junction clamp a few lines further down — which carries an
//! explicit rationale and states its confluence cost — this refusal is an
//! unannotated consequence of doing window arithmetic in axis space. A fix
//! converts the span to transcript coordinates for the window computation and
//! back afterwards, the way the projector already does; the region split is a
//! rendering concern, not a geometric one.
//!
//! # What is measured here
//!
//! One 20-base transcript, one 8-base block at transcript positions 5..=12
//! (`AATGCACA`), one payload (`TGTGCATT`, its exact reverse complement). Only
//! the CDS annotation moves, so the block, the bases and the edit are held
//! constant and the boundary is the sole variable:
//!
//! | CDS | the same block reads | lone `delins` re-derives to `inv` |
//! |---|---|---|
//! | `1..=20` | `c.5_12`  | yes |
//! | `1..=8`  | `c.5_*4`  | **no** |
//! | `8..=20` | `c.-3_5`  | **no** |
//!
//! # Status
//!
//! **Unfixed.** [`a_delins_crossing_a_cds_boundary_converges_with_its_inv_spelling`]
//! is `#[ignore]`d and RED against #1536 — an ignored red guard is a recorded
//! defect; a weakened green one is a lie. Run it with
//! `cargo nextest run --features dev -E 'test(issue_1536)' --run-ignored all`.
//!
//! The control — the same block entirely inside the CDS — is **not** ignored,
//! so the half that works today is gated.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// 20 bases, with the 8-base block `AATGCACA` at transcript positions 5..=12.
/// `TGTGCATT` is its exact whole-block reverse complement, so `delins` of that
/// payload denotes precisely `inv` of that span — nothing about the pair
/// depends on where the CDS is.
const CORE: &str = "GCTGAATGCACATCCGTAGC";

fn provider(cds_start: u64, cds_end: u64) -> MockProvider {
    let mut provider = MockProvider::new();
    let length = CORE.len() as u64;
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        CORE.to_string(),
        Some(cds_start),
        Some(cds_end),
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

fn normalize(provider: &MockProvider, input: &str, direction: ShuffleDirection) -> String {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
    Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::default().with_direction(direction),
    )
    .normalize(&variant)
    .unwrap_or_else(|e| panic!("{input} must normalize: {e}"))
    .to_string()
}

/// The three placements of one block: `(cds_start, cds_end, delins, inv)`.
///
/// The `inv` spelling is the same variant written the other way. Confluence
/// means the two reach one output; the defect is that across a boundary they
/// do not.
const PLACEMENTS: [(u64, u64, &str, &str); 3] = [
    // Entirely inside the CDS — the control, which converges today.
    (
        1,
        20,
        "NM_TEST.1:c.5_12delinsTGTGCATT",
        "NM_TEST.1:c.5_12inv",
    ),
    // Across the stop codon.
    (
        1,
        8,
        "NM_TEST.1:c.5_*4delinsTGTGCATT",
        "NM_TEST.1:c.5_*4inv",
    ),
    // Across the start codon.
    (
        8,
        20,
        "NM_TEST.1:c.-3_5delinsTGTGCATT",
        "NM_TEST.1:c.-3_5inv",
    ),
];

/// The half that works: inside the CDS, the two spellings converge.
///
/// Not ignored. It is the control the whole issue rests on — "the identical
/// block placed entirely inside the CDS converges correctly" — so if it ever
/// stops holding, the defect below has been mis-stated and the record needs
/// rewriting rather than the fix landing.
#[test]
fn a_delins_inside_the_cds_converges_with_its_inv_spelling() {
    let (cds_start, cds_end, delins, inv) = PLACEMENTS[0];
    let provider = provider(cds_start, cds_end);
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let from_delins = normalize(&provider, delins, direction);
        assert_eq!(
            from_delins,
            normalize(&provider, inv, direction),
            "{delins} and {inv} denote one variant inside the CDS and must converge"
        );
        // Named rather than merely equal: the issue's claim is not only that
        // they agree but that the `delins` is *re-derived from the sequence*
        // and typed `inv`. Two spellings that both froze would satisfy an
        // equality check and prove nothing.
        assert_eq!(
            from_delins, inv,
            "{delins} must be re-derived and typed as an inversion"
        );
    }
}

/// #1536: across a CDS boundary they do not.
///
/// Measured on `320e98dc` with the geometry in the module docs — the same
/// block, the same payload, only the CDS annotation moved:
///
/// ```text
/// CDS 1..=20   c.5_12delinsTGTGCATT   -> c.5_12inv     (converges)
/// CDS 1..=8    c.5_*4delinsTGTGCATT   -> unchanged     (does not)
/// CDS 8..=20   c.-3_5delinsTGTGCATT   -> unchanged     (does not)
/// ```
///
/// So `c.5_*4inv` and `c.5_*4delinsTGTGCATT` are both fixed points: two
/// answers for one variant, discriminated by nothing but where the stop codon
/// falls.
///
/// Ignored because #1536 is unfixed, not because the assertion is wrong.
/// Un-ignore it with the fix, in the same change.
#[test]
#[ignore = "#1536: a delins crossing a CDS boundary is never re-derived — unfixed, this guard is red"]
fn a_delins_crossing_a_cds_boundary_converges_with_its_inv_spelling() {
    let mut divergent = Vec::new();
    for (cds_start, cds_end, delins, inv) in PLACEMENTS.iter().skip(1) {
        let provider = provider(*cds_start, *cds_end);
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            let from_delins = normalize(&provider, delins, direction);
            let from_inv = normalize(&provider, inv, direction);
            if from_delins != from_inv {
                divergent.push(format!(
                    "  CDS {cds_start}..={cds_end}: {delins} -> {from_delins}, but {inv} -> {from_inv}"
                ));
            }
        }
    }
    assert!(
        divergent.is_empty(),
        "#1536: {} boundary-crossing spelling pair(s) reached two fixed points for one \
         variant:\n{}",
        divergent.len(),
        divergent.join("\n")
    );
}

/// The premise, asserted separately so a fix cannot be credited to a change of
/// geometry.
///
/// If the block ever stops being an exact whole-block reverse complement, the
/// guard above becomes vacuous — `delins` and `inv` would no longer denote one
/// variant and their divergence would be correct. Checked from the bases
/// themselves rather than trusted from the constant.
#[test]
fn the_payload_is_the_exact_reverse_complement_of_the_block() {
    let block = &CORE[4..12];
    let revcomp: String = block
        .chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => panic!("unexpected base {other}"),
        })
        .collect();
    assert_eq!(block, "AATGCACA");
    assert_eq!(revcomp, "TGTGCATT");
}
