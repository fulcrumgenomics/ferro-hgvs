//! A cis group that nets to an insertion 5' of a transcript's first base must
//! not be spelled `c.-1_1ins…` on a transcript that has no base there.
//!
//! # The defect
//!
//! `NM_TEST.1:c.[1dup;1T>A]` on a transcript whose `cds_start` is `1` came back
//! as `NM_TEST.1:c.-1_1insA`. Both members sit on the first base of the CDS,
//! which on that transcript is also the first base of the *sequence*, so the
//! combined haplotype is "insert `A` immediately 5' of base 1" — a position
//! HGVS cannot anchor an insertion at, and one this transcript does not have.
//!
//! # Why it is quiet
//!
//! `c.-1` is a perfectly legal 5'UTR spelling, so the output parses,
//! re-parses, and is a fixed point. It is wrong only against the reference it
//! was derived from, and only for transcripts whose CDS starts at base 1.
//! `FERRO_ASSERT_IN_BOUNDS` is the one seam oracle that can see it, and it runs
//! in the `Test oracle` job alone.
//!
//! # Where it came from
//!
//! `merge::collapse_overlapping_cis_edits` builds the collapsed group's anchor
//! in dense window arithmetic, so a group netting to a pure insertion at the
//! window's 5' edge produces the insertion-shaped anchor `(del_start,
//! del_start - 1)` — i.e. axis endpoint `0`. #1327 gave that function a `past_end`
//! refusal for the mirror-image case at the 3' end and left the 5' side open;
//! #1772 then renamed the endpoint `0` to `-1` (`name_on_zeroless_axis`) so it
//! would stop rendering as `c.?`/`n.0`, and said in its own doc comment that the
//! residue — what such a group should be *called* — was still open. This is that
//! residue.
//!
//! Nothing downstream could repair it, because every later pass is gated on the
//! axis's **positive body region**: once the collapse has emitted a member
//! starting at `c.-1`, `collect_canonical_edits` refuses the group and
//! `canonicalize_from_sequence` — which already answers this shape correctly on
//! the genomic axis, via `boundary_delins_anchor` — never runs.
//!
//! # What is asserted, and what makes each assertion non-vacuous
//!
//! The load-bearing check is **not** a string comparison. Each member of the
//! output is resolved through [`hgvs_to_spdi`], which is the coordinate mapper
//! rather than the normalizer, and which declines a `c.`/`n.` position that
//! resolves to a non-positive transcript base. A test that only compared
//! strings would pass against any future spelling that happened to differ from
//! the recorded one, including another out-of-range one.
//!
//! Three controls keep the guard honest, because "every coordinate resolves" is
//! satisfiable by doing nothing useful:
//!
//! * [`the_in_bounds_instrument_tells_the_two_transcripts_apart`] pins that the
//!   instrument answers *differently* for the two references — `c.-1_1insA` is
//!   refused against a CDS starting at base 1 and accepted against one with a
//!   real 5'UTR. Without it, a change that made the mapper accept `c.-1`
//!   everywhere would turn every other test here green while making the defect
//!   worse.
//! * [`a_transcript_with_a_real_five_prime_utr_still_spells_it_at_c_minus_one`]
//!   is the over-fix control: banning `c.-1` outright would satisfy the
//!   in-bounds checks and destroy a correct answer.
//! * [`an_interior_pair_on_the_same_transcript_still_merges`] and
//!   [`the_group_is_still_merged_into_one_member`] are the under-fix controls:
//!   refusing to combine anything at all also produces coordinates that all
//!   resolve.
//!
//! # The form, and its authority
//!
//! `c.1delinsAT` is not a free choice. `DNA/insertion.md:95-101` requires an
//! insertion to name two adjacent flanking positions, so a payload resting 5' of
//! base 1 has no `ins` spelling at all; the single-position `delins` is the same
//! identity `insertion_to_boundary_delins` already applies at every other
//! sequence bound — the contig start on `g.` (#1205), the mitochondrial origin
//! (#1217) and the transcript bounds on `n.`/`r.` (#1202). Deriving it from the
//! resulting sequence rather than preserving the input's spelling is
//! `rulings[canonical-form-choice-when-both-legal]`.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};

/// The transcript sequence every case here is drawn against.
///
/// Its first base is `T` and the rest is the `ACGT` period-4 rotation, so no
/// base equals either of its neighbours and nothing in these cases can shift
/// along a run — the only movement on show is the one the collapse performs.
const CORE: &str = "TACGTACGTACGTACGTACG";

/// A codon-complete CDS end, leaving five 3'UTR bases. Nothing here reaches the
/// 3' end; it is fixed so the two providers differ in `cds_start` alone.
const CDS_END: u64 = 15;

/// `cds_start = 1`: the CDS begins at the transcript's first base, so the
/// transcript has **no** 5'UTR and `c.-1` names nothing.
const NO_FIVE_PRIME_UTR: u64 = 1;

/// `cds_start = 4`: three 5'UTR bases, so `c.-1` is transcript base 3 and the
/// same insertion spelling is correct.
const WITH_FIVE_PRIME_UTR: u64 = 4;

fn coding_provider(cds_start: u64) -> MockProvider {
    SyntheticBuilder::cds(CORE, cds_start, CDS_END, Strand::Plus).build()
}

fn noncoding_provider() -> MockProvider {
    SyntheticBuilder::noncoding(CORE, Strand::Plus).build()
}

/// The members of `description`, whether or not it is bracketed.
fn members(description: &str) -> Vec<HgvsVariant> {
    match parse_hgvs(description).unwrap_or_else(|e| panic!("`{description}` parses: {e}")) {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    }
}

fn normalized(provider: &MockProvider, input: &str) -> String {
    Normalizer::new(provider.clone())
        .normalize(&parse_hgvs(input).unwrap_or_else(|e| panic!("`{input}` parses: {e}")))
        .unwrap_or_else(|e| panic!("`{input}` normalizes: {e}"))
        .to_string()
}

/// The members of `description` that name a coordinate `provider`'s transcript
/// does not have, each with the mapper's own reason.
///
/// The question is put to [`hgvs_to_spdi`] — the coordinate mapper — and not to
/// the normalizer, so nothing here can agree with an output merely because
/// normalization produced it. The mapper declines a `c.`/`n.` position that
/// resolves to a non-positive transcript base, which is exactly the property
/// this issue is about.
fn unresolvable_members(provider: &MockProvider, description: &str) -> Vec<String> {
    members(description)
        .iter()
        .filter_map(|member| match hgvs_to_spdi(member, provider) {
            Ok(_) => None,
            Err(e) => Some(format!("{member}: {e}")),
        })
        .collect()
}

/// The bases `description` denotes over [`CORE`], by splicing each member's
/// SPDI triple into the transcript sequence.
///
/// `Err` for a description the mapper cannot resolve or whose members overlap;
/// a description that denotes no sequence is a failure to report, never a
/// comparison to skip.
fn denoted(provider: &MockProvider, description: &str) -> Result<String, String> {
    let mut triples = Vec::new();
    for member in members(description) {
        triples.push(hgvs_to_spdi(&member, provider).map_err(|e| format!("{member}: {e}"))?);
    }
    // 3'-to-5', longer deletion first at a tied position, so a zero-width
    // insertion sharing an interbase with a deletion does not read as an
    // overrun. Same order `common::synthetic`'s applier uses.
    triples.sort_by_key(|t| std::cmp::Reverse((t.position, t.deletion.len())));
    let mut edited = CORE.as_bytes().to_vec();
    let mut claimed = CORE.len();
    for triple in &triples {
        let start =
            usize::try_from(triple.position).map_err(|_| "position overflows".to_string())?;
        let end = start + triple.deletion.len();
        if end > claimed {
            return Err(format!("members overlap at interbase {start}"));
        }
        edited.splice(start..end, triple.insertion.bytes());
        claimed = start;
    }
    String::from_utf8(edited).map_err(|e| e.to_string())
}

/// The reproducer from the issue: both members on the first coding base of a
/// transcript whose CDS starts at base 1.
const REPRODUCER: &str = "NM_TEST.1:c.[1dup;1T>A]";

/// The out-of-range string the reproducer produced before the fix. Named so the
/// controls below can state what they are calibrated against.
const OUT_OF_RANGE_OUTPUT: &str = "NM_TEST.1:c.-1_1insA";

#[test]
fn the_allele_at_the_first_coding_base_names_only_coordinates_the_transcript_has() {
    let provider = coding_provider(NO_FIVE_PRIME_UTR);
    let output = normalized(&provider, REPRODUCER);
    let unresolvable = unresolvable_members(&provider, &output);
    assert!(
        unresolvable.is_empty(),
        "`{REPRODUCER}` -> `{output}` names a coordinate NM_TEST.1 does not \
         have (cds_start = {NO_FIVE_PRIME_UTR}, so there is no base 5' of c.1): \
         {unresolvable:?}",
    );
}

#[test]
fn the_output_denotes_the_bases_the_input_did() {
    let provider = coding_provider(NO_FIVE_PRIME_UTR);
    let output = normalized(&provider, REPRODUCER);
    let from_input = denoted(&provider, REPRODUCER).expect("the input denotes a sequence");
    let from_output = denoted(&provider, &output)
        .unwrap_or_else(|e| panic!("`{REPRODUCER}` -> `{output}` denotes no sequence: {e}"));
    assert_eq!(
        from_output, from_input,
        "`{REPRODUCER}` -> `{output}` changed the sequence",
    );
}

#[test]
fn the_canonical_form_is_the_single_position_delins_at_the_first_base() {
    let provider = coding_provider(NO_FIVE_PRIME_UTR);
    assert_eq!(
        normalized(&provider, REPRODUCER),
        "NM_TEST.1:c.1delinsAT",
        "an insertion resting 5' of base 1 has no two-anchor `ins` spelling \
         (DNA/insertion.md:95-101), so it is the single-position delins the \
         genomic axis already emits for the same shape",
    );
}

#[test]
fn the_group_is_still_merged_into_one_member() {
    // The under-fix control for the string above: declining to combine also
    // produces an output whose every coordinate resolves, and would leave the
    // two members as they were.
    let provider = coding_provider(NO_FIVE_PRIME_UTR);
    let output = normalized(&provider, REPRODUCER);
    assert_eq!(
        members(&output).len(),
        1,
        "`{REPRODUCER}` -> `{output}`: the two members denote one contiguous \
         change and must still collapse to one",
    );
}

#[test]
fn the_non_coding_axis_has_the_same_shape_and_the_same_answer() {
    // `n.` has no negative zone at all, so `n.-1` is worse than `c.-1`: it is
    // refused at strict parse (#1751) as well as being off the sequence.
    let provider = noncoding_provider();
    let input = "NR_TEST.1:n.[1dup;1T>A]";
    let output = normalized(&provider, input);
    let unresolvable = unresolvable_members(&provider, &output);
    assert!(
        unresolvable.is_empty(),
        "`{input}` -> `{output}` names a coordinate NR_TEST.1 does not have: \
         {unresolvable:?}",
    );
    assert_eq!(output, "NR_TEST.1:n.1delinsAT");
}

#[test]
fn a_transcript_with_a_real_five_prime_utr_still_spells_it_at_c_minus_one() {
    // The over-fix control. Same shape, same first coding base — but here
    // `c.-1` is transcript base 3, so the insertion has two real anchors and
    // `c.-1_1insA` is the right answer. A fix that refused the collapse
    // whenever an endpoint rendered negative would break this.
    let provider = coding_provider(WITH_FIVE_PRIME_UTR);
    let input = "NM_TEST.1:c.[1dup;1G>A]";
    let output = normalized(&provider, input);
    assert_eq!(output, OUT_OF_RANGE_OUTPUT);
    assert!(
        unresolvable_members(&provider, &output).is_empty(),
        "`{input}` -> `{output}`: with cds_start = {WITH_FIVE_PRIME_UTR} every \
         coordinate here exists",
    );
}

#[test]
fn an_interior_pair_on_the_same_transcript_still_merges() {
    // The other under-fix control: a guard written on the collapse's output
    // rather than on its 5' endpoint would switch this off too.
    let provider = coding_provider(NO_FIVE_PRIME_UTR);
    let input = "NM_TEST.1:c.[2dup;2A>C]";
    let output = normalized(&provider, input);
    assert_eq!(output, "NM_TEST.1:c.1_2insC");
    assert!(unresolvable_members(&provider, &output).is_empty());
}

#[test]
fn the_in_bounds_instrument_tells_the_two_transcripts_apart() {
    // Calibration. Every assertion above is only as good as the mapper's
    // ability to answer this question differently for the two references — if
    // it accepted `c.-1` unconditionally the guards would all go green while
    // the defect got worse, and if it refused `c.-1` unconditionally the
    // over-fix control above would be measuring nothing.
    let no_utr = coding_provider(NO_FIVE_PRIME_UTR);
    let with_utr = coding_provider(WITH_FIVE_PRIME_UTR);
    assert_eq!(
        unresolvable_members(&no_utr, OUT_OF_RANGE_OUTPUT).len(),
        1,
        "`{OUT_OF_RANGE_OUTPUT}` must be refused against a CDS starting at \
         base 1 — that refusal is what every other test in this file reads",
    );
    assert!(
        unresolvable_members(&with_utr, OUT_OF_RANGE_OUTPUT).is_empty(),
        "`{OUT_OF_RANGE_OUTPUT}` must be accepted against a transcript with a \
         real 5'UTR, or the instrument is a blanket ban rather than a bounds \
         check",
    );
}
