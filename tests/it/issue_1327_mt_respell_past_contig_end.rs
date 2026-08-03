//! Issue #1327 — a junction re-spelling could place an edit one base past the
//! end of a circular mitochondrial contig.
//!
//! `respell_at_gap` writes its member over the gap `[junction, junction + 1]`.
//! When a member 3'-shifts to rest on the contig's final base that gap is
//! `[16569, 16570]` on rCRS, and 16570 does not exist. Both callers admit
//! `CisKind::Mt`, and both were reachable:
//!
//! ```text
//! NC_012920.1:m.[16566del;16567dup]            -> m.16569_16570=
//! NC_012920.1:m.[16566_16567insA;16567_16568insA] -> m.16569_16570insAA
//! ```
//!
//! ## Why the existing guards do not catch it
//!
//! `respell_at_gap`'s `landed` check only confirms the member reads back at the
//! gap it was *told* to use, which an out-of-range gap does. And
//! `FERRO_ASSERT_REPARSE` is blind here: `parse_hgvs` does not know the
//! contig's length, so `m.16569_16570insAA` re-parses cleanly. The defect is
//! therefore silent in both oracles, which is why it needs its own assertions.
//!
//! ## What the repair emits instead, and why not the other two options
//!
//! **Not the wrapped form.** The obvious circular answer — spell the gap
//! `m.16569_1` — is not available for an insertion. SVD-WG006 authorises the
//! reversed `<high>_<low>` range only for deletions and duplications
//! ("deletions/duplications", line 33, with the `NC_012920.1:m.16563_13del` and
//! `J01749.1:o.4344_197dup` examples); `insertion.md` is silent, so the general
//! 5'→3' rule applies and ferro's parser rejects `m.16569_1insA` (#129). Every
//! re-spelling here produces an insertion, so there is no legal wrapped form.
//!
//! **Not a refusal either**, though the issue offers it. Measured: refusing
//! left `coalesce_members_at_one_junction`'s two members unmerged, and each had
//! already been canonicalised to `m.16569dup` — an allele claiming one position
//! twice, which is the #1286 defect that pass exists to remove, and which
//! re-parses so no oracle sees it. That trades an out-of-range spelling for an
//! overlapping one.
//!
//! **A boundary `delins`.** Inserting `A'` immediately 3' of the last base *is*
//! deleting that base and inserting `ref[last] ++ A'` — the identity every
//! other clamp in the crate already uses (`insertion_to_boundary_delins`,
//! shared with #1170 / #387 / #1202 / #1205 / #1217). In range, valid on a
//! circular or linear reference alike, and denoting exactly what the
//! out-of-range insertion denoted.
//!
//! A second producer of the same overrun — `collapse_overlapping_cis_edits`
//! naming `del_start` one past the end when a group nets to a boundary
//! insertion — is fixed alongside it. There refusing *is* right, because it
//! leaves the members in the in-range spelling they already had.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

const MT: &str = "NC_012920.1";
const MT_LENGTH: usize = 16569;

/// An rCRS-length circular contig ending in a 4-base run, so a junction member
/// can 3'-shift onto the final base.
fn mt_provider(tail: &str) -> MockProvider {
    let mut sequence = "ACGT".repeat(4000);
    sequence.push_str(&"C".repeat(MT_LENGTH - sequence.len() - tail.len()));
    sequence.push_str(tail);
    assert_eq!(sequence.len(), MT_LENGTH, "contig must be rCRS length");
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(MT, sequence);
    provider
}

fn normalize(input: &str, tail: &str) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::new(mt_provider(tail))
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string()
}

/// Every endpoint the output names must exist on the contig.
///
/// Checked by parsing the rendered output and reading its coordinates back,
/// rather than by string-matching `"16570"`, so an overrun by any amount at any
/// endpoint is caught.
///
/// Deliberately over-approximate: it scans *every* integer after the `:`, so a
/// repeat count (`A[6]`) or an inserted length is checked as though it were a
/// position. That is sound for these fixtures — the counts are single digits
/// against a 16569-base contig — and it fails safe, since the risk is a
/// spurious failure rather than a missed overrun. If a future case here emits a
/// count above the contig length, narrow this to the location field rather than
/// loosening the bound.
fn assert_within_contig(input: &str, output: &str) {
    let reparsed = parse_hgvs(output).unwrap_or_else(|e| {
        panic!("`{input}` normalized to `{output}`, which will not re-parse: {e}")
    });
    // `parse_hgvs` does not length-check, so the coordinates have to be read
    // out of the rendered string; every integer in the position field must be
    // a real 1-based position on the contig.
    let positions = output
        .rsplit_once(':')
        .map(|(_, rest)| rest)
        .unwrap_or(output);
    for run in positions.split(|c: char| !c.is_ascii_digit()) {
        if run.is_empty() {
            continue;
        }
        let Ok(position) = run.parse::<usize>() else {
            continue;
        };
        assert!(
            (1..=MT_LENGTH).contains(&position),
            "`{input}` normalized to `{output}` ({reparsed}), which names position \
             {position} on a {MT_LENGTH}-base contig"
        );
    }
}

/// The `respell_colliding_duplications` caller: a deletion and a duplication
/// that both 3'-shift through the terminal run and collide on the last base.
#[test]
fn a_collision_repaired_at_the_last_base_stays_on_the_contig() {
    for input in [
        "NC_012920.1:m.[16566del;16567dup]",
        "NC_012920.1:m.[16566del;16568dup]",
        "NC_012920.1:m.[16567dup;16566del]",
    ] {
        let output = normalize(input, "AAAA");
        assert_within_contig(input, &output);
        // The pair nets to no change at all, so an identity is the right
        // answer; what matters is that it is spelled in range. Pinned as well
        // as bounded, because "stays on the contig" alone would also be
        // satisfied by silently dropping a member.
        assert_eq!(
            output, "NC_012920.1:m.16568=",
            "`{input}` must settle on an in-range identity"
        );
    }
}

/// The `coalesce_members_at_one_junction` caller (#1289), which places its
/// merged insertion through the same helper.
#[test]
fn a_junction_merge_at_the_last_base_stays_on_the_contig() {
    for input in [
        "NC_012920.1:m.[16566_16567insA;16567_16568insA]",
        "NC_012920.1:m.[16565_16566insA;16566_16567insA]",
    ] {
        let output = normalize(input, "AAAA");
        assert_within_contig(input, &output);
        // Two `A`s added to the terminal 4-`A` tract. The boundary `delins`
        // this repair emits is re-canonicalised to the repeat form by the next
        // pass, which is the spelling the same edit gets anywhere else in the
        // contig — the terminus is no longer a special case in the output.
        //
        // Pinning this specifically guards the regression that declining
        // caused: refusing left the two members unmerged as
        // `m.[16569dup;16569dup]`, one position claimed twice (#1286).
        assert_eq!(
            output, "NC_012920.1:m.16566_16569A[6]",
            "`{input}` must merge into one in-range member"
        );
        assert!(
            !output.contains(';'),
            "`{input}` must not come back as a multi-member allele; got `{output}`"
        );
    }
}

/// The wrapped insertion spelling is not an option, so this pins the reason the
/// repair declines rather than producing it.
///
/// If `m.16569_1insA` ever becomes parseable, this test fails and the decision
/// above should be revisited — the wrapped form would then be strictly better
/// than declining.
#[test]
fn the_wrapped_insertion_spelling_is_not_valid_hgvs() {
    assert!(
        parse_hgvs("NC_012920.1:m.16569_1insA").is_err(),
        "SVD-WG006 authorises the reversed range for del/dup only; a wrapped \
         `ins` must stay rejected (#129)"
    );
    // The del/dup forms it *does* authorise are still accepted, so the
    // assertion above is about the edit type and not about wraparound at large.
    assert!(parse_hgvs("NC_012920.1:m.16569_1del").is_ok());
    assert!(parse_hgvs("NC_012920.1:m.16569_1dup").is_ok());
}

/// Declining must not cost the in-range repairs. A collision one base earlier
/// resolves exactly as it did before.
#[test]
fn an_in_range_collision_is_still_repaired() {
    // Terminal run of `T`s with the members in a `C` run further 5', so the
    // shuffle settles well inside the contig.
    let input = "NC_012920.1:m.[16000del;16001dup]";
    let output = normalize(input, "TTTT");
    assert_within_contig(input, &output);
    // Pinned, not merely "parseable". A refusal, a dropped member, or some other
    // in-range spelling all satisfy a parse check, and every one of those would
    // be the regression this control exists to catch — the point is that
    // declining at the contig end costs the in-range repair nothing, which only
    // a pinned output can state.
    assert_eq!(
        output, "NC_012920.1:m.16000T>C",
        "`{input}` must still collapse to the single substitution it did before"
    );
    assert!(
        parse_hgvs(&output).is_ok(),
        "an in-range repair must still produce a parseable result"
    );
}

/// The same overrun on a **transcript** axis, where the last base is `c.*N`.
///
/// Not a mitochondrial special case: any sequence has a last base, and a
/// duplication resting on it lands its copy past the end whatever the axis is
/// called. This is here because the first cut of the fix refused instead of
/// re-spelling, and refusing on this axis left #1284's collision unrepaired and
/// the output unparseable — the re-parse oracle caught it. Pinned so the
/// transcript half cannot silently regress to a refusal again.
#[test]
fn the_same_overrun_on_a_transcript_axis_is_repaired() {
    use crate::common::synthetic::SyntheticBuilder;
    use ferro_hgvs::reference::transcript::Strand;

    // 35 bases, CDS 13..24, so the 3'UTR is `c.*1`..`c.*11` and `c.*11` is the
    // transcript's final base. The terminal `TTTTAA` gives the members a run to
    // shuffle through.
    const CORE: &str = "GCAAAGCGCGCGATGAAACCCTAAGGCATTTTTAA";
    let input = "NM_TEST.1:c.[*10del;*11dup]";
    let variant = parse_hgvs(input).expect("input must parse");
    let output = Normalizer::new(SyntheticBuilder::cds(CORE, 13, 24, Strand::Plus).build())
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string();
    parse_hgvs(&output).unwrap_or_else(|e| {
        panic!("`{input}` normalized to `{output}`, which will not re-parse: {e}")
    });
    assert!(
        !output.contains("*12"),
        "`{input}` normalized to `{output}`, which names a position past the transcript end"
    );
    // Pinned, not merely bounded. `*10` and `*11` are both `A`, so deleting the
    // first and duplicating the second cancel exactly and the allele denotes the
    // reference — the single identity member #1284 settles on.
    //
    // A bounds-and-parseability check is satisfied by a refusal, a dropped
    // member, or any other in-range spelling, and that is how #1344 hid here: an
    // earlier cut of this fix re-spelled the member as `c.*11delinsAA`, giving
    // two coincident edits on `c.*11` — rejected by ferro's own strict mode as
    // `OverlapConflictingEdits / W5002`, and a stable fixed point, so a permanent
    // wrong answer that passed both of the weaker checks.
    assert_eq!(
        output, "NM_TEST.1:c.*11=",
        "`{input}` must settle on the identity, not an overlapping pair"
    );
}

/// A lone duplication at the final base has no sibling to collide with and must
/// pass through untouched — the control that shows the guard is scoped to the
/// repair and is not clamping ordinary terminal edits.
#[test]
fn a_lone_terminal_duplication_is_unchanged() {
    assert_eq!(
        normalize("NC_012920.1:m.16566dup", "AAAA"),
        "NC_012920.1:m.16569dup"
    );
}
