//! Tests using a real ferro-prepared reference.
//!
//! The reference is located through `FERRO_MANIFEST`, falling back to
//! `benchmark-output/manifest.json` — the convention every reference-aware
//! module here uses. Without one these guards stand down; under
//! `FERRO_REQUIRE_MANIFEST` standing down is a failure, whether the manifest is
//! absent or merely does not carry the transcript a guard names.
//!
//! **That promotion bought execution, not verification, and #1858 is the
//! difference.** Three of these guards used to carry no assertion at all — they
//! were investigation scripts that printed a finding and returned successfully,
//! so the condition each existed to detect was reported by *passing* — and
//! `test_compare_with_mutalyzer_sequence` kept both of its assertions inside a
//! bounds check. This module now asserts, and what each guard asserts is chosen
//! by what the reference can actually decide.
//!
//! # The one property that needs no adjudication
//!
//! Which of two legal spellings ferro ships is a rule 6 choice (`README.md`'s
//! ruleset) and needs an operator to settle. Whether the shipped spelling still
//! denotes the input's bases is rule 1, and a failure there is a bug under
//! every reading of the recommendations. So [`assert_denotes_the_same_bases`]
//! carries the load wherever the spelling question is open, and no guard here
//! pins a string that nobody has justified.
//!
//! # `unwrap_or(1)` is removed because it *can* fabricate, not because it did
//!
//! Every guard used to read its CDS origin as `transcript.cds_start
//! .unwrap_or(1)`. That substitutes an origin when the reference serves none,
//! and the substituted coordinates print exactly like real ones — a guard
//! reporting `c.1324` would be naming a base the description does not name,
//! with nothing in its output saying so. The fallback is gone; an unserved CDS
//! now routes through [`cds_start_or_skip`], which stands the guard down and
//! records it.
//!
//! **On the blessed reference that branch is dead, and this is the correction
//! of a claim this module previously shipped as measured fact.** It stated that
//! `NM_033517` is served as bare sequence with no cdot record —
//! `cds_start: None`, `cds_end: None`, zero exons. Re-measured 2026-08-14
//! against both prepared volumes, all four accessions this module names come
//! back annotated:
//!
//! ```text
//! NM_033517.1/.2/.3  cds_start=Some(1)   cds_end=Some(5157)  exons=22  len=7113
//! NM_001408491.1     cds_start=Some(296) cds_end=Some(2365)  exons=23  len=3748
//! NM_001282424.3     cds_start=Some(169) cds_end=Some(3060)  exons=25  len=3791
//! NM_001394148.2     cds_start=Some(137) cds_end=Some(328)   exons=3   len=570
//! ```
//!
//! `NM_033517`'s cdot record carries `start_codon: 0`, so `cds_start` is 1 and
//! a `c.N` on it resolves to transcript offset `N` — the old fallback happened
//! to agree with the reference here rather than fabricate against it. Only the
//! length (7113) survives from the withdrawn reading. Removing the fallback is
//! still right, and [`cds_start_or_skip`] is still the right shape; what is
//! withdrawn is the claim that anything was observed going through it.
//!
//! Note `cds_end` is the one figure above that is **cdot-version-dependent**,
//! and nothing here keys on it. Read out of both cdot builds directly:
//! `data_v0.2.32` — the version the prepared manifest pins — collapses this
//! accession's 39-base alignment gap (transcript 1302..1342, between exons 10
//! and 11) into the stop codon and records `stop_codon: 5157`, while
//! `data_v0.2.34` records `stop_codon: 5196`, which is what RefSeq/GenBank
//! annotate (`CDS 1..5196`). That discrepancy is the open sibling of
//! `rulings[c-and-n-positions-are-flat-transcript-offsets]` and is not this
//! module's business.
//!
//! **Both versions carry the record**, with 22 exons and `start_codon: 0`, so
//! `cds_start` is 1 either way and no cdot in play yields the record-less state
//! the withdrawn reading described. The correction above therefore does not
//! depend on which cdot the reference was prepared from.

use ferro_hgvs::hgvs::HgvsVariant;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::spdi::{compare_denoted_sequences, DenotedSequenceComparison};
use ferro_hgvs::{parse_hgvs, MultiFastaProvider, Normalizer, ReferenceProvider};
use std::path::Path;
use std::sync::Arc;

/// This module's name, as it appears in a `FERRO_REQUIRE_MANIFEST` failure.
const MODULE: &str = "real_data_normalization_tests";

/// The prepared manifest, by the same convention every other reference-aware
/// module here uses: `FERRO_MANIFEST` is authoritative when set, with
/// `benchmark-output/manifest.json` as the well-known fallback.
///
/// This module used to consult **only** the relative path, which made it
/// orphaned in a stronger sense than its siblings: an operator following the
/// documented local recipe (`export FERRO_MANIFEST=…`) did not arm it either,
/// so it was dead even for the developer trying to reproduce a nightly result.
fn manifest_path() -> Option<std::path::PathBuf> {
    if let Ok(path) = std::env::var("FERRO_MANIFEST") {
        let path = std::path::PathBuf::from(path);
        return path.exists().then_some(path);
    }
    let fallback = Path::new("benchmark-output/manifest.json");
    fallback.exists().then(|| fallback.to_path_buf())
}

/// The provider, or `None` when there is no manifest to build one from.
///
/// A manifest that is **present but will not load** panics rather than
/// returning `None`, which is the other half of the convention `manifest_path`
/// above adopts — `biocommons_normalize_tests::provider` and
/// `multi_member_cis_axis::provider` each say so in their own words. It matters
/// more here than the wording suggests: `None` reaches
/// `common::manifest::absent`, whose text states that no manifest was found and
/// tells the reader to prepare one. For a manifest that *is* there and is
/// broken that is a wrong diagnosis, and under `FERRO_REQUIRE_MANIFEST` it is a
/// wrong diagnosis on the one job whose whole purpose is to notice that
/// `ferro prepare` did not leave a usable reference behind.
fn get_provider() -> Option<MultiFastaProvider> {
    let path = manifest_path()?;
    Some(
        MultiFastaProvider::from_manifest(&path)
            .unwrap_or_else(|e| panic!("from_manifest({}) failed: {e}", path.display())),
    )
}

/// The transcript a guard names, or `None` after recording the stand-down.
///
/// A manifest that is present but does not carry the accession leaves the guard
/// asserting nothing, which is the same coverage loss an absent manifest causes
/// and is promoted by the same `FERRO_REQUIRE_MANIFEST` — see
/// `ferro_hgvs::conformance::manifest_gate::unserved` for why one flag and two
/// messages rather than two flags.
///
/// Each call site keeps its own `return`, so what a reader sees is still an
/// ordinary early exit.
fn transcript_or_skip(provider: &MultiFastaProvider, accession: &str) -> Option<Arc<Transcript>> {
    match provider.get_transcript(accession) {
        Ok(transcript) => Some(transcript),
        Err(error) => {
            crate::common::manifest::unserved(MODULE, accession, &error.to_string());
            None
        }
    }
}

/// The transcript's CDS origin, or `None` after recording the stand-down.
///
/// **Never substitute a default here.** `cds_start.unwrap_or(1)` turns "the
/// reference cannot resolve this accession's `c.` axis" into "every `c.N` is
/// transcript offset `N`" — real bases, at coordinates the description does not
/// name, printed identically to a resolved reading. A guard that publishes a
/// ferro-versus-Mutalyzer verdict off that has measured something other than
/// what it reports, and nothing in its output distinguishes the two.
///
/// **No accession this module names currently reaches the `None` branch** — see
/// the module header for the measured `cds_start` of all four. So this is a
/// safety net rather than a repair of an observed failure, and it must not be
/// cited as evidence that one occurred. It earns its place from the shape of
/// the failure it removes, not from a sighting: a stand-down is loud and a
/// substituted origin is silent, and #1858 is the record of what silence costs.
/// Since #1870, `Normalizer::normalize` refuses a `c.` description against an
/// unresolvable CDS rather than answering, so the branch would now be reached
/// before the normalizer is asked
/// (`rulings[c-description-against-an-unresolvable-cds-is-refused]`).
///
/// A served sequence with no served CDS is the same coverage loss as an
/// unserved transcript one level up, so it is promoted by the same
/// `FERRO_REQUIRE_MANIFEST` and through the same
/// `ferro_hgvs::conformance::manifest_gate::unserved`.
fn cds_start_or_skip(transcript: &Transcript, accession: &str) -> Option<u64> {
    match transcript.cds_start {
        Some(cds_start) => Some(cds_start),
        None => {
            crate::common::manifest::unserved(
                MODULE,
                &format!("{accession} CDS annotation"),
                "the reference serves this accession's bases but no cdot record, so it has no \
                 CDS origin and a `c.` position does not resolve",
            );
            None
        }
    }
}

/// The transcript's cached bases, or `None` after recording the stand-down.
fn bases_or_skip<'a>(transcript: &'a Transcript, accession: &str) -> Option<&'a [u8]> {
    match transcript.sequence.as_deref() {
        Some(sequence) => Some(sequence.as_bytes()),
        None => {
            crate::common::manifest::unserved(
                MODULE,
                &format!("{accession} cached bases"),
                "the reference serves this accession's metadata but not its sequence",
            );
            None
        }
    }
}

/// Fail unless normalization left the bases the description denotes untouched.
///
/// Asserts on the **bases** rather than on the output string, so it keeps its
/// meaning if a canonical spelling is ever revised: the property is that
/// normalization does not change what the description means, which is rule 1 of
/// `README.md`'s ruleset and is not open to a rule 6 choice. That makes it the
/// assertion the guards below can carry without an operator settling which of
/// two conformant spellings ferro should ship.
///
/// A `NotComparable` verdict is a **panic**, not a pass. This module exists
/// because a guard that declines to check reads exactly like a guard that
/// checked and agreed (#1858), and folding the skips in with the agreements
/// would rebuild that hole inside the very helper meant to close it.
fn assert_denotes_the_same_bases(
    provider: &MultiFastaProvider,
    input: &HgvsVariant,
    output: &HgvsVariant,
) {
    match compare_denoted_sequences(input, output, provider) {
        DenotedSequenceComparison::Agree => {}
        DenotedSequenceComparison::Differ {
            accession,
            start,
            reference,
            from_input,
            from_output,
        } => panic!(
            "normalization changed the sequence the description denotes\n  input     {input}\n  \
             output    {output}\n  over {accession} from {start}\n  reference   {reference}\n  \
             from input  {from_input}\n  from output {from_output}"
        ),
        other => panic!(
            "expected a base-level verdict, got {other:?}\n  input  {input}\n  output {output}"
        ),
    }
}

#[test]
fn test_nm_001408491_c517dela_should_shift() {
    // This test reproduces the issue found in mutalyzer comparison:
    // NM_001408491.1:c.517delA should normalize to c.518del (3' rule)
    // because c.517 and c.518 are both 'A'
    //
    // Mutalyzer output: NM_001408491.1:c.518del
    // Ferro output (before fix): NM_001408491.1:c.517del

    let Some(provider) = get_provider() else {
        crate::common::manifest::absent(MODULE);
        return;
    };

    // First, let's get the transcript and examine the sequence
    let Some(transcript) = transcript_or_skip(&provider, "NM_001408491.1") else {
        return;
    };

    println!("=== Transcript NM_001408491.1 ===");
    println!("Sequence length: {}", transcript.sequence_length());
    println!("CDS start (1-based): {:?}", transcript.cds_start);
    println!("CDS end: {:?}", transcript.cds_end);
    println!("Exon count: {}", transcript.exons.len());
    for (i, exon) in transcript.exons.iter().enumerate() {
        println!(
            "  Exon {}: {} - {} (number: {})",
            i + 1,
            exon.start,
            exon.end,
            exon.number
        );
    }

    // Calculate transcript positions for c.517 and c.518
    let Some(cds_start) = cds_start_or_skip(&transcript, "NM_001408491.1") else {
        return;
    };
    let tx_pos_517 = cds_start + 517 - 1; // c.517 -> tx position
    let tx_pos_518 = cds_start + 518 - 1; // c.518 -> tx position

    println!("\nCoordinate mapping:");
    println!("  CDS start: {}", cds_start);
    println!(
        "  c.517 -> tx position {} (0-based: {})",
        tx_pos_517,
        tx_pos_517 - 1
    );
    println!(
        "  c.518 -> tx position {} (0-based: {})",
        tx_pos_518,
        tx_pos_518 - 1
    );

    // Get the bases at these positions
    let Some(seq) = bases_or_skip(&transcript, "NM_001408491.1") else {
        return;
    };
    let idx_517 = (tx_pos_517 - 1) as usize; // 0-based
    let idx_518 = (tx_pos_518 - 1) as usize; // 0-based

    // Asserted rather than tested, for the reason given on the same bound in
    // `test_compare_with_mutalyzer_sequence`: skipping the premise check when
    // the premise is unavailable is how a guard passes without checking.
    assert!(
        idx_517 < seq.len() && idx_518 < seq.len(),
        "c.517/c.518 are past the end of NM_001408491.1 ({} bases): the reference no longer \
         carries the coordinates this guard is about",
        seq.len()
    );

    println!("\nSequence around c.517-518:");
    // Show context
    let start = idx_517.saturating_sub(3);
    let end = (idx_518 + 4).min(seq.len());
    let context: String = seq[start..end].iter().map(|&b| b as char).collect();
    println!("  ... {} ...", context);
    println!(
        "  c.517 (tx pos {}, 0-based {}): {}",
        tx_pos_517, idx_517, seq[idx_517] as char
    );
    println!(
        "  c.518 (tx pos {}, 0-based {}): {}",
        tx_pos_518, idx_518, seq[idx_518] as char
    );

    // The premise the 3' shift rests on: the two bases must be equal, or
    // moving the deletion between them would change what it denotes. Checked
    // here rather than printed, so a reference that stopped satisfying it
    // could not leave the assertion below passing for the wrong reason.
    assert_eq!(
        seq[idx_517], seq[idx_518],
        "c.517 ({}) and c.518 ({}) differ, so no 3' shift between them is sequence-preserving \
         and the expectation below no longer follows from the reference",
        seq[idx_517] as char, seq[idx_518] as char
    );

    // Find which exon contains c.517
    println!("\nExon containing c.517 (tx pos {}):", tx_pos_517);
    for exon in &transcript.exons {
        if tx_pos_517 >= exon.start && tx_pos_517 <= exon.end {
            println!("  Exon {} (tx {}-{})", exon.number, exon.start, exon.end);
            // Check if c.518 is also in this exon
            if tx_pos_518 <= exon.end {
                println!("  c.518 is also in this exon - no boundary issues");
            } else {
                println!("  WARNING: c.518 is NOT in this exon - boundary may prevent shift!");
            }
        }
    }

    // Now test normalization
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_001408491.1:c.517delA").unwrap();
    let result = normalizer.normalize(&variant).unwrap();
    let output = format!("{}", result);

    println!("\n=== Normalization Result ===");
    println!("Input:  NM_001408491.1:c.517delA");
    println!("Output: {}", output);

    // The test expectation
    assert!(
        output.contains("c.518del"),
        "Expected c.517delA to normalize to c.518del (3' rule), got: {}",
        output
    );

    // …and that the re-spelling is sequence-preserving. The assertion above
    // fixes the string; this one fixes its meaning, so a future change that
    // moved the deletion to c.518 by some route that also changed the bases
    // could not satisfy the guard.
    assert_denotes_the_same_bases(normalizer.provider(), &variant, &result);
}

/// The guard's original question — "ferro outputs `c.1324del`, Mutalyzer
/// outputs `c.1326del`, one of these is wrong about 3' shifting" — **is
/// answered: they agree, and the disagreement the comment recorded no longer
/// occurs.**
///
/// Measured 2026-08-14 against both prepared volumes, `NM_033517.1:c.1324del`
/// normalizes to **`NM_033517.1:c.1326del`** — the Mutalyzer answer. It is a
/// real 3' shift and the reference says why: `c.1324`, `c.1325` and `c.1326`
/// are all `C`, `c.1327` is `G`, so the deleted base moves to the 3' end of a
/// `CCC` run and stops there. Four further deletions on this accession shift in
/// the same run — `c.125del` -> `c.127del`, `c.160del` -> `c.162del`,
/// `c.164del` -> `c.167del`, `c.517del` -> `c.518del` — so the 3' rule is
/// plainly reaching this transcript.
///
/// **An earlier version of this comment said the opposite, and the retraction
/// is left in place on purpose.** It claimed the reference serves this
/// accession as bare FASTA with `cds_start: None`, `cds_end: None` and zero
/// exons, and that `c.1324del` therefore came back unchanged as a passthrough
/// rather than a 3'-rule decision. None of that reproduces: the module header
/// carries the measured annotation for all four accessions, `cds_start` is
/// `Some(1)`, and every input listed above shifts. So there is no CDS-less
/// state here, no fabricated coordinate, and no ferro-versus-Mutalyzer
/// disagreement — and correspondingly no open question for an operator about
/// whether `ferro prepare` should serve `NM_033517` with a cdot record. It
/// does.
///
/// **Why the guard still asserts sequence preservation rather than pinning
/// `c.1326del`.** Nothing about the measurement forbids the pin — it would be
/// stable across the cdot versions in play, since `cds_start` is 1 under both
/// and the shift is decided by bases rather than by the `cds_end` those
/// versions disagree on. It is not added because this module's policy (see the
/// header) is that a guard here carries the rule 1 property, which needs no
/// adjudication, and leaves the choice of spelling to whoever adjudicates it;
/// `assert_denotes_the_same_bases` fails for any shift out of the `CCC` run,
/// which is how a 3'-rule defect actually presents. Pinning the string is a
/// separate, defensible change and should be argued as one.
///
/// `#1619` is the neighbouring work on this accession, but it runs on
/// `MockProvider` and a committed fixture, so it says nothing about what the
/// prepared reference serves.
#[test]
fn test_potential_bug_deletion_shift_nm033517() {
    let Some(provider) = get_provider() else {
        crate::common::manifest::absent(MODULE);
        return;
    };

    let Some(transcript) = transcript_or_skip(&provider, "NM_033517.1") else {
        return;
    };

    // Resolves to 1 on the blessed reference — this does not stand down, and an
    // earlier version of this comment claiming it did is withdrawn (see the doc
    // comment above). The call stays because the guard must not substitute an
    // origin the reference did not give it, whatever today's reference gives.
    let Some(cds_start) = cds_start_or_skip(&transcript, "NM_033517.1") else {
        return;
    };
    let Some(seq) = bases_or_skip(&transcript, "NM_033517.1") else {
        return;
    };

    // Check c.1324, c.1325, c.1326
    for pos in 1324..=1328 {
        let tx_pos = cds_start + pos - 1;
        let idx = (tx_pos - 1) as usize;
        if idx < seq.len() {
            println!(
                "c.{} (tx {}, idx {}): {}",
                pos, tx_pos, idx, seq[idx] as char
            );
        }
    }

    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_033517.1:c.1324del").unwrap();
    let result = normalizer.normalize(&variant).unwrap();
    println!("\nInput:  NM_033517.1:c.1324del");
    println!("Output: {}", result);

    // What is assertable without settling either open question above: whichever
    // position ferro names, it must denote the bases the input denotes. A 3'
    // shift within a run is sequence-preserving by construction, so this passes
    // for `c.1324del` and for `c.1326del` alike and takes no side between them
    // — while a shift that landed outside the run, which is the way a 3'-rule
    // defect actually goes wrong, fails it.
    assert_denotes_the_same_bases(normalizer.provider(), &variant, &result);
}

/// Whether ferro may re-spell `NM_001282424.3:c.2139_2140delinsTATGCA`.
///
/// **The guard's original premise is withdrawn.** It read "HGVS spec: delins
/// should NOT be 3' shifted like del/dup", which grades a description by its
/// *edit type* — a property of the input's spelling. `README.md` rule 3 says
/// every rule is evaluated over the resulting sequence "never over the input's
/// spelling", and `canonical-form-choice-when-both-legal` (decided 2026-08-07)
/// says ferro derives from the resulting sequence and emits what falls out. So
/// "delins does not shift" is not a rule this project holds, and the guard
/// cannot be repaired by asserting it.
///
/// **What it was really written to catch is still real, and is now asserted.**
/// The feared output, `c.2140_2141delinsTATGCA`, does not denote the input's
/// bases. Over the measured reference (`c.2138` = A, `c.2139` = G, `c.2140` =
/// C, `c.2141` = T) the input yields `A TATGCA T` while that spelling yields
/// `A G TATGCA T` — a different sequence, so it would be a rule 1 violation
/// rather than a debatable choice. [`assert_denotes_the_same_bases`] is exactly
/// the test that separates the two, and it needs no ruling to do it.
///
/// Measured 2026-08-13, ferro emits neither form: it returns
/// `c.[2138_2139insTAT;2140_2141insA]`, a two-insertion split that leaves the
/// reference's `GC` at `c.2139_2140` in place and is sequence-preserving. That
/// is what re-derivation produces here, and it is the form the ledger predicts
/// — this is a **net insertion** (6 payload bases for a 2-base span), and
/// `delins-merge-vs-individual-gap-two-or-more` scopes the merge ruling to the
/// net-deletion case, leaving `general.md:34`'s individual description
/// governing. The exact spelling is deliberately **not** pinned here: it is a
/// rule 6 choice, and this module is not where that gets decided.
///
/// The old body ended by printing "POTENTIAL BUG: Ferro shifted delins from
/// 2139_2140 to 2140_2141" whenever the output contained `2140_2141`. That
/// substring now matches the `insA` member of a correct split, so the warning
/// fired on conforming output — a detector that had gone wrong in both
/// directions at once, reporting its finding by passing and reporting the wrong
/// finding as well.
#[test]
fn test_potential_bug_delins_shift() {
    let Some(provider) = get_provider() else {
        crate::common::manifest::absent(MODULE);
        return;
    };

    let Some(transcript) = transcript_or_skip(&provider, "NM_001282424.3") else {
        return;
    };

    let Some(cds_start) = cds_start_or_skip(&transcript, "NM_001282424.3") else {
        return;
    };
    let Some(seq) = bases_or_skip(&transcript, "NM_001282424.3") else {
        return;
    };

    // Check sequence around c.2139-2141
    println!("Sequence around c.2139-2141:");
    for pos in 2137..=2145 {
        let tx_pos = cds_start + pos - 1;
        let idx = (tx_pos - 1) as usize;
        if idx < seq.len() {
            println!(
                "c.{} (tx {}, idx {}): {}",
                pos, tx_pos, idx, seq[idx] as char
            );
        }
    }

    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_001282424.3:c.2139_2140delinsTATGCA").unwrap();
    let result = normalizer.normalize(&variant).unwrap();
    println!("\nInput:  NM_001282424.3:c.2139_2140delinsTATGCA");
    println!("Output: {}", result);

    // The assertion the printed warning should always have been: not "did the
    // span move" but "does the output still mean what the input meant". It
    // fires on `c.2140_2141delinsTATGCA` — the shift the guard was named for —
    // and stays silent for any re-spelling that denotes the same bases.
    assert_denotes_the_same_bases(normalizer.provider(), &variant, &result);
}

/// **This input cannot answer the question the guard was written to ask, and
/// that is now asserted rather than left as a note.**
///
/// The question was "does the 3' rule reach the 5'UTR (negative coordinates)?",
/// posed because ferro once emitted `c.-55_-46dup` where Mutalyzer emitted
/// `c.-56_-47dup`. It is a real question. This input does not discriminate it,
/// because over the measured reference the duplicated unit cannot rotate 3' at
/// all: a duplication of `[a, b]` may be re-spelled `[a+1, b+1]` only when
/// `seq[a] == seq[b+1]`, and here `c.-56` is `T` while `c.-46` is `G`. Both
/// readings of the 5'UTR question therefore predict the same output, so a green
/// run is evidence about neither. The old body said the answer "depends on HGVS
/// spec interpretation" and asserted nothing; the non-discrimination is
/// checkable, so it is checked.
///
/// **What the guard does still cover, and why it is worth keeping.** The 1-base
/// difference it was filed for was never a 3'-rule ruling — it was an off-by-one
/// in the 5'UTR CDS-to-transcript mapping, fixed in #97. So this is a
/// coordinate-mapping regression guard, and the shift it must not see is
/// arithmetic rather than adjudicated. The expected string is not invented
/// here: `normalize_tests::clinvar_normalization` already pins this exact input
/// to this exact output hermetically, and this guard re-checks it against the
/// real prepared reference — a different provider, a different coordinate path,
/// the same answer.
#[test]
fn test_5utr_duplication_shifting() {
    let Some(provider) = get_provider() else {
        crate::common::manifest::absent(MODULE);
        return;
    };

    let Some(transcript) = transcript_or_skip(&provider, "NM_001394148.2") else {
        return;
    };

    let Some(cds_start) = cds_start_or_skip(&transcript, "NM_001394148.2") else {
        return;
    };
    let Some(seq) = bases_or_skip(&transcript, "NM_001394148.2") else {
        return;
    };

    // Check sequence around c.-56 to c.-45
    println!("CDS start: {}", cds_start);
    println!("Sequence around c.-56 to c.-45:");
    for offset in -60i64..=-40 {
        // c.-N means position cds_start - N in 1-based
        let tx_pos = (cds_start as i64 + offset) as u64;
        if tx_pos >= 1 {
            let idx = (tx_pos - 1) as usize;
            if idx < seq.len() {
                println!(
                    "c.{} (tx {}, idx {}): {}",
                    offset, tx_pos, idx, seq[idx] as char
                );
            }
        }
    }

    // `c.-N` resolves to transcript position `cds_start + (-N)`, so the
    // duplicated unit c.-56_-47 occupies these two 0-based bounds.
    let idx_first = (cds_start as i64 - 56 - 1) as usize;
    let idx_last = (cds_start as i64 - 47 - 1) as usize;
    assert!(
        idx_last + 1 < seq.len(),
        "c.-47 and its 3' neighbour are past the end of NM_001394148.2 ({} bases)",
        seq.len()
    );

    // The non-discrimination, asserted. If this ever fails, the reference has
    // changed such that a 3' rotation *is* available, this input starts telling
    // the two readings apart, and the expectation below must be re-derived
    // rather than kept — so the failure is the signal to reopen the question,
    // not a defect to route around.
    assert_ne!(
        seq[idx_first],
        seq[idx_last + 1],
        "c.-56 ({}) and c.-46 ({}) are equal, so the duplicated unit CAN be re-spelled one base \
         3' and this input now discriminates whether the 3' rule reaches the 5'UTR — a question \
         this guard has never answered. Re-derive the expectation; do not update it to match.",
        seq[idx_first] as char,
        seq[idx_last + 1] as char
    );

    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_001394148.2:c.-56_-47dup").unwrap();
    let result = normalizer.normalize(&variant).unwrap();
    println!("\nInput:  NM_001394148.2:c.-56_-47dup");
    println!("Output: {}", result);

    // No rotation exists, so the input is already the only spelling of this
    // duplication — which makes this a pin on the 5'UTR coordinate arithmetic
    // (#97), not on a normalization preference. It cross-checks the hermetic
    // pin in `normalize_tests::clinvar_normalization` against the real
    // reference.
    assert_eq!(
        result.to_string(),
        "NM_001394148.2:c.-56_-47dup",
        "the 5'UTR duplication moved, but the assertion above says no sequence-preserving \
         re-spelling exists — so this is the #97 CDS-to-transcript off-by-one, not a 3'-rule \
         choice"
    );

    assert_denotes_the_same_bases(normalizer.provider(), &variant, &result);
}

#[test]
fn test_compare_with_mutalyzer_sequence() {
    // Compare what we have in cdot vs what mutalyzer reports
    let Some(provider) = get_provider() else {
        crate::common::manifest::absent(MODULE);
        return;
    };

    let Some(transcript) = transcript_or_skip(&provider, "NM_001408491.1") else {
        return;
    };

    // From mutalyzer normalizer output for NM_001408491.1:c.517delA:
    // - Output: NM_001408491.1:c.518del
    // - This means the deleted base shifts from c.517 to c.518
    // - This is only correct if c.517 == c.518 (both are the same base)

    let Some(cds_start) = cds_start_or_skip(&transcript, "NM_001408491.1") else {
        return;
    };
    let tx_pos_517 = cds_start + 517 - 1;
    let tx_pos_518 = cds_start + 518 - 1;

    let Some(seq) = bases_or_skip(&transcript, "NM_001408491.1") else {
        return;
    };
    let idx_517 = (tx_pos_517 - 1) as usize;
    let idx_518 = (tx_pos_518 - 1) as usize;

    // This was an `if`, and the two assertions below were its body — so the one
    // condition it exists to detect, a coordinate the reference cannot serve,
    // used to leave the whole test asserting nothing and passing (#1858).
    //
    // On the blessed reference the bound is never false, and measurably so:
    // `NM_001408491.1` is 3748 bases with `cds_start` 296, putting these two at
    // 0-based 811 and 812 — inside by a factor of four. So the change is not
    // fixing a live failure; it is removing the one way this guard could go
    // quiet, which is precisely the shape that cannot be seen from a run's
    // duration or its pass/fail.
    assert!(
        idx_517 < seq.len() && idx_518 < seq.len(),
        "c.517 (0-based {idx_517}) and c.518 (0-based {idx_518}) must both be inside \
         NM_001408491.1, which is {} bases: a reference that no longer serves them cannot \
         support the comparison this guard makes, and must not pass it by skipping",
        seq.len()
    );

    let base_517 = seq[idx_517] as char;
    let base_518 = seq[idx_518] as char;

    println!("cdot sequence at c.517: {}", base_517);
    println!("cdot sequence at c.518: {}", base_518);

    // For the mutalyzer result to be correct, both must be 'A'
    assert_eq!(
        base_517, 'A',
        "Expected c.517 to be 'A' per mutalyzer result"
    );
    assert_eq!(
        base_518, 'A',
        "Expected c.518 to be 'A' for 3' shift to be valid"
    );
}
