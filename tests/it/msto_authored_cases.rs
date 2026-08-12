//! Authored test cases harvested out of the issue tracker into runnable code.
//!
//! # Why this module exists
//!
//! A large body of concrete HGVS cases exists only as prose in issue bodies:
//! spec `✓`/`✗` pairs, diagnostic probe tables, worked examples. Filed is not
//! the same as guarded — an issue records that something was once measured, but
//! nothing re-checks it, so behaviour drifts silently and the case has to be
//! re-derived by hand every time someone needs it. Several of the rows below
//! turned out to have been **fixed years ago** and several are **still open**,
//! and there was no way to tell which without running them.
//!
//! # These are characterizations, not aspirations
//!
//! Every row asserts **what ferro does today**, so this module is green on
//! `main` and stays green. Rows are classified:
//!
//! * [`Verdict::Guard`] — current behaviour matches what the issue (or the spec
//!   line it cites) says it should be. Changing it is a regression.
//! * [`Verdict::Gap`] — current behaviour differs from what the issue argues
//!   for. The assertion still pins today's answer; the comment names the
//!   expected one. **A `Gap` flipping to `Guard` is progress** — update the row
//!   and the census, deliberately, in the PR that fixes it.
//!
//! Pinning the wrong answer on purpose is the point: it makes the wrong answer
//! *visible in a diff* the day somebody changes it, which prose in a closed
//! issue cannot do.
//!
//! # Provenance
//!
//! Every row is traceable to a specific issue. Nothing here is invented. The
//! `TEMPLATE` contig is reproduced verbatim from the #1419/#1420/#1421
//! reproducers; the #179 contig is reconstructed (see [`issue_179_reference`]).
//!
//! # Overlap with `reported_confluence_pairs`
//!
//! The nine headline **pairs** of #1419/#1420/#1421 are guarded elsewhere as a
//! confluence axis. This module deliberately does not repeat them: it harvests
//! the material those pairs left behind — #1421's five-row diagnostic table,
//! which is about a *single* spelling's verdict rather than a pair, and which
//! is the evidence that the split decision keys on the inserted block rather
//! than on the rule.

use ferro_hgvs::reference::mock::JsonProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};
use std::io::Write;

/// Whether current behaviour matches what the source issue argues for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Verdict {
    /// ferro agrees with the issue / the spec line it cites.
    Guard,
    /// ferro disagrees; the row pins today's answer and names the wanted one.
    Gap,
}

// ---------------------------------------------------------------------------
// #1421 — does a net-insertion spanning delins split across unchanged interior
// bases? The answer depends on the inserted block, which is the defect.
// ---------------------------------------------------------------------------

/// The 125 nt contig from the #1419/#1420/#1421 reproducers, verbatim.
///
/// Synthetic: it names no real sequence.
const TEMPLATE: &str = concat!(
    "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCG",
    "TACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA"
);

/// #1421's "Diagnostic note" table, reproduced exactly.
///
/// Every row is the **same locus** (29–33, interior `AC` at 30–31 unchanged)
/// with only the inserted block varied. `delins.md:17` asks for the separated
/// form in all five — the separation is two nucleotides, and `delins.md:18`'s
/// one-nucleotide exception cannot apply to a `g.` description. ferro split
/// three and retained two.
///
/// **All five now split.** What retained the last two was
/// `separations_are_meaningful`'s raise to `RAISED_PIECE_SEPARATION`, which
/// disbelieves a one-base separation inside a large net change. That raise is
/// `DNA/delins.md:44-47`'s payload-coincidence carve-out, and the operator
/// ruling `delins-payload-coincidence-carve-out-is-coding-dna-scoped` scopes it
/// to the coding DNA axis. Every row here is `TEMPLATE:g.`, so the carve-out
/// does not reach them and `delins.md:17` gets the answer this table says it
/// asks for.
///
/// `(label, input, current output, splits)`
const ISSUE_1421_INSERT_BLOCKS: &[(&str, &str, &str, bool)] = &[
    (
        "net+2",
        "TEMPLATE:g.29_33delinsAACTGTG",
        "TEMPLATE:g.[29C>A;32_33delinsTGTG]",
        true,
    ),
    (
        "net+3",
        "TEMPLATE:g.29_33delinsAACTGCTG",
        "TEMPLATE:g.[29C>A;32_33delinsTGCTG]",
        true,
    ),
    (
        "net+6-split",
        "TEMPLATE:g.29_33delinsAACTGCATGTG",
        "TEMPLATE:g.[29C>A;32_33delinsTGCATGTG]",
        true,
    ),
    // The headline case of #1421. Same locus, same net length as the row above,
    // and it used to reach the OPPOSITE verdict — which was the whole point of
    // the table, and is what scoping the carve-out to `c.` removed. The two
    // labels are kept as `-retained` so the table still reads against the
    // issue's own text; what they retained is history, not behaviour.
    (
        "net+6-retained",
        "TEMPLATE:g.29_33delinsAACACATACTG",
        "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        true,
    ),
    (
        "net+7-retained",
        "TEMPLATE:g.29_33delinsAACTGCATGCTG",
        "TEMPLATE:g.[29C>A;32_33delinsTGCATGCTG]",
        true,
    ),
];

/// How many of the five split today. `delins.md:17` wants all five, and all
/// five now do.
///
/// **This may only ever go up.** A drop means a spelling that once obeyed the
/// separated-description rule stopped. At 5 of 5 the only way left to move it is
/// down, which is exactly the regression this constant exists to catch.
const ISSUE_1421_ROWS_THAT_SPLIT: usize = 5;

fn provider_for(sequence: &str, accession: &str) -> JsonProvider {
    let n = sequence.len() as u64;
    let json = serde_json::json!({
        "version": "1.0",
        "transcripts": [{
            "id": format!("{accession}-gene.1"),
            "chromosome": accession,
            "strand": "+",
            "sequence": sequence,
            "cds_start": 1,
            "cds_end": n - (n % 3),
            "genomic_start": 1,
            "genomic_end": n,
            "protein_id": format!("{accession}-gene.1"),
            "exons": [{
                "number": 1, "start": 1, "end": n,
                "genomic_start": 1, "genomic_end": n
            }]
        }],
        "genomic_sequences": { accession: sequence }
    });
    let mut f = tempfile::NamedTempFile::new().expect("tempfile");
    write!(f, "{json}").expect("write json");
    JsonProvider::from_json(f.path()).expect("load reference")
}

fn normalize_with(normalizer: &Normalizer<JsonProvider>, input: &str) -> String {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"));
    format!("{normalized}")
}

/// Each row of #1421's diagnostic table still produces what it produced when
/// the table was written.
#[test]
fn issue_1421_insert_block_table_is_unchanged() {
    let normalizer = Normalizer::new(provider_for(TEMPLATE, "TEMPLATE"));
    for (label, input, expected, splits) in ISSUE_1421_INSERT_BLOCKS {
        let got = normalize_with(&normalizer, input);
        assert_eq!(&got, expected, "#1421 row `{label}`: {input}");
        // The recorded flag is checked against behaviour rather than left as
        // prose, so a row whose `expected` is edited without its `splits` — or
        // the reverse — fails here instead of silently disagreeing with the
        // census below.
        assert_eq!(
            got.contains(';'),
            *splits,
            "#1421 row `{label}`: recorded `splits` = {splits} disagrees with `{got}`"
        );
    }
}

/// The defect itself, stated as an assertion rather than as prose: two rows at
/// one locus, identical net length, opposite split verdicts.
///
/// Whether a spanning `delins` is separated must be a function of which
/// positions changed — not of which bases the payload happens to carry. While
/// this test passes, it is not.
#[test]
fn issue_1421_split_decision_no_longer_depends_on_the_inserted_block() {
    let normalizer = Normalizer::new(provider_for(TEMPLATE, "TEMPLATE"));
    let split = normalize_with(&normalizer, "TEMPLATE:g.29_33delinsAACTGCATGTG");
    let retained = normalize_with(&normalizer, "TEMPLATE:g.29_33delinsAACACATACTG");

    // #1421's complaint, in one assertion: two inputs at one locus with one net
    // length reached opposite verdicts, so the verdict was a property of the
    // payload's bases rather than of the locus. On a frameless axis it no longer
    // is — the payload-coincidence carve-out that read those bases is scoped to
    // `c.` (`delins-payload-coincidence-carve-out-is-coding-dna-scoped`), so
    // `general.md:34` decides both the same way.
    //
    // Asserted as a pair rather than as two independent `contains(';')` checks:
    // what the issue is about is the two rows AGREEING, and two separate checks
    // would both keep passing if a later change split them apart again in the
    // other direction.
    assert_eq!(
        (split.contains(';'), retained.contains(';')),
        (true, true),
        "the two net+6 rows must reach the same verdict; got `{split}` and `{retained}`"
    );
}

/// Census. Free to rise, never to fall silently.
#[test]
fn issue_1421_split_census_holds() {
    let normalizer = Normalizer::new(provider_for(TEMPLATE, "TEMPLATE"));
    let splits = ISSUE_1421_INSERT_BLOCKS
        .iter()
        .filter(|(_, input, _, _)| normalize_with(&normalizer, input).contains(';'))
        .count();
    assert_eq!(
        splits, ISSUE_1421_ROWS_THAT_SPLIT,
        "the number of #1421 rows obeying delins.md:17 moved; if it rose, raise \
         ISSUE_1421_ROWS_THAT_SPLIT and update the affected rows"
    );
}

// ---------------------------------------------------------------------------
// #87 — v21 spec-compliance ✓/✗ pairs. Parser-level, so no reference needed.
// ---------------------------------------------------------------------------

/// What the parser does with a string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Outcome {
    /// Parses, and displays as the given string (which may differ from the
    /// input — a silent rewrite is itself behaviour worth pinning).
    Parses(&'static str),
    /// Rejected.
    Rejected,
}

/// `(input, outcome, verdict, spec line, note)`.
///
/// #87 was closed as subsumed by #81/#91/#92. Running its checklist shows most
/// items did land — with spec-citing rejection messages — and three did not.
/// Those three are marked [`Verdict::Gap`] and are the reason this table is
/// worth keeping rather than deleting along with the issue.
const ISSUE_87_SPEC_CASES: &[(&str, Outcome, Verdict, &str, &str)] = &[
    // --- Insertion form: two adjacent flanking positions -------------------
    (
        "NM_000000.1:c.123_124insG",
        Outcome::Parses("NM_000000.1:c.123_124insG"),
        Verdict::Guard,
        "DNA/insertion.md:15",
        "the ✓ form",
    ),
    (
        "NC_000000.1:g.32867861_32867862insT",
        Outcome::Parses("NC_000000.1:g.32867861_32867862insT"),
        Verdict::Guard,
        "DNA/insertion.md:15",
        "the ✓ form",
    ),
    (
        "NM_000000.1:c.123insG",
        Outcome::Rejected,
        Verdict::Guard,
        "DNA/insertion.md:96-101",
        "single-position insertion; #87 recorded this as parsing and then panicking \
         in shuffle.rs on an underflow. Now rejected at parse with a spec citation.",
    ),
    (
        "NM_000000.1:c.123_125insATG",
        Outcome::Rejected,
        Verdict::Guard,
        "DNA/insertion.md:15",
        "non-adjacent anchor; was accepted when #87 was filed",
    ),
    // --- Alleles §54: compact homozygous forms ------------------------------
    (
        "NM_000000.1:c.2376[G>C];[G>C]",
        Outcome::Rejected,
        Verdict::Guard,
        "DNA/alleles.md:54",
        "compact homozygous shorthand is not allowed",
    ),
    (
        "NM_000000.1:c.2376G>C[];[]",
        Outcome::Rejected,
        Verdict::Guard,
        "DNA/alleles.md:54",
        "the even shorter shorthand, likewise not allowed",
    ),
    (
        "NM_000000.1:c.[2376G>C];[2376G>C]",
        Outcome::Parses("NM_000000.1:c.[2376G>C];[2376G>C]"),
        Verdict::Guard,
        "DNA/alleles.md:54",
        "the ✓ long form",
    ),
    // --- Alleles §60-61: bare-`=` allele ------------------------------------
    (
        "NM_000000.1:c.[2376G>C];[=]",
        Outcome::Parses("NM_000000.1:c.[2376G>C];[=]"),
        Verdict::Guard,
        "DNA/alleles.md:60-61",
        "whole-reference-WT semantic; #87 recorded this as NOT parsing. It does now.",
    ),
    (
        "NM_000000.1:c.[2376G>C];[2376=]",
        Outcome::Parses("NM_000000.1:c.[2376G>C];[2376=]"),
        Verdict::Guard,
        "DNA/alleles.md:60-61",
        "position-specific `=`, a distinct semantic from the bare form above",
    ),
    // --- Alleles §19/§49: redundant `=` in both alleles ---------------------
    (
        "NM_000000.1:c.[2376G>C;3103=];[2376=;3103del]",
        Outcome::Rejected,
        Verdict::Guard,
        "DNA/alleles.md:19,49",
        "the ✗ form. Rejected — though by the grammar declining `=` in cis with a \
         non-`=` edit, not by a named rule, so it is compliant here by construction \
         rather than by decision.",
    ),
    // --- Alleles §106: c.0 accession placement ------------------------------
    (
        "LRG_199t1:c.[76A>C];[0]",
        Outcome::Parses("[LRG_199t1:c.76A>C];[0]"),
        Verdict::Gap,
        "DNA/alleles.md:106",
        "STILL OPEN. The spec writes the accession OUTSIDE the bracket group — \
         `LRG_199t1:c.[76A>C];[0]` — but Display hoists it inside the first bracket. \
         Round-trips, but not to the spec's own spelling.",
    ),
    // --- Inversion: span > 1 nt ---------------------------------------------
    (
        "NC_000000.1:g.234inv",
        Outcome::Rejected,
        Verdict::Guard,
        "DNA/inversion.md:16",
        "one-nucleotide inversion; was accepted when #87 was filed",
    ),
    (
        "NM_000000.1:c.20inv",
        Outcome::Rejected,
        Verdict::Guard,
        "DNA/inversion.md:16",
        "same, on the coding axis",
    ),
    (
        "NC_000000.1:g.234_235inv",
        Outcome::Parses("NC_000000.1:g.234_235inv"),
        Verdict::Guard,
        "DNA/inversion.md:16",
        "the ✓ form",
    ),
    (
        "NM_000000.1:c.20_22inv",
        Outcome::Parses("NM_000000.1:c.20_22inv"),
        Verdict::Guard,
        "DNA/inversion.md:16",
        "the ✓ form",
    ),
    // --- Translocations: delins, not `::` ------------------------------------
    (
        "NC_000002.12:g.pter_8247756::NC_000011.10:g.15825273_cen_qter",
        Outcome::Rejected,
        Verdict::Guard,
        "DNA/complex.md:60-72",
        "`::` between two genomic accessions in a translocation shape",
    ),
    (
        "NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825272]",
        Outcome::Parses("NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825272]"),
        Verdict::Guard,
        "DNA/complex.md:60-72",
        "the ✓ delins form",
    ),
    // --- Range completeness --------------------------------------------------
    (
        "NM_000000.1:c.123-65_-50",
        Outcome::Parses("NM_000000.1:c.123-65_-50"),
        Verdict::Gap,
        "checklist.md:26",
        "STILL OPEN. The spec calls this incomplete — the second endpoint must \
         repeat its anchor (`c.123-65_123-50`). ferro accepts it, and the shorthand \
         means something else entirely: a 5' UTR position.",
    ),
    (
        "NM_000000.1:c.123-65_123-50",
        Outcome::Parses("NM_000000.1:c.123-65_123-50"),
        Verdict::Guard,
        "checklist.md:26",
        "the ✓ fully-specified form",
    ),
    // --- RNA uncertainty forms -----------------------------------------------
    (
        "NM_000000.1:r.0",
        Outcome::Parses("NM_000000.1:r.0"),
        Verdict::Guard,
        "uncertain.md:184",
        "no transcript",
    ),
    (
        "NM_000000.1:r.0?",
        Outcome::Rejected,
        Verdict::Gap,
        "uncertain.md:184",
        "STILL OPEN. `r.0?` (possibly no transcript) is a distinct spec form from \
         `r.0` (no transcript) and does not parse.",
    ),
    (
        "NM_000000.1:r.=",
        Outcome::Parses("NM_000000.1:r.="),
        Verdict::Guard,
        "uncertain.md:187",
        "observed no change",
    ),
    (
        "NM_000000.1:r.(=)",
        Outcome::Parses("NM_000000.1:r.(=)"),
        Verdict::Guard,
        "uncertain.md:187",
        "predicted no change; #87 recorded this as NOT parsing. It does now.",
    ),
    (
        "NM_000000.1:r.spl",
        Outcome::Parses("NM_000000.1:r.spl"),
        Verdict::Guard,
        "uncertain.md:190-196",
        "#87 recorded this as NOT parsing. It does now.",
    ),
    (
        "NM_000000.1:r.spl?",
        Outcome::Parses("NM_000000.1:r.spl?"),
        Verdict::Guard,
        "uncertain.md:190-196",
        "#87 recorded this as NOT parsing. It does now.",
    ),
    // --- Repeats --------------------------------------------------------------
    (
        "NC_000000.1:g.112036755_112036823CTG[9]TTG[1]CTG[13]",
        Outcome::Parses("NC_000000.1:g.112036755_112036823CTG[9]TTG[1]CTG[13]"),
        Verdict::Guard,
        "DNA/repeated.md:18,85",
        "composite mixed-unit repeat",
    ),
    (
        "NM_000000.1:c.1210-33_1210-6GT[11]T[6]",
        Outcome::Parses("NM_000000.1:c.1210-33_1210-6GT[11]T[6]"),
        Verdict::Guard,
        "DNA/repeated.md:81",
        "intronic composite repeat",
    ),
    (
        "NM_000000.1:r.123cug[23]",
        Outcome::Parses("NM_000000.1:r.123cug[23]"),
        Verdict::Guard,
        "RNA/repeated.md:20-21",
        "sequence form — must be a pure repeat",
    ),
    (
        "NM_000000.1:r.123_125[23]",
        Outcome::Parses("NM_000000.1:r.123_125[23]"),
        Verdict::Guard,
        "RNA/repeated.md:20-21",
        "range form — may be interrupted. Distinct semantic from the row above; \
         both round-trip, so the two spellings stay distinguishable.",
    ),
    // --- Protein forms ---------------------------------------------------------
    (
        "NP_000000.1:p.[(Arg727Ser;Cys1334Trp)]",
        Outcome::Parses("NP_000000.1:p.[(Arg727Ser;Cys1334Trp)]"),
        Verdict::Guard,
        "protein/alleles.md:16,32-34",
        "predicted parens INSIDE the brackets; #87 recorded this as NOT parsing. \
         It does now.",
    ),
    (
        "NP_000000.1:p.([Arg727Ser;Cys1334Trp])",
        Outcome::Rejected,
        Verdict::Guard,
        "protein/alleles.md:16,32-34",
        "the ✗ inverted nesting",
    ),
    (
        "NP_000000.1:p.Trp10X",
        Outcome::Parses("NP_000000.1:p.Trp10Xaa"),
        Verdict::Guard,
        "protein/substitution.md:103-104",
        "#87 left this as an open policy choice between silent rewrite, warn, and \
         reject. The shipped answer is the silent rewrite — pinned here so the \
         choice is visible rather than implicit.",
    ),
    (
        "NP_000000.1:p.Gln151Thrfs*9",
        Outcome::Parses("NP_000000.1:p.Gln151ThrfsTer9"),
        Verdict::Guard,
        "protein/frameshift.md:19",
        "`*` normalized to `Ter` on display",
    ),
    (
        "NP_000000.1:p.His150Hisfs*10",
        Outcome::Rejected,
        Verdict::Guard,
        "protein/frameshift.md:47-49",
        "a frameshift anchored at an UNCHANGED residue; #87's headline ✗ example, \
         now rejected with a spec citation",
    ),
    (
        "NP_000000.1:p.Tyr4Ter",
        Outcome::Parses("NP_000000.1:p.Tyr4Ter"),
        Verdict::Guard,
        "protein/frameshift.md:22",
        "immediate stop is a nonsense substitution — the ✓ form",
    ),
    (
        "NP_000000.1:p.Tyr4TerfsTer1",
        Outcome::Rejected,
        Verdict::Guard,
        "protein/frameshift.md:22",
        "the same variant written as a frameshift — the ✗ form",
    ),
    (
        "NP_000000.1:p.=",
        Outcome::Parses("NP_000000.1:p.="),
        Verdict::Guard,
        "uncertain.md",
        "observed no change",
    ),
    (
        "NP_000000.1:p.(=)",
        Outcome::Parses("NP_000000.1:p.(=)"),
        Verdict::Guard,
        "uncertain.md",
        "predicted no change",
    ),
];

/// How many #87 rows still disagree with the spec line they cite.
///
/// **This may only ever go down.** A rise means a compliance item regressed.
const ISSUE_87_OPEN_GAPS: usize = 3;

/// Every #87 checklist row still parses (or is rejected) exactly as it does
/// today, and displays exactly as it does today.
#[test]
fn issue_87_spec_checklist_is_unchanged() {
    for (input, outcome, _, spec, note) in ISSUE_87_SPEC_CASES {
        match (parse_hgvs(input), outcome) {
            (Ok(variant), Outcome::Parses(display)) => assert_eq!(
                &format!("{variant}"),
                display,
                "`{input}` ({spec}) parsed but displayed differently — {note}"
            ),
            (Ok(variant), Outcome::Rejected) => panic!(
                "`{input}` ({spec}) is expected to be rejected but parsed as \
                 `{variant}` — if this is the fix, flip the row to Parses and \
                 lower ISSUE_87_OPEN_GAPS. {note}"
            ),
            (Err(_), Outcome::Rejected) => {}
            (Err(e), Outcome::Parses(display)) => panic!(
                "`{input}` ({spec}) is expected to parse as `{display}` but was \
                 rejected: {e} — {note}"
            ),
        }
    }
}

/// The three items #87 raised that are still not compliant.
///
/// Named individually rather than only counted, so a fix names itself in the
/// diff instead of just moving a number.
#[test]
fn issue_87_open_gaps_are_still_open() {
    // alleles.md:106 — the accession migrates inside the bracket group.
    let hoisted = parse_hgvs("LRG_199t1:c.[76A>C];[0]").expect("parses");
    assert_eq!(format!("{hoisted}"), "[LRG_199t1:c.76A>C];[0]");

    // checklist.md:26 — an incomplete intronic range is accepted.
    assert!(parse_hgvs("NM_000000.1:c.123-65_-50").is_ok());

    // uncertain.md:184 — `r.0?` is not supported.
    assert!(parse_hgvs("NM_000000.1:r.0?").is_err());
}

/// Census over the table, so the count cannot drift away from the rows.
#[test]
fn issue_87_gap_census_holds() {
    let gaps = ISSUE_87_SPEC_CASES
        .iter()
        .filter(|(_, _, verdict, _, _)| *verdict == Verdict::Gap)
        .count();
    assert_eq!(
        gaps, ISSUE_87_OPEN_GAPS,
        "the number of open #87 compliance gaps moved; if one was fixed, flip its \
         Verdict to Guard and lower ISSUE_87_OPEN_GAPS"
    );
}

// ---------------------------------------------------------------------------
// #179 — inversions whose interior bases are self-complementary.
// ---------------------------------------------------------------------------

/// #179's reference, reconstructed.
///
/// The issue states the bases at each locus it uses but never published the
/// contig, so this rebuilds one: a repeating filler with #179's stated bases
/// written in at #179's stated positions.
///
/// **The reconstruction is corroborated, not assumed.** Examples 1 and 2 below
/// reproduce the `Ferro:` output recorded in the issue *byte for byte*, which
/// they could not do if the filler were affecting those loci. Example 3 does
/// not — its behaviour has moved since #179 was filed, which is recorded rather
/// than hidden.
fn issue_179_reference() -> String {
    let mut seq: Vec<u8> = "CTAGGTCA".bytes().cycle().take(1200).collect();
    // (1-based position, base) — every base #179 states.
    let pins: &[(usize, u8)] = &[
        // Example 3: 1003-1005 = AGA (deleted), 1006-1008 = AGG (inverted to CCT).
        (1003, b'A'),
        (1004, b'G'),
        (1005, b'A'),
        (1006, b'A'),
        (1007, b'G'),
        (1008, b'G'),
        // Example 1: 1040-1043 = TTAC (revcomp GTAA), interior 1041-1042 `TA`
        // self-complementary; plus the flanking edits at 1009 and 1051-1052.
        (1009, b'T'),
        (1040, b'T'),
        (1041, b'T'),
        (1042, b'A'),
        (1043, b'C'),
        (1051, b'T'),
        (1052, b'A'),
        // Example 2: 1018-1021 = TTAC, interior 1019-1020 `TA`; plus 1017, 1037, 1154.
        (1017, b'C'),
        (1018, b'T'),
        (1019, b'T'),
        (1020, b'A'),
        (1021, b'C'),
        (1037, b'C'),
        (1154, b'A'),
    ];
    for (pos, base) in pins {
        seq[pos - 1] = *base;
    }
    String::from_utf8(seq).expect("ascii")
}

/// `(label, input, current output, verdict, what #179 wanted)`
const ISSUE_179_CASES: &[(&str, &str, &str, Verdict, &str)] = &[
    (
        "ex1",
        "REF:g.[1009T>G;1040T>G;1043C>A;1051T>G;1052A>G]",
        "REF:g.[1009T>G;1040T>G;1043C>A;1051_1052delinsGG]",
        Verdict::Gap,
        "REF:g.[1009T>G;1040_1043inv;1051_1052delinsGG] — 1040-1043 TTAC -> GTAA is \
         revcomp; the interior TA is self-complementary so the two substitutions are \
         separated and never merge into a delins for the inv check to see. Matches \
         the `Ferro:` line in #179 exactly.",
    ),
    (
        "ex2",
        "REF:g.[1017C>G;1018T>G;1021C>A;1037C>A;1154A>C]",
        "REF:g.[1017_1018delinsGG;1021C>A;1037C>A;1154A>C]",
        Verdict::Gap,
        "REF:g.[1017C>G;1018_1021inv;1037C>A;1154A>C] — same shape one base over, and \
         the adjacent pair merges the WRONG way: 1017+1018 coalesce, which consumes \
         the base the inversion needed. Matches the `Ferro:` line in #179 exactly.",
    ),
    (
        "ex3",
        "REF:g.[1003del;1004del;1005del;1008_1009insCCT]",
        "REF:g.1005_1008delinsGCCT",
        Verdict::Gap,
        "REF:g.[1005A>G;1006_1008inv] — a net inversion assembled from deletions plus \
         an insertion. NOTE: behaviour has MOVED since #179 was filed, which reported \
         `[1004_1005del;1006del;1008_1009insCCT]`. It is now a single spanning delins \
         — a different answer, still not the inversion.",
    ),
];

/// How many #179 rows still disagree with what the issue argued for.
///
/// All three, today. Mirrors [`ISSUE_87_OPEN_GAPS`] so that the `Verdict` column
/// of [`ISSUE_179_CASES`] is *read* rather than merely documented — a row flipped
/// to `Guard` without a matching change here is a failure, not a silent edit.
const ISSUE_179_OPEN_GAPS: usize = 3;

/// #179's examples still normalize as they do today.
///
/// All three are gaps. That is expected and is not a bug being tolerated
/// silently: `inv` recognition across a self-complementary interior is
/// **discretionary**, not mandated — the mandatory merge threshold is one
/// intervening nucleotide and these have two, and SVD-WG010, which would have
/// widened it, was rejected. So the value here is the record, not a demand:
/// these are the shapes a discretionary inv pass would have to catch, with
/// today's answers pinned so such a pass can be measured against them.
#[test]
fn issue_179_self_complementary_inversions_are_unchanged() {
    let normalizer = Normalizer::new(provider_for(&issue_179_reference(), "REF"));
    for (label, input, expected, _, wanted) in ISSUE_179_CASES {
        assert_eq!(
            &normalize_with(&normalizer, input),
            expected,
            "#179 {label}: {input}\n  #179 argued for: {wanted}"
        );
    }
}

/// The reconstruction's premise: the bases #179 states are the bases it gets.
///
/// If this fails, every row above is measuring the wrong locus.
#[test]
fn issue_179_reference_has_the_stated_bases() {
    let seq = issue_179_reference();
    assert_eq!(&seq[1002..1005], "AGA", "1003-1005");
    assert_eq!(&seq[1005..1008], "AGG", "1006-1008");
    assert_eq!(&seq[1017..1021], "TTAC", "1018-1021");
    assert_eq!(&seq[1039..1043], "TTAC", "1040-1043");
}

/// None of #179's three examples is recognized as an inversion today.
///
/// Split out from the table so that the day one *is*, the failure names the
/// achievement rather than looking like a broken string comparison.
#[test]
fn issue_179_no_example_is_recognized_as_an_inversion() {
    let normalizer = Normalizer::new(provider_for(&issue_179_reference(), "REF"));
    for (label, input, _, _, _) in ISSUE_179_CASES {
        let got = normalize_with(&normalizer, input);
        assert!(
            !got.contains("inv"),
            "#179 {label} now yields an inversion (`{got}`) — if this is a \
             deliberate discretionary inv pass, update the row and this test"
        );
    }
}

/// Census over the #179 table, so the count cannot drift away from the rows.
#[test]
fn issue_179_gap_census_holds() {
    let gaps = ISSUE_179_CASES
        .iter()
        .filter(|(_, _, _, verdict, _)| *verdict == Verdict::Gap)
        .count();
    assert_eq!(
        gaps, ISSUE_179_OPEN_GAPS,
        "the number of open #179 gaps moved; if a discretionary inv pass closed \
         one, flip its Verdict to Guard and lower ISSUE_179_OPEN_GAPS"
    );
}
