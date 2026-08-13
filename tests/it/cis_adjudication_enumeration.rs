//! Every `cis_confluence_adjudication` record, turned from a point assertion
//! into a property over an **enumerated description space**.
//!
//! Each record over there pins one or two spellings. Here the complete bounded
//! space of spellings that denote the *same variant* is generated, and the whole
//! space is normalized and bucketed by output. For a record whose ruling is
//! `decided` **and implemented**, the property is that the space collapses to one
//! string. For a record that is a scope carve-out or an open disagreement, the
//! property is the measured partition of the space, pinned as an observation with
//! the record that owns the residue named on it.
//!
//! # What this proves, and — more importantly — what it does NOT
//!
//! **It proves only that one variant's spellings converge.** It says *nothing*
//! about different variants staying apart: a normalizer that mapped every input
//! to the single string `"x"` would pass every convergence assertion in this
//! file. Convergence is one half of canonicalization and this module measures
//! only that half. **Do not read a green run here as more than "the space did not
//! fragment".**
//!
//! One guard rail is in place and it is worth being precise about what it buys.
//! [`assert_outputs_are_sequence_preserving`] applies every distinct output back
//! to the reference and requires it to denote what its inputs denoted — the
//! `FERRO_ASSERT_SEQUENCE` question, asked over exactly these outputs, and asked
//! on **every** partitioner arm rather than only the shipping one. That rules out
//! the degenerate normalizer collapsing this class onto a string that means
//! something else. It does **not** rule out two *distinct* variants being mapped
//! together, because nothing here ever normalizes a second variant. The
//! separating half genuinely lives elsewhere.
//!
//! Two further limits, stated so a green run is not over-read:
//!
//!  - Every expectation below is an **observation of ferro's output** unless its
//!    doc comment cites a ruling record. Where a record is cited, the record is
//!    the authority and the assertion is a guard against regressing away from it.
//!  - The space is **bounded**, in two named dimensions (below). A spelling
//!    outside those dimensions is not covered, and the audit in
//!    [`the_enumerated_space_can_vary_what_it_claims_to_vary`] names what this
//!    fixture provably cannot reach.
//!
//! # The two dimensions the space varies
//!
//! **(a) Spanning single-member descriptions.** Every `delins`/`del` whose span
//! covers all of the change — start at or before the first changed column, end at
//! or after the last — with the payload solved from the resulting sequence. On a
//! `c.` axis this reaches spans that end in the 3'UTR (`c.{i}_*1`), which is
//! where #1536's shape lives.
//!
//! **(b) Multi-member spellings with each member at every legal shift.** For each
//! seed partition, each member is independently re-positioned over a window and
//! its payload re-solved, keeping exactly those positions at which the member
//! *applied alone* denotes what the seed member denotes. The cartesian product of
//! those per-member shift classes is taken, non-monotone and overlapping
//! combinations are dropped, and what survives is verified against the apply
//! oracle. A **shift preserves the deleted length**, so this dimension does not
//! re-partition; re-partitioning is dimension (a)'s job for the one-member case
//! and is not enumerated for the multi-member case.
//!
//! Each surviving member is rendered in every spelling HGVS admits for it — a
//! deletion with and without its stated bases, a single-column replacement as
//! both a substitution and a `delins`, an insertion also as a `dup` where the
//! payload copies its 5' flank — so the space varies *typing* as well as
//! position.
//!
//! **Nothing generated here is trusted.** Every spelling is round-tripped through
//! `spdi::hgvs_to_spdi` and applied to the reference, and a spelling that does not
//! denote the class's own sequence is a hard failure rather than a silent drop
//! ([`Enumeration::build`]).
//!
//! # Arms
//!
//! `FERRO_PARTITION` is read once per process, so a single test process sees one
//! partitioner. These tests assert the **shipping** (`live`) arm, which is what
//! CI runs; under any other arm they enumerate and count exactly as before but
//! assert only the arm-independent facts (the space's size, and that every member
//! of it denotes one sequence). That is deliberate: the size assertions are what
//! make a structurally-collapsed enumeration fail loudly rather than pass as an
//! empty green.
//!
//! Running the other arms is one command per arm and the output is the `report`
//! line plus its buckets:
//!
//! ```text
//! FERRO_PARTITION=canonical cargo nextest run --features dev --test it \
//!   -E 'test(cis_adjudication_enumeration)' --no-capture
//! ```
//!
//! **One thing that fell out of doing so, recorded because it is a question for
//! an operator rather than a defect to fix here.** Re-measure it rather than
//! quoting this paragraph. Every class in this module collapses to a single
//! output under `canonical` and under `canonical-coalesced`, including the two
//! that are divergent under `live` — but for two of those classes the single
//! output the candidate arms pick is a *third* form, named by no decided record:
//! the adjacent-members class lands on the authored `[dup;del]` partition rather
//! than on the spanning `delins` that
//! `delins-adjacent-members-when-both-consume-reference` adjudicates correct, and
//! the spanning-versus-split class lands on a partition matching neither of the
//! two fixed points `cis_confluence_adjudication` pins. Converging a class is not
//! the same as converging it on the adjudicated form, and only the enumeration
//! makes the difference visible.

use std::collections::{BTreeMap, BTreeSet};

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::{MockProvider, ReferenceProvider};
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

// ---------------------------------------------------------------------------
// Fixtures — the adjudication module's own, verbatim
// ---------------------------------------------------------------------------

/// `cis_confluence_adjudication::CORE`, copied so this module's classes are that
/// module's classes rather than newly invented ones.
///
/// 64 bases, `CDS 1..=63`. That geometry matters and is audited in
/// [`the_enumerated_space_can_vary_what_it_claims_to_vary`]: there is **one**
/// 3'UTR base and **no** 5'UTR.
const CORE: &str = "TTTTTTTTTAATATATTTTAATATAATTAAAAAAATAATTTTTATAAATATATTATTTTAAAAA";

const TX_ACCESSION: &str = "NM_TEST.1";
const TX_CONTIG: &str = "chr_synth";
const PAD_OFFSET: usize = 256;
const CDS_START: u64 = 1;
const CDS_END: u64 = 63;

/// `cis_confluence_adjudication::GENOMIC_WINDOW` — `NC_000001.11:g.1001002_1001016`.
const GENOMIC_WINDOW: &str = "ATGAGGGGCCACTGT";
const GENOMIC_WINDOW_START: usize = 1_001_002;
const GENOMIC_CONTIG: &str = "NC_000001.11";

/// How far either side of the genomic case's own fifteen bases the enumeration
/// window reaches. Shifts are expected to stay inside the fifteen, and the flanks
/// give an equivalent that escaped them somewhere to be *found* rather than
/// clipped at the record's own edge.
///
/// **Read what the closure check does and does not buy here.**
/// [`Enumeration::build`] asserts the shift class is closed within
/// [`SHIFT_SEARCH`], which bounds the search *radius*. This flank is smaller than
/// that radius, so for a member in the 3' half of the window `member_shifts`'
/// upper bound is clamped by `bases.len()` before the radius is reached — and a
/// class that continued past the window edge would therefore not trip the
/// assertion. It does not arise on this fixture (every shift class here is one or
/// two positions wide, far inside the flank), but do not read the closure check as
/// a guarantee about the *window*. Widening the flank is not free either:
/// dimension (a) enumerates spans out to the window's ends, so it would move
/// [`SEPARATION_GENOMIC_SPELLINGS`].
const GENOMIC_FLANK: usize = 8;

fn transcript_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let tx_len = CORE.len() as u64;
    let g_start = PAD_OFFSET as u64 + 1;
    let g_end = PAD_OFFSET as u64 + tx_len;
    let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
    let transcript = Transcript::new(
        TX_ACCESSION.to_string(),
        Some("SYNTH".to_string()),
        Strand::Plus,
        CORE.to_string(),
        Some(CDS_START),
        Some(CDS_END),
        vec![exon],
        Some(TX_CONTIG.to_string()),
        Some(g_start),
        Some(g_end),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    provider.add_genomic_sequence(TX_CONTIG, format!("{pad}{CORE}{pad}"));
    provider.add_transcript(transcript);
    provider
}

fn genomic_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let mut sequence = "ACGT".repeat(GENOMIC_WINDOW_START.div_ceil(4)).into_bytes();
    sequence.truncate(GENOMIC_WINDOW_START - 1);
    sequence.extend_from_slice(GENOMIC_WINDOW.as_bytes());
    sequence.extend_from_slice(b"ACGTACGTAC");
    provider.add_genomic_sequence(
        GENOMIC_CONTIG,
        String::from_utf8(sequence).expect("ascii filler"),
    );
    provider
}

// ---------------------------------------------------------------------------
// The window an enumeration runs over
// ---------------------------------------------------------------------------

/// Which axis a [`Window`] renders coordinates on.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum AxisKind {
    /// `c.`, with a CDS end past which coordinates are spelled `*n`.
    Coding { cds_end: usize },
    /// `g.`, with the axis coordinate of `bases[0]`.
    Genomic { first: usize },
}

/// The reference slice an enumeration works over, plus how to spell positions in
/// it.
struct Window {
    /// The bases. Index 0 is the first base of the window.
    bases: String,
    /// The 0-based offset of `bases[0]` in the coordinate space
    /// `hgvs_to_spdi` reports triples in — 0 for a transcript axis, the
    /// contig-absolute offset of the window for a genomic one.
    spdi_origin: usize,
    axis: AxisKind,
    /// Everything before the coordinate, e.g. `NM_TEST.1:c.`.
    prefix: String,
}

impl Window {
    fn coding() -> Self {
        Self {
            bases: CORE.to_string(),
            spdi_origin: 0,
            axis: AxisKind::Coding {
                cds_end: CDS_END as usize,
            },
            prefix: format!("{TX_ACCESSION}:c."),
        }
    }

    fn genomic(provider: &MockProvider) -> Self {
        let first = GENOMIC_WINDOW_START - GENOMIC_FLANK;
        let last = GENOMIC_WINDOW_START + GENOMIC_WINDOW.len() - 1 + GENOMIC_FLANK;
        // `get_genomic_sequence` takes a 0-based half-open range, so the window
        // spelled `first..=last` in axis coordinates is `[first - 1, last)`.
        let bases = provider
            .get_genomic_sequence(GENOMIC_CONTIG, (first - 1) as u64, last as u64)
            .expect("the genomic window is served");
        assert!(
            bases.contains(GENOMIC_WINDOW),
            "the widened window must still contain the record's own fifteen bases"
        );
        Self {
            bases,
            spdi_origin: first - 1,
            axis: AxisKind::Genomic { first },
            prefix: format!("{GENOMIC_CONTIG}:g."),
        }
    }

    /// Spell a 1-based index into [`Self::bases`] as an axis coordinate.
    fn coordinate(&self, one_based: usize) -> String {
        match self.axis {
            AxisKind::Coding { cds_end } => {
                if one_based <= cds_end {
                    one_based.to_string()
                } else {
                    format!("*{}", one_based - cds_end)
                }
            }
            AxisKind::Genomic { first } => (first + one_based - 1).to_string(),
        }
    }
}

// ---------------------------------------------------------------------------
// Members, shifts, and rendering
// ---------------------------------------------------------------------------

/// One member of a description, as an SPDI-shaped triple on the window: a
/// 0-based start, the reference bases it consumes, and the bases it leaves.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct Member {
    start: usize,
    del: String,
    ins: String,
}

/// The sequence `member` denotes when applied to `bases` **on its own**.
///
/// `None` when it does not fit, or when its stated deleted bases disagree with
/// the reference — both of which make it not a spelling of anything here.
fn apply_alone(bases: &str, member: &Member) -> Option<String> {
    let end = member.start.checked_add(member.del.len())?;
    if end > bases.len() {
        return None;
    }
    if bases[member.start..end] != member.del {
        return None;
    }
    Some(format!(
        "{}{}{}",
        &bases[..member.start],
        member.ins,
        &bases[end..]
    ))
}

/// Apply a whole member list, which must be non-overlapping and in ascending
/// order.
fn apply_all(bases: &str, members: &[Member]) -> Option<String> {
    let mut cursor = 0usize;
    let mut out = String::with_capacity(bases.len());
    for member in members {
        if member.start < cursor {
            return None;
        }
        let end = member.start.checked_add(member.del.len())?;
        if end > bases.len() || bases[member.start..end] != member.del {
            return None;
        }
        out.push_str(&bases[cursor..member.start]);
        out.push_str(&member.ins);
        cursor = end;
    }
    out.push_str(&bases[cursor..]);
    Some(out)
}

/// How far either side of a seed member's own position the shift search looks.
///
/// Larger than any homopolymer run in either fixture. [`Enumeration::build`]
/// asserts the resulting shift class is **closed** — no equivalent sits at the
/// edge of the search — so a window too small to hold the whole class is a
/// failure rather than a quietly truncated space.
const SHIFT_SEARCH: usize = 14;

/// Every position at which `member` denotes exactly what it denotes where it is.
///
/// The deleted length is held fixed (that is what makes this a *shift* rather
/// than a re-partition) and the payload is **solved** at each candidate position
/// from the sequence the seed member denotes, rather than guessed — so an
/// insertion's rotation inside a tandem tract falls out instead of having to be
/// enumerated.
fn member_shifts(bases: &str, member: &Member) -> Vec<Member> {
    let Some(target) = apply_alone(bases, member) else {
        return Vec::new();
    };
    let deleted = member.del.len();
    let inserted = member.ins.len();
    let lo = member.start.saturating_sub(SHIFT_SEARCH);
    let hi = (member.start + SHIFT_SEARCH).min(bases.len());
    let mut found = Vec::new();
    for start in lo..=hi {
        if start + deleted > bases.len() || start + inserted > target.len() {
            continue;
        }
        let candidate = Member {
            start,
            del: bases[start..start + deleted].to_string(),
            ins: target[start..start + inserted].to_string(),
        };
        if apply_alone(bases, &candidate).as_deref() == Some(target.as_str()) {
            found.push(candidate);
        }
    }
    found
}

/// Every HGVS spelling of one member.
///
/// More than one where HGVS admits more than one: a deletion with and without its
/// stated bases, a single-column replacement as a substitution and as a `delins`,
/// an insertion also as a `dup` when its payload copies the bases immediately 5'
/// of it.
fn render_member(window: &Window, member: &Member) -> Vec<String> {
    let start = member.start + 1; // 1-based
    let end = member.start + member.del.len(); // 1-based inclusive
    let (del, ins) = (member.del.as_str(), member.ins.as_str());
    let mut spellings = Vec::new();

    match (del.is_empty(), ins.is_empty()) {
        (true, true) => {}
        (true, false) => {
            // A pure insertion sits between two bases, so it needs a base on
            // each side.
            if member.start == 0 || member.start >= window.bases.len() {
                return spellings;
            }
            let left = window.coordinate(member.start);
            let right = window.coordinate(member.start + 1);
            spellings.push(format!("{left}_{right}ins{ins}"));
            if member.start >= ins.len()
                && &window.bases[member.start - ins.len()..member.start] == ins
            {
                let from = window.coordinate(member.start - ins.len() + 1);
                let to = window.coordinate(member.start);
                spellings.push(if ins.len() == 1 {
                    format!("{from}dup")
                } else {
                    format!("{from}_{to}dup")
                });
            }
        }
        (false, true) => {
            let span = if del.len() == 1 {
                window.coordinate(start)
            } else {
                format!("{}_{}", window.coordinate(start), window.coordinate(end))
            };
            spellings.push(format!("{span}del"));
            spellings.push(format!("{span}del{del}"));
        }
        (false, false) => {
            if del.len() == 1 && ins.len() == 1 {
                let at = window.coordinate(start);
                spellings.push(format!("{at}{del}>{ins}"));
                spellings.push(format!("{at}delins{ins}"));
            } else {
                let span = if del.len() == 1 {
                    window.coordinate(start)
                } else {
                    format!("{}_{}", window.coordinate(start), window.coordinate(end))
                };
                spellings.push(format!("{span}delins{ins}"));
                spellings.push(format!("{span}del{del}ins{ins}"));
            }
        }
    }
    spellings
}

/// Assemble rendered members into a whole description.
fn assemble(window: &Window, rendered: &[String]) -> String {
    if rendered.len() == 1 {
        format!("{}{}", window.prefix, rendered[0])
    } else {
        format!("{}[{}]", window.prefix, rendered.join(";"))
    }
}

// ---------------------------------------------------------------------------
// The enumeration itself
// ---------------------------------------------------------------------------

/// A ceiling on one class's space, so a combinatorial blow-up is a failure rather
/// than a slow test. Asserted, never silently applied.
const MAX_SPELLINGS: usize = 200_000;

struct Enumeration {
    /// Every generated spelling, deduplicated and sorted.
    spellings: Vec<String>,
    /// How many spellings came from dimension (a).
    spanning: usize,
    /// How many distinct spans dimension (a) accepted — the positional size of
    /// that dimension, before the per-member spelling choices multiply it.
    spans: usize,
    /// The sequence every one of them denotes.
    denoted: String,
}

impl Enumeration {
    /// Build the space for one class.
    ///
    /// `seeds` are spellings of the class taken from the adjudication module;
    /// they supply both the target sequence and the partitions dimension (b)
    /// shifts. Every seed must denote the same sequence — that is checked here
    /// rather than assumed, because a class whose seeds disagree is not a class.
    fn build(window: &Window, provider: &MockProvider, seeds: &[&str]) -> Self {
        let denoted = denotes(window, provider, seeds[0])
            .unwrap_or_else(|| panic!("seed {} must apply", seeds[0]));
        for seed in seeds {
            assert_eq!(
                denotes(window, provider, seed).as_deref(),
                Some(denoted.as_str()),
                "{seed} does not denote what {} does — these are not one class",
                seeds[0]
            );
        }

        let mut spellings: BTreeSet<String> = BTreeSet::new();

        // -- dimension (a): every span covering the change ------------------
        let (first, last) = changed_columns(&window.bases, &denoted);
        let trailing = window.bases.len() - last;
        let mut spanning = 0usize;
        let mut spans = 0usize;
        for start in 0..=first {
            for end in last..=window.bases.len() {
                let member = Member {
                    start,
                    del: window.bases[start..end].to_string(),
                    ins: format!(
                        "{}{}{}",
                        &window.bases[start..first],
                        &denoted[first..denoted.len() - trailing],
                        &window.bases[last..end]
                    ),
                };
                if apply_alone(&window.bases, &member).as_deref() != Some(denoted.as_str()) {
                    continue;
                }
                spans += 1;
                for rendered in render_member(window, &member) {
                    if spellings.insert(assemble(window, std::slice::from_ref(&rendered))) {
                        spanning += 1;
                    }
                }
            }
        }

        // -- dimension (b): each seed partition, each member shifted --------
        for seed in seeds {
            let members = seed_members(window, provider, seed);
            if members.len() < 2 {
                continue;
            }
            let classes: Vec<Vec<Member>> = members
                .iter()
                .map(|member| {
                    let shifts = member_shifts(&window.bases, member);
                    assert!(
                        !shifts.is_empty(),
                        "{seed}: member {member:?} has an empty shift class"
                    );
                    // Closure: a shift sitting at the edge of the search window
                    // means the class is larger than the window and the space is
                    // silently truncated.
                    for shift in &shifts {
                        assert!(
                            shift.start.abs_diff(member.start) < SHIFT_SEARCH,
                            "{seed}: member {member:?} has an equivalent at the edge of the \
                             shift search — SHIFT_SEARCH is too small and this space is truncated"
                        );
                    }
                    shifts
                })
                .collect();

            let mut combinations: Vec<Vec<Member>> = vec![Vec::new()];
            for class in &classes {
                let mut next = Vec::new();
                for prefix in &combinations {
                    for candidate in class {
                        let fits = prefix.last().is_none_or(|previous: &Member| {
                            previous.start + previous.del.len() <= candidate.start
                                && !(previous.del.is_empty()
                                    && candidate.del.is_empty()
                                    && previous.start == candidate.start)
                        });
                        if !fits {
                            continue;
                        }
                        let mut extended = prefix.clone();
                        extended.push(candidate.clone());
                        next.push(extended);
                    }
                }
                combinations = next;
                assert!(
                    combinations.len() <= MAX_SPELLINGS,
                    "{seed}: the shift product exceeded {MAX_SPELLINGS}"
                );
            }

            for combination in &combinations {
                if apply_all(&window.bases, combination).as_deref() != Some(denoted.as_str()) {
                    continue;
                }
                let rendered: Vec<Vec<String>> = combination
                    .iter()
                    .map(|member| render_member(window, member))
                    .collect();
                if rendered.iter().any(Vec::is_empty) {
                    continue;
                }
                for spelling in cartesian(&rendered) {
                    spellings.insert(assemble(window, &spelling));
                }
                assert!(
                    spellings.len() <= MAX_SPELLINGS,
                    "{seed}: the rendered space exceeded {MAX_SPELLINGS}"
                );
            }
        }

        // -- nothing generated is trusted -----------------------------------
        for spelling in &spellings {
            assert_eq!(
                denotes(window, provider, spelling).as_deref(),
                Some(denoted.as_str()),
                "{spelling} was generated into this class but does not denote its sequence"
            );
        }

        Self {
            spellings: spellings.into_iter().collect(),
            spanning,
            spans,
            denoted,
        }
    }

    /// Normalize the whole space and bucket it by output string.
    fn outputs(&self, provider: &MockProvider) -> BTreeMap<String, Vec<String>> {
        let normalizer = Normalizer::with_config(
            provider.clone(),
            NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
        );
        let mut buckets: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for spelling in &self.spellings {
            let variant = parse_hgvs(spelling).unwrap_or_else(|e| panic!("{spelling}: {e}"));
            let output = match normalizer.normalize(&variant) {
                Ok(normalized) => normalized.to_string(),
                Err(error) => format!("<error: {error}>"),
            };
            buckets.entry(output).or_default().push(spelling.clone());
        }
        buckets
    }

    fn len(&self) -> usize {
        self.spellings.len()
    }

    /// How many spellings in the space contain `fragment`.
    ///
    /// Used to assert **reachability**: a dimension nobody can show reaching a
    /// shape is a dimension that is not varying it, and a zero from such a space
    /// is structural rather than informative.
    fn containing(&self, fragment: &str) -> usize {
        self.spellings
            .iter()
            .filter(|spelling| spelling.contains(fragment))
            .count()
    }
}

/// The half-open column range `[first, last)` in which `bases` and `denoted`
/// differ, after trimming the common prefix and suffix.
fn changed_columns(bases: &str, denoted: &str) -> (usize, usize) {
    let reference = bases.as_bytes();
    let result = denoted.as_bytes();
    let mut first = 0usize;
    while first < reference.len() && first < result.len() && reference[first] == result[first] {
        first += 1;
    }
    let mut tail = 0usize;
    while tail < reference.len() - first
        && tail < result.len() - first
        && reference[reference.len() - 1 - tail] == result[result.len() - 1 - tail]
    {
        tail += 1;
    }
    (first, reference.len() - tail)
}

/// The members of a spelling, as triples on the window, via `hgvs_to_spdi`.
fn seed_members(window: &Window, provider: &MockProvider, spelling: &str) -> Vec<Member> {
    let members: Vec<HgvsVariant> = match parse_hgvs(spelling).expect("parses") {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let mut triples: Vec<Member> = members
        .iter()
        .map(|member| {
            let triple = hgvs_to_spdi(member, provider).expect("applies");
            let absolute = usize::try_from(triple.position).expect("non-negative");
            Member {
                start: absolute - window.spdi_origin,
                del: triple.deletion.clone(),
                ins: triple.insertion.clone(),
            }
        })
        .collect();
    triples.sort();
    triples
}

/// What a spelling denotes on the window, through `hgvs_to_spdi` and an SPDI
/// splice — the apply oracle, independent of the normalizer by construction.
///
/// `None` when the spelling's members overlap or fall outside the window, which
/// is how a generated candidate that is not actually a spelling of this class is
/// detected rather than assumed away.
fn denotes(window: &Window, provider: &MockProvider, spelling: &str) -> Option<String> {
    let variant = parse_hgvs(spelling).ok()?;
    let members: Vec<HgvsVariant> = match variant {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let mut triples: Vec<Member> = Vec::new();
    for member in &members {
        let triple = hgvs_to_spdi(member, provider).ok()?;
        let absolute = usize::try_from(triple.position).ok()?;
        triples.push(Member {
            start: absolute.checked_sub(window.spdi_origin)?,
            del: triple.deletion.clone(),
            ins: triple.insertion.clone(),
        });
    }
    triples.sort();
    apply_all(&window.bases, &triples)
}

/// The cartesian product of per-member rendering choices.
fn cartesian(choices: &[Vec<String>]) -> Vec<Vec<String>> {
    let mut product: Vec<Vec<String>> = vec![Vec::new()];
    for options in choices {
        let mut next = Vec::with_capacity(product.len() * options.len());
        for prefix in &product {
            for option in options {
                let mut extended = prefix.clone();
                extended.push(option.clone());
                next.push(extended);
            }
        }
        product = next;
    }
    product
}

// ---------------------------------------------------------------------------
// Arm reporting
// ---------------------------------------------------------------------------

/// Whether this process is running the shipping partitioner.
///
/// `FERRO_PARTITION` is read once per process by the library, so a test cannot
/// choose its own arm. The convergence expectations below are the `live` arm's;
/// under any other arm the enumeration still runs and its size is still asserted,
/// but the per-arm output partition is printed rather than asserted.
fn shipping_arm() -> bool {
    match std::env::var("FERRO_PARTITION") {
        Err(_) => true,
        Ok(value) => value.trim().is_empty() || value.trim() == "live",
    }
}

/// Print the output partition, so a bake-off run under another arm reports a
/// number instead of nothing.
fn report(class: &str, enumeration: &Enumeration, buckets: &BTreeMap<String, Vec<String>>) {
    let arm = std::env::var("FERRO_PARTITION").unwrap_or_else(|_| "live (unset)".to_string());
    println!(
        "[{class}] arm={arm} spellings={} spanning={} spans={} outputs={}",
        enumeration.len(),
        enumeration.spanning,
        enumeration.spans,
        buckets.len()
    );
    for (output, inputs) in buckets {
        println!("    {:>5}  {output}", inputs.len());
        if inputs.len() <= 12 {
            for input in inputs {
                println!("             <- {input}");
            }
        }
    }
}

/// Every **output** must denote the class's own sequence.
///
/// This is the half of canonicalization the convergence assertions do not reach.
/// Convergence alone is satisfied by a normalizer that maps everything to one
/// string; this refuses an output whose bases differ from the input's, which is
/// the `FERRO_ASSERT_SEQUENCE` question asked over exactly the outputs this
/// module produced. It does **not** establish that two *different* variants stay
/// apart — no assertion in this file does, and the module docs say so.
fn assert_outputs_are_sequence_preserving(
    class: &str,
    window: &Window,
    provider: &MockProvider,
    buckets: &BTreeMap<String, Vec<String>>,
    denoted: &str,
) {
    for output in buckets.keys() {
        assert_eq!(
            denotes(window, provider, output).as_deref(),
            Some(denoted),
            "[{class}] the output {output} does not denote the sequence its inputs did"
        );
    }
}

/// The one bucket the class must collapse to, or a failure naming the strays.
fn assert_converges(class: &str, buckets: &BTreeMap<String, Vec<String>>, expected: &str) {
    let total: usize = buckets.values().map(Vec::len).sum();
    let converged = buckets.get(expected).map_or(0, Vec::len);
    assert_eq!(
        converged,
        total,
        "[{class}] {converged}/{total} converged on {expected}; the strays are {:?}",
        buckets
            .iter()
            .filter(|(output, _)| output.as_str() != expected)
            .map(|(output, inputs)| (output.clone(), inputs.len(), inputs.first().cloned()))
            .collect::<Vec<_>>()
    );
}

/// The measured partition, as `output -> count`, pinned as an observation.
fn assert_partition(
    class: &str,
    buckets: &BTreeMap<String, Vec<String>>,
    expected: &[(&str, usize)],
) {
    let measured: Vec<(String, usize)> = buckets
        .iter()
        .map(|(output, inputs)| (output.clone(), inputs.len()))
        .collect();
    let want: Vec<(String, usize)> = expected
        .iter()
        .map(|(output, count)| ((*output).to_string(), *count))
        .collect();
    assert_eq!(measured, want, "[{class}] the output partition moved");
}

// ---------------------------------------------------------------------------
// Class 1 — record 1 on the `c.` axis
// ---------------------------------------------------------------------------

/// **Adjudicated-correct, as a property.**
///
/// `cis_confluence_adjudication::the_separation_two_members_present_is_not_a_property_of_the_variant`
/// pins three spellings of this variant. This enumerates the whole bounded space
/// of its spellings and asserts every one of them reaches the decided form.
///
/// **Authority.** Ruling record
/// `separation-is-a-property-of-the-spelling-not-of-the-variant`, `decided`: the
/// separation is read off the partition re-derived from the resulting sequence,
/// never off the input's spelling. A rule read off the *spelling* would partition
/// this space; a rule read off the *sequence* cannot, whatever the spelling. That
/// is exactly the property, and it is why the enumeration is the right instrument
/// for this record rather than a third pinned string.
#[test]
fn the_separation_class_converges_over_its_whole_description_space() {
    let provider = transcript_provider();
    let window = Window::coding();
    let enumeration = Enumeration::build(
        &window,
        &provider,
        &[
            "NM_TEST.1:c.[9del;10del;13del]",
            "NM_TEST.1:c.[9del;11del;13del]",
            "NM_TEST.1:c.9_13delinsAT",
            "NM_TEST.1:c.[9_10del;13del]",
        ],
    );
    let buckets = enumeration.outputs(&provider);
    report("separation/c.", &enumeration, &buckets);
    assert_outputs_are_sequence_preserving(
        "separation/c.",
        &window,
        &provider,
        &buckets,
        &enumeration.denoted,
    );

    // Reachability, not decoration: dimension (a) runs the span's 3' end all the
    // way to the transcript's last base, which is the single 3'UTR base `c.*1`.
    // Those are the spellings whose CDS-boundary crossing #1536 was about, so a
    // convergence claim here is worth only as much as their presence.
    assert_eq!(
        enumeration.containing("_*1"),
        CODING_STAR_ONE_SPELLINGS,
        "the CDS-boundary-crossing spellings left the space — any convergence \
         reported without them is structural"
    );

    assert_eq!(
        enumeration.len(),
        SEPARATION_CODING_SPELLINGS,
        "the enumerated space changed size — see the structural audit before re-pinning"
    );
    if shipping_arm() {
        assert_converges("separation/c.", &buckets, "NM_TEST.1:c.[9_10del;13del]");
    }
}

/// The size of class 1's space. Pinned so a change that shrinks the enumeration —
/// the way a corpus generator's `0 moved` goes quiet — fails instead of passing
/// faster.
const SEPARATION_CODING_SPELLINGS: usize = 1_084;

/// How many of class 1's spellings end their span on `c.*1` — the one 3'UTR base
/// this transcript has, and so the whole of the CDS-boundary-crossing shape the
/// space can reach.
const CODING_STAR_ONE_SPELLINGS: usize = 18;

// ---------------------------------------------------------------------------
// Class 2 — record 1 on the `g.` axis, the coordinates the ruling is written on
// ---------------------------------------------------------------------------

/// The same record on real coordinates and on an axis with no reading frame, so
/// `general.md:35`'s codon exception cannot be what carries it.
///
/// **Authority.** The `OPERATOR RULING, 2026-08-10` paragraph of
/// `separation-is-a-property-of-the-spelling-not-of-the-variant`: form A,
/// `g.[1001009_1001010del;1001013del]`.
#[test]
fn the_separation_class_converges_over_its_whole_description_space_on_real_coordinates() {
    let provider = genomic_provider();
    let window = Window::genomic(&provider);
    let enumeration = Enumeration::build(
        &window,
        &provider,
        &[
            "NC_000001.11:g.[1001009del;1001010del;1001013del]",
            "NC_000001.11:g.[1001009del;1001011del;1001013del]",
            "NC_000001.11:g.1001009_1001013delinsCA",
            "NC_000001.11:g.[1001009_1001010del;1001013del]",
        ],
    );
    let buckets = enumeration.outputs(&provider);
    report("separation/g.", &enumeration, &buckets);
    assert_outputs_are_sequence_preserving(
        "separation/g.",
        &window,
        &provider,
        &buckets,
        &enumeration.denoted,
    );

    assert_eq!(enumeration.len(), SEPARATION_GENOMIC_SPELLINGS);
    if shipping_arm() {
        assert_converges(
            "separation/g.",
            &buckets,
            "NC_000001.11:g.[1001009_1001010del;1001013del]",
        );
    }
}

const SEPARATION_GENOMIC_SPELLINGS: usize = 452;

// ---------------------------------------------------------------------------
// Class 3 — record 2, the adjacent-members ruling
// ---------------------------------------------------------------------------

/// **Observation, not adjudication, and the residue is named.**
///
/// `cis_confluence_adjudication::two_adjacent_members_that_both_consume_reference_are_one_delins`
/// pins the spanning spelling as adjudicated-correct against
/// `DNA/substitution.md:32`, and pins the other spelling as still keeping its own
/// partition — so the class is divergent on `main` by design. This enumerates the
/// space and pins the partition, so the *shape* of the residue is visible rather
/// than two of its points.
///
/// **What the enumeration adds over the two pinned points**, and it is not a
/// detail: the split is by *provenance*, exactly. Every spelling reached through
/// dimension (a) — a single spanning member, however written and wherever its
/// span starts and ends — lands on the spanning form; every spelling reached
/// through dimension (b) — the authored partition, at any shift and in any
/// typing — lands on the authored form. Not one crosses. That is the mechanism
/// `separation-is-a-property-of-the-spelling-not-of-the-variant` names — ferro
/// preserving the partition it was handed — measured as a clean bipartition of a
/// whole space rather than inferred from two rows.
#[test]
fn the_adjacent_members_class_partitions_as_measured() {
    let provider = transcript_provider();
    let window = Window::coding();
    let enumeration = Enumeration::build(
        &window,
        &provider,
        &["NM_TEST.1:c.10_13delinsTAAT", "NM_TEST.1:c.[9dup;13del]"],
    );
    let buckets = enumeration.outputs(&provider);
    report("adjacent-members", &enumeration, &buckets);
    assert_outputs_are_sequence_preserving(
        "adjacent-members",
        &window,
        &provider,
        &buckets,
        &enumeration.denoted,
    );

    assert_eq!(enumeration.len(), ADJACENT_MEMBERS_SPELLINGS);
    if shipping_arm() {
        assert_partition(
            "adjacent-members",
            &buckets,
            &[
                ("NM_TEST.1:c.10_13delinsTAAT", 1_040),
                ("NM_TEST.1:c.[9dup;13del]", 36),
            ],
        );
    }
}

const ADJACENT_MEMBERS_SPELLINGS: usize = 1_076;

// ---------------------------------------------------------------------------
// Class 4 — record 3, the `dup`-flush-against-`del` carve-out
// ---------------------------------------------------------------------------

/// **Out of scope, not right.** `cis_confluence_adjudication::a_dup_flush_against_a_del_is_left_alone`
/// pins the limit of the adjacency ruling: merging a `dup` into a flush neighbour
/// would destroy the duplication `DNA/duplication.md:18` requires, so the ruling
/// does not reach it. The partition below is therefore an observation of where a
/// deliberately unruled shape lands, pinned so it cannot drift unnoticed.
#[test]
fn the_dup_flush_carve_out_class_partitions_as_measured() {
    let provider = transcript_provider();
    let window = Window::coding();
    let enumeration = Enumeration::build(
        &window,
        &provider,
        &[
            "NM_TEST.1:c.[9T>A;13_20dup;24_31del;34_35insCACCAAAA]",
            "NM_TEST.1:c.[9T>A;16_23dup;24_31del;34_35insCACCAAAA]",
        ],
    );
    let buckets = enumeration.outputs(&provider);
    report("dup-flush-carve-out", &enumeration, &buckets);
    assert_outputs_are_sequence_preserving(
        "dup-flush-carve-out",
        &window,
        &provider,
        &buckets,
        &enumeration.denoted,
    );

    assert_eq!(enumeration.len(), DUP_FLUSH_SPELLINGS);
    if shipping_arm() {
        assert_partition(
            "dup-flush-carve-out",
            &buckets,
            &[
                ("NM_TEST.1:c.9_30delinsAAATATATTTTAATATTTTAATAAAACACC", 630),
                ("NM_TEST.1:c.[9T>A;16_23dup;24T>A;27_30delinsCACC]", 112),
                ("NM_TEST.1:c.[9T>A;16_23dup;24_31del;34_35insCACCAAAA]", 388),
                (
                    "NM_TEST.1:c.[9T>A;25_26inv;28_30delinsAAT;34_35insCACCAAAA]",
                    28,
                ),
                ("NM_TEST.1:c.[9T>A;25_30delinsTTTAATAAAACACC]", 8),
            ],
        );
    }
}

const DUP_FLUSH_SPELLINGS: usize = 1_166;

// ---------------------------------------------------------------------------
// Class 5 — record 4, the spanning `delins` against its aligned split
// ---------------------------------------------------------------------------

/// **The open one, and the one the instrument was built on.**
///
/// `cis_confluence_adjudication::a_spanning_delins_and_its_aligned_split_are_two_fixed_points`
/// pins two fixed points and names the ruling the residue waits on
/// (`delins-merge-vs-individual-gap-two-or-more`, whose scope on an input that
/// arrives already split is what is open). Nothing here decides that; the
/// partition is pinned as an observation.
///
/// **The #1536 shape is in this space, and on this tree it is no longer a
/// straggler.** The enumeration was built because an earlier run of it, on an
/// older tree, isolated the `c.{i}_*1` spellings — a lone `delins` whose span
/// ends in the 3'UTR — as the only members of the space that failed to join the
/// rest under the sequence-first `canonical` arm. Those spellings are still
/// generated here, and their count is asserted below so their *absence* can never
/// be mistaken for their convergence. On this tree they separate out on no arm:
/// #1536 is fixed and `issue_1536_cds_boundary_delins` is the guard that keeps it
/// so. What is left is the class's own residue, which is a partition question and
/// not a boundary one.
#[test]
fn the_spanning_versus_split_class_partitions_as_measured() {
    let provider = transcript_provider();
    let window = Window::coding();
    let enumeration = Enumeration::build(
        &window,
        &provider,
        &["NM_TEST.1:c.9_17delinsA", "NM_TEST.1:c.[9_12del;14_17del]"],
    );
    let buckets = enumeration.outputs(&provider);
    report("spanning-vs-split", &enumeration, &buckets);
    assert_outputs_are_sequence_preserving(
        "spanning-vs-split",
        &window,
        &provider,
        &buckets,
        &enumeration.denoted,
    );

    // The same reachability check the separation class makes, on the class the
    // instrument was built for: `c.{i}_*1` is #1536's shape and it must be in
    // the space for anything this test says about it to mean anything.
    assert_eq!(
        enumeration.containing("_*1"),
        SPANNING_VS_SPLIT_STAR_ONE_SPELLINGS,
        "the CDS-boundary-crossing spellings left the space"
    );

    assert_eq!(enumeration.len(), SPANNING_VS_SPLIT_SPELLINGS);
    if shipping_arm() {
        assert_partition(
            "spanning-vs-split",
            &buckets,
            &[
                ("NM_TEST.1:c.9_17delinsA", 864),
                ("NM_TEST.1:c.[9_12del;15_18del]", 16),
            ],
        );
    }
}

const SPANNING_VS_SPLIT_SPELLINGS: usize = 880;

/// Class 5's CDS-boundary-crossing spellings — see
/// [`CODING_STAR_ONE_SPELLINGS`].
const SPANNING_VS_SPLIT_STAR_ONE_SPELLINGS: usize = 18;

// ---------------------------------------------------------------------------
// The structural audit
// ---------------------------------------------------------------------------

/// **Name the property the enumeration keys on, and check the fixture can vary
/// it.** A bounded space bounded in the wrong dimension is this repository's
/// recurring failure — member geometry (#1456), scale (#1460), transcript
/// geometry (#1478) — each invisible until the one before it was fixed.
///
/// This test asserts what these fixtures **can** reach, and states what they
/// provably cannot, so a reader does not take "every spelling converged" for more
/// coverage than it is.
#[test]
fn the_enumerated_space_can_vary_what_it_claims_to_vary() {
    // CAN vary: a span that ends in the 3'UTR. There is exactly one 3'UTR base,
    // `c.*1`, and dimension (a) reaches it.
    assert_eq!(CORE.len(), 64, "the coding fixture is 64 bases");
    assert_eq!(CDS_END, 63, "…of which 63 are CDS, leaving one 3'UTR base");
    let window = Window::coding();
    assert_eq!(window.coordinate(64), "*1");
    assert_eq!(window.coordinate(63), "63");

    // CANNOT vary: a span that reaches the 5'UTR. `CDS_START` is 1, so there is
    // no `c.-n` coordinate on this transcript at all and no enumeration over it
    // can produce one. Any zero this module reports for a 5'UTR-crossing shape is
    // STRUCTURAL.
    assert_eq!(CDS_START, 1, "no 5'UTR exists on this fixture");

    // CANNOT vary: an exon junction. One exon, so no spelling in any of these
    // spaces crosses one, and no intronic (`c.n+m`) coordinate is reachable.
    let provider = transcript_provider();
    let transcript = provider
        .get_transcript(TX_ACCESSION)
        .expect("the fixture transcript");
    assert_eq!(
        transcript.exons.len(),
        1,
        "one exon: junction-crossing is structurally unreachable here"
    );

    // CANNOT vary: scale. Every block is far below `MAX_SPLIT_BLOCK` (1024), so
    // the length short-circuit that produces two of `spec_conformance_axis`'s
    // residual guard violations is never reached from this module.
    assert!(
        CORE.len() < 1024,
        "the whole transcript is shorter than MAX_SPLIT_BLOCK, so the length cap is unreachable"
    );

    // CAN vary: the genomic window is real coordinates on a frameless axis, and
    // it is widened past the record's own fifteen bases so the shift-closure
    // check has somewhere to fail.
    let genomic = genomic_provider();
    let genomic_window = Window::genomic(&genomic);
    assert_eq!(
        genomic_window.bases.len(),
        GENOMIC_WINDOW.len() + 2 * GENOMIC_FLANK
    );
}
