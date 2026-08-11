//! An exhaustive, spec-derived conformance corpus, enumerated deterministically.
//!
//! # The measurement gap this closes
//!
//! Conformance is currently measured by the generated spec fixture, which
//! harvests HGVS strings out of the recommendations: 934 rows, of which **52 are
//! multi-member cis alleles** and **40 carry a spec-stated answer**. That ceiling
//! is set by the spec, not by us — "matches the form the spec publishes" needs
//! the spec to have published one. So a change to the multi-member partitioner
//! measures `0 gained, 0 lost`, which is measurement blindness rather than
//! neutrality.
//!
//! The way past it is to stop needing a published answer. Generate a variant,
//! then generate *many spellings of the same variant* by re-partitioning it and
//! shifting each partition within its ambiguous run. Every spelling in such a
//! family denotes the same sequence **by construction**, which makes four
//! properties checkable with no oracle at all:
//!
//! 1. **Validity** — the output must parse, and must not violate an absolute
//!    prohibition. Two members claiming the same territory denote no sequence.
//! 2. **Confluence** — every spelling in a family must reach one output.
//! 3. **Idempotency** — each output must be its own fixed point.
//! 4. **Sequence preservation** — each output must denote the input's bases,
//!    verified through [`hgvs_to_spdi`] rather than through the normalizer.
//!
//! # Why this lives in the library
//!
//! For the same reason the rest of `conformance/` does: so `examples/` and
//! `tests/` share one definition. It buys something extra here — the corpus is
//! **enumerated at test time** rather than read from a generated artifact, so
//! every case runs on every `cargo test` with no fixture to stage, no CI artifact
//! to plumb, and no committed data file (which the project's "generate test data
//! programmatically" rule forbids anyway).
//!
//! # Multi-member is deliberately overweighted
//!
//! Real corpora contain 592 multi-member cis alleles in 9,949,738 rows
//! (0.006%, `examples/harvest_multi_member_cis.rs`). They are simultaneously the
//! rarest shape in nature and the most intricate in the code: partitioning,
//! separation, member ordering, overlap detection, coalescing and typing
//! interact *only* there. So this corpus samples them at roughly four orders of
//! magnitude above their natural rate — see [`SpecCorpus::multi_member_rows`],
//! which the axis test pins.
//!
//! # The three recorded blindnesses, and where each is varied
//!
//! This repository has been bitten three times by a corpus that could not vary
//! the property under test, each invisible to the one before. Every generator
//! added here must state where it varies all three:
//!
//! | blindness | what could not vary | varied by |
//! |---|---|---|
//! | #1456 | member **geometry** — all families paired members on disjoint territory, so no conflicting allele existed | [`Geometry`], including four geometries that denote *no* sequence and must be refused |
//! | #1460 | **scale** — 20-mer cores, so nothing reached `MAX_SPLIT_BLOCK` (1024) | [`BLOCK_LADDER`] and [`SEPARATION_LADDER`], which straddle every reachable threshold |
//! | #1478 | **transcript geometry** — `Exon::with_genomic(1, 1, tx_len, …)` builds one exon with `CDS_START = 1`, so no row crossed a junction or sat in the 5'UTR | [`RefShape::CodingMultiExon`] on both strands, with [`Region`] placing rows in the 5'UTR, across each junction, and in the 3'UTR |
//!
//! [`RefShape::CodingSingleExon`] is retained on purpose as the **control**: it
//! is exactly the shape #1478 names, so a divergence that appears only on the
//! multi-exon shapes is attributable.

use std::collections::BTreeMap;
use std::fmt;

use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use crate::reference::MockProvider;
use crate::spdi::hgvs_to_spdi;
use crate::{parse_hgvs, HgvsVariant};

// ---------------------------------------------------------------------------
// Synthetic reference conventions
// ---------------------------------------------------------------------------

/// Bases of period-4 `ACGT` padding either side of a genomic core.
///
/// Shared with `tests/it/common/synthetic.rs` and
/// `examples/generate_cis_confluence_corpus.rs`: the base immediately 5' of the
/// core is `T` and the one immediately 3' of it is `A`, so a core that neither
/// starts with `A` nor ends with `T` cannot extend the pad's own rotation and a
/// repeat tract is bounded to the core.
pub const PAD_OFFSET: usize = 256;

/// The genomic contig a `g.` row is drawn against.
pub const GENOMIC_CONTIG: &str = "NC_TEST.1";
/// The coding transcript a `c.` row is drawn against.
pub const CODING_ACCESSION: &str = "NM_TEST.1";
/// The non-coding transcript an `n.` row is drawn against.
pub const NONCODING_ACCESSION: &str = "NR_TEST.1";
/// The contig every synthetic transcript is placed on.
///
/// Spelled as an accession rather than `chr_synth` so the genomic-wrapper form
/// `NC_SYNTH.1(NM_TEST.1):c.20+2del` is expressible. That matters because
/// `checklist.md:20`'s note makes the wrapper **mandatory** for an intronic
/// position: "NM reference sequences … can only be used to describe variants in
/// introns using a `c.` prefix when a genomic reference sequence is given on
/// which the coding DNA reference sequence is annotated".
pub const TRANSCRIPT_CONTIG: &str = "NC_SYNTH.1";

/// Genomic bases between consecutive exons of a synthetic multi-exon
/// transcript.
///
/// Wide enough that an intronic offset ladder up to `±5` stays inside the
/// intron, and wide enough that a 3'-shift inside one exon cannot reach the
/// next.
const INTRON_LEN: usize = 60;

/// How far either way [`equivalent_placements`] searches for an equivalent
/// placement of one member.
///
/// A brute-force search rather than a re-derivation of the shift rules: the
/// question asked is exactly "does this denote the same sequence", and answering
/// it with the normalizer's own rules would make the corpus agree with the
/// normalizer by construction.
const SHIFT_SEARCH: isize = 24;

// ---------------------------------------------------------------------------
// The declared shape space
// ---------------------------------------------------------------------------

/// Block lengths the scale stratum enumerates, straddling every threshold the
/// canonicalizer keys on that is reachable inside a bounded corpus.
///
/// `CANONICAL_PAD` (128), `MAX_TIE_BREAK_SWEEP` (256), `MAX_SPLIT_BLOCK` (1024)
/// and `MAX_CANONICAL_WINDOW` (4096) each get a value below, on and above them,
/// because every one of those is an inclusive-or-exclusive boundary whose side
/// matters. `MAX_SHIFT_TRACT` (32768) and `MAX_APPLY_WINDOW` (100000) are
/// **declared out of bounds** rather than silently absent: see
/// [`CorpusBounds::extended_scale`].
pub const BLOCK_LADDER: &[usize] = &[
    1, 2, 4, 8, 120, 128, 136, 248, 256, 264, 1016, 1024, 1032, 4088, 4096, 4104,
];

/// Separations the scale stratum enumerates. Separation drives the *window*
/// rather than the block: the canonicalizer fetches the members' hull widened by
/// `CANONICAL_PAD` either side and refuses past `MAX_CANONICAL_WINDOW`, so a
/// large separation exercises the refusal-and-fall-back path that a large block
/// does not.
pub const SEPARATION_LADDER: &[usize] = &[
    0, 1, 2, 3, 5, 8, 120, 128, 136, 1016, 1024, 1032, 3832, 3968, 4104,
];

/// Block lengths added by [`CorpusBounds::extended_scale`], crossing
/// `MAX_SHIFT_TRACT` (32768) and approaching `MAX_APPLY_WINDOW` (100000).
///
/// Off by default and stated as a bound rather than omitted: one cell here costs
/// a ~200 kB synthetic reference and, measured, two orders of magnitude more
/// normalization time than the whole rest of the stratum, so including it by
/// default would make the corpus's runtime a property of three cells.
pub const EXTENDED_BLOCK_LADDER: &[usize] = &[32_760, 32_768, 32_776, 65_536];

/// Separations between consecutive members in the dense multi-member strata.
///
/// `0` is flush adjacency — members sharing a block with nothing between them —
/// and `8` is far enough that no partitioner should merge them. `1` and `2` are
/// the two values `general.md:34` and `general.md:35` disagree about.
pub const DENSE_SEPARATIONS: &[usize] = &[0, 1, 2, 3, 5, 8];

/// Member footprint / payload sizes in the dense strata.
pub const DENSE_PAYLOADS: &[usize] = &[1, 2, 4];

/// Member counts the dense multi-member strata enumerate.
///
/// Two is the smallest allele that reaches the partitioner at all
/// (`canonicalize_from_sequence` is gated on `members.len() > 1`); four is where
/// the bulk corpora run out of evidence entirely — the real-data harvest found
/// seven three-member and three four-member rows in 9.9M.
pub const MEMBER_COUNTS: &[usize] = &[2, 3, 4];

/// The bounds of one corpus. Every pinned number is measured over these, so
/// changing one re-rolls the census.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorpusBounds {
    /// Reference-sequence seeds. Two cores per seed (an `AT` and an `ACGT`
    /// alphabet). Prefix-stable: a smaller count is a strict prefix of a larger
    /// one, so a reduced run enumerates a strict subset of the full run's cases.
    pub seeds: u32,
    /// Include [`EXTENDED_BLOCK_LADDER`].
    pub extended_scale: bool,
}

impl Default for CorpusBounds {
    fn default() -> Self {
        Self {
            seeds: 1,
            extended_scale: false,
        }
    }
}

// ---------------------------------------------------------------------------
// Row taxonomy
// ---------------------------------------------------------------------------

/// Which properties a row can be asked about.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum RowKind {
    /// Two or more spellings verified to denote one sequence. All four
    /// properties apply.
    Family,
    /// One spelling whose denoted sequence is not recoverable by the corpus's
    /// oracle — an intronic offset, which SPDI cannot express. **Validity and
    /// idempotency only**; confluence and sequence preservation are reported
    /// `VACUOUS` for these rather than silently counted as passes.
    Single,
    /// A description that denotes no single sequence. The property is that the
    /// implementation **refuses** it.
    Conflict,
    /// A description the recommendations prohibit outright — a genomic offset, a
    /// hyphen range, an `X` base, an intronic position on a bare transcript
    /// accession. The property is again that the implementation **refuses** it,
    /// and the row carries the clause and its [`Strength`].
    Prohibited,
}

/// How strongly the recommendations state a prohibition.
///
/// The split exists because "is not allowed" and "can only be used … when" are
/// not the same claim, and this repository's own `CLAUDE.md` records that
/// **uppercase RFC 2119 keywords appear exactly once outside `style.md`** — so
/// keyword strength cannot rank clauses and the wording has to be quoted instead.
/// The axis test pins the two counts separately and asserts on neither: a
/// [`Strength::Conditional`] acceptance is a finding to adjudicate, not a
/// regression to fail.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Strength {
    /// The spec says "is not allowed", "is not correct", "is invalid", or
    /// "MUST NOT" in as many words.
    Absolute,
    /// The spec states a condition or a preference from which a prohibition
    /// follows, without using prohibitive words.
    Conditional,
}

impl Strength {
    /// Stable label, used in censuses.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Absolute => "absolute",
            Self::Conditional => "conditional",
        }
    }
}

impl fmt::Display for RowKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = match self {
            Self::Family => "family",
            Self::Single => "single",
            Self::Conflict => "conflict",
            Self::Prohibited => "prohibited",
        };
        f.write_str(text)
    }
}

/// How a row's members are laid out relative to one another (#1456).
///
/// The first three denote a sequence and produce [`RowKind::Family`] rows; the
/// last four denote none and produce [`RowKind::Conflict`] rows. Before #1456
/// every generated family was [`Geometry::Disjoint`], so a conflicting allele
/// could not exist in the corpus and `0 of 18,432` was reported three times as
/// though it were evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Geometry {
    /// One or more unchanged reference bases between consecutive footprints.
    Disjoint,
    /// Zero unchanged bases between consecutive footprints.
    FlushAdjacent,
    /// A pure insertion at the interbase immediately 5' of another member's
    /// footprint. Legal — the footprints do not intersect — but the two members
    /// share an endpoint, which is where an off-by-one in overlap detection
    /// lives.
    CoincidentEndpoint,
    /// One member's footprint lies strictly inside another's.
    Nested,
    /// Two footprints partially intersect.
    Overlapping,
    /// Two pure insertions at one interbase, whose joint denotation is
    /// undefined because they have no order.
    CoincidentInsertions,
    /// A deletion of a span and a duplication drawn from inside it —
    /// `general.md:58`'s "descriptions removing part of a reference sequence and
    /// replacing it with part of the same sequence are not allowed".
    SelfReplacement,
}

impl Geometry {
    /// Whether a design with this geometry denotes a single sequence.
    #[must_use]
    pub fn denotes_a_sequence(self) -> bool {
        matches!(
            self,
            Self::Disjoint | Self::FlushAdjacent | Self::CoincidentEndpoint
        )
    }

    /// Stable label, used in row ids and censuses.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Disjoint => "disjoint",
            Self::FlushAdjacent => "flush",
            Self::CoincidentEndpoint => "coincident-endpoint",
            Self::Nested => "nested",
            Self::Overlapping => "overlapping",
            Self::CoincidentInsertions => "coincident-insertions",
            Self::SelfReplacement => "self-replacement",
        }
    }
}

/// The syntactic mechanism a row's members are combined by.
///
/// **Enumerated by mechanism, never by scanning for `[`.** A bracket scan is
/// wrong in both directions, measured over the DNA slice of the recommendations:
/// of 14 bracket-bearing entries only 9 use the allele-membership mechanism —
/// repeat-count `[n]` and composite insertion payload `ins[a;b]` are
/// **single-member** shapes that merely look multi-member — while five genuinely
/// multi-member rules carry no brackets at all, the unknown-phase `(;)` operator
/// (`DNA/alleles.md:20`) among them. Keying on `[` would therefore inflate the
/// multi-member share with single-member shapes *and* leave `(;)` ungenerated
/// while the report showed the axis covered. That is the same blindness class as
/// #1456/#1460/#1478.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Mechanism {
    /// One variant, no combining operator.
    Lone,
    /// `[a;b]` — two or more variants in cis on one allele
    /// (`DNA/alleles.md:16`).
    Cis,
    /// `[a](;)[b]` / `a(;)b` — phase unknown (`DNA/alleles.md:20`).
    UnknownPhase,
    /// `[a];[b]` — in trans, on different chromosomes
    /// (`DNA/alleles.md:17`).
    Trans,
    /// `ins[a;b]` — a composite insertion payload. **Single-member.**
    CompositePayload,
    /// `<span><unit>[n]` — a repeat count. **Single-member.**
    RepeatCount,
}

impl Mechanism {
    /// Whether the mechanism combines two or more *allele members*.
    ///
    /// `false` for [`Self::CompositePayload`] and [`Self::RepeatCount`], which is
    /// the whole point of the type.
    #[must_use]
    pub fn combines_members(self) -> bool {
        matches!(self, Self::Cis | Self::UnknownPhase | Self::Trans)
    }

    /// Stable label, used in row ids and censuses.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Lone => "lone",
            Self::Cis => "cis",
            Self::UnknownPhase => "unknown-phase",
            Self::Trans => "trans",
            Self::CompositePayload => "composite-payload",
            Self::RepeatCount => "repeat-count",
        }
    }
}

/// Where in a transcript's structure a row sits (#1478).
///
/// Only [`Region::Anywhere`] is meaningful on a genomic axis. The rest are the
/// placements a single-exon `CDS_START = 1` transcript makes structurally
/// impossible.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Region {
    /// Genomic, or a transcript row placed without regard to structure.
    Anywhere,
    /// Wholly inside the 5'UTR — `c.-n` positions.
    Utr5,
    /// Straddling `c.-1` / `c.1`.
    CdsStart,
    /// Interior of the CDS, inside one exon.
    MidCds,
    /// Straddling the exon 1 / exon 2 junction.
    ExonJunction1,
    /// Straddling the exon 2 / exon 3 junction.
    ExonJunction2,
    /// Straddling the last CDS base and the first 3'UTR base.
    CdsEnd,
    /// Wholly inside the 3'UTR — `c.*n` positions.
    Utr3,
    /// An intronic offset position — `c.n+m` / `c.n-m`.
    Intronic,
}

impl Region {
    /// Stable label, used in row ids and censuses.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Anywhere => "anywhere",
            Self::Utr5 => "utr5",
            Self::CdsStart => "cds-start",
            Self::MidCds => "mid-cds",
            Self::ExonJunction1 => "junction-1",
            Self::ExonJunction2 => "junction-2",
            Self::CdsEnd => "cds-end",
            Self::Utr3 => "utr3",
            Self::Intronic => "intronic",
        }
    }
}

/// Which synthetic reference a row is drawn against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefShape {
    /// A padded synthetic contig addressed with `g.`.
    Genomic,
    /// One exon, `CDS_START == 1` — the exact shape #1478 names, retained as the
    /// control so a multi-exon-only divergence is attributable.
    CodingSingleExon,
    /// Three exons with a real 5'UTR, CDS and 3'UTR, on the given strand.
    CodingMultiExon(Strand),
    /// Three exons and no CDS, addressed with `n.`, on the given strand.
    NonCodingMultiExon(Strand),
}

impl RefShape {
    /// Every shape the corpus draws against, in a deterministic order.
    #[must_use]
    pub fn all() -> Vec<Self> {
        vec![
            Self::Genomic,
            Self::CodingSingleExon,
            Self::CodingMultiExon(Strand::Plus),
            Self::CodingMultiExon(Strand::Minus),
            Self::NonCodingMultiExon(Strand::Plus),
            Self::NonCodingMultiExon(Strand::Minus),
        ]
    }

    /// The shapes that carry transcript structure, so a [`Region`] means
    /// something.
    #[must_use]
    pub fn structured() -> Vec<Self> {
        Self::all()
            .into_iter()
            .filter(|shape| !matches!(shape, Self::Genomic))
            .collect()
    }

    /// Stable label, used in row ids and censuses.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Genomic => "g",
            Self::CodingSingleExon => "c1",
            Self::CodingMultiExon(Strand::Minus) => "c3m",
            Self::CodingMultiExon(_) => "c3p",
            Self::NonCodingMultiExon(Strand::Minus) => "n3m",
            Self::NonCodingMultiExon(_) => "n3p",
        }
    }

    /// The HGVS coordinate prefix.
    #[must_use]
    pub fn prefix(self) -> &'static str {
        match self {
            Self::Genomic => "g",
            Self::CodingSingleExon | Self::CodingMultiExon(_) => "c",
            Self::NonCodingMultiExon(_) => "n",
        }
    }

    /// Whether the shape has more than one exon, so a junction exists.
    #[must_use]
    pub fn is_multi_exon(self) -> bool {
        matches!(self, Self::CodingMultiExon(_) | Self::NonCodingMultiExon(_))
    }
}

/// One corpus row: a synthetic reference, what it denotes, and the spellings.
#[derive(Debug, Clone)]
pub struct Row {
    /// Deterministic identifier, derived from the design parameters. Unique
    /// within a corpus and independent of iteration order, so a failing row can
    /// be named in a committed regression test.
    pub id: String,
    /// Which properties apply.
    pub kind: RowKind,
    /// Which stratum enumerated it.
    pub stratum: &'static str,
    /// The synthetic reference.
    pub shape: RefShape,
    /// The variable part of the reference sequence.
    pub core: String,
    /// The sequence after the design is applied, on the axis's own reference —
    /// the row's ground truth. `None` for [`RowKind::Single`] and
    /// [`RowKind::Conflict`].
    pub denoted: Option<String>,
    /// Distinct spellings. At least two for a [`RowKind::Family`], exactly one
    /// otherwise.
    pub spellings: Vec<String>,
    /// Members in the authored design.
    pub members: usize,
    /// Unchanged reference bases between consecutive members.
    pub separation: usize,
    /// Reference bases the union of the members spans.
    pub block_len: usize,
    /// How the members are combined syntactically. Never inferred from
    /// brackets; see [`Mechanism`].
    pub mechanism: Mechanism,
    /// The clause a [`RowKind::Prohibited`] row cites, and how strongly it is
    /// stated.
    pub prohibition: Option<(&'static str, Strength)>,
    /// Guards this row would catch a violation of: behaviour that would
    /// implement a **rejected** consultation proposal.
    ///
    /// `general.md:36-39` announces a "new recommendation" that "two variants
    /// separated by less than two nucleotides should be described as a
    /// 'delins'". That is SVD-WG010, opened 2021-05-15 and **closed rejected**
    /// 2021-07-31 (`consultation/SVD-WG010.md:5-8`) — the note in the
    /// recommendations is stale and announces a change that did not happen. A
    /// separation floor of two on a **frameless** axis is that rejected proposal,
    /// and ferro has shipped it before: `coalesce_coding_frame_separation`
    /// emitted SVD-WG010's own worked example. So the shape is generated with the
    /// guard attached rather than left to be re-discovered.
    pub negative_guards: Vec<&'static str>,
    /// Member layout (#1456).
    pub geometry: Geometry,
    /// Transcript placement (#1478).
    pub region: Region,
    /// Thresholds this row's block or window straddles (#1460).
    pub scale_bands: Vec<&'static str>,
    /// Rule tags this row exercises. Joined against the committed rule
    /// inventory, which is what turns "we generate a lot" into "we generate
    /// *this clause*".
    pub rules: Vec<&'static str>,
}

impl Row {
    /// The design as authored — the first spelling, by construction of
    /// [`candidate_spellings`], which emits the authored form before any
    /// respelling.
    ///
    /// Load-bearing for the negative guard: the *spanning `delins`* respelling of
    /// a two-member design is itself a single member, so a guard evaluated over
    /// every spelling would report the corpus's own candidate as a merge and
    /// count 575 violations where the question was only ever about what the
    /// normalizer does to the authored pair.
    #[must_use]
    pub fn authored_spelling(&self) -> &str {
        self.spellings
            .first()
            .map_or("", std::string::String::as_str)
    }

    /// Whether the row has more than one authored **allele member**.
    ///
    /// Both halves are required: a composite insertion payload or a repeat count
    /// has one member however many `;`-separated fragments it holds, and a
    /// bracket-free `(;)` row has two however few brackets it holds.
    #[must_use]
    pub fn is_multi_member(&self) -> bool {
        self.members > 1 && self.mechanism.combines_members()
    }

    /// The provider the row's coordinates resolve against, and the sequence they
    /// address.
    ///
    /// Rebuilt from [`Self::core`] and [`Self::shape`] rather than stored, so a
    /// corpus of 14,000 rows does not carry 14,000 copies of a padded contig.
    /// Rows sharing a `(shape, core)` share a provider, which is what makes the
    /// axis test's cost linear in rows rather than in bases.
    #[must_use]
    pub fn frame(&self) -> Frame {
        Frame::build(self.shape, &self.core)
    }

    /// The key rows sharing one provider agree on.
    #[must_use]
    pub fn frame_key(&self) -> (RefShape, &str) {
        (self.shape, self.core.as_str())
    }
}

// ---------------------------------------------------------------------------
// Frames: a synthetic reference, plus the coordinate arithmetic for it
// ---------------------------------------------------------------------------

/// A materialized synthetic reference: the provider, the sequence HGVS
/// coordinates address, and the mapping from an index into that sequence to an
/// HGVS position label.
///
/// The served sequence is the padded contig for a `g.` row and the *transcript*
/// for a `c.`/`n.` row, because that is the frame `hgvs_to_spdi` reports
/// positions in for each axis (0-based, verified by
/// [`tests::spdi_positions_are_zero_based_on_the_served_sequence`]).
#[derive(Clone)]
pub struct Frame {
    shape: RefShape,
    accession: &'static str,
    served: String,
    /// 0-based offset of the core within [`Self::served`].
    core_offset: usize,
    /// 1-based inclusive transcript CDS bounds, for a coding shape.
    cds: Option<(usize, usize)>,
    /// 1-based inclusive transcript bounds of each exon, in transcript order.
    exons: Vec<(usize, usize)>,
    provider: MockProvider,
}

impl Frame {
    /// Build the reference for `shape` over `core`.
    ///
    /// # Panics
    ///
    /// If `core` is shorter than 24 bases, which no stratum generates: the
    /// multi-exon layouts need three non-empty exons plus a 5'UTR and a 3'UTR.
    #[must_use]
    pub fn build(shape: RefShape, core: &str) -> Self {
        assert!(
            core.len() >= 24,
            "a synthetic core must be at least 24 bases, got {}",
            core.len()
        );
        match shape {
            RefShape::Genomic => {
                let served = padded(core);
                let mut provider = MockProvider::new();
                provider.add_genomic_sequence(GENOMIC_CONTIG, served.clone());
                Self {
                    shape,
                    accession: GENOMIC_CONTIG,
                    served,
                    core_offset: PAD_OFFSET,
                    cds: None,
                    exons: Vec::new(),
                    provider,
                }
            }
            RefShape::CodingSingleExon => {
                // The #1478 control: one exon, `CDS_START == 1`, so `c.p` is the
                // core's 1-based position `p` and no row can cross a junction or
                // reach the 5'UTR.
                let cds = (1usize, core.len() - 1);
                let exons = vec![(1usize, core.len())];
                let provider = transcript_provider(
                    CODING_ACCESSION,
                    Strand::Plus,
                    core,
                    Some(cds),
                    &[(1, core.len())],
                );
                Self {
                    shape,
                    accession: CODING_ACCESSION,
                    served: core.to_string(),
                    core_offset: 0,
                    cds: Some(cds),
                    exons,
                    provider,
                }
            }
            RefShape::CodingMultiExon(strand) | RefShape::NonCodingMultiExon(strand) => {
                let exons = three_exon_layout(core.len());
                let coding = matches!(shape, RefShape::CodingMultiExon(_));
                let cds = coding.then(|| cds_layout(core.len()));
                let accession = if coding {
                    CODING_ACCESSION
                } else {
                    NONCODING_ACCESSION
                };
                let provider = transcript_provider(accession, strand, core, cds, &exons);
                Self {
                    shape,
                    accession,
                    served: core.to_string(),
                    core_offset: 0,
                    cds,
                    exons,
                    provider,
                }
            }
        }
    }

    /// The provider.
    #[must_use]
    pub fn provider(&self) -> &MockProvider {
        &self.provider
    }

    /// The sequence the row's HGVS coordinates address.
    #[must_use]
    pub fn served(&self) -> &str {
        &self.served
    }

    /// 0-based offset of the core within [`Self::served`].
    #[must_use]
    pub fn core_offset(&self) -> usize {
        self.core_offset
    }

    /// The accession spellings are written against.
    #[must_use]
    pub fn accession(&self) -> &'static str {
        self.accession
    }

    /// The accession an **intronic** position must be written against:
    /// `NC_SYNTH.1(NM_TEST.1)`.
    ///
    /// `checklist.md:20`'s note makes the genomic wrapper mandatory there — "NM
    /// reference sequences cover mature transcripts and **do not contain** intron
    /// and gene flanking sequences, and can only be used to describe variants in
    /// introns using a `c.` prefix when a genomic reference sequence is given on
    /// which the coding DNA reference sequence is annotated". `None` for a
    /// genomic frame, which has no introns and admits no offsets at all.
    #[must_use]
    pub fn wrapped_accession(&self) -> Option<String> {
        self.shape
            .is_multi_exon()
            .then(|| format!("{TRANSCRIPT_CONTIG}({})", self.accession))
    }

    /// The HGVS position label for a 0-based index into [`Self::served`].
    ///
    /// For a coding shape this is where `-n` / `n` / `*n` is decided, which is
    /// the whole point of [`Region::Utr5`] and [`Region::Utr3`] being reachable
    /// at all.
    #[must_use]
    pub fn label(&self, index: usize) -> String {
        let one_based = index + 1;
        match self.cds {
            None => one_based.to_string(),
            Some((cds_start, cds_end)) => {
                if one_based < cds_start {
                    format!("-{}", cds_start - one_based)
                } else if one_based <= cds_end {
                    (one_based - cds_start + 1).to_string()
                } else {
                    format!("*{}", one_based - cds_end)
                }
            }
        }
    }

    /// An intronic position label: the offset `delta` from the exonic base at
    /// `index`, spelled `+`/`-` per HGVS.
    ///
    /// Returns `None` when the shape has no introns, or when `index` is not an
    /// exon boundary in the direction asked for — an intronic offset is only
    /// meaningful hanging off the last base of an exon (`+`) or the first base of
    /// the next (`-`).
    #[must_use]
    pub fn intronic_label(&self, index: usize, delta: isize) -> Option<String> {
        if !self.shape.is_multi_exon() || delta == 0 {
            return None;
        }
        let one_based = index + 1;
        let boundary_ok = self.exons.iter().any(|&(start, end)| {
            (delta > 0 && one_based == end && end != self.served.len())
                || (delta < 0 && one_based == start && start != 1)
        });
        if !boundary_ok {
            return None;
        }
        let base = self.label(index);
        Some(if delta > 0 {
            format!("{base}+{delta}")
        } else {
            format!("{base}{delta}")
        })
    }

    /// A 0-based served index inside `region`, or `None` when the region does not
    /// exist on this shape.
    ///
    /// `width` is how many bases the design needs from that index onward; a
    /// region too narrow for the design yields `None` rather than a silently
    /// clamped placement.
    #[must_use]
    pub fn region_start(&self, region: Region, width: usize) -> Option<usize> {
        let len = self.served.len();
        let (cds_start, cds_end) = self.cds.unwrap_or((1, len));
        let fits = |start: usize| (start + width <= len).then_some(start);
        match region {
            // Eight bases into the **core**, so a 5'-shifting member has
            // somewhere to travel that is still inside the core.
            //
            // `core_offset` is not decoration: a genomic frame's served sequence
            // begins with 256 bases of period-4 `ACGT` padding, so an offset
            // measured from the served sequence would place every genomic row
            // inside a perfect tandem repeat rather than in the drawn core. That
            // was the first revision's behaviour and it silently turned the whole
            // genomic half of the corpus into a repeat-tract measurement.
            Region::Anywhere | Region::MidCds => {
                let start = if self.shape.is_multi_exon() {
                    let mid = self.exons.get(1).map_or((1, len), |&e| e);
                    mid.0 + 3 - 1
                } else {
                    self.core_offset + 8
                };
                fits(start)
            }
            Region::Utr5 => {
                if cds_start < 4 {
                    return None;
                }
                fits(1).filter(|&start| start + width < cds_start - 1)
            }
            Region::CdsStart => {
                // Straddle `c.-1` / `c.1`: the design's first base is the last
                // 5'UTR base.
                (cds_start >= 2).then_some(())?;
                fits(cds_start.saturating_sub(2))
            }
            Region::ExonJunction1 | Region::ExonJunction2 => {
                let which = if region == Region::ExonJunction1 {
                    0
                } else {
                    1
                };
                let exon = *self.exons.get(which)?;
                // Straddle the junction: start one base before the exon's last
                // base, so a design of width ≥ 2 crosses into the next exon.
                (width >= 2).then_some(())?;
                (exon.1 >= 2).then_some(())?;
                fits(exon.1 - 1)
            }
            Region::CdsEnd => {
                self.cds?;
                (cds_end >= 2).then_some(())?;
                fits(cds_end - 1)
            }
            Region::Utr3 => {
                self.cds?;
                fits(cds_end + 1).filter(|&start| start + width <= len)
            }
            // An intronic row is not placed by served index; see
            // [`Self::intronic_label`].
            Region::Intronic => None,
        }
    }

    /// The exon boundaries, 1-based inclusive transcript coordinates.
    #[must_use]
    pub fn exons(&self) -> &[(usize, usize)] {
        &self.exons
    }

    /// The 1-based inclusive CDS bounds, for a coding shape.
    #[must_use]
    pub fn cds(&self) -> Option<(usize, usize)> {
        self.cds
    }
}

/// Wrap `core` in [`PAD_OFFSET`] bases of period-4 `ACGT` on each side.
#[must_use]
pub fn padded(core: &str) -> String {
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    format!("{pad}{core}{pad}")
}

/// Three exons over a transcript of `tx_len` bases, 1-based inclusive.
///
/// Proportional rather than fixed so the scale stratum's long cores still get a
/// junction in the middle of the design's reach.
///
/// `pub(crate)` so [`crate::conformance::synthetic_protein`] draws its exon
/// layout from the same function rather than a second copy of the arithmetic.
pub(crate) fn three_exon_layout(tx_len: usize) -> Vec<(usize, usize)> {
    let first = (tx_len / 3).max(8);
    let second = (tx_len / 3).max(8);
    vec![
        (1, first),
        (first + 1, first + second),
        (first + second + 1, tx_len),
    ]
}

/// CDS bounds leaving a real 5'UTR and 3'UTR, 1-based inclusive.
///
/// The 5'UTR is at least 4 bases so `c.-1`..`c.-4` exist, and the 3'UTR at least
/// 4 so `c.*1`..`c.*4` do.
fn cds_layout(tx_len: usize) -> (usize, usize) {
    let utr5 = (tx_len / 8).clamp(4, tx_len / 3);
    let utr3 = (tx_len / 8).clamp(4, tx_len / 3);
    (utr5 + 1, tx_len - utr3)
}

/// Reverse complement, uppercase DNA.
fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// Build a provider holding one synthetic transcript and the contig it sits on.
///
/// The contig interleaves the exon blocks with [`INTRON_LEN`] intronic bases. On
/// the minus strand the exon blocks are reverse-complemented and laid out in
/// reverse exon order, so exon 1 occupies the highest genomic coordinates — the
/// arrangement `tests/it/issue_214_repeat_unit_divides.rs` and
/// `tests/it/coverage_gap_tests.rs` build by hand.
///
/// `pub(crate)` so [`crate::conformance::synthetic_protein`] builds its
/// transcript through this one function: a protein frame whose contig, exon
/// records and strand handling came from a second implementation would not be
/// the same molecule the corpus's `c.` rows are drawn against.
pub(crate) fn transcript_provider(
    accession: &'static str,
    strand: Strand,
    tx: &str,
    cds: Option<(usize, usize)>,
    exons: &[(usize, usize)],
) -> MockProvider {
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    let intron = "GATTACA".repeat(INTRON_LEN / 7 + 1);
    let intron = &intron[..INTRON_LEN];

    // Exon blocks in genomic order.
    let mut blocks: Vec<String> = exons
        .iter()
        .map(|&(start, end)| tx[start - 1..end].to_string())
        .collect();
    if strand == Strand::Minus {
        blocks = blocks.iter().rev().map(|b| reverse_complement(b)).collect();
    }

    let mut contig = String::with_capacity(2 * PAD_OFFSET + tx.len() + 2 * INTRON_LEN);
    contig.push_str(&pad);
    // 0-based genomic offsets of each block, in genomic order.
    let mut block_starts = Vec::with_capacity(blocks.len());
    for (index, block) in blocks.iter().enumerate() {
        if index > 0 {
            contig.push_str(intron);
        }
        block_starts.push(contig.len());
        contig.push_str(block);
    }
    contig.push_str(&pad);

    // Back to transcript order, and 1-based.
    let mut exon_records = Vec::with_capacity(exons.len());
    for (index, &(tx_start, tx_end)) in exons.iter().enumerate() {
        let genomic_index = if strand == Strand::Minus {
            exons.len() - 1 - index
        } else {
            index
        };
        let g_start = block_starts[genomic_index] as u64 + 1;
        let g_end = g_start + (tx_end - tx_start) as u64;
        let number = u32::try_from(index + 1).unwrap_or(u32::MAX);
        exon_records.push(Exon::with_genomic(
            number,
            tx_start as u64,
            tx_end as u64,
            g_start,
            g_end,
        ));
    }

    let span_start = *block_starts.first().expect("at least one exon") as u64 + 1;
    let span_end = contig.len() as u64 - PAD_OFFSET as u64;
    let transcript = Transcript::new(
        accession.to_string(),
        Some("SYNTH".to_string()),
        strand,
        tx.to_string(),
        cds.map(|(start, _)| start as u64),
        cds.map(|(_, end)| end as u64),
        exon_records,
        Some(TRANSCRIPT_CONTIG.to_string()),
        Some(span_start),
        Some(span_end),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );

    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(TRANSCRIPT_CONTIG, contig);
    provider.add_transcript(transcript);
    provider
}

// ---------------------------------------------------------------------------
// Edit kinds and members
// ---------------------------------------------------------------------------

/// The edit types a member can take.
///
/// `Repeat` is only emitted where the reference genuinely holds a tandem array,
/// which is why it has its own stratum: `hgvs_to_spdi` refuses to expand a
/// repeat whose unit is not spelled out, and refuses a spelled unit that does not
/// match the span, so a repeat member drawn against arbitrary sequence would be
/// dropped rather than measured.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Kind {
    /// `del`.
    Del,
    /// `ins`.
    Ins,
    /// `delins`.
    Delins,
    /// A single-base substitution.
    Sub,
    /// `dup`.
    Dup,
    /// `inv`.
    Inv,
    /// A tandem repeat count, `<span><unit>[<n>]`.
    Repeat,
}

/// The kinds the dense strata pair exhaustively. `Repeat` is excluded; see
/// [`Kind`].
pub const PAIRED_KINDS: &[Kind] = &[
    Kind::Del,
    Kind::Ins,
    Kind::Delins,
    Kind::Sub,
    Kind::Dup,
    Kind::Inv,
];

impl Kind {
    /// Stable label, used in row ids and censuses.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Del => "del",
            Self::Ins => "ins",
            Self::Delins => "delins",
            Self::Sub => "sub",
            Self::Dup => "dup",
            Self::Inv => "inv",
            Self::Repeat => "repeat",
        }
    }

    /// Reference bases the member occupies for a design of `payload` size.
    fn span(self, payload: usize) -> usize {
        match self {
            Self::Ins => 0,
            Self::Sub => 1,
            Self::Del | Self::Delins | Self::Dup | Self::Inv | Self::Repeat => payload,
        }
    }
}

/// One member of a design, in 0-based served coordinates.
#[derive(Debug, Clone)]
struct Member {
    kind: Kind,
    /// Start of the reference footprint. For [`Kind::Ins`] this is the interbase
    /// the payload goes into and the footprint is empty.
    start: usize,
    /// Reference bases occupied. Zero for [`Kind::Ins`].
    span: usize,
    /// Inserted bases, for the kinds that state them.
    payload: String,
    /// Repeat unit and copy count, for [`Kind::Repeat`].
    repeat: Option<(String, usize)>,
}

/// Deterministic replacement bases that differ from the reference base for base,
/// so a substitution's alt is never its ref and a `delins` is never an identity.
fn payload_bases(served: &str, at: usize, size: usize) -> String {
    let bytes = served.as_bytes();
    (0..size)
        .map(|k| match bytes.get(at + k).copied().unwrap_or(b'A') {
            b'A' => 'C',
            b'C' => 'G',
            b'G' => 'T',
            _ => 'A',
        })
        .collect()
}

// ---------------------------------------------------------------------------
// The ground-truth applier — independent of the normalizer
// ---------------------------------------------------------------------------

/// Each member's `(position, deletion, insertion)` in the served frame, in
/// authored order. `None` when the description does not parse or a member has no
/// SPDI triple.
///
/// The whole corpus rests on this being normalizer-independent: `hgvs_to_spdi`
/// resolves coordinates and reads reference bases, and does not consult the
/// shuffler.
#[must_use]
pub fn triples_of(
    provider: &MockProvider,
    descriptor: &str,
) -> Option<Vec<(usize, String, String)>> {
    let members: Vec<HgvsVariant> = match parse_hgvs(descriptor).ok()? {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let mut triples = Vec::with_capacity(members.len());
    for member in &members {
        let triple = hgvs_to_spdi(member, provider).ok()?;
        triples.push((
            usize::try_from(triple.position).ok()?,
            triple.deletion.clone(),
            triple.insertion.clone(),
        ));
    }
    Some(triples)
}

/// Apply `triples` to `reference`, or decline.
///
/// The three rules that make this an oracle rather than a formatter, carried over
/// from `tests/it/common/cis_apply_oracle.rs::apply_reason`:
///
/// - a 3'→5' walk with a `claimed` cursor, so an overlapping description is
///   declined rather than double-spliced;
/// - a longer-deletion-first tie-break, without which a zero-width member flush
///   against a deletion reads as an overlap;
/// - rejection of two pure insertions at one interbase, whose joint denotation is
///   undefined.
#[must_use]
pub fn apply_triples(reference: &str, triples: &[(usize, String, String)]) -> Option<String> {
    let mut ordered: Vec<&(usize, String, String)> = triples.iter().collect();
    ordered.sort_by_key(|t| (std::cmp::Reverse(t.0), std::cmp::Reverse(t.1.len())));
    let bytes = reference.as_bytes();
    let mut edited = bytes.to_vec();
    let mut claimed = reference.len();
    let mut insertion_at: Option<usize> = None;
    for (position, deletion, insertion) in ordered {
        let end = position.checked_add(deletion.len())?;
        if end > reference.len() || end > claimed {
            return None;
        }
        if deletion.is_empty() && insertion_at == Some(*position) {
            return None;
        }
        if !bytes[*position..end].eq_ignore_ascii_case(deletion.as_bytes()) {
            return None;
        }
        edited.splice(*position..end, insertion.bytes());
        if deletion.is_empty() {
            insertion_at = Some(*position);
        }
        claimed = *position;
    }
    String::from_utf8(edited).ok()
}

/// What the oracle can say about `descriptor`.
///
/// The distinction between [`Self::NoSequence`] and [`Self::Inexpressible`] is
/// load-bearing and collapsing it produced a wrong headline once already: an
/// intronic `c.` position has no SPDI triple at all ("SPDI is positional and has
/// no offset notation"), so a single `Option<String>` reported 381 outputs as
/// "denoting no sequence" when the great majority were outputs that had *left the
/// transcript* — a different and separately-citable defect
/// (`general.md:44`, `checklist.md:20`), not two members claiming one territory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Denotation {
    /// The description denotes this sequence.
    Sequence(String),
    /// `parse_hgvs` rejected it.
    Unparseable,
    /// A member has no SPDI triple on this axis — an intronic offset, or a repeat
    /// whose unit is not spelled. A limit of the instrument, **not** a verdict on
    /// the description. For an output whose input was exonic it does carry
    /// information: the description has left the transcript.
    Inexpressible,
    /// Every member has a triple, but applying them overlaps or is ordered
    /// undefinably, so the description denotes no single sequence.
    NoSequence,
}

/// What `descriptor` denotes on `reference`, with the reason when it denotes
/// nothing.
#[must_use]
pub fn denotation_of(provider: &MockProvider, reference: &str, descriptor: &str) -> Denotation {
    if parse_hgvs(descriptor).is_err() {
        return Denotation::Unparseable;
    }
    let Some(triples) = triples_of(provider, descriptor) else {
        return Denotation::Inexpressible;
    };
    match apply_triples(reference, &triples) {
        Some(sequence) => Denotation::Sequence(sequence),
        None => Denotation::NoSequence,
    }
}

/// The sequence `descriptor` denotes on `reference`, or `None` when it denotes
/// none. A convenience over [`denotation_of`] for callers that only need the
/// yes/no.
#[must_use]
pub fn denoted_by(provider: &MockProvider, reference: &str, descriptor: &str) -> Option<String> {
    match denotation_of(provider, reference, descriptor) {
        Denotation::Sequence(sequence) => Some(sequence),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Equivalent placements
// ---------------------------------------------------------------------------

/// The extreme equivalent placements of one member's edit, furthest 5' and
/// furthest 3' of where it is written.
///
/// A placement is equivalent when applying it *alone* yields the sequence the
/// original does alone. Brute force over offsets within [`SHIFT_SEARCH`], for the
/// reason given on that constant. An unambiguous member yields none.
fn equivalent_placements(
    served: &str,
    position: usize,
    deletion: &str,
    insertion: &str,
) -> Vec<(usize, String, String)> {
    let Some(target) = splice(served, position, deletion.len(), insertion) else {
        return Vec::new();
    };
    let mut lowest: Option<(usize, String, String)> = None;
    let mut highest: Option<(usize, String, String)> = None;
    for delta in -SHIFT_SEARCH..=SHIFT_SEARCH {
        if delta == 0 {
            continue;
        }
        let Ok(candidate) = usize::try_from(position as isize + delta) else {
            continue;
        };
        let end = candidate + deletion.len();
        if end > served.len() {
            continue;
        }
        if !target.starts_with(&served[..candidate]) || !target.ends_with(&served[end..]) {
            continue;
        }
        let tail = served.len() - end;
        if target.len() < candidate + tail {
            continue;
        }
        let new_deletion = served[candidate..end].to_string();
        let new_insertion = target[candidate..target.len() - tail].to_string();
        if splice(served, candidate, new_deletion.len(), &new_insertion).as_deref() != Some(&target)
        {
            continue;
        }
        let placement = (candidate, new_deletion, new_insertion);
        if delta < 0 {
            if lowest.is_none() {
                lowest = Some(placement);
            }
        } else {
            highest = Some(placement);
        }
    }
    lowest.into_iter().chain(highest).collect()
}

/// `sequence[..at] + insertion + sequence[at + deleted..]`, or `None` when the
/// span runs off the end.
fn splice(sequence: &str, at: usize, deleted: usize, insertion: &str) -> Option<String> {
    let end = at.checked_add(deleted)?;
    if end > sequence.len() {
        return None;
    }
    Some(format!(
        "{}{insertion}{}",
        &sequence[..at],
        &sequence[end..]
    ))
}

// ---------------------------------------------------------------------------
// Sequences
// ---------------------------------------------------------------------------

/// Deterministic cores, two per seed (an `AT` and an `ACGT` alphabet).
///
/// The same xorshift64 as `tests/it/common/cis_apply_oracle.rs::sweep_sequences`
/// and `examples/dump_normalized_corpus.rs::corpus_sequences`, with the draw
/// length lifted into a parameter. Prefix-stable on both axes: a smaller `seeds`
/// is a strict prefix of a larger one, and a shorter `length` is a strict prefix
/// of each core, because the stream is re-seeded per `(seed, alphabet)` and
/// consumed one base at a time.
///
/// Prefix stability is load-bearing, not cosmetic: it is what makes a reduced
/// run a strict *subset* of a full run's cases, so a zero measured at a prefix
/// cannot be non-zero at the full corpus.
#[must_use]
pub fn corpus_cores(seeds: u32, length: usize) -> Vec<String> {
    let mut cores = Vec::with_capacity(2 * seeds as usize);
    for seed in 0..seeds {
        for alphabet in [b"AT".as_slice(), b"ACGT".as_slice()] {
            let mut state = u64::from(seed).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
            cores.push(
                (0..length)
                    .map(|_| {
                        state ^= state << 13;
                        state ^= state >> 7;
                        state ^= state << 17;
                        alphabet[(state % alphabet.len() as u64) as usize] as char
                    })
                    .collect(),
            );
        }
    }
    cores
}

/// Core length for the dense strata.
///
/// Long enough for four members of payload 4 separated by eight unchanged bases,
/// plus shift room either side, plus a three-exon layout with a 5'UTR and a
/// 3'UTR at each end.
pub const DENSE_CORE_LEN: usize = 96;

/// A core holding a tandem array of `copies` copies of a `unit_len`-base unit,
/// flanked by non-matching sequence so the array's extent is unambiguous.
///
/// The flanks are drawn from a different xorshift stream than the unit, and the
/// unit is rotated so it cannot be a homopolymer for `unit_len > 1`.
#[must_use]
pub fn repeat_core(unit_len: usize, copies: usize, seed: u32) -> Option<(String, String, usize)> {
    if unit_len == 0 || copies < 2 {
        return None;
    }
    let unit: String = "ACGTAC"
        .chars()
        .skip(seed as usize % 3)
        .take(unit_len)
        .collect();
    if unit.len() != unit_len {
        return None;
    }
    // A homopolymer unit of length > 1 is really a shorter unit repeated, which
    // makes the array's unit ambiguous; reject rather than measure it.
    if unit_len > 1
        && unit
            .chars()
            .all(|c| c == unit.chars().next().unwrap_or('A'))
    {
        return None;
    }
    let flank_len = 24usize;
    let flanks = corpus_cores(seed + 7, flank_len * 2);
    let flank = flanks.into_iter().nth(1)?;
    let array = unit.repeat(copies);
    // Break any accidental continuation of the array into the flanks.
    let left: String = flank[..flank_len].to_string();
    let right: String = flank[flank_len..].to_string();
    let core = format!("{left}{array}{right}");
    Some((core, unit, left.len()))
}

// ---------------------------------------------------------------------------
// The corpus
// ---------------------------------------------------------------------------

/// Why a design produced no row.
///
/// Every variant is a *legitimate* outcome of a deliberately-adversarial
/// enumeration, which is precisely why they are counted rather than skipped: a
/// generator whose designs all collapsed into one of these would produce an empty
/// corpus that reads as "nothing to find".
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DropReason {
    /// The design ran past the end of the reference, or a region was too narrow
    /// for it.
    OutOfRange,
    /// The authored spelling denotes no single sequence, in a stratum that
    /// expected one.
    NoDenotedSequence,
    /// The members cancel out and the design denotes the reference, so every
    /// spelling of it is a family about no variant.
    DenotesTheReference,
    /// Fewer than two spellings survived the ground-truth filter, so confluence
    /// is not decidable.
    Singleton,
    /// An `inv` whose span is its own reverse complement, which denotes nothing.
    PalindromicInversion,
    /// A shape the axis cannot express — an intronic offset where the reference
    /// has no introns, or a region a genomic frame does not have.
    UnavailableOnThisAxis,
    /// A geometry that does not conflict, reached from the conflict stratum.
    ///
    /// Recorded rather than skipped because the alternative was a bare
    /// `_ => continue`, which produced nothing and said nothing. That arm was
    /// unreachable while `GEOMETRIES` listed exactly the four conflicting
    /// geometries — but it is the arm someone hits by adding a fifth to that
    /// list, and the failure it would have produced is the silent one this
    /// module exists to make impossible: a design that vanishes without
    /// appearing in `designs_considered` or `dropped_by_reason`.
    NotAConflictingGeometry,
}

impl DropReason {
    /// Stable label, grouped in the census and in the ledger's `dropped_by_reason`.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::OutOfRange => "out of range",
            Self::NoDenotedSequence => "no denoted sequence",
            Self::DenotesTheReference => "denotes the reference",
            Self::Singleton => "fewer than two spellings",
            Self::PalindromicInversion => "palindromic inversion",
            Self::UnavailableOnThisAxis => "shape unavailable on this axis",
            Self::NotAConflictingGeometry => "not a conflicting geometry",
        }
    }
}

impl fmt::Display for DropReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.label())
    }
}

/// One enumerated design: either a row, or the reason it produced none.
///
/// Returned as a `Vec` of these rather than pre-folded so a consumer can route
/// the population through [`CaptureLedger`](crate::conformance::completeness::CaptureLedger)
/// — the *last* fallible step before an artifact is written, per that module's
/// contract.
pub type Attempt = Result<Row, (String, DropReason)>;

/// A corpus, folded from [`enumerate`].
#[derive(Debug, Clone)]
pub struct SpecCorpus {
    /// The bounds it was enumerated over.
    pub bounds: CorpusBounds,
    /// Designs the enumeration proposed, before any filtering.
    pub designs_considered: usize,
    /// Drops, by reason.
    pub drops: BTreeMap<&'static str, usize>,
    /// The rows.
    pub rows: Vec<Row>,
}

impl SpecCorpus {
    /// Fold [`enumerate`]'s attempts.
    #[must_use]
    pub fn from_attempts(bounds: CorpusBounds, attempts: Vec<Attempt>) -> Self {
        let designs_considered = attempts.len();
        let mut drops: BTreeMap<&'static str, usize> = BTreeMap::new();
        let mut rows = Vec::new();
        for attempt in attempts {
            match attempt {
                Ok(row) => rows.push(row),
                Err((_, reason)) => *drops.entry(reason.label()).or_default() += 1,
            }
        }
        rows.sort_by(|a, b| a.id.cmp(&b.id));
        Self {
            bounds,
            designs_considered,
            drops,
            rows,
        }
    }

    /// Rows with more than one authored member — the priority axis.
    #[must_use]
    pub fn multi_member_rows(&self) -> usize {
        self.rows.iter().filter(|row| row.is_multi_member()).count()
    }

    /// Spellings across every row, which is the number of normalizations a
    /// consumer will perform.
    #[must_use]
    pub fn spellings(&self) -> usize {
        self.rows.iter().map(|row| row.spellings.len()).sum()
    }

    /// Rows of each kind.
    #[must_use]
    pub fn by_kind(&self) -> BTreeMap<RowKind, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.kind).or_default() += 1;
        }
        counts
    }

    /// Rows per stratum.
    #[must_use]
    pub fn by_stratum(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.stratum).or_default() += 1;
        }
        counts
    }

    /// Rows per combining mechanism, and how many of each the corpus counts as
    /// genuinely multi-member.
    ///
    /// Reported per mechanism rather than as one share because the two
    /// look-alike mechanisms — [`Mechanism::CompositePayload`] and
    /// [`Mechanism::RepeatCount`] — are exactly what a bracket scan would
    /// miscount, so a reader has to be able to see them separately.
    #[must_use]
    pub fn by_mechanism(&self) -> BTreeMap<&'static str, RuleCoverage> {
        let mut counts: BTreeMap<&'static str, RuleCoverage> = BTreeMap::new();
        for row in &self.rows {
            let entry = counts.entry(row.mechanism.label()).or_default();
            entry.rows += 1;
            if row.is_multi_member() {
                entry.multi_member_rows += 1;
            }
        }
        counts
    }

    /// [`RowKind::Prohibited`] rows per [`Strength`].
    #[must_use]
    pub fn by_prohibition_strength(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            if let Some((_, strength)) = row.prohibition {
                *counts.entry(strength.label()).or_default() += 1;
            }
        }
        counts
    }

    /// Rows carrying each negative guard.
    #[must_use]
    pub fn by_negative_guard(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            for guard in &row.negative_guards {
                *counts.entry(*guard).or_default() += 1;
            }
        }
        counts
    }

    /// Rows per member geometry (#1456).
    #[must_use]
    pub fn by_geometry(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.geometry.label()).or_default() += 1;
        }
        counts
    }

    /// Rows per transcript region (#1478).
    #[must_use]
    pub fn by_region(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.region.label()).or_default() += 1;
        }
        counts
    }

    /// Rows per reference shape.
    #[must_use]
    pub fn by_shape(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.shape.label()).or_default() += 1;
        }
        counts
    }

    /// Rows per scale band (#1460). A row with no band crossed is counted under
    /// `below-all`.
    #[must_use]
    pub fn by_scale_band(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            if row.scale_bands.is_empty() {
                *counts.entry("below-all").or_default() += 1;
            }
            for band in &row.scale_bands {
                *counts.entry(*band).or_default() += 1;
            }
        }
        counts
    }

    /// For each rule tag: how many rows exercise it, and how many of those are
    /// multi-member.
    ///
    /// The second half is the point. A rule only ever exercised on a
    /// single-member row is a gap even though it reads as covered, and
    /// multi-member is where partitioning, separation, ordering, overlap
    /// detection, coalescing and typing interact.
    #[must_use]
    pub fn by_rule(&self) -> BTreeMap<&'static str, RuleCoverage> {
        let mut counts: BTreeMap<&'static str, RuleCoverage> = BTreeMap::new();
        for row in &self.rows {
            for rule in &row.rules {
                let entry = counts.entry(*rule).or_default();
                entry.rows += 1;
                if row.is_multi_member() {
                    entry.multi_member_rows += 1;
                }
            }
        }
        counts
    }
}

/// How thoroughly one rule tag is exercised.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RuleCoverage {
    /// Rows carrying the tag.
    pub rows: usize,
    /// Of those, rows with more than one authored member.
    pub multi_member_rows: usize,
}

// ---------------------------------------------------------------------------
// Enumeration
// ---------------------------------------------------------------------------

/// Enumerate every design in the declared shape space, deterministically.
///
/// The order is fixed by the loop nesting and does not depend on hashing, so two
/// runs at the same [`CorpusBounds`] produce byte-identical output.
#[must_use]
pub fn enumerate(bounds: &CorpusBounds) -> Vec<Attempt> {
    let mut attempts = Vec::new();
    enumerate_members(bounds, &mut attempts);
    enumerate_regions(bounds, &mut attempts);
    enumerate_scale(bounds, &mut attempts);
    enumerate_repeats(bounds, &mut attempts);
    enumerate_conflicts(bounds, &mut attempts);
    enumerate_intronic(bounds, &mut attempts);
    enumerate_mechanisms(bounds, &mut attempts);
    enumerate_prohibited(bounds, &mut attempts);
    enumerate_ambiguity(bounds, &mut attempts);
    attempts
}

/// Build a corpus at `bounds`.
#[must_use]
pub fn corpus(bounds: &CorpusBounds) -> SpecCorpus {
    SpecCorpus::from_attempts(bounds.clone(), enumerate(bounds))
}

/// The rule tags a design exercises, derived from its shape.
///
/// Kept in one place so the join against the committed rule inventory has a
/// single definition, and so a reviewer can see the whole mapping at once
/// instead of hunting it across six enumerators.
fn rules_for(
    kinds: &[Kind],
    geometry: Geometry,
    separation: usize,
    region: Region,
    shape: RefShape,
    members: usize,
) -> Vec<&'static str> {
    let mut rules = vec!["general.md:41-three-prime-rule"];
    if members > 1 {
        rules.push("general.md:80-allele-semicolon-separator");
        match separation {
            0 => rules.push("DNA/delins.md:17-adjacent-members-merge"),
            1 => {
                rules.push("general.md:34-separation-split");
                rules.push("general.md:35-codon-exception");
            }
            _ => rules.push("general.md:34-separation-split"),
        }
    }
    // Only geometries that denote NO sequence claim `general.md:58`.
    //
    // `CoincidentEndpoint` used to claim it too, on a separate branch above.
    // That was a false coverage claim: the variant's own doc calls it "Legal —
    // the footprints do not intersect", `denotes_a_sequence()` returns true for
    // it, and `:58` prohibits "removing part of a reference sequence and
    // replacing it with part of the same sequence" — which a pure insertion
    // beside another member's footprint does not do. Because rule tags are what
    // the committed inventory joins against, tagging a clause on rows that do
    // not exercise it made coverage of `:58` read as satisfied by rows that
    // demonstrate nothing about it. The flattering direction again.
    //
    // `:58` is still covered: `SelfReplacement` is not in `denotes_a_sequence`,
    // so its rows reach this branch, and those rows ARE the clause's own
    // example shape (`NM_004006.2:c.[762_768del;767_774dup]`).
    if !geometry.denotes_a_sequence() {
        rules.push("general.md:58-no-self-replacement");
        rules.push("checklist.md:6-most-offended");
    }
    for kind in kinds {
        rules.push(match kind {
            Kind::Del => "DNA/deletion.md:5-del-format",
            Kind::Ins => "DNA/insertion.md:5-ins-flanking-format",
            Kind::Delins => "DNA/delins.md:5-delins-format",
            Kind::Sub => "DNA/substitution.md:5-sub-format",
            Kind::Dup => "DNA/duplication.md:5-dup-format",
            Kind::Inv => "DNA/inversion.md:5-inv-format",
            Kind::Repeat => "DNA/repeated.md:5-repeat-format",
        });
    }
    if kinds.contains(&Kind::Dup) && kinds.contains(&Kind::Ins) {
        rules.push("general.md:57-dup-outranks-ins");
    }
    if kinds.contains(&Kind::Inv) {
        rules.push("general.md:56-type-priority");
    }
    match region {
        Region::Utr5 | Region::CdsStart => rules.push("checklist.md:17-cds-numbering-from-atg"),
        Region::Utr3 | Region::CdsEnd => rules.push("checklist.md:17-cds-numbering-from-atg"),
        Region::ExonJunction1 | Region::ExonJunction2 => {
            rules.push("general.md:44-junction-exception");
        }
        Region::Intronic => rules.push("checklist.md:24-intronic-positions-not-exon-numbers"),
        Region::Anywhere | Region::MidCds => {}
    }
    if shape == RefShape::Genomic {
        rules.push("checklist.md:16-genomic-has-no-offsets");
    }
    rules.sort_unstable();
    rules.dedup();
    rules
}

/// Which thresholds a block of `block_len` reference bases and a window of
/// `window` bases straddle (#1460).
fn scale_bands(block_len: usize, window: usize) -> Vec<&'static str> {
    let mut bands = Vec::new();
    for (threshold, label) in [
        (128usize, "canonical-pad-128"),
        (256, "tie-break-sweep-256"),
        (1024, "split-block-1024"),
        (4096, "canonical-window-4096"),
        (32_768, "shift-tract-32768"),
    ] {
        if block_len >= threshold {
            bands.push(label);
        }
    }
    if window >= 4096 {
        bands.push("window-past-canonical-4096");
    }
    bands.sort_unstable();
    bands.dedup();
    bands
}

/// A design resolved against one reference, ready to be turned into a row.
struct Design<'a> {
    id: String,
    stratum: &'static str,
    frame: &'a Frame,
    core: &'a str,
    members: Vec<Member>,
    separation: usize,
    geometry: Geometry,
    region: Region,
    mechanism: Mechanism,
    /// Whether the design's *shape* is the rejected-SVD-WG010 shape. Still
    /// subject to the irreducibility check in [`build_family`]: see
    /// [`negative_guards_for`].
    negative_guard_candidate: bool,
    rules: Vec<&'static str>,
}

/// Lay out `kinds.len()` members from `start`, each separated from the next by
/// `separation` unchanged reference bases.
///
/// Returns the *reason* rather than a bare `None`, because the two reasons are
/// not interchangeable: an earlier revision attributed every failure containing
/// an `inv` to a palindrome, so a corpus with no `inv` rows at all would have
/// looked like a corpus full of palindromes. Misattributed accounting is the same
/// defect as no accounting.
fn layout(
    frame: &Frame,
    start: usize,
    kinds: &[Kind],
    payload: usize,
    separation: usize,
) -> Result<Vec<Member>, DropReason> {
    let served = frame.served();
    let mut members = Vec::with_capacity(kinds.len());
    let mut cursor = start;
    for &kind in kinds {
        let span = kind.span(payload);
        if cursor == 0 {
            // An insertion needs a base 5' of its interbase to name, and a
            // footprint at index 0 has no 5' shift room.
            return Err(DropReason::OutOfRange);
        }
        let highest = if span == 0 { cursor } else { cursor + span - 1 };
        if highest >= served.len() {
            return Err(DropReason::OutOfRange);
        }
        let bases = match kind {
            Kind::Del | Kind::Dup | Kind::Inv | Kind::Repeat => String::new(),
            Kind::Ins | Kind::Delins => payload_bases(served, cursor, payload.max(1)),
            Kind::Sub => payload_bases(served, cursor, 1),
        };
        if kind == Kind::Inv {
            // `inversion.md:5` requires "**more than one nucleotide**", and a
            // span that is its own reverse complement denotes no change at all.
            if span < 2 {
                return Err(DropReason::OutOfRange);
            }
            let span_bases = &served[cursor..cursor + span];
            if reverse_complement(span_bases) == span_bases {
                return Err(DropReason::PalindromicInversion);
            }
        }
        members.push(Member {
            kind,
            start: cursor,
            span,
            payload: bases,
            repeat: None,
        });
        cursor = cursor + span + separation;
    }
    Ok(members)
}

/// Render one member as HGVS, on the frame's axis.
fn render_member(frame: &Frame, member: &Member) -> Option<String> {
    let served = frame.served();
    let label = |index: usize| frame.label(index);
    let last = member.start + member.span.saturating_sub(1);
    Some(match member.kind {
        Kind::Del if member.span == 1 => format!("{}del", label(member.start)),
        Kind::Del => format!("{}_{}del", label(member.start), label(last)),
        Kind::Dup if member.span == 1 => format!("{}dup", label(member.start)),
        Kind::Dup => format!("{}_{}dup", label(member.start), label(last)),
        Kind::Inv => format!("{}_{}inv", label(member.start), label(last)),
        Kind::Sub => format!(
            "{}{}>{}",
            label(member.start),
            &served[member.start..member.start + 1],
            member.payload
        ),
        Kind::Delins if member.span == 1 => {
            format!("{}delins{}", label(member.start), member.payload)
        }
        Kind::Delins => format!(
            "{}_{}delins{}",
            label(member.start),
            label(last),
            member.payload
        ),
        Kind::Ins => format!(
            "{}_{}ins{}",
            label(member.start - 1),
            label(member.start),
            member.payload
        ),
        Kind::Repeat => {
            let (unit, copies) = member.repeat.as_ref()?;
            format!(
                "{}_{}{unit}[{copies}]",
                label(member.start),
                label(member.start + unit.len() - 1)
            )
        }
    })
}

/// `accession:prefix.[m1;m2;…]`, or `accession:prefix.m` for one member.
fn render_allele(frame: &Frame, rendered: &[String]) -> String {
    render_allele_as(
        frame.accession(),
        frame.shape.prefix(),
        Mechanism::Cis,
        rendered,
    )
}

/// Render an allele under a chosen accession and combining mechanism.
///
/// Split out because two shapes need a different accession from the frame's own:
/// the genomic-wrapper form `NC_(NM_):c.…`, which `checklist.md:20` makes
/// mandatory for an intronic position, and the bare form that is the prohibited
/// counterpart of it.
fn render_allele_as(
    accession: &str,
    prefix: &str,
    mechanism: Mechanism,
    rendered: &[String],
) -> String {
    if rendered.len() == 1 {
        return format!("{accession}:{prefix}.{}", rendered[0]);
    }
    match mechanism {
        // `[a];[b]` — different chromosomes (`DNA/alleles.md:17`).
        Mechanism::Trans => format!(
            "{accession}:{prefix}.{}",
            rendered
                .iter()
                .map(|member| format!("[{member}]"))
                .collect::<Vec<_>>()
                .join(";")
        ),
        // `a(;)b` — phase unknown, and `DNA/alleles.md:20` says "i.e. without
        // using `[ ]`" in as many words.
        Mechanism::UnknownPhase => {
            format!("{accession}:{prefix}.{}", rendered.join("(;)"))
        }
        _ => format!("{accession}:{prefix}.[{}]", rendered.join(";")),
    }
}

/// Render an arbitrary `(position, deletion, insertion)` triple in served
/// coordinates, choosing the plainest spelling its shape admits.
fn render_triple(
    frame: &Frame,
    position: usize,
    deletion: &str,
    insertion: &str,
) -> Option<String> {
    let label = |index: usize| frame.label(index);
    match (deletion.is_empty(), insertion.is_empty()) {
        (true, true) => None,
        (true, false) => {
            if position == 0 {
                return None;
            }
            Some(format!(
                "{}_{}ins{insertion}",
                label(position - 1),
                label(position)
            ))
        }
        (false, true) if deletion.len() == 1 => Some(format!("{}del", label(position))),
        (false, true) => Some(format!(
            "{}_{}del",
            label(position),
            label(position + deletion.len() - 1)
        )),
        (false, false) if deletion.len() == 1 => {
            Some(format!("{}delins{insertion}", label(position)))
        }
        (false, false) => Some(format!(
            "{}_{}delins{insertion}",
            label(position),
            label(position + deletion.len() - 1)
        )),
    }
}

/// Every spelling that might denote the design, before the ground-truth filter.
///
/// Six families, which between them are what makes confluence measurable without
/// a published answer:
///
/// 1. the design as authored;
/// 2. the spanning `delins` over the union block, which always exists for a
///    block with reference bases;
/// 3. each member independently moved to either end of its own ambiguous run;
/// 4. the members authored in the opposite order — order-independence is its own
///    property and is the cheapest respelling there is;
/// 5. each member with a footprint of two or more split into two *touching*
///    members, which is a different **partition** of the same variant rather
///    than a different typing of one member;
/// 6. the union block retyped as an `inv` or a `dup` where the alt block admits
///    it, since `general.md:56`'s priority list makes those the preferred forms
///    and a family that cannot reach them cannot measure whether ferro does.
fn candidate_spellings(
    design: &Design<'_>,
    triples: &[(usize, String, String)],
    block: (usize, usize),
    alt_block: &str,
) -> Vec<String> {
    let frame = design.frame;
    let served = frame.served();
    let (lo, hi) = block;
    let Some(rendered) = design
        .members
        .iter()
        .map(|member| render_member(frame, member))
        .collect::<Option<Vec<String>>>()
    else {
        return Vec::new();
    };

    let mut out = vec![render_allele(frame, &rendered)];

    // 2. The spanning form.
    if let Some(spanning) = render_triple(frame, lo, &served[lo..hi], alt_block) {
        out.push(render_allele(frame, std::slice::from_ref(&spanning)));
    }

    // 3. Each member at either end of its own run.
    for (index, triple) in triples.iter().enumerate() {
        for shifted in equivalent_placements(served, triple.0, &triple.1, &triple.2) {
            let Some(text) = render_triple(frame, shifted.0, &shifted.1, &shifted.2) else {
                continue;
            };
            let mut variant = rendered.clone();
            variant[index] = text;
            out.push(render_allele(frame, &variant));
        }
    }

    // 4. The opposite authored order.
    if rendered.len() > 1 {
        let mut reversed = rendered.clone();
        reversed.reverse();
        out.push(render_allele(frame, &reversed));
    }

    // 5. One member re-partitioned into two touching members.
    for (index, triple) in triples.iter().enumerate() {
        let (position, deletion, insertion) = triple;
        if deletion.len() < 2 {
            continue;
        }
        let cut = deletion.len() / 2;
        let insert_cut = insertion.len().min(cut);
        let Some(left) =
            render_triple(frame, *position, &deletion[..cut], &insertion[..insert_cut])
        else {
            continue;
        };
        let Some(right) = render_triple(
            frame,
            position + cut,
            &deletion[cut..],
            &insertion[insert_cut..],
        ) else {
            continue;
        };
        let mut variant = Vec::with_capacity(rendered.len() + 1);
        variant.extend_from_slice(&rendered[..index]);
        variant.push(left);
        variant.push(right);
        variant.extend_from_slice(&rendered[index + 1..]);
        out.push(render_allele(frame, &variant));
    }

    // 6. Retyping the union block.
    let reference_block = &served[lo..hi];
    if !reference_block.is_empty() {
        if alt_block == reverse_complement(reference_block) && reference_block.len() >= 2 {
            let inv = format!("{}_{}inv", frame.label(lo), frame.label(hi - 1));
            out.push(render_allele(frame, std::slice::from_ref(&inv)));
        }
        if alt_block == format!("{reference_block}{reference_block}") {
            let dup = if reference_block.len() == 1 {
                format!("{}dup", frame.label(lo))
            } else {
                format!("{}_{}dup", frame.label(lo), frame.label(hi - 1))
            };
            out.push(render_allele(frame, std::slice::from_ref(&dup)));
        }
    }

    out
}

/// Turn a design into a [`RowKind::Family`] row, or say why it did not.
fn build_family(design: Design<'_>) -> Attempt {
    let frame = design.frame;
    let served = frame.served();
    let rendered: Option<Vec<String>> = design
        .members
        .iter()
        .map(|member| render_member(frame, member))
        .collect();
    let Some(rendered) = rendered else {
        return Err((design.id, DropReason::NoDenotedSequence));
    };
    let authored = render_allele(frame, &rendered);

    let Some(triples) = triples_of(frame.provider(), &authored) else {
        return Err((design.id, DropReason::NoDenotedSequence));
    };
    let Some(denoted) = apply_triples(served, &triples) else {
        return Err((design.id, DropReason::NoDenotedSequence));
    };
    if denoted == served {
        return Err((design.id, DropReason::DenotesTheReference));
    }

    // The union block, in served coordinates.
    let lo = triples.iter().map(|t| t.0).min().unwrap_or(0);
    let hi = triples.iter().map(|t| t.0 + t.1.len()).max().unwrap_or(0);
    if hi > served.len() || lo > hi {
        return Err((design.id, DropReason::NoDenotedSequence));
    }
    let suffix = served.len() - hi;
    if denoted.len() < lo + suffix
        || !denoted.starts_with(&served[..lo])
        || !denoted.ends_with(&served[hi..])
    {
        return Err((design.id, DropReason::NoDenotedSequence));
    }
    let alt_block = denoted[lo..denoted.len() - suffix].to_string();

    let mut spellings = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for candidate in candidate_spellings(&design, &triples, (lo, hi), &alt_block) {
        if !seen.insert(candidate.clone()) {
            continue;
        }
        if denoted_by(frame.provider(), served, &candidate).as_deref() == Some(denoted.as_str()) {
            spellings.push(candidate);
        }
    }
    if spellings.len() < 2 {
        return Err((design.id, DropReason::Singleton));
    }

    // The negative guard survives only when the intervening base is
    // irreducible: see `is_svd_wg010_shape`.
    let negative_guards = if design.negative_guard_candidate
        && triples.len() == 2
        && triples
            .iter()
            .all(|t| equivalent_placements(served, t.0, &t.1, &t.2).is_empty())
    {
        vec![SVD_WG010_GUARD]
    } else {
        Vec::new()
    };

    let window = hi - lo + 2 * 128;
    Ok(Row {
        id: design.id,
        kind: RowKind::Family,
        stratum: design.stratum,
        shape: frame.shape,
        core: design.core.to_string(),
        denoted: Some(denoted),
        spellings,
        members: design.members.len(),
        separation: design.separation,
        block_len: hi - lo,
        mechanism: design.mechanism,
        prohibition: None,
        negative_guards,
        geometry: design.geometry,
        region: design.region,
        scale_bands: scale_bands(hi - lo, window),
        rules: design.rules,
    })
}

/// Stratum 1: dense multi-member designs — every ordered pair of edit kinds for
/// two members, and rotations for three and four.
///
/// This is the priority axis and is where the corpus spends most of its cells.
fn enumerate_members(bounds: &CorpusBounds, out: &mut Vec<Attempt>) {
    for (core_index, core) in corpus_cores(bounds.seeds, DENSE_CORE_LEN)
        .into_iter()
        .enumerate()
    {
        for shape in RefShape::all() {
            let frame = Frame::build(shape, &core);
            // Two members: every ordered pair of kinds, exhaustively.
            for &first in PAIRED_KINDS {
                for &second in PAIRED_KINDS {
                    for &payload in DENSE_PAYLOADS {
                        for &separation in DENSE_SEPARATIONS {
                            let geometry = if separation == 0 {
                                Geometry::FlushAdjacent
                            } else {
                                Geometry::Disjoint
                            };
                            let kinds = [first, second];
                            let id = format!(
                                "s{core_index:02}-{}-pair-{}-{}-p{payload}-sep{separation}",
                                shape.label(),
                                first.label(),
                                second.label()
                            );
                            push_design(
                                out,
                                &frame,
                                &core,
                                id,
                                "members-pairs",
                                &kinds,
                                payload,
                                separation,
                                geometry,
                                Region::MidCds,
                            );
                        }
                    }
                }
            }
            // Two members sharing an endpoint: a pure insertion at the interbase
            // immediately 5' of the second member's footprint. Legal, and the
            // shape an off-by-one in overlap detection lives in.
            for &second in PAIRED_KINDS {
                for &payload in DENSE_PAYLOADS {
                    let kinds = [Kind::Ins, second];
                    let id = format!(
                        "s{core_index:02}-{}-endpoint-ins-{}-p{payload}",
                        shape.label(),
                        second.label()
                    );
                    push_design(
                        out,
                        &frame,
                        &core,
                        id,
                        "members-endpoint",
                        &kinds,
                        payload,
                        0,
                        Geometry::CoincidentEndpoint,
                        Region::MidCds,
                    );
                }
            }
            // Three and four members: rotations and uniforms over the paired
            // kinds. The ordered-pair enumeration above does not scale to three
            // members (216 cells per parameter combination), and a rotation
            // still covers every ordered *adjacent* pair.
            for &count in &MEMBER_COUNTS[1..] {
                for pattern in kind_patterns() {
                    for &payload in DENSE_PAYLOADS {
                        for &separation in DENSE_SEPARATIONS {
                            let kinds: Vec<Kind> =
                                (0..count).map(|index| pattern.kind(index)).collect();
                            let geometry = if separation == 0 {
                                Geometry::FlushAdjacent
                            } else {
                                Geometry::Disjoint
                            };
                            let id = format!(
                                "s{core_index:02}-{}-m{count}-{}-p{payload}-sep{separation}",
                                shape.label(),
                                pattern.label()
                            );
                            push_design(
                                out,
                                &frame,
                                &core,
                                id,
                                "members-rotations",
                                &kinds,
                                payload,
                                separation,
                                geometry,
                                Region::MidCds,
                            );
                        }
                    }
                }
            }
        }
    }
}

/// How a design with three or more members assigns edit kinds.
#[derive(Debug, Clone, Copy)]
enum KindPattern {
    /// Member `j` gets `PAIRED_KINDS[(r + j) % 6]`, so over `r` every ordered
    /// adjacent pair of kinds is covered.
    Rotate(usize),
    /// Every member gets the same kind, which is what a systematic design
    /// walking one window actually looks like.
    Uniform(Kind),
}

impl KindPattern {
    fn kind(self, member: usize) -> Kind {
        match self {
            Self::Rotate(r) => PAIRED_KINDS[(r + member) % PAIRED_KINDS.len()],
            Self::Uniform(kind) => kind,
        }
    }

    fn label(self) -> String {
        match self {
            Self::Rotate(r) => format!("rot{r}"),
            Self::Uniform(kind) => format!("all-{}", kind.label()),
        }
    }
}

fn kind_patterns() -> Vec<KindPattern> {
    let mut patterns: Vec<KindPattern> = (0..PAIRED_KINDS.len()).map(KindPattern::Rotate).collect();
    patterns.extend(PAIRED_KINDS.iter().copied().map(KindPattern::Uniform));
    patterns
}

/// Lay a design out and push the resulting attempt.
#[allow(clippy::too_many_arguments)]
fn push_design(
    out: &mut Vec<Attempt>,
    frame: &Frame,
    core: &str,
    id: String,
    stratum: &'static str,
    kinds: &[Kind],
    payload: usize,
    separation: usize,
    geometry: Geometry,
    region: Region,
) {
    let width: usize = kinds
        .iter()
        .map(|kind| kind.span(payload) + separation)
        .sum::<usize>()
        + 2;
    let Some(start) = frame.region_start(region, width) else {
        out.push(Err((id, DropReason::UnavailableOnThisAxis)));
        return;
    };
    let members = match layout(frame, start, kinds, payload, separation) {
        Ok(members) => members,
        Err(reason) => {
            out.push(Err((id, reason)));
            return;
        }
    };
    let mut rules = rules_for(
        kinds,
        geometry,
        separation,
        region,
        frame.shape,
        kinds.len(),
    );
    let mechanism = if kinds.len() > 1 {
        Mechanism::Cis
    } else {
        Mechanism::Lone
    };
    let negative_guard_candidate = is_svd_wg010_shape(frame.shape, kinds, separation);
    if negative_guard_candidate {
        rules.push("rejected:svd-wg010-frameless-floor-two");
        rules.sort_unstable();
        rules.dedup();
    }
    out.push(build_family(Design {
        id,
        stratum,
        frame,
        core,
        members,
        separation,
        geometry,
        region,
        mechanism,
        negative_guard_candidate,
        rules,
    }));
}

/// The one negative guard: the rejected-SVD-WG010 shape.
///
/// A design on a **frameless** axis (`g.`, or `n.` on a non-coding transcript)
/// whose two members each consume reference and sit exactly one unchanged base
/// apart. `general.md:34` says such a pair "should be described individually and
/// **not** as a 'delins'"; merging it is rejected SVD-WG010
/// (`consultation/SVD-WG010.md:8`, "The proposal has been **rejected**"), whose
/// stated ground was that needing to know "whether the two variants are in a
/// coding sequence and affecting one amino acid" is undesirable — the condition
/// `general.md:35` actually imposes.
///
/// `general.md:35`'s codon exception cannot reach a frameless axis at all, since
/// it is conditioned on "together affecting one amino acid". So a merge here is
/// not the exception applied broadly; it is the rejected proposal.
///
/// # This is only half the test — the other half is irreducibility
///
/// A separation of one is not necessarily a separation of one. Inside a run, both
/// members can shift, and two deletions separated by one `T` in a `T`-tract are
/// the *same variant* as two adjacent deletions — merging them is the 3' rule,
/// not SVD-WG010. The first revision of this guard reported 575 violations,
/// almost all of that shape. [`build_family`] therefore keeps the guard only when
/// **neither** member has any equivalent placement, i.e. the intervening base
/// cannot be shifted away. Conservative on purpose: a false negative here costs
/// coverage, a false positive costs a wrong published number.
fn is_svd_wg010_shape(shape: RefShape, kinds: &[Kind], separation: usize) -> bool {
    let frameless = matches!(shape, RefShape::Genomic | RefShape::NonCodingMultiExon(_));
    let all_consume_reference = kinds.iter().all(|kind| *kind != Kind::Ins);
    frameless && separation == 1 && kinds.len() == 2 && all_consume_reference
}

/// The guard label a design earns once irreducibility is confirmed.
const SVD_WG010_GUARD: &str = "svd-wg010-frameless-separation-floor-of-two";

/// Stratum 2: transcript geometry (#1478) — the placements a single-exon
/// `CDS_START = 1` transcript makes structurally impossible.
fn enumerate_regions(bounds: &CorpusBounds, out: &mut Vec<Attempt>) {
    const REGIONS: &[Region] = &[
        Region::Utr5,
        Region::CdsStart,
        Region::MidCds,
        Region::ExonJunction1,
        Region::ExonJunction2,
        Region::CdsEnd,
        Region::Utr3,
    ];
    const PAIRS: &[[Kind; 2]] = &[
        [Kind::Del, Kind::Del],
        [Kind::Del, Kind::Ins],
        [Kind::Dup, Kind::Dup],
        [Kind::Sub, Kind::Sub],
        [Kind::Delins, Kind::Del],
    ];
    for (core_index, core) in corpus_cores(bounds.seeds, DENSE_CORE_LEN)
        .into_iter()
        .enumerate()
    {
        for shape in RefShape::structured() {
            let frame = Frame::build(shape, &core);
            for &region in REGIONS {
                for kinds in PAIRS {
                    for &payload in &[1usize, 2] {
                        for &separation in &[0usize, 1, 2] {
                            let geometry = if separation == 0 {
                                Geometry::FlushAdjacent
                            } else {
                                Geometry::Disjoint
                            };
                            let id = format!(
                                "s{core_index:02}-{}-{}-{}-{}-p{payload}-sep{separation}",
                                shape.label(),
                                region.label(),
                                kinds[0].label(),
                                kinds[1].label()
                            );
                            push_design(
                                out, &frame, &core, id, "regions", kinds, payload, separation,
                                geometry, region,
                            );
                        }
                    }
                }
            }
        }
    }
}

/// Stratum 3: scale (#1460) — block lengths and separations that straddle every
/// threshold the canonicalizer keys on.
///
/// Enumerated one ladder at a time against a fixed spine rather than as a cross
/// product: the two ladders are 16 and 15 values, and the 4096-base cells are
/// measured at ~35-77 ms per normalization, so a cross product would make the
/// corpus's runtime a property of a few hundred cells.
fn enumerate_scale(bounds: &CorpusBounds, out: &mut Vec<Attempt>) {
    const PATTERNS: &[[Kind; 2]] = &[
        [Kind::Del, Kind::Del],
        [Kind::Inv, Kind::Del],
        [Kind::Delins, Kind::Del],
    ];
    let mut blocks: Vec<usize> = BLOCK_LADDER.to_vec();
    if bounds.extended_scale {
        blocks.extend_from_slice(EXTENDED_BLOCK_LADDER);
    }

    for shape in [RefShape::Genomic, RefShape::CodingMultiExon(Strand::Plus)] {
        // Block ladder at a fixed separation of 1 — the value `general.md:34`
        // and `general.md:35` disagree about, so every cell is also a
        // separation-rule cell.
        for &payload in &blocks {
            let core_len = payload * 3 + 128;
            let core = corpus_cores(1, core_len).into_iter().next_back();
            let Some(core) = core else { continue };
            let frame = Frame::build(shape, &core);
            for kinds in PATTERNS {
                let id = format!(
                    "scale-{}-block{payload}-{}-{}",
                    shape.label(),
                    kinds[0].label(),
                    kinds[1].label()
                );
                push_design(
                    out,
                    &frame,
                    &core,
                    id,
                    "scale-block",
                    kinds,
                    payload,
                    1,
                    Geometry::Disjoint,
                    Region::MidCds,
                );
            }
        }
        // Separation ladder at a fixed payload of 2 — separation drives the
        // canonicalizer's *window*, which is bounded separately from the block.
        for &separation in SEPARATION_LADDER {
            let core_len = separation * 2 + 256;
            let core = corpus_cores(1, core_len).into_iter().next_back();
            let Some(core) = core else { continue };
            let frame = Frame::build(shape, &core);
            for kinds in PATTERNS {
                let id = format!(
                    "scale-{}-sep{separation}-{}-{}",
                    shape.label(),
                    kinds[0].label(),
                    kinds[1].label()
                );
                push_design(
                    out,
                    &frame,
                    &core,
                    id,
                    "scale-separation",
                    kinds,
                    2,
                    separation,
                    Geometry::Disjoint,
                    Region::MidCds,
                );
            }
        }
    }
}

/// Stratum 4: tandem repeats, over cores engineered to hold a real array.
///
/// Separate from the dense strata because `hgvs_to_spdi` refuses a repeat whose
/// unit is not spelled out, and refuses a spelled unit that does not match the
/// span — so a repeat member drawn against arbitrary sequence is dropped rather
/// than measured, and a stratum that mixed the two would report a repeat
/// coverage it does not have.
fn enumerate_repeats(bounds: &CorpusBounds, out: &mut Vec<Attempt>) {
    for seed in 0..bounds.seeds.max(1) {
        for unit_len in [1usize, 2, 3, 4] {
            for copies in [3usize, 5, 8] {
                let Some((core, unit, array_start)) = repeat_core(unit_len, copies, seed) else {
                    continue;
                };
                for shape in [RefShape::Genomic, RefShape::CodingMultiExon(Strand::Plus)] {
                    let frame = Frame::build(shape, &core);
                    let offset = frame.core_offset() + array_start;
                    for new_copies in [copies + 1, copies.saturating_sub(1)] {
                        if new_copies < 2 || new_copies == copies {
                            continue;
                        }
                        let id = format!(
                            "repeat-{}-u{unit_len}-c{copies}-to{new_copies}-s{seed}",
                            shape.label()
                        );
                        let member = Member {
                            kind: Kind::Repeat,
                            start: offset,
                            span: unit_len * copies,
                            payload: String::new(),
                            repeat: Some((unit.clone(), new_copies)),
                        };
                        // A second, disjoint member 3' of the array, so the row
                        // is multi-member: the repeat/sibling interaction is
                        // where the tract logic swallowed a sibling.
                        let sibling_start = offset + unit_len * copies + 2;
                        if sibling_start + 1 >= frame.served().len() {
                            out.push(Err((id, DropReason::OutOfRange)));
                            continue;
                        }
                        let sibling = Member {
                            kind: Kind::Del,
                            start: sibling_start,
                            span: 1,
                            payload: String::new(),
                            repeat: None,
                        };
                        let rules = rules_for(
                            &[Kind::Repeat, Kind::Del],
                            Geometry::Disjoint,
                            2,
                            Region::MidCds,
                            shape,
                            2,
                        );
                        out.push(build_family(Design {
                            id,
                            stratum: "repeats",
                            frame: &frame,
                            core: &core,
                            members: vec![member, sibling],
                            separation: 2,
                            geometry: Geometry::Disjoint,
                            region: Region::MidCds,
                            mechanism: Mechanism::Cis,
                            negative_guard_candidate: false,
                            rules,
                        }));
                    }
                }
            }
        }
    }
}

/// Stratum 5: conflicting alleles (#1456) — descriptions that denote no single
/// sequence, whose property is that they are **refused**.
///
/// Before #1456 the designed corpus could not build one of these at all, so
/// three separate changes reported `0 of 18,432` as though it were evidence of
/// neutrality.
fn enumerate_conflicts(bounds: &CorpusBounds, out: &mut Vec<Attempt>) {
    const GEOMETRIES: &[Geometry] = &[
        Geometry::Nested,
        Geometry::Overlapping,
        Geometry::CoincidentInsertions,
        Geometry::SelfReplacement,
    ];
    for (core_index, core) in corpus_cores(bounds.seeds, DENSE_CORE_LEN)
        .into_iter()
        .enumerate()
    {
        for shape in RefShape::all() {
            let frame = Frame::build(shape, &core);
            for &geometry in GEOMETRIES {
                for &payload in &[2usize, 4] {
                    let id = format!(
                        "conflict-s{core_index:02}-{}-{}-p{payload}",
                        shape.label(),
                        geometry.label()
                    );
                    let width = payload * 3 + 4;
                    let Some(start) = frame.region_start(Region::MidCds, width) else {
                        out.push(Err((id, DropReason::UnavailableOnThisAxis)));
                        continue;
                    };
                    let members = match geometry {
                        // A wide member with a narrower one strictly inside it.
                        Geometry::Nested => vec![
                            simple_member(Kind::Del, start, payload * 2),
                            simple_member(Kind::Del, start + 1, payload.max(1)),
                        ],
                        // Footprints that partially intersect.
                        Geometry::Overlapping => vec![
                            simple_member(Kind::Del, start, payload * 2),
                            simple_member(Kind::Del, start + payload, payload * 2),
                        ],
                        // Two pure insertions at one interbase have no order.
                        Geometry::CoincidentInsertions => vec![
                            Member {
                                kind: Kind::Ins,
                                start,
                                span: 0,
                                payload: payload_bases(frame.served(), start, payload),
                                repeat: None,
                            },
                            Member {
                                kind: Kind::Ins,
                                start,
                                span: 0,
                                payload: payload_bases(frame.served(), start + 1, payload),
                                repeat: None,
                            },
                        ],
                        // `general.md:58`: a description that removes part of the
                        // reference and replaces it with part of the same
                        // sequence. The spec's own example is
                        // `NM_004006.2:c.[762_768del;767_774dup]`.
                        Geometry::SelfReplacement => vec![
                            simple_member(Kind::Del, start, payload * 2),
                            simple_member(Kind::Dup, start + payload, payload * 2),
                        ],
                        // Named rather than left to a `_` arm, so this match is
                        // EXHAUSTIVE over `Geometry`. Two failure modes close
                        // together: adding a fifth entry to `GEOMETRIES` now
                        // records a drop instead of producing nothing silently,
                        // and adding a variant to `Geometry` itself stops
                        // compiling here instead of being absorbed by a
                        // catch-all. The three below do not conflict — they are
                        // the disjoint/adjacent/shared-endpoint geometries the
                        // other strata enumerate — so reaching them from the
                        // conflict stratum is a mistake worth surfacing.
                        Geometry::Disjoint
                        | Geometry::FlushAdjacent
                        | Geometry::CoincidentEndpoint => {
                            out.push(Err((id, DropReason::NotAConflictingGeometry)));
                            continue;
                        }
                    };
                    let last = members
                        .iter()
                        .map(|m| m.start + m.span.max(1))
                        .max()
                        .unwrap_or(0);
                    if last >= frame.served().len() {
                        out.push(Err((id, DropReason::OutOfRange)));
                        continue;
                    }
                    let rendered: Option<Vec<String>> = members
                        .iter()
                        .map(|member| render_member(&frame, member))
                        .collect();
                    let Some(rendered) = rendered else {
                        out.push(Err((id, DropReason::NoDenotedSequence)));
                        continue;
                    };
                    let spelling = render_allele(&frame, &rendered);
                    // The row is only a conflict if it really denotes nothing.
                    // A design the applier accepts would be a *family*, and
                    // counting it here would make the refusal census a claim
                    // about the generator rather than about the implementation.
                    if denoted_by(frame.provider(), frame.served(), &spelling).is_some() {
                        out.push(Err((id, DropReason::DenotesTheReference)));
                        continue;
                    }
                    let rules = rules_for(
                        &[Kind::Del],
                        geometry,
                        0,
                        Region::MidCds,
                        shape,
                        members.len(),
                    );
                    out.push(Ok(Row {
                        id,
                        kind: RowKind::Conflict,
                        stratum: "conflicts",
                        shape,
                        core: core.clone(),
                        denoted: None,
                        spellings: vec![spelling],
                        members: members.len(),
                        separation: 0,
                        block_len: payload * 3,
                        mechanism: Mechanism::Cis,
                        prohibition: Some((
                            "general.md:58-no-self-replacement",
                            Strength::Absolute,
                        )),
                        negative_guards: Vec::new(),
                        geometry,
                        region: Region::MidCds,
                        scale_bands: Vec::new(),
                        rules,
                    }));
                }
            }
        }
    }
}

/// A member with no payload and no repeat — `del`, `dup` and `inv`.
fn simple_member(kind: Kind, start: usize, span: usize) -> Member {
    Member {
        kind,
        start,
        span,
        payload: String::new(),
        repeat: None,
    }
}

/// Stratum 6: intronic offsets — validity and idempotency only.
///
/// `hgvs_to_spdi` declines an intronic `c.` position ("SPDI is positional and has
/// no offset notation"), so the corpus has **no sequence oracle** for these rows
/// and cannot state what they denote. Rather than drop them — the shape is a
/// documented source of defects, and a single-exon `CDS_START = 1` transcript
/// cannot express one at all (#1478) — they are carried as [`RowKind::Single`]
/// rows and the two properties that need no oracle are measured. The axis test
/// reports the other two as `VACUOUS` for this stratum rather than counting them
/// as passes.
fn enumerate_intronic(bounds: &CorpusBounds, out: &mut Vec<Attempt>) {
    for (core_index, core) in corpus_cores(bounds.seeds, DENSE_CORE_LEN)
        .into_iter()
        .enumerate()
    {
        for shape in RefShape::structured() {
            if !shape.is_multi_exon() {
                continue;
            }
            let frame = Frame::build(shape, &core);
            let boundaries: Vec<(usize, isize)> = frame
                .exons()
                .iter()
                .enumerate()
                .flat_map(|(index, &(start, end))| {
                    let mut sites = Vec::new();
                    if index + 1 < frame.exons().len() {
                        sites.push((end - 1, 1isize));
                        sites.push((end - 1, 2));
                        sites.push((end - 1, 5));
                    }
                    if index > 0 {
                        sites.push((start - 1, -1));
                        sites.push((start - 1, -2));
                    }
                    sites
                })
                .collect();
            for (index, delta) in boundaries {
                for kind in [Kind::Del, Kind::Dup, Kind::Sub, Kind::Ins] {
                    for members in [1usize, 2] {
                        let id = format!(
                            "intronic-s{core_index:02}-{}-i{index}{delta:+}-{}-m{members}",
                            shape.label(),
                            kind.label()
                        );
                        let Some(label) = frame.intronic_label(index, delta) else {
                            out.push(Err((id, DropReason::UnavailableOnThisAxis)));
                            continue;
                        };
                        let member = match kind {
                            Kind::Del => format!("{label}del"),
                            Kind::Dup => format!("{label}dup"),
                            Kind::Sub => {
                                // The intronic base is unknown to the corpus (it
                                // lives on the contig, not the transcript), so a
                                // substitution states an unverifiable reference
                                // base. Spell it as a `delins` instead, which
                                // states none.
                                format!("{label}delinsTT")
                            }
                            _ => format!("{label}_{label}insTT"),
                        };
                        // An insertion needs two flanking positions; spelling
                        // both as the same intronic position is not legal HGVS,
                        // so an intronic insertion is written against the
                        // adjacent offset instead.
                        let member = if kind == Kind::Ins {
                            let Some(next) = frame.intronic_label(
                                index,
                                if delta > 0 { delta + 1 } else { delta - 1 },
                            ) else {
                                out.push(Err((id, DropReason::UnavailableOnThisAxis)));
                                continue;
                            };
                            if delta > 0 {
                                format!("{label}_{next}insTT")
                            } else {
                                format!("{next}_{label}insTT")
                            }
                        } else {
                            member
                        };
                        let rendered = if members == 1 {
                            vec![member]
                        } else {
                            // A second, exonic member, so the row exercises the
                            // interaction rather than the offset alone.
                            let Some(exonic) = frame
                                .region_start(Region::MidCds, 2)
                                .map(|start| format!("{}del", frame.label(start)))
                            else {
                                out.push(Err((id, DropReason::UnavailableOnThisAxis)));
                                continue;
                            };
                            if delta > 0 {
                                vec![exonic, member]
                            } else {
                                vec![member, exonic]
                            }
                        };
                        // Written against the genomic wrapper, not the bare
                        // transcript accession: `checklist.md:20` makes the
                        // wrapper mandatory for an intronic position, so a bare
                        // `NM_…:c.20+2del` is an *invalid input* rather than a
                        // test case. The bare form is generated separately, as a
                        // `RowKind::Prohibited` row whose property is refusal.
                        let Some(accession) = frame.wrapped_accession() else {
                            out.push(Err((id, DropReason::UnavailableOnThisAxis)));
                            continue;
                        };
                        let spelling =
                            render_allele_as(&accession, shape.prefix(), Mechanism::Cis, &rendered);
                        if parse_hgvs(&spelling).is_err() {
                            out.push(Err((id, DropReason::UnavailableOnThisAxis)));
                            continue;
                        }
                        let rules = rules_for(
                            &[kind],
                            Geometry::Disjoint,
                            2,
                            Region::Intronic,
                            shape,
                            members,
                        );
                        out.push(Ok(Row {
                            id,
                            kind: RowKind::Single,
                            stratum: "intronic",
                            shape,
                            core: core.clone(),
                            denoted: None,
                            spellings: vec![spelling],
                            members,
                            separation: 0,
                            block_len: 1,
                            mechanism: if members > 1 {
                                Mechanism::Cis
                            } else {
                                Mechanism::Lone
                            },
                            prohibition: None,
                            negative_guards: Vec::new(),
                            geometry: Geometry::Disjoint,
                            region: Region::Intronic,
                            scale_bands: Vec::new(),
                            rules,
                        }));
                    }
                }
            }
        }
    }
}

/// Stratum 7: the combining mechanisms, enumerated by mechanism rather than by
/// bracket (see [`Mechanism`]).
///
/// Each mechanism gets its own rows so the coverage table can say which are
/// exercised. Two of the six are deliberately **single-member** shapes that look
/// multi-member, and they are here so the census shows them being counted
/// correctly rather than leaving that claim to a reader's trust.
fn enumerate_mechanisms(bounds: &CorpusBounds, out: &mut Vec<Attempt>) {
    const MECHANISMS: &[Mechanism] = &[
        Mechanism::Cis,
        Mechanism::UnknownPhase,
        Mechanism::Trans,
        Mechanism::CompositePayload,
        Mechanism::RepeatCount,
    ];
    for (core_index, core) in corpus_cores(bounds.seeds, DENSE_CORE_LEN)
        .into_iter()
        .enumerate()
    {
        for shape in RefShape::all() {
            let frame = Frame::build(shape, &core);
            for &mechanism in MECHANISMS {
                for &separation in &[1usize, 3] {
                    let id = format!(
                        "mech-s{core_index:02}-{}-{}-sep{separation}",
                        shape.label(),
                        mechanism.label()
                    );
                    let Some(start) = frame.region_start(Region::MidCds, separation + 12) else {
                        out.push(Err((id, DropReason::UnavailableOnThisAxis)));
                        continue;
                    };
                    let served = frame.served();
                    let members = match mechanism {
                        Mechanism::Cis | Mechanism::UnknownPhase | Mechanism::Trans => vec![
                            format!("{}del", frame.label(start)),
                            format!("{}del", frame.label(start + 1 + separation)),
                        ],
                        // `ins[a;b]`: one member whose payload has two fragments.
                        Mechanism::CompositePayload => vec![format!(
                            "{}_{}ins[{};{}]",
                            frame.label(start),
                            frame.label(start + 1),
                            payload_bases(served, start, 2),
                            payload_bases(served, start + 1, 2)
                        )],
                        // `<span><unit>[n]`: one member with a repeat count. The
                        // unit is read out of the reference so the span matches.
                        Mechanism::RepeatCount => {
                            let unit = &served[start..start + 1];
                            vec![format!("{}{unit}[3]", frame.label(start))]
                        }
                        Mechanism::Lone => continue,
                    };
                    let spelling =
                        render_allele_as(frame.accession(), shape.prefix(), mechanism, &members);
                    if parse_hgvs(&spelling).is_err() {
                        // A mechanism the parser declines on this axis is a
                        // *finding*, not a silent absence: it is emitted as a
                        // prohibited row so the census reports it.
                        out.push(Ok(Row {
                            id,
                            kind: RowKind::Prohibited,
                            stratum: "mechanisms",
                            shape,
                            core: core.clone(),
                            denoted: None,
                            spellings: vec![spelling],
                            members: members.len(),
                            separation,
                            block_len: 1,
                            mechanism,
                            prohibition: Some((
                                "DNA/alleles.md:20-unknown-phase-without-brackets",
                                Strength::Conditional,
                            )),
                            negative_guards: Vec::new(),
                            geometry: Geometry::Disjoint,
                            region: Region::MidCds,
                            scale_bands: Vec::new(),
                            rules: vec!["DNA/alleles.md:20-unknown-phase-without-brackets"],
                        }));
                        continue;
                    }
                    // A trans or unknown-phase allele does not denote one
                    // sequence — its two members are on different molecules, or
                    // it is not known whether they are — so there is no
                    // sequence oracle and no confluence family. Validity and
                    // idempotency are what remain measurable.
                    let denotes_one_sequence = mechanism == Mechanism::Cis
                        || mechanism == Mechanism::CompositePayload
                        || mechanism == Mechanism::RepeatCount;
                    let denoted = denotes_one_sequence
                        .then(|| denoted_by(frame.provider(), served, &spelling))
                        .flatten();
                    let rules = vec![match mechanism {
                        Mechanism::Cis => "DNA/alleles.md:16-cis-semicolon",
                        Mechanism::UnknownPhase => {
                            "DNA/alleles.md:20-unknown-phase-without-brackets"
                        }
                        Mechanism::Trans => "DNA/alleles.md:17-trans-bracket-pairs",
                        Mechanism::CompositePayload => "general.md:79-brackets-composite-insertion",
                        Mechanism::RepeatCount => "DNA/repeated.md:5-repeat-format",
                        Mechanism::Lone => continue,
                    }];
                    out.push(Ok(Row {
                        id,
                        kind: RowKind::Single,
                        stratum: "mechanisms",
                        shape,
                        core: core.clone(),
                        denoted,
                        spellings: vec![spelling],
                        members: members.len(),
                        separation,
                        block_len: 1,
                        mechanism,
                        prohibition: None,
                        negative_guards: Vec::new(),
                        geometry: Geometry::Disjoint,
                        region: Region::MidCds,
                        scale_bands: Vec::new(),
                        rules,
                    }));
                }
            }
        }
    }
}

/// One prohibited shape: how to spell it, the clause, and how strongly the
/// clause is stated.
struct Prohibition {
    slug: &'static str,
    rule: &'static str,
    strength: Strength,
    /// Renders the offending **member** — no accession and no coordinate prefix —
    /// so the same prohibition can be emitted alone and inside a cis allele. A
    /// prohibited member next to a legal sibling is where a hygiene check that
    /// runs per description rather than per member stops firing.
    ///
    /// `None` when the shape is not expressible on that axis, which is itself the
    /// point for some rows: a genomic frame has no intron to hang an offset off.
    render: fn(&Frame) -> Option<String>,
    /// Whether the shape is single-member by nature, so no cis-paired row is
    /// emitted for it.
    lone_only: bool,
}

/// Stratum 8: shapes the recommendations prohibit outright.
///
/// The property is refusal. Every row carries its clause and its [`Strength`],
/// and the axis test pins the two acceptance counts separately: an
/// [`Strength::Absolute`] acceptance is a conformance defect, a
/// [`Strength::Conditional`] one is a finding to adjudicate.
fn enumerate_prohibited(bounds: &CorpusBounds, out: &mut Vec<Attempt>) {
    const PROHIBITIONS: &[Prohibition] = &[
        // `checklist.md:16`: "genomic (`g.`) reference sequences start with
        // nucleotide 1 and can not have nucleotides with additions like a `+`,
        // `-`, or `*`."
        Prohibition {
            slug: "genomic-plus-offset",
            rule: "checklist.md:16-genomic-has-no-offsets",
            strength: Strength::Absolute,
            render: |frame| {
                (frame.shape == RefShape::Genomic).then(|| format!("{}+2del", PAD_OFFSET + 10))
            },
            lone_only: false,
        },
        Prohibition {
            slug: "genomic-star-offset",
            rule: "checklist.md:16-genomic-has-no-offsets",
            strength: Strength::Absolute,
            render: |frame| (frame.shape == RefShape::Genomic).then(|| "*10del".to_string()),
            lone_only: false,
        },
        // `checklist.md:45`: "Not correct is `c.12-14del`, this describes a
        // deletion of nucleotide -14 in the intron directly 5' of nucleotide
        // `c.12`." On a genomic axis a hyphen cannot be an intronic offset
        // either, so the shape has no legal reading at all.
        Prohibition {
            slug: "hyphen-range",
            rule: "checklist.md:45-range-is-underscore",
            strength: Strength::Absolute,
            render: |frame| {
                (frame.shape == RefShape::Genomic)
                    .then(|| format!("{}-{}del", PAD_OFFSET + 10, PAD_OFFSET + 12))
            },
            lone_only: false,
        },
        // `checklist.md:49`: "Descriptions like `g.123del3` are not allowed".
        Prohibition {
            slug: "sized-deletion-suffix",
            rule: "checklist.md:49-deletion-names-both-endpoints",
            strength: Strength::Absolute,
            render: |frame| Some(format!("{}del3", frame.label(9))),
            lone_only: false,
        },
        // `checklist.md:31`: "The format `c.52insT` is **ambiguous**, and not
        // allowed."
        Prohibition {
            slug: "single-anchor-insertion",
            rule: "checklist.md:30-insertion-needs-two-anchors",
            strength: Strength::Absolute,
            render: |frame| Some(format!("{}insT", frame.label(9))),
            lone_only: false,
        },
        // `checklist.md:33`: "Describing a variant as `c.5439_5430ins6` is not
        // allowed, the inserted sequence … should be specified."
        Prohibition {
            slug: "sized-insertion-payload",
            rule: "checklist.md:32-insertion-states-its-sequence",
            strength: Strength::Absolute,
            render: |frame| Some(format!("{}_{}ins6", frame.label(9), frame.label(10))),
            lone_only: false,
        },
        // `background/standards.md:39` footnotes `X` and `-` as "used in
        // alignment only", so neither is a base a description may state.
        Prohibition {
            slug: "alignment-only-base-x",
            rule: "standards.md:39-alignment-only-symbols",
            strength: Strength::Conditional,
            render: |frame| Some(format!("{}delinsX", frame.label(9))),
            lone_only: false,
        },
        // `general.md:96`: "spaces are *not* permitted in any HGVS description".
        Prohibition {
            slug: "internal-space",
            rule: "general.md:96-no-spaces",
            strength: Strength::Absolute,
            render: |frame| Some(format!("{}_{} del", frame.label(9), frame.label(11))),
            lone_only: false,
        },
        // `checklist.md:20` NOTE: an NM_ "can only be used to describe variants
        // in introns using a `c.` prefix when a genomic reference sequence is
        // given on which the coding DNA reference sequence is annotated".
        // Conditional, not absolute: the clause states a condition rather than
        // using prohibitive words.
        Prohibition {
            slug: "bare-transcript-intronic",
            rule: "checklist.md:20-intron-needs-a-genomic-wrapper",
            strength: Strength::Conditional,
            render: |frame| {
                let exon = frame.exons().first().copied()?;
                (frame.exons().len() > 1).then(|| format!("{}+2del", frame.label(exon.1 - 1)))
            },
            lone_only: false,
        },
        // `checklist.md:26`: "the format `c.123-65_-50` is not, it is
        // **incomplete**."
        Prohibition {
            slug: "incomplete-intronic-range",
            rule: "checklist.md:26-intronic-range-is-complete",
            strength: Strength::Absolute,
            render: |frame| {
                let exon = frame.exons().first().copied()?;
                (frame.exons().len() > 1).then(|| format!("{}+2_+5del", frame.label(exon.1 - 1)))
            },
            lone_only: false,
        },
    ];

    for (core_index, core) in corpus_cores(bounds.seeds, DENSE_CORE_LEN)
        .into_iter()
        .enumerate()
    {
        for shape in RefShape::all() {
            let frame = Frame::build(shape, &core);
            for prohibition in PROHIBITIONS {
                let pairings: &[Mechanism] = if prohibition.lone_only {
                    &[Mechanism::Lone]
                } else {
                    &[Mechanism::Lone, Mechanism::Cis]
                };
                for &mechanism in pairings {
                    let id = format!(
                        "prohibited-s{core_index:02}-{}-{}-{}",
                        shape.label(),
                        prohibition.slug,
                        mechanism.label()
                    );
                    let Some(offender) = (prohibition.render)(&frame) else {
                        out.push(Err((id, DropReason::UnavailableOnThisAxis)));
                        continue;
                    };
                    // A legal sibling three bases 3' of the offender, so the row
                    // is a cis allele whose *other* member is fine. A hygiene
                    // check that runs once per description rather than once per
                    // member stops firing exactly here.
                    let members = match mechanism {
                        Mechanism::Cis => {
                            let Some(start) = frame.region_start(Region::MidCds, 4) else {
                                out.push(Err((id, DropReason::UnavailableOnThisAxis)));
                                continue;
                            };
                            vec![offender, format!("{}del", frame.label(start))]
                        }
                        _ => vec![offender],
                    };
                    let spelling =
                        render_allele_as(frame.accession(), shape.prefix(), mechanism, &members);
                    out.push(Ok(Row {
                        id,
                        kind: RowKind::Prohibited,
                        stratum: "prohibited",
                        shape,
                        core: core.clone(),
                        denoted: None,
                        spellings: vec![spelling],
                        members: members.len(),
                        separation: 0,
                        block_len: 1,
                        mechanism,
                        prohibition: Some((prohibition.rule, prohibition.strength)),
                        negative_guards: Vec::new(),
                        geometry: Geometry::Disjoint,
                        region: Region::Anywhere,
                        scale_bands: Vec::new(),
                        rules: vec![prohibition.rule],
                    }));
                }
            }
        }
    }
}

/// The 15 nucleotide symbols the recommendations admit, per
/// `background/standards.md`.
///
/// `X` and `-` appear in that table and are footnoted at `:39` as "used in
/// alignment only", so they are **excluded** here and generated instead as
/// prohibited shapes. Getting that boundary wrong in the other direction —
/// emitting `X` as content — would make every row carrying it an invalid input
/// whose result means nothing.
pub const DNA_SYMBOLS: &[char] = &[
    'A', 'C', 'G', 'T', 'B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V', 'W', 'Y',
];

/// Stratum 9: IUPAC ambiguity codes in **inserted** payloads.
///
/// Only on the inserted side. An ambiguous *reference* base is not something a
/// synthetic reference can hold and still be applied against, so the corpus does
/// not claim to cover it — see the rule inventory's not-generatable list.
fn enumerate_ambiguity(bounds: &CorpusBounds, out: &mut Vec<Attempt>) {
    for (core_index, core) in corpus_cores(bounds.seeds, DENSE_CORE_LEN)
        .into_iter()
        .enumerate()
    {
        for shape in [RefShape::Genomic, RefShape::CodingMultiExon(Strand::Plus)] {
            let frame = Frame::build(shape, &core);
            // Skip the four unambiguous symbols: they are what every other
            // stratum already draws from.
            for &symbol in &DNA_SYMBOLS[4..] {
                let id = format!("ambiguity-s{core_index:02}-{}-{symbol}", shape.label());
                let Some(start) = frame.region_start(Region::MidCds, 8) else {
                    out.push(Err((id, DropReason::UnavailableOnThisAxis)));
                    continue;
                };
                let members = vec![
                    format!("{}delins{symbol}{symbol}", frame.label(start)),
                    format!(
                        "{}_{}ins{symbol}",
                        frame.label(start + 2),
                        frame.label(start + 3)
                    ),
                ];
                let spelling = render_allele(&frame, &members);
                let Some(denoted) = denoted_by(frame.provider(), frame.served(), &spelling) else {
                    out.push(Err((id, DropReason::NoDenotedSequence)));
                    continue;
                };
                out.push(Ok(Row {
                    id,
                    kind: RowKind::Single,
                    stratum: "ambiguity",
                    shape,
                    core: core.clone(),
                    denoted: Some(denoted),
                    spellings: vec![spelling],
                    members: 2,
                    separation: 1,
                    block_len: 4,
                    mechanism: Mechanism::Cis,
                    prohibition: None,
                    negative_guards: Vec::new(),
                    geometry: Geometry::Disjoint,
                    region: Region::MidCds,
                    scale_bands: Vec::new(),
                    rules: vec!["standards.md:15-iupac-nucleotide-symbols"],
                }));
            }
        }
    }
}

/// Rules that are about one member's own spelling, so no multi-member row can
/// exercise them any more thoroughly than a lone one.
///
/// Named rather than derived, because "this rule cannot be multi-member" is a
/// judgement and the list is short enough to review. Everything *not* here must
/// reach a multi-member row — see
/// [`tests::every_rule_tag_is_reachable_and_mostly_multi_member`].
pub const SINGLE_MEMBER_BY_NATURE: &[&str] = &[
    // A composite insertion payload is one member by definition; that is the
    // whole reason it is enumerated as its own mechanism.
    "general.md:79-brackets-composite-insertion",
];

#[cfg(test)]
mod tests {
    use super::*;

    /// The corpus must be reproducible from its bounds alone, and a smaller run
    /// must be a strict prefix of a larger one on both axes — seed count and
    /// draw length. Prefix stability is what makes a reduced run a strict subset
    /// of the full run's cases, so a zero at a prefix cannot be non-zero at the
    /// full corpus.
    #[test]
    fn cores_are_deterministic_and_prefix_stable() {
        assert_eq!(corpus_cores(2, 64), corpus_cores(2, 64));
        let few = corpus_cores(2, 64);
        let many = corpus_cores(6, 64);
        assert_eq!(few[..], many[..few.len()]);
        for (short, long) in corpus_cores(2, 20).iter().zip(corpus_cores(2, 64)) {
            assert_eq!(short.as_str(), &long[..20]);
        }
        // The same stream `sweep_sequences` draws, which is what "shares the
        // generator" has to mean to be worth claiming.
        assert_eq!(corpus_cores(1, 20)[0], "TTTTTTTTTAATATATTTTA");
        assert_eq!(corpus_cores(1, 20)[1], "CCCCCCCCTGACGTATCCTA");
    }

    /// `hgvs_to_spdi` reports positions 0-based on the sequence the axis serves —
    /// the padded contig for `g.`, the transcript for `c.`/`n.`. Everything the
    /// applier does rests on this, so it is pinned rather than assumed.
    #[test]
    fn spdi_positions_are_zero_based_on_the_served_sequence() {
        let core = corpus_cores(1, DENSE_CORE_LEN).remove(1);
        let genomic = Frame::build(RefShape::Genomic, &core);
        let at = PAD_OFFSET + 9;
        let spelling = format!("{GENOMIC_CONTIG}:g.{}del", at + 1);
        let triples = triples_of(genomic.provider(), &spelling).expect("a genomic triple");
        assert_eq!(triples[0].0, at);
        assert_eq!(triples[0].1, core[9..10]);

        let coding = Frame::build(RefShape::CodingMultiExon(Strand::Plus), &core);
        let (cds_start, _) = coding.cds().expect("a coding frame has a CDS");
        let spelling = format!("{CODING_ACCESSION}:c.1del");
        let triples = triples_of(coding.provider(), &spelling).expect("a coding triple");
        assert_eq!(triples[0].0, cds_start - 1);
    }

    /// A coding frame must reach `-n`, `n` and `*n` labels, because a
    /// single-exon `CDS_START = 1` transcript reaches only the middle one — that
    /// is #1478 exactly, and a regression to it would silently delete a third of
    /// the region stratum.
    #[test]
    fn a_coding_frame_labels_all_three_regions() {
        let core = corpus_cores(1, DENSE_CORE_LEN).remove(1);
        let frame = Frame::build(RefShape::CodingMultiExon(Strand::Plus), &core);
        let (cds_start, cds_end) = frame.cds().expect("a CDS");
        assert!(cds_start > 1, "the 5'UTR must be non-empty");
        assert!(cds_end < core.len(), "the 3'UTR must be non-empty");
        assert_eq!(frame.label(0), format!("-{}", cds_start - 1));
        assert_eq!(frame.label(cds_start - 1), "1");
        assert_eq!(frame.label(cds_end), "*1");
        assert_eq!(frame.exons().len(), 3, "three exons, so two junctions");

        // The control shape is deliberately the blind one.
        let control = Frame::build(RefShape::CodingSingleExon, &core);
        assert_eq!(control.exons().len(), 1);
        assert_eq!(control.cds().map(|c| c.0), Some(1));
        assert!(
            control.region_start(Region::Utr5, 2).is_none(),
            "the #1478 control shape must have no 5'UTR to place a row in"
        );
    }

    /// A minus-strand transcript's exon 1 must occupy the highest genomic
    /// coordinates, and its transcript sequence must still read 5'->3' in
    /// transcript order. Otherwise the frame is a plus-strand transcript wearing
    /// a label, and every minus-strand row measures the same thing twice.
    #[test]
    fn a_minus_strand_frame_reverses_the_genomic_layout() {
        let core = corpus_cores(1, DENSE_CORE_LEN).remove(1);
        let frame = Frame::build(RefShape::CodingMultiExon(Strand::Minus), &core);
        assert_eq!(frame.served(), core.as_str());
        // The projection has to agree with the transcript sequence: a `c.` row
        // in exon 2 must read the base the transcript holds, not its complement.
        let (cds_start, _) = frame.cds().expect("a CDS");
        let spelling = format!("{CODING_ACCESSION}:c.1del");
        let triples = triples_of(frame.provider(), &spelling).expect("a triple");
        assert_eq!(triples[0].1, core[cds_start - 1..cds_start]);
    }

    /// The applier declines exactly what denotes no sequence, and nothing else.
    /// It is the whole ground truth, so its two non-obvious rules get their own
    /// assertions: a zero-width member flush against a deletion is not an
    /// overlap, and two pure insertions at one interbase are declined.
    #[test]
    fn the_applier_declines_only_what_denotes_no_sequence() {
        let reference = "AAAACCCCGGGG";
        assert_eq!(
            apply_triples(
                reference,
                &[
                    (4, "CCCC".to_string(), String::new()),
                    (4, String::new(), "TT".to_string()),
                ],
            )
            .as_deref(),
            Some("AAAATTGGGG")
        );
        assert_eq!(
            apply_triples(
                reference,
                &[
                    (4, "CCCC".to_string(), String::new()),
                    (6, "CC".to_string(), String::new()),
                ],
            ),
            None
        );
        assert_eq!(
            apply_triples(
                reference,
                &[
                    (4, String::new(), "TT".to_string()),
                    (4, String::new(), "GG".to_string()),
                ],
            ),
            None
        );
        assert_eq!(
            apply_triples(reference, &[(4, "GGGG".to_string(), String::new())]),
            None
        );
    }

    /// A deletion inside a run has an ambiguous placement and both ends must be
    /// found; one outside a run has none.
    #[test]
    fn equivalent_placements_find_both_ends_of_a_run() {
        //          0123456789
        let core = "GTAAAAATCG";
        let placements = equivalent_placements(core, 4, "A", "");
        let offsets: Vec<usize> = placements.iter().map(|p| p.0).collect();
        assert_eq!(offsets, vec![2, 6], "the run spans core[2..7]");
        assert!(equivalent_placements(core, 0, "G", "").is_empty());
    }

    /// Every emitted family satisfies the contract the axis test relies on: at
    /// least two distinct spellings, each of which really does denote the row's
    /// sequence when applied independently of the normalizer.
    #[test]
    fn every_family_is_a_real_confluence_family() {
        let built = corpus(&CorpusBounds::default());
        assert!(
            built.rows.len() > 2_000,
            "the corpus collapsed to {} rows",
            built.rows.len()
        );
        let mut checked = 0usize;
        for row in &built.rows {
            if row.kind != RowKind::Family {
                continue;
            }
            let frame = row.frame();
            let expected = row.denoted.as_deref().expect("a family denotes a sequence");
            assert!(row.spellings.len() >= 2, "{} is a singleton", row.id);
            for spelling in &row.spellings {
                let applied = denoted_by(frame.provider(), frame.served(), spelling)
                    .unwrap_or_else(|| panic!("{spelling} does not apply"));
                assert_eq!(applied, expected, "{spelling} denotes a different sequence");
            }
            checked += 1;
            if checked >= 400 {
                break;
            }
        }
        assert!(checked > 0, "no families to check — the corpus is vacuous");
    }

    /// A conflict row must really denote nothing, or the refusal census would be
    /// a claim about the generator rather than about the implementation.
    #[test]
    fn every_conflict_row_denotes_no_sequence() {
        let built = corpus(&CorpusBounds::default());
        let conflicts: Vec<&Row> = built
            .rows
            .iter()
            .filter(|row| row.kind == RowKind::Conflict)
            .collect();
        assert!(
            conflicts.len() >= 20,
            "only {} conflict rows — #1456 is rebuilt",
            conflicts.len()
        );
        for row in conflicts {
            let frame = row.frame();
            assert_eq!(
                denoted_by(frame.provider(), frame.served(), &row.spellings[0]),
                None,
                "{} denotes a sequence, so it is not a conflict",
                row.id
            );
        }
    }

    /// Each of the three recorded blindnesses must be varied, and the assertion
    /// is on the corpus rather than on the generator's intent: #1456, #1460 and
    /// #1478 were each invisible to the one before it precisely because the
    /// generator believed it covered them.
    #[test]
    fn the_three_recorded_blindnesses_are_varied() {
        let built = corpus(&CorpusBounds::default());

        // #1456 — member geometry, including geometries that denote nothing.
        let geometries = built.by_geometry();
        for geometry in [
            Geometry::Disjoint,
            Geometry::FlushAdjacent,
            Geometry::CoincidentEndpoint,
            Geometry::Nested,
            Geometry::Overlapping,
            Geometry::CoincidentInsertions,
            Geometry::SelfReplacement,
        ] {
            assert!(
                geometries.get(geometry.label()).copied().unwrap_or(0) > 0,
                "no rows with geometry {} — #1456 is rebuilt",
                geometry.label()
            );
        }

        // #1460 — scale. A row must actually reach each band, not merely be
        // enumerated at a length that would.
        let bands = built.by_scale_band();
        for band in [
            "canonical-pad-128",
            "tie-break-sweep-256",
            "split-block-1024",
            "canonical-window-4096",
            "window-past-canonical-4096",
        ] {
            assert!(
                bands.get(band).copied().unwrap_or(0) > 0,
                "no rows reach {band} — #1460 is rebuilt"
            );
        }

        // #1478 — transcript geometry.
        let regions = built.by_region();
        for region in [
            Region::Utr5,
            Region::CdsStart,
            Region::MidCds,
            Region::ExonJunction1,
            Region::ExonJunction2,
            Region::CdsEnd,
            Region::Utr3,
            Region::Intronic,
        ] {
            assert!(
                regions.get(region.label()).copied().unwrap_or(0) > 0,
                "no rows in region {} — #1478 is rebuilt",
                region.label()
            );
        }
        let shapes = built.by_shape();
        for shape in RefShape::all() {
            assert!(
                shapes.get(shape.label()).copied().unwrap_or(0) > 0,
                "no rows on reference shape {}",
                shape.label()
            );
        }
    }

    /// Multi-member is the priority axis, so its share is asserted rather than
    /// hoped for. Real corpora are 0.006% multi-member; this must be
    /// overwhelmingly the opposite.
    #[test]
    fn multi_member_rows_dominate_the_corpus() {
        let built = corpus(&CorpusBounds::default());
        let share = built.multi_member_rows() as f64 / built.rows.len() as f64;
        assert!(
            share > 0.9,
            "multi-member share is {share:.4}; real corpora are 0.00006 and this corpus \
             exists to invert that"
        );
    }

    /// Every rule tag a design can carry must be reachable, and the census must
    /// say whether it is reached in a multi-member context — a rule only ever
    /// exercised single-member is a gap even when it reads as covered.
    #[test]
    fn every_rule_tag_is_reachable_and_mostly_multi_member() {
        let built = corpus(&CorpusBounds::default());
        let coverage = built.by_rule();
        assert!(!coverage.is_empty(), "no rule tags — coverage is vacuous");
        for (rule, seen) in &coverage {
            assert!(seen.rows > 0, "{rule} is tagged but reached by no row");
        }
        // A rule only ever exercised on a single-member row is a gap even when
        // it reads as covered, because partitioning, separation, ordering,
        // overlap detection, coalescing and typing interact only in a
        // multi-member allele. So the exceptions are named rather than tolerated:
        // this list is what "single-member by nature" means, and anything else
        // showing up here is a coverage gap to close.
        let single_member_only: Vec<&str> = coverage
            .iter()
            .filter(|(_, seen)| seen.multi_member_rows == 0)
            .map(|(rule, _)| *rule)
            .filter(|rule| !SINGLE_MEMBER_BY_NATURE.contains(rule))
            .collect();
        assert!(
            single_member_only.is_empty(),
            "these rules are never exercised in a multi-member row, and are not in \
             SINGLE_MEMBER_BY_NATURE: {single_member_only:?}"
        );
    }
}
