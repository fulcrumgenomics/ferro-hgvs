//! The HGVS recommendations' worked examples, W1–W95, as executable guards.
//!
//! # What this file is
//!
//! `specs/2026-08-08-hgvs-rule-catalog.md` catalogues 162 rules from the
//! vendored spec checkout and, in its Appendix A, the **95 input → published
//! answer pairs** the spec prints for itself (W1–W95). Those pairs are the only
//! rows in the corpus that need no interpretation: the spec names an input and
//! prints, in the same sentence, the description to use instead. A divergence
//! is therefore a defect or a spec self-conflict — never a matter of taste.
//!
//! This module executes every W-row that *can* be executed hermetically, and
//! **records every row that cannot, with its reason**, so the size of what is
//! not covered is itself pinned ([`UNEXECUTABLE`], [`the_roster_covers_every_one_of_the_95_worked_examples`]).
//!
//! It is a companion to `tests/it/spec_worked_examples.rs`, which already runs a
//! hand-curated twelve rows (W33–W36, W65) through all three shipped error
//! modes against real reference bases. Those rows are **not duplicated** here;
//! they are listed in [`CROSS_REFERENCED`] and a test asserts they stay out of
//! this file's tables.
//!
//! # Four hermetic fixtures, no skip path
//!
//! Every row runs with no `FERRO_MANIFEST`, no network, and no environment
//! variable, so nothing here can pass by being skipped:
//!
//! - **[`Fixture::Slice`]** — the committed reference window beside
//!   `spec_worked_examples`'s `cases.json` (`LRG_199t1`, `NM_024312.4`), i.e.
//!   **real _DMD_ bases under the accession the spec itself cites**. Nearly
//!   every `LRG_199t1`/`NM_004006.2` worked example lives on this transcript,
//!   and [`the_slice_carries_the_bases_the_spec_quotes`] re-derives the spec's
//!   own quoted context from it before any row is trusted.
//! - **[`Fixture::Genomic`]**, **[`Fixture::Coding`]**, **[`Fixture::Protein`]** —
//!   synthetic contigs, transcripts and protein sequences built from the
//!   sequence the spec *prints beside the example* (`TGT`+`GC`+`CA`,
//!   `MetLeuTrpTrpGlu`, …), under deliberately synthetic accessions
//!   (`SPEC_W02.1`, `NM_SPECW16.1`, `NP_SPECW29.1`). Never a real accession
//!   over invented bases.
//!
//! Where the accession had to change, the row carries [`Rebase::Accession`] and
//! [`every_rebased_answer_differs_from_the_published_string_only_in_its_accession`]
//! proves the answer is the spec's string with nothing but the accession
//! swapped. Where the spec's example is a coordinate on a real chromosome and
//! the fixture is a short synthetic contig, the row carries [`Rebase::Local`]
//! and the same test checks the edit token is unchanged; the *position* is then
//! pinned by [`PARTITIONS`], which is what the row is really about. Where the
//! spec states its answer as prose rather than as a description ("a change of
//! the `T` at position 4"), the row carries [`Rebase::Prose`] and a
//! [`PARTITIONS`] entry is **required**, since the position is all there is.
//!
//! # Assertion policy
//!
//! - Exact equality on the **whole** output string. Never `contains`, never a
//!   prefix.
//! - Rows about a partition additionally pin **member count and every member's
//!   span** ([`PARTITIONS`]), read back off the *parsed* output — a rendered
//!   string can look right over a wrong partition.
//! - **The forbidden half is asserted too.** [`NEGATIVE`] names the exact string
//!   ferro must not emit for a legal input; [`REJECTED`] names the spellings the
//!   parser must refuse. This is the half that was missing, and the half that
//!   catches a manufactured separation.
//! - No `#[ignore]`, and no test that can pass vacuously: a row whose fixture
//!   cannot exhibit the property it tests is in [`UNEXECUTABLE`], not in a table.
//!
//! # Where ferro does not produce the spec's answer
//!
//! Three dispositions, kept apart on purpose:
//!
//! - [`DIVERGENT`] — pinned divergences. Each records ferro's actual output, the
//!   spec's, and *why* the difference is not a defect (a spec self-conflict, or
//!   a form the spec licenses in more than one spelling). Pinned in the style of
//!   `KNOWN_DIVERGENT_INPUTS` in `hgvs_spec_normalization_tests.rs`: a fix that
//!   makes one of these conform **fails**, rather than rotting unnoticed.
//! - [`PARSE_GAPS`] — a form the spec publishes as *correct* that ferro's parser
//!   refuses. A defect, pinned so that closing it fails loudly.
//! - [`FERRO_IS_WRONG`] — the spec is explicit, ferro disagrees, and there is no
//!   conflicting clause. **[`the_protein_axis_must_split_a_separation_one_delins`]
//!   is red on purpose.** See its doc comment; it is not a regression, and it
//!   must not be "fixed" by pinning ferro's answer.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use ferro_hgvs::conformance::reference_window::{WindowFixture, WindowProvider};
use ferro_hgvs::conformance::spec_worked_examples::WINDOWS_PATH;
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::reference::provider::ReferenceProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// The vendored spec checkout's `docs/` root, relative to the crate.
const SPEC_DOCS: &str = "assets/hgvs-nomenclature/docs";

/// How many worked examples Appendix A of the rule catalog records.
const WORKED_EXAMPLE_COUNT: usize = 95;

// ---------------------------------------------------------------------------
// Row types
// ---------------------------------------------------------------------------

/// Which hermetic reference a row is served by.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Fixture {
    /// The committed spec worked-example window: real `LRG_199t1` /
    /// `NM_024312.4` bases.
    Slice,
    /// A synthetic genomic contig carrying the sequence the spec prints.
    Genomic,
    /// A synthetic transcript (`c.`/`r.`/`n.`) carrying the printed sequence.
    Coding,
    /// A synthetic protein sequence carrying the printed residues.
    Protein,
}

/// How far the row's expected string had to move off the spec's literal text.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Rebase {
    /// The expected string is the spec's, byte for byte.
    Same,
    /// The spec's string with only the accession replaced by the hermetic
    /// fixture's. Checked by
    /// [`every_rebased_answer_differs_from_the_published_string_only_in_its_accession`].
    Accession,
    /// The spec states the sequence context but numbers it on a real
    /// chromosome; the fixture is that context as its own short contig, so the
    /// coordinates are local. The edit token must still match, and the position
    /// is pinned by [`PARTITIONS`].
    Local,
    /// The spec gives its answer as prose ("a change of the `T` at position 4")
    /// rather than as a description, so there is no string to compare. The
    /// answer's *position* is what the row asserts, and [`PARTITIONS`] pins it.
    Prose,
}

/// A worked example ferro satisfies.
struct Worked {
    /// The catalog's row id (`W1` … `W95`).
    id: &'static str,
    /// `path:line`, relative to `assets/hgvs-nomenclature/docs/`.
    clause: &'static str,
    /// A byte-exact substring of that line, checked by
    /// [`every_row_quotes_the_spec_verbatim`].
    quote: &'static str,
    /// What is normalized.
    input: &'static str,
    /// The spec's published answer, verbatim.
    published: &'static str,
    /// What ferro must emit — `published`, rebased onto the hermetic fixture.
    answer: &'static str,
    rebase: Rebase,
    fixture: Fixture,
    /// Why the row is here and what it proves.
    why: &'static str,
}

/// A legal input together with a description the spec forbids for it.
struct Negative {
    id: &'static str,
    clause: &'static str,
    quote: &'static str,
    input: &'static str,
    /// The exact string ferro must **not** emit for `input`.
    forbidden: &'static str,
    fixture: Fixture,
    why: &'static str,
}

/// A spelling the spec bans outright, which the parser must refuse.
struct Rejected {
    id: &'static str,
    clause: &'static str,
    quote: &'static str,
    /// The forbidden description.
    input: &'static str,
    why: &'static str,
}

/// The partition a description asserts: how many members, and where each sits.
struct Partition {
    id: &'static str,
    clause: &'static str,
    quote: &'static str,
    input: &'static str,
    /// One `(start, end)` per member, in rendered order.
    members: &'static [(i64, i64)],
    fixture: Fixture,
    why: &'static str,
}

/// A row where ferro's answer is not the spec's, deliberately.
struct Divergence {
    id: &'static str,
    clause: &'static str,
    quote: &'static str,
    input: &'static str,
    /// The spec's answer, rebased onto the fixture.
    spec_answer: &'static str,
    /// What ferro actually emits, pinned.
    ferro_answer: &'static str,
    fixture: Fixture,
    /// The classification, and the second clause where one is in tension.
    why: &'static str,
}

/// A description the spec publishes as **correct** that ferro cannot parse.
struct ParseGap {
    id: &'static str,
    clause: &'static str,
    quote: &'static str,
    input: &'static str,
    why: &'static str,
}

/// A row this file does not execute, with the reason it cannot.
struct Unexecutable {
    id: &'static str,
    reason: &'static str,
}

// ---------------------------------------------------------------------------
// The rows ferro satisfies
// ---------------------------------------------------------------------------

/// Every worked example whose published answer ferro reproduces exactly.
const WORKED: &[Worked] = &[
    // -- type choice: consecutive changes are one delins ---------------------
    Worked {
        id: "W1",
        clause: "recommendations/DNA/substitution.md:32",
        quote: "changes involving two or more consecutive nucleotides are described as deletion-insertion (delins)",
        input: "LRG_199t1:c.[79G>T;80C>T]",
        published: "LRG_199t1:c.79_80delinsTT",
        answer: "LRG_199t1:c.79_80delinsTT",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "separation 0 always merges — the uncontested end of the merge ladder",
    },
    Worked {
        id: "W2",
        clause: "recommendations/DNA/delins.md:77",
        quote: "should be described as `g.4_5delinsTG`, i.e. a deletion/insertion (delins).",
        input: "SPEC_W02.1:g.[4G>T;5C>G]",
        published: "g.4_5delinsTG",
        answer: "SPEC_W02.1:g.4_5delinsTG",
        rebase: Rebase::Accession,
        fixture: Fixture::Genomic,
        why: "the same merge on a genomic axis, where no codon can be invoked",
    },
    Worked {
        id: "W3",
        clause: "recommendations/DNA/substitution.md:90",
        quote: "should be described as `NG_012232.1:g.12_13delinsTG`",
        input: "SPEC_W03.1:g.[12G>T;13C>G]",
        published: "NG_012232.1:g.12_13delinsTG",
        answer: "SPEC_W03.1:g.12_13delinsTG",
        rebase: Rebase::Accession,
        fixture: Fixture::Genomic,
        why: "cis phase merges …",
    },
    Worked {
        id: "W3",
        clause: "recommendations/DNA/substitution.md:91",
        quote: "When phase information is not available, the variant should be described as `NG_012232.1:g.12G>T(;)13C>G`",
        input: "SPEC_W03.1:g.12G>T(;)13C>G",
        published: "NG_012232.1:g.12G>T(;)13C>G",
        answer: "SPEC_W03.1:g.12G>T(;)13C>G",
        rebase: Rebase::Accession,
        fixture: Fixture::Genomic,
        why: "… and unknown phase does not. Same two nucleotides, same bases: the \
              merge is licensed by the asserted phase, not by the sequence. A \
              normalizer that re-derived the partition from the resulting bases \
              could not tell these two apart",
    },
    Worked {
        id: "W4",
        clause: "recommendations/RNA/delins.md:66",
        quote: "should be described as `r.5_6delinsug`, i.e. a deletion/insertion (delins).",
        input: "NM_SPECW04.1:r.[5g>u;6c>g]",
        published: "r.5_6delinsug",
        answer: "NM_SPECW04.1:r.5_6delinsug",
        rebase: Rebase::Accession,
        fixture: Fixture::Coding,
        why: "the RNA axis merges consecutive substitutions like the DNA axis",
    },
    Worked {
        id: "W5",
        clause: "recommendations/protein/substitution.md:98",
        quote: "should be described as `NP_003997.1:p.Trp24_Val25delinsCysArg`",
        input: "NP_SPECW05.1:p.[Trp24Cys;Val25Arg]",
        published: "NP_003997.1:p.Trp24_Val25delinsCysArg",
        answer: "NP_SPECW05.1:p.Trp24_Val25delinsCysArg",
        rebase: Rebase::Accession,
        fixture: Fixture::Protein,
        why: "and so does the protein axis, for adjacent residues \
              (protein/substitution.md:23). Its separation-1 sibling is the red \
              row below — the two together locate the gap exactly",
    },
    // -- type choice: preferred spellings ------------------------------------
    Worked {
        id: "W6",
        clause: "recommendations/DNA/delins.md:29",
        quote: "the recommendation is not to describe the variant as",
        input: "LRG_199t1:c.4661delAinsTC",
        published: "LRG_199t1:c.4661delinsTC",
        answer: "LRG_199t1:c.4661delinsTC",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "the deleted base is redundant and is dropped (`delins.md:28` gives \
              this c. form as the g. example's coding equivalent)",
    },
    Worked {
        id: "W7",
        clause: "recommendations/DNA/delins.md:34",
        quote: "the recommendation is not to describe the variant as",
        input: "LRG_199t1:c.6775_6777delGAGinsC",
        published: "NM_004006.2:c.6775_6777delinsC",
        answer: "LRG_199t1:c.6775_6777delinsC",
        rebase: Rebase::Accession,
        fixture: Fixture::Slice,
        why: "the same rule over a range rather than one nucleotide",
    },
    Worked {
        id: "W8",
        clause: "recommendations/DNA/duplication.md:34",
        quote: "it is **not** allowed to describe the variant as `c.19_20insT`",
        input: "LRG_199t1:c.19_20insT",
        published: "NM_004006.2:c.20dup",
        answer: "LRG_199t1:c.20dup",
        rebase: Rebase::Accession,
        fixture: Fixture::Slice,
        why: "prioritisation (general.md:57): a duplicating insertion is a dup",
    },
    Worked {
        id: "W9",
        clause: "recommendations/DNA/duplication.md:50",
        quote: "the recommendation is not to describe the variant as",
        input: "LRG_199t1:c.20_23dupTAGA",
        published: "NM_004006.2:c.20_23dup",
        answer: "LRG_199t1:c.20_23dup",
        rebase: Rebase::Accession,
        fixture: Fixture::Slice,
        why: "the duplicated sequence is redundant and is dropped",
    },
    Worked {
        id: "W10",
        clause: "recommendations/RNA/duplication.md:30",
        quote: "it is **not** allowed to describe the variant as",
        input: "NM_SPECW10.1:r.6_7insu",
        published: "r.7dup",
        answer: "NM_SPECW10.1:r.7dup",
        rebase: Rebase::Accession,
        fixture: Fixture::Coding,
        why: "prioritisation holds on the RNA axis too",
    },
    Worked {
        id: "W15",
        clause: "recommendations/DNA/inversion.md:28",
        quote: "changing `..CA`<code class=\"sub\">TCAG</code>`CCT..` to `..CA`<code class=\"sub\">CTGA</code>`CCT..`",
        input: "SPEC_W15.1:g.3_6delinsCTGA",
        published: "NC_000023.10:g.32361330_32361333inv",
        answer: "SPEC_W15.1:g.3_6inv",
        rebase: Rebase::Local,
        fixture: Fixture::Genomic,
        why: "a whole-block reverse complement is an `inv`, not a `delins` \
              (inversion.md:5, and the ruling `inversion-vs-two-delins-76-83`)",
    },
    Worked {
        id: "W38",
        clause: "recommendations/DNA/other.md:28",
        quote: "the description <code class=\"invalid\">c.123C>C</code> is not allowed.",
        input: "LRG_199t1:c.123C>C",
        published: "LRG_199t1:c.123=",
        answer: "LRG_199t1:c.123=",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "an unchanged nucleotide is `=`, never a self-substitution",
    },
    // -- merge vs split ------------------------------------------------------
    Worked {
        id: "W16",
        clause: "recommendations/DNA/insertion.md:117",
        quote: "the correct description is `c.9_32dup`",
        input: "NM_SPECW16.1:c.1_24dup",
        published: "c.9_32dup",
        answer: "NM_SPECW16.1:c.9_32dup",
        rebase: Rebase::Accession,
        fixture: Fixture::Coding,
        why: "a 24-nt block shifted 8 nt by the 3' rule — the claimed `c.1_24dup` \
              is the same variant spelled 5'-most",
    },
    Worked {
        id: "W43",
        clause: "recommendations/DNA/delins.md:88",
        quote: "The correct description of this variant is `NM_007294.3:c.2077delinsATA`.",
        input: "NM_SPECW43.1:c.[2077G>A;2077_2078insTA]",
        published: "NM_007294.3:c.2077delinsATA",
        answer: "NM_SPECW43.1:c.2077delinsATA",
        rebase: Rebase::Accession,
        fixture: Fixture::Coding,
        why: "abutting members merge unconditionally — `delins.md:89` records that \
              the sentence permitting the two-member spelling was removed by the \
              committee",
    },
    Worked {
        id: "W44",
        clause: "recommendations/RNA/delins.md:70",
        quote: "The correct description of this variant is `NM_007294.3:r.2077delinsaua`.",
        input: "NM_SPECW43.1:r.[2077g>a;2077_2078insua]",
        published: "NM_007294.3:r.2077delinsaua",
        answer: "NM_SPECW43.1:r.2077delinsaua",
        rebase: Rebase::Accession,
        fixture: Fixture::Coding,
        why: "the RNA restatement of W43, on the same synthetic transcript",
    },
    Worked {
        id: "W45",
        clause: "recommendations/DNA/substitution.md:37",
        quote: "two variants separated by one nucleotide, together affecting one amino acid, should be described as a \"delins\"",
        input: "LRG_199t1:c.[145C>T;147C>G]",
        published: "LRG_199t1:c.145_147delinsTGG",
        answer: "LRG_199t1:c.145_147delinsTGG",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "the canonical separation-1 codon case (ruling \
              `delins-codon-carve-out-gap-one`); `c.145`–`c.147` is codon 49",
    },
    Worked {
        id: "W46",
        clause: "recommendations/general.md:35",
        quote: "**exception**: two variants separated by one nucleotide, together affecting one amino acid, should be described as a \"delins\".",
        input: "LRG_199t1:c.[235A>T;237G>T]",
        published: "c.235_237delinsTAT",
        answer: "LRG_199t1:c.235_237delinsTAT",
        rebase: Rebase::Accession,
        fixture: Fixture::Slice,
        why: "the general-page statement of the same exception, on its own example",
    },
    Worked {
        id: "W47",
        clause: "recommendations/RNA/delins.md:18",
        quote: "should be described as a \"delins\" (e.g., `r.142_144delinsugg`",
        input: "NM_SPECW47.1:r.[142c>u;144a>g]",
        published: "r.142_144delinsugg",
        answer: "NM_SPECW47.1:r.142_144delinsugg",
        rebase: Rebase::Accession,
        fixture: Fixture::Coding,
        why: "the codon exception on the RNA axis; `r.142`–`r.144` is codon 48. \
              (`RNA/delins.md:41` prefers the split when either variant is a known \
              polymorphism — a population fact no reference carries, so it is not \
              executable; catalogued as conflict C4)",
    },
    Worked {
        id: "W49",
        clause: "recommendations/DNA/delins.md:47",
        quote: "**The \"delins\" format is recommended**",
        input: "LRG_199t1:c.850_901delinsTTCCTCGATGCCTG",
        published: "LRG_199t1:c.850_901delinsTTCCTCGATGCCTG",
        answer: "LRG_199t1:c.850_901delinsTTCCTCGATGCCTG",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "the spanning delins is kept, not decomposed into the aligning split \
              `:46` names as its alternative — see \
              `the_spanning_delins_and_its_published_split_are_two_fixed_points`",
    },
    Worked {
        id: "W51",
        clause: "recommendations/DNA/alleles.md:19",
        quote: "`LRG_199t1:c.[2376G>C];[3103del]` is correct",
        input: "LRG_199t1:c.[2376G>C];[3103del]",
        published: "LRG_199t1:c.[2376G>C];[3103del]",
        answer: "LRG_199t1:c.[2376G>C];[3103del]",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "two variants in trans stay two members on two alleles",
    },
    Worked {
        id: "W52",
        clause: "recommendations/DNA/alleles.md:60",
        quote: "the other allele (chromosome) contains the reference sequence, `c.2376=`",
        input: "LRG_199t1:c.[2376G>C];[2376=]",
        published: "NM_004006.2:c.[2376G>C];[2376=]",
        answer: "LRG_199t1:c.[2376G>C];[2376=]",
        rebase: Rebase::Accession,
        fixture: Fixture::Slice,
        why: "a positioned `=` on the second allele is preserved, and is not the \
              same claim as a bare `[=]` (`alleles.md:61`)",
    },
    Worked {
        id: "W55",
        clause: "recommendations/protein/alleles.md:90",
        quote: "this should be described as `p.[(Phe233Leu;Cys690Trp)]`",
        input: "NP_SPECW48.1:p.[(Phe233Leu;Cys690Trp)]",
        published: "p.[(Phe233Leu;Cys690Trp)]",
        answer: "NP_SPECW48.1:p.[(Phe233Leu;Cys690Trp)]",
        rebase: Rebase::Accession,
        fixture: Fixture::Protein,
        why: "prediction parentheses sit inside the allele brackets, per member",
    },
    Worked {
        id: "W56",
        clause: "consultation/SVD-WG010.md:45",
        quote: "separated by two nucleotides and described as separate variants, not as \"delins\"",
        input: "LRG_199t1:c.235_238delinsTAGT",
        published: "LRG_199t1:c.[235A>T;238G>T]",
        answer: "LRG_199t1:c.[235A>T;238G>T]",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "separation 2 splits — and this is the *converging* direction: fed the \
              merged spelling, ferro produces the spec's split. SVD-WG010, which \
              would have merged it, was rejected (`SVD-WG010.md:5`), so \
              `general.md:34` governs",
    },
    Worked {
        id: "W56",
        clause: "recommendations/general.md:34",
        quote: "two variants separated by one or more nucleotides should be described individually and **not** as a \"delins\".",
        input: "LRG_199t1:c.[235A>T;238G>T]",
        published: "LRG_199t1:c.[235A>T;238G>T]",
        answer: "LRG_199t1:c.[235A>T;238G>T]",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "and the split spelling is a fixed point, so the pair is confluent",
    },
    Worked {
        id: "W58",
        clause: "recommendations/general.md:34",
        quote: "two variants separated by one or more nucleotides should be described individually and **not** as a \"delins\".",
        input: "LRG_199t1:c.[992_1002del;1004T>C]",
        published: "LRG_199t1:c.[992_1002del;1004T>C]",
        answer: "LRG_199t1:c.[992_1002del;1004T>C]",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "SVD-WG010:51 would have merged a del and a substitution separated by \
              one nucleotide; it was rejected. The in-force exception \
              (`general.md:35`) is about two *variants* affecting one amino acid, \
              and ferro restricts it to the sub/unchanged/sub shape (ruling \
              `codon-carve-out-shape-restriction`), so the pair stays split",
    },
    Worked {
        id: "W59",
        clause: "recommendations/protein/delins.md:21",
        quote: "two variants separated by one or more amino acids should be described individually and not as a \"delins\".",
        input: "LRG_199t1:c.[2962T>G;2970A>T]",
        published: "LRG_199t1:c.[2962T>G;2970A>T]",
        answer: "LRG_199t1:c.[2962T>G;2970A>T]",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "the DNA half of the rejected SVD-WG010:53 example: separation 7 \
              stays split. Its protein half is the red row below",
    },
    // -- placement: the 3' rule and its exception ---------------------------
    Worked {
        id: "W60",
        clause: "background/glossary.md:10",
        quote: "HGVS describes this as a change of the `T` at position 4 (not the `T` at position 2 or 3)",
        input: "SPEC_W60.1:g.2del",
        published: "a change of the `T` at position 4",
        answer: "SPEC_W60.1:g.4del",
        rebase: Rebase::Prose,
        fixture: Fixture::Genomic,
        why: "the smallest 3'-rule example in the corpus, on its own bases \
              (`ATT`+`T`+`G`). The spec states the answer as a position rather \
              than as a description, so the position — not a string — is what \
              carries it, and it is pinned in PARTITIONS",
    },
    Worked {
        id: "W29",
        clause: "recommendations/protein/deletion.md:107",
        quote: "the description `p.Ser4del` is not correct",
        input: "NP_SPECW29.1:p.Ser4del",
        published: "p.Ser5del",
        answer: "NP_SPECW29.1:p.Ser5del",
        rebase: Rebase::Accession,
        fixture: Fixture::Protein,
        why: "the 3' rule on the protein axis: in a Ser run the most C-terminal \
              residue is the one described as deleted, and the DNA-level fact \
              that codon 4 was the one removed must not be used",
    },
    Worked {
        id: "W30",
        clause: "recommendations/protein/duplication.md:98",
        quote: "the description `p.Ser4dup` is not correct",
        input: "NP_SPECW29.1:p.Ser4dup",
        published: "p.Ser5dup",
        answer: "NP_SPECW29.1:p.Ser5dup",
        rebase: Rebase::Accession,
        fixture: Fixture::Protein,
        why: "the duplication half of the same run",
    },
    Worked {
        id: "W31",
        clause: "recommendations/protein/deletion.md:39",
        quote: "for deletions in single amino acid stretches or tandem repeats, the most C-terminal residue is arbitrarily assigned",
        input: "NP_SPECW31.1:p.Trp3del",
        published: "NP_003997.2:p.Trp4del",
        answer: "NP_SPECW31.1:p.Trp4del",
        rebase: Rebase::Accession,
        fixture: Fixture::Protein,
        why: "`MetLeuTrp`+`Trp`+`Glu`: the deletion is assigned to Trp4, not Trp3",
    },
    Worked {
        id: "W32",
        clause: "recommendations/protein/duplication.md:40",
        quote: "for duplications in single amino acid stretches or tandem repeats, the most C-terminal residue is arbitrarily assigned",
        input: "NP_SPECW31.1:p.Trp4dup",
        published: "NP_003997.2:p.Trp4dup",
        answer: "NP_SPECW31.1:p.Trp4dup",
        rebase: Rebase::Accession,
        fixture: Fixture::Protein,
        why: "and the C-terminal-most spelling is already a fixed point",
    },
    Worked {
        id: "W61",
        clause: "recommendations/DNA/deletion.md:34",
        quote: "the last A of the 8 nucleotide A-stretch running from position `c.5690` to `c.5697`",
        input: "LRG_199t1:c.5690del",
        published: "NM_004006.2:c.5697del",
        answer: "LRG_199t1:c.5697del",
        rebase: Rebase::Accession,
        fixture: Fixture::Slice,
        why: "the 3' rule inside a mononucleotide run, on the real _DMD_ A-stretch",
    },
    Worked {
        id: "W61",
        clause: "recommendations/DNA/duplication.md:39",
        quote: "the last `A` of the 8 nucleotide A-stretch running from position `c.5690` to `c.5697`",
        input: "LRG_199t1:c.5690dup",
        published: "NM_004006.2:c.5697dup",
        answer: "LRG_199t1:c.5697dup",
        rebase: Rebase::Accession,
        fixture: Fixture::Slice,
        why: "the duplication half of the same run",
    },
    Worked {
        id: "W64",
        clause: "background/numbering.md:25",
        quote: "the variant is described as `LRG_199t1:c.3921del`",
        input: "LRG_199t1:c.3921del",
        published: "LRG_199t1:c.3921del",
        answer: "LRG_199t1:c.3921del",
        rebase: Rebase::Same,
        fixture: Fixture::Slice,
        why: "the exon/exon junction exception to the 3' rule: the spec's own \
              answer is a fixed point and is not shifted into the next exon. \
              (Whether `c.3922del` should converge *back* onto it is the open \
              ruling `exon-junction-dup-converge-from-the-far-side`, whose `dup` \
              half is pinned in `spec_worked_examples.rs`)",
    },
];

// ---------------------------------------------------------------------------
// The forbidden half
// ---------------------------------------------------------------------------

/// Legal inputs, each paired with the exact description the spec forbids.
///
/// A positive assertion already implies these, but only as long as the positive
/// answer stays what it is. Naming the forbidden string makes the failure say
/// *which rule* was broken, and it is the half that catches a normalizer
/// manufacturing a separation or a merge it was not licensed to make.
const NEGATIVE: &[Negative] = &[
    Negative {
        id: "W1",
        clause: "recommendations/DNA/substitution.md:32",
        quote: "so the description <code class=\"invalid\">c.[79G>T;80C>T]</code> is not correct",
        input: "LRG_199t1:c.[79G>T;80C>T]",
        forbidden: "LRG_199t1:c.[79G>T;80C>T]",
        fixture: Fixture::Slice,
        why: "two consecutive substitutions must not survive as two members",
    },
    Negative {
        id: "W8",
        clause: "recommendations/DNA/duplication.md:34",
        quote: "it is **not** allowed to describe the variant as `c.19_20insT`",
        input: "LRG_199t1:c.19_20insT",
        forbidden: "LRG_199t1:c.19_20insT",
        fixture: Fixture::Slice,
        why: "a duplicating insertion must not stay an insertion",
    },
    Negative {
        id: "W8",
        clause: "recommendations/DNA/duplication.md:35",
        quote: "the recommendation is not to describe the variant as <code class=\"invalid\">NM_004006.2:c.20dupT</code>",
        input: "LRG_199t1:c.20dupT",
        forbidden: "LRG_199t1:c.20dupT",
        fixture: Fixture::Slice,
        why: "and the duplicated base must not be restated",
    },
    Negative {
        id: "W10",
        clause: "recommendations/RNA/duplication.md:30",
        quote: "it is **not** allowed to describe the variant as",
        input: "NM_SPECW10.1:r.6_7insu",
        forbidden: "NM_SPECW10.1:r.6_7insu",
        fixture: Fixture::Coding,
        why: "the RNA restatement of the same prohibition",
    },
    Negative {
        id: "W13",
        clause: "recommendations/DNA/inversion.md:56",
        quote: "`TTCG` is only the **complement** of `AAGC`.",
        input: "SPEC_W13.1:g.1_4delinsTTCG",
        forbidden: "SPEC_W13.1:g.1_4inv",
        fixture: Fixture::Genomic,
        why: "a complement is not an inversion. This is the negative side of the \
              `delins -> inv` canonicalisation W15 exercises: it must fire on the \
              reverse complement and on nothing else",
    },
    Negative {
        id: "W14",
        clause: "recommendations/DNA/inversion.md:62",
        quote: "`CGAA` is only the **reverse** of `AAGC`.",
        input: "SPEC_W13.1:g.1_4delinsCGAA",
        forbidden: "SPEC_W13.1:g.1_4inv",
        fixture: Fixture::Genomic,
        why: "nor is a reverse",
    },
    Negative {
        id: "W16",
        clause: "recommendations/DNA/insertion.md:118",
        quote: "seems correct but neglects the **3'rule**",
        input: "NM_SPECW16.1:c.1_24dup",
        forbidden: "NM_SPECW16.1:c.1_24dup",
        fixture: Fixture::Coding,
        why: "the 5'-most spelling of a shiftable block must not survive",
    },
    Negative {
        id: "W45",
        clause: "recommendations/DNA/delins.md:42",
        quote: "so the description <code class=\"invalid\">c.[145C>T;147C>G]</code> is not correct",
        input: "LRG_199t1:c.[145C>T;147C>G]",
        forbidden: "LRG_199t1:c.[145C>T;147C>G]",
        fixture: Fixture::Slice,
        why: "the separation-1 codon pair must not stay split. (This NOTE sits \
              under an unrelated example, `c.9002_9009delinsTTT` — catalog \
              conflict C8, a spec editing slip; the rule it states is the one \
              `substitution.md:37` attaches to the right example)",
    },
    Negative {
        id: "W56",
        clause: "consultation/SVD-WG010.md:45",
        quote: "described as separate variants, not as \"delins\" (`c.235_238delinsTAGT`)",
        input: "LRG_199t1:c.[235A>T;238G>T]",
        forbidden: "LRG_199t1:c.235_238delinsTAGT",
        fixture: Fixture::Slice,
        why: "separation 2 must not merge. SVD-WG010 proposed exactly this merge \
              and was rejected; `general.md:36-39` still announces it as \
              forthcoming (catalog conflict C2), so a normalizer written to the \
              general page's NOTE would fail here",
    },
    Negative {
        id: "W57",
        clause: "consultation/SVD-WG010.md:16",
        quote: "is described as `g.3_4delinsGGT`, not as `g.[3C>G;7dup]`",
        input: "SPEC_W57.1:g.[3C>G;7dup]",
        forbidden: "SPEC_W57.1:g.3_4delinsGGT",
        fixture: Fixture::Genomic,
        why: "the rejected proposal's flagship example, read the way the rejection \
              leaves it. SVD-WG010:15 is the only place in the corpus that would \
              measure separation over the whole shift interval rather than at the \
              canonical position; it is not in force, so a substitution at g.3 and \
              a dup at g.7 stay two members",
    },
    Negative {
        id: "W58",
        clause: "consultation/SVD-WG010.md:51",
        quote: "separated by fewer than two nucleotides and described as a \"delins\" variant",
        input: "LRG_199t1:c.[992_1002del;1004T>C]",
        forbidden: "LRG_199t1:c.992_1004delinsAC",
        fixture: Fixture::Slice,
        why: "the rejected proposal's answer must not be produced",
    },
    Negative {
        id: "W60",
        clause: "background/glossary.md:10",
        quote: "not the `T` at position 2 or 3",
        input: "SPEC_W60.1:g.2del",
        forbidden: "SPEC_W60.1:g.2del",
        fixture: Fixture::Genomic,
        why: "the spec names the two wrong positions explicitly; this pins one of \
              them, and the partition row pins the right one",
    },
    Negative {
        id: "W61",
        clause: "recommendations/DNA/deletion.md:30",
        quote: "the recommendation is not to describe the variant as",
        input: "LRG_199t1:c.5697delA",
        forbidden: "LRG_199t1:c.5697delA",
        fixture: Fixture::Slice,
        why: "the deleted base must not be restated",
    },
    Negative {
        id: "W64",
        clause: "background/numbering.md:25",
        quote: "**not as** `c.3922del`",
        input: "LRG_199t1:c.3921del",
        forbidden: "LRG_199t1:c.3922del",
        fixture: Fixture::Slice,
        why: "the junction exception must hold: a 3' shift here would cross into \
              the next exon and, projected back, move 2,790 bp on the genome",
    },
];

/// Spellings the spec bans outright: ferro must refuse them at the parser.
///
/// A prohibition ferro merely *normalizes away* still accepts the string; these
/// rows assert the stronger thing, which is what the spec's wording ("is not
/// allowed", "can not be used") states.
const REJECTED: &[Rejected] = &[
    Rejected {
        id: "W1",
        clause: "recommendations/DNA/substitution.md:33",
        quote: "this change can not be described as a substitution like",
        input: "LRG_199t1:c.79_80GC>TT",
        why: "a substitution replaces one nucleotide by one other; a multi-base \
              reference in a `>` is not a description ferro will read",
    },
    Rejected {
        id: "W2",
        clause: "recommendations/DNA/delins.md:73",
        quote: "Can I describe a `GC` to `TG` variant as a di-nucleotide substitution",
        input: "SPEC_W02.1:g.4GC>TG",
        why: "the genomic form of the same prohibition",
    },
    Rejected {
        id: "W4",
        clause: "recommendations/RNA/delins.md:65",
        quote: "By definition, a substitution changes **one** nucleotide into **one** other nucleotide",
        input: "NM_SPECW04.1:r.4gc>ug",
        why: "and the RNA form",
    },
    Rejected {
        id: "W5",
        clause: "recommendations/protein/substitution.md:96",
        quote: "No, this is not allowed.",
        input: "NP_SPECW05.1:p.TrpVal24CysArg",
        why: "and the protein form — a two-residue substitution is a delins",
    },
    Rejected {
        id: "W42",
        clause: "recommendations/DNA/duplication.md:92",
        quote: "the variant is not described using <code class=\"invalid\">dupins</code>, a format not used in HGVS nomenclature",
        input: "SPEC_W42.1:g.123_234dupinv",
        why: "an inverted duplication is an insertion carrying an `inv` payload \
              (`duplication.md:20`), never a `dup` with a suffix. The positive \
              half of W42 — whether the payload is spelled as a coordinate range \
              or as literal bases — is the same under-specification pinned at W11",
    },
    Rejected {
        id: "W50",
        clause: "recommendations/general.md:58",
        quote: "descriptions removing part of a reference sequence and replacing it with part of the same sequence are not allowed",
        input: "LRG_199t1:c.[762_768del;767_774dup]",
        why: "the spec's own forbidden example, and the strongest wording in the \
              prioritisation family. Note ferro reads `:58` narrowly — as barring \
              members that *overlap* on the reference — which is what lets it keep \
              the disjoint aligning split `delins.md:46` publishes (catalog \
              conflict C1). Both halves are asserted, so a future widening of \
              either reading breaks a test rather than passing silently",
    },
    Rejected {
        id: "W55",
        clause: "recommendations/protein/alleles.md:34",
        quote: "the description <code class=\"invalid\">p.([Ser68Arg;Asn594del])</code> is not correct",
        input: "NP_SPECW48.1:p.([Phe233Leu;Cys690Trp])",
        why: "prediction parentheses outside the allele brackets are not a legal \
              spelling",
    },
];

// ---------------------------------------------------------------------------
// Partitions
// ---------------------------------------------------------------------------

/// The partition each description asserts, pinned member by member.
///
/// Read off the **parsed** output rather than the rendered text: an output can
/// render plausibly over a wrong member count, and the member count is precisely
/// what the merge/split rules decide.
const PARTITIONS: &[Partition] = &[
    Partition {
        id: "W1",
        clause: "recommendations/DNA/delins.md:16",
        quote: "changes involving two or more consecutive nucleotides are described as deletion/insertion (delins) variants.",
        input: "LRG_199t1:c.[79G>T;80C>T]",
        members: &[(79, 80)],
        fixture: Fixture::Slice,
        why: "two members in, one member out, spanning both nucleotides",
    },
    Partition {
        id: "W2",
        clause: "recommendations/DNA/delins.md:77",
        quote: "should be described as `g.4_5delinsTG`",
        input: "SPEC_W02.1:g.[4G>T;5C>G]",
        members: &[(4, 5)],
        fixture: Fixture::Genomic,
        why: "same shape on the genomic axis",
    },
    Partition {
        id: "W3",
        clause: "recommendations/DNA/substitution.md:91",
        quote: "When phase information is not available",
        input: "SPEC_W03.1:g.12G>T(;)13C>G",
        members: &[(12, 12), (13, 13)],
        fixture: Fixture::Genomic,
        why: "unknown phase keeps two members on consecutive nucleotides — the one \
              shape `delins.md:16` does not reach, because the two changes are not \
              asserted to be on one molecule",
    },
    Partition {
        id: "W43",
        clause: "recommendations/DNA/delins.md:88",
        quote: "The correct description of this variant is `NM_007294.3:c.2077delinsATA`.",
        input: "NM_SPECW43.1:c.[2077G>A;2077_2078insTA]",
        members: &[(2077, 2077)],
        fixture: Fixture::Coding,
        why: "an abutting substitution and insertion collapse to one member",
    },
    Partition {
        id: "W45",
        clause: "recommendations/DNA/substitution.md:37",
        quote: "should be described as a \"delins\"",
        input: "LRG_199t1:c.[145C>T;147C>G]",
        members: &[(145, 147)],
        fixture: Fixture::Slice,
        why: "the merged member spans the unchanged c.146 between them",
    },
    Partition {
        id: "W46",
        clause: "recommendations/general.md:35",
        quote: "together affecting one amino acid, should be described as a \"delins\"",
        input: "LRG_199t1:c.[235A>T;237G>T]",
        members: &[(235, 237)],
        fixture: Fixture::Slice,
        why: "same, one codon further along the same transcript",
    },
    Partition {
        id: "W49",
        clause: "recommendations/DNA/delins.md:46",
        quote: "giving an alternative description like `c.[850_869del;874_881del;887_897del;901_902insG]`",
        input: "LRG_199t1:c.[850_869del;874_881del;887_897del;901_902insG]",
        members: &[(850, 869), (874, 881), (887, 897), (901, 902)],
        fixture: Fixture::Slice,
        why: "the split the spec names as legal survives with all four members and \
              their gaps of 4, 5 and 3 nt intact — nothing is merged into the hull",
    },
    Partition {
        id: "W49",
        clause: "recommendations/DNA/delins.md:47",
        quote: "**The \"delins\" format is recommended**",
        input: "LRG_199t1:c.850_901delinsTTCCTCGATGCCTG",
        members: &[(850, 901)],
        fixture: Fixture::Slice,
        why: "and the recommended spanning form stays one member",
    },
    Partition {
        id: "W51",
        clause: "recommendations/DNA/alleles.md:17",
        quote: "when two variants are identified in a gene that are **on different chromosomes** (in trans)",
        input: "LRG_199t1:c.[2376G>C];[3103del]",
        members: &[(2376, 2376), (3103, 3103)],
        fixture: Fixture::Slice,
        why: "a trans pair is two members on two alleles and never merges, whatever \
              their separation",
    },
    Partition {
        id: "W56",
        clause: "consultation/SVD-WG010.md:45",
        quote: "separated by two nucleotides and described as separate variants",
        input: "LRG_199t1:c.235_238delinsTAGT",
        members: &[(235, 235), (238, 238)],
        fixture: Fixture::Slice,
        why: "the split is where the spec puts it: at 235 and 238, not at the hull. \
              A splitter that merely re-derived a minimal partition could reach the \
              same rendered string by a different route; the spans say it did not",
    },
    Partition {
        id: "W57",
        clause: "recommendations/general.md:34",
        quote: "should be described individually and **not** as a \"delins\"",
        input: "SPEC_W57.1:g.[3C>G;7dup]",
        members: &[(3, 3), (7, 7)],
        fixture: Fixture::Genomic,
        why: "the dup stays at its 3'-most position (the T-run is g.5–g.7) and the \
              substitution stays its own member",
    },
    Partition {
        id: "W58",
        clause: "recommendations/general.md:34",
        quote: "two variants separated by one or more nucleotides should be described individually",
        input: "LRG_199t1:c.[992_1002del;1004T>C]",
        members: &[(992, 1002), (1004, 1004)],
        fixture: Fixture::Slice,
        why: "an 11-nt deletion and a substitution one nucleotide apart",
    },
    Partition {
        id: "W59",
        clause: "recommendations/general.md:34",
        quote: "should be described individually and **not** as a \"delins\"",
        input: "LRG_199t1:c.[2962T>G;2970A>T]",
        members: &[(2962, 2962), (2970, 2970)],
        fixture: Fixture::Slice,
        why: "separation 7 on the DNA axis",
    },
    Partition {
        id: "W60",
        clause: "background/glossary.md:10",
        quote: "a change of the `T` at position 4",
        input: "SPEC_W60.1:g.2del",
        members: &[(4, 4)],
        fixture: Fixture::Genomic,
        why: "the 3'-rule position the spec names, asserted as a coordinate rather \
              than as a rendered string",
    },
    Partition {
        id: "W61",
        clause: "recommendations/DNA/deletion.md:34",
        quote: "the last A of the 8 nucleotide A-stretch",
        input: "LRG_199t1:c.5690del",
        members: &[(5697, 5697)],
        fixture: Fixture::Slice,
        why: "the run is c.5690–c.5697 and the deletion lands on its 3' end",
    },
    Partition {
        id: "W64",
        clause: "background/numbering.md:23",
        quote: "the 3' rule is not applied when there is a deletion/duplication around exon/exon junctions",
        input: "LRG_199t1:c.3921del",
        members: &[(3921, 3921)],
        fixture: Fixture::Slice,
        why: "the exception holds the member at the exon's last nucleotide",
    },
];

// ---------------------------------------------------------------------------
// Recorded divergences
// ---------------------------------------------------------------------------

/// Rows where ferro's answer is not the spec's, and that is not a defect.
///
/// Pinned both ways: ferro must still produce `ferro_answer`, and it must still
/// **not** produce `spec_answer`. If a change makes one of these conform, the
/// test fails and the row must be promoted into [`WORKED`] deliberately.
const DIVERGENT: &[Divergence] = &[
    Divergence {
        id: "W11",
        clause: "recommendations/DNA/insertion.md:111",
        quote: "The variant should be described as an insertion; `g.17_18ins5_16`.",
        input: "SPEC_W11.1:g.17_18ins5_16",
        spec_answer: "SPEC_W11.1:g.17_18ins5_16",
        ferro_answer: "SPEC_W11.1:g.17_18insATCGATCGATCG",
        fixture: Fixture::Genomic,
        why: "**Spelling the spec leaves open.** `insertion.md:22` licenses the \
              inserted sequence as literal bases *or* as an in-reference \
              coordinate range, and `consultation/open-issues.md:72-78` records \
              that the recommendations \"do not specify when to use which of these \
              formats ... This is undesired\" (catalog conflict C11). Ferro \
              resolves the range against the reference and emits the literal, \
              which is the same variant in the other licensed spelling. What the \
              spec is actually settling here — that this is an insertion and not \
              a `dup`, because the copy is not directly 3'-flanking \
              (`duplication.md:17`) — ferro gets right, and that is asserted \
              below rather than left implied",
    },
    Divergence {
        id: "W53",
        clause: "recommendations/DNA/alleles.md:106",
        quote: "which can be described as `LRG_199t1:c.[76A>C];[0]`",
        input: "LRG_199t1:c.[76A>C];[0]",
        spec_answer: "LRG_199t1:c.[76A>C];[0]",
        ferro_answer: "[LRG_199t1:c.76A>C];[0]",
        fixture: Fixture::Slice,
        why: "**Render placement, not content.** Ferro puts the accession inside \
              the first bracket where the spec puts it in front; the members, \
              their phase and their positions are identical. The same deliberate \
              choice is already pinned as class 2 of `KNOWN_DIVERGENT_INPUTS` in \
              `hgvs_spec_normalization_tests.rs`; recorded here so the trans-plus- \
              null-allele shape is covered by a worked example rather than only by \
              the harvested fixture",
    },
];

/// Descriptions the spec publishes as **correct** that ferro cannot parse.
///
/// A defect, not a divergence: there is no second reading under which refusing
/// the spec's own published answer is right. Pinned as a gap so that closing it
/// fails here and the row can be promoted into [`WORKED`].
const PARSE_GAPS: &[ParseGap] = &[ParseGap {
    id: "W12",
    clause: "recommendations/protein/duplication.md:19",
    quote: "duplication may only be used when the additional copy is **directly C-terminal**",
    input: "NP_SPECW48.1:p.His7_Gln8insGly4_Ser6",
    why: "`protein/duplication.md:87-91` and `protein/insertion.md:76-80` both \
          publish `p.His7_Gln8insGly4_Ser6` as the correct description of a \
          non-tandem extra copy — the protein analogue of W11. Ferro's parser \
          rejects the residue-range payload (`Unexpected trailing characters: \
          '4_Ser6'`), so the row cannot be executed as a normalization at all",
}];

/// Rows where the spec is explicit, ferro disagrees, and no clause competes.
///
/// Executed by [`the_protein_axis_must_split_a_separation_one_delins`], which is
/// **red on purpose**.
const FERRO_IS_WRONG: &[Worked] = &[
    Worked {
        id: "W48",
        clause: "recommendations/protein/delins.md:64",
        quote: "the variant is not described as `p.Ser44_Trp46delinsArgLeuArg`.",
        input: "NP_SPECW48.1:p.Ser44_Trp46delinsArgLeuArg",
        published: "p.[Ser44Arg;Trp46Arg]",
        answer: "NP_SPECW48.1:p.[Ser44Arg;Trp46Arg]",
        rebase: Rebase::Accession,
        fixture: Fixture::Protein,
        why: "the spec's single clearest split-side worked example",
    },
    Worked {
        id: "W59",
        clause: "recommendations/protein/delins.md:21",
        quote: "two variants separated by one or more amino acids should be described individually and not as a \"delins\".",
        input: "NP_SPECW48.1:p.Ser988_Gln990delinsAlaLeuHis",
        published: "p.[Ser988Ala;Gln990His]",
        answer: "NP_SPECW48.1:p.[Ser988Ala;Gln990His]",
        rebase: Rebase::Accession,
        fixture: Fixture::Protein,
        why: "the same shape at a second locus, which is what makes it systematic \
              rather than a one-off",
    },
];

// ---------------------------------------------------------------------------
// The rows this file does not execute
// ---------------------------------------------------------------------------

/// Already executed, in all three shipped error modes, by
/// `tests/it/spec_worked_examples.rs`. Deliberately not duplicated here.
const CROSS_REFERENCED: &[(&str, &str)] = &[
    (
        "W33",
        "spec_worked_examples: NM_024312.4:c.2686A[10] -> c.2692_2693dup",
    ),
    (
        "W34",
        "spec_worked_examples: NM_024312.4:c.1738TA[6] -> c.1741_1742insTATATATA",
    ),
    (
        "W35",
        "spec_worked_examples: NM_024312.4:c.-6_-3G[6] is valid in the 5'UTR",
    ),
    (
        "W36",
        "spec_worked_examples: the RNA analogues, r.2686a[10] / r.1738ua[6] / r.-6_-3g[6]",
    ),
    (
        "W65",
        "spec_worked_examples: LRG_199t1:c.3921dup, and the c.3922dup far-side pair",
    ),
];

/// Worked examples that cannot be executed hermetically, each with its reason.
///
/// Pinned by [`the_unexecutable_census_is_pinned`]. Shrinking this list is the
/// point; growing it silently is the failure mode this table exists to prevent.
const UNEXECUTABLE: &[Unexecutable] = &[
    // -- needs a derived protein consequence -------------------------------
    Unexecutable { id: "W17", reason: "protein consequence of a DNA edit: needs translation of NM_004006-like CDS through a real transcript/protein pair" },
    Unexecutable { id: "W18", reason: "frameshift consequence (p.Glu5ValfsTer5) derived from c.6_13dup: needs translation" },
    Unexecutable { id: "W19", reason: "frameshift first-new-residue rule (p.Gln151Thrfs*9): needs the shifted reading frame" },
    Unexecutable { id: "W20", reason: "dup inside the stop codon: the spec gives three different answers for this geometry (catalog conflict C6), so there is no single published answer to assert" },
    Unexecutable { id: "W21", reason: "an upstream ATG created by c.-23A>T: needs 5'UTR translation" },
    Unexecutable { id: "W22", reason: "start-loss activating an upstream initiation site: the prose says delins and the worked answer is an insertion (catalog conflict C7)" },
    Unexecutable { id: "W23", reason: "protein consequence of NM_004380.2:c.139delinsTCATCATGAGCTG: needs translation" },
    Unexecutable { id: "W24", reason: "protein consequence of NM_004380.2:c.138_139ins...: needs translation" },
    Unexecutable { id: "W25", reason: "protein consequence of c.9_10insGGGTAG: needs translation" },
    Unexecutable { id: "W26", reason: "protein consequence of NM_080877.2:c.1733_1735delinsTTT: needs translation" },
    Unexecutable { id: "W27", reason: "protein consequence of NM_000222.3:c.1676_1684del: needs translation" },
    Unexecutable { id: "W28", reason: "protein consequence of an inversion (NM_003079.4:c.374_395inv): needs translation" },
    Unexecutable { id: "W37", reason: "the FMR1 mixed GGC/GGA repeat: needs the real NM_002024.5 5'UTR bases" },
    Unexecutable { id: "W39", reason: "keyed on population variability ('when the repeat is variable in the population'), which no reference sequence carries" },
    Unexecutable { id: "W40", reason: "same, on the RNA axis" },
    Unexecutable { id: "W41", reason: "an L1 insertion whose payload cites an external accession and a target-site duplication; no hermetic stand-in for the inserted element" },
    Unexecutable { id: "W54", reason: "two competing partitions of one evidence set joined by `^`, with the protein consequence of c.472A>G/c.473A>T: needs translation" },
    Unexecutable { id: "W62", reason: "the g./c. pair for a minus-strand transcript: needs a real genome plus cdot alignment, not a transcript slice" },
    Unexecutable { id: "W63", reason: "same, at NM_004006.1/NC_000023.10 coordinates" },
    Unexecutable { id: "W66", reason: "the exon/intron border case (c.1704+1del): the committed slice serves no intronic bases at that junction, so a pass would be vacuous" },
    Unexecutable { id: "W67", reason: "the intron/exon border case (c.1813del): same" },
    Unexecutable { id: "W68", reason: "LRG_2t1:r.1034_1036del: that accession is not in any committed fixture" },
    Unexecutable { id: "W69", reason: "the CFTR deltaF508 pair: needs real NM_000492.3 bases" },
    Unexecutable { id: "W70", reason: "the same variant's full DNA/RNA/protein stack: needs translation" },
    Unexecutable { id: "W71", reason: "LRG_232t1:c.1365_1373del projected to protein: needs translation" },
    Unexecutable { id: "W72", reason: "a DMD exon-5 deletion projected to protein: needs translation and real exon structure" },
    Unexecutable { id: "W73", reason: "the ATXN7 CAG->AGC unit rotation: needs real NM_000333.3 bases" },
    Unexecutable { id: "W74", reason: "the FMR1 CGG->GGC rotation: needs real bases" },
    Unexecutable { id: "W75", reason: "the HTT CAG->AGC rotation: needs real bases" },
    Unexecutable { id: "W76", reason: "the CFTR intron-9 GT/T tract: needs real intronic bases" },
    Unexecutable { id: "W77", reason: "the DAB1 repeat spelled differently on g. and c.: needs a real minus-strand genome/transcript pair" },
    Unexecutable { id: "W78", reason: "origin-spanning del/dup on a circular molecule: needs full-length NC_012920.1 / J01749.1 contigs, and rebasing the coordinates onto a short contig would fabricate the published answer" },
    Unexecutable { id: "W79", reason: "the pre-SVD-WG006 interim spelling m.[1del;16569del]: a superseded suggestion, not a published answer" },
    Unexecutable { id: "W80", reason: "the 3' rule applied to an MLPA probe's ligation centre: an assay artefact, not a variant description" },
    Unexecutable { id: "W81", reason: "the published answer is 'use genomic coordinates', not a string" },
    Unexecutable { id: "W82", reason: "a numbering scheme (c.-15+1 ... c.-14-1), with no input/output pair" },
    Unexecutable { id: "W83", reason: "same, for an intron immediately after the stop codon" },
    Unexecutable { id: "W84", reason: "a pericentric inversion plus an abutting substitution at megabase coordinates; a short synthetic contig cannot carry it and rebasing would fabricate the answer" },
    Unexecutable { id: "W85", reason: "same, with a 275 bp deletion and an anchorless insG" },
    Unexecutable { id: "W86", reason: "a cross-axis g./c./r./p. quadruple: needs projection through a real reference" },
    Unexecutable { id: "W87", reason: "an RNA consequence (r.(3277_3432del)) predicted from a coding substitution: needs splicing prediction" },
    Unexecutable { id: "W88", reason: "alternative RNA products plus their protein consequences: needs projection and translation" },
    Unexecutable { id: "W89", reason: "same shape at c.93G>T" },
    Unexecutable { id: "W90", reason: "an RNA insertion whose payload cites intronic c. coordinates: needs a genomic reference" },
    Unexecutable { id: "W91", reason: "same shape for c.831+2T>A, whose RNA answer inserts intronic sequence" },
    Unexecutable { id: "W92", reason: "SVD-WG011 is an OPEN proposal, not in force" },
    Unexecutable { id: "W93", reason: "SVD-WG011 again, for c.1704+2T>A: an OPEN proposal, not in force" },
    Unexecutable { id: "W94", reason: "SVD-WG011 again, for c.8669-24_8669-19del: an OPEN proposal, not in force" },
    Unexecutable { id: "W95", reason: "a predicted silent protein consequence of LRG_199t1:c.123C>T: needs translation" },
];

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

fn repo_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

/// The committed spec worked-example reference window, loaded once.
fn slice_provider() -> &'static WindowProvider {
    static SLICE: OnceLock<WindowProvider> = OnceLock::new();
    SLICE.get_or_init(|| {
        WindowFixture::from_json_path(&repo_path(WINDOWS_PATH))
            .expect("load the committed spec worked-example reference slice")
            .to_provider()
    })
}

fn transcript(id: &str, sequence: &str, cds: Option<(u64, u64)>) -> Transcript {
    let (cds_start, cds_end) = match cds {
        Some((start, end)) => (Some(start), Some(end)),
        None => (None, None),
    };
    Transcript::new(
        id.to_string(),
        Some("SPEC".to_string()),
        Strand::Plus,
        sequence.to_string(),
        cds_start,
        cds_end,
        vec![Exon::new(1, 1, sequence.len() as u64)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    )
}

/// Cyclic `ACGT` filler with literal blocks planted at 1-based positions.
///
/// Used where a worked example names a coordinate far from any sequence it
/// prints (`c.2077`, `r.142`): the printed context goes at the printed
/// coordinate and the filler only has to be *something*. The filler is not a
/// homopolymer or a tandem repeat, so it cannot create a shift the example does
/// not intend.
fn planted(length: usize, blocks: &[(usize, &str)]) -> String {
    let mut bases: Vec<u8> = (0..length).map(|i| b"ACGT"[(i * 7 + i / 4) % 4]).collect();
    for (start, block) in blocks {
        for (offset, base) in block.bytes().enumerate() {
            bases[start - 1 + offset] = base;
        }
    }
    String::from_utf8(bases).expect("planted sequences are ASCII")
}

/// A protein of `length` residues, `Ala` everywhere except the planted ones.
fn protein(length: usize, residues: &[(usize, u8)]) -> String {
    let mut sequence = vec![b'A'; length];
    for (position, residue) in residues {
        sequence[position - 1] = *residue;
    }
    String::from_utf8(sequence).expect("protein sequences are ASCII")
}

/// Every synthetic fixture, built once.
///
/// Each sequence is the one the spec prints beside its example, under an
/// accession that cannot be mistaken for a real one.
fn synthetic_provider() -> &'static MockProvider {
    static SYNTHETIC: OnceLock<MockProvider> = OnceLock::new();
    SYNTHETIC.get_or_init(|| {
        let mut provider = MockProvider::new();

        // W2  — `TGT` + `GC` + `CA` (DNA/delins.md:77): g.4_5 is the `GC`.
        provider.add_genomic_sequence("SPEC_W02.1", "TGTGCCA");
        // W3  — `..GAA` + `GC` + `CAG..` (DNA/substitution.md:90): g.12_13 is
        //       the `GC`, g.9_11 the `GAA`, g.14_16 the `CAG`.
        provider.add_genomic_sequence("SPEC_W03.1", "TTAGCTTAGAAGCCAGTTAC");
        // W13/W14 — `AAGC` at g.1_4 (DNA/inversion.md:53-62).
        provider.add_genomic_sequence("SPEC_W13.1", "AAGCTTGACCTG");
        // W15 — `..CA` + `TCAG` + `CCT..` (DNA/inversion.md:28): g.3_6 inverts.
        provider.add_genomic_sequence("SPEC_W15.1", "CATCAGCCT");
        // W60 — `ATT` + `T` + `G` (glossary.md:10): the T-run is g.2_4.
        provider.add_genomic_sequence("SPEC_W60.1", "ATTTG");
        // W57 — `AG` + `C` + `GTTTAGC` (SVD-WG010.md:16): the T-run is g.5_7.
        provider.add_genomic_sequence("SPEC_W57.1", "AGCGTTTAGC");
        // W11 — `ATCG` + `ATCGATCGATCG` + `A` + `GGGTCCC` (DNA/insertion.md:109).
        provider.add_genomic_sequence("SPEC_W11.1", "ATCGATCGATCGATCGAGGGTCCC");
        // W42 — a plain contig long enough to carry the published coordinates.
        provider.add_genomic_sequence("SPEC_W42.1", planted(300, &[]));

        // W4  — `augu` + `gc` + `ca` (RNA/delins.md:66), served as DNA bases and
        //       rendered on the r. axis.
        provider.add_transcript(transcript("NM_SPECW04.1", "ATGTGCCAGGGGCCCAAA", None));
        // W10 — `acuuacu` + `u` + `gcc` (RNA/duplication.md:29).
        provider.add_transcript(transcript("NM_SPECW10.1", "ACTTACTGCCAGGGTTT", None));
        // W16 — `cagc` + `ATGGAGCC` + `GGCGGCGGGGAGCAGC` + `ATGGAGCC` + `TTCG`
        //       (DNA/insertion.md:117), CDS starting at the first `ATG`.
        provider.add_transcript(transcript(
            "NM_SPECW16.1",
            "CAGCATGGAGCCGGCGGCGGGGAGCAGCATGGAGCCTTCGAAAAAAAAAAAAAAAAAAAA",
            Some((5, 58)),
        ));
        // W43/W44 — BRCA1 `c.2074`–`c.2080` is `..CAT` + `G` + `ACA..`
        //            (DNA/delins.md:86, RNA/delins.md:68).
        provider.add_transcript(transcript(
            "NM_SPECW43.1",
            &planted(2400, &[(2074, "CATGACA")]),
            Some((1, 2400)),
        ));
        // W47 — `r.142`–`r.144` is `cga`, one codon (RNA/delins.md:18); with the
        //       CDS starting at position 1, r.142 opens codon 48.
        provider.add_transcript(transcript(
            "NM_SPECW47.1",
            &planted(600, &[(142, "CGA")]),
            Some((1, 600)),
        ));

        // W5  — `TrpVal` at p.24_25 (protein/substitution.md:98).
        provider.add_protein("NP_SPECW05.1", protein(40, &[(24, b'W'), (25, b'V')]));
        // W48/W55/W59 — `Ser44`, `Leu45`, `Trp46` and `Ser988`, `Leu989`,
        //               `Gln990` (protein/delins.md:62-64, SVD-WG010.md:53-54),
        //               plus `Phe233`/`Cys690` for W55.
        provider.add_protein(
            "NP_SPECW48.1",
            protein(
                1200,
                &[
                    (44, b'S'),
                    (45, b'L'),
                    (46, b'W'),
                    (233, b'F'),
                    (690, b'C'),
                    (988, b'S'),
                    (989, b'L'),
                    (990, b'Q'),
                ],
            ),
        );
        // W29/W30 — `MetTrpSerSerSerHisAsp..` (general.md:159,
        //           protein/duplication.md:97).
        provider.add_protein("NP_SPECW29.1", "MWSSSHDKLMNPQRSTVWY");
        // W31/W32 — `MetLeuTrpTrpGlu` (protein/deletion.md:38).
        provider.add_protein("NP_SPECW31.1", "MLWWEKLMNPQRSTVWY");

        provider
    })
}

// ---------------------------------------------------------------------------
// Running a row
// ---------------------------------------------------------------------------

fn normalize_with<P: ReferenceProvider + Clone>(
    provider: &P,
    input: &str,
) -> Result<String, String> {
    let variant = parse_hgvs(input).map_err(|e| format!("parse error: {e}"))?;
    Normalizer::new(provider.clone())
        .normalize(&variant)
        .map(|n| n.to_string())
        .map_err(|e| format!("normalize error: {e}"))
}

/// Normalize `input` through the fixture `fixture` names.
fn run(fixture: Fixture, input: &str) -> Result<String, String> {
    match fixture {
        Fixture::Slice => normalize_with(slice_provider(), input),
        _ => normalize_with(synthetic_provider(), input),
    }
}

/// Normalize and re-parse, returning the parsed normal form.
fn normalized_variant(fixture: Fixture, input: &str) -> HgvsVariant {
    let rendered = run(fixture, input).unwrap_or_else(|e| panic!("{input}: {e}"));
    parse_hgvs(&rendered).unwrap_or_else(|e| {
        panic!("{input} normalized to {rendered}, which does not re-parse: {e}")
    })
}

/// The `(start, end)` of one member, on that member's own axis.
fn span_of(variant: &HgvsVariant) -> Option<(i64, i64)> {
    match variant {
        HgvsVariant::Genome(v) => {
            let location = &v.loc_edit.location;
            Some((
                location.start.inner()?.base as i64,
                location.end.inner()?.base as i64,
            ))
        }
        HgvsVariant::Mt(v) => {
            let location = &v.loc_edit.location;
            Some((
                location.start.inner()?.base as i64,
                location.end.inner()?.base as i64,
            ))
        }
        HgvsVariant::Cds(v) => {
            let location = &v.loc_edit.location;
            Some((location.start.inner()?.base, location.end.inner()?.base))
        }
        HgvsVariant::Tx(v) => {
            let location = &v.loc_edit.location;
            Some((location.start.inner()?.base, location.end.inner()?.base))
        }
        HgvsVariant::Rna(v) => {
            let location = &v.loc_edit.location;
            Some((location.start.inner()?.base, location.end.inner()?.base))
        }
        HgvsVariant::Protein(v) => {
            let location = &v.loc_edit.location;
            Some((
                location.start.inner()?.number as i64,
                location.end.inner()?.number as i64,
            ))
        }
        _ => None,
    }
}

/// Every member's span, in rendered order. A bare (single-member) description
/// yields one span.
fn member_spans(variant: &HgvsVariant) -> Vec<(i64, i64)> {
    match variant {
        HgvsVariant::Allele(allele) => allele.variants.iter().flat_map(member_spans).collect(),
        other => span_of(other).into_iter().collect(),
    }
}

/// Everything before the first `:`, i.e. the accession, or `""` when there is
/// none (the spec often prints a bare `c.`/`p.` fragment).
fn without_accession(description: &str) -> &str {
    match description.split_once(':') {
        Some((_, rest)) => rest,
        None => description,
    }
}

/// The trailing edit token of a description: `del`, `dup`, `inv`, `delinsTG`, …
fn edit_token(description: &str) -> String {
    description
        .chars()
        .rev()
        .take_while(|c| c.is_ascii_alphabetic())
        .collect::<Vec<_>>()
        .into_iter()
        .rev()
        .collect()
}

// ---------------------------------------------------------------------------
// Tests — the citations
// ---------------------------------------------------------------------------

/// Every row's quote is still on the line it cites.
///
/// A pinned answer with no authority behind it is a change detector, not a
/// record. Verifying the quote byte-for-byte is what keeps the citation honest
/// across a spec-submodule bump: a bare line number resolves against any file
/// long enough to have that line.
#[test]
fn every_row_quotes_the_spec_verbatim() {
    let docs = repo_path(SPEC_DOCS);
    assert!(
        docs.join("recommendations/general.md").is_file(),
        "the vendored HGVS spec checkout at {SPEC_DOCS} is empty. Initialise it:\n    \
         git submodule update --init assets/hgvs-nomenclature"
    );

    let mut citations: Vec<(&str, &str, &str)> = Vec::new();
    for row in WORKED.iter().chain(FERRO_IS_WRONG) {
        citations.push((row.id, row.clause, row.quote));
    }
    for row in NEGATIVE {
        citations.push((row.id, row.clause, row.quote));
    }
    for row in REJECTED {
        citations.push((row.id, row.clause, row.quote));
    }
    for row in PARTITIONS {
        citations.push((row.id, row.clause, row.quote));
    }
    for row in DIVERGENT {
        citations.push((row.id, row.clause, row.quote));
    }
    for row in PARSE_GAPS {
        citations.push((row.id, row.clause, row.quote));
    }

    let mut broken = Vec::new();
    for (id, clause, quote) in &citations {
        let (path, line_number) = clause
            .rsplit_once(':')
            .unwrap_or_else(|| panic!("{id}: citation {clause} is not `path:line`"));
        let line_number: usize = line_number
            .parse()
            .unwrap_or_else(|_| panic!("{id}: citation {clause} has no line number"));
        let text = match std::fs::read_to_string(docs.join(path)) {
            Ok(text) => text,
            Err(e) => {
                broken.push(format!("{id}: cannot read {path}: {e}"));
                continue;
            }
        };
        match text.lines().nth(line_number - 1) {
            None => broken.push(format!("{id}: {path} has no line {line_number}")),
            Some(line) if !line.contains(quote) => broken.push(format!(
                "{id}: {clause} no longer contains its quote\n    quote: {quote}\n    line:  {}",
                line.trim()
            )),
            Some(_) => {}
        }
    }
    assert!(
        broken.is_empty(),
        "{} citation(s) no longer quote the spec:\n{}",
        broken.len(),
        broken.join("\n")
    );

    assert!(
        citations.len() >= 60,
        "the citation set shrank to {} — rows are being dropped rather than moved",
        citations.len()
    );
}

/// Every row explains itself.
#[test]
fn every_row_records_why_it_is_here() {
    for row in WORKED.iter().chain(FERRO_IS_WRONG) {
        assert!(!row.why.trim().is_empty(), "{} records no reason", row.id);
    }
    for row in UNEXECUTABLE {
        assert!(
            row.reason.len() > 20,
            "{} gives no usable reason for being unexecutable",
            row.id
        );
    }
}

// ---------------------------------------------------------------------------
// Tests — the fixtures are what the spec says they are
// ---------------------------------------------------------------------------

/// The committed slice really carries the bases the spec quotes.
///
/// Every `LRG_199t1` row below is only as good as this: if the slice's bases
/// were not the ones the spec prints beside its example, a row could pass while
/// describing a different variant. So the spec's own quoted contexts are
/// re-derived from the fixture before anything else runs.
#[test]
fn the_slice_carries_the_bases_the_spec_quotes() {
    let fixture = WindowFixture::from_json_path(&repo_path(WINDOWS_PATH))
        .expect("load the committed spec worked-example reference slice");
    let tx = fixture
        .transcripts
        .iter()
        .find(|t| t.id == "LRG_199t1")
        .expect("LRG_199t1 in the committed slice");
    let sequence = tx.sequence.as_deref().expect("LRG_199t1 bases");
    let cds_start = tx.cds_start.expect("LRG_199t1 CDS start") as usize;
    // `c.n` is transcript position `cds_start + n - 1`, 1-based.
    let context =
        |from: usize, to: usize| -> &str { &sequence[cds_start + from - 2..cds_start + to - 1] };

    // W1: c.79_80 is the `GC` replaced by `TT`.
    assert_eq!(context(79, 80), "GC", "W1");
    // W45: c.145_147 is `CGC` (delins.md:38 names it), replaced by `TGG`.
    assert_eq!(context(145, 147), "CGC", "W45");
    // W46/W56: c.235 is `A`, c.237 and c.238 are `G`.
    assert_eq!(context(235, 238), "AAGG", "W46/W56");
    // W8/W9: `AGAAG` + `TAGA` + `GG` around c.20.
    assert_eq!(context(15, 25), "AGAAGTAGAGG", "W8/W9");
    // W61: the 8-nt A-stretch from c.5690 to c.5697, inside `ATTGAAAAAAAA` + `TTA`.
    assert_eq!(context(5686, 5700), "ATTGAAAAAAAATTA", "W61");
    // W7: `..GGAA` + `GAG` + `TTGC..` at c.6775_6777.
    assert_eq!(context(6771, 6780), "GGAAGAGTTG", "W7");
    // W38: c.123 is a `C` (other.md:27).
    assert_eq!(context(123, 123), "C", "W38");
    // W51/W52: c.2376 is a `G`; W59: c.2962 is `T` and c.2970 is `A`.
    assert_eq!(context(2376, 2376), "G", "W51/W52");
    assert_eq!(context(2962, 2962), "T", "W59");
    assert_eq!(context(2970, 2970), "A", "W59");
    // W64: the exon boundary at c.3921, whose neighbour c.3922 carries the same
    // base — which is why the junction exception exists at all.
    assert_eq!(context(3920, 3922), "ATT", "W64");
}

/// The synthetic fixtures carry the sequence their worked example prints.
#[test]
fn the_synthetic_fixtures_carry_the_sequence_the_spec_prints() {
    let provider = synthetic_provider();
    let genomic = |accession: &str, from: u64, to: u64| -> String {
        provider
            .get_sequence(accession, from - 1, to)
            .unwrap_or_else(|e| panic!("{accession}:{from}-{to}: {e}"))
            .to_uppercase()
    };
    // W2: `TGT` + `GC` + `CA`.
    assert_eq!(genomic("SPEC_W02.1", 1, 7), "TGTGCCA", "W2");
    // W3: `..GAA` + `GC` + `CAG..` with the `GC` at g.12_13.
    assert_eq!(genomic("SPEC_W03.1", 9, 16), "GAAGCCAG", "W3");
    // W13/W14: `AAGC`.
    assert_eq!(genomic("SPEC_W13.1", 1, 4), "AAGC", "W13/W14");
    // W15: `..CA` + `TCAG` + `CCT..`; note `TCAG` reverse-complements to `CTGA`.
    assert_eq!(genomic("SPEC_W15.1", 1, 9), "CATCAGCCT", "W15");
    // W60: `ATT` + `T` + `G`.
    assert_eq!(genomic("SPEC_W60.1", 1, 5), "ATTTG", "W60");
    // W57: `AG` + `C` + `GTTTAGC`.
    assert_eq!(genomic("SPEC_W57.1", 1, 10), "AGCGTTTAGC", "W57");
    // W11: the 12-mer at g.5_16 is repeated at g.1_4 and is *not* immediately 5'
    // of the insertion point g.17_18 — the whole point of the example.
    assert_eq!(genomic("SPEC_W11.1", 5, 16), "ATCGATCGATCG", "W11");
    assert_eq!(genomic("SPEC_W11.1", 17, 17), "A", "W11");
}

// ---------------------------------------------------------------------------
// Tests — the answers
// ---------------------------------------------------------------------------

/// Every worked example produces the spec's published answer, exactly.
///
/// Reported as one batch rather than row by row: when a normalizer change moves
/// these, *which* rows moved together is the diagnosis.
#[test]
fn every_worked_example_produces_the_spec_published_answer() {
    let mut failures = Vec::new();
    for row in WORKED {
        match run(row.fixture, row.input) {
            Ok(actual) if actual == row.answer => {}
            Ok(actual) => failures.push(format!(
                "  [{}] {}\n    expected: {}\n    actual:   {}\n    {} — {}",
                row.id, row.input, row.answer, actual, row.clause, row.why
            )),
            Err(e) => failures.push(format!(
                "  [{}] {}\n    expected: {}\n    error:    {}\n    {}",
                row.id, row.input, row.answer, e, row.clause
            )),
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} worked examples no longer produce the spec's answer:\n{}",
        failures.len(),
        WORKED.len(),
        failures.join("\n")
    );
}

/// A rebased answer is the spec's string with only the accession changed.
///
/// The hermetic fixtures force some rows off the spec's literal accession. This
/// is what stops that concession from widening into "close enough": for a
/// [`Rebase::Accession`] row the two strings must be identical after the
/// accession, and for a [`Rebase::Local`] row — where the spec numbers its
/// example on a real chromosome — the edit must at least be the same edit, with
/// the position pinned separately by [`PARTITIONS`].
#[test]
fn every_rebased_answer_differs_from_the_published_string_only_in_its_accession() {
    for row in WORKED.iter().chain(FERRO_IS_WRONG) {
        match row.rebase {
            Rebase::Same => assert_eq!(
                row.answer, row.published,
                "[{}] claims Rebase::Same but the strings differ",
                row.id
            ),
            Rebase::Accession => {
                assert_ne!(
                    row.answer, row.published,
                    "[{}] claims Rebase::Accession but nothing was rebased",
                    row.id
                );
                assert_eq!(
                    without_accession(row.answer),
                    without_accession(row.published),
                    "[{}] rebased more than the accession",
                    row.id
                );
            }
            Rebase::Local => {
                let (answer, published) = (edit_token(row.answer), edit_token(row.published));
                assert!(
                    !answer.is_empty() && published.ends_with(&answer),
                    "[{}] is Rebase::Local but its edit token ({answer}) is not the \
                     spec's ({published}); a local rebase may move coordinates, never the edit",
                    row.id
                );
            }
            Rebase::Prose => {
                assert!(
                    parse_hgvs(row.published).is_err(),
                    "[{}] claims the spec states its answer as prose, but {} is a \
                     description — compare against it directly",
                    row.id,
                    row.published
                );
                assert!(
                    PARTITIONS.iter().any(|p| p.id == row.id),
                    "[{}] is Rebase::Prose, so its position is the whole content of \
                     the row and must be pinned in PARTITIONS",
                    row.id
                );
            }
        }
    }
}

/// Every answer is its own normal form.
///
/// Cheap, and worth having: each output is itself a legal description, so
/// re-normalizing it must be a no-op. A rule that reaches the right answer by
/// one shift too many and one back satisfies the test above and fails this.
#[test]
fn every_answer_is_a_fixed_point() {
    let mut inputs: Vec<(&str, Fixture, &str)> = Vec::new();
    for row in WORKED {
        inputs.push((row.id, row.fixture, row.input));
    }
    for row in NEGATIVE {
        inputs.push((row.id, row.fixture, row.input));
    }
    for row in PARTITIONS {
        inputs.push((row.id, row.fixture, row.input));
    }
    for row in DIVERGENT {
        inputs.push((row.id, row.fixture, row.input));
    }
    for (id, fixture, input) in inputs {
        let Ok(once) = run(fixture, input) else {
            continue;
        };
        let twice = run(fixture, &once);
        assert_eq!(
            twice.as_deref(),
            Ok(once.as_str()),
            "[{id}] {input} normalized to {once}, which is not a fixed point"
        );
    }
}

// ---------------------------------------------------------------------------
// Tests — the forbidden half
// ---------------------------------------------------------------------------

/// The description the spec forbids is never the one ferro emits.
#[test]
fn the_forbidden_description_is_never_what_ferro_emits() {
    let mut violations = Vec::new();
    for row in NEGATIVE {
        match run(row.fixture, row.input) {
            Ok(actual) if actual == row.forbidden => violations.push(format!(
                "  [{}] {} -> {}\n    forbidden by {} — {}",
                row.id, row.input, actual, row.clause, row.why
            )),
            Ok(_) => {}
            Err(e) => violations.push(format!(
                "  [{}] {} failed to normalize: {e}",
                row.id, row.input
            )),
        }
    }
    assert!(
        violations.is_empty(),
        "{} forbidden description(s) were produced:\n{}",
        violations.len(),
        violations.join("\n")
    );
}

/// The spellings the spec bans are refused by the parser.
#[test]
fn the_banned_spellings_are_refused_by_the_parser() {
    let mut accepted = Vec::new();
    for row in REJECTED {
        if let Ok(variant) = parse_hgvs(row.input) {
            accepted.push(format!(
                "  [{}] {} parsed as {variant}\n    banned by {} — {}",
                row.id, row.input, row.clause, row.why
            ));
        }
    }
    assert!(
        accepted.is_empty(),
        "{} banned spelling(s) are accepted by the parser:\n{}",
        accepted.len(),
        accepted.join("\n")
    );
}

// ---------------------------------------------------------------------------
// Tests — partitions
// ---------------------------------------------------------------------------

/// Each description's normal form asserts the partition the spec asserts.
///
/// Member count first, then every member's span. The string-level tests above
/// can be satisfied by a right answer reached over a wrong partition — most
/// obviously where a member's rendered coordinates coincide with a different
/// member's, which is exactly the shape `general.md:34` and `delins.md:16`
/// legislate.
#[test]
fn the_normal_forms_assert_the_partitions_the_spec_asserts() {
    let mut failures = Vec::new();
    for row in PARTITIONS {
        let normalized = normalized_variant(row.fixture, row.input);
        let spans = member_spans(&normalized);
        if spans.len() != row.members.len() {
            failures.push(format!(
                "  [{}] {} -> {normalized}\n    expected {} member(s), found {}: {:?}\n    {} — {}",
                row.id,
                row.input,
                row.members.len(),
                spans.len(),
                spans,
                row.clause,
                row.why
            ));
            continue;
        }
        if spans != row.members {
            failures.push(format!(
                "  [{}] {} -> {normalized}\n    expected spans {:?}, found {:?}\n    {} — {}",
                row.id, row.input, row.members, spans, row.clause, row.why
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} partitions moved:\n{}",
        failures.len(),
        PARTITIONS.len(),
        failures.join("\n")
    );
}

/// W11: the copy is not directly 3'-flanking, so the answer is an insertion.
///
/// The rendered payload is a recorded divergence (see [`DIVERGENT`]), which
/// would let the *substance* of `duplication.md:17` go unasserted. It does not:
/// this pins that the member is an insertion at `g.17_18` and not a `dup`,
/// independently of how its payload is spelled.
#[test]
fn a_copy_that_is_not_directly_three_prime_flanking_is_an_insertion() {
    let normalized = normalized_variant(Fixture::Genomic, "SPEC_W11.1:g.17_18ins5_16");
    assert_eq!(
        member_spans(&normalized),
        vec![(17, 18)],
        "the insertion must stay anchored between g.17 and g.18"
    );
    let rendered = normalized.to_string();
    assert!(
        !rendered.contains("dup"),
        "`duplication.md:17`: a duplication may only be used when the additional \
         copy is directly 3'-flanking, and here `g.17` intervenes. Got: {rendered}"
    );
}

// ---------------------------------------------------------------------------
// Tests — recorded divergences
// ---------------------------------------------------------------------------

/// The recorded divergences still diverge, exactly as recorded.
///
/// Both directions are asserted. A change that makes one of these produce the
/// spec's answer fails here — which is the point: promoting a divergence to a
/// conformance is a decision, and it should cost a deliberate edit rather than
/// happening as a side effect.
#[test]
fn the_recorded_divergences_are_still_exactly_what_was_recorded() {
    let mut moved = Vec::new();
    for row in DIVERGENT {
        match run(row.fixture, row.input) {
            Ok(actual) if actual == row.ferro_answer => {}
            Ok(actual) if actual == row.spec_answer => moved.push(format!(
                "  [{}] {} now produces the spec's answer ({}). That may be right — \
                 promote it into WORKED and say so.\n    {} — {}",
                row.id, row.input, row.spec_answer, row.clause, row.why
            )),
            Ok(actual) => moved.push(format!(
                "  [{}] {}\n    pinned: {}\n    actual: {}\n    {} — {}",
                row.id, row.input, row.ferro_answer, actual, row.clause, row.why
            )),
            Err(e) => moved.push(format!("  [{}] {} errored: {e}", row.id, row.input)),
        }
    }
    assert!(
        moved.is_empty(),
        "{} recorded divergence(s) moved:\n{}",
        moved.len(),
        moved.join("\n")
    );
}

/// The parse gaps are still gaps.
///
/// Pinned as a failure so that closing one is loud. `PARSE_GAPS` is a defect
/// list, not an exemption list.
#[test]
fn the_recorded_parse_gaps_are_still_gaps() {
    for row in PARSE_GAPS {
        assert!(
            parse_hgvs(row.input).is_err(),
            "[{}] {} now parses. The gap is closed — move the row into WORKED and \
             assert its answer.\n    {} — {}",
            row.id,
            row.input,
            row.clause,
            row.why
        );
    }
}

// ---------------------------------------------------------------------------
// Tests — non-confluence the spec's own examples expose
// ---------------------------------------------------------------------------

/// `DNA/delins.md:44-47`'s two descriptions of one variant are two fixed points.
///
/// Stated as an observation, not a requirement. `:44` gives the spanning
/// `c.850_901delinsTTCCTCGATGCCTG`, `:46` names
/// `c.[850_869del;874_881del;887_897del;901_902insG]` as an "alternative
/// description" of the *same* variant, and `:47` recommends the delins. Ferro
/// preserves both, so the pair is non-confluent — the spec's own worked example
/// of the problem #1235 is about.
///
/// Converging them is a representation change with a downstream cost, and which
/// form wins is the undecided ruling `canonical-form-choice-when-both-legal`. So
/// this test asserts only that the situation is what the record says it is; if a
/// change *does* converge them, it fails and points at the ruling.
#[test]
fn the_spanning_delins_and_its_published_split_are_two_fixed_points() {
    let delins = run(Fixture::Slice, "LRG_199t1:c.850_901delinsTTCCTCGATGCCTG")
        .expect("the spanning delins normalizes");
    let split = run(
        Fixture::Slice,
        "LRG_199t1:c.[850_869del;874_881del;887_897del;901_902insG]",
    )
    .expect("the published alternative normalizes");
    assert_eq!(delins, "LRG_199t1:c.850_901delinsTTCCTCGATGCCTG");
    assert_eq!(
        split,
        "LRG_199t1:c.[850_869del;874_881del;887_897del;901_902insG]"
    );
    assert_ne!(
        delins, split,
        "the spec's own non-confluent pair converged. That may well be right — but \
         it is a representation change and an open question (ruling \
         `canonical-form-choice-when-both-legal`), so it must be decided rather \
         than land as a side effect."
    );
}

/// SVD-WG010's example is likewise two fixed points, for the same reason.
///
/// `SVD-WG010.md:16` prints one variant twice — `g.3_4delinsGGT` and
/// `g.[3C>G;7dup]` — and would have made the first canonical. The proposal was
/// rejected (`:5`), so `general.md:34` governs and the split is the conformant
/// form; ferro produces the split from the split and the delins from the delins.
///
/// This is the *sequence-identical* case: both spellings denote `AGGGTTTTAGC`.
/// It is recorded rather than asserted-away because converging them means
/// choosing which of two legal forms ships — the same open ruling as above.
#[test]
fn the_two_svd_wg010_spellings_of_one_variant_are_two_fixed_points() {
    let split = run(Fixture::Genomic, "SPEC_W57.1:g.[3C>G;7dup]").expect("the split normalizes");
    let delins = run(Fixture::Genomic, "SPEC_W57.1:g.3_4delinsGGT").expect("the delins normalizes");
    assert_eq!(split, "SPEC_W57.1:g.[3C>G;7dup]");
    assert_eq!(delins, "SPEC_W57.1:g.3_4delinsGGT");
    assert_ne!(
        split, delins,
        "the SVD-WG010 pair converged; see the ruling \
         `canonical-form-choice-when-both-legal` before accepting it"
    );
}

// ---------------------------------------------------------------------------
// Tests — the defect
// ---------------------------------------------------------------------------

/// **RED ON PURPOSE.** The protein axis does not split a separation-1 delins.
///
/// `protein/delins.md:21` — "two variants separated by one or more amino acids
/// should be described individually and not as a 'delins'" — and its worked
/// example `:62-64` publishes `p.[Ser44Arg;Trp46Arg]` and states, in the same
/// breath, that "the variant is not described as
/// `p.Ser44_Trp46delinsArgLeuArg`". The protein axis has **no codon exception**:
/// unlike `general.md:35` on the DNA/RNA axes, nothing licenses merging across
/// one unchanged residue. So the merged spelling is not a legal description and
/// must be split.
///
/// Ferro preserves it. Both directions were measured on a hermetic protein whose
/// residue 45 is `Leu`, so the two spellings denote one protein sequence:
///
/// ```text
///   p.[Ser44Arg;Trp46Arg]           -> p.[Ser44Arg;Trp46Arg]           (right)
///   p.Ser44_Trp46delinsArgLeuArg    -> p.Ser44_Trp46delinsArgLeuArg    (wrong)
/// ```
///
/// It is not a rendering quirk and not a one-off: the same shape at
/// `p.Ser988_Gln990delinsAlaLeuHis` (SVD-WG010.md:53-54, whose merge was
/// *rejected*) behaves identically. The adjacent case is handled correctly —
/// `p.[Trp24Cys;Val25Arg]` merges to `p.Trp24_Val25delinsCysArg` (W5) — so the
/// merge direction works on this axis and only the split direction is missing.
/// Under the partition model that split is one of the two licensed moves:
/// an equal-length member whose interior holds an unchanged run reaching it
/// (`general.md:34`, read backwards). The DNA axis performs it; the protein axis
/// does not.
///
/// **Do not make this test pass by pinning ferro's answer.** There is no
/// competing clause to cite, so this is not a `DIVERGENT` row. It is either
/// fixed in `src/` or it stays red.
#[test]
fn the_protein_axis_must_split_a_separation_one_delins() {
    let mut failures = Vec::new();
    for row in FERRO_IS_WRONG {
        match run(row.fixture, row.input) {
            Ok(actual) if actual == row.answer => {}
            Ok(actual) => failures.push(format!(
                "  [{}] {}\n    spec ({}): {}\n    ferro:   {}\n    {}",
                row.id, row.input, row.clause, row.answer, actual, row.why
            )),
            Err(e) => failures.push(format!("  [{}] {} errored: {e}", row.id, row.input)),
        }
    }
    assert!(
        failures.is_empty(),
        "{} protein worked example(s) are not what the spec publishes. This is a \
         known, unfixed defect — see this test's doc comment — not a regression, \
         and not something to pin:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

// ---------------------------------------------------------------------------
// Tests — the census
// ---------------------------------------------------------------------------

/// Every id named anywhere in this file, in one set.
fn covered_ids() -> BTreeSet<&'static str> {
    let mut ids = BTreeSet::new();
    for row in WORKED.iter().chain(FERRO_IS_WRONG) {
        ids.insert(row.id);
    }
    for row in NEGATIVE {
        ids.insert(row.id);
    }
    for row in REJECTED {
        ids.insert(row.id);
    }
    for row in PARTITIONS {
        ids.insert(row.id);
    }
    for row in DIVERGENT {
        ids.insert(row.id);
    }
    for row in PARSE_GAPS {
        ids.insert(row.id);
    }
    ids
}

fn unexecutable_ids() -> BTreeSet<&'static str> {
    UNEXECUTABLE.iter().map(|row| row.id).collect()
}

fn cross_referenced_ids() -> BTreeSet<&'static str> {
    CROSS_REFERENCED.iter().map(|(id, _)| *id).collect()
}

/// All 95 worked examples are accounted for, one way or another.
///
/// This is the test that makes the coverage claim checkable. Without it, the
/// honest way to raise the pass rate is to delete the hard rows — and nothing
/// would notice.
#[test]
fn the_roster_covers_every_one_of_the_95_worked_examples() {
    let executed = covered_ids();
    let unexecutable = unexecutable_ids();
    let cross_referenced = cross_referenced_ids();

    let overlap: Vec<&str> = executed.intersection(&unexecutable).copied().collect();
    assert!(
        overlap.is_empty(),
        "these rows are both executed and listed as unexecutable: {overlap:?}"
    );
    let duplicated: Vec<&str> = executed.intersection(&cross_referenced).copied().collect();
    assert!(
        duplicated.is_empty(),
        "these rows are already covered by tests/it/spec_worked_examples.rs and must \
         not be duplicated here: {duplicated:?}"
    );

    let mut all: BTreeSet<&str> = executed;
    all.extend(&unexecutable);
    all.extend(&cross_referenced);

    let expected: BTreeSet<String> = (1..=WORKED_EXAMPLE_COUNT)
        .map(|n| format!("W{n}"))
        .collect();
    let seen: BTreeSet<String> = all.iter().map(|id| (*id).to_string()).collect();
    let missing: Vec<&String> = expected.difference(&seen).collect();
    let unknown: Vec<&String> = seen.difference(&expected).collect();
    assert!(
        missing.is_empty() && unknown.is_empty(),
        "the roster does not account for all {WORKED_EXAMPLE_COUNT} worked examples.\n  \
         missing: {missing:?}\n  not in the catalog: {unknown:?}"
    );
}

/// The size of what is *not* covered is pinned, in both directions.
///
/// Growing [`UNEXECUTABLE`] is how coverage silently rots; shrinking it without
/// adding rows is how a hard row gets dropped. Either direction fails here.
#[test]
fn the_unexecutable_census_is_pinned() {
    assert_eq!(
        UNEXECUTABLE.len(),
        49,
        "the set of worked examples this file cannot execute changed size. Adding \
         one means claiming a row is not executable hermetically — an argument, not \
         a table edit. Removing one means it is now executed, so add it to a table."
    );
    assert_eq!(
        CROSS_REFERENCED.len(),
        5,
        "the set of rows delegated to tests/it/spec_worked_examples.rs changed"
    );
    assert_eq!(
        covered_ids().len(),
        WORKED_EXAMPLE_COUNT - UNEXECUTABLE.len() - CROSS_REFERENCED.len(),
        "the executed / unexecutable / delegated counts no longer partition the 95 rows"
    );
    assert_eq!(
        (DIVERGENT.len(), PARSE_GAPS.len(), FERRO_IS_WRONG.len()),
        (2, 1, 2),
        "the set of rows ferro does not satisfy changed. Each of these three tables \
         is an adjudication — a divergence with a competing clause, a parser defect, \
         or a plain defect — so a change here is a ruling, not a fixture edit."
    );
}
