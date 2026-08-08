//! Issue #1040: probes for `inv` over-recognition (reverse-complement window
//! carved out of a `delins`), pinning ferro's post-#1036 (v0.8.1) behavior.
//!
//! Background
//! ----------
//! #1040 reports that ferro splits a *contiguous* multi-nucleotide change to
//! expose an internal reverse-complement sub-window as an `inv`, where the
//! whole run should stay a single `delins`. That mechanism was the #1034 bug,
//! fixed by #1036: `rules::decompose_delins` now types a maximal contiguous
//! mismatch run as `inv` only when the *entire* run is a reverse complement
//! (`DNA/inversion.md`; `general.md:34`). A reverse-complement *sub-run* of a
//! longer contiguous change is no longer carved out.
//!
//! The literal anonymized example in #1040 is a single 10-nt contiguous
//! substitution run whose whole-run reverse complement does NOT equal the alt,
//! so on v0.8.1 it already normalizes to one `delins` — the spec-correct form
//! mutalyzer produces. This suite pins that, and probes the surrounding
//! boundary of the `inv` typer so any *future* over-recognition regression is
//! caught and so the true residual divergence surface is documented.
//!
//! Two `inv`-producing code paths exist, and both are probed here:
//!   * `decompose_delins` — length-preserving compound runs (the #1040 path).
//!   * `canonicalize_delins` step 4d — a single span; types `inv` only when the
//!     *shared-affix-trimmed* span is a whole reverse complement. Its trimmed
//!     ref and alt must be equal length, so a genuinely length-changing delins
//!     can never become an `inv`, and `shorten_inversion`'s outer-pair peeling
//!     is a provable no-op post-trim (a peelable outer pair would already have
//!     been trimmed as a shared affix). Length-changing revcomp probes are
//!     therefore expected to stay `delins`.
//!
//! Spec basis for the pinned outputs (`DNA/inversion.md`):
//!   * an inversion is "more than one nucleotide" that is the reverse
//!     complement — so a 2-nt reverse complement is a *valid* `inv`;
//!   * "two variants separated by one or more nucleotides should be described
//!     individually and not as a delins" — the delins exception is codon-only
//!     (two variants separated by one nucleotide affecting one amino acid), so
//!     on a genomic axis reverse-complement runs separated by an unchanged base
//!     are individual `inv`s, not one merged `delins`.
//!
//! Every assertion below is ferro's actual v0.8.1 output and is HGVS-spec
//! compliant. The remaining ferro-vs-mutalyzer `inv` divergence candidate is
//! the 2-nt inversion (P3/P4/P5): ferro emits `inv` per the spec's ">1 nt"
//! rule, whereas some tools render a 2-nt reverse complement as a `delins`.
//! That is a spec-vs-tool typing choice, not the #1040 sub-run carving, which
//! is gone.
//!
//! 3'-end code path
//! ----------------
//! The related #1041 investigation found that a variant within `window_size`
//! (100 bp) of a short contig's 3' end took a *different* code path on v0.8.1:
//! the fetch window overran EOF, the read errored, and the variant fell into
//! the minimal-notation fallback that runs neither 3' shifting nor the
//! `delins -> inv/sub/dup` canonicalization (fixed for #1041 by #1042's window
//! clamp). Because the real ferro-vs-mutalyzer comparison sites variants on
//! short `NG_` contigs, the near-3'-end path is the one that actually runs
//! there — so the last section pins the #1040 example on a short contig too.
//! Near the 3' end the fallback simply leaves a `delins` as a `delins`, which
//! is the correct answer for #1040 (the fallback cannot *carve* an `inv`); the
//! result is a single `delins` on both v0.8.1 and the #1042 branch.

use crate::common::synthetic::{assert_padded_preserving, normalize_to_string, SyntheticBuilder};
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};
use std::io::Write;

const G: &str = "NC_TEST.1";
const C: &str = "NM_TEST.1";

/// Build a genomic provider whose core sits at 1-based position 257 and
/// normalize `input` against it.
fn g(core: &str, input: &str) -> String {
    normalize_to_string(SyntheticBuilder::genomic(core).build(), input)
}

/// Build a CDS provider (whole core is the CDS) and normalize `input`.
fn c(core: &str, input: &str) -> String {
    let cds_end = core.len() as u64;
    normalize_to_string(
        SyntheticBuilder::cds(core, 1, cds_end, Strand::Plus).build(),
        input,
    )
}

// ---------------------------------------------------------------------------
// The #1040 report: a contiguous run must stay a single `delins` (FIXED).
// ---------------------------------------------------------------------------

#[test]
fn issue_1040_literal_ten_nt_contiguous_run_stays_delins() {
    // ref[257..266] = CAGTGACTAG (10 explicit substitutions, all mismatch).
    //   revcomp(CAGTGACTAG) = CTAGTCACTG != TGTCACGACT  => NOT a whole-run inv.
    //   Interior sub-window 259_262 (GTGA -> TCAC) IS a revcomp but must not be
    //   carved out of the contiguous run.
    // Spec-correct (and what mutalyzer produces): one delins.
    assert_eq!(
        g(
            "CAGTGACTAG",
            &format!(
                "{G}:g.[257C>T;258A>G;259G>T;260T>C;261G>A;262A>C;263C>G;264T>A;265A>C;266G>T]"
            ),
        ),
        format!("{G}:g.257_266delinsTGTCACGACT"),
    );
}

#[test]
fn issue_1040_same_run_as_delins_input_stays_delins() {
    // The identical change expressed as a delins rather than an allele of subs.
    assert_eq!(
        g("CAGTGACTAG", &format!("{G}:g.257_266delinsTGTCACGACT")),
        format!("{G}:g.257_266delinsTGTCACGACT"),
    );
}

#[test]
fn issue_1034_minimal_contiguous_run_stays_delins() {
    // ref CTG -> alt ACA. revcomp(CTG)=CAG != ACA; interior TG->CA is a revcomp
    // sub-run that must not be carved. The #1034 minimal case, both input forms.
    assert_eq!(
        g("CTG", &format!("{G}:g.257_259delinsACA")),
        format!("{G}:g.257_259delinsACA"),
    );
    assert_eq!(
        g("CTG", &format!("{G}:g.[257C>A;258T>C;259G>A]")),
        format!("{G}:g.257_259delinsACA"),
    );
}

// ---------------------------------------------------------------------------
// Guard: a genuine whole-run reverse complement STILL emits `inv`.
// (Under-recognizing these would be the opposite defect.)
// ---------------------------------------------------------------------------

#[test]
fn whole_run_reverse_complement_emits_inv() {
    // 3-nt whole-run revcomp: GCT -> AGC (= revcomp(GCT)).
    assert_eq!(
        g("GCT", &format!("{G}:g.257_259delinsAGC")),
        format!("{G}:g.257_259inv"),
    );
    // 6-nt whole-run revcomp: ACGTGC -> GCACGT (= revcomp(ACGTGC)).
    assert_eq!(
        g("ACGTGC", &format!("{G}:g.257_262delinsGCACGT")),
        format!("{G}:g.257_262inv"),
    );
}

// ---------------------------------------------------------------------------
// The sequence-first derivation reaches the same `inv`, from every spelling.
//
// `whole_run_reverse_complement_emits_inv` above pins the *per-member* answer:
// one `delins` member, typed by `rules::canonicalize_delins`. The three tests
// below pin what `merge::canonicalize_from_sequence` derives from `(reference,
// denoted sequence)` — the path a multi-member cis allele takes, and the one
// that used to refuse an inversion outright (the since-removed
// `needs_unsupported_form`).
//
// They use `assert_padded_preserving` rather than this file's `g` helper: it
// projects input and output through `hgvs_to_spdi` and compares the bases each
// denotes, so a re-typing that quietly changed the sequence fails here rather
// than merely disagreeing with a pinned string.
//
// Spec basis, and it is not prioritisation. `inversion.md:5` defines an
// inversion by the content of the whole span and says nothing about interior
// columns coinciding, so a 5 nt span whose replacement is its reverse
// complement qualifies. `delins.md:5` defines a delins as a replacement "which
// is not a substitution or inversion", so the *spanning* alternative is not
// merely lower-priority, it is excluded. `general.md:56`'s prioritised list is
// deliberately not cited: it ranks single-variant type labels for one span, it
// omits `delins`, and it puts substitution above inversion.
// ---------------------------------------------------------------------------

#[test]
fn a_derived_whole_block_inversion_renders_as_inv() {
    // GTTAA -> TTAAC is a 5 nt whole-span reverse complement whose 2nd and 4th
    // columns coincide, so a position-wise partition finds three substitutions
    // in it and the inversion disappears. Spelled as those three substitutions
    // the input is a genuine multi-member allele, so it is the derivation — not
    // the per-member typer — that has to put the inversion back together.
    assert_eq!(
        assert_padded_preserving("GTTAA", &format!("{G}:g.[257G>T;259T>A;261A>C]")),
        format!("{G}:g.257_261inv"),
    );
}

#[test]
fn every_spelling_of_a_derived_whole_block_inversion_converges_on_inv() {
    // The same 5 nt inversion written every way its changed columns allow: the
    // spanning description, the three substitutions, and each grouping of those
    // columns into `delins` runs. Two spellings of one allele must normalize to
    // one string (#1235).
    //
    // The three-substitution spelling is what pins the mechanism. It weighs 3
    // changed columns against the derivation's 5, so the changed-columns bound
    // in `canonicalize_from_sequence` refuses it whenever the widening to the
    // whole span happens *before* that bound — while the spanning spelling
    // weighs 5 and is accepted. One variant, two answers. That is why
    // `coalesce_whole_block_inversion` runs after the bound rather than inside
    // `partition_block`.
    //
    // `g.257_261delinsTTAAC` is included even though `delins.md:5` makes it an
    // incorrect spelling: ferro accepts it as input and must land on the same
    // normal form as the spellings that are correct.
    for input in [
        format!("{G}:g.257_261inv"),
        format!("{G}:g.257_261delinsTTAAC"),
        format!("{G}:g.[257G>T;259T>A;261A>C]"),
        format!("{G}:g.[257_259delinsTTA;261A>C]"),
        format!("{G}:g.[257G>T;259_261delinsAAC]"),
    ] {
        assert_eq!(
            assert_padded_preserving("GTTAA", &input),
            format!("{G}:g.257_261inv"),
            "spelling `{input}` did not converge on the inv",
        );
    }
    // ACGTGC -> GCACGT, the 6 nt case `whole_run_reverse_complement_emits_inv`
    // already pins from the `delins` side, now from its substitution spelling.
    assert_eq!(
        assert_padded_preserving("ACGTGC", &format!("{G}:g.[257A>G;259G>A;260T>C;262C>T]")),
        format!("{G}:g.257_262inv"),
    );
}

#[test]
fn a_whole_span_reverse_complement_is_not_merged_across_a_multi_base_separation() {
    // The control for the two tests above, and the reason they pin `inv` rather
    // than merely "one member". AAGCTA -> TAGCTT is *also* a whole-span reverse
    // complement, but only its first and last bases change and four unchanged
    // bases lie between them.
    //
    // `general.md:34` describes two variants separated by one or more
    // nucleotides individually; its letter names only `delins`, but its
    // rationale (`delins.md:81-84` — the variants may have been reported, or may
    // occur, individually) reaches any single spanning description. At four
    // unchanged bases these are two independent changes rather than interior
    // columns of one reverse-complement relation, so
    // `coalesce_whole_block_inversion`'s `every_separation_is_a_single_base`
    // gate must refuse them.
    //
    // Without this case the two tests above would pass just as well against a
    // rule that merged *any* whole-span reverse complement — the shape, not the
    // thing they name.
    let expected = format!("{G}:g.[257A>T;262A>T]");
    assert_eq!(
        assert_padded_preserving("AAGCTA", &format!("{G}:g.[257A>T;262A>T]")),
        expected,
    );
    assert_eq!(
        assert_padded_preserving("AAGCTA", &format!("{G}:g.257_262delinsTAGCTT")),
        expected,
    );
}

#[test]
fn a_dense_inversion_is_recognised_across_multi_base_separations() {
    // The counterpart of the control above, and the #1461 regression.
    //
    // `AACAACCGCGTG -> CACGCGGTTGTT` is a whole-span reverse complement whose
    // changed columns fall in three runs separated by **two** unchanged bases
    // each, so `every_separation_is_a_single_base` refuses it. Before #1461 that
    // was the end of it and the span came back as separate members.
    //
    // Separation cannot decide between this and the control above, and that is
    // the whole point of the pair. Measured on the real case #1461 reports —
    // `NC_000013.10:g.100809575_100810031inv`, 457 nt — the shipped rule emitted
    // **109** members whose widest gap is **5**, while the control's single gap
    // is **4**. The inversion that must be recognised has the *wider* gaps, so
    // every threshold on separation either refuses #1461 or admits the control.
    //
    // Density does separate them, which is `changed_columns_dominate_the_span`
    // and hence this test: 8 of 12 columns change here, against 2 of 6 in the
    // control.
    //
    // It separates them as a *length-correlated proxy*, though, and not because
    // an unchanged fraction says whether a reverse complement is structural.
    // Coincidences come in mirror pairs, so for a random block the unchanged
    // count is `2 * Binomial(floor(n/2), 1/4)` — mean `n/4` at every `n`.
    // Against that distribution the control's 66.7% is only +1.67 sd, which
    // about **one in six** genuine 6 nt inversions reaches; #1461's 32.8% is
    // +2.75 sd, which 0.45% reach. The separation comes from `n`, not from the
    // fraction. `changed_columns_dominate_the_span`'s doc works this through.
    let core = "AACAACCGCGTG";
    let expected = format!("{G}:g.257_268inv");
    for input in [
        format!("{G}:g.257_268inv"),
        format!("{G}:g.257_268delinsCACGCGGTTGTT"),
        format!("{G}:g.[257_257delinsC;260_265delinsGCGGTT;268G>T]"),
    ] {
        assert_eq!(
            assert_padded_preserving(core, &input),
            expected,
            "spelling `{input}` did not converge on the inv",
        );
    }
}

#[test]
fn a_derived_whole_block_inversion_outranks_the_codon_exception() {
    // The `c.`-axis half of the two tests above. It matters because every other
    // case in this section is genomic, where `general.md:35`'s codon exception
    // cannot fire at all — so this is the only coverage of what happens when the
    // whole-block inversion rule and the codon rule reach for the same pieces.
    //
    // `GTTAA -> TTAAC` is a whole-span reverse complement, and it changes c.1,
    // c.3 and c.5. c.1 and c.3 are two variants separated by one nucleotide
    // inside codon 1: precisely the `[Sub@p; Identity@p+1; Sub@p+2]` triplet the
    // codon exception merges into a `delins`. So both rules are live on the same
    // bases, which is what this row is for.
    //
    // **All four spellings give one answer, `c.1_5inv`.** A revision inside the
    // partition-model branch asserted two answers here — the loop was split in
    // two, with the divergence itself pinned — on the reasoning that the
    // one-member and three-member spellings assert different partitions. That is
    // true of the *partition* and false of the *type*: see the second half of
    // this test for why `coalesce_whole_block_inversion` re-types the hull, and
    // for the `2fa5bec0` evidence that `c.1_5inv` is the shipped baseline.

    // ---- The two rules, and which one wins.
    //
    // `general.md:34` has nothing to say about a LONE member, so the one-member
    // spellings keep the whole span as one piece; the equal-length split move
    // cannot cut them either, because it cuts only at an interior unchanged run
    // reaching the axis floor and the unchanged columns here (c.2, c.4) are
    // single bases against a floor of 2 on a reading-frame axis
    // (`general.md:35`).
    //
    // The MULTI-member spellings do get their members merged by the codon half of
    // `general.md:35`: c.1 and c.3 are both in codon 1, so they become
    // `1_3delinsTTA`, while c.3 and c.5 sit in codons 1 and 2 and so stay
    // individual — "two variants separated by one or more nucleotides should be
    // described individually and **not** as a \"delins\"".
    //
    // `coalesce_whole_block_inversion` then runs over whichever pieces resulted
    // and outranks the codon exception — the property this test is named for and
    // the one `merge::canonicalize_from_sequence` cites this test to pin.
    // Stubbing `apply_coding_codon_exception` instead leaves these rows green
    // (while reddening `issue_165_delins_sub_only_decompose`'s codon tests, so
    // the stub is doing something).
    //
    // **RE-ADJUDICATED 2026-08-08: all four spellings reach the `inv`, and the
    // split this row briefly asserted is withdrawn.** The paragraph above stops
    // one step short. The codon exception does merge c.1 with c.3 into
    // `1_3delinsTTA`, and c.3/c.5 do stay individual — but
    // `coalesce_whole_block_inversion` then runs, reads the hull of those two
    // pieces (they are one unchanged base apart, so its
    // `every_separation_is_a_single_base` arm admits them), finds `GTTAA` against
    // `TTAAC` to be a whole-span reverse complement, and types the hull `inv`.
    //
    // That ordering is not incidental, it is the ruling. `inversion-vs-two-delins-76-83`
    // is `decided` and governs `DNA/inversion.md:5` — "a sequence change where,
    // compared to a reference sequence, **more than one nucleotide** replacing the
    // original sequence is the reverse complement of the original sequence" —
    // which is a property of the *whole span* and not of any decomposition of it.
    // `DNA/delins.md:5` then forces the type: a delins is a replacement "**and
    // which is not** a substitution or inversion", so over this span `delins` is
    // not a description the spec admits at all. And the spec does not decompose
    // its own inversions: `DNA/inversion.md:33-34` publishes
    // `NM_004006.2:c.4145_4160inv` — on a CODING axis, with two three-base
    // unchanged interior runs — as one `inv`.
    //
    // **This restores the shipped baseline exactly.** `git show
    // 2fa5bec0:tests/it/issue_1040_inv_overrecognition_probes.rs` (the v0.13.1
    // release) asserted one loop over all four spellings, all expecting
    // `c.1_5inv`. The two-group split was introduced inside the partition-model
    // branch and measured against an intermediate tree in which the changed-column
    // weight bound was still vetoing every reading-frame merge, so the `inv` was
    // unreachable and the codon-merged form looked final. So this row is not a
    // representation change against the release — it is the reversal of one that
    // never shipped, and it restores four-way confluence.
    //
    // #1230's split is untouched and is asserted next door by
    // `a_whole_span_reverse_complement_is_not_merged_across_a_multi_base_separation`:
    // that case has a MULTI-base separation, so the same arm refuses it.
    for input in [
        format!("{C}:c.1_5inv"),
        format!("{C}:c.1_5delinsTTAAC"),
        format!("{C}:c.[1G>T;3T>A;5A>C]"),
        format!("{C}:c.[1_3delinsTTA;5A>C]"),
    ] {
        assert_eq!(
            c("GTTAA", &input),
            format!("{C}:c.1_5inv"),
            "spelling `{input}` did not converge on the inv",
        );
    }

    // The negative half, replacing the `assert_ne!` that used to assert the two
    // groups were deliberately distinct — that assertion pinned a divergence the
    // released version never had, so inverting it would keep the same mistake.
    // What needs guarding instead is the codon-merged form: it is what a
    // regression in `coalesce_whole_block_inversion`'s ranking would emit, it is a
    // legitimate `delins` on its own sub-span (`GTT -> TTA` is not a reverse
    // complement — revcomp(GTT) is AAC), so nothing downstream would reject it,
    // and this is the only assertion that would notice.
    for input in [format!("{C}:c.1_5inv"), format!("{C}:c.[1G>T;3T>A;5A>C]")] {
        assert_ne!(
            c("GTTAA", &input),
            format!("{C}:c.[1_3delinsTTA;5A>C]"),
            "`{input}` reached the codon-merged form; `coalesce_whole_block_inversion` \
             must outrank `apply_coding_codon_exception` on these bases \
             (`DNA/inversion.md:5`, ruling `inversion-vs-two-delins-76-83`)",
        );
    }

    // `revcomp_runs_in_distinct_codons_are_individual_on_cds_axis` above pins the
    // other side of the same axis — the case where neither rule may fire.
}

// ---------------------------------------------------------------------------
// Boundary probes: `canonicalize_delins` step 4d + `shorten_inversion`.
// ---------------------------------------------------------------------------

#[test]
fn shared_affix_trimming_reveals_inner_inv() {
    // ACGAGT -> ACTCGT: shared prefix "AC", shared suffix "GT"; the trimmed
    // core GA -> TC is a whole revcomp, so the inv is the trimmed 259_260 span.
    // This exercises the trimming path; the inv is spec-correct (whole trimmed
    // span is the reverse complement).
    assert_eq!(
        g("ACGAGT", &format!("{G}:g.257_262delinsACTCGT")),
        format!("{G}:g.259_260inv"),
    );
}

#[test]
fn length_changing_revcomp_delins_stays_delins() {
    // A revcomp-looking change whose ref/alt differ in length can NOT be typed
    // `inv` (4d requires equal trimmed lengths). ref GCT (257_259) -> alt AG:
    //   revcomp(GCT) = AGC, and the alt "AG" is a truncated (length-changing)
    //   reverse complement. No shared affix reduces it to a sub/del/ins/dup, so
    //   it stays a delins rather than being over-typed as an `inv`.
    assert_eq!(
        g("GCT", &format!("{G}:g.257_259delinsAG")),
        format!("{G}:g.257_259delinsAG"),
    );
}

// ---------------------------------------------------------------------------
// 2-nt inversions: the remaining ferro-vs-tool divergence candidate.
// Per DNA/inversion.md an inversion is ">1 nt", so a 2-nt reverse complement
// is a valid `inv`. ferro emits `inv`; some tools render it as a `delins`.
// ---------------------------------------------------------------------------

#[test]
fn two_nt_reverse_complement_emits_inv() {
    // TG -> CA (= revcomp(TG)).
    assert_eq!(
        g("TG", &format!("{G}:g.257_258delinsCA")),
        format!("{G}:g.257_258inv"),
    );
    // GG -> CC (= revcomp(GG)).
    assert_eq!(
        g("GG", &format!("{G}:g.257_258delinsCC")),
        format!("{G}:g.257_258inv"),
    );
}

// ---------------------------------------------------------------------------
// Reverse-complement runs separated by an unchanged base: individual `inv`s on
// a genomic axis (general.md:34 — the delins exception is codon-only).
// ---------------------------------------------------------------------------

#[test]
fn revcomp_runs_separated_by_identity_are_individual_on_genomic_axis() {
    // ref TGACA -> alt CAATG: 257_258 (TG->CA) and 260_261 (CA->TG) are each a
    // revcomp, separated by the unchanged base 259 (A). On a genomic axis they
    // are described individually, not merged into one delins.
    assert_eq!(
        g("TGACA", &format!("{G}:g.257_261delinsCAATG")),
        format!("{G}:g.[257_258inv;260_261inv]"),
    );
}

#[test]
fn revcomp_runs_in_distinct_codons_are_individual_on_cds_axis() {
    // Same change on a CDS axis: the two invs fall in *different* codons
    // (1_2 in codon 1, 4_5 in codon 2), so the codon-frame delins exception
    // (two variants separated by one nt affecting one amino acid) does not
    // apply and they stay individual.
    assert_eq!(
        c("TGACA", &format!("{C}:c.1_5delinsCAATG")),
        format!("{C}:c.[1_2inv;4_5inv]"),
    );
}

// ---------------------------------------------------------------------------
// #1040 near a short contig's 3' end — the path the real comparison exercises.
// ---------------------------------------------------------------------------

/// Genome-capable provider: one `contig` of `len` cyclic ACGT bytes with
/// `payload` written 1-based at `pos1` (mirrors the #1041 repro harness).
fn provider_len(contig: &str, len: usize, pos1: usize, payload: &str) -> JsonProvider {
    let mut bases: Vec<u8> = "ACGT".bytes().cycle().take(len).collect();
    for (i, b) in payload.bytes().enumerate() {
        bases[pos1 - 1 + i] = b;
    }
    let seq = String::from_utf8(bases).unwrap();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { contig: seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

fn norm(p: &JsonProvider, input: &str) -> String {
    Normalizer::new(p.clone())
        .normalize(&parse_hgvs(input).unwrap())
        .unwrap()
        .to_string()
}

#[test]
fn issue_1040_contiguous_run_near_3prime_end_stays_delins() {
    // ref[200..209] = CAGTGACTAG on a 250 bp contig: the variant end (209) plus
    // window_size (100) = 309 overruns the contig, so this takes the near-3'-end
    // path that the real comparison exercises. It must still be a single delins,
    // not a carved `[..delins;..inv;..delins]`. Holds on v0.8.1 and post-#1042.
    let near = provider_len("c", 250, 200, "CAGTGACTAG");
    assert_eq!(
        norm(&near, "c:g.200_209delinsTGTCACGACT"),
        "c:g.200_209delinsTGTCACGACT",
    );
}

#[test]
fn issue_1040_compound_allele_near_3prime_end_stays_delins() {
    // The reported input shape (a compound of adjacent substitutions), sited
    // near the 3' end. It collapses to the same single delins, not a carve.
    let near = provider_len("c", 250, 200, "CAGTGACTAG");
    assert_eq!(
        norm(
            &near,
            "c:g.[200C>T;201A>G;202G>T;203T>C;204G>A;205A>C;206C>G;207T>A;208A>C;209G>T]",
        ),
        "c:g.200_209delinsTGTCACGACT",
    );
}
