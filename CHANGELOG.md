# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.14.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/ferro-hgvs-v0.13.1...ferro-hgvs-v0.14.0) - 2026-08-16

### Representation changes

- *(normalize)* accept a duplication and an insertion at one junction ([#2008](https://github.com/fulcrumgenomics/ferro-hgvs/pull/2008))
  > a cis allele naming a duplication (or an
  > equivalent growing repeat) and an insertion at one shared junction now
  > settles on a single insertion whose payload is the copied bases followed
  > by the inserted ones — `c.[5_6dup;6_7insA]` becomes `c.6_7insAAA`, and
  > `g.[20_25GT[4];25_26insA]` becomes `g.25_26insGTA`. The effect depends
  > on the error mode and both halves are a migration. Under
  > `ErrorMode::Strict` the shape was refused with W5002, so those are
  > previously-rejected descriptions that now normalize and no
  > previously-accepted description changes form. Under the library default,
  > which is lenient (`NormalizeConfig::default()` uses
  > `ErrorConfig::lenient()`, `src/normalize/config.rs:129`), W5002 was only
  > a warning, so the same allele was already accepted and its normalized
  > string moves: measured on this branch against its merge-base,
  > `NM_TEST.1:c.[5_6dup;6_7insA]` normalized to
  > `NM_TEST.1:c.[5_6dup;6_7insA]` (preserved as authored) and now
  > normalizes to `NM_TEST.1:c.6_7insAAA`. A default-config consumer that
  > keys stored data on the normalized string therefore has a real migration
  > for this shape. Both cis-confluence censuses are byte-identical. Two
  > insertions at one junction remain refused. A same-junction pair whose
  > members were listed in the two possible orders previously produced two
  > different diagnostic strings and now produces one description.
- *(normalize)* re-present a contig-start insertion instead of refusing it ([#1969](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1969))
  > a pure insertion that 5'-shuffles onto interbase
  > 0 of a window beginning at the accession's first base is now
  > re-presented at the leftmost nameable interbase (`g.1_2insTA`) or folded
  > into a base-1-anchored delins, instead of refusing with "5' of the
  > window's first base".
  > Previously-refused inputs therefore start producing a description; no
  > description ferro already emitted is spelled differently, because the
  > escape fires only where the old path errored. A lone insertion genuinely
  > before base 1, with no piece to absorb it, still refuses.
  > Bounded to `w_lo == 1` — at any interior window the answer remains to
  > widen the 5' flank, not to relabel.
  >
  > ---------
- *(rulings)* type a whole-span reverse complement as inv, uniformly ([#1839](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1839))
  > 0 rows move in this PR. The normalizer is
  > untouched, the corpus is untouched, and every pinned string and census
  > constant is unchanged — including `CYP21A2_TARGET` and the inversion
  > sweep's `cases.tsv`, both of which this ruling decides *against* and
  > which are deliberately left green until the implementation moves them.
  > What is declared is the DECISION, and it is a real one: a whole-span
  > reverse complement will be typed `inv` uniformly, so #1541's
  > `NM_000500.9:c.[710T>A;713T>A]` becomes `c.710_713inv`, the rows the
  > `NM_004006.2` inversion sweep still repartitions return to their
  > authored `inv` — that population is `CENSUS_REPARTITIONED` = **92** on
  > the default arm at `origin/main` (93 on `live`), **not** the 155 this
  > record was first costed on and not only the 25 `delins`-bearing ones
  > #1575 costed — and #1230's checked-in guards are flipped. The
  > real-corpus measurement is **owed by the implementing change and is not
  > supplied here** — no A/B was run, because running one requires the
  > normalizer change this PR deliberately does not make. Disclosed in the
  > release the decision was made in, not only in the release that ships the
  > moves.
- allow repeat/complex variants to be converted to SPDI ([#1967](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1967))
  > two further insertion payloads now resolve
  > through `hgvs_to_spdi` instead of being refused — an exact tandem repeat
  > expanded under the existing `MAX_REPEAT_EXPANSION_BASES` cap, and a
  > bracketed compound insert whose parts mix literals, same-reference spans
  > and exact repeats, each part resolved by the leaf that resolves it
  > alone. A part whose bases are undetermined (uncertain or range count,
  > intronic CDS offsets, external accession) still declines the whole
  > insert, with one exception that predates this PR and is not closed by
- allow ref span insertions when converting to SPDI ([#1966](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1966))
  > a same-reference position-range insert
  > (`g.100_101ins50_57`, its `inv` and `delins` spellings, and the
  > single-part bracketed form) now resolves to its literal bases through
  > `hgvs_to_spdi` instead of being refused.
  > No normalized HGVS string moves — nothing under `src/normalize/` reaches
  > `src/spdi/convert.rs`, and the normalizer already expanded these
  > payloads through `rules::expand_inserted_sequence`; what changes is
  > which descriptions are convertible, not what a convertible one converts
  > to. On the `c.` and `r.` axes the payload span is resolved through the
  > same coordinate mapping as the location, so a coding-axis range insert
  > names the bases it says it does.
  > 18 of the 9,949,738 rows across the four release-asset corpora carry the
  > shape. Every affected input was previously rejected, so no consumer
  > holds a stored old value; the visible consequence is those rows leaving
  > `unkeyable` for a `group_by_spdi_key` bucket.
  >
  > ---------
- *(normalize)* define cis-member overlap once, on write footprints ([#1752](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1752))
  > 108 of 86,398 corpus rows move — 102 in
  > `repeat_beside_a_sibling` and 6 in `junction_interior_to_span`. Every
  > moved row is a member-conflicting input whose previous output was a
  > re-derivation that denoted different bases; the description is now left
  > exactly as authored, so all 108 become fixed points. **0** were fixed
  > points before, so no stored string migrates. **This supersedes an
  > earlier `6 of 86,020, all respell`**, which was measured on a corpus
  > that could build only GROWING repeats: the `repeat_beside_a_sibling`
  > family emitted `{unit}[3]` over a 6-base tract and never `{unit}[1]`, so
  > `repeat_footprint`'s shrinking branch — the branch this PR adds — was
  > structurally unreachable and its 102 rows were invisible. The family now
  > emits both directions and
  > `the_repeat_sibling_family_varies_the_repeat_direction` fails if it
  > stops.
- *(reference)* serve the synthesized transcript's exon table on its own axis ([#1783](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1783))
  > on a cdot-SYNTHESIZED transcript (one no
  > transcript FASTA record serves) whose cdot exon table is not contiguous
  > from 0 — a head offset, an internal hole, or an exon emitting a
  > different number of bases than its declared transcript span — the genome
  > alignment is now withheld: every exon's `genomic_start`/`genomic_end`
  > and the transcript's own `genomic_start`/`genomic_end` are `None`, so
  > `CoordinateMapper::tx_to_genomic` and every consumer above it
  > (normalize, VCF conversion, equivalence checking, W-code correction) now
  > REFUSE a c./n. position that previously received a genomic coordinate.
  > Measured on the prepared reference, 15 such transcripts are served, on
  > 14 distinct accessions — 4 on the GRCh38 primary path and 11 on GRCh37,
  > among them NM_001354366.1, NM_001387112.1, NM_004007.2 and NM_020979.5.
  > Nothing else moves: 482,352 of 482,519 GRCh38 and 197,404 of 197,846
  > GRCh37 cdot records have a contiguous zero-based table and keep their
  > genome axis unchanged, no FASTA-backed transcript is affected, and for
  > the 15 the served bases, the (gap-collapsed) exon table and the CDS
  > bounds are all still served. The genome coordinate for those accessions
  > remains available on the route stated on the deposited transcript's own
  > axis — `CdotTranscript::tx_to_genome`, used by the projector's
  > non-coding arm, `data/mapping.rs` and the convert service — which is
  > untouched. That route is why the withholding is a fix rather than a
- *(normalize)* decline a reversed range rather than reading it as a linear window ([#1923](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1923))
  > 6 rows move, all from a wrong description to a
  > decline. Measured by census, not by corpus. Over every parseable
  > reversed form x twelve spans on a 4000-base contig, 24 rows were
  > answered before this change: 10 panicked, 8 already declined, and 6
  > returned a wrong description - one dup answering g.4001_4000dup on a
  > 4000-base contig (past its own end), and five inv rows answering `=`, an
  > identity, for an inversion. Those 6 now return the description as
  > authored, which is what the other 8 already did. The 10 panicking rows
  > are not counted as moved, having had no previous answer to move from.
  > The synthetic corpus measures 0 rows moved and that figure is
  > STRUCTURAL, not evidence: dump_normalized_corpus.rs builds no reversed
  > range anywhere, so it cannot observe this change in principle - the same
  > blindness as #1460 and #1478. Report a structural zero as structural.
- *(convert)* refuse a CDS position that names no base instead of computing on it ([#1747](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1747))
  > Three parsed descriptions stop converting; no
  > normalized
  > output string moves. `NM_001234.1:c.pterdel`, `c.qterdel` and `c.cendel`
  > were
  > answered by `hgvs_to_spdi` with the triple `NM_001234.1:4:A:` — the same
  > triple it gives the ordinary `c.1del`, and the same one for both ends of
  > a
  > chromosome — and are now refused. **Previously accepted**, so a
  > migration for
  > anyone converting a terminus marker on a bare transcript reference: the
  > answer
  >   they hold is `c.1`'s, not the marker's. Direction is away from the
  > already-shipped string and toward a refusal, because there is no correct
  > string to move to — `pter`/`qter`/`cen` is a genomic landmark and names
  > no
  > position on a transcript's own CDS axis; where it resolves at all it
  > resolves
  > through `project_cds_terminus_to_parent` against a named genomic parent,
  > which
  >   is a different frame and is untouched. Not a corpus change:
  > `Normalizer::normalize` is byte-identical before and after over
  > `c.pterdel`,
  > `c.qterdel`, `c.cendel`, `c.pter_qterdel`, `c.1del`, `c.4del`, `c.*1del`
  > and
  > `c.-1del` on `JsonProvider::with_test_data`, measured on both arms with
  > the
  > same command, so the movement is confined to the SPDI conversion layer.
  > No
  > watched directory is touched — of the six (`src/normalize/`,
  > `src/hgvs/`,
  > `src/spdi/`, `src/project/`, `src/reference/`, `src/error_handling/`)
  > this
  > PR's four changed files match none, being `src/convert/mapper.rs` plus
  > three
  > under `tests/it/` — so this declaration is voluntary, and is made
  > because the
  > check gates on the directory while the effect is real. Two further arms
  > are
  > refused and are reachable only by a Rust caller constructing a `CdsPos`:
  > base
  > 0 with `utr3` set (`c.*0`, which `parse_hgvs` refuses) and the `c.?`
  > sentinel,
  > which `hgvs_to_spdi` already refused a layer up at `parser/variant.rs`.
  >   Neither is reachable from the Python bindings any more — #1741 added
  > `reject_zero_base`, which `c_to_g`, `c_to_p`, `c_to_n` and `n_to_c` all
  > call
  > before converting, so the Python reach this PR originally claimed for
  > those
  > two arms closed at the rebase and is withdrawn. An extreme coordinate
  > moves
  > from a panic or a silent wrap to a refusal. To reproduce output from
  > before
  > this change, resolve a terminus marker against a named genomic parent
  > rather
  >   than against a bare transcript.
- *(normalize)* clamp the 3' rule at an exon junction from both sides ([#1845](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1845))
  > 13 of 85,642 rows move (0.0%) on the designed
  > shape-family corpus — `coincident_bounds` 9, `repeat_expansion` 4 — and
  > **all 13 were previously not fixed points**, so 0 previously-accepted
  > rows move and there is no consumer migration. Quote that rate as of the
  > affected families, not repo-wide. Independently, over
  > `spec_conformance_axis`'s 58,552 outputs exactly four family instances
  > change, all `junction-1-del-del`: `s01-c3p-…-sep{0,1}` fall one arity
  > each as the far-side spelling `c.21_22del` converges onto `c.19_20del`,
  > and `s01-c3m-…-sep{0,1}` re-spell that same output as
  > `NC_SYNTH.1(NM_TEST.1):c.18_20+2A[3]`. No family loses or gains
  > convergence, `converged` is unchanged at 11,016, no failure counter
  > rises, and the 5' census does not move at all. The spec's own row
  > `LRG_199t1:c.3922dup` → `c.3921dup` moves its projected genomic form
  > ~2,790 bp (that figure is the ruling's, relayed, not re-derived here).
- *(normalize)* report REFSEQ_MISMATCH at a coordinate, not a window offset ([#1906](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1906))
  > The emitted WARNING SET changes on cis alleles;
  > no normalized output moves (0 rows respell).
  > This is NOT a diagnostic-only change, and the previous `none` verdict
  > was wrong.
  > `ValidationResult::valid` is untouched and the new `frame` argument is
  > read at exactly one place — where the
  > warning string is built — so no normalization decision and no normalized
  > string moves. What moves is the
  > set of warnings a consumer receives, plus the text of the ones that
  > survive.
  > (1) WARNING SET. The `REFSEQ_MISMATCH` `position` string is the sole key
  > of the cis-allele duplicate
  > suppression (`src/normalize/mod.rs:4160-4182`), so repairing the
  > coordinate changes how many warnings
  > survive. Measured on `NC_000001.11:g.[201C>T;300G>A]` over a 400-base
  > ACGT-cycling contig: 4 warnings
  > before (`["100-100","100-100","201-201","300-300"]`), 2 after
  > (`["201-201","300-300"]`) — pinned by
  > `each_genomic_cis_member_is_reported_exactly_once`. A consumer reading
  > `ferro normalize -f json`'s
  > `warnings` array, or the Python `NormalizationWarning` list, sees two
  > entries where it saw four.
  > (2) TEXT, on five surfaces beyond `Display`, four of them
- *(reference)* model the LRG <diff> list so indel-bearing parents project exactly ([#1847](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1847))
  > 51 of 5176 LRG projection rows change status, all
  > `unavailable` -> `ok`
  > (2764 -> 2815 ok, 147 -> 96 unavailable under --error-mode lenient; 43
  > rows and 2489 -> 2532
  > under the CLI-default strict mode, the 8-row difference being inputs
  > strict refuses before a
  > placement is consulted). 0 previously-answered row changes value, and 0
  > regress. Measured
  > against origin/main AFTER #1830 merged: #1830 narrowed indel-bearing LRG
  > placements to their
  > affine-exact prefix and declined the rest; this models the alignment
  > instead, so those rows are
  > answered rather than declined. Direction: TOWARD the reference oracle --
  > all 51 recovered
  > coordinates were adjudicated against the frozen LRG_<N>g records with an
  > unshiftable
  > substitution probe (11-base window of the record against the chromosome
  > coordinate cdot alone
  > resolves): 51 agree, 0 disagree. 49 of the 51 sat at a base the gapless
  > affine got demonstrably
  > wrong (drift -8..+23); the other 2 sat on a cancelling indel pair. The
  > earlier "55 move / 0
  > status" figure was measured against the pre-#1830 base and no longer
  > describes this merge.
- *(hgvs)* stop dropping and leaking in-band position markers ([#1762](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1762))
  > `n.<b>+?`, `n.<b>-?` and `r.<b>+?` now render as
  > written
  > instead of leaking the in-band sentinel — `n.5+?` previously
  > round-tripped as
  > `n.5+9223372036854775807`. 3 outputs move. No committed corpus contains
  > an
  >   input that reaches them: the real corpora carry 709 unknown-offset
  > descriptions and every one is on `c.`/`g.`, the two axes that already
  > rendered
  >   `?`, and `n.*N` occurs 0 times in 103,762 `n.`-axis rows. Separately,
  > `tx_to_cds` and the projection fan-out seed now refuse `n.*N` where they
  > previously answered — an `Ok` becomes an `Err`, so no string moves, but
  > it is
  > a behaviour change and is disclosed here rather than left to be
  > discovered.
- *(normalize)* make a cis partition independent of the shuffle direction ([#1840](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1840))
  > 54 of 85,642 corpus rows move, 44 merge / 10
  > respell. All 54 fall in the two inversion families (`inv_member` 38,
  > `delins_hiding_an_inversion` 16); the other 19 families move 0. 19 of
  > the 54 are on the shipped 3' direction, 35 on the internal 5' oracle.
  > NOT a migration: all 54 were previously not fixed points, so no stored
  > string was stable for any of them, and 38 become fixed points. Measured
  > with examples/dump_normalized_corpus.rs against origin/main@2d8b490b
  > with only src/normalize/merge.rs swapped. The previously declared 150
  > was measured on the pre-#1835 base and is withdrawn.
- *(project)* stop renaming the n. axis onto an LRG transcript with a different frame ([#1844](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1844))
  > 2 of 681 projected non-coding axis rows respell.
  > Previously-accepted output; both moved rows previously named bases they
  > did not denote, so this is a correction rather than churn.
- *(normalize)* refuse a c./r. position on base 0 instead of answering it as -1 ([#1799](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1799))
  > Direct-construction only; a structural zero over
  > any parsed corpus. A `c.`/`r.` position built on base 0 was answered as
  > this axis's own `-1` — `NM_TEST.1:c.?del` normalized to `c.-1del` — and
  > is now left as authored. No parseable input can carry the shape: the
  > parser refuses `c.0` (#269), a bare `c.?` maps to
  > `UncertainBoundary::Single(Mu::Unknown)` whose `inner()` is `None` and
  > so exits before the conversion, and with #1777 in the base the cis
  > collapse no longer manufactures a base-0 anchor. Measured: the full
  > suite is unchanged at 10,442 passing, and the four seam oracles show 17
  > failures on this branch against the same 17 on its base.
- *(reference)* narrow an LRG placement to its affine-exact prefix ([#1830](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1830))
  > 51 of 5176 LRG projection rows move, all
  > previously-accepted and all to a refusal; 0 accepted rows change their
  > string. Measured over the 1294 LRG records of the prepared reference; 49
  > of the 51 named a demonstrably wrong base (drift -8..+23) and 0 were
  > below their record's affine horizon.
- *(normalize)* refuse to clamp a repeat onto bases that are not its unit ([#1836](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1836))
  > 98 rows respell, of 98 measured over the #1597
  > corpus. All 98 previously denoted a **different sequence than their own
  > input**, so every moved string was wrong — this is a correction, not a
  > re-spelling of correct output. Separately, 0 of 85,642 rows move in the
  > committed corpus harness; that one is a **structural** zero and is
  > evidenced as such below.
- *(spdi)* refuse a genomic offset inside a complex boundary, and an end that names no coordinate ([#1807](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1807))
  > Previously-accepted input is now refused on the
  > SPDI conversion path; no normalized output moves.
  > `hgvs_to_spdi`/`hgvs_to_spdi_simple` now decline a `g.`/`m.`/`o.`
  > description whose complex `(a_b)` boundary carries a `+`/`-` offset or a
  > `pter`/`qter`/`cen`, and one whose end boundary names no single
  > coordinate (`(a_b)` or `?`); each previously converted by dropping the
  > boundary and collapsing the end onto the start. `normalize` does not
  > route through `hgvs_to_spdi`, and 0 of 85,642 `dump_normalized_corpus`
  > rows can carry any of these shapes — a **structural** zero, not a
  > measured one: the generator spells every position with `format!` from a
  > plain integer, so it cannot vary the property this change keys on.
- *(hgvs)* refuse a bare-transcript intronic position at the strict parse stage ([#1872](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1872))
  > 1,031,655 of 5,163,595 parseable `c.`/`n.` corpus
  > rows change their strict-mode PARSE verdict, from accepted to refused
  > (`W4007`). **No normalized output moves — 0 rows respell**, and nothing
  > under `src/normalize/` changes what it emits; what moves is the
  > accept/reject boundary of `parse_hgvs_with_config(_,
  > ErrorConfig::strict())` for a bare-transcript intronic description
  > (`NM_…:c.20+2del`, `NR_…:n.20+2del`, and inside an allele). Those inputs
  > were already **rejected** in strict mode one stage later, at normalize,
  > by a rung that short-circuits before normalization — so a caller that
  > parses and then normalizes strictly sees the same verdict, reached
  > earlier, with `FerroError::Parse` in place of
  > `FerroError::IntronicVariant`. The genuinely new refusal is for a strict
  > caller that parses **without** normalizing, and that is the migration
  > this discloses: it is a previously-accepted input, not a
  > previously-rejected one. Lenient parse newly emits a `W4007` warning
  > where it was silent; silent mode is unchanged; the CLI and the Mutalyzer
  > conformance runner both use the config-less `parse_hgvs` and are
  > unaffected.
- *(normalize)* refuse a c. description whose reference resolves no CDS ([#1882](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1882))
  > 67,576 of 225,662 served `c.`-addressable
  > accessions (~30%) move from answered to REFUSED — 1,242 curated `NM_`
  > and 66,334 `XM_` predicted models, measured through ferro over the
  > prepared GRCh38 reference, not extrapolated. These inputs were
  > previously ACCEPTED (they returned a value, `status: ok`, `changed:
  > false`, `warnings: []`, in every mode), so this is a real migration for
  > any consumer keying on the normalized string. The sharpest
  > consumer-facing fact is deduplication: on an affected accession
  > `normalize` was the IDENTITY function, so all five of
  > `NM_000546.3:c.528del`..`c.532del` — one deletion from one `CCCCC` run,
  > one variant — came back as five distinct "normalized" strings and a
  > caller deduplicating by normalized string got five variants where there
  > is one. The clinically weighty `NM_` rows are overwhelmingly SUPERSEDED
  > transcript versions (601 of 1,242 have a higher served version that
  > resolves a CDS), i.e. a legacy-input exposure — the version pinned in an
  > older report, a database export or a lab's calling config, not the
  > version a fresh pipeline picks. 27 of the ACMG SF v3.2 rows named in
  > #1870 reproduce as refused, `TP53 NM_000546.3` and `BRCA1 NM_007294.2`
  > among them; that is a legacy-input exposure and not "ferro silently
  > mis-normalizes TP53". Nothing moves on an accession whose CDS the
  > reference resolves. Committed-corpus movement: 10 rows on the
  > spec-enumeration corpus (all `NM_004006.1`, divergence budget unchanged)
  > and 6 rows on the real-corpus multi-member cis census (unwindowed ->
  > declined, three transcripts).
- *(project)* refuse an `n.*N` endpoint on the genome-pivot path ([#1781](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1781))
  > outputs move for two shapes. An `n.*N` input on
  > the genome-pivot path previously produced a projection — measured as
  > `g.1012del` / `c.*1del` / `n.13del` on a coding transcript, and a `g.`
  > coordinate identical to the corresponding `n.<base>` on a non-coding one
  > — and now returns `UnsupportedProjection`. The same now holds for an
  > `r.*N` input on a **non-coding** transcript, which previously returned
  > the same genomic coordinate as the corresponding `r.<base>` (measured:
  > both `NC_000001.11(NCGENE):g.1004del`) and now returns
  > `UnsupportedProjection`. **0 corpus rows move for `n.*N`**: `n.*N`
  > occurs 0 times in the four committed corpora (`clinvar_hgvs_500k`,
  > `clinvar_hgvs_unique`, `cmrg_genes_exhaustive`,
  > `paraphase_genes_exhaustive`, 103,762 `n.`-axis rows), and that zero is
  > real rather than structural — the same scan finds the sibling `n.-N`
  > five times. The corpus incidence of `r.*N` was not separately measured,
  > so that zero is claimed for `n.*N` only. No stored string changes; what
  > changes is that a description ferro could not answer correctly is now
  > declined instead of answered wrongly.
- *(spdi)* refuse a c./n./r. end boundary that names no coordinate ([#1823](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1823))
  > previously-accepted input is newly refused on the
  > `hgvs_to_spdi` / `hgvs_to_spdi_simple` path, so this is a real migration
  > for any consumer converting a `c.`/`n.`/`r.` description whose end
  > boundary is a range `(a_b)` or `?`. 0 of 85,642 corpus rows move, and
  > that zero is **structural** — see *Representation stability* below for
  > the three checks that establish it. `normalize` output is unaffected: it
  > does not route through `hgvs_to_spdi`.
- *(normalize)* [**breaking**] make canonical-coalesced the default partition rule ([#1835](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1835))
  > this is the largest representation change the
  > project has made. On the 11,272-class designed cis corpus every
  > divergent class but one converges: 2,910 classes change their normalized
  > form at 3' and 2,905 at 5', with `sequence_changed` 0 on both sides, so
  > every move is a re-spelling and none changes denoted bases. Direction of
  > movement is overwhelmingly merge — a payload-coincidence split becomes
  > the spanning `delins` on the `c.` axis, and members that `live` left as
  > the input's own spelling are re-derived from the resulting sequence. On
  > the manifest-backed conformance axes a further 45 rows move, all one
  > mechanism and all net insertions where the split form is canonical. A
  > consumer keying on the exact normalized string must expect
  > previously-accepted inputs to return a different, still-conformant
  > description, and must pin `FERRO_PARTITION=live` to reproduce
  > pre-v0.15.0 output. The real-corpus row count over
  > ClinVar/CMRG/Paraphase is NOT supplied here and remains owed; it is
  > tracked as #1885. The earlier reason given for its absence — that the
  > unresolved inversion family would contaminate it — is stale now that
  > #1706 has merged, so the honest statement is that it is derivable and
  > has not been run.
  >
  > **Editorial correction (#1886) — not the trailer's words.** The trailer names
  > `pre-v0.15.0` as the output `FERRO_PARTITION=live` reproduces. The release carrying
  > this change is the one heading this section, so read that as *output from before
  > this change* — every release up to and including v0.13.1. `release-plz.toml` sets
  > `features_always_increment_minor`, so which release a `feat:` lands in is not
  > knowable when its trailer is written; naming one is now against house style, and
  > `CONTRIBUTING.md` asks for the version-free phrasing instead.
  >
  > **Editorial correction (#1886) — not the trailer's words.** The trailer discloses
  > confluence only. The flip also moves **2 rows of 932** on the HGVS spec corpus,
  > measured per arm at `45820926` against a prepared GRCh38 reference, 3' direction,
  > with the two rows #1846 cannot terminate on excluded: `live` 646, `shadow` 646,
  > `canonical` 642, **`canonical-coalesced` — the new default — 644**. Two rows stop
  > matching the spec's stated form and none starts, and both are named here because
  > the bare count reads as a regression and neither row is one. `DNA/delins.md`
  > publishes one variant as two corpus rows — the spanning
  > `c.850_901delinsTTCCTCGATGCCTG` that `:47` recommends, and the split `:46` calls
  > its "alternative description" — and both echoed themselves before; the new
  > default converges the split onto the spanning form, which is what any change that
  > converges a published pair must cost. The other row,
  > `LRG_199t1:c.992_1004delinsAC`, is harvested from `consultation/SVD-WG010.md` —
  > the **rejected** proposal, whose answer `tests/it/spec_worked_example_rules.rs`
  > already pins as one ferro must not produce — and the new default no longer echoes
  > it either. Read this as neither more nor less conformant: the measurement supports
  > neither claim. 3' only, because #1879 removes the 5' direction from the public
  > surface in this same release.
  >
  > **Editorial correction (#1885) — not the trailer's words.** The real-corpus row
  > count this trailer defers is now measured. Over all four bulk corpora in full —
  > **9,949,738** stored ClinVar, CMRG and Paraphase expressions, each normalized
  > through both arms of one build — **11,491 normalized strings move (0.115%)**:
  > 11,476 respellings whose SPDI key is unchanged, 15 whose denotation SPDI cannot
  > key, and **0 that denote different bases**. No row starts or stops normalizing.
  > A further 9,418 rows (0.095%) are excluded by a deterministic coordinate-span
  > bound rather than found unchanged — the #1846 cost defect, which cannot
  > terminate on them — and are recorded as excluded rather than folded into the
  > total. The measurement keys on which partitioner cuts the block, a property
  > 86.5% of the corpus reaches, so the zero above is a measured zero and not a
  > structural one.
- *(normalize)* anchor a 5'-edge cis insertion on -1, not the axis's absent 0 ([#1777](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1777))
  > 1 shape moves. A cis allele that nets to a pure
  > insertion at the 5' edge of the window on the `n.` axis respells from
  > `n.0_1insA` to `n.-1_1insA` — measured on `NM_TEST.1:n.[1G>A;1dup]`. The
  > previous output names a position that does not exist and `parse_hgvs`
  > rejects it, so no consumer can be holding it as a parsed value; this
  > replaces an unreadable description with a readable one. The `c.` and
  > `r.` equivalents do **not** move (`c.-1dup` and `r.-1dup` on both sides,
  > measured), because their malformed intermediate was already repaired
  > downstream by the `base == 0` conversion arms. No genomic-axis row
  > moves, and no row away from the axis origin moves.
- *(reference)* never trade a cdot exon table for a supplemental record that carries nothing ([#1724](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1724))
  > 50 accessions move from declining ("Transcript
  > has no CDS") to served on the c./n. axes. Previously-erroring inputs now
  > produce output, so this is a real change for a consumer that keyed on
  > the decline.
- *(normalize)* recognise inversions post hoc over blocks and runs of pieces ([#1706](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1706))
  > 62 of 2,075 rows move. An authored inversion that
  > was shredded into an allele comes back as the single `inv`;
  > previously-accepted inputs, so a real migration for a consumer keying on
  > the string. Denotation is unchanged on all 62, verified through the SPDI
  > applier with a firing negative control. Re-derived by set difference
  > against `origin/main` at `cc8407bc`, not carried forward from any
  > earlier measurement: the same 62 rows. Denotation is re-stated from the
  > `1ea75334` run rather than re-derived here; see the table below, which
  > marks which figures were re-run at this base and which were not.
- *(spdi)* emit a transcript-axis SPDI on the transcript, not the compound accession ([#1821](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1821))
  > SPDI output only (`hgvs_to_spdi` /
  > `canonical_spdi`), and only for compound `NG_xxx(NM_yyy)` references on
  > the c./n./r. axes — the `sequence` field becomes the inner transcript
  > and the del/ins bases it carries become the transcript's instead of the
  > genomic parent's. No HGVS description string changes: normalize and
  > project output is byte-identical, and the g./m./o. axes are untouched.
  >
  > ```
  > NG_008939.1(NM_000532.5):n.192_197del
  >   before ->  NG_008939.1(NM_000532.5):191:CACACA:      (parent's bases at a transcript offset)
  >   after  ->  NM_000532.5:191:ACGCCG:                   (= what the bare NM_000532.5 spelling already gave)
  > ```
  >
  > The affected rows' normalized outputs are byte-identical (`c.101del` →
  > `c.102del`, and so on for the other six). The genomic axes keep
  > `Display` deliberately: there the coordinates *are* the parent's, so
  > stripping the wrapper would name the transcript for a genomic offset —
  > the same defect with the frames swapped.
- *(normalize)* carry junction-exit provenance per leaf, from the gate that makes it ([#1726](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1726))
  > 1 shape moves, and it is a repair — a
  > manufactured intronic offset on a SECOND bare transcript accession now
  > carries the genomic reference `checklist.md:20` requires, where #1704
  > left it bare. 0 rows of any committed corpus or fixture move: the 3'
  > census reads `outputs_leaving_the_transcript: 0` /
  > `outputs_intronic_under_a_genomic_wrapper: 371`, byte-identical to
  > #1704. That zero is STRUCTURAL and is a claim about the corpus, not
  > about the change — the spec corpus's junction stratum builds
  > single-accession descriptions and so cannot emit a two-accession allele
  > at all; the moved shape is constructed explicitly in
  > `defect_371_transcript_exit::a_second_bare_accession_is_repaired_too`.
  > Every other output, including the mixed-allele residue, is
  > byte-identical to #1704.
- *(normalize)* test the span the codon-frame merge authorises, not one edge ([#1720](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1720))
  > 6 of the 1,221 coding separation-one corpus rows
  > move, and 0 of the 592 real multi-member cis alleles harvested from
  > 9,949,738 ClinVar/CMRG/Paraphase rows. Every moved row is a merge whose
  > `delins` span crossed a codon boundary; the merge is withdrawn and the
  > two variants are described individually per `general.md:34`. The
  > real-corpus zero is a property of that corpus rather than of the change
  > — the reproducers are real transcripts (`NM_004006.2` / `LRG_199t1`),
  > the shape just does not occur in ClinVar.
- *(convert)* resolve a c./n. position on the flat transcript axis ([#1735](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1735))
  > moves representation on transcripts whose cdot
  > exon table has a transcript-coordinate gap (58 of 474,818 GRCh38
  > multi-exon builds, 23-2718 bases; 159 of 190,754 on GRCh37). On such an
  > accession every c./n. position, and therefore every SPDI triple and
  > every applied base derived from one, moves onto the flat transcript
  > offset the accession's own numbering names. The movement is NOT the gap
- render a junction-crossing output against its genomic reference ([#1704](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1704)) ([#1708](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1708))
  > 371 rows respell (389 under
  > `canonical-coalesced`), all accession-only; 4 of 2,075 real-transcript
  > inversion-sweep rows and 6 spec-enumeration rows move, 2 of those from
  > "no genomic representation" to a rendered `g.`. No coordinate moves
  > anywhere.
- *(hgvs)* refuse an alignment-only symbol on the r. axis ([#1721](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1721))
  > yes — descriptions previously accepted and
  > re-emitted are now refused. Every `r.` description whose insert reaches
  > `InsertedSequence::Named` and states a non-leading lower-case `x` moves
  > from re-emitted to refused; 0 respell, 0 merge, 0 split. 0 corpus rows
  > move, and that zero is STRUCTURAL rather than reassuring —
  > `RefShape::all()` builds no `r.` reference at all, so 0 of 58,552
  > spellings are on this axis and no corpus counter could move in either
  > direction; the population is established by construction instead, which
  > is what `REACHING_SHAPES` is. Every affected description states a masked
  > base and therefore denotes no sequence, so nothing that currently
  > produces a legal description changes. The refusal is newly reachable by
  > lenient and silent callers, which previously got the string back; strict
  > callers now fail one stage earlier, at parse. The generated spec fixture
  > and enumeration are byte-identical against a baseline generated from
  > `origin/main`'s sources — 934 and 2172 rows, 0 added, 0 removed, 0 field
  > differences.
- *(normalize)* require a gap-bearing insert before the delins.md:44-47 merge ([#1698](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1698))
  > 4,616 of 8,855,661 compared rows split on the
  > unstable `FERRO_PARTITION=canonical-coalesced` evaluation arm across
  > ClinVar 500k, ClinVar unique, CMRG and Paraphase, plus 4 CMRG rows that
  > become clean `W4007` declines. The shipped `live` arm never runs this
  > pass and 0 rows move there — measured on the same four corpora against
  > base `ce933533`, and re-confirmed at `2f4e3bb9` by regenerating both
  > spec artifacts at base and head and comparing them byte-identical. 0
  > rows change the bases they denote, verified with
  > `Normalizer::canonical_spdi` over all 4,002 unique changed pairs against
  > a negative control that fired 3,050 times.
- *(hgvs)* refuse an alignment-only symbol, at the stage the ruling names ([#1684](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1684))
  > 24 rows move (3' and 5' alike), all from
  > re-emitted to refused; 0 respell, 0 merge, 0 split. Every affected row's
  > output is invalid today — a description stating `X` denotes no sequence
  > — so nothing that currently produces a legal description changes. The
  > refusal is newly reachable by lenient and silent callers, which
  > previously got a string back; strict callers now fail one stage earlier,
  > at parse. The generated spec fixture and enumeration are byte-identical
  > (934 and 2172 rows, zero status changes).
- *(normalize)* grow the fetch window so a long homopolymer converges ([#1697](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1697))
  > two classes move. (1) A `g.` or `m.` description
  > whose 3' shuffle runs past `window_size` (100 bases) now converges on
  > the tract's 3' end instead of advancing 100 bases per normalize call.
  > (2) Past the 64 KiB growth cap the shift is refused but the
  > reference-aware edit type is no longer refused with it, so a description
  > that is really a duplication, substitution or inversion is now spelled
  > as one — `g.65_66insA` -> `g.67dup`, `g.65delinsG` -> `g.65A>G` — where
  > before it was echoed in the authored edit type. Repeat notation stays
  > refused past the cap, because its `<first>_<last>unit[N]` extent is
  > unreadable there. Measured 0 of 85,642 corpus rows moved, but that zero
  > is STRUCTURAL and is not evidence of safety — the corpus's largest shift
  > is 14 bases and no member is 64 KiB long, so it cannot reach either
  > changed path.
- *(normalize)* let the splitter express two deletions, not just two insertions ([#1649](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1649))
  > 9 confluence censuses move and none regresses;
  > `converged` rises 8027->8361 (3') and 8023->8367 (5') in
  > cis_confluence_axis, 9141->9402 and 8944->9228 in spec_conformance_axis,
  > and by 167/167/172/172 across cis_confluence_nr_axis's n./r. x 3'/5'
  > pins, while multi_member_cis_axis's respelling_converged goes 395->401
  > with 0 rows lost - that last one measured against a local prepared
  > reference and NOT pre-merge coverage, since multi_member_cis_axis is
  > manifest-gated and skips silently on PR CI; only the nightly runs it.
  > In every one of those censuses the three split-arity deltas sum EXACTLY
  > to the converged delta, so each moved family converged outright rather
  > than dropping an arity.
  > Measured by dumping the full divergence row-id lists on origin/main and
  > on this branch and diffing the sets, not by reading the nets: 0 families
  > lose convergence, 0 rise in arity, and 0 become divergent that were not,
  > on any corpus in either direction.
  > Moving the other way, 3 spec-enumeration rows go from projection-pinned
  > to projection-splits-single-member (1168->1165 and 9->12) - the three
  > projections project-{c,n,r} of NM_004006.3:r.123_127delinsag,
  > sequence-preserving via spdi::compare_denoted_sequences, ruled a fix on
  > n., unruled on r., and one deferred residue on c.
  > Six previously-accepted pinned strings also move: the 3' and 5' outputs
  > of 1419-r1/span, 1419-r2/span and 1419-r3/span in
  > reported_partition_verdicts, which converge each pair onto its /cis
  > sibling's form - a real migration for any consumer storing those.
- *(normalize)* carry the reference copies a re-phased repeat absorbs ([#1678](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1678))
  > 13 shapes move, all on an anchored
  > single-position repeat
  >   whose spelled unit tiles a different tract than the rotation the tract
  > maximization picks. Measured against origin/main (d5f26fcb) on the
  > branch's own
  > synthetic references (core base 1 = g.257), 3'-shuffle, in three
  > classes.
  > A. SEEDED OVERLAP -- and here THE EMITTED EDIT KIND CHANGES, not only
  > the count.
  > The spelled unit's tract SPANS the anchor, so it seeds the maximization,
  > and a
  > strictly longer overlapping rotation displaces it; the count is then
  > re-based by
  > the copies the wider window swallows. Because a repeat whose count sits
  > one unit
  > either side of the reference tract's count is rendered as `del` or `dup`
  > rather
  > than as a repeat, shifting the count by one also shifts WHICH count
  > renders as
  > WHICH edit. On core ATGTGTGTGA anchored at g.259 with unit GT -- a
  > 3-copy GT
  > tract sitting inside a 4-copy TG tract -- every count 1..8 moves and
  > four of them
  >   change kind:
  >     g.259GT[1]  g.258_265TG[1]  -> g.258_265TG[2]
  >     g.259GT[2]  g.258_265TG[2]  -> g.264_265del    repeat -> del
  > g.259GT[3] g.264_265del -> g.259GT[3] del -> repeat, left alone
  >     g.259GT[4]  g.259GT[4]      -> g.264_265dup    repeat -> dup
  >     g.259GT[5]  g.264_265dup    -> g.258_265TG[6]  dup    -> repeat
  >     g.259GT[6]  g.258_265TG[6]  -> g.258_265TG[7]
  >     g.259GT[7]  g.258_265TG[7]  -> g.258_265TG[8]
  >     g.259GT[8]  g.258_265TG[8]  -> g.258_265TG[9]
  > So a consumer holding a stored `repeat` for one of these inputs can read
  > back a
  > `del` or a `dup`, and one holding a stored `del` or `dup` can read back
  > a repeat.
  > This class is NOT merely a count change, and "the count now carries the
  > copies
  > the wider window swallows" does not predict it. The family is unbounded
  > in the
  >   count; 1..8 is what is measured and asserted.
  > B. OVERLAPPING but NOT seeded -- the spelled unit's tract ends AT the
  > anchor, so
  > the span filter rejects it as a seed. The count carries the copies the
  > wider
  >   window swallows and the edit kind is unchanged --
  >     AAGTGTTA      g.262TG[6]  g.259_262GT[6] -> g.259_262GT[7]
  > C. DISJOINT tracts -- the re-phase is declined and the count stands
  > against the
  > tract it was written for, because relocating the window denotes
  > different bases
  >   whatever the count --
  >     AATGTGTGGTAA  g.265TG[6]  g.265_266GT[6] -> g.259_264TG[6]
  >     AATGTGTGGTAA  g.265TG[2]  g.265_266dup   -> g.263_264del
  >     AATGTGTGGTAA  g.265TG[1]  g.265TG[1]     -> g.259_264TG[1]
  >     GGACCAGG      g.261AC[5]  g.261_262CA[5] -> g.259_260AC[5]
  >   Every new output denotes its input's bases, checked with
  >   `spdi::compare_denoted_sequences` through `hgvs_to_spdi`, not with the
  > normalizer. Eleven of the thirteen previously denoted a DIFFERENT
  > sequence rather
  > than merely a different spelling -- 7 of the 8 in A and 4 of the 5 in B
  > and C
  > together (an earlier revision of this trailer said 3 of 5; re-measured,
  > it is 4,
  > because GGACCAGG g.261AC[5] shipped a same-LENGTH, different-BASES
  > output) -- so
  > a stored string of those shapes changes and the new one is the correct
  > one.
  > Unmoved and asserted: a pure re-phase (`repeated.md:97`'s CFTR
  > `g.258TG[13]`), an
  > expansion (`g.258GT[17]`), and the arm where the spelled unit tiles
  > nowhere at
  > the anchor (`g.259TG[6]`), which is pinned as a limit rather than fixed
  > --
  > `hgvs_to_spdi` reads that input as denoting no sequence, so there is no
  > baseline
  >   to conserve. Not measured against a real corpus:
  > `examples/dump_normalized_corpus.rs` cannot build a repeat whose spelled
  > unit is
  > out of phase with a different tract, so a zero from it would be
  > structural rather
  >   than evidence.
  >
  > Closes #1618
- *(error_handling)* keep the repeat unit when the W3013 range disagrees with it ([#1661](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1661))
  > 1 pinned (row, mode) result moves —
  > `NM_024312.4:r.-6_-3g[6]` under `--error-mode lenient` and `--error-mode
  > silent`, from `r.-6_-3[6]` to `r.-6_-3g[6]`. Previously **accepted**, so
  > it is a real migration for a consumer storing that string — but the
  > previously stored string denoted a different variant (24 nt against 6
  > nt), so the migration is a correction rather than a re-spelling. The
  > affected family is exactly the W3013 repair path where the range's span
  > differs from the unit's length; the equal-length arm is unchanged and
  > pinned. Measured by the full suite (one pinned row moved, nothing else)
  > and by the four-oracle A/B above. `examples/dump_normalized_corpus.rs`
  > is **structurally incapable** of measuring this and was not run for it:
  > it parses with `parse_hgvs`, which applies no `ErrorConfig`, so no
  > corpus row ever reaches the preprocessor and a `0 moved` from it would
  > be a fact about the harness.
- *(project)* decide the codon-frame delins exception per rendered axis ([#1672](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1672))
  > 3 rows move, all on the transcript axes of a
  > codon delins
  > authored on the genomic axis, whose partition was decided without a
  > frame and
  >   then inherited. `REF:g.4_6delinsATA` projects as `c.4_6delinsATA` (was
  >   `c.[4C>A;6G>A]`), `n.4_6delinsATA` (was `n.[4C>A;6G>A]`) and
  > `r.(4_6delinsaua)` (was `r.[(4c>a;6g>a)]`). Toward the form
  > `DNA/delins.md:42`
  > asks for on an axis that declares a reading frame; the `r.` row merges
  > on
  > `RNA/delins.md:18`'s own authority. The `g.` axis does NOT move -- it
  > stays
  > `g.[4C>A;6G>A]`, which is #79's scope, and the projection census
  > constants are
  > byte-identical to main (`ProjectionSplitsSingleMember` 9,
  > `ProjectionPinned`
  > 1168), so the two `LRG_199t1` spec rows this PR argues about are
  > unchanged
  > rather than reversed. Previously-accepted inputs, so a migration on
  > paper for
  > anyone storing a projected `c.`/`n.`/`r.` string of this shape.
  > Normalization
  >   itself is untouched.
- *(normalize)* decline a protein delins whose affix trim lands on the initiation codon ([#1670](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1670))
  > a protein `delins` whose affix trim leaves a
  > single-residue residual at residue 1 replacing the initiator `Met` is
  > now left as authored instead of being rewritten to a start-loss
  > substitution. The old output (`p.Met1_Ala3delinsValAlaAla` ->
  > `p.Met1Val`) **does not parse**, so no consumer can be holding a legal
  > string that changed. 0 rows move across the four committed corpora.
  >
  > Closes #1607
- *(normalize)* scope the delins payload-coincidence carve-out to coding DNA ([#1682](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1682))
  > The shipped default is unaffected —
  > `FERRO_PARTITION` unset is `live`, which never runs this pass; measured
  > at **0 rows moved** over 785,461 normalized ClinVar and Paraphase rows.
  > What moves is the opt-in `FERRO_PARTITION=canonical-coalesced` arm,
  > which stops merging across unchanged bases on axes with no reading
- *(normalize)* coalesce protein cis adjacency created by the 3'-shift, and add the corpus axis that found it ([#1614](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1614))
  > 171 of 85642 corpus rows move (0.2%), all in
  > `protein_cis_separated`. The forms: a two-member protein bracket allele
  > whose members became adjacent under the 3'-shift,
  > `p.[Gly16Ala;Gly17del]`, becomes the merged delins
  > `p.Gly16_Gly17delinsAla`. Direction: TOWARD the already-shipped form,
  > not away from it — those 171 rows were not fixed points, so a second
  > normalization pass already produced the merged form; the change moves
  > pass one onto the answer pass two was giving, and toward the spec's
  > stated form, since `protein/substitution.md:32` marks the split spelling
  > `class="invalid"` (the decided
  > `delins-adjacent-members-when-both-consume-reference` record scopes that
  > to members which both consume reference bases, as a substitution and a
  > single-residue deletion do). 0 of the moved rows were previously
  > accepted forms: all 171 were not fixed points, so nothing downstream
  > holds them and the move is free. Measured by dumping this branch before
  > and after the fix over the identical 85,642-row corpus.
- *(normalize)* canonicalize a member whose span crosses a CDS boundary ([#1651](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1651))
  > 1826 of 78298 corpus rows move (2.3%) — 1747
  > respell and 79 change member count. Re-measured on this branch rebased
  > onto bfcc1802, not carried over from the base it was opened against.
  > 1460 of the moved rows were previously accepted forms; the other 366
  > were not fixed points, so moving those is free. Every move is a
  > straddling member being re-typed or affix-trimmed against the bases it
  > already denoted; no move changes what a description denotes.
- *(equivalence)* expose a groupable SPDI key for bucketing by denoted bases ([#1609](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1609))
  > SPDI keys move on a soft-masked reference.
  > `canonical_spdi`'s emitted `deletion` and `insertion` are now folded
  > through `apply_alphabet`, so a lowercase repeat-masked contig keys
  > `2:ATTAC:GGCTA` where it previously keyed `2:attac:GGCTA`, and an
  > `r.`-axis provider serving uracil now folds to the RefSeq DNA spelling
  > wherever the fold reaches (see above for the three shapes it does not,
  > none reachable from a RefSeq-spelled provider). No normalized HGVS
  > string changes, and no key derived from an uppercase reference changes.
  > The affected consumer is any that keys off a soft-masked FASTA, for whom
  > the previous keys split buckets that denote identical bases — so this is
  > a fix to the type's contract, migrated by re-deriving keys once.
- *(normalize)* split an equal-length protein delins whose interior is unchanged ([#1606](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1606))
  > Protein axis only — an equal-length `p.` delins
  > whose interior residue is unchanged now splits into its members, a split
  > and never a merge. Measured blast radius is zero on all three corpora
  > available, each for a stated reason (see Blast radius); the gate cannot
  > fire against the prepared reference at all, since it carries no protein
  > FASTAs.
- *(normalize)* gate the coding delins carve-out on the amino-acid precondition ([#1599](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1599))
  > 1581 of 78298 corpus rows move, all merge ->
  > split; 124 previously accepted (a real migration), 1457 previously not a
  > fixed point. On cis_confluence_axis both census pins are re-blessed
  > upward (3' 8026->8027, 5' 8021->8023). On the harsher
  > spec_conformance_axis corpus the move is mixed and is re-blessed in both
  > directions here: 3' non_idempotent_outputs 7->4 (three rows fixed), 3'
  > converged 9140->9139 (six classes lose convergence against five gaining
  > it, all six the spanning-delins respelling, promoted to
  > spec_corpus_regressions::the_codon_gate_splits_a_spanning_delins_its_own_members_do_not),
  > 3' split_two 2435->2436, 5' split_two 2696->2698 and split_three
  > 204->202. Three mutalyzer conformance rows leave agreement and carry a
  > spec_citation.
- *(rulings)* decide the exon-junction duplication and codon carve-out records ([#1623](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1623))
  > 0 rows move in this PR — no ferro output changes
  > here, and the corpus is untouched. What is declared is the DECISION:
  > `exon-junction-dup-converge-from-the-far-side` is now decided for
  > convergence, so every stored `c.3922dup`-shaped description will change
  > to `c.3921dup` and its projected genomic form will shift 2,790 bp once
  > #1621 implements it; and `codon-carve-out-shape-restriction` is decided
  > for widening, whose implementation #1599 measures at 1,581 of 78,298
  > corpus rows. Recorded here so both decisions are disclosed in the
  > release they were made in, not only in the release that ships the moves.
- *(project)* describe the protein consequence of a net-frameshift cis allele ([#1590](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1590))
  > 62 of the 402 protein-projectable `:c.` rows in
  > the committed multi-member cis fixture (15.4%) move off `p.?` onto a
  > concrete consequence — 59 onto a `p.(…fs…)`, and 3 onto the nonsense
  > substitution `p.(…Ter)` that the immediate-stop rule requires in place
  > of an illegal `fsTer1`. Nothing moves the other way and no
  > already-concrete form changes, so this is a gain of information, not a
  > re-spelling. Those same 3 rows are the only ones whose `is_frameshift`
  > moves (true -> false): the flag is now read off the built consequence,
  > so it agrees with the single-member spelling of the same sequence change
  > instead of with member arity. One further population is reachable but
  > empty on this corpus (0 rows): a net-in-frame combined product on a
  > reference CDS with no terminal stop, previously `p.?` and now described
  > — again matching what the single-member spelling has always produced.
  > The g./c./n. axes and `normalize` output are untouched.
- *(normalize)* demote a tract-wide repeat that a reducing delins grew ([#1602](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1602))
  > the wrong-sequence outputs this fixes move — each
  > previously denoted different bases than its input, so there is no
  > correct stored string to migrate from. The committed corpus measures 0
  > changed of 78,298 rows compared, which is a STRUCTURAL zero for the
  > reasons given above, and not evidence of safety.
- *(normalize)* bound the junction of a member that reduces to an insertion ([#1598](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1598))
  > 371 rows move of 5,629,002 swept; all 371 were
  > previously incorrect — 321 denoted the wrong bases outright and 50
  > denoted no sequence at all because their members overlapped. This is a
  > correctness repair, not a respelling: no row that already denoted its
  > input's bases changes.
- *(error-handling)* only read codon-aligned windows in the amino-acid case detector ([#1594](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1594))
  > 82 previously-rejected inputs now parse; no
  > normalized output moves. Measured over 368,870 ClinVar protein
  > descriptions in both strict and lenient mode — nothing previously
  > accepted is now rejected and no previously-accepted output changed.
- *(parser)* validate ring segments and sup inners in every post-parse check via a shared leaf walker ([#1578](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1578)) ([#1583](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1583))
  > zero normalized outputs move; previously-accepted
  > invalid ring/sup spellings across six validators now reject,
  > lenient/silent modes reject a ring del<size> they previously passed
  > verbatim, and member/segment repair-suggestion text changes.
  >
  > **NB** Same disclaimer as #1576 that I'm new to both Rust & the
  > material. This all "feels right" but I could be off base in a few
  > different dimensions.
- *(normalize)* widen the axis gate to c./n./r. ([#1484](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1484))
  > 2 027 rows move of 5 761 302 (0.035 %), measured
  > against `origin/main` with both sides built locally and row counts
  > validated; 2 016 gain a member, 0 lose, 11 same arity. Arity change, not
  > respelling, so a consumer keying on the normalized string sees the shape
  > of their data change — a real migration. By axis: 1 827 `c.`, 126 `n.`,
  > 74 `g.` — the `g.` rows are pre-existing `general.md:34` violations the
  > axis gate never touched, fixed by the axis-neutral member audit.
  > Confluence improves: 1 373 classes of 11 272 from the gate, plus 23 (3')
  > / 18 (5') from the member audit, with `declined`, `underdetermined` and
  > `sequence_changed` all 0. Separation violations fall by a net 874 (−14.7
  > %) before the member audit and further after it. One genuine
  > canonical-form move in the tail (`NM_000342.4`, 5-member → 4-member);
  > the other five tail rows are confluence fixes.

### Other

- bound each archive run, and name the phase a stall stopped in ([#2010](https://github.com/fulcrumgenomics/ferro-hgvs/pull/2010))
- *(normalize)* narrow the alignment DAG's storage to the band ([#2006](https://github.com/fulcrumgenomics/ferro-hgvs/pull/2006))
- *(ci)* take the last two heavy censuses off the debug shards ([#2001](https://github.com/fulcrumgenomics/ferro-hgvs/pull/2001))
- *(rulings)* correct the axis record's enforcement pin from 8 to 3 ([#1999](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1999))
- *(rulings)* re-ground the minimal-alignment rule as a house choice ([#1980](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1980))
- *(rulings)* correct the four equal-length rows' ground in the weight-bound record ([#1956](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1956))
- *(ci)* re-measure the denoted-sequence oracle; #1690 is no longer a blocker ([#1995](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1995))
- *(normalize)* denote each row under its own accession, not a hardcoded one ([#1997](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1997))
- *(scripts)* stop reading a fenced trailer as a representation declaration ([#1984](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1984))
- *(normalize)* five allocation and formatting wins on the normalization hot path ([#1981](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1981))
- *(spec)* bump hgvs-nomenclature to 565b973 and repair the moved citations ([#1793](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1793))
- *(corpus)* derive the corpus header's family counts instead of restating them ([#1976](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1976))
- *(ci)* build the soak archive from a driver package, not the whole `it` crate ([#1977](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1977))
- *(normalize)* band the sequence-first alignment DP ([#1973](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1973))
- *(ledger)* qualify what the brute-force dominator oracle corroborates ([#1985](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1985))
- [**breaking**] remove unused convert API (cds_to_transcript_pos, strand, … ([#1681](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1681))
- audit each normalization stage as governed, ruled, or undecided ([#1855](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1855))
- *(ci)* move the three slow censuses onto the optimized archive ([#1960](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1960))
- *(reference)* derive the intronic band predicates from the one splice ladder ([#1825](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1825))
- add sequence_normalize() helper ([#1965](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1965))
- *(normalize)* narrow the alignment cost grid to u16, and make divergence a benchmark axis ([#1952](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1952))
- *(spdi)* pin the flush-insertion rule against a deletion and an inversion ([#1968](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1968))
- *(rulings)* give the ledger a home for a clause-free house choice ([#1953](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1953))
- *(ledger)* record the measured scope of the minimal-alignment ruling ([#1958](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1958))
- *(python)* report the #1455 GIL margin on every run, and retry once more ([#1898](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1898))
- give the five real-data normalization guards assertions that can fail ([#1867](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1867))
- *(normalize)* record the minimal-alignment cost model gap, and correct an alignment count ([#1949](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1949))
- enumerate minimal alignments, so the unchanged-base rule can be computed ([#1944](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1944))
- *(corpus)* straddle the gate that fires, and import the bound instead of restating it ([#1938](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1938))
- *(prepare)* make prepare re-runnable into a directory it prepared ([#1962](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1962))
- *(normalize)* cite #1899 in the weight-bound ruling record ([#1920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1920))
- state the "assert the property, never let a count be the property" doctrine ([#1927](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1927))
- withdraw a stale flip instruction and correct a backwards record citation ([#1951](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1951))
- *(tests)* record that one stray symbol reclassifies the whole inserted run ([#1936](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1936))
- *(readme)* equal block length does not make the column correspondence unique ([#1937](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1937))
- *(ci)* delete the dead nextest ci profile that killed the censuses ([#1948](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1948))
- label the coordinate-basis docs by layer, one splice ladder per question ([#1742](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1742))
- *(rulings)* judge the guard scan on its reach, not on what the tree holds ([#1945](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1945))
- make the applier's overlap verdict order-independent, and let reanchor widen a pure-deletion window ([#1837](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1837))
- *(python)* validate the numeric conversion entry points at the boundary ([#1741](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1741))
- make a present-but-incomplete reference fail, and delete a test whose corpus never existed ([#1861](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1861))
- *(rulings)* make a record state what enforces it, and reject silence ([#1883](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1883))
- *(reference)* decline an intron distance the offset cannot measure ([#1918](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1918))
- *(data)* populate MappingInfo::transcript_id on the g. to c. path ([#1902](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1902))
- *(changelog)* correct #1835's flip disclosure, and make a merged one correctable ([#1897](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1897))
- *(normalize)* measure the canonical-coalesced flip over the four real corpora ([#1900](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1900))
- *(normalize)* [**breaking**] delete the input-relative weight bound and derive the split cap ([#1899](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1899))
- Add more test coverage based on mutatino testing ([#1676](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1676))
- *(rulings)* check that a record's guard citation resolves ([#1880](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1880))
- watch src/reference/ and src/error_handling/ for representation changes ([#1876](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1876))
- *(normalize)* relabel #1421's six rows SpecExplicit on delins.md:17 ([#1869](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1869))
- *(benchmark)* repair a hollow supplemental record from the FASTA, and refuse a corrupt one ([#1732](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1732))
- *(data)* refuse an inverted-CDS record instead of underflowing on it ([#1877](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1877))
- *(rulings)* scope two ruling claims to what was actually measured ([#1838](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1838))
- *(normalize)* [**breaking**] remove the 5' shuffle direction from the public surface ([#1879](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1879))
- *(ci)* refuse a Representation-Change trailer no consumer can read ([#1852](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1852))
- run the 45 orphaned reference-aware guards, and fail on a missing manifest ([#1851](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1851))
- publish the adjudication ledger as a generated normative document ([#1850](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1850))
- *(normalize)* name the partition rule the derivation surface cuts with ([#1848](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1848))
- *(effect)* [**breaking**] decline an unknown offset from the two splice classifiers that could not ([#1843](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1843))
- *(spec)* measure spec conformance per FERRO_PARTITION arm over a real reference ([#1842](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1842))
- *(examples)* refuse a FERRO_PARTITION naming no arm from every cargo target ([#1822](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1822))
- *(effect)* decline the unknown-offset sentinel instead of classifying it as a splice site ([#1806](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1806))
- *(project)* record why the g. axis declines on a bare transcript input ([#1812](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1812))
- *(hgvs)* refuse a `+`/`-` offset on a `g.`/`m.`/`o.` position ([#1792](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1792))
- *(normalize)* drop the subsumed separation disjunct from the inversion gate ([#1787](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1787))
- *(convert)* pin that the sequence and genome frames disagree on a gapped transcript ([#1776](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1776))
- *(rulings)* decide what delins.md:47 reaches — an inserted sequence, not any coincidence ([#1760](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1760))
- *(conformance)* a codon-designed protein corpus, and the protein axis's first census ([#1818](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1818))
- *(harvest)* un-ignore #1541's convergence guard and re-ground its form pin ([#1817](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1817))
- *(conformance)* tighten the protein empty-projection budget to what main produces ([#1803](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1803))
- *(normalize)* pin the codon exception at a non-zero window offset ([#1810](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1810))
- enumerate each cis adjudication class's description space ([#1814](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1814))
- *(cli)* report a declined description in hgvs-to-vcf and keep going ([#1811](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1811))
- record the selection-wide FERRO_ASSERT_SEQUENCE measurement test-oracle asked for ([#1813](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1813))
- *(equivalence)* drive the cross-axis rung over constructed forms, not a diverging corpus pair ([#1809](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1809))
- *(rulings)* duplication.md:18 ranks the dup label, not the partition ([#1801](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1801))
- make the reported pair model a closed three-state enumeration ([#1805](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1805))
- *(data)* read cdot's raw genomic cds_end as exclusive, not inclusive ([#1786](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1786))
- *(benchmark)* read internal cdot transcript coordinates as 0-based ([#1740](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1740))
- *(test)* parallelize the three remaining slow censuses, and cut the per-row cost of the two biggest ([#1782](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1782))
- *(cdot)* propagate cds_start_incomplete through Transcript to CdotTranscript ([#1770](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1770))
- *(benchmark)* make the supplemental FASTA the single authority for presence, and upsert by accession ([#1791](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1791))
- *(projector)* pin the non-coding axis non-idempotency residue ([#1788](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1788))
- *(conformance)* count coding-axis merges across two or more unchanged nucleotides ([#1710](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1710))
- *(convert)* assert the cdot 0-based tx basis without a prepared reference ([#1779](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1779))
- correct the #466 filing-status claim in a decided ruling record ([#1784](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1784))
- *(hgvs)* refuse `n.*N` at parse in every mode; `n.-N` stays strict-only ([#1751](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1751))
- *(vcf)* refuse an unresolvable genomic position instead of dropping it ([#1734](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1734))
- generate and publish the n./r. cis confluence corpus, and guard the wiring ([#1775](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1775))
- *(conformance)* refuse a zero denominator on a mutalyzer axis, and measure the corpus by member arity ([#1718](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1718))
- *(normalize)* scope the compensating-gap coalesce to the coding DNA axis ([#1727](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1727))
- *(rulings)* record that a derivation may not be bounded by the input's spelling ([#1733](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1733))
- *(rulings)* correct a false cost claim in the payload-coincidence axis-scope record ([#1696](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1696))
- host bulk test fixtures as release assets instead of Git LFS ([#1750](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1750))
- *(rulings)* record general.md:34 as a preference, not a prohibition ([#1725](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1725))
- *(normalize)* make each frame-derived rule declare its axis scope ([#1688](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1688))
- *(normalize)* pin that normalization is axis-preserving and frame-free axes ignore annotation ([#1686](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1686))
- *(rulings)* refuse a ruling that speaks to a molecule it does not cite ([#1685](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1685))
- *(ci)* make the local seam-oracle run reproduce CI's selection ([#1689](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1689))
- *(equivalence)* a cross-axis rung, and an undecidable one ([#1652](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1652))
- *(python)* expose from_sequences, to_sequences and re-anchoring ([#1695](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1695))
- *(normalize)* re-anchor a window to a bound the variant must not leave ([#1694](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1694))
- *(normalize)* derive a conformant HGVS description from a reference/alternate pair ([#1693](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1693))
- *(ci)* derive the reporter's provenance sentence from the triggering event ([#1705](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1705))
- *(error-handling)* [**breaking**] remove the unwireable mode-level warning predicates, and stamp the mode censuses are measured under ([#1674](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1674))
- *(fixtures)* restore transcript exon ends and guard the Exon contract ([#1673](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1673))
- *(normalize)* pin the stranded identity member, measured and non-confluent ([#1668](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1668))
- *(spdi)* refuse a genomic pter/qter/cen instead of flattening it to base 0 ([#1662](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1662))
- *(readme)* state that rule 3 governs rule 2's evaluation basis ([#1660](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1660))
- *(normalize)* refuse an extreme coordinate instead of wrapping the window ([#1658](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1658))
- *(rulings)* two adjudication records, and make four census ratchets state their trap ([#1648](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1648))
- *(rulings)* scope the delins merge ruling by direction, close the ins carve-out ([#1680](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1680))
- describe the citation-quote guard as the substring check it is ([#1683](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1683))
- *(changelog)* attach each Representation-Change trailer to its changelog entry ([#1666](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1666))
- *(parser)* say what parse_hgvs actually applies — the grammar, and no ErrorConfig ([#1667](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1667))
- *(rulings)* decide the absolute-prohibition enforcement stage as mode-dependent ([#1634](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1634))
- *(ci)* let the failure reporter resolve the repository it files against ([#1663](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1663))
- enforce the adjudication tables against the ledger, and correct three records ([#1671](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1671))
- *(ci)* read the count, not the denominator, in a quantified-zero trailer ([#1657](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1657))
- *(cis)* correct the confluence headline and check it against the pin ([#1653](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1653))
- *(protein)* assert the reference stop decodes to `AminoAcid::Ter` ([#1656](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1656))
- *(normalize)* make a partitioner decline observable, and stop the bake-off switch aborting a shipping process ([#1639](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1639))
- *(normalize)* a fourth seam oracle that compares the denoted sequence ([#1622](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1622))
- *(cis)* census the n. and r. confluence axes, and pin what #1484 did there ([#1646](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1646))
- *(normalize)* merge a split whose boundaries the alignment manufactured ([#1644](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1644))
- *(normalize)* let the inversion coalesce reach a sequence-derived partition ([#1640](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1640))
- *(adjudication)* correct two stale claims and pin what two respell tests denote ([#1645](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1645))
- *(conformance)* pin what the output denotes, not only how it is spelled ([#1633](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1633))
- *(corpus)* verify a row against its reference sequence, not its label ([#1625](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1625))
- *(spdi)* refuse a genomic offset instead of dropping it ([#1641](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1641))
- *(rulings)* decide the separation record for the re-derived form ([#1620](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1620))
- *(parser)* require ring `::` segments to be ordered deletions ([#1595](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1595))
- *(rulings)* resolve the ledger's own record-to-record citations ([#1611](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1611))
- *(rulings)* make the precedence record a pointer to the README ruleset ([#1604](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1604))
- *(conformance)* an exhaustive spec-derived corpus, and the burn-down it measures ([#1585](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1585))
- *(normalize)* pin the two worked examples that explain the weight bound ([#1591](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1591))
- [**breaking**] correct three stale claims, and disambiguate the two NormalizeConfig types ([#1596](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1596))
- *(readme)* state the project's normalization ruleset ([#1605](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1605))
- *(cis)* pin the converged rows, not just their count, on the cis axis census ([#1581](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1581))
- *(normalize)* let normalizer warnings reach the caller ([#1580](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1580))
- *(normalize)* refuse an unknown FERRO_PARTITION instead of serving live ([#1582](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1582))
- *(equivalence)* decline a ring no rung can evaluate; rule that self-cancellation does not cross `::` (#1578 items 3-4) ([#1589](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1589))
- *(inversion)* gate 2,075 authored inversions on real bases, hermetically ([#1554](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1554))
- *(changelog)* make the audit agree with the checker, and anchor its range to the merge ref ([#1593](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1593))
- *(changelog)* audit the Representation changes section against its trailers ([#1560](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1560))
- *(parser)* cover the last three allele-recursion arms ([#1587](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1587))
- Adding tests from mutation testing findings ([#1576](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1576))
- *(changelog)* read a decline as a decline when text follows the trailer ([#1574](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1574))
- *(conformance)* pin the adjudications that keep getting re-litigated ([#1579](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1579))
- *(ci)* 62% off the critical path on a re-run, 43% off runner time, cache 13.8→7.0 GiB ([#1564](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1564))
- *(changelog)* state what v0.13.1 moved, in the released section ([#1561](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1561))
- *(conformance)* refuse to write an artifact built from a partial generator pass ([#1551](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1551))

### Added

- **`from_sequences` — derive an HGVS description from a reference/alternate pair.** For callers
  whose input is bases rather than a description: a window out of a BAM, a VCF row, an aligner's
  output. It reads no reference sequence, so the output is a pure function of its arguments,
  and it delivers the two normalization rules that are always achievable — **conformant** and
  **deterministic** — leaving **recommended form** and **confluent** to `normalize`, which holds
  the reference. Available as `ferro_hgvs::from_sequences` / `from_sequences_detailed` in Rust and
  `ferro_hgvs.from_sequences` / `from_sequences_detailed` in Python.

  It derives `g.` descriptions, and `m.` on the two rCRS mitochondrial accessions
  (`NC_012920`, `NC_001807`), which HGVS requires the `m.` axis for. Every other accession class is
  refused with a message naming the class.

  Because the derivation never sees a spelling, two spellings of one variant reach one
  description: 5 636 cis confluence corpus classes with no divergence — the corpus's **genomic
  half**, which is all a `g.`-only surface can reach, its 5 636 `c.` classes being drawn against a
  transcript accession this refuses by design — and nine of the nine externally-reported
  #1419 / #1420 / #1421 confluence pairs in both shuffle directions.

  Nothing here changes any existing output — `normalize` is untouched, and a caller who has a
  description and wants it normalized still needs it to converge. This is a new entry point for a
  different input.

- **`Normalizer::to_sequences`** — the inverse, turning any HGVS description into that pair, so a
  caller already holding descriptions needs no new plumbing to reach the derivation. Pads **both**
  sides by `pad` — so the window is `span + 2 * pad` — because `dup` typing reads the reference
  bases immediately 5' of an insertion point (`DNA/duplication.md:18`). The window is returned
  upper-cased, so a soft-masked region does not come back as a mixed-case pair.

- **`Normalizer::from_sequences`** — the same derivation against a held reference, which lets it
  additionally range-check the interval and optionally `normalize` the result.

- **`Normalizer::sequence_normalize`** — the "one canonical description per variant" round trip, as
  a single call. It takes a parsed description to the bases it denotes (`to_sequences`), re-derives
  a description from those bases alone (`from_sequences`), and — while a member still rests on a
  window edge that can still move — doubles the pad and retries, so two spellings of one variant
  reach one description decided by the observed bases rather than by how either was written. The
  widening loop reads the two per-side flags apart, so a placement pinned to the sequence's own 5'
  or 3' terminus is recognised as settled rather than chased; an interior tract wider than the
  widest window it will read is declined, naming the unsettled side, rather than answered over a
  truncated window. Available as `Normalizer::sequence_normalize` in Rust and
  `Normalizer.sequence_normalize` in Python.

- **Re-anchoring, for callers whose variant must stay inside a region.** `SequencePair::trim_to`
  narrows a window to given bounds, trimming matching bases only and needing no reference;
  `Normalizer::reanchor` moves both edges, padding from the reference where it must widen. Bounds
  are 1-based inclusive and optional per edge. Both refuse rather than clamp — including a bound
  outside the sequence, since a window silently pulled back to the contig would hide a bug upstream
  of the call.

  **`reanchor` moves a window's edges; it does not relocate the window.** The window asked for must
  overlap the pair's own, because the changed bases exist only in the pair — a disjoint request is
  refused, with a message naming `reanchor` and both windows. Its bases come back upper-cased, as
  `to_sequences`' do, so widening a soft-masked window cannot splice provider bases onto caller
  bases and return a mixed-case pair; `trim_to` fetches nothing and leaves case alone. Case is not
  a disagreement in either: a soft-masked reference against an upper-case alternate trims normally.

  This is for a bound that is a *requirement* (a target region, an amplicon, a tiling window). It
  is not the way to make heterogeneous raw pairs agree in general: `normalize = true` and a
  `to_sequences` round trip both already do that and reach the reference-anchored placement, which
  can shift further than any fixed window allows.

- **`SequencePair::new`**, so a caller holding bases and no description can build one. The type is
  `#[non_exhaustive]` and was previously only ever returned, which put both re-anchoring entry
  points out of reach of exactly the caller they are for. It validates through the same check
  `from_sequences` uses rather than a second copy, so a pair that constructs is a pair that
  derives.

- Both sequences are **case-folded** and the alphabet is DNA. A soft-masked (lower-case) window
  derives exactly as its upper-case twin does, which keeps `to_sequences` -> `from_sequences` an
  inverse over a masked region; `U` is refused on either side rather than admitted into a `g.`/`m.`
  description.

- `FromSequencesOptions::with_direction` and `with_max_grid_cells`. The struct is
  `#[non_exhaustive]`, so a struct expression is forbidden outside the crate and
  `Default::default()` was the only reachable constructor — the two documented knobs could not be
  set by any downstream caller.

### Fixed

- **Reference-span insertions are now resolved rather than refused.** An insertion whose bases are
  named by position in the same reference — `g.100_101ins50_57`, and its inverted form
  `g.100_101ins50_57inv` — was rejected as "insertion sequence is not a literal sequence" whenever
  a provider-backed path (`sequence_normalize`, `to_spdi`, `canonical_spdi`, apply-to-reference)
  needed its bases. It now reads those bases from the reference exactly as `del`/`dup`/`delins`
  read their omitted bases, so `sequence_normalize` re-derives the literal `g.100_101insCGTA…`
  form. `delins` with a same-reference range insert is resolved on the same path. The payload is
  read on the **axis the description is written on**: a `c.` (and a coding `r.`) payload gets the
  same `cds_start + N - 1` shift the location already gets, so `NM_000532.5:c.156_157ins180_188`
  resolves to `insGCGAGGAAA` — the bases `normalize` names — rather than to the transcript-offset
  span 36 bases 5' of it. `n.`/`g.`/`m.`/`o.` payloads are transcript/genomic offsets already and
  are read unshifted.
- **Exact tandem-repeat insertions are now expanded rather than refused.** An insertion whose bases
  are named by a spelled unit and an exact copy count — `g.100_101insC[4]`, `g.100_101insCG[3]` —
  was rejected as "not a literal sequence" whenever a provider-backed path needed its bases. It now
  expands the unit the named number of times, exactly as the short-form repeat edit does, so
  `sequence_normalize` re-derives the literal `g.100_101insCCCC` / `g.100_101insCGCGCG` form (the
  same `MAX_REPEAT_EXPANSION_BASES` cap applies). `delins` with an exact repeat insert is resolved
  on the same path. Uncertain or range counts (`insC[10_15]`, `insC[?]`) remain undetermined and
  still decline.
- **Compound-insert brackets mixing shapes are now resolved rather than refused.** A bracketed
  insertion whose parts combine a literal with a same-reference coordinate span or an exact tandem
  repeat — `g.185_201delins[T;213_271]`, `g.100_101ins[C[3];A]`, `g.100_101ins[213_271inv;A]` —
  was rejected as "Unsupported variant type … cannot apply to reference" whenever a provider-backed
  path needed its bases, because only single-part brackets were read out. Each part is now resolved
  through the same leaf that resolves it on its own — a span reads from the reference (and
  reverse-complements for the `…inv` form), an exact repeat expands its unit — and the parts are
  concatenated in written order into one contiguous inserted sequence, so a compound insert encodes
  to SPDI just as a plain literal does. This applies on both the `ins` and `delins` arms. A part
  whose bases are genuinely undetermined — an uncertain or range repeat count, a CDS position range
  carrying intronic offsets, or an external reference SPDI cannot dereference — still declines the
  whole insert.
- **A `from_sequences` derivation whose insertion rolls to the contig's first base is now
  re-anchored rather than refused.** When the observed window begins at position 1 — a variant
  anchored at `1A>C`/`1A>T`, say — the 5'-shuffle can roll a pure insertion onto interbase 0, whose
  only HGVS spelling names position 0, a base that does not exist. Because the window already starts
  at the accession's first base, `sequence_normalize` cannot widen the 5' flank to escape it, so
  these derivations were dropped with "the derivation places an inserted payload immediately 5' of
  the window's first base". There are now two escapes, one per shape:
  - **The insertion crosses the terminal base unchanged** (an ambiguous run reaching the terminus,
    base 1 itself not changed): it is presented at the leftmost *nameable* interbase — `1_2ins…`,
    the payload rotated one base — the terminal analogue of the 3'-most rule.
  - **Base 1 is itself rewritten** (a dense substitution stack, the primer/adapter-artifact shape):
    the leading insertion is folded rightward into a single span anchored at base 1, so
    `g.[1A>C;2T>A;3G>T]` derives to `g.1_3inv` (the span is a reverse complement) and
    `g.[1G>C;2T>A;3G>T]` to `g.1_3delinsCAT`. This is the span form Mutalyzer emits for the same
    inputs (`g.1_5delinsCAAGA`, `g.1_3inv`), reached without needing an interbase 5' of position 1.

  A pure insertion genuinely *before* base 1 with no piece to fold into still refuses, since it
  names a sequence no in-window placement does.

## [0.13.1](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.13.0...v0.13.1) - 2026-08-08

### Representation changes

**3 stored strings move — 3 of 500,004 ClinVar rows (0.0006%).** Two are respelled
and one merges back to the form the submitter wrote. These inputs were **previously
accepted**, so a consumer who normalized them under v0.13.0 holds the old string:
a real migration, if a very small one. The merge, end to end:

```
NM_000083.3:c.2461_2464delinsCTCC
  v0.13.0 ->  c.[2461_2463delinsCTC;2464G>C]
  v0.13.1 ->  NM_000083.3:c.2461_2464delinsCTCC   (unchanged — the submitted spelling)
```

Both entries below are the same fix in different directions, and only the first
moves anything. [#1535] moves **0 rows over 5,761,302 real expressions** (ClinVar
500,004 / Paraphase 435,235 / CMRG 4,826,063), and that zero is *structural* rather
than reassuring: those corpora contain no instance of the shape it converts, so the
measurement could not have found a move. Where the shape does occur the split
spelling moves onto `inv`.

One deliberate cost, recorded because it is the property this project ranks above
stability: **confluence drops by 6 classes of 11,272** in [#1537], by operator
ruling. Those classes agreed only by way of a spec-illegal intermediate the fix
removes; `sequence_changed` stays 0, so nothing converged on a wrong sequence, and
the surviving form is the better-formed one. What was lost is agreement, not
correctness. The same change fixes 58 of 59 known `delins.md:16` separation
violations.

[#1535]: https://github.com/fulcrumgenomics/ferro-hgvs/pull/1535
[#1537]: https://github.com/fulcrumgenomics/ferro-hgvs/pull/1537

- *(normalize)* never split a delins into members on consecutive nucleotides ([#1537](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1537))
  > 3 rows of 500,004 move (0.0006%) — 2 respell, 1
  > merge toward the submitted form. Confluence drops by 6 classes of
  > 11,272, deliberately and by operator ruling: those classes agreed only
  > through a spec-illegal intermediate, `sequence_changed` stays 0, and the
  > surviving form is the better-formed one. 58 separation-0 violations of
  > `delins.md:16` are fixed.
- *(normalize)* type the whole-block inversion by what it competes with ([#1535](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1535))
  > 0 rows move over 5,761,302 real expressions
  > (ClinVar 500,004 / Paraphase 435,235 / CMRG 4,826,063), diffed row by
  > row between this branch and its base on both binaries. That zero is
  > structural for the direction that could move: a cis-allele of >=2
  > `delins` members forming a whole-block reverse complement has **0**
  > instances in those corpora, so the corpus cannot build the input this
  > change converts. It is not vacuous either — 2,639 `inv` rows (897 of
  > them genomic) and 38,190 multi-base spanning `delins` rows do reach the
  > typing path and are byte-identical across the change, so the gate
  > demonstrably does not over-fire. Where the shape does occur, the split
  > spelling moves onto `inv`; the `inv` and spanning-`delins` spellings do
  > not move on `c.`, and on a frameless axis all three converge on `inv`.

### Other

- *(changelog)* read a decline that gives a reason as a decline ([#1558](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1558))
- *(confluence)* record two adjudications from the divergent cis classes ([#1549](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1549))
- *(normalize)* guard that emitted members obey general.md:34 ([#1539](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1539)) ([#1540](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1540))
- *(spec)* run the spec's own worked examples against real bases, and record #1536 ([#1544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1544))
- *(conformance)* gate the harvested unguarded cases on every PR ([#1545](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1545))
- *(normalize)* validate a cis-allele member's stated reference bases before a merge can consume them ([#1547](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1547))
- *(rulings)* decide the precedence order and the canonical-form choice ([#1548](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1548))
- *(corpus)* add a separated reverse-complement family and three real oracles ([#1546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1546))
- *(rulings)* restore the delins.md:47 decision dropped by #1532's squash ([#1538](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1538))
- *(spec-fixture)* close the review-sweep follow-ups from #1519/#1520/#1521 ([#1533](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1533))
- record every adjudication as a test in the same change ([#1531](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1531))
- record the unrecorded adjudications and correct two ruling rationales ([#1532](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1532))
- *(claude)* document the required Representation-Change trailer ([#1534](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1534))
- *(changelog)* stop grouping a declined change as a representation change ([#1527](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1527))
- require a Representation-Change declaration on output-bearing changes ([#1523](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1523))

## [0.13.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.12.0...v0.13.0) - 2026-08-07

### Representation changes

Normalization output is a shipped key: a downstream consumer stores read counts
against the normalized HGVS string, so a form that moves between releases
re-buckets what they have already stored.

No commit this cycle carried a `Representation-Change:` trailer, so this section
was reconstructed after the fact by normalizing three real corpora through
v0.12.0 and v0.13.0 against the same prepared reference and diffing row by row —
5,761,302 expressions (ClinVar 500k, Paraphase, CMRG).

| outcome | rows | what it means for a consumer |
|---|---|---|
| identical on both releases | 4,490,581 | nothing to do |
| **moved** | **577** (0.0100%) | **re-map these keys** |
| newly normalizes | 326,404 | free — no key existed before |
| now rejected | 2 | key disappears; the input was malformed |
| failed on both releases | 943,738 | no key on either release |

**577 stored strings move — about 1 in 10,000.** The direction is mixed: 360
collapse to fewer allele members, 205 expand into more, and 12 are respelled at
the same arity. Two opposing mechanisms are responsible — the `delins.md:44-47`
payload-alignment merge ([#1470]) pulls toward one delins, and the sequence-first
cis derivation ([#1392]) pulls toward split members. Per corpus: ClinVar 23,
Paraphase 36, CMRG 518.

**326,404 inputs newly normalize, and every one is an LRG accession** — [#1501]
(LRG accessions are versionless, so W3001 is no longer demanded) and [#1490]
(LRG genomic bases are served from the LRG record). These previously returned a
parse error, so no consumer holds a stored string for them and adopting them
costs nothing.

**2 inputs are now rejected rather than echoed back**, both the same class: two
coincident cis-allele edits, for which the spec defines no canonical form
(`OverlapConflictingEdits` / W5002). `NM_033028.4:r.[118_261del;118_373del]` and
`NM_003803.3:c.[64_66delGTGinsCTC;65T>G]`. v0.12.0 accepted both silently and
returned a string that denoted nothing.

**The 943,738 that fail on both releases are unchanged in kind**, and are a
property of the corpora rather than of this release:

| reason | rows | share |
|---|---|---|
| intronic variant — declined by design | 860,832 | 91.2% |
| parse error (spelling ferro rejects) | 49,162 | 5.2% |
| reference mismatch (stated base differs from the reference) | 32,302 | 3.4% |
| invalid coordinates / spec-undefined | 1,442 | 0.2% |

The intronic decline is longstanding — the same error variant exists at
`src/error.rs:429` in v0.12.0 and v0.13.0 — and is not new here.

For contrast, the synthetic shape-family corpus (`dump_normalized_corpus`,
78,028 rows) moves 12,530 rows, 2,201 of them previously stable. That corpus is
deliberately enriched for the shapes that churn and its rate is not a
library-wide figure; the real-corpus numbers above are the ones that describe
consumer impact.

[#1392]: https://github.com/fulcrumgenomics/ferro-hgvs/pull/1392
[#1470]: https://github.com/fulcrumgenomics/ferro-hgvs/pull/1470
[#1490]: https://github.com/fulcrumgenomics/ferro-hgvs/pull/1490
[#1501]: https://github.com/fulcrumgenomics/ferro-hgvs/pull/1501

### Added

- *(dev)* add --verify-spdi to the normalized-corpus harness ([#1515](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1515))
- *(normalize)* record provenance when reported members are coalesced ([#1486](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1486))
- *(normalize)* coalesce payload-alignment splits back into one delins ([#1470](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1470))
- *(normalize)* add the member-count-minimal partitioner behind a flag ([#1464](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1464))
- *(spdi)* reference-aware SPDI and apply-to-reference for an encoding-invariant key ([#1374](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1374))
- *(normalize)* make the sequence-first derivation authoritative for cis alleles ([#1392](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1392))
- *(effect)* [**breaking**] drop the dead inversion validator that contradicted the live one ([#1359](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1359))

### Fixed

- *(normalize)* apply a junction barrier the translate cannot express ([#1516](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1516))
- *(spdi)* name the caller's anchor when a repeat names no tract ([#1496](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1496))
- *(normalize)* type a split-created span inv on the first pass ([#1493](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1493))
- *(normalize)* see cis members that span a region boundary ([#1506](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1506))
- *(reference)* length a record-less placed parent from coordinate one ([#1502](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1502))
- *(normalize)* report cis overlaps across a region boundary ([#1511](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1511))
- *(spdi)* normalise the reference window before searching a repeat tract ([#1489](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1489))
- *(reference)* serve LRG genomic bases from the LRG record, not the placement ([#1490](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1490))
- *(errors)* LRG accessions are versionless, so stop demanding W3001 ([#1501](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1501))
- *(normalize)* fold case when counting a tandem repeat tract ([#1495](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1495))
- *(normalize)* recognise a dense inversion across multi-base separations ([#1485](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1485))
- *(normalize)* keep the sequence-first derivation inside one exon ([#1480](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1480))
- *(normalize)* repair non-coding `r.` members instead of reverting them ([#1465](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1465))
- *(reference)* a zero-length sequence read is empty, not out of bounds ([#1479](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1479))
- *(normalize)* report span edits whose spans intersect without coinciding ([#1451](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1451))
- *(normalize)* key a dup on what it writes in the interior-junction branch ([#1446](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1446))
- *(normalize)* order a conflicting cis allele whose members were rewritten ([#1438](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1438))
- *(normalize)* keep a real overlap conflict visible in the output ([#1445](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1445))
- *(spdi)* read a single-position repeat anchor as the whole tract ([#1439](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1439))
- *(normalize)* cap the split only where there is a gap to place ([#1403](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1403))
- *(normalize)* demote a repeat that outgrew its tract, on both branches of the guard ([#1399](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1399))
- *(normalize)* shift derived pieces in the direction that was requested ([#1432](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1432))
- *(normalize)* let a 5'-shifting junction insertion finish in one pass ([#1428](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1428))
- *(normalize)* bound a sibling-crossing shift across a region boundary ([#1425](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1425))
- *(normalize)* let a junction insertion finish its shift in one pass ([#1427](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1427))
- *(normalize)* drop a cancelled identity member a sibling partly covers ([#1423](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1423))
- *(normalize)* clamp a derived insertion at the CDS/3'UTR junction ([#1417](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1417))
- *(normalize)* take the overlap verdict from the input, not the repair ([#1411](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1411))
- *(normalize)* refuse to lower a repeat past a sibling inside its tract ([#1407](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1407))
- *(normalize)* drop a cancelled cis member where the merge creates it ([#1380](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1380))
- *(normalize)* ask the codon gate about the window the repeat is emitted over ([#1393](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1393))
- *(normalize)* apply the strict-mode ladder on both public exits ([#1385](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1385))
- *(parse)* let a UTR-marker ins payload reach the message written for it ([#1377](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1377))
- *(spdi)* name the span and the way out when bases are unspelled ([#1390](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1390))
- *(normalize)* guard a net deletion's split by rule, not by a length bound ([#1373](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1373))
- *(normalize)* stop an out-of-phase rotation from killing a valid repeat ([#1368](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1368))
- *(normalize)* validate the reference base an identity member states ([#1367](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1367))
- *(spdi)* [**breaking**] give an identity edit the span it claims ([#1363](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1363))
- *(normalize)* stop a terminal duplication respelling past the contig end ([#1365](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1365))

### Other

- *(msto)* pin the reported family's expected outputs and link its guards ([#1520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1520))
- *(normalize)* correct the coincidence statistics behind the density gate ([#1521](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1521))
- *(spec)* assert equivalence classes converge and record contested spec readings ([#1519](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1519))
- compile the default-feature unit tests so a gated call cannot ship ([#1510](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1510))
- retire the `parallel` feature gate ([#1509](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1509))
- Repairs 'cargo test' ([#1504](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1504))
- document FERRO_PARTITION as an unstable evaluation switch ([#1483](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1483))
- *(repr)* add a multi-exon coding axis to the representation corpus ([#1481](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1481))
- *(python)* release the GIL across reference reads, and stop deep-copying the reference per call ([#1475](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1475))
- *(sweep)* exercise the non-coding `r.` axis in the transcript cis sweep ([#1476](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1476))
- pin the reported partition family's per-spelling verdicts and spec arguments ([#1474](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1474))
- *(cis)* measure confluence over designed multi-member cis alleles ([#1458](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1458))
- *(repr)* reach #1454's shape, and close the corpus over its output vocabulary ([#1468](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1468))
- *(repr)* draw a block past the split cap in the representation corpus ([#1467](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1467))
- *(corpus)* re-enumerate the msto issue catalog and let it see Python guards ([#1466](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1466))
- harvest authored spec cases from the issue tracker into runnable code ([#1457](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1457))
- *(repr)* build conflicting cis alleles in the representation corpus ([#1459](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1459))
- correct two normalizer comments that name the wrong mechanism ([#1447](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1447))
- *(repr)* commit the normalized-output dump/diff harness ([#1442](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1442))
- surface representation changes on the release PR and changelog ([#1433](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1433))
- describe batch workloads generically in binding and test comments ([#1434](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1434))
- *(cis)* sweep three-member alleles, with the causes counted apart ([#1402](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1402))
- *(cis)* sweep the transcript axes, and find a CDS-boundary defect ([#1400](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1400))
- *(ci)* pin every sweep_seeds consumer to SWEEP_FILTER ([#1410](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1410))
- *(reference)* cache decoded sequence blocks under a byte budget ([#1371](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1371))
- *(cis)* draw a corpus prefix by default, the full sweep in CI ([#1395](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1395))
- *(spdi)* count the UnspelledBases arms as six, not five ([#1412](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1412))
- *(spec)* accept either spec-legal stop spelling as conformant ([#1381](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1381))
- *(review)* stop implying three suites cover both normalization exits ([#1387](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1387))
- lint the release profile so the oracles' debug-only gating is enforced ([#1386](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1386))
- *(normalize)* repoint the circular-normalizer deferral at its live tracker ([#1375](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1375))
- stream the library and Python batch APIs in bounded chunks ([#1372](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1372))
- *(normalize)* assert at the seam that no coordinate runs off its sequence ([#1366](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1366))
- record that cargo build --features python cannot link ([#1364](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1364))
- *(normalize)* record the 5' shuffle confluence gap ([#1361](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1361))

## [0.12.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.11.0...v0.12.0) - 2026-08-03

### Added

- *(sequence)* [**breaking**] finish the complement unification, and fix an RNA inversion false alarm ([#1356](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1356))
- *(normalize)* sequence-first cis-allele splitter, and the measurements that bound it ([#1341](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1341))
- *(normalize)* add a tsv output format reporting which variants changed ([#1242](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1242))

### Fixed

- *(normalize)* bound a 5'-shifting member at a sibling's insertion junction ([#1357](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1357))
- *(normalize)* stop an out-of-phase rotation from killing a valid dup ([#1354](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1354))
- *(normalize)* bound a cis junction's 5' shift at its siblings ([#1350](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1350))
- *(extractor)* stop three delins shapes from classifying as inversions ([#1346](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1346))
- *(normalize)* stop refusing a derivation for disagreeing with its input ([#1345](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1345))
- *(normalize)* keep a junction re-spelling inside the sequence ([#1339](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1339))
- *(normalize)* answer the del→repeat codon gate on the tract, not the input span ([#1340](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1340))
- *(normalize)* repair del+dup collisions on the CDS-relative axes ([#1338](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1338))
- *(normalize)* unify the complement helpers and case-fold the delins path ([#1336](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1336))
- *(liftover)* return an error for position 0 instead of panicking ([#1335](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1335))
- *(normalize)* order protein cis members by junction, like the nucleotide axes ([#1333](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1333))
- *(project)* decline axes ferro cannot render, instead of emitting them ([#1332](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1332))
- *(normalize)* describe a payload that lands before the first base ([#1331](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1331))
- *(normalize)* repair del+dup collisions on the non-coding axis ([#1315](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1315))
- *(normalize)* order a shared junction's payloads by where they came from ([#1326](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1326))
- *(normalize)* let a duplication cover the identity member inside it ([#1324](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1324))
- *(normalize)* re-spell a duplication that spans a sibling's junction ([#1322](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1322))
- *(normalize)* combine two repeats that grew one tract ([#1319](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1319))
- *(normalize)* test commuting against the payload a member lands with ([#1317](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1317))
- *(normalize)* bound a commuting sibling at the junction it settles on ([#1313](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1313))
- *(normalize)* drop a cancelled member left inside the one that absorbed it ([#1309](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1309))
- *(normalize)* read a junction barrier's payload from its own snapshot ([#1306](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1306))
- *(normalize)* order two members sharing a span by where they add sequence ([#1305](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1305))
- *(normalize)* let a repeat claim the bases under its tract ([#1299](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1299))
- *(normalize)* keep a moved member's junction payload in phase ([#1298](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1298))
- *(normalize)* stop a junction shifting past another junction ([#1293](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1293))
- *(normalize)* stop a repeat's tract swallowing a sibling's junction ([#1291](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1291))
- *(normalize)* merge cis members that settle on one junction ([#1289](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1289))
- *(normalize)* make a duplication a barrier to a sibling's shift ([#1288](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1288))
- *(normalize)* describe a one-base inversion residue as a substitution ([#1250](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1250))
- *(python)* make the native enums honour the contract their stubs declare ([#1253](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1253))
- *(python)* return changed_positions as list[int], not bytes ([#1252](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1252))
- *(normalize)* cis member ordering and a repositioned member's clamp ([#1285](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1285))
- *(normalize)* canonicalize cis alleles on the transcript axes ([#1243](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1243))
- *(normalize)* stop a cis member's normalization crossing its siblings ([#1259](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1259))
- *(equivalence)* decline overlapping SPDI triples instead of panicking
- *(normalize)* canonicalize cis alleles from the resulting sequence ([#1237](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1237))

### Other

- *(normalize)* stop citing prioritisation for a rule it does not state ([#1348](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1348))
- *(normalize)* widen the sibling-crossing sweep to dup/ins siblings and a tandem tract ([#1337](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1337))
- *(spec)* give the generated spec fixtures a committed oracle ([#1330](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1330))
- *(review)* add a CodeRabbit review checklist for recurring miss patterns ([#1278](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1278))
- *(coderabbit)* repoint stale test paths, cover untuned modules, enforce liveness ([#1275](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1275))
- *(normalize)* use the shared cis-allele apply oracle in three suites ([#1329](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1329))
- *(normalize)* strengthen the cis-allele test oracles to their stated contracts ([#1310](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1310))
- *(test)* draw confluence cores by index instead of a string union ([#1300](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1300))
- *(normalize)* drop the indel model's adjacent-gap-insertion filter ([#1302](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1302))
- *(normalize)* lock the dup-junction overlap contract end to end ([#1277](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1277))
- *(normalize)* property-test cis-allele confluence ([#1238](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1238))
- *(dev)* make the spec generators binaries so cargo locates them ([#1228](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1228))

## [0.11.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.10.1...v0.11.0) - 2026-07-27

### Added

- *(error-handling)* diagnose a bare `ins` as W3027 InsertionWithoutInsertedSequence ([#1166](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1166))
- *(prepare)* record the cdot data release per artifact ([#1164](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1164))
- *(reference)* per-artifact schema handshake + recorded cdot data version (#1001 item 3) ([#1154](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1154))
- *(reference)* stamp & verify a reference content identity (#1001 items 2+4) ([#1153](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1153))

### Fixed

- *(test)* pass the target directory to the generator-example build ([#1225](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1225))
- *(test)* build the generator examples for the running profile ([#1224](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1224))
- *(error-handling)* make every enforced warning code reach the error configuration ([#1220](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1220))
- *(project)* surface the engine's own reason for an unavailable axis ([#1221](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1221))
- *(normalize)* clamp an insertion saturated at a mitochondrial terminus ([#1222](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1222))
- *(normalize)* clamp an insertion saturated at a contig bound on the g. axis ([#1215](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1215))
- *(normalize)* fall back to dup, not ins, where the codon gate refuses a repeat ([#1213](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1213))
- *(project)* refuse out-of-transcript positions and surface projection warnings ([#1193](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1193))
- *(test)* stop spec_generator_preconditions running a stale example binary ([#1219](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1219))
- *(cli)* pass the error configuration to the normalizer ([#1191](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1191))
- *(normalize)* answer the repeat codon gate for the tract, not the input span ([#1212](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1212))
- *(normalize)* shift a c.-axis insertion past cds_end in one pass ([#1211](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1211))
- *(normalize)* clamp a saturated insertion at the CDS bounds on the r. axis ([#1208](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1208))
- *(normalize)* complete the 3'-shift across cds_end in one pass for delins and dup ([#1189](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1189))
- *(normalize)* clamp boundary-saturated insertions at the transcript bounds ([#1203](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1203))
- *(project)* read r.-axis inputs as CDS-relative, not transcript-relative ([#1179](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1179))
- *(normalize)* apply the codon-frame gate on the r. axis ([#1194](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1194))
- *(normalize)* run the ins[...] expansion on the r. axis ([#1188](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1188))
- *(normalize)* handle an inv suffix on a cross-reference insertion range ([#1187](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1187))
- *(normalize)* maximize the repeat tract instead of taking the first literal hit ([#1176](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1176))
- *(normalize)* direction-aware copy selection in normalize_repeat; fuzz unit[N] inputs ([#1172](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1172))
- *(normalize)* 5' boundary + repeat/dup canonicalization; enable the general 5' idempotency proptest ([#1170](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1170))
- *(dev)* make the spec generators usable on a fresh worktree ([#1201](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1201))
- *(normalize)* rotate the inserted sequence when 5'-shifting a multi-base insertion ([#1169](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1169))
- *(normalize)* make the #418 CDS-start boundary clamp a stable fixed point ([#1167](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1167))
- compare resulting sequence in EquivalenceChecker for complex indels ([#1158](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1158)) ([#1161](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1161))
- *(reference)* enforce the canonical-overrides schema version at load ([#1155](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1155))
- 3'-shift delins that reduces to a deletion or duplication ([#1157](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1157)) ([#1160](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1160))
- *(parser)* reject an insertion with no inserted sequence (bare `ins`) ([#1152](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1152))

### Other

- *(reference)* construct the FASTA index through a single PreparedIndex type ([#1186](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1186))
- *(reference)* stop building FASTA-index state that is immediately discarded ([#1178](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1178))
- *(normalize)* require an explicit ErrorConfig at every entry point ([#1223](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1223))
- *(io)* buffer serde_json reads from files ([#1171](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1171))
- *(compliance)* pin cross-document HGVS rules (round-trip identity + contiguous changes) ([#1190](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1190))
- *(conformance)* accept the sized-suffix parse-rejection taxonomy divergence ([#1173](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1173))
- *(normalize)* fold is_coding + cds_span into one CodonGate ([#1214](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1214))
- shard the test suite and gate every PR on a 1M-case idempotency soak ([#1216](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1216))
- *(normalize)* fuzz the n./r. axes and inv/con edits for idempotency ([#1180](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1180))
- *(normalize)* extend the idempotency oracle to the projector path ([#1174](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1174))
- *(normalize)* make the direction-invariance test and proptest flanks honest ([#1175](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1175))
- [**breaking**] mark the remaining additive-change-prone public types non_exhaustive ([#1168](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1168))
- [**breaking**] make the public equivalence API robust to additive changes ([#1165](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1165))

## [0.10.1](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.10.0...v0.10.1) - 2026-07-23

### Fixed

- *(lenient)* stop the one-letter AA expansion at the next member's accession ([#1150](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1150))
- *(render)* apply the protein render style to an allele's members ([#1149](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1149))
- *(parser)* accept a specification on custom and assembly references ([#1148](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1148))
- *(hgvs)* single-source the gene-symbol selector's Display rule ([#1141](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1141))
- *(render)* drop the gene-symbol selector on every protein allele ([#1144](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1144))
- *(project)* name the codons an all-silent change rewrote ([#1140](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1140))
- *(normalize)* describe a self-cancelling del+ins merge as an identity ([#1136](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1136))

### Other

- *(hgvs)* pin that an assembly reference keeps its gene-symbol selector ([#1147](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1147))

## [0.10.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.9.1...v0.10.0) - 2026-07-22

### Fixed

- *(normalize)* scope the remaining protein coalesce declines to the affected members ([#1134](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1134))
- *(normalize)* coalesce a protein cis run whose to-Ter member is last ([#1133](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1133))
- *(normalize)* fall back to named endpoints when a provider lacks the accession ([#1132](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1132))
- *(normalize)* scope the protein coalesce to-Ter decline to the affected run ([#1127](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1127))
- *(normalize)* minimize a named-flank pure insertion reference-free ([#1128](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1128))
- *(project)* retain a genomic input's accession as c./n. context ([#1115](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1115))
- *(normalize)* canonicalize protein delins to its minimal form (reference-free) ([#1122](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1122))
- *(normalize)* coalesce mixed sub+del protein cis runs into one delins ([#1121](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1121))
- *(parser)* accept the '*' stop-codon glyph in strict mode ([#1118](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1118))
- *(normalize)* sort protein cis members before the delins coalesce ([#1117](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1117))
- *(tooling)* make mutalyzer refresh merge fully non-destructive ([#1110](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1110))
- *(normalize)* sort cis-allele members before merge so merges are input-order-independent ([#1106](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1106))
- *(normalize)* render cis-allele members in genomic order ([#1101](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1101))
- *(service)* validate HGVS format via the parser, not a regex gate ([#1104](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1104))
- *(parser)* reject a frameshift anchored at an unchanged residue ([#1096](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1096))
- *(normalize)* canonicalize consecutive-residue protein cis substitutions to delins ([#1095](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1095))
- *(tests)* make the failure-expectations blesser merge instead of regenerate ([#1093](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1093))
- *(project)* retain the input's genomic context on the c./n. axis ([#1090](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1090))
- *(project)* decline the r. axis for intronic spans; emit genuine n. for r. inputs ([#1089](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1089))
- *(project)* decline unknown intronic offsets instead of overflowing ([#1088](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1088))
- *(project)* describe a decomposed delins as one protein consequence ([#1083](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1083))
- *(spec-fixture)* drop prose position references from the harvester and repoint the tracker ([#1082](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1082))
- *(parser)* reject six classes of spec-forbidden variant description ([#1091](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1091))

### Other

- *(audit)* scan corrections.rs and ignore #[cfg(test)] in the emission-site check ([#1124](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1124))
- *(project)* cover cross-axis r.→n. reframe in PR CI ([#1112](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1112))
- *(conformance)* type spec-enumeration row vocabularies ([#1113](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1113))
- *(spec)* gate every enumeration status against budget or allowlist ([#1108](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1108))
- *(spec)* exhaustive non-redundant HGVS spec test enumeration ([#1085](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1085))
- *(conformance)* grow the hermetic differential protein corpus ([#1084](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1084))
- *(normalize)* guard RefSeqMismatch=Reject across all substitution syntaxes ([#1100](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1100))

## [0.9.1](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.9.0...v0.9.1) - 2026-07-20

### Fixed

- *(project)* wrong protein consequence for same-codon in-cis substitutions ([#1077](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1077))

### Other

- point Zenodo DOI badge at the concept DOI ([#1075](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1075))

## [0.9.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.8.1...v0.9.0) - 2026-07-19

### Added

- opt-in protein p. render style (Ter vs * / 3- vs 1-letter) ([#1056](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1056))
- *(reference)* reject unknown manifest fields at load, fail loud on a misread reference ([#1055](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1055))

### Fixed

- render stop-spanning deletions/delins as frameshift when a sense residue changes ([#1073](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1073))
- derive projection is_frameshift from the protein consequence ([#1072](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1072))
- phase-aware frameshift aggregation with combined cis protein consequence ([#1071](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1071))
- report an allele's coordinate axis from HgvsVariant ([#1064](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1064))
- *(normalize)* preserve the uncertain/predicted (...) wrapper through normalization ([#1063](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1063)) ([#1065](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1065))
- *(hgvs)* make allele compaction gene-agnostic on a transcript reference ([#1062](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1062))
- *(hgvs)* classify species-qualified Ensembl transcripts uniformly ([#1061](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1061))
- *(normalize)* validate substitution reference bases and expose error_config to Python ([#1052](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1052)) ([#1060](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1060))
- *(hgvs)* drop the gene-symbol selector on transcript-reference Display ([#1054](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1054))
- *(convert-gff)* fail fast when the FASTA cannot cover an exon, instead of masking it ([#1053](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1053))
- *(normalize)* clamp the mitochondrial fetch window to contig length ([#1048](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1048))
- *(dev)* resolve spec submodule commit_sha under git hook env ([#1047](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1047))
- *(normalize)* clamp genome fetch window to contig length ([#1042](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1042))

### Other

- build generate_spec_fixture example into the outer target dir ([#1074](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1074))
- bump minor version for feature commits on 0.x releases ([#1066](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1066))
- *(conformance)* update issue_506 expectation for #1054 transcript-reference selector drop ([#1067](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1067))
- *(normalize)* probe inv over-recognition boundary for issue #1040 ([#1043](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1043))

## [0.8.1](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.8.0...v0.8.1) - 2026-07-16

### Added

- *(python)* expose convert_gff and build_transcript, mirroring prepare_reference_data ([#1037](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1037))

### Fixed

- *(normalize)* only type a maximal contiguous run as inv, not a sub-run ([#1036](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1036))

## [0.8.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.7.1...v0.8.0) - 2026-07-15

### Added

- *(reference,cli)* validate transcripts.json version and check it standalone ([#1031](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1031))
- *(normalize)* warn-and-degrade when genomic data is unavailable ([#1015](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1015))
- *(cli)* genome-capable transcripts.json via --emit-genomic-sequences ([#1029](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1029))
- *(python)* value equality on projections/effects, idempotency tests, typed-Raises docs ([#1022](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1022))
- *(python)* surface normalization warnings on projections ([#1021](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1021))
- *(python)* batch error-isolation on the projector (return_exceptions) ([#1020](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1020))
- *(python)* typed exception hierarchy, py.typed, and structured error codes ([#1019](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1019))
- *(python)* surface reference capability and warn on reduced-capability build ([#1014](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1014))

### Fixed

- *(python)* reject an unrecognized direction argument instead of defaulting to 3' ([#1017](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1017))

### Other

- make public error/projection API robust to additive changes ([#1033](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1033))
- *(pre-commit)* fix the broken mypy hook (bump to v2.2.0, check by directory) ([#1023](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1023))
- *(reference)* rename MockProvider to JsonProvider ([#1030](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1030))
- *(readme)* add PyPI version and Python-versions badges ([#1010](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1010))

## [0.7.1](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.7.0...v0.7.1) - 2026-07-09

### Fixed

- *(ci)* publish wheels directly to PyPI and drop the TestPyPI gate that had failed on every prior release, so Python wheels are published to PyPI for the first time ([#1008](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1008))

## [0.7.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.6.0...v0.7.0) - 2026-07-09

### Added

- *(reference)* validate manifest schema/version at load, fail loud on an incompatible reference ([#1003](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1003))
- *(mosaic)* parse predicted-wrapper and whole-entity-LHS =/ forms ([#992](https://github.com/fulcrumgenomics/ferro-hgvs/pull/992))
- *(protein)* parse insXaa[n] and insTer<n>/ins*<n> inserted sequences ([#991](https://github.com/fulcrumgenomics/ferro-hgvs/pull/991))
- *(python)* forward the remaining MultiFasta-overridden methods to the Python provider ([#977](https://github.com/fulcrumgenomics/ferro-hgvs/pull/977))
- *(python)* forward legacy gene-selector resolution to the Python provider ([#967](https://github.com/fulcrumgenomics/ferro-hgvs/pull/967))
- feat(project)+test(conformance): source-scope enumeration & provision Ensembl reference ([#933](https://github.com/fulcrumgenomics/ferro-hgvs/pull/933)) ([#942](https://github.com/fulcrumgenomics/ferro-hgvs/pull/942))
- *(normalize)* synthesize the transcript selector on a bare-NG_ c. input ([#923](https://github.com/fulcrumgenomics/ferro-hgvs/pull/923)) ([#932](https://github.com/fulcrumgenomics/ferro-hgvs/pull/932))
- *(python)* normalize project_to_genomic output by default; raw via normalize=False ([#888](https://github.com/fulcrumgenomics/ferro-hgvs/pull/888))
- *(project)* NG_/LRG in-place c./n.→g. projection for compound alleles and UTR coordinates ([#881](https://github.com/fulcrumgenomics/ferro-hgvs/pull/881))
- *(project)* predict C-terminal extension / stop-loss for CDS→3'UTR-spanning deletion ([#880](https://github.com/fulcrumgenomics/ferro-hgvs/pull/880))
- *(reference)* provision NM_002001.2 conformance rows ([#802](https://github.com/fulcrumgenomics/ferro-hgvs/pull/802)) ([#849](https://github.com/fulcrumgenomics/ferro-hgvs/pull/849))
- *(prepare)* version-aware transcript backfill ([#842](https://github.com/fulcrumgenomics/ferro-hgvs/pull/842)) ([#848](https://github.com/fulcrumgenomics/ferro-hgvs/pull/848))
- *(service/effect)* real amino-acid resolution and junction-based NMD ([#837](https://github.com/fulcrumgenomics/ferro-hgvs/pull/837))
- *(vcf)* RNA→VCF conversion and complex-allele decomposition ([#822](https://github.com/fulcrumgenomics/ferro-hgvs/pull/822))
- *(vcf)* parse per-sample FORMAT genotype values ([#823](https://github.com/fulcrumgenomics/ferro-hgvs/pull/823))
- *(project)* translate downstream variants on non-AUG-initiation transcripts ([#780](https://github.com/fulcrumgenomics/ferro-hgvs/pull/780)) ([#796](https://github.com/fulcrumgenomics/ferro-hgvs/pull/796))
- *(reference)* provision exon→genome structure for cdot-absent old transcript versions (NM_003002.2) ([#790](https://github.com/fulcrumgenomics/ferro-hgvs/pull/790)) ([#795](https://github.com/fulcrumgenomics/ferro-hgvs/pull/795))
- *(project)* resolve legacy GENE_vNNN selectors on NG_(...):c./n./r. projection ([#784](https://github.com/fulcrumgenomics/ferro-hgvs/pull/784))
- *(project)* predict p.? for an initiation-codon member of a cis compound allele ([#778](https://github.com/fulcrumgenomics/ferro-hgvs/pull/778))
- *(normalize)* expand r. (RNA) cross-reference insert/delins payloads ([#777](https://github.com/fulcrumgenomics/ferro-hgvs/pull/777))
- *(normalize)* expand axis-less cross-reference insert payloads ([#769](https://github.com/fulcrumgenomics/ferro-hgvs/pull/769))
- *(project)* predict protein consequence for whole-exon deletions ([#761](https://github.com/fulcrumgenomics/ferro-hgvs/pull/761))
- *(conformance)* corpus-driven NG_ placement derivation ([#744](https://github.com/fulcrumgenomics/ferro-hgvs/pull/744)) ([#747](https://github.com/fulcrumgenomics/ferro-hgvs/pull/747))
- *(conformance)* hermetic genomic-axis PR gate ([#737](https://github.com/fulcrumgenomics/ferro-hgvs/pull/737)) ([#746](https://github.com/fulcrumgenomics/ferro-hgvs/pull/746))
- *(prepare)* produce derived_refseqgene_placements via --derive-ng-placements ([#740](https://github.com/fulcrumgenomics/ferro-hgvs/pull/740)) ([#743](https://github.com/fulcrumgenomics/ferro-hgvs/pull/743))
- *(project)* emit the g. axis for NG/LRG/NC-parent c.pter/c.qter inputs ([#537](https://github.com/fulcrumgenomics/ferro-hgvs/pull/537)) ([#739](https://github.com/fulcrumgenomics/ferro-hgvs/pull/739))
- *(conformance)* hermetic protein-axis PR gate ([#719](https://github.com/fulcrumgenomics/ferro-hgvs/pull/719)) ([#733](https://github.com/fulcrumgenomics/ferro-hgvs/pull/733))
- *(conformance)* hermetic snapshot-backed normalized-axis PR gate (#719 I4) ([#729](https://github.com/fulcrumgenomics/ferro-hgvs/pull/729))
- *(reference)* provision derived NG_ placements — load + produce ([#728](https://github.com/fulcrumgenomics/ferro-hgvs/pull/728)) ([#732](https://github.com/fulcrumgenomics/ferro-hgvs/pull/732))
- *(reference)* derive version-independent NG_ genomic placement — core ([#728](https://github.com/fulcrumgenomics/ferro-hgvs/pull/728)) ([#731](https://github.com/fulcrumgenomics/ferro-hgvs/pull/731))
- *(conformance)* pinned version-exact transcript reference snapshot (#719 I2) ([#727](https://github.com/fulcrumgenomics/ferro-hgvs/pull/727))
- *(normalize,project)* resolve legacy gene-model selectors on genomic references ([#709](https://github.com/fulcrumgenomics/ferro-hgvs/pull/709))
- *(reference)* cross-build version-exact cdot fallback for c. normalization ([#720](https://github.com/fulcrumgenomics/ferro-hgvs/pull/720))
- *(project)* curate project_*_all enumeration — collapse superseded versions, prefer curated transcripts ([#710](https://github.com/fulcrumgenomics/ferro-hgvs/pull/710))
- *(project,python)* --assembly genome-build override for build-agnostic NG_/LRG_ inputs ([#723](https://github.com/fulcrumgenomics/ferro-hgvs/pull/723))
- *(reference,project)* bundle GRCh37 RefSeqGene alignment + thread build through projection ([#721](https://github.com/fulcrumgenomics/ferro-hgvs/pull/721))
- *(reference)* build-aware NG_/LRG_ genomic placement (GRCh37 + GRCh38) foundation ([#711](https://github.com/fulcrumgenomics/ferro-hgvs/pull/711))
- *(conformance)* add reference-gap report tool ([#722](https://github.com/fulcrumgenomics/ferro-hgvs/pull/722))
- *(cli,project)* ferro project --axis {g,c,n,p,r} multi-axis projection surface ([#712](https://github.com/fulcrumgenomics/ferro-hgvs/pull/712))
- *(conformance)* add corpus accession inventory ([#707](https://github.com/fulcrumgenomics/ferro-hgvs/pull/707))
- *(conformance)* add reference_unavailable disposition ([#706](https://github.com/fulcrumgenomics/ferro-hgvs/pull/706))
- *(project)* add c.→r. RNA consequence prediction surface ([#701](https://github.com/fulcrumgenomics/ferro-hgvs/pull/701))
- *(normalize)* apply the 3′ rule across exon/intron boundaries ([#670](https://github.com/fulcrumgenomics/ferro-hgvs/pull/670)) ([#700](https://github.com/fulcrumgenomics/ferro-hgvs/pull/700))
- *(normalize)* emit the NP_ protein-accession selector for p. on genomic references ([#502](https://github.com/fulcrumgenomics/ferro-hgvs/pull/502)) ([#676](https://github.com/fulcrumgenomics/ferro-hgvs/pull/676))
- *(prepare)* add Ensembl reference support (ENST/ENSG/ENSP) ([#677](https://github.com/fulcrumgenomics/ferro-hgvs/pull/677))
- *(check)* validate CDS start-codon consistency in `ferro check` ([#629](https://github.com/fulcrumgenomics/ferro-hgvs/pull/629)) ([#674](https://github.com/fulcrumgenomics/ferro-hgvs/pull/674))
- *(normalize)* collapse overlapping cis-allele edits into a single delins ([#667](https://github.com/fulcrumgenomics/ferro-hgvs/pull/667))
- *(project)* extend the direct c.→p. path to bare n./r. transcript inputs ([#506](https://github.com/fulcrumgenomics/ferro-hgvs/pull/506)) ([#661](https://github.com/fulcrumgenomics/ferro-hgvs/pull/661))
- *(project)* enumerate + frame c./n./r. inputs with NG_/LRG_ context ([#646](https://github.com/fulcrumgenomics/ferro-hgvs/pull/646)) ([#647](https://github.com/fulcrumgenomics/ferro-hgvs/pull/647))
- *(project)* project NG_/LRG_ genomic inputs onto transcripts ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#642](https://github.com/fulcrumgenomics/ferro-hgvs/pull/642))
- *(project)* re-anchor NG_/LRG_-parent projections into the parent frame ([#480](https://github.com/fulcrumgenomics/ferro-hgvs/pull/480)) ([#616](https://github.com/fulcrumgenomics/ferro-hgvs/pull/616))
- *(benchmark)* perf-matrix engine + measured performance tables ([#617](https://github.com/fulcrumgenomics/ferro-hgvs/pull/617))
- *(parser)* parse and/or between bracketed alleles ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#563](https://github.com/fulcrumgenomics/ferro-hgvs/pull/563))
- *(parser)* parse predicted compound cis allele [(a;b)] ([#545](https://github.com/fulcrumgenomics/ferro-hgvs/pull/545)) ([#569](https://github.com/fulcrumgenomics/ferro-hgvs/pull/569))
- *(parser)* parse trans alleles with a mosaic (;) tail ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#559](https://github.com/fulcrumgenomics/ferro-hgvs/pull/559))
- *(perf-table)* render README performance tables from a results JSON ([#604](https://github.com/fulcrumgenomics/ferro-hgvs/pull/604))
- *(reference)* load + index the GRCh37 cdot so multi-build manifests project GRCh37 inputs ([#605](https://github.com/fulcrumgenomics/ferro-hgvs/pull/605))
- *(error)* emit AlignmentGap (E3007) for variants in a transcript-genome CIGAR indel (Cycle 1c-iii) ([#603](https://github.com/fulcrumgenomics/ferro-hgvs/pull/603))
- *(parser)* parse N-padded deletion over an uncertain range ([#545](https://github.com/fulcrumgenomics/ferro-hgvs/pull/545)) ([#568](https://github.com/fulcrumgenomics/ferro-hgvs/pull/568))
- *(parser)* parse comma products allele [a,b] ([#545](https://github.com/fulcrumgenomics/ferro-hgvs/pull/545)) ([#571](https://github.com/fulcrumgenomics/ferro-hgvs/pull/571))
- *(parser)* parse positionless breakpoint insertion in composite alleles ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#584](https://github.com/fulcrumgenomics/ferro-hgvs/pull/584))
- *(tool-support)* single source of truth for the tool-comparison tables ([#590](https://github.com/fulcrumgenomics/ferro-hgvs/pull/590))
- *(parser)* parse the sup supernumerary marker ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#596](https://github.com/fulcrumgenomics/ferro-hgvs/pull/596))
- *(parser)* accept pter/qter endpoints in inserted position-ranges ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#597](https://github.com/fulcrumgenomics/ferro-hgvs/pull/597))
- *(parser)* parse same-chromosome ring `::` deletion-join ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#595](https://github.com/fulcrumgenomics/ferro-hgvs/pull/595))
- *(parser)* support inverted insertion of an uncertain-boundary range ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#576](https://github.com/fulcrumgenomics/ferro-hgvs/pull/576))
- *(parser)* accept all-identity trans alleles via cross-spell check ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#558](https://github.com/fulcrumgenomics/ferro-hgvs/pull/558))
- *(project)* add the c↔n (non-coding transcript) axis to VariantProjection ([#592](https://github.com/fulcrumgenomics/ferro-hgvs/pull/592))
- *(parser)* parse shared-position repeat trans alleles ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#555](https://github.com/fulcrumgenomics/ferro-hgvs/pull/555))
- *(parser)* parse predicted-wrapper members in protein alleles ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#552](https://github.com/fulcrumgenomics/ferro-hgvs/pull/552))
- *(parser)* parse residue-level and/or `^` substitutions ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#550](https://github.com/fulcrumgenomics/ferro-hgvs/pull/550))
- *(parser)* parse frameshift with alternative new residues ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#557](https://github.com/fulcrumgenomics/ferro-hgvs/pull/557))
- *(parser)* parse trans alleles whose members are cis groups ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#549](https://github.com/fulcrumgenomics/ferro-hgvs/pull/549))
- *(parser)* extend mosaic `=/` to protein + fix RNA case ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#548](https://github.com/fulcrumgenomics/ferro-hgvs/pull/548))
- *(parser)* parse the and/or `^` operator ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#547](https://github.com/fulcrumgenomics/ferro-hgvs/pull/547))
- *(conformance)* generate failure-patterns.md from cases.json; retire hand-maintained counts ([#509](https://github.com/fulcrumgenomics/ferro-hgvs/pull/509)) ([#538](https://github.com/fulcrumgenomics/ferro-hgvs/pull/538))
- *(reference)* reconcile cdot for synthesis-only transcripts via exon-sum length ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#532](https://github.com/fulcrumgenomics/ferro-hgvs/pull/532))
- *(reference)* ingest protein FASTAs + get_protein_sequence + translation check ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#530](https://github.com/fulcrumgenomics/ferro-hgvs/pull/530))
- *(reference)* apply canonical overrides in MultiFastaProvider + reconcile cdot ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#528](https://github.com/fulcrumgenomics/ferro-hgvs/pull/528))
- *(prepare)* canonical-overrides acquisition — ferro prepare --validate-canonical ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#527](https://github.com/fulcrumgenomics/ferro-hgvs/pull/527))
- *(normalize)* resolve bare-transcript c. pter/qter telomere markers (#488 Phase 2) ([#533](https://github.com/fulcrumgenomics/ferro-hgvs/pull/533))
- *(reference)* version-exact correction of transcript metadata ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#524](https://github.com/fulcrumgenomics/ferro-hgvs/pull/524))
- *(reference)* offline structural validation of transcript records ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#521](https://github.com/fulcrumgenomics/ferro-hgvs/pull/521))
- *(normalize)* resolve genomic pter/qter telomere markers to concrete coordinates ([#488](https://github.com/fulcrumgenomics/ferro-hgvs/pull/488)) ([#526](https://github.com/fulcrumgenomics/ferro-hgvs/pull/526))
- *(project)* direct c.→p. path for bare-NM_ inputs; demote 4 bare-NP protein rows ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#508](https://github.com/fulcrumgenomics/ferro-hgvs/pull/508))
- *(reference)* strict exact-version transcript resolution (#478 pillar 3) ([#490](https://github.com/fulcrumgenomics/ferro-hgvs/pull/490))
- *(error-handling)* reject standalone single-member allele brackets uniformly (W3026, closes #493) ([#496](https://github.com/fulcrumgenomics/ferro-hgvs/pull/496))
- *(error-handling)* add W3023/W3024/W3025 for dup<N>/dup<seq>/del<seq> soft-prohibition forms (closes #460) ([#462](https://github.com/fulcrumgenomics/ferro-hgvs/pull/462))
- *(parser)* accept whole-entity edits in protein bracket members (closes #468) ([#475](https://github.com/fulcrumgenomics/ferro-hgvs/pull/475))
- *(normalize)* expand ins[ACC:g.A_B] cross-reference to flat literal (closes #422) ([#437](https://github.com/fulcrumgenomics/ferro-hgvs/pull/437))
- *(normalize)* canonicalize p.delins → p.dup and surface W3022 InitiatorMetCanonicalization (closes #92) ([#407](https://github.com/fulcrumgenomics/ferro-hgvs/pull/407))
- contig-length-aware math + wraparound boundary policies for m./o. variants (closes #399) ([#411](https://github.com/fulcrumgenomics/ferro-hgvs/pull/411))
- *(normalize,project,python,validate)* four misc gaps from #395 (items 1, 4, 5, 6) ([#421](https://github.com/fulcrumgenomics/ferro-hgvs/pull/421))
- *(parser)* spec-compliant wraparound ranges on m./o. circular refs ([#380](https://github.com/fulcrumgenomics/ferro-hgvs/pull/380))
- *(normalize)* protein 3'-shifting for deletions and duplications ([#377](https://github.com/fulcrumgenomics/ferro-hgvs/pull/377))
- *(normalize)* expand ins[...] to flat literal sequence (closes #333) ([#365](https://github.com/fulcrumgenomics/ferro-hgvs/pull/365))
- *(normalize)* structured info-code surface mirroring W#### warnings ([#373](https://github.com/fulcrumgenomics/ferro-hgvs/pull/373))
- *(project)* VariantProjector accepts c./n./r. inputs ([#379](https://github.com/fulcrumgenomics/ferro-hgvs/pull/379))
- *(error_handling)* mutalyzer ↔ ferro error-code mapping table (closes #329) ([#362](https://github.com/fulcrumgenomics/ferro-hgvs/pull/362))
- *(error-handling)* W5003 VariantExceedsReference, strict-mode rejection (closes #355) ([#363](https://github.com/fulcrumgenomics/ferro-hgvs/pull/363))
- *(test/mutalyzer-normalize)* accepted_divergence + spec_citation (closes #335) ([#359](https://github.com/fulcrumgenomics/ferro-hgvs/pull/359))
- *(parser)* canonicalize r. thymine input to u with W3020 soft-warning (closes #282) ([#299](https://github.com/fulcrumgenomics/ferro-hgvs/pull/299))
- *(parser)* recognize non-human Ensembl accession prefixes ([#320](https://github.com/fulcrumgenomics/ferro-hgvs/pull/320))
- *(project)* project_to_genomic for c./n./r. inputs (closes #327) ([#358](https://github.com/fulcrumgenomics/ferro-hgvs/pull/358))
- *(parser)* RNA/TX trans-allele predicted-edit wrapper in compound brackets (closes #287) ([#309](https://github.com/fulcrumgenomics/ferro-hgvs/pull/309))
- *(parser)* whole-entity predicted forms for m./o./n. variants (closes #288) ([#308](https://github.com/fulcrumgenomics/ferro-hgvs/pull/308))
- *(normalize)* MultiRepeat partial validation for non-Exact counts (closes #279) ([#296](https://github.com/fulcrumgenomics/ferro-hgvs/pull/296))
- *(annotation)* populate Transcript.protein_id from GFF/GTF attributes ([#321](https://github.com/fulcrumgenomics/ferro-hgvs/pull/321))
- *(parser)* whole-entity predicted forms `c.(=)`, `c.(?)`, `r.(=)`, `r.(?)`, `r.(0)` (closes #245) ([#246](https://github.com/fulcrumgenomics/ferro-hgvs/pull/246))
- *(parser)* predicted-edit wrapper inside compound brackets (closes #243) ([#244](https://github.com/fulcrumgenomics/ferro-hgvs/pull/244))
- *(parser)* accept ?con<src> and ?copy<N> on unknown position (closes #286) ([#303](https://github.com/fulcrumgenomics/ferro-hgvs/pull/303))
- *(parser)* bracket-aware split for chimeric + non-compact mosaic forms (closes #216) ([#217](https://github.com/fulcrumgenomics/ferro-hgvs/pull/217))
- *(normalize)* validate repeat unit divides reference span and matches bases (closes #214) ([#215](https://github.com/fulcrumgenomics/ferro-hgvs/pull/215))

### Fixed

- *(ci)* satisfy Rust 1.97 clippy lints and pin the toolchain ([#1007](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1007))
- *(normalize)* preserve co-located insertion overlaps instead of a non-idempotent merge ([#1005](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1005))
- *(normalize)* coalesce shift-created cis adjacency into a single delins ([#1002](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1002))
- *(project)* decline c./p./r. on 5′-incomplete-CDS (cds_start_NF) transcripts ([#972](https://github.com/fulcrumgenomics/ferro-hgvs/pull/972)) ([#994](https://github.com/fulcrumgenomics/ferro-hgvs/pull/994))
- *(reference)* resolve transcript_genome_span for deferred Ensembl transcripts ([#995](https://github.com/fulcrumgenomics/ferro-hgvs/pull/995)) ([#996](https://github.com/fulcrumgenomics/ferro-hgvs/pull/996))
- *(normalize)* render plain past-CDS coordinates in canonical c.*N form (#920/#336) ([#987](https://github.com/fulcrumgenomics/ferro-hgvs/pull/987))
- *(spdi)* bound exonic c./n./r. positions past the 3' end; resolve coding r.N as CDS-relative ([#989](https://github.com/fulcrumgenomics/ferro-hgvs/pull/989))
- *(convert)* CDS<->tx is transcript-native across intra-exon CIGAR insertions ([#944](https://github.com/fulcrumgenomics/ferro-hgvs/pull/944)) ([#988](https://github.com/fulcrumgenomics/ferro-hgvs/pull/988))
- *(normalize)* canonicalize bracketed inserted-repeat counts >= 2 ([#920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/920)) ([#983](https://github.com/fulcrumgenomics/ferro-hgvs/pull/983))
- *(parser)* validate compound-reference accession pairings (outer=gene/genomic) ([#981](https://github.com/fulcrumgenomics/ferro-hgvs/pull/981))
- *(effect)* stop_retained, codon-based start_lost, and deterministic most_severe ([#978](https://github.com/fulcrumgenomics/ferro-hgvs/pull/978))
- *(spdi)* reject r.*N/c.*N positions past the transcript 3' end ([#970](https://github.com/fulcrumgenomics/ferro-hgvs/pull/970))
- *(spdi)* resolve r.-N (5'UTR) through the exon-aware mapper so it agrees with c.-N ([#969](https://github.com/fulcrumgenomics/ferro-hgvs/pull/969))
- *(vcf)* decline unknown-position/edit c./n./g. variants instead of panicking ([#943](https://github.com/fulcrumgenomics/ferro-hgvs/pull/943)) ([#947](https://github.com/fulcrumgenomics/ferro-hgvs/pull/947))
- *(project)* make c.→p. consequence decline observable, not silent ([#956](https://github.com/fulcrumgenomics/ferro-hgvs/pull/956))
- *(parser)* resolve ENSG(ENST) compound refs; require versions on Ensembl accessions ([#933](https://github.com/fulcrumgenomics/ferro-hgvs/pull/933)) ([#938](https://github.com/fulcrumgenomics/ferro-hgvs/pull/938))
- *(normalize)* reconcile MultiRepeat needs_normalization with the shuffle dispatch ([#958](https://github.com/fulcrumgenomics/ferro-hgvs/pull/958))
- *(spdi)* resolve r.*N through the exon-aware mapper so it agrees with c.*N ([#944](https://github.com/fulcrumgenomics/ferro-hgvs/pull/944)) ([#950](https://github.com/fulcrumgenomics/ferro-hgvs/pull/950))
- *(parser)* correct self-cancelling-allele guard docs and fixture classification ([#959](https://github.com/fulcrumgenomics/ferro-hgvs/pull/959))
- *(prepare)* make --dry-run a true no-op preview ([#939](https://github.com/fulcrumgenomics/ferro-hgvs/pull/939)) ([#940](https://github.com/fulcrumgenomics/ferro-hgvs/pull/940))
- *(normalize)* collapse overlapping cis del+ins allele members on the transcript axis ([#920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/920)) ([#929](https://github.com/fulcrumgenomics/ferro-hgvs/pull/929))
- *(normalize)* 3'-shift CDS and whole-CDS-spanning del/dup into the UTR per the HGVS 3'-rule ([#918](https://github.com/fulcrumgenomics/ferro-hgvs/pull/918)) ([#935](https://github.com/fulcrumgenomics/ferro-hgvs/pull/935))
- *(normalize)* canonicalize degenerate bracketed inserted-repeat counts ([#920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/920)) ([#926](https://github.com/fulcrumgenomics/ferro-hgvs/pull/926))
- *(project)* render a codon-straddling single-codon in-frame deletion as del, not empty delins ([#931](https://github.com/fulcrumgenomics/ferro-hgvs/pull/931))
- *(benchmark)* surface the failure reason when ferro normalize yields zero successes ([#916](https://github.com/fulcrumgenomics/ferro-hgvs/pull/916))
- *(project)* render premature-stop protein consequences per spec, not a C-terminal delins ([#911](https://github.com/fulcrumgenomics/ferro-hgvs/pull/911)) ([#913](https://github.com/fulcrumgenomics/ferro-hgvs/pull/913))
- *(reference)* skip hidden/AppleDouble files when scanning FASTA directories ([#915](https://github.com/fulcrumgenomics/ferro-hgvs/pull/915))
- *(project)* cis-sort minus-strand compound-allele members on the project_variant genomic axis ([#894](https://github.com/fulcrumgenomics/ferro-hgvs/pull/894)) ([#898](https://github.com/fulcrumgenomics/ferro-hgvs/pull/898))
- *(project)* de-anchor bare NG_/LRG_ genomic input on the single-transcript projection paths ([#879](https://github.com/fulcrumgenomics/ferro-hgvs/pull/879)) ([#896](https://github.com/fulcrumgenomics/ferro-hgvs/pull/896))
- *(project)* surface the genomic axis for LRG inputs and pter/qter termini on project_variant (#886, #887) ([#893](https://github.com/fulcrumgenomics/ferro-hgvs/pull/893))
- *(normalize)* keep out-of-phase insertions as ins, not spurious dup/repeat ([#882](https://github.com/fulcrumgenomics/ferro-hgvs/pull/882)) ([#892](https://github.com/fulcrumgenomics/ferro-hgvs/pull/892))
- *(normalize)* 3'-align insertion→duplication; guarantee normalization idempotency ([#883](https://github.com/fulcrumgenomics/ferro-hgvs/pull/883))
- *(project)* preserve LRG_t/LRG_p namespace on projection output ([#860](https://github.com/fulcrumgenomics/ferro-hgvs/pull/860)) ([#874](https://github.com/fulcrumgenomics/ferro-hgvs/pull/874))
- render in-frame single-codon deletion in a residue run as del, not empty ins ([#850](https://github.com/fulcrumgenomics/ferro-hgvs/pull/850))
- *(project)* spec-canonical minus-strand genomic repeat projection ([#852](https://github.com/fulcrumgenomics/ferro-hgvs/pull/852)) ([#869](https://github.com/fulcrumgenomics/ferro-hgvs/pull/869))
- *(parser)* collapse single-element bracketed literal ins/delins to a plain literal ([#863](https://github.com/fulcrumgenomics/ferro-hgvs/pull/863))
- *(project)* canonicalize protein consequences (delins→sub, tandem dup, adjacent in-cis delins) ([#862](https://github.com/fulcrumgenomics/ferro-hgvs/pull/862))
- *(project)* build-scope the predicted-r. allele transcript fetch ([#843](https://github.com/fulcrumgenomics/ferro-hgvs/pull/843)) ([#845](https://github.com/fulcrumgenomics/ferro-hgvs/pull/845))
- *(normalize)* treat LRG transcript references as bare for EINTRONIC ([#834](https://github.com/fulcrumgenomics/ferro-hgvs/pull/834)) ([#844](https://github.com/fulcrumgenomics/ferro-hgvs/pull/844))
- *(project)* genomic projection of c.* variants in the 3'UTR poly-A region (NM_003002.2; #790 residual) ([#839](https://github.com/fulcrumgenomics/ferro-hgvs/pull/839))
- *(project)* wire rna_description axis and frame predicted r. from input ([#838](https://github.com/fulcrumgenomics/ferro-hgvs/pull/838))
- *(vcf)* drop per-sample genotype when splitting multi-allelic records ([#841](https://github.com/fulcrumgenomics/ferro-hgvs/pull/841))
- *(reference)* apply cdot CIGAR during transcript synthesis, decline on insertions ([#831](https://github.com/fulcrumgenomics/ferro-hgvs/pull/831))
- *(project)* force initiator Met at residue 1 for non-AUG ref-protein translation ([#836](https://github.com/fulcrumgenomics/ferro-hgvs/pull/836))
- *(reference)* suppress derived-tx injection for the #790 producer ([#835](https://github.com/fulcrumgenomics/ferro-hgvs/pull/835))
- *(reference)* cover canonical-sequence reingestion through the provider ([#791](https://github.com/fulcrumgenomics/ferro-hgvs/pull/791)) ([#832](https://github.com/fulcrumgenomics/ferro-hgvs/pull/832))
- *(vcf)* correct genomic anchor-base recovery and consolidate NC_→chrom mapping ([#821](https://github.com/fulcrumgenomics/ferro-hgvs/pull/821))
- *(vcf)* emit RefSeq NC_ accessions on the genomic HGVS axis ([#820](https://github.com/fulcrumgenomics/ferro-hgvs/pull/820))
- *(equivalence/service)* make unsupported variant classes explicit ([#816](https://github.com/fulcrumgenomics/ferro-hgvs/pull/816))
- *(data/cdot)* drop frequently-wrong NM_→NP_ protein-accession inference ([#815](https://github.com/fulcrumgenomics/ferro-hgvs/pull/815))
- *(vcf)* classify coding indels by typed NaEdit and emit Frameshift ([#817](https://github.com/fulcrumgenomics/ferro-hgvs/pull/817))
- *(error)* give TranscriptVersionNotExact a distinct error code (E2004) ([#814](https://github.com/fulcrumgenomics/ferro-hgvs/pull/814))
- *(normalize)* resolve legacy GENE_v001 selector against the NG_ parent's hosted transcript ([#792](https://github.com/fulcrumgenomics/ferro-hgvs/pull/792)) ([#793](https://github.com/fulcrumgenomics/ferro-hgvs/pull/793))
- decline silent transcript version substitution in normalization and projection ([#787](https://github.com/fulcrumgenomics/ferro-hgvs/pull/787))
- *(project)* report p.? for an initiation-codon variant on a non-AUG transcript ([#772](https://github.com/fulcrumgenomics/ferro-hgvs/pull/772))
- *(data)* anchor minus-strand intronic offsets to the correct exon end ([#766](https://github.com/fulcrumgenomics/ferro-hgvs/pull/766))
- *(convert)* map c.*N+offset 3'UTR positions past a transcript terminus ([#765](https://github.com/fulcrumgenomics/ferro-hgvs/pull/765))
- *(normalize)* resolve the transcript for a genomic-context c. cross-reference payload ([#654](https://github.com/fulcrumgenomics/ferro-hgvs/pull/654)) ([#756](https://github.com/fulcrumgenomics/ferro-hgvs/pull/756))
- *(normalize)* reject overlapping insertions in cis alleles (#486 EOVERLAP) ([#749](https://github.com/fulcrumgenomics/ferro-hgvs/pull/749))
- *(data)* load cdot exons in the HGVS coordinate convention ([#748](https://github.com/fulcrumgenomics/ferro-hgvs/pull/748))
- *(project)* re-normalize projected genomic variants in their own frame ([#741](https://github.com/fulcrumgenomics/ferro-hgvs/pull/741))
- *(normalize)* normalize r. edits carrying the RNA base `u` ([#736](https://github.com/fulcrumgenomics/ferro-hgvs/pull/736)) ([#738](https://github.com/fulcrumgenomics/ferro-hgvs/pull/738))
- *(normalize)* extend the exon/intron 3′ rule to n. and r. variants ([#704](https://github.com/fulcrumgenomics/ferro-hgvs/pull/704)) ([#734](https://github.com/fulcrumgenomics/ferro-hgvs/pull/734))
- *(parser)* accept a parenthesized gene-selector in a cross-reference insert ([#730](https://github.com/fulcrumgenomics/ferro-hgvs/pull/730))
- *(reference)* forward dropped methods through Arc/Box provider wrappers ([#726](https://github.com/fulcrumgenomics/ferro-hgvs/pull/726))
- *(project)* resolve range-reference insertions against the NG_/LRG_ parent before de-anchoring ([#708](https://github.com/fulcrumgenomics/ferro-hgvs/pull/708))
- *(reference)* resolve cdot version-exact in normalize, never substitute a sibling version ([#717](https://github.com/fulcrumgenomics/ferro-hgvs/pull/717))
- *(project)* re-anchor the genomic axis on the single-transcript projection path ([#702](https://github.com/fulcrumgenomics/ferro-hgvs/pull/702)) ([#703](https://github.com/fulcrumgenomics/ferro-hgvs/pull/703))
- *(normalize)* correct 1-based↔0-based off-by-one in genomic-shuffle fetch ([#690](https://github.com/fulcrumgenomics/ferro-hgvs/pull/690))
- *(normalize)* keep W4007 on warn-only intronic-bare resolve-failure ([#686](https://github.com/fulcrumgenomics/ferro-hgvs/pull/686))
- *(project)* error instead of emitting invalid HGVS when an NG_/LRG_ re-anchor declines ([#655](https://github.com/fulcrumgenomics/ferro-hgvs/pull/655)) ([#662](https://github.com/fulcrumgenomics/ferro-hgvs/pull/662))
- *(project)* sequence-aware g.↔c. projection across unmodeled exon indels ([#644](https://github.com/fulcrumgenomics/ferro-hgvs/pull/644)) ([#668](https://github.com/fulcrumgenomics/ferro-hgvs/pull/668))
- *(reference)* don't abort LRG mapping scan on a malformed earlier `<mapping>` ([#684](https://github.com/fulcrumgenomics/ferro-hgvs/pull/684))
- *(check)* warn when `refseqgene_alignments` file is missing ([#683](https://github.com/fulcrumgenomics/ferro-hgvs/pull/683))
- *(python)* add missing ErrorType discriminants 41/42 to the type stub ([#658](https://github.com/fulcrumgenomics/ferro-hgvs/pull/658))
- *(normalize)* spec-compliant repeat normalization — 3' rotation + flank absorption ([#649](https://github.com/fulcrumgenomics/ferro-hgvs/pull/649))
- *(reference)* resolve bare LRG_N genomic accession to its LRG_Ng record ([#638](https://github.com/fulcrumgenomics/ferro-hgvs/pull/638))
- *(benchmark)* make `setup uta` reliable ([#623](https://github.com/fulcrumgenomics/ferro-hgvs/pull/623)) ([#632](https://github.com/fulcrumgenomics/ferro-hgvs/pull/632))
- *(project)* decline protein prediction when the CDS frame is inconsistent ([#625](https://github.com/fulcrumgenomics/ferro-hgvs/pull/625)) ([#628](https://github.com/fulcrumgenomics/ferro-hgvs/pull/628))
- *(validate)* reject intronic offsets on a bare transcript reference (#486 EINTRONIC) ([#577](https://github.com/fulcrumgenomics/ferro-hgvs/pull/577))
- *(project)* route stop-disrupting indels to a C-terminal extension ([#615](https://github.com/fulcrumgenomics/ferro-hgvs/pull/615)) ([#621](https://github.com/fulcrumgenomics/ferro-hgvs/pull/621))
- *(benchmark)* time normalize tools by subprocess-internal timer ([#609](https://github.com/fulcrumgenomics/ferro-hgvs/pull/609)) ([#612](https://github.com/fulcrumgenomics/ferro-hgvs/pull/612))
- *(validate)* reject del/dup/delins numeric length not matching the position span (#486 length mismatch) ([#565](https://github.com/fulcrumgenomics/ferro-hgvs/pull/565))
- *(tool-support)* assign footnote markers in canonical order ([#602](https://github.com/fulcrumgenomics/ferro-hgvs/pull/602))
- *(reference)* report the true secondary-build transcript count ([#613](https://github.com/fulcrumgenomics/ferro-hgvs/pull/613))
- *(normalize)* size the intronic shuffle window to the enclosing intron ([#573](https://github.com/fulcrumgenomics/ferro-hgvs/pull/573)) ([#575](https://github.com/fulcrumgenomics/ferro-hgvs/pull/575))
- make HashMap-ordered output deterministic (VCF INFO + benchmark reports) ([#594](https://github.com/fulcrumgenomics/ferro-hgvs/pull/594))
- *(parser)* bring fast path to full parity with the generic parser (ncRNA, U base, positions, panic) ([#560](https://github.com/fulcrumgenomics/ferro-hgvs/pull/560))
- *(reference)* make cdot base→version fallback deterministic ([#583](https://github.com/fulcrumgenomics/ferro-hgvs/pull/583))
- *(conformance)* map ECOORDINATESYSTEMMISMATCH to ferro's Parse error ([#486](https://github.com/fulcrumgenomics/ferro-hgvs/pull/486)) ([#579](https://github.com/fulcrumgenomics/ferro-hgvs/pull/579))
- *(tool-support)* address CodeRabbit findings on the matrix generator ([#601](https://github.com/fulcrumgenomics/ferro-hgvs/pull/601))
- *(bench)* make the intronic projection fixture self-consistent ([#600](https://github.com/fulcrumgenomics/ferro-hgvs/pull/600))
- *(project)* emit nonsense (not fsTer1) for an immediate-stop frameshift ([#589](https://github.com/fulcrumgenomics/ferro-hgvs/pull/589))
- *(project)* scan the 3'UTR for a frameshift's new stop (fsTer{K}, not fsTer?) ([#587](https://github.com/fulcrumgenomics/ferro-hgvs/pull/587))
- *(build)* repair test code that broke under non-`dev` feature flags ([#580](https://github.com/fulcrumgenomics/ferro-hgvs/pull/580))
- *(validate)* reject delins/del/dup/repeat reference-sequence mismatches in strict mode ([#486](https://github.com/fulcrumgenomics/ferro-hgvs/pull/486)) ([#556](https://github.com/fulcrumgenomics/ferro-hgvs/pull/556))
- *(validate)* reject the RNA base U in a DNA-context edit (#486 ENODNA) ([#542](https://github.com/fulcrumgenomics/ferro-hgvs/pull/542))
- *(normalize)* coerce g. on a mitochondrial accession to m. (#487 mito group) ([#541](https://github.com/fulcrumgenomics/ferro-hgvs/pull/541))
- *(validate)* reject c./g./m. on a non-coding RNA reference (#486 coord-system mismatch) ([#543](https://github.com/fulcrumgenomics/ferro-hgvs/pull/543))
- *(normalize)* number coding-transcript r. positions CDS-relative ([#469](https://github.com/fulcrumgenomics/ferro-hgvs/pull/469)) ([#539](https://github.com/fulcrumgenomics/ferro-hgvs/pull/539))
- *(normalize)* spec-correct refusal of NG/LRG-parent c. pter/qter flank coordinates (#488 Phase 2b) ([#534](https://github.com/fulcrumgenomics/ferro-hgvs/pull/534))
- *(project)* mid-codon in-frame insertions render as delins, not a clean ins ([#511](https://github.com/fulcrumgenomics/ferro-hgvs/pull/511)) ([#517](https://github.com/fulcrumgenomics/ferro-hgvs/pull/517))
- *(normalize)* guard special/offset genome positions to avoid pterdel overflow panic ([#488](https://github.com/fulcrumgenomics/ferro-hgvs/pull/488)) ([#518](https://github.com/fulcrumgenomics/ferro-hgvs/pull/518))
- *(reference)* decline protein prediction for non-version-exact transcripts ([#505](https://github.com/fulcrumgenomics/ferro-hgvs/pull/505)) ([#519](https://github.com/fulcrumgenomics/ferro-hgvs/pull/519))
- *(project)* extend stop-loss substitution prediction with extTer notation ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#514](https://github.com/fulcrumgenomics/ferro-hgvs/pull/514))
- *(project)* predict p.(Met1?) for CDS-boundary edits reaching the initiation codon ([#504](https://github.com/fulcrumgenomics/ferro-hgvs/pull/504)) ([#513](https://github.com/fulcrumgenomics/ferro-hgvs/pull/513))
- *(hgvs)* distinguish frameshift fsTer? from short fs; predict fsTer? for no-stop frameshifts ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#515](https://github.com/fulcrumgenomics/ferro-hgvs/pull/515))
- *(project)* report p.(Met1?) for initiation-codon variants ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#512](https://github.com/fulcrumgenomics/ferro-hgvs/pull/512))
- *(project)* never emit descending positions in in-frame protein delins ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#510](https://github.com/fulcrumgenomics/ferro-hgvs/pull/510))
- *(normalize)* avoid subtract-overflow panic in minus-strand intronic shuffle ([#497](https://github.com/fulcrumgenomics/ferro-hgvs/pull/497))
- *(test/mutalyzer-normalize)* spec-correct bare-NP protein refs (7 cases) ([#495](https://github.com/fulcrumgenomics/ferro-hgvs/pull/495))
- *(reference)* flag cdot base-synthesis divergence risk (closes #471) ([#482](https://github.com/fulcrumgenomics/ferro-hgvs/pull/482))
- *(error-handling)* make detect_swapped_positions axis-aware for m./o. wraparound (closes #467) ([#474](https://github.com/fulcrumgenomics/ferro-hgvs/pull/474))
- *(cdot)* plumb transcript protein accession from cdot JSON ([#483](https://github.com/fulcrumgenomics/ferro-hgvs/pull/483))
- *(test/mutalyzer-normalize)* accept infos codes ferro does not model (22→11) ([#481](https://github.com/fulcrumgenomics/ferro-hgvs/pull/481))
- *(parser)* depth-aware separator scan in expanded allele entrypoints ([#476](https://github.com/fulcrumgenomics/ferro-hgvs/pull/476))
- *(error-handling)* extend W3009/W3010 corrector to uncertain-stop forms ([#461](https://github.com/fulcrumgenomics/ferro-hgvs/pull/461))
- *(parser)* reject single-position insertion per insertion.md:95-101 ([#446](https://github.com/fulcrumgenomics/ferro-hgvs/pull/446)) ([#450](https://github.com/fulcrumgenomics/ferro-hgvs/pull/450))
- *(edit)* drop brackets on single-payload InsertedSequence::Complex ([#443](https://github.com/fulcrumgenomics/ferro-hgvs/pull/443))
- *(service)* extend frameshift classifier span_len to Deletion + Duplication (closes #427) ([#435](https://github.com/fulcrumgenomics/ferro-hgvs/pull/435))
- *(parser)* reject `dupins<seq>` per duplication.md:92 ([#445](https://github.com/fulcrumgenomics/ferro-hgvs/pull/445)) ([#451](https://github.com/fulcrumgenomics/ferro-hgvs/pull/451))
- *(error-handling)* provider-aware W3016 for mixed-shape intronic-endpoint ranges (closes #429) ([#436](https://github.com/fulcrumgenomics/ferro-hgvs/pull/436))
- *(error-handling)* W3016 detects numeric-length-suffix disagreement on del/dup/inv (closes #439) ([#441](https://github.com/fulcrumgenomics/ferro-hgvs/pull/441))
- *(service)* inversion is in-frame at DNA level (closes #438) ([#440](https://github.com/fulcrumgenomics/ferro-hgvs/pull/440))
- *(validate)* tighten validate_multirepeat_tract for nonzero-min middle (closes #428) ([#434](https://github.com/fulcrumgenomics/ferro-hgvs/pull/434))
- *(service,validation,spdi)* three delins consistency bugs (closes #394) ([#419](https://github.com/fulcrumgenomics/ferro-hgvs/pull/419))
- *(normalize)* canon boundary-spanning ins + side-aware ins→dup anchor (closes #402) ([#413](https://github.com/fulcrumgenomics/ferro-hgvs/pull/413))
- *(normalize)* close untracked biocommons-burndown divergences (closes #418) ([#420](https://github.com/fulcrumgenomics/ferro-hgvs/pull/420))
- *(normalize)* direction-aware tie-break in insertion_to_duplication (closes #403) ([#408](https://github.com/fulcrumgenomics/ferro-hgvs/pull/408))
- fix(spdi)+(error-handling)+(test): six HGVS↔SPDI / W3016 / VCF audit gaps (closes #390) ([#406](https://github.com/fulcrumgenomics/ferro-hgvs/pull/406))
- *(project)* build-aware cdot lookups + widen fan-out path for c./n./r. (closes #389) ([#398](https://github.com/fulcrumgenomics/ferro-hgvs/pull/398))
- *(normalize)* extend W4004 PositionPastEnd to intronic offsets (closes #392) ([#414](https://github.com/fulcrumgenomics/ferro-hgvs/pull/414))
- *(normalize)* preserve spanning duplications in CDS-start canon clamp (closes #401) ([#410](https://github.com/fulcrumgenomics/ferro-hgvs/pull/410))
- *(normalize)* clamp 3'-canonicalization at CDS-end for c.-axis inputs (closes #387) ([#388](https://github.com/fulcrumgenomics/ferro-hgvs/pull/388))
- *(normalize)* reject c. positions past CDS-end / transcript-end (W4004 PositionPastEnd, closes #336) ([#342](https://github.com/fulcrumgenomics/ferro-hgvs/pull/342))
- *(normalize)* clamp 5'-canonicalization at CDS-start for c.-axis inputs (closes #383) ([#385](https://github.com/fulcrumgenomics/ferro-hgvs/pull/385))
- *(normalize)* apply exon-junction exception to n. and r. axes ([#374](https://github.com/fulcrumgenomics/ferro-hgvs/pull/374))
- *(parser)* propagate structured chunk-level diagnostics in slash-form alleles (closes #375) ([#376](https://github.com/fulcrumgenomics/ferro-hgvs/pull/376))
- *(parser)* extend E3006 self-cancelling detection to *-region and intronic ranges (closes #172) ([#371](https://github.com/fulcrumgenomics/ferro-hgvs/pull/371))
- *(diagnostics)* preserve source span on E3006 SelfCancellingAllele (closes #171) ([#369](https://github.com/fulcrumgenomics/ferro-hgvs/pull/369))
- *(preprocessor)* accept lowercase IUPAC for r. multi-base subs (closes #170) ([#367](https://github.com/fulcrumgenomics/ferro-hgvs/pull/367))
- *(normalize)* intronic + boundary normalization with NG/NC-parent inputs (closes #332) ([#364](https://github.com/fulcrumgenomics/ferro-hgvs/pull/364))
- *(parser,error-handling)* wire W3019 NonSpecMosaicForm for nested + [a/b] mosaic (closes #281) ([#298](https://github.com/fulcrumgenomics/ferro-hgvs/pull/298))
- *(normalize)* clamp 3'-rule shuffle at CDS↔UTR axis boundary (closes #337) ([#343](https://github.com/fulcrumgenomics/ferro-hgvs/pull/343))
- *(normalize)* strip explicit deleted-sequence from delins (closes #338) ([#344](https://github.com/fulcrumgenomics/ferro-hgvs/pull/344))
- *(parser,error-handling)* mito heteroplasmy prose diagnostic + W3017/W3018 SVAs (closes #278) ([#295](https://github.com/fulcrumgenomics/ferro-hgvs/pull/295))
- *(normalize)* bail gracefully when ref window mismatches HGVS span (closes #339) ([#345](https://github.com/fulcrumgenomics/ferro-hgvs/pull/345))
- *(parser,display)* preserve parens on inner repeat-count range insN[(150_180)] (closes #285) ([#306](https://github.com/fulcrumgenomics/ferro-hgvs/pull/306))
- *(normalize)* set RefSeqMismatch.corrected to reflect actual correction (closes #280) ([#297](https://github.com/fulcrumgenomics/ferro-hgvs/pull/297))
- *(normalize)* direction-aware ins→dup for homopolymers (closes #340 subgroup (i)) ([#346](https://github.com/fulcrumgenomics/ferro-hgvs/pull/346))
- *(reference)* synthesize transcript bases via cdot exon-alignment when FASTA version is missing (closes #331) ([#341](https://github.com/fulcrumgenomics/ferro-hgvs/pull/341))
- *(parser)* accept p.(0) whole-entity predicted no-protein (closes #289) ([#307](https://github.com/fulcrumgenomics/ferro-hgvs/pull/307))
- *(parser)* chain-expand mito = arms to inherit position, not synthesize m.1= (closes #284) ([#302](https://github.com/fulcrumgenomics/ferro-hgvs/pull/302))
- *(parser)* route [0] inside protein compact trans-allele to NoProtein (closes #277) ([#294](https://github.com/fulcrumgenomics/ferro-hgvs/pull/294))
- *(display)* canonicalize r. emission T→u in repeat units and bases (closes #276) ([#293](https://github.com/fulcrumgenomics/ferro-hgvs/pull/293))
- *(normalize)* extend codon-frame exception to r., 3+ chains, and sub+del (closes #275) ([#292](https://github.com/fulcrumgenomics/ferro-hgvs/pull/292))
- *(mock)* disambiguate genomic vs transcript lookups when ids collide (closes #311) ([#312](https://github.com/fulcrumgenomics/ferro-hgvs/pull/312))
- *(parser,error-handling)* wire W3021 ProteinBracketedAaInsertion (closes #290) ([#305](https://github.com/fulcrumgenomics/ferro-hgvs/pull/305))
- *(fasta)* route known-contig lookups in FastaProvider through the FASTA path (closes #315) ([#318](https://github.com/fulcrumgenomics/ferro-hgvs/pull/318))
- *(fasta)* require version-boundary equality in MmapFastaProvider transcript lookup (closes #314) ([#317](https://github.com/fulcrumgenomics/ferro-hgvs/pull/317))
- *(projector)* protein prediction for non-RefSeq transcripts + drop (GENE) selector from p. Display ([#313](https://github.com/fulcrumgenomics/ferro-hgvs/pull/313))
- *(validation)* wire W3016 LengthMismatch soft-validation warning (closes #81 K1 follow-up) ([#272](https://github.com/fulcrumgenomics/ferro-hgvs/pull/272))
- *(error-handling)* extend W4001 SwappedPositions to offset and *N markers (closes #81 L2 follow-up to #264) ([#271](https://github.com/fulcrumgenomics/ferro-hgvs/pull/271))
- *(spdi)* preserve inv and m. coord system on SPDI -> HGVS (closes #81 K1 follow-up) ([#274](https://github.com/fulcrumgenomics/ferro-hgvs/pull/274))
- *(error-handling)* wire W1004 MixedCaseEditType warning into preprocessor (closes #81 L2 follow-up) ([#273](https://github.com/fulcrumgenomics/ferro-hgvs/pull/273))
- *(error-handling)* wire W4001 SwappedPositions warning into preprocessor (closes #81 L2 remaining) ([#264](https://github.com/fulcrumgenomics/ferro-hgvs/pull/264))
- *(parser,display)* predicted-edit wrapper `c.(<pos><edit>)` on every NA coord (closes #241) ([#242](https://github.com/fulcrumgenomics/ferro-hgvs/pull/242))
- *(parser,display)* single uncertain position range c.(a_b)<edit> (closes #237) ([#238](https://github.com/fulcrumgenomics/ferro-hgvs/pull/238))
- *(protein)* extension Display emits extTerN to match Frameshift Ter canonicalization (closes #224) ([#225](https://github.com/fulcrumgenomics/ferro-hgvs/pull/225))
- *(parser)* accept `?<edit>` on g./r./n. and `?A>G` / `?=` (closes #239) ([#240](https://github.com/fulcrumgenomics/ferro-hgvs/pull/240))
- *(normalize)* read dup bases from reference (single-pass idempotency for dup with mismatched stated-ref) (closes #219) ([#230](https://github.com/fulcrumgenomics/ferro-hgvs/pull/230))
- *(vcf)* classify protein consequences by enum variant, not substring (closes #228) ([#229](https://github.com/fulcrumgenomics/ferro-hgvs/pull/229))
- *(normalize)* wire window-based normalization for m. variants ([#210](https://github.com/fulcrumgenomics/ferro-hgvs/pull/210)) ([#213](https://github.com/fulcrumgenomics/ferro-hgvs/pull/213))
- *(normalize)* exempt introns and UTR from codon-frame gate ([#209](https://github.com/fulcrumgenomics/ferro-hgvs/pull/209)) ([#211](https://github.com/fulcrumgenomics/ferro-hgvs/pull/211))
- *(parser)* depth-aware close-bracket lookup in nested HGVS brackets (closes #207) ([#208](https://github.com/fulcrumgenomics/ferro-hgvs/pull/208))
- *(normalize)* decompose delins to multi-sub canonical form when sub > delins applies (closes #165) ([#206](https://github.com/fulcrumgenomics/ferro-hgvs/pull/206))

### Other

- *(benchmark)* refresh measured performance tables for 0.7.0 ([#1006](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1006))
- *(spec-fixture)* classify the remaining spec-audit parse gaps ([#993](https://github.com/fulcrumgenomics/ferro-hgvs/pull/993))
- *(spec-fixture)* reassemble spot-split variants, drop harvest noise ([#990](https://github.com/fulcrumgenomics/ferro-hgvs/pull/990))
- *(conformance)* regenerate stale mutalyzer failure-patterns summary ([#998](https://github.com/fulcrumgenomics/ferro-hgvs/pull/998))
- *(conformance)* reclassify NM_017668.3 c.CAA[5] repeat as accepted_divergence, not known_bug ([#920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/920)) ([#986](https://github.com/fulcrumgenomics/ferro-hgvs/pull/986))
- *(conformance)* reclassify c.235_237delinsTAT as spec_citation, not known_bug ([#920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/920)) ([#984](https://github.com/fulcrumgenomics/ferro-hgvs/pull/984))
- *(conformance)* re-bless mutalyzer baseline-failures — all ledgers empty (closes #326) ([#982](https://github.com/fulcrumgenomics/ferro-hgvs/pull/982))
- *(reference)* defer Ensembl cdot load so RefSeq startup skips the merge ([#965](https://github.com/fulcrumgenomics/ferro-hgvs/pull/965))
- *(conformance)* document the dormant proj-near-disc + NR2E3 rows ([#946](https://github.com/fulcrumgenomics/ferro-hgvs/pull/946)) ([#949](https://github.com/fulcrumgenomics/ferro-hgvs/pull/949))
- *(conformance)* make the whole-CDS-deletion p.0?/p.(Met1?) divergence terminal ([#945](https://github.com/fulcrumgenomics/ferro-hgvs/pull/945)) ([#948](https://github.com/fulcrumgenomics/ferro-hgvs/pull/948))
- *(conformance)* reclassify the 8 NG_007485.1 noncoding rows as spec_citation ([#921](https://github.com/fulcrumgenomics/ferro-hgvs/pull/921)) ([#934](https://github.com/fulcrumgenomics/ferro-hgvs/pull/934))
- *(normalize)* correct misleading circular o. normalization comments ([#957](https://github.com/fulcrumgenomics/ferro-hgvs/pull/957))
- *(benchmark)* refresh comparator pins and record verified provenance ([#890](https://github.com/fulcrumgenomics/ferro-hgvs/pull/890)) ([#941](https://github.com/fulcrumgenomics/ferro-hgvs/pull/941))
- *(conformance)* reclassify copy-range delins rows as spec_citation ([#922](https://github.com/fulcrumgenomics/ferro-hgvs/pull/922)) ([#937](https://github.com/fulcrumgenomics/ferro-hgvs/pull/937))
- *(conformance)* re-home stale-tracker annotations to live issues ([#912](https://github.com/fulcrumgenomics/ferro-hgvs/pull/912)) ([#924](https://github.com/fulcrumgenomics/ferro-hgvs/pull/924))
- *(conformance)* disposition residual infos rows; reclassify NM_017668.3 exon-junction rows as spec_citation (#908, #918) ([#919](https://github.com/fulcrumgenomics/ferro-hgvs/pull/919))
- *(conformance)* disposition the 3 NG_008939 stop-insertion coding_protein rows (#911/#913 follow-up) ([#927](https://github.com/fulcrumgenomics/ferro-hgvs/pull/927))
- *(conformance)* reconcile the redundant cis-= allele axis inconsistency (#912 action 4) ([#928](https://github.com/fulcrumgenomics/ferro-hgvs/pull/928))
- *(conformance)* accept hgvs-rs-projection form-currency divergences ([#925](https://github.com/fulcrumgenomics/ferro-hgvs/pull/925)) ([#930](https://github.com/fulcrumgenomics/ferro-hgvs/pull/930))
- *(tool-support)* label the support matrix for the ferro 0.7.0 release ([#914](https://github.com/fulcrumgenomics/ferro-hgvs/pull/914))
- *(conformance)* harden the reference-identity signature ([#905](https://github.com/fulcrumgenomics/ferro-hgvs/pull/905)) ([#910](https://github.com/fulcrumgenomics/ferro-hgvs/pull/910))
- *(conformance)* accepted-divergence annotation sweep for the mutalyzer parity corpus ([#861](https://github.com/fulcrumgenomics/ferro-hgvs/pull/861)) ([#909](https://github.com/fulcrumgenomics/ferro-hgvs/pull/909))
- *(conformance)* correct insCAT protein_description; demote the 5'UTR-insertion divergence ([#891](https://github.com/fulcrumgenomics/ferro-hgvs/pull/891))
- *(benchmark)* record comparator provenance in the mutalyzer-normalize corpus header ([#890](https://github.com/fulcrumgenomics/ferro-hgvs/pull/890))
- *(conformance)* let accepted_rejection disposition an empty-projection Err ([#903](https://github.com/fulcrumgenomics/ferro-hgvs/pull/903))
- *(conformance)* demote stale known_bug #487 for LRG_303t1 out-of-phase insertion ([#906](https://github.com/fulcrumgenomics/ferro-hgvs/pull/906))
- *(project)* untangle project_to_genomic_nc — extract the coordinate map and #785 gate ([#868](https://github.com/fulcrumgenomics/ferro-hgvs/pull/868)) ([#902](https://github.com/fulcrumgenomics/ferro-hgvs/pull/902))
- *(conformance)* bless genomic-axis rejections via multi-axis accepted_rejection (#870 follow-up)
- *(conformance)* route genomic axis through the user-facing project_variant path ([#870](https://github.com/fulcrumgenomics/ferro-hgvs/pull/870))
- *(normalize)* make normalize_repeat the single source of truth for ins→repeat canonicalization ([#866](https://github.com/fulcrumgenomics/ferro-hgvs/pull/866)) ([#897](https://github.com/fulcrumgenomics/ferro-hgvs/pull/897))
- remove machine-local paths from committed example and runbook ([#889](https://github.com/fulcrumgenomics/ferro-hgvs/pull/889))
- *(conformance)* converge NG_-annotated transcript version on RefSeqGene selectors ([#859](https://github.com/fulcrumgenomics/ferro-hgvs/pull/859)) ([#871](https://github.com/fulcrumgenomics/ferro-hgvs/pull/871))
- *(conformance)* accept ferro's decline of NM_017668.3 coords outside NG_009299.1 coverage (#853, #865) ([#877](https://github.com/fulcrumgenomics/ferro-hgvs/pull/877))
- *(nightly)* derive ng_hosted_transcripts + NG_ placements in the nightly reference ([#875](https://github.com/fulcrumgenomics/ferro-hgvs/pull/875)) ([#876](https://github.com/fulcrumgenomics/ferro-hgvs/pull/876))
- *(conformance)* cite ferro's spec-correct inserted-inversion resolution ([#854](https://github.com/fulcrumgenomics/ferro-hgvs/pull/854)) ([#878](https://github.com/fulcrumgenomics/ferro-hgvs/pull/878))
- *(conformance)* accept ferro's spec-correct refusal of non-standard/absent NG_ selectors ([#858](https://github.com/fulcrumgenomics/ferro-hgvs/pull/858)) ([#873](https://github.com/fulcrumgenomics/ferro-hgvs/pull/873))
- *(conformance)* rebless rna_description ledger and protein empty-projection pin ([#799](https://github.com/fulcrumgenomics/ferro-hgvs/pull/799)) ([#846](https://github.com/fulcrumgenomics/ferro-hgvs/pull/846))
- *(normalize)* resolve the cigar-aware normalization backlog ([#811](https://github.com/fulcrumgenomics/ferro-hgvs/pull/811)) ([#833](https://github.com/fulcrumgenomics/ferro-hgvs/pull/833))
- *(conformance)* allow multiple per-axis spec_citation annotations on a Case ([#830](https://github.com/fulcrumgenomics/ferro-hgvs/pull/830))
- *(conformance)* disposition bare-NP protein-framing rows as spec_overridden ([#826](https://github.com/fulcrumgenomics/ferro-hgvs/pull/826))
- *(spdi)* run SPDI roundtrip + rsID tests offline via synthetic fixture ([#825](https://github.com/fulcrumgenomics/ferro-hgvs/pull/825))
- *(python_helpers)* remove dead edit_type_from_debug {:?}-string matcher ([#829](https://github.com/fulcrumgenomics/ferro-hgvs/pull/829))
- *(conformance)* cover exon-junction rows in the genomic-axis gate ([#751](https://github.com/fulcrumgenomics/ferro-hgvs/pull/751)) ([#789](https://github.com/fulcrumgenomics/ferro-hgvs/pull/789))
- *(reference)* derive accession→assembly from NCBI assembly_report ([#716](https://github.com/fulcrumgenomics/ferro-hgvs/pull/716)) ([#788](https://github.com/fulcrumgenomics/ferro-hgvs/pull/788))
- *(conformance)* decrement protein empty-projection budget after #778 ([#786](https://github.com/fulcrumgenomics/ferro-hgvs/pull/786))
- *(project)* predict VSIR reverse-strand whole-exon frameshift ([#774](https://github.com/fulcrumgenomics/ferro-hgvs/pull/774))
- *(conformance)* re-bless empty-projection budget on canonical reference ([#775](https://github.com/fulcrumgenomics/ferro-hgvs/pull/775))
- *(conformance)* disposition pter/qter-in-c flank payloads; green the normalized axis ([#770](https://github.com/fulcrumgenomics/ferro-hgvs/pull/770))
- *(conformance)* accept coding_protein transcript-set enumeration divergences ([#768](https://github.com/fulcrumgenomics/ferro-hgvs/pull/768))
- *(conformance)* pin the empty-projection budget to its reference ([#767](https://github.com/fulcrumgenomics/ferro-hgvs/pull/767))
- *(conformance)* restore clobbered genomic #745 + add accepted_rejection ([#654](https://github.com/fulcrumgenomics/ferro-hgvs/pull/654)) ([#757](https://github.com/fulcrumgenomics/ferro-hgvs/pull/757))
- *(conformance)* disposition the normalized-axis Ok-mismatch divergences ([#654](https://github.com/fulcrumgenomics/ferro-hgvs/pull/654)) ([#755](https://github.com/fulcrumgenomics/ferro-hgvs/pull/755))
- *(conformance)* forward dropped methods through the ArcProvider wrapper (#726 follow-up) ([#754](https://github.com/fulcrumgenomics/ferro-hgvs/pull/754))
- *(conformance)* green errors axis — bucket accepted divergences ([#486](https://github.com/fulcrumgenomics/ferro-hgvs/pull/486)) ([#752](https://github.com/fulcrumgenomics/ferro-hgvs/pull/752))
- *(conformance)* reclassify #745 homopolymer repeat-contraction as accepted divergence ([#750](https://github.com/fulcrumgenomics/ferro-hgvs/pull/750))
- *(reference)* NCBI assembly_report parser + data-driven build inference ([#716](https://github.com/fulcrumgenomics/ferro-hgvs/pull/716)) ([#735](https://github.com/fulcrumgenomics/ferro-hgvs/pull/735))
- consolidate ~230 integration-test binaries into one crate ([#725](https://github.com/fulcrumgenomics/ferro-hgvs/pull/725))
- *(annotate-vcf)* annotate VCF records in place instead of cloning ([#698](https://github.com/fulcrumgenomics/ferro-hgvs/pull/698))
- *(annotate-vcf)* parallelize file-input record annotation behind -j/--workers ([#696](https://github.com/fulcrumgenomics/ferro-hgvs/pull/696))
- *(parse)* parallelize the batch CLI loop behind -j/--workers ([#688](https://github.com/fulcrumgenomics/ferro-hgvs/pull/688))
- *(vcf)* remove dead batch processors from src/vcf/batch.rs ([#699](https://github.com/fulcrumgenomics/ferro-hgvs/pull/699))
- *(annotate)* interval-index TranscriptDb::get_by_region ([#697](https://github.com/fulcrumgenomics/ferro-hgvs/pull/697))
- *(annotate-vcf)* buffer the output writer ([#695](https://github.com/fulcrumgenomics/ferro-hgvs/pull/695))
- *(python)* parallelize the batch API and release the GIL ([#691](https://github.com/fulcrumgenomics/ferro-hgvs/pull/691))
- *(normalize)* parallelize the batch CLI loop behind -j/--workers ([#687](https://github.com/fulcrumgenomics/ferro-hgvs/pull/687))
- *(pre-commit)* bump ruff hook to v0.15.12 to match CI ([#705](https://github.com/fulcrumgenomics/ferro-hgvs/pull/705))
- *(normalize)* keep `IntronicVariant.variant` machine-readable ([#685](https://github.com/fulcrumgenomics/ferro-hgvs/pull/685))
- *(conformance)* add an output-quality signal to the mutalyzer axis harness ([#651](https://github.com/fulcrumgenomics/ferro-hgvs/pull/651)) ([#675](https://github.com/fulcrumgenomics/ferro-hgvs/pull/675))
- *(conformance)* classify remaining #487 3′-shift audit rows ([#673](https://github.com/fulcrumgenomics/ferro-hgvs/pull/673))
- *(parser)* skip Unicode trim when input has no surrounding whitespace ([#659](https://github.com/fulcrumgenomics/ferro-hgvs/pull/659))
- *(readme)* honest benchmark framing — offline best-case + single-threaded parse ([#660](https://github.com/fulcrumgenomics/ferro-hgvs/pull/660))
- *(normalize)* pre-size runtime cdot maps when rebuilding from the archive ([#666](https://github.com/fulcrumgenomics/ferro-hgvs/pull/666))
- *(normalize)* use FxHashMap for the runtime cdot reference maps ([#665](https://github.com/fulcrumgenomics/ferro-hgvs/pull/665))
- *(normalize)* avoid per-variant allocations and walks in input preprocessing ([#664](https://github.com/fulcrumgenomics/ferro-hgvs/pull/664))
- *(normalize)* skip reference teardown at process exit ([#663](https://github.com/fulcrumgenomics/ferro-hgvs/pull/663))
- *(parser)* fast-path identity edits (protein and nucleotide) ([#657](https://github.com/fulcrumgenomics/ferro-hgvs/pull/657))
- *(parser)* fast-path protein missense substitutions ([#650](https://github.com/fulcrumgenomics/ferro-hgvs/pull/650))
- *(parser)* fast-path intronic and UTR coding substitutions ([#643](https://github.com/fulcrumgenomics/ferro-hgvs/pull/643))
- *(benchmark)* refresh measured tables with startup-excluded normalize timing ([#648](https://github.com/fulcrumgenomics/ferro-hgvs/pull/648))
- *(parser)* fast-path plain deletions and duplications ([#641](https://github.com/fulcrumgenomics/ferro-hgvs/pull/641))
- *(conformance)* annotate coding non-codon-aligned repeats as spec-overrides ([#640](https://github.com/fulcrumgenomics/ferro-hgvs/pull/640))
- *(parser)* eliminate per-parse accession allocations (~27% faster parsing) ([#639](https://github.com/fulcrumgenomics/ferro-hgvs/pull/639))
- *(conformance)* annotate + hermetic-gate the biocommons normalize corpus (closes #325) ([#636](https://github.com/fulcrumgenomics/ferro-hgvs/pull/636))
- *(reference)* lazily load the GRCh37 secondary cdot until first use ([#635](https://github.com/fulcrumgenomics/ferro-hgvs/pull/635))
- remove committed machine-specific manifest path from benchmark tests ([#634](https://github.com/fulcrumgenomics/ferro-hgvs/pull/634))
- *(conformance)* disposition the residual signal — fully classify the hgvs-rs dashboard (burn-down PR 3) ([#620](https://github.com/fulcrumgenomics/ferro-hgvs/pull/620))
- add CodeRabbit config (assertive profile + path instructions) ([#627](https://github.com/fulcrumgenomics/ferro-hgvs/pull/627))
- *(conformance)* quarantine hgvs-rs data-source divergences so red = real (burn-down PR 2) ([#614](https://github.com/fulcrumgenomics/ferro-hgvs/pull/614))
- *(conformance)* version-insensitive comparator + selection-coverage for the hgvs-rs corpus ([#607](https://github.com/fulcrumgenomics/ferro-hgvs/pull/607))
- fix broken DOI badge in README ([#574](https://github.com/fulcrumgenomics/ferro-hgvs/pull/574))
- *(conformance)* accept compound-allele SHUFFLE_APPLIED as a tracked divergence ([#499](https://github.com/fulcrumgenomics/ferro-hgvs/pull/499)) ([#572](https://github.com/fulcrumgenomics/ferro-hgvs/pull/572))
- *(conformance)* reject position-past-end on the errors axis (#486 EOUTOFBOUNDARY) ([#567](https://github.com/fulcrumgenomics/ferro-hgvs/pull/567))
- *(reference)* precompute has_genomic_data instead of scanning the index ([#610](https://github.com/fulcrumgenomics/ferro-hgvs/pull/610))
- *(cli)* format the normalized variant once per output line ([#611](https://github.com/fulcrumgenomics/ferro-hgvs/pull/611))
- *(conformance)* reject ISCN-only / superseded `::` rearrangement forms ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#608](https://github.com/fulcrumgenomics/ferro-hgvs/pull/608))
- perf-regression gate (cdot load-source structural test + nightly startup budget) ([#606](https://github.com/fulcrumgenomics/ferro-hgvs/pull/606))
- *(normalize)* resolve intronic transcripts by accession, not a cloned variant ([#582](https://github.com/fulcrumgenomics/ferro-hgvs/pull/582))
- *(reference)* reuse FASTA file handles with positioned reads (~1.9x normalize) ([#570](https://github.com/fulcrumgenomics/ferro-hgvs/pull/570))
- *(parser)* route the default parse_hgvs through the fast path (~1.7x on ClinVar) ([#562](https://github.com/fulcrumgenomics/ferro-hgvs/pull/562))
- *(reference)* memoize resolved transcripts in MultiFastaProvider ([#578](https://github.com/fulcrumgenomics/ferro-hgvs/pull/578))
- *(reference)* faster, deterministic FASTA index build (FxHashMap + highest-version fallback) ([#588](https://github.com/fulcrumgenomics/ferro-hgvs/pull/588))
- *(project)* add c→p, c↔n, ref-AA cache, and batch projection microbenchmarks (Cycle 1c-ii) ([#599](https://github.com/fulcrumgenomics/ferro-hgvs/pull/599))
- *(conformance)* import the hgvs-rs projection corpus as a third oracle (Cycle 1c-i) ([#598](https://github.com/fulcrumgenomics/ferro-hgvs/pull/598))
- *(reference)* binary cache for the GFF/GTF TranscriptDb (annotate-vcf ~3s → ~0.6s) ([#593](https://github.com/fulcrumgenomics/ferro-hgvs/pull/593))
- *(reference)* rkyv archive for the cdot cache (zero-copy load, ~0.37s off startup) ([#591](https://github.com/fulcrumgenomics/ferro-hgvs/pull/591))
- *(reference)* version + self-heal the cdot bincode cache (fix silent JSON fallback) ([#585](https://github.com/fulcrumgenomics/ferro-hgvs/pull/585))
- *(tests)* stop committing the generated spec-normalization fixture ([#586](https://github.com/fulcrumgenomics/ferro-hgvs/pull/586))
- *(normalize)* two no-functional-change normalize hot-path wins (shuffle-info + protein-length) ([#553](https://github.com/fulcrumgenomics/ferro-hgvs/pull/553))
- *(bench)* broaden microbenchmark suite + tag-sweep runner (perf-effort foundation) ([#554](https://github.com/fulcrumgenomics/ferro-hgvs/pull/554))
- fix nightly reference-aware tests (`ferro prepare --output-dir`) + add nightly badge ([#551](https://github.com/fulcrumgenomics/ferro-hgvs/pull/551))
- *(biocommons)* mirror the improvement conformance disposition ([#503](https://github.com/fulcrumgenomics/ferro-hgvs/pull/503)) ([#540](https://github.com/fulcrumgenomics/ferro-hgvs/pull/540))
- *(mutalyzer)* demote c.41A>C — substitution is a missense, not a frameshift ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#516](https://github.com/fulcrumgenomics/ferro-hgvs/pull/516))
- *(mutalyzer)* add improvement disposition; demote NG_ gene-symbol selector divergences ([#500](https://github.com/fulcrumgenomics/ferro-hgvs/pull/500)) ([#501](https://github.com/fulcrumgenomics/ferro-hgvs/pull/501))
- *(mutalyzer)* known_bug disposition + XPASS detection ([#478](https://github.com/fulcrumgenomics/ferro-hgvs/pull/478)) ([#489](https://github.com/fulcrumgenomics/ferro-hgvs/pull/489))
- *(biocommons)* annotation + XPASS conformance harness (#478 pillars 1-2) ([#484](https://github.com/fulcrumgenomics/ferro-hgvs/pull/484))
- *(mutalyzer-normalize)* demote ins[...] cases fixed by #333 (49 rows) ([#479](https://github.com/fulcrumgenomics/ferro-hgvs/pull/479))
- *(test)* prune stale NM_000051.3:c.1_2insCA from biocommons baseline-failures ([#477](https://github.com/fulcrumgenomics/ferro-hgvs/pull/477))
- trans-allele dedup, mixed phase reject, RNA whole-entity in brackets (closes #396) ([#425](https://github.com/fulcrumgenomics/ferro-hgvs/pull/425))
- *(protein)* pin stop-codon glyph canonicalization ([#453](https://github.com/fulcrumgenomics/ferro-hgvs/pull/453)) ([#454](https://github.com/fulcrumgenomics/ferro-hgvs/pull/454))
- *(protein)* pin construct-boundary canonicalization probes (F13) ([#458](https://github.com/fulcrumgenomics/ferro-hgvs/pull/458))
- *(allele)* pin allele-bracket grammar corner cases (F12) ([#457](https://github.com/fulcrumgenomics/ferro-hgvs/pull/457))
- *(parser)* pin strict-reject mode for W3011 DelSizeSuffix ([#447](https://github.com/fulcrumgenomics/ferro-hgvs/pull/447)) ([#459](https://github.com/fulcrumgenomics/ferro-hgvs/pull/459))
- *(parser)* fix red main — drop obsolete dupins-accept probe (contradicts #445) ([#463](https://github.com/fulcrumgenomics/ferro-hgvs/pull/463))
- *(parser)* pin input-hygiene & spec-rejected parser behavior (F1+F2) ([#444](https://github.com/fulcrumgenomics/ferro-hgvs/pull/444))
- *(spec)* pin miscellaneous coverage corners (F6+F7+F10-11+F14-18) ([#456](https://github.com/fulcrumgenomics/ferro-hgvs/pull/456))
- omnibus follow-up: loader + diagnostics + ci polish (closes #397) ([#426](https://github.com/fulcrumgenomics/ferro-hgvs/pull/426))
- *(parser)* pin strict-whitespace parse mode ([#449](https://github.com/fulcrumgenomics/ferro-hgvs/pull/449)) ([#455](https://github.com/fulcrumgenomics/ferro-hgvs/pull/455))
- *(rna)* add r.↔c. consistency + RNA splicing marker coverage ([#442](https://github.com/fulcrumgenomics/ferro-hgvs/pull/442))
- *(parser)* drop dead 0? clause in protein trans-allele shorthand (closes #424) ([#431](https://github.com/fulcrumgenomics/ferro-hgvs/pull/431))
- *(benchmark)* wire supplemental_fasta into manifest after patterns fetch ([#417](https://github.com/fulcrumgenomics/ferro-hgvs/pull/417))
- *(test)* resync biocommons-normalize baseline-failures with manifest-mode xfail ([#405](https://github.com/fulcrumgenomics/ferro-hgvs/pull/405))
- *(normalize)* 3'-rule completeness — repeat-depth audit + #343 deferred symmetry tests (closes #391) ([#409](https://github.com/fulcrumgenomics/ferro-hgvs/pull/409))
- *(test)* sync biocommons-normalize baseline-failures with current xfail ([#386](https://github.com/fulcrumgenomics/ferro-hgvs/pull/386))
- *(normalize)* pin delins → dup canonicalization over downstream tract (closes #382) ([#384](https://github.com/fulcrumgenomics/ferro-hgvs/pull/384))
- *(python)* migrate PyO3 bindings to abi3 stable ABI ([#381](https://github.com/fulcrumgenomics/ferro-hgvs/pull/381))
- *(project)* rebuild VariantProjector hot paths ([#361](https://github.com/fulcrumgenomics/ferro-hgvs/pull/361))
- *(error-handling)* retire W4002 PositionZero duplicate identity (closes #269) ([#370](https://github.com/fulcrumgenomics/ferro-hgvs/pull/370))
- *(parser,display)* pin c./r. predicted-edit wrapper Display parity (closes #300) ([#372](https://github.com/fulcrumgenomics/ferro-hgvs/pull/372))
- *(biocommons)* import normalize test corpus (88 cases, stacked on #323) ([#324](https://github.com/fulcrumgenomics/ferro-hgvs/pull/324))
- update stale Display test expectations to HGVS spec-canonical forms ([#366](https://github.com/fulcrumgenomics/ferro-hgvs/pull/366))
- *(contributing)* defer CHANGELOG.md to release-plz ([#360](https://github.com/fulcrumgenomics/ferro-hgvs/pull/360))
- *(reference)* impl Default for Transcript to enable spread-update fixtures ([#322](https://github.com/fulcrumgenomics/ferro-hgvs/pull/322))
- *(mutalyzer)* import normalize test corpus (320 cases, 8 axes) ([#323](https://github.com/fulcrumgenomics/ferro-hgvs/pull/323))
- *(normalize)* pin r.-axis (transcript-1-relative) convention (closes #291) ([#304](https://github.com/fulcrumgenomics/ferro-hgvs/pull/304))
- *(convert)* pin c. → r. reference-driven conversion across edit shapes (closes #283) ([#301](https://github.com/fulcrumgenomics/ferro-hgvs/pull/301))
- *(audit)* HGVS ↔ SPDI round-trip coverage matrix (closes #81 K1 remaining) ([#260](https://github.com/fulcrumgenomics/ferro-hgvs/pull/260))
- *(audit)* HGVS ↔ VCF conversion surface (closes #81 K2) ([#262](https://github.com/fulcrumgenomics/ferro-hgvs/pull/262))
- *(audit)* minus-strand intronic position ordering at parse/Display (closes #81 J3) ([#258](https://github.com/fulcrumgenomics/ferro-hgvs/pull/258))
- *(audit)* n. upstream/downstream * and - markers canonical output (closes #81 J2) ([#256](https://github.com/fulcrumgenomics/ferro-hgvs/pull/256))
- *(audit)* LRG references full round-trip across edit types (closes #81 I1) ([#248](https://github.com/fulcrumgenomics/ferro-hgvs/pull/248))
- *(audit)* boundary-spanning del/delins ranges across UTR/CDS (closes #81 J1) ([#254](https://github.com/fulcrumgenomics/ferro-hgvs/pull/254))
- *(audit)* chromosome alias input + canonical output (closes #81 I4) ([#252](https://github.com/fulcrumgenomics/ferro-hgvs/pull/252))
- *(audit)* versioned vs unversioned accession policy (closes #81 I2) ([#250](https://github.com/fulcrumgenomics/ferro-hgvs/pull/250))
- *(mito)* extended heteroplasmy / prose-shape audit and policy doc (closes #235) ([#236](https://github.com/fulcrumgenomics/ferro-hgvs/pull/236))
- *(rna)* pin r./c. canonicalization parity across edit shapes (closes #233) ([#234](https://github.com/fulcrumgenomics/ferro-hgvs/pull/234))
- *(rna)* dedicated lowercase enforcement coverage across edit types (closes #231) ([#232](https://github.com/fulcrumgenomics/ferro-hgvs/pull/232))
- protein canonicalization audit umbrella for #81 § D — D1 / D3 / D4 / D8 (closes #226) ([#227](https://github.com/fulcrumgenomics/ferro-hgvs/pull/227))
- 3+ variants in one allele bracket where merging does not apply (closes #221) ([#222](https://github.com/fulcrumgenomics/ferro-hgvs/pull/222))
- dedicated coverage for mixed-accession compound alleles (closes #218) ([#220](https://github.com/fulcrumgenomics/ferro-hgvs/pull/220))

- *(parser)* accept `?con<src>` and `?copy<N>` at unknown position, consistent with `?del`/`?dup`/`?ins`/`?inv` (closes #286)
- *(fasta)* require version-boundary equality in `MmapFastaProvider::get_transcript`'s unversioned-prefix fallback (closes #314)
- *(fasta)* route `FastaProvider::get_sequence` for known contigs through the FASTA path so a transcript registered with a chromosome-colliding id no longer wins over the genomic index (closes #315)

### Fixed

- *(project)* decline rather than emit invalid HGVS when an `NG_`/`LRG_`-parented projection cannot be re-anchored into the parent's own frame: with no chromosomal placement (or an endpoint outside the placed span / an uncertain boundary), `project_to_genomic` now returns `UnsupportedProjection`/`InvalidCoordinates` instead of stamping a chromosome (`NC_`) coordinate under the parent accession. The cross-isoform enumerate path (`project_variant_all`, #646) degrades gracefully — dropping only the unframable genomic axis while keeping the parent-framed coding/protein forms (closes #655).
- *(reference)* make cdot base→version fallback deterministic; previously, when a cdot file contained multiple versions of the same base accession (e.g. `NM_000088.3` and `NM_000088.4`), `base_to_versioned` could resolve to any version across runs due to `HashMap` iteration order, silently shifting CDS coordinates for callers requesting an absent version (closes #583).
- *(edit)* drop brackets on single-payload `InsertedSequence::Complex` `Display` so e.g. `delins[78185355_78199419inv]` round-trips as `delins78185355_78199419inv`, matching HGVS v21 (`DNA/insertion.md:22`, `DNA/inversion.md:39`).
- *(normalize)* extend HGVS codon-frame exception beyond the two-`c.`-SNV case (closes #275, follow-up to #79 / #104): the spec's "two variants separated by one nucleotide, together affecting one amino acid" carve-out now also covers (1) `r.` coding regions, (2) chains of 3+ SNVs where a strict-adjacency merge leaves `prev_a` as a multi-base delins, and (3) `sub`+`del` and `del`+`sub` pairs separated by one unchanged nucleotide.

## [0.6.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.5.0...v0.6.0) - 2026-05-16

### Added

- *(project)* compound alleles + project_all + project_normalized (PR 3: closes #200) ([#204](https://github.com/fulcrumgenomics/ferro-hgvs/pull/204))
- *(cli)* add build-transcript subcommand for FASTA + CDS construction (closes #184) ([#199](https://github.com/fulcrumgenomics/ferro-hgvs/pull/199))
- *(project)* protein consequence for CDS indels (PR 2: indel p. nomenclature) ([#203](https://github.com/fulcrumgenomics/ferro-hgvs/pull/203))
- *(project)* variant-level projection g. -> c./p. (PR 1: foundation + substitutions) ([#202](https://github.com/fulcrumgenomics/ferro-hgvs/pull/202))
- *(loader)* Phase 4 - FASTA validation + complete migration off shims ([#190](https://github.com/fulcrumgenomics/ferro-hgvs/pull/190)) ([#196](https://github.com/fulcrumgenomics/ferro-hgvs/pull/196))
- *(loader,cli)* Phase 3 - diagnostics registry + CLI flags + cross-format tests ([#190](https://github.com/fulcrumgenomics/ferro-hgvs/pull/190)) ([#195](https://github.com/fulcrumgenomics/ferro-hgvs/pull/195))
- *(loader)* Phase 2 - CDS phase/start_codon/stop_codon, UTR-merge, MANE wiring ([#190](https://github.com/fulcrumgenomics/ferro-hgvs/pull/190)) ([#194](https://github.com/fulcrumgenomics/ferro-hgvs/pull/194))

### Fixed

- *(loader)* rewrite GFF/GTF loader pipeline, close #183 ([#191](https://github.com/fulcrumgenomics/ferro-hgvs/pull/191))
- *(mock)* accept version/genome_build metadata keys in from_json (closes #185) ([#198](https://github.com/fulcrumgenomics/ferro-hgvs/pull/198))
- *(normalize)* group consecutive sub-flanks of inv-split into delins ([#182](https://github.com/fulcrumgenomics/ferro-hgvs/pull/182)) ([#188](https://github.com/fulcrumgenomics/ferro-hgvs/pull/188))
- *(normalize)* apply 3' rule to allele-merged del/dup/ins ([#180](https://github.com/fulcrumgenomics/ferro-hgvs/pull/180)) ([#187](https://github.com/fulcrumgenomics/ferro-hgvs/pull/187))

### Other

- *(loader)* Phase 5 - swap parser internals for noodles ([#190](https://github.com/fulcrumgenomics/ferro-hgvs/pull/190)) ([#197](https://github.com/fulcrumgenomics/ferro-hgvs/pull/197))
- *(spec)* refresh vendored hgvs-nomenclature spec snapshot ([#193](https://github.com/fulcrumgenomics/ferro-hgvs/pull/193))

### Added

- *(loader)* unified `load_annotations` entry point with `LoaderConfig`/`LoaderReport`; auto-detects GFF3 vs GTF by extension and content (#191, #194, #195)
- *(loader)* exon-derivation ladder closes single-exon GFF3 without `exon` records (#183 fix) (#191)
- *(loader)* phase-aware CDS bounds with `start_codon` / `stop_codon` precedence; GFF3 input now extends `cds_end` to include the stop codon, matching ferro's downstream convention (fixes a latent off-by-one-codon bug in GFF3 protein conversion) (#194)
- *(loader)* UTR-merge ladder step: synthesizes exons from UTR + CDS when no explicit `exon` records (#194)
- *(loader)* `gene_symbol`, `mane_status`, `refseq_match`, `ensembl_match` extracted from attributes (#194)
- *(loader)* optional FASTA-aware validation: CDS length mod 3 and start codon (ATG/CTG/GTG/TTG) checks when `--fasta` is supplied
- *(error_handling)* 13 loader diagnostic codes (E-LOAD-*, W-LOAD-*) registered for `ferro explain` (#195)
- *(cli)* `ferro convert-gff --strict / --silent / --no-validate-fasta / --diagnostics-json` flags (#195)
- *(strand)* `Strand::Unknown` variant for GFF3 `.` and `?` strand values; transcripts with unknown strand are dropped at load with `E-LOAD-103` (#191)

### Changed

- *(strand)* `Strand`, `LoaderConfig`, and `LoaderReport` are marked `#[non_exhaustive]` for forward compatibility (#191)
- *(loader)* `load_gff3` and `load_gtf` have been removed; use `load_annotations` instead

## [0.5.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.4.1...v0.5.0) - 2026-05-12

### Added

- *(parser)* spec compact-mosaic form across sub/del/dup ([#133](https://github.com/fulcrumgenomics/ferro-hgvs/pull/133)) ([#153](https://github.com/fulcrumgenomics/ferro-hgvs/pull/153))
- *(normalize)* canonicalize con to delins per SVD-WG009 (#81 H1) ([#142](https://github.com/fulcrumgenomics/ferro-hgvs/pull/142))
- *(rna)* support r.spl, r.spl?, r.(spl), r.(spl?) splicing markers (#81 E2) ([#134](https://github.com/fulcrumgenomics/ferro-hgvs/pull/134))

### Fixed

- *(spdi)* emit inversion and repeat as SPDI delins ([#159](https://github.com/fulcrumgenomics/ferro-hgvs/pull/159))
- *(error-handling)* error codes for spec-mandated input rejections ([#115](https://github.com/fulcrumgenomics/ferro-hgvs/pull/115)) ([#152](https://github.com/fulcrumgenomics/ferro-hgvs/pull/152))
- *(spdi)* reference-aware HGVS→SPDI for del/dup/delins ([#158](https://github.com/fulcrumgenomics/ferro-hgvs/pull/158))
- *(parser)* preserve explicit deleted sequence in delins round-trip ([#154](https://github.com/fulcrumgenomics/ferro-hgvs/pull/154))
- *(error-handling)* soft-warn non-canonical input forms at parse time ([#127](https://github.com/fulcrumgenomics/ferro-hgvs/pull/127)) ([#151](https://github.com/fulcrumgenomics/ferro-hgvs/pull/151))
- *(ci)* drop redundant `#[allow(dead_code)]` on `mod common;` ([#178](https://github.com/fulcrumgenomics/ferro-hgvs/pull/178))
- *(spdi)* accept c./n./r./m. variants in HGVS->SPDI conversion ([#157](https://github.com/fulcrumgenomics/ferro-hgvs/pull/157))
- *(normalize)* apply 3' rule across phase-mismatched cyclic rotations for ins ([#132](https://github.com/fulcrumgenomics/ferro-hgvs/pull/132)) ([#155](https://github.com/fulcrumgenomics/ferro-hgvs/pull/155))
- *(error-handling)* surface deprecated stop-codon and frameshift forms as soft-warns ([#125](https://github.com/fulcrumgenomics/ferro-hgvs/pull/125)) ([#150](https://github.com/fulcrumgenomics/ferro-hgvs/pull/150))
- *(error-handling)* wire W1001/W1002/W3001 soft-validation warnings ([#124](https://github.com/fulcrumgenomics/ferro-hgvs/pull/124)) ([#149](https://github.com/fulcrumgenomics/ferro-hgvs/pull/149))
- *(spdi)* recover dup form on SPDI→HGVS for duplicated insertions ([#156](https://github.com/fulcrumgenomics/ferro-hgvs/pull/156))
- *(parser)* extend unknown-phase (;) support to g./n./m./o./p. coord systems ([#123](https://github.com/fulcrumgenomics/ferro-hgvs/pull/123)) ([#148](https://github.com/fulcrumgenomics/ferro-hgvs/pull/148))
- *(normalize)* detect cis-allele edits with coincident bounds (#81 A8) ([#147](https://github.com/fulcrumgenomics/ferro-hgvs/pull/147))
- *(parser)* soft-warn embedded whitespace and zero-width chars ([#128](https://github.com/fulcrumgenomics/ferro-hgvs/pull/128)) ([#145](https://github.com/fulcrumgenomics/ferro-hgvs/pull/145))
- *(normalize)* preserve r.*N UTR flag and translate via cds_end (closes #163) ([#164](https://github.com/fulcrumgenomics/ferro-hgvs/pull/164))
- *(parser)* correct LRG accession variant-type inference and compound-ref handling ([#122](https://github.com/fulcrumgenomics/ferro-hgvs/pull/122)) ([#141](https://github.com/fulcrumgenomics/ferro-hgvs/pull/141))
- *(normalize)* recognize revcomp inv sub-spans within delins ([#166](https://github.com/fulcrumgenomics/ferro-hgvs/pull/166))
- *(normalize)* apply 3'-rule to merged cis-allele deletions (closes #161) ([#162](https://github.com/fulcrumgenomics/ferro-hgvs/pull/162))
- *(normalize)* rewrite revcomp delins as inversion (#81 A2) ([#109](https://github.com/fulcrumgenomics/ferro-hgvs/pull/109))
- *(normalize)* enforce HGVS c. codon-frame exception for repeat notation (#81 B1) ([#110](https://github.com/fulcrumgenomics/ferro-hgvs/pull/110))
- *(normalize)* rewrite empty-insert delins as del (#81 A3) ([#113](https://github.com/fulcrumgenomics/ferro-hgvs/pull/113))
- *(normalize)* degenerate substitution (ref==alt) -> identity (#81 A4) ([#111](https://github.com/fulcrumgenomics/ferro-hgvs/pull/111))

### Other

- *(parser)* triage failure expectations (#174 phase 2) ([#176](https://github.com/fulcrumgenomics/ferro-hgvs/pull/176))
- *(parser)* per-input failure-expectations framework (#174 phase 1) ([#175](https://github.com/fulcrumgenomics/ferro-hgvs/pull/175))
- drop ferro_version from HGVS spec fixture ([#177](https://github.com/fulcrumgenomics/ferro-hgvs/pull/177))
- *(test)* cut Test job runtime ~398s → ~20s ([#173](https://github.com/fulcrumgenomics/ferro-hgvs/pull/173))
- *(test)* correct misleading comment on trans-allele expanded-form test ([#167](https://github.com/fulcrumgenomics/ferro-hgvs/pull/167))
- *(normalize)* tag r. positive bases as Region::Rna (closes #168) ([#169](https://github.com/fulcrumgenomics/ferro-hgvs/pull/169))
- *(allele)* pin trans-phase round-trip across all coord systems and merge-barrier (#81 C1) ([#146](https://github.com/fulcrumgenomics/ferro-hgvs/pull/146))
- *(compound)* pin cross-reference / cross-coord compound round-trip and merge-barrier (#81 H2) ([#143](https://github.com/fulcrumgenomics/ferro-hgvs/pull/143))
- *(error-handling)* audit error codes against HGVS spec sections (#81 L1) ([#137](https://github.com/fulcrumgenomics/ferro-hgvs/pull/137))
- *(parser)* pin gene-selector round-trip end-to-end with Display preservation (#81 I3) ([#135](https://github.com/fulcrumgenomics/ferro-hgvs/pull/135))
- *(mito)* audit heteroplasmy notation; tracking #133 (#81 F2) ([#139](https://github.com/fulcrumgenomics/ferro-hgvs/pull/139))
- *(normalize)* pin RNA path + edge cases for A9 substitution-after-trim ([#114](https://github.com/fulcrumgenomics/ferro-hgvs/pull/114))
- *(protein)* pin p.? unknown-effect round-trip across allele forms and edge cases (#81 D7) ([#136](https://github.com/fulcrumgenomics/ferro-hgvs/pull/136))
- *(protein)* pin p.0 no-product round-trip and adjacent guards (#81 D6) ([#130](https://github.com/fulcrumgenomics/ferro-hgvs/pull/130))
- *(protein)* pin silent `=` round-trip across allele forms and edge cases (#81 D5) ([#131](https://github.com/fulcrumgenomics/ferro-hgvs/pull/131))
- *(mito)* audit m. coord-system parse + wraparound behavior; tracking #129 (#81 F1) ([#138](https://github.com/fulcrumgenomics/ferro-hgvs/pull/138))
- pin normalize() against HGVS v21.0 spec fixture (closes #84) ([#105](https://github.com/fulcrumgenomics/ferro-hgvs/pull/105))

### Added

- *(normalize)* detect cis-allele edits with coincident reference bounds (e.g. `g.[100G>A;100A>C]`) and emit an advisory `OVERLAP_CONFLICTING_EDITS` (`W5002`) warning. Variant output is preserved unchanged. Addresses [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A8.
- *(test)* `expected_warnings` field on the HGVS v21.0 spec-fixture row schema, pinning the warning set ferro emits per row.

### Fixed

- *(normalize)* apply 3' rule across phase-mismatched cyclic rotations for single-copy ins. When the inserted alt is a non-zero cyclic rotation of an adjacent reference repeat unit (e.g. `g.X_(X+1)insTG` against a `GT[3]` tract), shuffle's first-base check (`alt[0] == ref[ins_point]`) failed and the variant never moved. The new helper `insertion_to_duplication` mirrors `insertion_to_repeat`'s rotation iteration for the 1-copy case, so the variant now canonicalizes to a `dup` at the most-3' rotation-aligned position. Closes [#132](https://github.com/fulcrumgenomics/ferro-hgvs/issues/132).
- *(normalize)* recognize a reverse-complement sub-span within a `delins` (synthesized by cis-allele merge OR user-typed) and emit the spec-canonical `inv`, splitting the surrounding span into `[…; inv; …]`. The HGVS edit-priority rule places `inv` above `delins` in the priority order (`general.md:56`) and defines `delins` as the residual when no higher-priority form applies (`delins.md`). For example, `g.[1150T>G;1151C>A;1152C>G]` (over `TCC`) now normalizes to `g.[1150_1151inv;1152C>G]` instead of `g.1150_1152delinsGAG`; `g.[1092G>C;1093G>C]` (over `GG`) normalizes to `g.1092_1093inv` instead of `g.1092_1093delinsCC`. The same rule fires for a user-typed `g.1150_1152delinsGAG`, since the canonical form depends on `(ref, position, alt)` and not on input shape. Applies across all NA coord systems: `g.`, `m.`, `c.` (CDS-proper positions), `n.`, `r.` (T/U-equivalent comparison so `r.` alts with `U` align with transcript ref bytes that contain `T`). Sub-only decomposition (rewriting a non-inv multi-sub `delinsXY` to `[X>...; Y>...]`) is intentionally left out of scope and is a separate spec interpretation question. The codon-frame `c.` merge from issue [#79](https://github.com/fulcrumgenomics/ferro-hgvs/issues/79) is preserved automatically — the synthesized middle base in a codon-frame-merged delins makes a length-2 inv across the middle mathematically impossible. Closes issue [#160](https://github.com/fulcrumgenomics/ferro-hgvs/issues/160).
- *(normalize)* rewrite degenerate substitutions (ref == alt, e.g. `c.100A>A`) to identity (`=`) per HGVS v21 spec, which marks `c.X>X` as "not allowed" (`docs/recommendations/DNA/other.md`). The rule is purely syntactic on the edit's stated bases, so it fires in both the full-normalization path and the no-reference canonicalization path — `c.123C>C` rewrites to `c.123=` regardless of provider availability. ([#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) A4)
- *(normalize)* rewrite a `delins` whose inserted sequence is empty as a deletion, per the HGVS spec requirement that an empty insert is semantically a deletion and must be rendered as `del`. The rewritten deletion is then 3'-shifted under the standard del rule. Issue [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A3.
- *(normalize)* rewrite a `delins` whose inserted sequence is the reverse complement of the deleted reference as an inversion, per the HGVS spec definition of `inv`. The complementary-outer-bases shortening rule applies to the result so that a `delins`-encoded inversion produces the same canonical output as a directly-encoded `inv`. Issue [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A2.
- *(normalize)* canonicalize `delins` to the minimal HGVS form by trimming any shared prefix/suffix between the inserted sequence and the deleted reference, then reclassifying the residual edit as substitution / deletion / insertion / inversion / smaller `delins` per the sub > del > inv > dup > ins priority. For example, `c.1_4delinsAAGC` against ref `ATGC` collapses to `c.2T>A`; `c.1_4delinsAC` against ref `ATGC` collapses to `c.2_3del`; `c.1_3delinsATCG` against ref `ATG` collapses to `c.2_3insC`. Extension to issue [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A2.
- *(normalize)* apply inversion shortening to direct `inv` inputs. `NaEdit::Inversion` was missing from `needs_normalization()`, so direct inversions bypassed `normalize_na_edit` entirely and the `shorten_inversion()` / identity-collapse logic was never exercised. Direct `inv` variants now also emit minimal notation after shortening (no stale explicit `sequence`/`length` from the input).
- *(normalize)* HGVS spec compliance ([#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) B1): repeat-notation rewrites in `c.` (coding DNA) context now enforce the spec's codon-frame exception — repeat notation requires `unit_len % 3 == 0`. Previously `insertion_to_repeat`, `deletion_to_repeat`, `duplication_to_repeat`, and `normalize_repeat` could emit `c.X_YA[N]`, `c.X_YAT[N]`, etc. for non-codon-aligned units, violating the spec (`docs/recommendations/DNA/repeated.md`). Variants in `c.` with non-codon-aligned units now retain the spec-prescribed alternative form: `dup` for 1 added copy, `ins<literal>` for ≥2 added copies, and plain `del` for contractions of ≥2 unit copies.

### Changed (public API)

- `insertion_to_repeat`, `duplication_to_repeat`, and `normalize_repeat` in `ferro_hgvs::normalize::rules` gain an `is_coding: bool` parameter to drive the codon-frame gate. (`deletion_to_repeat` is `pub(crate)` and gains the same parameter as an internal change.)
- `RepeatNormResult::Insertion { start, end, sequence }` and `DupToRepeatResult::GatedInsertion { start, end, sequence }` variants added so the rule layer can route gated rewrites to the spec-canonical literal `ins` form.
- *(normalize)* `NormalizationWarning` is now a sum type (`RefSeqMismatch` / `OverlapConflict`). Each warning code carries only the fields relevant to it. Read sites migrate from `.code` field access to `.code()` method and pattern-matching on the variant.

### Changed

- Internal: the four `delins_is_*` boolean helpers in `normalize::rules` are unified into one `canonicalize_delins()` function returning a `DelinsCanonical` enum, expressing HGVS edit-priority (sub > del > inv > dup > ins) in a single decision tree. The unreachable second delins arm in `normalize_na_edit` has been removed.

## [0.4.1](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.4.0...v0.4.1) - 2026-05-04

### Added

- *(python)* load reference data via Normalizer.from_manifest and extended from_json ([#86](https://github.com/fulcrumgenomics/ferro-hgvs/pull/86))
- *(normalize)* merge consecutive edits in cis alleles per HGVS spec ([#80](https://github.com/fulcrumgenomics/ferro-hgvs/pull/80))

### Fixed

- *(normalize)* codon-frame exception for c. SNV pairs separated by one nucleotide ([#79](https://github.com/fulcrumgenomics/ferro-hgvs/pull/79)) ([#104](https://github.com/fulcrumgenomics/ferro-hgvs/pull/104))
- *(normalize)* 5'UTR CDS↔tx off-by-one collapsed UTR del to c.? ([#97](https://github.com/fulcrumgenomics/ferro-hgvs/pull/97)) ([#102](https://github.com/fulcrumgenomics/ferro-hgvs/pull/102))
- *(normalize)* merge same-region UTR adjacency in cis alleles ([#89](https://github.com/fulcrumgenomics/ferro-hgvs/pull/89)) ([#103](https://github.com/fulcrumgenomics/ferro-hgvs/pull/103))
- *(normalize)* minus-strand intronic ref-base orientation ([#98](https://github.com/fulcrumgenomics/ferro-hgvs/pull/98)) ([#100](https://github.com/fulcrumgenomics/ferro-hgvs/pull/100))
- rewrite delins as identity when insert matches reference ([#78](https://github.com/fulcrumgenomics/ferro-hgvs/pull/78))
- rewrite single-base delins as substitution per HGVS priority ([#77](https://github.com/fulcrumgenomics/ferro-hgvs/pull/77))
- Emit single-variant alleles in bare spec form ([#76](https://github.com/fulcrumgenomics/ferro-hgvs/pull/76))

### Other

- unpin CI Rust toolchain and fix 1.95 clippy lints ([#108](https://github.com/fulcrumgenomics/ferro-hgvs/pull/108))
- dup 3'-shift coverage matrix (#81 A6) ([#107](https://github.com/fulcrumgenomics/ferro-hgvs/pull/107))
- del 3'-shift coverage matrix + tandem-repeat del canonical-form fix (#81 A5, B2) ([#106](https://github.com/fulcrumgenomics/ferro-hgvs/pull/106))
- tighten coverage_gap_tests assertions; restore intronic-ins coverage ([#94](https://github.com/fulcrumgenomics/ferro-hgvs/pull/94)) ([#99](https://github.com/fulcrumgenomics/ferro-hgvs/pull/99))
- *(coverage)* restore intronic insertion tests dropped in #93 ([#95](https://github.com/fulcrumgenomics/ferro-hgvs/pull/95))
- ins 3'-shift coverage matrix + tandem-repeat ins canonical-form fix (#81 A1, A7) ([#93](https://github.com/fulcrumgenomics/ferro-hgvs/pull/93))

### Added

- *(test)* HGVS v21.0 spec normalization regression fixture pinning ferro's current `normalize()` output for every variant string in the spec (823 rows). Each row carries a status (`preserved` / `diverges` / `correctly-rejected` / `false-acceptance` / `parse-error` / `needs-reference`) derived from `(parse_ok, normalize_ok, spec_expected, current == spec_expected)`. `spec_expected: null` is the sentinel for "spec rejects this," set automatically for inputs the spec marks via `<code class="invalid">…</code>` (35 occurrences in v21.0) and overridable by hand. Bare illustrative fragments like `c.1083A>C` get a default accession prepended per coord system so they parse — the prefixed form is recorded in a separate `input_prefixed` field. The fixture is generated by `examples/generate_spec_fixture.rs` from the vendored spec at `assets/hgvs-nomenclature/` (git submodule pinned to tag `21.0.0`). Hand overrides live in `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json` and accept `status`, `spec_expected` (string / null / absent), `input_prefixed`, `requires_reference`, and `todo`. CI runs the generator in `--check` mode to enforce byte-identical regeneration. Closes [#84](https://github.com/fulcrumgenomics/ferro-hgvs/issues/84); companion to [#83](https://github.com/fulcrumgenomics/ferro-hgvs/issues/83).
- 3'-shift coverage matrix for insertion variants
  (`tests/ins_shift_matrix.rs`): 84 rstest cases across all 7 nucleotide
  coord-system / strand combinations × 8 shuffle scenarios, with a
  shared `SyntheticBuilder` fixture helper at `tests/common/synthetic.rs`
  (reusable by future del / dup / repeat-notation matrices). Issue #81
  item A7.
- Tightened `tests/coverage_gap_tests.rs` from canary-shaped
  (`contains(...)`, `is_ok() || is_err()`) assertions to strict
  `assert_eq!` against the exact normalizer output. Builds on PR #95
  (intronic-insertion restoration) by extending the same tightening
  pass to the remaining gap-test buckets, and adds cross-references
  and procedural detail to the restored intronic-insertion comments.
  Several locked outputs are suspected-buggy and carry `FIXME(#NN)`
  comments pointing to follow-up tracking issues — most converge on
  #98 (minus-strand intronic reference-base orientation, identified
  as a likely common root cause behind A1 / A5 / A6 misfires on
  minus-strand intronic positions and the PR #78 identity-rewrite
  missing on minus strand), with #96 (wrong-strand repeat unit
  emission on minus-strand multi-base dup) and #97 (`c.?del` collapse
  on minus-strand 5'UTR del) as separate sub-symptoms. True
  boundary-spanning panic-canary tests are intentionally left as
  `is_ok() || is_err()` per the issue's scope. Issue #94.
- *(normalize)* Merge consecutive sub-variants in cis alleles into a single delins per HGVS spec. `g.[1000G>A;1001A>C]` now normalizes to `g.1000_1001delinsAC`; `g.[1000del;1001del]` to `g.1000_1001del`. Covers `g./c./n./r./m.` coordinate systems, sub/del/delins/ins edit combinations, chains, and same-boundary insertion pairs. Non-adjacent variants, intronic/UTR boundaries, uncertain edits, and non-`Literal` insertion payloads are barriers and pass through unchanged. The codon-frame exception (one-nt gap within a codon) is tracked separately in [#79](https://github.com/fulcrumgenomics/ferro-hgvs/issues/79). ([#80](https://github.com/fulcrumgenomics/ferro-hgvs/pull/80))

### Fixed

- *(normalize)* Cis-allele consecutive-edit merging now collapses
  adjacent sub-variants *within* a UTR / upstream / downstream region,
  not only within the CDS / transcript body. The PR #80 implementation
  rejected every `is_5utr() / is_3utr() / is_upstream() / is_downstream()`
  position outright, so inputs like `c.[-2A>G;-1C>T]` (both 5'UTR) and
  `c.[*1A>G;*2C>T]` (both 3'UTR) round-tripped unchanged even though
  HGVS allows ranges within those regions (`c.-2_-1`, `c.*1_*2`). The
  fix replaces the `Option<u64>` position keys with a `Region`-tagged
  `(Region, i64)` axis (covering CDS, 5'UTR, 3'UTR, transcript-body,
  upstream, downstream, and genomic / mitochondrial); merge eligibility
  becomes "same region + integer adjacency on the region's axis".
  `build_cds_merged` / `build_tx_merged` / `build_rna_merged` consume
  the region tag to reconstruct the right `CdsPos` / `TxPos` / `RnaPos`
  shape (negative base for 5'UTR / upstream, `utr3` / `downstream` flag
  for 3'UTR / downstream). Cross-region pairs (`c.[-1A>G;1A>T]` 5'UTR↔CDS,
  `c.[40C>T;*1A>G]` CDS↔3'UTR, …) still correctly do not merge.
  Issue #89.
- *(normalize)* CDS ↔ transcript coordinate mappings now respect the
  HGVS no-c.0 numbering rule for 5'UTR positions. The forward mapping
  (`Normalizer::cds_to_tx_pos`, `convert::coding::cds_to_transcript_pos`)
  previously computed `tx = cds_start + base - 1` for negative `base`,
  which double-counted the gap between c.-1 and c.1 and emitted a tx
  position one base 5' of the true location. The inverse mapping
  (`Normalizer::tx_to_cds_pos`) had the mirror bug: tx positions one
  base before `cds_start` mapped to `base = 0`, which `CdsPos::Display`
  renders as `c.?` (`CDS_BASE_UNKNOWN`). The most visible symptom was a
  5'UTR single-base deletion on a minus-strand transcript collapsing
  to `c.?del` instead of resolving to a real position (e.g. `c.-1del`
  after 3'-shifting within a UTR homopolymer). Forward and inverse are
  now `tx = cds_start + base` and `base = tx - cds_start` for negative
  / pre-cds_start positions, matching the spec's "c.-1 is one base 5'
  of c.1" rule and the existing exon-aware mapper at
  `convert::mapper::cds_to_tx` / `tx_to_cds`. Issue #97.
- *(normalize)* Cis-allele consecutive-edit merging now also collapses
  the codon-frame exception case: two `c.` exonic SNVs in the CDS
  proper, separated by exactly one nucleotide, that fall within the
  same codon merge into a single delins with the unchanged middle
  reference base preserved verbatim — per HGVS spec
  (`c.[145C>T;147C>G]` → `c.145_147delinsTGG`, where the middle base
  is the reference at `c.146`). Eligibility is narrow: same accession,
  both endpoints in `Region::Cds`, gap-of-one on the axis, both prev
  and next are single-base SUB anchors, and `(prev-1)/3 == (next-1)/3`
  (same codon). The unchanged middle base is fetched via the
  `ReferenceProvider` threaded into `merge_consecutive_edits`; if no
  transcript is registered or the position is out of range, the merge
  is silently declined and the variants pass through unchanged.
  Codon-frame–merged delins continue to participate in the
  strictly-consecutive walk, so a third SNV one base 3' of the pair
  (`c.[10A>G;12A>C;13A>T]`) still folds into the running delins.
  Cross-codon (`c.[3G>T;5A>C]`), gap-of-two (`c.[10A>G;13A>C]`), and
  non-CDS pairs (`g.`, UTR, `n.`, `r.`) all correctly do not merge.
  Issue #79.
- *(normalize)* Minus-strand intronic normalization now reads the
  reference window in transcript-view orientation. `normalize_intronic_cds`
  and `normalize_intronic_tx` previously passed the genomic-strand bytes
  fetched from `get_genomic_sequence` directly into `normalize_na_edit`
  alongside the variant's transcript-view edit alt; on minus-strand
  transcripts the two were mis-oriented, defeating every rule that
  compared the alt against the local reference window. The fix
  reverse-complements the genomic window and flips the relative
  positions / shuffle boundaries on minus strand before normalization,
  then maps the resulting positions back to the genomic frame for the
  CDS / tx coordinate conversion. As a single root-cause fix this
  resurfaces #81 A1 / A5 / A6 canonicalization, the PR #78 delins-as-
  identity rewrite, and the transcript-view repeat-unit letter on
  minus-strand intronic positions — all of which had been observed to
  misfire in #94's locking pass. Issue #98.
- Insertions that add ≥2 copies of a multi-base tandem repeat unit now
  emit repeat notation (`unit[N+k]`) instead of a duplication of the
  inserted sequence, per HGVS spec ("when more than one additional copies
  are inserted directly 3' of the original copy, the change is indicated
  using the format for Repeated sequences"). Single-unit additions remain
  `dup`. Issue #81 item A7.
- `MockProvider::get_sequence` now falls through to genomic contig
  lookup when the id is not a transcript, matching `FastaProvider`'s
  behavior. This unblocks 3'-shift normalization for genomic test
  fixtures that register only `add_genomic_sequence`. Issue #81 item A7.

## [0.4.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.3.0...v0.4.0) - 2026-04-30

### Added

- *(python)* add poethepoet task runner for Python dev workflow ([#53](https://github.com/fulcrumgenomics/ferro-hgvs/pull/53))
- Add variant accessor properties to Python bindings ([#49](https://github.com/fulcrumgenomics/ferro-hgvs/pull/49))

### Fixed

- accept gene selectors on non-RefSeq accessions ([#70](https://github.com/fulcrumgenomics/ferro-hgvs/pull/70))
- Use HGVS spec compact form for allele Display ([#48](https://github.com/fulcrumgenomics/ferro-hgvs/pull/48))

### Other

- *(readme)* vendor Fulcrum logo and use absolute URLs ([#71](https://github.com/fulcrumgenomics/ferro-hgvs/pull/71))
- bump pyo3 0.23 → 0.28 and add Python 3.14 wheels ([#57](https://github.com/fulcrumgenomics/ferro-hgvs/pull/57))
- *(python)* use dependency-groups, stricter mypy, and --locked in CI ([#52](https://github.com/fulcrumgenomics/ferro-hgvs/pull/52))
- publish Python wheels to PyPI via Trusted Publishing ([#58](https://github.com/fulcrumgenomics/ferro-hgvs/pull/58))
- *(python)* drop Python 3.8 and 3.9 support ([#55](https://github.com/fulcrumgenomics/ferro-hgvs/pull/55))
- *(python)* add uv lockfile for reproducible dev environment ([#50](https://github.com/fulcrumgenomics/ferro-hgvs/pull/50))
- build Python wheels and attach to GitHub Releases ([#39](https://github.com/fulcrumgenomics/ferro-hgvs/pull/39))
- *(prepare)* Modularize ReferenceManifest ([#44](https://github.com/fulcrumgenomics/ferro-hgvs/pull/44))
- switch reqwest from native-tls to rustls-tls ([#38](https://github.com/fulcrumgenomics/ferro-hgvs/pull/38))

### Changed

- Allele `Display` now emits HGVS spec-correct compact form (`ACC:c.[edit1;edit2]`) when sub-variants share an accession and coordinate type, instead of the expanded form (`[ACC:c.edit1;ACC:c.edit2]`). Mixed-accession alleles and alleles containing the per-variant unknown form (`c.?`, `r.?`, etc.) still emit the expanded form. Downstream consumers parsing the previous expanded output (including Python `str(variant)`) will see the new format. ([#48](https://github.com/fulcrumgenomics/ferro-hgvs/pull/48))

### Fixed

- Prevent panic when `NullAllele`/`UnknownAllele` are used as sub-variants in an allele ([#48](https://github.com/fulcrumgenomics/ferro-hgvs/pull/48))

## [0.3.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.2.0...v0.3.0) - 2026-03-30

### Added

- Bincode serialization for fast cdot loading ([#26](https://github.com/fulcrumgenomics/ferro-hgvs/pull/26))
- Compound reference syntax `NC_*(NM_*):c.…` support ([#21](https://github.com/fulcrumgenomics/ferro-hgvs/pull/21))
- GRCh37 cdot transcript metadata download in `ferro prepare` ([#24](https://github.com/fulcrumgenomics/ferro-hgvs/pull/24))

### Fixed

- Genomic-to-coding coordinate conversion ([#25](https://github.com/fulcrumgenomics/ferro-hgvs/pull/25))
- Coding-order positions after genomic-space normalization ([#20](https://github.com/fulcrumgenomics/ferro-hgvs/pull/20))
- Mutalyzer configuration across web service and benchmarks ([#19](https://github.com/fulcrumgenomics/ferro-hgvs/pull/19))
- Serde defaults for config structs and deploy workflow ([#27](https://github.com/fulcrumgenomics/ferro-hgvs/pull/27))
- Missing `in_memory` field in HgvsRsConfig construction ([#28](https://github.com/fulcrumgenomics/ferro-hgvs/pull/28))

### Changed

- Replaced vendored hgvs-rs with published crate v0.20.1 ([#22](https://github.com/fulcrumgenomics/ferro-hgvs/pull/22))

## [0.2.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.1.0...v0.2.0) - 2026-03-23

### Fixed

- derive Python package version from Cargo.toml via maturin ([#14](https://github.com/fulcrumgenomics/ferro-hgvs/pull/14))
- allow valid HGVS characters in web service input validation ([#13](https://github.com/fulcrumgenomics/ferro-hgvs/pull/13))
- work-stealing for mutalyzer normalize shards ([#4](https://github.com/fulcrumgenomics/ferro-hgvs/pull/4))
- CIGAR-aware CDS-to-transcript coordinate mapping ([#7](https://github.com/fulcrumgenomics/ferro-hgvs/pull/7))
- prevent integer overflow and handle edge cases in normalization ([#6](https://github.com/fulcrumgenomics/ferro-hgvs/pull/6))
- prevent panics on malformed HGVS patterns during normalization ([#5](https://github.com/fulcrumgenomics/ferro-hgvs/pull/5))

### Other

- Add Bioconda and Zenodo badges ([#3](https://github.com/fulcrumgenomics/ferro-hgvs/pull/3))
- refresh test fixtures and data extraction scripts ([#8](https://github.com/fulcrumgenomics/ferro-hgvs/pull/8))
- release v0.1.0

## [0.1.0](https://github.com/fulcrumgenomics/ferro-hgvs/releases/tag/v0.1.0) - 2026-02-19

### Other

- Initial commit
