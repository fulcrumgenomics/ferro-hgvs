//! The four properties, over the exhaustive spec-derived conformance corpus.
//!
//! # What this measures, and why it can measure it without an oracle
//!
//! Conformance used to be measured by the generated spec fixture: 934 harvested
//! rows, **52** of them multi-member cis alleles and **40** carrying a
//! spec-stated answer. That ceiling is the spec's, not ours — "matches the form
//! the spec publishes" needs the spec to have published one — so a change to the
//! multi-member partitioner measured `0 gained, 0 lost`. That is measurement
//! blindness, not neutrality.
//!
//! `ferro_hgvs::conformance::spec_corpus` gets past it by generating *many
//! spellings of one variant*: re-partition it, and shift each partition within
//! its ambiguous run. Every spelling in such a family denotes the same sequence
//! **by construction**, which makes four properties checkable with no published
//! answer at all:
//!
//! | rank | property | how it is checked |
//! |---|---|---|
//! | 1 | **validity** | the output re-parses, denotes a sequence, and violates no absolute prohibition |
//! | 2 | **confluence** | every spelling of one variant reaches ONE output |
//! | — | **idempotency** | each output is its own fixed point |
//! | — | **sequence preservation** | each output denotes the input's bases, via `hgvs_to_spdi` |
//!
//! Rank 1 and rank 2 are the operator's precedence order
//! (`rulings[adjudication-precedence-order]`); the other two are invariants
//! rather than ranked choices.
//!
//! # This is a census, not a gate
//!
//! None of these properties holds today and a permanently red suite is worth
//! nothing, so every figure below is a **pinned baseline**: each failure figure
//! may only ever go DOWN and `converged` only ever UP. A change that moves one
//! must re-bless the pin in the same commit and say so in the PR.
//!
//! # The pins are NOT an acceptable state — five counters are rank-1 defects
//!
//! Pinning a non-zero figure records a defect; it does not excuse one. These five
//! are conformance defects that this instrument found on its first run, and none
//! is fixed here (a corpus PR that also changed the normalizer would make its own
//! baseline unmeasurable).
//!
//! **The base they are measured against is `main`.** Every figure below was
//! measured on the branch that introduced this module, which added no
//! `src/normalize/` change of its own — so the census records `main`'s
//! behaviour, not that PR's. The corpus is new; the defects it counts are not.
//!
//! **That paragraph is the CORPUS branch speaking, and it no longer holds.**
//! Five changes have re-blessed figures here since, and there are five
//! RE-BLESSED sections below — one per change, in the order they landed. Read
//! all five before quoting any figure as a `main` baseline:
//!
//! 1. the amino-acid precondition on `coalesce_coding_frame_separation` (#1599),
//!    which moved five figures, one of them a net regression in `converged`;
//! 2. the span-preserving re-typing carve-out for a member that straddles a CDS
//!    boundary (#1536), which moved seven and regressed none;
//! 3. the two-deletion alignment in the payload splitter (#1649), which moved
//!    eight and regressed none — the largest confluence move recorded here;
//! 4. the `standards.md:39` refusal (#1627), which moved two REFUSAL/validity
//!    figures, both down, identically in both directions;
//! 5. the codon-frame merge's own span (#1716), which moved six rank-2 figures,
//!    left `converged` alone in both directions, and raised `split_three` at 3'.
//!
//! **The second was measured on the first, not composed with it.** The two
//! changes' affected row sets are disjoint — verified by diffing the row ids, not
//! assumed — so the deltas happen to be additive here, but a census counts
//! *classes* rather than row-deltas and a composed figure would not have been a
//! measurement. Every figure in the second section came from running the corpus
//! on this branch rebased onto #1599.
//!
//! | counter | 3' | 5' | what it is |
//! |---|---|---|---|
//! | `outputs_leaving_the_transcript` | ~~371~~ **0** | 0 | a `c.`/`n.` output naming an INTRONIC position its input did not, **on a bare accession**. `checklist.md:20` — a bare `NM_` cannot express an intronic position at all. **Closed by #1704**; see RE-BLESSED (3 of 3) |
//! | `outputs_intronic_under_a_genomic_wrapper` | **371** | 0 | the same 371 rows, now rendered as `NC_…(NM_…):c.…` — CONFORMANT, and not a defect counter |
//! | `outputs_denoting_no_sequence` | 10 | 18 | two members of the OUTPUT claim one territory, so it denotes nothing — rank-1 invalid |
//! | `sequence_changed` | 4 | 0 | the output denotes different bases from the input; a member was dropped |
//! | `non_idempotent_outputs` | ~~7~~ **4** | 4 | the output is not its own fixed point. Re-blessed: see the note below |
//! | `conflicts_accepted` | 72 | 72 | a conflicting allele was normalized instead of refused — nested, partially overlapping, and two insertions at one interbase |
//! | `prohibited_absolute_accepted` | 32 | 32 | a shape the spec calls "not allowed" was accepted |
//! | `prohibition_violating_outputs` | ~~32~~ **8** | ~~32~~ **8** | and then EMITTED, so the prohibition is not enforced on output either. Re-blessed by #1627: see the fourth section below |
//!
//! `guard_violations` is **0 of 210 guarded rows**, and
//! [`RESIDUAL_GUARD_VIOLATIONS`] is empty. The denominator is asserted non-zero,
//! because `0 of 0` is what a rebuilt #1456 looks like — and the instrument's
//! ability to read non-zero is not inferred either: it read **18** on this
//! branch's first commit and **2** on its second, so the zero below is a
//! measurement and not a structural silence.
//!
//! **It was 0, then 18, then 2, and is now 0 again**, and the middle figures are
//! the ones to understand. Deleting the input-relative weight bound
//! (`rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`) let the
//! partitioner re-derive multi-member spellings instead of handing them back —
//! and the partitioner was applying `DNA/delins.md:44-47`'s payload-coincidence
//! carve-out on **every** axis, so 18 authored frameless pairs at separation one
//! came back merged. Scoping that carve-out to `c.`, which
//! `rulings[delins-payload-coincidence-carve-out-is-coding-dna-scoped]` had
//! already decided and which the non-shipping `CanonicalCoalesced` arm already
//! honoured, closed 16 of them.
//!
//! **The last 2 were `MAX_SPLIT_BLOCK`, and the fix is the same clause one hop
//! earlier.** Both are 2 065-base blocks, over the 1 024 cap, so
//! `partition_block` returns the whole block on length *before any rule about
//! the derived pieces runs*. What comes back is therefore not a finding that the
//! block holds one variant — it is the absence of a finding, and emitting it
//! asserts "no separation here" about a block nothing examined. Off `c.`
//! `general.md:34` governs that assertion and nothing licenses it, so
//! `canonicalize_from_sequence` now **declines** such a block and the per-member
//! pipeline answers. On `c.` the carve-out is in reach and `:47` recommends the
//! spanning `delins`, so nothing there moves. See
//! `merge::CoincidenceCarveOut::block_may_be_asserted_unexamined` and
//! `tests/it/issue_1616_unexamined_block_stays_split.rs`.
//!
//! **What that costs, disclosed rather than netted away.** A block past the cap
//! cannot be split, so its spanning spelling and its split spelling can no
//! longer reach one string: `converged` falls **18** at 3' (10 995 -> 10 977) and
//! **14** at 5' (10 789 -> 10 775), with `split_two` rising by exactly the same
//! amounts. Confluence is ferro's own policy (README rule 3) and `general.md:34`
//! is a spec preference clause (rule 2), which outranks it — so this is the
//! ruled trade and not an accident. Against `origin/main` the direction is still
//! strongly positive: 9 402 -> 10 977 (+1 575) and 9 228 -> 10 775 (+1 547).
//! **Raising the cap would buy those classes back and is still open**; it is a
//! cost question with no clause behind it, and it is deliberately not answered
//! here.
//!
//! **These are a rule-2 preference miss, not a rule-1 violation.**
//! `rulings[separation-rule-force-modal-or-negation]` holds that
//! `general.md:34`'s modal grades the whole clause and its "and not" names the
//! excluded alternative rather than prohibiting it — the two halves are
//! complements, so a prohibition would make the individual form mandatory and
//! leave "should" doing no work. So this counter is a best-effort preference
//! figure. It is still pinned at its true value and still tripwired, because a
//! preference clause outranks maintainer judgment.
//!
//! The 3'/5' asymmetry is itself a finding: the junction-crossing class is
//! **371 to 0**, and 3' is the default direction. #1704 makes those 371 outputs
//! conformant without moving a single coordinate, so the asymmetry survives the
//! fix intact — it is now carried by `outputs_intronic_under_a_genomic_wrapper`,
//! and it is still a claim about the code rather than about the corpus.
//!
//! # RE-BLESSED — #1616, and the row that retired with it
//!
//! `spec_corpus_regressions::the_codon_gate_splits_a_spanning_delins_its_own_members_do_not`
//! is **deleted** by this change, and the reason belongs here because it is the
//! same corpus. It held the `m4-all-del-p1-sep2` design's two spellings apart:
//!
//! | spelling | before | now |
//! |---|---|---|
//! | spanning `c.24_33delinsAATTTA` | `c.[24T>A;26_29del;32_33delinsTA]` | **unchanged** |
//! | members `c.[24del;27del;30del;33del]` | `c.[24del;30_33delinsA]` | **`c.[24T>A;26_29del;32_33delinsTA]`** |
//!
//! Measured both strands, both outputs idempotent. Only the member side moved,
//! onto the string the spanning side already printed — a convergence, so that
//! design contributes to `converged` rather than to `split_two` in the censuses
//! below, and the row was the last survivor of a six-row family (#1649 disposed
//! of four the same way). The clause is
//! `rulings[canonical-form-choice-when-both-legal]`.
//!
//! **Not a structural zero.** `spec_corpus::DENSE_SEPARATIONS` is still
//! `[0, 1, 2, 3, 5, 8]`, so the generator still builds the `sep2` design the
//! deleted row's doc warned could quietly stop existing (#1456/#1460/#1478).
//!
//! # RE-BLESSED (1 of 4) — #1599, the amino-acid precondition
//! # RE-BLESSED (1 of 5) — #1599, the amino-acid precondition
//!
//! Everything above and in [`THREE_PRIME`]'s own doc was written on the corpus
//! PR, which touched no normalizer file. **That is no longer the base.** #1599
//! adds `general.md:35`'s amino-acid conjunct to
//! `coalesce_coding_frame_separation`, so the pass now declines a merge whose
//! span crosses a codon boundary, and five figures moved:
//!
//! | figure | was | now | direction |
//! |---|---|---|---|
//! | 3' `non_idempotent_outputs` | 7 | **4** | improves — the three `scale-c3p-sep{120,128,136}` rows became fixed points |
//! | 3' `converged` | 9,140 | **9,139** | **regresses by one, net** — 6 classes lost convergence and 5 gained it |
//! | 3' `split_two` | 2,435 | **2,436** | the same net one, on the other side of the ledger |
//! | 5' `split_two` | 2,696 | **2,698** | improves — two 3-output classes became 2-output |
//! | 5' `split_three` | 204 | **202** | the same two classes |
//!
//! **`converged` going DOWN is the disclosure this section exists for**, and the
//! net figure hides the shape of it. Measured by dumping both directions' full
//! divergence lists on `origin/main` and on this branch and diffing the row ids:
//!
//! - **6 rows lost convergence**, three designs on each of the plus and minus
//!   coding multi-exon shapes: `s00-c3{p,m}-m3-all-del-p1-sep3`,
//!   `s00-c3{p,m}-m4-all-del-p1-sep2` and `s00-c3{p,m}-pair-del-del-p1-sep8`.
//!   In every one it is the **spanning-`delins` respelling** that moved; the
//!   authored multi-member spellings are unchanged. They are promoted to
//!   `spec_corpus_regressions::the_codon_gate_splits_a_spanning_delins_its_own_members_do_not`.
//! - **5 rows gained convergence**: `s01-c1-pair-del-del-p2-sep5`,
//!   `s01-c3{p,m}-m3-all-del-p1-sep2` (whose spanning `delins` used to be a
//!   fixed point disagreeing with its own members) and
//!   `scale-c3p-sep{120,128}-del-del` (which used to need a second pass to reach
//!   the answer their members reach).
//!
//! Both counters therefore move for one reason — the pass no longer merges
//! across a codon boundary — and neither direction of the move is a
//! representation *choice* left open: `general.md:35`'s second conjunct is
//! either met or it is not. `outputs_leaving_the_transcript` stayed **371** across
//! #1599, which was the pre-existing baseline and not a regression of its — it is
//! #1704, three sections down, that finally moves it.
//!
//! # RE-BLESSED (2 of 5) — #1536, the cross-axis re-typing carve-out
//!
//! This branch adds a third carve-out from the #350 cross-axis bail in
//! `normalize_cds`, so a `delins`/`inv` whose span straddles a CDS boundary is
//! re-typed against the bases it denotes instead of being returned unprocessed.
//! Seven figures moved and **none regressed**:
//!
//! | figure | was (post-#1599) | now | direction |
//! |---|---|---|---|
//! | 3' `converged` | 9,139 | **9,141** | improves by 2 |
//! | 3' `split_two` | 2,436 | **2,440** | net +4, from higher arities |
//! | 3' `split_three` | 250 | **248** | improves |
//! | 3' `split_more` | 57 | **53** | improves |
//! | 5' `split_two` | 2,698 | **2,706** | net +8, from higher arities |
//! | 5' `split_three` | 202 | **200** | improves |
//! | 5' `split_more` | 38 | **32** | improves |
//!
//! **The `was` column is measured, not carried over.** It was re-derived by
//! running this branch's corpus against `main`'s `src/normalize/mod.rs` — which
//! reproduced the committed post-#1599 pins exactly, `rc=0` — so the two sides of
//! the table come from one corpus and one machine.
//!
//! **`split_two` rising is not a regression here**, and the net counters cannot
//! settle that on their own: a family entering `split_two` from `split_three` and
//! a family leaving `converged` for `split_two` look identical in the totals. So
//! the full divergence lists were dumped on both sides and the row ids diffed:
//!
//! - **2 families newly converge** (3' only): `s00-c3{p,m}-cds-end-del-ins-p1-sep2`,
//!   each from arity 2. At 5' the same two rows only drop 3 -> 2.
//! - **10 families drop arity at 3'** and **14 at 5'**, 10 of them the same row
//!   ids in both directions: `s00-c3{p,m}-cds-end-` `del-del-p1-sep2` (3 -> 2),
//!   `del-del-p2-sep1` (4 -> 3), `del-ins-p2-sep1` (4 -> 3), `sub-sub-p1-sep0`
//!   (3 -> 2) and `sub-sub-p2-sep0` (3 -> 2). The 5'-only four are
//!   `s00-c3{p,m}-cds-end-del-ins-p1-sep2` (3 -> 2) and
//!   `s00-c3{p,m}-cds-end-del-ins-p2-sep2` (4 -> 3).
//! - **0 families lose convergence, and 0 rise in arity, in either direction.**
//!
//! So there is **nothing to promote** to `spec_corpus_regressions.rs` for this
//! change — which is a measured claim about a diffed row set, not the absence of
//! a search. Contrast #1599 directly above, which did have six rows to promote.
//!
//! **Every moved row is `s00-c3{p,m}-cds-end-*`.** That is the corpus speaking
//! about itself, not a statement that the fix only reaches the CDS-end boundary:
//! this corpus builds boundary designs at `cds_end` and none at `cds_start`, so
//! the 5'UTR/CDS half of the carve-out is **unmeasured here**. It is covered
//! instead by `issue_1536_cds_boundary_delins`, whose `PLACEMENTS` table crosses
//! `cds_start` (`c.-3_5`) and whose 33-placement walk crosses both.
//!
//! `non_idempotent_outputs` stays **4** at 3' and **4** at 5' — measured, and the
//! recursion step in `normalize_cds` is what keeps it there; without that step it
//! measures **6** at 3', the two extra rows being
//! `s00-c3{p,m}-cds-end-del-del-p2-sep1`. Every rank-1 counter is unchanged.
//!
//! # RE-BLESSED (3 of 5) — #1649, the two-deletion alignment
//!
//! `merge`'s payload splitter could express *insertion, retained reference,
//! insertion* but not *deletion, retained reference, deletion*, so a payload
//! whose alignment against its span carries two gaps was forced onto a coarser
//! partition than the sequence supports. Letting it express the second shape
//! moves eight figures, all rank 2, and **none regresses**:
//!
//! | figure | was (post-#1536) | now | direction |
//! |---|---|---|---|
//! | 3' `converged` | 9,141 | **9,402** | improves by 261 |
//! | 3' `split_two` | 2,440 | **2,225** | improves by 215 |
//! | 3' `split_three` | 248 | **223** | improves by 25 |
//! | 3' `split_more` | 53 | **32** | improves by 21 |
//! | 5' `converged` | 8,944 | **9,228** | improves by 284 |
//! | 5' `split_two` | 2,706 | **2,461** | improves by 245 |
//! | 5' `split_three` | 200 | **167** | improves by 33 |
//! | 5' `split_more` | 32 | **26** | improves by 6 |
//!
//! This is the largest confluence move recorded in this module, which is exactly
//! why the net counters are not allowed to carry it. **A family entering
//! `split_two` from `split_three` and a family leaving `converged` for
//! `split_two` look identical in the totals** — so the full divergence row-id
//! lists were dumped on `origin/main` and on this branch (the `divergences` cap
//! in [`measure`] raised locally for the measurement, then restored) and the sets
//! diffed:
//!
//! - **0 families lose convergence**, at 3' or at 5'.
//! - **0 families rise in arity**, and **0** families become divergent that were
//!   not divergent on `main` — the branch's divergent set is a strict *subset* of
//!   `main`'s in both directions (3': 2,741 -> 2,480; 5': 2,938 -> 2,654).
//! - **0 still-divergent families change arity at all.** Every one of the 261 /
//!   284 moved families goes straight to `converged`, which is why the three
//!   `split_*` deltas sum exactly to the `converged` delta in each direction
//!   (215 + 25 + 21 = 261; 245 + 33 + 6 = 284). They come from arity 2 / 3 / 4 / 5
//!   in the counts 215 / 25 / 16 / 5 at 3' and 245 / 33 / 6 / 0 at 5'.
//!
//! The moved set is broad rather than local — 16 corpus shapes at each direction,
//! led by `pair-del-*` (134 at 3', 150 at 5'), `m3-all-*`, `m4-all-*` and
//! `mid-cds-*` — which is what a partitioner change looks like as against
//! #1536's, whose every moved row was `s00-c3{p,m}-cds-end-*`.
//!
//! **Every rank-1 counter is unchanged**, in both directions, as are
//! `non_idempotent_outputs` (4 / 4), `sequence_changed` (4 / 0), all three
//! refusal counters and `guard_violations`.
//!
//! **Two of the six rows #1599 promoted come back**, and they are removed from
//! `spec_corpus_regressions::the_codon_gate_splits_a_spanning_delins_its_own_members_do_not`
//! per that test's own instruction: `s00-c3{p,m}-m3-all-del-p1-sep3` (arity 2 at
//! 3', 3 at 5') and `s00-c3{p,m}-pair-del-del-p1-sep8` (arity 2 at both) are now
//! converged. `s00-c3{p,m}-m4-all-del-p1-sep2` is untouched and stays pinned
//! there — its spanning `delins` still stops where `general.md:34` puts it.
//!
//! # RE-BLESSED (4 of 5) — #1627, the alignment-only symbol
//!
//! `background/standards.md:39` footnotes the table's daggered `X` and `-` as
//! "used in alignment only", and the decided
//! `rulings[alignment-only-symbol-in-a-description]` rules that neither may
//! appear in a description. Ferro accepted `X` in every mode and re-emitted it.
//! It now refuses: strict at parse, lenient and silent at normalize, per the
//! decided `rulings[absolute-prohibition-enforcement-stage]`.
//!
//! **Two figures moved, both DOWN, identically in both directions:**
//!
//! | figure | was | now | direction |
//! |---|---|---|---|
//! | 3' and 5' `prohibition_violating_outputs` | 32 | **8** | improves by 24 |
//! | 3' and 5' `prohibited_conditional_accepted` | 40 | **16** | improves by 24 |
//!
//! **It is one population counted twice, not two.** All 24 rows are
//! `standards.md:39-alignment-only-symbols`, and they were the only entry in
//! either counter that this clause contributed:
//! `the_absolute_prohibitions_ferro_accepts_are_three_clauses_and_no_others`
//! and `every_prohibition_violating_output_is_a_re_emitted_prohibited_input`
//! pin the per-clause decomposition of each, and both drop exactly that one
//! entry. The residues are the shapes this change does not reach —
//! `checklist.md:20`'s 16 conditional rows (correct as they stand; see
//! `rulings[bare-transcript-intronic-position]`) and `checklist.md:16`'s 8
//! violating outputs (that is #1628).
//!
//! **`prohibited_absolute_accepted` deliberately does NOT move**, and its
//! staying at 32 is the useful half of the measurement. The corpus grades an
//! `X` row `Strength::Conditional` — the footnote states a scope rather than
//! using prohibitive words — while this axis counts an `X` *output* as an
//! absolute violation. The decided ruling records that the grading is moot
//! either way, because `general.md:48` admits only IUPAC-IUBMB symbols. So the
//! two counters that moved are the two the clause was actually in, and the
//! surviving 32 absolute acceptances are `checklist.md:32`'s `ins6` (24) plus
//! `checklist.md:16`/`:45`'s genomic offset and hyphen range (4 + 4) — the
//! other half of #1627 and #1628 respectively, neither touched here.
//!
//! **Every rank-2, idempotency and sequence-preservation figure is unchanged**,
//! in both directions, which is what a refusal of un-denotable inputs should
//! look like: the refused rows contributed no confluence family and no
//! sequence, so nothing they leave behind can move. `converged`, `split_*`,
//! `non_idempotent_outputs`, `sequence_changed`, `conflicts_accepted`,
//! `outputs_denoting_no_sequence`, `outputs_leaving_the_transcript` and
//! `guard_violations` all hold at their **#1649** values — this section was
//! first written against #1536's and re-measured on top of #1649 when the
//! branch was rebased, which is the only reason the sentence names a revision
//! at all: "unchanged" is a claim about a baseline, and the baseline moved.
//!
//! **Nothing to promote to `spec_corpus_regressions.rs`** — no row newly fails.
//! The row that stopped failing has an adjudicated-correct guard instead:
//! `corpus_prohibited_inputs::an_alignment_only_symbol_is_refused_in_every_mode_for_both_x_and_dash`,
//! which covers the embedded shapes (`delinsACGTX`, `delinsXACGT`,
//! `delinsACXGT`) the corpus does not generate.
//!
//! # RE-BLESSED (5 of 5) — #1716, the codon-frame merge's own span
//!
//! `merge_consecutive_edits`' per-member codon-frame predicate asked
//! `same_codon(prev_a.end, next_start)` — its left anchor's **right** edge —
//! while authorising a `delins` over `prev_a.start ..= next.end`. A codon is
//! three positions, so any anchor wider than one base makes that span four or
//! more and `general.md:35`'s "together affecting one amino acid" cannot hold;
//! `general.md:34` / `DNA/delins.md:17` then govern and the members stay
//! individual. Asking the span instead is the fix. (Sibling passes always did:
//! `apply_coding_codon_exception` and `coalesce_coding_frame_separation`.)
//!
//! **Six figures move, all rank 2, and `converged` is unchanged in BOTH
//! directions:**
//!
//! | figure | was (post-#1627) | now | direction |
//! |---|---|---|---|
//! | 3' `converged` | 9,402 | **9,402** | unchanged |
//! | 3' `split_two` | 2,225 | **2,223** | −2, both to `split_three` |
//! | 3' `split_three` | 223 | **226** | **+3 — rises** |
//! | 3' `split_more` | 32 | **31** | −1, to `split_three` |
//! | 5' `converged` | 9,228 | **9,228** | unchanged |
//! | 5' `split_two` | 2,461 | **2,462** | +1, from `split_more` |
//! | 5' `split_three` | 167 | **168** | +1, from `split_more` |
//! | 5' `split_more` | 26 | **24** | −2 |
//!
//! **Exactly three family instances move per direction, and they are named**,
//! because the net counters cannot tell an arity rise from an arity fall when
//! both happen. The full divergence row-id lists were dumped either side (the
//! `divergences` cap in [`measure`] raised locally, then restored) and diffed:
//!
//! - **0 families lose convergence, 0 gain it, and 0 become divergent that were
//!   not** — the divergent *set* is byte-identical in both directions (3': 2,480
//!   families; 5': 2,654). Only three members of it change arity.
//! - 3': `s00-c1-m4-all-del-p2-sep1` 4 -> 3 (improves);
//!   `s00-c3p-m4-all-del-p1-sep2` and `s00-c3m-m4-all-del-p1-sep2` 2 -> 3.
//! - 5': `s00-c1-m4-all-del-p2-sep1` 3 -> 2 and both `m4-all-del-p1-sep2` rows
//!   4 -> 3 (all three improve).
//!
//! **The two rows that rise are already promoted**, which is why this section
//! adds nothing to `spec_corpus_regressions.rs`:
//! `the_codon_gate_splits_a_spanning_delins_its_own_members_do_not` has pinned
//! `s00-c3{p,m}-m4-all-del-p1-sep2` since #1599. What changed there is its
//! *members'* answer, and it changed toward that test's own clause — the old
//! `c.[24del;30_33delinsA]` merged across `c.30..c.33`, codon 10 into codon 11.
//! The pin is updated in this commit with the reason on the row.
//!
//! **Every rank-1 counter is unchanged**, in both directions — `outputs`,
//! `declined`, `unparseable_outputs`, `outputs_denoting_no_sequence`,
//! `outputs_leaving_the_transcript`, `prohibition_violating_outputs` — as are
//! `non_idempotent_outputs` (4 / 4), `sequence_changed` (4 / 0), all three
//! refusal counters and `guard_violations`.
//!
//! # CORRECTED — what the transcript-leaving class actually violates
//!
//! An earlier form of the table above said `general.md:44` "exempts
//! junction-adjacent deletions/duplications from the 3' rule", implying the shift
//! into the intron is itself forbidden. **It is not, and citing it that way
//! points at an exception the spec withholds.** `general.md:44` defers to
//! `background/numbering.md`, and that file scopes the exception twice:
//!
//! - `:23` — "the 3' rule is not applied when there is a deletion/duplication
//!   around **exon/exon** junctions with identical nucleotides flanking the
//!   junction, **where shifting the variant 3' would place it in the next
//!   exon**."
//! - `:26` — "**NOTE**: this exception **does not apply** to a
//!   deletion/duplication around **exon/intron and intron/exon** junctions with
//!   identical nucleotides flanking the junction".
//!
//! The 371 rows are exactly the exon→intron case `:26` excludes, so the shift is
//! legal and ferro already blocks the exon→exon case the exception does cover.
//! What is invalid is the **accession**: `checklist.md:20` forbids a bare `NM_`
//! from naming an intronic position. So converging these rows on the clamped
//! in-exon answer would implement an exception the spec explicitly withholds —
//! and would silently revert #670. **Settled by #1704: re-parenting onto a genomic
//! wrapper, not refusal** — see the RE-BLESSED section below and
//! `Normalizer::reparent_junction_exit`.
//!
//! Two further corrections from the same investigation, both measured:
//!
//! - **The 3'/5' split is a claim about the CODE, not the corpus.** All three
//!   copies of the #670 gate (`normalize_cds`, `normalize_tx`, `normalize_rna`)
//!   are guarded on `ThreePrime` with no 5' mirror. Sweeping 16 placements with a
//!   matching preceding intron, the 5' direction never leaves the exon.
//! - **The strand is a confound, not a finding.** "Entirely minus-strand" is an
//!   artifact of the fixture: the provider writes a literal `GATTACA` intron and
//!   reverse-complements only the exon blocks, so against `AT` cores `G` never
//!   continues a run and `AA` always does. Holding the transcript-direction
//!   intron fixed makes both strands agree. See `defect_371_transcript_exit.rs`,
//!   whose `junction_provider` takes the intron as a parameter for this reason.
//!
//! # RE-BLESSED (3 of 3) — #1704, the junction exit's accession
//!
//! `outputs_leaving_the_transcript` goes **371 -> 0** at 3', and every one of
//! those 371 rows is now counted by the new
//! `outputs_intronic_under_a_genomic_wrapper`, at **371**. That identity is the
//! whole disclosure: **no row stopped reaching the check**, and no coordinate
//! moved. The fix renders a junction-crossing output as
//! `NC_SYNTH.1(NM_TEST.1):c.20+2del` instead of `NM_TEST.1:c.20+2del`.
//!
//! | figure | was | now | direction |
//! |---|---|---|---|
//! | 3' `outputs_leaving_the_transcript` | 371 | **0** | rank-1 class CLOSED |
//! | 3' `outputs_intronic_under_a_genomic_wrapper` | — | **371** | new counter; conformant, not a defect |
//! | 5' both of the above | 0 | **0** | unchanged — no 5' mirror of the #670 gate exists |
//!
//! **Every other counter is byte-identical in both directions**, measured rather
//! than argued: 58,552 outputs / 0 declined / 0 unparseable / 10 (3') and 18 (5')
//! denoting no sequence / 32 prohibition-violating / 9,141 converged / 2,440 /
//! 248 / 53 / 4 non-idempotent / 4 sequence-changed / 72 / 32 / 40 / 0 guard
//! violations. So nothing was traded for this.
//!
//! **The counter had to be re-scoped, and that is a finding rather than a
//! convenience.** The old counter was every `Denotation::Inexpressible` output,
//! and re-parenting moves that not at all: SPDI is positional and has no offset
//! notation, so `hgvs_to_spdi` declines an intronic position whatever accession
//! carries it. Run against the pre-#1704 census, the fix reports **371 before and
//! 371 after** — a change that made all 371 outputs conformant, invisible to the
//! instrument. The proxy and the clause had silently come apart, and the census
//! now asks `checklist.md:20`'s own question (is the accession bare?) directly
//! from the AST. See `names_bare_transcript_intronic` below, including why it is
//! written there rather than borrowed from `src/normalize/`.
//!
//! **Under `FERRO_PARTITION=canonical-coalesced` the same fix reads 389 -> 0**,
//! with 389 on the wrapper counter. That arm has no pin of its own and fails on
//! an unrelated pre-existing counter (`guard_violations` 0 -> 6), which is the
//! coalesced partitioner's own divergence and is unmoved by #1704 — verified by
//! re-running that arm with `src/normalize/` reverted.
//!
//! # The coding-axis merge counter — an INSTRUMENT, and the guard beside it is not one
//!
//! `guard_violations` above is a **negative guard**: its rows are the shape
//! rejected SVD-WG010 would have merged, so a violation there is a verdict.
//! `coding_axis_separation_two_or_more_merges` is not that. It counts a
//! population — coding-axis multi-member alleles whose members sit **two or
//! more** unchanged nucleotides apart and which normalization merged into fewer
//! members than were authored — and makes no claim that merging them is wrong.
//!
//! **Why that population is worth counting, and which records govern it, are
//! stated once**, on `spec_corpus::is_coding_axis_separation_two_or_more_shape`,
//! and are not repeated here. In short: `general.md:34` / `DNA/delins.md:17` are
//! what make it interesting, and three ruling records — two `decided`, one not —
//! govern different parts of it, so a rise is not one thing and neither is the
//! zero. Do not argue from the two clauses without reading those records.
//!
//! **It exists because `guard_violations` provably cannot count these rows.**
//! `is_svd_wg010_shape` admits `Genomic | NonCodingMultiExon` at a separation of
//! exactly one with exactly two members; this counter's rows are
//! `CodingSingleExon | CodingMultiExon` at two or more, with two or more
//! members. The two domains are disjoint twice over, so a merge on the coding
//! axis leaves `guard_violations` at zero however many rows it moves — which is
//! why the guard is left exactly as it is and this is a second counter beside
//! it, not a widening of it. That disjointness is pinned as a property of the
//! corpus by `the_coding_axis_merge_population_is_what_it_says_it_is`.
//!
//! **Pinned at 0 / 0, measured over 997 rows in each direction.** The zero is a
//! real negative result, not an absence. Re-measured on `origin/main`
//! `e98fa77e`, which is **after #1698** — see the paragraph below. (The table
//! was first derived at `1ea75334` and reproduced cell-for-cell after the
//! rebase; `git diff 1ea75334 origin/main -- src/` is empty, so the arms could
//! not have moved.)
//!
//! | evaluation arm (`FERRO_PARTITION`) | 3' | 5' |
//! |---|---|---|
//! | unset — the shipped rule | **0** of 997 | **0** of 997 |
//! | `shadow` | 0 | 0 |
//! | `canonical` | 1 | 1 |
//! | `canonical-coalesced` | **3** | **3** |
//!
//! The 3 are one row at a separation of 3 and two at a separation of 2, all of
//! them two-member designs collapsing to one member — every one a row
//! `guard_violations` cannot see. The `canonical` arm's single residual is
//! `s01-c1-pair-inv-del-p4-sep3` (`c.[9_12inv;16_19del]` -> `c.9_19delinsGTCAGTA`),
//! which is not the coalesce pass and is also one of the coalesced arm's 3.
//!
//! **That table read 122 in the `canonical-coalesced` arm when this change was
//! written, and #1698 took it to 3.** `fix(normalize): require a gap-bearing
//! insert before the delins.md:44-47 merge` landed on `merge.rs` between the
//! branch point and the rebase, removing 119 of the 122 — a 97.5% collapse in
//! the figure this change was originally motivated by, with **no test failing**,
//! because the figures live in prose. That is the hazard this module is supposed
//! to be about, realised on the module itself; it is recorded here rather than
//! quietly restated. Attribution was measured, not assumed: re-running all four
//! arms with the corrected numerator (below) reproduces 0 / 0 / 1 / 3, so the
//! movement is #1698's and none of it is the numerator's.
//!
//! **The scale argument this counter was introduced with does not survive that,
//! and the disjointness argument does.** At 122 rows the case was "the canonical
//! arm's largest blocking family is invisible to this axis". At 3 it is not a
//! large family. What is unchanged is that `guard_violations` *cannot* count
//! them at any size — a structural fact about `is_svd_wg010_shape`'s domain, not
//! a claim about how many rows there happen to be — and that a counter is how
//! that stays true after the next change to `merge.rs`, which is exactly what
//! #1698 demonstrated.
//!
//! **The counter was watched moving, in both directions, before being trusted.**
//! Stubbing `merge::coalesce_payload_alignment_split` to decline unconditionally
//! takes the `canonical-coalesced` arm to the `canonical` figure exactly.
//! Removing the arm gate so the pass runs on the shipped rule takes the shipping
//! arm off zero. (Both figures were measured pre-#1698 — 122 -> 1, and 0 -> 110
//! (3') / 113 (5') — so quote the direction, not the numbers.) The shipped
//! zero is the pass not running, not the instrument failing to look. Since
//! neither mutation ships, `the_coding_axis_merge_counter_can_observe_a_merge`
//! is the positive control that keeps a blinded numerator from passing CI.
//!
//! # The honest-zero discipline is enforced, not described
//!
//! Every property whose denominator is zero fails loudly as **VACUOUS** rather
//! than passing. A corpus that stopped generating conflict rows would otherwise
//! report `0 accepted` — indistinguishable from perfect refusal, and exactly what
//! #1456 did three times (`0 of 18,432`, reported as evidence of neutrality when
//! the corpus simply could not build a conflicting allele).
//!
//! # No fixture, deliberately
//!
//! The corpus is enumerated **at test time** from the library. There is no
//! generated artifact to stage, no CI cache entry to add, and no committed data
//! file — so this runs on every `cargo test` and every CI shard with no wiring at
//! all. `examples/generate_spec_conformance_corpus.rs` exists for the other half
//! of the job: checking the committed rule inventory's citations against the spec
//! checkout, which needs the submodule and therefore belongs in CI's
//! `spec-fixtures` job.
//!
//! # Tier 2: a failing row is promoted, not left in the corpus
//!
//! A corpus row can silently stop being generated — that is the #1456/#1460/#1478
//! class exactly. So every row that fails here is promoted to a **named committed
//! test** in `tests/it/spec_corpus_regressions.rs`, with the governing
//! `file.md:line`. Corpus for breadth; named tests for anything that ever went
//! red.

use std::collections::{BTreeMap, BTreeSet};

use rayon::prelude::*;

use ferro_hgvs::conformance::spec_corpus::{
    corpus, denotation_of, denoted_by, CorpusBounds, Denotation, Row, RowKind, SpecCorpus, Strength,
};
use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::hgvs::interval::{Interval, UncertainBoundary};
use ferro_hgvs::hgvs::location::{CdsPos, TxPos};
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

// ---------------------------------------------------------------------------
// The pinned baselines
// ---------------------------------------------------------------------------

/// The corpus's own shape, pinned separately from the properties measured over
/// it.
///
/// A generator change that halved the corpus would otherwise improve every rate
/// below while measuring less — the failure mode `#1460` is, where a corpus that
/// could not reach `MAX_SPLIT_BLOCK` reported a length-gated change as neutral.
const SHAPE: CorpusShape = CorpusShape {
    designs_considered: 16_650,
    rows: 12_946,
    spellings: 58_552,
    family_rows: 11_882,
    single_rows: 804,
    conflict_rows: 96,
    prohibited_rows: 164,
    multi_member_rows: 12_496,
    guarded_rows: 210,
    coding_axis_separation_two_or_more_rows: 997,
};

/// The 3'-direction census, pinned. See the module docs for which way each
/// figure is allowed to move.
///
/// # Measured on `main`, and it measures `main`
///
/// This branch adds a corpus, an axis test and a regression suite; it touches no
/// file under `src/normalize/`, `src/hgvs/`, `src/spdi/` or `src/project/`. So
/// every figure here is `main`'s behaviour observed for the first time, not
/// behaviour this change introduced. Re-measured against `main` @ 5b54b2b6, and
/// identical row-for-row to the corpus's first measurement on `main` @ 35de96c8.
///
/// A pin measured somewhere the branch does not reproduce is red on every run
/// and can therefore detect nothing — the instrument would be permanently blind
/// in exactly the way it exists to prevent. That is the whole reason these are
/// pinned at the measured values rather than at zero.
///
/// The confluence figure is the headline and it is **not** good: this corpus is
/// deliberately far harsher than `cis_confluence_axis`'s, because it varies
/// member geometry, transcript geometry and scale rather than holding all three
/// fixed. Read the three asserted-zero counters before the headline, per the
/// module docs.
///
/// **The "measures `main`" heading above is now historical** — three normalizer
/// changes have re-blessed figures here (#1599, then #1536, then #1649), and the
/// four confluence figures below are measured against the third rather than
/// against the corpus PR's `main`. The module docs carry a RE-BLESSED section per
/// change, naming each figure, which way it moved, and the row ids behind it.
pub(crate) const THREE_PRIME: Census = Census {
    // -- provenance --
    //
    // LENIENT, because `NormalizeConfig::default()` substitutes
    // `ErrorConfig::lenient()`. Every figure below that is downstream of
    // normalization is therefore a lenient-mode figure;
    // `corpus_prohibited_inputs.rs` re-measures the refusal counters in strict
    // mode and gets different numbers. `attempted` and `unparseable_outputs`
    // are the exception and carry no mode at all — see the field's own docs.
    measured_under: ErrorMode::Lenient,
    // -- validity (rank 1) --
    outputs: 58_552,
    declined: 0,
    unparseable_outputs: 0,
    outputs_denoting_no_sequence: 10,
    // #1704: was **371**, and every one of those rows is now on the line below.
    // The class did not shrink — it stopped being a violation, because the
    // junction-crossing answer is now rendered against the genomic reference the
    // crossing itself resolved. 371 = 0 + 371 is the accounting identity that
    // makes this a fix rather than a re-partition.
    outputs_leaving_the_transcript: 0,
    outputs_intronic_under_a_genomic_wrapper: 371,
    // Re-blessed DOWN by #1627: `standards.md:39`'s 24 `X` rows are refused
    // rather than re-emitted. The residual 8 are `checklist.md:16`'s genomic
    // offsets, which is #1628. See the module docs' fourth RE-BLESSED section.
    //
    // Independent of the #1704 split directly above: #1627 refuses rows on the
    // alignment-symbol ground, #1704 re-parents rows on the junction-crossing
    // ground, and the two touch disjoint row sets.
    prohibition_violating_outputs: 8,
    // -- confluence (rank 2) --
    //
    // Re-blessed by #1649's two-deletion alignment, on top of #1536's re-bless
    // and #1599's before it — **measured on this branch rebased onto both**, not
    // composed from the three changes' deltas. The row-level diff is in the
    // module docs' third RE-BLESSED section.
    //
    // 261 families newly converge and **zero** lose convergence, rise in arity,
    // or become divergent that were not — measured by diffing the full
    // divergence row-id lists either side, not inferred from the net counters
    // (which cannot tell net from gross). Every moved family goes straight to
    // `converged`, which is why the three deltas below sum exactly to
    // `converged`'s: 215 + 25 + 21 = 261.
    //
    // **Re-blessed by #1616, and this is the branch's headline figure.**
    // Deleting the input-relative weight bound
    // (`rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`) lets a
    // multi-member spelling be re-derived from the sequence instead of being
    // handed back, so families whose spellings previously reached their own
    // inputs now reach one answer: converged 9 402 -> 10 995, +1 593. The three
    // split counters fall by 1 548 + 42 + 3 = 1 593, exactly the gain, so every
    // family that moved went straight to `converged` and none left the corpus.
    //
    // **Then 18 of them are given back**, by the `general.md:34` gate that takes
    // `guard_violations` to zero below: a block past `MAX_SPLIT_BLOCK` cannot be
    // split, so off `c.` the derivation declines and the spanning and split
    // spellings of those 18 classes stay two strings. `split_two` rises by
    // exactly 18 and nothing else moves. Net against `origin/main`: +1 575.
    // Rule 2 outranks rule 3, so a preference clause is paid for in confluence
    // and not the other way round.
    converged: 10_977,
    split_two: 695,
    split_three: 181,
    split_more: 29,
    underdetermined: 0,
    // -- idempotency --
    //
    // Re-blessed DOWN: the three `scale-c3p-sep{120,128,136}-del-del` rows now
    // settle in one pass. The four `cds-end` families remain.
    non_idempotent_outputs: 4,
    // -- sequence preservation --
    sequence_changed: 4,
    // -- refusal --
    conflicts_accepted: 72,
    prohibited_absolute_accepted: 32,
    // Re-blessed DOWN by #1627, same 24 rows as `prohibition_violating_outputs`
    // above. The residual 16 are `checklist.md:20`'s bare-transcript intronic
    // rows, which lenient is CORRECT to accept — see
    // `rulings[bare-transcript-intronic-position]`.
    prohibited_conditional_accepted: 16,
    // -- negative guards --
    //
    // **Zero, and measured rather than structural.** The same instrument read
    // 18 on this branch's first commit and 2 on its second, so it can report
    // non-zero over this corpus and this denominator; see the module docs.
    //
    // 18 while the shipping partitioner applied `DNA/delins.md:44-47`'s
    // payload-coincidence carve-out on every axis; 16 closed when that carve-out
    // was scoped to `c.` per
    // `rulings[delins-payload-coincidence-carve-out-is-coding-dna-scoped]`, and
    // the last 2 — 2 065-base blocks over `MAX_SPLIT_BLOCK` (1024) — closed when
    // `canonicalize_from_sequence` stopped emitting the whole block the length
    // short-circuit hands back **unexamined**. `MAX_SPLIT_BLOCK` itself is
    // untouched: raising it is a cost question and is still open.
    guard_violations: 0,
    // The denominator, pinned so the zero above cannot be a structural silence.
    // Equals `CorpusShape::guarded_rows` (210) exactly, which is the check that
    // every guarded row's authored spelling re-parsed and was actually decided.
    guard_evaluations: 210,
    // -- instruments (#1710), kept across the #1616 rebase --
    coding_axis_separation_two_or_more_rows: 997,
    coding_axis_separation_two_or_more_merges: 0,
};

/// The rows that still merge a frameless separation of one, named.
///
/// **Empty, and kept rather than deleted.** The list is the mechanism that stops
/// a residual rotting into a silent permanent exemption:
/// `the_residual_guard_violations_are_exactly_these_rows` asserts the set both
/// ways, so a row entering FAILS and a listed row that starts conforming FAILS
/// too. Empty, it is the strongest form of that assertion — *no* row may merge a
/// frameless separation of one — and it is where the next regression will land
/// with its id already printed.
///
/// It last held `scale-g-block1032-delins-del` and `scale-g-block1032-inv-del`,
/// whose 2 065-base blocks are the only ones in the corpus over
/// `merge::MAX_SPLIT_BLOCK`; both closed when `canonicalize_from_sequence`
/// stopped emitting the block that cap hands back unexamined. They are also why
/// no smaller-scale corpus can see this class at all (#1460), so do not read a
/// zero here as covering a corpus that cannot build the shape.
pub(crate) const RESIDUAL_GUARD_VIOLATIONS: &[&str] = &[];

/// The 5'-direction census, pinned.
///
/// `--direction 5prime` is a supported public option and confluence is a
/// property of the normalizer rather than of one shuffle direction, so it is
/// measured in full. Two directions landing on materially different numbers
/// would mean a fix was treating a symptom of the shuffle rather than the
/// partitioner.
///
/// Measured on the same base as [`THREE_PRIME`], and re-blessed by the same
/// five changes — see the module docs' five RE-BLESSED sections.
///
/// The two directions agreeing on every rank-1 counter except
/// `outputs_leaving_the_transcript` is itself the cross-check the two-direction
/// measurement exists for: a figure that appears at 3' and not at 5' is a claim
/// about the code (all three copies of the #670 junction gate are guarded on
/// `ThreePrime` with no 5' mirror), not about one shuffle direction's luck.
pub(crate) const FIVE_PRIME: Census = Census {
    // Lenient, as at 3' — see [`THREE_PRIME`].
    measured_under: ErrorMode::Lenient,
    outputs: 58_552,
    declined: 0,
    unparseable_outputs: 0,
    outputs_denoting_no_sequence: 18,
    outputs_leaving_the_transcript: 0,
    // Zero on both sides of the split at 5', and for the same structural reason
    // the sibling counter was zero here before #1704: all three copies of the
    // #670 junction gate are guarded on `ThreePrime` with no 5' mirror, so the 5'
    // direction never produces an intronic output at all — wrapped or bare. Read
    // as the negative control for the 3' pin, not as a second measurement of it.
    outputs_intronic_under_a_genomic_wrapper: 0,
    // Re-blessed DOWN by #1627, by the same 24 rows as at 3'. Validity does not
    // depend on shuffle direction, so the two directions moving identically is
    // the cross-check rather than a coincidence.
    prohibition_violating_outputs: 8,
    // Unmoved by #1599 or by #1536 — no 5' family changed which side of the
    // converged/split line it sat on in either. **#1649 moves it for the first
    // time**, by 284, and the same three measured zeros hold as at 3': no family
    // loses convergence, rises in arity, or becomes divergent that was not.
    // #1627 does not move it either: a refused row contributes no family.
    // **Re-blessed by #1616**, the mirror of [`THREE_PRIME`]'s move: 9 228 ->
    // 10 789, +1 561, with `split 2`/`3`/`4+` falling by 1 516 + 43 + 2 = 1 561.
    // The two directions gain by different amounts (1 593 against 1 561) for the
    // reason every earlier re-bless in this file records: the partition is
    // decided the same way in both, and the follow-on shuffle that it enables
    // runs 3' or 5', so only some landings coincide with the sibling spelling's.
    converged: 10_775,
    // Re-blessed by #1649, rank-2 only, same as [`THREE_PRIME`] — and measured on
    // the rebased branch rather than composed. Every moved family goes straight
    // to `converged`, so these three deltas sum exactly to `converged`'s:
    // 245 + 33 + 6 = 284.
    //
    // The two directions moving by different amounts is expected rather than
    // suspicious, for the reason #1536's re-bless already recorded: the partition
    // is decided the same way in both, but the *follow-on* shuffle it enables
    // runs 3' or 5'. Here 5' moves the larger set (284 against 261), the reverse
    // of #1536.
    split_two: 959,
    split_three: 124,
    split_more: 24,
    underdetermined: 0,
    non_idempotent_outputs: 4,
    sequence_changed: 0,
    conflicts_accepted: 72,
    prohibited_absolute_accepted: 32,
    // Re-blessed DOWN by #1627, by the same 24 rows as at 3'.
    prohibited_conditional_accepted: 16,
    // Zero, as at 3'. The guard is a property of the partition rather than of a
    // shuffle direction, so the two directions agreeing here is the cross-check
    // — and both read 18, then 2, then 0 through the same two commits.
    guard_violations: 0,
    // As at 3', and for the same reason. The two directions agreeing is the
    // cross-check: a direction that stopped evaluating would show up here rather
    // than as a quietly unchanged zero above.
    guard_evaluations: 210,
    // -- instruments (#1710), kept across the #1616 rebase --
    coding_axis_separation_two_or_more_rows: 997,
    coding_axis_separation_two_or_more_merges: 0,
};

/// The corpus's shape, independent of any property measured over it.
#[derive(Debug, Default, PartialEq, Eq)]
struct CorpusShape {
    designs_considered: usize,
    rows: usize,
    spellings: usize,
    family_rows: usize,
    single_rows: usize,
    conflict_rows: usize,
    prohibited_rows: usize,
    multi_member_rows: usize,
    guarded_rows: usize,
    /// The coding-axis merge instrument's population — see
    /// [`Census::coding_axis_separation_two_or_more_merges`]. Pinned here as well as in the census
    /// because the two answer different questions: this is how many rows the
    /// *generator* built, and the census figure is how many the *measurement*
    /// reached.
    coding_axis_separation_two_or_more_rows: usize,
}

impl CorpusShape {
    fn of(built: &SpecCorpus) -> Self {
        let by_kind = built.by_kind();
        Self {
            designs_considered: built.designs_considered,
            rows: built.rows.len(),
            spellings: built.spellings(),
            family_rows: by_kind.get(&RowKind::Family).copied().unwrap_or(0),
            single_rows: by_kind.get(&RowKind::Single).copied().unwrap_or(0),
            conflict_rows: by_kind.get(&RowKind::Conflict).copied().unwrap_or(0),
            prohibited_rows: by_kind.get(&RowKind::Prohibited).copied().unwrap_or(0),
            multi_member_rows: built.multi_member_rows(),
            guarded_rows: built.by_negative_guard().values().sum(),
            coding_axis_separation_two_or_more_rows: built
                .coding_axis_separation_two_or_more_rows(),
        }
    }
}

/// What one direction's run found. Every field is pinned.
///
/// `pub(crate)`, with `pub(crate)` fields, so sibling modules can reference the
/// pins instead of re-typing them. `corpus_prohibited_inputs.rs` decomposes
/// three of these figures and used to hard-code them; a literal copy of a pin
/// goes stale silently at the next re-bless, and its failure message then
/// accuses the implementation of a defect that is really a stale constant.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct Census {
    /// The error mode the **normalizer** behind every figure below was built
    /// with (#1629).
    ///
    /// **It does not reach the parse half.** Every spelling is read with the
    /// bare `parse_hgvs`, which constructs no `InputPreprocessor` and so applies
    /// no `ErrorConfig` at all — not lenient's repairs, not strict's refusals
    /// (#1632, pinned by `issue_1632_parse_entry_applies_no_mode.rs`). So
    /// `unparseable_outputs` and `attempted` are mode-*less* figures, and only
    /// the counters downstream of normalization are lenient-mode figures. "No
    /// mode" is a third thing rather than a synonym for strict, which is exactly
    /// why it has to be said instead of inferred from this field.
    ///
    /// **Not decoration.** `corpus_prohibited_inputs.rs` re-measures three of
    /// the four refusal counters in strict mode and two of them move
    /// (`prohibited_absolute_accepted`'s 32 stay accepted) —
    /// so a census compared against one taken under a different mode is a
    /// category error, and until this field existed nothing in the artifact or
    /// the pin said which mode either was. The decided ruling
    /// `bare-transcript-intronic-position` has to end with the sentence "so that
    /// counter is a lenient-mode figure" for exactly this reason.
    ///
    /// Set from [`measurement_config`], never written by hand, and carried in
    /// the same `assert_eq!(census, pinned)` as every other field — so changing
    /// what the axis measures under fails here, naming the mode, instead of
    /// silently re-basing seventeen numbers.
    ///
    /// [`Default`] gives `Strict` (the enum's own default) rather than the
    /// lenient value the axis actually uses, which is deliberate: an accumulator
    /// that was never stamped must not read as a valid claim.
    pub(crate) measured_under: ErrorMode,
    /// Spellings **attempted**, across every row — incremented before
    /// `parse_hgvs` runs, so it counts every spelling the sweep reached.
    ///
    /// Deliberately not "spellings normalized": a spelling that fails to parse,
    /// or that normalization declines, is counted here **and** in
    /// [`Census::declined`]. The two are therefore NOT disjoint buckets, and the
    /// `report` line's `"{} outputs, {} declined"` must not be read as though
    /// they were — subtracting one from the other is wrong the moment `declined`
    /// leaves zero. Both pins hold `declined: 0` today, so the two readings
    /// coincide and nothing is currently mis-stated; this doc is what stops a
    /// future re-bless with a non-zero `declined` being read wrongly.
    ///
    /// Counting attempts rather than successes is the useful choice: it makes
    /// `outputs == SHAPE.spellings` a cross-check between the shape pin and the
    /// census pin, which a success-only counter could not provide.
    pub(crate) outputs: usize,
    /// Spellings normalization declined or panicked on, in a stratum where an
    /// output was expected. Counted apart from a divergence, so one cannot hide
    /// the other.
    pub(crate) declined: usize,
    /// Outputs `parse_hgvs` rejects. **Asserted at zero** — rank-1 invalid.
    pub(crate) unparseable_outputs: usize,
    /// Outputs whose members claim overlapping territory, so the description
    /// denotes no sequence. Rank-1 invalid.
    pub(crate) outputs_denoting_no_sequence: usize,
    /// Outputs naming an intronic offset the input did not name, **on a bare
    /// transcript accession**. Rank-1 invalid, on `checklist.md:20`: an `NM_`
    /// "can only be used to describe variants in introns using a `c.` prefix when
    /// a genomic reference sequence is given", so the description is not
    /// expressible against the accession it carries.
    ///
    /// **The coordinate is not the defect, and the counter was re-scoped in
    /// #1704 to say so.** It used to be every `Denotation::Inexpressible`
    /// output — i.e. every intronic output, wrapper or not — which conflated a
    /// spec violation with a limit of the oracle. `general.md:44`'s exemption of
    /// "deletions/duplications around exon/exon junctions" was read as requiring
    /// the shift to stop at the junction; `background/numbering.md:26` explicitly
    /// **withholds** that exemption from exon/intron junctions, so the shift is
    /// required and only the accession was ever wrong. The conformant outputs the
    /// old counter also swept up are now [`Self::outputs_intronic_under_a_genomic_wrapper`].
    pub(crate) outputs_leaving_the_transcript: usize,
    /// Outputs naming an intronic offset under the genomic wrapper
    /// `checklist.md:20` asks for (`NC_…(NM_…):c.10+2del`). **Conformant** — this
    /// is not a defect counter, and it is not asserted at zero.
    ///
    /// It exists so the class stays *accounted for* rather than merely stopping
    /// being counted. #1704 moved every one of the 371 rows the sibling counter
    /// held from the defect side to this one, and the two figures summing to the
    /// old 371 is what shows the fix closed the class rather than shrinking the
    /// population that reached the check.
    ///
    /// The oracle cannot verify these outputs' bases — SPDI is positional and has
    /// no offset notation, so `hgvs_to_spdi` declines an intronic position
    /// whatever accession carries it. That was equally true before, and is a
    /// property of the instrument rather than of the descriptions:
    /// `defect_371_transcript_exit` is where the coordinates themselves are
    /// pinned, string by string.
    pub(crate) outputs_intronic_under_a_genomic_wrapper: usize,
    /// Outputs violating an absolute textual prohibition.
    ///
    /// **Ratchet-pinned at 8, not asserted at zero** — this doc once said the
    /// latter while both censuses pinned a non-zero figure, which reads as
    /// "ferro emits no prohibited output" to anyone who does not go and check
    /// the pin. It emits 8, and `corpus_prohibited_inputs` decomposes them:
    /// 8 × `checklist.md:16`, every one a re-emission of an input that already
    /// violated the clause rather than a violation ferro manufactured.
    ///
    /// It was 32 until #1627, whose 24 × `standards.md:39` rows are now refused
    /// rather than re-emitted. Quote the pin, not this prose, if the two ever
    /// disagree again — a count in a doc comment is the thing that goes stale.
    pub(crate) prohibition_violating_outputs: usize,
    /// Families whose spellings all reached ONE output. The goal is for this to
    /// equal `family_rows`.
    pub(crate) converged: usize,
    /// Families reaching exactly two distinct outputs.
    pub(crate) split_two: usize,
    /// … exactly three.
    pub(crate) split_three: usize,
    /// … four or more.
    pub(crate) split_more: usize,
    /// Families left with fewer than two spellings that normalized at all, so
    /// confluence is not decidable. Counted rather than skipped: a rise here
    /// would hollow the axis out while every other figure improved.
    pub(crate) underdetermined: usize,
    /// Outputs that are not their own fixed point.
    pub(crate) non_idempotent_outputs: usize,
    /// Outputs whose applied sequence differs from their row's.
    ///
    /// **Ratchet-pinned, and it is not zero in both directions**: 4 at 3' and 0
    /// at 5'. The zero is the 5' census alone, so a doc claiming the counter is
    /// asserted at zero describes half the axis. The four are the CDS-end flush
    /// pairs `spec_corpus_regressions` dissects, and they are the most serious
    /// class the axis measures — a normalization that changes the sequence its
    /// input denoted.
    pub(crate) sequence_changed: usize,
    /// Conflicting alleles — two members claiming one territory — the
    /// implementation accepted instead of refusing.
    pub(crate) conflicts_accepted: usize,
    /// Shapes the spec prohibits in as many words ("is not allowed", "is not
    /// correct") that the implementation accepted.
    pub(crate) prohibited_absolute_accepted: usize,
    /// Shapes a conditional clause prohibits that the implementation accepted. A
    /// finding to adjudicate rather than a regression.
    pub(crate) prohibited_conditional_accepted: usize,
    /// Rows whose output implements a REJECTED consultation proposal — a
    /// separation floor of two on a frameless axis, which is SVD-WG010
    /// (`consultation/SVD-WG010.md:8`, "The proposal has been **rejected**").
    pub(crate) guard_violations: usize,
    /// How many times the negative guard above was actually **evaluated**.
    ///
    /// The denominator of [`Census::guard_violations`], and pinned separately
    /// from it on purpose. `CorpusShape::guarded_rows` says how many rows *carry*
    /// a guard; this says how many reached the point where it could fire, which
    /// is a strictly narrower claim — the guard runs only on a row's authored
    /// spelling and only when that spelling's output re-parses. Without it a run
    /// in which every guarded output became unparseable would report
    /// `guard_violations: 0` with nothing having been checked, which is exactly
    /// the shape of a structural zero read as a real one.
    ///
    /// `unparseable_outputs` being pinned at 0 covers most of that, but it does
    /// so from two counters and a reader has to compose them. This asserts it
    /// directly, at the site.
    pub(crate) guard_evaluations: usize,
    /// The denominator [`Self::coding_axis_separation_two_or_more_merges`] is
    /// `n of` — coding-axis rows separated by two or more unchanged nucleotides
    /// that the measurement actually reached and evaluated.
    ///
    /// Pinned, and asserted non-zero. A separate figure from
    /// `SHAPE.coding_axis_separation_two_or_more_rows` on purpose: that one says the
    /// generator built the rows, this one says the run *normalized* them. A
    /// measurement whose every eligible row declined would report `0` merges
    /// with the generator-side denominator still healthy, which is the shape
    /// of vacuity this axis exists to refuse.
    pub(crate) coding_axis_separation_two_or_more_rows: usize,
    /// Coding-axis multi-member alleles whose authored spelling normalized to
    /// **fewer members than were authored** — i.e. the members were merged
    /// across the two or more unchanged nucleotides between them.
    ///
    /// # An instrument. It counts a population; it does not judge it
    ///
    /// This counter reports how often a merge of this shape happens. It does
    /// **not** assert the merge is wrong. Contrast [`Self::guard_violations`],
    /// which counts a shape whose merge implements a **rejected** proposal and
    /// is therefore a verdict. A count of `n` here is a measurement of how large
    /// the merged population is, not a claim that `n` rows are defective.
    ///
    /// # Everything else about this population is stated ONCE, and not here
    ///
    /// Why the floor is two; what the sub-floor population contains and which
    /// part of it is a known defect; which three ruling records govern which
    /// part of it, and therefore what a rise and a zero each mean; and why this
    /// is a second counter rather than a widened guard — all of that is on
    /// `spec_corpus::is_coding_axis_separation_two_or_more_shape`.
    ///
    /// It is deliberately not restated here. The first revision of this change
    /// carried that argument in seven places, and two of the copies ended up
    /// contradicting each other 406 lines apart in this file — one saying ferro
    /// implements the codon exception and it accounts for the sub-floor merges,
    /// the other retracting exactly that. A pointer cannot drift.
    ///
    /// # The numerator, and what it can and cannot see
    ///
    /// [`coding_axis_merge_observed`] is the test: `members_of(output).len() <
    /// row.members`. It was `< 2`, which on the 313 of 997 flagged rows carrying
    /// three or four members could only ever see a collapse **all the way** to
    /// one member — a 4-member row merged to 2 counted as zero. Below the floor
    /// that blindness is 25 of 41 merges.
    ///
    /// It still cannot see a merge whose output is unparseable, since the count
    /// is taken on a parsed output; `unparseable_outputs` is the counter for
    /// that and is pinned at 0.
    ///
    /// # Pinned at zero in both directions, and that zero is a measurement
    ///
    /// The pass that merges this population under the canonical-coalesced
    /// evaluation arm — `merge::coalesce_payload_alignment_split`, reached only
    /// when `payload_coalesce_applies` sees `PartitionRule::CanonicalCoalesced`
    /// — does not run on the shipped rule, so the shipping arm merges none of
    /// these rows. That is a mechanical reason rather than an empirical
    /// coincidence, and the counter is separately shown able to move: see the
    /// module docs' arm table and
    /// [`the_coding_axis_merge_counter_can_observe_a_merge`].
    pub(crate) coding_axis_separation_two_or_more_merges: usize,
}

// ---------------------------------------------------------------------------
// Absolute prohibitions, checked on the output
// ---------------------------------------------------------------------------

/// Whether `output` violates a prohibition the spec states in as many words, and
/// which.
///
/// Deliberately a **closed, small** list of textual checks rather than a general
/// validator: each entry cites a clause that uses prohibitive words, so a
/// non-zero count is a conformance defect and not an adjudication. The parser
/// already rejects most of these on *input*; the question here is whether the
/// normalizer can emit one.
fn violated_prohibition(output: &str) -> Option<&'static str> {
    // `general.md:96`: "spaces are *not* permitted in any HGVS description".
    if output.contains(' ') {
        return Some("general.md:96 — spaces are not permitted in any HGVS description");
    }
    // `background/standards.md:39` footnotes `X` and `-` as "used in alignment
    // only", so neither is a nucleotide a description may state. `X` cannot
    // appear anywhere in a legal DNA/RNA description.
    if output.contains('X') || output.contains('x') {
        return Some("standards.md:39 — `X` is an alignment-only symbol, not a nucleotide");
    }
    // `checklist.md:16`: a genomic reference "can not have nucleotides with
    // additions like a `+`, `-`, or `*`". Checked only on the description body,
    // since an accession may legitimately hold neither.
    if let Some(body) = output.split_once(":g.").map(|(_, body)| body) {
        if body.contains('+') || body.contains('*') || body.contains('-') {
            return Some("checklist.md:16 — a `g.` description admits no `+`/`-`/`*` offset");
        }
    }
    None
}

/// Members of a parsed description, for the negative guard and the applier.
fn members_of(variant: &HgvsVariant) -> Vec<HgvsVariant> {
    match variant {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single.clone()],
    }
}

/// The coding-axis merge instrument's **numerator**: did normalization collapse
/// `row`'s authored members into fewer members?
///
/// # It is `< row.members`, not `< 2`, and the difference is 31.4% of the corpus
///
/// The sibling `guard_violations` test may write `< 2` because
/// `is_svd_wg010_shape` pins `kinds.len() == 2`, so for its rows "fewer than two
/// members" and "merged" are the same statement. This population is wider —
/// `MEMBER_COUNTS = [2, 3, 4]` — and there `< 2` means "collapsed **all the way**
/// to one member". A 4-member row merged to 2 is two merges across two or more
/// unchanged nucleotides and would have been counted as **zero**. Measured over
/// the flagged rows: `{2: 684, 3: 159, 4: 154}`, so **313 of 997 (31.4%)** were
/// blind to a partial merge under the narrower test.
///
/// That is the property the doc comments claim ("the members were merged across
/// the two or more unchanged nucleotides between them"), so measuring anything
/// narrower would be an assertion weaker than the stated contract.
///
/// Extracted as a named function rather than inlined so
/// `the_coding_axis_merge_counter_can_observe_a_merge` can exercise it directly
/// — see that test for why a zero pin needs a positive control.
fn coding_axis_merge_observed(row: &Row, output: &HgvsVariant) -> bool {
    members_of(output).len() < row.members
}

/// Whether either endpoint of `interval` names an intronic offset.
///
/// Both `Single` and `Range` (uncertain-breakpoint) boundaries are inspected, and
/// an unknown (`?`) offset counts — it is still an intronic position, just one
/// whose magnitude is unspecified.
fn interval_is_intronic<T>(
    interval: &Interval<T>,
    is_intronic: impl Fn(&T) -> bool + Copy,
) -> bool {
    let boundary = |b: &UncertainBoundary<T>| match b {
        UncertainBoundary::Single(mu) => mu.inner().is_some_and(is_intronic),
        UncertainBoundary::Range { start, end } => {
            start.inner().is_some_and(is_intronic) || end.inner().is_some_and(is_intronic)
        }
    };
    boundary(&interval.start) || boundary(&interval.end)
}

/// Whether `output` names an intronic position under a **bare** transcript
/// accession — the form `checklist.md:20` refuses.
///
/// > a (non-)coding DNA reference sequence "can only be used to describe variants
/// > in introns using a `c.` prefix when a genomic reference sequence is given"
///
/// so `NC_000023.10(NM_004006.2):c.94-2A>G` is the conformant spelling and
/// `NM_004006.2:c.94-2A>G` is not.
///
/// **Written here from the AST rather than borrowed from `src/normalize/`.** The
/// normalizer has its own reading of this clause
/// (`intronic_on_bare_transcript_warning`, the W4007 rung), and calling it from
/// the census would make the instrument agree with the code by construction —
/// the census exists to judge ferro from outside. The two readings agreeing is a
/// result; sharing one implementation would make it an identity.
fn names_bare_transcript_intronic(output: &HgvsVariant) -> bool {
    members_of(output).iter().any(|member| {
        let (accession, intronic) = match member {
            HgvsVariant::Cds(v) => (
                &v.accession,
                interval_is_intronic(&v.loc_edit.location, CdsPos::is_intronic),
            ),
            HgvsVariant::Tx(v) => (
                &v.accession,
                interval_is_intronic(&v.loc_edit.location, TxPos::is_intronic),
            ),
            _ => return false,
        };
        intronic && accession.genomic_context.is_none()
    })
}

// ---------------------------------------------------------------------------
// Measurement
// ---------------------------------------------------------------------------

/// One divergence, for the failure message when a pin moves.
struct Divergence {
    id: String,
    outputs: Vec<String>,
}

/// A finding worth naming in the report: a validity failure, an accepted
/// prohibition, or a negative-guard violation.
struct Finding {
    id: String,
    what: String,
}

struct Measured {
    census: Census,
    divergences: Vec<Divergence>,
    findings: Vec<Finding>,
}

/// The built corpus, shared across both directions.
///
/// Enumerating it costs ~0.3 s, so the two direction tests each build their own
/// rather than sharing a `OnceLock` — which keeps them independent and lets
/// nextest run them in parallel.
fn built() -> SpecCorpus {
    corpus(&CorpusBounds::default())
}

/// Rows in an order that groups every row sharing a synthetic reference, so one
/// provider is built per reference rather than per row.
fn grouped(built: &SpecCorpus) -> Vec<&Row> {
    let mut rows: Vec<&Row> = built.rows.iter().collect();
    rows.sort_by(|a, b| {
        a.shape
            .label()
            .cmp(b.shape.label())
            .then_with(|| a.core.cmp(&b.core))
            .then_with(|| a.id.cmp(&b.id))
    });
    rows
}

/// The normalization configuration every figure in the census is measured under.
///
/// Extracted so the mode stamped into [`Census::measured_under`] is read off the
/// very config the normalizer is built from, rather than restated beside it
/// (#1629). A restated mode is a claim nothing checks, which is the defect the
/// stamp exists to close.
///
/// `NormalizeConfig::default()` substitutes `ErrorConfig::lenient()` — note that
/// `ErrorConfig::default()` and `ErrorMode`'s own `#[default]` are both *strict*,
/// so "the default" is not a description of this measurement.
fn measurement_config(direction: ShuffleDirection) -> NormalizeConfig {
    NormalizeConfig::default().with_direction(direction)
}

/// The contiguous runs of `rows` that share one synthetic reference.
///
/// [`grouped`] sorts by `(shape, core, id)`, so every row sharing a reference is
/// already adjacent — this only names the boundaries. It replaces the running
/// `frame` cache the serial loop carried: same "one provider per reference, not
/// per row" property, expressed as an explicit partition so [`measure`] can hand
/// each group to its own task.
fn reference_groups<'a>(rows: &'a [&'a Row]) -> Vec<&'a [&'a Row]> {
    let mut groups = Vec::new();
    let mut start = 0usize;
    while start < rows.len() {
        let mut end = start + 1;
        while end < rows.len()
            && rows[end].shape == rows[start].shape
            && rows[end].core == rows[start].core
        {
            end += 1;
        }
        groups.push(&rows[start..end]);
        start = end;
    }
    groups
}

/// Census one reference group. Builds its own frame and normalizer, reads
/// nothing outside its own slice, and so is safe to run alongside its siblings.
fn measure_group(group: &[&Row], direction: ShuffleDirection) -> Measured {
    let mut census = Census {
        measured_under: measurement_config(direction).error_config.mode,
        ..Census::default()
    };
    let mut divergences: Vec<Divergence> = Vec::new();
    let mut findings: Vec<Finding> = Vec::new();

    let active = group[0].frame();
    let normalizer =
        Normalizer::with_config(active.provider().clone(), measurement_config(direction));
    let (active, normalizer) = (&active, &normalizer);

    for row in group {
        let mut outputs: BTreeSet<String> = BTreeSet::new();
        let mut normalized_spellings = 0usize;
        for spelling in &row.spellings {
            census.outputs += 1;
            let parsed = match parse_hgvs(spelling) {
                Ok(parsed) => parsed,
                Err(_) => {
                    // A row whose own spelling does not parse is expected only
                    // where refusal IS the property.
                    if !matches!(row.kind, RowKind::Conflict | RowKind::Prohibited) {
                        census.declined += 1;
                    }
                    continue;
                }
            };
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                normalizer.normalize(&parsed)
            }));
            let output = match result {
                Ok(Ok(value)) => value.to_string(),
                _ => {
                    if !matches!(row.kind, RowKind::Conflict | RowKind::Prohibited) {
                        census.declined += 1;
                    }
                    continue;
                }
            };
            normalized_spellings += 1;

            // Refusal is the property for these two kinds, so reaching an output
            // at all is the finding.
            match row.kind {
                RowKind::Conflict => {
                    census.conflicts_accepted += 1;
                    findings.push(Finding {
                        id: row.id.clone(),
                        what: format!(
                            "conflicting allele accepted, geometry {}: {spelling} -> {output}",
                            row.geometry.label()
                        ),
                    });
                }
                RowKind::Prohibited => {
                    let (clause, strength) =
                        row.prohibition.expect("a prohibited row names its clause");
                    match strength {
                        Strength::Absolute => census.prohibited_absolute_accepted += 1,
                        Strength::Conditional => census.prohibited_conditional_accepted += 1,
                    }
                    findings.push(Finding {
                        id: row.id.clone(),
                        what: format!(
                            "{} prohibition accepted ({clause}): {spelling} -> {output}",
                            strength.label()
                        ),
                    });
                }
                RowKind::Family | RowKind::Single => {}
            }

            // --- rank 1: validity ---------------------------------------------
            let reparsed = parse_hgvs(&output);
            if reparsed.is_err() {
                census.unparseable_outputs += 1;
                findings.push(Finding {
                    id: row.id.clone(),
                    what: format!("output does not re-parse: {spelling} -> {output}"),
                });
            }
            if let Some(violation) = violated_prohibition(&output) {
                census.prohibition_violating_outputs += 1;
                findings.push(Finding {
                    id: row.id.clone(),
                    what: format!("output violates {violation}: {spelling} -> {output}"),
                });
            }

            // --- sequence preservation, and the overlap half of validity ------
            //
            // Three outcomes, kept apart deliberately. Collapsing
            // `Inexpressible` into `NoSequence` published a wrong headline once:
            // an intronic output has no SPDI triple, so 381 outputs were reported
            // as "two members claiming one territory" when most had simply left
            // the transcript — a different defect with different clauses.
            if let Some(expected) = row.denoted.as_deref() {
                match denotation_of(active.provider(), active.served(), &output) {
                    Denotation::NoSequence => {
                        census.outputs_denoting_no_sequence += 1;
                        findings.push(Finding {
                            id: row.id.clone(),
                            what: format!(
                                "output denotes no sequence (members claim one territory): \
                                 {spelling} -> {output}"
                            ),
                        });
                    }
                    // The oracle cannot express the output, and the input was
                    // exonic — so the description has left the transcript. Two
                    // very different things live under that one symptom, and
                    // #1704 is the reason they are now told apart:
                    //
                    // - a **bare** transcript accession naming an intronic
                    //   position is what `checklist.md:20` refuses, and is a
                    //   rank-1 output defect;
                    // - the same coordinate under the genomic wrapper that
                    //   clause asks for (`NC_…(NM_…):c.10+2del`) is CONFORMANT.
                    //   The oracle still cannot reach it — SPDI is positional and
                    //   has no offset notation, so `hgvs_to_spdi` declines the
                    //   offset whatever accession carries it — but that is a
                    //   limit of the instrument, which `Denotation` says of
                    //   itself in as many words.
                    //
                    // Keeping one counter over both would have made the rank-1
                    // figure unfixable: re-parenting closes the defect and moves
                    // `Denotation` not at all, so the census would have reported
                    // 371 before and after a change that made every one of those
                    // 371 outputs conformant.
                    Denotation::Inexpressible => {
                        // An output that does not re-parse cannot be shown to
                        // carry the wrapper, so it stays on the DEFECT side. A
                        // row must never fall into the conformant bucket for
                        // want of evidence.
                        if reparsed
                            .as_ref()
                            .ok()
                            .is_none_or(names_bare_transcript_intronic)
                        {
                            census.outputs_leaving_the_transcript += 1;
                            findings.push(Finding {
                                id: row.id.clone(),
                                what: format!(
                                    "output left the transcript (an intronic position the input \
                                     did not name, on a BARE transcript accession): \
                                     {spelling} -> {output}"
                                ),
                            });
                        } else {
                            census.outputs_intronic_under_a_genomic_wrapper += 1;
                        }
                    }
                    // Already counted as `unparseable_outputs` above.
                    Denotation::Unparseable => {}
                    Denotation::Sequence(applied) if applied != expected => {
                        census.sequence_changed += 1;
                        findings.push(Finding {
                            id: row.id.clone(),
                            what: format!("output denotes different bases: {spelling} -> {output}"),
                        });
                    }
                    Denotation::Sequence(_) => {}
                }
            }

            // --- idempotency --------------------------------------------------
            if let Ok(parsed_output) = reparsed {
                let again = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    normalizer.normalize(&parsed_output)
                }));
                let stable = matches!(&again, Ok(Ok(value)) if value.to_string() == output);
                if !stable {
                    census.non_idempotent_outputs += 1;
                    findings.push(Finding {
                        id: row.id.clone(),
                        what: format!(
                            "output is not a fixed point: {output} -> {}",
                            match &again {
                                Ok(Ok(value)) => value.to_string(),
                                Ok(Err(error)) => format!("<declined: {error}>"),
                                Err(_) => "<panicked>".to_string(),
                            }
                        ),
                    });
                }
            }

            // --- negative guards ---------------------------------------------
            // Only on the AUTHORED spelling: the spanning-`delins` respelling of
            // a two-member design is a single member already, so evaluating this
            // over every spelling would count the corpus's own candidate as a
            // merge. See `Row::authored_spelling`.
            if !row.negative_guards.is_empty() && spelling == row.authored_spelling() {
                if let Ok(parsed_output) = parse_hgvs(&output) {
                    // The guard is about to be decided on this output, so this is
                    // the denominator `guard_violations` is a numerator over.
                    census.guard_evaluations += 1;
                    if members_of(&parsed_output).len() < 2 {
                        census.guard_violations += 1;
                        findings.push(Finding {
                            id: row.id.clone(),
                            what: format!(
                                "merged a frameless separation of one into a single member, \
                                 which is rejected SVD-WG010: {spelling} -> {output}"
                            ),
                        });
                    }
                }
            }

            // --- the coding-axis merge instrument -----------------------------
            // Read on the AUTHORED spelling only, for the same reason the guard
            // above is: the corpus offers a spanning-`delins` respelling of its
            // own designs, and that respelling is a single member before
            // normalization ever runs.
            //
            // Not a verdict — see `Census::coding_axis_separation_two_or_more_merges`.
            // The denominator is counted here rather than off the corpus so
            // `0 of 0` cannot pass as a result even if every eligible row
            // stopped normalizing.
            if row.coding_axis_separation_two_or_more && spelling == row.authored_spelling() {
                census.coding_axis_separation_two_or_more_rows += 1;
                if let Ok(parsed_output) = parse_hgvs(&output) {
                    if coding_axis_merge_observed(row, &parsed_output) {
                        census.coding_axis_separation_two_or_more_merges += 1;
                        findings.push(Finding {
                            id: row.id.clone(),
                            what: format!(
                                "coding-axis allele merged across {} unchanged nucleotide(s), \
                                 {} authored members -> {} (`general.md:34` / \
                                 `DNA/delins.md:17` speak to this population, and the three \
                                 ruling records named on \
                                 `spec_corpus::is_coding_axis_separation_two_or_more_shape` \
                                 govern parts of it — counted, not adjudicated): \
                                 {spelling} -> {output}",
                                row.separation,
                                row.members,
                                members_of(&parsed_output).len()
                            ),
                        });
                    }
                }
            }

            outputs.insert(output);
        }

        // --- rank 2: confluence, over families only --------------------------
        if row.kind == RowKind::Family {
            if normalized_spellings < 2 {
                census.underdetermined += 1;
            } else {
                match outputs.len() {
                    1 => census.converged += 1,
                    2 => census.split_two += 1,
                    3 => census.split_three += 1,
                    _ => census.split_more += 1,
                }
                if outputs.len() > 1 && divergences.len() < MAX_DIVERGENCES {
                    divergences.push(Divergence {
                        id: row.id.clone(),
                        outputs: outputs.into_iter().collect(),
                    });
                }
            }
        }
    }

    Measured {
        census,
        divergences,
        findings,
    }
}

/// How many divergent families the failure message names. Applied per group and
/// again to the concatenation, which is what keeps the parallel result
/// byte-identical to a serial one — see [`measure`].
const MAX_DIVERGENCES: usize = 12;

impl Census {
    /// Fold another census in.
    ///
    /// Every field below is a count over a partition of the rows, so summing is
    /// exact rather than approximate — that is what makes [`measure`]'s
    /// per-group split safe.
    ///
    /// **The destructuring is the point, not style.** It carries no `..`, so
    /// adding a field to [`Census`] makes this function fail to COMPILE until
    /// the field is folded. An earlier version of this listed the fields by hand
    /// and #1710 added two underneath it; they were silently never summed, and
    /// the census reported `coding_axis_separation_two_or_more_merges` as
    /// **0 of 0** — a denominator of zero, which is the flattering direction
    /// this file's own VACUOUS guard exists to catch. It did catch it. This
    /// makes the next one a build error instead.
    ///
    /// [`Census::measured_under`] is bound and deliberately NOT summed: it is a
    /// mode stamp rather than a counter, and every group is built from the same
    /// [`measurement_config`], so the accumulator's own stamp already carries
    /// it. Binding it keeps it inside the exhaustiveness check.
    fn absorb(&mut self, other: &Census) {
        let Census {
            measured_under: _,
            outputs,
            declined,
            unparseable_outputs,
            outputs_denoting_no_sequence,
            outputs_leaving_the_transcript,
            outputs_intronic_under_a_genomic_wrapper,
            prohibition_violating_outputs,
            converged,
            split_two,
            split_three,
            split_more,
            underdetermined,
            non_idempotent_outputs,
            sequence_changed,
            conflicts_accepted,
            prohibited_absolute_accepted,
            prohibited_conditional_accepted,
            guard_violations,
            guard_evaluations,
            coding_axis_separation_two_or_more_rows,
            coding_axis_separation_two_or_more_merges,
        } = other;
        self.outputs += outputs;
        self.declined += declined;
        self.unparseable_outputs += unparseable_outputs;
        self.outputs_denoting_no_sequence += outputs_denoting_no_sequence;
        self.outputs_leaving_the_transcript += outputs_leaving_the_transcript;
        self.outputs_intronic_under_a_genomic_wrapper += outputs_intronic_under_a_genomic_wrapper;
        self.prohibition_violating_outputs += prohibition_violating_outputs;
        self.converged += converged;
        self.split_two += split_two;
        self.split_three += split_three;
        self.split_more += split_more;
        self.underdetermined += underdetermined;
        self.non_idempotent_outputs += non_idempotent_outputs;
        self.sequence_changed += sequence_changed;
        self.conflicts_accepted += conflicts_accepted;
        self.prohibited_absolute_accepted += prohibited_absolute_accepted;
        self.prohibited_conditional_accepted += prohibited_conditional_accepted;
        self.guard_violations += guard_violations;
        self.guard_evaluations += guard_evaluations;
        self.coding_axis_separation_two_or_more_rows += coding_axis_separation_two_or_more_rows;
        self.coding_axis_separation_two_or_more_merges += coding_axis_separation_two_or_more_merges;
    }
}

/// Census one direction, one reference group per rayon task.
///
/// **The result is byte-identical to the serial walk, not merely equivalent**,
/// which is the only basis on which a pinned census may be parallelized at all:
///
/// * every counter is a count over a partition of the rows, so [`Census::absorb`]
///   over the groups is exact and order-free;
/// * `par_iter().collect()` preserves *input* order, so the per-group results
///   fold in corpus order however they finish — which matters for `findings`,
///   an uncapped list whose order the failure message reproduces;
/// * `divergences` is "the first [`MAX_DIVERGENCES`] divergent families in
///   corpus order". Each group keeps its own first `MAX_DIVERGENCES`, and
///   concatenating in corpus order and truncating again yields exactly that set
///   — the serial loop's global cap can only have kept the same ones.
///
/// Nothing is shared across tasks: [`measure_group`] builds its own `Frame` and
/// `Normalizer` (the serial loop already rebuilt both per reference — that is
/// what the `frame` cache was), and the only borrow that crosses is the
/// immutable corpus slice.
fn measure(direction: ShuffleDirection) -> Measured {
    let built = built();
    let rows = grouped(&built);
    let groups = reference_groups(&rows);

    let per_group: Vec<Measured> = groups
        .par_iter()
        .map(|group| measure_group(group, direction))
        .collect();

    let mut measured = Measured {
        census: Census {
            measured_under: measurement_config(direction).error_config.mode,
            ..Census::default()
        },
        divergences: Vec::new(),
        findings: Vec::new(),
    };
    for group in per_group {
        measured.census.absorb(&group.census);
        measured.divergences.extend(group.divergences);
        measured.findings.extend(group.findings);
    }
    measured.divergences.truncate(MAX_DIVERGENCES);
    measured
}

fn report(label: &str, measured: &Measured) -> String {
    let census = &measured.census;
    let mut out = format!(
        "spec conformance axis ({label})\n  \
         NORMALIZED UNDER: error mode `{}` — every figure downstream of normalization is \
         a figure for THAT mode; the parse half applies no mode at all (#1632)\n  \
         VALIDITY (rank 1): {} outputs, {} declined, {} unparseable, {} denoting no sequence, \
         {} leaving the transcript on a bare accession, {} intronic under a genomic wrapper \
         (conformant), {} violating an absolute prohibition\n  \
         CONFLUENCE (rank 2): converged {}, split 2 {}, split 3 {}, split 4+ {}, \
         underdetermined {}\n  \
         IDEMPOTENCY: {} outputs are not their own fixed point\n  \
         SEQUENCE PRESERVATION: {} outputs denote different bases\n  \
         REFUSAL: {} conflicting alleles accepted, {} absolute prohibitions accepted, \
         {} conditional prohibitions accepted\n  \
         NEGATIVE GUARDS: {} outputs implement rejected SVD-WG010\n  \
         CODING-AXIS MERGE (an instrument, not a verdict): {} of {} coding-axis alleles \
         separated by two or more unchanged nucleotides normalized to fewer members than \
         were authored (a full or a partial merge)\n",
        census.measured_under,
        census.outputs,
        census.declined,
        census.unparseable_outputs,
        census.outputs_denoting_no_sequence,
        census.outputs_leaving_the_transcript,
        census.outputs_intronic_under_a_genomic_wrapper,
        census.prohibition_violating_outputs,
        census.converged,
        census.split_two,
        census.split_three,
        census.split_more,
        census.underdetermined,
        census.non_idempotent_outputs,
        census.sequence_changed,
        census.conflicts_accepted,
        census.prohibited_absolute_accepted,
        census.prohibited_conditional_accepted,
        census.guard_violations,
        census.coding_axis_separation_two_or_more_merges,
        census.coding_axis_separation_two_or_more_rows,
    );
    for divergence in &measured.divergences {
        out.push_str(&format!(
            "  DIVERGES {} -> {:?}\n",
            divergence.id, divergence.outputs
        ));
    }
    // Findings are grouped, because an accepted prohibition fires once per row
    // and there are hundreds of rows per clause.
    let mut grouped: BTreeMap<String, (usize, String)> = BTreeMap::new();
    for finding in &measured.findings {
        let key = finding
            .what
            .split_once(':')
            .map_or_else(|| finding.what.clone(), |(head, _)| head.to_string());
        let entry = grouped
            .entry(key)
            .or_insert_with(|| (0, format!("{} | {}", finding.id, finding.what)));
        entry.0 += 1;
    }
    for (kind, (count, example)) in grouped {
        out.push_str(&format!("  FINDING x{count} [{kind}] e.g. {example}\n"));
    }
    out
}

/// Assert one direction's census against its pin, printing the measured numbers
/// either way so a moved pin can be re-blessed from the test output.
fn assert_census(direction: ShuffleDirection, label: &str, pinned: &Census) {
    let measured = measure(direction);
    println!("{}", report(label, &measured));
    let census = &measured.census;

    // The honest-zero discipline: a property whose denominator is zero is
    // VACUOUS, not passing.
    //
    // Asserted FIRST, ahead of the residual set below. On an empty corpus every
    // downstream property is vacuously satisfied except the residual one, which
    // fails with "the negative-guard residual moved" — a message that sends the
    // reader after a representation change when what actually happened is that
    // nothing ran. The vacuity check names that, so it must get there first.
    assert!(
        census.outputs > 0,
        "{label}: VACUOUS — the corpus produced no outputs at all"
    );
    assert!(
        census.coding_axis_separation_two_or_more_rows > 0,
        "{label}: VACUOUS — no coding-axis allele separated by two or more unchanged \
         nucleotides was evaluated, so `coding_axis_separation_two_or_more_merges` measures \
         nothing and its zero is a claim about the corpus rather than about the \
         implementation.\n{}",
        report(label, &measured)
    );

    // The negative-guard residual, as a SET and not only as a count.
    //
    // The FERRO_IS_WRONG shape: this fails when a listed row starts conforming,
    // not only when an unlisted one stops. A residual that can only be asserted
    // downwards becomes a permanent exemption, because the day it is fixed
    // nothing says so and the list stays.
    //
    // Asserted here rather than in a test of its own because `measure` is the
    // expensive call in this file — a standalone test would run it a third and
    // fourth time for no new information, and it already takes ~80 s per
    // direction.
    //
    // Deduplicated, because this is a claim about *membership* and
    // `RESIDUAL_GUARD_VIOLATIONS` is a list of distinct row ids. The guard at the
    // `negative_guards` block above fires once per authored spelling, so a row
    // whose `spellings` carried its authored spelling twice would contribute the
    // same id twice and fail this on arity rather than on membership — which
    // reads as a moved residual and is not one. The multiplicity is not
    // discarded, it is asserted separately just below.
    let mut violating: Vec<&str> = measured
        .findings
        .iter()
        .filter(|finding| finding.what.starts_with("merged a frameless separation"))
        .map(|finding| finding.id.as_str())
        .collect();
    violating.sort_unstable();
    let fired = violating.len();
    violating.dedup();
    let mut residual: Vec<&str> = RESIDUAL_GUARD_VIOLATIONS.to_vec();
    residual.sort_unstable();
    assert_eq!(
        violating, residual,
        "{label}: the negative-guard residual moved. A row that STOPPED \
         violating is a fix — delete it from RESIDUAL_GUARD_VIOLATIONS and lower \
         `guard_violations` in the same commit, so the residual cannot sit here \
         as a permanent exemption. A row that STARTED violating is a \
         representation regression on `general.md:34`."
    );

    // Membership and count must agree, or the pinned `guard_violations` below is
    // counting something the set assertion above cannot see. They can only
    // diverge by a row firing twice, which means a duplicated authored spelling
    // in the corpus — a corpus defect, and one the dedup above would otherwise
    // absorb silently.
    assert_eq!(
        fired,
        violating.len(),
        "{label}: a guarded row fired more than once, so `guard_violations` \
         ({}) counts spellings while RESIDUAL_GUARD_VIOLATIONS counts rows. \
         Look for a duplicated authored spelling in the corpus, not for a \
         representation change.",
        census.guard_violations
    );

    // The zero above is only a finding if the guard ran. Asserted at the
    // denominator, not inferred from `unparseable_outputs` elsewhere.
    assert!(
        census.guard_evaluations > 0,
        "{label}: VACUOUS — the rejected-SVD-WG010 guard was never evaluated, so \
         `guard_violations: {}` is a structural silence rather than a measurement",
        census.guard_violations
    );

    // The rank-1 counters are pinned rather than asserted at zero, because they
    // are not zero on `main` — see the module docs' table. What is asserted here
    // is that they never RISE, checked ahead of the full census so the message
    // names the property rather than dumping a struct diff.
    for (measured_value, pinned_value, what) in [
        (
            census.unparseable_outputs,
            pinned.unparseable_outputs,
            "outputs that do not re-parse",
        ),
        (
            census.outputs_denoting_no_sequence,
            pinned.outputs_denoting_no_sequence,
            "outputs that denote no sequence (two members claiming one territory)",
        ),
        (
            census.outputs_leaving_the_transcript,
            pinned.outputs_leaving_the_transcript,
            "outputs that left the transcript (an intronic position the input did not name, \
             on a bare transcript accession — checklist.md:20)",
        ),
        (
            census.sequence_changed,
            pinned.sequence_changed,
            "outputs that denote different bases from their input",
        ),
        (
            census.non_idempotent_outputs,
            pinned.non_idempotent_outputs,
            "outputs that are not their own fixed point",
        ),
        (
            census.conflicts_accepted,
            pinned.conflicts_accepted,
            "conflicting alleles accepted instead of refused",
        ),
        (
            census.prohibited_absolute_accepted,
            pinned.prohibited_absolute_accepted,
            "shapes the spec calls \"not allowed\" that were accepted",
        ),
        // The two the loop was missing. Both were still caught by the full
        // census assertion below, but only as a struct diff — which is the one
        // thing this loop exists to avoid, since a reader then has to work out
        // which of twenty fields moved and whether it moved the wrong way.
        (
            census.prohibition_violating_outputs,
            pinned.prohibition_violating_outputs,
            "outputs that violate an absolute textual prohibition",
        ),
        (
            census.prohibited_conditional_accepted,
            pinned.prohibited_conditional_accepted,
            "shapes a conditional clause prohibits that were accepted",
        ),
        (
            census.guard_violations,
            pinned.guard_violations,
            "outputs implementing REJECTED SVD-WG010",
        ),
    ] {
        assert!(
            measured_value <= pinned_value,
            "{label}: {what} rose from {pinned_value} to {measured_value}. This is a rank-1 \
             conformance regression, not a representation choice.\n{}",
            report(label, &measured)
        );
    }

    // The instrument is checked apart from the loop above, because that loop's
    // message calls a rise a rank-1 conformance regression and this counter
    // makes no such claim. A rise here is a population getting bigger — a
    // finding to explain and, if it is the right behaviour, to re-pin.
    assert!(
        census.coding_axis_separation_two_or_more_merges
            <= pinned.coding_axis_separation_two_or_more_merges,
        "{label}: coding-axis alleles merged across two or more unchanged nucleotides rose \
         from {} to {} (of {} evaluated). This counter is an INSTRUMENT, not a ruling: \
         `general.md:34` / `DNA/delins.md:17` say such members are described individually, and \
         `general.md:35`'s one-amino-acid exception cannot reach a separation of two by its own \
         stated conjunct — but the ruling ledger governs part of this population and a rise is \
         not automatically a defect. Read the three records named on \
         `spec_corpus::is_coding_axis_separation_two_or_more_shape` BEFORE deciding what this \
         means: within the payload-coincidence shape on this axis a rise is a decided ruling \
         arriving, and outside it `general.md:34` still governs. So this is a population that \
         grew, and it needs explaining and re-pinning rather than assuming either verdict.\n{}",
        pinned.coding_axis_separation_two_or_more_merges,
        census.coding_axis_separation_two_or_more_merges,
        census.coding_axis_separation_two_or_more_rows,
        report(label, &measured)
    );

    assert_eq!(
        census,
        pinned,
        "{label}: the conformance census moved. Every failure figure must only ever go DOWN and \
         `converged` only ever UP; if this change lowers one, re-bless the pin in the same commit \
         and say so in the PR, and promote any newly-failing row to a named test in \
         `spec_corpus_regressions.rs`. Measured:\n{}",
        report(label, &measured)
    );
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// The corpus's own shape, before any property is read out of it.
///
/// #1460's lesson mechanized: a generator change that shrank the corpus would
/// otherwise improve every rate below while measuring less.
#[test]
fn the_corpus_has_the_shape_its_censuses_are_measured_over() {
    let built = built();
    let measured = CorpusShape::of(&built);
    println!("corpus shape: {measured:#?}");
    assert_eq!(
        measured, SHAPE,
        "the corpus shape moved; every census below is measured over it, so re-bless both \
         together and say what changed"
    );

    // Each property must have a non-zero denominator, or its census is a claim
    // about the corpus rather than about the implementation.
    assert!(measured.family_rows > 0, "VACUOUS: no family rows");
    assert!(measured.conflict_rows > 0, "VACUOUS: no conflict rows");
    assert!(
        measured.prohibited_rows > 0,
        "VACUOUS: no prohibited-shape rows"
    );
    assert!(
        measured.guarded_rows > 0,
        "VACUOUS: no rows carry a negative guard, so the rejected-SVD-WG010 check measures nothing"
    );
    assert!(
        measured.single_rows > 0,
        "VACUOUS: no single-spelling rows, so intronic and trans/unknown-phase idempotency \
         measures nothing"
    );
    assert!(
        measured.coding_axis_separation_two_or_more_rows > 0,
        "VACUOUS: no coding-axis rows separated by two or more unchanged nucleotides, so the \
         coding-axis merge counter measures nothing"
    );

    // Multi-member is the priority axis. Real corpora are 592 of 9,949,738.
    let share = measured.multi_member_rows as f64 / measured.rows as f64;
    assert!(
        share > 0.9,
        "multi-member share is {share:.4}; this corpus exists to invert the natural 0.00006"
    );
}

/// A **positive control** for the coding-axis merge counter: it must be able to
/// observe a merge.
///
/// # Why a zero pin needs one, and why nothing else supplies it
///
/// `coding_axis_separation_two_or_more_merges` is pinned at **0**, so every
/// other test that mentions it is satisfied by a numerator that can never
/// increment. Measured as a mutation matrix over the four tests that claim to
/// cover the counter — the two censuses, the corpus-shape test and
/// `spec_corpus`'s population test — hard-wiring
/// [`coding_axis_merge_observed`] to `false` leaves **all four green**. "Pinned
/// at 0" and "wired to 0" were the same observation, which is exactly what this
/// module's own thesis says an instrument must not be: *a counter nobody has
/// seen move is not an instrument*.
///
/// The other five mutations in that matrix are caught: predicate always-true,
/// predicate always-false (as VACUOUS), the separation floor weakened to one,
/// dropping the irreducibility conjunct, and the numerator forced always-true.
/// Only the always-**false** direction was blind, because it is the direction a
/// zero cannot distinguish.
///
/// # What it asserts
///
/// The predicate is exercised directly against real corpus rows and parsed
/// outputs, so it covers `members_of`, the parse path and the member-count
/// comparison — everything between normalization and the counter:
///
/// - an authored multi-member spelling is **not** a merge (no false positive);
/// - a single-member output **is** (the plain merge);
/// - a 3-or-more-member row collapsed to two members **is** — the partial merge
///   the old `< 2` test could not see on 313 of the 997 flagged rows.
#[test]
fn the_coding_axis_merge_counter_can_observe_a_merge() {
    let built = built();
    let pair = built
        .rows
        .iter()
        .find(|row| row.coding_axis_separation_two_or_more && row.members == 2)
        .expect(
            "VACUOUS: the corpus builds no two-member coding-axis row separated by two or more \
             unchanged nucleotides, so this control controls nothing",
        );

    let authored = parse_hgvs(pair.authored_spelling()).expect("the authored spelling parses");
    assert!(
        !coding_axis_merge_observed(pair, &authored),
        "{}: its own authored spelling ({}) counted as a merge, so the counter has a false \
         positive and its zero would be unreadable in the other direction",
        pair.id,
        pair.authored_spelling()
    );

    // The merged form is the row's OWN spanning-`delins` respelling, not a
    // hand-written literal: the corpus offers one for every family it builds,
    // and it is a single member before normalization ever runs (which is why
    // `measure` reads this counter on the authored spelling only). Deriving it
    // from the row keeps the control anchored to the population it controls.
    let (merged_spelling, merged) = pair
        .spellings
        .iter()
        .filter_map(|spelling| Some((spelling, parse_hgvs(spelling).ok()?)))
        .find(|(_, parsed)| members_of(parsed).len() == 1)
        .expect("VACUOUS: the row offers no single-member respelling to merge INTO");
    assert!(
        coding_axis_merge_observed(pair, &merged),
        "{}: its own single-member respelling ({merged_spelling}) was NOT counted as a merge. \
         The counter is blind, and its pinned 0 is a property of this predicate rather than of \
         the normalizer",
        pair.id
    );

    // The partial half. `pair`'s two-member authored spelling stands in for a
    // partial merge of `wide`, which authored three or more — 2 < wide.members
    // is exactly the property, and both sides are corpus rows rather than
    // literals. Re-spelling `wide` itself with a member dropped would mean
    // reassembling an allele out of `Display` output, which is string surgery
    // with its own failure mode and would not test anything further: the
    // predicate reads member counts and nothing else.
    let wide = built
        .rows
        .iter()
        .find(|row| row.coding_axis_separation_two_or_more && row.members >= 3)
        .expect(
            "VACUOUS: the corpus builds no 3-or-more-member coding-axis row, so the partial-merge \
             half of this control controls nothing",
        );
    assert!(
        coding_axis_merge_observed(wide, &authored),
        "{}: {} authored members collapsing to {} was NOT counted. That is the 313 of 997 rows \
         (31.4%) `coding_axis_merge_observed` exists to stop losing",
        wide.id,
        wide.members,
        members_of(&authored).len()
    );
}

/// Every family really is a family: two or more distinct spellings, each of
/// which denotes the row's sequence when applied independently of the
/// normalizer.
///
/// Without this the confluence figure could be a statement about a generator bug
/// — a "family" whose members denote different variants would diverge for a
/// reason that has nothing to do with the normalizer.
#[test]
fn every_family_denotes_one_sequence_independently_of_the_normalizer() {
    let built = built();
    let mut checked = 0usize;
    // Every family, not a sample. An earlier revision sampled, which is the wrong
    // trade here: the strata are enumerated in order, so any prefix or stride
    // risks skipping `scale` and `repeats` entirely — and those are the two whose
    // ground truth is hardest to get right. Measured at well under a second,
    // because the applier is a splice and a byte compare.
    for row in grouped(&built) {
        if row.kind != RowKind::Family {
            continue;
        }
        let frame = row.frame();
        let expected = row.denoted.as_deref().expect("a family denotes a sequence");
        assert!(row.spellings.len() >= 2, "{} is a singleton", row.id);
        assert_eq!(
            row.spellings.iter().collect::<BTreeSet<_>>().len(),
            row.spellings.len(),
            "{} repeats a spelling",
            row.id
        );
        for spelling in &row.spellings {
            assert_eq!(
                denoted_by(frame.provider(), frame.served(), spelling).as_deref(),
                Some(expected),
                "{}: {spelling} does not denote the row's sequence",
                row.id
            );
        }
        checked += 1;
    }
    assert!(checked > 0, "VACUOUS: no families were checked");
}

#[test]
fn three_prime_conformance_census() {
    assert_census(ShuffleDirection::ThreePrime, "3prime", &THREE_PRIME);
}

#[test]
fn five_prime_conformance_census() {
    assert_census(ShuffleDirection::FivePrime, "5prime", &FIVE_PRIME);
}
