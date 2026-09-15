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
//! | `prohibited_absolute_accepted` | ~~32~~ ~~24~~ ~~8~~ **0** | ~~32~~ ~~24~~ ~~8~~ **0** | a shape the spec calls "not allowed" was accepted. Re-blessed by #1628 (genomic offsets) and then to zero by #1789 (`ins6`) |
//! | `prohibition_violating_outputs` | ~~32~~ ~~8~~ **0** | ~~32~~ ~~8~~ **0** | and then EMITTED, so the prohibition is not enforced on output either. Re-blessed by #1627 and then to zero by #1628: see the fourth section below |
//!
//! `guard_violations` was **0 of 210 guarded rows** through #1899, when the
//! payload-coincidence carve-out was scoped to `c.` alone and every guarded row
//! (frameless `g.`/`n.` axes) sat outside it. **#2155 supersedes that scope to
//! all DNA axes** (`c./g./m./n.`, per `AxisFrame::is_dna`), which puts the
//! guard's own domain inside the widened carve-out's reach for the first time —
//! see the dedicated section below, "`guard_violations` is reframed", for the
//! measured count and why a non-zero pin here is now a disclosed deviation
//! rather than a fault. `r.` is unaffected and stays out, both structurally (the
//! corpus builds no `r.` row at all) and by a targeted hermetic test. The
//! denominator is still asserted non-zero, because `0 of 0` is what a rebuilt
//! #1456 looks like.
//!
//! The 3'/5' asymmetry is itself a finding: the junction-crossing class is
//! **371 to 0**, and 3' is the default direction. #1704 makes those 371 outputs
//! conformant without moving a single coordinate, so the asymmetry survives the
//! fix intact — it is now carried by `outputs_intronic_under_a_genomic_wrapper`,
//! and it is still a claim about the code rather than about the corpus.
//!
//! # RE-BLESSED (1 of 6) — #1599, the amino-acid precondition
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
//! # RE-BLESSED (2 of 6) — #1536, the cross-axis re-typing carve-out
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
//! # RE-BLESSED (3 of 6) — #1649, the two-deletion alignment
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
//! lists were dumped on `origin/main` and on this branch (at the time, by raising
//! the `divergences` cap in `measure` locally and then restoring it; since #2063
//! `measure` returns every divergent family and only `report` caps, so that step
//! is no longer needed) and the sets diffed:
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
//! # RE-BLESSED (4 of 6) — #1627, the alignment-only symbol
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
//! | 3' and 5' `prohibition_violating_outputs` | 32 | **8** | improves by 24 (#1627); #1628 then takes it to **0** |
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
//! violating outputs — **which #1628 has since closed, taking that counter to
//! 0**; the paragraph below is the state as of #1627 and is kept because the
//! reasoning about the two counters is what makes the two re-blessings
//! separable.
//!
//! **`prohibited_absolute_accepted` did NOT move for #1627**, and its staying
//! at 32 was the useful half of that measurement. The corpus grades an
//! `X` row `Strength::Conditional` — the footnote states a scope rather than
//! using prohibitive words — while this axis counts an `X` *output* as an
//! absolute violation. The decided ruling records that the grading is moot
//! either way, because `general.md:48` admits only IUPAC-IUBMB symbols. So the
//! two counters #1627 moved are the two the clause was actually in, and the
//! then-surviving 32 absolute acceptances were `checklist.md:33`'s `ins6` (24)
//! plus `checklist.md:16`/`:45`'s genomic offset and hyphen range (4 + 4).
//! **#1628 took the second group** and **#1789 takes the first**, so the counter
//! now reads 0 and no absolute-prohibition residue remains.
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
//! # RE-BLESSED (5 of 6) — #1716, the codon-frame merge's own span
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
//! both happen. The full divergence row-id lists were dumped either side (again
//! by raising the then-existing `divergences` cap locally — see above, it is gone)
//! and diffed:
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
//! # RE-BLESSED (6 of 6) — #1621, the exon-junction clamp from the far side
//!
//! `general.md:44` exempts deletions and duplications around exon/exon
//! junctions from the 3' rule on `c.`/`n.`, and ferro applied that exemption in
//! one direction only: a description approaching a junction stopped at it, one
//! already spelled past it stayed put. The operator ruling
//! `exon-junction-dup-converge-from-the-far-side` reads the exemption as a clamp
//! binding in **both** directions, so a description pinned at its exon's first
//! base whose equivalence run continues back over the junction is pulled onto
//! the near side.
//!
//! **Three figures move. `converged` is unchanged, no failure counter rises,
//! and both rank-2 moves are arity FALLS:**
//!
//! | figure | was (this branch's base, `96305ded`) | now | direction |
//! |---|---|---|---|
//! | 3' `outputs_intronic_under_a_genomic_wrapper` | 389 | **391** | +2 — conformant, not a defect counter |
//! | 3' `converged` | 11,016 | **11,016** | unchanged |
//! | 3' `split_two` | 673 | **674** | +1, from `split_three` |
//! | 3' `split_three` | 164 | **164** | unchanged — one in, one out |
//! | 3' `split_more` | 29 | **28** | −1, to `split_three` |
//!
//! **Every absolute above was restated on the #1835 rebase; every DELTA is the
//! one first measured.** The table read `371 -> 373`, `9,402` unchanged,
//! `2,223 -> 2,224`, `226` unchanged and `31 -> 30` when it was written, against
//! a base predating the partition flip (#1835), which moved 1,614 families into
//! `converged` (`split_two` −1,550, `split_three` −62, `split_more` −2) — and
//! #1840 and #1704 moved the intronic counter. So `converged` was **9,402** in
//! the prose and **11,016** in the constant three hundred lines below, which is
//! the shape this repository's `CLAUDE.md` names: a derived figure quoted
//! against a base nobody restates. Quote the deltas, which survived the rebase
//! unchanged, and re-read the absolutes off [`THREE_PRIME`] rather than from
//! here.
//!
//! **The 5' census does not move at all**, which is a check on the change
//! rather than a coincidence: the clamp is gated on
//! `ShuffleDirection::ThreePrime`, because the 3' rule is what `general.md:44`
//! carves an exception out of and a 5'/VCF shuffle has no such rule to clamp.
//!
//! **Exactly four family instances change, all in the `junction-1-del-del`
//! stratum, and they are named** — the net counters cannot tell an arity rise
//! from a fall. The full divergence row-id lists were dumped either side (the
//! `MAX_DIVERGENCES` cap raised locally, then restored) and diffed:
//!
//! - **0 families lose convergence, 0 gain it, and 0 become divergent that were
//!   not** — the divergent *set* is identical either side, 2,480 families.
//! - `s01-c3p-junction-1-del-del-p1-sep0` 4 -> 3 and `…-sep1` 3 -> 2: in both,
//!   the output that disappears is `NM_TEST.1:c.21_22del`, the far-side spelling,
//!   which now converges onto `c.19_20del` — the arity fall IS the ruling.
//! - `s01-c3m-junction-1-del-del-p1-sep{0,1}` keep their arity and re-spell that
//!   same output as `NC_SYNTH.1(NM_TEST.1):c.18_20+2A[3]`. Those are the +2 on
//!   the intronic counter: on the minus strand the pulled-back description rests
//!   on the junction's near side, where #670's junction-crossing continuation
//!   carries it into the intron. Both clauses apply and they compose —
//!   `general.md:44` forbids the exon/EXON shift, the exception 3' rule permits
//!   the exon/INTRON one — so the row is conformant, not a new violation.
//!
//! **Nothing is promoted to `spec_corpus_regressions.rs`**, because nothing
//! newly fails: every moved figure improves or is conformant. The ruling's own
//! row is pinned in `spec_worked_examples.rs`
//! (`the_decided_target_is_convergence_on_the_near_side`, no longer `#[ignore]`d)
//! and the rest of its scope — deletions, the `n.` axis, and both directions of
//! the clamp — in `issue_1621_exon_junction_far_side.rs`.
//!
//! **Every rank-1 counter is unchanged** — `outputs`, `declined`,
//! `unparseable_outputs`, `outputs_denoting_no_sequence`,
//! `outputs_leaving_the_transcript`, `prohibition_violating_outputs` — as are
//! `non_idempotent_outputs` (4), `sequence_changed` (4), all three refusal
//! counters and `guard_violations`.
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
//! # RE-BLESSED (5 of 5) — #1789, the size-count insertion
//!
//! `checklist.md:33` says a description like `c.5439_5430ins6` "is not allowed,
//! the inserted sequence (for `ins6`, e.g., `TGCCAT`) should be specified", and
//! `DNA/insertion.md:22` / `RNA/insertion.md:20` enumerate the payload forms
//! without a bare count among them. Ferro accepted `ins6` in every mode and
//! re-emitted it. It now refuses: strict at parse, lenient and silent at
//! normalize, per the decided `rulings[absolute-prohibition-enforcement-stage]`
//! — which names `checklist.md:33` in its own clause list.
//!
//! **One figure moved, DOWN, identically in both directions:**
//!
//! | figure | was | now | direction |
//! |---|---|---|---|
//! | 3' and 5' `prohibited_absolute_accepted` | 32 | **8** | improves by 24 |
//!
//! **`prohibition_violating_outputs` deliberately does NOT move**, and its
//! staying at 8 is this measurement's useful half. That counter was already 8
//! after #1627, and its per-clause decomposition — pinned by
//! `every_prohibition_violating_output_is_a_re_emitted_prohibited_input` — is
//! entirely `checklist.md:16`. The `ins6` rows were never in it: this axis's
//! output-side predicate is a text check for a `g.` offset and a hyphen range,
//! and `ins6` is neither. So the input counter is the only one the clause was
//! in, which is the mirror image of #1627's split and is why both are worth
//! stating.
//!
//! **Every rank-1, rank-2, idempotency and sequence-preservation figure is
//! unchanged**, in both directions — `outputs`, `declined`, `converged`,
//! `split_*`, `non_idempotent_outputs`, `sequence_changed`,
//! `outputs_denoting_no_sequence`, `conflicts_accepted` and `guard_violations`
//! all hold. That is what refusing un-denotable inputs should look like: the 24
//! refused rows contributed no confluence family and no sequence, so nothing
//! they leave behind can move.
//!
//! **Nothing to promote to `spec_corpus_regressions.rs`** — no row newly fails.
//! Two rows there *stopped* failing and are re-blessed in place; the guard that
//! replaces them is `issue_1789_insertion_size_count.rs`, which also covers the
//! `ins(6)` spelling and the `r.` axis, neither of which the corpus generates.
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
//! | `live` — the rule shipped to v0.14.0 | 0 of 997 | 0 of 997 |
//! | `shadow` | 0 | 0 |
//! | `canonical` | 1 | 1 |
//! | unset = `canonical-coalesced` — **the shipped rule since #1835** | **3** of 997 | **3** of 997 |
//!
//! The 3 are one row at a separation of 3 and two at a separation of 2, all of
//! them two-member designs collapsing to one member — every one a row
//! `guard_violations` cannot see.
//!
//! # The 3 are named and adjudicated, because a count is not a verdict (#1835)
//!
//! The counter's own text says a rise "needs explaining and re-pinning rather
//! than assuming either verdict", and that within the payload-coincidence shape a
//! rise is a decided ruling arriving while outside it `general.md:34` still
//! governs. So the rows were identified rather than counted, by instrumenting the
//! census to print each merging row and re-running under three arms:
//!
//! ```text
//! s01-c1-pair-inv-del-p4-sep3     c.[9_12inv;16_19del]  -> c.9_19delinsGTCAGTA
//! s01-c3p-pair-del-sub-p4-sep2    c.[24_27del;30C>G]    -> c.24_30delinsTTG
//! s01-c3m-pair-del-sub-p4-sep2    c.[24_27del;30C>G]    -> c.24_30delinsTTG
//! ```
//!
//! **The `sep3` row is NOT the coalesce pass**, and the table above already says
//! so — it is the `canonical` arm's single residual and it produces the identical
//! string on both candidate arms. What that means is that the canonical
//! partitioner returns this block as **one piece**: the minimal alignment of the
//! resulting sequence has no unchanged interior column at all. `general.md:34`
//! speaks about "two variants separated by one or more nucleotides", and
//! `separation-is-a-property-of-the-spelling-not-of-the-variant` (decided) reads
//! that separation off the partition re-derived from the resulting sequence, not
//! off the authored spelling. There is no separation in the derived partition, so
//! the clause has nothing to govern here. The row is counted because the counter
//! keys on the *authored* member count, which is what makes it an instrument.
//!
//! **The two `sep2` rows ARE the coalesce pass, and they are the decided ruling
//! arriving.** This is the case where reading only the authored spelling gets the
//! answer backwards, so it is worth stating fully. Authored, the members are
//! `24_27del` (a pure deletion) and `30C>G` (one base for one base) — neither
//! supplies bases while consuming a different number, so on the authored spelling
//! `delins-recommendation-reach-when-the-input-arrives-split` would keep `:47`
//! away and the pair would stay split. But that record keys on "some member of
//! the split **derived from the sequence**", and the derived split is different:
//!
//! ```text
//! FERRO_PARTITION=canonical   c.[24_27del;30C>G] -> c.[24_26del;29_30delinsG]
//! ```
//!
//! `29_30delinsG` consumes two reference bases to place one. It is gap-bearing,
//! so an inserted sequence re-aligned, so `:47` reaches the split and recommends
//! the span — which is what the coalesced arm then emits. The record anticipates
//! exactly this misreading, in the paragraph headed "THE READING THAT LOOKS RIGHT
//! AND IS MEASURABLY WRONG": the predicate keys on the gap, which the derived
//! member has, not on the rendering, which the authored one does not.
//!
//! The other two scopes hold as well: both rows are on the `c.` axis
//! (`delins-payload-coincidence-carve-out-is-coding-dna-scoped`), and the block
//! replaces 7 reference bases with 3 — a net deletion, which is the direction
//! `delins-merge-vs-individual-gap-two-or-more` is scoped to.
//!
//! So all three are accounted for and none is a `general.md:34` deviation. The
//! figure is re-pinned at 3, not waived.
//!
//! # #1610 adds 2, re-pinning 3 -> 5 — named and adjudicated on the same terms
//!
//! Re-derived by the same method (the census instrumented to print every merging
//! row, re-run armed against the prepared reference), not composed from the
//! previous figure:
//!
//! ```text
//! s01-c3m-pair-delins-del-p4-sep3  c.[24_27delinsGTCA;31_34del] -> c.24_34delinsGTCATTC
//! s01-c3p-pair-delins-del-p4-sep3  c.[24_27delinsGTCA;31_34del] -> c.24_34delinsGTCATTC
//! ```
//!
//! One variant, both shuffle directions, so the pair is one row's 3'/5' images
//! rather than two independent findings. The three rows above are unmoved.
//!
//! **These are `split_is_a_placed_gap_coincidence`'s rows**, and they are the
//! same shape as the `sep3` row already named above: the counter keys on the
//! **authored** separation, which is 3, while the rule keys on the separation of
//! the partition re-derived from the resulting sequence, which is 1
//! (`separation-is-a-property-of-the-spelling-not-of-the-variant`, decided). The
//! block is a net deletion whose derived split is one gap-bearing `delins` plus
//! the residual placed gap, one unchanged base apart, so
//! `unequal-length-block-a-placed-gap-is-not-a-separation` (decided) keeps it
//! whole.
//!
//! All three scopes hold: the rows are `c.`
//! (`delins-payload-coincidence-carve-out-is-coding-dna-scoped` — and this rule
//! is gated on that carve-out, which at the time was `c.`-only, which is why the
//! figure rose **only** here and `guard_violations` was 0 at the time — see
//! "`guard_violations` is reframed" below for what #2155 later did to that
//! carve-out's scope and to this counter's sibling), the block loses length
//! (`delins-merge-vs-individual-gap-two-or-more`'s direction scope), and some
//! derived member supplies bases
//! (`delins-recommendation-reach-when-the-input-arrives-split`).
//!
//! So these two are within the population
//! `coding-axis-merges-are-a-disclosed-general-34-deviation` (decided
//! 2026-08-13) governs, and that record's own instruction applies: a rise above
//! its pinned figure must **name the clause that carried it**. The clause is
//! `DNA/delins.md:47`, via the record named above. Re-pinned at 5, not waived.
//!
//! **#1616 is concurrently re-pinning this same counter to 8**, under a separate
//! operator ruling, on the branch that deletes the input-relative weight bound.
//! Two open PRs moving one counter to different values is a semantic conflict no
//! textual check sees — the constants live on different lines, so `merge-tree`
//! reports no conflict between the two branches. Whichever lands second must
//! re-measure this figure on the merged tree rather than keeping its own; the
//! merge queue is what forces that, since it rebuilds and re-runs against the
//! projected `main`.
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
//! # `guard_violations` is reframed: #2155 supersedes the axis scope, and this IS the merge
//!
//! Everything above this section, and every historical note inside
//! [`THREE_PRIME`]/[`FIVE_PRIME`] that credits a gate with holding
//! `guard_violations` at 0, was written against
//! `delins-payload-coincidence-carve-out-is-coding-dna-scoped` as it stood
//! through #1899: `DNA/delins.md:44-47`'s payload-coincidence carve-out reached
//! **only** `c.`, so every guarded row here — which is drawn from the
//! **frameless** axes, `g.` and non-coding multi-exon `n.` (`is_svd_wg010_shape`'s
//! domain, `Genomic | NonCodingMultiExon`) — sat entirely outside the carve-out's
//! reach. That is why `guard_violations` could be pinned at 0 as a genuine
//! negative result: the mechanism that would have produced the rejected
//! SVD-WG010 shape was, by construction, never in a position to run on these
//! rows.
//!
//! **#2155 supersedes that ruling to all DNA axes.** `CoincidenceCarveOut::for_axis`
//! now gates on `AxisFrame::is_dna` (`c./g./m./n.`) instead of the removed
//! coding-DNA-only predicate (`c.` alone) — see `merge.rs`'s doc comment on
//! `CoincidenceCarveOut` for the mechanism. The frameless axes this guard's rows
//! live on are squarely inside `is_dna`'s domain, so the carve-out now reaches
//! exactly the rows it used to be excluded from, and the guard fires:
//!
//! | direction | was (through #1899/#1610) | now (#2155) |
//! |---|---|---|
//! | 3' `guard_violations` | 0 | **12** |
//! | 5' `guard_violations` | 0 | **12** |
//!
//! Measured over the same 210 guarded rows in both directions (the denominator
//! is unmoved — the corpus's own shape, `guarded_rows: 210`, does not change).
//! Every other counter in both censuses is byte-identical to its pre-#2155
//! value, including `coding_axis_separation_two_or_more_merges` (5 / 5, unmoved
//! — its domain is `CodingSingleExon | CodingMultiExon`, disjoint from this
//! guard's `Genomic | NonCodingMultiExon`, so widening the carve-out to reach
//! `g.`/`m.`/`n.` cannot move a counter whose rows are all `c.`).
//!
//! **This is NOT a repeat of the #1899 defect, and the pin is reframed rather
//! than merely re-blessed.** The #1899 paragraph above describes an *accident*:
//! a rule shipped axis-blind and re-admitted, on `g.`/`n.`, exactly the shape a
//! sibling gate had just been scoped to exclude — a rank-1 regression, fixed by
//! adding the missing gate. What #2155 does is different in kind: it is a
//! **deliberate, operator-adjudicated widening** of the carve-out's own scope,
//! decided on the same terms as the original scoping (`delins.md:47`'s
//! simplicity ground is axis-neutral; the "incorrect protein-level predictions"
//! ground that originally justified stopping at `c.` has nothing to bite on
//! outside a reading frame either way). So the twelve rows per direction are not
//! a gate someone forgot to wire — they are the intended, disclosed consequence
//! of the ruling's own text, superseded in the same commit that widened the code.
//!
//! **Read as a rule-2 deviation, exactly as
//! `coding-axis-merges-are-a-disclosed-general-34-deviation` already reads the
//! `c.`-only case.** `general.md:33` / `DNA/delins.md:17` say a frameless
//! separation of one nucleotide is described individually; merging it anyway is
//! what the *rejected* SVD-WG010 proposal would have required
//! (`consultation/SVD-WG010.md:8`). `separation-rule-force-modal-or-negation`
//! (decided) is what makes that a deviation to **disclose and pin with a
//! tripwire**, under README rule 6's "should" reading, rather than a rule-7 bug
//! that blocks a release outright — the same classification the coding-axis
//! ruling already uses. So this pin is not "0 of 210, a fault we do not have";
//! it is "12 of 210, a disclosed house choice", and it carries the same
//! obligation that ruling states in as many words: **a rise OR a fall from 12
//! must name the clause or ruling that carried it**, not be waived or absorbed
//! silently. A fall could mean the carve-out narrowed again, or that a row's
//! irreducibility changed; a rise could mean the carve-out reached further
//! still. Neither is assumed — re-derive it the way the `coding_axis_*` rises
//! were re-derived above, by instrumenting and naming the rows.
//!
//! **`r.` stays excluded, on the same jurisdictional ground the ruling always
//! stated** — "a `DNA/` clause has no jurisdiction over the RNA axis regardless
//! of how widely the DNA scope itself is drawn" (`merge.rs`, `CoincidenceCarveOut`
//! doc comment) — and `AxisFrame::is_dna` says so exhaustively:
//! `CisKind::Rna => false`, with no wildcard arm, so a future `CisKind` variant
//! is a compile error to decide rather than a silent default. Two things confirm
//! it rather than merely asserting it:
//!
//! * **Structurally**, this corpus cannot even ask the question: `RefShape::all()`
//!   has no RNA shape (`Genomic`, `CodingSingleExon`, `CodingMultiExon`,
//!   `NonCodingMultiExon` — four variants, none of them `r.`), so 0 of the
//!   58,552 spellings measured above are on the `r.` axis. That zero is
//!   structural, exactly as `rna-axis-alignment-only-symbol-reach` records for
//!   the same corpus's blindness to `r.` elsewhere — it is not evidence that
//!   `r.` behaves correctly, only that this corpus cannot show it misbehaving.
//! * **Directly**, `the_rna_axis_stays_split_under_the_widened_carve_out` below
//!   takes a `NonCodingMultiExon` corpus row that actually merges under the
//!   widened carve-out — a net-deletion pair separated by exactly one unchanged
//!   base, where the payload happens to align with the reference — and
//!   normalizes its own authored spelling re-prefixed `r.` on the SAME
//!   provider. The `n.` spelling merges to one member; the `r.` twin, same
//!   transcript and same positions, stays a multi-member allele. One input
//!   pair, one provider, so the axis is the only thing that differs.
//!
//! # A sibling confluence gain, same cause, different population — not a guard row
//!
//! `guard_violations` is not the only counter this widening moves.
//! `s00-n3m-m3-rot5-p4-sep1` and its `n3p` twin — a **three**-member
//! non-coding design, `n.[36_39inv;41_44del;45_46insCCCA]` — diverged at 5'
//! before #2155: four of their six equivalent spellings normalized to the
//! fully authored split, the other two to a partially coalesced
//! `n.[36_39inv;42_45delinsCCCA]`, because the coalescing pass that resolves
//! that disagreement was still `c.`-only. Widening it to `n.` makes all six
//! agree on the coalesced (still two-member) form, moving [`FIVE_PRIME`]'s
//! `converged` 10,803 -> 10,805 and `split_two` 957 -> 955. **This is not a
//! guard-violation row** — `is_svd_wg010_shape` requires exactly two authored
//! members and this design has three, so it sits outside that guard's domain
//! by construction (the same disjointness `is_coding_axis_separation_two_or_more_shape`
//! already documents for its own counter). It is an ordinary rank-2 confluence
//! gain in the direction the ratchet welcomes; 3' is unaffected, verified by
//! diffing the full divergence-id list either side of the scoping change. See
//! [`FIVE_PRIME`]'s own note on `split_two` for the measurement.
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

use std::collections::BTreeSet;

use ferro_hgvs::conformance::census::{
    built, coding_axis_merge_observed, grouped, measure, members_of, report, Census, CorpusShape,
    Equivalence,
};
use ferro_hgvs::conformance::spec_corpus::{denoted_by, RowKind};
use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::{parse_hgvs, ShuffleDirection};

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
    // #1835: 371 -> 389. Eighteen more outputs are rendered against the genomic
    // reference their junction crossing resolved, because the re-derived
    // partition places a member across an exon junction where the previous one
    // did not. The class is conformant by construction — that is what #1704's
    // re-parenting made it — and its sibling `outputs_leaving_the_transcript`
    // stays at 0, which is the counter that would move if any of the eighteen
    // were being emitted on a bare transcript instead
    // (`bare-transcript-intronic-position`, decided).
    //
    // #1621 adds **2** on top of that: the two
    // `s01-c3m-junction-1-del-del-p1-sep{0,1}` rows whose far-side spelling is
    // now pulled back onto the junction's near side, where #670's
    // junction-crossing continuation then carries it into the intron as
    // `NC_SYNTH.1(NM_TEST.1):c.18_20+2A[3]`. Both clauses apply and they
    // compose: `general.md:44` forbids the exon/EXON shift, the exception 3'
    // rule permits the exon/INTRON one. Conformant, and not a defect counter —
    // see the module docs' sixth RE-BLESSED section. Re-measured on the rebase
    // onto #1835 rather than composed from the two branches' deltas: 391.
    outputs_intronic_under_a_genomic_wrapper: 391,
    // Re-blessed DOWN to ZERO. #1627 refused `standards.md:39`'s 24 `X` rows
    // rather than re-emitting them; #1628 refused the residual 8 —
    // `checklist.md:16`'s `+` offset (4) and `checklist.md:45`'s hyphen range
    // (4) — so ferro now emits no description that violates a prohibition the
    // spec states in as many words. See the module docs' fourth RE-BLESSED
    // section, and the empty clause map in
    // `corpus_prohibited_inputs::every_prohibition_violating_output_is_a_re_emitted_prohibited_input`.
    //
    // Independent of the #1704 split directly above: #1627 and #1628 refuse rows
    // on input-conformance grounds, #1704 re-parents rows on the
    // junction-crossing ground, and the two touch disjoint row sets.
    prohibition_violating_outputs: 0,
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
    // #1716 leaves `converged` alone and moves three family instances between
    // the `split_*` buckets — see the module docs' fifth RE-BLESSED section for
    // the three row ids and the diffed divergence sets.
    // # #1835 — the partition default flip, and the largest confluence move this
    // census has taken
    //
    // `converged` 9 402 -> 11 016, with `split_two` -1 550, `split_three` -62 and
    // `split_more` -2. Those three sum to exactly 1 614, the gain in `converged`,
    // so **every class that left a divergent bucket went straight to
    // convergence** and none merely dropped an arity. No divergence figure rises.
    //
    // Read it with the rank-1 counters directly above and below, which is what
    // makes the gain safe to accept rather than merely large:
    // `declined`, `unparseable_outputs`, `outputs_leaving_the_transcript` and
    // `prohibition_violating_outputs` are all still **0**, `sequence_changed`
    // and `non_idempotent_outputs` are unmoved at 4, and
    // `outputs_denoting_no_sequence` is unmoved at 10. So nothing converged by
    // losing a member, changing a base or ceasing to be a fixed point.
    //
    // # #1621 — on top of that, two families move DOWN one arity each
    //
    // `converged` is untouched — `split_more` -1 and `split_three` +1 from
    // `s01-c3p-junction-1-del-del-p1-sep0` (4 -> 3), then `split_three` -1 and
    // `split_two` +1 from `…-sep1` (3 -> 2), which is why only two of the four
    // counters below move off #1835's values. Measured by diffing the full
    // divergence row-id lists either side (the `MAX_DIVERGENCES` cap raised
    // locally, then restored), not inferred from the net counters: **no family
    // lost or gained convergence**. Re-measured on the rebase onto #1835 rather
    // than composed from the two branches' deltas: on this base the two moves
    // land in `split_two` +1 and `split_more` -1, leaving `split_three` where
    // #1835 left it, and the divergent-bucket deltas sum to zero against
    // `converged`'s zero gain.
    //
    // # #1610 — the predicted fall DID NOT HAPPEN, and the rule now takes the axis gate
    //
    // This block previously re-pinned `converged` 11 016 -> 10 959 and
    // `split_two` 673 -> 730, attributing the whole fall to
    // `canonicalize_from_sequence`'s input-relative weight bound with the
    // prediction that "disabling that bound returns `converged` to exactly
    // 11 016". **#1899 deleted that bound**, and the prediction is confirmed:
    // re-measured on this branch rebased onto it, `converged` and every
    // `split_*` bucket read exactly the values pinned here, so the re-pins are
    // reverted rather than carried forward. This change costs no confluence.
    //
    // What arrived with the bound's deletion was a rank-1 cost: axis-blind,
    // `split_is_a_placed_gap_coincidence` took `guard_violations` 0 -> 5 (the
    // rejected SVD-WG010 shape on the FRAMELESS `g.`/`n.` axes). #1899 had also
    // scoped the SIBLING collapse in `partition_block` to the coding DNA axis
    // (`CoincidenceCarveOut::may_disbelieve_a_separation`, per
    // `rulings[delins-payload-coincidence-carve-out-is-coding-dna-scoped]`), and
    // this rule sat directly after that gate re-admitting on `g.`/`n.` exactly
    // what it excludes. The rule is now gated on the same carve-out, which
    // returns `guard_violations` to 0 and leaves every counter in this block
    // where #1621 left it.
    //
    // **#2155 SUPERSEDES the carve-out scope this paragraph rests on.** The
    // ruling cited above is no longer `c.`-only — it now reaches every DNA axis,
    // `g.`/`n.` included — so the gate this paragraph credits with returning
    // `guard_violations` to 0 no longer excludes the frameless axes at all. See
    // "`guard_violations` is reframed" below for the re-measured count and why
    // this is a disclosed deviation rather than a repeat of the #1899 defect.
    //
    // # #1650 — the 3'UTR extension of the sequence-first body, on top of #1835
    //           and #1610
    //
    // Families move out of the `split_*` buckets and into `converged`, and the
    // three decreases sum exactly to `converged`'s increase. **That identity is
    // a NET statement and is quoted as one** — the module's own warning is that
    // the counters cannot tell net from gross. What is checked directly is the
    // reported divergence set, whose family ids are byte-identical either side
    // (all `s01-c1-m3-all-del-*`, none of them a `cds-end` row), so nothing
    // that was divergent before is divergent for a different reason now. The
    // four figures below are **re-measured on each rebase** rather than composed
    // from the branches' deltas — and the move is the same 40 families it was on
    // the pre-flip base, with the same 22 + 12 + 6 breakdown, so the flip, #1610
    // and this change reach disjoint populations.
    //
    // # #1816 — the CDS/3'UTR straddle anchor (per-endpoint `Anchor`)
    //
    // `merge::render_on_its_own_region` no longer refuses a member whose two
    // endpoints fall in different zones: `Anchor` now carries a region per
    // endpoint, so a piece running from `c.<cds_len-1>` to `c.*2` is *built* as
    // the spanning `c.N_*M` member the `g.` axis already produces, and the two
    // spellings of such a variant converge. Families move out of `split_two`
    // (-20) and `split_three` (-24) and into `converged` (+44); `split_more` is
    // unchanged and the two decreases sum exactly to the increase. The move is
    // meaning-preserving: `sequence_changed` holds at 4 and
    // `outputs_denoting_no_sequence` at 10 (verified against the pre-#1816 base,
    // identical), `coding_axis_separation_two_or_more_merges` holds at 19 (no new
    // `general.md:34` deviation), and every moved family is a CDS/3'UTR straddle.
    //
    // # #1816 — the 5'UTR straddle fold (`ExtendedBody`'s `FivePrimeUtr` arm)
    //
    // The mirror of the 3'UTR half above, on the CDS-start seam. `ExtendedBody`'s
    // fold/unfold now express a `c.-N` endpoint on the extended axis (flat `0` for
    // `c.-1`, descending below), so `canonicalize_from_sequence` reaches a group
    // crossing `cds_start` exactly as it reaches one crossing `cds_end`, and the
    // window low clamp is released for that shape. Families move out of `split_two`
    // (-72), `split_three` (-24) and `split_more` (-12) and into `converged`
    // (+108); the three decreases sum exactly to the increase. Meaning-preserving:
    // `sequence_changed` holds at 4 and `outputs_denoting_no_sequence` at 10,
    // `coding_axis_separation_two_or_more_merges` holds at 19 (no new
    // `general.md:34` deviation), and every moved family is a CDS/5'UTR straddle.
    // (Both shuffle-direction censuses move, since a 5'-seam family is measured
    // under each — see `FIVE_PRIME` for the 5'-direction figures.)
    converged: 11_208,
    split_two: 560,
    split_three: 104,
    split_more: 10,
    underdetermined: 0,
    // -- idempotency --
    //
    // Re-blessed DOWN, **to zero**, by #1650. The three
    // `scale-c3p-sep{120,128,136}-del-del` rows had already gone; the four
    // `cds-end` families were the whole residual and they now settle in one
    // pass, because the derivation can finally read a member on either side of
    // `cds_end` and so answers in one step what it previously answered in two.
    // `defect_non_idempotent_outputs` flips with this number, as its own failure
    // messages instruct.
    non_idempotent_outputs: 0,
    // -- sequence preservation --
    sequence_changed: 4,
    // -- refusal --
    conflicts_accepted: 72,
    // Re-blessed DOWN to zero. #1628 took the 8 genomic-offset rows (4 ×
    // `checklist.md:16-genomic-has-no-offsets`, 4 ×
    // `checklist.md:45-range-is-underscore`) and #1789 takes the 24
    // `checklist.md:33` `ins6` rows. Strict refuses each at parse (`W4009`,
    // `W3029`); lenient and silent accept the input and then fail at normalize,
    // which is why they leave the LENIENT counter this figure measures.
    prohibited_absolute_accepted: 0,
    // Re-blessed DOWN by #1627, same 24 rows as `prohibition_violating_outputs`
    // above. The residual 16 are `checklist.md:20`'s bare-transcript intronic
    // rows, which lenient is CORRECT to accept — see
    // `rulings[bare-transcript-intronic-position]`.
    prohibited_conditional_accepted: 16,
    // -- negative guards --
    //
    // #2155: 0 -> 12. `delins-payload-coincidence-carve-out-is-coding-dna-scoped`
    // is superseded from `c.`-only to all DNA axes (`c./g./m./n.`), which puts
    // this guard's own domain (the frameless `g.`/non-coding-multi-exon `n.`
    // axes) inside the widened carve-out's reach for the first time. This is a
    // DISCLOSED house deviation (README rule 6), not a fault: see the module
    // docs' "`guard_violations` is reframed" section. Any further change to this
    // figure, up or down, must name the clause or ruling that carried it —
    // exactly the obligation `coding-axis-merges-are-a-disclosed-general-34-deviation`
    // already states for its own counter below.
    guard_violations: 12,
    // -- instruments --
    coding_axis_separation_two_or_more_rows: 997,
    // #1835: 0 -> 3. The three rows are named and adjudicated individually in
    // the module docs' "The 3 are named and adjudicated" section; two are the
    // `delins.md:47` pass reaching a derived split that carries a gap-bearing
    // member, and one is a block the canonical partitioner returns whole, so it
    // has no derived separation for `general.md:34` to govern.
    //
    // #1610: 3 -> 5. The two added rows are `s01-c3{m,p}-pair-delins-del-p4-sep3`
    // — one variant in both shuffle directions — named and adjudicated in the
    // module docs' "#1610 adds 2" section. They are
    // `split_is_a_placed_gap_coincidence`'s rows; that rule is gated to the
    // coding DNA axis, which is why this instrument moves and `guard_violations`
    // directly above is now 12 (#2155 moved it 0 -> 12); it was 0 when this note
    // was first written.
    //
    // #1703: 5 -> 19. The 14 added rows are whole-span reverse complements the
    // decided `rulings[whole-span-reverse-complement-types-as-inv]` types as one
    // `inv`, uniformly, regardless of interior coincidence — 8 authored as a
    // triple of `delins` members (`m3-all-delins`, separation 2) and 6 as a pair
    // (`pair-delins-delins`, separations up to 8), every one length-preserving so
    // its resulting block equals `revcomp` of the reference block. They are NOT
    // `general.md:34` deviations: a whole-span inversion is ONE variant, so `:34`
    // ("two variants separated by...") has no antecedent (the ruling's GROUND 2),
    // and `sequence_changed` / `converged` below are unchanged in both directions,
    // so nothing lost confluence and no output denotes different bases. The pre-
    // existing 5 are the two net-deletion families (`pair-del-sub` x2,
    // `pair-inv-del` x3), which are shorter than their span and so cannot be a
    // whole-span reverse complement — this pass does not reach them.
    //
    // #2155 does not move this counter: its domain
    // (`CodingSingleExon | CodingMultiExon`) is disjoint from the widened
    // carve-out's non-coding reach (`g.`/`m.`/`n.`), so it cannot touch a counter
    // whose every row is `c.`. RE-MEASURE after rebase to confirm 19 holds.
    coding_axis_separation_two_or_more_merges: 19,
};

/// The 5'-direction census, pinned.
///
/// Confluence is a property of the normalizer rather than of one shuffle
/// direction, so it is measured in full. Two directions landing on materially
/// different numbers would mean a fix was treating a symptom of the shuffle
/// rather than the partitioner.
///
/// The 5' direction is **no longer a public option** — see `README.md` rule 6
/// and `tests/it/five_prime_public_surface_removed.rs`. It is measured in full
/// for a different reason now: the 5' arm is ferro's differential oracle over its own 3' output — the instrument that found #1542, where 7 of 8 `FERRO_PARTITION` x direction configurations agreed and only the shipped `live`/3' arm diverged.
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
    // Re-blessed DOWN by #1627, by the same 24 rows as at 3', and then to ZERO
    // by #1628, by the same 8. Validity does not depend on shuffle direction, so
    // the two directions moving identically is the cross-check rather than a
    // coincidence.
    prohibition_violating_outputs: 0,
    // Unmoved by #1599 or by #1536 — no 5' family changed which side of the
    // converged/split line it sat on in either. **#1649 moves it for the first
    // time**, by 284, and the same three measured zeros hold as at 3': no family
    // loses convergence, rises in arity, or becomes divergent that was not.
    // #1627 does not move it either: a refused row contributes no family.
    // #1542 (PR #1840) moves it +2 — see the note on `split_two` below.
    //
    // #1610 does not move it: see the `#1610` block on [`THREE_PRIME`]'s
    // `converged` for the measurement and for why the rule takes the axis gate.
    //
    // Re-blessed by #1650, the 5' half of the 3' re-bless above: families move
    // into `converged` and the decreases sum to it exactly. Read as net, for the
    // same reason stated there, and re-measured on each rebase rather than
    // composed — the move is the same 48 families it was on the pre-flip base,
    // with the same 22 + 26 + 0 breakdown.
    //
    // #1816 (the CDS/3'UTR straddle anchor — see the `#1816` block on
    // [`THREE_PRIME`]'s `converged`) moves it a further +46: the straddle families
    // converge in the 5'-shuffle direction too. The decreases are `split_two` -22,
    // `split_three` -22, `split_more` -2, summing exactly to +46, and the move is
    // meaning-preserving (`sequence_changed` holds at 0,
    // `coding_axis_separation_two_or_more_merges` at 19).
    //
    // #1816's 5'UTR fold (`ExtendedBody`'s `FivePrimeUtr` arm — see the same
    // `#1816` block on [`THREE_PRIME`]) moves it a further +102: the CDS/5'UTR
    // straddle families converge in the 5'-shuffle direction. The decreases are
    // `split_two` -64, `split_three` -26, `split_more` -12, summing exactly to
    // +102, and the move is meaning-preserving (`sequence_changed` holds at 0,
    // `coding_axis_separation_two_or_more_merges` at 19).
    converged: 10_953,
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
    //
    // #1716 moves the same three family instances as at 3', all three
    // *improving* on this side; `converged` is unchanged. Fifth RE-BLESSED
    // section.
    // #1835: `converged` 9 228 -> 10 753, `split_two` -1 481, `split_three` -44,
    // `split_more` unchanged — again summing exactly to the gain, so every moved
    // class converged outright. The 5' direction gains 1 525 against the 3'
    // direction's 1 614; the two moving by different amounts is the asymmetry
    // every entry above records, and both move the same way.
    //
    // #1542 (PR #1840): `converged` 10 753 -> 10 755, `split_two` 981 -> 979,
    // `split_three` and `split_more` unchanged — the two deltas cancel, so both
    // moved families converged outright and none rose in arity. **Every other
    // field in this census is unmoved**, so nothing traded a rank-2 gain for a
    // rank-1 loss.
    //
    // This is the FIRST entry to move 5' and leave 3' alone: `THREE_PRIME` is
    // untouched by the same change. That asymmetry is the point rather than a
    // side effect — the change makes the *member count* independent of the
    // shuffle direction, which it does by having each direction adopt the other's
    // partition where that one merges and this one did not. On these two families
    // the 3' partition was already the merged one, so only 5' had anywhere to
    // move; convergence is what agreeing with 3' looks like from here.
    //
    // every entry above records, and both move the same way. #1650 then moves
    // these on top of that, re-measured on the rebase: -22 and -26, with
    // `split_more` untouched.
    //
    // #2155: `converged` 10 803 -> 10 805, `split_two` 957 -> 955 — the SAME
    // widened carve-out that produces `guard_violations`, on a DIFFERENT and
    // disjoint population. `s00-n3m-m3-rot5-p4-sep1` and its `n3p` twin
    // (`n.[36_39inv;41_44del;45_46insCCCA]`, a THREE-member design) diverged
    // before this change: four of their six spellings normalized to the fully
    // authored split, the other two to a partially coalesced
    // `n.[36_39inv;42_45delinsCCCA]`, because the coalescing pass that resolves
    // that difference was still `c.`-only. Widening it to `n.` makes all six
    // spellings agree on the coalesced form. **This is not a guard-violation
    // row**: `is_svd_wg010_shape` requires exactly two authored members and
    // this design has three, so it is outside that guard's domain by
    // construction — it is an ordinary rank-2 confluence gain, in the
    // direction the ratchet welcomes, and needed no separate adjudication
    // beyond naming the cause. Verified by diffing the full divergence-id list
    // against the pre-#2155 scoping (the removed coding-DNA-only predicate, at
    // all three `merge.rs` call sites `CoincidenceCarveOut::for_axis`,
    // `payload_coalesce_applies` and `compensating_gap_coalesce_applies`
    // gate): those two ids are the only ones removed, none added, and 3' is
    // byte-identical on the same diff (its `THREE_PRIME` pins are untouched by
    // this change).
    // #1816 (3'UTR anchor) moves all three DOWN — `split_two` -22, `split_three`
    // -22, `split_more` -2 — as the CDS/3'UTR straddle families converge; the
    // three sum to `converged`'s +46 above. Its 5'UTR fold then moves all three
    // DOWN again — `split_two` -64, `split_three` -26, `split_more` -12 — as the
    // CDS/5'UTR straddle families converge; those three sum to `converged`'s +102.
    split_two: 869,
    split_three: 50,
    split_more: 10,
    underdetermined: 0,
    // Re-blessed DOWN to zero by #1650, the same four `cds-end` families as at 3'.
    non_idempotent_outputs: 0,
    sequence_changed: 0,
    conflicts_accepted: 72,
    // Re-blessed DOWN to zero by #1628 and #1789, by the same 32 rows as at 3'.
    prohibited_absolute_accepted: 0,
    // Re-blessed DOWN by #1627, by the same 24 rows as at 3'.
    prohibited_conditional_accepted: 16,
    // #2155: 0 -> 12, identically to 3' — see [`THREE_PRIME`]'s note on this
    // field. The two directions moving by the same amount is the cross-check
    // this pair of censuses exists for: the widened carve-out is not a property
    // of one shuffle direction.
    guard_violations: 12,
    coding_axis_separation_two_or_more_rows: 997,
    // #1835: 0 -> 3. The three rows are named and adjudicated individually in
    // the module docs' "The 3 are named and adjudicated" section; two are the
    // `delins.md:47` pass reaching a derived split that carries a gap-bearing
    // member, and one is a block the canonical partitioner returns whole, so it
    // has no derived separation for `general.md:34` to govern.
    //
    // #1610: 3 -> 5. The two added rows are `s01-c3{m,p}-pair-delins-del-p4-sep3`
    // — one variant in both shuffle directions — named and adjudicated in the
    // module docs' "#1610 adds 2" section. They are
    // `split_is_a_placed_gap_coincidence`'s rows; that rule is gated to the
    // coding DNA axis, which is why this instrument moves and `guard_violations`
    // directly above is now 12 (#2155 moved it 0 -> 12); it was 0 when this note
    // was first written.
    //
    // #1703: 5 -> 19. The 14 added rows are whole-span reverse complements the
    // decided `rulings[whole-span-reverse-complement-types-as-inv]` types as one
    // `inv`, uniformly, regardless of interior coincidence — 8 authored as a
    // triple of `delins` members (`m3-all-delins`, separation 2) and 6 as a pair
    // (`pair-delins-delins`, separations up to 8), every one length-preserving so
    // its resulting block equals `revcomp` of the reference block. They are NOT
    // `general.md:34` deviations: a whole-span inversion is ONE variant, so `:34`
    // ("two variants separated by...") has no antecedent (the ruling's GROUND 2),
    // and `sequence_changed` / `converged` below are unchanged in both directions,
    // so nothing lost confluence and no output denotes different bases. The pre-
    // existing 5 are the two net-deletion families (`pair-del-sub` x2,
    // `pair-inv-del` x3), which are shorter than their span and so cannot be a
    // whole-span reverse complement — this pass does not reach them.
    //
    // #2155 does not move this counter — see [`THREE_PRIME`]'s identical note.
    // RE-MEASURE after rebase to confirm 19 holds.
    coding_axis_separation_two_or_more_merges: 19,
};

/// Assert one direction's census against its pin, printing the measured numbers
/// either way so a moved pin can be re-blessed from the test output.
fn assert_census(direction: ShuffleDirection, label: &str, pinned: &Census) {
    let measured = measure(direction, Equivalence::ExactString);
    println!("{}", report(label, &measured));
    let census = &measured.census;

    // The honest-zero discipline: a property whose denominator is zero is
    // VACUOUS, not passing.
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
        // `guard_violations` is deliberately NOT in this loop as of #2155 — see
        // the dedicated assert below. Its rise is no longer automatically "a
        // rank-1 conformance regression, not a representation choice"; folding
        // it into this loop would keep asserting a framing that is no longer
        // true.
    ] {
        assert!(
            measured_value <= pinned_value,
            "{label}: {what} rose from {pinned_value} to {measured_value}. This is a rank-1 \
             conformance regression, not a representation choice.\n{}",
            report(label, &measured)
        );
    }

    // `guard_violations` is a DISCLOSED deviation as of #2155, not a fault: the
    // frameless payload-coincidence carve-out this guard measures is now a
    // named house choice under
    // `delins-payload-coincidence-carve-out-is-coding-dna-scoped` (superseded
    // from `c.`-only to all DNA axes), not "ferro does not implement rejected
    // SVD-WG010" — see the module docs' "`guard_violations` is reframed"
    // section. So neither direction gets the generic loop's "rank-1
    // regression" framing: a RISE is not automatically a fault, and — unlike
    // every counter in that loop — a FALL is not automatically welcome either,
    // since it could mean the carve-out narrowed in a way nobody decided. The
    // obligation is symmetric: any move away from the pin, either way, must
    // name the clause or ruling that carried it and re-pin in the same commit.
    // `assert_eq!` both directions rather than `<=`, deliberately.
    assert_eq!(
        census.guard_violations,
        pinned.guard_violations,
        "{label}: outputs merging the rejected-SVD-WG010 shape on a frameless DNA axis moved \
         from {} to {}. Read this WITH the module docs' \"`guard_violations` is reframed\" \
         section before re-pinning either way: a rise or a fall must name the clause or ruling \
         that carried it, and neither is to be waived or assumed correct from the number \
         alone.\n{}",
        pinned.guard_violations,
        census.guard_violations,
        report(label, &measured)
    );

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
/// # Why a zero pin needed one, and why nothing else supplied it
///
/// Written when `coding_axis_separation_two_or_more_merges` was pinned at
/// **0**, so every other test that mentioned it was satisfied by a numerator
/// that could never increment. Measured as a mutation matrix over the four
/// tests that claim to cover the counter — the two censuses, the corpus-shape
/// test and `spec_corpus`'s population test — hard-wiring
/// [`coding_axis_merge_observed`] to `false` left **all four green**. "Pinned
/// at 0" and "wired to 0" were the same observation, which is exactly what this
/// module's own thesis says an instrument must not be: *a counter nobody has
/// seen move is not an instrument*.
///
/// **The pin is now 3** (#1835), which closes that particular blindness on its
/// own — a counter wired to `false` would read 0 against a pin of 3 and fail
/// both censuses' `assert_eq!`. This control is kept regardless: it is the only
/// test that exercises the predicate *directly* rather than through a census
/// total, and the blindness returns the moment the pin returns to 0.
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

/// One payload-coincidence shape `guard_violations` fires on (#2155), taken
/// directly from a corpus row rather than hand-built — `RefShape::all()` has
/// no `r.` shape, so the corpus itself cannot ask this question, but a
/// `RefShape::NonCodingMultiExon` row's provider is a real transcript, and the
/// same transcript serves `r.` and `n.` identically when it carries no CDS
/// (`cis_confluence_nr_axis.rs` establishes that equivalence over a much
/// larger corpus: with no CDS, `r.` numbering is the plain transcript offset,
/// same as `n.`). So the row's own authored `n.` spelling, re-prefixed `r.`,
/// names the exact same edit on the exact same provider — the only thing that
/// differs between the two assertions below is the axis.
///
/// See the module docs' "`guard_violations` is reframed" section for the wider
/// argument; this is its direct half.
#[test]
fn the_rna_axis_stays_split_under_the_widened_carve_out() {
    use ferro_hgvs::conformance::spec_corpus::RefShape;
    use ferro_hgvs::{HgvsVariant, NormalizeConfig, Normalizer};

    /// A parsed description's member count: an allele's own variants, or one.
    /// Reads no reference, mirroring `case_harvest.rs`'s `member_count`.
    fn member_count(description: &str) -> usize {
        match parse_hgvs(description)
            .unwrap_or_else(|e| panic!("{description} does not parse: {e}"))
        {
            HgvsVariant::Allele(allele) => allele.variants.len(),
            _ => 1,
        }
    }

    // Candidates only say a merge WOULD be a guard violation, not that the
    // normalizer actually merges them — `s00-n3m-junction-1-sub-sub-p1-sep1`
    // (two substitutions) is a candidate that stays split, because neither
    // member is gap-bearing. So the row that actually demonstrates the guard is
    // found by normalizing each candidate and taking the first that merges,
    // exactly mirroring what the census's own `guard_violations` counter does.
    let built = built();
    let (row, n_output, normalizer) = built
        .rows
        .iter()
        .filter(|row| {
            matches!(row.shape, RefShape::NonCodingMultiExon(_)) && !row.negative_guards.is_empty()
        })
        .find_map(|row| {
            let frame = row.frame();
            let normalizer =
                Normalizer::with_config(frame.provider().clone(), NormalizeConfig::default());
            let n_input = parse_hgvs(row.authored_spelling()).ok()?;
            let n_output = normalizer.normalize(&n_input).ok()?.to_string();
            (member_count(&n_output) == 1).then_some((row, n_output, normalizer))
        })
        .expect(
            "VACUOUS: no NonCodingMultiExon negative-guard candidate actually merges, so this \
             test has no real transcript to build the r. twin from",
        );

    let (accession, n_body) = row
        .authored_spelling()
        .split_once(":n.")
        .unwrap_or_else(|| {
            panic!(
                "{}: authored spelling is not n. — {}",
                row.id,
                row.authored_spelling()
            )
        });

    let r_description = format!("{accession}:r.{n_body}");
    let r_input = parse_hgvs(&r_description)
        .unwrap_or_else(|e| panic!("{r_description}: r. twin does not parse: {e}"));
    let r_output = normalizer
        .normalize(&r_input)
        .expect("r. normalizes")
        .to_string();
    assert!(
        member_count(&r_output) > 1,
        "{}: its n. twin ({}) merged to {n_output}, but the r. twin of the same shape \
         ({r_description}) must stay split — `DNA/delins.md:47` has no jurisdiction over the \
         RNA axis regardless of how widely the DNA scope is drawn; got {r_output}",
        row.id,
        row.authored_spelling()
    );
}

#[test]
fn three_prime_conformance_census() {
    assert_census(ShuffleDirection::ThreePrime, "3prime", &THREE_PRIME);
}

#[test]
fn five_prime_conformance_census() {
    assert_census(ShuffleDirection::FivePrime, "5prime", &FIVE_PRIME);
}
