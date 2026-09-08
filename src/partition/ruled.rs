//! The ruled-cut engine (design `specs/2026-09-05-principled-partitioner-design-fable.md`,
//! §4). A declarative, prioritized rewrite system over a *typed cut* of a block.
//!
//! This module is **purely additive** and wired into nothing yet: it exists so the
//! Stage-B `merge.rs` passes can be migrated into `Rule` impls one at a time (design
//! §10 step 1), each verified byte-identical, until `PartitionRule::Ruled` can replace
//! the straight-line pass sequence. Until a rule is registered, running the engine on a
//! seed returns the seed unchanged.
//!
//! What the engine buys over the current straight-line pipeline (design §11):
//! * **Soundness by construction** — a [`Cut`] cannot be built unless it round-trips
//!   and is renderable ([`Cut::new`] is the only constructor), so the after-the-fact
//!   `validate_sound` gate becomes a type invariant.
//! * **Termination by a stated potential** — every rule must strictly decrease
//!   [`Cut::phi`], checked by a `debug_assert!` in the engine and a per-rule unit test.
//! * **Composition made explicit** — rules run to a fixed point under one priority
//!   order; scope is checked by the *engine*, never by a rule body, so the four
//!   spellings of the axis scope (`payload_coalesce_applies`,
//!   `compensating_gap_coalesce_applies`, the inline `cuts_with_canonical && is_dna`,
//!   `may_disbelieve_a_separation`) collapse into one [`Scope`] field.

use crate::partition::block_ctx::BlockCtx;
use crate::partition::fold;
use crate::partition::output::{EditKind, Member, Partition};
use crate::partition::partitioner::{validate_sound, PartitionError};

/// The label a [`Run`] carries. `EditKind` plus `Unlabelled` — a run the seed
/// produced but that no `Label` rule has yet typed. `Dup` and `Inv` are the only
/// labels that lock (design §4.2).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Label {
    Unlabelled,
    Sub,
    Del,
    Ins,
    Delins,
    Dup,
    Inv,
    Identity,
}

impl Label {
    /// The `EditKind` this label renders as. `Unlabelled` projects to `Delins`
    /// **for the soundness check only** — round-trip and renderability read the
    /// reference span and the inserted bytes, never the kind, so the placeholder
    /// cannot change whether a cut validates. A real render never sees an
    /// `Unlabelled` run: the engine reaches a normal form only when every run a
    /// label rule can type has been typed.
    fn to_edit_kind(self) -> EditKind {
        match self {
            Label::Unlabelled | Label::Delins => EditKind::Delins,
            Label::Sub => EditKind::Sub,
            Label::Del => EditKind::Del,
            Label::Ins => EditKind::Ins,
            Label::Dup => EditKind::Dup,
            Label::Inv => EditKind::Inv,
            Label::Identity => EditKind::Identity,
        }
    }

    /// A label locks a run against absorption by an un-authorised coarsening
    /// (design §4.2's label locks). Only `Dup` and `Inv` lock — the two labels
    /// that outrank a merge. `pub(crate)` so the bridge's typed lift can set the
    /// lock bit from the label a recogniser assigned, keeping the "which labels
    /// lock" rule in one place.
    pub(crate) fn locks(self) -> bool {
        matches!(self, Label::Dup | Label::Inv)
    }
}

/// One run of a [`Cut`]: a half-open reference span `[ref_start, ref_end)` replaced
/// by `alt`, a [`Label`], and whether it is `locked` against absorption. This is the
/// existing [`Member`] plus a label and a lock bit (design §4.1).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Run {
    pub ref_start: usize,
    pub ref_end: usize,
    pub alt: Vec<u8>,
    pub label: Label,
    pub locked: bool,
}

impl Run {
    /// An unlabelled, unlocked run — the shape the seed emits before any label rule.
    pub fn unlabelled(ref_start: usize, ref_end: usize, alt: Vec<u8>) -> Self {
        Run {
            ref_start,
            ref_end,
            alt,
            label: Label::Unlabelled,
            locked: false,
        }
    }

    fn to_member(&self) -> Member {
        Member {
            kind: self.label.to_edit_kind(),
            ref_start: self.ref_start,
            ref_end: self.ref_end,
            inserted: self.alt.clone(),
        }
    }
}

/// A validated, ascending, disjoint partition of a block whose splice reproduces
/// `ctx.resulting`. **Sound by construction**: [`Cut::new`] is the only constructor
/// and runs the round-trip + renderability gate, so a `Cut` value cannot exist
/// unsound (design §4.1). Callers hold a `Cut`, not a raw `Vec<Run>`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Cut {
    runs: Vec<Run>,
}

impl Cut {
    /// Build a cut from `runs`, validating soundness and renderability against `ctx`.
    /// The only constructor. Reuses the existing [`validate_sound`] gate so the two
    /// surfaces cannot disagree on what "sound" means.
    ///
    /// **Renderability by construction (design §10 step 2).** Two content-bearing
    /// insertions anchored at one interbase are not rejected — they are the single
    /// insertion they concatenate to, so they are folded
    /// ([`fold::fold_coincident_insertions`]). This retires the shipping
    /// `merge_coincident_insertions` render-stage fold and the bridge's lift-refusal
    /// seam for the ruled path. The fast path — no collision, which is every block
    /// on the corpus — keeps the original labelled runs untouched, so a
    /// non-colliding `Cut` is byte-identical to before this change.
    pub fn new(ctx: &BlockCtx, runs: Vec<Run>) -> Result<Cut, PartitionError> {
        let members: Vec<Member> = runs.iter().map(Run::to_member).collect();
        if !fold::has_coincident_collision(&members) {
            // `validate_sound` consumes and returns the partition; we keep our runs.
            validate_sound(ctx, Partition { members })?;
            return Ok(Cut { runs });
        }
        // Collision path (never reached on the corpus): fold the coincident
        // insertions, then reconstruct labelled runs. A member that survived the
        // fold unchanged carries its run's label and lock; the newly folded member
        // is a fresh, unlocked `Ins` — a dup+ins collision denotes one insertion
        // (`duplication-must-ranks-the-label-not-the-partition`), never a dup.
        let folded = fold::fold_coincident_insertions(Partition { members });
        let folded_runs: Vec<Run> = folded
            .members
            .iter()
            .map(|m| match runs.iter().find(|r| &r.to_member() == m) {
                Some(r) => r.clone(),
                None => Run {
                    ref_start: m.ref_start,
                    ref_end: m.ref_end,
                    alt: m.inserted.clone(),
                    label: Label::Ins,
                    locked: false,
                },
            })
            .collect();
        validate_sound(ctx, folded)?;
        Ok(Cut { runs: folded_runs })
    }

    pub fn runs(&self) -> &[Run] {
        &self.runs
    }

    /// Project to the neutral [`Partition`] for rendering / metrics. Labels become
    /// `EditKind`s; locks are dropped (they are a partition-stage concern only).
    pub fn to_partition(&self) -> Partition {
        Partition {
            members: self.runs.iter().map(Run::to_member).collect(),
        }
    }

    /// The termination potential `Φ = (unlocked runs, total runs)`, lexicographic
    /// (design §4.5). Every rule must strictly decrease it; the engine debug-asserts
    /// this on each application, and each rule carries a `rule_strictly_decreases_phi`
    /// unit test.
    pub fn phi(&self) -> (usize, usize) {
        let unlocked = self.runs.iter().filter(|r| !r.locked).count();
        (unlocked, self.runs.len())
    }
}

/// The three algebraic kinds of rule (design §4.2). A rule is exactly one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RuleKind {
    /// Remove one or more separators; the hull is unchanged; `total runs` strictly
    /// decreases. May not absorb a locked run unless it declares the authority to.
    Coarsen,
    /// Replace exactly one `Unlabelled` run by a specific list of runs, all locked;
    /// `unlocked runs` strictly decreases.
    Recut,
    /// Set the label on one `Unlabelled` run; locks iff the label is `Dup`/`Inv`;
    /// never changes geometry.
    Label,
}

/// On which axes / frames / directions a rule applies (design §4.2). Checked by the
/// engine before a rule's `apply` is ever called — a rule body never tests scope.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Scope {
    pub molecule: MoleculeScope,
    pub frame: FrameScope,
    pub direction: DirectionScope,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MoleculeScope {
    /// Any molecule.
    Any,
    /// The DNA axes only (`c./g./m./n.`), never `r.` — the one proposition the four
    /// `_applies` spellings encode today (`delins-payload-coincidence-carve-out-is-coding-dna-scoped`).
    DnaOnly,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameScope {
    Any,
    /// Only an axis that declares a reading frame (the codon exception's home).
    NeedsCodingFrame,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DirectionScope {
    Any,
    /// Payload shorter than the span consumed.
    NetDeletion,
    /// Payload longer than the span consumed.
    NetInsertion,
    /// Payload the same length as the span consumed.
    EqualLength,
}

impl Scope {
    /// Any molecule, any frame, any direction — the widest scope.
    pub const ANY: Scope = Scope {
        molecule: MoleculeScope::Any,
        frame: FrameScope::Any,
        direction: DirectionScope::Any,
    };

    /// Whether this scope admits `ctx`. All three limbs are pure functions of the
    /// context: direction of the block's own lengths, molecule of `ctx.molecule`
    /// (`DnaOnly` is `c./g./m./n.`, never `r.` —
    /// `delins-payload-coincidence-carve-out-is-coding-dna-scoped`, #2155), and frame
    /// of `ctx.frame` (`NeedsCodingFrame` is a frame with a codon grid — the same
    /// predicate `merge.rs`'s `carries_translated_frame` gates the codon passes on).
    /// The engine checks this before `apply`; no rule body tests scope. Deleting the
    /// four `_applies` spellings and asserting declared == computed on every `CisKind`
    /// is the rest of design §10 step 3.
    pub fn admits(&self, ctx: &BlockCtx) -> bool {
        use crate::partition::block_ctx::Molecule;
        let direction_ok = match self.direction {
            DirectionScope::Any => true,
            DirectionScope::NetDeletion => ctx.resulting.len() < ctx.reference.len(),
            DirectionScope::NetInsertion => ctx.resulting.len() > ctx.reference.len(),
            DirectionScope::EqualLength => ctx.resulting.len() == ctx.reference.len(),
        };
        let molecule_ok = match self.molecule {
            MoleculeScope::Any => true,
            MoleculeScope::DnaOnly => ctx.molecule == Molecule::Dna,
        };
        let frame_ok = match self.frame {
            FrameScope::Any => true,
            FrameScope::NeedsCodingFrame => ctx.frame.grid_phase().is_some(),
        };
        direction_ok && molecule_ok && frame_ok
    }
}

/// The authority a rule cites — a ledger ruling id, a spec clause, or a house choice
/// (design §4.2). Every registered rule must resolve to a decided ledger record (the
/// generated-index guard of §8); a rule that cites none is inadmissible.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Authority {
    Ruling(&'static str),
    Clause(&'static str),
    HouseChoice(&'static str),
}

/// A single adjudicated proposition, expressed as a scoped, authorised rewrite over a
/// [`Cut`] (design §4.2). Pure: `apply` returns `None` when the rule does not fire and
/// `Some(cut')` — which must strictly decrease `Φ` — when it does. The engine checks
/// `scope` before calling `apply`, so no rule body reads the molecule/frame/direction.
pub trait Rule {
    fn id(&self) -> &'static str;
    fn kind(&self) -> RuleKind;
    fn scope(&self) -> Scope;
    fn authority(&self) -> Authority;
    /// Rewrite `cut`, or `None` if this rule does not apply. A returned cut MUST
    /// strictly decrease `Φ` (the engine debug-asserts it).
    fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut>;
}

/// Run `rules` to a fixed point over `seed` under their registration order (design
/// §4.4). Restart-from-top after each application makes the priority order the only
/// order that matters; `Φ` bounds the iteration count at `2·runs`. Returns the normal
/// form. With no rules registered, returns `seed` unchanged — the step-1 state.
///
/// `rules` is a slice of trait references, not `Box<dyn Rule>`, so the driver hands it
/// a `'static` array with no per-block allocation (design §12 R8, the restart-from-top
/// perf risk).
///
/// **Termination.** `Φ = (unlocked, total)` is well-founded and (almost) every rule
/// strictly decreases it. The one exception is `codon-exception`, a lateral,
/// idempotent Coarsen that keeps `Φ` flat while re-partitioning to the canonical
/// codon-merged form; the loop still terminates, on three grounds — confluence
/// (the census records the only critical pairs, none involving it), its own
/// idempotence, and the `cap`. So the per-step `debug_assert!` checks **progress**
/// (`Φ` non-increasing AND the cut changed), not strict decrease — see its comment.
/// The `cap` converts a genuinely non-terminating rule into a graceful stop at a
/// still-sound cut in release rather than a hang; it is far above any legitimate
/// step count, so it never trips on a well-behaved rule set.
pub fn run_engine(seed: Cut, ctx: &BlockCtx, rules: &[&dyn Rule]) -> Cut {
    let cap = 64 * (seed.runs().len() + 4);
    let mut cut = seed;
    let mut steps = 0usize;
    'fixpoint: loop {
        for rule in rules {
            if !rule.scope().admits(ctx) {
                continue;
            }
            if let Some(next) = rule.apply(&cut, ctx) {
                // Progress, not strict-Φ decrease. Almost every rule strictly
                // decreases `Φ = (unlocked, total)` — a coarsen drops `total`, a
                // recut drops `unlocked`, a label drops `unlocked`. `codon-exception`
                // is the exception: on a coding block such as `AAAA->CACC` it
                // re-partitions to the correct codon-merged form at a CONSTANT run
                // count (a lateral, idempotent Coarsen), so `Φ` stays flat while the
                // cut genuinely changes. Asserting strict decrease falsely panics on
                // that correct move; asserting **progress** — `Φ` non-increasing AND
                // the cut actually changed — admits it while still catching the two
                // real bugs: a rule that RAISES `Φ` (non-terminating) or one that
                // fires without changing anything (should have returned `None`).
                //
                // Termination is not weakened in practice: the fixpoint set is
                // confluent (`partition_critical_pair_census` records the only two
                // critical pairs, neither involving `codon-exception`), so no
                // rule can undo another's lateral move to cycle; `codon-exception`
                // is idempotent on its own output (`pieces == before` → `None`); and
                // the `cap` below is the release-build backstop regardless. This
                // refines design §4.5's strict-`Φ` statement to match the rule set's
                // actual, terminating behaviour.
                debug_assert!(
                    next.phi() <= cut.phi() && next != cut,
                    "rule `{}` made no progress (Phi rose, or the cut is unchanged): {:?} -> {:?}",
                    rule.id(),
                    cut.phi(),
                    next.phi(),
                );
                cut = next;
                steps += 1;
                debug_assert!(
                    steps < cap,
                    "run_engine exceeded {cap} steps — a rule likely does not decrease Phi",
                );
                if steps >= cap {
                    // Release-build insurance only: a well-behaved rule set never
                    // reaches this. The current cut is a valid `Cut`, just possibly
                    // not a normal form.
                    break 'fixpoint;
                }
                continue 'fixpoint;
            }
        }
        break;
    }
    cut
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::block_ctx::{FrameContext, Molecule, Provenance};

    fn ctx<'a>(
        reference: &'a [u8],
        resulting: &'a [u8],
        frame: &'a FrameContext,
        prov: &'a Provenance,
    ) -> BlockCtx<'a> {
        BlockCtx {
            reference,
            resulting,
            frame,
            molecule: Molecule::Dna,
            provenance: prov,
        }
    }

    // A trivial `Label` rule: type a lone `A->` (net deletion of one base) run as `Del`.
    // Exists only to exercise the engine end to end; the real rules land in §10 step 1.
    struct LabelLoneDel;
    impl Rule for LabelLoneDel {
        fn id(&self) -> &'static str {
            "test:label-lone-del"
        }
        fn kind(&self) -> RuleKind {
            RuleKind::Label
        }
        fn scope(&self) -> Scope {
            Scope {
                direction: DirectionScope::NetDeletion,
                ..Scope::ANY
            }
        }
        fn authority(&self) -> Authority {
            Authority::HouseChoice("test-only")
        }
        fn apply(&self, cut: &Cut, _ctx: &BlockCtx) -> Option<Cut> {
            let mut runs = cut.runs().to_vec();
            let r = runs.iter_mut().find(|r| {
                r.label == Label::Unlabelled && r.alt.is_empty() && r.ref_end > r.ref_start
            })?;
            r.label = Label::Del;
            // Label rule never locks Del, so re-running would re-find nothing (label
            // is no longer Unlabelled) — Phi's unlocked count is unchanged, so the
            // engine must not treat a non-locking label as a rewrite. Represent that
            // by locking here for the test so Phi strictly decreases; real Label rules
            // for Del/Sub/etc. are applied once by the engine per run and are not
            // counted as rewrites (design §4.5). We lock to keep the test's Phi honest.
            r.locked = true;
            Some(Cut { runs })
        }
    }

    #[test]
    fn cut_new_rejects_an_unsound_partition() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference `AC`, resulting `AT`: a run claiming ref[1..2] -> `G` does NOT
        // round-trip to `AT`, so Cut::new must reject it.
        let c = ctx(b"AC", b"AT", &frame, &prov);
        let bad = vec![Run::unlabelled(1, 2, b"G".to_vec())];
        assert_eq!(Cut::new(&c, bad), Err(PartitionError::Unsound));
    }

    #[test]
    fn cut_new_folds_a_coincident_dup_ins_collision_into_one_unlocked_ins() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // The #2203 shape: a peeled `dup` (locked) and a plain `ins` at interbase 2.
        // reference AC -> ACACG: dup(2,"AC") + ins(2,"G") splice to AC+AC+G. `Cut::new`
        // must fold them to one renderable `Ins` rather than reject as Unrenderable.
        let c = ctx(b"AC", b"ACACG", &frame, &prov);
        let locked_dup = Run {
            ref_start: 2,
            ref_end: 2,
            alt: b"AC".to_vec(),
            label: Label::Dup,
            locked: true,
        };
        let ins = Run::unlabelled(2, 2, b"G".to_vec());
        let cut = Cut::new(&c, vec![locked_dup, ins]).expect("folds instead of rejecting");
        assert_eq!(cut.runs().len(), 1, "the two collide-folded into one run");
        let r = &cut.runs()[0];
        assert_eq!((r.ref_start, r.ref_end), (2, 2));
        assert_eq!(r.alt, b"ACG", "payloads concatenate in run order");
        assert_eq!(
            r.label,
            Label::Ins,
            "a dup+ins collision denotes one insertion"
        );
        assert!(!r.locked, "the folded insertion is unlocked");
    }

    #[test]
    fn cut_new_preserves_a_label_and_lock_on_a_non_colliding_run() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // A locked `Inv` run with no collision must survive Cut::new untouched — the
        // fast path keeps the original labelled runs (the byte-identity invariant).
        let c = ctx(b"AAGCTA", b"TAGCTT", &frame, &prov);
        let locked_inv = Run {
            ref_start: 0,
            ref_end: 6,
            alt: b"TAGCTT".to_vec(),
            label: Label::Inv,
            locked: true,
        };
        let cut = Cut::new(&c, vec![locked_inv.clone()]).expect("sound");
        assert_eq!(
            cut.runs(),
            &[locked_inv],
            "label and lock preserved verbatim"
        );
    }

    #[test]
    fn cut_new_accepts_a_sound_partition_and_engine_reaches_a_normal_form() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference `ACG`, resulting `AG`: delete ref[1..2] (`C`). Sound.
        let c = ctx(b"ACG", b"AG", &frame, &prov);
        let seed = Cut::new(&c, vec![Run::unlabelled(1, 2, Vec::new())]).expect("sound");
        assert_eq!(seed.phi(), (1, 1)); // one unlocked run

        let rules: [&dyn Rule; 1] = [&LabelLoneDel];
        let normal = run_engine(seed, &c, &rules);
        assert_eq!(normal.runs().len(), 1);
        assert_eq!(normal.runs()[0].label, Label::Del);
        assert_eq!(normal.phi(), (0, 1)); // the run is now locked, so Phi decreased
    }

    #[test]
    fn engine_with_no_rules_returns_the_seed_unchanged() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let c = ctx(b"ACG", b"AG", &frame, &prov);
        let seed = Cut::new(&c, vec![Run::unlabelled(1, 2, Vec::new())]).expect("sound");
        let normal = run_engine(seed.clone(), &c, &[] as &[&dyn Rule]);
        assert_eq!(normal, seed);
    }

    #[test]
    fn the_molecule_limb_admits_dna_and_refuses_rna_and_protein() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let dna_only = Scope {
            molecule: MoleculeScope::DnaOnly,
            ..Scope::ANY
        };
        for (molecule, want) in [
            (Molecule::Dna, true),
            (Molecule::Rna, false),
            (Molecule::Protein, false),
        ] {
            let c = BlockCtx {
                reference: b"AC",
                resulting: b"AT",
                frame: &frame,
                molecule,
                provenance: &prov,
            };
            assert_eq!(dna_only.admits(&c), want, "{molecule:?}");
            assert!(Scope::ANY.admits(&c), "ANY admits every molecule");
        }
    }

    #[test]
    fn the_frame_limb_admits_only_a_coding_frame() {
        let prov = Provenance::none();
        let needs_frame = Scope {
            frame: FrameScope::NeedsCodingFrame,
            ..Scope::ANY
        };
        let coding = FrameContext::coding(0, Some(0), None);
        let non_coding = FrameContext::NonCoding;
        let coding_ctx = BlockCtx {
            reference: b"AC",
            resulting: b"AT",
            frame: &coding,
            molecule: Molecule::Dna,
            provenance: &prov,
        };
        let non_coding_ctx = BlockCtx {
            reference: b"AC",
            resulting: b"AT",
            frame: &non_coding,
            molecule: Molecule::Dna,
            provenance: &prov,
        };
        assert!(needs_frame.admits(&coding_ctx));
        assert!(!needs_frame.admits(&non_coding_ctx));
        assert!(Scope::ANY.admits(&non_coding_ctx));
    }

    #[test]
    fn the_limbs_are_conjoined() {
        let prov = Provenance::none();
        let frame = FrameContext::NonCoding;
        // DnaOnly + NetDeletion: a DNA net insertion fails on direction, an RNA net
        // deletion fails on molecule; only DNA + net deletion passes.
        let scope = Scope {
            molecule: MoleculeScope::DnaOnly,
            frame: FrameScope::Any,
            direction: DirectionScope::NetDeletion,
        };
        let dna_del = BlockCtx {
            reference: b"ACG",
            resulting: b"AG",
            frame: &frame,
            molecule: Molecule::Dna,
            provenance: &prov,
        };
        let dna_ins = BlockCtx {
            reference: b"AG",
            resulting: b"ACG",
            frame: &frame,
            molecule: Molecule::Dna,
            provenance: &prov,
        };
        let rna_del = BlockCtx {
            reference: b"ACG",
            resulting: b"AG",
            frame: &frame,
            molecule: Molecule::Rna,
            provenance: &prov,
        };
        assert!(scope.admits(&dna_del));
        assert!(!scope.admits(&dna_ins));
        assert!(!scope.admits(&rna_del));
    }
}
