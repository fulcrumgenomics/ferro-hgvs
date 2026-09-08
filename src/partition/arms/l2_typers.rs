//! The Dial-A type-selection arms not already in `arms.rs` (design §5.2, the
//! locked #1–9). `arms.rs` carries #1 (recommended-form ladder), #6 (position-
//! wise) and #8 (always-one-delins control); this module adds the other six:
//!
//! * **#2 `LadderCostTiebreakL2`** — the ladder with a rendered-cost residue
//!   tiebreak (`adjudication-precedence-order` made executable): spec-priority
//!   first, cost only where priority leaves a genuine residue.
//! * **#3 `CoalesceMinimalL2`** — the shipped coalesce-on-minimal back half:
//!   minimal-edit typing, then coalesce each maximal run of consecutive changed
//!   columns into one `delins` (`delins.md:44-47`).
//! * **#4 `MinimalEdit3PrimeL2`** — raw canonical minimal-edit typing, indel
//!   placed 3′-most, no coalescing (the fragmentation-bug contrast).
//! * **#5 `RenderedCostMinL2`** — cost-primary (Mutalyzer-style): the min
//!   rendered-cost typing among the candidates, spec priority ignored.
//! * **#7 `PriorityGreedyL2`** — priority-only greedy, to isolate `general.md:55`
//!   substitution-first ordering: an equal-length block is typed as individual
//!   substitutions (substitution ranks above `delins`), not one spanning `delins`.
//! * **#9 `SpdiCanonicalL2`** — SPDI/sequence-first convention: trim common
//!   flanks, then one `del`/`ins`/`delins`; never `sub`/`dup`/`inv`.
//!
//! Each is a pure function of `(reference, resulting, frame)` and shares the
//! Dial-B merge-back via `apply_dial_b`, so it composes in the grid exactly like
//! the arms in `arms.rs`. Several of these coincide on a pre-trimmed single-edit
//! segment; the design expects that and reads it off the agreement matrix rather
//! than assuming it — they distinguish under a coarse L1 (wide segments).

use crate::partition::arm::{DialBConfig, L1Segmenter, L2Typer, MergeContext, Monolithic};
use crate::partition::arms::arms::LadderL2;
use crate::partition::arms::dial_b::apply_dial_b;
use crate::partition::arms::l1_align::MaximalSplitL1;
use crate::partition::block_ctx::{FrameContext, Molecule, Provenance};
use crate::partition::output::{EditKind, Member, Partition, Segment};

/// Largest changed block the O(n·m) aligner will process; over it the typers fall
/// back to one spanning member. Matches the L1 cutters' cap.
const MAX_BLOCK: usize = 1024;

/// One aligned column of a member-level alignment.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Op {
    Match,
    Sub,
    Del,
    Ins,
}

fn revcomp(b: &[u8]) -> Vec<u8> {
    b.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => *other,
        })
        .collect()
}

/// Rendered-cost of a member list (mirrors `metrics::rendered_cost`): per
/// non-identity member, a base of 4 plus its inserted bytes.
fn rendered_cost(members: &[Member]) -> usize {
    members
        .iter()
        .filter(|m| m.kind != EditKind::Identity)
        .map(|m| 4 + m.inserted.len())
        .sum()
}

/// Needleman–Wunsch (unit substitution and indel costs) over two blocks,
/// returning the column ops. The traceback is fixed and deterministic — the
/// property the confluence gate needs — and biased to place gaps 3′-most by
/// preferring, at equal cost, a diagonal step, then the gap that defers the match
/// (insertion before deletion on the alt-longer side and vice versa is handled by
/// the same preference order). Returns `None` over `MAX_BLOCK`.
fn align_ops(r: &[u8], a: &[u8]) -> Option<Vec<Op>> {
    if r.len() > MAX_BLOCK || a.len() > MAX_BLOCK {
        return None;
    }
    let (n, m) = (r.len(), a.len());
    // cost[i][j] = edit distance of r[i..] vs a[j..], filled from the back so the
    // forward traceback places gaps 3′-most.
    let mut cost = vec![vec![0u32; m + 1]; n + 1];
    for i in (0..n).rev() {
        cost[i][m] = (n - i) as u32;
    }
    for j in (0..m).rev() {
        cost[n][j] = (m - j) as u32;
    }
    for i in (0..n).rev() {
        for j in (0..m).rev() {
            let diag = cost[i + 1][j + 1] + if r[i] == a[j] { 0 } else { 1 };
            let del = cost[i + 1][j] + 1;
            let ins = cost[i][j + 1] + 1;
            cost[i][j] = diag.min(del).min(ins);
        }
    }
    let mut ops = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < n || j < m {
        if i < n && j < m {
            let diag = cost[i + 1][j + 1] + if r[i] == a[j] { 0 } else { 1 };
            let del = cost[i + 1][j] + 1;
            let ins = cost[i][j + 1] + 1;
            let best = diag.min(del).min(ins);
            // Prefer diagonal, then deletion, then insertion — a fixed order.
            if diag == best {
                ops.push(if r[i] == a[j] { Op::Match } else { Op::Sub });
                i += 1;
                j += 1;
            } else if del == best {
                ops.push(Op::Del);
                i += 1;
            } else {
                ops.push(Op::Ins);
                j += 1;
            }
        } else if i < n {
            ops.push(Op::Del);
            i += 1;
        } else {
            ops.push(Op::Ins);
            j += 1;
        }
    }
    Some(ops)
}

fn shared_prefix(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}
fn shared_suffix(a: &[u8], b: &[u8], prefix: usize) -> usize {
    let max = a.len().min(b.len()) - prefix;
    a.iter()
        .rev()
        .zip(b.iter().rev())
        .take(max)
        .take_while(|(x, y)| x == y)
        .count()
}

/// One spanning member over `[rs, re)` with `alt` — `Del` if `alt` empty,
/// `Identity` if both empty, else `Delins`.
fn spanning(rs: usize, re: usize, alt: Vec<u8>) -> Member {
    let kind = if alt.is_empty() {
        if re == rs {
            EditKind::Identity
        } else {
            EditKind::Del
        }
    } else {
        EditKind::Delins
    };
    Member {
        kind,
        ref_start: rs,
        ref_end: re,
        inserted: alt,
    }
}

/// Reference / alt block bytes and the reference start for a segment. Typing
/// reads alt bytes via the block-relative index (`a[ai]`), so no res-side
/// absolute coordinate is needed.
struct Blocks<'a> {
    r: &'a [u8],
    a: &'a [u8],
    ref_start: usize,
}

fn blocks<'a>(seg: &Segment, reference: &'a [u8], resulting: &'a [u8]) -> Blocks<'a> {
    Blocks {
        r: &reference[seg.ref_start..seg.ref_end],
        a: &resulting[seg.res_start..seg.res_end],
        ref_start: seg.ref_start,
    }
}

/// Type a segment from aligned ops, coalescing per `coalesce`:
/// * `false` (#4): subs individual, `Del`/`Ins` runs coalesced, matches skipped.
/// * `true`  (#3): each maximal run of consecutive changed columns → one `delins`.
fn type_from_ops(b: &Blocks, ops: &[Op], coalesce: bool) -> Vec<Member> {
    let mut members = Vec::new();
    let (mut ri, mut ai) = (0usize, 0usize);
    let mut k = 0usize;
    while k < ops.len() {
        match ops[k] {
            Op::Match => {
                ri += 1;
                ai += 1;
                k += 1;
            }
            _ if coalesce => {
                // Maximal run of non-Match ops → one delins.
                let run_rs = ri;
                let mut alt = Vec::new();
                while k < ops.len() && ops[k] != Op::Match {
                    match ops[k] {
                        Op::Sub => {
                            alt.push(b.a[ai]);
                            ri += 1;
                            ai += 1;
                        }
                        Op::Del => ri += 1,
                        Op::Ins => {
                            alt.push(b.a[ai]);
                            ai += 1;
                        }
                        Op::Match => unreachable!(),
                    }
                    k += 1;
                }
                members.push(spanning(b.ref_start + run_rs, b.ref_start + ri, alt));
            }
            Op::Sub => {
                members.push(Member {
                    kind: EditKind::Sub,
                    ref_start: b.ref_start + ri,
                    ref_end: b.ref_start + ri + 1,
                    inserted: vec![b.a[ai]],
                });
                ri += 1;
                ai += 1;
                k += 1;
            }
            Op::Del => {
                let start = ri;
                while k < ops.len() && ops[k] == Op::Del {
                    ri += 1;
                    k += 1;
                }
                members.push(Member {
                    kind: EditKind::Del,
                    ref_start: b.ref_start + start,
                    ref_end: b.ref_start + ri,
                    inserted: Vec::new(),
                });
            }
            Op::Ins => {
                let mut alt = Vec::new();
                while k < ops.len() && ops[k] == Op::Ins {
                    alt.push(b.a[ai]);
                    ai += 1;
                    k += 1;
                }
                members.push(Member {
                    kind: EditKind::Ins,
                    ref_start: b.ref_start + ri,
                    ref_end: b.ref_start + ri,
                    inserted: alt,
                });
            }
        }
    }
    members
}

/// Whole-segment spanning delins (the over-`MAX_BLOCK` fallback).
fn whole_segment(seg: &Segment, resulting: &[u8]) -> Vec<Member> {
    vec![spanning(
        seg.ref_start,
        seg.ref_end,
        resulting[seg.res_start..seg.res_end].to_vec(),
    )]
}

fn merge_back_dial_b(
    m: Vec<Member>,
    reference: &[u8],
    resulting: &[u8],
    ctx: &MergeContext,
    cfg: &DialBConfig,
) -> Vec<Member> {
    apply_dial_b(m, reference, resulting, ctx, cfg)
}

// ---- #4 minimal-edit, 3′-most, no coalesce --------------------------------

/// Dial A #4: raw minimal-edit typing with the indel placed 3′-most. Subs are
/// individual and `Del`/`Ins` runs are coalesced, but changed columns separated
/// by a match are NOT joined — the fragmentation-bug contrast.
pub struct MinimalEdit3PrimeL2;
impl L2Typer for MinimalEdit3PrimeL2 {
    fn name(&self) -> &str {
        "minimal-edit-3prime"
    }
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        _f: &FrameContext,
    ) -> Vec<Member> {
        let b = blocks(seg, reference, resulting);
        match align_ops(b.r, b.a) {
            Some(ops) => type_from_ops(&b, &ops, false),
            None => whole_segment(seg, resulting),
        }
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member> {
        merge_back_dial_b(m, reference, resulting, ctx, cfg)
    }
}

// ---- #3 coalesce-on-minimal -----------------------------------------------

/// Dial A #3: minimal-edit typing, then coalesce each maximal run of consecutive
/// changed columns into one `delins` (the shipped coalesce back half).
pub struct CoalesceMinimalL2;
impl L2Typer for CoalesceMinimalL2 {
    fn name(&self) -> &str {
        "coalesce-minimal"
    }
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        _f: &FrameContext,
    ) -> Vec<Member> {
        let b = blocks(seg, reference, resulting);
        match align_ops(b.r, b.a) {
            Some(ops) => type_from_ops(&b, &ops, true),
            None => whole_segment(seg, resulting),
        }
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member> {
        merge_back_dial_b(m, reference, resulting, ctx, cfg)
    }
}

// ---- #9 SPDI/sequence-first canonical typing ------------------------------

/// Dial A #9: the SPDI convention — trim common flanks, emit one `del`/`ins`/
/// `delins` over the residual. Never types `sub`/`dup`/`inv` (a convention
/// cross-check against the spec-shaped arms).
pub struct SpdiCanonicalL2;
impl L2Typer for SpdiCanonicalL2 {
    fn name(&self) -> &str {
        "spdi-canonical"
    }
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        _f: &FrameContext,
    ) -> Vec<Member> {
        let b = blocks(seg, reference, resulting);
        let p = shared_prefix(b.r, b.a);
        let s = shared_suffix(b.r, b.a, p);
        let rs = seg.ref_start + p;
        let re = seg.ref_end - s;
        let alt = resulting[seg.res_start + p..seg.res_end - s].to_vec();
        if rs == re && alt.is_empty() {
            return vec![Member {
                kind: EditKind::Identity,
                ref_start: rs,
                ref_end: re,
                inserted: Vec::new(),
            }];
        }
        vec![spanning(rs, re, alt)]
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member> {
        merge_back_dial_b(m, reference, resulting, ctx, cfg)
    }
}

// ---- #7 priority-only greedy ----------------------------------------------

/// Dial A #7: priority-only. Types the whole segment by the highest-priority
/// operator that fits, to isolate `general.md:55` substitution-first ordering:
/// an equal-length block becomes individual **substitutions** (substitution ranks
/// above `delins`), where the ladder (#1) would emit one spanning `delins`.
pub struct PriorityGreedyL2;
impl L2Typer for PriorityGreedyL2 {
    fn name(&self) -> &str {
        "priority-greedy"
    }
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        _f: &FrameContext,
    ) -> Vec<Member> {
        let b = blocks(seg, reference, resulting);
        let ref_len = b.r.len();
        let alt = b.a.to_vec();
        if ref_len == 0 && alt.is_empty() {
            return vec![Member {
                kind: EditKind::Identity,
                ref_start: seg.ref_start,
                ref_end: seg.ref_end,
                inserted: Vec::new(),
            }];
        }
        // dup (>=2, copies immediate 5' reference) ranks first.
        if ref_len == 0 {
            let n = alt.len();
            if n >= 2
                && seg.ref_start >= n
                && reference[seg.ref_start - n..seg.ref_start] == alt[..]
            {
                return vec![Member {
                    kind: EditKind::Dup,
                    ref_start: seg.ref_start,
                    ref_end: seg.ref_start,
                    inserted: alt,
                }];
            }
            return vec![Member {
                kind: EditKind::Ins,
                ref_start: seg.ref_start,
                ref_end: seg.ref_start,
                inserted: alt,
            }];
        }
        if alt.is_empty() {
            return vec![Member {
                kind: EditKind::Del,
                ref_start: seg.ref_start,
                ref_end: seg.ref_end,
                inserted: Vec::new(),
            }];
        }
        // inv over the whole span (>=2) ranks above delins.
        if alt.len() == ref_len && ref_len >= 2 && alt == revcomp(b.r) {
            return vec![Member {
                kind: EditKind::Inv,
                ref_start: seg.ref_start,
                ref_end: seg.ref_end,
                inserted: Vec::new(),
            }];
        }
        // Substitution-first: an equal-length block is individual subs.
        if alt.len() == ref_len {
            let mut members = Vec::new();
            for (i, (rb, ab)) in b.r.iter().zip(&alt).enumerate() {
                if rb != ab {
                    members.push(Member {
                        kind: EditKind::Sub,
                        ref_start: seg.ref_start + i,
                        ref_end: seg.ref_start + i + 1,
                        inserted: vec![*ab],
                    });
                }
            }
            return members;
        }
        // Unequal length with no higher-priority whole-span type → one delins.
        vec![spanning(seg.ref_start, seg.ref_end, alt)]
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member> {
        merge_back_dial_b(m, reference, resulting, ctx, cfg)
    }
}

// ---- #10 structural typer (position-wise split for equal length) ----------

/// Dial A #10: type a changed block WITHOUT an internal edit-alignment, then type
/// each piece by its structure on the recommended ladder.
///
/// The rule is the `unchanged-is-read-over-every-minimal-alignment` / #2174 reading:
/// * A whole-span exact reverse complement (>=2) is one `inv`, checked FIRST so an
///   interior coincidence cannot split it (`inversion.md:5`).
/// * An **equal-length** block is split ONLY at columns unchanged at their own
///   coordinate (position-wise fixed points); each maximal run of position-wise-changed
///   columns is typed structurally. A minimal-edit alignment's interior matches are
///   shift artifacts, not separators, so a contiguous equal-length change (`CAG->AGA`)
///   is one member, never del+ins.
/// * An **unequal-length** block (net indel) is one unit — there is no position-wise
///   correspondence to read a separator from, so a payload coincidence cannot manufacture
///   a cut (`delins.md:47`, #1610).
///
/// Each piece is typed: 1-base change -> `sub`; pure deletion -> `del`; pure insertion
/// copying its immediate 5′ reference -> `dup` (>=2, `duplication.md:18`) else `ins`;
/// equal-length exact revcomp (>=2) -> `inv`; else `delins`. It never spells consecutive
/// changes as individual subs (`substitution.md:15`). The earlier revision aligned with
/// Needleman–Wunsch and coalesced its non-match runs, which fragmented contiguous changes
/// exactly the way the over-splitting typers do (measured on the hard-forms grader); this
/// revision removes the alignment.
pub struct CoalesceStructuralL2;

impl CoalesceStructuralL2 {
    /// Structurally type one maximal changed run over absolute ref `[rs, re)` yielding
    /// `alt`. `reference` is needed only for the dup 5′-copy and inv reverse-complement
    /// tests. Mirrors the member conventions of #7 (`Dup`/`Ins` zero-width with the
    /// payload in `inserted`; `Inv` over the span with empty `inserted`).
    fn type_run(reference: &[u8], rs: usize, re: usize, alt: Vec<u8>) -> Member {
        let span_len = re - rs;
        if span_len == 0 {
            let n = alt.len();
            // n >= 1: a single-nucleotide insertion copying its immediate 5′ reference is a
            // `dup` (`duplication.md:5` "one or more", the `c.20dup` example forbidding
            // `c.19_20insT`, and the `:18` MUST). Operator ruling 2026-08-22 (Side A).
            if n >= 1 && rs >= n && reference[rs - n..rs] == alt[..] {
                return Member {
                    kind: EditKind::Dup,
                    ref_start: rs,
                    ref_end: rs,
                    inserted: alt,
                };
            }
            return Member {
                kind: EditKind::Ins,
                ref_start: rs,
                ref_end: rs,
                inserted: alt,
            };
        }
        if alt.is_empty() {
            return Member {
                kind: EditKind::Del,
                ref_start: rs,
                ref_end: re,
                inserted: Vec::new(),
            };
        }
        if span_len == 1 && alt.len() == 1 {
            return Member {
                kind: EditKind::Sub,
                ref_start: rs,
                ref_end: re,
                inserted: alt,
            };
        }
        if alt.len() == span_len && span_len >= 2 && alt == revcomp(&reference[rs..re]) {
            return Member {
                kind: EditKind::Inv,
                ref_start: rs,
                ref_end: re,
                inserted: Vec::new(),
            };
        }
        spanning(rs, re, alt)
    }
}

impl L2Typer for CoalesceStructuralL2 {
    fn name(&self) -> &str {
        "coalesce-structural"
    }
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        _f: &FrameContext,
    ) -> Vec<Member> {
        let b = blocks(seg, reference, resulting);
        if b.r.is_empty() && b.a.is_empty() {
            return Vec::new();
        }
        // Whole-segment inversion FIRST: an exact reverse complement must be typed `inv`
        // before any split, so an interior coincidence (which the position-wise pass
        // below would treat as a wall) cannot shatter it (`inversion.md:5`,
        // `whole-span-reverse-complement-types-as-inv`).
        if b.r.len() == b.a.len() && b.r.len() >= 2 && b.a == revcomp(b.r) {
            return vec![Member {
                kind: EditKind::Inv,
                ref_start: seg.ref_start,
                ref_end: seg.ref_end,
                inserted: Vec::new(),
            }];
        }
        if b.r.len() != b.a.len() {
            // Net indel: no position-wise correspondence, so one structural unit — a
            // payload coincidence cannot manufacture a separator (#1610).
            return vec![Self::type_run(
                reference,
                seg.ref_start,
                seg.ref_end,
                b.a.to_vec(),
            )];
        }
        // Equal length: split ONLY at position-wise fixed points; type each run.
        let mut members = Vec::new();
        let mut i = 0usize;
        while i < b.r.len() {
            if b.r[i] == b.a[i] {
                i += 1;
                continue;
            }
            let s = i;
            while i < b.r.len() && b.r[i] != b.a[i] {
                i += 1;
            }
            members.push(Self::type_run(
                reference,
                seg.ref_start + s,
                seg.ref_start + i,
                b.a[s..i].to_vec(),
            ));
        }
        members
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member> {
        merge_back_dial_b(m, reference, resulting, ctx, cfg)
    }
}

// ---- #2 ladder + rendered-cost tiebreak -----------------------------------

/// Dial A #2: the recommended-form ladder with a rendered-cost residue tiebreak.
/// Spec priority decides first (the ladder's single-member type); where the ladder
/// and the substitution-first reading are BOTH spec-legal (an equal-length block),
/// the lower rendered-cost form wins. This is the arm the design flags as the
/// likely real target — `adjudication-precedence-order` with cost as the residue.
pub struct LadderCostTiebreakL2;
impl L2Typer for LadderCostTiebreakL2 {
    fn name(&self) -> &str {
        "ladder-cost-tiebreak"
    }
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        f: &FrameContext,
    ) -> Vec<Member> {
        let ladder = LadderL2.type_segment(seg, reference, resulting, f);
        // Only an equal-length multi-mismatch block has a spec-legal alternative
        // (individual subs). Elsewhere the ladder's choice stands.
        let b = blocks(seg, reference, resulting);
        if b.r.len() == b.a.len() && b.r.len() >= 2 && b.r != b.a {
            let subs = PriorityGreedyL2.type_segment(seg, reference, resulting, f);
            if rendered_cost(&subs) < rendered_cost(&ladder) {
                return subs;
            }
        }
        ladder
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member> {
        merge_back_dial_b(m, reference, resulting, ctx, cfg)
    }
}

// ---- #5 rendered-cost minimization ----------------------------------------

/// Dial A #5: cost-primary (Mutalyzer-style). Among the candidate typings —
/// spanning delins (#9-shape), coalesce-on-minimal (#3), minimal-edit (#4) and
/// the ladder (#1) — pick the one with the least rendered cost, spec priority
/// ignored. Ties break toward fewer members, then a fixed candidate order.
pub struct RenderedCostMinL2;
impl L2Typer for RenderedCostMinL2 {
    fn name(&self) -> &str {
        "rendered-cost-min"
    }
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        f: &FrameContext,
    ) -> Vec<Member> {
        let candidates = vec![
            SpdiCanonicalL2.type_segment(seg, reference, resulting, f),
            CoalesceMinimalL2.type_segment(seg, reference, resulting, f),
            MinimalEdit3PrimeL2.type_segment(seg, reference, resulting, f),
            LadderL2.type_segment(seg, reference, resulting, f),
        ];
        candidates
            .into_iter()
            .min_by_key(|c| (rendered_cost(c), c.len()))
            .unwrap_or_else(|| whole_segment(seg, resulting))
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member> {
        merge_back_dial_b(m, reference, resulting, ctx, cfg)
    }
}

/// Merge maximal runs of atomic segments that are contiguous in BOTH reference
/// and resulting (no unchanged base between them) into one segment. The breaks are
/// exactly the unchanged bases — the separation-rule walls.
fn coalesce_atomic_segments(segs: &[Segment]) -> Vec<Segment> {
    let mut out: Vec<Segment> = Vec::new();
    for seg in segs {
        match out.last_mut() {
            Some(prev) if prev.ref_end == seg.ref_start && prev.res_end == seg.res_start => {
                prev.ref_end = seg.ref_end;
                prev.res_end = seg.res_end;
            }
            _ => out.push(seg.clone()),
        }
    }
    out
}

/// The merge-only canonicaliser (design experiment): start from the MAXIMAL split
/// (atomic edits), coalesce each maximal changed run, ladder-type it, then apply
/// the Dial-B carve-outs. It never splits — every step is a merge — so it cannot
/// hit the split-then-remerge failure class of ferro's `merge.rs`. Monolithic
/// because it re-derives the changed runs from the sequences, which makes it
/// L1-invariant by construction: the whole point of the architecture is that the
/// L1 dial disappears and the partitioner IS the merge policy.
pub struct MaxSplitMerge;
impl Monolithic for MaxSplitMerge {
    fn name(&self) -> &str {
        "max-split-merge"
    }
    // Bakeoff-experimental arm; not a production canonical-coalesce rule. Its
    // `cuts_with_canonical` value is not asserted by any test (Task 2 R2).
    fn cuts_with_canonical(&self) -> bool {
        false
    }
    fn partition(&self, reference: &[u8], resulting: &[u8], frame: &FrameContext) -> Partition {
        // `MaximalSplitL1` is axis-blind (finest sound split), so the molecule is
        // immaterial here; this monolithic's whole premise is that L1 disappears.
        let atoms =
            MaximalSplitL1.segment(reference, resulting, Molecule::Dna, &Provenance::none());
        let runs = coalesce_atomic_segments(&atoms);
        let ladder = LadderL2;
        let mut members = Vec::new();
        for seg in &runs {
            members.extend(ladder.type_segment(seg, reference, resulting, frame));
        }
        // The `Monolithic` trait carries no molecule/provenance; Dial-B used to be
        // molecule-blind here, which is operationally the DNA reading. Pin it so
        // this arm's output is byte-identical across the threading (residual:
        // monolithic arms stay molecule-blind — see the module note in `arm.rs`).
        let none = Provenance::none();
        let ctx = MergeContext {
            frame,
            molecule: Molecule::Dna,
            provenance: &none,
        };
        let members = apply_dial_b(
            members,
            reference,
            resulting,
            &ctx,
            &DialBConfig::ledger_current(),
        );
        Partition { members }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::metrics::{RefApplier, SequenceApplier};

    fn seg(rs: usize, re: usize, ss: usize, se: usize) -> Segment {
        Segment {
            ref_start: rs,
            ref_end: re,
            res_start: ss,
            res_end: se,
        }
    }

    // Round-trip a typer's output over a whole-span segment.
    fn roundtrip(t: &dyn L2Typer, reference: &[u8], resulting: &[u8]) -> Vec<Member> {
        let s = seg(0, reference.len(), 0, resulting.len());
        let members = t.type_segment(&s, reference, resulting, &FrameContext::NonCoding);
        let p = Partition {
            members: members.clone(),
        };
        assert_eq!(
            RefApplier.apply(reference, &p).as_deref(),
            Some(resulting),
            "typer {} must round-trip",
            t.name()
        );
        members
    }

    #[test]
    fn all_six_typers_round_trip() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"ACGTACGT", b"ATGTACGT"),     // one sub
            (b"ACGTACGT", b"ACGTACG"),      // trailing del
            (b"ACGTACGT", b"ACGTTACGT"),    // insertion
            (b"AAAAAA", b"CAGAAA"),         // two subs (equal length)
            (b"CGCG", b"AAC"),              // net-deletion delins (#1610 shape)
            (b"ACGT", b"TGCA"),             // whole-span revcomp (inv candidate)
            (b"ACGTACGTAC", b"ACgtACGTAC"), // no-op-ish (lowercase won't match; treated as sub)
        ];
        let typers: Vec<Box<dyn L2Typer>> = vec![
            Box::new(MinimalEdit3PrimeL2),
            Box::new(CoalesceMinimalL2),
            Box::new(SpdiCanonicalL2),
            Box::new(PriorityGreedyL2),
            Box::new(LadderCostTiebreakL2),
            Box::new(RenderedCostMinL2),
        ];
        for (r, a) in cases {
            for t in &typers {
                roundtrip(t.as_ref(), r, a);
            }
        }
    }

    #[test]
    fn priority_greedy_types_equal_length_as_subs() {
        // AAAAAA -> CAGAAA: two subs at 0 and 2. #7 emits two Subs; the ladder (#1)
        // would emit one delins over [0,3).
        let m = roundtrip(&PriorityGreedyL2, b"AAAAAA", b"CAGAAA");
        assert!(m.iter().all(|x| x.kind == EditKind::Sub));
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn spdi_never_types_sub_or_inv() {
        // Whole-span revcomp ACGT->TGCA: SPDI trims nothing (no common flank) and
        // emits one delins, never inv.
        let m = roundtrip(&SpdiCanonicalL2, b"ACGT", b"TGCA");
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].kind, EditKind::Delins);
    }

    #[test]
    fn cost_min_beats_or_ties_the_ladder() {
        // The min-cost candidate can never cost more than the ladder alone.
        let s = seg(0, 6, 0, 6);
        let ladder = LadderL2.type_segment(&s, b"AAAAAA", b"CAGAAA", &FrameContext::NonCoding);
        let costmin =
            RenderedCostMinL2.type_segment(&s, b"AAAAAA", b"CAGAAA", &FrameContext::NonCoding);
        assert!(rendered_cost(&costmin) <= rendered_cost(&ladder));
    }
}
