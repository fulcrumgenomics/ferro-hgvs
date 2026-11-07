//! `OperatorExtractL2` — the Opus-5 round-1 independent design for the peel frontier.
//!
//! Thesis (deliberately the OPPOSITE architectural bet from `anchor_peel.rs`, which
//! puts the peel in L1): the member-beside-a-member peel is a **TYPER (Dial-A)
//! decision, not an L1 wall decision.** An L1 wall is "where is unchanged content"
//! (a molecule-only, ruling-grounded fact); a dup/inv is "what KIND is this piece"
//! (a typing fact). Keeping the two apart lets the typer fall back to the
//! coalesce-structural (#2174) reading whenever no structured member is significant
//! — so it inherits B's measured 100% on Delins / Inv / multi-Sub instead of
//! fragmenting them, which is the price the L1-located peel pays (measured: peel L1
//! regresses the Delins bucket 100% -> ~42%).
//!
//! Given ONE block from a coarse L1 (`trim-only`), `decompose` runs a fixed
//! structural precedence, peeling a structured member and recursing on the residual:
//!
//!   0. whole-span revcomp (>=2)            -> one Inv  (`whole-span-reverse-complement-types-as-inv`)
//!   1. INV anchor (revcomp run)            -> peel Inv, recurse both flanks, iff a
//!      #2155 FLANK GUARD passes (>=1 unchanged/edge flank); cost then breaks the
//!      remaining peel-vs-no-peel tie (Sub+Inv, Del+Inv, Ins+Inv)
//!   2. pure insertion (R==0)               -> Dup if copies 5' flank else Ins
//!   3. pure deletion (A==0)                -> Del
//!   4. single sub (R==A==1)                -> Sub
//!   5. DUP anchor on a net insertion       -> peel Dup (5' or 3' copy), recurse
//!      (dup-beside-a-member, #2175 / #2194)
//!   6. equal-length block                  -> split at position-wise fixed points,
//!      type each run (sub / inv / delins) — the coalesce-structural reading
//!   7. STRUCTURAL LEAF                      -> split at the forced-unchanged walls,
//!      one member per segment, cost NEVER consulted (`general.md:33`)
//!
//! The design is `adjudication-precedence-order` as architecture: SPEC PRIORITY FIRST,
//! cost ONLY on the genuine residue. Every gate is a structural (ruling-grounded) test;
//! cost is the LAST tiebreak among spec-legal spellings and is confined to the INV
//! anchor's peel-vs-no-peel choice — the leaf is fully structural.
//!   * INV anchor: peel is ELIGIBLE only with a clean (unchanged/edge) flank — a
//!     revcomp changed on both sides is a payload coincidence that collapses into the
//!     spanning delins, a structural rule cost cannot express (a structured member
//!     spells 0 bases, so cost always wants to peel). Among eligible readings, cost
//!     (weighted when `weighted`) breaks the tie: a REAL inversion NAMES A RANGE and
//!     spells zero literal bases, so it is cheaper than the delins spelling the same
//!     bases; a chance revcomp is not. A MAJORITY-IDENTITY veto (M1) additionally
//!     refuses a peel whose span is >= half positionally unchanged — a chance seed that
//!     annexed unchanged bases to a block edge, not an inversion. This recovers mdl's
//!     inv-neighbour strength (Del+Inv ~93, Ins+Inv ~92, Sub+Inv ~96 on the 17k corpus)
//!     without mdl's collapse of multi-member forms.
//!   * DUP anchor: an EXACT tandem copy (k copies of ONE unit) of the adjacent reference
//!     flank (`duplication-must-ranks-the-label-not-the-partition`).
//!
//! Molecule-only by construction: nothing reads the frame. The coding C1 codon merge
//! stays downstream in Dial-B (`apply_dial_b`). The molecule-dependent collapse-vs-split
//! of a payload-coincidence delins is deliberately NOT decided here — it is a downstream
//! DNA/RNA axis layer's job (`delins-payload-coincidence-carve-out-is-coding-dna-
//! scoped`). Which way that layer moves depends on the block: an UNEQUAL-length leaf
//! (step 7) is split here and the DNA axis would COLLAPSE it; an EQUAL-length run
//! (step 6) is emitted as one delins here and the RNA axis would SPLIT it. This typer
//! makes the molecule-blind structural call and leaves that adjustment to the axis.

use crate::partition::arm::{DialBConfig, L1Segmenter, L2Typer, MergeContext};
use crate::partition::arms::arms::AllAlignmentSplitL1;
use crate::partition::arms::dial_b::apply_dial_b;
use crate::partition::block_ctx::{FrameContext, Molecule, Provenance};
use crate::partition::output::{EditKind, Member, Segment};

/// Largest block the O(n*m) LCS will process; over it, one spanning member.
const MAX_BLOCK: usize = 1024;
const MAX_DEPTH: usize = 64;

fn revcomp(b: &[u8]) -> Vec<u8> {
    b.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b'u',
            b'u' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => *other,
        })
        .collect()
}

/// Watson-Crick complement of a single base (DNA + RNA), identity otherwise.
fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b'u',
        b'u' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

/// Longest common substring of `r` and `a` as `(r_off, a_off, len)`, deterministic
/// (first maximal in scan order), `None` if nothing matches.
fn lcs_substring(r: &[u8], a: &[u8]) -> Option<(usize, usize, usize)> {
    if r.is_empty() || a.is_empty() {
        return None;
    }
    let mut best = (0usize, 0usize, 0usize); // (len, r_off, a_off)
    let mut prev = vec![0usize; a.len() + 1];
    let mut cur = vec![0usize; a.len() + 1];
    for i in 1..=r.len() {
        cur[0] = 0;
        for j in 1..=a.len() {
            cur[j] = if r[i - 1] == a[j - 1] {
                prev[j - 1] + 1
            } else {
                0
            };
            if cur[j] > best.0 {
                best = (cur[j], i - cur[j], j - cur[j]);
            }
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    if best.0 == 0 {
        None
    } else {
        Some((best.1, best.2, best.0))
    }
}

/// The peel typer. Tunable: `inv_min` (min inv-anchor seed span), `dup_min` (min
/// tandem-copy length for a peeled dup).
pub struct OperatorExtractL2 {
    /// Per-member structural cost `K` in the rendered-cost peel test.
    member_cost: usize,
    inv_min: usize,
    dup_min: usize,
    /// When true, the net-insertion dup peel uses the "wholly-consumed" completion
    /// condition (multi-unit, 1-base allowed) grafted from `anchor_peel` — a net
    /// insertion is peeled as dup copies only if the copies consume ALL of it,
    /// leaving an equal-length residual. This reaches the 1-base tandem dup beside a
    /// member that `dup_min = 2` misses (measured: it is why the length-floored peel
    /// sat at Sub+Dup 52 vs the completion-condition peel's 88), while the completion
    /// condition itself stops a chance 1-base flank match from fragmenting an `ins`.
    dup_complete: bool,
    /// When true, the rendered-cost peel test charges each SPELLED base -log2(p) for
    /// its frequency p in the local window, instead of a flat 1 (step 2b). This gives
    /// the structural inv gate teeth: under flat cost an inv (0 spelled bases) beats
    /// almost any spelling of the same span, so a coincidental short revcomp in
    /// low-complexity (A/T-rich) content is peeled (the octamer break — measured: at
    /// K=2 a chance 3-base GTT->AAC inv undercuts one delins). Weighting makes those
    /// A/T bases cheap to spell, so the chance inv no longer wins — collapse becomes
    /// K-robust and the K-Pareto dissolves.
    weighted: bool,
    name: String,
}

impl OperatorExtractL2 {
    pub fn new(member_cost: usize, inv_min: usize, dup_min: usize, name: &str) -> Self {
        // Minor 3: `dup_min = 0` would let the floored peel match a zero-length copy
        // (`[] == []`) and peel empty dups until MAX_DEPTH. No caller passes 0; guard it.
        debug_assert!(dup_min >= 1, "dup_min must be >= 1");
        Self {
            member_cost,
            inv_min,
            dup_min,
            dup_complete: false,
            weighted: false,
            name: name.to_string(),
        }
    }

    /// As [`Self::new`] but with the wholly-consumed tandem dup peel. `dup_min` is set
    /// to 1 so a genuine 1-base tandem is reachable via the completion branch; the
    /// non-completing fall-through floors itself at 2 (see M2 in `decompose`), so this
    /// 1 does NOT admit a chance 1-base single peel.
    pub fn v2(member_cost: usize, inv_min: usize, name: &str) -> Self {
        Self {
            member_cost,
            inv_min,
            dup_min: 1,
            dup_complete: true,
            weighted: false,
            name: name.to_string(),
        }
    }

    /// The step-2b synthesis: v2's dup peel + information-weighted cost.
    pub fn weighted(member_cost: usize, inv_min: usize, name: &str) -> Self {
        Self {
            member_cost,
            inv_min,
            dup_min: 1,
            dup_complete: true,
            weighted: true,
            name: name.to_string(),
        }
    }

    /// The peel-vs-leave cost of `members` over the composition window `ctx`.
    /// Flat `rendered_cost` unless `weighted`, then information-weighted.
    fn score(&self, members: &[Member], ctx: &[u8]) -> f64 {
        if self.weighted {
            weighted_cost(members, self.member_cost, ctx)
        } else {
            rendered_cost(members, self.member_cost) as f64
        }
    }
}

/// Information-weighted rendered cost: like [`rendered_cost`], but each SPELLED base
/// costs `-log2(p)` for its frequency `p` in the window `ctx` (Laplace-smoothed),
/// instead of a flat 1. The per-member token `k` stays flat. In a low-complexity
/// window a spelled base is cheap, so a coincidental short revcomp/tandem no longer
/// undercuts spelling the bases — which is what gives the structural inv gate teeth.
fn weighted_cost(members: &[Member], k: usize, ctx: &[u8]) -> f64 {
    let idx = |b: u8| match b {
        b'A' | b'a' => 0usize,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        _ => 3,
    };
    let mut counts = [1.0f64; 4]; // Laplace pseudocounts, so an absent base is finite
    for &b in ctx {
        counts[idx(b)] += 1.0;
    }
    let total: f64 = counts.iter().sum();
    let base_bits = |b: u8| -(counts[idx(b)] / total).log2();
    members
        .iter()
        .filter(|m| m.kind != EditKind::Identity)
        .map(|m| {
            let spelled: f64 = match m.kind {
                EditKind::Sub | EditKind::Ins | EditKind::Delins => {
                    m.inserted.iter().map(|&b| base_bits(b)).sum()
                }
                _ => 0.0,
            };
            k as f64 + spelled
        })
        .sum()
}

/// A pure-insertion (or dup) leaf at `rs`: Dup if the whole insert copies the
/// immediately-5' reference (>=1 base), else Ins. Mirrors `LadderL2`/`type_run`.
fn insertion_leaf(reference: &[u8], rs: usize, ins: &[u8], out: &mut Vec<Member>) {
    let n = ins.len();
    if n == 0 {
        return;
    }
    if rs >= n && reference[rs - n..rs] == *ins {
        out.push(Member {
            kind: EditKind::Dup,
            ref_start: rs,
            ref_end: rs,
            inserted: ins.to_vec(),
        });
    } else {
        out.push(Member {
            kind: EditKind::Ins,
            ref_start: rs,
            ref_end: rs,
            inserted: ins.to_vec(),
        });
    }
}

/// The #2175 / P2 "one reference-tandem dup beside one substitution" reading of a PURE
/// INSERTION of `ins` at reference offset `rs` (`r_len == 0`). Reclaims the single byte
/// `b = reference[rs-1]` the shared-flank trim consumed (it equals `resulting[ss-1]` by
/// construction) and asks whether `[b] ++ ins` is `k` whole copies of an `l`-byte
/// (`l >= 2`) tandem unit that copies the reference motif immediately 5' of `rs-1`, with
/// exactly the last inserted byte left over as a substitution (`ins[n-1] != b`). Returns
/// `(unit, copies)` for the dup, or `None` to fall through to a flat `insertion_leaf`.
///
/// Grounds: `duplication-must-ranks-the-label-not-the-partition` (#2175 extension) + P2
/// (the preference is spelling-independent, so it applies even though a flat `ins`
/// spelling also exists). Anti-over-peel is categorical, never a probability floor
/// (`inverted-duplication-is-derived-as-ins-range-inv` rejects a flat length minimum):
/// `l >= 2` is #2175's own discriminator (a 1-base anchor is the #2174 solid-run
/// collapse, not a dup); the reclaim depth is fixed at exactly one byte, encoding
/// #2175's "exactly one substitution" as a hard structural cap; and the wholly-consumed
/// periodicity test mirrors the `dup_complete` completion discipline. Finding #7
/// (bakeoff `c7-g-2175`).
fn dup_beside_sub_unit(reference: &[u8], rs: usize, ins: &[u8]) -> Option<(Vec<u8>, usize)> {
    let n = ins.len();
    if n == 0 || rs < 1 {
        return None;
    }
    let b = reference[rs - 1];
    // The genuine substitution: the last inserted byte must DIFFER from the reclaimed
    // reference byte. If it matched, this is a plain tandem, not a dup-beside-sub.
    if ins[n - 1] == b {
        return None;
    }
    // motif = [b] ++ ins (length n + 1): k copies of an l-byte unit, then the sub byte.
    let mut motif = Vec::with_capacity(n + 1);
    motif.push(b);
    motif.extend_from_slice(ins);
    // Smallest period l >= 2 that tiles n, whose l-byte unit copies the reference motif
    // immediately 5' of `rs-1`, and under which motif[0..n] is l-periodic.
    for l in 2..=n {
        if !n.is_multiple_of(l) || rs - 1 < l {
            continue;
        }
        let unit = &reference[rs - 1 - l..rs - 1];
        if &motif[..l] != unit {
            continue;
        }
        if (0..n).all(|j| motif[j] == motif[j % l]) {
            return Some((unit.to_vec(), n / l));
        }
    }
    None
}

/// Rendered description cost: each non-identity member costs `k` (a position range:
/// `pos` or `pos_pos`) plus the LITERAL bases it must spell. `inv`, `del`, and `dup`
/// name a range and spell zero bases (a `dup`'s payload is a back-reference, not
/// literal); `sub`/`ins`/`delins` spell their inserted bases. This is what makes a
/// real inversion or duplication cheaper than the delins that would spell the same
/// bases, and a chance revcomp not. `k` is the peel-aggressiveness knob: a larger `k`
/// penalises extra members and so peels more conservatively.
fn rendered_cost(members: &[Member], k: usize) -> usize {
    members
        .iter()
        .filter(|m| m.kind != EditKind::Identity)
        .map(|m| {
            k + match m.kind {
                EditKind::Sub | EditKind::Ins | EditKind::Delins => m.inserted.len(),
                EditKind::Del | EditKind::Inv | EditKind::Dup | EditKind::Identity => 0,
            }
        })
        .sum()
}

/// Structurally type one maximal changed run over absolute ref `[rs, re)` -> `alt`
/// (the coalesce-structural leaf). `sub`/`del`/`ins`/`dup`/`inv`/`delins`.
/// Find a cheap wall to split an over-`MAX_BLOCK` block on: the midpoint of the
/// longest run of position-wise-identical bases, over the shorter of the two lengths.
///
/// Returns an offset `mid` (`1 <= mid < min(r.len(), a.len())`) that sits on an
/// unchanged base (`r[mid] == a[mid]`), so splitting both `r` and `a` there is clean —
/// no member straddles the cut — and both sides are non-empty (progress is guaranteed).
/// Returns `None` when no interior match run of length `>= 2` exists, i.e. a genuinely
/// dense oversized block the caller must fall back to typing as one member.
///
/// Position-wise (not an alignment) on purpose: it is O(n) and avoids the very O(n*m)
/// LCS the `MAX_BLOCK` cap exists to skip. It finds the long identical interior that a
/// pair of far-apart point changes leaves behind; an indel-shifted block still yields
/// its longest unshifted run, which is enough to make the block smaller and recurse.
fn longest_unchanged_wall(r: &[u8], a: &[u8]) -> Option<usize> {
    let n = r.len().min(a.len());
    if n < 2 {
        return None;
    }
    let (mut best_start, mut best_len) = (0usize, 0usize);
    let (mut run_start, mut run_len) = (0usize, 0usize);
    for i in 0..n {
        if r[i] == a[i] {
            if run_len == 0 {
                run_start = i;
            }
            run_len += 1;
            if run_len > best_len {
                best_len = run_len;
                best_start = run_start;
            }
        } else {
            run_len = 0;
        }
    }
    if best_len < 2 {
        return None;
    }
    let mid = best_start + best_len / 2;
    // mid lies strictly inside the block on both axes, so both halves are non-empty.
    if mid == 0 || mid >= n {
        return None;
    }
    Some(mid)
}

fn type_run(reference: &[u8], rs: usize, re: usize, alt: &[u8], out: &mut Vec<Member>) {
    let span = re - rs;
    if span == 0 {
        insertion_leaf(reference, rs, alt, out);
        return;
    }
    if alt.is_empty() {
        out.push(Member {
            kind: EditKind::Del,
            ref_start: rs,
            ref_end: re,
            inserted: Vec::new(),
        });
        return;
    }
    if span == 1 && alt.len() == 1 {
        out.push(Member {
            kind: EditKind::Sub,
            ref_start: rs,
            ref_end: re,
            inserted: alt.to_vec(),
        });
        return;
    }
    if alt.len() == span && span >= 2 && *alt == *revcomp(&reference[rs..re]) {
        out.push(Member {
            kind: EditKind::Inv,
            ref_start: rs,
            ref_end: re,
            inserted: Vec::new(),
        });
        return;
    }
    out.push(Member {
        kind: EditKind::Delins,
        ref_start: rs,
        ref_end: re,
        inserted: alt.to_vec(),
    });
}

impl OperatorExtractL2 {
    #[allow(clippy::too_many_arguments)]
    fn decompose(
        &self,
        reference: &[u8],
        resulting: &[u8],
        rs0: usize,
        re0: usize,
        ss0: usize,
        se0: usize,
        depth: usize,
        allow_inv: bool,
        out: &mut Vec<Member>,
    ) {
        let (mut rs, mut re, mut ss, mut se) = (rs0, re0, ss0, se0);
        // Trim shared flanks of this sub-block.
        while rs < re && ss < se && reference[rs] == resulting[ss] {
            rs += 1;
            ss += 1;
        }
        while rs < re && ss < se && reference[re - 1] == resulting[se - 1] {
            re -= 1;
            se -= 1;
        }
        let (r_len, a_len) = (re - rs, se - ss);
        if r_len == 0 && a_len == 0 {
            return;
        }
        let r = &reference[rs..re];
        let a = &resulting[ss..se];

        if depth >= MAX_DEPTH {
            type_run(reference, rs, re, a, out);
            return;
        }
        if r_len > MAX_BLOCK || a_len > MAX_BLOCK {
            // A block too large for the O(n*m) LCS must NOT be ASSERTED as one
            // spanning delins — that reports a merge nothing verified. This is the
            // shipping-side ruling `derivation-may-not-be-bounded-by-the-inputs-
            // spelling` at the bakeoff layer: a block nothing examined is the ABSENCE
            // of a finding, not a finding of no separation. So wall at the longest
            // unchanged run (cheap, position-wise; the cut sits on an unchanged base,
            // so no member straddles it) and recurse, examining each half. Only a
            // genuinely dense oversized block with no interior match run of length
            // >= 2 (essentially never in real cis data — two point changes leave a
            // long identical interior) falls through to one spanning member.
            //
            // Without this, two independent variants separated by > MAX_BLOCK bases
            // collapse into one spanning delins (issue: op-extract over-merge). The
            // >= 2-base wall floor keeps a payload-coincidence block whose interior
            // matches only at single positions (the octamer) whole and <= MAX_BLOCK,
            // so its collapse is unaffected.
            if let Some(mid) = longest_unchanged_wall(r, a) {
                self.decompose(
                    reference,
                    resulting,
                    rs,
                    rs + mid,
                    ss,
                    ss + mid,
                    depth + 1,
                    allow_inv,
                    out,
                );
                self.decompose(
                    reference,
                    resulting,
                    rs + mid,
                    re,
                    ss + mid,
                    se,
                    depth + 1,
                    allow_inv,
                    out,
                );
                return;
            }
            type_run(reference, rs, re, a, out);
            return;
        }

        // 0. Whole-span inversion.
        if r_len >= 2 && r_len == a_len && *a == *revcomp(r) {
            out.push(Member {
                kind: EditKind::Inv,
                ref_start: rs,
                ref_end: re,
                inserted: Vec::new(),
            });
            return;
        }

        // 1. INV anchor: the longest revcomp run, peeled BEFORE the pure-ins/del
        //    leaves so an inversion beside an indel is found. The accept criterion is
        //    RENDERED-COST, not a chance floor: peel the inversion iff the derivation
        //    that names it spells strictly fewer literal bases (each member costs a
        //    fixed K; inv/del/dup name a RANGE and spell zero bases; ins/delins/sub
        //    spell theirs). This is the load-bearing distinction the significance
        //    floor could not make — a REAL inversion of length L replaces a delins
        //    that would spell L bases with an inv that spells none, so peeling wins;
        //    a CHANCE revcomp inside a delins spells the same bases either way plus an
        //    extra structural token, so it loses. It is `mdl`'s cost objective made
        //    local and paired with a structural fallback, which is why it recovers
        //    mdl's inv-neighbour strength without mdl's collapse of multi-member forms.
        //    (`whole-span-reverse-complement-types-as-inv` supplies rule 0's
        //    unconditional whole-span case; here we adjudicate a SUB-span inversion.)
        if allow_inv && r_len >= 2 {
            // Minor 6: allocate the revcomp only on the inv-anchor path, not for
            // every no-inv alternative decomposition.
            let rc = revcomp(r);
            // Minor 2 (NOTE, not fixed): only the SINGLE longest revcomp run is tried as
            // the inv anchor. If it fails the flank guard or the majority-identity veto,
            // a shorter genuine revcomp run with a clean flank is never considered. This
            // is rare (the longest run is almost always the real inversion when there is
            // one), and a multi-run search is a behaviour change that needs its own
            // measurement, so it is documented rather than done here.
            if let Some((x, u, n)) = lcs_substring(&rc, a) {
                if n >= self.inv_min {
                    // revcomp(r)[x..x+n] == a[u..u+n]; ref span [r_len-x-n, r_len-x].
                    let mut ref_lo = rs + (r_len - x - n);
                    let mut ref_hi = rs + (r_len - x);
                    let mut res_lo = ss + u;
                    let mut res_hi = ss + u + n;
                    // Extend outward while the complement holds — recovers an inverted
                    // edge a shared-flank trim nibbled. Each step keeps >=1 base inside
                    // the block so a chance seam pair cannot convert unchanged content.
                    while ref_lo > rs0
                        && res_hi < se0
                        && (ref_lo > rs || res_hi < se)
                        && resulting[res_hi] == comp(reference[ref_lo - 1])
                    {
                        ref_lo -= 1;
                        res_hi += 1;
                    }
                    while ref_hi < re0
                        && res_lo > ss0
                        && (ref_hi < re || res_lo > ss)
                        && resulting[res_lo - 1] == comp(reference[ref_hi])
                    {
                        ref_hi += 1;
                        res_lo -= 1;
                    }
                    // #2155 FLANK GUARD (structural eligibility, adjudicated BEFORE
                    // cost). A sub-span inversion may be peeled only when at least one
                    // of its flanks is unchanged content or a block edge. An inversion
                    // with CHANGED content on BOTH sides is a payload coincidence
                    // interior to a changed region: it collapses into the spanning
                    // delins. `whole-span-reverse-complement-types-as-inv` types a WHOLE
                    // span's exact revcomp as `inv`, and a doubly-interior revcomp is
                    // not a whole span — so it is not that ruling's `inv`. This is a
                    // STRUCTURAL rule cost cannot express: a structured member always
                    // spells zero literal bases, so cost ALWAYS prefers to peel it
                    // (the octamer break — a chance 3-base revcomp sandwiched between
                    // two delins undercuts the spanning delins at any K). The guard,
                    // not cost, keeps such a coincidence collapsed. A flank is "clean"
                    // iff it is empty (edge) or byte-identical (unchanged); an indel or
                    // substitution in it makes it changed.
                    let flank5_clean = reference[rs0..ref_lo] == resulting[ss0..res_lo];
                    let flank3_clean = reference[ref_hi..re0] == resulting[res_hi..se0];
                    // M1: MAJORITY-IDENTITY VETO. The flank guard's edge arm is a hole:
                    // a chance revcomp SEED (inv_min = 2) in low-complexity content can
                    // extend outward via the loops above, ANNEX unchanged bases, and
                    // reach a block edge — so one flank slice is empty and the guard
                    // passes, while the peeled span is mostly positionally UNCHANGED.
                    // That is a mostly-unchanged region wearing an `inv` label, not an
                    // inversion: its matched columns are unchanged under every minimal
                    // alignment (`unchanged-is-read-over-every-minimal-alignment`), so
                    // `general.md:33` requires the members individually — cost cannot
                    // refuse (the inv spells 0 bases AND drops member count, so it wins
                    // at any K). A REAL inversion is mostly-CHANGED (a genuine revcomp
                    // has only ~25% chance interior matches). So refuse the peel when
                    // >= half the span's columns match positionally. Measured: lifts
                    // Sub+Sub+Sub 35% -> 100% with the aggregate flat, and the blocks
                    // gained sit in decided-ruling territory. The 1/2 threshold is house
                    // policy (operator's call), not a clause. `span` is the inv length;
                    // ref and res spans are equal by construction (an exact revcomp).
                    let span = ref_hi - ref_lo;
                    let positional_matches = (0..span)
                        .filter(|&i| reference[ref_lo + i] == resulting[res_lo + i])
                        .count();
                    let majority_identity = positional_matches * 2 >= span;
                    if (flank5_clean || flank3_clean) && !majority_identity {
                        // Eligible. Peel candidate (inv allowed inside the flanks too).
                        let mut cand = Vec::new();
                        self.decompose(
                            reference,
                            resulting,
                            rs0,
                            ref_lo,
                            ss0,
                            res_lo,
                            depth + 1,
                            true,
                            &mut cand,
                        );
                        cand.push(Member {
                            kind: EditKind::Inv,
                            ref_start: ref_lo,
                            ref_end: ref_hi,
                            inserted: Vec::new(),
                        });
                        self.decompose(
                            reference,
                            resulting,
                            ref_hi,
                            re0,
                            res_hi,
                            se0,
                            depth + 1,
                            true,
                            &mut cand,
                        );
                        // No-inv alternative: the SAME block typed without any inv peel.
                        // Cost breaks this remaining peel-vs-no-peel tie among the
                        // spec-legal readings (the inv buckets measure ~92-96 on the 17k
                        // corpus); it is the ONLY place the typer consults cost.
                        let mut alt = Vec::new();
                        self.decompose(
                            reference,
                            resulting,
                            rs0,
                            re0,
                            ss0,
                            se0,
                            depth + 1,
                            false,
                            &mut alt,
                        );
                        let ctx = &resulting[ss0..se0];
                        if self.score(&cand, ctx) < self.score(&alt, ctx) {
                            out.extend(cand);
                        } else {
                            out.extend(alt);
                        }
                        return;
                    }
                    // Not eligible — fall through to the no-inv leaves (steps 2-7).
                }
            }
        }

        // 2/2a/3/4. Pure leaves.
        if r_len == 0 {
            // 2a. #2175 / P2 dup-beside-sub: a pure insertion whose bases re-align as a
            //     tandem copy of the immediate 5' reference plus one substituted byte
            //     renders as [dup;sub], not a flat ins. Gated on `dup_complete` (reuses
            //     that family's wholly-consumed completion discipline). Finding #7.
            if self.dup_complete {
                if let Some((unit, copies)) = dup_beside_sub_unit(reference, rs, a) {
                    for _ in 0..copies {
                        out.push(Member {
                            kind: EditKind::Dup,
                            ref_start: rs - 1,
                            ref_end: rs - 1,
                            inserted: unit.clone(),
                        });
                    }
                    out.push(Member {
                        kind: EditKind::Sub,
                        ref_start: rs - 1,
                        ref_end: rs,
                        inserted: vec![a[a.len() - 1]],
                    });
                    return;
                }
            }
            insertion_leaf(reference, rs, a, out);
            return;
        }
        if a_len == 0 {
            out.push(Member {
                kind: EditKind::Del,
                ref_start: rs,
                ref_end: re,
                inserted: Vec::new(),
            });
            return;
        }
        if r_len == 1 && a_len == 1 {
            out.push(Member {
                kind: EditKind::Sub,
                ref_start: rs,
                ref_end: re,
                inserted: a.to_vec(),
            });
            return;
        }

        // 5. DUP anchor on a net insertion: a maximal exact tandem copy of the
        //    adjacent reference flank, peeled as a zero-width dup, recurse residual.
        if a_len > r_len {
            // 2a graft: wholly-consumed dup peel (multi-unit, 1-base allowed). Peel
            // copies of the immediate 5' reference until the net insertion is fully
            // consumed; accept ONLY on complete consumption (equal-length residual),
            // which is what stops a chance 1-base flank from fragmenting an ins.
            if self.dup_complete {
                // 5' edge tandem (M2). The net insertion is peeled as a dup only when
                // it is k >= 1 whole copies of a SINGLE unit that copies the immediate
                // 5' reference flank — the genuine tandem-dup shape
                // (`duplication-must-ranks-the-label-not-the-partition`, #2175). The
                // old code matched a greedy sequence of DIFFERENT flank-prefix lengths,
                // which let a "parade of chance 1-base copies" consume a net insertion
                // in low-complexity content (e.g. `AATATATTT -> TAATTATAATTTT` fragmented
                // into `Ins,Dup,Dup,Dup`). Requiring one repeated unit keeps genuine
                // homopolymer/tandem dups (a 1-base unit repeated is a real homopolymer
                // expansion) while rejecting the parade. The canonical unit is the
                // SHORTEST period, so a single non-repetitive copy resolves to one dup.
                let net = a_len - r_len;
                let mut tandem_unit = 0usize;
                for l in 1..=net {
                    if net % l != 0 || rs < l {
                        continue;
                    }
                    // The unit must copy the immediate 5' flank ...
                    if reference[rs - l..rs] != resulting[ss..ss + l] {
                        continue;
                    }
                    // ... and the whole insertion must be l-periodic (k copies of it).
                    if (0..net).all(|j| resulting[ss + j] == resulting[ss + (j % l)]) {
                        tandem_unit = l;
                        break;
                    }
                }
                if tandem_unit > 0 {
                    let unit = resulting[ss..ss + tandem_unit].to_vec();
                    for _ in 0..(net / tandem_unit) {
                        out.push(Member {
                            kind: EditKind::Dup,
                            ref_start: rs,
                            ref_end: rs,
                            inserted: unit.clone(),
                        });
                    }
                    self.decompose(
                        reference,
                        resulting,
                        rs,
                        re,
                        ss + net,
                        se,
                        depth + 1,
                        allow_inv,
                        out,
                    );
                    return;
                }
                // 3' edge, single-unit shifted form (mirror of `anchor_peel`): a suffix
                // copying reference[re..re+u] 3'-shifts through those bases; the dup's
                // canonical position is `re + u`, payload read from the matched suffix.
                let mut end = se;
                let mut net3 = a_len - r_len;
                let mut units3 = 0usize;
                let mut last_u = 0usize;
                while net3 > 0 {
                    let maxu = net3.min(reference.len() - re).min(end - ss);
                    let mut found = 0usize;
                    for u in (1..=maxu).rev() {
                        if reference[re..re + u] == resulting[end - u..end] {
                            found = u;
                            break;
                        }
                    }
                    if found == 0 {
                        break;
                    }
                    units3 += 1;
                    last_u = found;
                    end -= found;
                    net3 -= found;
                }
                if net3 == 0 && units3 == 1 && end + last_u <= se {
                    self.decompose(
                        reference,
                        resulting,
                        rs,
                        re,
                        ss,
                        end,
                        depth + 1,
                        allow_inv,
                        out,
                    );
                    out.push(Member {
                        kind: EditKind::Dup,
                        ref_start: re + last_u,
                        ref_end: re + last_u,
                        inserted: resulting[end..end + last_u].to_vec(),
                    });
                    return;
                }
                // Neither edge completed — fall through to the floored single-peel.
            }
            // Floored single-peel. In `dup_complete` mode the floor is raised to >= 2
            // (M2): a genuine 1-base tandem is already caught by the completion branch
            // above, so a 1-base peel HERE — where completion FAILED — is a chance flank
            // match that would fragment an insertion, exactly the parade this fix closes.
            let dup_floor = if self.dup_complete {
                self.dup_min.max(2)
            } else {
                self.dup_min
            };
            // 5' prefix copy: a[0..k] == reference[rs-k..rs].
            let mut k5 = 0usize;
            let maxk = rs.min(a_len);
            for k in (dup_floor..=maxk).rev() {
                if reference[rs - k..rs] == a[..k] {
                    k5 = k;
                    break;
                }
            }
            // 3' suffix copy: a[a_len-k..] == reference[re..re+k].
            let mut k3 = 0usize;
            let maxk3 = (reference.len() - re).min(a_len);
            for k in (dup_floor..=maxk3).rev() {
                if reference[re..re + k] == a[a_len - k..] {
                    k3 = k;
                    break;
                }
            }
            // Prefer the longer copy; tie -> 5' (the earlier-in-reference tandem).
            if k5 >= k3 && k5 >= dup_floor {
                out.push(Member {
                    kind: EditKind::Dup,
                    ref_start: rs,
                    ref_end: rs,
                    inserted: a[..k5].to_vec(),
                });
                self.decompose(
                    reference,
                    resulting,
                    rs,
                    re,
                    ss + k5,
                    se,
                    depth + 1,
                    allow_inv,
                    out,
                );
                return;
            }
            if k3 >= dup_floor {
                out.push(Member {
                    kind: EditKind::Dup,
                    ref_start: re,
                    ref_end: re,
                    inserted: a[a_len - k3..].to_vec(),
                });
                self.decompose(
                    reference,
                    resulting,
                    rs,
                    re,
                    ss,
                    se - k3,
                    depth + 1,
                    allow_inv,
                    out,
                );
                return;
            }
        }

        // 6. Equal-length block: split at position-wise fixed points; type each run.
        //    A position-wise-unchanged column in an equal-length block is a genuine
        //    wall; a minimal-edit interior match (the #2174 rotation) is not, and this
        //    keeps a contiguous equal-length change as one run.
        if r_len == a_len {
            let mut i = 0usize;
            while i < r_len {
                if r[i] == a[i] {
                    i += 1;
                    continue;
                }
                let st = i;
                while i < r_len && r[i] != a[i] {
                    i += 1;
                }
                type_run(reference, rs + st, rs + i, &a[st..i], out);
            }
            return;
        }

        // 7. STRUCTURAL LEAF (cost NEVER consulted). Split at the FORCED-UNCHANGED
        //    walls (`unchanged-is-read-over-every-minimal-alignment`) and keep the
        //    split — one member per segment. A forced-unchanged column is a genuine
        //    separation under `general.md:33` (unchanged under EVERY minimal
        //    alignment), so splitting there is what the spec REQUIRES, not a choice to
        //    adjudicate. A single delins has no such interior wall (`segs.len() <= 1`)
        //    and types as one delins.
        //
        //    This deliberately drops the old cost test that could MERGE across a
        //    genuine wall into one spanning delins. That merge is the DNA
        //    payload-coincidence collapse (`delins-payload-coincidence-carve-out-is-
        //    coding-dna-scoped`, `unequal-length-block-a-placed-gap-is-not-a-
        //    separation`) — c.-axis-scoped, molecule-dependent, and therefore the job
        //    of the downstream DNA/RNA axis layer, NOT of this molecule-blind typer.
        //    The spec-priority reading here is: keep members individual; let the axis
        //    collapse where a molecule licenses it. Removing cost from the leaf is what
        //    recovers Delins (never mis-merged) and multi-Sub (positional walls never
        //    cost-merged into a delins).
        // Minor 7: `Molecule::Dna` is inert here — `AllAlignmentSplitL1` ignores the
        // molecule (its walls are the forced-unchanged columns, a molecule-blind fact);
        // it is not an axis leak. And the ∅ provenance is likewise inert here.
        let segs = AllAlignmentSplitL1.segment(r, a, Molecule::Dna, &Provenance::none());
        if segs.len() <= 1 {
            type_run(reference, rs, re, a, out);
            return;
        }
        // Minor 1: recurse `decompose` per segment rather than `type_run`, so each
        // segment is RE-TRIMMED (a forced-unchanged wall can leave a non-forced matching
        // edge base inside a segment) and gets the pure-leaf refinements of steps 2-4
        // (a Sub/Ins/Del rather than a spanning Delins). `allow_inv = false` keeps the
        // leaf inv-free — the inv anchor already had its chance at this block, and a
        // genuine whole-segment revcomp is still caught by step 0 / step 6 / `type_run`.
        // Segments are strictly smaller when `segs.len() >= 2`, and depth increments, so
        // recursion terminates.
        for s in &segs {
            self.decompose(
                reference,
                resulting,
                rs + s.ref_start,
                rs + s.ref_end,
                ss + s.res_start,
                ss + s.res_end,
                depth + 1,
                false,
                out,
            );
        }
    }
}

impl L2Typer for OperatorExtractL2 {
    fn name(&self) -> &str {
        &self.name
    }
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        _f: &FrameContext,
    ) -> Vec<Member> {
        let mut out = Vec::new();
        self.decompose(
            reference,
            resulting,
            seg.ref_start,
            seg.ref_end,
            seg.res_start,
            seg.res_end,
            0,
            true,
            &mut out,
        );
        out
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member> {
        apply_dial_b(m, reference, resulting, ctx, cfg)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Run the typer over a whole block and return its members.
    fn members_of(t: &OperatorExtractL2, r: &[u8], a: &[u8]) -> Vec<Member> {
        let mut out = Vec::new();
        t.decompose(r, a, 0, r.len(), 0, a.len(), 0, true, &mut out);
        out
    }

    /// Reconstruct the resulting sequence from the members (mirrors `RefApplier`),
    /// so every test can assert the partition is VALID.
    fn reconstruct(reference: &[u8], members: &[Member]) -> Vec<u8> {
        let mut ms = members.to_vec();
        ms.sort_by_key(|m| (m.ref_start, m.ref_end));
        let mut out = Vec::new();
        let mut cursor = 0usize;
        for m in &ms {
            out.extend_from_slice(&reference[cursor..m.ref_start]);
            match m.kind {
                EditKind::Identity => out.extend_from_slice(&reference[m.ref_start..m.ref_end]),
                EditKind::Del => {}
                EditKind::Inv => out.extend(revcomp(&reference[m.ref_start..m.ref_end])),
                _ => out.extend_from_slice(&m.inserted),
            }
            cursor = m.ref_end;
        }
        out.extend_from_slice(&reference[cursor..]);
        out
    }

    fn count(members: &[Member], k: EditKind) -> usize {
        members.iter().filter(|m| m.kind == k).count()
    }

    /// #2155: a 3-base revcomp (`GTT`->`AAC`) sandwiched between two delins — changed
    /// on BOTH sides. The flank guard must keep it collapsed, not peel it as an inv.
    #[test]
    fn octamer_payload_coincidence_is_not_peeled_as_inv() {
        let t = OperatorExtractL2::v2(4, 2, "t");
        let (r, a): (&[u8], &[u8]) = (b"CTTAGTTA", b"AAACAAAC");
        let m = members_of(&t, r, a);
        assert_eq!(reconstruct(r, &m), a, "must reconstruct");
        assert_eq!(
            count(&m, EditKind::Inv),
            0,
            "no inv on a doubly-interior revcomp: {m:?}"
        );
    }

    /// A genuine inversion at a block edge (`ACG`->`CGT`, unchanged `TT` 3' flank):
    /// mostly-changed span, one clean flank — accepted as one inv.
    #[test]
    fn a_genuine_edge_inversion_is_peeled() {
        let t = OperatorExtractL2::v2(4, 2, "t");
        let (r, a): (&[u8], &[u8]) = (b"ACGTT", b"CGTTT");
        let m = members_of(&t, r, a);
        assert_eq!(reconstruct(r, &m), a);
        assert_eq!(
            count(&m, EditKind::Inv),
            1,
            "genuine edge inv accepted: {m:?}"
        );
    }

    /// M1: three isolated subs (cols 0/5/9). A chance 2-base revcomp seed can extend to
    /// the block edge and satisfy the guard's edge arm; the majority-identity veto
    /// refuses it (the span is mostly unchanged). Must be three subs, no inv.
    #[test]
    fn m1_chance_seed_extended_to_edge_is_not_an_inv() {
        let t = OperatorExtractL2::v2(4, 2, "t");
        let (r, a): (&[u8], &[u8]) = (b"AAATTAATAA", b"TAATTTATAT");
        let m = members_of(&t, r, a);
        assert_eq!(reconstruct(r, &m), a);
        assert_eq!(count(&m, EditKind::Inv), 0, "no chance inv: {m:?}");
        assert_eq!(count(&m, EditKind::Sub), 3, "three individual subs: {m:?}");
    }

    /// M2: a net insertion in low-complexity content must NOT be consumed by a "parade
    /// of chance 1-base copies" of DIFFERENT units. At most one dup survives (a genuine
    /// tandem is one repeated unit); the remainder is an ins/delins leaf.
    #[test]
    fn m2_parade_of_chance_1base_dups_is_refused() {
        let t = OperatorExtractL2::v2(4, 2, "t");
        let (r, a): (&[u8], &[u8]) = (b"AATATATTT", b"TAATTATAATTTT");
        let m = members_of(&t, r, a);
        assert_eq!(reconstruct(r, &m), a);
        assert!(
            count(&m, EditKind::Dup) <= 1,
            "no parade of chance dups: {m:?}"
        );
    }

    /// M2: a REAL 1-base tandem (`A` duplicated) beside a sub still peels as a Dup —
    /// the fix rejects the parade without losing genuine tandems.
    #[test]
    fn m2_genuine_1base_tandem_dup_is_kept() {
        let t = OperatorExtractL2::v2(4, 2, "t");
        let (r, a): (&[u8], &[u8]) = (b"CAG", b"CAAT");
        let m = members_of(&t, r, a);
        assert_eq!(reconstruct(r, &m), a);
        assert_eq!(
            count(&m, EditKind::Dup),
            1,
            "genuine tandem dup kept: {m:?}"
        );
    }
}
