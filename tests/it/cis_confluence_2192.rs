//! #2192 / #2194 confluence corpus — generated spellings, one canonical output.
//!
//! The fragmentation family (#2155 / #2174 / #2175 / #2194) is about a contiguous
//! run, or a run beside a net-imbalanced cis member, coalescing to one `delins`
//! instead of fragmenting. The property that matters is **confluence**: every legal
//! spelling of one `(reference, resulting)` pair must normalize to the SAME form.
//!
//! This builds synthetic multi-member patterns from typed members whose equivalent
//! spellings are known by construction (a run as `delins`, as position-wise
//! substitutions, or as `del`+`ins`; a tandem `dup` as `dup` or as `ins`), takes
//! the cross-product of per-member spellings and member orderings (cis members are
//! order-independent), and **filters every candidate through the independent
//! `apply` oracle** — so a candidate that does not genuinely denote the pair is
//! discarded rather than counted as a confluence failure. The survivors (20–100+
//! per pattern) must all normalize to one output, in both shuffle directions.
//!
//! It complements `fragmentation_corpus` (which pins the recommended FORM) by
//! attacking the CONVERGENCE of many spellings onto whatever that form is, without
//! pinning the string — a spelling-diversity test, not an output pin.

use crate::common::cis_apply_oracle::{apply, normalize_in};
use ferro_hgvs::ShuffleDirection;
use std::collections::BTreeSet;

const PAD: &str = "CACGTACGTC";
/// >= 2 unchanged bases so `general.md:33` keeps flanking members individual.
const SPACER: &str = "TCGGATCAG";

/// A member of a synthetic pattern. Each knows the reference bases it sits on, the
/// bases it contributes to the resulting sequence, and every equivalent HGVS body
/// (sans the `TEMPLATE:g.` prefix) for its own locus.
#[derive(Clone, Copy)]
enum Member {
    /// Equal-length replacement `ref_seg -> alt_seg` (a run). Spellings: one
    /// spanning `delins`, the position-wise per-column substitutions, and — when it
    /// is a single-rotation — a `del`+`ins`. All are apply-filtered.
    Run {
        ref_seg: &'static str,
        alt_seg: &'static str,
    },
    /// A tandem duplication of `unit` (>= 2 bases): the reference carries `unit`
    /// once, the result twice. Spellings: `dup` and the equivalent `ins`.
    Dup { unit: &'static str },
    /// A plain deletion of `seg`.
    Del { seg: &'static str },
    /// A single substitution `from -> to`.
    Sub { from: char, to: char },
}

impl Member {
    fn ref_seg(&self) -> String {
        match self {
            Member::Run { ref_seg, .. } => (*ref_seg).to_string(),
            Member::Dup { unit } => (*unit).to_string(),
            Member::Del { seg } => (*seg).to_string(),
            Member::Sub { from, .. } => from.to_string(),
        }
    }

    /// Bases this member contributes to the resulting sequence.
    fn alt_seg(&self) -> String {
        match self {
            Member::Run { alt_seg, .. } => (*alt_seg).to_string(),
            Member::Dup { unit } => format!("{unit}{unit}"),
            Member::Del { .. } => String::new(),
            Member::Sub { to, .. } => to.to_string(),
        }
    }

    /// Every equivalent HGVS body for this member placed with 1-based reference
    /// start `p`. Candidates only — the caller apply-filters, so a body that does
    /// not denote the pair is dropped rather than trusted.
    fn spellings(&self, p: usize) -> Vec<String> {
        match self {
            Member::Run { ref_seg, alt_seg } => {
                let len = ref_seg.len();
                let end = p + len - 1;
                let mut out = vec![format!("{p}_{end}delins{alt_seg}")];
                // Position-wise substitutions: one `X>Y` per differing column. Always
                // a legal spelling of an equal-length block (replace each column).
                let subs: Vec<String> = ref_seg
                    .bytes()
                    .zip(alt_seg.bytes())
                    .enumerate()
                    .filter(|(_, (r, a))| r != a)
                    .map(|(i, (r, a))| format!("{}{}>{}", p + i, r as char, a as char))
                    .collect();
                if subs.len() >= 2 {
                    out.push(subs.join(";"));
                }
                // A single left/right rotation reads as one `del` plus one `ins`.
                if len >= 2 {
                    let last = alt_seg.as_bytes()[len - 1] as char;
                    out.push(format!("{p}del;{end}_{}ins{last}", end + 1));
                    let first = alt_seg.as_bytes()[0] as char;
                    out.push(format!("{p}_{p}ins{first};{end}del"));
                }
                out
            }
            Member::Dup { unit } => {
                let len = unit.len();
                let end = p + len - 1;
                vec![
                    format!("{p}_{end}dup"),
                    format!("{end}_{}ins{unit}", end + 1),
                ]
            }
            Member::Del { seg } => {
                let len = seg.len();
                let end = p + len - 1;
                if len == 1 {
                    vec![format!("{p}del")]
                } else {
                    vec![format!("{p}_{end}del")]
                }
            }
            Member::Sub { from, to } => vec![format!("{p}{from}>{to}")],
        }
    }
}

/// A pattern: members laid out left to right in a padded reference, separated by an
/// unchanged spacer. Returns `(reference, resulting, one_based_start_per_member)`.
fn build(members: &[Member]) -> (String, String, Vec<usize>) {
    let mut reference = String::from(PAD);
    let mut resulting = String::from(PAD);
    let mut starts = Vec::new();
    for (i, m) in members.iter().enumerate() {
        if i > 0 {
            reference.push_str(SPACER);
            resulting.push_str(SPACER);
        }
        starts.push(reference.len() + 1); // 1-based reference coordinate
        reference.push_str(&m.ref_seg());
        resulting.push_str(&m.alt_seg());
    }
    reference.push_str(PAD);
    resulting.push_str(PAD);
    (reference, resulting, starts)
}

/// Bounded permutations of `0..n` (member orderings), capped so a wide pattern does
/// not blow up. Cis members are order-independent, so every ordering is a spelling.
fn orderings(n: usize, cap: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut cur: Vec<usize> = (0..n).collect();
    // Identity first, then a bounded set of rotations + swaps — enough diversity
    // without factorial blow-up.
    out.push(cur.clone());
    for shift in 1..n {
        cur.rotate_left(1);
        out.push(cur.clone());
        if out.len() >= cap {
            return out;
        }
        let _ = shift;
    }
    // Adjacent-swap variants of the identity.
    for i in 0..n.saturating_sub(1) {
        let mut v: Vec<usize> = (0..n).collect();
        v.swap(i, i + 1);
        out.push(v);
        if out.len() >= cap {
            break;
        }
    }
    out
}

/// Every apply-verified spelling of a pattern: cross-product of per-member spellings
/// and member orderings, wrapped as an allele, kept only if it denotes the pair.
fn spellings_of(
    members: &[Member],
    reference: &str,
    resulting: &str,
    starts: &[usize],
) -> Vec<String> {
    let per_member: Vec<Vec<String>> = members
        .iter()
        .zip(starts)
        .map(|(m, &p)| m.spellings(p))
        .collect();
    // Cross-product of one spelling choice per member.
    let mut choices: Vec<Vec<String>> = vec![Vec::new()];
    for opts in &per_member {
        let mut next = Vec::new();
        for prefix in &choices {
            for o in opts {
                let mut c = prefix.clone();
                c.push(o.clone());
                next.push(c);
            }
        }
        choices = next;
        if choices.len() > 4096 {
            choices.truncate(4096);
        }
    }
    let mut spellings: BTreeSet<String> = BTreeSet::new();
    for choice in &choices {
        for ord in orderings(members.len(), 8) {
            let bodies: Vec<&str> = ord.iter().map(|&i| choice[i].as_str()).collect();
            // Each body may itself carry `;` (a run's del+ins / subs), so joining
            // yields a flat member list — exactly what a cis allele is.
            let inner = bodies.join(";");
            let desc = if members.len() == 1 && !inner.contains(';') {
                format!("TEMPLATE:g.{inner}")
            } else {
                format!("TEMPLATE:g.[{inner}]")
            };
            if apply(reference, &desc).as_deref() == Some(resulting) {
                spellings.insert(desc);
            }
        }
    }
    spellings.into_iter().collect()
}

/// The synthetic #2192/#2194 corpus: net-imbalanced members (`dup`, `del`) beside
/// equal-length runs and substitutions, in the orderings the family stresses.
fn corpus() -> Vec<Vec<Member>> {
    let run4 = Member::Run {
        ref_seg: "GACT",
        alt_seg: "ACTG",
    };
    let run3 = Member::Run {
        ref_seg: "ATA",
        alt_seg: "TAC",
    };
    let dup = Member::Dup { unit: "AG" };
    let del = Member::Del { seg: "T" };
    let del2 = Member::Del { seg: "GT" };
    let sub = Member::Sub { from: 'C', to: 'A' };
    vec![
        // Single runs (the #2174 base case).
        vec![run4],
        vec![run3],
        // Run beside a substitution (#2192 Ex1/Ex3 shape).
        vec![sub, run4],
        vec![run4, sub],
        // Dup 5' of a run (#2194 — the fix under test).
        vec![dup, run4],
        vec![dup, run3, sub],
        // Del 5' of a run (the broader net-imbalanced-5'-member class).
        vec![del, run4],
        vec![del2, run4, sub],
        // Runs on both sides of a net-imbalanced member.
        vec![run4, dup, run3],
        vec![run3, del, run4],
        // Denser mixes.
        vec![sub, run4, dup, run3],
        vec![dup, sub, run4, del, run3],
    ]
}

fn assert_confluent(members: &[Member], direction: ShuffleDirection) {
    let (reference, resulting, starts) = build(members);
    let spellings = spellings_of(members, &reference, &resulting, &starts);
    // Non-vacuity: the generator must have produced a real spread of spellings, or a
    // silent zero would pass while testing nothing (the "a corpus zero is a claim
    // about the corpus" trap).
    // A lone run still yields three genuinely different spellings (spanning
    // `delins`, position-wise substitutions, and `del`+`ins`); a multi-member
    // pattern yields many more (per-member forms x orderings).
    let floor = if members.len() == 1 { 3 } else { 6 };
    assert!(
        spellings.len() >= floor,
        "only {} apply-verified spellings for {members:?}-shaped pattern (floor {floor}) — \
         the generator is not exercising the confluence it claims to",
        spellings.len(),
    );
    // Every spelling normalizes to ONE canonical form.
    let mut outputs: BTreeSet<String> = BTreeSet::new();
    let mut first: Option<(String, String)> = None;
    for s in &spellings {
        let out = normalize_in(&reference, s, direction);
        // Meaning preservation: the normalized output still denotes the pair.
        assert_eq!(
            apply(&reference, &out).as_deref(),
            Some(resulting.as_str()),
            "{direction:?}: normalized output does not denote the pair\n  spelling: {s}\n  output:   {out}"
        );
        if first.is_none() {
            first = Some((s.clone(), out.clone()));
        }
        outputs.insert(out);
    }
    assert_eq!(
        outputs.len(),
        1,
        "{direction:?}: {} spellings did NOT converge — {} distinct outputs:\n{}",
        spellings.len(),
        outputs.len(),
        outputs
            .iter()
            .take(8)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n"),
    );
}

// Debug is needed for the assertion messages above.
impl std::fmt::Debug for Member {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Member::Run { ref_seg, alt_seg } => write!(f, "Run({ref_seg}->{alt_seg})"),
            Member::Dup { unit } => write!(f, "Dup({unit})"),
            Member::Del { seg } => write!(f, "Del({seg})"),
            Member::Sub { from, to } => write!(f, "Sub({from}>{to})"),
        }
    }
}

#[test]
fn generated_spellings_converge_three_prime() {
    for members in corpus() {
        assert_confluent(&members, ShuffleDirection::ThreePrime);
    }
}

#[test]
fn generated_spellings_converge_five_prime() {
    for members in corpus() {
        assert_confluent(&members, ShuffleDirection::FivePrime);
    }
}
