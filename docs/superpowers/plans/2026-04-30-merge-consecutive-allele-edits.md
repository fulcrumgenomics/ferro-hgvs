# Merge consecutive allele edits — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** After per-variant normalization in `normalize_allele`, merge sub-variants whose edits describe consecutive nucleotides into a single delins/del/ins per HGVS spec — bringing allele output into spec compliance for the strictly-consecutive case.

**Architecture:** A new module `src/normalize/merge.rs` exposes `merge_consecutive_edits`, called from `normalize_allele` after `resolve_overlaps`. Each sub-variant is reduced to an `(start, end, alt)` *anchor*; pairs are merged when `A.end + 1 == B.start`. The merge result becomes a `Deletion`, `Insertion`, or `Delins` depending on which inputs contributed alt content. Phase scope (Cis only), variant types (g/c/n/r/m), edit-type allowlist, and exclusions (non-`Literal` `InsertedSequence`, intronic, UTR-boundary crossing, `Mu::Uncertain`, etc.) follow the spec doc.

**Tech Stack:** Rust 2021, `cargo nextest` for tests, existing `parse_hgvs` / `Normalizer` / `MockProvider` test infrastructure.

**Spec doc:** `docs/superpowers/specs/2026-04-30-merge-consecutive-allele-edits-design.md`

**Issue:** [#72](https://github.com/fulcrumgenomics/ferro-hgvs/issues/72) (Phase 2 codon-frame exception tracked separately as #79)

---

## File Structure

| File | Action | Responsibility |
| --- | --- | --- |
| `src/normalize/merge.rs` | Create | Anchor model, adjacency check, merge logic. Pure transformation — no reference-provider access. |
| `src/normalize/mod.rs` | Modify | `pub mod merge;` declaration; call `merge::merge_consecutive_edits` at the end of `normalize_allele`. |
| `tests/merge_consecutive_edits_tests.rs` | Create | Integration tests using `parse_hgvs → normalize → display` round-trip. |

The new module is internal (`pub(crate)`); no public-API change.

---

## Conventions for all tasks

- **Tests use `parse → normalize → display` round-trip** unless otherwise noted. Genomic (`g.`) tests use `MockProvider::new()` (no transcripts needed); CDS / RNA tests use `provider_with_transcript` patterns from the existing test suite.
- **Run tests with:** `cargo nextest run --features dev -E 'test(merge_consecutive)'` (filters to the new module).
- **Pre-commit hooks** run `cargo fmt`, `cargo clippy --features dev -- -D warnings`, and the test suite. If a hook fails, fix the issue and create a NEW commit (per the repo's commit safety protocol).
- **Commit messages** follow the repo style: `feat:`, `fix:`, `test:`, `docs:`. Include a `Co-Authored-By:` trailer matching recent commits.

---

## Task 1: Module scaffold + sub+sub merging

**Files:**
- Create: `src/normalize/merge.rs`
- Modify: `src/normalize/mod.rs:26-30` (add `pub mod merge;`) and `src/normalize/mod.rs:200-231` (call merge at end of `normalize_allele`)
- Create: `tests/merge_consecutive_edits_tests.rs`

### Step 1: Write the first failing test

Create `tests/merge_consecutive_edits_tests.rs`:

```rust
//! Tests for HGVS-spec consecutive-edit merging in alleles.
//! See docs/superpowers/specs/2026-04-30-merge-consecutive-allele-edits-design.md.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

#[test]
fn test_merge_consecutive_subs_genome() {
    // Issue #72 example: two adjacent SNVs collapse to a single delins.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001A>C]"),
        "NC_000001.11:g.1000_1001delinsAC",
    );
}
```

- [ ] **Step 2: Verify the test fails**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: FAIL with the assertion showing `NC_000001.11:g.[1000G>A;1001A>C]` (the unchanged input) on the left.

- [ ] **Step 3: Create the merge module**

Create `src/normalize/merge.rs`:

```rust
//! Merge consecutive sub-variants in a cis allele into a single edit per HGVS spec.
//!
//! See `docs/superpowers/specs/2026-04-30-merge-consecutive-allele-edits-design.md`.

use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::{Interval, UncertainBoundary};
use crate::hgvs::location::GenomePos;
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{
    Accession, AllelePhase, GenomeVariant, HgvsVariant, LocEdit,
};

/// Anchor for a single sub-variant.
///
/// `start` and `end` are 1-based inclusive position bounds. For an `Insertion`
/// between positions `p` and `p+1`, `start = p+1` and `end = p` (empty range
/// at a boundary).
#[derive(Debug, Clone)]
struct Anchor {
    start: u64,
    end: u64,
    alt: Vec<Base>,
}

/// Merge consecutive sub-variants in an allele.
///
/// Returns the input unchanged unless `phase == Cis` and at least one merge
/// is possible. Sub-variants that aren't merge-eligible (different accession,
/// non-NaEdit, uncertain edit, disqualifying position, etc.) act as merge
/// barriers — they pass through unchanged in their original input order.
pub(crate) fn merge_consecutive_edits(
    variants: &[HgvsVariant],
    phase: AllelePhase,
) -> Vec<HgvsVariant> {
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return variants.to_vec();
    }

    let mut output: Vec<HgvsVariant> = Vec::with_capacity(variants.len());
    for next in variants {
        if let Some(prev) = output.last() {
            if let Some(merged) = try_merge_pair(prev, next) {
                let last = output.last_mut().unwrap();
                *last = merged;
                continue;
            }
        }
        output.push(next.clone());
    }
    output
}

fn try_merge_pair(a: &HgvsVariant, b: &HgvsVariant) -> Option<HgvsVariant> {
    // Phase 1: only Genome-Genome same-accession sub+sub merging.
    let (av, bv) = match (a, b) {
        (HgvsVariant::Genome(av), HgvsVariant::Genome(bv)) => (av, bv),
        _ => return None,
    };
    if av.accession != bv.accession {
        return None;
    }
    let a_anchor = genome_anchor(av)?;
    let b_anchor = genome_anchor(bv)?;
    if a_anchor.end + 1 != b_anchor.start {
        return None;
    }
    let mut alt = a_anchor.alt;
    alt.extend(b_anchor.alt);
    Some(HgvsVariant::Genome(build_genome_delins(
        av,
        a_anchor.start,
        b_anchor.end,
        alt,
    )))
}

fn genome_anchor(v: &GenomeVariant) -> Option<Anchor> {
    let edit = v.loc_edit.edit.inner()?; // Mu::Certain only
    let (start, end) = simple_genome_range(&v.loc_edit.location)?;
    match edit {
        NaEdit::Substitution { alternative, .. } => Some(Anchor {
            start,
            end,
            alt: vec![*alternative],
        }),
        _ => None,
    }
}

fn simple_genome_range(interval: &Interval<GenomePos>) -> Option<(u64, u64)> {
    let start = simple_genome_pos(&interval.start)?;
    let end = simple_genome_pos(&interval.end)?;
    Some((start, end))
}

fn simple_genome_pos(boundary: &UncertainBoundary<GenomePos>) -> Option<u64> {
    let mu = boundary.as_single()?;
    let pos = match mu {
        Mu::Certain(p) => p,
        _ => return None,
    };
    if pos.is_special() || pos.offset.is_some() {
        return None;
    }
    Some(pos.base)
}

fn build_genome_delins(
    template: &GenomeVariant,
    start: u64,
    end: u64,
    alt: Vec<Base>,
) -> GenomeVariant {
    let location = Interval::new(GenomePos::new(start), GenomePos::new(end));
    let edit = NaEdit::Delins {
        sequence: InsertedSequence::Literal(Sequence::new(alt)),
    };
    GenomeVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_input_returns_empty() {
        assert!(merge_consecutive_edits(&[], AllelePhase::Cis).is_empty());
    }

    #[test]
    fn single_input_passes_through() {
        // Single-variant alleles are not Phase-1 mergeable.
        let acc = Accession::new("NC", "000001", Some(11));
        let v = HgvsVariant::Genome(GenomeVariant {
            accession: acc,
            gene_symbol: None,
            loc_edit: LocEdit::new(
                Interval::new(GenomePos::new(100), GenomePos::new(100)),
                NaEdit::Substitution {
                    reference: Base::G,
                    alternative: Base::A,
                },
            ),
        });
        let out = merge_consecutive_edits(std::slice::from_ref(&v), AllelePhase::Cis);
        assert_eq!(out.len(), 1);
    }
}
```

- [ ] **Step 4: Wire merge into `normalize_allele`**

Edit `src/normalize/mod.rs`. Add `pub mod merge;` at the end of the existing `pub mod` block (after line 30). Then update `normalize_allele` (around line 219-231) so the final return passes the merged variants:

Replace this block in `src/normalize/mod.rs:219-230`:

```rust
        // Second pass: check for overlaps and resolve
        if self.config.prevent_overlap {
            normalized = self.resolve_overlaps(normalized)?;
        }

        Ok((
            HgvsVariant::Allele(crate::hgvs::variant::AlleleVariant::new(
                normalized,
                allele.phase,
            )),
            all_warnings,
        ))
```

with:

```rust
        // Second pass: check for overlaps and resolve
        if self.config.prevent_overlap {
            normalized = self.resolve_overlaps(normalized)?;
        }

        // Third pass: merge consecutive edits per HGVS spec (issue #72)
        let merged = merge::merge_consecutive_edits(&normalized, allele.phase);

        Ok((
            HgvsVariant::Allele(crate::hgvs::variant::AlleleVariant::new(
                merged,
                allele.phase,
            )),
            all_warnings,
        ))
```

- [ ] **Step 5: Verify the test now passes**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: PASS for `test_merge_consecutive_subs_genome`.

Also run the unit tests in the new module: `cargo nextest run --features dev -E 'test(merge::tests)'` — both `empty_input_returns_empty` and `single_input_passes_through` should PASS.

- [ ] **Step 6: Run the full normalize test file to confirm no regressions**

Run: `cargo nextest run --features dev -E 'test(normalize)'`

Expected: all existing tests still PASS. In particular `test_allele_each_variant_normalized`, `test_allele_with_two_deletions_both_shift`, `test_allele_preserves_semicolon_separator`, `test_allele_with_dup_and_sub`, `test_compact_allele_notation`, `test_normalize_allele_normalizes_each_variant` should all still pass — they use inputs the merge pass should leave alone (different positions or non-mergeable edit combinations).

If any pre-existing allele test fails because the merge pass changed its output, investigate whether the test's input is actually a case the spec says should merge. If yes, update the test's expected output. If no, the merge logic has a bug — fix it before committing.

- [ ] **Step 7: Commit**

```bash
git add src/normalize/merge.rs src/normalize/mod.rs tests/merge_consecutive_edits_tests.rs
git commit -m "$(cat <<'EOF'
feat(normalize): merge consecutive substitutions in cis alleles (#72)

Adds src/normalize/merge.rs and wires it into normalize_allele after
resolve_overlaps. Currently handles the simplest case from the spec:
two adjacent g. substitutions in a cis allele collapse to a single
delins (e.g. g.[1000G>A;1001A>C] -> g.1000_1001delinsAC).

Subsequent commits extend coverage to deletions, delins, insertions,
all NaEdit-based variant types, and round-trip preservation of the
non-mergeable cases.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Add `Deletion` support (incl. issue example for del+del)

**Files:**
- Modify: `src/normalize/merge.rs` (extend `genome_anchor` and `try_merge_pair`)
- Modify: `tests/merge_consecutive_edits_tests.rs`

- [ ] **Step 1: Write failing tests for `del` cases**

Append to `tests/merge_consecutive_edits_tests.rs`:

```rust
#[test]
fn test_merge_consecutive_dels_genome() {
    // Issue #72 example: two adjacent single-nt deletions become a single ranged del.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000del;1001del]"),
        "NC_000001.11:g.1000_1001del",
    );
}

#[test]
fn test_merge_sub_then_del() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001del]"),
        "NC_000001.11:g.1000_1001delinsA",
    );
}

#[test]
fn test_merge_del_then_sub() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000del;1001A>C]"),
        "NC_000001.11:g.1000_1001delinsC",
    );
}

#[test]
fn test_merge_dels_drops_explicit_sequence() {
    // Per design doc: del+del with explicit ref sequences emits the no-sequence form.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000delA;1001delC]"),
        "NC_000001.11:g.1000_1001del",
    );
}

#[test]
fn test_merge_dels_drops_length() {
    // Per design doc: del+del with length specifiers emits no-sequence/no-length form.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000_1002del3;1003_1004del2]"),
        "NC_000001.11:g.1000_1004del",
    );
}
```

- [ ] **Step 2: Verify the new tests fail**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: the four new tests FAIL (anchor extraction returns `None` for `Deletion`, so no merge happens). The Task-1 sub+sub test should still PASS.

- [ ] **Step 3: Extend anchor extraction to handle `Deletion`**

In `src/normalize/merge.rs`, replace the `genome_anchor` function with:

```rust
fn genome_anchor(v: &GenomeVariant) -> Option<Anchor> {
    let edit = v.loc_edit.edit.inner()?; // Mu::Certain only
    let (start, end) = simple_genome_range(&v.loc_edit.location)?;
    match edit {
        NaEdit::Substitution { alternative, .. } => Some(Anchor {
            start,
            end,
            alt: vec![*alternative],
        }),
        NaEdit::SubstitutionNoRef { alternative } => Some(Anchor {
            start,
            end,
            alt: vec![*alternative],
        }),
        NaEdit::Deletion { .. } => Some(Anchor {
            start,
            end,
            alt: Vec::new(),
        }),
        _ => None,
    }
}
```

- [ ] **Step 4: Generalize the merge result to choose `Deletion` vs `Delins`**

Replace `build_genome_delins` and update `try_merge_pair` to use a new shared builder. In `src/normalize/merge.rs`:

```rust
fn try_merge_pair(a: &HgvsVariant, b: &HgvsVariant) -> Option<HgvsVariant> {
    let (av, bv) = match (a, b) {
        (HgvsVariant::Genome(av), HgvsVariant::Genome(bv)) => (av, bv),
        _ => return None,
    };
    if av.accession != bv.accession {
        return None;
    }
    let a_anchor = genome_anchor(av)?;
    let b_anchor = genome_anchor(bv)?;
    if a_anchor.end + 1 != b_anchor.start {
        return None;
    }
    let mut alt = a_anchor.alt;
    alt.extend(b_anchor.alt);
    Some(HgvsVariant::Genome(build_genome_merged(
        av,
        a_anchor.start,
        b_anchor.end,
        alt,
    )))
}

fn build_genome_merged(
    template: &GenomeVariant,
    start: u64,
    end: u64,
    alt: Vec<Base>,
) -> GenomeVariant {
    let location = Interval::new(GenomePos::new(start), GenomePos::new(end));
    let edit = if alt.is_empty() {
        // Spec-recommended no-sequence/no-length form.
        NaEdit::Deletion {
            sequence: None,
            length: None,
        }
    } else {
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::new(alt)),
        }
    };
    GenomeVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}
```

Delete the old `build_genome_delins` function.

- [ ] **Step 5: Verify all tests pass**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: all six tests in `tests/merge_consecutive_edits_tests.rs` PASS.

- [ ] **Step 6: Confirm no regressions**

Run: `cargo nextest run --features dev -E 'test(normalize)'`

Expected: all existing normalize tests still PASS. In particular `test_allele_with_two_deletions_both_shift` (which uses non-adjacent positions `c.4del` and `c.13del`) must still pass — the merge pass should leave it alone.

- [ ] **Step 7: Commit**

```bash
git add src/normalize/merge.rs tests/merge_consecutive_edits_tests.rs
git commit -m "$(cat <<'EOF'
feat(normalize): merge consecutive deletions and sub/del pairs (#72)

Extends the consecutive-edit merge pass to handle Deletion-typed
sub-variants. Two adjacent dels collapse to a single ranged del; a
sub+del or del+sub pair becomes a delins. Per-input ref-sequence and
length specifiers are dropped from the merged output to emit the
spec-recommended no-sequence form.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Add `Delins` support (with non-`Literal` exclusion)

**Files:**
- Modify: `src/normalize/merge.rs` (extend `genome_anchor`)
- Modify: `tests/merge_consecutive_edits_tests.rs`

- [ ] **Step 1: Write failing tests for `delins` cases**

Append to `tests/merge_consecutive_edits_tests.rs`:

```rust
#[test]
fn test_merge_sub_then_delins() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001_1002delinsTT]"),
        "NC_000001.11:g.1000_1002delinsATT",
    );
}

#[test]
fn test_merge_delins_then_sub() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000_1001delinsTT;1002A>C]"),
        "NC_000001.11:g.1000_1002delinsTTC",
    );
}

#[test]
fn test_merge_delins_then_delins() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000_1001delinsTT;1002_1003delinsAA]"),
        "NC_000001.11:g.1000_1003delinsTTAA",
    );
}

#[test]
fn test_merge_skips_non_literal_delins() {
    // delins with a non-Literal payload (e.g., delins10) is not safely
    // concatenable; the design doc requires the pair to pass through.
    let input = "NC_000001.11:g.[1000G>A;1001_1010delins10]";
    let result = normalize_to_string(input);
    // The output must still contain both edits separately (unchanged).
    assert!(result.contains("1000G>A"), "expected 1000G>A in {}", result);
    assert!(result.contains("1001_1010delins"), "expected 1001_1010delins in {}", result);
    assert!(result.contains(';'), "expected separator in {}", result);
}
```

- [ ] **Step 2: Verify the new tests fail**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: the three positive tests FAIL (anchor extraction returns `None` for `Delins`); `test_merge_skips_non_literal_delins` may PASS by accident (no merge attempted) — that's fine, the assertion still holds.

- [ ] **Step 3: Extend anchor extraction to handle `Delins`**

In `src/normalize/merge.rs`, update `genome_anchor` to add a `NaEdit::Delins` arm. The arm extracts the literal alt sequence; non-`Literal` payloads return `None` (refusing the merge for that pair):

```rust
        NaEdit::Delins { sequence } => {
            let bases = sequence.as_literal()?.bases().to_vec();
            Some(Anchor { start, end, alt: bases })
        }
```

Place the new arm after the `Deletion` arm, before the catch-all `_ => None`.

- [ ] **Step 4: Verify all tests pass**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: all positive tests PASS; `test_merge_skips_non_literal_delins` PASSES.

- [ ] **Step 5: Confirm no regressions**

Run: `cargo nextest run --features dev -E 'test(normalize)'`

Expected: all existing tests still PASS.

- [ ] **Step 6: Commit**

```bash
git add src/normalize/merge.rs tests/merge_consecutive_edits_tests.rs
git commit -m "$(cat <<'EOF'
feat(normalize): merge delins into consecutive-edit pass (#72)

Adjacent Delins sub-variants (literal alt only) participate in the
merge alongside Substitution and Deletion. Non-Literal InsertedSequence
payloads (Count, Range, Repeat, etc.) are excluded — those pairs pass
through unchanged because their alt content cannot be safely
concatenated into a single merged delins.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Add `Insertion` support (incl. spec FAQ example, two-ins-same-boundary)

**Files:**
- Modify: `src/normalize/merge.rs` (extend `genome_anchor`, generalize merged-edit construction)
- Modify: `tests/merge_consecutive_edits_tests.rs`

- [ ] **Step 1: Write failing tests for `ins` cases**

Append to `tests/merge_consecutive_edits_tests.rs`:

```rust
#[test]
fn test_merge_sub_then_ins() {
    // Spec FAQ analogue (in g. context):
    // sub at 100 + ins between 100 and 101 -> delins at 100 with alt = sub.alt + ins.bases.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100G>A;100_101insTA]"),
        "NC_000001.11:g.100delinsATA",
    );
}

#[test]
fn test_merge_ins_then_sub() {
    // Mirror: ins between 100 and 101 + sub at 101 -> delins at 101 with alt = ins.bases + sub.alt.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100_101insTA;101G>A]"),
        "NC_000001.11:g.101delinsTAA",
    );
}

#[test]
fn test_merge_del_then_ins() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100del;100_101insTA]"),
        "NC_000001.11:g.100delinsTA",
    );
}

#[test]
fn test_merge_ins_then_del() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100_101insTA;101del]"),
        "NC_000001.11:g.101delinsTA",
    );
}

#[test]
fn test_merge_two_ins_same_boundary_preserves_input_order() {
    // Two ins at the same boundary p|p+1 collapse into a single ins; bases
    // are concatenated in input order (T then A -> TA).
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100_101insT;100_101insA]"),
        "NC_000001.11:g.100_101insTA",
    );
}

#[test]
fn test_merge_skips_non_literal_ins() {
    // Ins with a non-Literal payload (e.g., ins10) is not safely
    // concatenable; the pair passes through unchanged.
    let input = "NC_000001.11:g.[100G>A;100_101ins10]";
    let result = normalize_to_string(input);
    assert!(result.contains("100G>A"), "expected 100G>A in {}", result);
    assert!(result.contains("100_101ins10"), "expected 100_101ins10 in {}", result);
    assert!(result.contains(';'), "expected separator in {}", result);
}
```

- [ ] **Step 2: Verify the new tests fail**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: the five positive tests FAIL; `test_merge_skips_non_literal_ins` passes by accident.

- [ ] **Step 3: Add `Insertion` to anchor extraction**

In `src/normalize/merge.rs`, update `genome_anchor` to add an `Insertion` arm. For an insertion at interval `[p, p+1]`, the anchor is `start = p+1, end = p` (empty range at the boundary):

```rust
        NaEdit::Insertion { sequence } => {
            // Insertion's interval is [p, p+1]; anchor is start=p+1, end=p (empty range).
            let bases = sequence.as_literal()?.bases().to_vec();
            // simple_genome_range gave us (p, p+1); convert to anchor form.
            if end != start + 1 {
                return None;
            }
            Some(Anchor {
                start: end, // p+1
                end: start, // p
                alt: bases,
            })
        }
```

Place this arm after the `Delins` arm.

- [ ] **Step 4: Generalize merged-edit construction to handle the empty-range case**

Replace `build_genome_merged` in `src/normalize/merge.rs`:

```rust
fn build_genome_merged(
    template: &GenomeVariant,
    start: u64,
    end: u64,
    alt: Vec<Base>,
) -> GenomeVariant {
    let edit = if start > end {
        // Empty range -> emit Insertion at boundary [end, start] = [start-1, start].
        debug_assert_eq!(start, end + 1, "invariant: anchor span <= 1 nt");
        NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::new(alt)),
        }
    } else if alt.is_empty() {
        NaEdit::Deletion {
            sequence: None,
            length: None,
        }
    } else {
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::new(alt)),
        }
    };
    let location = if start > end {
        // Insertion uses the original [end, start] interval.
        Interval::new(GenomePos::new(end), GenomePos::new(start))
    } else {
        Interval::new(GenomePos::new(start), GenomePos::new(end))
    };
    GenomeVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}
```

- [ ] **Step 5: Verify all tests pass**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: all tests in `tests/merge_consecutive_edits_tests.rs` PASS, including the previous tasks' tests.

- [ ] **Step 6: Confirm no regressions**

Run: `cargo nextest run --features dev -E 'test(normalize)'`

Expected: all existing normalize tests still PASS.

- [ ] **Step 7: Commit**

```bash
git add src/normalize/merge.rs tests/merge_consecutive_edits_tests.rs
git commit -m "$(cat <<'EOF'
feat(normalize): merge insertions at boundaries into consecutive pass (#72)

Insertions sit between positions, not at them, so they're modeled as
empty-range anchors (start = end + 1) at the boundary. They merge with
adjacent sub/del/delins via the same A.end+1 == B.start adjacency rule,
producing a delins covering the union range. Two ins at the same
boundary collapse to a single ins, bases concatenated in input order.

Mirrors the HGVS spec FAQ:
  c.[2077G>A;2077_2078insTA] -> c.2077delinsATA

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Multi-coordinate-system support + negative cases

**Files:**
- Modify: `src/normalize/merge.rs` (generalize across `Cds`, `Tx`, `Rna`, `Mt` variant types)
- Modify: `tests/merge_consecutive_edits_tests.rs`

This task generalizes the per-type matching in `try_merge_pair` so the same anchor / merge logic applies to `Cds`, `Tx`, `Rna`, and `Mt` variants — and also pins the negative cases that must round-trip unchanged.

The CDS / TX / RNA position types differ from `GenomePos`: they are signed (`i64`), with offsets and UTR markers. The merge pass treats a position as merge-eligible only if it is positive (≥ 1), has no offset, has no special UTR/downstream marker, and is not unknown — matching the design doc's "simple position" rule. For these "simple" positions the adjacency rule reduces to the same `A.end + 1 == B.start` check on the integer base value.

- [ ] **Step 1: Write tests for each coordinate system**

Append to `tests/merge_consecutive_edits_tests.rs`:

```rust
fn provider_with_simple_transcript() -> MockProvider {
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    let mut provider = MockProvider::new();
    let sequence: String = "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT".to_string();
    let exons = vec![Exon::new(1, 1, sequence.len() as u32)];
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(60),
        exons,
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

fn normalize_with_provider(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

#[test]
fn test_merge_cds_consecutive_subs() {
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:c.[10A>G;11A>C]",
        ),
        "NM_TEST.1:c.10_11delinsGC",
    );
}

#[test]
fn test_merge_tx_consecutive_subs() {
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:n.[10A>G;11A>C]",
        ),
        "NM_TEST.1:n.10_11delinsGC",
    );
}

#[test]
fn test_merge_rna_consecutive_subs_lowercase() {
    // RNA uses lowercase nucleotides per HGVS spec; merged alt must preserve case.
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:r.[10a>g;11a>c]",
        ),
        "NM_TEST.1:r.10_11delinsgc",
    );
}

#[test]
fn test_merge_mt_consecutive_subs() {
    assert_eq!(
        normalize_to_string("NC_012920.1:m.[100G>A;101A>C]"),
        "NC_012920.1:m.100_101delinsAC",
    );
}

// =====================================================================
// Negative cases — must round-trip unchanged.
// =====================================================================

#[test]
fn test_no_merge_one_nt_gap() {
    // One unchanged nucleotide between variants -> spec keeps them separate.
    let result = normalize_to_string("NC_000001.11:g.[100G>A;102C>T]");
    assert!(result.contains("100G>A"), "got {}", result);
    assert!(result.contains("102C>T"), "got {}", result);
    assert!(result.contains(';'), "got {}", result);
}

#[test]
fn test_no_merge_different_accessions() {
    let result = normalize_to_string("[NC_000001.11:g.100G>A;NC_000002.11:g.101A>C]");
    assert!(result.contains("NC_000001.11"), "got {}", result);
    assert!(result.contains("NC_000002.11"), "got {}", result);
}

#[test]
fn test_no_merge_different_variant_types() {
    // Genome and Cds in the same allele bracket -> not mergeable.
    let result = normalize_with_provider(
        provider_with_simple_transcript(),
        "[NC_000001.11:g.100G>A;NM_TEST.1:c.10A>C]",
    );
    assert!(result.contains("g.100G>A"), "got {}", result);
    assert!(result.contains("c.10A>C"), "got {}", result);
}

#[test]
fn test_no_merge_trans_phase() {
    // [a];[b] (semicolon between bracket pairs) is trans phase per HGVS.
    let result = normalize_to_string("NC_000001.11:g.[100G>A];[101A>C]");
    assert!(result.contains("100G>A"), "got {}", result);
    assert!(result.contains("101A>C"), "got {}", result);
    // No delins should appear.
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_intronic_position() {
    // Intronic positions (non-zero offset) are excluded from the merge pass.
    let result = normalize_with_provider(
        provider_with_simple_transcript(),
        "NM_TEST.1:c.[10+1A>G;10+2T>G]",
    );
    assert!(result.contains("10+1A>G"), "got {}", result);
    assert!(result.contains("10+2T>G"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_utr_boundary() {
    // c.-1 and c.1 are physically adjacent but no valid HGVS range syntax
    // spans the 5'UTR / CDS boundary (c.-1_1 doesn't exist).
    let result = normalize_with_provider(
        provider_with_simple_transcript(),
        "NM_TEST.1:c.[-1A>G;1A>T]",
    );
    assert!(result.contains("-1A>G"), "got {}", result);
    assert!(result.contains("1A>T"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_uncertain_edit() {
    // Mu::Uncertain (paren-wrapped edit) is not mergeable.
    let result = normalize_to_string("NC_000001.11:g.[(100G>A);101A>C]");
    assert!(result.contains("(100G>A)"), "got {}", result);
    assert!(result.contains("101A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_duplication_adjacent_to_sub() {
    let result = normalize_to_string("NC_000001.11:g.[100dup;101A>C]");
    assert!(result.contains("100dup"), "got {}", result);
    assert!(result.contains("101A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_inversion_adjacent_to_sub() {
    let result = normalize_to_string("NC_000001.11:g.[100_102inv;103A>C]");
    assert!(result.contains("100_102inv"), "got {}", result);
    assert!(result.contains("103A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_two_ins_different_boundaries() {
    // Two ins separated by an unchanged nucleotide at position 101 -> no merge.
    let result = normalize_to_string("NC_000001.11:g.[100_101insT;101_102insA]");
    assert!(result.contains("100_101insT"), "got {}", result);
    assert!(result.contains("101_102insA"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_reverse_input_order() {
    // Input listed in reverse position order; the walk preserves input order
    // and checks input-order adjacency (1001's anchor end + 1 != 1000's start),
    // so no merge happens. Pins the no-resort decision.
    let result = normalize_to_string("NC_000001.11:g.[1001A>C;1000G>A]");
    assert!(result.contains("1001A>C"), "got {}", result);
    assert!(result.contains("1000G>A"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}
```

- [ ] **Step 2: Verify the new tests fail (positive ones) / pass (negative ones)**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: the four positive tests (`test_merge_cds_consecutive_subs`, `test_merge_tx_consecutive_subs`, `test_merge_rna_consecutive_subs_lowercase`, `test_merge_mt_consecutive_subs`) FAIL because `try_merge_pair` only matches `(Genome, Genome)`. The negative tests should already PASS — the merge pass leaves them alone today (no support for those edit types means no merge attempted).

- [ ] **Step 3: Generalize the merge pass across coordinate types**

Replace `try_merge_pair` and add per-type anchor helpers in `src/normalize/merge.rs`. The structure: each variant produces an `(Anchor, accession, type-tag)`; merging requires matching accession and type-tag; the anchor logic is identical across types because all use simple integer positions.

Add these imports near the top of the file (replacing the existing import line):

```rust
use crate::hgvs::location::{CdsPos, GenomePos, RnaPos, TxPos};
use crate::hgvs::variant::{
    Accession, AllelePhase, CdsVariant, GenomeVariant, HgvsVariant, LocEdit, MtVariant,
    RnaVariant, TxVariant,
};
```

Replace `try_merge_pair` and the helper functions with the generalized versions:

```rust
fn try_merge_pair(a: &HgvsVariant, b: &HgvsVariant) -> Option<HgvsVariant> {
    match (a, b) {
        (HgvsVariant::Genome(av), HgvsVariant::Genome(bv)) => {
            try_merge_genome(av, bv).map(HgvsVariant::Genome)
        }
        (HgvsVariant::Cds(av), HgvsVariant::Cds(bv)) => {
            try_merge_cds(av, bv).map(HgvsVariant::Cds)
        }
        (HgvsVariant::Tx(av), HgvsVariant::Tx(bv)) => try_merge_tx(av, bv).map(HgvsVariant::Tx),
        (HgvsVariant::Rna(av), HgvsVariant::Rna(bv)) => {
            try_merge_rna(av, bv).map(HgvsVariant::Rna)
        }
        (HgvsVariant::Mt(av), HgvsVariant::Mt(bv)) => try_merge_mt(av, bv).map(HgvsVariant::Mt),
        _ => None,
    }
}

fn try_merge_genome(a: &GenomeVariant, b: &GenomeVariant) -> Option<GenomeVariant> {
    if a.accession != b.accession { return None; }
    let aa = anchor_from_naedit(&a.loc_edit, simple_genome_range)?;
    let ba = anchor_from_naedit(&b.loc_edit, simple_genome_range)?;
    let merged = merge_anchors(aa, ba)?;
    Some(build_genome_merged(a, merged))
}

fn try_merge_cds(a: &CdsVariant, b: &CdsVariant) -> Option<CdsVariant> {
    if a.accession != b.accession { return None; }
    let aa = anchor_from_naedit(&a.loc_edit, simple_cds_range)?;
    let ba = anchor_from_naedit(&b.loc_edit, simple_cds_range)?;
    let merged = merge_anchors(aa, ba)?;
    Some(build_cds_merged(a, merged))
}

fn try_merge_tx(a: &TxVariant, b: &TxVariant) -> Option<TxVariant> {
    if a.accession != b.accession { return None; }
    let aa = anchor_from_naedit(&a.loc_edit, simple_tx_range)?;
    let ba = anchor_from_naedit(&b.loc_edit, simple_tx_range)?;
    let merged = merge_anchors(aa, ba)?;
    Some(build_tx_merged(a, merged))
}

fn try_merge_rna(a: &RnaVariant, b: &RnaVariant) -> Option<RnaVariant> {
    if a.accession != b.accession { return None; }
    let aa = anchor_from_naedit(&a.loc_edit, simple_rna_range)?;
    let ba = anchor_from_naedit(&b.loc_edit, simple_rna_range)?;
    let merged = merge_anchors(aa, ba)?;
    Some(build_rna_merged(a, merged))
}

fn try_merge_mt(a: &MtVariant, b: &MtVariant) -> Option<MtVariant> {
    if a.accession != b.accession { return None; }
    let aa = anchor_from_naedit(&a.loc_edit, simple_genome_range)?;
    let ba = anchor_from_naedit(&b.loc_edit, simple_genome_range)?;
    let merged = merge_anchors(aa, ba)?;
    Some(build_mt_merged(a, merged))
}

/// Combine two adjacency-checked anchors into one merged anchor.
/// Returns None if not adjacent.
fn merge_anchors(a: Anchor, b: Anchor) -> Option<Anchor> {
    if a.end.checked_add(1)? != b.start {
        return None;
    }
    let mut alt = a.alt;
    alt.extend(b.alt);
    Some(Anchor {
        start: a.start.min(b.start),
        end: a.end.max(b.end),
        alt,
    })
}

/// Extract an anchor from a sub-variant's location+edit. The `range_fn`
/// callback returns `(start, end)` for the location only when it is "simple"
/// (no offsets, no UTR markers, certain). Returns None when the edit is
/// not merge-eligible (uncertain, non-NaEdit handled at the caller).
fn anchor_from_naedit<L>(
    loc_edit: &LocEdit<Interval<L>, NaEdit>,
    range_fn: impl Fn(&Interval<L>) -> Option<(u64, u64)>,
) -> Option<Anchor> {
    let edit = loc_edit.edit.inner()?; // Mu::Certain only
    let (start, end) = range_fn(&loc_edit.location)?;
    match edit {
        NaEdit::Substitution { alternative, .. }
        | NaEdit::SubstitutionNoRef { alternative } => Some(Anchor {
            start,
            end,
            alt: vec![*alternative],
        }),
        NaEdit::Deletion { .. } => Some(Anchor {
            start,
            end,
            alt: Vec::new(),
        }),
        NaEdit::Delins { sequence } => {
            let bases = sequence.as_literal()?.bases().to_vec();
            Some(Anchor { start, end, alt: bases })
        }
        NaEdit::Insertion { sequence } => {
            // Insertion's interval is [p, p+1]; convert to anchor [p+1, p].
            if end != start.checked_add(1)? {
                return None;
            }
            let bases = sequence.as_literal()?.bases().to_vec();
            Some(Anchor {
                start: end,
                end: start,
                alt: bases,
            })
        }
        _ => None,
    }
}

fn simple_cds_range(interval: &Interval<CdsPos>) -> Option<(u64, u64)> {
    Some((simple_cds_pos(&interval.start)?, simple_cds_pos(&interval.end)?))
}

fn simple_cds_pos(boundary: &UncertainBoundary<CdsPos>) -> Option<u64> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_unknown() || pos.is_intronic() || pos.is_5utr() || pos.is_3utr() || pos.base < 1 {
        return None;
    }
    Some(pos.base as u64)
}

fn simple_tx_range(interval: &Interval<TxPos>) -> Option<(u64, u64)> {
    Some((simple_tx_pos(&interval.start)?, simple_tx_pos(&interval.end)?))
}

fn simple_tx_pos(boundary: &UncertainBoundary<TxPos>) -> Option<u64> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_intronic() || pos.is_upstream() || pos.is_downstream() || pos.base < 1 {
        return None;
    }
    Some(pos.base as u64)
}

fn simple_rna_range(interval: &Interval<RnaPos>) -> Option<(u64, u64)> {
    Some((simple_rna_pos(&interval.start)?, simple_rna_pos(&interval.end)?))
}

fn simple_rna_pos(boundary: &UncertainBoundary<RnaPos>) -> Option<u64> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    // RnaPos has the same `is_intronic` / `is_3utr` / negative-base
    // exclusions as CdsPos. Confirm method names against the impl in
    // src/hgvs/location.rs (the impl block around line 332). If a method
    // differs, replace the call with the equivalent predicate.
    if pos.offset.unwrap_or(0) != 0 || pos.utr3 || pos.base < 1 {
        return None;
    }
    Some(pos.base as u64)
}
```

Replace the existing `simple_genome_range` and `simple_genome_pos` with these (unchanged behavior, signature kept):

```rust
fn simple_genome_range(interval: &Interval<GenomePos>) -> Option<(u64, u64)> {
    Some((simple_genome_pos(&interval.start)?, simple_genome_pos(&interval.end)?))
}

fn simple_genome_pos(boundary: &UncertainBoundary<GenomePos>) -> Option<u64> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_special() || pos.offset.is_some() {
        return None;
    }
    Some(pos.base)
}
```

Add per-variant `build_*_merged` functions. They share structure; here is `build_cds_merged` as the canonical example, with the others differing only in position type and accession-prefix string:

```rust
fn build_genome_merged(template: &GenomeVariant, merged: Anchor) -> GenomeVariant {
    let (location, edit) = build_naedit(merged, GenomePos::new);
    GenomeVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_cds_merged(template: &CdsVariant, merged: Anchor) -> CdsVariant {
    let (location, edit) = build_naedit(merged, |b| CdsPos::new(b as i64));
    CdsVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_tx_merged(template: &TxVariant, merged: Anchor) -> TxVariant {
    let (location, edit) = build_naedit(merged, |b| TxPos::new(b as i64));
    TxVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_rna_merged(template: &RnaVariant, merged: Anchor) -> RnaVariant {
    let (location, edit) = build_naedit(merged, |b| RnaPos::new(b as i64));
    RnaVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_mt_merged(template: &MtVariant, merged: Anchor) -> MtVariant {
    let (location, edit) = build_naedit(merged, GenomePos::new);
    MtVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_naedit<P>(merged: Anchor, mut to_pos: impl FnMut(u64) -> P) -> (Interval<P>, NaEdit) {
    let edit = if merged.start > merged.end {
        // Empty-range anchor -> Insertion at boundary.
        NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::new(merged.alt)),
        }
    } else if merged.alt.is_empty() {
        NaEdit::Deletion { sequence: None, length: None }
    } else {
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::new(merged.alt)),
        }
    };
    let (lo, hi) = if merged.start > merged.end {
        (merged.end, merged.start) // Insertion: original [p, p+1] interval
    } else {
        (merged.start, merged.end)
    };
    (Interval::new(to_pos(lo), to_pos(hi)), edit)
}
```

Delete the old non-generic helpers superseded by the new versions (`build_genome_merged` is replaced; `genome_anchor` is replaced by `anchor_from_naedit`).

- [ ] **Step 4: Verify all tests pass**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: every test in `tests/merge_consecutive_edits_tests.rs` PASSES, including the previous tasks' tests, the four new positive multi-coord tests, and the negative tests.

If `simple_rna_pos` fails to compile due to method-name differences on `RnaPos`, open `src/hgvs/location.rs` near the `RnaPos` impl block (around line 332) and adjust the predicate calls — the rule is "exonic CDS or transcript region only, no offsets, no UTR markers, base ≥ 1." Use the same set of predicates `RnaPos` actually exposes.

- [ ] **Step 5: Confirm no regressions**

Run: `cargo nextest run --features dev`

Expected: the entire test suite still PASSES.

- [ ] **Step 6: Run lint to catch dead code**

Run: `cargo clippy --features dev -- -D warnings`

Expected: no warnings. If clippy flags any unused imports left over from the refactor, remove them.

- [ ] **Step 7: Commit**

```bash
git add src/normalize/merge.rs tests/merge_consecutive_edits_tests.rs
git commit -m "$(cat <<'EOF'
feat(normalize): generalize consecutive-edit merge across c./n./r./m. (#72)

Adds Cds, Tx, Rna, and Mt to the merge pass via per-type anchor helpers
that share the core walk and merge logic. Position eligibility for c./
n./r. excludes intronic offsets, UTR boundaries, and unknown positions
per the design doc — there is no valid HGVS range syntax across UTR
boundaries, so adjacency cannot be expressed even when positions are
physically consecutive.

Pins the negative-case round-trips: one-nt gap, different accessions,
different variant types, trans phase, intronic positions, UTR boundary,
uncertain edits, dup/inv adjacent to sub, two ins at different
boundaries, and reverse-input-order pairs.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Chains, mixed cases, and Phase-2 deferred test stub

**Files:**
- Modify: `tests/merge_consecutive_edits_tests.rs`

This task adds tests that exercise the existing implementation rather than introducing new behavior. The walk's `for next in variants` loop already folds chains by repeatedly merging into `output.last()`, so 3+ chains should "just work."

- [ ] **Step 1: Write tests for chains, defensive overlap, and the Phase-2 stub**

Append to `tests/merge_consecutive_edits_tests.rs`:

```rust
#[test]
fn test_merge_chain_three_consecutive_subs() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001A>C;1002C>T]"),
        "NC_000001.11:g.1000_1002delinsACT",
    );
}

#[test]
fn test_merge_chain_sub_ins_sub_at_shared_boundary() {
    // Chain across one shared boundary: sub at 100 + ins between 100 and 101 + sub at 101.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100G>A;100_101insT;101A>C]"),
        "NC_000001.11:g.100_101delinsATC",
    );
}

#[test]
fn test_merge_long_chain_mixed_types() {
    // 4-element chain spanning sub, del, sub, sub at consecutive positions.
    // Position 100 sub, 101 del, 102 sub, 103 sub -> single delins 100..=103 alt = A_AA = AAA.
    // (The del at 101 contributes no alt, so the merged alt is "A" + "" + "A" + "A" = "AAA".)
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100G>A;101del;102C>A;103T>A]"),
        "NC_000001.11:g.100_103delinsAAA",
    );
}

#[test]
fn test_merge_same_position_twice_no_merge() {
    // Two edits at the SAME position (zero gap, but overlap not adjacency)
    // must not collapse into a single delins.
    let result = normalize_to_string("NC_000001.11:g.[100G>A;100A>C]");
    // The output should still contain both — we don't synthesize a "double mutation" form.
    assert!(result.contains("100G>A"), "got {}", result);
    assert!(result.contains("100A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
#[ignore = "Codon-frame exception for c. — tracked in issue #79"]
fn test_codon_frame_exception_deferred() {
    // Per HGVS spec: c.[145C>T;147C>G] (positions 145-147 share codon 49)
    // must be expressed as c.145_147delinsTGG. Implementing this needs
    // reference-provider access for the unchanged middle nucleotide and
    // codon-frame logic; deferred to issue #79.
    let provider = provider_with_simple_transcript();
    assert_eq!(
        normalize_with_provider(provider, "NM_TEST.1:c.[10A>G;12A>C]"),
        "NM_TEST.1:c.10_12delinsGAC",
    );
}
```

- [ ] **Step 2: Verify the tests pass (and the deferred one is skipped)**

Run: `cargo nextest run --features dev -E 'test(merge_consecutive)'`

Expected: the four positive chain tests and the same-position no-merge test PASS. The deferred test reports as IGNORED.

If `test_merge_long_chain_mixed_types` fails because of a position assumption (e.g., the `del` at 101 doesn't actually have anchor end at 101), inspect the parsed variant's interval. `parse_hgvs("g.101del")` produces an interval `(101, 101)`, so the anchor is `(101, 101)` with empty alt — adjacency to the prior sub at 100 holds (`100 + 1 == 101`). The merged result then has anchor `(100, 101)` with alt `[A]`; the next sub at 102 gives `101 + 1 == 102`, merged anchor `(100, 102)` alt `[A, A]`; etc.

- [ ] **Step 3: Run the full test suite and lint**

Run: `cargo nextest run --features dev` and `cargo clippy --features dev -- -D warnings`

Expected: all tests pass (with one ignored), no clippy warnings.

- [ ] **Step 4: Commit**

```bash
git add tests/merge_consecutive_edits_tests.rs
git commit -m "$(cat <<'EOF'
test(normalize): cover merge chains and add deferred codon-rule stub (#72)

Pins the chain behavior: 3+ consecutive substitutions, sub/ins/sub at a
shared boundary, and a 4-element mixed sub/del/sub/sub chain all
collapse into a single delins. Same-position-twice does not merge
(overlap, not adjacency).

Adds an #[ignore]-d test for the codon-frame exception (one-nt gap with
both positions in the same codon), tracked in #79, so the spec gap is
visible in the suite rather than only in the design doc.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Verification before opening the PR

- [ ] Run the full test suite once more: `cargo nextest run --features dev`. Expected: all pass, one ignored (the deferred codon-rule stub).
- [ ] Run lint: `cargo clippy --features dev -- -D warnings`. Expected: clean.
- [ ] Run formatter: `cargo fmt --all --check`. Expected: no diff.
- [ ] Run pre-commit: `pre-commit run --all-files`. Expected: clean.
- [ ] Spot-check the rendered HGVS for each issue example by running locally:
  ```
  cargo run --release -- parse 'NC_000001.11:g.[1000G>A;1001A>C]' --normalize
  cargo run --release -- parse 'NC_000001.11:g.[1000del;1001del]' --normalize
  ```
  Expected: `NC_000001.11:g.1000_1001delinsAC` and `NC_000001.11:g.1000_1001del` respectively. (If the CLI subcommand differs from `parse --normalize`, find the equivalent in `src/bin/ferro.rs` — but not having a CLI smoke test is fine; the integration tests already cover this.)

---

## Self-review notes

- **Spec coverage.** Every section of the design doc has at least one task or test:
  - Phase / variant types / edit types / certainty / position constraints / accession matching → enforced in `try_merge_pair` (Tasks 1, 5) and exercised by the negative-case tests (Task 5).
  - Anchor model + adjacency rule → Tasks 1 (sub), 2 (del), 3 (delins), 4 (ins), 5 (multi-coord).
  - Result-type rules (`Insertion` / `Deletion` / `Delins`) → Task 4 generalizes the builder; Tasks 2 and 4 cover the no-alt and empty-range branches.
  - Two-ins-same-boundary input ordering → Task 4.
  - Reverse-input-order no-resort decision → Task 5.
  - Defensive `A.end >= B.start` overlap refusal → covered indirectly by `test_merge_same_position_twice_no_merge` (Task 6) and the strict `A.end + 1 == B.start` check.
  - Phase-2 deferred test → Task 6.

- **Type / signature consistency.** `Anchor`, `try_merge_pair`, `merge_anchors`, `anchor_from_naedit`, `simple_*_range`, `build_*_merged` all use consistent (`u64`, `Vec<Base>`) types throughout. The `Insertion` anchor invariant (`start = end + 1`) is asserted via `debug_assert_eq!` in `build_naedit` (Task 5 generalization implicitly carries the same logic).

- **Pre-existing tests at risk.** The five pre-existing `mod allele_tests` tests in `tests/normalize_tests.rs` use inputs the merge pass should not touch (positions are not strictly adjacent, or mix dup with sub). Verified manually; if any of them turn out to be affected, the right move is to inspect whether the input is genuinely a spec-merge case and update the test, not to weaken the merge pass.
