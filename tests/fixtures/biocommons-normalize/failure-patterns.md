# Biocommons normalize — failure pattern triage

Companion to `cases.json` and `baseline-failures/normalized.txt`. Groups the
45 currently-known FAILs from `axis_normalized` (run against the
manifest at `/Volumes/scratch-00001/.../ferro/manifest.json`) into root-cause
patterns. Each pattern carries the HGVS v21.0 spec section that arbitrates the
case, the disposition (`ferro-bug` / `upstream-bug` / `accepted-divergence` /
`spec-ambiguous`), and linked ferro-hgvs issues / PRs.

Spec snapshot: `assets/hgvs-nomenclature/` at `6f85311989e76ead95d3547828f97ebaa3802e35` (v21.0.0).

## Pattern A — `NM_001166478.1` is missing from the reference bundle (9 cases × 2 directions, with dedupe ≈ 13 lines)

ferro's `MultiFastaProvider::from_manifest` does not bundle any version of
`NM_001166478.1` (gene DEFB133); the version-fallback path in
`resolve_name` only fires when the input has no version. cdot has the
exon alignment to `NC_000006.12` but ferro doesn't currently use that as a
transcript-base fallback. ferro returns the input unchanged for all
shift-dependent cases on this accession.

| Disposition | `ferro-bug` (infrastructure, not a normalize-logic bug) |
|---|---|
| Spec verdict | biocommons is correct (§general.md "3' rule"). |
| Linked | None directly. Closest: `tests/real_data_normalization_tests.rs` proves the manifest-loader works for other accessions. |
| Recommendation | New issue: *"MultiFastaProvider: use cdot exon-alignment as fallback for transcript base fetch when transcript FASTA is missing"*. |

Cases (input + direction + cross):

```text
NM_001166478.1:c.2_7delinsTTTAGA       (3prime+cross,  5prime+cross)
NM_001166478.1:c.30_31insT             (3prime+cross,  5prime+cross)
NM_001166478.1:c.31del                 (3prime+cross)
NM_001166478.1:c.34del                 (5prime+cross)
NM_001166478.1:c.35_36dup              (3prime+cross)
NM_001166478.1:c.35_36insT             (3prime+cross,  5prime+cross)
NM_001166478.1:c.36_37insTC            (3prime+cross,  5prime+cross)
NM_001166478.1:c.36_37insTCTCTC        (3prime+cross,  5prime+cross)
NM_001166478.1:c.59delG                (3prime+cross,  3prime+no-cross)
NM_001166478.1:c.59_61del              (3prime+no-cross — also Pattern E)
NM_001166478.1:c.61delG                (5prime+cross,  5prime+no-cross)
```

## Pattern B — Past-CDS-end positions (`NM_001001656.1`, CDS=945) (6 cases)

| Input | bio | ferro |
|---|---|---|
| `c.946G>C` (cross=false, cross=true) | rejects | accepts |
| `c.946dup` (cross=false, cross=true) | rejects | accepts (shifts past end to `c.*2dup`) |
| `c.935_946del` (cross=false, cross=true) | rejects | accepts (shifts past end to `c.935_*1del`) |

| Disposition | `ferro-bug` |
|---|---|
| Spec verdict | biocommons is correct (§general.md — positions must lie within the reference). |
| Linked | Not directly tracked. #266 (CLOSED, length-mismatch) covers a different validation gap; this needs a new SVA code (e.g. `W4xxx PositionPastEnd`). |
| Recommendation | New issue: *"normalize/validate: reject `c.`/`n.` positions past the CDS-end / transcript-end of the resolved reference"*. |

## Pattern C — Past-CDS-start positions / wrong-side shift (`NM_001001656.1:c.1del/dup`) (4 cases × 2 cross)

| Input | bio | ferro |
|---|---|---|
| `c.1del` (5prime+cross, 5prime+no-cross) | `c.1del` | `c.-2del` |
| `c.1dup` (5prime+cross, 5prime+no-cross) | `c.1dup` | `c.-2dup` |

ferro shifts a CDS-start variant into the 5' UTR under 5'-direction; biocommons
clamps. The CLI default is `cross_boundaries=false` which would have prevented
the cross — but these cases assert `cross=false` and ferro still crosses.

| Disposition | `ferro-bug` (boundary semantics) |
|---|---|
| Spec verdict | biocommons is correct (§general.md §"boundaries"; 3' rule does not authorize crossing the CDS-start with `cross_boundaries=false`). |
| Linked | Possibly related to #97 (CLOSED, "minus-strand 5'UTR single-base del collapses"). Different accession family (NM_001001656 is plus-strand). |
| Recommendation | New issue: *"normalize: cross_boundaries=false must clamp at CDS-start for 5'-direction shift"*. |

## Pattern D — `top-level delins` not canonicalized to `inv` (1 case)

| Input | bio | ferro |
|---|---|---|
| `NC_000009.11:g.36233991_36233992delCAinsAC` | `…delinsAC` | `…delCAinsAC` (unchanged) |

ferro keeps the explicit deleted bases instead of dropping them per
§DNA/delins.md: *"the recommendation is not to describe the variant as
`g.32386323delTinsGA`"*.

| Disposition | `ferro-bug` (canonicalization) |
|---|---|
| Spec verdict | biocommons is correct. |
| Linked | #75 (CLOSED, "Normalize same-base single-position delins as identity"), #73 (CLOSED, "Apply HGVS edit-type priority to simplify single-base delins"). Neither closes this case. |
| Recommendation | New issue: *"normalize: strip explicit deleted-sequence from delins where the alt is the canonical form"*. |

## Pattern E — Boundary-spanning del with `cross_boundaries=false` (1 case)

| Input | bio | ferro |
|---|---|---|
| `NM_001166478.1:c.59_61del` (cross=false) | rejects | accepts (pass-through; also Pattern A) |

| Disposition | `accepted-divergence` |
|---|---|
| Spec verdict | Per #253 audit: "HGVS Nomenclature v21 does not explicitly forbid these forms." Spec-compatible to accept. |
| Linked | #253 (CLOSED) — ferro deliberately accepts. |
| Recommendation | No action. Document as a known biocommons-deviation. |

## Pattern F — Ref-length mismatch (`NM_000059.3:c.7790delAAG`) (1 case)

1-position interval, 3-base deletion sequence. Biocommons rejects with
`HGVSError`.

| Disposition | `ferro-bug` (fix in flight) |
|---|---|
| Spec verdict | biocommons is correct (§general.md / §DNA/deletion.md: ref length must match position range length). |
| Linked | **#266 (CLOSED) + PR #272 (MERGED on main as a97a771)** wire `W3016 LengthMismatch`. The fix is now on `origin/main` but this branch was rebuilt off PR #323 which branched before #272 merged. After this branch rebases on top of merged main, the case should flip to `correctly-rejected` in strict mode. |
| Recommendation | Defer; re-check after rebasing onto post-#272 main. |

## Pattern G — Ferro panic on large delins (1 case, `NG_032871.1:g.32476_53457delinsAATTAAGGTATA`)

```text
assertion `left == right` failed: ref_bytes length must match HGVS interval span
  left: 15539
 right: 20982
note: thrown from src/normalize/mod.rs:2779
```

20,982 bp deleted, 12 bp inserted. ferro fetches a 15,539-base window
(probably the genomic alignment span) and then asserts that the byte length
equals the HGVS interval span (20,982) — they differ because the alignment
collapses gaps. The variant gets two upstream copies for the same input
(biocommons python copy-paste), so this panic surfaces twice in
`baseline-failures/normalized.txt`.

| Disposition | `ferro-bug` (panic in normalize path) |
|---|---|
| Spec verdict | Ambiguous on whether to reject biocommons-style (spec is silent on size limits, but the variant is well-formed). The panic itself is a ferro bug regardless. |
| Linked | Not directly tracked. |
| Recommendation | New issue: *"normalize: panic on `ref_bytes length must match HGVS interval span` for delins spanning alignment gaps in NG_*"*. |

## Pattern H — `issue_293` regression (`NG_029146.1:g.6494delG`) (1 case)

| Disposition | `ferro-bug` candidate (depends on reference mismatch resolution) |
|---|---|
| Spec verdict | Implementation-specific — biocommons' rejection in #293 was a Python attribute-assignment error, not a spec rule. Whether ferro should reject depends on whether position 6494 of NG_029146.1 is actually a G; if not, this is a RefSeqMismatch case (W3001). |
| Linked | None. |
| Recommendation | Investigate reference base at NG_029146.1:6494; if mismatch, expected behavior is W3001 (warn or reject depending on error_mode). |

## Pattern I — 3'-rule shift miss (various accessions in FASTA bundle) (~14 cases)

Cases where ferro DOES have reference data but the shift result differs
from biocommons. Examples:

| Input | bio | ferro |
|---|---|---|
| `NM_000051.3:c.14_15insT` (5prime+cross) | `c.14dup` | `c.15dup` |
| `NM_000051.3:c.-2_-1insCA` (3prime+no-cross) | `c.-1_1insAC` | `c.-1_1dup` |
| `NM_000051.3:c.-4_-3insA` (5prime+cross) | `c.-4dup` | `c.-3dup` |
| `NM_000051.3:c.1_2insCA` (5prime+cross) | `c.-1_1dup` | `c.-3_-2dup` |
| `NM_000051.3:c.9171_*1insA` (5prime+cross) | `c.9171dup` | `c.9170dup` |
| `NM_000051.3:c.*4_*5insT` (5prime+cross) | `c.*3dup` | `c.*4dup` |
| `NC_000001.10:g.1647893delinsCTTTCTT` (3prime+cross) | `g.1647895_1647900dup` | `g.1647893_1647894insTTTCTT` |
| `NC_000006.11:g.49917121_49917122insGA` (5prime+cross) | `g.49917121_49917122dup` | `g.49917122_49917123dup` |
| `NC_000006.11:g.49917122_49917123insA` (5prime+cross) | `g.49917123dup` | `g.49917127dup` |
| `NM_212556.2:c.1400_1401insAC` (3prime+cross, 3prime+no-cross) | `c.1401delinsACA` | `c.1401_*1insCA` |
| `NM_212556.2:c.1_2insCA` (5prime+cross) | `c.1delinsACA` | `c.-2_-1dup` |
| `NM_212556.2:c.1delinsCA` (5prime+cross, 3prime+cross) | `c.1delinsCA` | `c.-1_1insC` |
| `NM_212556.2:c.2_3insCAT` (5prime+cross) | `c.1delinsATCA` | `c.-1_1insCAT` |

The 5'-direction cases generally shift one or two positions further than
biocommons (ferro is more aggressive); the `NM_212556.2` cluster crosses
CDS-start when biocommons clamps (related to Pattern C). The
`NC_000001.10` case is a delins→dup canonicalization miss where ferro keeps
the insertion form.

| Disposition | `ferro-bug` (subgroups) |
|---|---|
| Spec verdict | biocommons is correct on the 3'-rule outcomes (§general.md). Some cross-CDS cases overlap with Pattern C. |
| Linked | #161 (CLOSED, "3'-shuffle terminator under-shifts deletions by one position") — adjacent area. #180 (CLOSED, "3' shifting in allele normalization") — adjacent. |
| Recommendation | Split into per-subgroup issues: (i) 5'-shuffle over-shift by one in single-base homopolymers, (ii) delins → dup canonicalization on multi-base inserts (`NC_000001.10`), (iii) intron/UTR boundary clamp (`NM_000051.3:c.*4_*5insT` etc.) — overlaps with Pattern C. |

## Summary

| Pattern | Cases | Disposition |
|---|---|---|
| A — NM_001166478.1 missing from bundle | 13 | ferro-bug (infrastructure) |
| B — past-CDS-end positions | 6 | ferro-bug |
| C — 5'-direction crosses CDS-start with cross=false | 4 | ferro-bug |
| D — top-level delins not canonicalized | 1 | ferro-bug |
| E — boundary-spanning del with cross=false | 1 | accepted-divergence (#253) |
| F — ref-length mismatch | 1 | ferro-bug (in flight, PR #272) |
| G — panic on large delins | 1 | ferro-bug (panic) |
| H — issue #293 regression | 1 | needs investigation |
| I — 3'-rule shift miss (subgroups) | ~14 | ferro-bug (multi-issue) |
| **Total** | **42–45** | (some inputs span multiple patterns) |

Plus 43 cases currently PASS — the regression baseline.

## Open issues to file

1. *MultiFastaProvider: use cdot exon-alignment as fallback for transcript base fetch* (Pattern A)
2. *normalize/validate: reject `c.`/`n.` positions past CDS-end / transcript-end* (Pattern B)
3. *normalize: cross_boundaries=false must clamp at CDS-start for 5'-direction shift* (Pattern C)
4. *normalize: strip explicit deleted-sequence from delins when alt is the canonical form* (Pattern D)
5. *normalize: panic on `ref_bytes length must match HGVS interval span` for delins spanning alignment gaps* (Pattern G)
6. Subgroup issues for Pattern I (see Recommendation)

After this branch rebases on top of post-#272 main, Pattern F should
demote from `ferro-bug` to `expected-pass` and the BRCA2 entry should
fall out of `baseline-failures/normalized.txt`.
