<!-- GENERATED — do not edit; regenerate via `cargo run --features dev --example generate_conformance_summary` -->

# biocommons-normalize conformance summary

Generated from `cases.json` — do not edit by hand; regenerate with the example above. Every row below is a tracked disposition. The live divergence set against full reference data is emitted only by the manifest run and is never committed (it is non-hermetic).

## Root-cause clusters

### Boundary-spanning del/delins accepted (5'UTR/CDS, CDS/3'UTR)

Spec: `recommendations/general.md, DNA/deletion.md`

Per the #253 audit, HGVS v21 does not forbid a boundary-spanning del with cross_boundaries=false; ferro accepts and is spec-compatible (e.g. NM_001166478.1:c.59_61del). A terminal accepted divergence from biocommons, seeded from a manifest-backed run (#325).

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NM_001166478.1:c.59_61del` | normalized | accepted_divergence | `NM_001166478.1:c.59_61del` | — |

### Spanning-dup canonicalization across the c.-1/c.1 boundary

Spec: `recommendations/DNA/duplication.md (edit-type priority)`

Per #404, an insertion spanning the CDS-start boundary canonicalizes to a spanning dup by edit-type priority (e.g. NM_000051.3:c.-2_-1insCA). A terminal accepted divergence, seeded from a manifest-backed run (#325).

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NM_000051.3:c.-2_-1insCA` | normalized | accepted_divergence | `NM_000051.3:c.-1_1dup` | — |

## Disposition tallies

| axis | accepted_divergence | known_bug | improvement | reference_unavailable | spec_citation |
|---|---:|---:|---:|---:|---:|
| normalized | 2 | 0 | 0 | 0 | 0 |
