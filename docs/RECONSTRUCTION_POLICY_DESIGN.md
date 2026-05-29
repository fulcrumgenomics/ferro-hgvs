# Retire silent cdot-genome reconstruction — deposited-or-error transcript sourcing

Parent context: [#325](https://github.com/fulcrumgenomics/ferro-hgvs/issues/325) burn-down, [#471](https://github.com/fulcrumgenomics/ferro-hgvs/issues/471) (cdot-reconstruction mismatch), [#400](https://github.com/fulcrumgenomics/ferro-hgvs/issues/400) (Pattern J).

Builds on landed work: [#341](https://github.com/fulcrumgenomics/ferro-hgvs/pull/341)/#331 (the fallback this proposes to demote), #417 (deposited bases via supplemental FASTA), PR #482 (synthesis-risk diagnostic), PR #490 (strict exact-version resolution).

## Principle

ferro must not silently fabricate transcript **sequence** from the genome. The genome is the wrong source for transcript bases whenever the deposited transcript differs from it — by substitutions (invisible to cdot's `M`/`I`/`D` CIGAR), indels (which corrupt the synthesized length and therefore *all* downstream transcript coordinates), or non-reference haplotypes. No reconstruction logic can recover bases the genome does not contain.

The correct long-term model:

1. **Deposited transcript bases are the single source of truth** — the RefSeq/supplemental FASTA bundle (#417).
2. **cdot is for coordinates only** — exon alignment / genome↔transcript mapping — never for manufacturing sequence.
3. **When deposited bases are unavailable, fail loudly** — a structured error, not a fabricated answer. Genome reconstruction, if kept at all, is an explicit, off-by-default, clearly-labelled best-effort mode.
4. **Close the gap at the source** — `ferro prepare` should fetch deposited bases for every transcript the corpus/users need. Reconstruction exists only because the bundle was incomplete; the fix is to complete the bundle, not perfect the reconstruction.

## Why this reverses a prior decision (and must be explicit)

The cdot-genome synthesis fallback was added deliberately (#341/#331) to *produce answers* for transcripts missing from the FASTA bundle — "a possibly-wrong answer beats no answer." This RFC argues the opposite for a normalizer whose output is taken as authoritative: a silently-wrong answer (especially indel-driven coordinate corruption, which `axis_normalized` then misreads as a ferro divergence) is worse than a clean error. Because this changes ferro's public contract, the default-flip is called out below as an explicit, gated decision — not folded silently into a refactor.

## Current behaviour (problem statement)

`MultiFastaProvider::synthesize_transcript_from_cdot` pulls genome bases over each exon's genomic span `[e0,e1)` and concatenates (revcomp per-exon on the minus strand), **ignoring the CIGAR entirely**. Per exon `genome_len = e1−e0 = M+D` while `tx_len = e3−e2 = M+I`, so:

- **Substitutions** → emitted as the genome base (wrong base, no signal; this is #400/#473).
- **Indels** → wrong exon length → every downstream transcript coordinate is shifted (silent corruption).
- PR #482 *warns* on these but still returns the corrupt transcript; PR #490 adds an opt-in strict getter that refuses non-exact versions, but reconstruction is still the silent default for missing transcripts.

## Phased migration (each phase is independently shippable and the first phases are non-breaking)

### Phase 1 — resolution-policy seam (non-breaking)
Introduce an explicit policy on the provider/config, e.g.

```rust
enum TranscriptSourcePolicy {
    /// Deposited bases only; structured error if unavailable. (target default)
    DepositedOnly,
    /// Deposited bases, else genome reconstruction with a loud warning. (current behaviour)
    ReconstructBestEffort,
}
```

Default it to `ReconstructBestEffort` (today's behaviour — **no functional change**). `DepositedOnly` reuses PR #490's exact-resolution path plus a deposited-supplemental check, returning the existing/added structured error otherwise. Thread the policy through `MultiFastaProvider` construction (and `from_manifest`).

### Phase 2 — make reconstruction safe while it still runs
In `ReconstructBestEffort`, upgrade #482 from *warn* to *refuse* on **detectable corruption**: any exon with CIGAR `I`/`D` ops **or** a span-length mismatch (`e1−e0 ≠ e3−e2`) returns an error instead of emitting a length-corrupted transcript. (Span-length detection catches no-CIGAR net indels that the CIGAR-only check misses.) Substitutions remain undetectable from genome+cdot and stay a documented limitation. This is the "Option D" from the #471 discussion and is valuable on its own.

### Phase 3 — flip the default *(the gated decision)*
Switch the default to `DepositedOnly`. Transcripts missing deposited bases now error rather than reconstruct. The gate is quantitative: **bundle completeness = (distinct corpus transcripts whose deposited bases resolve from the FASTA bundle) ÷ (distinct corpus transcripts referenced by the gating suites), measured over a `DepositedOnly` run of the #478 conformance corpora in CI, must be 100%** — i.e. zero transcripts in the gating suites fall back to reconstruction or error for a missing deposit (see Phase 4). (100%, not a 99%-style threshold, because any falling-back transcript becomes a hard error under the new default and would break the gate.) Flipping also requires sign-off that the contract change is acceptable. Keep `ReconstructBestEffort` available behind explicit opt-in for callers who still want best-effort.

### Phase 4 — close the data gap upstream
Extend `ferro prepare` / the manifest pipeline to fetch deposited bases (supplemental FASTA, #417 mechanism) for every transcript the corpora and known users need — so `DepositedOnly` is the *complete*, not merely *correct*, answer. This is what actually demotes the #400/Pattern J rows (deposited `NM_001166478.1` makes ferro match biocommons).

## Interaction with the conformance harness (#478)
Under `DepositedOnly`, a transcript missing from the bundle yields a clean error, which the #478 annotation model records as a `known_bug`/`accepted_divergence` with a clear cause ("transcript not in bundle") rather than a phantom base-divergence. PR #490's strict resolver + this policy together are what let `axis_normalized` become a trustworthy CI gate.

## Risks / back-compat
- **Behaviour change at Phase 3** — callers relying on best-effort reconstruction must opt in. Mitigated by keeping the mode available and flipping only after the bundle is complete.
- **Bundle completeness is operational** — Phase 3's value depends on Phase 4; sequence them accordingly.
- **Substitutions remain unrecoverable** — out of scope for any genome-based approach; only deposited bases fix them. Documented, not solved here.

## Open decision for the maintainer
Phases 1–2 are safe to land now (additive + a corruption-refusal that only turns silent-wrong into loud-error). **Phase 3 (flip the default) is the contract-changing call** — proceed when the bundle is complete and the error-on-missing behaviour is agreed.
