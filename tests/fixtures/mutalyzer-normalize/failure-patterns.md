<!-- GENERATED — do not edit; regenerate via `cargo run --features dev --example generate_conformance_summary` -->

# mutalyzer-normalize conformance summary

Generated from `cases.json` — do not edit by hand; regenerate with the example above. Every row below is a tracked disposition. The live divergence set against full reference data is emitted only by the manifest run and is never committed (it is non-hermetic).

## Root-cause clusters

### RefSeqGene transcript selector (gene-symbol → NM_)

Spec: `background/refseq.md L38-42, L138-139`

On an `NG_/NC_/LRG_` reference ferro preserves the input's gene-symbol selector; the transcript-accession (NM_) form mutalyzer emits is spec-preferred. Valid HGVS, convergence tracked in #500. (The NM_(GENE) policy from #121 is unaffected.)

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC` | normalized | improvement | — | #500 |
| `NG_008939.1(PCCB_v001):c.156_157ins[GTCCTGTGCTC;ATTATCTGGC]` | normalized | improvement | — | #500 |
| `NG_008939.1(PCCB_v001):c.156_157ins[GTCCTGTGCTCATTATCTGGC]` | normalized | improvement | — | #500 |
| `NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC` | normalized | improvement | — | #500 |
| `NG_008939.1(PCCB_v001):c.156_161delins[GTCCTGTGCT;CATTATCTGGC]` | normalized | improvement | — | #500 |
| `NG_008939.1(PCCB_v001):c.156_161delins[GTCCTGTGCTCATTATCTGGC]` | normalized | improvement | — | #500 |
| `NG_008939.1:c.155_157del3` | normalized | improvement | — | #500 |
| `NG_008939.1:c.155_157delAAC` | normalized | improvement | — | #500 |
| `NG_008939.1:c.274_275inv` | normalized | improvement | — | #500 |
| `NG_012337.1(10683):c.274G>T` | normalized | improvement | — | #500 |
| `NG_012337.1(SDHD):c.274G>T` | normalized | improvement | — | #500 |
| `NG_012337.1(SDHD_v001):c.274G>T` | normalized | improvement | — | #500 |
| `NG_012337.1(TIMM8B):c.12_15delCAGC` | normalized | improvement | — | #500 |
| `NG_012337.1(TIMM8B):c.12_15delCAGCinsTTT` | normalized | improvement | — | #500 |
| `NG_012337.1(TIMM8B):c.12_15dupCAGC` | normalized | improvement | — | #500 |
| `NG_012337.1(TIMM8B):c.12delC` | normalized | improvement | — | #500 |
| `NG_012337.1(TIMM8B):c.12delCinsAT` | normalized | improvement | — | #500 |
| `NG_012337.1(TIMM8B):c.12dupC` | normalized | improvement | — | #500 |
| `NG_012337.1(TIMM8B_v001):c.12_13insGATC` | normalized | improvement | — | #500 |
| `NG_012337.1(TIMM8B_v001):c.12_13ins[GATC]` | normalized | improvement | — | #500 |
| `NG_012337.1(TIMM8B_v001):c.12_13ins[TTT;GATC]` | normalized | improvement | — | #500 |
| `NG_012772.1(BRCA2_v001):c.622_672del` | normalized | improvement | — | #500 |
| `NG_012772.1(BRCA2_v001):c.622_674del` | normalized | improvement | — | #500 |
| `NG_012772.1(BRCA2_v001):c.632_681del` | normalized | improvement | — | #500 |

### Bare NP_ protein reference

Spec: `recommendations/protein/*`

A p. description references the protein accession alone (NP_); mutalyzer's genomic-context(NP_) wrapper has no spec counterpart, so ferro's bare-NP_ output is the spec-correct value.

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_007485.1(NM_000077.4):c.161_162delTGinsATCCC` | protein_description | spec_citation | — | — |
| `NG_007485.1(NM_000077.4):c.161_162delTGins[ATCCC]` | protein_description | spec_citation | — | — |
| `NG_007485.1(NM_000077.4):c.161_162delinsATCCC` | protein_description | spec_citation | — | — |
| `NG_007485.1(NM_000077.4):c.161_162delins[ATCCC]` | protein_description | spec_citation | — | — |
| `NG_007485.1(NM_000077.4):c.161_162insATC` | protein_description | spec_citation | — | — |
| `NG_007485.1(NM_000077.4):c.161_162ins[ATC]` | protein_description | spec_citation | — | — |
| `NG_007485.1(NM_000077.4):c.161_163del` | protein_description | spec_citation | — | — |
| `NG_007485.1(NM_058195.3):c.141_142del` | protein_description | spec_citation | — | — |
| `NG_008835.1(NM_022153.2):c.568del` | protein_description | spec_citation | — | — |
| `NM_000143.3:c.1531T>G` | protein_description | spec_citation | — | — |
| `NM_000193.2:c.108_109del2insG` | protein_description | spec_citation | — | — |
| `NM_001199.3:c.2188dup` | protein_description | spec_citation | — | — |
| `NM_003002.2:c.273del` | protein_description | spec_citation | — | — |
| `NM_003002.2:c.274del` | protein_description | spec_citation | — | — |
| `NM_003002.4:c.169_170insATA` | protein_description | spec_citation | — | — |
| `NM_003002.4:c.206_210delins190_220inv` | protein_description | spec_citation | — | — |

### Protein initiation codon → p.(Met1?)

Spec: `recommendations/protein/substitution.md:51, deletion.md:62`

A variant reaching the initiation codon (CDS 1-3) has an unpredictable consequence, reported p.(Met1?). ferro derives it from the canonical normalized variant (#504, #512).

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NM_000143.3:c.-1_1insCAT` | protein_description | spec_citation | — | — |
| `NM_000143.3:c.1_2insCAT` | protein_description | spec_citation | — | — |

### A substitution consequence is never a frameshift

Spec: `recommendations/protein/substitution.md:5, frameshift.md:18`

A single-base coding substitution preserves length and frame; a p.…fs… consequence for a substitution is invalid against a fixed reference. ferro emits the in-frame missense against the canonical RefSeq transcript.

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_009299.1(NM_017668.3):c.41A>C` | protein_description | spec_citation | — | — |

### Gene-symbol selector policy (#121)

Spec: `background/refseq.md (NM_(GENE) form)`

ferro's terminal #121 policy preserves the NM_(GENE) gene-symbol selector; both forms are spec-allowed (an accepted divergence, not a convergence target).

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_012337.3(NM_003002.4):p.(Asp92Tyr)` | normalized | accepted_divergence | — | — |

### Ungrouped

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_012337.1(NM_003002.2):c.pterdel` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_003002.2):c.qterdel` | normalized | spec_citation | — | — |

## Disposition tallies

| axis | accepted_divergence | known_bug | improvement | spec_citation |
|---|---:|---:|---:|---:|
| normalized | 1 | 0 | 24 | 2 |
| protein_description | 0 | 0 | 0 | 19 |
