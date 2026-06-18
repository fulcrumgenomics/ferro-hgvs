<!-- GENERATED — do not edit; regenerate via `cargo run --features dev --example generate_conformance_summary` -->

# mutalyzer-normalize conformance summary

Generated from `cases.json` — do not edit by hand; regenerate with the example above. Every row below is a tracked disposition. The live divergence set against full reference data is emitted only by the manifest run and is never committed (it is non-hermetic).

## Root-cause clusters

### RefSeqGene transcript selector (gene-symbol → NM_)

Spec: `background/refseq.md L38-42, L138-139`

On an `NG_/NC_/LRG_` reference ferro preserves the input's gene-symbol selector; the transcript-accession (NM_) form mutalyzer emits is spec-preferred. Valid HGVS, convergence tracked in #500. (The NM_(GENE) policy from #121 is unaffected.)

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
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

_No seeded member rows yet (manifest-gated seeding, #325)._

### Compound-allele SHUFFLE_APPLIED (#499)

Spec: `recommendations/DNA (3' shifting); mutalyzer info model`

ferro emits a per-member SHUFFLE_APPLIED for a genuine 3' shift in a compound allele; mutalyzer reports only the allele-level ISORTEDVARIANTS (it sorts/collapses members, ferro preserves order). Both spec-defensible; ferro's shift info is correct. Accepted divergence per #499 (precedent #481).

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_009299.1(NM_002474.3):c.[310del;295G>A]` | infos | accepted_divergence | — | — |
| `NG_012337.1:g.[104_105insA;105_106insC;105del]` | infos | accepted_divergence | — | — |
| `NG_012337.1:g.[105_106insC;105del;104_105insA]` | infos | accepted_divergence | — | — |
| `NG_012337.1:g.[105del;104_105insA;105_106insC]` | infos | accepted_divergence | — | — |

### Coding repeat: codon-multiple exception

Spec: `recommendations/DNA/repeated.md L21-23`

On a coding DNA reference (c.), repeat notation unit[N] is permitted only for repeat units whose length is a multiple of 3 (cannot affect the reading frame); other expansions/contractions must be dup/ins/del. ferro emits the spec-mandated dup/ins/del for non-codon-aligned units in coding c.; mutalyzer's unit[N] form is invalid HGVS, so ferro's output is the spec-correct value.

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_012337.1(NM_012459.2):c.10CA[5]` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.3GC[5]` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.6C[4]` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.6delinsCCCC` | normalized | spec_citation | — | — |

### 3' rule across exon/intron boundaries (#670)

Spec: `background/numbering.md#DNAc (exception 3' rule NOTE); recommendations/DNA/deletion.md L51-64`

ferro shuffles genomic-context c./n. variants in spliced transcript space (cross_boundaries hardcoded false), which has no intronic bases, so it under-applies the 3' rule at exon/intron and intron/exon junctions. The spec applies the 3' rule across these junctions (suppressing it only at exon/exon junctions). The fix requires 3'-shifting in genomic coordinate space, which depends on the c↔g projection keystone (#642/#647/#616). Tracked in #670.

_No seeded member rows yet (manifest-gated seeding, #325)._

### Ungrouped

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `ENSG00000184937.16(ENST00000452863.10):c.9del` | normalized | reference_unavailable | — | #671 |
| `ENSG00000204370.13(ENST00000375549.8):c.100del` | normalized | reference_unavailable | — | #671 |
| `ENST00000375549.7:c.100del` | errors | accepted_divergence | — | — |
| `ENST00000375549.8:c.100del` | normalized | reference_unavailable | — | #671 |
| `ENST00000375549:c.100del` | normalized | reference_unavailable | — | #671 |
| `ENST00000452863.10:c.9del` | normalized | reference_unavailable | — | #671 |
| `LRG_199t1:c.11LRG_199t1:c.11` | errors | accepted_divergence | — | — |
| `LRG_1t1:52_153CAG[21]CAA[1]CAG[1]CCG[1]CCA[1]CCG[7]CCT[2]` | errors | accepted_divergence | — | — |
| `NG_007485.1(1787):n.204_205insATC` | errors | accepted_divergence | — | — |
| `NG_007485.1(CDKN2A):n.204_205insATC` | errors | accepted_divergence | — | — |
| `NG_007485.1(CDKN2A_v001):n.204_205insATC` | errors | accepted_divergence | — | — |
| `NG_008835.1(NM_001168390.2):c.*3186_*3187del` | genomic | accepted_divergence | — | — |
| `NG_008835.1(NM_022124.6):c.3304_3305del` | genomic | accepted_divergence | — | — |
| `NG_009930.1(NM_001099625.2):c.1010` | errors | accepted_divergence | — | — |
| `NG_009930.1(NM_001099625.2):n.110del` | errors | accepted_divergence | — | — |
| `NG_012337.1(DUMMYACCNO_9999.9):c.12_13insGATC` | errors | accepted_divergence | — | — |
| `NG_012337.1(NM_003002.2):c.pterdel` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_003002.2):c.qterdel` | normalized | spec_citation | — | — |
| `NG_012337.1:g.[4_5insT;4insA;4_5insA;5del]` | errors | accepted_divergence | — | — |
| `NG_012337.1:g.[4_5insT;4insA;4_5insA]` | errors | accepted_divergence | — | — |
| `NG_012337.1:g.[4del;4_5insT;4insA;4_5insA;5del]` | errors | accepted_divergence | — | — |
| `NG_012337.1:g.[4del;4_5insT;4insA;4_5insA]` | errors | accepted_divergence | — | — |
| `NG_012337.3(NM_003002):274+400G>T` | errors | accepted_divergence | — | — |
| `NG_012337.3(NM_003002.4):r.274+10G>T` | errors | accepted_divergence | — | — |
| `NM_003002.2:c.[100del;200_201insNM_003002.2:274+20]` | errors | accepted_divergence | — | — |

## Disposition tallies

| axis | accepted_divergence | known_bug | improvement | reference_unavailable | spec_citation |
|---|---:|---:|---:|---:|---:|
| errors | 16 | 0 | 0 | 0 | 0 |
| genomic | 2 | 0 | 0 | 0 | 0 |
| infos | 4 | 0 | 0 | 0 | 0 |
| normalized | 0 | 0 | 15 | 5 | 6 |
| protein_description | 0 | 0 | 0 | 0 | 19 |
