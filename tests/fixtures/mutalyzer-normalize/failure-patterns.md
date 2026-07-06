<!-- GENERATED — do not edit; regenerate via `cargo run --features dev --example generate_conformance_summary` -->

# mutalyzer-normalize conformance summary

Generated from `cases.json` — do not edit by hand; regenerate with the example above. Every row below is a tracked disposition. The live divergence set against full reference data is emitted only by the manifest run and is never committed (it is non-hermetic).

## Root-cause clusters

### RefSeqGene transcript selector (gene-symbol → NM_)

Spec: `background/refseq.md L38-42, L138-139`

On an `NG_/NC_/LRG_` reference ferro preserves the input's gene-symbol selector; the transcript-accession (NM_) form mutalyzer emits is spec-preferred. Valid HGVS, convergence tracked in #500. (The NM_(GENE) policy from #121 is unaffected.)

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_008939.1:c.155_157del3` | normalized | improvement | — | #923 |
| `NG_008939.1:c.155_157delAAC` | normalized | improvement | — | #923 |
| `NG_008939.1:c.274_275inv` | normalized | improvement | — | #923 |
| `NG_012337.1(10683):c.274G>T` | normalized | improvement | — | #923 |

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
| `NG_008835.1(NM_022153.2):c.677-20_704+62del` | protein_description | spec_citation | — | — |
| `NG_008835.1(NM_022153.2):c.677-21_704+62del` | protein_description | spec_citation | — | — |
| `NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC` | protein_description | spec_citation | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins[GTCCTGTGCTC;ATTATCTGGC]` | protein_description | spec_citation | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins[GTCCTGTGCTCATTATCTGGC]` | protein_description | spec_citation | — | — |
| `NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC` | protein_description | spec_citation | — | — |
| `NG_008939.1(PCCB_v001):c.156_161delins[GTCCTGTGCT;CATTATCTGGC]` | protein_description | spec_citation | — | — |
| `NG_008939.1(PCCB_v001):c.156_161delins[GTCCTGTGCTCATTATCTGGC]` | protein_description | spec_citation | — | — |
| `NG_009299.1(NM_017668.3):c.[41A>C;250del]` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_003002.2):c.169T>A` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_003002.2):c.[274G>T;278A>G]` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.129G>T` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.130G>A` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.15_16insTGGAAGGTGGCGAGCCT` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.15_17dup` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.274delinsNG_012337.3(NM_003002.4):c.200_203` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.297_*1del` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.3GC[5]` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.4_5insGTA` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.4_6delinsGTA` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.5_6delinsTAG` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.6C[4]` | protein_description | spec_citation | — | — |
| `NG_012337.1(TIMM8B):c.12_15dupCAGC` | protein_description | spec_citation | — | — |
| `NG_012337.1(TIMM8B):c.12dupC` | protein_description | spec_citation | — | — |
| `NG_012337.1(TIMM8B_v001):c.12_13insGATC` | protein_description | spec_citation | — | — |
| `NG_012337.1(TIMM8B_v001):c.12_13ins[GATC]` | protein_description | spec_citation | — | — |
| `NG_012337.1(TIMM8B_v001):c.12_13ins[TTT;GATC]` | protein_description | spec_citation | — | — |
| `NG_012337.3(NM_003002.4):c.274delinsNG_012337.1(NM_012459.2):c.200_203` | protein_description | spec_citation | — | — |
| `NG_012772.1(BRCA2_v001):c.622_672del` | protein_description | spec_citation | — | — |
| `NG_012772.1(BRCA2_v001):c.622_674del` | protein_description | spec_citation | — | — |
| `NG_012772.1(BRCA2_v001):c.632_681del` | protein_description | spec_citation | — | — |
| `NM_000143.3:c.1531T>G` | protein_description | spec_citation | — | — |
| `NM_000193.2:c.108_109del2insG` | protein_description | spec_citation | — | — |
| `NM_000193.2:c.1388G>C` | protein_description | spec_citation | — | — |
| `NM_000193.2:c.1388_1389insC` | protein_description | spec_citation | — | — |
| `NM_001199.3:c.2188dup` | protein_description | spec_citation | — | — |
| `NM_003002.2:c.273del` | protein_description | spec_citation | — | — |
| `NM_003002.2:c.274del` | protein_description | spec_citation | — | — |
| `NM_003002.4:c.169_170insATA` | protein_description | spec_citation | — | — |
| `NM_003002.4:c.206_210delins190_220inv` | protein_description | spec_citation | — | — |
| `NM_003002.4:n.206_210del` | protein_description | spec_citation | — | — |
| `NM_024426.4:c.1107A>G` | protein_description | spec_citation | — | — |
| `NM_024426.4:c.1C>A` | protein_description | spec_citation | — | — |
| `NM_024426.4:c.1C>G` | protein_description | spec_citation | — | — |
| `NM_024426.4:c.1_4delinsATGA` | protein_description | spec_citation | — | — |

### Protein initiation codon → p.(Met1?)

Spec: `recommendations/protein/substitution.md:51, deletion.md:62`

A variant reaching the initiation codon (CDS 1-3) has an unpredictable consequence, reported p.(Met1?). ferro derives it from the canonical normalized variant (#504, #512).

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
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
| `NG_012337.1(NM_003002.2):c.[274+20C>T;400_401insNM_003002.4:100_102]` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.10CA[5]` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.3GC[5]` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.6C[4]` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.6delinsCCCC` | normalized | spec_citation | — | — |

### 3' rule across exon/intron boundaries (#670)

Spec: `background/numbering.md#DNAc (exception 3' rule NOTE); recommendations/DNA/deletion.md L51-64`

ferro shuffles genomic-context c./n. variants in spliced transcript space (cross_boundaries hardcoded false), which has no intronic bases, so it under-applies the 3' rule at exon/intron and intron/exon junctions. The spec applies the 3' rule across these junctions (suppressing it only at exon/exon junctions). The fix requires 3'-shifting in genomic coordinate space, which depends on the c↔g projection keystone (#642/#647/#616). Tracked in #670.

_No seeded member rows yet (manifest-gated seeding, #325)._

### coding_protein_descriptions transcript-set enumeration

Spec: `background/refseq.md L256, L259, L264`

On `coding_protein_descriptions` (expected ⊆ got, version- and gene-suffix-insensitive) ferro enumerates the latest curated version of every transcript cdot reports overlapping the locus (#710: collapse superseded versions, prefer curated NM_/NR_ over predicted XM_/XR_). Mutalyzer instead lists the neighbouring overlapping transcripts and partitions the input gene's own consequence onto its separate `protein_description` field — an info-model choice, not spec-mandated. `background/refseq.md` L256/L259 pose 'which transcript / which reference sequence should one use' as open questions and L264 endorses the latest build, so ferro's complete-current-curated-set enumeration is spec-valid; both representations are correct. Accepted divergence, not a convergence target — tightening would regress the #710 curation and cannot fix the rows anyway (the matcher is version-insensitive). Reference-coverage gaps where cdot lacks a transcript base mutalyzer used are tracked separately in #645.

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_007485.1(NM_000077.4):c.161_162delTGinsATCCC` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_007485.1(NM_000077.4):c.161_162delTGins[ATCCC]` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_007485.1(NM_000077.4):c.161_162delinsATCCC` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_007485.1(NM_000077.4):c.161_162delins[ATCCC]` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_007485.1(NM_000077.4):c.161_162insATC` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_007485.1(NM_000077.4):c.161_162ins[ATC]` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_007485.1(NM_000077.4):c.161_163del` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_007485.1(NM_058195.3):c.141_142del` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_007485.1:g.5350_5352del3` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_007485.1:g.5350_5352delCCA` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_008835.1(NM_022153.2):c.568del` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_009299.1(NM_017668.3):c.41A>C` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_009299.1(NM_017668.3):c.[41A>C;250del]` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_012337.1(NM_012459.2):c.129G>T` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_012337.1(NM_012459.2):c.130G>A` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_012337.1(TIMM8B):c.12_15dupCAGC` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_012337.1(TIMM8B):c.12dupC` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_012337.1:g.13124_13125del` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_012337.1:g.4812_4813insTAC` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_012337.1:g.7125G>T` | coding_protein_descriptions | accepted_divergence | — | — |

### Transcript 5'/3' flank not c.-numberable (pter/qter)

Spec: `background/refseq.md L45-46`

A coding/non-coding DNA reference does not contain the gene's 5'/3' flanking sequences and can not be used to describe variants there (refseq.md L45-46). pter/qter denote those flanks, so they are not c.-numberable — as a coordinate or inside a cross-reference/delins payload. ferro declines; mutalyzer extrapolates the flank and resolves the bases. ferro's refusal is spec-correct (#758/#537).

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `LRG_24:g.5525_5532delinsNM_003002.2:c.*835_qter` | normalized | accepted_divergence | — | — |
| `LRG_24:g.5525_5532delinsNM_003002.2:c.pter_-51` | normalized | accepted_divergence | — | — |
| `LRG_24t1:c.126_133delins[NM_003002.2:c.pter_-51;NM_003002.2:c.*835_qter]` | genomic | accepted_divergence | — | — |
| `LRG_24t1:c.126_133delins[NM_003002.2:c.pter_-51;NM_003002.2:c.*835_qter]` | normalized | accepted_divergence | — | — |
| `LRG_24t1:c.pter_qterdelins[pter_qter;NM_003002.2:c.*835_qter]` | genomic | accepted_divergence | — | — |
| `LRG_24t1:c.pter_qterdelins[pter_qter;NM_003002.2:c.*835_qter]` | normalized | accepted_divergence | — | — |
| `LRG_24t1:c.pter_qterdelinspter_qter` | genomic | accepted_divergence | — | — |
| `LRG_24t1:c.pter_qterdelinspter_qter` | normalized | accepted_divergence | — | — |

### RNA predicted allele bracket placement

Spec: `recommendations/RNA/alleles.md`

A predicted (uncertain) cis RNA allele places the '( )' predicted wrapper inside the allele bracket (r.[(a;b)], e.g. LRG_199t1:r.[(578c>u;1339a>g;1680del)]); mutalyzer wraps the whole bracket (r.([a;b])), which has no spec counterpart. ferro's inside-bracket form is the spec-correct value (#693).

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_012337.1(NM_003002.2):c.[274G>T;278A>G]` | rna_description | spec_citation | — | — |

### n. projection onto an overlapping transcript outside its exonic span

Spec: `background/numbering.md#DNAn (3' downstream offset)`

The n. corpus rows expect a projection onto the overlapping single-exon non-coding transcript NR_024274.1 (CDKN2A-AS1), where the variant lies 3' downstream of its single exon and mutalyzer addresses it as n.616+offset. ferro's transcript fan-out enumerates only transcripts the variant overlaps, so it emits valid n. forms on the overlapping NM_ isoforms but declines NR_024274.1 ('variant does not overlap transcript'). Spec tension: numbering.md L43-44 bars describing variants beyond a bare transcript's boundaries, so a bare NR_024274.1:n.616+offset would be invalid; but mutalyzer (and the corpus) frame it under its NG_007485.1 genomic parent (NG_007485.1(NR_024274.1):n.616+offset), the same genomic-context construction numbering.md L65 endorses (e.g. NG_012232.1(NM_004006.2):r.186+1_186+4), where the parent supplies the flanking sequence. Under that framing the form is valid and ferro should produce it -> known_bug, not accepted_divergence. Cross-isoform enumeration scope tracked by #646 (core landed in #647); ferro should extend enumeration/addressing to the downstream-of-last-exon case when the genomic parent is present. XPASS-guarded: if #646 lands and ferro starts matching, the annotation fails loudly and the row demotes.

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_007485.1(NM_000077.4):c.161_162delTGinsATCCC` | noncoding | known_bug | — | #921 |
| `NG_007485.1(NM_000077.4):c.161_162delTGins[ATCCC]` | noncoding | known_bug | — | #921 |
| `NG_007485.1(NM_000077.4):c.161_162delinsATCCC` | noncoding | known_bug | — | #921 |
| `NG_007485.1(NM_000077.4):c.161_162delins[ATCCC]` | noncoding | known_bug | — | #921 |
| `NG_007485.1(NM_000077.4):c.161_162insATC` | noncoding | known_bug | — | #921 |
| `NG_007485.1(NM_000077.4):c.161_162ins[ATC]` | noncoding | known_bug | — | #921 |
| `NG_007485.1(NM_058195.3):c.141_142del` | noncoding | known_bug | — | #921 |
| `NG_007485.1:g.5479_5480insACT` | noncoding | known_bug | — | #921 |

### Inserted inversion (reverse complement) not applied by mutalyzer

Spec: `HGVS §Inserted inversion (reverse complement)`

An `<range>inv` segment in an ins/delins payload inserts the reverse complement (inverted copy) of the range's bases (DNA/insertion.md L52-62; DNA/inversion.md L5, L48-49, L53-62). ferro applies the inversion and normalizes the delins; mutalyzer leaves the segment forward (inversion not applied). ferro's output is the spec-correct value. #854.

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NG_008939.1:g.5207_5212delins[4300_4309inv;4310_4320]` | genomic | spec_citation | — | — |

### Ungrouped

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `ENSG00000184937.16(ENST00000452863.10):c.9del` | normalized | reference_unavailable | — | #933 |
| `ENSG00000204370.13(ENST00000375549.8):c.100del` | normalized | reference_unavailable | — | #933 |
| `ENST00000375549.7:c.100del` | errors | accepted_divergence | — | — |
| `ENST00000375549.8:c.100del` | normalized | reference_unavailable | — | #933 |
| `ENST00000375549:c.100del` | normalized | reference_unavailable | — | #933 |
| `ENST00000452863.10:c.9del` | normalized | reference_unavailable | — | #933 |
| `LRG_199t1:c.11LRG_199t1:c.11` | errors | accepted_divergence | — | — |
| `LRG_199t1:c.11LRG_199t1:c.11[10]` | normalized | accepted_divergence | — | — |
| `LRG_199t1:c.11NG_012337.3(NM_003002.4):c.14[10]` | normalized | accepted_divergence | — | — |
| `LRG_199t1:c.235_237delinsTAT` | normalized | known_bug | — | #920 |
| `LRG_1t1:52_153CAG[21]CAA[1]CAG[1]CCG[1]CCA[1]CCG[7]CCT[2]` | errors | accepted_divergence | — | — |
| `LRG_303:g.4_5insGT[5]` | normalized | known_bug | — | #920 |
| `LRG_303:g.6932_6933insAGCAACGTGATCGCCTCCCTCACCT[2]` | normalized | known_bug | — | #920 |
| `NG_007485.1(1787):n.204_205insATC` | errors | accepted_divergence | — | — |
| `NG_007485.1(CDKN2A):n.204_205insATC` | errors | accepted_divergence | — | — |
| `NG_007485.1(CDKN2A_v001):n.204_205insATC` | errors | accepted_divergence | — | — |
| `NG_007485.1(NM_058195.3):n.204_205insATC` | normalized | accepted_divergence | — | — |
| `NG_007485.1:g.5479_5480insACT` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_008835.1(NM_001168390.2):c.*3186_*3187del` | genomic | accepted_divergence | — | — |
| `NG_008835.1(NM_001168390.2):c.*3186_*3187del` | normalized | spec_citation | — | — |
| `NG_008835.1(NM_022124.6):c.3304_3305del` | genomic | accepted_divergence | — | — |
| `NG_008835.1(NM_022153.2):c.82del` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_008835.1:g.320802_320803del` | genomic | accepted_divergence | — | — |
| `NG_008835.1:g.320802_320803del` | normalized | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins180_188` | genomic | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins180_188` | normalized | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins180_188` | protein_description | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins180_188inv` | genomic | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins180_188inv` | normalized | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins180_188inv` | protein_description | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins[180_188]` | genomic | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins[180_188]` | normalized | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins[180_188]` | protein_description | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins[180_188inv]` | genomic | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins[180_188inv]` | normalized | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_157ins[180_188inv]` | protein_description | accepted_divergence | — | — |
| `NG_008939.1(PCCB_v001):c.156_161delins180_188` | genomic | known_bug | — | #922 |
| `NG_008939.1(PCCB_v001):c.156_161delins180_188` | normalized | known_bug | — | #922 |
| `NG_008939.1(PCCB_v001):c.156_161delins180_188inv` | genomic | known_bug | — | #922 |
| `NG_008939.1(PCCB_v001):c.156_161delins180_188inv` | normalized | known_bug | — | #922 |
| `NG_008939.1(PCCB_v001):c.156_161delins[180_188]` | genomic | known_bug | — | #922 |
| `NG_008939.1(PCCB_v001):c.156_161delins[180_188]` | normalized | known_bug | — | #922 |
| `NG_008939.1(PCCB_v001):c.156_161delins[180_188inv]` | genomic | known_bug | — | #922 |
| `NG_008939.1(PCCB_v001):c.156_161delins[180_188inv]` | normalized | known_bug | — | #922 |
| `NG_008939.1:c.155_157del3` | genomic | accepted_divergence | — | — |
| `NG_008939.1:c.155_157del3` | protein_description | accepted_divergence | — | — |
| `NG_008939.1:c.155_157delAAC` | genomic | accepted_divergence | — | — |
| `NG_008939.1:c.155_157delAAC` | protein_description | accepted_divergence | — | — |
| `NG_008939.1:c.274_275inv` | genomic | accepted_divergence | — | — |
| `NG_008939.1:c.274_275inv` | protein_description | accepted_divergence | — | — |
| `NG_008939.1:g.5207_5208ins4300_4320inv` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_008939.1:g.5207_5208ins[4300_4320inv]` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_008939.1:g.5207_5212delins[GTCCTGTGCT;4310_4320inv]` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_009299.1(NM_002474.3):c.[310del;295G>A]` | normalized | spec_citation | — | — |
| `NG_009299.1(NM_017668.3):c.33_35CAA[5]` | genomic | accepted_divergence | — | — |
| `NG_009299.1(NM_017668.3):c.33_35CAA[5]` | normalized | known_bug | — | #920 |
| `NG_009299.1(NM_017668.3):c.41A>C` | genomic | accepted_divergence | — | — |
| `NG_009299.1(NM_017668.3):c.[250del;41>CA]` | infos | accepted_divergence | — | — |
| `NG_009299.1(NM_017668.3):c.[250del;41>CA]` | normalized | accepted_divergence | — | — |
| `NG_009299.1(NM_017668.3):c.[250del;41A>C]` | infos | spec_citation | — | — |
| `NG_009299.1(NM_017668.3):c.[250del;41A>C]` | normalized | spec_citation | — | — |
| `NG_009299.1(NM_017668.3):c.[41>CA;250del]` | genomic | accepted_divergence | — | — |
| `NG_009299.1(NM_017668.3):c.[41>CA;250del]` | infos | accepted_divergence | — | — |
| `NG_009299.1(NM_017668.3):c.[41>CA;250del]` | normalized | accepted_divergence | — | — |
| `NG_009299.1(NM_017668.3):c.[41A>C;250del]` | genomic | accepted_divergence | — | — |
| `NG_009299.1(NM_017668.3):c.[41A>C;250del]` | normalized | spec_citation | — | — |
| `NG_009930.1(NM_001099625.2):c.1010` | errors | accepted_divergence | — | — |
| `NG_009930.1(NM_001099625.2):n.110del` | errors | accepted_divergence | — | — |
| `NG_012337.1(10683):c.274G>T` | genomic | accepted_divergence | — | — |
| `NG_012337.1(DUMMYACCNO_9999.9):c.12_13insGATC` | errors | accepted_divergence | — | — |
| `NG_012337.1(NM_003002.2):c.1A>G` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_003002.2):c.1A>T` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_003002.2):c.1_4delinsTTGA` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_003002.2):c.276C>T` | protein_description | spec_citation | — | — |
| `NG_012337.1(NM_003002.2):c.[274;600]` | genomic | accepted_divergence | — | — |
| `NG_012337.1(NM_003002.2):c.[274;600]` | normalized | known_bug | — | #920 |
| `NG_012337.1(NM_003002.2):c.[274=;275del]` | genomic | accepted_divergence | — | — |
| `NG_012337.1(NM_003002.2):c.[274=;275del]` | normalized | accepted_divergence | — | — |
| `NG_012337.1(NM_003002.2):c.[274del;275=]` | genomic | accepted_divergence | — | — |
| `NG_012337.1(NM_003002.2):c.[274del;275=]` | normalized | accepted_divergence | — | — |
| `NG_012337.1(NM_003002.2):c.pterdel` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_003002.2):c.qterdel` | normalized | spec_citation | — | — |
| `NG_012337.1(NM_012459.2):c.-1_*1del` | protein_description | accepted_divergence | — | — |
| `NG_012337.1(NM_012459.2):c.15_16insTGGAAGGTGGCGAGCCT` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_012337.1(NM_012459.2):c.15_16insTGGAAGGTGGCGAGCCT` | genomic | accepted_divergence | — | — |
| `NG_012337.1(NM_012459.2):c.15_17dup` | coding_protein_descriptions | accepted_divergence | — | — |
| `NG_012337.1:274` | normalized | accepted_divergence | — | — |
| `NG_012337.1:g.4delins7_31` | normalized | accepted_divergence | — | — |
| `NG_012337.1:g.4delins7_50` | normalized | accepted_divergence | — | — |
| `NG_012337.1:g.[4_5insT;4insA;4_5insA;5del]` | errors | accepted_divergence | — | — |
| `NG_012337.1:g.[4_5insT;4insA;4_5insA]` | errors | accepted_divergence | — | — |
| `NG_012337.1:g.[4del;4_5insT;4insA;4_5insA;5del]` | errors | accepted_divergence | — | — |
| `NG_012337.1:g.[4del;4_5insT;4insA;4_5insA]` | errors | accepted_divergence | — | — |
| `NG_012337.3(NM_003002):274+400G>T` | errors | accepted_divergence | — | — |
| `NG_012337.3(NM_003002.4):c.274delinsNG_009930.1(NM_001099625.2):c.1010` | errors | accepted_divergence | — | — |
| `NG_012337.3(NM_003002.4):r.274+10G>T` | errors | accepted_divergence | — | — |
| `NG_017013.2:g.17496_17497insAGCTGCTCAGATAGCGA` | normalized | accepted_divergence | — | — |
| `NG_029724.1(NM_004321.7):c.101del` | normalized | accepted_divergence | — | — |
| `NM_000143.3:c.-1_1insCAT` | protein_description | accepted_divergence | — | — |
| `NM_002001.2:c.1_3delinsATG` | normalized | accepted_divergence | — | — |
| `NM_003002.2:c.273del` | infos | accepted_divergence | — | — |
| `NM_003002.2:c.[100del;200_201insNM_003002.2:274+20]` | errors | accepted_divergence | — | — |
| `NM_003002.2:c.[100del;200_201insNM_003002.2:274+20]` | infos | accepted_divergence | — | — |
| `NM_003002.4:c.206_210delins190_220inv` | normalized | accepted_divergence | — | — |
| `NM_003002.4:n.206_210del` | normalized | accepted_divergence | — | — |

## Disposition tallies

| axis | accepted_divergence | known_bug | improvement | reference_unavailable | spec_citation |
|---|---:|---:|---:|---:|---:|
| coding_protein_descriptions | 27 | 0 | 0 | 0 | 0 |
| errors | 17 | 0 | 0 | 0 | 0 |
| genomic | 22 | 4 | 0 | 0 | 1 |
| infos | 8 | 0 | 0 | 0 | 1 |
| noncoding | 0 | 8 | 0 | 0 | 0 |
| normalized | 25 | 9 | 4 | 5 | 11 |
| protein_description | 9 | 0 | 0 | 0 | 60 |
| rna_description | 0 | 0 | 0 | 0 | 1 |
