# HGVS spec-compliant variant examples (extracted from hgvs-nomenclature.org/stable)

Source: https://hgvs-nomenclature.org/stable/ — every page under
`/recommendations/`, `/background/`, and `/versions/21.0/` was crawled.

For each variant string the spec shows on a page, this file records:

1. **Compliant** — examples the spec presents as the *correct* form. Normalization
   should leave these unchanged.
2. **Pairs (incorrect → correct)** — examples the spec shows as wrong, paired
   with the recommended replacement. Grouped by *kind of transformation*, since
   only some kinds are the responsibility of the normalizer.
3. **Standalone non-compliant** — the spec calls these out as wrong, but does
   not give a paired correct form. Useful as "must fail to parse" or "must be
   flagged" cases, but not as normalization-roundtrip cases.

Notes on completeness and caveats:

- The spec uses many partial strings (e.g. `c.76A>G` without a reference
  prefix). These are kept as-is — for unit tests against the parser/normalizer
  some will need a reference prefix added.
- The spec also shows tabular *syntax* placeholders that look like variants but
  aren't (e.g. `g.123_456` as a numbering example). Those are not included
  unless they're flagged in context as a real example.
- Pages crawled: `recommendations/{general,summary,uncertain,style,grammar,
  checklist}`, `recommendations/DNA/{substitution,deletion,duplication,
  insertion,inversion,delins,repeated,alleles,complex}`, `recommendations/RNA/
  {substitution,deletion,duplication,insertion,inversion,delins,repeated,
  alleles}`, `recommendations/protein/{substitution,deletion,duplication,
  insertion,delins,repeated,frameshift,extension,alleles}`,
  `background/{numbering,glossary,simple,refseq,standards,basics,history}`,
  `versions/21.0/`.

---

## 1. Compliant examples

### 1.1 DNA substitution

- `NC_000023.10:g.33038255C>A`
- `NG_012232.1(NM_004006.2):c.93+1G>T`
- `LRG_199t1:c.54G>H` — IUPAC ambiguity code
- `NM_004006.2:c.123=` — no change at tested position
- `LRG_199t1:c.94-23_188+33=` — no variants in range
- `LRG_199t1:c.85=/T>C` — mosaic
- `NM_004006.2:c.85=//T>C` — chimeric
- `r.76a>g`
- `c.(76A>G)` — DNA-level inferred from RNA data
- `NC_000023.10:g.33357783G>A` — promoter region
- `NC_000023.10(NM_004006.1):c.-128354C>T` — promoter via genomic ref
- `NC_000023.10(NM_000109.3):c.-401C>T` — promoter via genomic ref
- `L01538.1:g.1407C>T`
- `NC_000023.11:g.32849790T>A`
- `NC_000023.11:g.32389644G>A`
- `NC_000023.9:g.32317682G>A`
- `NC_000023.10:g.32407761G>A`
- `NG_012232.1:g.954966C>T`
- `LRG_199:g.954966C>T`
- `NM_004006.2:c.4375C>T`
- `NR_002196.1:n.601G>T`
- `NP_003997.1:p.Arg1459Ter`
- `NP_003997.1:p.Arg1459*`
- `NC_000023.9:g.32290917C>T`
- `NM_004006.3:c.5234G>A`
- `NC_000023.10(NM_004006.3):c.357+1G>A`
- `LRG_199t1:c.357+1G>A`
- `LRG_199t1:c.11T>G`
- `NC_012920.1:m.3243A>G`
- `NC_012920.1(MT-TL1):m.3243A>G`
- `NC_012920.1(MT-TL1):n.14A>G`
- `NC_012920.1:m.3460G>A`
- `NC_012920.1(MT-ND1):m.3460G>A`
- `c.78T>C`
- `c.-78G>A`
- `c.*78T>A`
- `c.78+45T>G`
- `c.79-45G>T`
- `c.-106+2T>A`
- `c.*639-1G>A`
- `c.88+2T>G`
- `c.89-1G>T`

### 1.2 DNA deletion

- `NC_000001.11:g.1234del`
- `NC_000001.11:g.1234_2345del`
- `NC_000023.11:g.33344591del`
- `NM_004006.2:c.5697del` — applies 3' rule on A-stretch
- `NC_000023.11:g.32343183del` — applies 3' rule on T-stretch
- `NC_000023.11:g.33344590_33344592del`
- `NC_000023.11(NM_004006.2):c.183_186+48del` — crosses exon/intron border
- `LRG_199t1:c.3921del` — exon/exon border (3' rule exception)
- `LRG_199t1:c.1704+1del` — exon/intron border
- `LRG_199t1:c.1813del` — intron/exon border
- `NC_000023.11(NM_004006.2):c.4072-1234_5155-246del` — multi-exon deletion
- `NC_000023.11(NM_004006.2):c.(4071+1_4072-1)_(5154+1_5155-1)del` — uncertain breakpoints
- `NC_000023.11(NM_004006.2):c.(3996_4196)_(5090_5284)del` — MLPA probe-based
- `NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del` — SNP-array
- `NC_000023.11:g.(?_31120496)_(33339477_?)del` — uncertain boundaries
- `NC_000023.11:g.33344590_33344592=/del` — mosaic
- `NC_000023.11:g.33344590_33344592=//del` — chimeric
- `NG_012232.1:g.123_128del`
- `NM_000492.3:c.1521_1523del`
- `NM_000492.3:c.1522_1524del`
- `NC_000023.10:g.32361300del`
- `NM_004006.1:c.5697del`
- `c.57_58del`
- `NC_000006.11:g.117198495_117198496del`
- `NC_000023.10:g.32459297del`
- `c.4375_4379del`

### 1.3 DNA duplication

- `NC_000001.11:g.1234dup`
- `NC_000001.11:g.1234_2345dup`
- `NM_004006.2:c.20dup`
- `NC_000023.10:g.33229410dup`
- `NM_004006.2:c.5697dup`
- `NC_000023.11:g.32343183dup`
- `NM_004006.2:c.20_23dup`
- `NC_000023.11:g.33211290_33211293dup`
- `NC_000023.11(NM_004006.2):c.260_264+48dup`
- `NC_000023.11:g.32844735_32844787dup`
- `NC_000023.11(NM_004006.2):c.3921dup` — exon/exon junction
- `NC_000023.11:g.32441180dup`
- `NC_000023.11(NM_004006.2):c.1704+1dup`
- `NC_000023.11(NM_004006.2):c.1813dup`
- `NC_000023.11(NM_004006.2):c.4072-1234_5155-246dup`
- `NC_000023.11(NM_004006.2):c.720_991dup`
- `NC_000023.11(NM_004006.2):c.(4071+1_4072-1)_(5154+1_5155-1)dup`
- `NC_000001.11(NM_206933.2):c.[675-542_1211-703dup;1211-703_1211-702insGTAAA]`
- `NC_000023.11:g.(32381076_32382698)_(32430031_32456357)[3]` — triplication
- `NC_000023.11(NM_004006.2):c.(4071+1_4072-1)_(5154+1_5155-1)[3]`
- `NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup`
- `NC_000023.11:g.(?_31120496)_(33339477_?)dup`
- `NC_000023.11:g.pter_qtersup`
- `NC_000023.11:g.pter_qter[2]`
- `NC_000023.11:g.33344590_33344592=/dup`
- `NC_000023.11:g.33344590_33344592=//dup`
- `c.4375_4385dup`

### 1.4 DNA insertion

- `NC_000001.11:g.1234_1235insACGT`
- `NC_000023.10:g.32867861_32867862insT`
- `NM_004006.2:c.169_170insA`
- `NC_000023.10:g.32862923_32862924insCCT`
- `LRG_199t1:c.240_241insAGG`
- `NM_004006.2:c.849_850ins858_895` — inserted seq from same ref
- `NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]`
- `NM_004006.2:c.419_420ins[T;401_419]`
- `LRG_199t1:c.419_420ins[T;450_470;AGGG]`
- `NC_000006.11:g.10791926_10791927ins[NC_000004.11:g.106370094_106370420;A[26]]`
- `NM_004006.2:c.849_850ins850_900inv` — inverted-duplication insertion
- `NM_004006.2:c.900_901ins850_900inv`
- `LRG_199t1:c.940_941ins[885_940inv;A;851_883inv]`
- `NM_004006.2:c.940_941ins[903_940inv;851_885inv]`
- `NM_004006.2:c.(222_226)insG`
- `NC_000004.11:g.(3076562_3076732)insN[12]`
- `NC_000023.10:g.32717298_32717299insN`
- `NM_004006.2:c.761_762insN`
- `NM_004006.2:c.761_762insNNNNN`
- `NM_004006.1:c.761_762insN[5]`
- `NC_000023.10:g.32717298_32717299insN[100]`
- `NC_000023.10:g.32717298_32717299insN[(80_120)]`
- `NC_000023.10:g.32717298_32717299insN[?]`
- `NC_000006.11:g.8897754_8897755ins[N[543];8897743_8897754]`
- `g.?_?ins[NC_000023.10:g.(12345_23456)_(34567_45678)]`
- `NC_000002.11:g.48031621_48031622ins[TAT;48026961_48027223;GGC]`
- `NC_000002.11:g.?_?ins[NC_000022.10:g.35788169_35788352]`
- `c.4375_4376insACCT`

### 1.5 DNA inversion

- `NC_000023.10:g.32361330_32361333inv`
- `NM_004006.2:c.5657_5660inv`
- `NM_004006.2:c.4145_4160inv`
- `NC_000023.10:g.111754331_111966764inv`

### 1.6 DNA deletion-insertion (delins)

- `NC_000001.11:g.123delinsAC`
- `NC_000001.11:g.123_129delinsAC`
- `NC_000023.11:g.32386323delinsGA`
- `NM_004006.2:c.6775_6777delinsC`
- `LRG_199t1:c.79_80delinsTT`
- `LRG_199t1:c.145_147delinsTGG`
- `LRG_199t1:c.9002_9009delinsTTT`
- `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG`
- `NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]`
- `NC_000022.10:g.42522624_42522669delins42536337_42536382`
- `NC_000012.11:g.6128892_6128954delins[NC_000022.10:g.17179029_17179091]`
- `NM_000797.3:c.812_829delins908_925`
- `NM_004006.2:c.812_829delinsN[12]`
- `NM_007294.3:c.2077delinsATA`
- `NM_004006.2:c.145_147delinsTGG`
- `NM_004006.2:c.2077delinsATA`
- `c.4375_4376delinsAGTT`

### 1.7 DNA repeats

- `NC_000014.8:g.123CAG[23]`
- `NC_000014.8:g.123_191CAG[19]CAA[4]`
- `NC_000014.8:g.101179660_101179695TG[14]`
- `NC_000014.8:g.[101179660_101179695TG[14]];[101179660_101179695TG[18]]`
- `NM_023035.2:c.6955_6993CAG[26]`
- `NC_000003.12:g.63912687_63912716AGC[13]`
- `NM_000333.3:c.89_118AGC[13]`
- `NC_000003.12:g.(63912602_63912844)insN[9]`
- `NM_000333.3:c.(4_246)insN[9]`
- `NC_000003.12:g.(63912602_63912844)delN[15]`
- `NM_000333.3:c.(4_246)delN[15]`
- `NM_002024.5:c.-128_-69GGC[10]GGA[1]GGC[9]GGA[1]GGC[10]`
- `NM_002024.5:c.-128_-69GGC[68]GGA[1]GGC[10]`
- `NM_002024.5:c.-128_-69GGM[108]`
- `NM_002024.5:c.(-144_-16)insN[(1800_2400)]`
- `LRG_763t1:c.52_153CAG[21]CAA[1]CAG[1]CCG[1]CCA[1]CCG[7]CCT[2]`
- `LRG_763t1:c.54_110GCA[23]`
- `NM_000492.3:c.1210-33_1210-6GT[11]T[6]`
- `NM_000492.3:c.1210-12_1210-6T[7]`
- `NC_000012.11:g.112036755_112036823CTG[9]TTG[1]CTG[13]`
- `NC_000001.11:g.57367047_57367121ATAAA[15]`
- `NM_021080.3:c.-136-75952_-136-75878ATTTT[15]`
- `c.1210-33_1210-6GT[(9_13)]T[(4_8)]`
- `c.1210-12_1210-6T[(5_9)]`
- `NC_000003.12:g.63912687AGC[(50_60)]`
- `NC_000003.12:g.63912687AGC[(60_?)]`

### 1.8 DNA alleles

- `NC_000001.11:g.[123G>A;345del]` — cis
- `NC_000001.11:g.[123G>A];[345del]` — trans
- `NC_000001.11:g.123G>A(;)345del` — uncertain phase
- `NM_004006.2:c.[2376G>C;3103del]`
- `NC_000023.10:g.[30683643A>G;33038273T>G]`
- `NM_004006.2:c.[2376G>C];[3103del]`
- `NM_004006.2:c.[2376G>C];[2376G>C]` — homozygous
- `NM_004006.2:c.[296T>G;476T>C;1083A>C];[296T>G;1083A>C]`
- `NM_004006.2:c.[2376G>C];[2376=]`
- `NM_004006.2:c.[2376del];[2376=]`
- `NM_004006.2:c.[2376_2399dup];[2376_2399=]`
- `NM_004006.2:c.[2376_2377insGT];[2376_2377=]`
- `NM_004006.2:c.[2376G>C];[?]`
- `NM_004006.2:c.2376G>C(;)3103del`
- `NM_004006.2:c.2376G>C(;)(2376G>C)`
- `NM_004006.2:c.[296T>G;476T>C];[476T>C](;)1083A>C`
- `NM_004006.2:c.[296T>G];[476T>C](;)1083G>C(;)1406del`
- `LRG_199t1:c.[76A>C];[76=]`
- `LRG_199t1:c.[76A>C];[0]` — male X-chromosome
- `[NM_004004.2:c.35del];[NM_006783.1:c.689_690insT]`
- `[GJB2:c.35del];[GJB6:c.689_690insT]`
- `[M59228.1:g.250G>C;AF209160.1:g.572CA[(11_21)];Z11861.1:g.61T>C;Z16803.1:g.114A[(18_22)]]`
- `[C;13;T;21]` — short haplotype
- `NM_004006.1:c.[837G>A;1704+51T>C;3734C>T;6438+2669T[(16_23)];6571C>T;7098+13212GT[(15_19)]]`
- `[G;C;C;18;T;17]`
- `OPRM1:c.[118A>G];[118A=]`
- `OPRM1:c.[118A=];[118A=]`
- `OPRM1:c.[118>G];[118>G]`
- `g.[123456A>G;345678G>C]`
- `g.[123456A>G];[345678G>C]`

### 1.9 DNA complex / large rearrangements

- `NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825272]`
- `NC_000011.10:g.pter_5825272delins[NC_000002.12:g.pter_8247756]`
- `NC_000002.12:g.17450009_qterdelins[NC_000011.10:g.108111987_qter]`
- `NC_000011.10:g.108111982_qterdelins[NC_000002.12:g.17450009_qter]`
- `NC_000003.12:g.pter_36969141delins[CATTTGTTCAAATTTAGTTCAAATGA;NC_000014.9:g.29745314_qterinv]`
- `NC_000014.9:g.29745314_qterdelins[NC_000003.12:g.36969141_pterinv]`
- `NC_000009.12:g.pter_26393001delins102425452_qterinv`
- `NC_000009.12:g.102425452_qterdelinspter_26393001inv`
- `NC_000003.12:g.158573187_qterdelins[NC_000008.11:g.(128534000_128546000)_qter]`
- `NC_000005.10:g.29658442delins[NC_000010.11:g.67539995_qterinv]`
- `NC_000006.12:g.[776788_93191545inv;93191546T>C]`
- `NC_000002.12:g.[32310435_32310710del;32310711_171827243inv;insG]`
- `NC_000023.11:g.(86200001_103700000)del`
- `NC_000022.11:g.pter_(12200001_14700000)del::(37600001_410000000)_qterdel`
- `NC_000008.11:g.(127300001_131500000)_(131500001_136400000)dup`
- `NC_000008.11:g.(131500001_136400000)ins(127300001_131500000)_(131500001_136400000)inv`
- `NC_000004.12:g.134850793_134850794ins[NC_000023.11:g.89555676_100352080inv]`
- `NC_000004.12:g.134850793_134850794ins[NC_000023.11:g.89555676_100352080]`
- `NC_000023.11:g.89555676_100352080del`
- `NC_000022.11:g.[pter_(12200001_14700000)del::(37600001_410000000)_qterdel]sup`

### 1.10 Methylation / epigenetic

- `g.12345678_12345901|gom`
- `g.12345678_12345901|lom`
- `g.12345678_12345901|met=`

### 1.11 RNA substitution

- `NM_004006.3:r.123c>g`
- `NM_004006.3:r.76a>c`
- `NM_004006.3:r.(1388g>a)` — predicted
- `NM_004006.3:r.123=`
- `NM_004006.1:r.-14a>c`
- `NM_004006.3:r.*41u>a`
- `NM_004006.3:r.[897u>g,832_960del]`
- `NM_004006.1:r.0` — no RNA detected
- `NM_004006.3:r.spl` — splicing affected
- `NM_004006.3:r.?` — uncertain
- `NM_004006.3:r.85=/u>c` — mosaic
- `NM_004006.3:r.85=//u>c` — chimeric
- `LRG_199t1:r.11u>g`
- `r.123a>u`
- `r.76a>u`
- `r.123c>u`
- `r.124a>u`
- `r.(76a>c)`
- `r.(=)` — no RNA changes expected
- `r.0?` — possibly no transcript
- `r.spl?`

### 1.12 RNA deletion

- `NM_004006.3:r.123_127del`
- `LRG_199t1:r.10del`
- `NM_004006.2:r.6_8del`
- `LRG_2t1:r.1034_1036del`
- `LRG_199t1:r.(4072_5145del)` — predicted, uncertain
- `LRG_199t1:r.=/6_8del` — mosaic
- `r.123_128del`

### 1.13 RNA duplication

- `NM_004006.3:r.123_345dup`
- `r.7dup`
- `r.6_8dup`
- `r.17_18ins5_16` — non-tandem duplicative event correctly written as insertion

### 1.14 RNA insertion

- `NM_004006.3:r.123_124insauc`
- `LRG_199t1:r.426_427insa`
- `LRG_199t1:r.756_757insuacu`
- `NM_004006.2:r.(222_226)insg`
- `NM_004006.2:r.549_550insn`
- `NM_004006.2:r.761_762insnnnnn`
- `LRG_199t1:r.1149_1150insn[100]`
- `NG_012232.1(NM_004006.2):r.2949_2950ins[2950-30_2950-12;2950-4_2950-1]`
- `r.2949_2950ins[2950-30_2950-12;uuag]`
- `LRG_199t1:r.186_187ins186+1_186+4`
- `NC_000023.10(NM_004006.2):r.186_187ins186+1_186+4`
- `NG_012232.1(NM_004006.2):r.186_187ins186+1_186+4`
- `NC_000023.10(NM_004006.2):r.357_358ins357+1_357+12`
- `NG_012232.1(NM_004006.2):r.357_358ins357+1_357+12`

### 1.15 RNA inversion

- `NM_004006.3:r.123_345inv`
- `r.177_180inv`
- `r.203_506inv`

### 1.16 RNA delins

- `NM_004006.3:r.123_127delinsag`
- `r.775delinsga`
- `r.775_777delinsc`
- `r.902_909delinsuuu`
- `r.142_144delinsugg`
- `NM_004006.2:r.2623_2803delins2804_2949`
- `r.415_1655delins[AC096506.5:g.409_1649]`
- `r.1401_1446delins[NR_002570.3:r.1513_1558]`
- `NM_007294.3:r.2077delinsaua`

### 1.17 RNA repeats

- `NM_004006.3:r.9495_9497[4]`
- `NM_004006.3:r.-110_-108[6]`
- `NM_004006.3:r.9495caa[4]`
- `NM_004006.3:r.-110gcu[6]`
- `r.-124_-123[14]`
- `r.-124_-123[14];[18]`
- `r.-128_-126[79]`
- `r.-128_-126[(600_800)]`
- `r.53agc[19]`
- `r.53_55[31]`

### 1.18 RNA alleles

- `NM_004006.3:r.[123c>a;345del]`
- `NM_004006.3:r.[123c>a];[345del]`
- `NM_004006.3:r.123c>a(;)345del`
- `LRG_199t1:r.[76a>u;103del]`
- `LRG_199t1:r.[(578c>u;1339a>g;1680del)]`
- `LRG_199t1:r.[76a>u];[103del]`
- `NM_004006.2:r.[76a>u];[76a>u]`
- `LRG_199t1:r.[76a>u];[76=]`
- `NM_004006.2:r.[76a>u];[?]`
- `LRG_199t1:r.[897u>g,832_960del]`
- `LRG_199t1:r.[76a>c];[76=]`
- `LRG_199t1:r.[76a>c];[0]`
- `r.[123a>u,122_154del]`

### 1.19 Protein substitution

- `NP_003997.1:p.Trp24Cys` — experimental
- `NP_003997.1:p.Trp24Ter`
- `NP_003997.1:p.W24*`
- `NP_003997.1:p.(Trp24Cys)` — predicted
- `LRG_199p1:p.Trp24Cys`
- `NP_003997.1:p.Cys188=`
- `LRG_199p1:p.0`
- `LRG_199p1:p.(Met1?)`
- `NP_003997.1:p.Leu2_Met124del`
- `p.Met1_Leu2insArgSerThrVal`
- `p.Met1ext-5`
- `NP_003997.1:p.?`
- `NP_003997.1:p.(Gly56Ala^Ser^Cys)` — uncertain identity
- `LRG_199p1:p.Trp24=/Cys` — mosaic
- `LRG_199t1:p.(Trp24=/Cys)` — predicted mosaic
- `NP_003997.1:p.Gln2366Lys`
- `NP_003997.1:p.Trp24_Val25delinsCysArg`
- `p.Arg782Xaa`
- `LRG_199p1:p.(Val25Gly)`
- `YP_003024026.1:p.(Ala52Thr)` — mitochondrial protein
- `p.(Ser123Arg)`
- `p.(Cys123Gly)`
- `p.(Ser42Cys)`
- `p.(Arg234=)`
- `p.Ser321Arg`
- `Lys76Asn`
- `p.Trp41*`
- `p.(Gly26Cys)`
- `p.(Leu54=)`
- `p.(Ile43Val)`
- `p.(Arg123Ser)`
- `p.(Arg1459Ter)`
- `p.Tyr4Ter`
- `p.Trp26Ter`
- `p.Trp26*`

### 1.20 Protein deletion / duplication / insertion

- `NP_003997.2:p.Val7del`
- `NP_003997.2:p.Lys23_Val25del`
- `NP_003997.2:p.(Val7del)`
- `NP_003997.2:p.Trp4del`
- `NP_003997.2:p.(Pro458_Gly460del)`
- `NP_003997.2:p.(His321_Glu383del)`
- `NP_003997.2:p.(Asp90_Val120del)` — 3' rule applied
- `NP_003997.2:p.(His321Leufs*3)`
- `NP_003997.1:p.Val7=/del` — mosaic
- `NP_003997.1:p.(Val7=/del)`
- `NP_000483.3:p.Phe508del`
- `p.Gly2_Met46del`
- `p.Ser5del`
- `NP_003997.2:p.Val7dup`
- `NP_003997.2:p.Lys23_Val25dup`
- `NP_003997.2:p.(Val7dup)`
- `NP_003997.2:p.Trp4dup`
- `NP_003997.2:p.(Pro458_Gly460dup)`
- `NP_003997.2:p.(His321_Glu383dup)`
- `NP_003997.2:p.(Asp90_Val120dup)`
- `NP_003997.2:p.(Asn444Lysfs*15)`
- `NP_003997.1:p.Val7=/dup`
- `NP_003997.1:p.(Val7=/dup)`
- `p.His4_Gln5insAla`
- `p.Lys2_Gly3insGlnSerLys`
- `p.(Met3_His4insGlyTer)`
- `NP_004371.2:p.(Pro46_Asn47insSerSerTer)`
- `p.Arg78_Gly79insXaa[23]`
- `NP_060250.2:p.Gln746_Lys747ins*63`
- `NP_003997.1:p.(Ser332_Ser333insXaa)`
- `NP_003997.1:p.(Val582_Asn583insXaa[5])`
- `NP_003997.1:p.(Val582_Asn583insXaaXaaXaaXaaXaa)`
- `p.His7_Gln8insGly4_Ser6` — non-tandem duplicative event written as insertion

### 1.21 Protein delins

- `NP_004371.2:p.Asn47delinsSerSerTer`
- `NP_004371.2:p.(Asn47delinsSerSerTer)`
- `NP_003070.3:p.Glu125_Ala132delinsGlyLeuHisArgPheIleValLeu`
- `NP_003070.3:p.(Glu125_Ala132delinsGlyLeuHisArgPheIleValLeu)`
- `p.Cys28delinsTrpVal`
- `p.Cys28_Lys29delinsTrp`
- `p.(Pro578_Lys579delinsLeuTer)`
- `NP_000213.1:p.(Val559_Glu561del)`
- `p.[Ser44Arg;Trp46Arg]` — 2 separate variants with intervening AA

### 1.22 Protein repeats

- `NP_0123456.1:p.Ala2[10]`
- `NP_0123456.1:p.Ala2[10];[11]`
- `NP_0123456.1:p.Arg65_Ser67[12]`
- `p.Ala2[10]`
- `p.Ala2[12]`
- `p.Gln18[23]`
- `p.(Gln18)[(70_80)]`

### 1.23 Protein frameshift

- `NP_0123456.1:p.Arg97ProfsTer23`
- `NP_0123456.1:p.Arg97fs`
- `p.(Arg123LysfsTer34)`
- `p.Arg97ProfsTer23`
- `p.Arg97Profs*23`
- `p.Arg97fs`
- `p.(Tyr4*)`
- `p.Glu5ValfsTer5`
- `p.Glu5fs`
- `p.Ile327Argfs*?`
- `p.Ile327fs`
- `p.Gln151Thrfs*9`
- `p.Arg456GlyfsTer17`
- `p.Arg456Glyfs*17`

### 1.24 Protein extension

- `NP_003997.2:p.Met1ext-5`
- `NP_003997.2:p.Ter110GlnextTer17`
- `p.(Ter110GlnextTer17)`
- `p.(*110Glnext*17)`
- `p.Ter327ArgextTer?`
- `p.*327Argext*?`
- `p.(Met1ext-8)`

### 1.25 Protein alleles

- `NP_003997.1:p.[Ser73Arg;Asn103del]` — cis
- `NP_003997.1:p.[Ser73Arg];[Asn103del]` — trans
- `NP_003997.1:p.(Ser73Arg)(;)(Asn103del)` — uncertain phase
- `NP_003997.1:p.[Ser68Arg;Asn594del]`
- `NP_003997.1:p.[(Ser68Arg;Asn594del)]` — predicted cis
- `NP_003997.1:p.[Ser68Arg];[Ser68Arg]` — homozygous
- `NP_003997.1:p.[(Ser68Arg)];[(Ser68Arg)]`
- `NP_003997.1:p.(Ser68Arg)(;)(Ser68Arg)`
- `NP_003997.1:p.[Ser68Arg];[Asn594del]`
- `NP_003997.1:p.[(Ser68Arg)];[?]`
- `NP_003997.1:p.[Ser68Arg];[Ser68=]`
- `NP_003997.1:p.(Ser68Arg)(;)(Asn594del)`
- `NP_003997.2:p.[(Asn158Asp)(;)(Asn158Ile)]^[(Asn158Val)]`
- `NP_003997.1:p.[Lys31Asn,Val25_Lys31del]` — two proteins from one allele
- `p.[Ser86Arg];[Ser86=]`
- `p.[Ser86Arg];[0]`

### 1.26 Uncertain / range

- `NC_000023.10:g.(33038277_33038278)C>T`
- `NC_000023.11:g.(31729716_31774235)_(32216847_32287541)del`
- `NC_000023.11:g.(31729663_31774080)_(32216973_32287624)del`
- `NC_000023.10:g.(32218983_32238146)_(32984039_33252615)del`
- `NC_000023.10:g.(?_32238146)_(32984039_?)del`
- `NC_000013.11:g.(19385993_19394916)_(25045592_25059364)del`
- `NC_000013.11:g.(?_19394916)_(25045592_?)del`
- `NC_000023.10:g.(32057077_32364657)_(32975163_33394206)del`
- `NC_000023.10:g.(?_32364657)_(32975163_?)del`
- `NC_000023.10:g.(32057077_32364657)_(32894352_33055973)del`
- `NC_000023.11:g.(31775822_31819974)_(32217064_32278336)del`
- `NC_000023.11:g.(31729716_31773911)_(32217064_32287541)dup`
- `NC_000023.11:g.(31729662_31774079)_(32216972_32287623)dup`
- `NC_000023.11:g.(31775822_31817965)_(32218461_32278336)dup`
- `g.(123456_234567)_(345678_456789)del`
- `g.(?_234567)_(345678_?)del`
- `c.(370A>C^372C>R)` — back-translation alternatives
- `NM_000517.4:c.424C>T^NM_000558.3:c.424C>T` — multi-gene mappable

### 1.27 Special protein cases

- `p.(Ala123_Pro131)Ter`
- `p.(Ala123_Pro131)fs`
- `p.(Ala123_Pro131)insXaa[4]`
- `p.Gly719(Ala^Ser)fsTer23`
- `p.(Gly23GlufsTer7^Gly23CysfsTer26)`

---

## 2. Pairs (incorrect → correct)

Grouped by transformation type. Some kinds (3' rule, ins→dup) are
**normalization** in the strict sense; others are **canonicalization** of
parser-level forms.

### 2.1 3' rule shifting (true normalization)

- `c.79_80GC>TT` → `NG_012232.1:g.12_13delinsTG` *(also a substitution→delins reclassification — see 2.4)*
- `c.3922del` → `LRG_199t1:c.3921del` — exon/exon junction 3' rule exception
- `c.1_24dup` → `c.9_32dup` — 3' rule pushes duplication to most 3' position
- `c.3922dup` → `NC_000023.11(NM_004006.2):c.3921dup` — exon/exon junction
- `p.Val89_Gln119del` → `p.(Asp90_Val120del)` — 3' rule (C-terminal)
- `p.Val89_Gln119dup` → `p.(Asp90_Val120dup)` — 3' rule (C-terminal)
- `p.Ser4del` → `p.Ser5del` — derive from protein sequence comparison; resulting position obeys 3' rule
- `p.Ser4dup` → (correct form not given on page) — same principle
- `r.1033_1035del` → `LRG_2t1:r.1034_1036del` — 3' rule application
- `NM_002024.5:c.-129CGG[79]` → `NM_002024.5:c.-128_-69GGC[68]GGA[1]GGC[10]` — single repeat unit cannot describe mixed reference

### 2.2 Drop redundant deleted/duplicated sequence

- `NC_000023.11:g.33344591delA` (or `delG`) → `NC_000023.11:g.33344591del`
- `NC_000023.11:g.33344590_33344592delGAT` (or `delTTA`) → `NC_000023.11:g.33344590_33344592del`
- `NM_004006.2:c.20dupT` → `NM_004006.2:c.20dup`
- `c.20_23dupTAGA` (or `dupTGGA`) → `NM_004006.2:c.20_23dup`
- `r.6_8deluug` → `NM_004006.2:r.6_8del` (per recommendation — sequence form discouraged)
- `r.6_8dupugc` → `r.6_8dup`
- `c.4375_4379delCGATT` → `c.4375_4379del`
- `c.4375_4385dupCGATTATTCCA` → `c.4375_4385dup`
- `NC_000023.11:g.32386323delTinsGA` → `NC_000023.11:g.32386323delinsGA`
- `NC_000023.11:g.32386323delCinsGA` → `NC_000023.11:g.32386323delinsGA` (also fixes wrong base)
- `NM_004006.2:c.6775_6777delGAGinsC` → `NM_004006.2:c.6775_6777delinsC`
- `c.4375_4376delCGinsAGTT` → `c.4375_4376delinsAGTT`

### 2.3 Range vs. count notation (`del6` → range)

- `NG_012232.1:g.123del6` → `NG_012232.1:g.123_128del`
- `g.123dup6` → `g.123_128dup` (page notes `g.124_129dup` is also wrong — must specify actual duplicated positions)
- `r.123del6` → `r.123_128del`
- `p.Arg45del6` → range form (not given explicitly on page)
- `c.5439_5430ins6` → use `insN[6]` style

### 2.4 Substitution → delins (consecutive nucleotides)

- `c.[79G>T;80C>T]` → `LRG_199t1:c.79_80delinsTT`
- `c.79GC>TT` → `NG_012232.1:g.12_13delinsTG`
- `NG_012232.1:g.12GC>TG` → `NG_012232.1:g.12_13delinsTG`
- `g.4GC>TG` → delins form (per spec `delinsTG` style)
- `c.[145C>T;147C>G]` → `NM_004006.2:c.145_147delinsTGG` — variants in same codon
- `NG_012232.1:g.12G>T(;)13C>G` → `NG_012232.1:g.12_13delinsTG` — phase-unknown consecutive
- `NM_007294.3:c.[2077G>A;2077_2078insTA]` → `NM_007294.3:c.2077delinsATA`
- `r.76_77aa>ug` (or `r.76aa>ug`) → `NM_004006.3:r.76_77delinsug`
- `r.4gc>ug` → delins form
- `NM_007294.3:r.[2077g>a;2077_2078insua]` → `NM_007294.3:r.2077delinsaua`
- `r.[142c>u;144a>g]` → `r.142_144delinsugg`
- `p.[Arg76Ser;Cys77Trp]` → `p.Arg76_Cys77delinsSerTrp` — adjacent AA changes
- `p.TrpVal24CysArg` → `NP_003997.1:p.Trp24_Val25delinsCysArg` — multi-AA must be delins
- `p.Ser44_Trp46delinsArgLeuArg` → `p.[Ser44Arg;Trp46Arg]` — non-adjacent variants written individually

### 2.5 Insertion → duplication (tandem case)

- `c.19_20insT` → `NM_004006.2:c.20dup`
- `g.456_457ins123_456` → use `dup` form
- `r.6_7insu` → use `dup` form
- `r.456_457ins123_456` → use `dup` form

### 2.6 Inverted-duplication notation

- `g.123_456dupinv` → `g.234_235ins123_234inv` (use `ins[..._...inv]`)
- `r.123_456dupinv` → `ins[..._...inv]` form
- `r.124_500delinsoAB053210.2:r.1289-365_1289-73` → `r.124_500delins[AB053210.2:r.1289-365_1289-73inv]` — `o` strand-flag deprecated, `inv` instead

### 2.7 IVS notation removed

- `c.IVS2+2T>G` → `c.88+2T>G`
- `c.IVS2-1G>T` → `c.89-1G>T`

### 2.8 Allele/phase notation

- `[c.76A>C+c.83G>C]` → `c.[76A>C;83G>C]` — `+` separator deprecated
- `c.[76A>C;c.83G>C]` → `c.[76A>C;83G>C]` — coord-type prefix only once
- `[r.76a>c+r.83g>c]` → `r.[76a>c;83g>c]`
- `[p.Ser73Arg+p.Asn103del]` → `p.[Ser73Arg;Asn103del]`
- `c.[76A>C]` (single allele) → `c.[76A>C];[76=]` — must specify both alleles
- `c.[76A>C];[]` → `c.[76A>C];[76=]`
- `r.[76a>c]` → `r.[76a>c];[=]`
- `p.[Ser73Arg];[]` → `p.[Ser73Arg];[Ser73=]`
- `p.([Phe233Leu;Cys690Trp])` → `p.[(Phe233Leu;Cys690Trp)]` — parentheses inside brackets

### 2.9 Range separator (hyphen → underscore)

- `c.12-14del` → `c.12_14del`
- `c.123-65_-50` → `c.123-65_123-50` — must specify both endpoints

### 2.10 Insertion missing flanking position

- `c.52insT` → `c.51_52insT`
- `g.123ˆ124insG` (or `g.123ˆ124G`) → underscore form (`g.123_124insG`)
- `r.123ˆ124insu` (or `r.123ˆ124u`) → underscore form
- `p.123ˆ124Ala` → underscore form
- `p.His4Gln5insAla` → `p.His4_Gln5insAla`

### 2.11 Reference-sequence form

- `NM_004006` → `NM_004006.3` — version mandatory
- `NM_004006.2:c.357+1G>A` → `NC_000023.10(NM_004006.3):c.357+1G>A` — intronic needs genomic context
- `LRG_199:c.357+1G>A` → `LRG_199t1:c.357+1G>A`

### 2.12 Reformulating nonsense / extension / frameshift

- `p.Trp24TerfsTer1` → `p.Tyr4Ter` — nonsense, not frameshift
- `p.Tyr4_Cys5insTerGluAsp` → `p.Tyr4Ter` — nonsense, not insertion
- `p.Cys5_Ser6delinsTerGluAsp` → `p.Tyr4Ter` — nonsense, not delins
- `p.Trp26_Arg1623del` → `p.Trp26Ter` — nonsense, not deletion
- `p.His150HisfsTer10` → `p.Gln151Thrfs*9` — frameshift starts at first AA changed
- `p.Met1Valext-4` → `p.Met1_Leu2insArgSerThrVal` — Met1 change uses delins/ins, not extension
- `p.(Met3_Ile3418delinsGly)` → `p.(Met3_His4insGlyTer)` — stop codon in insertion shouldn't replace downstream
- `insSerSerTerAlaPro` → `NP_004371.2:p.(Pro46_Asn47insSerSerTer)` — drop AAs after stop
- `NM_024312.4:c.2686A[10]` → `NM_024312.4:c.2692_2693dup` — coding DNA, non-multiple-of-3 unit
- `NM_024312.4:c.1738TA[6]` → `NM_024312.4:c.1741_1742insTATATA` — frame-preserving exception

### 2.13 Repeat notation

- `g.123CAG[23]` → `g.123_191CAG[23]` — entire range required
- `r.123_191cag[23]` → recommended form per consultation (not paired)
- `r.123cug[23]` → `r.123_125[23]` — position-based form preferred
- `r.-123ug[14]` → `r.-124_-123[14]`
- `r.-125_-123cug[4]` → drop redundant sequence

### 2.14 Polymorphism / silent / "/" notation

- `NM_004006.1:c.76A/G` → use substitution form `c.76A>G` (split into two separate variants)
- `r.76a/g` → split into substitutions
- `p.2366Gln/Lys` → `NP_003997.1:p.Gln2366Lys`
- `p.Leu54Leu` → `p.(Leu54=)` (silent variant must use `=`)
- `p.54L/L` → `p.(Leu54=)`

---

## 3. Standalone non-compliant (no fix given on page)

These should fail to parse, or should be flagged — they are not
normalization-roundtrip cases.

- `chr1:g.1000_1005del` — `chr1` not a valid reference sequence
- `c.IVS4-2A>G` — exon/intron numbering in positions
- `g.234inv` — single-nucleotide inversion (must be substitution)
- `r.234inv` — same
- `c.4072-?_5154+?dup` — incorrect uncertain notation when exon involvement is known
- `NC_000023.11(NM_004006.2):c.(?_-244)_(31+1_32-1)dup` — extends beyond transcript
- `NC_000023.11(NM_004006.2):c.(10086+1_10087-1)_(*2691_?)dup` — extends beyond transcript
- `NC_000023.11(NM_004006.2):c.(-205839_-62966)_(*21568_*61692)dup` — outside transcript
- `c.4072-1234_5155-246dupXXXXX` — size suffix not allowed
- `c.1813-1dup` — wrong intron position
- `r.EX17del` — exon-numbering form not allowed
- `r.109u=` / `r.4567_4569=` — no-change must include sequence ID
- `p.Arg45del6` — count form not allowed for protein
- `p.EX17del` — exon-numbering form not allowed
- `p.(Met1Val)` — Met1 substitution not permitted
- `p.Ala11del` — deprecated for variable repeats; use `p.Ala2[9]`
- `p.Ala10_Ala11dup` — deprecated for variable repeats; use `p.Ala2[12]`
- `p.Tyr4*` — alternative-only notation (recommended is `p.Tyr4Ter`)
- `LRG_199t1:c.[2376G>C;3103=];[2376=;3103del]` — should not list "no change" alongside variants
- `c.2376[G>C];[G>C]` — forbidden shorthand
- `c.2376G>C[];[]` — forbidden abbreviation
- `c.[2376G>C];[=]` — `=` ambiguous (means whole reference analyzed)
- `c.2376G>A(;)3103del` — no allele brackets with uncertain phase
- `LRG_199t1:r.76a>u(;)(76a>u)` — parens around second occurrence wrong
- `r.[76a>u];[=]` — incomplete
- `NM_004006.2:c.[762_768del;767_774dup]` — overlapping del/dup of same ref not allowed
- `p.His4insAla` — missing flanking position
- `r.123insg`, `r.-14insG` — missing flanking position
- `g.123insG` — same
- `c.23ins24` — incomplete/inadequate sequence specification
- `NC_000002.12:g.pter_8247756::NC_000011.10:g.15825273_cen_qter` — ISCN2016, replaced by delins
- `NC_000011.10:g.pter_15825272::NC_000002.12:g.8247757_cen_qter` — same

---

## 4. Suggested test partitioning

When drafting tests:

1. **Roundtrip set (compliant unchanged)**: every entry in section 1.
2. **Normalization set (incorrect → correct)**: section 2 entries —
   especially 2.1 (3' rule), 2.2 (drop redundant seq), 2.4 (subst→delins),
   2.5 (ins→dup). These are what most normalizers actually do.
3. **Canonicalization-only set**: 2.3, 2.7, 2.9, 2.10, 2.11, 2.13, 2.14 —
   purely syntactic; may be parser/lexer-level rather than semantic
   normalization.
4. **Reject set**: section 3.

Some rows include partial reference prefixes (e.g. `c.76A>G` without
`NM_xxx.x:`). For tests that need a real reference, prepend a fixed test
reference or use the same one given elsewhere in the same spec page.
