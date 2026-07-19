# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.9.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.8.1...v0.9.0) - 2026-07-19

### Added

- opt-in protein p. render style (Ter vs * / 3- vs 1-letter) ([#1056](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1056))
- *(reference)* reject unknown manifest fields at load, fail loud on a misread reference ([#1055](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1055))

### Fixed

- render stop-spanning deletions/delins as frameshift when a sense residue changes ([#1073](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1073))
- derive projection is_frameshift from the protein consequence ([#1072](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1072))
- phase-aware frameshift aggregation with combined cis protein consequence ([#1071](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1071))
- report an allele's coordinate axis from HgvsVariant ([#1064](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1064))
- *(normalize)* preserve the uncertain/predicted (...) wrapper through normalization ([#1063](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1063)) ([#1065](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1065))
- *(hgvs)* make allele compaction gene-agnostic on a transcript reference ([#1062](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1062))
- *(hgvs)* classify species-qualified Ensembl transcripts uniformly ([#1061](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1061))
- *(normalize)* validate substitution reference bases and expose error_config to Python ([#1052](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1052)) ([#1060](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1060))
- *(hgvs)* drop the gene-symbol selector on transcript-reference Display ([#1054](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1054))
- *(convert-gff)* fail fast when the FASTA cannot cover an exon, instead of masking it ([#1053](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1053))
- *(normalize)* clamp the mitochondrial fetch window to contig length ([#1048](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1048))
- *(dev)* resolve spec submodule commit_sha under git hook env ([#1047](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1047))
- *(normalize)* clamp genome fetch window to contig length ([#1042](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1042))

### Other

- build generate_spec_fixture example into the outer target dir ([#1074](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1074))
- bump minor version for feature commits on 0.x releases ([#1066](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1066))
- *(conformance)* update issue_506 expectation for #1054 transcript-reference selector drop ([#1067](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1067))
- *(normalize)* probe inv over-recognition boundary for issue #1040 ([#1043](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1043))

## [0.8.1](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.8.0...v0.8.1) - 2026-07-16

### Added

- *(python)* expose convert_gff and build_transcript, mirroring prepare_reference_data ([#1037](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1037))

### Fixed

- *(normalize)* only type a maximal contiguous run as inv, not a sub-run ([#1036](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1036))

## [0.8.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.7.1...v0.8.0) - 2026-07-15

### Added

- *(reference,cli)* validate transcripts.json version and check it standalone ([#1031](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1031))
- *(normalize)* warn-and-degrade when genomic data is unavailable ([#1015](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1015))
- *(cli)* genome-capable transcripts.json via --emit-genomic-sequences ([#1029](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1029))
- *(python)* value equality on projections/effects, idempotency tests, typed-Raises docs ([#1022](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1022))
- *(python)* surface normalization warnings on projections ([#1021](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1021))
- *(python)* batch error-isolation on the projector (return_exceptions) ([#1020](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1020))
- *(python)* typed exception hierarchy, py.typed, and structured error codes ([#1019](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1019))
- *(python)* surface reference capability and warn on reduced-capability build ([#1014](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1014))

### Fixed

- *(python)* reject an unrecognized direction argument instead of defaulting to 3' ([#1017](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1017))

### Other

- make public error/projection API robust to additive changes ([#1033](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1033))
- *(pre-commit)* fix the broken mypy hook (bump to v2.2.0, check by directory) ([#1023](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1023))
- *(reference)* rename MockProvider to JsonProvider ([#1030](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1030))
- *(readme)* add PyPI version and Python-versions badges ([#1010](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1010))

## [0.7.1](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.7.0...v0.7.1) - 2026-07-09

### Fixed

- *(ci)* publish wheels directly to PyPI and drop the TestPyPI gate that had failed on every prior release, so Python wheels are published to PyPI for the first time ([#1008](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1008))

## [0.7.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.6.0...v0.7.0) - 2026-07-09

### Added

- *(reference)* validate manifest schema/version at load, fail loud on an incompatible reference ([#1003](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1003))
- *(mosaic)* parse predicted-wrapper and whole-entity-LHS =/ forms ([#992](https://github.com/fulcrumgenomics/ferro-hgvs/pull/992))
- *(protein)* parse insXaa[n] and insTer<n>/ins*<n> inserted sequences ([#991](https://github.com/fulcrumgenomics/ferro-hgvs/pull/991))
- *(python)* forward the remaining MultiFasta-overridden methods to the Python provider ([#977](https://github.com/fulcrumgenomics/ferro-hgvs/pull/977))
- *(python)* forward legacy gene-selector resolution to the Python provider ([#967](https://github.com/fulcrumgenomics/ferro-hgvs/pull/967))
- feat(project)+test(conformance): source-scope enumeration & provision Ensembl reference ([#933](https://github.com/fulcrumgenomics/ferro-hgvs/pull/933)) ([#942](https://github.com/fulcrumgenomics/ferro-hgvs/pull/942))
- *(normalize)* synthesize the transcript selector on a bare-NG_ c. input ([#923](https://github.com/fulcrumgenomics/ferro-hgvs/pull/923)) ([#932](https://github.com/fulcrumgenomics/ferro-hgvs/pull/932))
- *(python)* normalize project_to_genomic output by default; raw via normalize=False ([#888](https://github.com/fulcrumgenomics/ferro-hgvs/pull/888))
- *(project)* NG_/LRG in-place c./n.→g. projection for compound alleles and UTR coordinates ([#881](https://github.com/fulcrumgenomics/ferro-hgvs/pull/881))
- *(project)* predict C-terminal extension / stop-loss for CDS→3'UTR-spanning deletion ([#880](https://github.com/fulcrumgenomics/ferro-hgvs/pull/880))
- *(reference)* provision NM_002001.2 conformance rows ([#802](https://github.com/fulcrumgenomics/ferro-hgvs/pull/802)) ([#849](https://github.com/fulcrumgenomics/ferro-hgvs/pull/849))
- *(prepare)* version-aware transcript backfill ([#842](https://github.com/fulcrumgenomics/ferro-hgvs/pull/842)) ([#848](https://github.com/fulcrumgenomics/ferro-hgvs/pull/848))
- *(service/effect)* real amino-acid resolution and junction-based NMD ([#837](https://github.com/fulcrumgenomics/ferro-hgvs/pull/837))
- *(vcf)* RNA→VCF conversion and complex-allele decomposition ([#822](https://github.com/fulcrumgenomics/ferro-hgvs/pull/822))
- *(vcf)* parse per-sample FORMAT genotype values ([#823](https://github.com/fulcrumgenomics/ferro-hgvs/pull/823))
- *(project)* translate downstream variants on non-AUG-initiation transcripts ([#780](https://github.com/fulcrumgenomics/ferro-hgvs/pull/780)) ([#796](https://github.com/fulcrumgenomics/ferro-hgvs/pull/796))
- *(reference)* provision exon→genome structure for cdot-absent old transcript versions (NM_003002.2) ([#790](https://github.com/fulcrumgenomics/ferro-hgvs/pull/790)) ([#795](https://github.com/fulcrumgenomics/ferro-hgvs/pull/795))
- *(project)* resolve legacy GENE_vNNN selectors on NG_(...):c./n./r. projection ([#784](https://github.com/fulcrumgenomics/ferro-hgvs/pull/784))
- *(project)* predict p.? for an initiation-codon member of a cis compound allele ([#778](https://github.com/fulcrumgenomics/ferro-hgvs/pull/778))
- *(normalize)* expand r. (RNA) cross-reference insert/delins payloads ([#777](https://github.com/fulcrumgenomics/ferro-hgvs/pull/777))
- *(normalize)* expand axis-less cross-reference insert payloads ([#769](https://github.com/fulcrumgenomics/ferro-hgvs/pull/769))
- *(project)* predict protein consequence for whole-exon deletions ([#761](https://github.com/fulcrumgenomics/ferro-hgvs/pull/761))
- *(conformance)* corpus-driven NG_ placement derivation ([#744](https://github.com/fulcrumgenomics/ferro-hgvs/pull/744)) ([#747](https://github.com/fulcrumgenomics/ferro-hgvs/pull/747))
- *(conformance)* hermetic genomic-axis PR gate ([#737](https://github.com/fulcrumgenomics/ferro-hgvs/pull/737)) ([#746](https://github.com/fulcrumgenomics/ferro-hgvs/pull/746))
- *(prepare)* produce derived_refseqgene_placements via --derive-ng-placements ([#740](https://github.com/fulcrumgenomics/ferro-hgvs/pull/740)) ([#743](https://github.com/fulcrumgenomics/ferro-hgvs/pull/743))
- *(project)* emit the g. axis for NG/LRG/NC-parent c.pter/c.qter inputs ([#537](https://github.com/fulcrumgenomics/ferro-hgvs/pull/537)) ([#739](https://github.com/fulcrumgenomics/ferro-hgvs/pull/739))
- *(conformance)* hermetic protein-axis PR gate ([#719](https://github.com/fulcrumgenomics/ferro-hgvs/pull/719)) ([#733](https://github.com/fulcrumgenomics/ferro-hgvs/pull/733))
- *(conformance)* hermetic snapshot-backed normalized-axis PR gate (#719 I4) ([#729](https://github.com/fulcrumgenomics/ferro-hgvs/pull/729))
- *(reference)* provision derived NG_ placements — load + produce ([#728](https://github.com/fulcrumgenomics/ferro-hgvs/pull/728)) ([#732](https://github.com/fulcrumgenomics/ferro-hgvs/pull/732))
- *(reference)* derive version-independent NG_ genomic placement — core ([#728](https://github.com/fulcrumgenomics/ferro-hgvs/pull/728)) ([#731](https://github.com/fulcrumgenomics/ferro-hgvs/pull/731))
- *(conformance)* pinned version-exact transcript reference snapshot (#719 I2) ([#727](https://github.com/fulcrumgenomics/ferro-hgvs/pull/727))
- *(normalize,project)* resolve legacy gene-model selectors on genomic references ([#709](https://github.com/fulcrumgenomics/ferro-hgvs/pull/709))
- *(reference)* cross-build version-exact cdot fallback for c. normalization ([#720](https://github.com/fulcrumgenomics/ferro-hgvs/pull/720))
- *(project)* curate project_*_all enumeration — collapse superseded versions, prefer curated transcripts ([#710](https://github.com/fulcrumgenomics/ferro-hgvs/pull/710))
- *(project,python)* --assembly genome-build override for build-agnostic NG_/LRG_ inputs ([#723](https://github.com/fulcrumgenomics/ferro-hgvs/pull/723))
- *(reference,project)* bundle GRCh37 RefSeqGene alignment + thread build through projection ([#721](https://github.com/fulcrumgenomics/ferro-hgvs/pull/721))
- *(reference)* build-aware NG_/LRG_ genomic placement (GRCh37 + GRCh38) foundation ([#711](https://github.com/fulcrumgenomics/ferro-hgvs/pull/711))
- *(conformance)* add reference-gap report tool ([#722](https://github.com/fulcrumgenomics/ferro-hgvs/pull/722))
- *(cli,project)* ferro project --axis {g,c,n,p,r} multi-axis projection surface ([#712](https://github.com/fulcrumgenomics/ferro-hgvs/pull/712))
- *(conformance)* add corpus accession inventory ([#707](https://github.com/fulcrumgenomics/ferro-hgvs/pull/707))
- *(conformance)* add reference_unavailable disposition ([#706](https://github.com/fulcrumgenomics/ferro-hgvs/pull/706))
- *(project)* add c.→r. RNA consequence prediction surface ([#701](https://github.com/fulcrumgenomics/ferro-hgvs/pull/701))
- *(normalize)* apply the 3′ rule across exon/intron boundaries ([#670](https://github.com/fulcrumgenomics/ferro-hgvs/pull/670)) ([#700](https://github.com/fulcrumgenomics/ferro-hgvs/pull/700))
- *(normalize)* emit the NP_ protein-accession selector for p. on genomic references ([#502](https://github.com/fulcrumgenomics/ferro-hgvs/pull/502)) ([#676](https://github.com/fulcrumgenomics/ferro-hgvs/pull/676))
- *(prepare)* add Ensembl reference support (ENST/ENSG/ENSP) ([#677](https://github.com/fulcrumgenomics/ferro-hgvs/pull/677))
- *(check)* validate CDS start-codon consistency in `ferro check` ([#629](https://github.com/fulcrumgenomics/ferro-hgvs/pull/629)) ([#674](https://github.com/fulcrumgenomics/ferro-hgvs/pull/674))
- *(normalize)* collapse overlapping cis-allele edits into a single delins ([#667](https://github.com/fulcrumgenomics/ferro-hgvs/pull/667))
- *(project)* extend the direct c.→p. path to bare n./r. transcript inputs ([#506](https://github.com/fulcrumgenomics/ferro-hgvs/pull/506)) ([#661](https://github.com/fulcrumgenomics/ferro-hgvs/pull/661))
- *(project)* enumerate + frame c./n./r. inputs with NG_/LRG_ context ([#646](https://github.com/fulcrumgenomics/ferro-hgvs/pull/646)) ([#647](https://github.com/fulcrumgenomics/ferro-hgvs/pull/647))
- *(project)* project NG_/LRG_ genomic inputs onto transcripts ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#642](https://github.com/fulcrumgenomics/ferro-hgvs/pull/642))
- *(project)* re-anchor NG_/LRG_-parent projections into the parent frame ([#480](https://github.com/fulcrumgenomics/ferro-hgvs/pull/480)) ([#616](https://github.com/fulcrumgenomics/ferro-hgvs/pull/616))
- *(benchmark)* perf-matrix engine + measured performance tables ([#617](https://github.com/fulcrumgenomics/ferro-hgvs/pull/617))
- *(parser)* parse and/or between bracketed alleles ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#563](https://github.com/fulcrumgenomics/ferro-hgvs/pull/563))
- *(parser)* parse predicted compound cis allele [(a;b)] ([#545](https://github.com/fulcrumgenomics/ferro-hgvs/pull/545)) ([#569](https://github.com/fulcrumgenomics/ferro-hgvs/pull/569))
- *(parser)* parse trans alleles with a mosaic (;) tail ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#559](https://github.com/fulcrumgenomics/ferro-hgvs/pull/559))
- *(perf-table)* render README performance tables from a results JSON ([#604](https://github.com/fulcrumgenomics/ferro-hgvs/pull/604))
- *(reference)* load + index the GRCh37 cdot so multi-build manifests project GRCh37 inputs ([#605](https://github.com/fulcrumgenomics/ferro-hgvs/pull/605))
- *(error)* emit AlignmentGap (E3007) for variants in a transcript-genome CIGAR indel (Cycle 1c-iii) ([#603](https://github.com/fulcrumgenomics/ferro-hgvs/pull/603))
- *(parser)* parse N-padded deletion over an uncertain range ([#545](https://github.com/fulcrumgenomics/ferro-hgvs/pull/545)) ([#568](https://github.com/fulcrumgenomics/ferro-hgvs/pull/568))
- *(parser)* parse comma products allele [a,b] ([#545](https://github.com/fulcrumgenomics/ferro-hgvs/pull/545)) ([#571](https://github.com/fulcrumgenomics/ferro-hgvs/pull/571))
- *(parser)* parse positionless breakpoint insertion in composite alleles ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#584](https://github.com/fulcrumgenomics/ferro-hgvs/pull/584))
- *(tool-support)* single source of truth for the tool-comparison tables ([#590](https://github.com/fulcrumgenomics/ferro-hgvs/pull/590))
- *(parser)* parse the sup supernumerary marker ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#596](https://github.com/fulcrumgenomics/ferro-hgvs/pull/596))
- *(parser)* accept pter/qter endpoints in inserted position-ranges ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#597](https://github.com/fulcrumgenomics/ferro-hgvs/pull/597))
- *(parser)* parse same-chromosome ring `::` deletion-join ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#595](https://github.com/fulcrumgenomics/ferro-hgvs/pull/595))
- *(parser)* support inverted insertion of an uncertain-boundary range ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#576](https://github.com/fulcrumgenomics/ferro-hgvs/pull/576))
- *(parser)* accept all-identity trans alleles via cross-spell check ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#558](https://github.com/fulcrumgenomics/ferro-hgvs/pull/558))
- *(project)* add the c↔n (non-coding transcript) axis to VariantProjection ([#592](https://github.com/fulcrumgenomics/ferro-hgvs/pull/592))
- *(parser)* parse shared-position repeat trans alleles ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#555](https://github.com/fulcrumgenomics/ferro-hgvs/pull/555))
- *(parser)* parse predicted-wrapper members in protein alleles ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#552](https://github.com/fulcrumgenomics/ferro-hgvs/pull/552))
- *(parser)* parse residue-level and/or `^` substitutions ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#550](https://github.com/fulcrumgenomics/ferro-hgvs/pull/550))
- *(parser)* parse frameshift with alternative new residues ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#557](https://github.com/fulcrumgenomics/ferro-hgvs/pull/557))
- *(parser)* parse trans alleles whose members are cis groups ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#549](https://github.com/fulcrumgenomics/ferro-hgvs/pull/549))
- *(parser)* extend mosaic `=/` to protein + fix RNA case ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#548](https://github.com/fulcrumgenomics/ferro-hgvs/pull/548))
- *(parser)* parse the and/or `^` operator ([#544](https://github.com/fulcrumgenomics/ferro-hgvs/pull/544)) ([#547](https://github.com/fulcrumgenomics/ferro-hgvs/pull/547))
- *(conformance)* generate failure-patterns.md from cases.json; retire hand-maintained counts ([#509](https://github.com/fulcrumgenomics/ferro-hgvs/pull/509)) ([#538](https://github.com/fulcrumgenomics/ferro-hgvs/pull/538))
- *(reference)* reconcile cdot for synthesis-only transcripts via exon-sum length ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#532](https://github.com/fulcrumgenomics/ferro-hgvs/pull/532))
- *(reference)* ingest protein FASTAs + get_protein_sequence + translation check ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#530](https://github.com/fulcrumgenomics/ferro-hgvs/pull/530))
- *(reference)* apply canonical overrides in MultiFastaProvider + reconcile cdot ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#528](https://github.com/fulcrumgenomics/ferro-hgvs/pull/528))
- *(prepare)* canonical-overrides acquisition — ferro prepare --validate-canonical ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#527](https://github.com/fulcrumgenomics/ferro-hgvs/pull/527))
- *(normalize)* resolve bare-transcript c. pter/qter telomere markers (#488 Phase 2) ([#533](https://github.com/fulcrumgenomics/ferro-hgvs/pull/533))
- *(reference)* version-exact correction of transcript metadata ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#524](https://github.com/fulcrumgenomics/ferro-hgvs/pull/524))
- *(reference)* offline structural validation of transcript records ([#520](https://github.com/fulcrumgenomics/ferro-hgvs/pull/520)) ([#521](https://github.com/fulcrumgenomics/ferro-hgvs/pull/521))
- *(normalize)* resolve genomic pter/qter telomere markers to concrete coordinates ([#488](https://github.com/fulcrumgenomics/ferro-hgvs/pull/488)) ([#526](https://github.com/fulcrumgenomics/ferro-hgvs/pull/526))
- *(project)* direct c.→p. path for bare-NM_ inputs; demote 4 bare-NP protein rows ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#508](https://github.com/fulcrumgenomics/ferro-hgvs/pull/508))
- *(reference)* strict exact-version transcript resolution (#478 pillar 3) ([#490](https://github.com/fulcrumgenomics/ferro-hgvs/pull/490))
- *(error-handling)* reject standalone single-member allele brackets uniformly (W3026, closes #493) ([#496](https://github.com/fulcrumgenomics/ferro-hgvs/pull/496))
- *(error-handling)* add W3023/W3024/W3025 for dup<N>/dup<seq>/del<seq> soft-prohibition forms (closes #460) ([#462](https://github.com/fulcrumgenomics/ferro-hgvs/pull/462))
- *(parser)* accept whole-entity edits in protein bracket members (closes #468) ([#475](https://github.com/fulcrumgenomics/ferro-hgvs/pull/475))
- *(normalize)* expand ins[ACC:g.A_B] cross-reference to flat literal (closes #422) ([#437](https://github.com/fulcrumgenomics/ferro-hgvs/pull/437))
- *(normalize)* canonicalize p.delins → p.dup and surface W3022 InitiatorMetCanonicalization (closes #92) ([#407](https://github.com/fulcrumgenomics/ferro-hgvs/pull/407))
- contig-length-aware math + wraparound boundary policies for m./o. variants (closes #399) ([#411](https://github.com/fulcrumgenomics/ferro-hgvs/pull/411))
- *(normalize,project,python,validate)* four misc gaps from #395 (items 1, 4, 5, 6) ([#421](https://github.com/fulcrumgenomics/ferro-hgvs/pull/421))
- *(parser)* spec-compliant wraparound ranges on m./o. circular refs ([#380](https://github.com/fulcrumgenomics/ferro-hgvs/pull/380))
- *(normalize)* protein 3'-shifting for deletions and duplications ([#377](https://github.com/fulcrumgenomics/ferro-hgvs/pull/377))
- *(normalize)* expand ins[...] to flat literal sequence (closes #333) ([#365](https://github.com/fulcrumgenomics/ferro-hgvs/pull/365))
- *(normalize)* structured info-code surface mirroring W#### warnings ([#373](https://github.com/fulcrumgenomics/ferro-hgvs/pull/373))
- *(project)* VariantProjector accepts c./n./r. inputs ([#379](https://github.com/fulcrumgenomics/ferro-hgvs/pull/379))
- *(error_handling)* mutalyzer ↔ ferro error-code mapping table (closes #329) ([#362](https://github.com/fulcrumgenomics/ferro-hgvs/pull/362))
- *(error-handling)* W5003 VariantExceedsReference, strict-mode rejection (closes #355) ([#363](https://github.com/fulcrumgenomics/ferro-hgvs/pull/363))
- *(test/mutalyzer-normalize)* accepted_divergence + spec_citation (closes #335) ([#359](https://github.com/fulcrumgenomics/ferro-hgvs/pull/359))
- *(parser)* canonicalize r. thymine input to u with W3020 soft-warning (closes #282) ([#299](https://github.com/fulcrumgenomics/ferro-hgvs/pull/299))
- *(parser)* recognize non-human Ensembl accession prefixes ([#320](https://github.com/fulcrumgenomics/ferro-hgvs/pull/320))
- *(project)* project_to_genomic for c./n./r. inputs (closes #327) ([#358](https://github.com/fulcrumgenomics/ferro-hgvs/pull/358))
- *(parser)* RNA/TX trans-allele predicted-edit wrapper in compound brackets (closes #287) ([#309](https://github.com/fulcrumgenomics/ferro-hgvs/pull/309))
- *(parser)* whole-entity predicted forms for m./o./n. variants (closes #288) ([#308](https://github.com/fulcrumgenomics/ferro-hgvs/pull/308))
- *(normalize)* MultiRepeat partial validation for non-Exact counts (closes #279) ([#296](https://github.com/fulcrumgenomics/ferro-hgvs/pull/296))
- *(annotation)* populate Transcript.protein_id from GFF/GTF attributes ([#321](https://github.com/fulcrumgenomics/ferro-hgvs/pull/321))
- *(parser)* whole-entity predicted forms `c.(=)`, `c.(?)`, `r.(=)`, `r.(?)`, `r.(0)` (closes #245) ([#246](https://github.com/fulcrumgenomics/ferro-hgvs/pull/246))
- *(parser)* predicted-edit wrapper inside compound brackets (closes #243) ([#244](https://github.com/fulcrumgenomics/ferro-hgvs/pull/244))
- *(parser)* accept ?con<src> and ?copy<N> on unknown position (closes #286) ([#303](https://github.com/fulcrumgenomics/ferro-hgvs/pull/303))
- *(parser)* bracket-aware split for chimeric + non-compact mosaic forms (closes #216) ([#217](https://github.com/fulcrumgenomics/ferro-hgvs/pull/217))
- *(normalize)* validate repeat unit divides reference span and matches bases (closes #214) ([#215](https://github.com/fulcrumgenomics/ferro-hgvs/pull/215))

### Fixed

- *(ci)* satisfy Rust 1.97 clippy lints and pin the toolchain ([#1007](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1007))
- *(normalize)* preserve co-located insertion overlaps instead of a non-idempotent merge ([#1005](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1005))
- *(normalize)* coalesce shift-created cis adjacency into a single delins ([#1002](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1002))
- *(project)* decline c./p./r. on 5′-incomplete-CDS (cds_start_NF) transcripts ([#972](https://github.com/fulcrumgenomics/ferro-hgvs/pull/972)) ([#994](https://github.com/fulcrumgenomics/ferro-hgvs/pull/994))
- *(reference)* resolve transcript_genome_span for deferred Ensembl transcripts ([#995](https://github.com/fulcrumgenomics/ferro-hgvs/pull/995)) ([#996](https://github.com/fulcrumgenomics/ferro-hgvs/pull/996))
- *(normalize)* render plain past-CDS coordinates in canonical c.*N form (#920/#336) ([#987](https://github.com/fulcrumgenomics/ferro-hgvs/pull/987))
- *(spdi)* bound exonic c./n./r. positions past the 3' end; resolve coding r.N as CDS-relative ([#989](https://github.com/fulcrumgenomics/ferro-hgvs/pull/989))
- *(convert)* CDS<->tx is transcript-native across intra-exon CIGAR insertions ([#944](https://github.com/fulcrumgenomics/ferro-hgvs/pull/944)) ([#988](https://github.com/fulcrumgenomics/ferro-hgvs/pull/988))
- *(normalize)* canonicalize bracketed inserted-repeat counts >= 2 ([#920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/920)) ([#983](https://github.com/fulcrumgenomics/ferro-hgvs/pull/983))
- *(parser)* validate compound-reference accession pairings (outer=gene/genomic) ([#981](https://github.com/fulcrumgenomics/ferro-hgvs/pull/981))
- *(effect)* stop_retained, codon-based start_lost, and deterministic most_severe ([#978](https://github.com/fulcrumgenomics/ferro-hgvs/pull/978))
- *(spdi)* reject r.*N/c.*N positions past the transcript 3' end ([#970](https://github.com/fulcrumgenomics/ferro-hgvs/pull/970))
- *(spdi)* resolve r.-N (5'UTR) through the exon-aware mapper so it agrees with c.-N ([#969](https://github.com/fulcrumgenomics/ferro-hgvs/pull/969))
- *(vcf)* decline unknown-position/edit c./n./g. variants instead of panicking ([#943](https://github.com/fulcrumgenomics/ferro-hgvs/pull/943)) ([#947](https://github.com/fulcrumgenomics/ferro-hgvs/pull/947))
- *(project)* make c.→p. consequence decline observable, not silent ([#956](https://github.com/fulcrumgenomics/ferro-hgvs/pull/956))
- *(parser)* resolve ENSG(ENST) compound refs; require versions on Ensembl accessions ([#933](https://github.com/fulcrumgenomics/ferro-hgvs/pull/933)) ([#938](https://github.com/fulcrumgenomics/ferro-hgvs/pull/938))
- *(normalize)* reconcile MultiRepeat needs_normalization with the shuffle dispatch ([#958](https://github.com/fulcrumgenomics/ferro-hgvs/pull/958))
- *(spdi)* resolve r.*N through the exon-aware mapper so it agrees with c.*N ([#944](https://github.com/fulcrumgenomics/ferro-hgvs/pull/944)) ([#950](https://github.com/fulcrumgenomics/ferro-hgvs/pull/950))
- *(parser)* correct self-cancelling-allele guard docs and fixture classification ([#959](https://github.com/fulcrumgenomics/ferro-hgvs/pull/959))
- *(prepare)* make --dry-run a true no-op preview ([#939](https://github.com/fulcrumgenomics/ferro-hgvs/pull/939)) ([#940](https://github.com/fulcrumgenomics/ferro-hgvs/pull/940))
- *(normalize)* collapse overlapping cis del+ins allele members on the transcript axis ([#920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/920)) ([#929](https://github.com/fulcrumgenomics/ferro-hgvs/pull/929))
- *(normalize)* 3'-shift CDS and whole-CDS-spanning del/dup into the UTR per the HGVS 3'-rule ([#918](https://github.com/fulcrumgenomics/ferro-hgvs/pull/918)) ([#935](https://github.com/fulcrumgenomics/ferro-hgvs/pull/935))
- *(normalize)* canonicalize degenerate bracketed inserted-repeat counts ([#920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/920)) ([#926](https://github.com/fulcrumgenomics/ferro-hgvs/pull/926))
- *(project)* render a codon-straddling single-codon in-frame deletion as del, not empty delins ([#931](https://github.com/fulcrumgenomics/ferro-hgvs/pull/931))
- *(benchmark)* surface the failure reason when ferro normalize yields zero successes ([#916](https://github.com/fulcrumgenomics/ferro-hgvs/pull/916))
- *(project)* render premature-stop protein consequences per spec, not a C-terminal delins ([#911](https://github.com/fulcrumgenomics/ferro-hgvs/pull/911)) ([#913](https://github.com/fulcrumgenomics/ferro-hgvs/pull/913))
- *(reference)* skip hidden/AppleDouble files when scanning FASTA directories ([#915](https://github.com/fulcrumgenomics/ferro-hgvs/pull/915))
- *(project)* cis-sort minus-strand compound-allele members on the project_variant genomic axis ([#894](https://github.com/fulcrumgenomics/ferro-hgvs/pull/894)) ([#898](https://github.com/fulcrumgenomics/ferro-hgvs/pull/898))
- *(project)* de-anchor bare NG_/LRG_ genomic input on the single-transcript projection paths ([#879](https://github.com/fulcrumgenomics/ferro-hgvs/pull/879)) ([#896](https://github.com/fulcrumgenomics/ferro-hgvs/pull/896))
- *(project)* surface the genomic axis for LRG inputs and pter/qter termini on project_variant (#886, #887) ([#893](https://github.com/fulcrumgenomics/ferro-hgvs/pull/893))
- *(normalize)* keep out-of-phase insertions as ins, not spurious dup/repeat ([#882](https://github.com/fulcrumgenomics/ferro-hgvs/pull/882)) ([#892](https://github.com/fulcrumgenomics/ferro-hgvs/pull/892))
- *(normalize)* 3'-align insertion→duplication; guarantee normalization idempotency ([#883](https://github.com/fulcrumgenomics/ferro-hgvs/pull/883))
- *(project)* preserve LRG_t/LRG_p namespace on projection output ([#860](https://github.com/fulcrumgenomics/ferro-hgvs/pull/860)) ([#874](https://github.com/fulcrumgenomics/ferro-hgvs/pull/874))
- render in-frame single-codon deletion in a residue run as del, not empty ins ([#850](https://github.com/fulcrumgenomics/ferro-hgvs/pull/850))
- *(project)* spec-canonical minus-strand genomic repeat projection ([#852](https://github.com/fulcrumgenomics/ferro-hgvs/pull/852)) ([#869](https://github.com/fulcrumgenomics/ferro-hgvs/pull/869))
- *(parser)* collapse single-element bracketed literal ins/delins to a plain literal ([#863](https://github.com/fulcrumgenomics/ferro-hgvs/pull/863))
- *(project)* canonicalize protein consequences (delins→sub, tandem dup, adjacent in-cis delins) ([#862](https://github.com/fulcrumgenomics/ferro-hgvs/pull/862))
- *(project)* build-scope the predicted-r. allele transcript fetch ([#843](https://github.com/fulcrumgenomics/ferro-hgvs/pull/843)) ([#845](https://github.com/fulcrumgenomics/ferro-hgvs/pull/845))
- *(normalize)* treat LRG transcript references as bare for EINTRONIC ([#834](https://github.com/fulcrumgenomics/ferro-hgvs/pull/834)) ([#844](https://github.com/fulcrumgenomics/ferro-hgvs/pull/844))
- *(project)* genomic projection of c.* variants in the 3'UTR poly-A region (NM_003002.2; #790 residual) ([#839](https://github.com/fulcrumgenomics/ferro-hgvs/pull/839))
- *(project)* wire rna_description axis and frame predicted r. from input ([#838](https://github.com/fulcrumgenomics/ferro-hgvs/pull/838))
- *(vcf)* drop per-sample genotype when splitting multi-allelic records ([#841](https://github.com/fulcrumgenomics/ferro-hgvs/pull/841))
- *(reference)* apply cdot CIGAR during transcript synthesis, decline on insertions ([#831](https://github.com/fulcrumgenomics/ferro-hgvs/pull/831))
- *(project)* force initiator Met at residue 1 for non-AUG ref-protein translation ([#836](https://github.com/fulcrumgenomics/ferro-hgvs/pull/836))
- *(reference)* suppress derived-tx injection for the #790 producer ([#835](https://github.com/fulcrumgenomics/ferro-hgvs/pull/835))
- *(reference)* cover canonical-sequence reingestion through the provider ([#791](https://github.com/fulcrumgenomics/ferro-hgvs/pull/791)) ([#832](https://github.com/fulcrumgenomics/ferro-hgvs/pull/832))
- *(vcf)* correct genomic anchor-base recovery and consolidate NC_→chrom mapping ([#821](https://github.com/fulcrumgenomics/ferro-hgvs/pull/821))
- *(vcf)* emit RefSeq NC_ accessions on the genomic HGVS axis ([#820](https://github.com/fulcrumgenomics/ferro-hgvs/pull/820))
- *(equivalence/service)* make unsupported variant classes explicit ([#816](https://github.com/fulcrumgenomics/ferro-hgvs/pull/816))
- *(data/cdot)* drop frequently-wrong NM_→NP_ protein-accession inference ([#815](https://github.com/fulcrumgenomics/ferro-hgvs/pull/815))
- *(vcf)* classify coding indels by typed NaEdit and emit Frameshift ([#817](https://github.com/fulcrumgenomics/ferro-hgvs/pull/817))
- *(error)* give TranscriptVersionNotExact a distinct error code (E2004) ([#814](https://github.com/fulcrumgenomics/ferro-hgvs/pull/814))
- *(normalize)* resolve legacy GENE_v001 selector against the NG_ parent's hosted transcript ([#792](https://github.com/fulcrumgenomics/ferro-hgvs/pull/792)) ([#793](https://github.com/fulcrumgenomics/ferro-hgvs/pull/793))
- decline silent transcript version substitution in normalization and projection ([#787](https://github.com/fulcrumgenomics/ferro-hgvs/pull/787))
- *(project)* report p.? for an initiation-codon variant on a non-AUG transcript ([#772](https://github.com/fulcrumgenomics/ferro-hgvs/pull/772))
- *(data)* anchor minus-strand intronic offsets to the correct exon end ([#766](https://github.com/fulcrumgenomics/ferro-hgvs/pull/766))
- *(convert)* map c.*N+offset 3'UTR positions past a transcript terminus ([#765](https://github.com/fulcrumgenomics/ferro-hgvs/pull/765))
- *(normalize)* resolve the transcript for a genomic-context c. cross-reference payload ([#654](https://github.com/fulcrumgenomics/ferro-hgvs/pull/654)) ([#756](https://github.com/fulcrumgenomics/ferro-hgvs/pull/756))
- *(normalize)* reject overlapping insertions in cis alleles (#486 EOVERLAP) ([#749](https://github.com/fulcrumgenomics/ferro-hgvs/pull/749))
- *(data)* load cdot exons in the HGVS coordinate convention ([#748](https://github.com/fulcrumgenomics/ferro-hgvs/pull/748))
- *(project)* re-normalize projected genomic variants in their own frame ([#741](https://github.com/fulcrumgenomics/ferro-hgvs/pull/741))
- *(normalize)* normalize r. edits carrying the RNA base `u` ([#736](https://github.com/fulcrumgenomics/ferro-hgvs/pull/736)) ([#738](https://github.com/fulcrumgenomics/ferro-hgvs/pull/738))
- *(normalize)* extend the exon/intron 3′ rule to n. and r. variants ([#704](https://github.com/fulcrumgenomics/ferro-hgvs/pull/704)) ([#734](https://github.com/fulcrumgenomics/ferro-hgvs/pull/734))
- *(parser)* accept a parenthesized gene-selector in a cross-reference insert ([#730](https://github.com/fulcrumgenomics/ferro-hgvs/pull/730))
- *(reference)* forward dropped methods through Arc/Box provider wrappers ([#726](https://github.com/fulcrumgenomics/ferro-hgvs/pull/726))
- *(project)* resolve range-reference insertions against the NG_/LRG_ parent before de-anchoring ([#708](https://github.com/fulcrumgenomics/ferro-hgvs/pull/708))
- *(reference)* resolve cdot version-exact in normalize, never substitute a sibling version ([#717](https://github.com/fulcrumgenomics/ferro-hgvs/pull/717))
- *(project)* re-anchor the genomic axis on the single-transcript projection path ([#702](https://github.com/fulcrumgenomics/ferro-hgvs/pull/702)) ([#703](https://github.com/fulcrumgenomics/ferro-hgvs/pull/703))
- *(normalize)* correct 1-based↔0-based off-by-one in genomic-shuffle fetch ([#690](https://github.com/fulcrumgenomics/ferro-hgvs/pull/690))
- *(normalize)* keep W4007 on warn-only intronic-bare resolve-failure ([#686](https://github.com/fulcrumgenomics/ferro-hgvs/pull/686))
- *(project)* error instead of emitting invalid HGVS when an NG_/LRG_ re-anchor declines ([#655](https://github.com/fulcrumgenomics/ferro-hgvs/pull/655)) ([#662](https://github.com/fulcrumgenomics/ferro-hgvs/pull/662))
- *(project)* sequence-aware g.↔c. projection across unmodeled exon indels ([#644](https://github.com/fulcrumgenomics/ferro-hgvs/pull/644)) ([#668](https://github.com/fulcrumgenomics/ferro-hgvs/pull/668))
- *(reference)* don't abort LRG mapping scan on a malformed earlier `<mapping>` ([#684](https://github.com/fulcrumgenomics/ferro-hgvs/pull/684))
- *(check)* warn when `refseqgene_alignments` file is missing ([#683](https://github.com/fulcrumgenomics/ferro-hgvs/pull/683))
- *(python)* add missing ErrorType discriminants 41/42 to the type stub ([#658](https://github.com/fulcrumgenomics/ferro-hgvs/pull/658))
- *(normalize)* spec-compliant repeat normalization — 3' rotation + flank absorption ([#649](https://github.com/fulcrumgenomics/ferro-hgvs/pull/649))
- *(reference)* resolve bare LRG_N genomic accession to its LRG_Ng record ([#638](https://github.com/fulcrumgenomics/ferro-hgvs/pull/638))
- *(benchmark)* make `setup uta` reliable ([#623](https://github.com/fulcrumgenomics/ferro-hgvs/pull/623)) ([#632](https://github.com/fulcrumgenomics/ferro-hgvs/pull/632))
- *(project)* decline protein prediction when the CDS frame is inconsistent ([#625](https://github.com/fulcrumgenomics/ferro-hgvs/pull/625)) ([#628](https://github.com/fulcrumgenomics/ferro-hgvs/pull/628))
- *(validate)* reject intronic offsets on a bare transcript reference (#486 EINTRONIC) ([#577](https://github.com/fulcrumgenomics/ferro-hgvs/pull/577))
- *(project)* route stop-disrupting indels to a C-terminal extension ([#615](https://github.com/fulcrumgenomics/ferro-hgvs/pull/615)) ([#621](https://github.com/fulcrumgenomics/ferro-hgvs/pull/621))
- *(benchmark)* time normalize tools by subprocess-internal timer ([#609](https://github.com/fulcrumgenomics/ferro-hgvs/pull/609)) ([#612](https://github.com/fulcrumgenomics/ferro-hgvs/pull/612))
- *(validate)* reject del/dup/delins numeric length not matching the position span (#486 length mismatch) ([#565](https://github.com/fulcrumgenomics/ferro-hgvs/pull/565))
- *(tool-support)* assign footnote markers in canonical order ([#602](https://github.com/fulcrumgenomics/ferro-hgvs/pull/602))
- *(reference)* report the true secondary-build transcript count ([#613](https://github.com/fulcrumgenomics/ferro-hgvs/pull/613))
- *(normalize)* size the intronic shuffle window to the enclosing intron ([#573](https://github.com/fulcrumgenomics/ferro-hgvs/pull/573)) ([#575](https://github.com/fulcrumgenomics/ferro-hgvs/pull/575))
- make HashMap-ordered output deterministic (VCF INFO + benchmark reports) ([#594](https://github.com/fulcrumgenomics/ferro-hgvs/pull/594))
- *(parser)* bring fast path to full parity with the generic parser (ncRNA, U base, positions, panic) ([#560](https://github.com/fulcrumgenomics/ferro-hgvs/pull/560))
- *(reference)* make cdot base→version fallback deterministic ([#583](https://github.com/fulcrumgenomics/ferro-hgvs/pull/583))
- *(conformance)* map ECOORDINATESYSTEMMISMATCH to ferro's Parse error ([#486](https://github.com/fulcrumgenomics/ferro-hgvs/pull/486)) ([#579](https://github.com/fulcrumgenomics/ferro-hgvs/pull/579))
- *(tool-support)* address CodeRabbit findings on the matrix generator ([#601](https://github.com/fulcrumgenomics/ferro-hgvs/pull/601))
- *(bench)* make the intronic projection fixture self-consistent ([#600](https://github.com/fulcrumgenomics/ferro-hgvs/pull/600))
- *(project)* emit nonsense (not fsTer1) for an immediate-stop frameshift ([#589](https://github.com/fulcrumgenomics/ferro-hgvs/pull/589))
- *(project)* scan the 3'UTR for a frameshift's new stop (fsTer{K}, not fsTer?) ([#587](https://github.com/fulcrumgenomics/ferro-hgvs/pull/587))
- *(build)* repair test code that broke under non-`dev` feature flags ([#580](https://github.com/fulcrumgenomics/ferro-hgvs/pull/580))
- *(validate)* reject delins/del/dup/repeat reference-sequence mismatches in strict mode ([#486](https://github.com/fulcrumgenomics/ferro-hgvs/pull/486)) ([#556](https://github.com/fulcrumgenomics/ferro-hgvs/pull/556))
- *(validate)* reject the RNA base U in a DNA-context edit (#486 ENODNA) ([#542](https://github.com/fulcrumgenomics/ferro-hgvs/pull/542))
- *(normalize)* coerce g. on a mitochondrial accession to m. (#487 mito group) ([#541](https://github.com/fulcrumgenomics/ferro-hgvs/pull/541))
- *(validate)* reject c./g./m. on a non-coding RNA reference (#486 coord-system mismatch) ([#543](https://github.com/fulcrumgenomics/ferro-hgvs/pull/543))
- *(normalize)* number coding-transcript r. positions CDS-relative ([#469](https://github.com/fulcrumgenomics/ferro-hgvs/pull/469)) ([#539](https://github.com/fulcrumgenomics/ferro-hgvs/pull/539))
- *(normalize)* spec-correct refusal of NG/LRG-parent c. pter/qter flank coordinates (#488 Phase 2b) ([#534](https://github.com/fulcrumgenomics/ferro-hgvs/pull/534))
- *(project)* mid-codon in-frame insertions render as delins, not a clean ins ([#511](https://github.com/fulcrumgenomics/ferro-hgvs/pull/511)) ([#517](https://github.com/fulcrumgenomics/ferro-hgvs/pull/517))
- *(normalize)* guard special/offset genome positions to avoid pterdel overflow panic ([#488](https://github.com/fulcrumgenomics/ferro-hgvs/pull/488)) ([#518](https://github.com/fulcrumgenomics/ferro-hgvs/pull/518))
- *(reference)* decline protein prediction for non-version-exact transcripts ([#505](https://github.com/fulcrumgenomics/ferro-hgvs/pull/505)) ([#519](https://github.com/fulcrumgenomics/ferro-hgvs/pull/519))
- *(project)* extend stop-loss substitution prediction with extTer notation ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#514](https://github.com/fulcrumgenomics/ferro-hgvs/pull/514))
- *(project)* predict p.(Met1?) for CDS-boundary edits reaching the initiation codon ([#504](https://github.com/fulcrumgenomics/ferro-hgvs/pull/504)) ([#513](https://github.com/fulcrumgenomics/ferro-hgvs/pull/513))
- *(hgvs)* distinguish frameshift fsTer? from short fs; predict fsTer? for no-stop frameshifts ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#515](https://github.com/fulcrumgenomics/ferro-hgvs/pull/515))
- *(project)* report p.(Met1?) for initiation-codon variants ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#512](https://github.com/fulcrumgenomics/ferro-hgvs/pull/512))
- *(project)* never emit descending positions in in-frame protein delins ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#510](https://github.com/fulcrumgenomics/ferro-hgvs/pull/510))
- *(normalize)* avoid subtract-overflow panic in minus-strand intronic shuffle ([#497](https://github.com/fulcrumgenomics/ferro-hgvs/pull/497))
- *(test/mutalyzer-normalize)* spec-correct bare-NP protein refs (7 cases) ([#495](https://github.com/fulcrumgenomics/ferro-hgvs/pull/495))
- *(reference)* flag cdot base-synthesis divergence risk (closes #471) ([#482](https://github.com/fulcrumgenomics/ferro-hgvs/pull/482))
- *(error-handling)* make detect_swapped_positions axis-aware for m./o. wraparound (closes #467) ([#474](https://github.com/fulcrumgenomics/ferro-hgvs/pull/474))
- *(cdot)* plumb transcript protein accession from cdot JSON ([#483](https://github.com/fulcrumgenomics/ferro-hgvs/pull/483))
- *(test/mutalyzer-normalize)* accept infos codes ferro does not model (22→11) ([#481](https://github.com/fulcrumgenomics/ferro-hgvs/pull/481))
- *(parser)* depth-aware separator scan in expanded allele entrypoints ([#476](https://github.com/fulcrumgenomics/ferro-hgvs/pull/476))
- *(error-handling)* extend W3009/W3010 corrector to uncertain-stop forms ([#461](https://github.com/fulcrumgenomics/ferro-hgvs/pull/461))
- *(parser)* reject single-position insertion per insertion.md:95-101 ([#446](https://github.com/fulcrumgenomics/ferro-hgvs/pull/446)) ([#450](https://github.com/fulcrumgenomics/ferro-hgvs/pull/450))
- *(edit)* drop brackets on single-payload InsertedSequence::Complex ([#443](https://github.com/fulcrumgenomics/ferro-hgvs/pull/443))
- *(service)* extend frameshift classifier span_len to Deletion + Duplication (closes #427) ([#435](https://github.com/fulcrumgenomics/ferro-hgvs/pull/435))
- *(parser)* reject `dupins<seq>` per duplication.md:92 ([#445](https://github.com/fulcrumgenomics/ferro-hgvs/pull/445)) ([#451](https://github.com/fulcrumgenomics/ferro-hgvs/pull/451))
- *(error-handling)* provider-aware W3016 for mixed-shape intronic-endpoint ranges (closes #429) ([#436](https://github.com/fulcrumgenomics/ferro-hgvs/pull/436))
- *(error-handling)* W3016 detects numeric-length-suffix disagreement on del/dup/inv (closes #439) ([#441](https://github.com/fulcrumgenomics/ferro-hgvs/pull/441))
- *(service)* inversion is in-frame at DNA level (closes #438) ([#440](https://github.com/fulcrumgenomics/ferro-hgvs/pull/440))
- *(validate)* tighten validate_multirepeat_tract for nonzero-min middle (closes #428) ([#434](https://github.com/fulcrumgenomics/ferro-hgvs/pull/434))
- *(service,validation,spdi)* three delins consistency bugs (closes #394) ([#419](https://github.com/fulcrumgenomics/ferro-hgvs/pull/419))
- *(normalize)* canon boundary-spanning ins + side-aware ins→dup anchor (closes #402) ([#413](https://github.com/fulcrumgenomics/ferro-hgvs/pull/413))
- *(normalize)* close untracked biocommons-burndown divergences (closes #418) ([#420](https://github.com/fulcrumgenomics/ferro-hgvs/pull/420))
- *(normalize)* direction-aware tie-break in insertion_to_duplication (closes #403) ([#408](https://github.com/fulcrumgenomics/ferro-hgvs/pull/408))
- fix(spdi)+(error-handling)+(test): six HGVS↔SPDI / W3016 / VCF audit gaps (closes #390) ([#406](https://github.com/fulcrumgenomics/ferro-hgvs/pull/406))
- *(project)* build-aware cdot lookups + widen fan-out path for c./n./r. (closes #389) ([#398](https://github.com/fulcrumgenomics/ferro-hgvs/pull/398))
- *(normalize)* extend W4004 PositionPastEnd to intronic offsets (closes #392) ([#414](https://github.com/fulcrumgenomics/ferro-hgvs/pull/414))
- *(normalize)* preserve spanning duplications in CDS-start canon clamp (closes #401) ([#410](https://github.com/fulcrumgenomics/ferro-hgvs/pull/410))
- *(normalize)* clamp 3'-canonicalization at CDS-end for c.-axis inputs (closes #387) ([#388](https://github.com/fulcrumgenomics/ferro-hgvs/pull/388))
- *(normalize)* reject c. positions past CDS-end / transcript-end (W4004 PositionPastEnd, closes #336) ([#342](https://github.com/fulcrumgenomics/ferro-hgvs/pull/342))
- *(normalize)* clamp 5'-canonicalization at CDS-start for c.-axis inputs (closes #383) ([#385](https://github.com/fulcrumgenomics/ferro-hgvs/pull/385))
- *(normalize)* apply exon-junction exception to n. and r. axes ([#374](https://github.com/fulcrumgenomics/ferro-hgvs/pull/374))
- *(parser)* propagate structured chunk-level diagnostics in slash-form alleles (closes #375) ([#376](https://github.com/fulcrumgenomics/ferro-hgvs/pull/376))
- *(parser)* extend E3006 self-cancelling detection to *-region and intronic ranges (closes #172) ([#371](https://github.com/fulcrumgenomics/ferro-hgvs/pull/371))
- *(diagnostics)* preserve source span on E3006 SelfCancellingAllele (closes #171) ([#369](https://github.com/fulcrumgenomics/ferro-hgvs/pull/369))
- *(preprocessor)* accept lowercase IUPAC for r. multi-base subs (closes #170) ([#367](https://github.com/fulcrumgenomics/ferro-hgvs/pull/367))
- *(normalize)* intronic + boundary normalization with NG/NC-parent inputs (closes #332) ([#364](https://github.com/fulcrumgenomics/ferro-hgvs/pull/364))
- *(parser,error-handling)* wire W3019 NonSpecMosaicForm for nested + [a/b] mosaic (closes #281) ([#298](https://github.com/fulcrumgenomics/ferro-hgvs/pull/298))
- *(normalize)* clamp 3'-rule shuffle at CDS↔UTR axis boundary (closes #337) ([#343](https://github.com/fulcrumgenomics/ferro-hgvs/pull/343))
- *(normalize)* strip explicit deleted-sequence from delins (closes #338) ([#344](https://github.com/fulcrumgenomics/ferro-hgvs/pull/344))
- *(parser,error-handling)* mito heteroplasmy prose diagnostic + W3017/W3018 SVAs (closes #278) ([#295](https://github.com/fulcrumgenomics/ferro-hgvs/pull/295))
- *(normalize)* bail gracefully when ref window mismatches HGVS span (closes #339) ([#345](https://github.com/fulcrumgenomics/ferro-hgvs/pull/345))
- *(parser,display)* preserve parens on inner repeat-count range insN[(150_180)] (closes #285) ([#306](https://github.com/fulcrumgenomics/ferro-hgvs/pull/306))
- *(normalize)* set RefSeqMismatch.corrected to reflect actual correction (closes #280) ([#297](https://github.com/fulcrumgenomics/ferro-hgvs/pull/297))
- *(normalize)* direction-aware ins→dup for homopolymers (closes #340 subgroup (i)) ([#346](https://github.com/fulcrumgenomics/ferro-hgvs/pull/346))
- *(reference)* synthesize transcript bases via cdot exon-alignment when FASTA version is missing (closes #331) ([#341](https://github.com/fulcrumgenomics/ferro-hgvs/pull/341))
- *(parser)* accept p.(0) whole-entity predicted no-protein (closes #289) ([#307](https://github.com/fulcrumgenomics/ferro-hgvs/pull/307))
- *(parser)* chain-expand mito = arms to inherit position, not synthesize m.1= (closes #284) ([#302](https://github.com/fulcrumgenomics/ferro-hgvs/pull/302))
- *(parser)* route [0] inside protein compact trans-allele to NoProtein (closes #277) ([#294](https://github.com/fulcrumgenomics/ferro-hgvs/pull/294))
- *(display)* canonicalize r. emission T→u in repeat units and bases (closes #276) ([#293](https://github.com/fulcrumgenomics/ferro-hgvs/pull/293))
- *(normalize)* extend codon-frame exception to r., 3+ chains, and sub+del (closes #275) ([#292](https://github.com/fulcrumgenomics/ferro-hgvs/pull/292))
- *(mock)* disambiguate genomic vs transcript lookups when ids collide (closes #311) ([#312](https://github.com/fulcrumgenomics/ferro-hgvs/pull/312))
- *(parser,error-handling)* wire W3021 ProteinBracketedAaInsertion (closes #290) ([#305](https://github.com/fulcrumgenomics/ferro-hgvs/pull/305))
- *(fasta)* route known-contig lookups in FastaProvider through the FASTA path (closes #315) ([#318](https://github.com/fulcrumgenomics/ferro-hgvs/pull/318))
- *(fasta)* require version-boundary equality in MmapFastaProvider transcript lookup (closes #314) ([#317](https://github.com/fulcrumgenomics/ferro-hgvs/pull/317))
- *(projector)* protein prediction for non-RefSeq transcripts + drop (GENE) selector from p. Display ([#313](https://github.com/fulcrumgenomics/ferro-hgvs/pull/313))
- *(validation)* wire W3016 LengthMismatch soft-validation warning (closes #81 K1 follow-up) ([#272](https://github.com/fulcrumgenomics/ferro-hgvs/pull/272))
- *(error-handling)* extend W4001 SwappedPositions to offset and *N markers (closes #81 L2 follow-up to #264) ([#271](https://github.com/fulcrumgenomics/ferro-hgvs/pull/271))
- *(spdi)* preserve inv and m. coord system on SPDI -> HGVS (closes #81 K1 follow-up) ([#274](https://github.com/fulcrumgenomics/ferro-hgvs/pull/274))
- *(error-handling)* wire W1004 MixedCaseEditType warning into preprocessor (closes #81 L2 follow-up) ([#273](https://github.com/fulcrumgenomics/ferro-hgvs/pull/273))
- *(error-handling)* wire W4001 SwappedPositions warning into preprocessor (closes #81 L2 remaining) ([#264](https://github.com/fulcrumgenomics/ferro-hgvs/pull/264))
- *(parser,display)* predicted-edit wrapper `c.(<pos><edit>)` on every NA coord (closes #241) ([#242](https://github.com/fulcrumgenomics/ferro-hgvs/pull/242))
- *(parser,display)* single uncertain position range c.(a_b)<edit> (closes #237) ([#238](https://github.com/fulcrumgenomics/ferro-hgvs/pull/238))
- *(protein)* extension Display emits extTerN to match Frameshift Ter canonicalization (closes #224) ([#225](https://github.com/fulcrumgenomics/ferro-hgvs/pull/225))
- *(parser)* accept `?<edit>` on g./r./n. and `?A>G` / `?=` (closes #239) ([#240](https://github.com/fulcrumgenomics/ferro-hgvs/pull/240))
- *(normalize)* read dup bases from reference (single-pass idempotency for dup with mismatched stated-ref) (closes #219) ([#230](https://github.com/fulcrumgenomics/ferro-hgvs/pull/230))
- *(vcf)* classify protein consequences by enum variant, not substring (closes #228) ([#229](https://github.com/fulcrumgenomics/ferro-hgvs/pull/229))
- *(normalize)* wire window-based normalization for m. variants ([#210](https://github.com/fulcrumgenomics/ferro-hgvs/pull/210)) ([#213](https://github.com/fulcrumgenomics/ferro-hgvs/pull/213))
- *(normalize)* exempt introns and UTR from codon-frame gate ([#209](https://github.com/fulcrumgenomics/ferro-hgvs/pull/209)) ([#211](https://github.com/fulcrumgenomics/ferro-hgvs/pull/211))
- *(parser)* depth-aware close-bracket lookup in nested HGVS brackets (closes #207) ([#208](https://github.com/fulcrumgenomics/ferro-hgvs/pull/208))
- *(normalize)* decompose delins to multi-sub canonical form when sub > delins applies (closes #165) ([#206](https://github.com/fulcrumgenomics/ferro-hgvs/pull/206))

### Other

- *(benchmark)* refresh measured performance tables for 0.7.0 ([#1006](https://github.com/fulcrumgenomics/ferro-hgvs/pull/1006))
- *(spec-fixture)* classify the remaining spec-audit parse gaps ([#993](https://github.com/fulcrumgenomics/ferro-hgvs/pull/993))
- *(spec-fixture)* reassemble spot-split variants, drop harvest noise ([#990](https://github.com/fulcrumgenomics/ferro-hgvs/pull/990))
- *(conformance)* regenerate stale mutalyzer failure-patterns summary ([#998](https://github.com/fulcrumgenomics/ferro-hgvs/pull/998))
- *(conformance)* reclassify NM_017668.3 c.CAA[5] repeat as accepted_divergence, not known_bug ([#920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/920)) ([#986](https://github.com/fulcrumgenomics/ferro-hgvs/pull/986))
- *(conformance)* reclassify c.235_237delinsTAT as spec_citation, not known_bug ([#920](https://github.com/fulcrumgenomics/ferro-hgvs/pull/920)) ([#984](https://github.com/fulcrumgenomics/ferro-hgvs/pull/984))
- *(conformance)* re-bless mutalyzer baseline-failures — all ledgers empty (closes #326) ([#982](https://github.com/fulcrumgenomics/ferro-hgvs/pull/982))
- *(reference)* defer Ensembl cdot load so RefSeq startup skips the merge ([#965](https://github.com/fulcrumgenomics/ferro-hgvs/pull/965))
- *(conformance)* document the dormant proj-near-disc + NR2E3 rows ([#946](https://github.com/fulcrumgenomics/ferro-hgvs/pull/946)) ([#949](https://github.com/fulcrumgenomics/ferro-hgvs/pull/949))
- *(conformance)* make the whole-CDS-deletion p.0?/p.(Met1?) divergence terminal ([#945](https://github.com/fulcrumgenomics/ferro-hgvs/pull/945)) ([#948](https://github.com/fulcrumgenomics/ferro-hgvs/pull/948))
- *(conformance)* reclassify the 8 NG_007485.1 noncoding rows as spec_citation ([#921](https://github.com/fulcrumgenomics/ferro-hgvs/pull/921)) ([#934](https://github.com/fulcrumgenomics/ferro-hgvs/pull/934))
- *(normalize)* correct misleading circular o. normalization comments ([#957](https://github.com/fulcrumgenomics/ferro-hgvs/pull/957))
- *(benchmark)* refresh comparator pins and record verified provenance ([#890](https://github.com/fulcrumgenomics/ferro-hgvs/pull/890)) ([#941](https://github.com/fulcrumgenomics/ferro-hgvs/pull/941))
- *(conformance)* reclassify copy-range delins rows as spec_citation ([#922](https://github.com/fulcrumgenomics/ferro-hgvs/pull/922)) ([#937](https://github.com/fulcrumgenomics/ferro-hgvs/pull/937))
- *(conformance)* re-home stale-tracker annotations to live issues ([#912](https://github.com/fulcrumgenomics/ferro-hgvs/pull/912)) ([#924](https://github.com/fulcrumgenomics/ferro-hgvs/pull/924))
- *(conformance)* disposition residual infos rows; reclassify NM_017668.3 exon-junction rows as spec_citation (#908, #918) ([#919](https://github.com/fulcrumgenomics/ferro-hgvs/pull/919))
- *(conformance)* disposition the 3 NG_008939 stop-insertion coding_protein rows (#911/#913 follow-up) ([#927](https://github.com/fulcrumgenomics/ferro-hgvs/pull/927))
- *(conformance)* reconcile the redundant cis-= allele axis inconsistency (#912 action 4) ([#928](https://github.com/fulcrumgenomics/ferro-hgvs/pull/928))
- *(conformance)* accept hgvs-rs-projection form-currency divergences ([#925](https://github.com/fulcrumgenomics/ferro-hgvs/pull/925)) ([#930](https://github.com/fulcrumgenomics/ferro-hgvs/pull/930))
- *(tool-support)* label the support matrix for the ferro 0.7.0 release ([#914](https://github.com/fulcrumgenomics/ferro-hgvs/pull/914))
- *(conformance)* harden the reference-identity signature ([#905](https://github.com/fulcrumgenomics/ferro-hgvs/pull/905)) ([#910](https://github.com/fulcrumgenomics/ferro-hgvs/pull/910))
- *(conformance)* accepted-divergence annotation sweep for the mutalyzer parity corpus ([#861](https://github.com/fulcrumgenomics/ferro-hgvs/pull/861)) ([#909](https://github.com/fulcrumgenomics/ferro-hgvs/pull/909))
- *(conformance)* correct insCAT protein_description; demote the 5'UTR-insertion divergence ([#891](https://github.com/fulcrumgenomics/ferro-hgvs/pull/891))
- *(benchmark)* record comparator provenance in the mutalyzer-normalize corpus header ([#890](https://github.com/fulcrumgenomics/ferro-hgvs/pull/890))
- *(conformance)* let accepted_rejection disposition an empty-projection Err ([#903](https://github.com/fulcrumgenomics/ferro-hgvs/pull/903))
- *(conformance)* demote stale known_bug #487 for LRG_303t1 out-of-phase insertion ([#906](https://github.com/fulcrumgenomics/ferro-hgvs/pull/906))
- *(project)* untangle project_to_genomic_nc — extract the coordinate map and #785 gate ([#868](https://github.com/fulcrumgenomics/ferro-hgvs/pull/868)) ([#902](https://github.com/fulcrumgenomics/ferro-hgvs/pull/902))
- *(conformance)* bless genomic-axis rejections via multi-axis accepted_rejection (#870 follow-up)
- *(conformance)* route genomic axis through the user-facing project_variant path ([#870](https://github.com/fulcrumgenomics/ferro-hgvs/pull/870))
- *(normalize)* make normalize_repeat the single source of truth for ins→repeat canonicalization ([#866](https://github.com/fulcrumgenomics/ferro-hgvs/pull/866)) ([#897](https://github.com/fulcrumgenomics/ferro-hgvs/pull/897))
- remove machine-local paths from committed example and runbook ([#889](https://github.com/fulcrumgenomics/ferro-hgvs/pull/889))
- *(conformance)* converge NG_-annotated transcript version on RefSeqGene selectors ([#859](https://github.com/fulcrumgenomics/ferro-hgvs/pull/859)) ([#871](https://github.com/fulcrumgenomics/ferro-hgvs/pull/871))
- *(conformance)* accept ferro's decline of NM_017668.3 coords outside NG_009299.1 coverage (#853, #865) ([#877](https://github.com/fulcrumgenomics/ferro-hgvs/pull/877))
- *(nightly)* derive ng_hosted_transcripts + NG_ placements in the nightly reference ([#875](https://github.com/fulcrumgenomics/ferro-hgvs/pull/875)) ([#876](https://github.com/fulcrumgenomics/ferro-hgvs/pull/876))
- *(conformance)* cite ferro's spec-correct inserted-inversion resolution ([#854](https://github.com/fulcrumgenomics/ferro-hgvs/pull/854)) ([#878](https://github.com/fulcrumgenomics/ferro-hgvs/pull/878))
- *(conformance)* accept ferro's spec-correct refusal of non-standard/absent NG_ selectors ([#858](https://github.com/fulcrumgenomics/ferro-hgvs/pull/858)) ([#873](https://github.com/fulcrumgenomics/ferro-hgvs/pull/873))
- *(conformance)* rebless rna_description ledger and protein empty-projection pin ([#799](https://github.com/fulcrumgenomics/ferro-hgvs/pull/799)) ([#846](https://github.com/fulcrumgenomics/ferro-hgvs/pull/846))
- *(normalize)* resolve the cigar-aware normalization backlog ([#811](https://github.com/fulcrumgenomics/ferro-hgvs/pull/811)) ([#833](https://github.com/fulcrumgenomics/ferro-hgvs/pull/833))
- *(conformance)* allow multiple per-axis spec_citation annotations on a Case ([#830](https://github.com/fulcrumgenomics/ferro-hgvs/pull/830))
- *(conformance)* disposition bare-NP protein-framing rows as spec_overridden ([#826](https://github.com/fulcrumgenomics/ferro-hgvs/pull/826))
- *(spdi)* run SPDI roundtrip + rsID tests offline via synthetic fixture ([#825](https://github.com/fulcrumgenomics/ferro-hgvs/pull/825))
- *(python_helpers)* remove dead edit_type_from_debug {:?}-string matcher ([#829](https://github.com/fulcrumgenomics/ferro-hgvs/pull/829))
- *(conformance)* cover exon-junction rows in the genomic-axis gate ([#751](https://github.com/fulcrumgenomics/ferro-hgvs/pull/751)) ([#789](https://github.com/fulcrumgenomics/ferro-hgvs/pull/789))
- *(reference)* derive accession→assembly from NCBI assembly_report ([#716](https://github.com/fulcrumgenomics/ferro-hgvs/pull/716)) ([#788](https://github.com/fulcrumgenomics/ferro-hgvs/pull/788))
- *(conformance)* decrement protein empty-projection budget after #778 ([#786](https://github.com/fulcrumgenomics/ferro-hgvs/pull/786))
- *(project)* predict VSIR reverse-strand whole-exon frameshift ([#774](https://github.com/fulcrumgenomics/ferro-hgvs/pull/774))
- *(conformance)* re-bless empty-projection budget on canonical reference ([#775](https://github.com/fulcrumgenomics/ferro-hgvs/pull/775))
- *(conformance)* disposition pter/qter-in-c flank payloads; green the normalized axis ([#770](https://github.com/fulcrumgenomics/ferro-hgvs/pull/770))
- *(conformance)* accept coding_protein transcript-set enumeration divergences ([#768](https://github.com/fulcrumgenomics/ferro-hgvs/pull/768))
- *(conformance)* pin the empty-projection budget to its reference ([#767](https://github.com/fulcrumgenomics/ferro-hgvs/pull/767))
- *(conformance)* restore clobbered genomic #745 + add accepted_rejection ([#654](https://github.com/fulcrumgenomics/ferro-hgvs/pull/654)) ([#757](https://github.com/fulcrumgenomics/ferro-hgvs/pull/757))
- *(conformance)* disposition the normalized-axis Ok-mismatch divergences ([#654](https://github.com/fulcrumgenomics/ferro-hgvs/pull/654)) ([#755](https://github.com/fulcrumgenomics/ferro-hgvs/pull/755))
- *(conformance)* forward dropped methods through the ArcProvider wrapper (#726 follow-up) ([#754](https://github.com/fulcrumgenomics/ferro-hgvs/pull/754))
- *(conformance)* green errors axis — bucket accepted divergences ([#486](https://github.com/fulcrumgenomics/ferro-hgvs/pull/486)) ([#752](https://github.com/fulcrumgenomics/ferro-hgvs/pull/752))
- *(conformance)* reclassify #745 homopolymer repeat-contraction as accepted divergence ([#750](https://github.com/fulcrumgenomics/ferro-hgvs/pull/750))
- *(reference)* NCBI assembly_report parser + data-driven build inference ([#716](https://github.com/fulcrumgenomics/ferro-hgvs/pull/716)) ([#735](https://github.com/fulcrumgenomics/ferro-hgvs/pull/735))
- consolidate ~230 integration-test binaries into one crate ([#725](https://github.com/fulcrumgenomics/ferro-hgvs/pull/725))
- *(annotate-vcf)* annotate VCF records in place instead of cloning ([#698](https://github.com/fulcrumgenomics/ferro-hgvs/pull/698))
- *(annotate-vcf)* parallelize file-input record annotation behind -j/--workers ([#696](https://github.com/fulcrumgenomics/ferro-hgvs/pull/696))
- *(parse)* parallelize the batch CLI loop behind -j/--workers ([#688](https://github.com/fulcrumgenomics/ferro-hgvs/pull/688))
- *(vcf)* remove dead batch processors from src/vcf/batch.rs ([#699](https://github.com/fulcrumgenomics/ferro-hgvs/pull/699))
- *(annotate)* interval-index TranscriptDb::get_by_region ([#697](https://github.com/fulcrumgenomics/ferro-hgvs/pull/697))
- *(annotate-vcf)* buffer the output writer ([#695](https://github.com/fulcrumgenomics/ferro-hgvs/pull/695))
- *(python)* parallelize the batch API and release the GIL ([#691](https://github.com/fulcrumgenomics/ferro-hgvs/pull/691))
- *(normalize)* parallelize the batch CLI loop behind -j/--workers ([#687](https://github.com/fulcrumgenomics/ferro-hgvs/pull/687))
- *(pre-commit)* bump ruff hook to v0.15.12 to match CI ([#705](https://github.com/fulcrumgenomics/ferro-hgvs/pull/705))
- *(normalize)* keep `IntronicVariant.variant` machine-readable ([#685](https://github.com/fulcrumgenomics/ferro-hgvs/pull/685))
- *(conformance)* add an output-quality signal to the mutalyzer axis harness ([#651](https://github.com/fulcrumgenomics/ferro-hgvs/pull/651)) ([#675](https://github.com/fulcrumgenomics/ferro-hgvs/pull/675))
- *(conformance)* classify remaining #487 3′-shift audit rows ([#673](https://github.com/fulcrumgenomics/ferro-hgvs/pull/673))
- *(parser)* skip Unicode trim when input has no surrounding whitespace ([#659](https://github.com/fulcrumgenomics/ferro-hgvs/pull/659))
- *(readme)* honest benchmark framing — offline best-case + single-threaded parse ([#660](https://github.com/fulcrumgenomics/ferro-hgvs/pull/660))
- *(normalize)* pre-size runtime cdot maps when rebuilding from the archive ([#666](https://github.com/fulcrumgenomics/ferro-hgvs/pull/666))
- *(normalize)* use FxHashMap for the runtime cdot reference maps ([#665](https://github.com/fulcrumgenomics/ferro-hgvs/pull/665))
- *(normalize)* avoid per-variant allocations and walks in input preprocessing ([#664](https://github.com/fulcrumgenomics/ferro-hgvs/pull/664))
- *(normalize)* skip reference teardown at process exit ([#663](https://github.com/fulcrumgenomics/ferro-hgvs/pull/663))
- *(parser)* fast-path identity edits (protein and nucleotide) ([#657](https://github.com/fulcrumgenomics/ferro-hgvs/pull/657))
- *(parser)* fast-path protein missense substitutions ([#650](https://github.com/fulcrumgenomics/ferro-hgvs/pull/650))
- *(parser)* fast-path intronic and UTR coding substitutions ([#643](https://github.com/fulcrumgenomics/ferro-hgvs/pull/643))
- *(benchmark)* refresh measured tables with startup-excluded normalize timing ([#648](https://github.com/fulcrumgenomics/ferro-hgvs/pull/648))
- *(parser)* fast-path plain deletions and duplications ([#641](https://github.com/fulcrumgenomics/ferro-hgvs/pull/641))
- *(conformance)* annotate coding non-codon-aligned repeats as spec-overrides ([#640](https://github.com/fulcrumgenomics/ferro-hgvs/pull/640))
- *(parser)* eliminate per-parse accession allocations (~27% faster parsing) ([#639](https://github.com/fulcrumgenomics/ferro-hgvs/pull/639))
- *(conformance)* annotate + hermetic-gate the biocommons normalize corpus (closes #325) ([#636](https://github.com/fulcrumgenomics/ferro-hgvs/pull/636))
- *(reference)* lazily load the GRCh37 secondary cdot until first use ([#635](https://github.com/fulcrumgenomics/ferro-hgvs/pull/635))
- remove committed machine-specific manifest path from benchmark tests ([#634](https://github.com/fulcrumgenomics/ferro-hgvs/pull/634))
- *(conformance)* disposition the residual signal — fully classify the hgvs-rs dashboard (burn-down PR 3) ([#620](https://github.com/fulcrumgenomics/ferro-hgvs/pull/620))
- add CodeRabbit config (assertive profile + path instructions) ([#627](https://github.com/fulcrumgenomics/ferro-hgvs/pull/627))
- *(conformance)* quarantine hgvs-rs data-source divergences so red = real (burn-down PR 2) ([#614](https://github.com/fulcrumgenomics/ferro-hgvs/pull/614))
- *(conformance)* version-insensitive comparator + selection-coverage for the hgvs-rs corpus ([#607](https://github.com/fulcrumgenomics/ferro-hgvs/pull/607))
- fix broken DOI badge in README ([#574](https://github.com/fulcrumgenomics/ferro-hgvs/pull/574))
- *(conformance)* accept compound-allele SHUFFLE_APPLIED as a tracked divergence ([#499](https://github.com/fulcrumgenomics/ferro-hgvs/pull/499)) ([#572](https://github.com/fulcrumgenomics/ferro-hgvs/pull/572))
- *(conformance)* reject position-past-end on the errors axis (#486 EOUTOFBOUNDARY) ([#567](https://github.com/fulcrumgenomics/ferro-hgvs/pull/567))
- *(reference)* precompute has_genomic_data instead of scanning the index ([#610](https://github.com/fulcrumgenomics/ferro-hgvs/pull/610))
- *(cli)* format the normalized variant once per output line ([#611](https://github.com/fulcrumgenomics/ferro-hgvs/pull/611))
- *(conformance)* reject ISCN-only / superseded `::` rearrangement forms ([#546](https://github.com/fulcrumgenomics/ferro-hgvs/pull/546)) ([#608](https://github.com/fulcrumgenomics/ferro-hgvs/pull/608))
- perf-regression gate (cdot load-source structural test + nightly startup budget) ([#606](https://github.com/fulcrumgenomics/ferro-hgvs/pull/606))
- *(normalize)* resolve intronic transcripts by accession, not a cloned variant ([#582](https://github.com/fulcrumgenomics/ferro-hgvs/pull/582))
- *(reference)* reuse FASTA file handles with positioned reads (~1.9x normalize) ([#570](https://github.com/fulcrumgenomics/ferro-hgvs/pull/570))
- *(parser)* route the default parse_hgvs through the fast path (~1.7x on ClinVar) ([#562](https://github.com/fulcrumgenomics/ferro-hgvs/pull/562))
- *(reference)* memoize resolved transcripts in MultiFastaProvider ([#578](https://github.com/fulcrumgenomics/ferro-hgvs/pull/578))
- *(reference)* faster, deterministic FASTA index build (FxHashMap + highest-version fallback) ([#588](https://github.com/fulcrumgenomics/ferro-hgvs/pull/588))
- *(project)* add c→p, c↔n, ref-AA cache, and batch projection microbenchmarks (Cycle 1c-ii) ([#599](https://github.com/fulcrumgenomics/ferro-hgvs/pull/599))
- *(conformance)* import the hgvs-rs projection corpus as a third oracle (Cycle 1c-i) ([#598](https://github.com/fulcrumgenomics/ferro-hgvs/pull/598))
- *(reference)* binary cache for the GFF/GTF TranscriptDb (annotate-vcf ~3s → ~0.6s) ([#593](https://github.com/fulcrumgenomics/ferro-hgvs/pull/593))
- *(reference)* rkyv archive for the cdot cache (zero-copy load, ~0.37s off startup) ([#591](https://github.com/fulcrumgenomics/ferro-hgvs/pull/591))
- *(reference)* version + self-heal the cdot bincode cache (fix silent JSON fallback) ([#585](https://github.com/fulcrumgenomics/ferro-hgvs/pull/585))
- *(tests)* stop committing the generated spec-normalization fixture ([#586](https://github.com/fulcrumgenomics/ferro-hgvs/pull/586))
- *(normalize)* two no-functional-change normalize hot-path wins (shuffle-info + protein-length) ([#553](https://github.com/fulcrumgenomics/ferro-hgvs/pull/553))
- *(bench)* broaden microbenchmark suite + tag-sweep runner (perf-effort foundation) ([#554](https://github.com/fulcrumgenomics/ferro-hgvs/pull/554))
- fix nightly reference-aware tests (`ferro prepare --output-dir`) + add nightly badge ([#551](https://github.com/fulcrumgenomics/ferro-hgvs/pull/551))
- *(biocommons)* mirror the improvement conformance disposition ([#503](https://github.com/fulcrumgenomics/ferro-hgvs/pull/503)) ([#540](https://github.com/fulcrumgenomics/ferro-hgvs/pull/540))
- *(mutalyzer)* demote c.41A>C — substitution is a missense, not a frameshift ([#498](https://github.com/fulcrumgenomics/ferro-hgvs/pull/498)) ([#516](https://github.com/fulcrumgenomics/ferro-hgvs/pull/516))
- *(mutalyzer)* add improvement disposition; demote NG_ gene-symbol selector divergences ([#500](https://github.com/fulcrumgenomics/ferro-hgvs/pull/500)) ([#501](https://github.com/fulcrumgenomics/ferro-hgvs/pull/501))
- *(mutalyzer)* known_bug disposition + XPASS detection ([#478](https://github.com/fulcrumgenomics/ferro-hgvs/pull/478)) ([#489](https://github.com/fulcrumgenomics/ferro-hgvs/pull/489))
- *(biocommons)* annotation + XPASS conformance harness (#478 pillars 1-2) ([#484](https://github.com/fulcrumgenomics/ferro-hgvs/pull/484))
- *(mutalyzer-normalize)* demote ins[...] cases fixed by #333 (49 rows) ([#479](https://github.com/fulcrumgenomics/ferro-hgvs/pull/479))
- *(test)* prune stale NM_000051.3:c.1_2insCA from biocommons baseline-failures ([#477](https://github.com/fulcrumgenomics/ferro-hgvs/pull/477))
- trans-allele dedup, mixed phase reject, RNA whole-entity in brackets (closes #396) ([#425](https://github.com/fulcrumgenomics/ferro-hgvs/pull/425))
- *(protein)* pin stop-codon glyph canonicalization ([#453](https://github.com/fulcrumgenomics/ferro-hgvs/pull/453)) ([#454](https://github.com/fulcrumgenomics/ferro-hgvs/pull/454))
- *(protein)* pin construct-boundary canonicalization probes (F13) ([#458](https://github.com/fulcrumgenomics/ferro-hgvs/pull/458))
- *(allele)* pin allele-bracket grammar corner cases (F12) ([#457](https://github.com/fulcrumgenomics/ferro-hgvs/pull/457))
- *(parser)* pin strict-reject mode for W3011 DelSizeSuffix ([#447](https://github.com/fulcrumgenomics/ferro-hgvs/pull/447)) ([#459](https://github.com/fulcrumgenomics/ferro-hgvs/pull/459))
- *(parser)* fix red main — drop obsolete dupins-accept probe (contradicts #445) ([#463](https://github.com/fulcrumgenomics/ferro-hgvs/pull/463))
- *(parser)* pin input-hygiene & spec-rejected parser behavior (F1+F2) ([#444](https://github.com/fulcrumgenomics/ferro-hgvs/pull/444))
- *(spec)* pin miscellaneous coverage corners (F6+F7+F10-11+F14-18) ([#456](https://github.com/fulcrumgenomics/ferro-hgvs/pull/456))
- omnibus follow-up: loader + diagnostics + ci polish (closes #397) ([#426](https://github.com/fulcrumgenomics/ferro-hgvs/pull/426))
- *(parser)* pin strict-whitespace parse mode ([#449](https://github.com/fulcrumgenomics/ferro-hgvs/pull/449)) ([#455](https://github.com/fulcrumgenomics/ferro-hgvs/pull/455))
- *(rna)* add r.↔c. consistency + RNA splicing marker coverage ([#442](https://github.com/fulcrumgenomics/ferro-hgvs/pull/442))
- *(parser)* drop dead 0? clause in protein trans-allele shorthand (closes #424) ([#431](https://github.com/fulcrumgenomics/ferro-hgvs/pull/431))
- *(benchmark)* wire supplemental_fasta into manifest after patterns fetch ([#417](https://github.com/fulcrumgenomics/ferro-hgvs/pull/417))
- *(test)* resync biocommons-normalize baseline-failures with manifest-mode xfail ([#405](https://github.com/fulcrumgenomics/ferro-hgvs/pull/405))
- *(normalize)* 3'-rule completeness — repeat-depth audit + #343 deferred symmetry tests (closes #391) ([#409](https://github.com/fulcrumgenomics/ferro-hgvs/pull/409))
- *(test)* sync biocommons-normalize baseline-failures with current xfail ([#386](https://github.com/fulcrumgenomics/ferro-hgvs/pull/386))
- *(normalize)* pin delins → dup canonicalization over downstream tract (closes #382) ([#384](https://github.com/fulcrumgenomics/ferro-hgvs/pull/384))
- *(python)* migrate PyO3 bindings to abi3 stable ABI ([#381](https://github.com/fulcrumgenomics/ferro-hgvs/pull/381))
- *(project)* rebuild VariantProjector hot paths ([#361](https://github.com/fulcrumgenomics/ferro-hgvs/pull/361))
- *(error-handling)* retire W4002 PositionZero duplicate identity (closes #269) ([#370](https://github.com/fulcrumgenomics/ferro-hgvs/pull/370))
- *(parser,display)* pin c./r. predicted-edit wrapper Display parity (closes #300) ([#372](https://github.com/fulcrumgenomics/ferro-hgvs/pull/372))
- *(biocommons)* import normalize test corpus (88 cases, stacked on #323) ([#324](https://github.com/fulcrumgenomics/ferro-hgvs/pull/324))
- update stale Display test expectations to HGVS spec-canonical forms ([#366](https://github.com/fulcrumgenomics/ferro-hgvs/pull/366))
- *(contributing)* defer CHANGELOG.md to release-plz ([#360](https://github.com/fulcrumgenomics/ferro-hgvs/pull/360))
- *(reference)* impl Default for Transcript to enable spread-update fixtures ([#322](https://github.com/fulcrumgenomics/ferro-hgvs/pull/322))
- *(mutalyzer)* import normalize test corpus (320 cases, 8 axes) ([#323](https://github.com/fulcrumgenomics/ferro-hgvs/pull/323))
- *(normalize)* pin r.-axis (transcript-1-relative) convention (closes #291) ([#304](https://github.com/fulcrumgenomics/ferro-hgvs/pull/304))
- *(convert)* pin c. → r. reference-driven conversion across edit shapes (closes #283) ([#301](https://github.com/fulcrumgenomics/ferro-hgvs/pull/301))
- *(audit)* HGVS ↔ SPDI round-trip coverage matrix (closes #81 K1 remaining) ([#260](https://github.com/fulcrumgenomics/ferro-hgvs/pull/260))
- *(audit)* HGVS ↔ VCF conversion surface (closes #81 K2) ([#262](https://github.com/fulcrumgenomics/ferro-hgvs/pull/262))
- *(audit)* minus-strand intronic position ordering at parse/Display (closes #81 J3) ([#258](https://github.com/fulcrumgenomics/ferro-hgvs/pull/258))
- *(audit)* n. upstream/downstream * and - markers canonical output (closes #81 J2) ([#256](https://github.com/fulcrumgenomics/ferro-hgvs/pull/256))
- *(audit)* LRG references full round-trip across edit types (closes #81 I1) ([#248](https://github.com/fulcrumgenomics/ferro-hgvs/pull/248))
- *(audit)* boundary-spanning del/delins ranges across UTR/CDS (closes #81 J1) ([#254](https://github.com/fulcrumgenomics/ferro-hgvs/pull/254))
- *(audit)* chromosome alias input + canonical output (closes #81 I4) ([#252](https://github.com/fulcrumgenomics/ferro-hgvs/pull/252))
- *(audit)* versioned vs unversioned accession policy (closes #81 I2) ([#250](https://github.com/fulcrumgenomics/ferro-hgvs/pull/250))
- *(mito)* extended heteroplasmy / prose-shape audit and policy doc (closes #235) ([#236](https://github.com/fulcrumgenomics/ferro-hgvs/pull/236))
- *(rna)* pin r./c. canonicalization parity across edit shapes (closes #233) ([#234](https://github.com/fulcrumgenomics/ferro-hgvs/pull/234))
- *(rna)* dedicated lowercase enforcement coverage across edit types (closes #231) ([#232](https://github.com/fulcrumgenomics/ferro-hgvs/pull/232))
- protein canonicalization audit umbrella for #81 § D — D1 / D3 / D4 / D8 (closes #226) ([#227](https://github.com/fulcrumgenomics/ferro-hgvs/pull/227))
- 3+ variants in one allele bracket where merging does not apply (closes #221) ([#222](https://github.com/fulcrumgenomics/ferro-hgvs/pull/222))
- dedicated coverage for mixed-accession compound alleles (closes #218) ([#220](https://github.com/fulcrumgenomics/ferro-hgvs/pull/220))

- *(parser)* accept `?con<src>` and `?copy<N>` at unknown position, consistent with `?del`/`?dup`/`?ins`/`?inv` (closes #286)
- *(fasta)* require version-boundary equality in `MmapFastaProvider::get_transcript`'s unversioned-prefix fallback (closes #314)
- *(fasta)* route `FastaProvider::get_sequence` for known contigs through the FASTA path so a transcript registered with a chromosome-colliding id no longer wins over the genomic index (closes #315)

### Fixed

- *(project)* decline rather than emit invalid HGVS when an `NG_`/`LRG_`-parented projection cannot be re-anchored into the parent's own frame: with no chromosomal placement (or an endpoint outside the placed span / an uncertain boundary), `project_to_genomic` now returns `UnsupportedProjection`/`InvalidCoordinates` instead of stamping a chromosome (`NC_`) coordinate under the parent accession. The cross-isoform enumerate path (`project_variant_all`, #646) degrades gracefully — dropping only the unframable genomic axis while keeping the parent-framed coding/protein forms (closes #655).
- *(reference)* make cdot base→version fallback deterministic; previously, when a cdot file contained multiple versions of the same base accession (e.g. `NM_000088.3` and `NM_000088.4`), `base_to_versioned` could resolve to any version across runs due to `HashMap` iteration order, silently shifting CDS coordinates for callers requesting an absent version (closes #583).
- *(edit)* drop brackets on single-payload `InsertedSequence::Complex` `Display` so e.g. `delins[78185355_78199419inv]` round-trips as `delins78185355_78199419inv`, matching HGVS v21 (`DNA/insertion.md:22`, `DNA/inversion.md:39`).
- *(normalize)* extend HGVS codon-frame exception beyond the two-`c.`-SNV case (closes #275, follow-up to #79 / #104): the spec's "two variants separated by one nucleotide, together affecting one amino acid" carve-out now also covers (1) `r.` coding regions, (2) chains of 3+ SNVs where a strict-adjacency merge leaves `prev_a` as a multi-base delins, and (3) `sub`+`del` and `del`+`sub` pairs separated by one unchanged nucleotide.

## [0.6.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.5.0...v0.6.0) - 2026-05-16

### Added

- *(project)* compound alleles + project_all + project_normalized (PR 3: closes #200) ([#204](https://github.com/fulcrumgenomics/ferro-hgvs/pull/204))
- *(cli)* add build-transcript subcommand for FASTA + CDS construction (closes #184) ([#199](https://github.com/fulcrumgenomics/ferro-hgvs/pull/199))
- *(project)* protein consequence for CDS indels (PR 2: indel p. nomenclature) ([#203](https://github.com/fulcrumgenomics/ferro-hgvs/pull/203))
- *(project)* variant-level projection g. -> c./p. (PR 1: foundation + substitutions) ([#202](https://github.com/fulcrumgenomics/ferro-hgvs/pull/202))
- *(loader)* Phase 4 - FASTA validation + complete migration off shims ([#190](https://github.com/fulcrumgenomics/ferro-hgvs/pull/190)) ([#196](https://github.com/fulcrumgenomics/ferro-hgvs/pull/196))
- *(loader,cli)* Phase 3 - diagnostics registry + CLI flags + cross-format tests ([#190](https://github.com/fulcrumgenomics/ferro-hgvs/pull/190)) ([#195](https://github.com/fulcrumgenomics/ferro-hgvs/pull/195))
- *(loader)* Phase 2 - CDS phase/start_codon/stop_codon, UTR-merge, MANE wiring ([#190](https://github.com/fulcrumgenomics/ferro-hgvs/pull/190)) ([#194](https://github.com/fulcrumgenomics/ferro-hgvs/pull/194))

### Fixed

- *(loader)* rewrite GFF/GTF loader pipeline, close #183 ([#191](https://github.com/fulcrumgenomics/ferro-hgvs/pull/191))
- *(mock)* accept version/genome_build metadata keys in from_json (closes #185) ([#198](https://github.com/fulcrumgenomics/ferro-hgvs/pull/198))
- *(normalize)* group consecutive sub-flanks of inv-split into delins ([#182](https://github.com/fulcrumgenomics/ferro-hgvs/pull/182)) ([#188](https://github.com/fulcrumgenomics/ferro-hgvs/pull/188))
- *(normalize)* apply 3' rule to allele-merged del/dup/ins ([#180](https://github.com/fulcrumgenomics/ferro-hgvs/pull/180)) ([#187](https://github.com/fulcrumgenomics/ferro-hgvs/pull/187))

### Other

- *(loader)* Phase 5 - swap parser internals for noodles ([#190](https://github.com/fulcrumgenomics/ferro-hgvs/pull/190)) ([#197](https://github.com/fulcrumgenomics/ferro-hgvs/pull/197))
- *(spec)* refresh vendored hgvs-nomenclature spec snapshot ([#193](https://github.com/fulcrumgenomics/ferro-hgvs/pull/193))

### Added

- *(loader)* unified `load_annotations` entry point with `LoaderConfig`/`LoaderReport`; auto-detects GFF3 vs GTF by extension and content (#191, #194, #195)
- *(loader)* exon-derivation ladder closes single-exon GFF3 without `exon` records (#183 fix) (#191)
- *(loader)* phase-aware CDS bounds with `start_codon` / `stop_codon` precedence; GFF3 input now extends `cds_end` to include the stop codon, matching ferro's downstream convention (fixes a latent off-by-one-codon bug in GFF3 protein conversion) (#194)
- *(loader)* UTR-merge ladder step: synthesizes exons from UTR + CDS when no explicit `exon` records (#194)
- *(loader)* `gene_symbol`, `mane_status`, `refseq_match`, `ensembl_match` extracted from attributes (#194)
- *(loader)* optional FASTA-aware validation: CDS length mod 3 and start codon (ATG/CTG/GTG/TTG) checks when `--fasta` is supplied
- *(error_handling)* 13 loader diagnostic codes (E-LOAD-*, W-LOAD-*) registered for `ferro explain` (#195)
- *(cli)* `ferro convert-gff --strict / --silent / --no-validate-fasta / --diagnostics-json` flags (#195)
- *(strand)* `Strand::Unknown` variant for GFF3 `.` and `?` strand values; transcripts with unknown strand are dropped at load with `E-LOAD-103` (#191)

### Changed

- *(strand)* `Strand`, `LoaderConfig`, and `LoaderReport` are marked `#[non_exhaustive]` for forward compatibility (#191)
- *(loader)* `load_gff3` and `load_gtf` have been removed; use `load_annotations` instead

## [0.5.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.4.1...v0.5.0) - 2026-05-12

### Added

- *(parser)* spec compact-mosaic form across sub/del/dup ([#133](https://github.com/fulcrumgenomics/ferro-hgvs/pull/133)) ([#153](https://github.com/fulcrumgenomics/ferro-hgvs/pull/153))
- *(normalize)* canonicalize con to delins per SVD-WG009 (#81 H1) ([#142](https://github.com/fulcrumgenomics/ferro-hgvs/pull/142))
- *(rna)* support r.spl, r.spl?, r.(spl), r.(spl?) splicing markers (#81 E2) ([#134](https://github.com/fulcrumgenomics/ferro-hgvs/pull/134))

### Fixed

- *(spdi)* emit inversion and repeat as SPDI delins ([#159](https://github.com/fulcrumgenomics/ferro-hgvs/pull/159))
- *(error-handling)* error codes for spec-mandated input rejections ([#115](https://github.com/fulcrumgenomics/ferro-hgvs/pull/115)) ([#152](https://github.com/fulcrumgenomics/ferro-hgvs/pull/152))
- *(spdi)* reference-aware HGVS→SPDI for del/dup/delins ([#158](https://github.com/fulcrumgenomics/ferro-hgvs/pull/158))
- *(parser)* preserve explicit deleted sequence in delins round-trip ([#154](https://github.com/fulcrumgenomics/ferro-hgvs/pull/154))
- *(error-handling)* soft-warn non-canonical input forms at parse time ([#127](https://github.com/fulcrumgenomics/ferro-hgvs/pull/127)) ([#151](https://github.com/fulcrumgenomics/ferro-hgvs/pull/151))
- *(ci)* drop redundant `#[allow(dead_code)]` on `mod common;` ([#178](https://github.com/fulcrumgenomics/ferro-hgvs/pull/178))
- *(spdi)* accept c./n./r./m. variants in HGVS->SPDI conversion ([#157](https://github.com/fulcrumgenomics/ferro-hgvs/pull/157))
- *(normalize)* apply 3' rule across phase-mismatched cyclic rotations for ins ([#132](https://github.com/fulcrumgenomics/ferro-hgvs/pull/132)) ([#155](https://github.com/fulcrumgenomics/ferro-hgvs/pull/155))
- *(error-handling)* surface deprecated stop-codon and frameshift forms as soft-warns ([#125](https://github.com/fulcrumgenomics/ferro-hgvs/pull/125)) ([#150](https://github.com/fulcrumgenomics/ferro-hgvs/pull/150))
- *(error-handling)* wire W1001/W1002/W3001 soft-validation warnings ([#124](https://github.com/fulcrumgenomics/ferro-hgvs/pull/124)) ([#149](https://github.com/fulcrumgenomics/ferro-hgvs/pull/149))
- *(spdi)* recover dup form on SPDI→HGVS for duplicated insertions ([#156](https://github.com/fulcrumgenomics/ferro-hgvs/pull/156))
- *(parser)* extend unknown-phase (;) support to g./n./m./o./p. coord systems ([#123](https://github.com/fulcrumgenomics/ferro-hgvs/pull/123)) ([#148](https://github.com/fulcrumgenomics/ferro-hgvs/pull/148))
- *(normalize)* detect cis-allele edits with coincident bounds (#81 A8) ([#147](https://github.com/fulcrumgenomics/ferro-hgvs/pull/147))
- *(parser)* soft-warn embedded whitespace and zero-width chars ([#128](https://github.com/fulcrumgenomics/ferro-hgvs/pull/128)) ([#145](https://github.com/fulcrumgenomics/ferro-hgvs/pull/145))
- *(normalize)* preserve r.*N UTR flag and translate via cds_end (closes #163) ([#164](https://github.com/fulcrumgenomics/ferro-hgvs/pull/164))
- *(parser)* correct LRG accession variant-type inference and compound-ref handling ([#122](https://github.com/fulcrumgenomics/ferro-hgvs/pull/122)) ([#141](https://github.com/fulcrumgenomics/ferro-hgvs/pull/141))
- *(normalize)* recognize revcomp inv sub-spans within delins ([#166](https://github.com/fulcrumgenomics/ferro-hgvs/pull/166))
- *(normalize)* apply 3'-rule to merged cis-allele deletions (closes #161) ([#162](https://github.com/fulcrumgenomics/ferro-hgvs/pull/162))
- *(normalize)* rewrite revcomp delins as inversion (#81 A2) ([#109](https://github.com/fulcrumgenomics/ferro-hgvs/pull/109))
- *(normalize)* enforce HGVS c. codon-frame exception for repeat notation (#81 B1) ([#110](https://github.com/fulcrumgenomics/ferro-hgvs/pull/110))
- *(normalize)* rewrite empty-insert delins as del (#81 A3) ([#113](https://github.com/fulcrumgenomics/ferro-hgvs/pull/113))
- *(normalize)* degenerate substitution (ref==alt) -> identity (#81 A4) ([#111](https://github.com/fulcrumgenomics/ferro-hgvs/pull/111))

### Other

- *(parser)* triage failure expectations (#174 phase 2) ([#176](https://github.com/fulcrumgenomics/ferro-hgvs/pull/176))
- *(parser)* per-input failure-expectations framework (#174 phase 1) ([#175](https://github.com/fulcrumgenomics/ferro-hgvs/pull/175))
- drop ferro_version from HGVS spec fixture ([#177](https://github.com/fulcrumgenomics/ferro-hgvs/pull/177))
- *(test)* cut Test job runtime ~398s → ~20s ([#173](https://github.com/fulcrumgenomics/ferro-hgvs/pull/173))
- *(test)* correct misleading comment on trans-allele expanded-form test ([#167](https://github.com/fulcrumgenomics/ferro-hgvs/pull/167))
- *(normalize)* tag r. positive bases as Region::Rna (closes #168) ([#169](https://github.com/fulcrumgenomics/ferro-hgvs/pull/169))
- *(allele)* pin trans-phase round-trip across all coord systems and merge-barrier (#81 C1) ([#146](https://github.com/fulcrumgenomics/ferro-hgvs/pull/146))
- *(compound)* pin cross-reference / cross-coord compound round-trip and merge-barrier (#81 H2) ([#143](https://github.com/fulcrumgenomics/ferro-hgvs/pull/143))
- *(error-handling)* audit error codes against HGVS spec sections (#81 L1) ([#137](https://github.com/fulcrumgenomics/ferro-hgvs/pull/137))
- *(parser)* pin gene-selector round-trip end-to-end with Display preservation (#81 I3) ([#135](https://github.com/fulcrumgenomics/ferro-hgvs/pull/135))
- *(mito)* audit heteroplasmy notation; tracking #133 (#81 F2) ([#139](https://github.com/fulcrumgenomics/ferro-hgvs/pull/139))
- *(normalize)* pin RNA path + edge cases for A9 substitution-after-trim ([#114](https://github.com/fulcrumgenomics/ferro-hgvs/pull/114))
- *(protein)* pin p.? unknown-effect round-trip across allele forms and edge cases (#81 D7) ([#136](https://github.com/fulcrumgenomics/ferro-hgvs/pull/136))
- *(protein)* pin p.0 no-product round-trip and adjacent guards (#81 D6) ([#130](https://github.com/fulcrumgenomics/ferro-hgvs/pull/130))
- *(protein)* pin silent `=` round-trip across allele forms and edge cases (#81 D5) ([#131](https://github.com/fulcrumgenomics/ferro-hgvs/pull/131))
- *(mito)* audit m. coord-system parse + wraparound behavior; tracking #129 (#81 F1) ([#138](https://github.com/fulcrumgenomics/ferro-hgvs/pull/138))
- pin normalize() against HGVS v21.0 spec fixture (closes #84) ([#105](https://github.com/fulcrumgenomics/ferro-hgvs/pull/105))

### Added

- *(normalize)* detect cis-allele edits with coincident reference bounds (e.g. `g.[100G>A;100A>C]`) and emit an advisory `OVERLAP_CONFLICTING_EDITS` (`W5002`) warning. Variant output is preserved unchanged. Addresses [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A8.
- *(test)* `expected_warnings` field on the HGVS v21.0 spec-fixture row schema, pinning the warning set ferro emits per row.

### Fixed

- *(normalize)* apply 3' rule across phase-mismatched cyclic rotations for single-copy ins. When the inserted alt is a non-zero cyclic rotation of an adjacent reference repeat unit (e.g. `g.X_(X+1)insTG` against a `GT[3]` tract), shuffle's first-base check (`alt[0] == ref[ins_point]`) failed and the variant never moved. The new helper `insertion_to_duplication` mirrors `insertion_to_repeat`'s rotation iteration for the 1-copy case, so the variant now canonicalizes to a `dup` at the most-3' rotation-aligned position. Closes [#132](https://github.com/fulcrumgenomics/ferro-hgvs/issues/132).
- *(normalize)* recognize a reverse-complement sub-span within a `delins` (synthesized by cis-allele merge OR user-typed) and emit the spec-canonical `inv`, splitting the surrounding span into `[…; inv; …]`. The HGVS edit-priority rule places `inv` above `delins` in the priority order (`general.md:56`) and defines `delins` as the residual when no higher-priority form applies (`delins.md`). For example, `g.[1150T>G;1151C>A;1152C>G]` (over `TCC`) now normalizes to `g.[1150_1151inv;1152C>G]` instead of `g.1150_1152delinsGAG`; `g.[1092G>C;1093G>C]` (over `GG`) normalizes to `g.1092_1093inv` instead of `g.1092_1093delinsCC`. The same rule fires for a user-typed `g.1150_1152delinsGAG`, since the canonical form depends on `(ref, position, alt)` and not on input shape. Applies across all NA coord systems: `g.`, `m.`, `c.` (CDS-proper positions), `n.`, `r.` (T/U-equivalent comparison so `r.` alts with `U` align with transcript ref bytes that contain `T`). Sub-only decomposition (rewriting a non-inv multi-sub `delinsXY` to `[X>...; Y>...]`) is intentionally left out of scope and is a separate spec interpretation question. The codon-frame `c.` merge from issue [#79](https://github.com/fulcrumgenomics/ferro-hgvs/issues/79) is preserved automatically — the synthesized middle base in a codon-frame-merged delins makes a length-2 inv across the middle mathematically impossible. Closes issue [#160](https://github.com/fulcrumgenomics/ferro-hgvs/issues/160).
- *(normalize)* rewrite degenerate substitutions (ref == alt, e.g. `c.100A>A`) to identity (`=`) per HGVS v21 spec, which marks `c.X>X` as "not allowed" (`docs/recommendations/DNA/other.md`). The rule is purely syntactic on the edit's stated bases, so it fires in both the full-normalization path and the no-reference canonicalization path — `c.123C>C` rewrites to `c.123=` regardless of provider availability. ([#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) A4)
- *(normalize)* rewrite a `delins` whose inserted sequence is empty as a deletion, per the HGVS spec requirement that an empty insert is semantically a deletion and must be rendered as `del`. The rewritten deletion is then 3'-shifted under the standard del rule. Issue [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A3.
- *(normalize)* rewrite a `delins` whose inserted sequence is the reverse complement of the deleted reference as an inversion, per the HGVS spec definition of `inv`. The complementary-outer-bases shortening rule applies to the result so that a `delins`-encoded inversion produces the same canonical output as a directly-encoded `inv`. Issue [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A2.
- *(normalize)* canonicalize `delins` to the minimal HGVS form by trimming any shared prefix/suffix between the inserted sequence and the deleted reference, then reclassifying the residual edit as substitution / deletion / insertion / inversion / smaller `delins` per the sub > del > inv > dup > ins priority. For example, `c.1_4delinsAAGC` against ref `ATGC` collapses to `c.2T>A`; `c.1_4delinsAC` against ref `ATGC` collapses to `c.2_3del`; `c.1_3delinsATCG` against ref `ATG` collapses to `c.2_3insC`. Extension to issue [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A2.
- *(normalize)* apply inversion shortening to direct `inv` inputs. `NaEdit::Inversion` was missing from `needs_normalization()`, so direct inversions bypassed `normalize_na_edit` entirely and the `shorten_inversion()` / identity-collapse logic was never exercised. Direct `inv` variants now also emit minimal notation after shortening (no stale explicit `sequence`/`length` from the input).
- *(normalize)* HGVS spec compliance ([#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) B1): repeat-notation rewrites in `c.` (coding DNA) context now enforce the spec's codon-frame exception — repeat notation requires `unit_len % 3 == 0`. Previously `insertion_to_repeat`, `deletion_to_repeat`, `duplication_to_repeat`, and `normalize_repeat` could emit `c.X_YA[N]`, `c.X_YAT[N]`, etc. for non-codon-aligned units, violating the spec (`docs/recommendations/DNA/repeated.md`). Variants in `c.` with non-codon-aligned units now retain the spec-prescribed alternative form: `dup` for 1 added copy, `ins<literal>` for ≥2 added copies, and plain `del` for contractions of ≥2 unit copies.

### Changed (public API)

- `insertion_to_repeat`, `duplication_to_repeat`, and `normalize_repeat` in `ferro_hgvs::normalize::rules` gain an `is_coding: bool` parameter to drive the codon-frame gate. (`deletion_to_repeat` is `pub(crate)` and gains the same parameter as an internal change.)
- `RepeatNormResult::Insertion { start, end, sequence }` and `DupToRepeatResult::GatedInsertion { start, end, sequence }` variants added so the rule layer can route gated rewrites to the spec-canonical literal `ins` form.
- *(normalize)* `NormalizationWarning` is now a sum type (`RefSeqMismatch` / `OverlapConflict`). Each warning code carries only the fields relevant to it. Read sites migrate from `.code` field access to `.code()` method and pattern-matching on the variant.

### Changed

- Internal: the four `delins_is_*` boolean helpers in `normalize::rules` are unified into one `canonicalize_delins()` function returning a `DelinsCanonical` enum, expressing HGVS edit-priority (sub > del > inv > dup > ins) in a single decision tree. The unreachable second delins arm in `normalize_na_edit` has been removed.

## [0.4.1](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.4.0...v0.4.1) - 2026-05-04

### Added

- *(python)* load reference data via Normalizer.from_manifest and extended from_json ([#86](https://github.com/fulcrumgenomics/ferro-hgvs/pull/86))
- *(normalize)* merge consecutive edits in cis alleles per HGVS spec ([#80](https://github.com/fulcrumgenomics/ferro-hgvs/pull/80))

### Fixed

- *(normalize)* codon-frame exception for c. SNV pairs separated by one nucleotide ([#79](https://github.com/fulcrumgenomics/ferro-hgvs/pull/79)) ([#104](https://github.com/fulcrumgenomics/ferro-hgvs/pull/104))
- *(normalize)* 5'UTR CDS↔tx off-by-one collapsed UTR del to c.? ([#97](https://github.com/fulcrumgenomics/ferro-hgvs/pull/97)) ([#102](https://github.com/fulcrumgenomics/ferro-hgvs/pull/102))
- *(normalize)* merge same-region UTR adjacency in cis alleles ([#89](https://github.com/fulcrumgenomics/ferro-hgvs/pull/89)) ([#103](https://github.com/fulcrumgenomics/ferro-hgvs/pull/103))
- *(normalize)* minus-strand intronic ref-base orientation ([#98](https://github.com/fulcrumgenomics/ferro-hgvs/pull/98)) ([#100](https://github.com/fulcrumgenomics/ferro-hgvs/pull/100))
- rewrite delins as identity when insert matches reference ([#78](https://github.com/fulcrumgenomics/ferro-hgvs/pull/78))
- rewrite single-base delins as substitution per HGVS priority ([#77](https://github.com/fulcrumgenomics/ferro-hgvs/pull/77))
- Emit single-variant alleles in bare spec form ([#76](https://github.com/fulcrumgenomics/ferro-hgvs/pull/76))

### Other

- unpin CI Rust toolchain and fix 1.95 clippy lints ([#108](https://github.com/fulcrumgenomics/ferro-hgvs/pull/108))
- dup 3'-shift coverage matrix (#81 A6) ([#107](https://github.com/fulcrumgenomics/ferro-hgvs/pull/107))
- del 3'-shift coverage matrix + tandem-repeat del canonical-form fix (#81 A5, B2) ([#106](https://github.com/fulcrumgenomics/ferro-hgvs/pull/106))
- tighten coverage_gap_tests assertions; restore intronic-ins coverage ([#94](https://github.com/fulcrumgenomics/ferro-hgvs/pull/94)) ([#99](https://github.com/fulcrumgenomics/ferro-hgvs/pull/99))
- *(coverage)* restore intronic insertion tests dropped in #93 ([#95](https://github.com/fulcrumgenomics/ferro-hgvs/pull/95))
- ins 3'-shift coverage matrix + tandem-repeat ins canonical-form fix (#81 A1, A7) ([#93](https://github.com/fulcrumgenomics/ferro-hgvs/pull/93))

### Added

- *(test)* HGVS v21.0 spec normalization regression fixture pinning ferro's current `normalize()` output for every variant string in the spec (823 rows). Each row carries a status (`preserved` / `diverges` / `correctly-rejected` / `false-acceptance` / `parse-error` / `needs-reference`) derived from `(parse_ok, normalize_ok, spec_expected, current == spec_expected)`. `spec_expected: null` is the sentinel for "spec rejects this," set automatically for inputs the spec marks via `<code class="invalid">…</code>` (35 occurrences in v21.0) and overridable by hand. Bare illustrative fragments like `c.1083A>C` get a default accession prepended per coord system so they parse — the prefixed form is recorded in a separate `input_prefixed` field. The fixture is generated by `examples/generate_spec_fixture.rs` from the vendored spec at `assets/hgvs-nomenclature/` (git submodule pinned to tag `21.0.0`). Hand overrides live in `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json` and accept `status`, `spec_expected` (string / null / absent), `input_prefixed`, `requires_reference`, and `todo`. CI runs the generator in `--check` mode to enforce byte-identical regeneration. Closes [#84](https://github.com/fulcrumgenomics/ferro-hgvs/issues/84); companion to [#83](https://github.com/fulcrumgenomics/ferro-hgvs/issues/83).
- 3'-shift coverage matrix for insertion variants
  (`tests/ins_shift_matrix.rs`): 84 rstest cases across all 7 nucleotide
  coord-system / strand combinations × 8 shuffle scenarios, with a
  shared `SyntheticBuilder` fixture helper at `tests/common/synthetic.rs`
  (reusable by future del / dup / repeat-notation matrices). Issue #81
  item A7.
- Tightened `tests/coverage_gap_tests.rs` from canary-shaped
  (`contains(...)`, `is_ok() || is_err()`) assertions to strict
  `assert_eq!` against the exact normalizer output. Builds on PR #95
  (intronic-insertion restoration) by extending the same tightening
  pass to the remaining gap-test buckets, and adds cross-references
  and procedural detail to the restored intronic-insertion comments.
  Several locked outputs are suspected-buggy and carry `FIXME(#NN)`
  comments pointing to follow-up tracking issues — most converge on
  #98 (minus-strand intronic reference-base orientation, identified
  as a likely common root cause behind A1 / A5 / A6 misfires on
  minus-strand intronic positions and the PR #78 identity-rewrite
  missing on minus strand), with #96 (wrong-strand repeat unit
  emission on minus-strand multi-base dup) and #97 (`c.?del` collapse
  on minus-strand 5'UTR del) as separate sub-symptoms. True
  boundary-spanning panic-canary tests are intentionally left as
  `is_ok() || is_err()` per the issue's scope. Issue #94.
- *(normalize)* Merge consecutive sub-variants in cis alleles into a single delins per HGVS spec. `g.[1000G>A;1001A>C]` now normalizes to `g.1000_1001delinsAC`; `g.[1000del;1001del]` to `g.1000_1001del`. Covers `g./c./n./r./m.` coordinate systems, sub/del/delins/ins edit combinations, chains, and same-boundary insertion pairs. Non-adjacent variants, intronic/UTR boundaries, uncertain edits, and non-`Literal` insertion payloads are barriers and pass through unchanged. The codon-frame exception (one-nt gap within a codon) is tracked separately in [#79](https://github.com/fulcrumgenomics/ferro-hgvs/issues/79). ([#80](https://github.com/fulcrumgenomics/ferro-hgvs/pull/80))

### Fixed

- *(normalize)* Cis-allele consecutive-edit merging now collapses
  adjacent sub-variants *within* a UTR / upstream / downstream region,
  not only within the CDS / transcript body. The PR #80 implementation
  rejected every `is_5utr() / is_3utr() / is_upstream() / is_downstream()`
  position outright, so inputs like `c.[-2A>G;-1C>T]` (both 5'UTR) and
  `c.[*1A>G;*2C>T]` (both 3'UTR) round-tripped unchanged even though
  HGVS allows ranges within those regions (`c.-2_-1`, `c.*1_*2`). The
  fix replaces the `Option<u64>` position keys with a `Region`-tagged
  `(Region, i64)` axis (covering CDS, 5'UTR, 3'UTR, transcript-body,
  upstream, downstream, and genomic / mitochondrial); merge eligibility
  becomes "same region + integer adjacency on the region's axis".
  `build_cds_merged` / `build_tx_merged` / `build_rna_merged` consume
  the region tag to reconstruct the right `CdsPos` / `TxPos` / `RnaPos`
  shape (negative base for 5'UTR / upstream, `utr3` / `downstream` flag
  for 3'UTR / downstream). Cross-region pairs (`c.[-1A>G;1A>T]` 5'UTR↔CDS,
  `c.[40C>T;*1A>G]` CDS↔3'UTR, …) still correctly do not merge.
  Issue #89.
- *(normalize)* CDS ↔ transcript coordinate mappings now respect the
  HGVS no-c.0 numbering rule for 5'UTR positions. The forward mapping
  (`Normalizer::cds_to_tx_pos`, `convert::coding::cds_to_transcript_pos`)
  previously computed `tx = cds_start + base - 1` for negative `base`,
  which double-counted the gap between c.-1 and c.1 and emitted a tx
  position one base 5' of the true location. The inverse mapping
  (`Normalizer::tx_to_cds_pos`) had the mirror bug: tx positions one
  base before `cds_start` mapped to `base = 0`, which `CdsPos::Display`
  renders as `c.?` (`CDS_BASE_UNKNOWN`). The most visible symptom was a
  5'UTR single-base deletion on a minus-strand transcript collapsing
  to `c.?del` instead of resolving to a real position (e.g. `c.-1del`
  after 3'-shifting within a UTR homopolymer). Forward and inverse are
  now `tx = cds_start + base` and `base = tx - cds_start` for negative
  / pre-cds_start positions, matching the spec's "c.-1 is one base 5'
  of c.1" rule and the existing exon-aware mapper at
  `convert::mapper::cds_to_tx` / `tx_to_cds`. Issue #97.
- *(normalize)* Cis-allele consecutive-edit merging now also collapses
  the codon-frame exception case: two `c.` exonic SNVs in the CDS
  proper, separated by exactly one nucleotide, that fall within the
  same codon merge into a single delins with the unchanged middle
  reference base preserved verbatim — per HGVS spec
  (`c.[145C>T;147C>G]` → `c.145_147delinsTGG`, where the middle base
  is the reference at `c.146`). Eligibility is narrow: same accession,
  both endpoints in `Region::Cds`, gap-of-one on the axis, both prev
  and next are single-base SUB anchors, and `(prev-1)/3 == (next-1)/3`
  (same codon). The unchanged middle base is fetched via the
  `ReferenceProvider` threaded into `merge_consecutive_edits`; if no
  transcript is registered or the position is out of range, the merge
  is silently declined and the variants pass through unchanged.
  Codon-frame–merged delins continue to participate in the
  strictly-consecutive walk, so a third SNV one base 3' of the pair
  (`c.[10A>G;12A>C;13A>T]`) still folds into the running delins.
  Cross-codon (`c.[3G>T;5A>C]`), gap-of-two (`c.[10A>G;13A>C]`), and
  non-CDS pairs (`g.`, UTR, `n.`, `r.`) all correctly do not merge.
  Issue #79.
- *(normalize)* Minus-strand intronic normalization now reads the
  reference window in transcript-view orientation. `normalize_intronic_cds`
  and `normalize_intronic_tx` previously passed the genomic-strand bytes
  fetched from `get_genomic_sequence` directly into `normalize_na_edit`
  alongside the variant's transcript-view edit alt; on minus-strand
  transcripts the two were mis-oriented, defeating every rule that
  compared the alt against the local reference window. The fix
  reverse-complements the genomic window and flips the relative
  positions / shuffle boundaries on minus strand before normalization,
  then maps the resulting positions back to the genomic frame for the
  CDS / tx coordinate conversion. As a single root-cause fix this
  resurfaces #81 A1 / A5 / A6 canonicalization, the PR #78 delins-as-
  identity rewrite, and the transcript-view repeat-unit letter on
  minus-strand intronic positions — all of which had been observed to
  misfire in #94's locking pass. Issue #98.
- Insertions that add ≥2 copies of a multi-base tandem repeat unit now
  emit repeat notation (`unit[N+k]`) instead of a duplication of the
  inserted sequence, per HGVS spec ("when more than one additional copies
  are inserted directly 3' of the original copy, the change is indicated
  using the format for Repeated sequences"). Single-unit additions remain
  `dup`. Issue #81 item A7.
- `MockProvider::get_sequence` now falls through to genomic contig
  lookup when the id is not a transcript, matching `FastaProvider`'s
  behavior. This unblocks 3'-shift normalization for genomic test
  fixtures that register only `add_genomic_sequence`. Issue #81 item A7.

## [0.4.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.3.0...v0.4.0) - 2026-04-30

### Added

- *(python)* add poethepoet task runner for Python dev workflow ([#53](https://github.com/fulcrumgenomics/ferro-hgvs/pull/53))
- Add variant accessor properties to Python bindings ([#49](https://github.com/fulcrumgenomics/ferro-hgvs/pull/49))

### Fixed

- accept gene selectors on non-RefSeq accessions ([#70](https://github.com/fulcrumgenomics/ferro-hgvs/pull/70))
- Use HGVS spec compact form for allele Display ([#48](https://github.com/fulcrumgenomics/ferro-hgvs/pull/48))

### Other

- *(readme)* vendor Fulcrum logo and use absolute URLs ([#71](https://github.com/fulcrumgenomics/ferro-hgvs/pull/71))
- bump pyo3 0.23 → 0.28 and add Python 3.14 wheels ([#57](https://github.com/fulcrumgenomics/ferro-hgvs/pull/57))
- *(python)* use dependency-groups, stricter mypy, and --locked in CI ([#52](https://github.com/fulcrumgenomics/ferro-hgvs/pull/52))
- publish Python wheels to PyPI via Trusted Publishing ([#58](https://github.com/fulcrumgenomics/ferro-hgvs/pull/58))
- *(python)* drop Python 3.8 and 3.9 support ([#55](https://github.com/fulcrumgenomics/ferro-hgvs/pull/55))
- *(python)* add uv lockfile for reproducible dev environment ([#50](https://github.com/fulcrumgenomics/ferro-hgvs/pull/50))
- build Python wheels and attach to GitHub Releases ([#39](https://github.com/fulcrumgenomics/ferro-hgvs/pull/39))
- *(prepare)* Modularize ReferenceManifest ([#44](https://github.com/fulcrumgenomics/ferro-hgvs/pull/44))
- switch reqwest from native-tls to rustls-tls ([#38](https://github.com/fulcrumgenomics/ferro-hgvs/pull/38))

### Changed

- Allele `Display` now emits HGVS spec-correct compact form (`ACC:c.[edit1;edit2]`) when sub-variants share an accession and coordinate type, instead of the expanded form (`[ACC:c.edit1;ACC:c.edit2]`). Mixed-accession alleles and alleles containing the per-variant unknown form (`c.?`, `r.?`, etc.) still emit the expanded form. Downstream consumers parsing the previous expanded output (including Python `str(variant)`) will see the new format. ([#48](https://github.com/fulcrumgenomics/ferro-hgvs/pull/48))

### Fixed

- Prevent panic when `NullAllele`/`UnknownAllele` are used as sub-variants in an allele ([#48](https://github.com/fulcrumgenomics/ferro-hgvs/pull/48))

## [0.3.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.2.0...v0.3.0) - 2026-03-30

### Added

- Bincode serialization for fast cdot loading ([#26](https://github.com/fulcrumgenomics/ferro-hgvs/pull/26))
- Compound reference syntax `NC_*(NM_*):c.…` support ([#21](https://github.com/fulcrumgenomics/ferro-hgvs/pull/21))
- GRCh37 cdot transcript metadata download in `ferro prepare` ([#24](https://github.com/fulcrumgenomics/ferro-hgvs/pull/24))

### Fixed

- Genomic-to-coding coordinate conversion ([#25](https://github.com/fulcrumgenomics/ferro-hgvs/pull/25))
- Coding-order positions after genomic-space normalization ([#20](https://github.com/fulcrumgenomics/ferro-hgvs/pull/20))
- Mutalyzer configuration across web service and benchmarks ([#19](https://github.com/fulcrumgenomics/ferro-hgvs/pull/19))
- Serde defaults for config structs and deploy workflow ([#27](https://github.com/fulcrumgenomics/ferro-hgvs/pull/27))
- Missing `in_memory` field in HgvsRsConfig construction ([#28](https://github.com/fulcrumgenomics/ferro-hgvs/pull/28))

### Changed

- Replaced vendored hgvs-rs with published crate v0.20.1 ([#22](https://github.com/fulcrumgenomics/ferro-hgvs/pull/22))

## [0.2.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.1.0...v0.2.0) - 2026-03-23

### Fixed

- derive Python package version from Cargo.toml via maturin ([#14](https://github.com/fulcrumgenomics/ferro-hgvs/pull/14))
- allow valid HGVS characters in web service input validation ([#13](https://github.com/fulcrumgenomics/ferro-hgvs/pull/13))
- work-stealing for mutalyzer normalize shards ([#4](https://github.com/fulcrumgenomics/ferro-hgvs/pull/4))
- CIGAR-aware CDS-to-transcript coordinate mapping ([#7](https://github.com/fulcrumgenomics/ferro-hgvs/pull/7))
- prevent integer overflow and handle edge cases in normalization ([#6](https://github.com/fulcrumgenomics/ferro-hgvs/pull/6))
- prevent panics on malformed HGVS patterns during normalization ([#5](https://github.com/fulcrumgenomics/ferro-hgvs/pull/5))

### Other

- Add Bioconda and Zenodo badges ([#3](https://github.com/fulcrumgenomics/ferro-hgvs/pull/3))
- refresh test fixtures and data extraction scripts ([#8](https://github.com/fulcrumgenomics/ferro-hgvs/pull/8))
- release v0.1.0

## [0.1.0](https://github.com/fulcrumgenomics/ferro-hgvs/releases/tag/v0.1.0) - 2026-02-19

### Other

- Initial commit
