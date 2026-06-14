//! HGVS variant types
//!
//! The main variant enum and supporting types for representing
//! complete HGVS variant descriptions.

use super::edit::{NaEdit, ProteinEdit};
use super::interval::{CdsInterval, GenomeInterval, ProtInterval, RnaInterval, TxInterval};
use super::uncertainty::Mu;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::sync::Arc;

/// Interned accession prefixes for common RefSeq and Ensembl types.
/// Using static `Arc<str>` avoids allocating new strings for common prefixes.
mod interned {
    use std::sync::{Arc, OnceLock};

    macro_rules! interned_prefix {
        ($name:ident, $value:literal) => {
            pub fn $name() -> Arc<str> {
                static INSTANCE: OnceLock<Arc<str>> = OnceLock::new();
                INSTANCE.get_or_init(|| Arc::from($value)).clone()
            }
        };
    }

    // RefSeq prefixes
    interned_prefix!(nc, "NC");
    interned_prefix!(ng, "NG");
    interned_prefix!(nt, "NT");
    interned_prefix!(nw, "NW");
    interned_prefix!(nm, "NM");
    interned_prefix!(nr, "NR");
    interned_prefix!(np, "NP");
    interned_prefix!(xm, "XM");
    interned_prefix!(xr, "XR");
    interned_prefix!(xp, "XP");

    // Ensembl prefixes
    interned_prefix!(enst, "ENST");
    interned_prefix!(ensg, "ENSG");
    interned_prefix!(ensp, "ENSP");
    interned_prefix!(ense, "ENSE");
    interned_prefix!(ensr, "ENSR");

    // LRG prefix
    interned_prefix!(lrg, "LRG");

    // Empty string for assembly refs
    interned_prefix!(empty, "");

    /// Get an interned prefix if available, otherwise return None
    #[inline]
    pub fn get_prefix(s: &str) -> Option<Arc<str>> {
        match s {
            "NC" => Some(nc()),
            "NG" => Some(ng()),
            "NT" => Some(nt()),
            "NW" => Some(nw()),
            "NM" => Some(nm()),
            "NR" => Some(nr()),
            "NP" => Some(np()),
            "XM" => Some(xm()),
            "XR" => Some(xr()),
            "XP" => Some(xp()),
            "ENST" => Some(enst()),
            "ENSG" => Some(ensg()),
            "ENSP" => Some(ensp()),
            "ENSE" => Some(ense()),
            "ENSR" => Some(ensr()),
            "LRG" => Some(lrg()),
            "" => Some(empty()),
            _ => None,
        }
    }
}

/// Accession number with optional version
///
/// Uses `Arc<str>` internally for cheap cloning - accession strings are frequently
/// cloned during parsing and normalization, and `Arc<str>` reduces this to a
/// simple reference count increment.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Accession {
    /// Accession prefix (e.g., "NC", "NM", "NP", "ENST")
    pub prefix: Arc<str>,
    /// Accession number
    pub number: Arc<str>,
    /// Version number (e.g., .3 in NM_000088.3)
    pub version: Option<u32>,
    /// Whether this is an Ensembl-style accession (no underscore separator)
    #[serde(default)]
    pub ensembl_style: bool,
    /// Assembly name (e.g., "GRCh37", "GRCh38") for assembly/chromosome notation
    #[serde(default)]
    pub assembly: Option<Arc<str>>,
    /// Chromosome name (e.g., "chr1", "chr23") for assembly/chromosome notation
    #[serde(default)]
    pub chromosome: Option<Arc<str>>,
    /// Genomic context accession for compound reference syntax
    /// e.g., NC_000013.11 in NC_000013.11(NM_004119.3):c.…
    #[serde(default)]
    pub genomic_context: Option<Box<Accession>>,
}

impl Accession {
    /// Intern a prefix string, using a static Arc if available
    #[inline]
    fn intern_prefix(prefix: impl AsRef<str> + Into<Arc<str>>) -> Arc<str> {
        // Look up the interned static *before* allocating: for the common known
        // prefixes (NC/NM/NP/ENST/…) this avoids allocating an `Arc<str>` on every
        // parse just to discard it for the static. Only unknown prefixes allocate.
        interned::get_prefix(prefix.as_ref()).unwrap_or_else(|| prefix.into())
    }

    pub fn new(
        prefix: impl AsRef<str> + Into<Arc<str>>,
        number: impl Into<Arc<str>>,
        version: Option<u32>,
    ) -> Self {
        let prefix = Self::intern_prefix(prefix);
        let ensembl_style = Self::is_ensembl_prefix(&prefix);
        // LRG accessions never carry a version per HGVS spec (SVD-WG008).
        // Drop any provided version so Display canonicalizes to the spec form.
        let version = if Self::is_lrg_prefix(&prefix) {
            None
        } else {
            version
        };
        Self {
            prefix,
            number: number.into(),
            version,
            ensembl_style,
            assembly: None,
            chromosome: None,
            genomic_context: None,
        }
    }

    /// Create with explicit Ensembl style setting
    pub fn with_style(
        prefix: impl AsRef<str> + Into<Arc<str>>,
        number: impl Into<Arc<str>>,
        version: Option<u32>,
        ensembl_style: bool,
    ) -> Self {
        let prefix = Self::intern_prefix(prefix);
        // LRG accessions never carry a version per HGVS spec (SVD-WG008).
        let version = if Self::is_lrg_prefix(&prefix) {
            None
        } else {
            version
        };
        Self {
            prefix,
            number: number.into(),
            version,
            ensembl_style,
            assembly: None,
            chromosome: None,
            genomic_context: None,
        }
    }

    /// Create an assembly/chromosome reference (e.g., GRCh37(chr23))
    pub fn from_assembly(assembly: impl Into<Arc<str>>, chromosome: impl Into<Arc<str>>) -> Self {
        Self {
            prefix: interned::empty(),
            number: interned::empty(),
            version: None,
            ensembl_style: false,
            assembly: Some(assembly.into()),
            chromosome: Some(chromosome.into()),
            genomic_context: None,
        }
    }

    /// Create a compound reference with genomic context
    /// e.g., NC_000013.11(NM_004119.3) where `context` is NC_000013.11
    /// and the primary accession is NM_004119.3
    pub fn with_genomic_context(mut self, context: Accession) -> Self {
        self.genomic_context = Some(Box::new(context));
        self
    }

    /// Check if this is an assembly/chromosome style reference
    pub fn is_assembly_ref(&self) -> bool {
        self.assembly.is_some() && self.chromosome.is_some()
    }

    /// Check if a prefix is an Ensembl-style prefix.
    ///
    /// Ensembl stable IDs have the shape `ENS<species_code><feature>` where:
    /// - `species_code` is 0–5 uppercase letters (`""` for human, `MUS` for
    ///   mouse, `RNO` for rat, `BTA` for cow, `DAR` for zebrafish, `CAF`
    ///   for dog, `SSC` for pig, etc.).
    /// - `feature` is one of `G` (gene), `T` (transcript), `P` (peptide),
    ///   `E` (exon), or `R` (regulatory).
    ///
    /// Examples that match: `ENST`, `ENSP`, `ENSMUST`, `ENSMUSG`,
    /// `ENSRNOT`, `ENSBTAG`, `ENSDARP`.
    ///
    /// Reference: <https://useast.ensembl.org/info/genome/stable_ids/index.html>.
    pub fn is_ensembl_prefix(prefix: &str) -> bool {
        // Must start with "ENS" and have at least one feature letter after it.
        let Some(body) = prefix.strip_prefix("ENS") else {
            return false;
        };
        if body.is_empty() {
            return false;
        }
        // body = species_code + feature. Per the Ensembl stable-ID convention
        // the species code is 0..=5 uppercase letters and the feature is 1
        // letter, so body must be at most 6 chars. This rejects arbitrarily
        // long uppercase prefixes that would otherwise pass the suffix check.
        if body.len() > 6 {
            return false;
        }
        // All characters in the body must be uppercase ASCII letters.
        if !body.chars().all(|c| c.is_ascii_uppercase()) {
            return false;
        }
        // The last letter is the feature suffix.
        matches!(
            body.as_bytes()[body.len() - 1],
            b'G' | b'T' | b'P' | b'E' | b'R'
        )
    }

    /// Check if this accession is an Ensembl accession
    pub fn is_ensembl(&self) -> bool {
        self.ensembl_style || Self::is_ensembl_prefix(&self.prefix)
    }

    /// Check if a prefix is the canonical LRG prefix (`"LRG"`, uppercase).
    ///
    /// Per HGVS spec the LRG prefix is case-sensitive (lowercase `lrg_` is not
    /// recognized as an LRG accession; it parses as a generic chromosome
    /// name). See `assets/hgvs-nomenclature/docs/background/refseq.md`.
    pub fn is_lrg_prefix(prefix: &str) -> bool {
        prefix == "LRG"
    }

    /// Check if this accession is an LRG accession.
    pub fn is_lrg(&self) -> bool {
        Self::is_lrg_prefix(&self.prefix)
    }

    /// Validate Ensembl accession format
    /// Returns true if valid, false if invalid
    pub fn validate_ensembl(&self) -> bool {
        if !self.is_ensembl() {
            return true; // Not an Ensembl ID, skip validation
        }
        // Ensembl IDs should have 11-15 digit numbers
        let digit_count = self.number.len();
        (11..=15).contains(&digit_count) && self.number.chars().all(|c| c.is_ascii_digit())
    }

    /// Infer the *primary* expected variant coordinate type from the accession
    /// prefix and (for LRG) the trailing discriminator inside `number`.
    ///
    /// This is informational, not a validator: it returns the canonical
    /// coordinate type for the accession class. The accession may legally
    /// appear with other coordinate types in the variant block — for
    /// example, `NM_*` always returns `"c"` even though `NM_*:n.…` and
    /// `NM_*:r.…` parse correctly, and `LRG_<N>t<M>` always returns `"c"`
    /// even though `LRG_<N>t<M>:n.…` (non-coding) and
    /// `LRG_<N>t<M>:r.…` (RNA) are valid per the HGVS spec.
    ///
    /// LRG accessions take three forms (see https://www.lrg-sequence.org/faq/
    /// and `assets/hgvs-nomenclature/docs/background/refseq.md`):
    /// - `LRG_<N>`        — the LRG genomic record itself; uses `g.` coordinates.
    /// - `LRG_<N>t<M>`    — transcript M of LRG_N; primary `c.` (coding).
    /// - `LRG_<N>p<M>`    — protein M of LRG_N; uses `p.` coordinates.
    ///
    /// The `t<M>` / `p<M>` discriminator is captured inside `self.number`
    /// (e.g. `prefix="LRG"`, `number="1t1"`), so the dispatch must inspect the
    /// trailing discriminator rather than only the prefix.
    pub fn inferred_variant_type(&self) -> Option<&'static str> {
        match &*self.prefix {
            // RefSeq genomic
            "NC" | "NG" | "NT" | "NW" => Some("g"),
            // RefSeq transcript/mRNA
            "NM" => Some("c"),
            // RefSeq non-coding RNA
            "NR" => Some("n"),
            // RefSeq protein
            "NP" => Some("p"),
            // Ensembl transcript
            "ENST" => Some("c"),
            // Ensembl gene (genomic)
            "ENSG" => Some("g"),
            // Ensembl protein
            "ENSP" => Some("p"),
            // LRG: dispatch on the t<M> / p<M> discriminator captured in `number`.
            // See `is_lrg_prefix` / `is_lrg` for the prefix-class predicate.
            p if Self::is_lrg_prefix(p) => Some(Self::lrg_inferred_variant_type(&self.number)),
            // UniProt protein (single letter prefix: O, P, Q, A-N, R-Z)
            p if p.len() == 1 && p.chars().next().is_some_and(|c| c.is_ascii_uppercase()) => {
                Some("p")
            }
            _ => None,
        }
    }

    /// Map an LRG `number` field (e.g. `"1"`, `"1t1"`, `"1p1"`) to its variant type.
    ///
    /// The discriminator is `t<digits>` for transcripts (coding, `c.`) or
    /// `p<digits>` for protein products (`p.`). A bare numeric `number`
    /// (the LRG record itself) maps to genomic (`g.`). Any malformed input
    /// falls back to `g`, matching the prior behavior for unrecognized shapes.
    fn lrg_inferred_variant_type(number: &str) -> &'static str {
        let bytes = number.as_bytes();
        // Find the discriminator letter, if any. A valid LRG number is one or
        // more digits, optionally followed by `t` or `p` and one or more digits.
        let disc_pos = bytes.iter().position(|&b| b == b't' || b == b'p');
        match disc_pos {
            None => "g",
            Some(pos) => {
                // Validate: digits before, digits after, no further non-digits.
                let before_ok = pos > 0 && bytes[..pos].iter().all(|b| b.is_ascii_digit());
                let after = &bytes[pos + 1..];
                let after_ok = !after.is_empty() && after.iter().all(|b| b.is_ascii_digit());
                if before_ok && after_ok {
                    if bytes[pos] == b't' {
                        "c"
                    } else {
                        "p"
                    }
                } else {
                    // Malformed discriminator — be conservative and treat as genomic.
                    "g"
                }
            }
        }
    }

    /// Check if this is a UniProt accession
    pub fn is_uniprot(&self) -> bool {
        // UniProt: single uppercase letter prefix, 5-character number, no version
        self.prefix.len() == 1
            && self
                .prefix
                .chars()
                .next()
                .is_some_and(|c| c.is_ascii_uppercase())
            && self.number.len() == 5
            && self.number.chars().all(|c| c.is_ascii_alphanumeric())
    }

    /// Returns true if this accession is a non-coding RNA transcript reference
    /// (`NR_` curated or `XR_` predicted). Such references have no CDS and are
    /// not genomic, so only `n.` and `r.` coordinate systems are valid on them
    /// (a `c.`/`g.`/`m.`/`o.` description is a coordinate-system mismatch — #486).
    pub fn is_noncoding_rna(&self) -> bool {
        matches!(&*self.prefix, "NR" | "XR")
    }

    /// Returns true if this is a known human mitochondrial reference accession.
    ///
    /// `NC_012920` is the GRCh38 rCRS mitochondrion; `NC_001807` is the older
    /// (deprecated but still seen) rCRS draft. HGVS requires the `m.` coordinate
    /// system for these accessions, so `normalize()` coerces a `g.` variant on
    /// one of them to `m.`.
    pub fn is_mitochondrial(&self) -> bool {
        let prefix: &str = self.prefix.as_ref();
        let number: &str = self.number.as_ref();
        matches!((prefix, number), ("NC", "012920") | ("NC", "001807"))
    }

    /// Full accession string without version
    pub fn base(&self) -> String {
        // Assembly/chromosome notation: GRCh37(chr23)
        if let (Some(assembly), Some(chromosome)) = (&self.assembly, &self.chromosome) {
            return format!("{}({})", assembly, chromosome);
        }
        let base = if self.ensembl_style {
            format!("{}{}", self.prefix, self.number)
        } else {
            format!("{}_{}", self.prefix, self.number)
        };
        // Compound reference: context(base)
        if let Some(ctx) = &self.genomic_context {
            format!("{}({})", ctx.full(), base)
        } else {
            base
        }
    }

    /// Full accession string with version if present
    pub fn full(&self) -> String {
        // Assembly/chromosome notation doesn't have versions
        if self.is_assembly_ref() {
            return self.base();
        }
        let full = match self.version {
            Some(v) => {
                if self.ensembl_style {
                    format!("{}{}.{}", self.prefix, self.number, v)
                } else {
                    format!("{}_{}.{}", self.prefix, self.number, v)
                }
            }
            None => {
                if self.ensembl_style {
                    format!("{}{}", self.prefix, self.number)
                } else {
                    format!("{}_{}", self.prefix, self.number)
                }
            }
        };
        // Compound reference: context(full)
        if let Some(ctx) = &self.genomic_context {
            format!("{}({})", ctx.full(), full)
        } else {
            full
        }
    }

    /// Returns the bare accession string for provider/transcript lookups,
    /// stripping any genomic context wrapper. For compound references like
    /// `NC_000013.11(NM_004119.3)`, this returns just `"NM_004119.3"`.
    pub fn transcript_accession(&self) -> String {
        // Assembly/chromosome notation (e.g. GRCh38(chr1)) uses base() directly
        if self.is_assembly_ref() {
            return self.base();
        }

        match self.version {
            Some(v) => {
                if self.ensembl_style {
                    format!("{}{}.{}", self.prefix, self.number, v)
                } else {
                    format!("{}_{}.{}", self.prefix, self.number, v)
                }
            }
            None => {
                if self.ensembl_style {
                    format!("{}{}", self.prefix, self.number)
                } else {
                    format!("{}_{}", self.prefix, self.number)
                }
            }
        }
    }
}

impl fmt::Display for Accession {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.full())
    }
}

/// Location and edit combined
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct LocEdit<L, E> {
    pub location: L,
    pub edit: Mu<E>,
}

impl<L, E> LocEdit<L, E> {
    pub fn new(location: L, edit: E) -> Self {
        Self {
            location,
            edit: Mu::Certain(edit),
        }
    }

    /// Create a predicted/uncertain loc_edit (edit wrapped in parentheses)
    /// Used for patterns like c.(9740C>A) where the variant is predicted
    pub fn new_predicted(location: L, edit: E) -> Self {
        Self {
            location,
            edit: Mu::Uncertain(edit),
        }
    }

    pub fn with_uncertainty(location: L, edit: Mu<E>) -> Self {
        Self { location, edit }
    }
}

impl<L: fmt::Display, E: fmt::Display> fmt::Display for LocEdit<L, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.location, self.edit)
    }
}

/// Main HGVS variant enum
///
/// The relationship/separator kind among the sub-descriptions of a single
/// allele bracket. Most variants are true *phase* (cis/trans/mosaic/chimeric),
/// but the enum also carries non-phase relationships expressed with the same
/// `[...]` machinery: `AndOr` (`^`) and `Products` (`,`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AllelePhase {
    /// Same chromosome (cis): [var1;var2]
    Cis,
    /// Different chromosomes (trans): [var1];[var2]
    Trans,
    /// Phase unknown
    Unknown,
    /// Mosaic - variant present in some cells: var1/var2
    /// (somatic variation within an individual)
    Mosaic,
    /// Chimeric - different cell lines: var1//var2
    /// (two distinct cell populations in one individual)
    Chimeric,
    /// And/or - alternatives joined by `^`: var1^var2
    /// ("variant A and/or variant B"; HGVS general.md)
    AndOr,
    /// Products - different transcripts/proteins derived from one allele,
    /// joined by `,`: `r.[a,b]` (HGVS general.md, RNA/splicing.md,
    /// RNA/alleles.md). Not a phase: the members are downstream products of a
    /// single DNA-level event (e.g. alternative splicing), not co-occurrence
    /// relationships. Only meaningful on the r./p. axes.
    Products,
}

impl fmt::Display for AllelePhase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AllelePhase::Cis => write!(f, "cis"),
            AllelePhase::Trans => write!(f, "trans"),
            AllelePhase::Unknown => write!(f, "unknown"),
            AllelePhase::Mosaic => write!(f, "mosaic"),
            AllelePhase::Chimeric => write!(f, "chimeric"),
            AllelePhase::AndOr => write!(f, "and/or"),
            AllelePhase::Products => write!(f, "products"),
        }
    }
}

impl AllelePhase {
    /// Returns true if this is a mosaic phase
    pub fn is_mosaic(&self) -> bool {
        matches!(self, AllelePhase::Mosaic)
    }

    /// Returns true if this is a chimeric phase
    pub fn is_chimeric(&self) -> bool {
        matches!(self, AllelePhase::Chimeric)
    }
}

/// Allele variant containing multiple variants
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct AlleleVariant {
    /// The variants in this allele
    pub variants: Vec<HgvsVariant>,
    /// Phase relationship
    pub phase: AllelePhase,
    /// Whether the whole allele group is wrapped in an uncertainty
    /// marker `(...)` (predicted). Currently used by the and/or (`^`)
    /// group, e.g. `c.(370A>C^372C>R)`. Defaults to `false`.
    #[serde(default)]
    pub uncertain: bool,
}

impl AlleleVariant {
    /// Create a new allele variant
    pub fn new(variants: Vec<HgvsVariant>, phase: AllelePhase) -> Self {
        Self {
            variants,
            phase,
            uncertain: false,
        }
    }

    /// Create a new allele variant wrapped in an uncertainty marker
    /// `(...)` (predicted), e.g. the and/or group `c.(370A>C^372C>R)`.
    pub fn new_uncertain(variants: Vec<HgvsVariant>, phase: AllelePhase) -> Self {
        Self {
            variants,
            phase,
            uncertain: true,
        }
    }

    /// Create a cis allele (same chromosome)
    pub fn cis(variants: Vec<HgvsVariant>) -> Self {
        Self::new(variants, AllelePhase::Cis)
    }

    /// Create a trans allele (different chromosomes)
    pub fn trans(variants: Vec<HgvsVariant>) -> Self {
        Self::new(variants, AllelePhase::Trans)
    }

    /// Create a mosaic (somatic variation, some cells have variant)
    pub fn mosaic(variants: Vec<HgvsVariant>) -> Self {
        Self::new(variants, AllelePhase::Mosaic)
    }

    /// Create a chimeric (two distinct cell populations)
    pub fn chimeric(variants: Vec<HgvsVariant>) -> Self {
        Self::new(variants, AllelePhase::Chimeric)
    }

    /// Create an unknown phase allele
    pub fn unknown_phase(variants: Vec<HgvsVariant>) -> Self {
        Self::new(variants, AllelePhase::Unknown)
    }

    /// Detect whether the allele contains a self-cancelling (`del` + `dup`)
    /// pair that the HGVS spec disallows ("descriptions removing part of a
    /// reference sequence replacing it with part of the same sequence are
    /// not allowed", `recommendations/general.md` line 47).
    ///
    /// Returns `Some((idx_del, idx_dup))` for the first offending pair found
    /// (deterministic by ascending pair index), or `None` otherwise.
    ///
    /// Detection covers any concrete `del` / `dup` pair on a definite
    /// position range on the same accession: bare integer positions,
    /// 5'-UTR (`c.-N`), 3'-UTR (`c.*N`), and intronic offsets
    /// (`c.X+N` / `c.X-N`) on all axes that carry those flags
    /// (`c.`, `n.`, `r.`). Positions are compared on a canonical
    /// `(region, base, offset)` lex-order so cross-region pairs are
    /// disjoint by construction and intronic offsets sort correctly
    /// within their anchor.
    ///
    /// Variants with uncertain boundaries (`?` markers, unbounded
    /// `<pos>?` endpoints) are still skipped — uncertainty makes
    /// overlap undecidable.
    ///
    /// A single accession can host multiple HGVS coordinate axes
    /// (e.g. `NM_*` carries both `c.` and `n.`/`r.`), so the detector
    /// also requires both variants to share an axis before comparing
    /// canonical ranges. A `c.` del and an `n.` dup on the same
    /// accession with the same numeric range are not reported.
    pub fn detect_self_cancelling_pair(variants: &[HgvsVariant]) -> Option<(usize, usize)> {
        // Indexed loops are intentional: we need both `i` and `j` to return
        // the pair indices to the caller.
        #[allow(clippy::needless_range_loop)]
        for i in 0..variants.len() {
            let Some((kind_i, axis_i, range_i, acc_i)) = self_cancelling_descriptor(&variants[i])
            else {
                continue;
            };
            for j in (i + 1)..variants.len() {
                let Some((kind_j, axis_j, range_j, acc_j)) =
                    self_cancelling_descriptor(&variants[j])
                else {
                    continue;
                };
                if acc_i.full() != acc_j.full() {
                    continue;
                }
                // Gate on coordinate-axis identity: a single accession
                // can carry multiple coordinate systems (e.g. NM_*
                // supports both `c.` and `n.`/`r.`), and numerically
                // identical ranges across axes describe disjoint
                // continua. Mixing them would yield false-positive E3006.
                if axis_i != axis_j {
                    continue;
                }
                let (del_idx, dup_idx) = match (kind_i, kind_j) {
                    (SelfCancellingEditKind::Del, SelfCancellingEditKind::Dup) => (i, j),
                    (SelfCancellingEditKind::Dup, SelfCancellingEditKind::Del) => (j, i),
                    _ => continue,
                };
                if ranges_overlap(range_i, range_j) {
                    return Some((del_idx, dup_idx));
                }
            }
        }
        None
    }

    /// Build a new allele after running the spec-mandated self-cancelling
    /// check. Returns `FerroError::Parse` with
    /// `ErrorCode::SelfCancellingAllele` if the allele violates HGVS
    /// `recommendations/general.md` line 47.
    ///
    /// `source_span` is the byte range of the offending allele in the
    /// caller's source string. Pass it through whenever it's known so that
    /// downstream tooling (LSP, structured-error consumers, the web
    /// service) can underline the precise location; pass `None` for
    /// programmatically constructed alleles with no source backing. When
    /// `Some(span)` is provided, the returned error's `pos` is set to
    /// `span.start`.
    pub fn try_new_validated(
        variants: Vec<HgvsVariant>,
        phase: AllelePhase,
        source_span: Option<crate::error::SourceSpan>,
    ) -> Result<Self, crate::error::FerroError> {
        if let Some((dl, du)) = Self::detect_self_cancelling_pair(&variants) {
            return Err(build_self_cancelling_error(dl, du, source_span, None));
        }
        Ok(Self::new(variants, phase))
    }
}

/// Build the E3006 `SelfCancellingAllele` error with the supplied source
/// span. Shared by the parser-side validator and the library-level
/// `AlleleVariant::try_new_validated` so the diagnostic payload stays in
/// lock-step. `source` is the surrounding input string, attached to the
/// `Diagnostic` for caret rendering when available.
pub(crate) fn build_self_cancelling_error(
    del_idx: usize,
    dup_idx: usize,
    source_span: Option<crate::error::SourceSpan>,
    source: Option<&str>,
) -> crate::error::FerroError {
    let pos = source_span.as_ref().map(|s| s.start).unwrap_or(0);
    let mut diagnostic = crate::error::Diagnostic::new()
        .with_code(crate::error::ErrorCode::SelfCancellingAllele)
        .with_hint(
            "HGVS does not allow describing both a deletion and a duplication \
             of overlapping reference positions in the same allele \
             (recommendations/general.md line 47); drop one edit or rewrite \
             as a single delins",
        );
    if let Some(span) = source_span {
        diagnostic = diagnostic.with_span(span);
    }
    if let Some(src) = source {
        diagnostic = diagnostic.with_source(src);
    }
    crate::error::FerroError::Parse {
        pos,
        msg: format!(
            "Self-cancelling allele: variants at index {} (del) and {} (dup) describe \
             overlapping reference positions",
            del_idx, dup_idx
        ),
        diagnostic: Some(Box::new(diagnostic)),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SelfCancellingEditKind {
    Del,
    Dup,
}

/// Coordinate axis of an `HgvsVariant` arm — used by the
/// self-cancelling detector to refuse cross-axis comparisons even when
/// two variants share an accession.
///
/// A single accession can host multiple HGVS coordinate systems
/// (`NM_*` carries `c.` plus `n.`/`r.`; `NC_*` can carry both `g.` and
/// `m.` on the same circular contig in some assemblies), and a numeric
/// range like `100_150` denotes a different stretch of sequence on each
/// axis. Gating on `acc.full()` alone is therefore not enough to
/// suppress false-positive E3006 verdicts; we also require
/// `axis_i == axis_j` before testing range overlap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SelfCancellingAxis {
    Cds,
    Genome,
    Tx,
    Rna,
    Mt,
    Circular,
}

/// Canonical comparable representation of a single position in
/// HGVS-style "region + base + offset" coordinates.
///
/// Lex-ordered as `(region, base, offset)`:
///
/// - `region` is `false` for CDS-proper (`c.123`), 5'UTR (`c.-3`), and
///   any axis without a `*` / `downstream` flag (genome, mito,
///   circular). It is `true` for 3'UTR (`c.*N`, `r.*N`) and tx
///   downstream positions (`n.*N`). Because `false < true`, every
///   CDS-or-5'UTR endpoint sorts before every 3'UTR endpoint, which
///   matches the HGVS coordinate continuum.
/// - `base` carries the signed numeric value as written. 5'UTR
///   positions are negative (`c.-3` → `base = -3`); CDS-proper and
///   3'UTR positions are positive (`c.123` → 123, `c.*5` → 5).
/// - `offset` is the intronic offset (`c.100+5` → 5, `c.100-5` → -5),
///   or 0 when not set. Two positions with the same anchor sort by
///   offset.
///
/// The genome/mito axes never carry a region flag and rarely an
/// offset, but we include both for uniformity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct SelfCancellingPoint {
    region: bool,
    base: i64,
    offset: i64,
}

impl SelfCancellingPoint {
    const fn new(region: bool, base: i64, offset: i64) -> Self {
        Self {
            region,
            base,
            offset,
        }
    }
}

type SelfCancellingRange = (SelfCancellingPoint, SelfCancellingPoint);

/// Extract `(edit-kind, axis, [start, end] canonical range, accession)`
/// if the variant is a candidate for self-cancelling-pair detection.
/// Returns `None` for any variant that is not a concrete `del` / `dup`
/// on a definite (non-uncertain) position range.
///
/// The axis is returned alongside the range because two variants on the
/// same accession but different coordinate systems (`c.` vs `n.`, etc.)
/// describe disjoint sequence continua and must not be compared.
fn self_cancelling_descriptor(
    v: &HgvsVariant,
) -> Option<(
    SelfCancellingEditKind,
    SelfCancellingAxis,
    SelfCancellingRange,
    &Accession,
)> {
    use crate::hgvs::edit::NaEdit;
    let (edit, axis, range, accession) = match v {
        HgvsVariant::Cds(inner) => (
            inner.loc_edit.edit.inner()?,
            SelfCancellingAxis::Cds,
            cds_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        HgvsVariant::Genome(inner) => (
            inner.loc_edit.edit.inner()?,
            SelfCancellingAxis::Genome,
            genome_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        HgvsVariant::Tx(inner) => (
            inner.loc_edit.edit.inner()?,
            SelfCancellingAxis::Tx,
            tx_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        HgvsVariant::Rna(inner) => (
            inner.loc_edit.edit.inner()?,
            SelfCancellingAxis::Rna,
            rna_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        HgvsVariant::Mt(inner) => (
            inner.loc_edit.edit.inner()?,
            SelfCancellingAxis::Mt,
            genome_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        HgvsVariant::Circular(inner) => (
            inner.loc_edit.edit.inner()?,
            SelfCancellingAxis::Circular,
            genome_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        _ => return None,
    };
    let kind = match edit {
        NaEdit::Deletion { .. } => SelfCancellingEditKind::Del,
        NaEdit::Duplication { .. } => SelfCancellingEditKind::Dup,
        _ => return None,
    };
    Some((kind, axis, range, accession))
}

/// Return `Some(offset)` if `raw` is a concrete intronic offset, or
/// `None` if it is one of the parser's "unknown offset" sentinels
/// (`+?` → `i64::MAX`, `-?` → `i64::MIN`; see
/// `parser::position::OFFSET_UNKNOWN_{POSITIVE,NEGATIVE}` and
/// `location::GENOME_OFFSET_UNKNOWN_{POSITIVE,NEGATIVE}`). The
/// self-cancelling detector skips any endpoint with an unknown offset
/// because overlap is undecidable.
fn certain_offset(raw: Option<i64>) -> Option<i64> {
    use crate::hgvs::location::{GENOME_OFFSET_UNKNOWN_NEGATIVE, GENOME_OFFSET_UNKNOWN_POSITIVE};
    use crate::hgvs::parser::position::{OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE};
    match raw {
        None => Some(0),
        Some(v)
            if v == OFFSET_UNKNOWN_POSITIVE
                || v == OFFSET_UNKNOWN_NEGATIVE
                || v == GENOME_OFFSET_UNKNOWN_POSITIVE
                || v == GENOME_OFFSET_UNKNOWN_NEGATIVE =>
        {
            None
        }
        Some(v) => Some(v),
    }
}

fn cds_simple_range(interval: &crate::hgvs::interval::CdsInterval) -> Option<SelfCancellingRange> {
    let s = interval.start.inner()?;
    let e = interval.end.inner()?;
    // `c.?` carries `base == CDS_BASE_UNKNOWN` even inside a `Mu::Certain`
    // wrapper (so `inner()` returns Some). Filter explicitly: overlap on
    // an unknown base is undecidable.
    if s.is_unknown() || e.is_unknown() {
        return None;
    }
    let s_off = certain_offset(s.offset)?;
    let e_off = certain_offset(e.offset)?;
    Some((
        SelfCancellingPoint::new(s.utr3, s.base, s_off),
        SelfCancellingPoint::new(e.utr3, e.base, e_off),
    ))
}

fn genome_simple_range(
    interval: &crate::hgvs::interval::GenomeInterval,
) -> Option<SelfCancellingRange> {
    let s = interval.start.inner()?;
    let e = interval.end.inner()?;
    let s_off = certain_offset(s.offset)?;
    let e_off = certain_offset(e.offset)?;
    Some((
        SelfCancellingPoint::new(false, s.base as i64, s_off),
        SelfCancellingPoint::new(false, e.base as i64, e_off),
    ))
}

fn tx_simple_range(interval: &crate::hgvs::interval::TxInterval) -> Option<SelfCancellingRange> {
    let s = interval.start.inner()?;
    let e = interval.end.inner()?;
    let s_off = certain_offset(s.offset)?;
    let e_off = certain_offset(e.offset)?;
    Some((
        SelfCancellingPoint::new(s.downstream, s.base, s_off),
        SelfCancellingPoint::new(e.downstream, e.base, e_off),
    ))
}

fn rna_simple_range(interval: &crate::hgvs::interval::RnaInterval) -> Option<SelfCancellingRange> {
    let s = interval.start.inner()?;
    let e = interval.end.inner()?;
    let s_off = certain_offset(s.offset)?;
    let e_off = certain_offset(e.offset)?;
    Some((
        SelfCancellingPoint::new(s.utr3, s.base, s_off),
        SelfCancellingPoint::new(e.utr3, e.base, e_off),
    ))
}

/// Returns `true` when ranges `a` and `b` overlap.
///
/// **Precondition**: both ranges are in canonical orientation
/// (`a.0 <= a.1`, `b.0 <= b.1`). The preprocessor's W4001 Phase 15a
/// (`correct_swapped_positions`) enforces this for parsed input;
/// programmatically constructed AlleleVariants must satisfy it too or
/// the verdict may be a false negative.
fn ranges_overlap(a: SelfCancellingRange, b: SelfCancellingRange) -> bool {
    let lo = a.0.max(b.0);
    let hi = a.1.min(b.1);
    lo <= hi
}

/// Set the count of a single-count `Repeat` edit on `v` to `new_count`,
/// returning `true` on success. Used to build the sibling alleles of a
/// shared-position repeat trans allele (`p.Ala2[10];[11]`): clone the first
/// variant and swap in each subsequent allele's count. Variants that are
/// not single-count repeats are left unchanged and return `false`.
pub(crate) fn set_repeat_count(
    v: &mut HgvsVariant,
    new_count: crate::hgvs::edit::RepeatCount,
) -> bool {
    fn na(edit: &mut Mu<NaEdit>, c: crate::hgvs::edit::RepeatCount) -> bool {
        match edit {
            Mu::Certain(NaEdit::Repeat {
                count,
                additional_counts,
                trailing,
                ..
            })
            | Mu::Uncertain(NaEdit::Repeat {
                count,
                additional_counts,
                trailing,
                ..
            }) if additional_counts.is_empty() && trailing.is_none() => {
                *count = c;
                true
            }
            _ => false,
        }
    }
    match v {
        HgvsVariant::Genome(x) => na(&mut x.loc_edit.edit, new_count),
        HgvsVariant::Cds(x) => na(&mut x.loc_edit.edit, new_count),
        HgvsVariant::Tx(x) => na(&mut x.loc_edit.edit, new_count),
        HgvsVariant::Rna(x) => na(&mut x.loc_edit.edit, new_count),
        HgvsVariant::Mt(x) => na(&mut x.loc_edit.edit, new_count),
        HgvsVariant::Circular(x) => na(&mut x.loc_edit.edit, new_count),
        HgvsVariant::Protein(x) => match &mut x.loc_edit.edit {
            Mu::Certain(ProteinEdit::Repeat { count, .. })
            | Mu::Uncertain(ProteinEdit::Repeat { count, .. }) => {
                *count = new_count;
                true
            }
            _ => false,
        },
        _ => false,
    }
}

/// The count of a single-count `Repeat` edit on `v`, if any. Used by the
/// shared-position repeat trans Display.
fn repeat_count_of(v: &HgvsVariant) -> Option<&crate::hgvs::edit::RepeatCount> {
    fn na(edit: &Mu<NaEdit>) -> Option<&crate::hgvs::edit::RepeatCount> {
        match edit.inner() {
            Some(NaEdit::Repeat {
                count,
                additional_counts,
                trailing,
                ..
            }) if additional_counts.is_empty() && trailing.is_none() => Some(count),
            _ => None,
        }
    }
    match v {
        HgvsVariant::Genome(x) => na(&x.loc_edit.edit),
        HgvsVariant::Cds(x) => na(&x.loc_edit.edit),
        HgvsVariant::Tx(x) => na(&x.loc_edit.edit),
        HgvsVariant::Rna(x) => na(&x.loc_edit.edit),
        HgvsVariant::Mt(x) => na(&x.loc_edit.edit),
        HgvsVariant::Circular(x) => na(&x.loc_edit.edit),
        HgvsVariant::Protein(x) => match x.loc_edit.edit.inner() {
            Some(ProteinEdit::Repeat { count, .. }) => Some(count),
            _ => None,
        },
        _ => None,
    }
}

/// True iff `variants` is a shared-position repeat trans allele: 2+ members
/// that are single-count repeats sharing accession, coord type, position,
/// and repeat unit — differing only in count. Such an allele renders in the
/// compact `<pos><unit>[c1];[c2]` form (HGVS repeated.md) rather than the
/// generic `[m1];[m2]`.
fn is_shared_repeat_trans(variants: &[HgvsVariant]) -> bool {
    if variants.len() < 2 || !HgvsVariant::all_share_accession_and_type(variants) {
        return false;
    }
    // The compact `<pos><unit>[c1];[c2]` form is the spec canonical for RNA
    // and protein (`recommendations/{protein,DNA}/repeated.md`); the DNA
    // axes (g./c./m./o.) spell each allele out in full (`[member];[member]`,
    // DNA/repeated.md:33). Gate the compact rendering to RNA + protein so
    // DNA repeat trans alleles keep their explicit spec form.
    if !variants
        .iter()
        .all(|v| matches!(v, HgvsVariant::Rna(_) | HgvsVariant::Protein(_)))
    {
        return false;
    }
    // Key = the full Display with the count normalized out; equal keys mean
    // the members differ only in their repeat count.
    let key = |v: &HgvsVariant| -> Option<String> {
        let mut probe = v.clone();
        set_repeat_count(&mut probe, crate::hgvs::edit::RepeatCount::Exact(0))
            .then(|| format!("{}", probe))
    };
    match key(&variants[0]) {
        Some(first) => variants
            .iter()
            .all(|v| key(v).as_deref() == Some(first.as_str())),
        None => false,
    }
}

/// Write the `ACC:<type>.` prefix for compact allele form.
///
/// Per HGVS spec, bracketed `?` is only valid as a whole-allele marker (`[?]`),
/// never mixed with concrete edits inside the same bracket. Compact form is
/// suppressed when any sub-variant carries the per-variant `?` (e.g. `c.?`).
fn use_compact_form(variants: &[HgvsVariant]) -> bool {
    !variants.is_empty()
        && HgvsVariant::all_share_accession_and_type(variants)
        && !variants.iter().any(HgvsVariant::is_loc_edit_unknown)
}

/// For a `Trans` allele, return the first concrete leaf sub-variant if the
/// compact form `ACC:c.[..];[..]` applies. Every concrete leaf — descending
/// one level into nested cis-group members — must share accession, coord
/// type, and gene selector; there must be at least one; and none may carry
/// the per-variant unknown form. `[0]`/`[?]` members contribute no leaves
/// and are permitted alongside concrete members (this is what keeps
/// `ACC:c.[X];[?]` / `;[0]` in compact form).
fn trans_compact_anchor(variants: &[HgvsVariant]) -> Option<&HgvsVariant> {
    let mut anchor: Option<&HgvsVariant> = None;
    for member in variants {
        let leaves: &[HgvsVariant] = match member {
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => &[],
            HgvsVariant::Allele(a) => &a.variants,
            other => std::slice::from_ref(other),
        };
        for leaf in leaves {
            if leaf.accession().is_none()
                || matches!(
                    leaf.variant_type(),
                    "allele" | "null" | "unknown" | "r::r" | "g::g"
                )
                || leaf.is_loc_edit_unknown()
            {
                return None;
            }
            match anchor {
                None => anchor = Some(leaf),
                Some(a) => {
                    if leaf.accession() != a.accession()
                        || leaf.variant_type() != a.variant_type()
                        || leaf.gene_symbol() != a.gene_symbol()
                    {
                        return None;
                    }
                }
            }
        }
    }
    anchor
}

/// Write one trans member wrapped in its bracket(s). A cis-group / leaf /
/// `[0]` / `[?]` member renders as `[<inner>]`. An unknown-phase sub-allele
/// member (the trans+mosaic mix `[B](;)C`) brackets only its first member and
/// appends the rest with `(;)` outside the bracket.
fn write_trans_member_bracketed(f: &mut fmt::Formatter<'_>, member: &HgvsVariant) -> fmt::Result {
    if let HgvsVariant::Allele(a) = member {
        if a.phase == AllelePhase::Unknown {
            for (i, sub) in a.variants.iter().enumerate() {
                if i == 0 {
                    write!(f, "[")?;
                    sub.fmt_loc_edit(f)?;
                    write!(f, "]")?;
                } else {
                    write!(f, "(;)")?;
                    sub.fmt_loc_edit(f)?;
                }
            }
            return Ok(());
        }
    }
    write!(f, "[")?;
    write_trans_member(f, member)?;
    write!(f, "]")
}

/// Write the bracket-inner content of one `Trans` allele member (without the
/// surrounding `[` `]`). A nested cis-group member renders its sub-variants'
/// loc-edits joined by `;`; `[0]`/`[?]` render their marker; a leaf renders
/// its loc-edit.
fn write_trans_member(f: &mut fmt::Formatter<'_>, member: &HgvsVariant) -> fmt::Result {
    match member {
        HgvsVariant::NullAllele => write!(f, "0"),
        HgvsVariant::UnknownAllele => write!(f, "?"),
        HgvsVariant::Allele(a) => {
            for (i, sub) in a.variants.iter().enumerate() {
                if i > 0 {
                    write!(f, ";")?;
                }
                sub.fmt_loc_edit(f)?;
            }
            Ok(())
        }
        other => other.fmt_loc_edit(f),
    }
}

/// The first leaf (non-allele) sub-variant that carries an accession,
/// descending through nested allele members. Used to derive the shared
/// accession + coord-type prefix of an and/or-of-alleles or to inherit it
/// onto bare bracketed operands.
pub(crate) fn first_leaf_with_accession(v: &HgvsVariant) -> Option<&HgvsVariant> {
    match v {
        HgvsVariant::Allele(a) => a.variants.iter().find_map(first_leaf_with_accession),
        HgvsVariant::NullAllele | HgvsVariant::UnknownAllele | HgvsVariant::RnaFusion(_) => None,
        other => {
            if other.accession().is_some() {
                Some(other)
            } else {
                None
            }
        }
    }
}

/// True iff every operand of an and/or allele is itself an `Allele` and all
/// their leaves share accession + coord type — the `ACC:p.[..]^[..]` form
/// where the shared prefix is written once.
fn andor_alleles_share_accession(variants: &[HgvsVariant]) -> bool {
    let mut anchor: Option<&HgvsVariant> = None;
    for member in variants {
        let HgvsVariant::Allele(_) = member else {
            return false;
        };
        let leaf = match first_leaf_with_accession(member) {
            Some(l) => l,
            None => return false,
        };
        match anchor {
            None => anchor = Some(leaf),
            Some(a) => {
                if leaf.accession() != a.accession() || leaf.variant_type() != a.variant_type() {
                    return false;
                }
            }
        }
    }
    anchor.is_some()
}

/// Write the bracket-inner content of one and/or allele operand: each
/// member's loc-edit, separated only between adjacent members — a
/// one-member allele emits only that member's loc-edit with no separator.
/// The separator is `(;)` for an unknown-phase group or `;` for cis.
fn write_andor_allele_inner(f: &mut fmt::Formatter<'_>, allele: &AlleleVariant) -> fmt::Result {
    let sep = if allele.phase == AllelePhase::Unknown {
        "(;)"
    } else {
        ";"
    };
    for (i, member) in allele.variants.iter().enumerate() {
        if i > 0 {
            write!(f, "{}", sep)?;
        }
        member.fmt_loc_edit(f)?;
    }
    Ok(())
}

fn write_compact_prefix(f: &mut fmt::Formatter<'_>, first: &HgvsVariant) -> fmt::Result {
    let accession = first
        .accession()
        .expect("compact form requires an accession; guarded by all_share_accession_and_type");
    write!(f, "{}", accession)?;
    if let Some(gene) = first.gene_symbol() {
        write!(f, "({})", gene)?;
    }
    write!(f, ":{}.", first.variant_type())
}

/// Write `accession(gene):` when `gene_symbol` is set, otherwise `accession:`.
///
/// Per HGVS Nomenclature, the gene-symbol selector is informational disambiguation.
/// ferro's round-trip policy is "preserve when present in input; do not synthesize
/// when absent" (#121). This helper centralizes the selector emission so all
/// per-kind `Display` impls follow the same rule.
fn write_accession_with_optional_gene(
    f: &mut fmt::Formatter<'_>,
    accession: &Accession,
    gene_symbol: Option<&str>,
) -> fmt::Result {
    write!(f, "{}", accession)?;
    if let Some(gene) = gene_symbol {
        write!(f, "({})", gene)?;
    }
    Ok(())
}

/// Format the position+edit portion of a genome loc-edit (without the
/// `accession:g.` prefix). Shared by [`GenomeVariant::fmt_loc_edit`] and
/// [`GenomeRing`]'s `Display` so a ring segment round-trips byte-identically
/// to a standalone genome variant.
///
/// Whole-entity identity (`g.=`) / unknown (`g.?`) and a positionless
/// breakpoint insertion (composite SV `[...;insG]`) skip the position; the
/// predicted form wraps position+edit in parens (#241); otherwise the plain
/// `LocEdit` Display is used.
fn fmt_genome_loc_edit(
    f: &mut fmt::Formatter<'_>,
    loc_edit: &LocEdit<GenomeInterval, NaEdit>,
) -> fmt::Result {
    // For whole-entity identity (g.=) or unknown (g.?), or a positionless
    // breakpoint insertion (composite SV `[...;insG]`), skip the position.
    if let Some(edit) = loc_edit.edit.inner() {
        if edit.is_whole_entity() || edit.is_positionless() {
            return write!(f, "{}", loc_edit.edit);
        }
    }
    // Predicted form: wrap position+edit in parens. #241.
    if loc_edit.edit.is_uncertain() {
        if let Some(edit) = loc_edit.edit.inner() {
            return write!(f, "({}{})", loc_edit.location, edit);
        }
    }
    write!(f, "{}", loc_edit)
}

/// True iff `variant` carries a position-bound (not whole-entity)
/// identity edit (`<pos>=` or `<pos><ref>=`). This is the LHS shape
/// of the HGVS spec compact mosaic / chimeric form.
fn is_position_identity_lhs(variant: &HgvsVariant) -> bool {
    let edit_is_pos_identity = |edit: &Mu<NaEdit>| {
        matches!(
            edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: false,
                ..
            })
        )
    };
    match variant {
        HgvsVariant::Genome(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Cds(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Tx(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Rna(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Mt(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Circular(v) => edit_is_pos_identity(&v.loc_edit.edit),
        // Protein position-bound identity (`p.Trp24=`) — the LHS of the
        // protein `=/` compact mosaic form.
        HgvsVariant::Protein(v) => matches!(
            v.loc_edit.edit.inner(),
            Some(crate::hgvs::edit::ProteinEdit::Identity {
                whole_protein: false,
                ..
            })
        ),
        _ => false,
    }
}

/// True iff `variant` carries a whole-entity identity edit (the
/// canonical shape produced by `create_identity_variant_from`, used
/// in the `var/=` shorthand). On Display this should render as bare
/// `=` rather than the full synthetic-position form.
fn is_whole_entity_identity_rhs(variant: &HgvsVariant) -> bool {
    let edit_is_whole = |edit: &Mu<NaEdit>| {
        matches!(
            edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: true,
                ..
            })
        )
    };
    match variant {
        HgvsVariant::Genome(v) => edit_is_whole(&v.loc_edit.edit),
        HgvsVariant::Cds(v) => edit_is_whole(&v.loc_edit.edit),
        HgvsVariant::Tx(v) => edit_is_whole(&v.loc_edit.edit),
        HgvsVariant::Rna(v) => edit_is_whole(&v.loc_edit.edit),
        HgvsVariant::Mt(v) => edit_is_whole(&v.loc_edit.edit),
        HgvsVariant::Circular(v) => edit_is_whole(&v.loc_edit.edit),
        _ => false,
    }
}

/// True iff `a` and `b` carry the same coord-system arm AND their
/// `loc_edit.location` values compare equal. Used to decide whether
/// the spec compact mosaic Display form applies to a 2-variant Allele.
fn intervals_match_for_compact_mosaic(a: &HgvsVariant, b: &HgvsVariant) -> bool {
    match (a, b) {
        (HgvsVariant::Genome(x), HgvsVariant::Genome(y)) => {
            x.loc_edit.location == y.loc_edit.location
        }
        (HgvsVariant::Cds(x), HgvsVariant::Cds(y)) => x.loc_edit.location == y.loc_edit.location,
        (HgvsVariant::Tx(x), HgvsVariant::Tx(y)) => x.loc_edit.location == y.loc_edit.location,
        (HgvsVariant::Rna(x), HgvsVariant::Rna(y)) => x.loc_edit.location == y.loc_edit.location,
        (HgvsVariant::Mt(x), HgvsVariant::Mt(y)) => x.loc_edit.location == y.loc_edit.location,
        (HgvsVariant::Circular(x), HgvsVariant::Circular(y)) => {
            x.loc_edit.location == y.loc_edit.location
        }
        (HgvsVariant::Protein(x), HgvsVariant::Protein(y)) => {
            x.loc_edit.location == y.loc_edit.location
        }
        _ => false,
    }
}

/// Shared Display path for `Mosaic` (`/`) and `Chimeric` (`//`) phases.
///
/// Three branches in priority order (mutually exclusive on each input):
///
/// 1. **Spec compact mosaic / chimeric** — LHS is a position-bound `=`
///    identity edit, RHS shares accession + coord type + interval.
///    Emits `<lhs-full>=/<rhs-bare-edit>`. Per
///    `recommendations/DNA/{substitution,deletion,duplication}.md`.
///
/// 2. **`var/=` shorthand cleanup** — LHS is a non-identity variant,
///    RHS is a whole-entity identity edit (the synthetic shape
///    produced by `create_identity_variant_from`), and the two share
///    accession + coord type. Emits `<lhs-full>/=` instead of
///    expanding the synthetic identity to its `<acc>:<type>.1=` form.
///
/// 3. **Long form** — fallback per-variant join.
fn write_mosaic_or_chimeric(
    f: &mut fmt::Formatter<'_>,
    variants: &[HgvsVariant],
    sep: &str,
) -> fmt::Result {
    if variants.len() == 2 && use_compact_form(variants) {
        let lhs = &variants[0];
        let rhs = &variants[1];

        // Branch 1: spec compact (LHS pos-identity, intervals match).
        if is_position_identity_lhs(lhs) && intervals_match_for_compact_mosaic(lhs, rhs) {
            write!(f, "{}", lhs)?;
            write!(f, "{}", sep)?;
            // Emit the RHS edit alone (no accession, no type prefix, no
            // position) — that is the spec compact RHS shape.
            return write_bare_edit(f, rhs);
        }

        // Branch 2: var/= cleanup. LHS must NOT be position-identity
        // (otherwise branch 1 owns it); RHS must be whole-entity Identity.
        if !is_position_identity_lhs(lhs) && is_whole_entity_identity_rhs(rhs) {
            write!(f, "{}", lhs)?;
            write!(f, "{}", sep)?;
            write!(f, "=")?;
            return Ok(());
        }
    }

    // Branch 3: long-form join.
    for (i, v) in variants.iter().enumerate() {
        if i > 0 {
            write!(f, "{}", sep)?;
        }
        write!(f, "{}", v)?;
    }
    Ok(())
}

/// Write the `NaEdit` portion of `variant` (no accession, no type
/// prefix, no position). Used by `write_mosaic_or_chimeric` for the
/// spec compact RHS.
///
/// Variants that are not NA-coord arms are written with their full
/// `Display` (callers should not invoke this on those — guarded by
/// `use_compact_form`).
fn write_bare_edit(f: &mut fmt::Formatter<'_>, variant: &HgvsVariant) -> fmt::Result {
    match variant {
        HgvsVariant::Genome(v) => write!(f, "{}", v.loc_edit.edit),
        HgvsVariant::Cds(v) => write!(f, "{}", v.loc_edit.edit),
        HgvsVariant::Tx(v) => write!(f, "{}", v.loc_edit.edit),
        // RNA renders nucleotides in lowercase per the HGVS spec. The bare
        // edit must go through `to_rna_string()` — the plain `NaEdit`
        // `Display` is uppercase and would corrupt the case here.
        HgvsVariant::Rna(v) => match v.loc_edit.edit.inner() {
            Some(edit) => write!(f, "{}", edit.to_rna_string()),
            None => write!(f, "{}", v.loc_edit.edit),
        },
        HgvsVariant::Mt(v) => write!(f, "{}", v.loc_edit.edit),
        HgvsVariant::Circular(v) => write!(f, "{}", v.loc_edit.edit),
        // Protein substitution Display renders just the alternative AA, and
        // del/dup render their keyword — exactly the bare-edit RHS shape.
        HgvsVariant::Protein(v) => write!(f, "{}", v.loc_edit.edit),
        _ => write!(f, "{}", variant),
    }
}

impl fmt::Display for AlleleVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Per HGVS spec, `[var]` denotes a multi-variant allele; a singleton
        // is just the inner variant. (Mosaic/chimeric already emit bare via
        // their separator-only loops; this unifies cis/trans/unknown.)
        // A valid predicted cis allele always has ≥2 members — enforced by the
        // ';' guard in `find_predicted_cis_bracket` — so this singleton
        // short-circuit never fires for `[(a;b)]` inputs.  The normalizer
        // also guards `allele.uncertain` before collapsing a cis allele to a
        // singleton, keeping this invariant intact through the normalize path.
        if self.variants.len() == 1 {
            return self.variants[0].fmt(f);
        }
        match self.phase {
            AllelePhase::Cis => {
                // A predicted/uncertain cis allele wraps the bracket content in
                // `(...)`: `[(a;b)]` (recommendations/uncertain.md).
                let (lp, rp) = if self.uncertain { ("(", ")") } else { ("", "") };
                if use_compact_form(&self.variants) {
                    // Compact form: ACC:g.[edit1;edit2] / ACC:p.[(edit1;edit2)]
                    write_compact_prefix(f, &self.variants[0])?;
                    write!(f, "[{}", lp)?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, ";")?;
                        }
                        v.fmt_loc_edit(f)?;
                    }
                    write!(f, "{}]", rp)
                } else {
                    // Expanded form: [ACC:g.edit1;ACC:g.edit2]
                    write!(f, "[{}", lp)?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, ";")?;
                        }
                        write!(f, "{}", v)?;
                    }
                    write!(f, "{}]", rp)
                }
            }
            AllelePhase::Products => {
                // `[a,b]` — different transcripts/proteins from one allele.
                // Mirrors the cis compact/expanded split (members may carry
                // distinct transcript accessions, so do not assume shared),
                // joining with `,`. A predicted group wraps the members in
                // `(...)` inside the bracket: `[(a,b)]`.
                let (lp, rp) = if self.uncertain { ("(", ")") } else { ("", "") };
                if use_compact_form(&self.variants) {
                    write_compact_prefix(f, &self.variants[0])?;
                    write!(f, "[{}", lp)?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, ",")?;
                        }
                        v.fmt_loc_edit(f)?;
                    }
                    write!(f, "{}]", rp)
                } else {
                    write!(f, "[{}", lp)?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, ",")?;
                        }
                        write!(f, "{}", v)?;
                    }
                    write!(f, "{}]", rp)
                }
            }
            AllelePhase::Trans => {
                // A trans member that is itself a nested cis group (`[a;b]`)
                // needs dedicated rendering; the existing flat-trans paths
                // (which `fmt_loc_edit` a single leaf per bracket) can't
                // express it. Only take the nested path when a nested member
                // is actually present, so flat trans Display is unchanged.
                let has_nested_group = self
                    .variants
                    .iter()
                    .any(|v| matches!(v, HgvsVariant::Allele(_)));
                if has_nested_group {
                    if let Some(anchor) = trans_compact_anchor(&self.variants) {
                        // Compact: ACC:c.[a;b];[c;d] — prefix once, then each
                        // member's bracket-inner content. An unknown-phase
                        // member renders as `[B](;)C` (mosaic tail outside the
                        // bracket).
                        write_compact_prefix(f, anchor)?;
                        for (i, v) in self.variants.iter().enumerate() {
                            if i > 0 {
                                write!(f, ";")?;
                            }
                            write_trans_member_bracketed(f, v)?;
                        }
                        Ok(())
                    } else {
                        // Mixed-reference nested trans (not a spec target):
                        // expanded, each member rendered in full.
                        for (i, v) in self.variants.iter().enumerate() {
                            if i > 0 {
                                write!(f, ";")?;
                            }
                            write_trans_member_bracketed(f, v)?;
                        }
                        Ok(())
                    }
                } else if is_shared_repeat_trans(&self.variants) {
                    // Shared-position repeat trans: `ACC:p.Ala2[10];[11]` —
                    // the position+unit is written once (first member in
                    // full), then each allele's count `;[cN]`. A shared-repeat
                    // group has no nested-cis members, so it is disjoint from
                    // the `has_nested_group` arm above.
                    write_compact_prefix(f, &self.variants[0])?;
                    self.variants[0].fmt_loc_edit(f)?;
                    for v in &self.variants[1..] {
                        write!(f, ";")?;
                        // `repeat_count_of` cannot be `None` here: every member
                        // is a single-count repeat (the same invariant
                        // `is_shared_repeat_trans` keyed on). The `if let`
                        // stays defensive — a missing count would otherwise
                        // emit a malformed bare `;`.
                        debug_assert!(
                            repeat_count_of(v).is_some(),
                            "is_shared_repeat_trans guarantees a single-count repeat per member"
                        );
                        if let Some(count) = repeat_count_of(v) {
                            write!(f, "{}", count)?;
                        }
                    }
                    Ok(())
                } else if use_compact_form(&self.variants) {
                    // Compact form: ACC:g.[edit1];[edit2]
                    write_compact_prefix(f, &self.variants[0])?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, ";")?;
                        }
                        write!(f, "[")?;
                        v.fmt_loc_edit(f)?;
                        write!(f, "]")?;
                    }
                    Ok(())
                } else {
                    // Expanded form: [ACC:g.edit1];[ACC:g.edit2]
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, ";")?;
                        }
                        write!(f, "[{}]", v)?;
                    }
                    Ok(())
                }
            }
            AllelePhase::Unknown => {
                // Unknown-phase compact form has no surrounding brackets, so a
                // singleton would print as `ACC:type.edit` and lose the allele
                // wrapper. Only use compact form when there are 2+ sub-variants.
                if self.variants.len() > 1 && use_compact_form(&self.variants) {
                    // Compact form: ACC:g.edit1(;)edit2
                    write_compact_prefix(f, &self.variants[0])?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, "(;)")?;
                        }
                        v.fmt_loc_edit(f)?;
                    }
                    Ok(())
                } else {
                    // Expanded form: [ACC:g.edit1(;)ACC:g.edit2]
                    write!(f, "[")?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, "(;)")?;
                        }
                        write!(f, "{}", v)?;
                    }
                    write!(f, "]")
                }
            }
            AllelePhase::Mosaic => write_mosaic_or_chimeric(f, &self.variants, "/"),
            AllelePhase::Chimeric => write_mosaic_or_chimeric(f, &self.variants, "//"),
            AllelePhase::AndOr => {
                if self.uncertain && use_compact_form(&self.variants) {
                    // Uncertainty-wrapped, shared accession + coord type:
                    // `ACC:c.(edit1^edit2)` — prefix once, then the
                    // `^`-joined bare edits inside the parens.
                    write_compact_prefix(f, &self.variants[0])?;
                    write!(f, "(")?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, "^")?;
                        }
                        v.fmt_loc_edit(f)?;
                    }
                    write!(f, ")")
                } else if self.uncertain {
                    // Uncertainty-wrapped but mixed references: wrap the
                    // full-description join in parens — `(ACC1:c.A^ACC2:c.B)`.
                    write!(f, "(")?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, "^")?;
                        }
                        write!(f, "{}", v)?;
                    }
                    write!(f, ")")
                } else if andor_alleles_share_accession(&self.variants) {
                    // And/or between bracketed alleles sharing an accession:
                    // `ACC:p.[a(;)b]^[c]` — prefix once, then each operand as
                    // a bracketed allele.
                    let anchor = first_leaf_with_accession(&self.variants[0])
                        .expect("andor_alleles_share_accession guarantees a leaf");
                    write_compact_prefix(f, anchor)?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, "^")?;
                        }
                        write!(f, "[")?;
                        if let HgvsVariant::Allele(a) = v {
                            write_andor_allele_inner(f, a)?;
                        } else {
                            unreachable!(
                                "andor_alleles_share_accession guarantees all operands are Alleles"
                            )
                        }
                        write!(f, "]")?;
                    }
                    Ok(())
                } else {
                    // And/or `^` join of full descriptions (possibly across
                    // references): render each in full — `ACC1:c.A^ACC2:c.B`.
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, "^")?;
                        }
                        write!(f, "{}", v)?;
                    }
                    Ok(())
                }
            }
        }
    }
}

/// Represents all types of HGVS variants:
/// - g. (genomic)
/// - c. (coding DNA)
/// - n. (non-coding transcript)
/// - r. (RNA)
/// - p. (protein)
/// - m. (mitochondrial)
/// - o. (circular DNA) - SVD-WG006
/// - RNA fusion (::) - SVD-WG007
/// - Ring chromosome (::) — HGVS DNA/complex.md:127
/// - sup (supernumerary chromosome marker) — HGVS DNA/complex.md:33-34
/// - Allele (compound variants)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum HgvsVariant {
    /// Genomic variant (g.)
    Genome(GenomeVariant),
    /// Coding DNA variant (c.)
    Cds(CdsVariant),
    /// Non-coding transcript variant (n.)
    Tx(TxVariant),
    /// RNA variant (r.)
    Rna(RnaVariant),
    /// Protein variant (p.)
    Protein(ProteinVariant),
    /// Mitochondrial variant (m.)
    Mt(MtVariant),
    /// Circular DNA variant (o.) - SVD-WG006
    Circular(CircularVariant),
    /// RNA Fusion variant (::) - SVD-WG007
    RnaFusion(RnaFusionVariant),
    /// Ring chromosome: two or more genome edits on a single accession joined
    /// by `::` break junctions (HGVS DNA/complex.md:127/130).
    GenomeRing(GenomeRing),
    /// Supernumerary chromosome marker `sup` (HGVS DNA/complex.md:33-34, the
    /// ISCN named extension for a supernumerary chromosome). Wraps the inner
    /// genome description and renders it with a trailing `sup`.
    ///
    /// Two inner shapes occur in the spec:
    /// - a supernumerary whole chromosome (`pter_qtersup`,
    ///   DNA/duplication.md:117), where the inner is a [`HgvsVariant::Genome`]
    ///   with a `pter_qter` interval and no edit; and
    /// - a supernumerary ring chromosome (`[<ring>]sup`,
    ///   DNA/complex.md:161), where the inner is a [`HgvsVariant::GenomeRing`].
    ///   Unlike a standalone ring (which renders without brackets), the
    ///   supernumerary ring is bracketed (`g.[<ring>]sup`) so the marker
    ///   applies to the whole ring; `Display` adds the brackets for this case.
    Supernumerary(Box<HgvsVariant>),
    /// Allele/compound variant
    Allele(AlleleVariant),
    /// Null allele marker [0] - indicates no variant on this chromosome (hemizygous)
    NullAllele,
    /// Unknown allele marker [?] - indicates unknown variant on second chromosome
    UnknownAllele,
}

impl HgvsVariant {
    /// Get the accession for this variant, if available.
    ///
    /// For allele variants, returns the accession of the first variant.
    /// For RNA fusion variants, returns the 5' partner accession.
    /// Returns `None` for NullAllele, UnknownAllele, or empty allele variants.
    pub fn accession(&self) -> Option<&Accession> {
        match self {
            HgvsVariant::Genome(v) => Some(&v.accession),
            HgvsVariant::Cds(v) => Some(&v.accession),
            HgvsVariant::Tx(v) => Some(&v.accession),
            HgvsVariant::Rna(v) => Some(&v.accession),
            HgvsVariant::Protein(v) => Some(&v.accession),
            HgvsVariant::Mt(v) => Some(&v.accession),
            HgvsVariant::Circular(v) => Some(&v.accession),
            HgvsVariant::RnaFusion(v) => Some(&v.five_prime.accession),
            HgvsVariant::GenomeRing(g) => Some(&g.accession),
            HgvsVariant::Supernumerary(inner) => inner.accession(),
            HgvsVariant::Allele(a) => a.variants.first().and_then(|v| v.accession()),
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => None,
        }
    }

    /// Get the gene-symbol selector for this variant, if one was present in the
    /// original input.
    ///
    /// Mirrors `accession()` and is used by `Display` to preserve the selector
    /// on round-trip (#121). Returns `None` for variant kinds that do not carry
    /// a single top-level selector (`Allele`, `RnaFusion`, `NullAllele`,
    /// `UnknownAllele`); selectors on those nested kinds are emitted by their
    /// own `Display` impls.
    pub fn gene_symbol(&self) -> Option<&str> {
        match self {
            HgvsVariant::Genome(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Cds(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Tx(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Rna(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Protein(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Mt(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Circular(v) => v.gene_symbol.as_deref(),
            HgvsVariant::GenomeRing(g) => g.gene_symbol.as_deref(),
            HgvsVariant::Supernumerary(inner) => inner.gene_symbol(),
            HgvsVariant::RnaFusion(_)
            | HgvsVariant::Allele(_)
            | HgvsVariant::NullAllele
            | HgvsVariant::UnknownAllele => None,
        }
    }

    /// Get the variant type string
    pub fn variant_type(&self) -> &'static str {
        match self {
            HgvsVariant::Genome(_) => "g",
            HgvsVariant::Cds(_) => "c",
            HgvsVariant::Tx(_) => "n",
            HgvsVariant::Rna(_) => "r",
            HgvsVariant::Protein(_) => "p",
            HgvsVariant::Mt(_) => "m",
            HgvsVariant::Circular(_) => "o",
            HgvsVariant::RnaFusion(_) => "r::r",
            HgvsVariant::GenomeRing(_) => "g::g",
            HgvsVariant::Supernumerary(inner) => inner.variant_type(),
            HgvsVariant::Allele(_) => "allele",
            HgvsVariant::NullAllele => "null",
            HgvsVariant::UnknownAllele => "unknown",
        }
    }

    /// Check if this is an allele (compound) variant
    pub fn is_allele(&self) -> bool {
        matches!(self, HgvsVariant::Allele(_))
    }

    /// Check if this is a null allele marker
    pub fn is_null_allele(&self) -> bool {
        matches!(self, HgvsVariant::NullAllele)
    }

    /// Check if this is an unknown allele marker
    pub fn is_unknown_allele(&self) -> bool {
        matches!(self, HgvsVariant::UnknownAllele)
    }

    /// Format just the position+edit portion (without accession and coordinate prefix).
    ///
    /// For `NM_000088.3:c.459A>G`, this writes `459A>G`. Used by `AlleleVariant::Display`
    /// to emit the spec-correct compact form (`ACC:c.[edit1;edit2]`).
    ///
    /// Variant types that have no simple loc/edit form (`Allele`, `RnaFusion`,
    /// `NullAllele`, `UnknownAllele`) fall back to their full `Display` output. Callers
    /// guard against this via `all_share_accession_and_type`, which excludes those
    /// types from the compact branch.
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HgvsVariant::Genome(v) => v.fmt_loc_edit(f),
            HgvsVariant::Cds(v) => v.fmt_loc_edit(f),
            HgvsVariant::Tx(v) => v.fmt_loc_edit(f),
            HgvsVariant::Rna(v) => v.fmt_loc_edit(f),
            HgvsVariant::Protein(v) => v.fmt_loc_edit(f),
            HgvsVariant::Mt(v) => v.fmt_loc_edit(f),
            HgvsVariant::Circular(v) => v.fmt_loc_edit(f),
            // These types have no compact form — fall back to full Display.
            HgvsVariant::RnaFusion(v) => write!(f, "{}", v),
            HgvsVariant::GenomeRing(g) => write!(f, "{}", g),
            HgvsVariant::Supernumerary(_) => write!(f, "{}", self),
            HgvsVariant::Allele(a) => write!(f, "{}", a),
            HgvsVariant::NullAllele => write!(f, "0"),
            HgvsVariant::UnknownAllele => write!(f, "?"),
        }
    }

    /// True if this variant's loc/edit portion is the per-variant unknown form
    /// (`g.?`, `c.?`, `n.?`, `r.?`, `p.?`, `m.?`, `o.?`).
    ///
    /// The HGVS spec uses bracketed `?` only as a whole-allele marker (`[?]`),
    /// never mixed with concrete edits inside the same bracket. Compact form is
    /// suppressed when any sub-variant matches this so the output isn't visually
    /// ambiguous with the spec-sanctioned `[?]` form.
    pub(crate) fn is_loc_edit_unknown(&self) -> bool {
        match self {
            HgvsVariant::Genome(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Cds(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Tx(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Rna(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Mt(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Circular(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Protein(v) => match &v.loc_edit.edit {
                Mu::Unknown => true,
                Mu::Certain(e) | Mu::Uncertain(e) => e.is_whole_protein_unknown(),
            },
            HgvsVariant::RnaFusion(_)
            | HgvsVariant::GenomeRing(_)
            | HgvsVariant::Supernumerary(_)
            | HgvsVariant::Allele(_)
            | HgvsVariant::NullAllele
            | HgvsVariant::UnknownAllele => false,
        }
    }

    /// Check if all variants in a slice share the same accession, coordinate type, and
    /// gene-symbol selector. Used to determine whether the compact allele form can be
    /// used without losing information.
    ///
    /// The compact form (`ACC(GENE):c.[edit1;edit2]`) can carry exactly one
    /// `(GENE)` selector on the prefix. Per the #121 round-trip policy
    /// ("preserve when present; do not synthesize"), the compact form is only
    /// safe when every sub-variant agrees on `gene_symbol`: both `Some(g)`
    /// (same `g`) and all-`None` are safe; any disagreement (different `g`,
    /// or `Some` vs `None`) falls back to the expanded form so each
    /// sub-variant emits its own selector via per-variant `Display`.
    pub(crate) fn all_share_accession_and_type(variants: &[HgvsVariant]) -> bool {
        let Some(first) = variants.first() else {
            return true;
        };
        let first_acc = first.accession();
        let first_type = first.variant_type();
        let first_gene = first.gene_symbol();

        // Don't use compact form for types that aren't simple coordinate-based variants,
        // or when the first variant has no accession (e.g. NullAllele, UnknownAllele)
        if first_acc.is_none()
            || matches!(first_type, "allele" | "null" | "unknown" | "r::r" | "g::g")
        {
            return false;
        }

        variants[1..].iter().all(|v| {
            v.variant_type() == first_type
                && v.accession() == first_acc
                && v.gene_symbol() == first_gene
        })
    }
}

/// True if a nucleotide-edit `Mu` wrapper is a certain, position-only identity
/// edit (`NaEdit::Identity { sequence: None, whole_entity: false }`) — the
/// editless "position reference" form produced when a genome interval is given
/// with no edit (e.g. `g.pter_qter`). Used by the supernumerary whole-chromosome
/// Display to drop the trailing `=` so `pter_qtersup` round-trips.
///
/// Note: this is stricter than [`is_position_identity_lhs`], which accepts any
/// position-bound identity (`whole_entity: false`) regardless of `sequence`.
/// This predicate additionally requires `sequence: None` (no `<ref>` written),
/// so `g.123A=` satisfies `is_position_identity_lhs` but not this — they are not
/// redundant.
fn is_position_only_identity(edit: &Mu<NaEdit>) -> bool {
    matches!(
        edit,
        Mu::Certain(NaEdit::Identity {
            sequence: None,
            whole_entity: false,
        })
    )
}

/// Wrap a valid supernumerary inner (an editless whole-chromosome genome
/// variant, or a [`GenomeRing`]) in [`HgvsVariant::Supernumerary`]. Returns
/// `None` for any other inner — `sup` only applies to those two shapes
/// (HGVS DNA/complex.md:33-34). This is the single constructor enforcing the
/// supernumerary invariant, mirroring [`GenomeRing::new`].
pub(crate) fn make_supernumerary(inner: HgvsVariant) -> Option<HgvsVariant> {
    let ok = matches!(&inner, HgvsVariant::GenomeRing(_))
        || matches!(&inner, HgvsVariant::Genome(g) if is_position_only_identity(&g.loc_edit.edit));
    ok.then(|| HgvsVariant::Supernumerary(Box::new(inner)))
}

/// True if a nucleotide-edit `Mu` wrapper represents the per-variant unknown form
/// (e.g. `c.?`, `r.?`).
fn is_na_edit_unknown(edit: &Mu<NaEdit>) -> bool {
    match edit {
        Mu::Unknown => true,
        Mu::Certain(e) | Mu::Uncertain(e) => e.is_whole_entity_unknown(),
    }
}

impl fmt::Display for HgvsVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HgvsVariant::Genome(v) => write!(f, "{}", v),
            HgvsVariant::Cds(v) => write!(f, "{}", v),
            HgvsVariant::Tx(v) => write!(f, "{}", v),
            HgvsVariant::Rna(v) => write!(f, "{}", v),
            HgvsVariant::Protein(v) => write!(f, "{}", v),
            HgvsVariant::Mt(v) => write!(f, "{}", v),
            HgvsVariant::Circular(v) => write!(f, "{}", v),
            HgvsVariant::RnaFusion(v) => write!(f, "{}", v),
            HgvsVariant::GenomeRing(g) => write!(f, "{}", g),
            HgvsVariant::Supernumerary(inner) => {
                // A supernumerary ring is bracketed (`g.[<ring>]sup`); a bare
                // supernumerary whole chromosome is not (`g.pter_qtersup`).
                match inner.as_ref() {
                    HgvsVariant::GenomeRing(g) => g.fmt_bracketed(f)?,
                    // A whole-chromosome supernumerary carries a position-only
                    // (editless) interval such as `pter_qter`; the spec form is
                    // `pter_qtersup`, with no trailing identity `=`. Render the
                    // location alone so `sup` attaches directly to it.
                    HgvsVariant::Genome(v) if is_position_only_identity(&v.loc_edit.edit) => {
                        write_accession_with_optional_gene(
                            f,
                            &v.accession,
                            v.gene_symbol.as_deref(),
                        )?;
                        write!(f, ":g.{}", v.loc_edit.location)?;
                    }
                    other => write!(f, "{}", other)?,
                }
                write!(f, "sup")
            }
            HgvsVariant::Allele(a) => write!(f, "{}", a),
            HgvsVariant::NullAllele => write!(f, "0"),
            HgvsVariant::UnknownAllele => write!(f, "?"),
        }
    }
}

/// Genomic variant (g.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GenomeVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<GenomeInterval, NaEdit>,
}

impl GenomeVariant {
    /// Format just the position+edit portion (without `accession:g.` prefix).
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_genome_loc_edit(f, &self.loc_edit)
    }
}

impl fmt::Display for GenomeVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":g.")?;
        self.fmt_loc_edit(f)
    }
}

/// A ring chromosome described as two or more genome edits on a single
/// accession joined by `::` break junctions (HGVS DNA/complex.md:127/130 —
/// the ISCN2020 meaning of `::`). Each segment is a normal genome loc-edit
/// (typically a `del`); `Display` joins them with `::`.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GenomeRing {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub segments: Vec<LocEdit<GenomeInterval, NaEdit>>,
}

impl GenomeRing {
    /// Build a ring from its segments. Requires at least two segments joined
    /// by `::` (a ring needs >= 2 break junctions on one accession); returns
    /// `None` otherwise, since a 0/1-segment ring would not round-trip.
    pub fn new(
        accession: Accession,
        gene_symbol: Option<String>,
        segments: Vec<LocEdit<GenomeInterval, NaEdit>>,
    ) -> Option<Self> {
        if segments.len() < 2 {
            return None;
        }
        Some(Self {
            accession,
            gene_symbol,
            segments,
        })
    }

    /// Format the `::`-joined ring segments (without the `accession:g.` prefix
    /// or any surrounding brackets). Shared by [`Display`] and
    /// [`GenomeRing::fmt_bracketed`] so the two renderings stay in lock-step.
    fn fmt_segments(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, segment) in self.segments.iter().enumerate() {
            if i > 0 {
                write!(f, "::")?;
            }
            fmt_genome_loc_edit(f, segment)?;
        }
        Ok(())
    }

    /// Format the ring with its `::`-joined segments wrapped in `[...]`, i.e.
    /// `accession:g.[seg1::seg2…]`. This is the shape used when the ring is the
    /// inner of a supernumerary ring chromosome (`g.[<ring>]sup`,
    /// DNA/complex.md:161). The standalone ring [`Display`] is unbracketed.
    fn fmt_bracketed(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":g.[")?;
        self.fmt_segments(f)?;
        write!(f, "]")
    }
}

impl fmt::Display for GenomeRing {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":g.")?;
        self.fmt_segments(f)
    }
}

/// Coding DNA variant (c.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CdsVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<CdsInterval, NaEdit>,
}

impl CdsVariant {
    /// Format just the position+edit portion (without `accession:c.` prefix).
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // For whole-entity identity (c.=) or unknown (c.?), skip the position
        if let Some(edit) = self.loc_edit.edit.inner() {
            if edit.is_whole_entity_identity() || edit.is_whole_entity_unknown() {
                return write!(f, "{}", self.loc_edit.edit);
            }
        }
        // Predicted form: wrap position+edit in parens (mirrors the
        // protein canonical form). #241.
        if self.loc_edit.edit.is_uncertain() {
            if let Some(edit) = self.loc_edit.edit.inner() {
                return write!(f, "({}{})", self.loc_edit.location, edit);
            }
        }
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for CdsVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":c.")?;
        self.fmt_loc_edit(f)
    }
}

/// Non-coding transcript variant (n.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct TxVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<TxInterval, NaEdit>,
}

impl TxVariant {
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // For whole-entity identity / unknown / no-product (n.=, n.?, n.0
        // and their predicted wrappers n.(=), n.(?), n.(0)), skip the
        // position. #288 — mirrors the same short-circuits on g./c./r..
        if let Some(edit) = self.loc_edit.edit.inner() {
            if edit.is_whole_entity_identity()
                || edit.is_whole_entity_unknown()
                || edit.is_no_product()
            {
                return write!(f, "{}", self.loc_edit.edit);
            }
        }
        // Predicted form: wrap position+edit in parens. #241.
        if self.loc_edit.edit.is_uncertain() {
            if let Some(edit) = self.loc_edit.edit.inner() {
                return write!(f, "({}{})", self.loc_edit.location, edit);
            }
        }
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for TxVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":n.")?;
        self.fmt_loc_edit(f)
    }
}

/// RNA variant (r.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RnaVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<RnaInterval, NaEdit>,
}

impl RnaVariant {
    /// Format just the position+edit portion (without `accession:r.` prefix).
    ///
    /// RNA uses lowercase nucleotides per HGVS spec.
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.loc_edit.edit {
            Mu::Certain(edit) => {
                // For whole-entity patterns (r.=, r.?, r.spl, r.0), skip the position
                if edit.is_whole_entity() {
                    write!(f, "{}", edit.to_rna_string())
                } else {
                    write!(f, "{}{}", self.loc_edit.location, edit.to_rna_string())
                }
            }
            Mu::Uncertain(edit) => {
                // For whole-entity patterns, skip the position
                if edit.is_whole_entity() {
                    write!(f, "({})", edit.to_rna_string())
                } else {
                    // Predicted form wraps position+edit in parens
                    // (mirrors the protein canonical form). #241.
                    write!(f, "({}{})", self.loc_edit.location, edit.to_rna_string())
                }
            }
            Mu::Unknown => write!(f, "?"),
        }
    }
}

impl fmt::Display for RnaVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":r.")?;
        self.fmt_loc_edit(f)
    }
}

/// Protein variant (p.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ProteinVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<ProtInterval, ProteinEdit>,
}

impl ProteinVariant {
    /// Format just the position+edit portion (without `accession:p.` prefix).
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // For whole-protein identity (p.= or p.(=)), no-protein (p.0), or whole-protein unknown (p.?), skip the position
        if let Some(edit) = self.loc_edit.edit.inner() {
            if edit.is_whole_protein_identity()
                || edit.is_no_protein()
                || edit.is_whole_protein_unknown()
            {
                return write!(f, "{}", self.loc_edit.edit);
            }
        }
        // For predicted protein changes (uncertain edit), wrap position+edit in parentheses
        // e.g., p.(Arg248Gln) instead of p.Arg248(Gln)
        if self.loc_edit.edit.is_uncertain() {
            if let Some(edit) = self.loc_edit.edit.inner() {
                return write!(f, "({}{})", self.loc_edit.location, edit);
            }
        }
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for ProteinVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Per HGVS syntax.yaml lines 119–128 the protein-variant grammar is
        // `sequence_identifier ":" coordinate_type "." aa_position alternate_base`
        // — no parenthesized selector position. The `accession(selector):c.`
        // form (syntax.yaml line 213) is reserved for genomic-to-transcript
        // projection on c./n., not for p. The struct still stores
        // `gene_symbol` so callers can read it programmatically and so the
        // parser (which is permissive about the legacy `NP_*(GENE):p.`
        // form) can round-trip the value; it is just not emitted here. See
        // #310.
        write!(f, "{}:p.", self.accession)?;
        self.fmt_loc_edit(f)
    }
}

/// Mitochondrial variant (m.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct MtVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<GenomeInterval, NaEdit>,
}

impl MtVariant {
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // For whole-entity identity / unknown (m.=, m.? and their
        // predicted wrappers m.(=), m.(?)), skip the position. #288 —
        // mirrors the genome short-circuit. There is no m.0 form
        // (DNA isn't "produced"), so no-product isn't accepted here.
        if let Some(edit) = self.loc_edit.edit.inner() {
            if edit.is_whole_entity_identity() || edit.is_whole_entity_unknown() {
                return write!(f, "{}", self.loc_edit.edit);
            }
        }
        // Predicted form: wrap position+edit in parens. #241.
        if self.loc_edit.edit.is_uncertain() {
            if let Some(edit) = self.loc_edit.edit.inner() {
                return write!(f, "({}{})", self.loc_edit.location, edit);
            }
        }
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for MtVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":m.")?;
        self.fmt_loc_edit(f)
    }
}

/// Circular DNA variant (o.) - SVD-WG006
///
/// Used for circular DNA sequences like plasmids, bacterial chromosomes,
/// and mitochondrial DNA where variants can span the origin.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CircularVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<GenomeInterval, NaEdit>,
}

impl CircularVariant {
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // For whole-entity identity / unknown (o.=, o.? and their
        // predicted wrappers o.(=), o.(?)), skip the position. #288 —
        // mirrors the genome short-circuit. There is no o.0 form.
        if let Some(edit) = self.loc_edit.edit.inner() {
            if edit.is_whole_entity_identity() || edit.is_whole_entity_unknown() {
                return write!(f, "{}", self.loc_edit.edit);
            }
        }
        // Predicted form: wrap position+edit in parens. #241.
        if self.loc_edit.edit.is_uncertain() {
            if let Some(edit) = self.loc_edit.edit.inner() {
                return write!(f, "({}{})", self.loc_edit.location, edit);
            }
        }
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for CircularVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":o.")?;
        self.fmt_loc_edit(f)
    }
}

/// RNA Fusion breakpoint - SVD-WG007
///
/// Represents one side of an RNA fusion with accession, gene symbol, and interval
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RnaFusionBreakpoint {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub interval: RnaInterval,
}

impl fmt::Display for RnaFusionBreakpoint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(ref gene) = self.gene_symbol {
            write!(f, "{}({}):r.{}", self.accession, gene, self.interval)
        } else {
            write!(f, "{}:r.{}", self.accession, self.interval)
        }
    }
}

/// RNA Fusion variant (::) - SVD-WG007
///
/// Represents a fusion transcript where two RNA sequences are joined.
/// Common in cancer genomics (BCR-ABL1, EML4-ALK, etc.)
///
/// Format: `accession:r.interval::accession:r.interval`
/// Example: `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924`
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RnaFusionVariant {
    /// 5' partner (upstream)
    pub five_prime: RnaFusionBreakpoint,
    /// 3' partner (downstream)
    pub three_prime: RnaFusionBreakpoint,
}

impl RnaFusionVariant {
    pub fn new(five_prime: RnaFusionBreakpoint, three_prime: RnaFusionBreakpoint) -> Self {
        Self {
            five_prime,
            three_prime,
        }
    }
}

impl fmt::Display for RnaFusionVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}::{}", self.five_prime, self.three_prime)
    }
}

/// Check if a variant causes a frameshift (net indel length not divisible by 3).
///
/// Returns false if the indel length is 0 or cannot be determined.
pub fn is_frameshift(variant: &HgvsVariant) -> bool {
    match crate::python_helpers::get_indel_length(variant) {
        Some(len) if len != 0 => len % 3 != 0,
        _ => false,
    }
}

/// Like [`is_frameshift`] but uses a provider so wraparound `m.`/`o.`
/// ranges compute the spec-correct circular indel length via
/// [`crate::python_helpers::get_indel_length_with_provider`].
pub fn is_frameshift_with_provider<P: crate::reference::provider::ReferenceProvider + ?Sized>(
    variant: &HgvsVariant,
    provider: &P,
) -> bool {
    match crate::python_helpers::get_indel_length_with_provider(variant, provider) {
        Some(len) if len != 0 => len % 3 != 0,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::Base;
    use crate::hgvs::location::GenomePos;

    #[test]
    fn test_accession_display() {
        let acc = Accession::new("NM", "000088", Some(3));
        assert_eq!(format!("{}", acc), "NM_000088.3");

        let acc_no_version = Accession::new("NM", "000088", None);
        assert_eq!(format!("{}", acc_no_version), "NM_000088");
    }

    #[test]
    fn test_genome_variant_display() {
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };
        assert_eq!(format!("{}", variant), "NC_000001.11:g.12345A>G");
    }

    #[test]
    fn test_cds_variant_display() {
        use crate::hgvs::location::CdsPos;

        let variant = CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(459)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        assert_eq!(format!("{}", variant), "NM_000088.3:c.459del");
    }

    #[test]
    fn test_allele_cis_display() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let allele = AlleleVariant::cis(vec![var1, var2]);
        let allele_variant = HgvsVariant::Allele(allele);

        assert_eq!(
            format!("{}", allele_variant),
            "NM_000088.3:c.[100A>G;200C>T]"
        );
    }

    #[test]
    fn test_allele_trans_display() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let allele = AlleleVariant::trans(vec![var1, var2]);
        let allele_variant = HgvsVariant::Allele(allele);

        assert_eq!(
            format!("{}", allele_variant),
            "NM_000088.3:c.[100A>G];[200C>T]"
        );
    }

    #[test]
    fn test_allele_cis_three_variants() {
        use crate::hgvs::location::CdsPos;

        let make_var = |pos, alt| {
            HgvsVariant::Cds(CdsVariant {
                accession: Accession::new("NM", "000088", Some(3)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(pos)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: alt,
                    },
                ),
            })
        };

        let allele = AlleleVariant::cis(vec![
            make_var(100, Base::G),
            make_var(200, Base::C),
            make_var(300, Base::T),
        ]);

        assert_eq!(
            format!("{}", HgvsVariant::Allele(allele)),
            "NM_000088.3:c.[100A>G;200A>C;300A>T]"
        );
    }

    #[test]
    fn test_allele_trans_three_variants() {
        use crate::hgvs::location::CdsPos;

        let make_var = |pos, alt| {
            HgvsVariant::Cds(CdsVariant {
                accession: Accession::new("NM", "000088", Some(3)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(pos)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: alt,
                    },
                ),
            })
        };

        let allele = AlleleVariant::trans(vec![
            make_var(100, Base::G),
            make_var(200, Base::C),
            make_var(300, Base::T),
        ]);

        assert_eq!(
            format!("{}", HgvsVariant::Allele(allele)),
            "NM_000088.3:c.[100A>G];[200A>C];[300A>T]"
        );
    }

    #[test]
    fn test_allele_unknown_phase_compact_form() {
        // 2+ unknown-phase variants sharing accession + coordinate type emit the
        // bracket-less compact form `ACC:c.edit1(;)edit2`.
        use crate::hgvs::location::CdsPos;

        let make_var = |pos, alt| {
            HgvsVariant::Cds(CdsVariant {
                accession: Accession::new("NM", "000088", Some(3)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(pos)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: alt,
                    },
                ),
            })
        };

        let pair =
            AlleleVariant::unknown_phase(vec![make_var(100, Base::G), make_var(200, Base::C)]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(pair)),
            "NM_000088.3:c.100A>G(;)200A>C"
        );

        let triple = AlleleVariant::unknown_phase(vec![
            make_var(100, Base::G),
            make_var(200, Base::C),
            make_var(300, Base::T),
        ]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(triple)),
            "NM_000088.3:c.100A>G(;)200A>C(;)300A>T"
        );
    }

    #[test]
    fn test_allele_mixed_coordinate_types_uses_expanded_form() {
        // Same accession, different coordinate types (c. and n.) → expanded form.
        // `all_share_accession_and_type` requires both fields to match.
        use crate::hgvs::location::{CdsPos, TxPos};

        let coding = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let noncoding = HgvsVariant::Tx(TxVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let cis = AlleleVariant::cis(vec![coding.clone(), noncoding.clone()]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis)),
            "[NM_000088.3:c.100A>G;NM_000088.3:n.200C>T]"
        );

        let trans = AlleleVariant::trans(vec![coding, noncoding]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(trans)),
            "[NM_000088.3:c.100A>G];[NM_000088.3:n.200C>T]"
        );
    }

    #[test]
    fn test_allele_singleton_emits_bare_form() {
        use crate::hgvs::location::CdsPos;

        let make_var = || {
            HgvsVariant::Cds(CdsVariant {
                accession: Accession::new("NM", "000088", Some(3)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(100)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::G,
                    },
                ),
            })
        };
        let bare = "NM_000088.3:c.100A>G";

        for allele in [
            AlleleVariant::cis(vec![make_var()]),
            AlleleVariant::trans(vec![make_var()]),
            AlleleVariant::unknown_phase(vec![make_var()]),
            AlleleVariant::mosaic(vec![make_var()]),
            AlleleVariant::chimeric(vec![make_var()]),
        ] {
            assert_eq!(format!("{}", HgvsVariant::Allele(allele)), bare);
        }
    }

    #[test]
    fn test_allele_mixed_accession_uses_expanded_form() {
        use crate::hgvs::location::CdsPos;

        // Different accessions → expanded form (no compact shorthand)
        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000099", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let cis = AlleleVariant::cis(vec![var1.clone(), var2.clone()]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis)),
            "[NM_000088.3:c.100A>G;NM_000099.1:c.200C>T]"
        );

        let trans = AlleleVariant::trans(vec![var1, var2]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(trans)),
            "[NM_000088.3:c.100A>G];[NM_000099.1:c.200C>T]"
        );
    }

    #[test]
    fn test_allele_with_unknown_sub_variant_uses_expanded_form() {
        // Per HGVS spec, bracketed `?` is reserved for the whole-allele marker (`[?]`).
        // An allele containing a `c.?` sub-variant must NOT collapse to `[?;100A>G]` —
        // that would visually conflict with the spec-sanctioned form.
        use crate::hgvs::location::CdsPos;

        let unknown = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::whole_entity_unknown(),
            ),
        });
        let concrete = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let trans = AlleleVariant::trans(vec![concrete, unknown]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(trans)),
            "[NM_000088.3:c.100A>G];[NM_000088.3:c.?]"
        );
    }

    #[test]
    fn test_allele_compact_form_requires_matching_gene_symbol() {
        // Hand-built cis allele: same accession, same coordinate type, but the
        // sub-variants disagree on gene_symbol. The compact form
        // (`ACC(GENE):c.[edit1;edit2]`) can only carry one selector on the
        // prefix, so collapsing two distinct selectors would silently drop
        // the second. Per the #121 round-trip policy ("preserve when present;
        // do not synthesize"), we fall back to the expanded form which
        // emits each sub-variant's selector verbatim.
        use crate::hgvs::location::CdsPos;

        let acc = Accession::new("NM", "000088", Some(3));
        let v_a = HgvsVariant::Cds(CdsVariant {
            accession: acc.clone(),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let v_b = HgvsVariant::Cds(CdsVariant {
            accession: acc.clone(),
            gene_symbol: Some("COL1A2".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let cis = AlleleVariant::cis(vec![v_a.clone(), v_b.clone()]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis)),
            "[NM_000088.3(COL1A1):c.100A>G;NM_000088.3(COL1A2):c.200C>T]",
            "mismatched gene symbols must use expanded form to preserve both"
        );

        let trans = AlleleVariant::trans(vec![v_a, v_b]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(trans)),
            "[NM_000088.3(COL1A1):c.100A>G];[NM_000088.3(COL1A2):c.200C>T]",
            "trans expanded form already carries per-variant selectors"
        );
    }

    #[test]
    fn test_allele_compact_form_requires_matching_gene_symbol_some_vs_none() {
        // First sub-variant has Some(gene), second has None. Compact form
        // would either lose the gene symbol (emit bare prefix) or synthesize
        // one for the bare sub-variant (emit prefix-with-gene). Both violate
        // the #121 policy; fall back to expanded form.
        use crate::hgvs::location::CdsPos;

        let acc = Accession::new("NM", "000088", Some(3));
        let v_with = HgvsVariant::Cds(CdsVariant {
            accession: acc.clone(),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let v_bare = HgvsVariant::Cds(CdsVariant {
            accession: acc.clone(),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let cis = AlleleVariant::cis(vec![v_with, v_bare]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis)),
            "[NM_000088.3(COL1A1):c.100A>G;NM_000088.3:c.200C>T]",
            "Some/None gene_symbol mismatch must use expanded form"
        );
    }

    #[test]
    fn test_allele_compact_form_keeps_compact_when_all_match() {
        // All sub-variants share Some(same gene): compact form is safe and
        // emits the selector once on the prefix (regression test for the
        // happy path so the new compact-form rule isn't over-constrained).
        use crate::hgvs::location::CdsPos;

        let acc = Accession::new("NM", "000088", Some(3));
        let mk = |pos: i64, alt: Base| {
            HgvsVariant::Cds(CdsVariant {
                accession: acc.clone(),
                gene_symbol: Some("COL1A1".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(pos)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: alt,
                    },
                ),
            })
        };
        let cis = AlleleVariant::cis(vec![mk(100, Base::G), mk(200, Base::T)]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis)),
            "NM_000088.3(COL1A1):c.[100A>G;200A>T]",
            "matching gene symbols must keep compact form"
        );

        // All-None: existing compact behavior must be preserved.
        let mk_none = |pos: i64, alt: Base| {
            HgvsVariant::Cds(CdsVariant {
                accession: acc.clone(),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(pos)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: alt,
                    },
                ),
            })
        };
        let cis_none = AlleleVariant::cis(vec![mk_none(100, Base::G), mk_none(200, Base::T)]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis_none)),
            "NM_000088.3:c.[100A>G;200A>T]",
            "all-None gene_symbol must use compact form unchanged"
        );
    }

    #[test]
    fn test_allele_null_variant_singleton_delegates() {
        // NullAllele has no accession; the singleton path must not invoke
        // `use_compact_form` (which would deref the missing accession).
        let null = HgvsVariant::NullAllele;
        let cis = AlleleVariant::cis(vec![null]);
        assert_eq!(format!("{}", HgvsVariant::Allele(cis)), "0");
    }

    #[test]
    fn test_allele_accession() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });

        let allele = AlleleVariant::cis(vec![var1]);
        let allele_variant = HgvsVariant::Allele(allele);

        assert_eq!(
            &*allele_variant
                .accession()
                .expect("Expected accession")
                .prefix,
            "NM"
        );
        assert!(allele_variant.is_allele());
    }

    #[test]
    fn test_allele_mosaic_display() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let allele = AlleleVariant::mosaic(vec![var1, var2]);
        let allele_variant = HgvsVariant::Allele(allele);

        // Mosaic uses single slash separator: var1/var2
        assert_eq!(
            format!("{}", allele_variant),
            "NM_000088.3:c.100A>G/NM_000088.3:c.200C>T"
        );
    }

    #[test]
    fn test_allele_chimeric_display() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let allele = AlleleVariant::chimeric(vec![var1, var2]);
        let allele_variant = HgvsVariant::Allele(allele);

        // Chimerism uses double slash separator: var1//var2
        assert_eq!(
            format!("{}", allele_variant),
            "NM_000088.3:c.100A>G//NM_000088.3:c.200C>T"
        );
    }

    #[test]
    fn test_allele_phase_variants() {
        // Verify all phase types are correctly identified
        assert!(!AllelePhase::Cis.is_mosaic());
        assert!(!AllelePhase::Trans.is_mosaic());
        assert!(!AllelePhase::Unknown.is_mosaic());
        assert!(AllelePhase::Mosaic.is_mosaic());
        assert!(!AllelePhase::Chimeric.is_mosaic());
        assert!(!AllelePhase::AndOr.is_mosaic());
        assert!(!AllelePhase::Products.is_mosaic());

        assert!(!AllelePhase::Cis.is_chimeric());
        assert!(!AllelePhase::Trans.is_chimeric());
        assert!(!AllelePhase::Unknown.is_chimeric());
        assert!(!AllelePhase::Mosaic.is_chimeric());
        assert!(AllelePhase::Chimeric.is_chimeric());
        assert!(!AllelePhase::AndOr.is_chimeric());
        assert!(!AllelePhase::Products.is_chimeric());
    }

    // ============== Accession Tests ==============

    #[test]
    fn test_accession_new() {
        let acc = Accession::new("NM", "000088", Some(3));
        assert_eq!(&*acc.prefix, "NM");
        assert_eq!(&*acc.number, "000088");
        assert_eq!(acc.version, Some(3));
        assert!(!acc.ensembl_style);
    }

    #[test]
    fn test_accession_with_style() {
        let acc = Accession::with_style("ENST", "00000012345", Some(1), true);
        assert!(acc.ensembl_style);
        assert_eq!(format!("{}", acc), "ENST00000012345.1");
    }

    #[test]
    fn test_accession_from_assembly() {
        let acc = Accession::from_assembly("GRCh37", "chr1");
        assert!(acc.is_assembly_ref());
        assert_eq!(format!("{}", acc), "GRCh37(chr1)");
    }

    #[test]
    fn test_accession_is_ensembl_prefix() {
        // Human prefixes (empty species code).
        assert!(Accession::is_ensembl_prefix("ENST"));
        assert!(Accession::is_ensembl_prefix("ENSG"));
        assert!(Accession::is_ensembl_prefix("ENSP"));
        assert!(Accession::is_ensembl_prefix("ENSE"));
        assert!(Accession::is_ensembl_prefix("ENSR"));

        // Non-human species prefixes — same five feature suffixes,
        // species code (`MUS`, `RNO`, `BTA`, `DAR`, `CAF`, `SSC`, …) sits
        // between `ENS` and the suffix.
        assert!(Accession::is_ensembl_prefix("ENSMUST"));
        assert!(Accession::is_ensembl_prefix("ENSMUSG"));
        assert!(Accession::is_ensembl_prefix("ENSMUSP"));
        assert!(Accession::is_ensembl_prefix("ENSMUSE"));
        assert!(Accession::is_ensembl_prefix("ENSMUSR"));
        assert!(Accession::is_ensembl_prefix("ENSRNOT"));
        assert!(Accession::is_ensembl_prefix("ENSBTAG"));
        assert!(Accession::is_ensembl_prefix("ENSDARP"));
        assert!(Accession::is_ensembl_prefix("ENSCAFT"));
        assert!(Accession::is_ensembl_prefix("ENSSSCG"));

        // Negative cases.
        assert!(!Accession::is_ensembl_prefix("NM"));
        assert!(!Accession::is_ensembl_prefix("NC"));
        assert!(!Accession::is_ensembl_prefix("ENS")); // no feature suffix
        assert!(!Accession::is_ensembl_prefix("ENSA")); // suffix not in {G,T,P,E,R}
        assert!(!Accession::is_ensembl_prefix("ENSmust")); // lowercase rejected
        assert!(!Accession::is_ensembl_prefix("ENSMUS1T")); // digit in species code
        assert!(!Accession::is_ensembl_prefix("")); // empty
        assert!(!Accession::is_ensembl_prefix("ANYTHING"));
        // Species-code length is bounded at 5 chars (Ensembl convention).
        assert!(Accession::is_ensembl_prefix("ENSAAAAAT")); // 5-letter species code + feature: accept
        assert!(!Accession::is_ensembl_prefix("ENSAAAAAAT")); // 6-letter species code + feature: reject
        assert!(!Accession::is_ensembl_prefix("ENSAAAAAAAAAT")); // grossly overlong: reject
    }

    #[test]
    fn test_accession_is_ensembl() {
        let ensembl = Accession::new("ENST", "00000012345", Some(1));
        assert!(ensembl.is_ensembl());

        let refseq = Accession::new("NM", "000088", Some(3));
        assert!(!refseq.is_ensembl());
    }

    #[test]
    fn test_accession_validate_ensembl() {
        // Valid Ensembl ID (11 digits)
        let valid = Accession::new("ENST", "00000012345", Some(1));
        assert!(valid.validate_ensembl());

        // Non-Ensembl IDs should always validate
        let refseq = Accession::new("NM", "000088", Some(3));
        assert!(refseq.validate_ensembl());
    }

    #[test]
    fn test_accession_inferred_variant_type() {
        assert_eq!(
            Accession::new("NC", "000001", Some(11)).inferred_variant_type(),
            Some("g")
        );
        assert_eq!(
            Accession::new("NG", "012345", Some(1)).inferred_variant_type(),
            Some("g")
        );
        assert_eq!(
            Accession::new("NM", "000088", Some(3)).inferred_variant_type(),
            Some("c")
        );
        assert_eq!(
            Accession::new("NR", "123456", Some(1)).inferred_variant_type(),
            Some("n")
        );
        assert_eq!(
            Accession::new("NP", "000079", Some(2)).inferred_variant_type(),
            Some("p")
        );
        assert_eq!(
            Accession::new("ENST", "00000012345", Some(1)).inferred_variant_type(),
            Some("c")
        );
        assert_eq!(
            Accession::new("ENSG", "00000012345", Some(1)).inferred_variant_type(),
            Some("g")
        );
        assert_eq!(
            Accession::new("ENSP", "00000012345", Some(1)).inferred_variant_type(),
            Some("p")
        );
        assert_eq!(
            Accession::new("LRG", "1", None).inferred_variant_type(),
            Some("g")
        );
        // LRG transcript discriminator (`t<M>`) → coding (c.)
        // The transcript discriminator is captured in `number` (e.g. "1t1"),
        // so the dispatch must inspect it rather than only the prefix.
        // See https://www.lrg-sequence.org/faq/ — `t<N>` is the transcript.
        assert_eq!(
            Accession::new("LRG", "1t1", None).inferred_variant_type(),
            Some("c")
        );
        assert_eq!(
            Accession::new("LRG", "292t1", None).inferred_variant_type(),
            Some("c")
        );
        // LRG protein discriminator (`p<M>`) → protein (p.)
        assert_eq!(
            Accession::new("LRG", "1p1", None).inferred_variant_type(),
            Some("p")
        );
        assert_eq!(
            Accession::new("LRG", "673p1", None).inferred_variant_type(),
            Some("p")
        );
        // UniProt single-letter prefix
        assert_eq!(
            Accession::new("P", "12345", None).inferred_variant_type(),
            Some("p")
        );
        // Unknown prefix
        assert_eq!(
            Accession::new("XX", "12345", None).inferred_variant_type(),
            None
        );
    }

    #[test]
    fn test_accession_is_uniprot() {
        // Valid UniProt
        let uniprot = Accession::new("P", "12345", None);
        assert!(uniprot.is_uniprot());

        // Not UniProt - prefix too long
        let not_uniprot = Accession::new("NM", "000088", Some(3));
        assert!(!not_uniprot.is_uniprot());

        // Not UniProt - number wrong length
        let wrong_length = Accession::new("P", "123", None);
        assert!(!wrong_length.is_uniprot());
    }

    #[test]
    fn test_accession_is_mitochondrial() {
        // GRCh38 rCRS mitochondrion
        assert!(Accession::new("NC", "012920", Some(1)).is_mitochondrial());

        // Deprecated rCRS draft
        assert!(Accession::new("NC", "001807", Some(1)).is_mitochondrial());

        // Version-agnostic — coercion keys off prefix+number, not version
        assert!(Accession::new("NC", "012920", None).is_mitochondrial());

        // Other NC accessions (e.g. chromosome 1) are not mitochondrial
        assert!(!Accession::new("NC", "000001", Some(11)).is_mitochondrial());

        // Same number, different prefix is not mitochondrial
        assert!(!Accession::new("NM", "012920", Some(1)).is_mitochondrial());
    }

    #[test]
    fn test_accession_is_lrg_prefix() {
        assert!(Accession::is_lrg_prefix("LRG"));
        assert!(!Accession::is_lrg_prefix("lrg"));
        assert!(!Accession::is_lrg_prefix("LRG_"));
        assert!(!Accession::is_lrg_prefix("NM"));
        assert!(!Accession::is_lrg_prefix(""));
    }

    #[test]
    fn test_accession_is_lrg() {
        assert!(Accession::new("LRG", "1", None).is_lrg());
        assert!(Accession::new("LRG", "1t1", None).is_lrg());
        assert!(Accession::new("LRG", "1p1", None).is_lrg());
        assert!(!Accession::new("NM", "000088", Some(3)).is_lrg());
        assert!(!Accession::new("NC", "000001", Some(11)).is_lrg());
    }

    #[test]
    fn test_accession_base() {
        let acc = Accession::new("NM", "000088", Some(3));
        assert_eq!(acc.base(), "NM_000088");

        let ensembl = Accession::with_style("ENST", "00000012345", Some(1), true);
        assert_eq!(ensembl.base(), "ENST00000012345");

        let assembly = Accession::from_assembly("GRCh38", "chrX");
        assert_eq!(assembly.base(), "GRCh38(chrX)");
    }

    #[test]
    fn test_accession_full() {
        let acc = Accession::new("NM", "000088", Some(3));
        assert_eq!(acc.full(), "NM_000088.3");

        let no_version = Accession::new("NM", "000088", None);
        assert_eq!(no_version.full(), "NM_000088");

        let ensembl = Accession::with_style("ENST", "00000012345", Some(1), true);
        assert_eq!(ensembl.full(), "ENST00000012345.1");
    }

    // ============== LocEdit Tests ==============

    #[test]
    fn test_loc_edit_new() {
        use crate::hgvs::location::GenomePos;

        let loc_edit: LocEdit<GenomeInterval, NaEdit> = LocEdit::new(
            GenomeInterval::point(GenomePos::new(100)),
            NaEdit::Deletion {
                sequence: None,
                length: None,
            },
        );
        assert!(matches!(loc_edit.edit, Mu::Certain(_)));
    }

    #[test]
    fn test_loc_edit_with_uncertainty() {
        use crate::hgvs::location::GenomePos;

        let loc_edit: LocEdit<GenomeInterval, NaEdit> = LocEdit::with_uncertainty(
            GenomeInterval::point(GenomePos::new(100)),
            Mu::Uncertain(NaEdit::Deletion {
                sequence: None,
                length: None,
            }),
        );
        assert!(matches!(loc_edit.edit, Mu::Uncertain(_)));
    }

    // ============== AllelePhase Tests ==============

    #[test]
    fn test_allele_phase_display() {
        assert_eq!(format!("{}", AllelePhase::Cis), "cis");
        assert_eq!(format!("{}", AllelePhase::Trans), "trans");
        assert_eq!(format!("{}", AllelePhase::Unknown), "unknown");
        assert_eq!(format!("{}", AllelePhase::Mosaic), "mosaic");
        assert_eq!(format!("{}", AllelePhase::Chimeric), "chimeric");
        assert_eq!(format!("{}", AllelePhase::AndOr), "and/or");
        assert_eq!(format!("{}", AllelePhase::Products), "products");
    }

    // ============== AlleleVariant Tests ==============

    #[test]
    fn test_allele_variant_unknown_phase() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let allele = AlleleVariant::unknown_phase(vec![var1, var2]);
        assert_eq!(allele.phase, AllelePhase::Unknown);
    }

    // ============== HgvsVariant Helper Tests ==============

    #[test]
    fn test_hgvs_variant_accession() {
        let variant = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let acc = variant.accession().expect("Expected accession");
        assert_eq!(&*acc.prefix, "NC");
        assert_eq!(&*acc.number, "000001");
    }

    #[test]
    fn test_hgvs_variant_is_allele() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        assert!(!var1.is_allele());

        let allele = AlleleVariant::cis(vec![var1.clone()]);
        let allele_variant = HgvsVariant::Allele(allele);
        assert!(allele_variant.is_allele());
    }

    #[test]
    fn test_hgvs_variant_protein_type() {
        use crate::hgvs::location::{AminoAcid, ProtPos};

        let protein = HgvsVariant::Protein(ProteinVariant {
            accession: Accession::new("NP", "000079", Some(2)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                ProtInterval::point(ProtPos::new(AminoAcid::Val, 600)),
                ProteinEdit::Substitution {
                    reference: AminoAcid::Val,
                    alternative: AminoAcid::Glu,
                },
            ),
        });
        assert_eq!(protein.variant_type(), "p");
        assert!(matches!(protein, HgvsVariant::Protein(_)));

        let genomic = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        assert_eq!(genomic.variant_type(), "g");
        assert!(!matches!(genomic, HgvsVariant::Protein(_)));
    }

    #[test]
    fn test_hgvs_variant_variant_type() {
        let genomic = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        assert_eq!(genomic.variant_type(), "g");

        use crate::hgvs::location::CdsPos;
        let cds = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        assert_eq!(cds.variant_type(), "c");
    }

    #[test]
    fn test_cds_variant_with_gene_symbol() {
        use crate::hgvs::location::CdsPos;

        let variant = CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };

        // Per #121: Display preserves the gene-symbol selector when set.
        assert_eq!(variant.gene_symbol, Some("COL1A1".to_string()));
        let display = format!("{}", variant);
        assert_eq!(display, "NM_000088.3(COL1A1):c.100A>G");
    }

    #[test]
    fn test_genome_variant_with_gene_symbol() {
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: Some("BRCA1".to_string()),
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };

        // Per #121: Display preserves the gene-symbol selector when set.
        assert_eq!(variant.gene_symbol, Some("BRCA1".to_string()));
        let display = format!("{}", variant);
        assert_eq!(display, "NC_000001.11(BRCA1):g.12345A>G");
    }

    #[test]
    fn test_transcript_accession_bare() {
        let acc = Accession::new("NM", "004119", Some(3));
        assert_eq!(acc.transcript_accession(), "NM_004119.3");
        assert_eq!(acc.full(), "NM_004119.3");
    }

    #[test]
    fn test_transcript_accession_with_genomic_context() {
        let inner = Accession::new("NM", "004119", Some(3));
        let outer = Accession::new("NC", "000013", Some(11));
        let compound = inner.with_genomic_context(outer);

        // full() includes the compound wrapper
        assert_eq!(compound.full(), "NC_000013.11(NM_004119.3)");
        // transcript_accession() returns just the bare inner accession
        assert_eq!(compound.transcript_accession(), "NM_004119.3");
    }

    #[test]
    fn test_transcript_accession_assembly_ref() {
        let acc = Accession::from_assembly("GRCh38", "chr1");
        assert_eq!(acc.transcript_accession(), "GRCh38(chr1)");
    }

    #[test]
    fn test_transcript_accession_ensembl_with_context() {
        let inner = Accession::with_style("ENST", "00000241453", Some(7), true);
        let outer = Accession::new("NC", "000013", Some(11));
        let compound = inner.with_genomic_context(outer);

        assert_eq!(compound.full(), "NC_000013.11(ENST00000241453.7)");
        assert_eq!(compound.transcript_accession(), "ENST00000241453.7");
    }

    // ----- Issue #115: self-cancelling allele detection (E3006) -----

    use crate::hgvs::parser::variant::parse_variant;

    #[test]
    fn test_self_cancelling_detect_spec_example() {
        // Spec example: c.[762_768del;767_774dup]
        let del = parse_variant("NM_004006.2:c.762_768del").unwrap();
        let dup = parse_variant("NM_004006.2:c.767_774dup").unwrap();
        let pair = AlleleVariant::detect_self_cancelling_pair(&[del, dup]);
        assert!(
            pair.is_some(),
            "spec-example overlapping del+dup must be detected"
        );
        let (dl, du) = pair.unwrap();
        assert_eq!(dl, 0, "first index is the del");
        assert_eq!(du, 1, "second index is the dup");
    }

    #[test]
    fn test_self_cancelling_no_overlap() {
        let del = parse_variant("NM_004006.2:c.100_110del").unwrap();
        let dup = parse_variant("NM_004006.2:c.500_510dup").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none());
    }

    #[test]
    fn test_self_cancelling_two_dels() {
        let a = parse_variant("NM_004006.2:c.100_110del").unwrap();
        let b = parse_variant("NM_004006.2:c.105_115del").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[a, b]).is_none());
    }

    #[test]
    fn test_self_cancelling_different_accessions() {
        let del = parse_variant("NM_004006.2:c.100_110del").unwrap();
        let dup = parse_variant("NM_000088.3:c.105_115dup").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none());
    }

    #[test]
    fn test_self_cancelling_dup_first_then_del() {
        // Order in the allele shouldn't matter; the dup-first variant is
        // still a self-cancelling pair.
        let dup = parse_variant("NM_004006.2:c.767_774dup").unwrap();
        let del = parse_variant("NM_004006.2:c.762_768del").unwrap();
        let pair = AlleleVariant::detect_self_cancelling_pair(&[dup, del]);
        assert!(pair.is_some());
        let (dl, du) = pair.unwrap();
        assert_eq!(dl, 1, "del is at index 1");
        assert_eq!(du, 0, "dup is at index 0");
    }

    #[test]
    fn test_self_cancelling_genomic_overlap() {
        let del = parse_variant("NC_000017.11:g.43044295_43044300del").unwrap();
        let dup = parse_variant("NC_000017.11:g.43044298_43044305dup").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_some());
    }

    // ----- Issue #172: extend E3006 to *-region and intronic ranges -----

    #[test]
    fn test_self_cancelling_cds_3utr_star_overlap() {
        // 3'UTR (`c.*N`) range; HGVS general.md line 47 forbids overlapping
        // del+dup just as much in *-region as in the CDS proper.
        let del = parse_variant("NM_004006.2:c.*100_*150del").unwrap();
        let dup = parse_variant("NM_004006.2:c.*145_*160dup").unwrap();
        let pair = AlleleVariant::detect_self_cancelling_pair(&[del, dup]);
        assert!(
            pair.is_some(),
            "overlapping 3'UTR del+dup must be detected (#172)"
        );
    }

    #[test]
    fn test_self_cancelling_cds_3utr_star_no_overlap() {
        // Non-overlapping 3'UTR pair: spec-permitted.
        let del = parse_variant("NM_004006.2:c.*100_*110del").unwrap();
        let dup = parse_variant("NM_004006.2:c.*200_*210dup").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none());
    }

    #[test]
    fn test_self_cancelling_cds_intronic_offset_overlap() {
        // Intronic offsets on a CDS-relative axis (`c.100+5`).
        let del = parse_variant("NM_004006.2:c.100+5_100+15del").unwrap();
        let dup = parse_variant("NM_004006.2:c.100+12_100+20dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_some(),
            "overlapping intronic del+dup must be detected (#172)"
        );
    }

    #[test]
    fn test_self_cancelling_cds_intronic_offset_no_overlap() {
        // Same anchor, non-overlapping intronic offsets.
        let del = parse_variant("NM_004006.2:c.100+5_100+10del").unwrap();
        let dup = parse_variant("NM_004006.2:c.100+15_100+20dup").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none());
    }

    #[test]
    fn test_self_cancelling_tx_intronic_overlap() {
        // Tx axis (`n.`) intronic case.
        let del = parse_variant("NR_001234.1:n.100+5_100+15del").unwrap();
        let dup = parse_variant("NR_001234.1:n.100+12_100+20dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_some(),
            "overlapping tx intronic del+dup must be detected (#172)"
        );
    }

    #[test]
    fn test_self_cancelling_cds_cross_region_no_overlap() {
        // 5'UTR del (negative base) and 3'UTR dup (`*` region) are on
        // disjoint regions even if their numeric bases coincide; must
        // NOT be reported as overlapping.
        let del = parse_variant("NM_004006.2:c.-50_-40del").unwrap();
        let dup = parse_variant("NM_004006.2:c.*40_*50dup").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none());
    }

    #[test]
    fn test_parse_self_cancelling_cds_3utr_rejected_via_parser() {
        // End-to-end: the parser-level validator must reject the same
        // overlapping 3'UTR pair with E3006.
        let result = parse_variant("NM_004006.2:c.[*100_*150del;*145_*160dup]");
        let err = result.expect_err("3'UTR overlapping del+dup must reject");
        assert_eq!(
            err.code(),
            Some(crate::error::ErrorCode::SelfCancellingAllele)
        );
    }

    #[test]
    fn test_parse_self_cancelling_cds_intronic_rejected_via_parser() {
        let result = parse_variant("NM_004006.2:c.[100+5_100+15del;100+12_100+20dup]");
        let err = result.expect_err("intronic overlapping del+dup must reject");
        assert_eq!(
            err.code(),
            Some(crate::error::ErrorCode::SelfCancellingAllele)
        );
    }

    #[test]
    fn test_self_cancelling_cds_5utr_overlap() {
        // 5'UTR positions (negative `base`) overlap by lex order on
        // `(false, base, 0)`; positive-result coverage for the cross-region
        // negative test above.
        let del = parse_variant("NM_004006.2:c.-50_-30del").unwrap();
        let dup = parse_variant("NM_004006.2:c.-35_-20dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_some(),
            "overlapping 5'UTR del+dup must be detected (#172)"
        );
    }

    #[test]
    fn test_self_cancelling_cds_negative_intronic_offset_overlap() {
        // `c.X-N` offsets are on the same intron as `c.(X-1)+M`, and
        // lex-order on `(false, X, -N)` correctly places them.
        let del = parse_variant("NM_004006.2:c.100-15_100-5del").unwrap();
        let dup = parse_variant("NM_004006.2:c.100-10_100-2dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_some(),
            "overlapping negative-intronic-offset del+dup must be detected"
        );
    }

    #[test]
    fn test_self_cancelling_cds_cross_exon_intronic_range_ok() {
        // Two non-overlapping ranges on opposite sides of the intron
        // midpoint do NOT collide.
        let del = parse_variant("NM_004006.2:c.99+1_99+10del").unwrap();
        let dup = parse_variant("NM_004006.2:c.100-10_100-1dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none(),
            "non-overlapping cross-exon intronic ranges must NOT collide"
        );
    }

    #[test]
    fn test_self_cancelling_rna_intronic_overlap() {
        // r. axis intronic overlap — same shape as c.X+N but with the RNA
        // coordinate system. RNA descriptions are canonically lowercase.
        let del = parse_variant("NM_004006.2:r.100+5_100+15del").unwrap();
        let dup = parse_variant("NM_004006.2:r.100+12_100+20dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_some(),
            "overlapping r. intronic del+dup must be detected (#172)"
        );
    }

    #[test]
    fn test_self_cancelling_tx_downstream_star_overlap() {
        // `n.*N` downstream positions are the tx-axis analogue of c.*N.
        let del = parse_variant("NR_001234.1:n.*100_*150del").unwrap();
        let dup = parse_variant("NR_001234.1:n.*145_*160dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_some(),
            "overlapping n. downstream del+dup must be detected (#172)"
        );
    }

    #[test]
    fn test_self_cancelling_unknown_base_skipped() {
        // `c.?` (unknown base) makes overlap undecidable; the detector
        // must skip the variant rather than emit a false positive.
        let del = parse_variant("NM_004006.2:c.?_200del").unwrap();
        let dup = parse_variant("NM_004006.2:c.150_160dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none(),
            "pair containing c.? must be skipped (undecidable)"
        );
    }

    #[test]
    fn test_self_cancelling_unknown_offset_sentinel_skipped() {
        // `c.X+?` parses to an offset sentinel (`i64::MAX`); treat it as
        // undecidable rather than letting the sentinel poison the lex
        // order.
        let del = parse_variant("NM_004006.2:c.100+?_300del").unwrap();
        let dup = parse_variant("NM_004006.2:c.150_160dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none(),
            "pair containing `+?` offset sentinel must be skipped"
        );
    }

    #[test]
    fn test_self_cancelling_mixed_axis_same_accession_no_false_positive() {
        // CodeRabbit #371: an accession (NM_*) carries both `c.` (CDS) and
        // `n.` / `r.` (tx / RNA) coordinate systems. A `c.` del and an
        // `n.` dup with the same numeric range live in disjoint coordinate
        // continua and MUST NOT be reported as self-cancelling, even
        // though they share `acc.full()`. The detector must gate on axis
        // identity before comparing canonical ranges.
        let del_cds = parse_variant("NM_004006.2:c.100_150del").unwrap();
        let dup_tx = parse_variant("NM_004006.2:n.100_150dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del_cds, dup_tx]).is_none(),
            "c. del and n. dup with matching numeric ranges must NOT collide \
             (different coordinate axes)"
        );

        // Cover the c./r. pairing too — same accession, different axis.
        let del_cds2 = parse_variant("NM_004006.2:c.100_150del").unwrap();
        let dup_rna = parse_variant("NM_004006.2:r.100_150dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del_cds2, dup_rna]).is_none(),
            "c. del and r. dup with matching numeric ranges must NOT collide \
             (different coordinate axes)"
        );

        // Cover the n./r. pairing — both `region`-bearing axes but
        // numerically different coordinate systems.
        let del_tx = parse_variant("NM_004006.2:n.100_150del").unwrap();
        let dup_rna2 = parse_variant("NM_004006.2:r.100_150dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del_tx, dup_rna2]).is_none(),
            "n. del and r. dup with matching numeric ranges must NOT collide \
             (different coordinate axes)"
        );
    }

    #[test]
    fn test_self_cancelling_genomic_vs_mito_no_false_positive() {
        // Genome and mito axes are both `region=false` in the canonical
        // representation; without an axis gate they would collide on
        // identical numeric ranges. Practical exposure is low (the same
        // accession won't appear as both g. and m. in real data), but
        // `detect_self_cancelling_pair` is a public API that operates on
        // a raw slice — synthesize the mismatch directly.
        let del = parse_variant("NC_000017.11:g.43044295_43044300del").unwrap();
        let dup = parse_variant("NC_000017.11:m.43044295_43044300dup").unwrap();
        assert!(
            AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none(),
            "g. and m. variants with matching numeric ranges must NOT collide"
        );
    }
}
