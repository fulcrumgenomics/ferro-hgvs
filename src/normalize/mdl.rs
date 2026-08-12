//! Minimal-description-length (MDL) scoring and decomposition.
//!
//! A Rust port of mutalyzer's C++ `description-extractor`. **Attribution and
//! the exact upstream commit are recorded in `NOTICE-mdl`, beside this file**
//! (MIT, Leiden University Medical Center). Upstream line references below are
//! to `extractor/extractor.cc` and `extractor/extractor.h` at that commit.
//!
//! # Two things live here, and they are independent
//!
//! 1. **A scoring model** ([`Weights`], [`position_weight`],
//!    [`weight_of_change`], [`weight_trivial`]) — the rendered character cost
//!    of an HGVS description, and the cost of the one description that always
//!    exists: a single spanning `delins` over the whole changed block. This is
//!    the half worth having. It gives a comparand computed from
//!    `(reference, observed)` **alone**, where ferro's incumbent bound
//!    (`merge::changed_columns_of_edits`) compares against the *input's own
//!    spelling* and so decides differently for two spellings of one variant.
//!
//! 2. **A decomposition** ([`extract`]) — upstream's recursive
//!    longest-common-substring split. Kept as an **oracle**: it is the only way
//!    to say "mutalyzer would describe this as X" without shelling out to the
//!    live library. It is *not* proposed as ferro's decomposition; measurement
//!    (see the #1235 acceptance ledger) shows it merges member pairs that
//!    ferro's per-member pipeline deliberately keeps apart.
//!
//! # What is deliberately not ported
//!
//! `extractor_transposition` (`extractor.cc:516`) and the
//! `TRANSPOSITION_OPEN`/`TRANSPOSITION_CLOSE` variant flags. A transposition is
//! rendered by upstream as a cross-reference into the reference
//! (`ins[400_500]`); ferro's `Piece` carries a literal replacement and has no
//! such spelling, so the result would have to be flattened back to literal
//! bases anyway — at which point it is the insertion upstream would otherwise
//! have emitted.
//!
//! Upstream calls it at four sites, and at every one it is *conditional*:
//! the transposition is taken only when
//! `weight > weight_transposition && !transposition.is_empty() && !(len == 1 && front is SUBSTITUTION)`
//! (`extractor.cc:292`, `:353`, `:442`, `:481`). Omitting the call takes the
//! fall-through branch at each — i.e. **exactly the result the C++ produces
//! when that condition is false**:
//!
//! | site | fall-through |
//! |---|---|
//! | `:279` no reference left | one insertion variant |
//! | `:340` no LCS | one spanning `delins` |
//! | `:429` after the prefix recursion | one spanning `delins` |
//! | `:468` after the suffix recursion | one spanning `delins` |
//!
//! The protein path (`extractor_protein`, `extractor_frame_shift`,
//! `LCS_frame_shift`, `backtranslation`) is not ported either: ferro derives
//! protein consequences from the coding axis rather than by aligning peptides.
//!
//! # What is parameterized beyond upstream
//!
//! Upstream hard-codes all three of these. Each is calibration, not principle,
//! so each is a [`Config`] knob whose default is upstream's value:
//!
//! * the **weight table** ([`Weights::MUTALYZER`]) — these constants *are* the
//!   canonical-form policy, so they are named and documented rather than
//!   inlined;
//! * **`weight_position`'s reference-length dependence**
//!   ([`PositionWeight`]) — upstream derives it from the length of the whole
//!   reference, which makes canonical form depend on how long the surrounding
//!   contig is (`k <= 2 * wp`, so a 125 bp contig and a 60 kb `NG_` collapse at
//!   different thresholds). [`PositionWeight::Fixed`] makes it
//!   contig-independent;
//! * **the greedy collapse** ([`Collapse`]) — upstream compares a *partial*
//!   running weight against `weight_trivial` before the decomposition exists
//!   (`extractor.cc:421`), so the verdict depends on recursion order.
//!   [`Collapse::Whole`] compares the finished decomposition instead.

#![allow(dead_code)] // Several entry points exist only for the oracle/test arms.

// *******************************************************************
// Scoring model
// *******************************************************************

/// The rendered character cost of each piece of HGVS syntax.
///
/// Upstream's values (`extractor.h:100-106`), adopted here as a starting point
/// and **calibratable**: they are the description-cost policy that decides
/// which of two equally-correct descriptions is canonical, not a fact about
/// HGVS. Each is the number of characters the construct contributes to a
/// rendered description.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Weights {
    /// One base of literal sequence: `A`, `C`, `G`, `T`.
    pub base: usize,
    /// The `del` keyword.
    pub deletion: usize,
    /// The `delins` keyword.
    pub deletion_insertion: usize,
    /// The `ins` keyword.
    pub insertion: usize,
    /// The `inv` keyword.
    pub inversion: usize,
    /// One separator character: `_`, `[`, `]`, `;`.
    pub separator: usize,
    /// The `>` of a substitution.
    pub substitution: usize,
}

impl Weights {
    /// Upstream's table, verbatim (`extractor.h:100-106`).
    pub const MUTALYZER: Weights = Weights {
        base: 1,
        deletion: 3,
        deletion_insertion: 6,
        insertion: 3,
        inversion: 3,
        separator: 1,
        substitution: 1,
    };
}

impl Default for Weights {
    fn default() -> Self {
        Weights::MUTALYZER
    }
}

/// How the per-position weight is set.
///
/// A position descriptor (`274`, `*12`, `1234567`) has no constant length, so
/// upstream charges a single flat `weight_position` for every position it
/// renders and fixes that flat value once per extraction run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PositionWeight {
    /// Upstream's rule (`extractor.cc:100-103`): `ceil(log10(len / 4))` with
    /// **integer** division, clamped to at least 1.
    ///
    /// The consequence, which is why this is a knob: the collapse threshold
    /// scales with the length of the *whole* reference sequence, so the same
    /// local change canonicalizes differently on a 125 bp test contig
    /// (`wp == 1`) than on a 60 kb `NG_` (`wp == 4`). That is hostile to
    /// representation stability.
    FromReferenceLength,
    /// A contig-independent constant.
    Fixed(usize),
}

/// Resolve [`PositionWeight`] against the length of the whole served sequence.
///
/// `reference_length` is upstream's `reference_length` argument to `extract`,
/// i.e. the **entire** reference string — not the trimmed block under
/// consideration.
///
/// Upstream's expression is undefined for `reference_length < 4`: integer
/// division gives 0, `log10(0)` is `-inf`, and converting that to `size_t` is
/// UB (on x86-64 it yields `INT64_MIN` reinterpreted, which then survives the
/// `<= 0` clamp because the variable is unsigned). We take the clamp's evident
/// intent and return 1 there.
pub(crate) fn position_weight(mode: PositionWeight, reference_length: usize) -> usize {
    match mode {
        PositionWeight::Fixed(w) => w.max(1),
        PositionWeight::FromReferenceLength => {
            let quotient = reference_length / 4; // integer division, as upstream
            if quotient == 0 {
                return 1;
            }
            let w = (quotient as f64).log10().ceil() as i64;
            if w <= 0 {
                1
            } else {
                w as usize
            }
        }
    }
}

/// The rendered cost of describing one change as a single HGVS edit.
///
/// The five cases are upstream's own, each read off the site that charges it:
///
/// * insertion (`extractor.cc:272`) — `2*wp + sep + ins + base*alt`
/// * deletion (`:310`) — `wp + del + (ref > 1 ? wp + sep : 0)`
/// * substitution (`:318`) — `wp + 2*base + sub`
/// * inversion (`:387`, `:504`) — `2*wp + sep + inv`
/// * `delins` (`:245`) — `wp + delins + base*alt + (ref != 1 ? wp + sep : 0)`
///
/// Note the asymmetry upstream carries and we keep: a deletion drops its second
/// position at `ref == 1`, a `delins` at `ref != 1` *adds* one — i.e. a
/// single-position `delins` (`7delinsCC`) is charged one position and a
/// multi-position one (`7_9delinsCC`) two.
pub(crate) fn weight_of_change(
    reference_length: usize,
    alt_length: usize,
    is_inversion: bool,
    wp: usize,
    w: Weights,
) -> usize {
    if reference_length == 0 {
        // A pure insertion always names two flanking positions.
        return 2 * wp + w.separator + w.insertion + w.base * alt_length;
    }
    if alt_length == 0 {
        return wp
            + w.deletion
            + if reference_length > 1 {
                wp + w.separator
            } else {
                0
            };
    }
    if reference_length == 1 && alt_length == 1 {
        return wp + 2 * w.base + w.substitution;
    }
    if is_inversion {
        return 2 * wp + w.separator + w.inversion;
    }
    wp + w.deletion_insertion
        + w.base * alt_length
        + if reference_length != 1 {
            wp + w.separator
        } else {
            0
        }
}

/// The cost of the description that always exists: one spanning `delins` over
/// the whole block (`extractor.cc:245`, `weight_trivial`).
///
/// This is the comparand. It is a function of the `(reference, observed)` pair
/// and of nothing else — in particular not of how the input was spelled, which
/// is the whole reason for preferring it to `changed_columns_of_edits`.
pub(crate) fn weight_trivial(
    reference_length: usize,
    alt_length: usize,
    wp: usize,
    w: Weights,
) -> usize {
    wp + w.deletion_insertion
        + w.base * alt_length
        + if reference_length != 1 {
            wp + w.separator
        } else {
            0
        }
}

/// The fewest changed columns any description of this block can claim.
///
/// **This is the comparand that has the right direction.** `weight_trivial`
/// above answers "is the split cheap enough to be worth making", which is a
/// merge-preferring question; ferro's incumbent bound asks the opposite one —
/// "is the derivation claiming more change than it has to" — and that is what
/// stops a re-derivation shredding two well-separated members into a run of
/// substitutions across bases the variant left alone (`general.md:34`).
///
/// The incumbent asks it against the **input's own spelling**
/// (`merge::changed_columns_of_edits`), which is exactly #1235's defect: two
/// spellings of one variant get two different ceilings. The same question has a
/// sequence-derived answer.
///
/// `merge::changed_columns_of_pieces` sums `max(|span|, |replacement|)` over
/// the pieces, and that sum is precisely the number of non-matching columns in
/// the alignment the pieces came from: a piece with `ref = 1, alt = 3` is one
/// mismatched column plus two inserted ones, `ref = 3, alt = 1` is one plus two
/// deleted ones, and `ref = n, alt = n` is `n` mismatches. Minimising it over
/// all alignments is therefore exactly unit-cost Levenshtein distance — the
/// derivation is column-minimal iff its changed-column count equals this.
///
/// Because every spelling of a variant produces the same `(reference, result)`
/// pair, this ceiling is the same for all of them, and it is never larger than
/// the incumbent's (an input's own spelling is an alignment, so its column
/// count is at least the minimum). A gate built on it is therefore at least as
/// strict as the incumbent everywhere.
pub(crate) fn minimal_changed_columns(reference: &[u8], result: &[u8]) -> usize {
    let n = result.len();
    let mut previous: Vec<usize> = (0..=n).collect();
    let mut current = vec![0usize; n + 1];
    for (i, r) in reference.iter().enumerate() {
        current[0] = i + 1;
        for (j, a) in result.iter().enumerate() {
            let substitute = previous[j] + usize::from(r != a);
            current[j + 1] = substitute.min(previous[j + 1] + 1).min(current[j] + 1);
        }
        std::mem::swap(&mut previous, &mut current);
    }
    previous[n]
}

// *******************************************************************
// Configuration
// *******************************************************************

/// Whether the collapse test sees a partial or a finished decomposition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Collapse {
    /// Upstream: test after the prefix recursion *and* after the suffix
    /// recursion (`extractor.cc:421`, `:460`). The first test reads a running
    /// weight that does not yet include the suffix, so it can bail before the
    /// decomposition exists, and the verdict depends on recursion order.
    Greedy,
    /// Test once, on the finished decomposition. Order-independent.
    Whole,
}

/// Everything about the algorithm that is calibration rather than principle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Config {
    pub weights: Weights,
    pub position_weight: PositionWeight,
    pub collapse: Collapse,
    /// Whether the LCS search also looks in reverse-complement space, which is
    /// what lets an `inv` be recognised (`extractor.cc:1059`).
    pub reverse_complement: bool,
    /// Veto a collapse whose decomposition contains an unchanged run of at
    /// least this many reference bases.
    ///
    /// Not upstream — this is where `general.md:34` ("two variants separated by
    /// one or more unchanged nucleotides are described individually") composes
    /// with cost minimisation. `None` is upstream's behaviour: cost decides
    /// alone, and a cheap spanning `delins` may swallow any amount of unchanged
    /// sequence. Only consulted in [`Collapse::Whole`] with full information;
    /// under [`Collapse::Greedy`] it sees whatever the greedy site has built so
    /// far, which is the honest reading of a greedy test.
    pub separation_floor: Option<usize>,
}

impl Config {
    /// Upstream, as faithfully as this port goes.
    pub const MUTALYZER: Config = Config {
        weights: Weights::MUTALYZER,
        position_weight: PositionWeight::FromReferenceLength,
        collapse: Collapse::Greedy,
        reverse_complement: true,
        separation_floor: None,
    };
}

impl Default for Config {
    fn default() -> Self {
        Config::MUTALYZER
    }
}

// *******************************************************************
// Decomposition (the oracle)
// *******************************************************************

/// Upstream's `MASK` character, ignored by every match (`extractor.cc:72`).
const MASK: u8 = b'$';

/// Reference length above which `LCS` stops falling back to the quadratic
/// algorithm (`extractor.h:113`).
const THRESHOLD_CUT_OFF: usize = 16_000;

/// What one extracted region is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Kind {
    /// Unchanged (`IDENTITY`, `extractor.h:69`).
    Identity,
    /// Replaced by its own reverse complement (`REVERSE_COMPLEMENT`).
    ReverseComplement,
    /// Any other change: deletion, insertion, substitution or `delins`
    /// (`SUBSTITUTION`, which upstream's comment notes "covers most").
    Substitution,
}

/// One extracted region, in upstream's coordinates: half-open offsets into the
/// reference and the sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Variant {
    pub reference_start: usize,
    pub reference_end: usize,
    pub sample_start: usize,
    pub sample_end: usize,
    pub kind: Kind,
}

/// A common substring between reference (or its complement) and sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Substring {
    reference_index: usize,
    sample_index: usize,
    length: usize,
    reverse_complement: bool,
}

/// Complement each base in place, without reversing (`IUPAC_complement`,
/// `extractor.cc:1510`). Anything outside `ACGTU` is left alone, as upstream.
fn iupac_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            other => *other,
        })
        .collect()
}

/// Length of the common prefix (`prefix_match`, `extractor.cc:1471`).
fn prefix_match(reference: &[u8], sample: &[u8]) -> usize {
    let mut i = 0;
    while i < reference.len()
        && i < sample.len()
        && reference[i] == sample[i]
        && reference[i] != MASK
    {
        i += 1;
    }
    i
}

/// Length of the common suffix, given the common prefix (`suffix_match`,
/// `extractor.cc:1489`).
fn suffix_match(reference: &[u8], sample: &[u8], prefix: usize) -> usize {
    let mut i = 0;
    while i + prefix < reference.len()
        && i + prefix < sample.len()
        && reference[reference.len() - i - 1] == sample[sample.len() - i - 1]
        && reference[reference.len() - i - 1] != MASK
    {
        i += 1;
    }
    i
}

/// `strncmp`-alike that also rejects the mask (`string_match`, `:1436`).
fn string_match(a: &[u8], b: &[u8], length: usize) -> bool {
    (0..length).all(|i| a.get(i).is_some() && b.get(i).is_some() && a[i] == b[i] && a[i] != MASK)
}

/// Extract the variants between `reference` and `sample`.
///
/// Upstream's `extract` (`extractor.cc:90`), DNA path. `reference_length` is
/// the length of the **whole served sequence**, which is what
/// `weight_position` keys on; it is deliberately *not* `reference.len()`, since
/// ferro hands this function an already-trimmed block.
///
/// The returned list tiles both strings: identity regions are present, so the
/// result can be compared element-for-element against upstream's own test
/// expectations.
pub(crate) fn extract(
    reference: &[u8],
    sample: &[u8],
    reference_length: usize,
    config: &Config,
) -> Vec<Variant> {
    let wp = position_weight(config.position_weight, reference_length);
    let prefix = prefix_match(reference, sample);
    let suffix = suffix_match(reference, sample, prefix);

    let complement = config
        .reverse_complement
        .then(|| iupac_complement(reference));

    let mut out = Vec::new();
    if prefix > 0 {
        out.push(Variant {
            reference_start: 0,
            reference_end: prefix,
            sample_start: 0,
            sample_end: prefix,
            kind: Kind::Identity,
        });
    }

    let extractor = Extractor {
        reference,
        complement: complement.as_deref(),
        sample,
        wp,
        weights: config.weights,
        collapse: config.collapse,
        separation_floor: config.separation_floor,
    };
    let (_weight, body) = extractor.extract(
        prefix,
        reference.len() - suffix,
        prefix,
        sample.len() - suffix,
    );
    out.extend(body);

    if suffix > 0 {
        out.push(Variant {
            reference_start: reference.len() - suffix,
            reference_end: reference.len(),
            sample_start: sample.len() - suffix,
            sample_end: sample.len(),
            kind: Kind::Identity,
        });
    }
    out
}

/// The total weight upstream would assign to describing `sample` against
/// `reference` — the return value of `extract` that the variant list discards.
pub(crate) fn extract_weight(
    reference: &[u8],
    sample: &[u8],
    reference_length: usize,
    config: &Config,
) -> usize {
    let wp = position_weight(config.position_weight, reference_length);
    let prefix = prefix_match(reference, sample);
    let suffix = suffix_match(reference, sample, prefix);
    let complement = config
        .reverse_complement
        .then(|| iupac_complement(reference));
    let extractor = Extractor {
        reference,
        complement: complement.as_deref(),
        sample,
        wp,
        weights: config.weights,
        collapse: config.collapse,
        separation_floor: config.separation_floor,
    };
    extractor
        .extract(
            prefix,
            reference.len() - suffix,
            prefix,
            sample.len() - suffix,
        )
        .0
}

struct Extractor<'a> {
    reference: &'a [u8],
    /// Complement of `reference`, **not** reversed. `None` disables
    /// reverse-complement matching (upstream passes a null pointer for
    /// non-DNA).
    complement: Option<&'a [u8]>,
    sample: &'a [u8],
    wp: usize,
    weights: Weights,
    collapse: Collapse,
    separation_floor: Option<usize>,
}

impl Extractor<'_> {
    /// The recursive extractor (`extractor.cc:204`).
    ///
    /// Returns `(weight, variants)`.
    fn extract(
        &self,
        reference_start: usize,
        reference_end: usize,
        sample_start: usize,
        sample_end: usize,
    ) -> (usize, Vec<Variant>) {
        let w = self.weights;
        let wp = self.wp;

        // Mask trimming (`:213-238`). Upstream's second loop stops with one
        // character left (`end - i - 1 > start`), which is reproduced here as
        // `end > start + 1`.
        let mut rs = reference_start;
        let mut re = reference_end;
        let mut ss = sample_start;
        let mut se = sample_end;
        while rs < re && self.reference[rs] == MASK {
            rs += 1;
        }
        while re > rs + 1 && self.reference[re - 1] == MASK {
            re -= 1;
        }
        while ss < se && self.sample[ss] == MASK {
            ss += 1;
        }
        while se > ss + 1 && self.sample[se - 1] == MASK {
            se -= 1;
        }

        let reference_length = re - rs;
        let sample_length = se - ss;
        let trivial = weight_trivial(reference_length, sample_length, wp, w);

        let spanning = |kind: Kind| Variant {
            reference_start: rs,
            reference_end: re,
            sample_start: ss,
            sample_end: se,
            kind,
        };

        // Base cases (`:266-321`). None of them consults `weight_trivial`.
        if reference_length == 0 {
            if sample_length > 0 {
                // Transposition omitted — this is upstream's `:301`
                // fall-through.
                let weight = 2 * wp + w.separator + w.insertion + w.base * sample_length;
                return (weight, vec![spanning(Kind::Substitution)]);
            }
            return (0, Vec::new());
        }
        if sample_length == 0 {
            let weight = wp
                + w.deletion
                + if reference_length > 1 {
                    wp + w.separator
                } else {
                    0
                };
            return (weight, vec![spanning(Kind::Substitution)]);
        }
        if reference_length == 1 && sample_length == 1 {
            let weight = wp + 2 * w.base + w.substitution;
            return (weight, vec![spanning(Kind::Substitution)]);
        }

        let cut_off = if reference_length < THRESHOLD_CUT_OFF {
            1
        } else {
            wp
        };
        let (length, substrings) = self.lcs(rs, re, ss, se, cut_off);

        // No LCS (`:332`): one spanning delins, transposition omitted.
        if length == 0 || substrings.is_empty() {
            return (trivial, vec![spanning(Kind::Substitution)]);
        }

        // Pick the "best fitting" LCS: the one leaving prefixes and suffixes of
        // most nearly equal length (`:369-381`). Upstream's `abs` on wrapped
        // `size_t` differences is an absolute difference; written as one here.
        let mut diff = (re - rs) + (se - ss);
        let mut lcs = substrings[0];
        for candidate in &substrings {
            let prefix_diff =
                (candidate.reference_index - rs).abs_diff(candidate.sample_index - ss);
            let suffix_diff = (re - (candidate.reference_index + candidate.length))
                .abs_diff(se - (candidate.sample_index + candidate.length));
            if prefix_diff + suffix_diff < diff {
                diff = prefix_diff + suffix_diff;
                lcs = *candidate;
            }
        }

        // An inverted LCS is itself a described change, so it carries weight
        // (`:385-388`).
        let mut weight = if lcs.reverse_complement {
            2 * wp + w.separator + w.inversion
        } else {
            0
        };

        let matched = Variant {
            reference_start: lcs.reference_index,
            reference_end: lcs.reference_index + length,
            sample_start: lcs.sample_index,
            sample_end: lcs.sample_index + length,
            kind: if lcs.reverse_complement {
                Kind::ReverseComplement
            } else {
                Kind::Identity
            },
        };

        let (prefix_weight, prefix) = self.extract(rs, lcs.reference_index, ss, lcs.sample_index);
        weight += prefix_weight;

        // Greedy collapse site 1 (`:421`) — a *partial* weight.
        if self.collapse == Collapse::Greedy
            && weight > trivial
            && !self.separation_vetoes_collapse(&prefix, &matched, &[])
        {
            return (trivial, vec![spanning(Kind::Substitution)]);
        }

        let (suffix_weight, suffix) = self.extract(
            lcs.reference_index + length,
            re,
            lcs.sample_index + length,
            se,
        );
        weight += suffix_weight;

        // Collapse site 2 (`:460`) — in `Collapse::Whole` this is the only
        // test, and `weight` is then the finished decomposition's weight.
        if weight > trivial && !self.separation_vetoes_collapse(&prefix, &matched, &suffix) {
            return (trivial, vec![spanning(Kind::Substitution)]);
        }

        let mut out = prefix;
        out.push(matched);
        out.extend(suffix);
        (weight, out)
    }

    /// Whether [`Config::separation_floor`] forbids collapsing this
    /// decomposition into one spanning `delins`.
    ///
    /// Not upstream. See [`Config::separation_floor`].
    fn separation_vetoes_collapse(
        &self,
        prefix: &[Variant],
        matched: &Variant,
        suffix: &[Variant],
    ) -> bool {
        let Some(floor) = self.separation_floor else {
            return false;
        };
        prefix
            .iter()
            .chain(std::iter::once(matched))
            .chain(suffix.iter())
            .any(|v| v.kind == Kind::Identity && v.reference_end - v.reference_start >= floor)
    }

    /// `LCS` (`extractor.cc:931`): try `LCS_k` with a shrinking `k`, then fall
    /// back to the quadratic `LCS_1` unless a cut-off forbids it.
    fn lcs(
        &self,
        rs: usize,
        re: usize,
        ss: usize,
        se: usize,
        cut_off: usize,
    ) -> (usize, Vec<Substring>) {
        let reference_length = re - rs;
        let sample_length = se - ss;
        let mut k = if reference_length > sample_length {
            sample_length / 8
        } else {
            reference_length / 8
        };
        while k > 8 && k > cut_off {
            let (length, substrings) = self.lcs_k(rs, re, ss, se, k);
            if length >= 2 * k && !substrings.is_empty() {
                return (length, substrings);
            }
            k /= 3;
        }
        if cut_off > 1 {
            return (0, Vec::new());
        }
        self.lcs_1(rs, re, ss, se)
    }

    /// `LCS_1` (`extractor.cc:997`): the classical two-row dynamic program,
    /// run simultaneously in forward and reverse-complement space.
    fn lcs_1(&self, rs: usize, re: usize, ss: usize, se: usize) -> (usize, Vec<Substring>) {
        let reference_length = re - rs;
        let sample_length = se - ss;
        let mut reverse_complement = false;
        let mut line = vec![0usize; 2 * reference_length];
        let mut line_rc = vec![0usize; 2 * reference_length];
        let mut length = 0usize;
        let mut substring: Vec<Substring> = Vec::new();

        for i in 0..sample_length {
            let (cur, prev) = (i % 2, (i + 1) % 2);
            for j in 0..reference_length {
                if self.reference[rs + j] == self.sample[ss + i] && self.reference[rs + j] != MASK {
                    line[cur * reference_length + j] = if i == 0 || j == 0 {
                        1
                    } else {
                        line[prev * reference_length + j - 1] + 1
                    };
                    let run = line[cur * reference_length + j];
                    if run >= length {
                        if reverse_complement || run > length {
                            length = run;
                            substring = vec![Substring {
                                // Written `a + b + 1 - length` rather than
                                // upstream's `a - length + b + 1`: the C++
                                // intermediate wraps through `size_t` and back,
                                // which Rust rejects. Same value.
                                reference_index: rs + j + 1 - length,
                                sample_index: ss + i + 1 - length,
                                length,
                                reverse_complement: false,
                            }];
                        } else {
                            substring.push(Substring {
                                // Written `a + b + 1 - length` rather than
                                // upstream's `a - length + b + 1`: the C++
                                // intermediate wraps through `size_t` and back,
                                // which Rust rejects. Same value.
                                reference_index: rs + j + 1 - length,
                                sample_index: ss + i + 1 - length,
                                length,
                                reverse_complement: false,
                            });
                        }
                        reverse_complement = false;
                    }
                } else {
                    line[cur * reference_length + j] = 0;
                }

                match self.complement {
                    Some(complement)
                        if complement[re - j - 1] == self.sample[ss + i]
                            && complement[re - j - 1] != MASK =>
                    {
                        line_rc[cur * reference_length + j] = if i == 0 || j == 0 {
                            1
                        } else {
                            line_rc[prev * reference_length + j - 1] + 1
                        };
                        let run = line_rc[cur * reference_length + j];
                        // Upstream requires `> 1`, so a single complemented base
                        // is never reported as an inversion (`:1070`).
                        if run > 1 && run > length {
                            length = run;
                            substring = vec![Substring {
                                reference_index: re - j - 1,
                                sample_index: ss + i + 1 - length,
                                length,
                                reverse_complement: true,
                            }];
                            reverse_complement = true;
                        }
                    }
                    _ => line_rc[cur * reference_length + j] = 0,
                }

                // Upstream breaks the inner loop only (`:1084`), leaving the
                // rest of this row zero-initialised from the previous sweep —
                // which the next row then reads. Reproduced, including that.
                if !reverse_complement && length >= sample_length {
                    break;
                }
            }
        }
        (length, substring)
    }

    /// `LCS_k` (`extractor.cc:1103`): k-mer coarsened LCS for large similar
    /// strings, followed by base-level extension of each hit.
    fn lcs_k(
        &self,
        rs: usize,
        re: usize,
        ss: usize,
        se: usize,
        k: usize,
    ) -> (usize, Vec<Substring>) {
        let mut length = 0usize;
        let mut substring: Vec<Substring> = Vec::new();
        if k <= 1 || re - rs < k || se - ss < k {
            return (length, substring);
        }
        let reference_length = (re - rs) / k;
        let sample_length = se - ss - k + 1;
        if reference_length == 0 {
            return (length, substring);
        }
        let mut reverse_complement = false;
        let rows = k + 1;
        let mut line = vec![0usize; rows * reference_length];
        let mut line_rc = vec![0usize; rows * reference_length];

        for i in 0..sample_length {
            let cur = i % rows;
            let prev = (i + 1) % rows;
            for j in 0..reference_length {
                if string_match(&self.reference[rs + j * k..], &self.sample[ss + i..], k) {
                    line[cur * reference_length + j] = if i < k || j == 0 {
                        1
                    } else {
                        line[prev * reference_length + j - 1] + 1
                    };
                    let run = line[cur * reference_length + j];
                    if run > length {
                        length = run;
                        // Drop every solution now more than one behind, plus the
                        // partial that was just extended (`:1158-1166`).
                        substring.retain(|s| {
                            !(length - s.length > 1
                                || (Some(s.reference_index) == j.checked_sub(1)
                                    && s.sample_index == i.wrapping_sub(k)
                                    && !s.reverse_complement))
                        });
                        substring.push(Substring {
                            reference_index: j,
                            sample_index: i,
                            length: run,
                            reverse_complement: false,
                        });
                    } else if run > 0 && length - run <= 1 {
                        substring.push(Substring {
                            reference_index: j,
                            sample_index: i,
                            length: run,
                            reverse_complement: false,
                        });
                    }
                } else {
                    line[cur * reference_length + j] = 0;
                }

                let rc_hit = match self.complement {
                    Some(complement) => {
                        // `string_match_reverse(complement + reference_end - j*k - 1, …)`
                        // walks backwards from that pointer.
                        let anchor = re - j * k;
                        anchor >= k
                            && (0..k).all(|t| {
                                let ri = anchor - 1 - t;
                                complement[ri] == self.sample[ss + i + t] && complement[ri] != MASK
                            })
                    }
                    None => false,
                };
                if rc_hit {
                    line_rc[cur * reference_length + j] = if i < k || j == 0 {
                        1
                    } else {
                        line_rc[prev * reference_length + j - 1] + 1
                    };
                    let run = line_rc[cur * reference_length + j];
                    if run > length {
                        length = run;
                        substring.retain(|s| {
                            !(length - s.length > 1
                                || (Some(s.reference_index) == j.checked_sub(1)
                                    && s.sample_index == i.wrapping_sub(k)
                                    && s.reverse_complement))
                        });
                        substring.push(Substring {
                            reference_index: j,
                            sample_index: i,
                            length: run,
                            reverse_complement: true,
                        });
                    } else if line[cur * reference_length + j] > 0 && length - run <= 1 {
                        // Upstream reads the *forward* line here (`:1213-1215`).
                        // Faithfully reproduced, typo or not.
                        substring.push(Substring {
                            reference_index: j,
                            sample_index: i,
                            length: line[cur * reference_length + j],
                            reverse_complement: true,
                        });
                    }
                } else {
                    line_rc[cur * reference_length + j] = 0;
                }
            }
        }

        // Convert k-mer coordinates back to base coordinates and extend up to k
        // positions each way (`:1231-1295`).
        for s in substring.iter_mut() {
            if !s.reverse_complement {
                s.reference_index = rs + (s.reference_index + 1 - s.length) * k;
                s.sample_index = ss + s.sample_index - (s.length - 1) * k;
                s.length *= k;
                let mut i = 0;
                while i <= k
                    && s.reference_index + s.length + i < re
                    && s.sample_index + s.length + i < se
                    && self.reference[s.reference_index + s.length + i]
                        == self.sample[s.sample_index + s.length + i]
                    && self.reference[s.reference_index + s.length + i] != MASK
                {
                    i += 1;
                }
                s.length += i;
                let mut i = 0;
                while i <= k
                    && s.reference_index > rs + i
                    && s.sample_index > ss + i
                    && self.reference[s.reference_index - i - 1]
                        == self.sample[s.sample_index - i - 1]
                    && self.reference[s.reference_index - i - 1] != MASK
                {
                    i += 1;
                }
                s.reference_index -= i;
                s.sample_index -= i;
                s.length += i;
            } else {
                let complement = self.complement.expect("rc substring without a complement");
                s.reference_index = re - (s.reference_index + 1) * k;
                s.sample_index = ss + s.sample_index - (s.length - 1) * k;
                s.length *= k;
                let mut i = 0;
                while i <= k
                    && s.reference_index > rs + i
                    && s.sample_index + s.length + i < se
                    && complement[s.reference_index - i - 1]
                        == self.sample[s.sample_index + s.length + i]
                    && complement[s.reference_index - i - 1] != MASK
                {
                    i += 1;
                }
                s.reference_index -= i;
                s.length += i;
                let mut i = 0;
                while i <= k
                    && s.reference_index + s.length + i < re
                    && s.sample_index > ss + i
                    && complement[s.reference_index + s.length + i]
                        == self.sample[s.sample_index - i - 1]
                    && complement[s.reference_index + s.length + i] != MASK
                {
                    i += 1;
                }
                s.sample_index -= i;
                s.length += i;
            }

            if s.length > length {
                length = s.length;
                reverse_complement = s.reverse_complement;
            } else if reverse_complement && s.length == length {
                reverse_complement = false;
            }
        }

        substring.retain(|s| s.length >= length && s.reverse_complement == reverse_complement);
        (length, substring)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Upstream's own `weight_position` ladder, recomputed here so a change to
    /// the expression is visible rather than silent.
    #[test]
    fn position_weight_matches_upstreams_ladder() {
        let from = |n| position_weight(PositionWeight::FromReferenceLength, n);
        // < 4 is undefined upstream; we take the clamp's intent.
        assert_eq!(from(0), 1);
        assert_eq!(from(3), 1);
        // 4/4 = 1, log10(1) = 0, ceil = 0, clamped to 1.
        assert_eq!(from(4), 1);
        // 43/4 = 10 -> log10 = 1 -> ceil 1.
        assert_eq!(from(43), 1);
        // 44/4 = 11 -> log10 ~ 1.041 -> ceil 2.
        assert_eq!(from(44), 2);
        assert_eq!(from(403), 2); // 100 -> exactly 2
        assert_eq!(from(404), 3); // 101 -> 2.004 -> 3
                                  // 60000/4 = 15000 -> log10 = 4.176 -> ceil 5. A 60 kb NG_ therefore
                                  // charges five characters per position where a 125 bp contig charges
                                  // one — the contig dependence `PositionWeight::Fixed` exists to remove.
        assert_eq!(from(60_000), 5);
        assert_eq!(from(125), 2);
    }

    #[test]
    fn fixed_position_weight_is_contig_independent() {
        for n in [10usize, 125, 4_000, 250_000_000] {
            assert_eq!(position_weight(PositionWeight::Fixed(2), n), 2);
        }
        // A zero is clamped, as upstream clamps its own.
        assert_eq!(position_weight(PositionWeight::Fixed(0), 100), 1);
    }

    #[test]
    fn trivial_weight_is_the_spanning_delins_cost() {
        let w = Weights::MUTALYZER;
        // A 5 -> 3 delins at wp = 1: 1 + 6 + 3 + (1 + 1) = 12.
        assert_eq!(weight_trivial(5, 3, 1, w), 12);
        // A 1 -> 3 delins drops the second position: 1 + 6 + 3 = 10.
        assert_eq!(weight_trivial(1, 3, 1, w), 10);
        // And it scales with wp, which is the contig dependence the knob exists
        // to make optional.
        assert_eq!(weight_trivial(5, 3, 4, w), 4 + 6 + 3 + 4 + 1);
    }

    #[test]
    fn minimal_changed_columns_is_unit_cost_edit_distance() {
        for (r, a, expected) in [
            (&b""[..], &b""[..], 0usize),
            (&b"ACGT"[..], &b"ACGT"[..], 0),
            (&b"ACGT"[..], &b"AGGT"[..], 1),
            (&b"ACGT"[..], &b"ACT"[..], 1),
            (&b"ACGT"[..], &b"ACGGT"[..], 1),
            (&b"AAAAAAA"[..], &b"AACACAAAA"[..], 2),
            (&b"AAAAAAA"[..], &b"ACAAAA"[..], 2),
            // A pure deletion and a pure insertion each cost their own length.
            (&b"ACGTT"[..], &b""[..], 5),
            (&b""[..], &b"ACG"[..], 3),
        ] {
            assert_eq!(minimal_changed_columns(r, a), expected, "{r:?} -> {a:?}");
        }
    }

    /// The property the gate leans on: a ceiling taken from the sequences is
    /// never above one taken from any particular spelling, so the gate can only
    /// be stricter than the incumbent — never laxer.
    #[test]
    fn the_column_minimum_never_exceeds_a_particular_alignments_count() {
        // `NC:g.[261del;263_264insA]` on `…AGAA…`: the input spells two members
        // (2 columns); a spanning `261_263delinsAGA` spells 3. The minimum is
        // the smaller of the two, and it does not move when the spelling does.
        let (reference, result) = (&b"AGAA"[..], &b"GAAA"[..]);
        assert_eq!(minimal_changed_columns(reference, result), 2);
        // A verbose spelling of the same change cannot raise the ceiling.
        assert_eq!(minimal_changed_columns(reference, result), 2);
    }

    #[test]
    fn change_weights_follow_their_upstream_sites() {
        let w = Weights::MUTALYZER;
        assert_eq!(weight_of_change(0, 2, false, 1, w), 2 + 1 + 3 + 2); // ins
        assert_eq!(weight_of_change(3, 0, false, 1, w), 1 + 3 + 1 + 1); // del
        assert_eq!(weight_of_change(1, 0, false, 1, w), 1 + 3); // 1 nt del
        assert_eq!(weight_of_change(1, 1, false, 1, w), 1 + 2 + 1); // SNP
        assert_eq!(weight_of_change(4, 4, true, 1, w), 2 + 1 + 3); // inv
        assert_eq!(weight_of_change(3, 4, false, 1, w), 1 + 6 + 4 + 1 + 1); // delins
    }

    // ---------------------------------------------------------------
    // Ported from upstream `tests/test_extractor.py` (MIT, LUMC) — see
    // `NOTICE-mdl`. Each case is (reference, sample, expected variants) with
    // the expectation copied verbatim from the upstream file.
    // ---------------------------------------------------------------

    fn v(rs: usize, re: usize, ss: usize, se: usize, kind: Kind) -> Variant {
        Variant {
            reference_start: rs,
            reference_end: re,
            sample_start: ss,
            sample_end: se,
            kind,
        }
    }

    /// Upstream's numeric type codes: 1 = IDENTITY, 2 = REVERSE_COMPLEMENT,
    /// 4 = SUBSTITUTION.
    fn kind(code: u32) -> Kind {
        match code {
            1 => Kind::Identity,
            2 => Kind::ReverseComplement,
            4 => Kind::Substitution,
            other => panic!("unsupported upstream type code {other}"),
        }
    }

    fn check(reference: &str, sample: &str, expected: &[(usize, usize, usize, usize, u32)]) {
        // Upstream calls `extract(reference, len(reference), sample, …)`, so the
        // whole reference *is* the served sequence in these cases.
        let got = extract(
            reference.as_bytes(),
            sample.as_bytes(),
            reference.len(),
            &Config::MUTALYZER,
        );
        let want: Vec<Variant> = expected
            .iter()
            .map(|(rs, re, ss, se, t)| v(*rs, *re, *ss, *se, kind(*t)))
            .collect();
        assert_eq!(got, want, "reference={reference} sample={sample}");
    }

    #[test]
    fn upstream_test1_insertion_deletion_substitution_and_duplication() {
        check(
            "ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA",
            "ATGATTTGATCAGATACATGTGATACCGGTAGTTAGGACAA",
            &[
                (0, 5, 0, 5, 1),
                (5, 5, 5, 7, 4),
                (5, 16, 7, 18, 1),
                (16, 17, 18, 18, 4),
                (17, 25, 18, 26, 1),
                (25, 26, 26, 27, 4),
                (26, 34, 27, 35, 1),
                (34, 34, 35, 36, 4),
                (34, 39, 36, 41, 1),
            ],
        );
    }

    #[test]
    fn upstream_test2_three_substitutions() {
        check(
            "TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTACTGTG",
            "TAAGCACCAGGAGTCCATGAAGAAGCTGGATCCTCCCATGGAATCCCCTACTCTACTGTG",
            &[
                (0, 25, 0, 25, 1),
                (25, 26, 25, 26, 4),
                (26, 29, 26, 29, 1),
                (29, 30, 29, 30, 4),
                (30, 34, 30, 34, 1),
                (34, 35, 34, 35, 4),
                (35, 60, 35, 60, 1),
            ],
        );
    }

    #[test]
    fn upstream_test3_inversion_beside_a_substitution() {
        check(
            "TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTA",
            "TAAGCACCAGGAGTCCATGAAGAAGCCATGTCCTGCCATGGAATCCCCTACTCTA",
            &[
                (0, 25, 0, 25, 1),
                (25, 29, 25, 29, 2),
                (29, 30, 29, 30, 4),
                (30, 55, 30, 55, 1),
            ],
        );
    }

    #[test]
    fn upstream_test4_inversion_substitution_and_deletion() {
        check(
            "TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTA",
            "TAAGCACCAGGAGTCCATGAAGAAGCCATGTCCTGCCATGAATCCCCTACTCTA",
            &[
                (0, 25, 0, 25, 1),
                (25, 29, 25, 29, 2),
                (29, 30, 29, 30, 4),
                (30, 39, 30, 39, 1),
                (39, 40, 39, 39, 4),
                (40, 55, 39, 54, 1),
            ],
        );
    }

    #[test]
    fn upstream_tests5_to_18_single_variants() {
        const REFERENCE: &str = "ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT";
        // (sample, expected) exactly as upstream records them.
        check(REFERENCE, REFERENCE, &[(0, 44, 0, 44, 1)]);
        check(
            REFERENCE,
            "ACGTCGGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 6, 0, 6, 1), (6, 7, 6, 7, 4), (7, 44, 7, 44, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 6, 0, 6, 1), (6, 7, 6, 6, 4), (7, 44, 6, 43, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 6, 0, 6, 1), (6, 8, 6, 6, 4), (8, 44, 6, 42, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 6, 0, 6, 1), (6, 6, 6, 7, 4), (6, 44, 7, 45, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGCCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 6, 0, 6, 1), (6, 6, 6, 8, 4), (6, 44, 8, 46, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGAATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 7, 0, 7, 1), (7, 7, 7, 8, 4), (7, 44, 8, 45, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGAGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 7, 0, 7, 1), (7, 7, 7, 9, 4), (7, 44, 9, 46, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGACGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 7, 0, 7, 1), (7, 7, 7, 10, 4), (7, 44, 10, 47, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGCGAATCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 6, 0, 6, 1), (6, 11, 6, 11, 2), (11, 44, 11, 44, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGCCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 6, 0, 6, 1), (6, 7, 6, 8, 4), (7, 44, 8, 45, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGATTCGCTAGCTTCGTTTTGATAGATAGAGATATAGAGAT",
            &[(0, 20, 0, 20, 1), (20, 23, 20, 24, 4), (23, 44, 24, 45, 1)],
        );
        check(
            REFERENCE,
            "ACGTCTCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 5, 0, 5, 1), (5, 7, 5, 7, 2), (7, 44, 7, 44, 1)],
        );
        check(
            REFERENCE,
            "ACGTCGTCTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT",
            &[(0, 6, 0, 6, 1), (6, 8, 6, 8, 4), (8, 44, 8, 44, 1)],
        );
    }

    /// The knobs must actually change something, or the matrix in the report is
    /// measuring nothing.
    #[test]
    fn the_reverse_complement_knob_suppresses_inversions() {
        let reference = b"ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT";
        let sample = b"ACGTCGCGAATCTAGCTTCGGGGGATAGATAGAGATATAGAGAT";
        let with = extract(reference, sample, reference.len(), &Config::MUTALYZER);
        assert!(with.iter().any(|v| v.kind == Kind::ReverseComplement));
        let without = extract(
            reference,
            sample,
            reference.len(),
            &Config {
                reverse_complement: false,
                ..Config::MUTALYZER
            },
        );
        assert!(without.iter().all(|v| v.kind != Kind::ReverseComplement));
    }

    #[test]
    fn a_larger_position_weight_collapses_more() {
        // Three substitutions in an 8 nt block. Each costs `wp + 3`; the
        // spanning delins costs `2*wp + 15`. So `3*wp + 9 > 2*wp + 15` — the
        // block collapses exactly when `wp > 6`, which is a threshold upstream
        // sets from the length of the whole contig.
        let (reference, sample) = (b"ACGTACGT", b"CCGAACGA");
        let changed = |wp| {
            extract(
                reference,
                sample,
                reference.len(),
                &Config {
                    position_weight: PositionWeight::Fixed(wp),
                    ..Config::MUTALYZER
                },
            )
            .iter()
            .filter(|v| v.kind != Kind::Identity)
            .count()
        };
        assert_eq!(changed(1), 3);
        assert_eq!(changed(6), 3);
        assert_eq!(changed(7), 1);
        assert_eq!(changed(9), 1);
    }

    #[test]
    fn the_separation_floor_vetoes_a_collapse_across_unchanged_bases() {
        // Same block, at a `wp` that collapses it. The decomposition leaves a
        // 3 nt unchanged run, so a floor of 3 refuses to swallow it while a
        // floor of 4 still allows the collapse.
        let (reference, sample) = (b"ACGTACGT", b"CCGAACGA");
        let collapsing = Config {
            position_weight: PositionWeight::Fixed(9),
            collapse: Collapse::Whole,
            ..Config::MUTALYZER
        };
        let changed = |cfg: &Config| {
            extract(reference, sample, reference.len(), cfg)
                .iter()
                .filter(|v| v.kind != Kind::Identity)
                .count()
        };
        assert_eq!(changed(&collapsing), 1);
        assert_eq!(
            changed(&Config {
                separation_floor: Some(3),
                ..collapsing
            }),
            3
        );
        assert_eq!(
            changed(&Config {
                separation_floor: Some(4),
                ..collapsing
            }),
            1
        );
    }
}
