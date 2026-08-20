//! Divergence estimation and alignment rescoring.
//!
//! Ports of three routines from RepeatMasker's `SearchResult.pm`:
//!
//! * [`kimura_divergence`] — `calcKimuraDivergence`
//! * [`k2p_gap_divergence`] — `calcK2PGapDivergence` (Sato's gap-aware K2P)
//! * [`rescore`] — `rescoreAlignment`, including CpG-modified scoring, the
//!   post-hoc xDrop pass, and Phil Green's complexity adjustment
//!
//! # Orientation
//!
//! Every routine here takes the alignment as a pair of equal-length **gapped**
//! strings and assumes RepeatMasker's convention:
//!
//! ```text
//!   subject = consensus / ancestral state
//!     query = genomic / derived state
//! ```
//!
//! CpG detection walks the *subject*; matrix lookups are `matrix[subj][query]`.
//! Swapping the two silently changes the answer, so the argument order is
//! `(query, subject)` everywhere and is checked by tests.
//!
//! # Undefined results
//!
//! The K2P log operand goes non-positive once substitutions saturate.
//! [`kimura_divergence`] returns `None` there; the Perl returns its `100.00`
//! initialiser, and `dfam-curator`'s port returns `NaN`.  Rather than pick a
//! winner, the value is an `Option<f64>` and
//! [`Divergence::or_repeatmasker_default`] reproduces the Perl exactly.
//!
//! [`k2p_gap_divergence`] behaves differently, and not as its structure
//! suggests: on saturation it substitutes the *literal* `1` for the inner log
//! term rather than skipping the calculation, so it returns a defined —
//! typically negative — value.  A fully saturated alignment yields `-50.0`, not
//! a sentinel.  Its `100.00` initialiser survives only when there are no
//! well-characterised bases at all, and is then multiplied by 100 on the way
//! out, so callers see `10000`.  That is the one case where its value is `None`
//! here.
//!
//! # CpG site counting is not consistent in the Perl
//!
//! [`kimura_divergence`] and [`k2p_gap_divergence`] count CpG sites only when
//! `div_cpg_mod` is on, because the counter lives inside that branch.
//! [`rescore`] counts them unconditionally.  Both behaviours are reproduced as
//! found, so the same alignment can report different `cpg_sites` depending on
//! which routine asked.

use crate::error::{Error, Result};
use crate::matrix::SubstMatrix;
use crate::seq;

// ── Substitution classification ───────────────────────────────────────────────

/// Transition (purine↔purine, pyrimidine↔pyrimidine) or transversion.
///
/// ```text
///        Purines
///      A --i-- G
///      |  \ /  |
///      v   v   v
///      |  / \  |
///      C --i-- T
///     Pyrimidines
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubstClass {
    Transition,
    Transversion,
}

/// How lowercase residues are treated when tabulating statistics.
///
/// `SearchResult.pm` looks substitutions up in uppercase-keyed hashes without
/// calling `uc()`, so a lowercase base matches nothing and the position drops
/// out of *both* the substitution counts and the well-characterised
/// denominator.  In practice its callers uppercase beforehand, which makes
/// [`Masking::Ignore`] the behaviour actually exercised — hence the default.
///
/// [`Masking::Lowercase`] is useful in its own right: soft-masking a region of
/// the query, the subject, or both is a way of saying "do not count these
/// positions as aligned".
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum Masking {
    /// Case is ignored; `acgt` counts exactly as `ACGT`.
    #[default]
    Ignore,
    /// Lowercase marks a masked position.  A column is skipped entirely when
    /// **either** side is lowercase.
    Lowercase,
}

impl Masking {
    /// True when this column should be excluded from every tally.
    #[inline]
    pub fn is_masked(self, query: u8, subj: u8) -> bool {
        self == Masking::Lowercase
            && (query.is_ascii_lowercase() || subj.is_ascii_lowercase())
    }
}

/// Classify a substitution, or `None` for identity, gaps, or IUPAC ambiguity.
///
/// Equivalent to RepeatMasker's `%mutType` lookup, which is keyed on the
/// concatenation `query . subject` and therefore symmetric.
#[inline]
pub fn classify_subst(a: u8, b: u8, masking: Masking) -> Option<SubstClass> {
    if masking.is_masked(a, b) {
        return None;
    }
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();
    if a == b || !seq::is_acgt(a) || !seq::is_acgt(b) {
        return None;
    }
    let a_purine = matches!(a, b'A' | b'G');
    let b_purine = matches!(b, b'A' | b'G');
    Some(if a_purine == b_purine {
        SubstClass::Transition
    } else {
        SubstClass::Transversion
    })
}

/// A pairing is *well characterised* when both sides are unambiguous `ACGT`.
///
/// Equivalent to RepeatMasker's `%wellCharacterizedBases` table.  IUPAC codes
/// and gaps are excluded, so this is the denominator for K2P.
#[inline]
pub fn is_well_characterized(query: u8, subj: u8, masking: Masking) -> bool {
    !masking.is_masked(query, subj) && seq::is_acgt(query) && seq::is_acgt(subj)
}

// ── Divergence results ────────────────────────────────────────────────────────

/// Outcome of a divergence calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct Divergence {
    /// Divergence as a percentage, or `None` where the formula is undefined
    /// (no well-characterised bases, or substitutions have saturated the log).
    pub value: Option<f64>,
    /// Transition count.  Fractional when CpG modification is enabled.
    pub transitions: f64,
    pub transversions: u32,
    /// Positions where both sides are unambiguous `ACGT`.
    pub well_characterized: u32,
    /// CpG dinucleotides observed in the subject.
    pub cpg_sites: u32,
    /// Columns with a gap on exactly one side.  Only populated by
    /// [`k2p_gap_divergence`].
    pub gap_len: u32,
}

impl Divergence {
    /// The value RepeatMasker's Perl would have returned, substituting its
    /// initialiser where the formula is undefined.
    ///
    /// `saturated_default` is `100.0` for `calcKimuraDivergence` and `10000.0`
    /// for `calcK2PGapDivergence` — the latter because the Perl multiplies its
    /// `100.00` initialiser by 100 on the way out.
    pub fn or_repeatmasker_default(&self, saturated_default: f64) -> f64 {
        self.value.unwrap_or(saturated_default)
    }
}

/// Kimura two-parameter formula, in percent.
///
/// `transitions` is `f64` to carry the fractional CpG adjustment.  Returns
/// `None` when the log operand is non-positive (saturation) or there are no
/// well-characterised bases.
pub fn k2p(transitions: f64, transversions: f64, well_characterized: u32) -> Option<f64> {
    if well_characterized < 1 {
        return None;
    }
    let n = well_characterized as f64;
    let p = transitions / n;
    let q = transversions / n;

    // Perl: logOperand = (1 - 2p - q) * (1 - 2q)**0.5
    // A negative (1 - 2q) yields NaN under **0.5, which fails the `> 0` test.
    let inner = 1.0 - 2.0 * q;
    if inner < 0.0 {
        return None;
    }
    let log_operand = (1.0 - 2.0 * p - q) * inner.sqrt();
    // Deliberately negated rather than `<= 0.0`: a NaN operand must also fail,
    // which is how the Perl's `if ($logOperand > 0)` behaves.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    if !(log_operand > 0.0) {
        return None;
    }
    Some((-0.5 * log_operand.ln()).abs() * 100.0)
}

// ── Kimura divergence ─────────────────────────────────────────────────────────

/// Port of `SearchResult.pm::calcKimuraDivergence`.
///
/// `query` and `subject` are equal-length gapped strings.  Columns where the
/// subject has a gap are skipped entirely (they are insertions in the query and
/// carry no ancestral state to diverge from).
///
/// # CpG modification
///
/// With `div_cpg_mod`, transitions inside a CpG dinucleotide are discounted:
/// two transitions across the `C` and the `G` count as one, a lone transition
/// counts as one tenth.  Transversions are unaffected.  In human, transitions
/// are ~15× more likely at CpG sites; in rodents they are less likely.
///
/// The one-position lag this needs (a transition at the `C` is only discounted
/// once the following `G` is seen) is implemented exactly as the Perl does, by
/// holding the credit in `prev_trans` and flushing it at the next non-CpG
/// position and at end of sequence.
pub fn kimura_divergence(
    query: &[u8],
    subject: &[u8],
    div_cpg_mod: bool,
    masking: Masking,
) -> Result<Divergence> {
    check_pair(query, subject)?;

    let mut transversions = 0u32;
    let mut transitions = 0f64;
    let mut cpg_sites = 0u32;
    let mut well = 0u32;
    let mut prev_subj = 0u8;
    let mut prev_trans = 0f64;

    for i in 0..subject.len() {
        let s = subject[i];
        let q = query[i];
        if seq::is_gap(s) {
            continue;
        }
        if is_well_characterized(q, s, masking) {
            well += 1;
        }

        let in_cpg = div_cpg_mod
            && prev_subj.eq_ignore_ascii_case(&b'C')
            && s.eq_ignore_ascii_case(&b'G');

        if in_cpg {
            cpg_sites += 1;
            match classify_subst(q, s, masking) {
                Some(SubstClass::Transition) => prev_trans += 1.0,
                Some(SubstClass::Transversion) => transversions += 1,
                None => {}
            }
            prev_trans = discount_cpg(prev_trans);
        } else {
            transitions += prev_trans;
            prev_trans = 0.0;
            match classify_subst(q, s, masking) {
                // Held back: the next position might complete a CpG.
                Some(SubstClass::Transition) => prev_trans = 1.0,
                Some(SubstClass::Transversion) => transversions += 1,
                None => {}
            }
        }
        prev_subj = s;
    }
    transitions += prev_trans;

    Ok(Divergence {
        value: k2p(transitions, transversions as f64, well),
        transitions,
        transversions,
        well_characterized: well,
        cpg_sites,
        gap_len: 0,
    })
}

/// The Perl's CpG discount ladder, reproduced with its exact equality tests.
///
/// `prev_trans == 2` → 1 (both CpG positions mutated, count as one).
/// `prev_trans == 1` → 0.1 (a lone transition, count as a tenth).
///
/// Note the comparisons are exact.  A value of `1.1` — reachable when a
/// discounted `0.1` is incremented inside a run of overlapping CpGs — matches
/// neither arm and passes through untouched.  That is what the Perl does, so it
/// is what this does.
#[inline]
fn discount_cpg(prev_trans: f64) -> f64 {
    if prev_trans == 2.0 {
        1.0
    } else if prev_trans == 1.0 {
        0.1
    } else {
        prev_trans
    }
}

/// Adjusted and unadjusted Kimura divergence from a single pass over the pair.
///
/// [`kimura_divergence`] answers one or the other; display code generally wants
/// both side by side (Linup prints them as adjacent columns), and computing them
/// together keeps the two consistent.
#[derive(Debug, Clone, PartialEq)]
pub struct KimuraStats {
    /// Classic K2P divergence, percent. `None` where the formula is undefined.
    pub kimura: Option<f64>,
    /// K2P with the CpG discount applied, percent.
    pub kimura_adjusted: Option<f64>,
    /// Whole transitions, no CpG discount.
    pub transitions: u32,
    /// Transitions after the CpG discount; fractional by construction.
    pub transitions_adjusted: f64,
    pub transversions: u32,
    /// Positions where both sides are unambiguous `ACGT`.
    pub well_characterized: u32,
    /// CpG dinucleotides observed in the subject.
    pub cpg_sites: u32,
}

/// Both divergences for one pair.
///
/// `subject` supplies the CpG context, matching `calcKimuraDivergence`, where
/// the look-back runs over the subject sequence.
pub fn kimura_stats(query: &[u8], subject: &[u8], masking: Masking) -> Result<KimuraStats> {
    let plain = kimura_divergence(query, subject, false, masking)?;
    let adj = kimura_divergence(query, subject, true, masking)?;
    Ok(KimuraStats {
        kimura: plain.value,
        kimura_adjusted: adj.value,
        transitions: plain.transitions as u32,
        transitions_adjusted: adj.transitions,
        transversions: plain.transversions,
        well_characterized: plain.well_characterized,
        cpg_sites: adj.cpg_sites,
    })
}

/// Mean divergence of `seqs` against `consensus`, percent.
///
/// Rows that yield no defined divergence are skipped rather than counted as
/// zero; returns 0.0 when none are usable.
pub fn mean_kimura(
    consensus: &[u8],
    seqs: &[&[u8]],
    adjusted: bool,
    masking: Masking,
) -> f64 {
    let mut total = 0.0;
    let mut n = 0u32;
    for s in seqs {
        if let Ok(d) = kimura_divergence(s, consensus, adjusted, masking) {
            if let Some(v) = d.value {
                total += v;
                n += 1;
            }
        }
    }
    if n == 0 { 0.0 } else { total / n as f64 }
}

// ── K2P-Gap divergence ────────────────────────────────────────────────────────

/// Port of `SearchResult.pm::calcK2PGapDivergence` — Sato's gap-aware K2P.
///
/// ```text
///   K = 3/4 * w * ln(w) - w/2 * ln[(S - P) * sqrt(S + P - Q)]
///
///   w = (2*alignLen - gapLen) / (2*alignLen)
///   S = (wellCharacterized - (transitions + transversions)) / alignLen
///   P = transitions   / alignLen
///   Q = transversions / alignLen
/// ```
///
/// Validated by RepeatMasker against results supplied by Sato.  Unlike
/// [`kimura_divergence`] this walks *all* columns, counting single-sided gaps
/// toward `gapLen` and skipping only double gaps.
///
/// The returned `value` is already multiplied by 100.
pub fn k2p_gap_divergence(
    query: &[u8],
    subject: &[u8],
    div_cpg_mod: bool,
    masking: Masking,
) -> Result<Divergence> {
    check_pair(query, subject)?;

    let align_len = subject.len() as f64;
    let site_count = align_len * 2.0;

    let mut transversions = 0u32;
    let mut transitions = 0f64;
    let mut cpg_sites = 0u32;
    let mut well = 0u32;
    let mut gap_len = 0u32;
    let mut prev_subj = 0u8;
    let mut prev_trans = 0f64;

    for i in 0..subject.len() {
        let s = subject[i];
        let q = query[i];
        let (sg, qg) = (seq::is_gap(s), seq::is_gap(q));
        if sg && qg {
            continue;
        }
        if sg || qg {
            gap_len += 1;
            continue;
        }
        if is_well_characterized(q, s, masking) {
            well += 1;
        }

        let in_cpg = div_cpg_mod
            && prev_subj.eq_ignore_ascii_case(&b'C')
            && s.eq_ignore_ascii_case(&b'G');

        if in_cpg {
            cpg_sites += 1;
            match classify_subst(q, s, masking) {
                Some(SubstClass::Transition) => prev_trans += 1.0,
                Some(SubstClass::Transversion) => transversions += 1,
                None => {}
            }
            prev_trans = discount_cpg(prev_trans);
        } else {
            transitions += prev_trans;
            prev_trans = 0.0;
            match classify_subst(q, s, masking) {
                Some(SubstClass::Transition) => prev_trans = 1.0,
                Some(SubstClass::Transversion) => transversions += 1,
                None => {}
            }
        }
        prev_subj = s;
    }
    transitions += prev_trans;

    let value = if well >= 1 && align_len > 0.0 {
        let p = transitions / align_len;
        let q_ = transversions as f64 / align_len;
        let s_ = (well as f64 - (transitions + transversions as f64)) / align_len;
        let w = (site_count - gap_len as f64) / site_count;

        // The Perl substitutes 1 (so ln contributes 0 after the outer
        // multiply... in fact it uses the literal 1, not ln(1)) when either
        // operand is out of domain.  Reproduced exactly.
        let log_term1 = if w > 0.0 { w.ln() } else { 1.0 };
        let inner = s_ + p - q_;
        let log_op = if inner >= 0.0 { (s_ - p) * inner.sqrt() } else { f64::NAN };
        let log_term2 = if log_op > 0.0 { log_op.ln() } else { 1.0 };

        let k = (0.75 * w * log_term1) - ((w / 2.0) * log_term2);
        Some(k * 100.0)
    } else {
        None
    };

    Ok(Divergence {
        value,
        transitions,
        transversions,
        well_characterized: well,
        cpg_sites,
        gap_len,
    })
}

// ── Rescoring ─────────────────────────────────────────────────────────────────

/// Parameters for [`rescore`] — the analogue of `rescoreAlignment`'s named args.
#[derive(Debug, Clone)]
pub struct RescoreParams<'a> {
    /// Scores are looked up as `matrix[subject][query]`.
    pub matrix: &'a SubstMatrix,

    /// Penalty for opening a gap, **including its first position**.  Negative.
    pub gap_open: i32,
    /// Penalty per insertion position beyond the first (a gap in the subject).
    pub ins_gap_extend: i32,
    /// Penalty per deletion position beyond the first (a gap in the query).
    pub del_gap_extend: i32,

    /// Score only transversions at CpG sites, substituting the matrix's `C/C`
    /// and `G/G` identity scores for transitions.
    pub score_cpg_mod: bool,
    /// Apply the CpG discount to the divergence calculation.
    pub div_cpg_mod: bool,
    /// Apply Phil Green's complexity-adjusted scoring.  Requires the matrix to
    /// carry `FREQS` (and hence a lambda).
    pub complexity_adjust: bool,
    /// Run the post-hoc xDrop pass, recording where HSPs would have been called.
    pub xdrop: Option<i32>,
    /// How lowercase residues are treated.  See [`Masking`].
    pub masking: Masking,
}

impl<'a> RescoreParams<'a> {
    /// Defaults matching a typical crossmatch invocation: `-gap_init -25`,
    /// `-gap_ext -5`, no CpG modification, no complexity adjustment.
    pub fn new(matrix: &'a SubstMatrix) -> Self {
        RescoreParams {
            matrix,
            gap_open: -25,
            ins_gap_extend: -5,
            del_gap_extend: -5,
            score_cpg_mod: false,
            div_cpg_mod: false,
            complexity_adjust: false,
            xdrop: None,
            masking: Masking::default(),
        }
    }

    /// Use a single extension penalty for both insertions and deletions —
    /// the Perl's `gapExtPenalty`, which overrides the separate values.
    pub fn with_gap_extend(mut self, ext: i32) -> Self {
        self.ins_gap_extend = ext;
        self.del_gap_extend = ext;
        self
    }
}

/// Everything `rescoreAlignment` returns, plus the gap tallies it computes
/// internally but discards.
#[derive(Debug, Clone)]
pub struct RescoreResult {
    /// Raw score, or the complexity-adjusted score when that was requested.
    pub score: i32,
    /// Raw score before complexity adjustment.
    pub raw_score: i32,
    /// Sum of matrix scores over aligned columns only — no gap penalties.
    pub ungapped_raw_score: i32,

    /// Kimura divergence and its components (CpG-modified iff `div_cpg_mod`).
    pub divergence: Divergence,

    /// Percent insertion: inserted bases as a fraction of subject bases.
    /// Clamped to 100.
    pub pct_insert: f64,
    /// Percent deletion: deleted bases as a fraction of query bases.
    /// Clamped to 100.
    pub pct_delete: f64,

    /// Cumulative score after each alignment column.
    pub position_scores: Vec<i32>,
    /// `(start, end)` column pairs where xDrop would have broken the alignment
    /// into HSPs.  Empty unless [`RescoreParams::xdrop`] was set.
    pub xdrop_fragments: Vec<(usize, usize)>,

    pub insertion_inits: u32,
    pub insertion_extns: u32,
    pub deletion_inits: u32,
    pub deletion_extns: u32,
}

impl RescoreResult {
    /// Total inserted bases (gap columns in the subject).
    pub fn insertions(&self) -> u32 {
        self.insertion_inits + self.insertion_extns
    }
    /// Total deleted bases (gap columns in the query).
    pub fn deletions(&self) -> u32 {
        self.deletion_inits + self.deletion_extns
    }
}

/// Port of `SearchResult.pm::rescoreAlignment`.
///
/// Rescores an existing alignment under a (possibly different) scoring system
/// without altering it.  `query` and `subject` are equal-length gapped strings.
///
/// Under the subject-is-consensus convention a gap in the **subject** is an
/// insertion in the genomic sequence, and a gap in the **query** is a deletion.
// The index loops over `position_scores` retroactively patch a range of already
// emitted cumulative scores.  They are transcribed from the Perl line for line
// — including the `<=` / `<` asymmetry between the two CpG branches — and an
// iterator rewrite would obscure exactly the detail that matters.
#[allow(clippy::needless_range_loop)]
pub fn rescore(query: &[u8], subject: &[u8], params: &RescoreParams<'_>) -> Result<RescoreResult> {
    check_pair(query, subject)?;
    let matrix = params.matrix;

    let c_idx = matrix
        .index_of(b'C')
        .ok_or_else(|| Error::Scoring("matrix alphabet lacks 'C'".into()))?;
    let g_idx = matrix
        .index_of(b'G')
        .ok_or_else(|| Error::Scoring("matrix alphabet lacks 'G'".into()))?;
    let c_score = matrix.score_idx(c_idx, c_idx);
    let g_score = matrix.score_idx(g_idx, g_idx);

    let mut score = 0i32;
    let mut ungapped_raw_score = 0i32;
    let mut position_scores: Vec<i32> = Vec::with_capacity(subject.len());
    let mut mat_counts = vec![0u64; matrix.size()];

    let (mut ins_inits, mut ins_extns) = (0u32, 0u32);
    let (mut del_inits, mut del_extns) = (0u32, 0u32);

    let mut transversions = 0u32;
    let mut transitions = 0f64;
    let mut cpg_sites = 0u32;
    let mut well = 0u32;

    let mut prev_subj = 0u8;
    // Matrix score charged at the previous subject position, and where its
    // cumulative entry landed — both needed to retroactively correct a
    // transition at the `C` of a CpG once the `G` is seen.
    let mut prev_score = -1i32;
    let mut prev_pos: isize = -1;
    let mut prev_trans = 0f64;

    for i in 0..subject.len() {
        let s = subject[i];
        let q = query[i];

        // ── Insertion: gap in the subject (consensus) ──────────────────────
        if seq::is_gap(s) {
            if i > 0 && seq::is_gap(subject[i - 1]) {
                score += params.ins_gap_extend;
                ins_extns += 1;
            } else {
                score += params.gap_open;
                ins_inits += 1;
            }
            position_scores.push(score);
            continue;
        }

        // ── Deletion: gap in the query (genomic) ───────────────────────────
        if seq::is_gap(q) {
            if i > 0 && seq::is_gap(query[i - 1]) {
                score += params.del_gap_extend;
                del_extns += 1;
            } else {
                score += params.gap_open;
                del_inits += 1;
            }
            position_scores.push(score);

            if prev_subj.eq_ignore_ascii_case(&b'C') && s.eq_ignore_ascii_case(&b'G') {
                cpg_sites += 1;
                if params.score_cpg_mod && prev_score < c_score && prev_trans != 0.0 {
                    let diff = c_score - prev_score;
                    // Inclusive of the just-pushed entry — unlike the aligned
                    // branch below, nothing here overwrites it afterwards.
                    for j in prev_pos.max(0) as usize..position_scores.len() {
                        position_scores[j] += diff;
                    }
                    score += diff;
                    ungapped_raw_score += diff;
                }
                if params.div_cpg_mod && prev_trans == 1.0 {
                    // The CpG straddles the gap:  query T-
                    //                           subject CG
                    // Apply the one-tenth rule and clear the held credit.
                    transitions += 0.1;
                    prev_trans = 0.0;
                }
            }
            prev_subj = s;
            // Deliberately the C/C identity score, not a real matrix lookup:
            // a `C/-` column must not be retroactively rescored as `C/C`.
            prev_score = c_score;
            prev_pos = position_scores.len() as isize - 1;
            continue;
        }

        // ── Aligned pair ───────────────────────────────────────────────────
        let mat_score = matrix.score(s, q).ok_or_else(|| {
            Error::Scoring(format!(
                "column {i}: pair {}/{} is outside the matrix alphabet {:?}",
                q as char,
                s as char,
                String::from_utf8_lossy(matrix.alphabet())
            ))
        })?;
        score += mat_score;
        position_scores.push(score);
        ungapped_raw_score += mat_score;

        if let Some(qi) = matrix.index_of(q) {
            mat_counts[qi] += 1;
        }
        if is_well_characterized(q, s, params.masking) {
            well += 1;
        }

        let mut in_cpg = false;
        if prev_subj.eq_ignore_ascii_case(&b'C') && s.eq_ignore_ascii_case(&b'G') {
            in_cpg = true;
            cpg_sites += 1;

            if params.score_cpg_mod {
                if prev_score < c_score && prev_trans != 0.0 {
                    let diff = c_score - prev_score;
                    // Exclusive of the last entry: it is rewritten below.
                    let last = position_scores.len() - 1;
                    for j in prev_pos.max(0) as usize..last {
                        position_scores[j] += diff;
                    }
                    score += diff;
                    ungapped_raw_score += diff;
                }
                let last = position_scores.len() - 1;
                let prev_cum = if last > 0 { position_scores[last - 1] } else { 0 };
                if classify_subst(q, s, params.masking) == Some(SubstClass::Transition) {
                    // Charge the G/G identity score instead of the transition.
                    position_scores[last] = prev_cum + g_score;
                    score += g_score - mat_score;
                    ungapped_raw_score += g_score - mat_score;
                } else {
                    position_scores[last] = prev_cum + mat_score;
                }
            }
        }

        if params.div_cpg_mod && in_cpg {
            match classify_subst(q, s, params.masking) {
                Some(SubstClass::Transition) => prev_trans += 1.0,
                Some(SubstClass::Transversion) => transversions += 1,
                None => {}
            }
            prev_trans = discount_cpg(prev_trans);
        } else {
            transitions += prev_trans;
            prev_trans = 0.0;
            match classify_subst(q, s, params.masking) {
                Some(SubstClass::Transition) => prev_trans = 1.0,
                Some(SubstClass::Transversion) => transversions += 1,
                None => {}
            }
        }

        prev_subj = s;
        prev_score = mat_score;
        prev_pos = position_scores.len() as isize - 1;
    }
    transitions += prev_trans;

    let xdrop_fragments = match params.xdrop {
        Some(x) => xdrop_pass(&position_scores, x),
        None => Vec::new(),
    };

    // Percent insert/delete are relative to the ungapped span of the *other*
    // sequence, matching the Perl's use of subjEnd-subjStart+1 / qryEnd-qryStart+1.
    let subj_bases = seq::ungapped_len(subject);
    let query_bases = seq::ungapped_len(query);
    let pct_insert = pct(ins_inits + ins_extns, subj_bases);
    let pct_delete = pct(del_inits + del_extns, query_bases);

    let divergence = Divergence {
        value: k2p(transitions, transversions as f64, well),
        transitions,
        transversions,
        well_characterized: well,
        cpg_sites,
        gap_len: 0,
    };

    let raw_score = score;
    let final_score = if params.complexity_adjust {
        complexity_adjust(score, &mat_counts, matrix)?
    } else {
        score
    };

    Ok(RescoreResult {
        score: final_score,
        raw_score,
        ungapped_raw_score,
        divergence,
        pct_insert,
        pct_delete,
        position_scores,
        xdrop_fragments,
        insertion_inits: ins_inits,
        insertion_extns: ins_extns,
        deletion_inits: del_inits,
        deletion_extns: del_extns,
    })
}

/// Phil Green's complexity-adjusted scoring, as cross_match computes it.
///
/// Discounts the raw score by the information content of the query bases
/// actually aligned, so a hit made of low-complexity sequence scores less than
/// the same raw score made of compositionally surprising sequence.
///
/// `mat_counts` is indexed parallel to the matrix alphabet and counts **query**
/// bases.  Symbols with zero background frequency, and symbols whose frequency
/// is exactly 1.0, are excluded — matching the Perl's `log(freq) != 0` guard.
pub fn complexity_adjust(
    raw_score: i32,
    mat_counts: &[u64],
    matrix: &SubstMatrix,
) -> Result<i32> {
    let lambda = matrix.lambda().ok_or_else(|| {
        Error::Scoring(format!(
            "complexity adjustment needs a lambda; matrix {:?} has no FREQS line",
            matrix.name()
        ))
    })?;
    if lambda <= 0.0 {
        return Err(Error::Scoring(format!("non-positive lambda {lambda}")));
    }
    let freqs = matrix.freqs();

    let mut t_factor = 0f64;
    let mut t_sum = 0f64;
    let mut t_counts = 0f64;

    for (i, &count) in mat_counts.iter().enumerate() {
        if count == 0 {
            continue;
        }
        let f = freqs.get(i).copied().unwrap_or(0.0);
        if f <= 0.0 || f.ln() == 0.0 {
            continue;
        }
        let c = count as f64;
        t_factor += c * c.ln();
        t_sum += c * f.ln();
        t_counts += c;
    }

    if t_counts != 0.0 {
        t_factor -= t_counts * t_counts.ln();
    }
    t_sum -= t_factor;

    // cross_match's rounding, not sprintf's: truncate `x + 0.999`.
    let adj = (raw_score as f64 + t_sum / lambda + 0.999).trunc();
    if !adj.is_finite() || adj < 0.0 {
        return Ok(0);
    }
    Ok(adj as i32)
}

/// Post-hoc xDrop pass over cumulative position scores.
///
/// Reports `(start, end)` column pairs where a real xDrop-limited aligner would
/// have emitted separate HSPs.  Trailing fragments shorter than six columns are
/// discarded, as in the Perl.
// `i` is both the loop index and a recorded fragment boundary, so the indexed
// form is the honest transcription of the Perl.
#[allow(clippy::needless_range_loop)]
pub fn xdrop_pass(position_scores: &[i32], xdrop: i32) -> Vec<(usize, usize)> {
    let mut fragments = Vec::new();
    let mut last_highest = 0i32;
    let mut last_highest_pos = 0usize;
    let mut start_point = 0usize;
    let mut subtract = 0i32;

    for i in 0..position_scores.len() {
        let mut adj = position_scores[i] - subtract;
        if adj < 0 {
            adj = 0;
            start_point = i;
            subtract = position_scores[i];
            last_highest = 0;
            last_highest_pos = i + 1;
        }
        if adj >= last_highest {
            last_highest = position_scores[i] - subtract;
            last_highest_pos = i;
            continue;
        }
        if last_highest - adj > xdrop {
            fragments.push((start_point, last_highest_pos));
            subtract = position_scores[i];
            start_point = i + 1;
            last_highest = 0;
            last_highest_pos = i + 1;
        }
    }

    if last_highest_pos > start_point && last_highest_pos - start_point > 5 {
        fragments.push((start_point, last_highest_pos));
    }
    fragments
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn check_pair(query: &[u8], subject: &[u8]) -> Result<()> {
    if query.len() != subject.len() {
        return Err(Error::Scoring(format!(
            "gapped strings differ in length: query {} vs subject {}",
            query.len(),
            subject.len()
        )));
    }
    Ok(())
}

fn pct(n: u32, denom: usize) -> f64 {
    if denom < 1 {
        return 100.0;
    }
    ((n as f64 * 100.0) / denom as f64).min(100.0)
}

#[cfg(test)]
mod tests {
    /// Lowercase masking must remove a column from *both* tallies, or the
    /// numerator and denominator disagree and divergence is over-reported.
    /// That inconsistency is exactly what `dfam-curator`'s kimura.rs has.
    #[test]
    fn lowercase_masking_excludes_a_column_from_every_tally() {
        // One transversion at position 1, everything else identical.
        let q = b"ACGTACGT";
        let s = b"AAGTACGT";

        let plain = kimura_divergence(q, s, false, Masking::Ignore).unwrap();
        assert_eq!(plain.well_characterized, 8);
        assert_eq!(plain.transversions, 1);

        // Soft-mask the substituted column on ONE side only.  Under
        // Masking::Lowercase the whole column drops out: the transversion is
        // not counted, and neither is the position.
        let masked_q = b"AcGTACGT";
        let m = kimura_divergence(masked_q, s, false, Masking::Lowercase).unwrap();
        assert_eq!(m.well_characterized, 7, "masked column must leave the denominator");
        assert_eq!(m.transversions, 0, "masked column must leave the numerator");

        // Ignore is case-blind: same answer as the uppercase input.
        let ignored = kimura_divergence(masked_q, s, false, Masking::Ignore).unwrap();
        assert_eq!(ignored.well_characterized, 8);
        assert_eq!(ignored.transversions, 1);
    }

    /// Masking applies when *either* side is lowercase, matching the Perl's
    /// failed hash lookup on `$q . $s`.
    #[test]
    fn masking_triggers_from_either_side() {
        assert!(Masking::Lowercase.is_masked(b'a', b'A'));
        assert!(Masking::Lowercase.is_masked(b'A', b'a'));
        assert!(!Masking::Lowercase.is_masked(b'A', b'A'));
        assert!(!Masking::Ignore.is_masked(b'a', b'a'));
    }

    use super::*;

    const M14P35G: &str = "\
FREQS A 0.325 C 0.175 G 0.175 T 0.325
  A   R   G   C   Y   T   K   M   S   W   N   X
  8   0 -10 -18 -19 -21 -15  -4 -14  -6  -1 -30
  3   3  12 -17 -18 -19  -9  -8  -8  -9  -1 -30
 -7   2  12 -16 -16 -17  -2 -11  -1 -12  -1 -30
-17 -16 -16  12   2  -7 -11  -2  -1 -12  -1 -30
-19 -18 -17  12   0   3  -8  -9  -8  -9  -1 -30
-21 -19 -18 -10   0   8  -4 -15 -14  -6  -1 -30
-14  -8  -2 -13  -8  -4  -3 -13  -8  -9  -1 -30
 -4  -8 -13  -2  -8 -14 -13  -3  -8  -9  -1 -30
-12  -7  -1  -1  -7 -12  -7  -7  -1 -12  -1 -30
 -6 -10 -14 -14 -10  -6 -10 -10 -14  -6  -1 -30
 -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
-30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
";

    fn matrix() -> SubstMatrix {
        SubstMatrix::parse(M14P35G).unwrap()
    }

    // ── Classification ────────────────────────────────────────────────────

    #[test]
    fn classification_matches_the_mutation_table() {
        use SubstClass::*;
        assert_eq!(classify_subst(b'C', b'T', Masking::Ignore), Some(Transition));
        assert_eq!(classify_subst(b'A', b'G', Masking::Ignore), Some(Transition));
        assert_eq!(classify_subst(b'G', b'T', Masking::Ignore), Some(Transversion));
        assert_eq!(classify_subst(b'A', b'C', Masking::Ignore), Some(Transversion));
        assert_eq!(classify_subst(b'A', b'A', Masking::Ignore), None);
        // IUPAC codes and gaps are unclassifiable.
        assert_eq!(classify_subst(b'A', b'R', Masking::Ignore), None);
        assert_eq!(classify_subst(b'-', b'A', Masking::Ignore), None);
    }

    #[test]
    fn well_characterized_requires_acgt_on_both_sides() {
        assert!(is_well_characterized(b'A', b'C', Masking::Ignore));
        assert!(is_well_characterized(b'a', b'c', Masking::Ignore));
        assert!(!is_well_characterized(b'N', b'C', Masking::Ignore));
        assert!(!is_well_characterized(b'A', b'-', Masking::Ignore));
    }

    // ── Kimura ────────────────────────────────────────────────────────────

    #[test]
    fn identical_sequences_have_zero_divergence() {
        let d = kimura_divergence(b"ACGT", b"ACGT", false, Masking::Ignore).unwrap();
        assert_eq!(d.transitions, 0.0);
        assert_eq!(d.transversions, 0);
        assert_eq!(d.well_characterized, 4);
        assert!(d.value.unwrap().abs() < 1e-12);
    }

    #[test]
    fn saturation_is_undefined_not_a_number() {
        // Every position a transition: p = 1, so 1 - 2p - q < 0.
        let d = kimura_divergence(b"GGGG", b"AAAA", false, Masking::Ignore).unwrap();
        assert_eq!(d.transitions, 4.0);
        assert_eq!(d.value, None);
        // RepeatMasker's Perl would have printed its initialiser here.
        assert_eq!(d.or_repeatmasker_default(100.0), 100.0);
    }

    #[test]
    fn subject_gaps_are_skipped_entirely() {
        // The subject gap column contributes nothing — not to the denominator,
        // not to the substitution counts, and it does not update the CpG state.
        let with_gap = kimura_divergence(b"ACGT", b"AC-T", false, Masking::Ignore).unwrap();
        assert_eq!(with_gap.well_characterized, 3);
        assert_eq!(with_gap.transitions, 0.0);
        assert_eq!(with_gap.transversions, 0);
    }

    #[test]
    fn cpg_both_positions_mutated_counts_as_one_transition() {
        // subject CG, query TA: C->T and G->A are both transitions.
        let d = kimura_divergence(b"ATAX", b"ACGX", true, Masking::Ignore).unwrap();
        assert_eq!(d.cpg_sites, 1);
        assert!((d.transitions - 1.0).abs() < 1e-12, "{}", d.transitions);
    }

    #[test]
    fn cpg_single_mutation_counts_as_a_tenth() {
        let d = kimura_divergence(b"ATGX", b"ACGX", true, Masking::Ignore).unwrap();
        assert_eq!(d.cpg_sites, 1);
        assert!((d.transitions - 0.1).abs() < 1e-12, "{}", d.transitions);
    }

    #[test]
    fn cpg_discount_is_off_unless_requested() {
        let off = kimura_divergence(b"ATAX", b"ACGX", false, Masking::Ignore).unwrap();
        assert_eq!(off.transitions, 2.0);
        assert_eq!(off.cpg_sites, 0);
    }

    #[test]
    fn non_cpg_transitions_are_not_discounted() {
        // Subject AG is not a CpG, so the A->G transition counts in full.
        let d = kimura_divergence(b"GG", b"AG", true, Masking::Ignore).unwrap();
        assert_eq!(d.transitions, 1.0);
        assert_eq!(d.cpg_sites, 0);
    }

    #[test]
    fn argument_order_is_query_then_subject() {
        // A subject CpG with a query transition is discounted; swapping the
        // arguments removes the CpG context and the discount with it.
        let right = kimura_divergence(b"ATGX", b"ACGX", true, Masking::Ignore).unwrap();
        let swapped = kimura_divergence(b"ACGX", b"ATGX", true, Masking::Ignore).unwrap();
        assert_eq!(right.cpg_sites, 1);
        assert_eq!(swapped.cpg_sites, 0);
        assert!(right.transitions < swapped.transitions);
    }

    // ── K2P-Gap ───────────────────────────────────────────────────────────

    #[test]
    fn k2p_gap_counts_single_sided_gaps_and_skips_double_gaps() {
        let d = k2p_gap_divergence(b"AC-GT", b"ACG-T", false, Masking::Ignore).unwrap();
        assert_eq!(d.gap_len, 2);
        assert_eq!(d.well_characterized, 3);

        let both = k2p_gap_divergence(b"AC--T", b"AC--T", false, Masking::Ignore).unwrap();
        assert_eq!(both.gap_len, 0);
        assert_eq!(both.well_characterized, 3);
    }

    #[test]
    fn k2p_gap_is_defined_for_a_clean_alignment() {
        let d = k2p_gap_divergence(b"ACGTACGTAC", b"ACGTACGTAC", false, Masking::Ignore).unwrap();
        assert!(d.value.is_some());
    }

    // ── Rescoring ─────────────────────────────────────────────────────────

    #[test]
    fn rescore_sums_matrix_scores_for_an_ungapped_alignment() {
        let m = matrix();
        let p = RescoreParams::new(&m);
        // Perfect match: A(8) C(12) G(12) T(8) = 40.
        let r = rescore(b"ACGT", b"ACGT", &p).unwrap();
        assert_eq!(r.score, 40);
        assert_eq!(r.ungapped_raw_score, 40);
        assert_eq!(r.position_scores, vec![8, 20, 32, 40]);
        assert_eq!(r.pct_insert, 0.0);
        assert_eq!(r.pct_delete, 0.0);
    }

    #[test]
    fn rescore_charges_gap_open_then_extension() {
        let m = matrix();
        let p = RescoreParams::new(&m); // open -25, extend -5
        // Query gap of length 3 => a deletion: -25, -5, -5.
        let r = rescore(b"A---T", b"ACGAT", &p).unwrap();
        assert_eq!(r.deletion_inits, 1);
        assert_eq!(r.deletion_extns, 2);
        assert_eq!(r.insertions(), 0);
        // A(8) + (-25) + (-5) + (-5) + T(8) = -19
        assert_eq!(r.score, -19);
        // Gap penalties never enter the ungapped score.
        assert_eq!(r.ungapped_raw_score, 16);
    }

    #[test]
    fn insertions_and_deletions_are_named_from_the_consensus_view() {
        let m = matrix();
        let p = RescoreParams::new(&m);
        // Gap in the subject (consensus) = insertion in the genomic query.
        let ins = rescore(b"ACGT", b"AC-T", &p).unwrap();
        assert_eq!(ins.insertions(), 1);
        assert_eq!(ins.deletions(), 0);
        // Gap in the query = deletion.
        let del = rescore(b"AC-T", b"ACGT", &p).unwrap();
        assert_eq!(del.insertions(), 0);
        assert_eq!(del.deletions(), 1);
    }

    #[test]
    fn pct_insert_and_delete_use_the_opposite_sequence_as_denominator() {
        let m = matrix();
        let p = RescoreParams::new(&m);
        // One inserted base; subject has 3 ungapped bases => 33.3%.
        let r = rescore(b"ACGT", b"AC-T", &p).unwrap();
        assert!((r.pct_insert - 100.0 / 3.0).abs() < 1e-9);
        assert_eq!(r.pct_delete, 0.0);
    }

    #[test]
    fn score_cpg_mod_replaces_a_cpg_transition_with_the_identity_score() {
        let m = matrix();
        let mut p = RescoreParams::new(&m);
        // Subject CG, query CA: the G->A column is a transition inside a CpG.
        let plain = rescore(b"CA", b"CG", &p).unwrap();
        p.score_cpg_mod = true;
        let modded = rescore(b"CA", b"CG", &p).unwrap();
        // G/A scores -7; with the modification it is charged G/G = 12 instead.
        assert_eq!(plain.score, 12 + -7);
        assert_eq!(modded.score, 12 + 12);
        assert!(modded.score > plain.score);
    }

    #[test]
    fn score_cpg_mod_leaves_transversions_alone() {
        let m = matrix();
        let mut p = RescoreParams::new(&m);
        // Subject CG, query CC: G->C is a transversion, so no substitution.
        let plain = rescore(b"CC", b"CG", &p).unwrap();
        p.score_cpg_mod = true;
        let modded = rescore(b"CC", b"CG", &p).unwrap();
        assert_eq!(plain.score, modded.score);
    }

    #[test]
    fn position_scores_stay_cumulative_after_cpg_correction() {
        let m = matrix();
        let p = RescoreParams { score_cpg_mod: true, ..RescoreParams::new(&m) };
        let r = rescore(b"ATA", b"ACG", &p).unwrap();
        // Whatever the corrections, the last cumulative entry must equal the score.
        assert_eq!(*r.position_scores.last().unwrap(), r.raw_score);
    }

    #[test]
    fn complexity_adjustment_penalises_low_complexity_hits() {
        let m = matrix();
        let plain = RescoreParams::new(&m);
        let adjusted = RescoreParams { complexity_adjust: true, ..RescoreParams::new(&m) };

        // A homopolymer run carries little information.
        let q = b"AAAAAAAAAAAAAAAAAAAA";
        let raw = rescore(q, q, &plain).unwrap();
        let adj = rescore(q, q, &adjusted).unwrap();
        assert!(adj.score < raw.score, "{} !< {}", adj.score, raw.score);
        assert_eq!(adj.raw_score, raw.score);
    }

    #[test]
    fn complexity_adjustment_barely_touches_a_compositionally_typical_hit() {
        let m = matrix();
        let adjusted = RescoreParams { complexity_adjust: true, ..RescoreParams::new(&m) };
        // Base composition close to the matrix's own background frequencies.
        let q = b"ATATATATACGCGATATATAT";
        let raw = rescore(q, q, &RescoreParams::new(&m)).unwrap();
        let adj = rescore(q, q, &adjusted).unwrap();
        let loss = raw.score - adj.score;
        assert!(loss >= 0);
        assert!(loss < raw.score / 2, "unexpected large penalty: {loss}");
    }

    #[test]
    fn complexity_adjustment_without_freqs_is_an_error() {
        let m = SubstMatrix::parse("  A   C   G   T\n1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1\n").unwrap();
        let p = RescoreParams { complexity_adjust: true, ..RescoreParams::new(&m) };
        assert!(rescore(b"ACGT", b"ACGT", &p).is_err());
    }

    #[test]
    fn rescore_rejects_symbols_outside_the_matrix_alphabet() {
        let m = matrix();
        let p = RescoreParams::new(&m);
        // '@' is in neither the alphabet nor the gap set.
        assert!(rescore(b"A@GT", b"ACGT", &p).is_err());
    }

    #[test]
    fn rescore_rejects_mismatched_lengths() {
        let m = matrix();
        let p = RescoreParams::new(&m);
        assert!(rescore(b"ACG", b"ACGT", &p).is_err());
    }

    // ── xDrop ─────────────────────────────────────────────────────────────

    #[test]
    fn xdrop_keeps_a_single_uninterrupted_run() {
        let scores: Vec<i32> = (1..=20).collect();
        let frags = xdrop_pass(&scores, 10);
        assert_eq!(frags, vec![(0, 19)]);
    }

    #[test]
    fn xdrop_splits_on_a_deep_valley() {
        // Climb to 100, sag by 50 (well past the drop threshold) while staying
        // positive, then climb again.
        let mut scores: Vec<i32> = (1..=20).map(|i| i * 5).collect(); // 5..100
        scores.extend((1..=10).map(|i| 100 - i * 5)); // 95..50
        scores.extend((1..=20).map(|i| 50 + i * 5)); // 55..150

        let frags = xdrop_pass(&scores, 25);
        assert_eq!(frags.len(), 2, "{frags:?}");
        // The first fragment ends at the peak, not partway down the slope.
        assert_eq!(frags[0], (0, 19));
        // The second picks up after the sag and runs to the end.
        assert_eq!(frags[1].1, scores.len() - 1);
    }

    #[test]
    fn xdrop_resets_instead_of_splitting_once_the_score_goes_negative() {
        // A characteristic of the Perl worth pinning: the `adjScore < 0` reset
        // is evaluated *before* the drop test, so an alignment that decays
        // through zero re-bases rather than emitting a fragment at its peak.
        let mut scores: Vec<i32> = (1..=20).collect(); // peak 20
        scores.extend((0..30).map(|i| 20 - i * 5)); // straight down through 0
        let frags = xdrop_pass(&scores, 25);
        assert!(
            !frags.iter().any(|&(s, _)| s == 0),
            "peak should have been re-based away, got {frags:?}"
        );
    }

    #[test]
    fn xdrop_discards_short_trailing_fragments() {
        // Only five columns of signal at the end — below the Perl's >5 cutoff.
        let scores = vec![1, 2, 3, 4, 5];
        assert!(xdrop_pass(&scores, 10).is_empty());
    }

    #[test]
    fn xdrop_is_skipped_unless_requested() {
        let m = matrix();
        let r = rescore(b"ACGT", b"ACGT", &RescoreParams::new(&m)).unwrap();
        assert!(r.xdrop_fragments.is_empty());
    }
}
