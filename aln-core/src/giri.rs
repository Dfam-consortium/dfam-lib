//! The original GIRI consensus caller.
//!
//! A port of `MultipleAlignment::getConsensus()` from `giri_cpp_lib` (libbio),
//! by way of the verified Rust in `dfam-curator`'s `giri.rs`. This is what
//! `acons`/`autocons` reach with `--orig`; the default is [`crate::consensus`].
//!
//! Module and function names mirror `dfam_curator::giri` so migrating that crate
//! onto `aln-core` stays mechanical.
//!
//! # How it differs from the Dfam caller
//!
//! * **12 symbols** (`A R G C Y T K M S W N X`) against the Dfam caller's 18 —
//!   it cannot emit `B/D/H/V/Z`.
//! * **Gap weights** `-5` base-vs-gap and `+2` gap-vs-gap, against `-6`/`+3`.
//! * **`N` penalty** `-5` against the Dfam matrix's softened `-2`.
//! * **A minimum-coverage gate** ([`get_consensus`]'s `min_non_gap_count`,
//!   `acons --min`) that the Dfam caller has no equivalent for.
//! * **Ties resolve the other way.** Candidates are scanned in *reverse* with a
//!   strict `>`, so gap beats `X` beats `N` beats … beats `A`. The Dfam caller
//!   prefers `N`.
//! * **CpG restoration is a separate pass.** `acons` restores CpG on this path
//!   only when a species is given (`--mam`), using a probabilistic model quite
//!   unlike the Dfam caller's deterministic bonus — see [`restore_cpg`].
//!
//! Both halves of the C++ are here: the consensus caller and the `--mam`
//! statistical CpG restoration ([`restore_cpg`]), ported from `dfam-curator`.
//!
//! [`fixed`] holds a corrected variant of the restoration — GIRI's `firstBase`
//! lookup table is inert and its `CG` mask neither skips gaps nor stops at the
//! motif end. It is kept separate rather than folded in because the faithful
//! version is what reproduces `acons`, and the Dfam CpG benchmarks score both.

/// Interior gap (GIRI `GAPCHAR`).
pub const GAP: u8 = b'-';

/// Flanking padding. GIRI uses `<`/`>`; this crate canonicalises on a space, and
/// [`crate::seq::from_giri_padding`] converts.
pub const PAD: u8 = b' ';

/// `DNACONMATRIX`'s alphabet, in its native order.
const ALPHABET: [u8; 12] = *b"ARGCYTKMSWNX";

/// Base vs gap. GIRI parses `GAP -25 -5` as `gapIni = -25; gapExt = -5`, and
/// only the extension participates in consensus scoring — the `-25` is an
/// *alignment* penalty that never reaches here.
const GAP_EXT_PENALTY: i32 = -5;

/// Gap vs gap: `-gapExt / 2` under C++ integer division, so `+2`.
const GAP_MATCH_SCORE: i32 = 2;

/// GIRI `DNACONMATRIX`; rows are the candidate, columns the observed base.
#[rustfmt::skip]
static DNACONMATRIX: [[i32; 12]; 12] = [
    //        A    R    G    C    Y    T    K    M    S    W    N    X
    /* A */ [  9,   0,  -8, -15, -16, -17, -13,  -3, -11,  -4,  -5,  -7],
    /* R */ [  2,   1,   1, -15, -15, -16,  -7,  -6,  -6,  -7,  -5,  -7],
    /* G */ [ -4,   3,  10, -14, -14, -15,  -2,  -9,  -2,  -9,  -5,  -7],
    /* C */ [-15, -14, -14,  10,   3,  -4,  -9,  -2,  -2,  -9,  -5,  -7],
    /* Y */ [-16, -15, -15,   1,   1,   2,  -6,  -7,  -6,  -7,  -5,  -7],
    /* T */ [-17, -16, -15,  -8,   0,   9,  -3, -13, -11,  -4,  -5,  -7],
    /* K */ [-11,  -6,  -2, -11,  -7,  -3,  -2, -11,  -6,  -7,  -5,  -7],
    /* M */ [ -3,  -7, -11,  -2,  -6, -11, -11,  -2,  -6,  -7,  -5,  -7],
    /* S */ [ -9,  -5,  -2,  -2,  -5,  -9,  -5,  -5,  -2,  -9,  -5,  -7],
    /* W */ [ -4,  -8, -11, -11,  -8,  -4,  -8,  -8, -11,  -4,  -5,  -7],
    /* N */ [ -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -1,  -7],
    /* X */ [ -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7],
];

static IDX: [u8; 256] = {
    let mut t = [255u8; 256];
    let mut i = 0;
    while i < ALPHABET.len() {
        let c = ALPHABET[i];
        t[c as usize] = i as u8;
        t[(c + 32) as usize] = i as u8;
        i += 1;
    }
    t
};

#[inline]
fn idx(b: u8) -> Option<usize> {
    let i = IDX[b as usize];
    (i != 255).then_some(i as usize)
}

/// Port of `ScoreMatrix::getScore(c1, c2)`.
///
/// The ordering matters and is faithful: padding short-circuits to 0 *before*
/// any gap handling. Unknown symbols raise `logic_error` in the C++, which
/// `getConsensus` catches and then leaves the running score untouched — so they
/// contribute 0 here.
#[inline]
pub fn get_score(candidate: u8, observed: u8) -> i32 {
    if candidate == PAD || observed == PAD {
        return 0;
    }
    if candidate == GAP || observed == GAP {
        if candidate != GAP || observed != GAP {
            return GAP_EXT_PENALTY;
        }
        return GAP_MATCH_SCORE;
    }
    match (idx(candidate), idx(observed)) {
        (Some(r), Some(c)) => DNACONMATRIX[r][c],
        _ => 0,
    }
}

/// Port of `MultipleAlignment::getConsensus()`.
///
/// Returns the **gapped** consensus, one byte per alignment column.
///
/// `min_non_gap_count` is `acons --min`: a column with fewer than this many
/// non-gap, non-padding residues is forced to a gap whatever the scores say.
/// `autocons` uses 1 for candidate scoring and `--min` (default 2) for
/// refinement, matching the C++'s two call sites.
///
/// Width is taken from the *first* sequence, as the C++'s
/// `alnWidth = (*beg).size()` does. Positions past the end of a shorter row are
/// treated as padding rather than read out of bounds.
// `x` indexes both the output and every input row at the same column; an
// iterator over `cons` would hide that correspondence, which is the recurrence.
#[allow(clippy::needless_range_loop)]
pub fn get_consensus(sequences: &[&[u8]], min_non_gap_count: usize) -> Vec<u8> {
    if sequences.is_empty() {
        return Vec::new();
    }
    let width = sequences[0].len();
    let mut cons = vec![GAP; width];

    // Candidates are the alphabet plus gap appended at index 12.  The C++ scans
    // in reverse with a strict `>`, so gap is tried first and wins any tie.
    let n_cand = ALPHABET.len() + 1;
    let cand_byte = |i: usize| if i == ALPHABET.len() { GAP } else { ALPHABET[i] };

    for x in 0..width {
        let mut max_score: i64 = -1_000_000;
        let mut max_char = GAP;

        for i in (0..n_cand).rev() {
            let cand = cand_byte(i);
            let mut score: i64 = 0;
            let mut count: usize = 0;

            for seq in sequences {
                let cell = seq.get(x).copied().unwrap_or(PAD);
                if cell != GAP && cell != PAD {
                    count += 1;
                }
                score += i64::from(get_score(cand, cell));
            }

            if score > max_score && count >= min_non_gap_count {
                max_score = score;
                max_char = cand;
            }
        }
        cons[x] = max_char;
    }
    cons
}

/// Reproduce `acons -t`: rewrite flanking gaps as padding.
///
/// GIRI does this at read time (`MultipleAlignment::read(is, trimGaps)`), so
/// flanking gaps score 0 and are excluded from the coverage count.
pub fn trim_flanking_gaps(seq: &[u8]) -> Vec<u8> {
    let mut out = seq.to_vec();
    for b in out.iter_mut() {
        if *b == GAP {
            *b = PAD
        } else {
            break
        }
    }
    for b in out.iter_mut().rev() {
        if *b == GAP {
            *b = PAD
        } else {
            break
        }
    }
    out
}

/// Strip gaps, as `acons` does before output (`keep_gaps = false`).
pub fn strip_gaps(cons: &[u8]) -> Vec<u8> {
    cons.iter().copied().filter(|&b| b != GAP).collect()
}

impl crate::msa::MultiAlign {
    /// Call a consensus with the GIRI caller — `acons --orig`.
    ///
    /// Excludes the reference row, matching `maln.erase(maln.begin())`.
    pub fn giri_consensus(&self, min_non_gap_count: usize) -> Vec<u8> {
        let rows: Vec<&[u8]> = self.sequences[1.min(self.sequences.len())..]
            .iter()
            .map(|r| r.seq.as_slice())
            .collect();
        get_consensus(&rows, min_non_gap_count)
    }
}


// ═══════════════════════════════════════════════════════════════════════════
// CpG restoration (`acons --mam`)
// ═══════════════════════════════════════════════════════════════════════════
//
// Port of `MultipleAlignment::restoreCpG()`.  Unlike the deterministic
// bonus model in `crate::consensus`, this is a species-aware statistical
// test: it estimates the alignment's background transition rate from the
// observed similarity, derives the expected rate of each CpG deamination
// product, and restores `CG` only when binomial tests on the observed
// counts clear a 0.05 threshold.
//
// `sim` is computed with the **alignment** matrix (`DNAALNMATRIX`), not the
// consensus matrix: `acons` never calls `PairwiseAlignment::loadScoreMatrix`,
// so the static `scoreMatrix` is default-constructed, and that constructor
// installs `DNAALNMATRIX`.

/// `DNAALNMATRIX` alphabet.  NOTE the order differs from [`ALPHABET`] —
/// `M W S K` here vs. `K M S W` in `DNACONMATRIX`.
const ALN_ALPHABET: [u8; 12] = *b"ARGCYTMWSKNX";

/// `DNAALNMATRIX`; rows = `c1` (the *sequence* base), cols = `c2` (the
/// *consensus* base).  The matrix is asymmetric (A→G = −5 but G→A = −7), so
/// argument order matters: GIRI calls `getScore(*botPtr, *topPtr)`.
#[rustfmt::skip]
static DNAALNMATRIX: [[i32; 12]; 12] = [
    //        A    R    G    C    Y    T    M    W    S    K    N    X
    /* A */ [  9,   3,  -5, -12, -12, -12,   0,   0,  -9,  -9,  -1,  -3],
    /* R */ [  3,   3,   3, -12, -12, -12,  -6,  -6,  -6,  -6,  -1,  -3],
    /* G */ [ -7,   3,   9, -12, -12, -12,  -9,  -9,   0,   0,  -1,  -3],
    /* C */ [-12, -12, -12,   9,   3,  -7,   0,  -9,   0,  -9,  -1,  -3],
    /* Y */ [-12, -12, -12,   3,   3,   3,  -6,  -6,  -6,  -6,  -1,  -3],
    /* T */ [-12, -12, -12,  -5,   3,   9,  -9,   0,  -9,   0,  -1,  -3],
    /* M */ [  0,  -6,  -9,   0,  -6,  -9,   0,  -6,  -6,  -9,  -1,  -3],
    /* W */ [  0,  -6,  -9,  -9,  -6,   0,  -6,   0,  -9,  -6,  -1,  -3],
    /* S */ [ -9,  -6,   0,   0,  -6,  -9,  -6,  -9,   0,  -6,  -1,  -3],
    /* K */ [ -9,  -6,   0,  -9,  -6,   0,  -9,  -6,  -6,   0,  -1,  -3],
    /* N */ [ -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -3],
    /* X */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3],
];

static ALN_IDX: [u8; 256] = {
    let mut t = [255u8; 256];
    let mut i = 0;
    while i < 12 {
        let c = ALN_ALPHABET[i];
        t[c as usize] = i as u8;
        t[(c + 32) as usize] = i as u8;
        i += 1;
    }
    t
};

#[inline]
fn aln_idx(b: u8) -> Option<usize> {
    let i = ALN_IDX[b as usize];
    if i == 255 { None } else { Some(i as usize) }
}

/// `getScore` for the alignment matrix.  `None` mirrors the C++ `logic_error`
/// that `accumulateAlignmentStats` catches (and responds to with `continue`).
#[inline]
fn aln_score(c1: u8, c2: u8) -> Option<i32> {
    match (aln_idx(c1), aln_idx(c2)) {
        (Some(r), Some(c)) => Some(DNAALNMATRIX[r][c]),
        _ => None,
    }
}

#[inline]
fn is_pad(b: u8) -> bool {
    // GIRI LPADDING/RPADDING; this crate additionally uses b' '.
    b == PAD || b == b'<' || b == b'>'
}

/// Port of `firstBase` — including a defect in the original that materially
/// changes behaviour.
///
/// The C++ is:
/// ```text
/// static const char index[] = "1000...";   // characters '0'/'1', not integers
/// if (*ptr >= '*' && *ptr <= 'z' && index[*ptr - '*']) break;
/// ```
/// `index[...]` evaluates to the *character* `'0'` (48) or `'1'` (49) — both
/// non-zero, hence always truthy.  The lookup table is therefore inert, and
/// the function really stops at the first byte in `'*'..='z'`.  That range
/// includes `-` (45) and `.` (46), so gaps are **not** skipped.
///
/// The visible consequence is that the `CG` mask in `accumulateAlignmentStats`
/// only matches a *contiguous* `CG`; a `C-G` spanning a gap does not match.
/// Reproducing this is required for `sim` — and therefore every `--mam`
/// decision — to agree with `acons`.
fn first_base(s: &[u8], mut i: usize) -> usize {
    while i < s.len() {
        let c = s[i];
        if (b'*'..=b'z').contains(&c) {
            return i;
        }
        i += 1;
    }
    i
}

/// Port of `accumulateAlignmentStats`, restricted to the match/mismatch
/// counters that feed `sim`.  `top` is the consensus, `bot` one sequence.
///
/// Faithfully reproduces two subtleties:
/// * a **run** of gaps counts as a single mismatch (not one per column);
/// * a masked `CG` motif in `bot` skips the motif **and the following base**,
///   because the C++ `continue` also fires the loop increment.
fn accumulate_stats(top: &[u8], bot: &[u8], m_cnt: &mut u64, mm_cnt: &mut u64) {
    const MASK: &[u8; 2] = b"CG";
    let mut c_gap = 0u32;
    let mut s_gap = 0u32;
    let (mut ti, mut bi) = (0usize, 0usize);

    while ti < top.len() && bi < bot.len() {
        let (t, b) = (top[ti], bot[bi]);
        if b == GAP && t == GAP {
            ti += 1; bi += 1; continue;
        }
        if is_pad(b) || is_pad(t) {
            ti += 1; bi += 1; continue;
        }

        // CG-motif mask, scanned along `bot` skipping non-bases.
        let mut skip = false;
        let mut s1 = bi;
        let mut s2 = 0usize;
        while s1 < bot.len() && s2 < MASK.len() {
            if bot[s1].to_ascii_uppercase() != MASK[s2] {
                break;
            }
            s1 = first_base(bot, s1 + 1);
            s2 += 1;
        }
        if s2 == MASK.len() {
            skip = true;
            while bi != s1 {
                ti += 1;
                bi += 1;
            }
        }
        if ti >= top.len() || bi >= bot.len() {
            break;
        }
        if skip {
            ti += 1; bi += 1; continue;
        }

        let (t, b) = (top[ti], bot[bi]);
        if b != t {
            if b == GAP { c_gap += 1; ti += 1; bi += 1; continue; }
            if t == GAP { s_gap += 1; ti += 1; bi += 1; continue; }
        } else if b == GAP {
            ti += 1; bi += 1; continue;
        }

        // C++ throws on unknown symbols; the catch does `continue`, leaving
        // the pending gap-run counters untouched.
        let Some(score) = aln_score(b, t) else {
            ti += 1; bi += 1; continue;
        };

        if c_gap > 0 { *mm_cnt += 1; c_gap = 0; }
        if s_gap > 0 { *mm_cnt += 1; s_gap = 0; }

        if aln_idx(t) == aln_idx(b) && score > 0 {
            *m_cnt += 1;
        } else {
            *mm_cnt += 1;
        }
        ti += 1;
        bi += 1;
    }
}

// ── Deamination probability model ────────────────────────────────────────────
// `sc` is the similarity; `p_ts` the per-site transition probability.

/// Per-site transition probability.  The `2/3` hard-codes a mammalian
/// ts/tv of 2 and is independent of any scoring matrix.
fn pts_mam(sim: f64) -> f64 { 2.0 / 3.0 * (1.0 - sim) }

fn p_ca_tg(p: f64, sc: f64) -> f64 { p * p + p * sc / 2.0 }
fn p_tg_ca(p: f64, sc: f64) -> f64 { p_ca_tg(p, sc) }
fn p_ca_ta(p: f64, sc: f64) -> f64 { p * sc + p * p / 2.0 }
fn p_tg_cg(_p: f64, sc: f64) -> f64 { sc * sc / 2.0 }
fn p_ca_cg(p: f64, sc: f64) -> f64 { p_tg_cg(p, sc) }
fn p_ta_yr(p: f64, sc: f64) -> f64 { sc * sc / 4.0 + p * sc }
fn p_tg_cr(p: f64, sc: f64) -> f64 { (p_tg_ca(p, sc) + p_tg_cg(p, sc)) / 2.0 }
fn p_ca_tr(p: f64, sc: f64) -> f64 { (p_ca_ta(p, sc) + p_ca_tg(p, sc)) / 2.0 }
fn p_tg_ya(p: f64, sc: f64) -> f64 { (p_ca_ta(p, sc) + p_tg_ca(p, sc)) / 2.0 }
fn p_ca_yg(p: f64, sc: f64) -> f64 { (p_ca_tg(p, sc) + p_ca_cg(p, sc)) / 2.0 }

// ── Exact binomial tails ─────────────────────────────────────────────────────
// Computed by direct summation in log space rather than via a statistics
// crate, so that boost/GSL last-bit differences cannot flip a decision at the
// p = 0.05 threshold.

fn binom_sum(n: usize, p: f64, lo: usize, hi: usize) -> f64 {
    if lo > hi {
        return 0.0;
    }
    if p <= 0.0 {
        return if lo == 0 { 1.0 } else { 0.0 };
    }
    if p >= 1.0 {
        return if hi >= n { 1.0 } else { 0.0 };
    }
    let lnp = p.ln();
    let ln1mp = (1.0 - p).ln();
    let mut ln_binom = 0.0f64; // ln C(n,0)
    let mut sum = 0.0f64;
    for k in 0..=n {
        if k > 0 {
            ln_binom += ((n - k + 1) as f64).ln() - (k as f64).ln();
        }
        if k >= lo && k <= hi {
            sum += (ln_binom + k as f64 * lnp + (n - k) as f64 * ln1mp).exp();
        }
    }
    sum
}

/// `P(X <= k)` — GSL `gsl_cdf_binomial_P`.
fn binom_cdf_le(k: usize, n: usize, p: f64) -> f64 {
    binom_sum(n, p, 0, k.min(n))
}

/// `P(X >= k)` — GSL `gsl_cdf_binomial_Q(k) + gsl_ran_binomial_pdf(k)`.
fn binom_sf_ge(k: usize, n: usize, p: f64) -> f64 {
    if k > n { 0.0 } else { binom_sum(n, p, k, n) }
}

/// Port of `MultipleAlignment::restoreCpG()` for `species != UNKNOWN`.
///
/// Operates on the **gapped** consensus in place, exactly where `acons` calls
/// it — after `getConsensus` and before gaps are stripped.
pub fn restore_cpg(cons: &mut [u8], sequences: &[&[u8]]) {
    for b in cons.iter_mut() {
        *b = b.to_ascii_uppercase();
    }

    let maln_size = sequences.len();
    let (mut m_cnt, mut mm_cnt) = (0u64, 0u64);
    for s in sequences {
        accumulate_stats(cons, s, &mut m_cnt, &mut mm_cnt);
    }
    // 0/0 yields NaN in C++ too; NaN then fails every `p > 0.0` guard below,
    // which skips the tests entirely — reproduced here by construction.
    let sim = m_cnt as f64 / (m_cnt + mm_cnt) as f64;
    let p_ts = pts_mam(sim);
    let dbg = std::env::var_os("CALLCONS_DEBUG").is_some();
    if dbg {
        eprintln!("SIM={:.6} (iMCnt={} iMMCnt={}) pTs={:.6}", sim, m_cnt, mm_cnt, p_ts);
    }

    let mut i = 0usize;
    while i < cons.len() {
        while i < cons.len() && !cons[i].is_ascii_alphabetic() {
            i += 1;
        }
        if i >= cons.len() {
            break;
        }
        let mut k = i + 1;
        while k < cons.len() && !cons[k].is_ascii_alphabetic() {
            k += 1;
        }
        if k >= cons.len() {
            break;
        }

        let (c1, c2) = (cons[i], cons[k]);
        let (mut p_tg, mut p_ca, mut p_ta) = (0.0f64, 0.0f64, 0.0f64);
        match (c1, c2) {
            (b'T', b'G') => p_ca = p_ca_tg(p_ts, sim),
            // NOTE: the C++ assigns both from pCA_TA for this doublet.
            (b'T', b'A') => { let v = p_ca_ta(p_ts, sim); p_tg = v; p_ca = v; }
            (b'T', b'R') => p_ca = p_ca_tr(p_ts, sim),
            (b'C', b'A') => p_tg = p_tg_ca(p_ts, sim),
            (b'C', b'R') => p_tg = p_tg_cr(p_ts, sim),
            (b'Y', b'G') => p_ca = p_ca_yg(p_ts, sim),
            (b'Y', b'A') => p_tg = p_tg_ya(p_ts, sim),
            (b'Y', b'R') => p_ta = p_ta_yr(p_ts, sim),
            _ => { i += 1; continue; }
        }

        // Counts are taken from the raw sequence bytes and compared against
        // uppercase literals, exactly as the C++ does — lowercase bases
        // therefore do NOT register as deamination products.
        let (mut n_tg, mut n_ca, mut n_ta, mut n_del) = (0usize, 0usize, 0usize, 0usize);
        for s in sequences {
            let a = s.get(i).copied().unwrap_or(PAD);
            let b = s.get(k).copied().unwrap_or(PAD);
            if a == b'T' && b == b'G' { n_tg += 1; continue; }
            if a == b'C' && b == b'A' { n_ca += 1; continue; }
            if a == b'T' && b == b'A' { n_ta += 1; continue; }
            if a == b'C' && b == b'G' { continue; }
            if a == GAP || b == GAP || is_pad(a) || is_pad(b) {
                n_del += 1;
            }
        }
        let n = maln_size.saturating_sub(n_del);
        let doublet = [c1, c2];

        if p_tg > 0.0 {
            let reject = if &doublet == b"CR" || &doublet == b"YA" {
                binom_cdf_le(n_tg, n, p_tg) < 0.05
            } else {
                binom_sf_ge(n_tg, n, p_tg) > 0.05
            };
            if reject { i += 1; continue; }
        }
        if p_ca > 0.0 {
            let reject = if &doublet == b"TR" || &doublet == b"YG" {
                binom_cdf_le(n_ca, n, p_ca) < 0.05
            } else {
                binom_sf_ge(n_ca, n, p_ca) > 0.05
            };
            if reject { i += 1; continue; }
        }
        if p_ta > 0.0 && binom_sf_ge(n_ta, n, p_ta) <= 0.05 {
            i += 1;
            continue;
        }

        cons[i] = b'C';
        cons[k] = b'G';
        // The C++ bumps sptr1 here *and* again via the loop increment.
        i += 2;
    }
}

pub mod fixed {
    use super::*;

    /// Fixed `firstBase`: skip non-base bytes, stop at the first nucleotide
    /// letter.  Contrast [`super::first_base`], whose lookup table is inert.
    fn first_base(s: &[u8], mut i: usize) -> usize {
        while i < s.len() && !s[i].is_ascii_alphabetic() {
            i += 1;
        }
        i
    }

    /// Fixed `accumulateAlignmentStats`: the `CG` mask now skips gaps (so
    /// `C-G` masks) and consumes only the motif span, not a trailing base.
    fn accumulate_stats(top: &[u8], bot: &[u8], m_cnt: &mut u64, mm_cnt: &mut u64) {
        const MASK: &[u8; 2] = b"CG";
        let mut c_gap = 0u32;
        let mut s_gap = 0u32;
        let (mut ti, mut bi) = (0usize, 0usize);

        while ti < top.len() && bi < bot.len() {
            let (t, b) = (top[ti], bot[bi]);
            if b == GAP && t == GAP {
                ti += 1; bi += 1; continue;
            }
            if is_pad(b) || is_pad(t) {
                ti += 1; bi += 1; continue;
            }

            // CG-motif mask along `bot`, correctly skipping gaps between bases.
            let mut s1 = bi;
            let mut s2 = 0usize;
            let mut last = bi; // column of the most recently matched motif base
            while s1 < bot.len() && s2 < MASK.len() {
                if bot[s1].to_ascii_uppercase() != MASK[s2] {
                    break;
                }
                last = s1;
                s1 = first_base(bot, s1 + 1);
                s2 += 1;
            }
            if s2 == MASK.len() {
                // Mask exactly [bi ..= last]: the C, any gap between, and the G.
                let end = last + 1;
                while bi < end {
                    ti += 1;
                    bi += 1;
                }
                continue;
            }

            let (t, b) = (top[ti], bot[bi]);
            if b != t {
                if b == GAP { c_gap += 1; ti += 1; bi += 1; continue; }
                if t == GAP { s_gap += 1; ti += 1; bi += 1; continue; }
            } else if b == GAP {
                ti += 1; bi += 1; continue;
            }

            let Some(score) = aln_score(b, t) else {
                ti += 1; bi += 1; continue;
            };

            if c_gap > 0 { *mm_cnt += 1; c_gap = 0; }
            if s_gap > 0 { *mm_cnt += 1; s_gap = 0; }

            if aln_idx(t) == aln_idx(b) && score > 0 {
                *m_cnt += 1;
            } else {
                *mm_cnt += 1;
            }
            ti += 1;
            bi += 1;
        }
    }

    /// `restoreCpG` using the fixed `sim`.  Body is otherwise identical to
    /// [`super::restore_cpg`] (same probability model and binomial tests).
    pub fn restore_cpg(cons: &mut [u8], sequences: &[&[u8]]) {
        for b in cons.iter_mut() {
            *b = b.to_ascii_uppercase();
        }

        let maln_size = sequences.len();
        let (mut m_cnt, mut mm_cnt) = (0u64, 0u64);
        for s in sequences {
            accumulate_stats(cons, s, &mut m_cnt, &mut mm_cnt);
        }
        let sim = m_cnt as f64 / (m_cnt + mm_cnt) as f64;
        let p_ts = pts_mam(sim);
        let dbg = std::env::var_os("CALLCONS_DEBUG").is_some();
        if dbg {
            eprintln!("SIM={:.6} (iMCnt={} iMMCnt={}) pTs={:.6} [FIXED]", sim, m_cnt, mm_cnt, p_ts);
        }

        let mut i = 0usize;
        while i < cons.len() {
            while i < cons.len() && !cons[i].is_ascii_alphabetic() {
                i += 1;
            }
            if i >= cons.len() {
                break;
            }
            let mut k = i + 1;
            while k < cons.len() && !cons[k].is_ascii_alphabetic() {
                k += 1;
            }
            if k >= cons.len() {
                break;
            }

            let (c1, c2) = (cons[i], cons[k]);
            let (mut p_tg, mut p_ca, mut p_ta) = (0.0f64, 0.0f64, 0.0f64);
            match (c1, c2) {
                (b'T', b'G') => p_ca = p_ca_tg(p_ts, sim),
                (b'T', b'A') => { let v = p_ca_ta(p_ts, sim); p_tg = v; p_ca = v; }
                (b'T', b'R') => p_ca = p_ca_tr(p_ts, sim),
                (b'C', b'A') => p_tg = p_tg_ca(p_ts, sim),
                (b'C', b'R') => p_tg = p_tg_cr(p_ts, sim),
                (b'Y', b'G') => p_ca = p_ca_yg(p_ts, sim),
                (b'Y', b'A') => p_tg = p_tg_ya(p_ts, sim),
                (b'Y', b'R') => p_ta = p_ta_yr(p_ts, sim),
                _ => { i += 1; continue; }
            }

            let (mut n_tg, mut n_ca, mut n_ta, mut n_del) = (0usize, 0usize, 0usize, 0usize);
            for s in sequences {
                let a = s.get(i).copied().unwrap_or(PAD);
                let b = s.get(k).copied().unwrap_or(PAD);
                if a == b'T' && b == b'G' { n_tg += 1; continue; }
                if a == b'C' && b == b'A' { n_ca += 1; continue; }
                if a == b'T' && b == b'A' { n_ta += 1; continue; }
                if a == b'C' && b == b'G' { continue; }
                if a == GAP || b == GAP || is_pad(a) || is_pad(b) {
                    n_del += 1;
                }
            }
            let n = maln_size.saturating_sub(n_del);
            let doublet = [c1, c2];

            if p_tg > 0.0 {
                let reject = if &doublet == b"CR" || &doublet == b"YA" {
                    binom_cdf_le(n_tg, n, p_tg) < 0.05
                } else {
                    binom_sf_ge(n_tg, n, p_tg) > 0.05
                };
                if reject { i += 1; continue; }
            }
            if p_ca > 0.0 {
                let reject = if &doublet == b"TR" || &doublet == b"YG" {
                    binom_cdf_le(n_ca, n, p_ca) < 0.05
                } else {
                    binom_sf_ge(n_ca, n, p_ca) > 0.05
                };
                if reject { i += 1; continue; }
            }
            if p_ta > 0.0 && binom_sf_ge(n_ta, n, p_ta) <= 0.05 {
                i += 1;
                continue;
            }

            cons[i] = b'C';
            cons[k] = b'G';
            i += 2;
        }
    }

    #[cfg(test)]
    mod tests {

        use super::*;

        #[test]
        fn fixed_first_base_skips_gaps_and_dots() {
            assert_eq!(first_base(b"C-G", 1), 2, "'-' is skipped");
            assert_eq!(first_base(b"C.G", 1), 2, "'.' is skipped");
            assert_eq!(first_base(b"C--G", 1), 3, "run of gaps skipped");
            assert_eq!(first_base(b"C--", 1), 3, "runs off the end when no base");
        }

        /// A `C-G` spanning a gap is masked by the fixed caller (it is a CpG),
        /// so it does not inflate the mismatch count the way the buggy port
        /// does — the two produce different `sim`.
        #[test]
        fn gapped_cpg_is_masked() {
            let top = b"C-G";
            let bot = b"C-G";
            let (mut m, mut mm) = (0u64, 0u64);
            accumulate_stats(top, bot, &mut m, &mut mm);
            // whole motif masked → nothing counted
            assert_eq!((m, mm), (0, 0));
        }
    }
}

#[cfg(test)]
mod tests {
    /// Soft-masked input must call the same consensus as uppercase. Replaces a
    /// test lost when dfam-curator's duplicate `giri.rs` was deleted — that file
    /// was untracked, so its cases could not be recovered.
    #[test]
    fn lowercase_input_is_accepted() {
        let upper: Vec<&[u8]> = vec![b"ACGTRYKM", b"ACGTRYKM", b"ACGTRYKM"];
        let lower: Vec<&[u8]> = vec![b"acgtrykm", b"acgtrykm", b"acgtrykm"];
        assert_eq!(
            get_consensus(&upper, 1),
            get_consensus(&lower, 1),
            "case must not change the call"
        );
    }

    /// A column split evenly between C and G resolves to the IUPAC ambiguity
    /// code S, not to either base.
    #[test]
    fn a_balanced_c_g_column_calls_s() {
        let rows: Vec<&[u8]> = vec![b"C", b"G"];
        assert_eq!(get_consensus(&rows, 1), b"S");
    }

    use super::*;

    #[test]
    fn score_ordering_puts_padding_before_gaps() {
        // Padding wins the short-circuit even against a gap.
        assert_eq!(get_score(PAD, GAP), 0);
        assert_eq!(get_score(GAP, PAD), 0);
        assert_eq!(get_score(GAP, GAP), 2);
        assert_eq!(get_score(GAP, b'A'), -5);
        assert_eq!(get_score(b'A', GAP), -5);
    }

    #[test]
    fn matrix_lookups_match_the_table() {
        assert_eq!(get_score(b'A', b'A'), 9);
        assert_eq!(get_score(b'C', b'C'), 10);
        assert_eq!(get_score(b'N', b'N'), -1);
        // Asymmetric: candidate G against observed A is not the transpose.
        assert_eq!(get_score(b'G', b'A'), -4);
        assert_eq!(get_score(b'A', b'G'), -8);
    }

    #[test]
    fn unknown_symbols_contribute_nothing() {
        // The C++ throws and getConsensus swallows it, leaving score unchanged.
        assert_eq!(get_score(b'A', b'@'), 0);
        assert_eq!(get_score(b'@', b'A'), 0);
    }

    #[test]
    fn calls_the_obvious_column() {
        let seqs: Vec<&[u8]> = vec![b"AAA", b"AAA", b"AAA"];
        assert_eq!(get_consensus(&seqs, 1), b"AAA");
    }

    /// The gate is independent of the candidate: `count` counts residues in the
    /// column, so when it falls short *no* candidate is accepted and the
    /// initial gap survives.
    #[test]
    fn the_coverage_gate_forces_a_gap() {
        // Two residues against one gap — `A` wins on score (9+9-5 = 13 against
        // the gap candidate's -5-5+2 = -8).
        let seqs: Vec<&[u8]> = vec![b"A", b"A", b"-"];
        assert_eq!(get_consensus(&seqs, 2), b"A", "two residues clear min=2");
        assert_eq!(get_consensus(&seqs, 3), b"-", "but not min=3");
    }

    /// Separately from the gate, a gap-dominated column is won by the gap
    /// candidate on score alone, because gap-vs-gap pays `+2`.
    #[test]
    fn a_gap_dominated_column_is_won_on_score_not_by_the_gate() {
        let seqs: Vec<&[u8]> = vec![b"A", b"-", b"-", b"-"];
        // Gap candidate: -5 + 2 + 2 + 2 = +1.  Candidate A: 9 - 5 - 5 - 5 = -6.
        assert_eq!(
            get_consensus(&seqs, 1),
            b"-",
            "the gap wins on score even though the gate passes at min=1"
        );
    }

    /// The load-bearing difference from the Dfam caller: reverse scan with a
    /// strict `>` means gap wins ties, where the Dfam caller prefers `N`.
    #[test]
    fn ties_resolve_towards_gap_not_n() {
        // An all-padding column scores 0 for every candidate, so it is a pure
        // tie — and the reverse scan hands it to the first candidate tried.
        let seqs: Vec<&[u8]> = vec![b" ", b" "];
        assert_eq!(
            get_consensus(&seqs, 0),
            b"-",
            "GIRI resolves a tie to gap; the Dfam caller would say N"
        );
    }

    #[test]
    fn width_comes_from_the_first_row() {
        // The C++ takes alnWidth from the first sequence; a longer row is
        // truncated and a shorter one is padded rather than read past its end.
        let seqs: Vec<&[u8]> = vec![b"AAA", b"AAAAAA"];
        assert_eq!(get_consensus(&seqs, 1).len(), 3);

        let seqs: Vec<&[u8]> = vec![b"AAAAAA", b"AAA"];
        let cons = get_consensus(&seqs, 1);
        assert_eq!(cons.len(), 6);
        assert_eq!(&cons[..3], b"AAA");
    }

    #[test]
    fn flanking_gaps_become_padding() {
        assert_eq!(trim_flanking_gaps(b"--ACGT--"), b"  ACGT  ");
        assert_eq!(trim_flanking_gaps(b"AC--GT"), b"AC--GT", "interior untouched");
    }

    #[test]
    fn strip_gaps_removes_only_gaps() {
        assert_eq!(strip_gaps(b"AC-GT"), b"ACGT");
        assert_eq!(strip_gaps(b" AC-GT "), b" ACGT ", "padding survives");
    }

    /// The two callers genuinely disagree, so a comparison run has to hold the
    /// choice constant.
    #[test]
    fn giri_and_dfam_callers_differ_on_the_same_input() {
        use crate::consensus::{build_consensus_from_sequences, ConsensusParams};
        // A CpG that has decayed: the Dfam caller restores CG, GIRI does not.
        let seqs: Vec<&[u8]> = vec![b"TG", b"TG", b"TG", b"TG", b"TG", b"CA", b"CA", b"CG"];
        let dfam = build_consensus_from_sequences(&seqs, &ConsensusParams::default());
        let giri = get_consensus(&seqs, 1);
        assert_eq!(dfam, b"CG");
        assert_eq!(giri, b"TG");
        assert_ne!(dfam, giri);
    }

    #[test]
    fn first_base_does_not_skip_gaps_or_dots() {
        assert_eq!(first_base(b"C-G", 1), 1, "'-' (45) is in range and stops the scan");
        assert_eq!(first_base(b"C.G", 1), 1, "'.' (46) likewise");
        assert_eq!(first_base(b"C9G", 1), 1, "digits likewise");
        // Only bytes outside '*'..='z' are skipped — e.g. space (32).
        assert_eq!(first_base(b"C G", 1), 2);
        assert_eq!(first_base(b"C  ", 1), 3, "runs off the end when none qualify");
    }

    #[test]
    fn ambiguous_doublets_use_lower_tail_and_restore() {
        let rows: Vec<&[u8]> = vec![
            b"ATGACAATAATAATA",
            b"AKGACMATMAKAAKM",
            b"AMGACSATSAMAAMS",
            b"ASGACSATSASAASS",
        ];
        let base = get_consensus(&rows, 1);
        assert_eq!(
            String::from_utf8(strip_gaps(&base)).unwrap(),
            "AYGACRATRAYAAYR",
            "columns must resolve to Y/R to reach the ambiguous doublets"
        );

        let mut cons = base.clone();
        restore_cpg(&mut cons, &rows);
        assert_eq!(
            String::from_utf8(strip_gaps(&cons)).unwrap(),
            "ACGACGACGACGACG",
            "YG, CR, TR, YA and YR should each restore to CG"
        );
    }

    #[test]
    fn mam_restoration_matches_acons() {
        let rows: Vec<&[u8]> = vec![
            b"CTGgCCAACCaGATATTCCACCGC", b"C-GTCCAACCCAATATTCCAaCAC",
            b"CTGTCCAACCCAATATTCCACCAC", b"CTGTCCAACCC-ATATTCGGCTGC",
            b"CTGTCGAACCCAATATTCCACCAC", b"CCGTCCAAGCCAATTGTCCACTGc",
            b"CTcTCCAACCCAATATTCCACTGC", b"--GTCCAACCCAATATTCC-----",
            b"CCATCCAACCCAATATTCCACTGC", b"-TGTCC-ACCCAATATT-CACCGC",
            b"CTGTCCAACCCAATATTCC-CTGC", b"CCGTCCAACCTAATATTCCAcTGC",
            b"CTGTCTAACCCGATATTCTGCCAC", b"CTGTCCAACCTGATATTCCACTGC",
            b"CTGTCCAtCCTGGTATTCCACTGC", b"C-GTCCAACCTGATATTCCACTGC",
            b"CTGTCCAACCCAATATTCCAC-GC", b"CTGTCCAACCCGATATTC-GCTGC",
            b"-TATCCAACCCAA-ATTCCACTGC", b"TTGTCCAACCCAATATTCCACTGC",
            b"CTGTCCAACCCA-TATTCCACTGC", b"CTATCGAAGC-AA-ATTCCACTGC",
            b"CTGACCAACCCAATATCCCACTG-", b"CTGGCCAACCCAATATTCgACTGC",
            b"ACATC-AACCCAATATTCC-----", b"CCATCCAA-CCAATATTCTACCAC",
        ];
        // Plain consensus (no CpG restoration).
        let base = get_consensus(&rows, 1);
        assert_eq!(
            String::from_utf8(strip_gaps(&base)).unwrap(),
            "CTGTCCAACCCAATATTCCACTGC"
        );

        // sim must come out as acons computes it.
        let mut m = 0u64;
        let mut mm = 0u64;
        for s in &rows {
            accumulate_stats(&base, s, &mut m, &mut mm);
        }
        assert_eq!((m, mm), (519, 63), "iMCnt/iMMCnt must match acons");

        // Only the 21–22 doublet is restored, not 1–2.
        let mut cons = base.clone();
        restore_cpg(&mut cons, &rows);
        assert_eq!(
            String::from_utf8(strip_gaps(&cons)).unwrap(),
            "CTGTCCAACCCAATATTCCACCGC"
        );
    }

}
