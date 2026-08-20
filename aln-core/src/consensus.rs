//! The Dfam consensus caller.
//!
//! A port of Perl `MultAln.pm::buildConsensusFromArray`, by way of the verified
//! Rust in `dfam-curator`'s `consensus.rs`. It is the default caller in `acons`
//! and `autocons` (the original GIRI caller is reached with `--orig`).
//!
//! Public names deliberately match `dfam_curator::consensus` so that migrating
//! that crate onto `aln-core` is a delete-and-re-export rather than a rewrite.
//!
//! # Two passes
//!
//! 1. **Per-column argmax.** Every alphabet symbol — gap included — is scored as
//!    a candidate against the column's observed counts under [`MATRIX`]; the
//!    best wins, with `N` preferred on ties.
//! 2. **CpG restoration.** For each adjacent pair of non-gap consensus
//!    positions, the score of the called dinucleotide is compared against the
//!    score of a hypothetical `CG` under a deamination model. If `CG` wins, both
//!    positions are overwritten. This is what recovers CpG sites that mutated to
//!    `TG`/`CA` in most copies.
//!
//! # Scoring matrix
//!
//! This is *not* the alignment matrix ([`crate::SubstMatrix`]). It is an 18×18
//! table over 17 IUPAC symbols plus gap, with gap-vs-base `-6` and gap-vs-gap
//! `+3`, used only for consensus calling.

/// Byte → matrix index. 255 marks a symbol outside the alphabet.
static IDX: [u8; 256] = {
    let mut t = [255u8; 256];
    // A=0 R=1 G=2 C=3 Y=4 T=5 K=6 M=7 S=8 W=9 N=10 X=11 Z=12 V=13 H=14 D=15 B=16 -=17
    const SYMBOLS: &[u8; 17] = b"ARGCYTKMSWNXZVHDB";
    let mut i = 0;
    while i < SYMBOLS.len() {
        let c = SYMBOLS[i];
        t[c as usize] = i as u8;
        t[(c + 32) as usize] = i as u8; // lower case
        i += 1;
    }
    t[b'-' as usize] = 17;
    t
};

/// 17 IUPAC symbols plus gap.
pub const ALPHA_LEN: usize = 18;

/// Index of the gap row/column.
pub const GAP_IDX: usize = 17;

/// Index of `N`, the tie-break winner.
const N_IDX: usize = 10;

/// Matrix index for a byte, or `None` if it is outside the alphabet.
#[inline]
pub fn alpha_idx(b: u8) -> Option<usize> {
    let i = IDX[b as usize];
    (i != 255).then_some(i as usize)
}

/// Canonical upper-case byte for a matrix index.
pub fn alpha_byte(idx: usize) -> u8 {
    const LUT: [u8; ALPHA_LEN] = *b"ARGCYTKMSWNXZVHDB-";
    LUT[idx]
}

/// The 18×18 consensus scoring matrix; rows are the candidate, columns the
/// observed base. Both indexed by [`alpha_idx`].
#[rustfmt::skip]
pub static MATRIX: [[i32; ALPHA_LEN]; ALPHA_LEN] = [
    //       A    R    G    C    Y    T    K    M    S    W    N    X    Z    V    H    D    B   [-]
    /* A */ [ 9,   0,  -8, -15, -16, -17, -13,  -3, -11,  -4,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* R */ [ 2,   1,   1, -15, -15, -16,  -7,  -6,  -6,  -7,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* G */ [-4,   3,  10, -14, -14, -15,  -2,  -9,  -2,  -9,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* C */ [-15, -14, -14,  10,   3,  -4,  -9,  -2,  -2,  -9,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* Y */ [-16, -15, -15,   1,   1,   2,  -6,  -7,  -6,  -7,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* T */ [-17, -16, -15,  -8,   0,   9,  -3, -13, -11,  -4,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* K */ [-11,  -6,  -2, -11,  -7,  -3,  -2, -11,  -6,  -7,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* M */ [ -3,  -7, -11,  -2,  -6, -11, -11,  -2,  -6,  -7,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* S */ [ -9,  -5,  -2,  -2,  -5,  -9,  -5,  -5,  -2,  -9,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* W */ [ -4,  -8, -11, -11,  -8,  -4,  -8,  -8, -11,  -4,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* N */ [ -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -1,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* X */ [ -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* Z */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -6],
    /* V */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -6],
    /* H */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -6],
    /* D */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -6],
    /* B */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -6],
    /* - */ [ -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,   3],
];

/// Score a candidate against an observed base; 0 for unrecognised input.
#[inline]
pub fn score(candidate: u8, observed: u8) -> i32 {
    match (alpha_idx(candidate), alpha_idx(observed)) {
        (Some(r), Some(c)) => MATRIX[r][c],
        _ => 0,
    }
}

/// Tuning for [`build_consensus_from_sequences`].
///
/// Defaults match Perl `MultAln.pm::buildConsensusFromArray`.
#[derive(Debug, Clone)]
pub struct ConsensusParams {
    /// Count the reference row in the profile (`inclRef` in the Perl).
    ///
    /// `autocons` drops the reference before calling, so this is normally false.
    pub include_reference: bool,

    /// Bonus when the observed dinucleotide is `TG` or `CA` — the one-step CpG
    /// deamination products. Default 12.
    pub cg_param: i32,

    /// Applied when the observed dinucleotide is `TA`, the two-step deamination
    /// product. Default -5, i.e. a penalty.
    pub ta_param: i32,

    /// Bonus when the observed dinucleotide is a transition pair that could have
    /// arisen from a CpG (`TC`/`TT` forward, `AA`/`GA` reverse). Default 2.
    pub cg_trans_param: i32,

    /// Run the CpG restoration pass. Default true.
    ///
    /// Turn it off to get the plain per-column call, which is what to compare
    /// against callers that do no restoration (GIRI `acons` without `--mam`).
    pub enable_cpg: bool,
}

impl Default for ConsensusParams {
    fn default() -> Self {
        ConsensusParams {
            include_reference: false,
            cg_param: 12,
            ta_param: -5,
            cg_trans_param: 2,
            enable_cpg: true,
        }
    }
}

/// Call a gapped consensus from gapped aligned rows.
///
/// Rows should share a width. Padding (`' '`) marks columns a row does not
/// reach and is excluded from the profile; interior gaps (`-`) are counted and
/// can win a sparsely-covered column.
///
/// Leading and trailing gaps are rewritten to padding first, matching the
/// Perl's `MultAln::consensus()`, which pads before the first base and does not
/// extend past the last.
// The argmax scans index two parallel tables by the same alphabet position;
// iterators would obscure that correspondence, which is the whole point.
#[allow(clippy::needless_range_loop)]
pub fn build_consensus_from_sequences(
    sequences: &[&[u8]],
    params: &ConsensusParams,
) -> Vec<u8> {
    if sequences.is_empty() {
        return Vec::new();
    }
    let width = sequences.iter().map(|s| s.len()).max().unwrap_or(0);

    let processed: Vec<Vec<u8>> = sequences
        .iter()
        .map(|&seq| {
            let first = seq
                .iter()
                .position(|b| b.is_ascii_alphabetic())
                .unwrap_or(seq.len());
            let last = seq
                .iter()
                .rposition(|b| b.is_ascii_alphabetic())
                .unwrap_or(0);
            let mut out = seq.to_vec();
            for b in &mut out[..first] {
                if *b == b'-' {
                    *b = b' ';
                }
            }
            if first < seq.len() {
                for b in &mut out[last + 1..] {
                    if *b == b'-' {
                        *b = b' ';
                    }
                }
            }
            out
        })
        .collect();
    let seqs: Vec<&[u8]> = processed.iter().map(|v| v.as_slice()).collect();

    // ── Pass 1: per-column argmax ─────────────────────────────────────────
    let mut profile = vec![[0u32; ALPHA_LEN]; width];
    for &seq in &seqs {
        for (col, &b) in seq.iter().enumerate() {
            if b == b' ' {
                continue;
            }
            if let Some(idx) = alpha_idx(b.to_ascii_uppercase()) {
                profile[col][idx] += 1;
            }
        }
    }

    let mut consensus: Vec<u8> = Vec::with_capacity(width);
    for col in &profile {
        let mut max_score = i64::MIN;
        let mut best_idx = N_IDX;
        let mut n_score = i64::MIN;

        for cand in 0..ALPHA_LEN {
            let mut s: i64 = 0;
            for obs in 0..ALPHA_LEN {
                let cnt = col[obs] as i64;
                if cnt > 0 {
                    s += cnt * MATRIX[cand][obs] as i64;
                }
            }
            if cand == N_IDX {
                n_score = s;
            }
            if s > max_score {
                max_score = s;
                best_idx = cand;
            }
        }
        // The Perl prefers N on a tie with any other candidate.
        if best_idx != N_IDX && n_score == max_score {
            best_idx = N_IDX;
        }
        consensus.push(alpha_byte(best_idx));
    }

    if !params.enable_cpg {
        return consensus;
    }

    // ── Pass 2: CpG restoration ───────────────────────────────────────────
    let c_idx = alpha_idx(b'C').unwrap();
    let g_idx = alpha_idx(b'G').unwrap();

    let mut i = 0usize;
    'outer: while i + 1 < consensus.len() {
        if consensus[i] == b'-' {
            i += 1;
            continue;
        }
        // The partner is the next non-gap consensus position, however far away.
        let mut k = i + 1;
        loop {
            if k >= consensus.len() {
                break 'outer;
            }
            if consensus[k] != b'-' {
                break;
            }
            k += 1;
        }

        let (Some(cl_idx), Some(cr_idx)) =
            (alpha_idx(consensus[i]), alpha_idx(consensus[k]))
        else {
            i += 1;
            continue;
        };

        let mut dn_score: i64 = 0;
        let mut cg_score: i64 = 0;

        for &seq in &seqs {
            if i >= seq.len() {
                continue;
            }
            let hl_raw = seq[i];
            if hl_raw == b' ' {
                continue;
            }
            let hr_raw = if k < seq.len() { seq[k] } else { b' ' };
            if hr_raw == b' ' {
                continue;
            }
            let hl = hl_raw.to_ascii_uppercase();
            let hr = hr_raw.to_ascii_uppercase();

            dn_score += MATRIX[cl_idx][alpha_idx(hl).unwrap_or(N_IDX)] as i64;
            dn_score += MATRIX[cr_idx][alpha_idx(hr).unwrap_or(N_IDX)] as i64;

            match (hl, hr) {
                // One-step deamination on either strand.
                (b'C', b'A') | (b'T', b'G') => cg_score += params.cg_param as i64,
                // Two-step: both positions mutated.
                (b'T', b'A') => cg_score += params.ta_param as i64,
                // Transition at the C, forward strand.
                (b'T', b'C') | (b'T', b'T') => {
                    cg_score += params.cg_trans_param as i64
                        + MATRIX[g_idx][alpha_idx(hr).unwrap_or(N_IDX)] as i64
                }
                // Transition at the G, reverse strand.
                (b'A', b'A') | (b'G', b'A') => {
                    cg_score += params.cg_trans_param as i64
                        + MATRIX[c_idx][alpha_idx(hl).unwrap_or(N_IDX)] as i64
                }
                _ => {
                    cg_score += MATRIX[c_idx][alpha_idx(hl).unwrap_or(N_IDX)] as i64;
                    cg_score += MATRIX[g_idx][alpha_idx(hr).unwrap_or(N_IDX)] as i64;
                }
            }
        }

        if cg_score > dn_score {
            consensus[i] = b'C';
            consensus[k] = b'G';
        }
        i += 1;
    }

    consensus
}

impl crate::msa::MultiAlign {
    /// Call a consensus over this alignment's rows.
    ///
    /// Honours [`ConsensusParams::include_reference`]; `autocons` leaves it
    /// false so only the aligned instances vote.
    pub fn consensus(&self, params: &ConsensusParams) -> Vec<u8> {
        let start = if params.include_reference { 0 } else { 1 };
        let rows: Vec<&[u8]> = self.sequences[start.min(self.sequences.len())..]
            .iter()
            .map(|r| r.seq.as_slice())
            .collect();
        build_consensus_from_sequences(&rows, params)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_conserved_cg_survives() {
        let seqs: Vec<&[u8]> = vec![b"CG", b"CG", b"CG"];
        assert_eq!(
            build_consensus_from_sequences(&seqs, &ConsensusParams::default()),
            b"CG"
        );
    }

    /// The dataset has to be chosen with care: because `C` and `G` have the
    /// highest self-scores in [`MATRIX`], a balanced mix of `TG` and `CA`
    /// already calls `CG` from the per-column pass alone, so it proves nothing
    /// about restoration.  Here the plain call gives `TG` and only the
    /// deamination model recovers `CG`.
    #[test]
    fn deaminated_cpg_is_restored() {
        let seqs: Vec<&[u8]> = vec![b"TG", b"TG", b"TG", b"TG", b"TG", b"CA", b"CA", b"CG"];

        let plain = ConsensusParams { enable_cpg: false, ..Default::default() };
        assert_eq!(
            build_consensus_from_sequences(&seqs, &plain),
            b"TG",
            "per-column call alone should land on TG"
        );

        assert_eq!(
            build_consensus_from_sequences(&seqs, &ConsensusParams::default()),
            b"CG",
            "restoration should recover the ancestral CpG"
        );
    }

    /// A balanced TG/CA mix — the shape one reaches for intuitively — is called
    /// `CG` by the per-column pass on its own.  Pinned so nobody mistakes such a
    /// dataset for evidence that restoration works.
    #[test]
    fn a_balanced_tg_ca_mix_needs_no_restoration() {
        let seqs: Vec<&[u8]> =
            vec![b"TG", b"TG", b"TG", b"CA", b"CA", b"CA", b"CG", b"CG"];
        let plain = ConsensusParams { enable_cpg: false, ..Default::default() };
        assert_eq!(build_consensus_from_sequences(&seqs, &plain), b"CG");
    }

    #[test]
    fn unrelated_columns_are_not_turned_into_cpg() {
        let seqs: Vec<&[u8]> = vec![b"AC", b"AC", b"AC", b"AC"];
        let cons = build_consensus_from_sequences(&seqs, &ConsensusParams::default());
        assert_eq!(cons, b"AC", "a well-supported AC must survive restoration");
    }

    #[test]
    fn restoration_looks_past_an_intervening_gap_column() {
        // The partner position is the next non-gap consensus column.
        let seqs: Vec<&[u8]> = vec![b"T-G", b"T-G", b"T-G", b"C-A", b"C-A", b"C-A"];
        let cons = build_consensus_from_sequences(&seqs, &ConsensusParams::default());
        assert_eq!(cons[1], b'-', "the all-gap column stays a gap");
        assert_eq!(cons[0], b'C');
        assert_eq!(cons[2], b'G');
    }

    #[test]
    fn padding_is_excluded_but_interior_gaps_vote() {
        // Two rows cover the whole width; two only the middle.  The flanking
        // columns must be called from the covering rows alone.
        let seqs: Vec<&[u8]> = vec![b"AAAA", b"AAAA", b"  AA", b"  AA"];
        let cons = build_consensus_from_sequences(&seqs, &ConsensusParams::default());
        assert_eq!(cons, b"AAAA");
    }

    #[test]
    fn a_gap_dominant_column_calls_a_gap() {
        let seqs: Vec<&[u8]> = vec![b"A-A", b"A-A", b"A-A", b"AGA"];
        let cons = build_consensus_from_sequences(&seqs, &ConsensusParams::default());
        assert_eq!(cons[1], b'-');
    }

    #[test]
    fn leading_and_trailing_gaps_become_padding() {
        // A row whose gaps flank its bases must not drag those columns to '-'.
        let seqs: Vec<&[u8]> = vec![b"--AA--", b"CCAACC", b"CCAACC"];
        let cons = build_consensus_from_sequences(&seqs, &ConsensusParams::default());
        assert_eq!(&cons[..2], b"CC", "flanking gaps should not vote");
        assert_eq!(&cons[4..], b"CC");
    }

    #[test]
    fn n_wins_ties() {
        // An empty profile leaves every candidate at 0; N must be chosen.
        let seqs: Vec<&[u8]> = vec![b" "];
        let cons = build_consensus_from_sequences(&seqs, &ConsensusParams::default());
        assert_eq!(cons, b"N");
    }

    #[test]
    fn score_lookup_matches_the_table() {
        assert_eq!(score(b'A', b'A'), 9);
        assert_eq!(score(b'C', b'C'), 10);
        assert_eq!(score(b'-', b'-'), 3);
        assert_eq!(score(b'-', b'A'), -6);
        assert_eq!(score(b'@', b'A'), 0);
    }
}
