//! A plain `O(mn)` affine-gap aligner — the correctness arbiter for the SIMD
//! backends.
//!
//! This is Gotoh's affine-gap recurrence written for legibility, not speed.  It
//! exists so that parasail, Farrar's striped SSE2 code and anything else added
//! later can be differentially tested against an implementation whose
//! recurrences can be read off the page.
//!
//! # What agreement means
//!
//! Compare **scores** exactly.  Do *not* require identical tracebacks: when two
//! paths tie, which one a backend reports depends on its tie-breaking order, and
//! striped implementations legitimately differ from a row-major scalar one.  A
//! traceback is correct if re-scoring it reproduces the reported score, which is
//! what [`aln_core::stats::rescore`] is for.
//!
//! # Recurrences
//!
//! Rows index the query, columns the subject.
//!
//! ```text
//!   E[i][j] = max(H[i-1][j] - open, E[i-1][j] - extend)   query consumed, gap in subject
//!   F[i][j] = max(H[i][j-1] - open, F[i][j-1] - extend)   subject consumed, gap in query
//!   H[i][j] = max( H[i-1][j-1] + s(subject[j-1], query[i-1]),
//!                  E[i][j], F[i][j]
//!                  [, 0 when local] )
//! ```
//!
//! `open` covers the gap's first position, so a length-`k` gap costs
//! `open + (k-1) * extend` — the same convention as
//! [`aln_core::stats::rescore`] and crossmatch.
//!
//! ## Why both gap states derive from `H`, not from a separate match state
//!
//! A three-state formulation that opens gaps only from the *match* state
//! forbids an insertion directly abutting a deletion.  That is a real variant,
//! but it is not the one parasail, Farrar, crossmatch or BLAST implement, and it
//! scores a minority of alignments 1–2 points lower.  Since this module's whole
//! job is to arbitrate the others, it follows the mainstream recurrence: `E` and
//! `F` are both opened from `H`, which already incorporates all three states.

use aln_core::align::{Alignment, EditOp, EditScript};
use aln_core::{Sequence, Strand, SubstMatrix};
use aln_engine::{AlignMode, AlignParams, AlignerCaps, EngineError, PairwiseAligner, Result};

/// Guard against pathological allocations: three traceback bytes per cell.
const MAX_CELLS: u64 = 1 << 28; // ~256M cells => ~768 MB of traceback

/// Sentinel for "unreachable", far enough from `i32::MIN` to survive one
/// subtraction of a gap penalty without wrapping.
const NEG: i32 = i32::MIN / 4;

/// How `H[i][j]` was reached.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(u8)]
enum HFrom {
    /// Local alignment restarted here.
    Zero,
    /// Diagonal — an aligned pair.
    Diag,
    /// Took the `E` (gap-in-subject) state at this cell.
    E,
    /// Took the `F` (gap-in-query) state at this cell.
    F,
}

/// How a gap state was reached: opened from `H`, or extended from itself.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(u8)]
enum GapFrom {
    Open,
    Extend,
}

/// Traceback state machine position.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum State {
    H,
    E,
    F,
}

/// The reference aligner.
pub struct ReferenceAligner {
    matrix: SubstMatrix,
    params: AlignParams,
    /// Index used for bytes outside the matrix alphabet, when the alphabet has
    /// an `N`.  `None` makes unknown bytes an error instead.
    unknown_idx: Option<usize>,
}

/// A subject encoded into matrix indices.
pub struct EncodedSubject {
    name: String,
    idx: Vec<usize>,
    len: usize,
}

impl ReferenceAligner {
    pub fn new(matrix: SubstMatrix, params: AlignParams) -> Result<Self> {
        params.validate()?;
        let unknown_idx = matrix.index_of(b'N');
        Ok(ReferenceAligner { matrix, params, unknown_idx })
    }

    pub fn matrix(&self) -> &SubstMatrix {
        &self.matrix
    }

    fn encode(&self, seq: &[u8], what: &str) -> Result<Vec<usize>> {
        seq.iter()
            .map(|&b| {
                self.matrix.index_of(b).or(self.unknown_idx).ok_or_else(|| {
                    EngineError::backend(
                        "reference",
                        format!(
                            "{what} contains {:?}, which is outside the matrix alphabet {:?} \
                             and the alphabet has no 'N' to fall back on",
                            b as char,
                            String::from_utf8_lossy(self.matrix.alphabet())
                        ),
                    )
                })
            })
            .collect()
    }

    /// `(free_query_ends, free_subject_ends)` for the configured mode.
    fn free_ends(&self) -> (bool, bool) {
        match self.params.mode {
            AlignMode::Local => (true, true),
            AlignMode::Global => (false, false),
            AlignMode::SemiGlobal { free_query_ends, free_subject_ends } => {
                (free_query_ends, free_subject_ends)
            }
        }
    }
}

impl PairwiseAligner for ReferenceAligner {
    type Profile = EncodedSubject;

    fn name(&self) -> &'static str {
        "reference"
    }

    fn caps(&self) -> AlignerCaps {
        AlignerCaps {
            name: "reference",
            modes: &[
                AlignMode::Local,
                AlignMode::Global,
                AlignMode::SemiGlobal { free_query_ends: true, free_subject_ends: true },
                AlignMode::SemiGlobal { free_query_ends: true, free_subject_ends: false },
                AlignMode::SemiGlobal { free_query_ends: false, free_subject_ends: true },
                AlignMode::SemiGlobal { free_query_ends: false, free_subject_ends: false },
            ],
            traceback: true,
            banded: false,
            simd: "scalar",
        }
    }

    fn params(&self) -> &AlignParams {
        &self.params
    }

    fn prepare_subject(&self, subject: &Sequence) -> Result<EncodedSubject> {
        Ok(EncodedSubject {
            name: subject.name.clone(),
            idx: self.encode(&subject.seq, "subject")?,
            len: subject.len(),
        })
    }

    fn align_prepared(
        &self,
        subject: &EncodedSubject,
        query: &Sequence,
    ) -> Result<Option<Alignment>> {
        if let Some(why) = self.caps().supports(&self.params) {
            return Err(EngineError::unsupported("reference", why));
        }
        let q_idx = self.encode(&query.seq, "query")?;
        let (m, n) = (q_idx.len(), subject.idx.len());
        if m == 0 || n == 0 {
            return Ok(None);
        }
        let cells = (m as u64 + 1) * (n as u64 + 1);
        if cells > MAX_CELLS {
            return Err(EngineError::backend(
                "reference",
                format!(
                    "{m} x {n} exceeds the {MAX_CELLS}-cell guard; use a banded \
                     or SIMD backend for sequences this size"
                ),
            ));
        }

        let dp = self.fill(&q_idx, &subject.idx);
        let Some((score, end_i, end_j)) = self.best_cell(&dp, m, n) else {
            return Ok(None);
        };
        if score < self.params.min_score {
            return Ok(None);
        }

        let (edits, start_i, start_j) = self.traceback(&dp, end_i, end_j);
        if edits.is_empty() {
            return Ok(None);
        }

        let mut a = Alignment::new(
            query.name.clone(),
            subject.name.clone(),
            start_i,
            start_j,
            Strand::Plus,
            score,
            edits,
        );
        a.query_len = Some(query.len());
        a.subj_len = Some(subject.len);
        Ok(Some(a))
    }
}

/// The filled dynamic-programming tables.
struct Dp {
    n: usize, // subject length; row stride is n + 1
    h: Vec<i32>,
    e: Vec<i32>,
    f: Vec<i32>,
    tb_h: Vec<HFrom>,
    tb_e: Vec<GapFrom>,
    tb_f: Vec<GapFrom>,
    local: bool,
}

impl Dp {
    #[inline]
    fn at(&self, i: usize, j: usize) -> usize {
        i * (self.n + 1) + j
    }
}

impl ReferenceAligner {
    fn fill(&self, q: &[usize], s: &[usize]) -> Dp {
        let (m, n) = (q.len(), s.len());
        let stride = n + 1;
        let size = (m + 1) * stride;
        let local = self.params.mode == AlignMode::Local;
        let (free_q, free_s) = self.free_ends();

        let open = self.params.gap_open as i32;
        let extend = self.params.gap_extend as i32;

        let mut dp = Dp {
            n,
            h: vec![NEG; size],
            e: vec![NEG; size],
            f: vec![NEG; size],
            tb_h: vec![HFrom::Zero; size],
            tb_e: vec![GapFrom::Open; size],
            tb_f: vec![GapFrom::Open; size],
            local,
        };

        dp.h[0] = 0;

        // Column 0: the query hangs off the start of the subject.  Reaching
        // H[i][0] means consuming i query bases against nothing, which is the
        // E state.  Free when the query's ends are free.
        for i in 1..=m {
            let k = dp.at(i, 0);
            let v = if free_q {
                0
            } else if i == 1 {
                -open
            } else {
                dp.e[dp.at(i - 1, 0)] - extend
            };
            dp.e[k] = v;
            dp.h[k] = v;
            dp.tb_e[k] = if i == 1 { GapFrom::Open } else { GapFrom::Extend };
            dp.tb_h[k] = HFrom::E;
        }
        // Row 0: the subject hangs off the start of the query — the F state.
        for j in 1..=n {
            let k = dp.at(0, j);
            let v = if free_s {
                0
            } else if j == 1 {
                -open
            } else {
                dp.f[dp.at(0, j - 1)] - extend
            };
            dp.f[k] = v;
            dp.h[k] = v;
            dp.tb_f[k] = if j == 1 { GapFrom::Open } else { GapFrom::Extend };
            dp.tb_h[k] = HFrom::F;
        }

        for i in 1..=m {
            for j in 1..=n {
                let k = dp.at(i, j);
                let kd = dp.at(i - 1, j - 1);
                let ku = dp.at(i - 1, j);
                let kl = dp.at(i, j - 1);

                // E — gap in the subject; the query advances.
                let open_e = dp.h[ku].saturating_sub(open);
                let ext_e = dp.e[ku].saturating_sub(extend);
                if ext_e > open_e {
                    dp.e[k] = ext_e;
                    dp.tb_e[k] = GapFrom::Extend;
                } else {
                    dp.e[k] = open_e;
                    dp.tb_e[k] = GapFrom::Open;
                }

                // F — gap in the query; the subject advances.
                let open_f = dp.h[kl].saturating_sub(open);
                let ext_f = dp.f[kl].saturating_sub(extend);
                if ext_f > open_f {
                    dp.f[k] = ext_f;
                    dp.tb_f[k] = GapFrom::Extend;
                } else {
                    dp.f[k] = open_f;
                    dp.tb_f[k] = GapFrom::Open;
                }

                // H — best of diagonal, E and F (and 0 when local).
                let sub = self.matrix.score_idx(s[j - 1], q[i - 1]);
                let diag = dp.h[kd].saturating_add(sub);
                let mut best = diag;
                let mut from = HFrom::Diag;
                if dp.e[k] > best {
                    best = dp.e[k];
                    from = HFrom::E;
                }
                if dp.f[k] > best {
                    best = dp.f[k];
                    from = HFrom::F;
                }
                if local && best < 0 {
                    best = 0;
                    from = HFrom::Zero;
                }
                dp.h[k] = best;
                dp.tb_h[k] = from;
            }
        }
        dp
    }

    /// Pick the cell to trace back from, honouring the mode's free ends.
    fn best_cell(&self, dp: &Dp, m: usize, n: usize) -> Option<(i32, usize, usize)> {
        let (free_q, free_s) = self.free_ends();

        if dp.local {
            let mut best = (NEG, 0usize, 0usize);
            for i in 1..=m {
                for j in 1..=n {
                    let v = dp.h[dp.at(i, j)];
                    if v > best.0 {
                        best = (v, i, j);
                    }
                }
            }
            return (best.0 > 0).then_some(best);
        }

        let mut best = (dp.h[dp.at(m, n)], m, n);
        // A free query end lets the alignment stop before the query is used up.
        if free_q {
            for i in 1..=m {
                let v = dp.h[dp.at(i, n)];
                if v > best.0 {
                    best = (v, i, n);
                }
            }
        }
        // A free subject end lets it stop before the subject is used up.
        if free_s {
            for j in 1..=n {
                let v = dp.h[dp.at(m, j)];
                if v > best.0 {
                    best = (v, m, j);
                }
            }
        }
        (best.0 > NEG).then_some(best)
    }

    /// Walk back from `(end_i, end_j)`, returning the edit script in display
    /// order and the 0-based start offsets.
    fn traceback(&self, dp: &Dp, end_i: usize, end_j: usize) -> (EditScript, usize, usize) {
        let (free_q, free_s) = self.free_ends();
        let mut rev: Vec<EditOp> = Vec::new();
        let (mut i, mut j) = (end_i, end_j);
        let mut state = State::H;

        loop {
            if i == 0 && j == 0 {
                break;
            }
            // A free end means the leftover prefix is simply outside the
            // alignment.  Stop rather than walking the zero-initialised edge,
            // which would emit it as a run of leading gaps.
            if i == 0 && free_s {
                break;
            }
            if j == 0 && free_q {
                break;
            }
            let k = dp.at(i, j);
            match state {
                State::H => {
                    if dp.local && dp.h[k] == 0 {
                        break;
                    }
                    match dp.tb_h[k] {
                        HFrom::Zero => break,
                        HFrom::Diag => {
                            if i == 0 || j == 0 {
                                break;
                            }
                            rev.push(EditOp::Sub);
                            i -= 1;
                            j -= 1;
                        }
                        // Switching state consumes no cell; the gap arm below
                        // emits and moves on the next iteration.
                        HFrom::E => state = State::E,
                        HFrom::F => state = State::F,
                    }
                }
                State::E => {
                    if i == 0 {
                        break;
                    }
                    let next = dp.tb_e[k];
                    rev.push(EditOp::GapInSubject);
                    i -= 1;
                    state = match next {
                        GapFrom::Extend => State::E,
                        GapFrom::Open => State::H,
                    };
                }
                State::F => {
                    if j == 0 {
                        break;
                    }
                    let next = dp.tb_f[k];
                    rev.push(EditOp::GapInQuery);
                    j -= 1;
                    state = match next {
                        GapFrom::Extend => State::F,
                        GapFrom::Open => State::H,
                    };
                }
            }
        }

        let mut script = EditScript::new();
        for op in rev.into_iter().rev() {
            script.push(op, 1);
        }
        (script, i, j)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use aln_core::stats::{rescore, RescoreParams};

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

    fn aligner(mode: AlignMode) -> ReferenceAligner {
        let p = AlignParams { mode, gap_open: 25, gap_extend: 5, min_score: 1, ..Default::default() };
        ReferenceAligner::new(matrix(), p).unwrap()
    }

    fn seq(name: &str, s: &[u8]) -> Sequence {
        Sequence::new(name, s.to_vec())
    }

    #[test]
    fn identical_sequences_score_the_diagonal_sum() {
        let a = aligner(AlignMode::Local);
        let q = seq("q", b"ACGT");
        let s = seq("s", b"ACGT");
        let aln = a.align(&q, &s).unwrap().unwrap();
        // A(8) C(12) G(12) T(8)
        assert_eq!(aln.score, 40);
        assert_eq!(aln.edits.to_cigar(), "4M");
        assert_eq!((aln.query_start, aln.query_end), (0, 4));
        assert_eq!((aln.subj_start, aln.subj_end), (0, 4));
    }

    #[test]
    fn local_alignment_trims_flanking_mismatches() {
        let a = aligner(AlignMode::Local);
        // A strong core with junk on both sides of the query.
        let q = seq("q", b"TTTTACGTACGTTTTT");
        let s = seq("s", b"ACGTACGT");
        let aln = a.align(&q, &s).unwrap().unwrap();
        assert_eq!(aln.query_start, 4);
        assert_eq!(aln.query_end, 12);
        assert_eq!((aln.subj_start, aln.subj_end), (0, 8));
        assert_eq!(aln.edits.to_cigar(), "8M");
    }

    #[test]
    fn traceback_rescores_to_the_reported_score() {
        // The property that actually matters when comparing backends.
        let a = aligner(AlignMode::Local);
        let q = seq("q", b"ACGTTTACGTACG");
        let s = seq("s", b"ACGTACGTACG");
        let aln = a.align(&q, &s).unwrap().unwrap();
        aln.validate().unwrap();

        let (gq, gs) = aln.gapped(&q.seq, &s.seq).unwrap();
        let m = matrix();
        let p = RescoreParams {
            gap_open: -25,
            ins_gap_extend: -5,
            del_gap_extend: -5,
            ..RescoreParams::new(&m)
        };
        let r = rescore(&gq, &gs, &p).unwrap();
        assert_eq!(r.score, aln.score, "cigar was {}", aln.edits.to_cigar());
    }

    #[test]
    fn an_affine_gap_costs_open_plus_extends() {
        let a = aligner(AlignMode::Global);
        // Query is missing three subject bases in the middle.
        let q = seq("q", b"ACGTACGT");
        let s = seq("s", b"ACGTAAAACGT");
        let aln = a.align(&q, &s).unwrap().unwrap();
        assert_eq!(aln.edits.subject_consumed(), 11);
        assert_eq!(aln.edits.query_consumed(), 8);

        let (gq, gs) = aln.gapped(&q.seq, &s.seq).unwrap();
        let m = matrix();
        let p = RescoreParams {
            gap_open: -25,
            ins_gap_extend: -5,
            del_gap_extend: -5,
            ..RescoreParams::new(&m)
        };
        assert_eq!(rescore(&gq, &gs, &p).unwrap().score, aln.score);
    }

    #[test]
    fn one_long_gap_beats_several_short_ones() {
        // With open=25 and extend=5 a single 4 bp gap (25+15=40) must be
        // preferred over two 2 bp gaps (2 * 30 = 60).  This is the property an
        // affine implementation has and a linear one does not.
        let a = aligner(AlignMode::Global);
        let q = seq("q", b"ACGTACGTACGT");
        let s = seq("s", b"ACGTACGTTTTTACGT");
        let aln = a.align(&q, &s).unwrap().unwrap();
        let gap_runs = aln
            .edits
            .ops
            .iter()
            .filter(|&&(op, _)| op == EditOp::GapInQuery)
            .count();
        assert_eq!(gap_runs, 1, "cigar was {}", aln.edits.to_cigar());
    }

    #[test]
    fn global_mode_consumes_both_sequences_entirely() {
        let a = aligner(AlignMode::Global);
        let q = seq("q", b"ACGTAC");
        let s = seq("s", b"ACGTGC");
        let aln = a.align(&q, &s).unwrap().unwrap();
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, q.len());
        assert_eq!(aln.subj_start, 0);
        assert_eq!(aln.subj_end, s.len());
    }

    #[test]
    fn semiglobal_free_subject_ends_fits_a_query_inside_a_long_subject() {
        let a = aligner(AlignMode::SemiGlobal {
            free_query_ends: false,
            free_subject_ends: true,
        });
        let q = seq("q", b"ACGTACGT");
        let s = seq("s", b"TTTTTTACGTACGTTTTTTT");
        let aln = a.align(&q, &s).unwrap().unwrap();
        // The whole query is used; the subject's flanks are free.
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, q.len());
        assert_eq!(aln.subj_start, 6);
        assert_eq!(aln.subj_end, 14);
    }

    #[test]
    fn min_score_suppresses_weak_alignments() {
        let p = AlignParams { mode: AlignMode::Local, min_score: 1000, ..Default::default() };
        let a = ReferenceAligner::new(matrix(), p).unwrap();
        let q = seq("q", b"ACGT");
        let s = seq("s", b"ACGT");
        assert!(a.align(&q, &s).unwrap().is_none());
    }

    #[test]
    fn empty_input_yields_no_alignment() {
        let a = aligner(AlignMode::Local);
        assert!(a.align(&seq("q", b""), &seq("s", b"ACGT")).unwrap().is_none());
        assert!(a.align(&seq("q", b"ACGT"), &seq("s", b"")).unwrap().is_none());
    }

    #[test]
    fn unknown_bytes_fall_back_to_the_n_row() {
        let a = aligner(AlignMode::Local);
        // '@' is not in the alphabet; the matrix has N, so it is scored as N.
        let q = seq("q", b"AC@T");
        let s = seq("s", b"ACGT");
        let aln = a.align(&q, &s).unwrap().unwrap();
        // A(8) + C(12) + G/N(-1) + T(8) = 27
        assert_eq!(aln.score, 27);
    }

    #[test]
    fn unknown_bytes_error_when_the_alphabet_has_no_n() {
        let m = SubstMatrix::parse("  A   C\n  1  -1\n -1   1\n").unwrap();
        let a = ReferenceAligner::new(m, AlignParams::default()).unwrap();
        assert!(a.align(&seq("q", b"AG"), &seq("s", b"AC")).is_err());
    }

    #[test]
    fn profile_reuse_gives_the_same_answer_as_one_shot() {
        let a = aligner(AlignMode::Local);
        let s = seq("s", b"ACGTACGTACGT");
        let prepared = a.prepare_subject(&s).unwrap();
        for q in [b"ACGTACGT".as_slice(), b"GTACGTAC".as_slice(), b"TTTT".as_slice()] {
            let query = seq("q", q);
            let via_profile = a.align_prepared(&prepared, &query).unwrap();
            let one_shot = a.align(&query, &s).unwrap();
            assert_eq!(via_profile.map(|x| x.score), one_shot.map(|x| x.score));
        }
    }

    #[test]
    fn alignment_coordinates_always_validate() {
        let a = aligner(AlignMode::Local);
        let pairs: &[(&[u8], &[u8])] = &[
            (b"ACGTACGTAC", b"ACGTACGTAC"),
            (b"ACGTTTTACGT", b"ACGTACGT"),
            (b"AAAACCCCGGGGTTTT", b"ACGTACGTACGT"),
            (b"GGGGGGGG", b"CCCCCCCC"),
        ];
        for (q, s) in pairs {
            let query = seq("q", q);
            let subject = seq("s", s);
            if let Some(aln) = a.align(&query, &subject).unwrap() {
                aln.validate().unwrap();
                let (gq, gs) = aln.gapped(&query.seq, &subject.seq).unwrap();
                assert_eq!(gq.len(), gs.len());
            }
        }
    }
}
