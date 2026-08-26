//! Using rmblast where a [`PairwiseAligner`] is expected.
//!
//! `autocons` is generic over [`PairwiseAligner`] — two sequences in, at most
//! one alignment out. rmblast is a [`SearchEngine`]:
//! it seeds, extends, and returns *many* HSPs per pair, on both strands. This
//! adapter bridges the two by taking the **single best-scoring HSP** for each
//! (subject, query) pair.
//!
//! # Why best-HSP rather than chaining
//!
//! A repeat instance interrupted by an insertion, or one that has picked up an
//! inversion, can come back as several HSPs. Two readings are possible:
//!
//! * **Best HSP** (what this does) — take the highest-scoring segment and drop
//!   the rest. Simple, and matches what a local aligner would have returned:
//!   Smith-Waterman also reports one maximal segment.
//! * **Chain compatible HSPs** — stitch together HSPs that are colinear and
//!   non-overlapping on both sequences, giving a longer composite alignment.
//!   Recovers more of an interrupted element, but needs a chaining policy
//!   (maximum gap, score for the join) that has no obvious right answer.
//!
//! Best-HSP is the conservative choice and keeps behaviour comparable with the
//! `PairwiseAligner` backends. Chaining is left for when there is evidence it is
//! needed.
//!
//! # What this costs
//!
//! rmblast's advantage is searching one query against a *whole database* at
//! once. Driving it one pair at a time throws that away: each call rebuilds a
//! query lookup table and scans a single subject. It is the right shape for
//! correctness and for slotting into `autocons` unchanged, but it is not how to
//! make rmblast fast — that would mean giving `autocons` a search-engine-shaped
//! path where one reference is searched against all instances in a single call.
//!
//! # Seeding can miss what full DP finds
//!
//! rmblast will return nothing for a pair whose similarity never produces a
//! seed. Short or highly diverged instances are the usual casualties. That is
//! inherent to seeded search, not a defect — but it means an `autocons` run on
//! this backend can legitimately place fewer sequences than one on parasail.

use aln_core::{Alignment, Sequence};
use aln_engine::engine::{SearchEngine, SearchParams, SeqSource};
use aln_engine::{AlignMode, AlignParams, AlignerCaps, EngineError, PairwiseAligner, Result};

use crate::{RmblastEngine, RmblastOptions};

const NAME: &str = "rmblast-pairwise";

/// A [`PairwiseAligner`] backed by rmblast, returning the best HSP per pair.
pub struct RmblastPairwise {
    engine: RmblastEngine,
    params: AlignParams,
}

/// The prepared side. rmblast has no reusable profile across subjects, so this
/// just carries the sequence.
pub struct PreparedSubject {
    seq: Sequence,
}

impl RmblastPairwise {
    /// Build from the same [`SearchParams`] the engine takes.
    ///
    /// `align` is used only for [`AlignParams::min_score`] and to reject modes
    /// rmblast cannot honour; the scoring itself comes from `search`.
    pub fn new(
        search: SearchParams,
        opts: RmblastOptions,
        align: AlignParams,
    ) -> Result<Self> {
        if align.mode != AlignMode::Local {
            return Err(EngineError::unsupported(
                NAME,
                format!("rmblast produces local alignments only; got {:?}", align.mode),
            ));
        }
        Ok(RmblastPairwise { engine: RmblastEngine::new(search, opts)?, params: align })
    }

    pub fn engine(&self) -> &RmblastEngine {
        &self.engine
    }
}

impl PairwiseAligner for RmblastPairwise {
    type Profile = PreparedSubject;

    fn name(&self) -> &'static str {
        NAME
    }

    fn caps(&self) -> AlignerCaps {
        AlignerCaps {
            name: NAME,
            modes: &[AlignMode::Local],
            traceback: true,
            banded: false,
            simd: "rmblast",
        }
    }

    fn params(&self) -> &AlignParams {
        &self.params
    }

    fn prepare_subject(&self, subject: &Sequence) -> Result<PreparedSubject> {
        Ok(PreparedSubject { seq: subject.clone() })
    }

    fn align_prepared(
        &self,
        subject: &PreparedSubject,
        query: &Sequence,
    ) -> Result<Option<Alignment>> {
        if query.is_empty() || subject.seq.is_empty() {
            return Ok(None);
        }
        let hits = self.engine.search(
            &SeqSource::Memory(vec![query.clone()]),
            &SeqSource::Memory(vec![subject.seq.clone()]),
        )?;
        // Best HSP wins; ties break towards the longer alignment so the choice
        // is deterministic rather than dependent on rmblast's emission order.
        Ok(hits
            .into_iter()
            .filter(|a| a.score >= self.params.min_score)
            .max_by_key(|a| (a.score, a.align_len())))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use aln_core::SubstMatrix;

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

    fn aligner(min_score: i32) -> RmblastPairwise {
        let search = SearchParams {
            matrix: Some(SubstMatrix::parse(M14P35G).unwrap()),
            gap_init: -25,
            ins_gap_ext: -5,
            del_gap_ext: -5,
            min_match: 7,
            min_score,
            mask_level: 101,
            ..Default::default()
        };
        let align = AlignParams { mode: AlignMode::Local, min_score, ..Default::default() };
        RmblastPairwise::new(search, RmblastOptions::default(), align).unwrap()
    }

    struct Rng(u64);
    impl Rng {
        fn new(s: u64) -> Self {
            Rng(s.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
        }
        fn below(&mut self, n: usize) -> usize {
            let mut x = self.0;
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            self.0 = x;
            (x.wrapping_mul(0x2545_F491_4F6C_DD1D) % n as u64) as usize
        }
    }
    fn rand_seq(r: &mut Rng, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[r.below(4)]).collect()
    }

    #[test]
    fn a_self_hit_comes_back_full_length() {
        let mut r = Rng::new(1);
        let s = rand_seq(&mut r, 400);
        let a = aligner(100);
        let q = Sequence::new("q", s.clone());
        let subj = Sequence::new("s", s);
        let prof = a.prepare_subject(&subj).unwrap();
        let aln = a.align_prepared(&prof, &q).unwrap().expect("expected a hit");
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, 400);
        aln.validate().unwrap();
    }

    #[test]
    fn only_the_best_hsp_is_returned() {
        // Two copies of the same block in the query: rmblast can report both,
        // the adapter must hand back exactly one.
        let mut r = Rng::new(2);
        let block = rand_seq(&mut r, 300);
        let mut query = block.clone();
        query.extend_from_slice(&rand_seq(&mut r, 50));
        query.extend_from_slice(&block);

        let a = aligner(100);
        let subj = Sequence::new("s", block);
        let prof = a.prepare_subject(&subj).unwrap();
        let aln = a
            .align_prepared(&prof, &Sequence::new("q", query))
            .unwrap()
            .expect("expected a hit");
        assert!(aln.query_span() <= 320, "should be one block, got {}", aln.query_span());
    }

    #[test]
    fn unrelated_sequences_produce_nothing() {
        let mut r = Rng::new(3);
        let a = aligner(200);
        let subj = Sequence::new("s", rand_seq(&mut r, 400));
        let prof = a.prepare_subject(&subj).unwrap();
        let q = Sequence::new("q", rand_seq(&mut r, 400));
        assert!(a.align_prepared(&prof, &q).unwrap().is_none());
    }

    #[test]
    fn empty_input_is_not_an_error() {
        let a = aligner(100);
        let prof = a.prepare_subject(&Sequence::new("s", b"ACGT".to_vec())).unwrap();
        assert!(a.align_prepared(&prof, &Sequence::new("q", vec![])).unwrap().is_none());
    }

    #[test]
    fn non_local_modes_are_rejected() {
        let search = SearchParams {
            matrix: Some(SubstMatrix::parse(M14P35G).unwrap()),
            ..Default::default()
        };
        let align = AlignParams { mode: AlignMode::Global, ..Default::default() };
        assert!(RmblastPairwise::new(search, RmblastOptions::default(), align).is_err());
    }
}
