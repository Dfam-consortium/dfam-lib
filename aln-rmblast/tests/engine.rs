//! Conventions this backend has to get right, pinned against real searches.
//!
//! The conversion from rmblast's `Hsp` to [`aln_core::Alignment`] is a field
//! copy, which is only safe because both sides use 0-based half-open offsets
//! with plus-strand subject coordinates. These tests prove that rather than
//! assume it — in particular that rmblast's BLASTNA sentinel bytes do not leak
//! into the reported offsets.

use aln_core::stats::{rescore, RescoreParams};
use aln_core::{Sequence, Strand, SubstMatrix};
use aln_engine::engine::{SearchEngine, SearchParams, SeqSource};
use aln_rmblast::{RmblastEngine, RmblastOptions};

const CROSSMATCH_14P35G: &str = "\
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
    SubstMatrix::parse(CROSSMATCH_14P35G).unwrap()
}

fn engine(min_score: i32) -> RmblastEngine {
    let params = SearchParams {
        matrix: Some(matrix()),
        gap_init: -25,
        ins_gap_ext: -5,
        del_gap_ext: -5,
        min_match: 7,
        min_score,
        mask_level: 101, // effectively disabled unless a test asks for it
        ..Default::default()
    };
    RmblastEngine::new(params, RmblastOptions::default()).unwrap()
}

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
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

fn random_seq(rng: &mut Rng, len: usize) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.below(4)]).collect()
}

fn mem(seqs: Vec<Sequence>) -> SeqSource {
    SeqSource::Memory(seqs)
}

/// A self-hit: the same 400 bp as query and subject.
#[test]
fn offsets_exclude_sentinels() {
    let mut rng = Rng::new(1);
    let s = random_seq(&mut rng, 400);
    let e = engine(100);

    let hits = e
        .search(
            &mem(vec![Sequence::new("q", s.clone())]),
            &mem(vec![Sequence::new("s", s.clone())]),
        )
        .unwrap();

    let full = hits
        .iter()
        .find(|a| a.strand == Strand::Plus)
        .expect("expected a plus-strand self-hit");

    // If the BLASTNA sentinel bytes leaked into the offsets these would be 1
    // and 401 rather than 0 and 400.
    assert_eq!(full.query_start, 0, "query start");
    assert_eq!(full.query_end, 400, "query end");
    assert_eq!(full.subj_start, 0, "subject start");
    assert_eq!(full.subj_end, 400, "subject end");
    assert_eq!(full.query_len, Some(400));
    assert_eq!(full.subj_len, Some(400));
    assert_eq!(full.edits.to_cigar(), "400M");
    full.validate().unwrap();
}

/// A hit embedded partway into a longer subject must land at the true offset.
///
/// Ignored in debug builds by an upstream bug — see
/// `examples/left_flank_panic.rs`. `rmblast-lib`'s `search/gapped.rs:184`
/// computes `b.as_ptr().add(n - 1 - first_b_index)` unconditionally in the
/// left-extension (`REVERSE`) pass, before the loop whose bound would make it
/// safe. When `first_b_index >= n` the `usize` subtraction underflows: it panics
/// under `debug_assertions` and wraps to an out-of-range pointer in release.
/// The pointer is never dereferenced, so release results are correct — hence
/// this runs, and passes, with `--release`.
#[cfg_attr(
    debug_assertions,
    ignore = "rmblast-lib gapped.rs:184 underflows on left extension in debug builds"
)]
#[test]
fn embedded_hit_lands_at_the_right_offset() {
    let mut rng = Rng::new(2);
    let left = random_seq(&mut rng, 250);
    let core = random_seq(&mut rng, 400);
    let right = random_seq(&mut rng, 180);

    let mut subject = left.clone();
    subject.extend_from_slice(&core);
    subject.extend_from_slice(&right);

    let e = engine(100);
    let hits = e
        .search(
            &mem(vec![Sequence::new("q", core.clone())]),
            &mem(vec![Sequence::new("s", subject.clone())]),
        )
        .unwrap();

    let best = hits.iter().max_by_key(|a| a.score).expect("no hit");
    assert_eq!(best.subj_start, left.len(), "subject offset");
    assert_eq!(best.subj_end, left.len() + core.len());
    assert_eq!(best.query_start, 0);
    assert_eq!(best.query_end, core.len());
}

/// On the minus strand the subject span must stay plus-strand and ascending,
/// and `Alignment::gapped` must reconstruct a sensible alignment from it.
///
/// Same upstream constraint as `embedded_hit_lands_at_the_right_offset`: the
/// subject carries a left flank, so a left extension is required.
#[cfg_attr(
    debug_assertions,
    ignore = "rmblast-lib gapped.rs:184 underflows on left extension in debug builds"
)]
#[test]
fn minus_strand_keeps_plus_strand_subject_coordinates() {
    let mut rng = Rng::new(3);
    let core = random_seq(&mut rng, 400);
    let rc: Vec<u8> = aln_core::seq::revcomp(&core);

    let mut subject = random_seq(&mut rng, 120);
    let offset = subject.len();
    subject.extend_from_slice(&rc);
    subject.extend_from_slice(&random_seq(&mut rng, 90));

    let e = engine(100);
    let hits = e
        .search(
            &mem(vec![Sequence::new("q", core.clone())]),
            &mem(vec![Sequence::new("s", subject.clone())]),
        )
        .unwrap();

    let best = hits.iter().max_by_key(|a| a.score).expect("no hit");
    assert_eq!(best.strand, Strand::Minus, "expected a minus-strand hit");
    assert!(best.subj_start < best.subj_end, "subject span must ascend");
    assert_eq!(best.subj_start, offset, "subject offset");
    assert_eq!(best.subj_end, offset + core.len());

    // The reconstructed pair must be a perfect match: gapped() reverse-
    // complements the subject slice for minus-strand hits.
    let (gq, gs) = best.gapped(&core, &subject).unwrap();
    assert_eq!(gq, gs, "minus-strand reconstruction should be identity here");
}

/// rmblast's reported score must survive re-scoring under `aln-core`'s model.
///
/// This is the end-to-end check that the matrix transpose, the gap penalties and
/// the edit-script mapping are all mutually consistent — a slip in any one of
/// them shows up here.
#[test]
fn reported_scores_survive_rescoring() {
    let mut rng = Rng::new(4);
    let base = random_seq(&mut rng, 500);

    // Diverge the query so the alignment carries real substitutions and indels.
    let mut query = Vec::with_capacity(base.len());
    for &b in &base {
        match rng.below(100) {
            0..=3 => {}                                  // deletion
            4..=6 => {
                query.push(b"ACGT"[rng.below(4)]);
                query.push(b);
            }                                             // insertion
            7..=16 => query.push(b"ACGT"[rng.below(4)]), // substitution
            _ => query.push(b),
        }
    }

    let e = engine(100);
    let hits = e
        .search(
            &mem(vec![Sequence::new("q", query.clone())]),
            &mem(vec![Sequence::new("s", base.clone())]),
        )
        .unwrap();
    assert!(!hits.is_empty(), "expected at least one hit");

    let m = matrix();
    let rp = RescoreParams {
        gap_open: -25,
        ins_gap_extend: -5,
        del_gap_extend: -5,
        ..RescoreParams::new(&m)
    };

    for a in &hits {
        let (gq, gs) = a.gapped(&query, &base).unwrap();
        let r = rescore(&gq, &gs, &rp).unwrap();
        assert_eq!(
            r.score, a.score,
            "rmblast reported {} but the traceback rescores to {} (cigar {})",
            a.score,
            r.score,
            a.edits.to_cigar()
        );
    }
}

/// Searching several subjects must attribute each hit to the right one.
#[test]
fn multiple_subjects_are_reported_separately() {
    let mut rng = Rng::new(5);
    let a = random_seq(&mut rng, 400);
    let b = random_seq(&mut rng, 400);

    let e = engine(100);
    let hits = e
        .search(
            &mem(vec![Sequence::new("q", a.clone())]),
            &mem(vec![
                Sequence::new("decoy", b.clone()),
                Sequence::new("target", a.clone()),
            ]),
        )
        .unwrap();

    let best = hits.iter().max_by_key(|x| x.score).expect("no hit");
    assert_eq!(best.subj_name, "target");
    assert_eq!(best.query_name, "q");
}

/// The statistics rmblast computes during the search are available without
/// recomputing them.
#[test]
fn stats_come_back_with_the_alignment() {
    let mut rng = Rng::new(6);
    let s = random_seq(&mut rng, 400);
    let e = engine(100);

    let pairs = e
        .search_with_stats(
            &mem(vec![Sequence::new("q", s.clone())]),
            &mem(vec![Sequence::new("s", s.clone())]),
        )
        .unwrap();

    let (aln, stats) = pairs
        .iter()
        .max_by_key(|(a, _)| a.score)
        .expect("no hit");
    assert_eq!(aln.edits.to_cigar(), "400M");
    assert_eq!(stats.mismatches, 0, "a self-hit has no mismatches");
    assert_eq!(stats.matches, 400);
    assert!(stats.kdiv.abs() < 1e-9, "kdiv should be 0 for a self-hit");
}

/// Gap penalties are stored signed in `SearchParams` and passed to rmblast as
/// positive magnitudes.
#[test]
fn gap_penalties_reach_rmblast_as_magnitudes() {
    let e = engine(150);
    // crossmatch gap_init -25 / gap_ext -5 becomes NCBI's 20 / 5: NCBI charges
    // open + k*extend where crossmatch charges open + (k-1)*extend.  This is the
    // same 20/5 pair RepeatMasker and dfam-curator pass on the command line.
    assert_eq!(e.rmblast_params().gap_open, 20);
    assert_eq!(e.rmblast_params().gap_extend, 5);
    assert_eq!(e.rmblast_params().min_raw_gapped_score, 150);
    // X-drops default to dfam-curator's derivation from the minimum score.
    assert_eq!(e.rmblast_params().xdrop_ungap, 300);
    assert_eq!(e.rmblast_params().xdrop_gap, 75);
    assert_eq!(e.rmblast_params().xdrop_gap_final, 150);
}

#[test]
fn a_missing_matrix_is_rejected() {
    let params = SearchParams { matrix: None, ..Default::default() };
    assert!(RmblastEngine::new(params, RmblastOptions::default()).is_err());
}

#[test]
fn asymmetric_gap_extension_is_rejected() {
    let params = SearchParams {
        matrix: Some(matrix()),
        ins_gap_ext: -5,
        del_gap_ext: -9,
        ..Default::default()
    };
    let err = match RmblastEngine::new(params, RmblastOptions::default()) {
        Err(e) => e,
        Ok(_) => panic!("asymmetric gap extension should have been rejected"),
    };
    assert!(err.to_string().contains("single gap-extension"), "{err}");
}

#[test]
fn a_prepared_blast_database_is_rejected_with_an_explanation() {
    let e = engine(150);
    let err = e
        .search(
            &mem(vec![Sequence::new("q", b"ACGT".to_vec())]),
            &SeqSource::BlastDb(std::path::PathBuf::from("/tmp/nonexistent")),
        )
        .unwrap_err();
    assert!(err.to_string().contains("external rmblastn binary"), "{err}");
    assert!(!e.accepts(&SeqSource::BlastDb(std::path::PathBuf::from("/tmp/x"))));
}
