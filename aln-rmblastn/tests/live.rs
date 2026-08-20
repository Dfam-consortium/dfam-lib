//! End-to-end check against a real `rmblastn`.
//!
//! Skipped when the binary is absent, so the suite still passes on machines
//! without RepeatMasker installed. Unit tests cover the parser; this covers
//! everything the parser cannot — flag construction, the makeblastdb step, and
//! whether the output we ask for is the output we get.

use aln_core::Sequence;
use aln_engine::engine::{SearchEngine, SearchParams, SeqSource};
use aln_rmblastn::{RmblastnEngine, RmblastnOptions};

/// RepeatMasker ships NCBI-format matrices; rmblastn will not accept the
/// crossmatch layout that `SubstMatrix` holds, so the live test points at a
/// real one and skips if it is not installed.
/// The rmblast install to shell out to. A machine can carry several BLAST+
/// versions; this one ships `rmblastn` and its matching `makeblastdb`.
const RMBLAST_BIN: &str = "/usr/local/rmblast/bin/rmblastn";

const NCBI_MATRIX: &str =
    "/usr/local/RepeatMasker-4.2.0-Dfam_3.9_RB/Matrices/ncbi/nt/20p43g.matrix";

fn have(exe: &str) -> bool {
    std::process::Command::new(exe)
        .arg("-version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

fn engine() -> RmblastnEngine {
    // gap_init -25 / ext -5 is the crossmatch pair RepeatMasker uses; the engine
    // converts it to NCBI's -gapopen 20 -gapextend 5.
    let p = SearchParams {
        gap_init: -25,
        ins_gap_ext: -5,
        del_gap_ext: -5,
        min_match: 7,
        min_score: 30,
        cores: Some(1),
        path_to_engine: Some(RMBLAST_BIN.into()),
        ..Default::default()
    };
    let opts = RmblastnOptions {
        matrix_path: Some(NCBI_MATRIX.into()),
        ..Default::default()
    };
    RmblastnEngine::new(p, opts).unwrap()
}

/// Two copies of one element, one carrying substitutions, must align.
#[test]
fn finds_a_planted_homology() {
    if !have(RMBLAST_BIN) || !std::path::Path::new(NCBI_MATRIX).exists() {
        eprintln!("skipping: {RMBLAST_BIN} or the NCBI matrix is unavailable");
        return;
    }
    let core = "ACGTTGCAAGGCTTACGGATCCGTTACAGGCATTACGGATCA".repeat(6);
    let mut mutated: Vec<u8> = core.clone().into_bytes();
    for i in (17..mutated.len()).step_by(23) {
        mutated[i] = if mutated[i] == b'A' { b'G' } else { b'A' };
    }

    let q = vec![Sequence::new("query", mutated)];
    let s = vec![Sequence::new("subject", core.into_bytes())];
    let hits = engine()
        .search(&SeqSource::Memory(q), &SeqSource::Memory(s))
        .expect("search should succeed");

    assert!(!hits.is_empty(), "planted homology should be found");
    let best = hits.iter().max_by_key(|a| a.score).unwrap();
    assert_eq!(best.query_name, "query");
    assert_eq!(best.subj_name, "subject");
    assert!(best.score > 100, "score was {}", best.score);
    // The traceback must survive the round trip, not just the coordinates.
    assert!(
        best.edits.align_len() > 100,
        "expected a long alignment, got {} columns",
        best.edits.align_len()
    );
    assert!(best.query_start < best.query_end);
    assert!(best.subj_start < best.subj_end);
}

/// `makeblastdb` must come from the same install as `rmblastn`, not from
/// whatever `PATH` happens to find first — mixing BLAST+ versions across the
/// database build and the search corrupts results silently.
#[test]
fn makeblastdb_is_taken_from_the_rmblastn_install() {
    if !have(RMBLAST_BIN) {
        eprintln!("skipping: {RMBLAST_BIN} unavailable");
        return;
    }
    let p = SearchParams { path_to_engine: Some(RMBLAST_BIN.into()), ..Default::default() };
    let e = RmblastnEngine::new(p, RmblastnOptions::default()).unwrap();
    let v = e.version();
    let (rmb, mkb) = v.split_once(" / ").expect("version reports both binaries");
    let ver = |s: &str| s.split_whitespace().last().unwrap_or("").to_string();
    assert_eq!(
        ver(rmb),
        ver(mkb),
        "rmblastn and makeblastdb versions differ: {v}"
    );
}

/// Handing the engine a parsed crossmatch matrix without a file must fail with
/// an explanation, not silently produce wrong scores.
#[test]
fn refuses_to_synthesise_an_ncbi_matrix() {
    let p = SearchParams {
        matrix: Some(aln_core::SubstMatrix::parse(
            "GAP -10 -2\n    A   C\n  A 3 -4\n  C -4 3\n",
        ).unwrap()),
        ..Default::default()
    };
    let err = match RmblastnEngine::new(p, RmblastnOptions::default()) {
        Err(e) => e.to_string(),
        Ok(_) => panic!("a matrix with no file must be rejected at construction"),
    };
    assert!(err.contains("transpose"), "unhelpful message: {err}");
}

/// A 2bit source is the one thing this engine cannot take; it must say so
/// rather than fail obscurely inside the child process.
#[test]
fn rejects_a_twobit_source() {
    let e = engine();
    assert!(!e.accepts(&SeqSource::TwoBit("x.2bit".into())));
    let err = e
        .search(
            &SeqSource::TwoBit("x.2bit".into()),
            &SeqSource::Memory(vec![]),
        )
        .unwrap_err()
        .to_string();
    assert!(err.contains("twoBitToFa"), "unhelpful message: {err}");
}
