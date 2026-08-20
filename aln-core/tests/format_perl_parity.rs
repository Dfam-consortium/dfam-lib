//! Byte-for-byte parity with RepeatMasker's `SearchResult.pm` writers.
//!
//! The fixtures under `tests/data/` were produced by the Perl itself — see
//! `scratch/cmdump.pl` in the commit that added them, reproduced below — and are
//! compared verbatim, trailing whitespace included. Column alignment is the
//! whole point of these formats, so anything less than an exact match is a bug.
//!
//! ```perl
//! use lib "/usr/local/RepeatMasker"; use SearchResult;
//! my $r = SearchResult->new(
//!     queryName => "chr1", subjName => "AluY",
//!     queryStart => 1, queryEnd => 20, queryRemaining => 0,
//!     subjStart  => 1, subjEnd  => 20, subjRemaining  => 0,
//!     orientation => "", score => 500,
//!     pctDiverge => 5.0, pctDelete => 0.0, pctInsert => 0.0,
//!     queryString => $q, subjString => $s);
//! $r->setId(1);
//! print $r->toStringFormatted(SearchResult::N_AlignWithQuerySeq);
//! print $r->toStringFormatted(SearchResult::OutFileFormat);
//! ```

use aln_core::fmt::{to_crossmatch, to_out_line, AlignmentMode};
use aln_core::{Alignment, SearchResult, Strand};

/// Build a result whose fields mirror the Perl object's, from a gapped pair.
#[allow(clippy::too_many_arguments)]
fn result(
    q_name: &str,
    s_name: &str,
    gapped_q: &[u8],
    gapped_s: &[u8],
    strand: Strand,
    score: i32,
    q_left: usize,
    s_left: usize,
    pct: (f64, f64, f64),
    id: Option<u32>,
    matrix: Option<&str>,
) -> (SearchResult, Vec<u8>, Vec<u8>) {
    let query: Vec<u8> = gapped_q.iter().copied().filter(|&b| b != b'-').collect();
    let subj_fwd: Vec<u8> = gapped_s.iter().copied().filter(|&b| b != b'-').collect();
    // `gapped()` reverse-complements the subject on the minus strand, so the
    // stored forward subject must be the reverse complement of the display row.
    let subject = if strand == Strand::Minus {
        aln_core::seq::revcomp(&subj_fwd)
    } else {
        subj_fwd.clone()
    };

    let mut a = Alignment::from_gapped(
        q_name, s_name, 0, 0, strand, score, gapped_q, gapped_s,
    )
    .unwrap();
    a.query_len = Some(query.len() + q_left);
    a.subj_len = Some(subject.len() + s_left);

    let mut r = SearchResult::new(a);
    r.pct_diverge = pct.0;
    r.pct_delete = pct.1;
    r.pct_insert = pct.2;
    r.id = id;
    r.matrix_name = matrix.map(|s| s.to_string());
    (r, query, subject)
}

fn assert_matches(expected: &str, actual: &str, label: &str) {
    if expected != actual {
        // Show the first differing line with visible trailing whitespace.
        let e: Vec<&str> = expected.lines().collect();
        let a: Vec<&str> = actual.lines().collect();
        for i in 0..e.len().max(a.len()) {
            let (el, al) = (e.get(i).copied().unwrap_or("<none>"), a.get(i).copied().unwrap_or("<none>"));
            if el != al {
                panic!(
                    "{label}: line {} differs\n  perl: {:?}\n  rust: {:?}\n\n\
                     full expected:\n{}\nfull actual:\n{}",
                    i + 1,
                    el,
                    al,
                    expected,
                    actual
                );
            }
        }
        panic!("{label}: outputs differ only in trailing newlines");
    }
}

#[test]
fn crossmatch_plus_strand() {
    let (r, q, s) = result(
        "chr1", "AluY",
        b"ACGTACGTACGTACGTACGT",
        b"ACGTACGAACGTACGTACGT",
        Strand::Plus, 500, 0, 0, (5.0, 0.0, 0.0), Some(1), None,
    );
    let out = to_crossmatch(&r, &q, &s, AlignmentMode::WithQuerySeq).unwrap();
    assert_matches(include_str!("data/cm_plus.txt"), &out, "cm_plus");
}

#[test]
fn crossmatch_wraps_at_fifty_columns() {
    let seq: Vec<u8> = b"ACGT".iter().cycle().take(120).copied().collect();
    let (r, q, s) = result(
        "q", "s", &seq, &seq, Strand::Plus, 1000, 0, 0, (0.0, 0.0, 0.0), None, None,
    );
    let out = to_crossmatch(&r, &q, &s, AlignmentMode::WithQuerySeq).unwrap();
    assert_matches(include_str!("data/cm_wrap.txt"), &out, "cm_wrap");
}

#[test]
fn crossmatch_with_gaps() {
    let (r, q, s) = result(
        "q", "s",
        b"ACGT--ACGTAAGT",
        b"ACGTTTACGTACGT",
        Strand::Plus, 200, 5, 3, (8.3, 0.0, 14.3), None, Some("14p35g"),
    );
    let out = to_crossmatch(&r, &q, &s, AlignmentMode::WithQuerySeq).unwrap();
    assert_matches(include_str!("data/cm_gaps.txt"), &out, "cm_gaps");
}

#[test]
fn crossmatch_minus_strand_with_query_seq() {
    let (r, q, s) = result(
        "q", "s",
        b"ACGTACGTACGTACGT",
        b"ACGTACGTACGTACGT",
        Strand::Minus, 400, 4, 9, (0.0, 0.0, 0.0), None, None,
    );
    let out = to_crossmatch(&r, &q, &s, AlignmentMode::WithQuerySeq).unwrap();
    assert_matches(include_str!("data/cm_minus_query.txt"), &out, "cm_minus_query");
}

#[test]
fn crossmatch_minus_strand_with_subj_seq() {
    let (r, q, s) = result(
        "q", "s",
        b"ACGTACGTACGTACGT",
        b"ACGTACGTACGTACGT",
        Strand::Minus, 400, 4, 9, (0.0, 0.0, 0.0), None, None,
    );
    let out = to_crossmatch(&r, &q, &s, AlignmentMode::WithSubjSeq).unwrap();
    assert_matches(include_str!("data/cm_minus_subj.txt"), &out, "cm_minus_subj");
}

#[test]
fn out_line_plus_strand() {
    let (r, _, _) = result(
        "chr1", "AluY",
        b"ACGTACGTACGTACGTACGT",
        b"ACGTACGAACGTACGTACGT",
        Strand::Plus, 500, 0, 0, (5.0, 0.0, 0.0), Some(1), None,
    );
    assert_matches(include_str!("data/out_plus.txt"), &to_out_line(&r), "out_plus");
}

#[test]
fn out_line_minus_strand() {
    let (r, _, _) = result(
        "q", "s",
        b"ACGTACGTACGTACGT",
        b"ACGTACGTACGTACGT",
        Strand::Minus, 400, 4, 9, (0.0, 0.0, 0.0), None, None,
    );
    assert_matches(include_str!("data/out_minus.txt"), &to_out_line(&r), "out_minus");
}
