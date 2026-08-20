//! CAF parity with `SearchResult.pm::_toCAF`, and a record of what
//! `_toCIGAR` actually emits.
//!
//! Fixtures in `tests/data/caf.txt` and `tests/data/cigar.txt` were produced by
//! the Perl, one `label<TAB>record` per line.
//!
//! CAF is reproduced byte for byte. The CIGAR *payload* deliberately is not —
//! `_toCIGAR` pairs each run's length with the following run's operator and
//! opens with a zero-length run, so ten aligned columns come out as `0M10M`.
//! The broken strings are asserted here so the divergence is documented and
//! visible rather than merely described in a comment; the shared record header
//! is still checked byte for byte.

use std::collections::HashMap;

use aln_core::fmt::{to_caf, to_cigar_record};
use aln_core::{Alignment, SearchResult, Strand};

/// `(label, gapped query, gapped subject, minus?, subject name, class, id)`
const CASES: &[(&str, &str, &str, bool, &str, &str, u32)] = &[
    ("simple", "ACGTACGTAC", "ACGTACGAAC", false, "AluY", "SINE/Alu", 1),
    ("gapq", "ACGT--GTAC", "ACGTTTGTAC", false, "AluY", "SINE/Alu", 2),
    ("gaps", "ACGTTTGTAC", "ACGT--GTAC", false, "AluY", "SINE/Alu", 3),
    ("both", "AAGACTT---A", "AAT--CTAATA", false, "L1", "LINE/L1", 4),
    ("minus", "ACGTACGTAC", "ACGTACGTAC", true, "AluY", "SINE/Alu", 5),
    ("hashname", "ACGTACGTAC", "ACGTACGTAC", false, "AluY#SINE/Alu", "", 6),
    ("trailgap", "ACGTACGT--", "ACGTACGTAC", false, "AluY", "SINE/Alu", 7),
];

fn fixtures(text: &str) -> HashMap<&str, &str> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let (label, rec) = l.split_once('\t').expect("fixture line needs a tab");
            (label, rec)
        })
        .collect()
}

fn build(case: &(&str, &str, &str, bool, &str, &str, u32)) -> (SearchResult, Vec<u8>, Vec<u8>) {
    let (_, gq, gs, minus, s_name, class, id) = *case;
    let strand = if minus { Strand::Minus } else { Strand::Plus };

    let query: Vec<u8> = gq.bytes().filter(|&b| b != b'-').collect();
    let subj_fwd: Vec<u8> = gs.bytes().filter(|&b| b != b'-').collect();
    let subject = if minus {
        aln_core::seq::revcomp(&subj_fwd)
    } else {
        subj_fwd.clone()
    };

    let mut a = Alignment::from_gapped(
        "chr1", s_name, 0, 0, strand, 500, gq.as_bytes(), gs.as_bytes(),
    )
    .unwrap();
    // The Perl driver used queryRemaining => 7, subjRemaining => 9.
    a.query_len = Some(query.len() + 7);
    a.subj_len = Some(subject.len() + 9);

    let mut r = SearchResult::new(a);
    r.pct_diverge = 5.0;
    r.pct_delete = 1.5;
    r.pct_insert = 2.5;
    r.id = Some(id);
    r.subj_class = if class.is_empty() {
        None
    } else {
        Some(class.to_string())
    };
    (r, query, subject)
}

#[test]
fn caf_matches_repeatmasker_byte_for_byte() {
    let expected = fixtures(include_str!("data/caf.txt"));
    for case in CASES {
        let (r, q, s) = build(case);
        let got = to_caf(&r, &q, &s).unwrap();
        assert_eq!(
            got, expected[case.0],
            "\ncase {}\n  perl: {}\n  rust: {}",
            case.0, expected[case.0], got
        );
    }
}

/// The record header is shared with CAF and must match exactly; only the
/// payload after the last comma differs.
#[test]
fn cigar_record_header_matches_repeatmasker() {
    let expected = fixtures(include_str!("data/cigar.txt"));
    for case in CASES {
        let (r, q, s) = build(case);
        let got = to_cigar_record(&r, &q, &s).unwrap();

        let header = |rec: &str| {
            let idx = rec.rfind(',').expect("no payload separator");
            rec[..=idx].to_string()
        };
        assert_eq!(
            header(&got),
            header(expected[case.0]),
            "case {}: record header differs",
            case.0
        );
    }
}

/// Pin the divergence: what the Perl emits, and what we emit instead.
///
/// If `_toCIGAR` is ever fixed upstream these expectations will drift apart and
/// this test will say so.
#[test]
fn cigar_payload_diverges_from_the_broken_perl_output() {
    let expected = fixtures(include_str!("data/cigar.txt"));
    let payload = |rec: &str| rec[rec.rfind(',').unwrap() + 1..].to_string();

    // What the Perl actually produces — note the leading zero-length run and
    // the length/operator mismatch.
    assert_eq!(payload(expected["simple"]), "0M10M");
    assert_eq!(payload(expected["both"]), "0M3D2M2I3M1M");
    assert_eq!(payload(expected["trailgap"]), "0M8I2I");

    // What we produce: the runs RepeatMasker's own doc comment documents.
    let check = |label: &str, want: &str| {
        let case = CASES.iter().find(|c| c.0 == label).unwrap();
        let (r, q, s) = build(case);
        let got = to_cigar_record(&r, &q, &s).unwrap();
        assert_eq!(payload(&got), want, "case {label}");
    };
    check("simple", "10M");
    // The exact example from _toCIGAR's doc comment.
    check("both", "3M2D2M3I1M");
    check("trailgap", "8M2I");
}

/// RepeatMasker's CIGAR calls a query gap `I`; SAM (and `EditScript::to_cigar`)
/// call it `D`. The two must stay mirror images of each other.
#[test]
fn cigar_orientation_is_inverted_relative_to_sam() {
    let case = CASES.iter().find(|c| c.0 == "gapq").unwrap();
    let (r, q, s) = build(case);

    let rm = to_cigar_record(&r, &q, &s).unwrap();
    let rm_payload = &rm[rm.rfind(',').unwrap() + 1..];
    let sam = r.alignment.edits.to_cigar();

    // Query gap run of 2: RepeatMasker says I, SAM says D.
    assert!(rm_payload.contains("2I"), "RepeatMasker payload: {rm_payload}");
    assert!(sam.contains("2D"), "SAM cigar: {sam}");

    let swapped: String = rm_payload
        .chars()
        .map(|c| match c {
            'I' => 'D',
            'D' => 'I',
            other => other,
        })
        .collect();
    assert_eq!(swapped, sam, "the two encodings should be mirror images");
}
