//! Worked examples for the boundaries listed in the migration plan.
//!
//! Copy one into the crate that owns the boundary, replacing the literal `Span`
//! with whatever the parser under test returns. The numbers here are the
//! contract; a parser that disagrees with them is wrong whatever its own
//! fixtures say.

use aln_coord::calibration::{check_all, Case};
use aln_coord::Span;

#[test]
fn smitten_identifiers() {
    // Smitten's grammar annotates every position "1-based ... fully closed",
    // and `convert_id` rejects a 0.  `Identifier::normalize` returns
    // `Range { start, end }` in that convention; dfam-stk-io converts here.
    check_all(&[
        Case {
            label: "chr1:101-200_+",
            span: Span::from_1b_closed(101, 200).unwrap(),
            expect_0b_half_open: (100, 200),
            expect_1b_closed: Some((101, 200)),
            expect_len: 100,
        },
        Case {
            label: "chr1:1-1_+ (first base)",
            span: Span::from_1b_closed(1, 1).unwrap(),
            expect_0b_half_open: (0, 1),
            expect_1b_closed: Some((1, 1)),
            expect_len: 1,
        },
    ]);
}

#[test]
fn blast_tabular_and_crossmatch() {
    // BLAST reports `sstart <= send` on both strands, 1-based closed, which is
    // also what crossmatch `.align` and RepeatMasker `.out` carry.
    check_all(&[Case {
        label: "qstart=1 qend=250",
        span: Span::from_1b_closed(1, 250).unwrap(),
        expect_0b_half_open: (0, 250),
        expect_1b_closed: Some((1, 250)),
        expect_len: 250,
    }]);
}

#[test]
fn parasail_result_accessors() {
    // `parasail_result_get_end_query` and `_end_ref` return the last aligned
    // offset, inclusive.  aln-parasail calls neither: it reads `beg_query` and
    // `beg_ref` from the CIGAR and derives the ends from the consumed counts.
    // There is no `from_0b_closed` constructor for that reason; anyone reaching
    // for those accessors has to add the 1 here, in the open.
    check_all(&[Case {
        label: "beg_ref=0 end_ref=249 (inclusive), +1 to close the span",
        span: Span::new(0, 249 + 1).unwrap(),
        expect_0b_half_open: (0, 250),
        expect_1b_closed: Some((1, 250)),
        expect_len: 250,
    }]);
}

#[test]
fn rmblast_hsp_needs_no_conversion() {
    // `rmblast-lib`'s Hsp is already the house convention; aln-rmblast copies
    // the fields across.  Pinning it here means a change upstream fails loudly.
    check_all(&[Case {
        label: "Hsp { query_start: 100, query_end: 200 }",
        span: Span::new(100, 200).unwrap(),
        expect_0b_half_open: (100, 200),
        expect_1b_closed: Some((101, 200)),
        expect_len: 100,
    }]);
}

/// The `(0, 0)` sentinel that `Option<Span>` replaced.
///
/// Parsers used to return `(0, 0)` for a name without coordinates and writers
/// tested for it.  Under 1-based coordinates 0 is unreachable, so the sentinel
/// worked.  Under 0-based it is the first base of the sequence, which is why
/// absence is now `None` and a 1-based 0 is a construction error.
#[test]
fn the_zero_sentinel_is_unrepresentable() {
    assert!(Span::from_1b_closed(0, 0).is_err());

    // What the sentinel has to become.
    let absent: Option<Span> = None;
    let present_at_first_base = Some(Span::new(0, 100).unwrap());
    assert_ne!(absent, present_at_first_base);
}
