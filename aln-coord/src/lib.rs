//! [`Span`]: a half-open interval on a sequence, with the convention in the
//! type instead of in a doc comment.
//!
//! # Why this crate exists
//!
//! Four conventions are in use across the Dfam and GIRI lineages:
//!
//! | convention | who uses it | `chr1` bases 101..200 |
//! |---|---|---|
//! | 0-based half-open | `rmblast-lib`, `aln_core::Alignment`, this crate | `100, 200` |
//! | 1-based fully closed | Smitten IDs, RepeatMasker `.out`, BLAST tabular, Stockholm | `101, 200` |
//! | 0-based fully closed | parasail's `end_query`/`end_ref` accessors | `100, 199` |
//! | 1-based half-open | nobody | n/a |
//!
//! Stored as a bare `pub u64` with the convention in a comment, the difference
//! between the first two is a silent off-by-one that survives every test whose
//! fixtures are self-consistent. `dfam-coord`'s `validate_sequences` catches
//! those after the fact: it fetches the genome and checks whether the sequence
//! matches. The `_halfopen` repair it applies is this bug, in Dfam's production
//! data.
//!
//! You cannot read a `Span` without naming the convention you want, so the
//! compiler rejects the mismatch instead.
//!
//! # The convention
//!
//! **0-based, half-open, forward-strand, `start <= end`.** An empty span is
//! legal and means zero bases. Absence is `Option<Span>`.
//!
//! Strand is absent on purpose. A `Span` is always the forward-strand span,
//! matching `aln_core::Alignment`, which keeps `strand` in its own field. Carry
//! strand here and a reverse-strand interval ends up stored as `start > end`,
//! which is what BLAST does on minus-strand hits and what the rest of the stack
//! forbids.
//!
//! # Interop
//!
//! Construct with the named constructor for whatever convention the input uses
//! ([`Span::from_1b_closed`]), read back with the named accessor for whatever
//! the output needs ([`Span::as_1b_closed`]). Nothing else should add or
//! subtract 1 from a coordinate.
//!
//! # Scope
//!
//! The surface is the minimum a real conversion needs. A method added later is
//! a compatible change; removing one breaks callers, and the repair has to move
//! dfam-lib, RepeatAfterMe and dfam-curator together. So a method arrives here
//! when a caller needs it, not before.
//!
//! Two things are missing on purpose. There is no 0-based-fully-closed
//! constructor: the only sources in that convention are parasail's result
//! accessors, which `aln-parasail` never calls, and RepeatAfterMe's `glocal`
//! and `library` internals, which have not moved onto `Span`. There are also no
//! interval operations (`contains`, `intersect`, `overlaps`), because the MSA
//! column code has not moved onto `Span` either.

use std::fmt;
use std::ops::Range;

// Genomes are memory-mapped whole here, so 32-bit targets were never supported.
// Stating it makes `range_usize` a total function rather than a lossy cast.
const _: () = assert!(
    std::mem::size_of::<usize>() >= std::mem::size_of::<u64>(),
    "aln-coord requires a 64-bit target"
);

// ── Errors ────────────────────────────────────────────────────────────────────

/// Why a pair of coordinates could not become a [`Span`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CoordError {
    /// `end` precedes `start`. A reverse-strand interval is still stored
    /// ascending; flip the strand, not the coordinates.
    Inverted { start: u64, end: u64 },

    /// A 1-based coordinate of 0. Usually the `(0, 0)` sentinel some parsers
    /// use for "no coordinates", which needs an `Option<Span>` instead.
    ZeroInOneBased { start: u64, end: u64 },
}

impl fmt::Display for CoordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CoordError::Inverted { start, end } => {
                write!(f, "span {start}..{end} is inverted: end precedes start")
            }
            CoordError::ZeroInOneBased { start, end } => write!(
                f,
                "1-based coordinates {start}-{end} include 0; 1-based positions begin at 1"
            ),
        }
    }
}

impl std::error::Error for CoordError {}

pub type Result<T> = std::result::Result<T, CoordError>;

// ── Span ──────────────────────────────────────────────────────────────────────

/// A 0-based, half-open, forward-strand interval: `[start, end)`.
///
/// Fields are private. Every way in and out names its convention.
///
/// `Default` is not derived on purpose. A default `Span` would be `0..0`, and a
/// zero coordinate standing in for "unknown" is the sentinel this type removes.
/// Use `Option<Span>`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Span {
    start: u64,
    end: u64,
}

impl Span {
    // ── Construction ──────────────────────────────────────────────────────

    /// From 0-based half-open coordinates, the unmarked house convention.
    pub fn new(start: u64, end: u64) -> Result<Self> {
        if end < start {
            return Err(CoordError::Inverted { start, end });
        }
        Ok(Span { start, end })
    }

    /// From 1-based fully-closed coordinates: Smitten identifiers, Stockholm,
    /// RepeatMasker `.out`, BLAST tabular, crossmatch `.align`.
    ///
    /// `(101, 200)` becomes `100..200`. The empty encoding `(n, n - 1)` is
    /// accepted and yields an empty span.
    pub fn from_1b_closed(start: u64, end: u64) -> Result<Self> {
        if start == 0 || end == 0 {
            return Err(CoordError::ZeroInOneBased { start, end });
        }
        if end + 1 < start {
            return Err(CoordError::Inverted { start, end });
        }
        Ok(Span { start: start - 1, end })
    }

    // ── Reading ───────────────────────────────────────────────────────────

    /// 0-based start.
    pub fn start(&self) -> u64 {
        self.start
    }

    /// 0-based end, exclusive.
    pub fn end(&self) -> u64 {
        self.end
    }

    /// 0-based half-open, for `rmblast-lib` and `aln_core::Alignment`.
    pub fn as_0b_half_open(&self) -> (u64, u64) {
        (self.start, self.end)
    }

    /// 1-based fully closed, for Smitten identifiers, `.out`, `.align`,
    /// Stockholm and BLAST tabular.
    ///
    /// `None` for an empty span, which has no representation in a closed
    /// convention. The arithmetic gives `(start + 1, start)`, and nothing
    /// downstream treats a descending pair as an error: Smitten reads
    /// `start > end` as minus strand under its V1 rule, so does the
    /// RepeatModeler fallback in `dfam-stk-io`, and dfam-curator's FASTA and
    /// clustal writers emit exactly that form for reverse-strand rows. An empty
    /// span at 50 written as `chr1:51-50` parses back as two bases at `50..51`
    /// on the opposite strand.
    ///
    /// `(0, 0)` is no better. 0 is unreachable as a 1-based coordinate, which is
    /// what makes it tempting, but Smitten's `convert_id` rejects any range
    /// containing one, so an identifier written that way cannot be read back. It
    /// also drops the position: an empty span at 50 and one at 900 flatten to
    /// the same marker.
    ///
    /// dfam-curator's `seq_label` already omits the coordinate suffix when it
    /// has no coordinates, returning the bare name. `None` maps onto that.
    ///
    /// `aln_core::Alignment::validate` rejects a zero-length alignment, so a
    /// writer working from a validated `Alignment` can say
    /// `.expect("validate rejects zero-length alignments")`.
    pub fn as_1b_closed(&self) -> Option<(u64, u64)> {
        (!self.is_empty()).then(|| (self.start + 1, self.end))
    }

    /// As a slice index: `&seq[span.range_usize()]`.
    ///
    /// The reason the house convention is half-open. Nothing adds or subtracts
    /// on this path, so a slice taken through it cannot be off by one.
    pub fn range_usize(&self) -> Range<usize> {
        self.start as usize..self.end as usize
    }

    // ── Measuring ─────────────────────────────────────────────────────────

    /// Length in bases.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Zero bases. Absent is a different thing, spelled `Option<Span>`.
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }
}

impl fmt::Display for Span {
    /// 0-based half-open, matching Rust range syntax so the convention is
    /// visible in every log line and assertion failure.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}..{}", self.start, self.end)
    }
}

// ── Calibration harness ───────────────────────────────────────────────────────

/// Table-driven checks for the boundaries where a coordinate changes convention.
///
/// Every parser and writer that crosses a convention gets one of these, with
/// literal numbers worked out by hand. `RepeatAfterMe`'s
/// `coordinates_match_stk2ranges_convention` is the pattern to copy. A test in
/// that shape fails when a parser's convention drifts from what its callers
/// assume, which no fixture-versus-fixture test can do.
///
/// ```
/// use aln_coord::{Span, calibration::{Case, check_all}};
///
/// // A Smitten identifier `chr1:101-200_+` names 100 bases.
/// check_all(&[Case {
///     label: "smitten chr1:101-200_+",
///     span: Span::from_1b_closed(101, 200).unwrap(),
///     expect_0b_half_open: (100, 200),
///     expect_1b_closed: Some((101, 200)),
///     expect_len: 100,
/// }]);
/// ```
pub mod calibration {
    use super::Span;

    /// One boundary, stated in every convention it has to survive.
    pub struct Case {
        /// Where the coordinate came from, quoted literally where possible.
        pub label: &'static str,
        pub span: Span,
        pub expect_0b_half_open: (u64, u64),
        pub expect_1b_closed: Option<(u64, u64)>,
        pub expect_len: u64,
    }

    /// Check every case, panicking with the label of the first that disagrees.
    pub fn check_all(cases: &[Case]) {
        for c in cases {
            assert_eq!(
                c.span.as_0b_half_open(),
                c.expect_0b_half_open,
                "{}: 0-based half-open",
                c.label
            );
            assert_eq!(
                c.span.as_1b_closed(),
                c.expect_1b_closed,
                "{}: 1-based fully closed",
                c.label
            );
            assert_eq!(c.span.len(), c.expect_len, "{}: length", c.label);

            // A span that survives a round trip through the 1-based form is one
            // an output file can carry without drift.
            if let Some((s, e)) = c.expect_1b_closed {
                assert_eq!(
                    Span::from_1b_closed(s, e).unwrap(),
                    c.span,
                    "{}: 1-based round trip",
                    c.label
                );
            }
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_based_closed_start_loses_one_and_end_does_not() {
        let s = Span::from_1b_closed(101, 200).unwrap();
        assert_eq!(s.as_0b_half_open(), (100, 200));
        assert_eq!(s.len(), 100);
        assert_eq!(s.as_1b_closed(), Some((101, 200)));
    }

    #[test]
    fn a_single_base_is_the_same_base_in_every_convention() {
        let s = Span::from_1b_closed(1, 1).unwrap();
        assert_eq!(s.as_0b_half_open(), (0, 1));
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn zero_is_rejected_as_a_one_based_coordinate() {
        // The `(0, 0)` sentinel in dfam-curator's clustal writer and
        // dfam-stk-io's name parser lands here rather than becoming base 0.
        assert_eq!(
            Span::from_1b_closed(0, 0),
            Err(CoordError::ZeroInOneBased { start: 0, end: 0 })
        );
    }

    #[test]
    fn zero_is_a_legal_zero_based_start() {
        assert_eq!(Span::new(0, 10).unwrap().len(), 10);
    }

    #[test]
    fn inverted_coordinates_are_rejected_rather_than_swapped() {
        // A minus-strand hit still stores its span ascending. Swapping here
        // would hide the caller's strand bug.
        assert_eq!(
            Span::new(200, 100),
            Err(CoordError::Inverted { start: 200, end: 100 })
        );
    }

    #[test]
    fn empty_spans_have_no_closed_form() {
        let e = Span::new(50, 50).unwrap();
        assert!(e.is_empty());
        assert_eq!(e.len(), 0);
        assert_eq!(e.as_1b_closed(), None);
    }

    #[test]
    fn the_closed_empty_encoding_would_read_back_as_a_minus_strand_range() {
        // `(start + 1, start)` is what the arithmetic gives for an empty span.
        // Smitten's V1 rule and the RepeatModeler fallback both read a
        // descending pair as minus strand, so downstream parsers accept it.
        let empty_at_50 = Span::new(50, 50).unwrap();
        assert_eq!(empty_at_50.as_1b_closed(), None);

        let what_a_v1_parser_sees = Span::from_1b_closed(50, 51).unwrap();
        assert_eq!(what_a_v1_parser_sees.len(), 2);
        assert_ne!(what_a_v1_parser_sees, empty_at_50);
    }

    #[test]
    fn the_one_based_empty_encoding_round_trips_to_an_empty_span() {
        let e = Span::from_1b_closed(51, 50).unwrap();
        assert!(e.is_empty());
        assert_eq!(e.as_0b_half_open(), (50, 50));
    }

    #[test]
    fn slicing_needs_no_arithmetic() {
        let seq = b"ACGTACGTAC";
        let s = Span::from_1b_closed(3, 5).unwrap();
        assert_eq!(&seq[s.range_usize()], b"GTA");
    }

    #[test]
    fn adjacent_spans_touch_without_overlapping() {
        // The property that closed coordinates cost you: `a.end == b.start`.
        let a = Span::new(0, 100).unwrap();
        let b = Span::new(100, 200).unwrap();
        assert_eq!(a.end(), b.start());
        assert_eq!(a.len() + b.len(), 200);
    }

    #[test]
    fn display_shows_the_convention() {
        assert_eq!(Span::from_1b_closed(101, 200).unwrap().to_string(), "100..200");
    }
}
