//! The pairwise alignment type and its run-length-encoded edit script.
//!
//! # Coordinate convention
//!
//! **All coordinates in this crate are 0-based, half-open, and expressed in the
//! forward orientation of their own source sequence.**  `start < end` always,
//! on both sides, regardless of strand.  This matches `rmblast-lib`'s `Hsp` and
//! BLAST's `sstart <= send`.
//!
//! The file formats around it (Smitten identifiers, Stockholm, RepeatMasker
//! `.out`, BLAST tabular) are 1-based and fully closed.  Conversion happens
//! *only* at I/O boundaries, through [`aln_coord::Span`]'s named constructors
//! and accessors. See [`Alignment::query_one_based`] and
//! [`Alignment::subject_one_based`].
//!
//! # Strand
//!
//! [`Alignment::strand`] is the strand of the **subject** relative to the query.
//! The query is always read forward.  When the strand is
//! [`Strand::Minus`] the alignment walks the query
//! left-to-right and the subject right-to-left, so the subject's forward-
//! coordinate span is still `subj_start .. subj_end` but its *aligned* bases are
//! the reverse complement of that span.
//!
//! # Why an edit script rather than gapped strings
//!
//! `autocons` aligns thousands of sequences all-against-all.  Materialising two
//! gapped strings per pairwise alignment (as GIRI's `PairwiseAlignment` does)
//! costs `O(n^2)` strings live at once.  The edit script is a handful of bytes
//! and expands to gapped strings on demand via [`Alignment::gapped`].

use crate::error::{Error, Result};
use crate::seq::{self, Strand, GAP};

// ── Edit script ───────────────────────────────────────────────────────────────

/// A single traceback operation.
///
/// Naming follows NCBI's `eGapAlign` conventions, matching `rmblast-lib` so
/// scripts can move between the crates without translation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EditOp {
    /// Aligned pair — match or mismatch.  Both sequences advance.
    Sub,
    /// Gap in the query; the subject advances alone.
    /// Alignment column reads `query='-'`, `subject=base`.
    ///
    /// With the subject-is-consensus convention this is a **deletion** in the
    /// genomic sequence relative to the consensus.
    GapInQuery,
    /// Gap in the subject; the query advances alone.
    /// Alignment column reads `query=base`, `subject='-'`.
    ///
    /// With the subject-is-consensus convention this is an **insertion** in the
    /// genomic sequence relative to the consensus.
    GapInSubject,
}

/// A run-length-encoded traceback.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct EditScript {
    pub ops: Vec<(EditOp, u32)>,
}

impl EditScript {
    pub fn new() -> Self {
        Self::default()
    }

    /// Append `count` of `op`, coalescing with the previous run.  A zero count
    /// is a no-op.
    pub fn push(&mut self, op: EditOp, count: u32) {
        if count == 0 {
            return;
        }
        if let Some(last) = self.ops.last_mut() {
            if last.0 == op {
                last.1 += count;
                return;
            }
        }
        self.ops.push((op, count));
    }

    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }

    /// Total number of alignment columns.
    pub fn align_len(&self) -> u32 {
        self.ops.iter().map(|&(_, n)| n).sum()
    }

    /// Query bases consumed (`Sub` + `GapInSubject`).
    pub fn query_consumed(&self) -> u32 {
        self.ops
            .iter()
            .filter(|&&(op, _)| op != EditOp::GapInQuery)
            .map(|&(_, n)| n)
            .sum()
    }

    /// Subject bases consumed (`Sub` + `GapInQuery`).
    pub fn subject_consumed(&self) -> u32 {
        self.ops
            .iter()
            .filter(|&&(op, _)| op != EditOp::GapInSubject)
            .map(|&(_, n)| n)
            .sum()
    }

    /// Reverse the run order — used to flip a left-extension script, and to
    /// re-orient a script built against a reverse-complemented subject.
    pub fn reverse(&mut self) {
        self.ops.reverse();
    }

    /// Iterate individual columns (un-run-length-encoded).
    pub fn iter_columns(&self) -> impl Iterator<Item = EditOp> + '_ {
        self.ops
            .iter()
            .flat_map(|&(op, n)| std::iter::repeat_n(op, n as usize))
    }

    /// Render as a CIGAR string using the SAM `M`/`I`/`D` opcodes, where `I` is
    /// an insertion **in the query** relative to the subject.
    pub fn to_cigar(&self) -> String {
        let mut s = String::new();
        for &(op, n) in &self.ops {
            let c = match op {
                EditOp::Sub => 'M',
                EditOp::GapInSubject => 'I',
                EditOp::GapInQuery => 'D',
            };
            s.push_str(&n.to_string());
            s.push(c);
        }
        s
    }

    /// Parse a SAM-style CIGAR string (`M`/`=`/`X` all map to [`EditOp::Sub`]).
    pub fn from_cigar(cigar: &str) -> Result<Self> {
        let mut out = EditScript::new();
        let mut n: u32 = 0;
        let mut saw_digit = false;
        for c in cigar.chars() {
            if let Some(d) = c.to_digit(10) {
                n = n
                    .checked_mul(10)
                    .and_then(|v| v.checked_add(d))
                    .ok_or_else(|| Error::Alignment(format!("CIGAR run too long: {cigar}")))?;
                saw_digit = true;
                continue;
            }
            if !saw_digit {
                return Err(Error::Alignment(format!(
                    "CIGAR operator {c:?} without a run length in {cigar:?}"
                )));
            }
            let op = match c {
                'M' | '=' | 'X' => EditOp::Sub,
                'I' => EditOp::GapInSubject,
                'D' => EditOp::GapInQuery,
                _ => {
                    return Err(Error::Alignment(format!(
                        "unsupported CIGAR operator {c:?} in {cigar:?}"
                    )))
                }
            };
            out.push(op, n);
            n = 0;
            saw_digit = false;
        }
        if saw_digit {
            return Err(Error::Alignment(format!(
                "CIGAR ends with a run length and no operator: {cigar:?}"
            )));
        }
        Ok(out)
    }

    /// Derive a script from a pair of equal-length gapped strings.
    ///
    /// Columns where both sides are gaps are dropped.  Padding (`' '`) is
    /// rejected — trim to the aligned region first.
    pub fn from_gapped(query: &[u8], subject: &[u8]) -> Result<Self> {
        if query.len() != subject.len() {
            return Err(Error::Alignment(format!(
                "gapped strings differ in length: {} vs {}",
                query.len(),
                subject.len()
            )));
        }
        let mut out = EditScript::new();
        for (i, (&q, &s)) in query.iter().zip(subject).enumerate() {
            if seq::is_pad(q) || seq::is_pad(s) {
                return Err(Error::Alignment(format!(
                    "padding character at column {i}; trim to the aligned region first"
                )));
            }
            let (qg, sg) = (seq::is_gap(q), seq::is_gap(s));
            match (qg, sg) {
                (false, false) => out.push(EditOp::Sub, 1),
                (true, false) => out.push(EditOp::GapInQuery, 1),
                (false, true) => out.push(EditOp::GapInSubject, 1),
                (true, true) => {} // an all-gap column carries no information
            }
        }
        Ok(out)
    }
}

// ── Alignment ─────────────────────────────────────────────────────────────────

/// One pairwise alignment: coordinates, strand, score, and traceback.
///
/// Sequence *bytes* are deliberately not stored — see [`Alignment::gapped`].
/// This keeps the type cheap enough to hold `O(n^2)` of them during an
/// all-against-all pass.
#[derive(Debug, Clone, PartialEq)]
pub struct Alignment {
    /// Query = genomic / derived sequence, under the RepeatMasker convention.
    pub query_name: String,
    /// Subject = consensus / ancestral sequence.
    pub subj_name: String,

    /// 0-based half-open span in the forward query.
    pub query_start: usize,
    pub query_end: usize,
    /// Full length of the query sequence, when known.
    pub query_len: Option<usize>,

    /// 0-based half-open span in the forward subject.
    pub subj_start: usize,
    pub subj_end: usize,
    /// Full length of the subject sequence, when known.
    pub subj_len: Option<usize>,

    /// Strand of the subject relative to the (always-forward) query.
    pub strand: Strand,

    /// Raw alignment score under the aligner's own scoring system.
    pub score: i32,

    /// Traceback, in display order (left to right).
    pub edits: EditScript,
}

impl Alignment {
    /// Construct with defaults for the optional fields.
    pub fn new(
        query_name: impl Into<String>,
        subj_name: impl Into<String>,
        query_start: usize,
        subj_start: usize,
        strand: Strand,
        score: i32,
        edits: EditScript,
    ) -> Self {
        let query_end = query_start + edits.query_consumed() as usize;
        let subj_end = subj_start + edits.subject_consumed() as usize;
        Alignment {
            query_name: query_name.into(),
            subj_name: subj_name.into(),
            query_start,
            query_end,
            query_len: None,
            subj_start,
            subj_end,
            subj_len: None,
            strand,
            score,
            edits,
        }
    }

    /// Aligned query length in bases (excludes gaps).
    pub fn query_span(&self) -> usize {
        self.query_end - self.query_start
    }

    /// Aligned subject length in bases (excludes gaps).
    pub fn subj_span(&self) -> usize {
        self.subj_end - self.subj_start
    }

    /// Number of alignment columns.
    pub fn align_len(&self) -> usize {
        self.edits.align_len() as usize
    }

    /// Query span as 1-based fully-closed coordinates, for `.out` / CAF /
    /// Stockholm / `dfam-curator` interop.
    pub fn query_one_based(&self) -> (u64, u64) {
        (self.query_start as u64 + 1, self.query_end as u64)
    }

    /// Subject span as 1-based fully-closed coordinates.
    ///
    /// Still forward-oriented and ascending even on the minus strand — that is
    /// BLAST's convention.  RepeatMasker's `.out` writer flips these itself for
    /// `C` orientation; do the flip there, not here.
    pub fn subject_one_based(&self) -> (u64, u64) {
        (self.subj_start as u64 + 1, self.subj_end as u64)
    }

    /// Bases of the query remaining beyond the alignment's end, if the query
    /// length is known.  RepeatMasker prints this as the parenthesised value.
    pub fn query_remaining(&self) -> Option<usize> {
        self.query_len.map(|l| l.saturating_sub(self.query_end))
    }

    /// Bases of the subject remaining beyond the alignment's end.
    pub fn subj_remaining(&self) -> Option<usize> {
        self.subj_len.map(|l| l.saturating_sub(self.subj_end))
    }

    /// Check that the edit script agrees with the recorded coordinates, that
    /// neither end runs past a known sequence length, and that both sides align
    /// at least one base.
    pub fn validate(&self) -> Result<()> {
        let q = self.edits.query_consumed() as usize;
        let s = self.edits.subject_consumed() as usize;
        if self.query_start + q != self.query_end {
            return Err(Error::Alignment(format!(
                "{}: edit script consumes {q} query bases but span is {}..{}",
                self.query_name, self.query_start, self.query_end
            )));
        }
        if self.subj_start + s != self.subj_end {
            return Err(Error::Alignment(format!(
                "{}: edit script consumes {s} subject bases but span is {}..{}",
                self.subj_name, self.subj_start, self.subj_end
            )));
        }
        if let Some(l) = self.query_len {
            if self.query_end > l {
                return Err(Error::Alignment(format!(
                    "{}: query end {} exceeds length {l}",
                    self.query_name, self.query_end
                )));
            }
        }
        if let Some(l) = self.subj_len {
            if self.subj_end > l {
                return Err(Error::Alignment(format!(
                    "{}: subject end {} exceeds length {l}",
                    self.subj_name, self.subj_end
                )));
            }
        }
        // A side that consumes no bases is not an alignment: an empty edit
        // script, or one made entirely of gaps against the other sequence.  No
        // aligner in this workspace emitted either when the check went in, and
        // the suite passed unchanged, so reaching this is a construction bug.
        //
        // It also lets the 1-based writers drop their empty case.  A closed
        // convention cannot express an empty range, so `Span::as_1b_closed`
        // returns `None`; a writer working from a validated `Alignment` can
        // `.expect()` that away instead of inventing a representation.
        if self.query_start == self.query_end || self.subj_start == self.subj_end {
            return Err(Error::Alignment(format!(
                "{} vs {}: zero-length alignment (query {}..{}, subject {}..{})",
                self.query_name,
                self.subj_name,
                self.query_start,
                self.query_end,
                self.subj_start,
                self.subj_end
            )));
        }
        Ok(())
    }

    /// Expand to a pair of equal-length gapped strings, `(query, subject)`.
    ///
    /// `query` and `subject` are the **full forward** source sequences; the
    /// spans are sliced out here.  On the minus strand the subject slice is
    /// reverse-complemented so the returned pair reads left-to-right in query
    /// order — the same presentation as crossmatch and BLAST.
    pub fn gapped(&self, query: &[u8], subject: &[u8]) -> Result<(Vec<u8>, Vec<u8>)> {
        if self.query_end > query.len() {
            return Err(Error::Alignment(format!(
                "query {:?} has {} bases but the alignment ends at {}",
                self.query_name,
                query.len(),
                self.query_end
            )));
        }
        if self.subj_end > subject.len() {
            return Err(Error::Alignment(format!(
                "subject {:?} has {} bases but the alignment ends at {}",
                self.subj_name,
                subject.len(),
                self.subj_end
            )));
        }

        let q_src = &query[self.query_start..self.query_end];
        let s_slice = &subject[self.subj_start..self.subj_end];
        let s_src: Vec<u8> = if self.strand.is_minus() {
            seq::revcomp(s_slice)
        } else {
            s_slice.to_vec()
        };

        let n = self.align_len();
        let mut q_out = Vec::with_capacity(n);
        let mut s_out = Vec::with_capacity(n);
        let (mut qi, mut si) = (0usize, 0usize);

        for (op, count) in self.edits.ops.iter().copied() {
            for _ in 0..count {
                match op {
                    EditOp::Sub => {
                        q_out.push(q_src[qi]);
                        s_out.push(s_src[si]);
                        qi += 1;
                        si += 1;
                    }
                    EditOp::GapInQuery => {
                        q_out.push(GAP);
                        s_out.push(s_src[si]);
                        si += 1;
                    }
                    EditOp::GapInSubject => {
                        q_out.push(q_src[qi]);
                        s_out.push(GAP);
                        qi += 1;
                    }
                }
            }
        }
        Ok((q_out, s_out))
    }

    /// Build from a pair of gapped strings plus coordinates — the inverse of
    /// [`gapped`](Self::gapped), for parsing crossmatch / BLAST / Stockholm.
    ///
    /// Every argument is an independent field of the record being parsed, so
    /// bundling them into a struct would only move the same list elsewhere.
    #[allow(clippy::too_many_arguments)]
    pub fn from_gapped(
        query_name: impl Into<String>,
        subj_name: impl Into<String>,
        query_start: usize,
        subj_start: usize,
        strand: Strand,
        score: i32,
        gapped_query: &[u8],
        gapped_subject: &[u8],
    ) -> Result<Self> {
        let edits = EditScript::from_gapped(gapped_query, gapped_subject)?;
        Ok(Alignment::new(
            query_name,
            subj_name,
            query_start,
            subj_start,
            strand,
            score,
            edits,
        ))
    }

    /// Counts of matches, mismatches and gap columns over the alignment.
    ///
    /// A "match" requires the two bytes to be equal after upper-casing; this is
    /// the simple identity notion, not the matrix-positive notion GIRI's
    /// `PairwiseAlignmentStats::iPCnt` uses.  For matrix-aware statistics use
    /// [`crate::stats::rescore`].
    pub fn identity_counts(&self, query: &[u8], subject: &[u8]) -> Result<IdentityCounts> {
        let (q, s) = self.gapped(query, subject)?;
        let mut c = IdentityCounts::default();
        for (&qb, &sb) in q.iter().zip(&s) {
            if seq::is_gap(qb) || seq::is_gap(sb) {
                c.gap_columns += 1;
            } else if qb.eq_ignore_ascii_case(&sb) {
                c.matches += 1;
            } else {
                c.mismatches += 1;
            }
        }
        Ok(c)
    }
}

/// Plain identity tallies over an alignment.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IdentityCounts {
    pub matches: u32,
    pub mismatches: u32,
    pub gap_columns: u32,
}

impl IdentityCounts {
    /// Fraction identical over aligned (non-gap) columns; `None` if there are none.
    pub fn identity(&self) -> Option<f64> {
        let n = self.matches + self.mismatches;
        (n > 0).then(|| self.matches as f64 / n as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn script(spec: &[(EditOp, u32)]) -> EditScript {
        let mut s = EditScript::new();
        for &(op, n) in spec {
            s.push(op, n);
        }
        s
    }

    #[test]
    fn push_coalesces_runs() {
        let s = script(&[(EditOp::Sub, 3), (EditOp::Sub, 2), (EditOp::GapInQuery, 1)]);
        assert_eq!(s.ops, vec![(EditOp::Sub, 5), (EditOp::GapInQuery, 1)]);
        assert_eq!(s.align_len(), 6);
    }

    #[test]
    fn consumption_accounting() {
        // 4 aligned, 2 deleted from query, 3 inserted into query.
        let s = script(&[
            (EditOp::Sub, 4),
            (EditOp::GapInQuery, 2),
            (EditOp::GapInSubject, 3),
        ]);
        assert_eq!(s.query_consumed(), 7); // Sub + GapInSubject
        assert_eq!(s.subject_consumed(), 6); // Sub + GapInQuery
        assert_eq!(s.align_len(), 9);
    }

    #[test]
    fn cigar_round_trips() {
        let s = script(&[
            (EditOp::Sub, 10),
            (EditOp::GapInSubject, 2),
            (EditOp::Sub, 5),
            (EditOp::GapInQuery, 1),
        ]);
        assert_eq!(s.to_cigar(), "10M2I5M1D");
        assert_eq!(EditScript::from_cigar("10M2I5M1D").unwrap(), s);
    }

    #[test]
    fn cigar_rejects_malformed_input() {
        assert!(EditScript::from_cigar("M").is_err());
        assert!(EditScript::from_cigar("10").is_err());
        assert!(EditScript::from_cigar("10Z").is_err());
    }

    #[test]
    fn from_gapped_drops_all_gap_columns() {
        //             both gap here ↓
        let q = b"ACGT--AC";
        let s = b"AC-T--AC";
        let e = EditScript::from_gapped(q, s).unwrap();
        // Columns: M M I M (skip) (skip) M M — the trailing runs coalesce.
        assert_eq!(e.to_cigar(), "2M1I3M");
        assert_eq!(e.align_len(), 6);
    }

    #[test]
    fn from_gapped_rejects_padding() {
        assert!(EditScript::from_gapped(b" ACGT", b"AACGT").is_err());
    }

    #[test]
    fn gapped_expands_a_forward_alignment() {
        // query   [2..7] = GTTAC
        // subject [3..9] = GTACGT
        // script: 2 aligned, 1 base deleted from the query, 3 aligned.
        let query = b"ACGTTACG";
        let subject = b"NNNGTACGT";
        let e = script(&[
            (EditOp::Sub, 2),
            (EditOp::GapInQuery, 1),
            (EditOp::Sub, 3),
        ]);
        let mut a = Alignment::new("q", "s", 2, 3, Strand::Plus, 100, e);
        a.query_len = Some(query.len());
        a.subj_len = Some(subject.len());
        a.validate().unwrap();
        assert_eq!((a.query_start, a.query_end), (2, 7));
        assert_eq!((a.subj_start, a.subj_end), (3, 9));

        let (gq, gs) = a.gapped(query, subject).unwrap();
        assert_eq!(gq.len(), gs.len());
        assert_eq!(&gq, b"GT-TAC");
        assert_eq!(&gs, b"GTACGT");
    }

    #[test]
    fn gapped_reverse_complements_the_subject_on_minus_strand() {
        // Subject forward span is CGT; on the minus strand the aligned bases
        // are its reverse complement, ACG.
        let query = b"ACG";
        let subject = b"NNCGTNN";
        let e = script(&[(EditOp::Sub, 3)]);
        let a = Alignment::new("q", "s", 0, 2, Strand::Minus, 30, e);
        a.validate().unwrap();

        let (gq, gs) = a.gapped(query, subject).unwrap();
        assert_eq!(&gq, b"ACG");
        assert_eq!(&gs, b"ACG"); // revcomp(CGT) == ACG — a perfect minus-strand match
    }

    #[test]
    fn coordinates_stay_forward_and_ascending_on_minus_strand() {
        let e = script(&[(EditOp::Sub, 3)]);
        let a = Alignment::new("q", "s", 0, 2, Strand::Minus, 30, e);
        assert!(a.subj_start < a.subj_end);
        assert_eq!(a.subject_one_based(), (3, 5));
        assert_eq!(a.query_one_based(), (1, 3));
    }

    #[test]
    fn validate_rejects_an_empty_edit_script() {
        // A construction bug, not something an aligner emits.
        let a = Alignment::new("q", "s", 0, 0, Strand::Plus, 0, EditScript::new());
        assert!(a.validate().is_err());
    }

    #[test]
    fn validate_rejects_an_alignment_that_is_all_gap_on_one_side() {
        // A non-empty script can still consume nothing on one side, which leaves
        // that span empty and unwritable in any closed convention.
        let e = script(&[(EditOp::GapInQuery, 3)]);
        let a = Alignment::new("q", "s", 0, 0, Strand::Plus, 0, e);
        assert_eq!(a.query_span(), 0);
        assert!(a.validate().is_err());
    }

    #[test]
    fn validate_accepts_a_single_aligned_base() {
        // The smallest thing that is still an alignment.
        let e = script(&[(EditOp::Sub, 1)]);
        let a = Alignment::new("q", "s", 0, 0, Strand::Plus, 1, e);
        a.validate().unwrap();
        assert_eq!(a.query_one_based(), (1, 1));
    }

    #[test]
    fn validate_catches_coordinate_drift() {
        let e = script(&[(EditOp::Sub, 3)]);
        let mut a = Alignment::new("q", "s", 0, 0, Strand::Plus, 0, e);
        a.query_end = 99; // corrupt it
        assert!(a.validate().is_err());
    }

    #[test]
    fn gapped_reports_a_too_short_sequence_rather_than_panicking() {
        let e = script(&[(EditOp::Sub, 10)]);
        let a = Alignment::new("q", "s", 0, 0, Strand::Plus, 0, e);
        let err = a.gapped(b"ACGT", b"ACGT").unwrap_err();
        assert!(matches!(err, Error::Alignment(_)));
    }

    #[test]
    fn round_trip_gapped_to_script_and_back() {
        let query = b"ACGTTACG";
        let subject = b"ACGATACG";
        let a = Alignment::from_gapped(
            "q", "s", 0, 0, Strand::Plus, 42,
            b"ACGTTACG", b"ACGATACG",
        )
        .unwrap();
        let (gq, gs) = a.gapped(query, subject).unwrap();
        assert_eq!(&gq, b"ACGTTACG");
        assert_eq!(&gs, b"ACGATACG");
    }

    #[test]
    fn identity_counts_ignore_case_and_tally_gaps() {
        let query = b"acgtt";
        let subject = b"ACGAT";
        let a = Alignment::from_gapped("q", "s", 0, 0, Strand::Plus, 0, b"acgtt", b"ACGAT").unwrap();
        let c = a.identity_counts(query, subject).unwrap();
        assert_eq!(c.matches, 4);
        assert_eq!(c.mismatches, 1);
        assert_eq!(c.gap_columns, 0);
        assert_eq!(c.identity(), Some(0.8));
    }

    #[test]
    fn remaining_bases_need_a_known_length() {
        let e = script(&[(EditOp::Sub, 3)]);
        let mut a = Alignment::new("q", "s", 0, 0, Strand::Plus, 0, e);
        assert_eq!(a.query_remaining(), None);
        a.query_len = Some(10);
        assert_eq!(a.query_remaining(), Some(7));
    }
}
