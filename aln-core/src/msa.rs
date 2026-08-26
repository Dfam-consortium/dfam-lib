//! The multiple alignment type, and assembly of one from pairwise alignments.
//!
//! [`MultiAlign`] and [`SequenceRow`] mirror `dfam_curator::alignment` field for
//! field, so migrating that crate onto this one is mechanical.  The one rename
//! is `Orientation` → [`Strand`], which is re-exported here under both names.
//!
//! # The insertion problem
//!
//! Merging N pairwise alignments against a shared reference into one MSA forces
//! a choice about insertions — query bases with no reference position to sit in.
//! The two lineages resolve it differently, and the difference is not cosmetic:
//!
//! * **`dfam-curator`** drops them (`blast.rs::hits_to_multialign`), so the MSA
//!   width always equals the reference length.  Cheap and stable, but insertion
//!   evidence is destroyed.
//! * **GIRI** grows the reference (`MultipleAlignment::adjustReference`),
//!   opening gap columns wide enough for the largest insertion at each point.
//!   Necessary for `acons`/`autocons`, where insertions must survive to be
//!   called into the consensus.
//!
//! [`InsertionPolicy`] makes the choice explicit at the call site rather than
//! burying it in whichever module happens to build the MSA.

use crate::error::{Error, Result};
use crate::seq::{self, Strand, GAP, PAD};
use aln_coord::Span;

/// Alias kept so `dfam-curator` code moves over without edits.
pub use crate::seq::Strand as Orientation;

// ── Rows ──────────────────────────────────────────────────────────────────────

/// One row of a multiple alignment: the reference (index 0) or an instance.
///
/// `seq` is stored gapped and padded to the full alignment width.  Interior
/// indels are `-`; positions the row does not reach are `' '` (see
/// [`crate::seq`] for why the two must stay distinct).
#[derive(Debug, Clone, PartialEq)]
pub struct SequenceRow {
    /// Identifier, e.g. `chr1:12345-12500(+)`.
    pub name: String,

    /// Gapped bytes; length equals the alignment width.
    pub seq: Vec<u8>,

    /// First non-padding column, 0-based.
    pub col_start: usize,
    /// One past the last non-padding column: `col_start..col_end` indexes the
    /// occupied part of `seq`. Both are 0 for a row of nothing but padding.
    pub col_end: usize,

    /// Where the row's residues sit in the original ungapped source sequence,
    /// on the forward strand; `orient` says which strand was aligned. `None`
    /// when the identifier carried no coordinates.
    pub span: Option<Span>,

    pub orient: Strand,

    /// Left-flanking genomic sequence, not part of the alignment.
    pub lf_seq: Option<Vec<u8>>,
    /// Right-flanking genomic sequence.
    pub rf_seq: Option<Vec<u8>>,

    /// GC fraction of the source region (0.0–1.0).
    pub gc_background: Option<f64>,

    /// Raw divergence from the reference.
    pub div: Option<f64>,
    /// Kimura divergence from the reference.
    pub kdiv: Option<f64>,
    /// Transition count against the reference.
    pub trans_i: Option<u32>,
    /// Transversion count against the reference.
    pub trans_v: Option<u32>,
    /// Divergence carried over from the source pairwise alignment.
    pub src_div: Option<f64>,
}

/// Half-open column range occupied by non-padding bytes; `(0, 0)` when the
/// row is all padding.
fn col_bounds(seq: &[u8]) -> (usize, usize) {
    match seq.iter().position(|&b| !seq::is_pad(b)) {
        Some(s) => (s, seq.iter().rposition(|&b| !seq::is_pad(b)).unwrap() + 1),
        None => (0, 0),
    }
}

impl SequenceRow {
    /// Build from a name and a gapped row, deriving the column bounds from
    /// padding.
    pub fn new(name: impl Into<String>, seq: Vec<u8>) -> Self {
        let (col_start, col_end) = col_bounds(&seq);
        SequenceRow {
            name: name.into(),
            seq,
            col_start,
            col_end,
            span: None,
            orient: Strand::Plus,
            lf_seq: None,
            rf_seq: None,
            gc_background: None,
            div: None,
            kdiv: None,
            trans_i: None,
            trans_v: None,
            src_div: None,
        }
    }

    /// Residue count, excluding gaps and padding.
    pub fn ungapped_len(&self) -> usize {
        seq::ungapped_len(&self.seq)
    }

    /// Recompute `col_start`/`col_end` after the row has been edited.
    pub fn refresh_bounds(&mut self) {
        let (s, e) = col_bounds(&self.seq);
        self.col_start = s;
        self.col_end = e;
    }
}

// ── Multiple alignment ────────────────────────────────────────────────────────

/// A multiple alignment: `sequences[0]` is the reference, the rest are instances.
///
/// Every row shares [`width`](Self::width).
#[derive(Debug, Clone, Default)]
pub struct MultiAlign {
    /// Reference first, then instances.
    pub sequences: Vec<SequenceRow>,
    width: usize,
    /// Column ranges flagged as low-quality, in gapped coordinates.
    pub low_quality_blocks: Vec<(usize, usize)>,
}

impl MultiAlign {
    pub fn new() -> Self {
        Self::default()
    }

    /// Assemble from a reference row and instance rows.
    ///
    /// Returns an error unless every row is exactly as wide as the reference —
    /// the invariant the rest of the type depends on.
    pub fn from_sequences(reference: SequenceRow, instances: Vec<SequenceRow>) -> Result<Self> {
        let width = reference.seq.len();
        for row in &instances {
            if row.seq.len() != width {
                return Err(Error::Msa(format!(
                    "row {:?} is {} columns wide; reference is {width}",
                    row.name,
                    row.seq.len()
                )));
            }
        }
        let mut sequences = Vec::with_capacity(1 + instances.len());
        sequences.push(reference);
        sequences.extend(instances);
        Ok(MultiAlign { sequences, width, low_quality_blocks: Vec::new() })
    }

    // ── Dimensions ────────────────────────────────────────────────────────

    /// Number of columns.
    pub fn width(&self) -> usize {
        self.width
    }

    /// Number of instances, excluding the reference.
    /// True when nothing alignable is left: no rows, zero width, or a reference
    /// with no instances beside it.
    ///
    /// [`slice_columns`](Self::slice_columns) and [`trim`](Self::trim) can reach
    /// this from ordinary input — trimming 3+3 off a reference with 6 ungapped
    /// bases empties it — and they deliberately keep the reference row rather
    /// than returning no rows at all. A zero-width reference is a state the
    /// caller can test for and warn about; an empty row list has lost the fact
    /// that there ever was a reference. Check this after any trim or slice that
    /// might be aggressive.
    pub fn is_degenerate(&self) -> bool {
        self.sequences.is_empty() || self.width == 0 || self.num_instances() == 0
    }

    pub fn num_instances(&self) -> usize {
        self.sequences.len().saturating_sub(1)
    }

    pub fn reference(&self) -> Option<&SequenceRow> {
        self.sequences.first()
    }

    pub fn reference_seq(&self) -> Option<&[u8]> {
        self.sequences.first().map(|s| s.seq.as_slice())
    }

    /// Instance at 0-based index `i` among instances (0 = first non-reference).
    pub fn instance(&self, i: usize) -> Option<&SequenceRow> {
        self.sequences.get(i + 1)
    }

    // ── Profile ───────────────────────────────────────────────────────────

    /// Per-column symbol counts, indexed by an arbitrary caller-supplied
    /// alphabet mapping.
    ///
    /// `index_of` maps an upper-cased byte to a column of the returned profile;
    /// return `None` to exclude a symbol.  Passing
    /// [`SubstMatrix::index_of`](crate::matrix::SubstMatrix::index_of) counts
    /// into the matrix's own alphabet, which is what the consensus callers want.
    ///
    /// Padding is never counted — a row that does not reach a column
    /// contributes nothing to it, which is what makes per-column coverage
    /// meaningful.
    pub fn build_profile<F>(
        &self,
        alphabet_len: usize,
        index_of: F,
        include_reference: bool,
    ) -> Vec<Vec<u32>>
    where
        F: Fn(u8) -> Option<usize>,
    {
        let mut profile = vec![vec![0u32; alphabet_len]; self.width];
        let start = if include_reference { 0 } else { 1 };
        for row in self.sequences.iter().skip(start) {
            for (col, &b) in row.seq.iter().enumerate() {
                if seq::is_pad(b) {
                    continue;
                }
                if let Some(i) = index_of(b.to_ascii_uppercase()) {
                    if i < alphabet_len {
                        profile[col][i] += 1;
                    }
                }
            }
        }
        profile
    }

    /// Number of instances spanning each column (padding excluded, gaps counted:
    /// a gap means the instance is present but deleted here).
    pub fn coverage(&self) -> Vec<u32> {
        let mut cov = vec![0u32; self.width];
        for row in self.sequences.iter().skip(1) {
            for (col, &b) in row.seq.iter().enumerate() {
                if !seq::is_pad(b) {
                    cov[col] += 1;
                }
            }
        }
        cov
    }

    // ── Editing ───────────────────────────────────────────────────────────

    /// Extract columns `col_start..col_end`, dropping rows that become empty.
    ///
    /// The reference is retained even if it empties, so index 0 keeps its
    /// meaning.
    pub fn slice_columns(&mut self, col_start: usize, col_end: usize) {
        let col_end = col_end.min(self.width);
        let col_start = col_start.min(col_end);
        for row in &mut self.sequences {
            row.seq = row.seq[col_start..col_end].to_vec();
            row.refresh_bounds();
        }
        self.width = col_end - col_start;
        // The reference (row 0) is kept unconditionally; instances that became
        // all-padding are dropped. Slicing to an empty range therefore leaves a
        // zero-width reference rather than no rows — see [`is_degenerate`].
        // dfam-curator's equivalent drops the reference too, which loses the
        // information that there was one; this contract was chosen over that.
        let mut first = true;
        self.sequences.retain(|row| {
            let keep = first || row.seq.iter().any(|&b| !seq::is_pad(b));
            first = false;
            keep
        });
    }

    /// Trim `left_bp` and `right_bp` *ungapped reference* positions from the ends.
    pub fn trim(&mut self, left_bp: usize, right_bp: usize) {
        if left_bp == 0 && right_bp == 0 {
            return;
        }
        let Some(ref_seq) = self.reference_seq().map(|s| s.to_vec()) else {
            return;
        };
        let left_col = ungapped_to_gapped_col(&ref_seq, left_bp);
        let right_col = if right_bp == 0 {
            self.width
        } else {
            let n = seq::ungapped_len(&ref_seq);
            ungapped_to_gapped_col(&ref_seq, n.saturating_sub(right_bp))
        };
        self.slice_columns(left_col, right_col);
    }

    /// Reverse-complement the whole alignment in place, flipping every row's
    /// strand and swapping its flanking sequences.
    pub fn reverse_complement(&mut self) {
        for row in &mut self.sequences {
            seq::revcomp_in_place(&mut row.seq);
            row.orient = row.orient.flip();
            std::mem::swap(&mut row.lf_seq, &mut row.rf_seq);
            if let Some(s) = row.lf_seq.as_mut() {
                seq::revcomp_in_place(s);
            }
            if let Some(s) = row.rf_seq.as_mut() {
                seq::revcomp_in_place(s);
            }
            row.refresh_bounds();
        }
        let w = self.width;
        for b in &mut self.low_quality_blocks {
            *b = (w.saturating_sub(b.1 + 1), w.saturating_sub(b.0 + 1));
        }
        self.low_quality_blocks.reverse();
    }
}

// ── Assembly ──────────────────────────────────────────────────────────────────

/// What to do with query bases that have no reference position.
///
/// The two `Grow*` policies place the same bases but can differ in **width**,
/// because they resolve differently where one member's insertion should share
/// columns with another's. Measured on a real family (round-5/family-4989),
/// `GrowPerSlot` opened 72 insertion columns against `GrowIncremental`'s 70 —
/// all five differing slots inside one 24 bp A/GA-rich window, where indel
/// placement is score-ambiguous anyway. Neither is provably more correct, which
/// is why both are available.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum InsertionPolicy {
    /// Discard them; the MSA width stays equal to the reference length.
    ///
    /// Matches `dfam_curator::blast::hits_to_multialign`.
    #[default]
    Drop,

    /// Open, at each inter-base slot, as many columns as the **largest single
    /// insertion** seen there across all members.
    ///
    /// Members are independent: two insertions at different reference positions
    /// never share columns, so a column implies genuine positional homology.
    /// The cost is a slightly wider alignment when the aligner places what is
    /// really one event at neighbouring positions in different members.
    GrowPerSlot,

    /// Merge each member into an accumulating reference, left to right, reusing
    /// gap columns already opened by earlier members.
    ///
    /// A faithful port of GIRI `MultipleAlignment::adjustReference`, and what
    /// `acons`/`autocons` produce. Because merging is greedy and positional, a
    /// later member's insertion can be packed into columns an earlier member
    /// opened nearby — narrower, and it coalesces one event that the aligner
    /// placed inconsistently, but it can also imply homology that is not there.
    ///
    /// **Order-sensitive by construction.** Feed members in a stable order.
    GrowIncremental,
}

impl InsertionPolicy {
    /// Former name of [`InsertionPolicy::GrowPerSlot`].
    #[deprecated(note = "renamed to GrowPerSlot now that GrowIncremental exists")]
    #[allow(non_upper_case_globals)]
    pub const GrowReference: InsertionPolicy = InsertionPolicy::GrowPerSlot;
}

/// One instance to be placed into the alignment.
///
/// The two gapped strings must be equal-length and describe the same pairwise
/// alignment — exactly what [`Alignment::gapped`](crate::align::Alignment::gapped)
/// returns, with the reference as the subject.
#[derive(Debug, Clone)]
pub struct MsaMember<'a> {
    pub name: &'a str,
    /// Gapped instance (query) row.
    pub gapped_query: &'a [u8],
    /// Gapped reference (subject) row.
    pub gapped_reference: &'a [u8],
    /// 0-based offset into the ungapped reference where this alignment starts.
    pub ref_start: usize,
    /// Forward-strand source coordinates of the instance, if known.
    pub span: Option<Span>,
    pub orient: Strand,
}

/// Merge pairwise alignments against a shared reference into one [`MultiAlign`].
///
/// `reference` is the **ungapped** reference sequence.  Under
/// [`InsertionPolicy::Drop`] the result is exactly `reference.len()` columns
/// wide; under [`InsertionPolicy::GrowPerSlot`] it is wider by the sum of the
/// largest insertion seen at each inter-base position.
pub fn assemble_msa(
    reference: &[u8],
    ref_name: &str,
    members: &[MsaMember<'_>],
    policy: InsertionPolicy,
) -> Result<MultiAlign> {
    for m in members {
        if m.gapped_query.len() != m.gapped_reference.len() {
            return Err(Error::Msa(format!(
                "member {:?}: gapped rows differ in length ({} vs {})",
                m.name,
                m.gapped_query.len(),
                m.gapped_reference.len()
            )));
        }
        let consumed = seq::ungapped_len(m.gapped_reference);
        if m.ref_start + consumed > reference.len() {
            return Err(Error::Msa(format!(
                "member {:?}: consumes reference {}..{} but the reference is {} long",
                m.name,
                m.ref_start,
                m.ref_start + consumed,
                reference.len()
            )));
        }
    }

    match policy {
        InsertionPolicy::Drop => assemble_drop(reference, ref_name, members),
        InsertionPolicy::GrowPerSlot => assemble_grow(reference, ref_name, members),
        InsertionPolicy::GrowIncremental => assemble_incremental(reference, ref_name, members),
    }
}

/// GIRI's incremental star merge — a port of `MultipleAlignment::
/// adjustReference` plus `adjustAlignment`.
///
/// Maintains a gapped reference that grows as members are merged in. For each
/// member, its gapped-reference row is walked against the accumulated
/// reference:
///
/// | accumulated | member's ref row | action |
/// |---|---|---|
/// | gap | gap | reuse the column; member takes its query base |
/// | gap | base | column belongs to an earlier member; this member gets a gap |
/// | base | gap | open a new column here, widening every row placed so far |
/// | base | base | ordinary aligned column |
///
/// Reusing a column in row 1 means a later member's insertion is packed into it,
/// which is exactly where this diverges from [`InsertionPolicy::GrowPerSlot`].
fn assemble_incremental(
    reference: &[u8],
    ref_name: &str,
    members: &[MsaMember<'_>],
) -> Result<MultiAlign> {
    // Accumulated gapped reference; starts ungapped.
    let mut acc: Vec<u8> = reference.to_vec();
    // Rows already merged, in accumulated-reference coordinates.
    let mut rows: Vec<(String, Vec<u8>)> = Vec::with_capacity(members.len());

    for m in members {
        // Walk to the ref_start-th *base* of the accumulated reference; gap
        // columns opened by earlier members do not count.
        let mut ai = 0usize;
        let mut seen = 0usize;
        while ai < acc.len() && seen < m.ref_start {
            if !seq::is_gap(acc[ai]) {
                seen += 1;
            }
            ai += 1;
        }
        // Leading padding, in current coordinates.
        let mut row: Vec<u8> = vec![PAD; ai];

        let mut qi = 0usize; // index into the member's gapped rows
        while qi < m.gapped_reference.len() {
            if ai >= acc.len() {
                break;
            }
            let gr = m.gapped_reference[qi];
            let gq = m.gapped_query[qi];
            let acc_gap = seq::is_gap(acc[ai]);
            let gr_gap = seq::is_gap(gr);

            if acc_gap && !gr_gap {
                // A column an earlier member opened; this member spans it but
                // contributes nothing.  The C++ advances only its reference
                // pointer here, leaving the member's row untouched.
                row.push(GAP);
                ai += 1;
                continue;
            }
            if !acc_gap && gr_gap {
                // The member inserts here and no column exists yet — open one,
                // widening every row already placed.
                acc.insert(ai, GAP);
                for (_, r) in rows.iter_mut() {
                    if ai <= r.len() {
                        // A row shorter than `ai` ends before this column.
                        let fill = if ai < r.len() { GAP } else { PAD };
                        r.insert(ai, fill);
                    }
                }
            }
            // Remaining cases advance both: an ordinary aligned column, or —
            // crucially — both sides gapped, where the C++'s `*sptr1 != *sptr2`
            // test is false and the existing column is *reused*.  Opening a
            // fresh column here instead is what makes the alignment too wide.
            row.push(gq);
            ai += 1;
            qi += 1;
        }
        rows.push((m.name.to_string(), row));
    }

    let width = acc.len();
    // Pad every row out to the final width.
    for (_, r) in rows.iter_mut() {
        if r.len() < width {
            r.resize(width, PAD);
        }
        r.truncate(width);
    }

    let ref_row = SequenceRow::new(ref_name, acc);
    let instances: Vec<SequenceRow> = rows
        .into_iter()
        .zip(members)
        .map(|((name, seq), m)| {
            let mut r = SequenceRow::new(name, seq);
            r.span = m.span;
            r.orient = m.orient;
            r
        })
        .collect();
    MultiAlign::from_sequences(ref_row, instances)
}

fn assemble_drop(
    reference: &[u8],
    ref_name: &str,
    members: &[MsaMember<'_>],
) -> Result<MultiAlign> {
    let width = reference.len();
    let ref_row = SequenceRow::new(ref_name, reference.to_vec());

    let mut rows = Vec::with_capacity(members.len());
    for m in members {
        let mut row = vec![PAD; width];
        let mut col = m.ref_start;
        for (&q, &r) in m.gapped_query.iter().zip(m.gapped_reference) {
            if seq::is_gap(r) {
                continue; // insertion — dropped
            }
            if col >= width {
                break;
            }
            row[col] = q;
            col += 1;
        }
        rows.push(finish_row(m, row));
    }
    MultiAlign::from_sequences(ref_row, rows)
}

fn assemble_grow(
    reference: &[u8],
    ref_name: &str,
    members: &[MsaMember<'_>],
) -> Result<MultiAlign> {
    let ref_len = reference.len();

    // Pass 1 — the widest insertion at each of the ref_len+1 inter-base slots.
    // Slot k sits immediately before reference base k.
    let mut max_insert = vec![0usize; ref_len + 1];
    for m in members {
        let mut ref_pos = m.ref_start;
        let mut run = 0usize;
        for (&_q, &r) in m.gapped_query.iter().zip(m.gapped_reference) {
            if seq::is_gap(r) {
                run += 1;
            } else {
                if run > 0 {
                    max_insert[ref_pos] = max_insert[ref_pos].max(run);
                    run = 0;
                }
                ref_pos += 1;
            }
        }
        if run > 0 {
            max_insert[ref_pos] = max_insert[ref_pos].max(run);
        }
    }

    // Pass 2 — lay out columns: [slot 0][base 0][slot 1][base 1] ... [slot n].
    let width = ref_len + max_insert.iter().sum::<usize>();
    let mut slot_col = vec![0usize; ref_len + 1]; // first column of each slot
    let mut base_col = vec![0usize; ref_len];
    {
        let mut col = 0usize;
        for i in 0..ref_len {
            slot_col[i] = col;
            col += max_insert[i];
            base_col[i] = col;
            col += 1;
        }
        slot_col[ref_len] = col;
        debug_assert_eq!(col + max_insert[ref_len], width);
    }

    // Reference row: bases at their columns, gaps through every insertion slot.
    let mut ref_seq = vec![GAP; width];
    for (i, &b) in reference.iter().enumerate() {
        ref_seq[base_col[i]] = b;
    }
    let ref_row = SequenceRow::new(ref_name, ref_seq);

    // Pass 3 — place each member.
    let mut rows = Vec::with_capacity(members.len());
    for m in members {
        let mut row = vec![PAD; width];
        let ref_consumed = seq::ungapped_len(m.gapped_reference);
        let ref_end = m.ref_start + ref_consumed;

        // Interior slots are part of this member's span, so unused insertion
        // columns there are gaps, not padding.  The slots at the two ends are
        // outside it and stay padded.
        for slot in (m.ref_start + 1)..ref_end {
            let from = slot_col[slot];
            row[from..from + max_insert[slot]].fill(GAP);
        }

        let mut ref_pos = m.ref_start;
        let mut ins_used = 0usize;
        for (&q, &r) in m.gapped_query.iter().zip(m.gapped_reference) {
            if seq::is_gap(r) {
                // Left-align inserted bases within their slot.
                let c = slot_col[ref_pos] + ins_used;
                if c < slot_col[ref_pos] + max_insert[ref_pos] {
                    row[c] = q;
                }
                ins_used += 1;
                continue;
            }
            ins_used = 0;
            if ref_pos < ref_len {
                row[base_col[ref_pos]] = q;
            }
            ref_pos += 1;
        }
        rows.push(finish_row(m, row));
    }

    MultiAlign::from_sequences(ref_row, rows)
}

fn finish_row(m: &MsaMember<'_>, row: Vec<u8>) -> SequenceRow {
    let mut r = SequenceRow::new(m.name, row);
    r.span = m.span;
    r.orient = m.orient;
    r
}

/// Column index in `seq` holding the `bp`-th ungapped residue, or `seq.len()`.
fn ungapped_to_gapped_col(seq_bytes: &[u8], bp: usize) -> usize {
    let mut n = 0usize;
    for (col, &b) in seq_bytes.iter().enumerate() {
        if !seq::is_structural(b) {
            if n == bp {
                return col;
            }
            n += 1;
        }
    }
    seq_bytes.len()
}

#[cfg(test)]
mod tests {
    /// Slicing to an empty range must leave a testable zero-width reference,
    /// not an empty alignment. Reachable from ordinary `trim` input, so the
    /// contract is pinned rather than left to whatever `retain` happens to do.
    #[test]
    fn an_emptied_alignment_keeps_a_zero_width_reference() {
        let r = SequenceRow::new("ref", b"-CGWN--T-Y---".to_vec());
        let i1 = SequenceRow::new("inst", b"  SCM-CKA-G  ".to_vec());
        let mut m = MultiAlign::from_sequences(r, vec![i1]).unwrap();

        m.slice_columns(4, 4); // empty range

        assert_eq!(m.sequences.len(), 1, "the reference row must survive");
        assert_eq!(m.width(), 0);
        assert!(m.sequences[0].seq.is_empty());
        assert!(m.is_degenerate(), "and the state must be detectable");
    }

    /// `trim` reaches the same state when the requested trim exceeds what the
    /// reference has to give: 6 ungapped bases, 3 off each end.
    #[test]
    fn over_trimming_is_degenerate_not_empty() {
        let r = SequenceRow::new("ref", b"-CGWN--T-Y---".to_vec());
        let i1 = SequenceRow::new("inst", b"  SCM-CKA-G  ".to_vec());
        let mut m = MultiAlign::from_sequences(r, vec![i1]).unwrap();

        m.trim(3, 3);

        assert!(m.is_degenerate());
        assert!(!m.sequences.is_empty(), "a caller can still see there was a reference");
    }

    use super::*;

    fn member<'a>(
        name: &'a str,
        q: &'a [u8],
        r: &'a [u8],
        ref_start: usize,
    ) -> MsaMember<'a> {
        MsaMember {
            name,
            gapped_query: q,
            gapped_reference: r,
            ref_start,
            span: Some(Span::new(0, seq::ungapped_len(q) as u64).unwrap()),
            orient: Strand::Plus,
        }
    }

    #[test]
    fn rows_must_match_the_reference_width() {
        let r = SequenceRow::new("ref", b"ACGT".to_vec());
        let bad = SequenceRow::new("x", b"ACG".to_vec());
        assert!(MultiAlign::from_sequences(r, vec![bad]).is_err());
    }

    #[test]
    fn drop_policy_keeps_the_reference_width() {
        // The instance carries a 2 bp insertion; under Drop it vanishes.
        let reference = b"ACGT";
        let m = member("q1", b"ACTTGT", b"AC--GT", 0);
        let msa = assemble_msa(reference, "ref", &[m], InsertionPolicy::Drop).unwrap();
        assert_eq!(msa.width(), 4);
        assert_eq!(msa.reference_seq().unwrap(), b"ACGT");
        assert_eq!(&msa.instance(0).unwrap().seq, b"ACGT");
    }

    #[test]
    fn grow_policy_opens_columns_for_the_insertion() {
        let reference = b"ACGT";
        let m = member("q1", b"ACTTGT", b"AC--GT", 0);
        let msa = assemble_msa(reference, "ref", &[m], InsertionPolicy::GrowPerSlot).unwrap();
        assert_eq!(msa.width(), 6);
        // The reference gains two gap columns between C and G.
        assert_eq!(msa.reference_seq().unwrap(), b"AC--GT");
        assert_eq!(&msa.instance(0).unwrap().seq, b"ACTTGT");
    }

    #[test]
    fn grow_policy_shares_one_slot_across_members_of_different_widths() {
        // Two instances insert at the same point: one 1 bp, one 3 bp.
        // The slot must be sized for the larger and the smaller padded with gaps.
        let reference = b"ACGT";
        let m1 = member("wide", b"ACTTTGT", b"AC---GT", 0);
        let m2 = member("narrow", b"ACAGT", b"AC-GT", 0);
        let msa =
            assemble_msa(reference, "ref", &[m1, m2], InsertionPolicy::GrowPerSlot).unwrap();

        assert_eq!(msa.width(), 7);
        assert_eq!(msa.reference_seq().unwrap(), b"AC---GT");
        assert_eq!(&msa.instance(0).unwrap().seq, b"ACTTTGT");
        // The 1 bp insertion is left-aligned in the 3-wide slot.
        assert_eq!(&msa.instance(1).unwrap().seq, b"ACA--GT");
    }

    #[test]
    fn partial_members_are_padded_not_gapped_outside_their_span() {
        // Instance covers reference positions 1..3 only.
        let reference = b"ACGTA";
        let m = member("part", b"CG", b"CG", 1);
        let msa = assemble_msa(reference, "ref", &[m], InsertionPolicy::GrowPerSlot).unwrap();
        let row = &msa.instance(0).unwrap().seq;
        assert_eq!(row, b" CG  ");
        assert_eq!(msa.instance(0).unwrap().col_start, 1);
        assert_eq!(msa.instance(0).unwrap().col_end, 3);
    }

    #[test]
    fn deletions_become_gaps_not_padding() {
        let reference = b"ACGT";
        let m = member("del", b"A-GT", b"ACGT", 0);
        let msa = assemble_msa(reference, "ref", &[m], InsertionPolicy::GrowPerSlot).unwrap();
        assert_eq!(&msa.instance(0).unwrap().seq, b"A-GT");
    }

    #[test]
    fn coverage_counts_gaps_but_not_padding() {
        let reference = b"ACGT";
        let full = member("full", b"ACGT", b"ACGT", 0);
        let deleted = member("del", b"A--T", b"ACGT", 0);
        let partial = member("part", b"AC", b"AC", 0);
        let msa = assemble_msa(
            reference,
            "ref",
            &[full, deleted, partial],
            InsertionPolicy::Drop,
        )
        .unwrap();
        // Columns 0-1 covered by all three; 2-3 by the two full-span rows.
        assert_eq!(msa.coverage(), vec![3, 3, 2, 2]);
    }

    #[test]
    fn assembly_rejects_a_member_running_past_the_reference() {
        let reference = b"ACGT";
        let m = member("over", b"ACGTAC", b"ACGTAC", 2);
        assert!(assemble_msa(reference, "ref", &[m], InsertionPolicy::Drop).is_err());
    }

    #[test]
    fn assembly_rejects_ragged_member_rows() {
        let reference = b"ACGT";
        let m = MsaMember {
            name: "bad",
            gapped_query: b"ACG",
            gapped_reference: b"ACGT",
            ref_start: 0,
            span: Some(Span::new(0, 3).unwrap()),
            orient: Strand::Plus,
        };
        assert!(assemble_msa(reference, "ref", &[m], InsertionPolicy::Drop).is_err());
    }

    #[test]
    fn profile_ignores_padding_and_can_include_the_reference() {
        let idx = |b: u8| match b {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        };
        let reference = b"AC";
        let m = member("x", b"AC", b"AC", 0);
        let msa = assemble_msa(reference, "ref", &[m], InsertionPolicy::Drop).unwrap();

        let without = msa.build_profile(4, idx, false);
        assert_eq!(without[0][0], 1); // one A from the instance

        let with = msa.build_profile(4, idx, true);
        assert_eq!(with[0][0], 2); // instance + reference
    }

    #[test]
    fn slice_columns_drops_emptied_instances_but_keeps_the_reference() {
        let reference = b"ACGT";
        let left = member("left", b"AC", b"AC", 0);
        let right = member("right", b"GT", b"GT", 2);
        let mut msa =
            assemble_msa(reference, "ref", &[left, right], InsertionPolicy::Drop).unwrap();
        msa.slice_columns(0, 2);
        assert_eq!(msa.width(), 2);
        assert_eq!(msa.num_instances(), 1);
        assert_eq!(msa.instance(0).unwrap().name, "left");
        assert_eq!(msa.reference_seq().unwrap(), b"AC");
    }

    #[test]
    fn reverse_complement_flips_rows_and_strands() {
        let reference = b"ACGT";
        let m = member("x", b"ACGT", b"ACGT", 0);
        let mut msa = assemble_msa(reference, "ref", &[m], InsertionPolicy::Drop).unwrap();
        msa.reverse_complement();
        assert_eq!(msa.reference_seq().unwrap(), b"ACGT"); // ACGT is its own RC
        assert_eq!(msa.instance(0).unwrap().orient, Strand::Minus);
    }

    #[test]
    fn trim_works_in_ungapped_reference_units() {
        let reference = b"ACGTAC";
        let m = member("x", b"ACGTAC", b"ACGTAC", 0);
        let mut msa = assemble_msa(reference, "ref", &[m], InsertionPolicy::Drop).unwrap();
        msa.trim(1, 2);
        assert_eq!(msa.reference_seq().unwrap(), b"CGT");
    }

    /// Whatever the merge strategy, every member's residues must survive.
    #[test]
    fn both_grow_policies_preserve_every_members_bases() {
        let reference = b"ACGTACGT";
        let m1 = member("a", b"ACGTTTACGT", b"ACGT--ACGT", 0);
        let m2 = member("b", b"ACGTAACGT", b"ACGT-ACGT", 0);
        let m3 = member("c", b"ACGTACGGGT", b"ACGTACG--T", 0);

        for policy in [InsertionPolicy::GrowPerSlot, InsertionPolicy::GrowIncremental] {
            let msa = assemble_msa(
                reference,
                "ref",
                &[m1.clone(), m2.clone(), m3.clone()],
                policy,
            )
            .unwrap();
            assert_eq!(msa.num_instances(), 3, "{policy:?}");
            for (i, expect) in [b"ACGTTTACGT".len(), b"ACGTAACGT".len(), b"ACGTACGGGT".len()]
                .iter()
                .enumerate()
            {
                assert_eq!(
                    msa.instance(i).unwrap().ungapped_len(),
                    *expect,
                    "{policy:?}: instance {i} lost bases"
                );
            }
            // The reference itself never loses residues either.
            assert_eq!(
                seq::ungapped_len(msa.reference_seq().unwrap()),
                reference.len(),
                "{policy:?}"
            );
        }
    }

    /// Two members inserting the same amount at the same slot need only one set
    /// of columns under either policy.
    #[test]
    fn a_shared_insertion_point_costs_one_set_of_columns() {
        let reference = b"ACGTACGT";
        let m1 = member("a", b"ACGTTTACGT", b"ACGT--ACGT", 0);
        let m2 = member("b", b"ACGTGGACGT", b"ACGT--ACGT", 0);

        for policy in [InsertionPolicy::GrowPerSlot, InsertionPolicy::GrowIncremental] {
            let msa =
                assemble_msa(reference, "ref", &[m1.clone(), m2.clone()], policy).unwrap();
            assert_eq!(msa.width(), 10, "{policy:?}: 8 reference + 2 insertion columns");
        }
    }

    /// The incremental merge never opens more columns than the per-slot maximum:
    /// it can only reuse columns an earlier member already opened.
    #[test]
    fn incremental_is_never_wider_than_per_slot() {
        let reference = b"ACGTACGTACGT";
        let cases: &[(&[u8], &[u8], usize)] = &[
            (b"ACGTTTACGTACGT", b"ACGT--ACGTACGT", 0),
            (b"ACGTACGAAAATACGT", b"ACGTACG----TACGT", 0),
            (b"ACGTACGTAAACGT", b"ACGTACGT--ACGT", 0),
            (b"CGTACGTACG", b"CGTACGTACG", 1),
        ];
        let members: Vec<MsaMember<'_>> = cases
            .iter()
            .enumerate()
            .map(|(i, (q, r, start))| MsaMember {
                name: match i {
                    0 => "m0",
                    1 => "m1",
                    2 => "m2",
                    _ => "m3",
                },
                gapped_query: q,
                gapped_reference: r,
                ref_start: *start,
                span: Some(Span::new(0, seq::ungapped_len(q) as u64).unwrap()),
                orient: Strand::Plus,
            })
            .collect();

        let per = assemble_msa(reference, "ref", &members, InsertionPolicy::GrowPerSlot)
            .unwrap();
        let inc = assemble_msa(reference, "ref", &members, InsertionPolicy::GrowIncremental)
            .unwrap();
        assert!(
            inc.width() <= per.width(),
            "incremental {} should not exceed per-slot {}",
            inc.width(),
            per.width()
        );
        // Both must still hold every base.
        for i in 0..members.len() {
            assert_eq!(
                inc.instance(i).unwrap().ungapped_len(),
                per.instance(i).unwrap().ungapped_len(),
                "instance {i}"
            );
        }
    }

    #[test]
    fn grow_and_drop_agree_when_there_are_no_insertions() {
        let reference = b"ACGTACGT";
        let m1 = member("a", b"ACGTACGT", b"ACGTACGT", 0);
        let m2 = member("b", b"GTAC", b"GTAC", 2);
        let grow = assemble_msa(
            reference, "ref",
            &[m1.clone(), m2.clone()],
            InsertionPolicy::GrowPerSlot,
        )
        .unwrap();
        let drop = assemble_msa(reference, "ref", &[m1, m2], InsertionPolicy::Drop).unwrap();
        assert_eq!(grow.width(), drop.width());
        for i in 0..grow.sequences.len() {
            assert_eq!(grow.sequences[i].seq, drop.sequences[i].seq);
        }
    }
}
