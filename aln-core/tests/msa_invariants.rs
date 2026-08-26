//! Property tests for [`MultiAlign`]: randomised alignments, invariants checked
//! after every mutating operation.
//!
//! Salvaged from `msa-difftest`, a differential harness that compared this type
//! against `dfam-curator`'s independently derived one. That comparison is
//! finished — the two agreed everywhere except `reverse_complement`, where
//! curator failed to refresh row bounds, and the retention contract for an
//! emptied alignment. The *generator* and the self-consistency check are what
//! generalise, and they need no second implementation to be useful.
//!
//! The invariant that found the bug: a row's `col_start`/`col_end` must be the first and
//! last non-padding column of its own `seq`. Any operation that moves columns
//! has to re-establish it.

use aln_core::msa::{MultiAlign, SequenceRow};
use aln_core::Strand;
use aln_coord::Span;

/// Deterministic xorshift so a failure can be replayed from its seed.
struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
    }
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n.max(1) as u64) as usize
    }
}

/// Build an alignment with the awkward shapes: ragged ends, interior gaps,
/// IUPAC codes, all-gap and all-padding rows, and single-row alignments.
///
/// Deliberately limited to what a `MultiAlign` legally holds — residues, `-`,
/// and ` `. `.` (Stockholm) and `<`/`>` (GIRI) are normalised away at the input
/// boundaries, so generating them tests undefined behaviour rather than a
/// contract.
fn make(rng: &mut Rng) -> MultiAlign {
    let width = 12 + rng.below(40);
    let nrows = 1 + rng.below(6);
    let alpha = b"ACGTACGTACGTRYKMSWN";
    let mut rows = Vec::new();
    for r in 0..nrows {
        let lpad = if r == 0 { 0 } else { rng.below(width / 3) };
        let rpad = if r == 0 { 0 } else { rng.below(width / 3) };
        let shape = rng.below(12);
        let mut seq = vec![b' '; width];
        if shape != 0 {
            for cell in seq.iter_mut().take(width - rpad).skip(lpad) {
                *cell = match rng.below(6) {
                    0 | 1 => b'-',
                    _ => alpha[rng.below(alpha.len())],
                };
            }
        }
        let mut row = SequenceRow::new(format!("row{r}"), seq);
        let s0 = rng.below(500) as u64;
        row.span = Some(Span::new(s0, s0 + 50).unwrap());
        row.orient = if r > 0 && rng.below(4) == 0 { Strand::Minus } else { Strand::Plus };
        rows.push(row);
    }
    let (first, rest) = rows.split_first().unwrap();
    MultiAlign::from_sequences(first.clone(), rest.to_vec()).expect("generated a valid alignment")
}

/// `col_start`/`col_end` must describe the row's own `seq`, and every row must be as
/// wide as the alignment.
fn check(m: &MultiAlign, op: &str, seed: u64) {
    let w = m.width();
    for (i, row) in m.sequences.iter().enumerate() {
        assert_eq!(
            row.seq.len(),
            w,
            "{op} (seed {seed}): row {i} is {} wide, alignment is {w}",
            row.seq.len()
        );
        let pad = |b: u8| matches!(b, b' ' | b'<' | b'>');
        let (want_s, want_e) = match row.seq.iter().position(|&b| !pad(b)) {
            Some(s) => (s, row.seq.iter().rposition(|&b| !pad(b)).unwrap() + 1),
            None => (0, 0),
        };
        assert_eq!(
            (row.col_start, row.col_end),
            (want_s, want_e),
            "{op} (seed {seed}): row {i} bounds {}..{} contradict its own seq {:?}",
            row.col_start,
            row.col_end,
            String::from_utf8_lossy(&row.seq)
        );
    }
}

#[test]
fn mutating_operations_preserve_row_invariants() {
    for n in 0..3000u64 {
        let seed = 20_260_811 + n;
        let mut rng = Rng::new(seed);

        let m = make(&mut rng);
        check(&m, "construct", seed);

        // reverse_complement moves every column; bounds must be re-derived.
        // Curator's equivalent skipped this and produced rows contradicting
        // their own sequence in ~4,000 cases.
        let mut a = make(&mut rng);
        a.reverse_complement();
        check(&a, "reverse_complement", seed);

        // Degenerate ranges are reachable from ordinary trim input, so they are
        // generated on purpose rather than avoided.
        let w = m.width();
        let lo = rng.below(w);
        let hi = match rng.below(6) {
            0 => lo,
            1 => w,
            _ => (lo + 1 + rng.below(w.saturating_sub(lo).max(1))).min(w),
        };
        let mut b = make(&mut rng);
        b.slice_columns(lo, hi);
        check(&b, "slice_columns", seed);

        let mut c = make(&mut rng);
        c.trim(rng.below(4), rng.below(4));
        check(&c, "trim", seed);
        assert!(
            !c.sequences.is_empty(),
            "trim (seed {seed}) must keep the reference row even when it empties it"
        );
    }
}

#[test]
fn reverse_complement_is_an_involution() {
    for n in 0..500u64 {
        let seed = 900_000 + n;
        let mut rng = Rng::new(seed);
        let before = make(&mut rng);
        let mut after = before.clone();
        after.reverse_complement();
        after.reverse_complement();
        for (x, y) in before.sequences.iter().zip(&after.sequences) {
            assert_eq!(x.seq, y.seq, "double reverse-complement changed the sequence (seed {seed})");
            assert_eq!(x.orient, y.orient, "orientation not restored (seed {seed})");
            assert_eq!((x.col_start, x.col_end), (y.col_start, y.col_end), "bounds not restored (seed {seed})");
        }
    }
}
