//! Sequence primitives: strand, IUPAC alphabet handling, and the gap/padding
//! character conventions shared by the whole stack.
//!
//! # Character conventions
//!
//! There are three distinct "not a base" characters and conflating them is the
//! most common source of bugs when moving data between the Dfam and GIRI
//! lineages:
//!
//! | char | name | meaning |
//! |------|------|---------|
//! | `-`  | [`GAP`] | an interior indel — the sequence *is* present here, aligned to nothing |
//! | ` `  | [`PAD`] | flanking padding — the sequence does not extend this far; missing data |
//! | `.`  | — | accepted on input as a synonym for `-` (Stockholm) |
//!
//! Dfam (`dfam-curator`) uses `' '` for padding.  GIRI `ScoreMatrix.hpp` instead
//! uses `'<'` (LPADDING) and `'>'` (RPADDING).  This crate canonicalises on the
//! Dfam form; [`from_giri_padding`] and [`to_giri_padding`] convert at the
//! boundary.

/// Interior gap character (GIRI `GAPCHAR`).
pub const GAP: u8 = b'-';

/// Flanking padding character — the sequence is absent, not deleted.
pub const PAD: u8 = b' ';

/// GIRI `LPADDING`.
pub const GIRI_LPAD: u8 = b'<';

/// GIRI `RPADDING`.
pub const GIRI_RPAD: u8 = b'>';

// ── Strand ────────────────────────────────────────────────────────────────────

/// Strand of a sequence relative to its source.
///
/// Named to match `dfam_curator::Orientation` so the two can be mapped 1:1.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Strand {
    #[default]
    Plus,
    Minus,
}

impl Strand {
    /// RepeatMasker `.out` / crossmatch orientation character: `+` or `C`.
    pub fn as_rm_char(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => 'C',
        }
    }

    /// BLAST `sstrand` token.
    pub fn as_blast_str(self) -> &'static str {
        match self {
            Strand::Plus => "plus",
            Strand::Minus => "minus",
        }
    }

    pub fn is_minus(self) -> bool {
        matches!(self, Strand::Minus)
    }

    #[must_use]
    pub fn flip(self) -> Strand {
        match self {
            Strand::Plus => Strand::Minus,
            Strand::Minus => Strand::Plus,
        }
    }
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
        })
    }
}

// ── IUPAC tables ──────────────────────────────────────────────────────────────

/// IUPAC complement, case-preserving, with `-`/` `/`<`/`>` passed through.
///
/// Unknown bytes map to `N` (upper) — matching the behaviour of both the Perl
/// and the `dfam-curator` implementations.
static COMPLEMENT: [u8; 256] = build_complement();

const fn build_complement() -> [u8; 256] {
    let mut t = [b'N'; 256];
    const PAIRS: &[(u8, u8)] = &[
        (b'A', b'T'), (b'T', b'A'), (b'G', b'C'), (b'C', b'G'),
        (b'R', b'Y'), (b'Y', b'R'), (b'K', b'M'), (b'M', b'K'),
        (b'S', b'S'), (b'W', b'W'), (b'B', b'V'), (b'V', b'B'),
        (b'D', b'H'), (b'H', b'D'), (b'N', b'N'), (b'X', b'X'),
        (b'U', b'A'), (b'Z', b'Z'),
    ];
    let mut i = 0;
    while i < PAIRS.len() {
        let (up, comp) = PAIRS[i];
        t[up as usize] = comp;
        // Lower-case input yields lower-case output (soft-masking survives).
        t[(up + 32) as usize] = comp + 32;
        i += 1;
    }
    // Non-base characters are structural and must survive complementation.
    t[GAP as usize] = GAP;
    t[PAD as usize] = PAD;
    t[b'.' as usize] = b'.';
    // GIRI padding swaps handedness under reverse-complement.
    t[GIRI_LPAD as usize] = GIRI_RPAD;
    t[GIRI_RPAD as usize] = GIRI_LPAD;
    t
}

/// Complement a single IUPAC byte, preserving case.
#[inline]
pub fn complement(b: u8) -> u8 {
    COMPLEMENT[b as usize]
}

/// Reverse-complement into a new buffer.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Reverse-complement in place.
pub fn revcomp_in_place(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = complement(*b);
    }
}

/// True for the four unambiguous DNA bases (either case).
///
/// This is the predicate behind RepeatMasker's `%wellCharacterizedBases` table:
/// a pairing is well-characterised iff *both* sides satisfy it.
#[inline]
pub fn is_acgt(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// True for `-`, `.`, ` `, `<`, `>` — anything that is not a residue.
#[inline]
pub fn is_structural(b: u8) -> bool {
    matches!(b, GAP | PAD | b'.' | GIRI_LPAD | GIRI_RPAD)
}

/// True for the interior gap character (accepts the Stockholm `.` synonym).
#[inline]
pub fn is_gap(b: u8) -> bool {
    b == GAP || b == b'.'
}

/// True for flanking padding in either the Dfam or GIRI convention.
#[inline]
pub fn is_pad(b: u8) -> bool {
    matches!(b, PAD | GIRI_LPAD | GIRI_RPAD)
}

// ── Padding conventions ───────────────────────────────────────────────────────

/// Rewrite GIRI `<`/`>` padding to the canonical `' '`.
pub fn from_giri_padding(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        if *b == GIRI_LPAD || *b == GIRI_RPAD {
            *b = PAD;
        }
    }
}

/// Rewrite canonical `' '` padding to GIRI `<`/`>`, deciding handedness by
/// position: padding before the first residue is `<`, after the last is `>`.
///
/// Interior spaces (which should not occur in well-formed data) are treated as
/// right padding, matching GIRI's read path.
pub fn to_giri_padding(seq: &mut [u8]) {
    let first = seq.iter().position(|&b| !is_pad(b));
    let Some(first) = first else {
        // Entirely padding — GIRI has no basis to pick a side; use left.
        seq.iter_mut().for_each(|b| *b = GIRI_LPAD);
        return;
    };
    for b in seq[..first].iter_mut() {
        *b = GIRI_LPAD;
    }
    for b in seq[first..].iter_mut() {
        if is_pad(*b) {
            *b = GIRI_RPAD;
        }
    }
}

/// Strip gaps and padding, yielding the ungapped residue sequence.
pub fn ungap(seq: &[u8]) -> Vec<u8> {
    seq.iter().copied().filter(|&b| !is_structural(b)).collect()
}

/// Count of residues (excludes gaps and padding).
pub fn ungapped_len(seq: &[u8]) -> usize {
    seq.iter().filter(|&&b| !is_structural(b)).count()
}

// ── Named sequence ────────────────────────────────────────────────────────────

/// An ungapped, named sequence — the input unit for the aligners.
///
/// Equivalent to GIRI `NamedSequence`, minus the sub-sequence machinery: this
/// crate carries coordinates on [`crate::align::Alignment`] instead of baking
/// them into the sequence type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub seq: Vec<u8>,
}

impl Sequence {
    pub fn new(name: impl Into<String>, seq: impl Into<Vec<u8>>) -> Self {
        Sequence { name: name.into(), seq: seq.into() }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Reverse-complement, leaving the name untouched.
    #[must_use]
    pub fn revcomp(&self) -> Sequence {
        Sequence { name: self.name.clone(), seq: revcomp(&self.seq) }
    }

    /// GC fraction over unambiguous bases; `None` if there are none.
    pub fn gc_fraction(&self) -> Option<f64> {
        let mut gc = 0usize;
        let mut n = 0usize;
        for &b in &self.seq {
            match b.to_ascii_uppercase() {
                b'G' | b'C' => { gc += 1; n += 1; }
                b'A' | b'T' => { n += 1; }
                _ => {}
            }
        }
        (n > 0).then(|| gc as f64 / n as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_preserves_case_and_structure() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'a'), b't');
        assert_eq!(complement(GAP), GAP);
        assert_eq!(complement(PAD), PAD);
        assert_eq!(complement(b'?'), b'N');
    }

    #[test]
    fn revcomp_roundtrips() {
        let s = b"ACGTNRYacgt";
        let once = revcomp(s);
        let twice = revcomp(&once);
        assert_eq!(&twice, s);
    }

    #[test]
    fn revcomp_keeps_gaps_in_place() {
        // A gap is a column, not a base — it must survive reversal as a gap.
        assert_eq!(revcomp(b"AC-GT"), b"AC-GT");
    }

    #[test]
    fn giri_padding_swaps_hands_under_revcomp() {
        assert_eq!(complement(GIRI_LPAD), GIRI_RPAD);
        assert_eq!(revcomp(b"<<ACGT>"), b"<ACGT>>");
    }

    #[test]
    fn padding_conversion_roundtrip() {
        let mut s = b"  ACGT  ".to_vec();
        to_giri_padding(&mut s);
        assert_eq!(&s, b"<<ACGT>>");
        from_giri_padding(&mut s);
        assert_eq!(&s, b"  ACGT  ");
    }

    #[test]
    fn ungap_drops_gaps_and_padding_but_not_n() {
        assert_eq!(ungap(b" AC-GN "), b"ACGN");
        assert_eq!(ungapped_len(b" AC-GN "), 4);
    }

    #[test]
    fn gc_fraction_ignores_ambiguity() {
        let s = Sequence::new("x", b"GCATNNNN".to_vec());
        assert_eq!(s.gc_fraction(), Some(0.5));
        assert_eq!(Sequence::new("y", b"NNNN".to_vec()).gc_fraction(), None);
    }
}
