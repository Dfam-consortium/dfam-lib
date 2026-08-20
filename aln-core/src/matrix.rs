//! Substitution matrices.
//!
//! One type covers the three formats in play:
//!
//! * **RepeatMasker / crossmatch** (`Matrices/crossmatch/*.matrix`) — an optional
//!   `FREQS` line, a single-character column header, then square rows.
//! * **GIRI** (`ScoreMatrix.hpp`) — the same layout plus an optional
//!   `GAP <ini> <ext>` line, and rows that may carry a leading row label.
//! * **NCBI** — same shape with row labels always present.
//!
//! GIRI's parser is already a superset of crossmatch's (it skips `FREQS` lines),
//! so a single parser reads all three.  Unlike GIRI we *retain* the frequencies,
//! because RepeatMasker needs them for lambda and for Phil Green's
//! complexity-adjusted scoring.
//!
//! # Orientation is not symmetric — read this before calling [`SubstMatrix::score`]
//!
//! Arian Smit's matrices were built with the assignment
//!
//! ```text
//!     query = genomic / derived state          (columns)
//!   subject = consensus / ancestral state      (rows)
//! ```
//!
//! and are **not symmetric**.  `SearchResult.pm::rescoreAlignment` looks scores
//! up as `matrix[subj_base][query_base]`.  [`SubstMatrix::score`] takes its
//! arguments in that order and is documented accordingly; getting it backwards
//! silently changes scores rather than failing.

use std::fmt;

use crate::error::{Error, Result};

/// Sentinel for "this byte is not in the alphabet".
const NO_INDEX: i16 = -1;

/// A square substitution matrix over a single-character alphabet.
#[derive(Debug, Clone)]
pub struct SubstMatrix {
    name: String,
    alphabet: Vec<u8>,
    /// ASCII byte -> alphabet index, or [`NO_INDEX`].  Case-insensitive.
    char2index: [i16; 256],
    /// Row-major `size * size` scores; `scores[subj * size + query]`.
    scores: Vec<i32>,
    /// Background frequencies parallel to `alphabet`; 0.0 where unspecified.
    freqs: Vec<f64>,
    /// Karlin-Altschul lambda, derived from `scores` + `freqs` when the
    /// frequencies form a distribution.  `None` when there is no `FREQS` line.
    lambda: Option<f64>,
    /// GIRI `GAP <ini> <ext>` line, if present.  Typically negative.
    gap_open: Option<i32>,
    gap_extend: Option<i32>,
}

impl SubstMatrix {
    // ── Accessors ─────────────────────────────────────────────────────────────

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = name.into();
    }

    /// The alphabet in matrix order (e.g. `A R G C Y T K M S W N X`).
    pub fn alphabet(&self) -> &[u8] {
        &self.alphabet
    }

    pub fn size(&self) -> usize {
        self.alphabet.len()
    }

    /// Alphabet index for an ASCII byte, case-insensitively; `None` if absent.
    #[inline]
    pub fn index_of(&self, b: u8) -> Option<usize> {
        let i = self.char2index[b as usize];
        (i != NO_INDEX).then_some(i as usize)
    }

    /// Score for a **subject (consensus) base against a query (genomic) base**.
    ///
    /// Argument order matters — see the module docs.  Returns `None` if either
    /// byte is outside the alphabet (gaps included: gaps are scored by the gap
    /// penalties, not by the matrix).
    #[inline]
    pub fn score(&self, subj: u8, query: u8) -> Option<i32> {
        let s = self.index_of(subj)?;
        let q = self.index_of(query)?;
        Some(self.scores[s * self.alphabet.len() + q])
    }

    /// Score by pre-resolved indices.  Panics if out of range.
    #[inline]
    pub fn score_idx(&self, subj: usize, query: usize) -> i32 {
        self.scores[subj * self.alphabet.len() + query]
    }

    /// Background frequencies parallel to [`alphabet`](Self::alphabet).
    /// All-zero when the source had no `FREQS` line.
    pub fn freqs(&self) -> &[f64] {
        &self.freqs
    }

    pub fn has_freqs(&self) -> bool {
        self.freqs.iter().any(|&f| f > 0.0)
    }

    /// Karlin-Altschul lambda; `None` without frequencies.
    ///
    /// Required by complexity-adjusted scoring
    /// ([`crate::stats::rescore`] with `complexity_adjust`).
    pub fn lambda(&self) -> Option<f64> {
        self.lambda
    }

    /// Gap-open penalty from a GIRI `GAP` line, if the source carried one.
    pub fn gap_open(&self) -> Option<i32> {
        self.gap_open
    }

    /// Gap-extend penalty from a GIRI `GAP` line, if the source carried one.
    pub fn gap_extend(&self) -> Option<i32> {
        self.gap_extend
    }

    /// GIRI's gap-vs-gap score, derived as `-gapExt / 2` under C++ integer
    /// division (so `gapExt = -5` yields `+2`).
    ///
    /// This participates in GIRI consensus scoring but never in alignment.
    pub fn gap_match_score(&self) -> Option<i32> {
        self.gap_extend.map(|gx| -gx / 2)
    }

    /// Minimum and maximum entries — handy for choosing a SIMD lane width.
    pub fn score_range(&self) -> (i32, i32) {
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for &v in &self.scores {
            lo = lo.min(v);
            hi = hi.max(v);
        }
        (lo, hi)
    }

    // ── Construction ──────────────────────────────────────────────────────────

    /// Build from an explicit alphabet and a row-major `size * size` score list.
    pub fn from_parts(
        name: impl Into<String>,
        alphabet: &[u8],
        scores: Vec<i32>,
        freqs: Option<Vec<f64>>,
    ) -> Result<Self> {
        let size = alphabet.len();
        if size == 0 {
            return Err(Error::Matrix("empty alphabet".into()));
        }
        if scores.len() != size * size {
            return Err(Error::Matrix(format!(
                "expected {} scores for a {}-symbol alphabet, got {}",
                size * size,
                size,
                scores.len()
            )));
        }
        let freqs = match freqs {
            Some(f) if f.len() != size => {
                return Err(Error::Matrix(format!(
                    "expected {} frequencies, got {}",
                    size,
                    f.len()
                )))
            }
            Some(f) => f,
            None => vec![0.0; size],
        };

        let mut char2index = [NO_INDEX; 256];
        for (i, &c) in alphabet.iter().enumerate() {
            char2index[c.to_ascii_uppercase() as usize] = i as i16;
            char2index[c.to_ascii_lowercase() as usize] = i as i16;
        }

        let mut m = SubstMatrix {
            name: name.into(),
            alphabet: alphabet.to_vec(),
            char2index,
            scores,
            freqs,
            lambda: None,
            gap_open: None,
            gap_extend: None,
        };
        m.lambda = m.calculate_lambda()?;
        Ok(m)
    }

    /// Parse a matrix from text in crossmatch / GIRI / NCBI layout.
    ///
    /// Recognised lines, in any order before the body:
    /// * `FREQS A 0.325 C 0.175 ...` — background frequencies (crossmatch).
    /// * `GAP -25 -5` — gap open / extend (GIRI).
    /// * `#` comments and blank lines are ignored.
    ///
    /// The first remaining line is the column header (single characters); every
    /// line after it is a score row, with an optional leading row label.
    pub fn parse(text: &str) -> Result<Self> {
        let mut freqs_by_char: Vec<(u8, f64)> = Vec::new();
        let mut gap_open = None;
        let mut gap_extend = None;
        let mut alphabet: Vec<u8> = Vec::new();
        let mut rows: Vec<Vec<i32>> = Vec::new();
        // Rows may be labelled and therefore out of order; track placement.
        let mut row_order: Vec<usize> = Vec::new();
        let mut next_unlabelled = 0usize;

        for raw in text.lines() {
            let line = raw.split('#').next().unwrap_or("").trim();
            if line.is_empty() {
                continue;
            }

            if let Some(rest) = line.strip_prefix("FREQS") {
                let toks: Vec<&str> = rest.split_whitespace().collect();
                for pair in toks.chunks(2) {
                    if pair.len() != 2 {
                        return Err(Error::Matrix(format!("malformed FREQS line: {raw}")));
                    }
                    let c = pair[0].as_bytes()[0];
                    let f: f64 = pair[1]
                        .parse()
                        .map_err(|_| Error::Matrix(format!("bad frequency {:?}", pair[1])))?;
                    freqs_by_char.push((c, f));
                }
                continue;
            }

            // GIRI matches "GAP" anywhere in the line; anchoring at the start is
            // stricter and still accepts every matrix we have seen.
            if let Some(rest) = line.strip_prefix("GAP") {
                let toks: Vec<&str> = rest.split_whitespace().collect();
                if toks.len() < 2 {
                    return Err(Error::Matrix(format!("malformed GAP line: {raw}")));
                }
                gap_open = Some(
                    toks[0]
                        .parse()
                        .map_err(|_| Error::Matrix(format!("bad gap open {:?}", toks[0])))?,
                );
                gap_extend = Some(
                    toks[1]
                        .parse()
                        .map_err(|_| Error::Matrix(format!("bad gap extend {:?}", toks[1])))?,
                );
                continue;
            }

            if alphabet.is_empty() {
                for tok in line.split_whitespace() {
                    if tok.len() != 1 {
                        return Err(Error::Matrix(format!(
                            "expected single-character column header, got {tok:?}"
                        )));
                    }
                    alphabet.push(tok.as_bytes()[0]);
                }
                if alphabet.is_empty() {
                    return Err(Error::Matrix("empty column header".into()));
                }
                continue;
            }

            // Body row, optionally prefixed with its own label.
            let mut toks: Vec<&str> = line.split_whitespace().collect();
            let mut target = next_unlabelled;
            let first_is_label = toks[0].len() == 1 && toks[0].as_bytes()[0].is_ascii_alphabetic();
            if first_is_label {
                let c = toks[0].as_bytes()[0].to_ascii_uppercase();
                target = alphabet
                    .iter()
                    .position(|&a| a.to_ascii_uppercase() == c)
                    .ok_or_else(|| {
                        Error::Matrix(format!("row label {:?} absent from column header", toks[0]))
                    })?;
                toks.remove(0);
            }
            next_unlabelled = target + 1;

            let vals: Result<Vec<i32>> = toks
                .iter()
                .map(|t| {
                    // Scores are integral in practice but GIRI reads them as
                    // doubles and truncates, so accept "12.0" too.
                    t.parse::<i32>().or_else(|_| {
                        t.parse::<f64>()
                            .map(|v| v as i32)
                            .map_err(|_| Error::Matrix(format!("bad score {t:?}")))
                    })
                })
                .collect();
            let vals = vals?;
            if vals.len() != alphabet.len() {
                return Err(Error::Matrix(format!(
                    "row {} has {} values, expected {}",
                    rows.len(),
                    vals.len(),
                    alphabet.len()
                )));
            }
            rows.push(vals);
            row_order.push(target);
        }

        if alphabet.is_empty() {
            return Err(Error::Matrix("no column header found".into()));
        }
        if rows.len() != alphabet.len() {
            return Err(Error::Matrix(format!(
                "expected {} rows, found {}",
                alphabet.len(),
                rows.len()
            )));
        }

        let size = alphabet.len();
        let mut scores = vec![0i32; size * size];
        for (row, &target) in rows.iter().zip(&row_order) {
            if target >= size {
                return Err(Error::Matrix(format!("row index {target} out of range")));
            }
            scores[target * size..(target + 1) * size].copy_from_slice(row);
        }

        let freqs = if freqs_by_char.is_empty() {
            None
        } else {
            let mut f = vec![0.0; size];
            for (c, v) in freqs_by_char {
                let idx = alphabet
                    .iter()
                    .position(|&a| a.eq_ignore_ascii_case(&c))
                    .ok_or_else(|| {
                        Error::Matrix(format!("FREQS symbol {:?} absent from alphabet", c as char))
                    })?;
                f[idx] = v;
            }
            Some(f)
        };

        let mut m = Self::from_parts(String::new(), &alphabet, scores, freqs)?;
        m.gap_open = gap_open;
        m.gap_extend = gap_extend;
        Ok(m)
    }

    /// Parse from a file, taking the file stem as the matrix name.
    pub fn from_file(path: impl AsRef<std::path::Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::Matrix(format!("cannot read {}: {e}", path.display())))?;
        let mut m = Self::parse(&text)?;
        if let Some(stem) = path.file_stem().and_then(|s| s.to_str()) {
            m.set_name(stem);
        }
        Ok(m)
    }

    // ── Lambda ────────────────────────────────────────────────────────────────

    /// Port of `Matrix.pm::_calculateLambda`.
    ///
    /// Solves `sum_ij f_i f_j exp(lambda * S_ij) = 1` by doubling to bracket the
    /// root then bisecting to a tolerance of 1e-5 — the same schedule as the
    /// Perl, so results agree to the printed precision.
    ///
    /// Returns `Ok(None)` when there are no frequencies to work from.
    fn calculate_lambda(&self) -> Result<Option<f64>> {
        if !self.has_freqs() {
            return Ok(None);
        }

        let check: f64 = {
            let mut c = 0.0;
            for &fi in &self.freqs {
                for &fj in &self.freqs {
                    if fi > 0.0 && fj > 0.0 {
                        c += fi * fj;
                    }
                }
            }
            c
        };
        if !(0.999..=1.001).contains(&check) {
            return Err(Error::Matrix(format!(
                "matrix frequencies sum to {check:.6} (squared); expected ~1.0"
            )));
        }

        let sum_at = |lambda: f64| -> f64 {
            let mut sum = 0.0;
            for (i, &fi) in self.freqs.iter().enumerate() {
                if fi <= 0.0 {
                    continue;
                }
                for (j, &fj) in self.freqs.iter().enumerate() {
                    if fj <= 0.0 {
                        continue;
                    }
                    sum += fi * fj * (lambda * self.score_idx(i, j) as f64).exp();
                }
            }
            sum
        };

        let mut lambda_lower = 0.0;
        let mut lambda = 0.5;
        // Double until the sum crosses 1.0 from below.
        loop {
            let sum = sum_at(lambda);
            if sum >= 1.0 {
                break;
            }
            lambda_lower = lambda;
            lambda *= 2.0;
            if lambda > 1e6 {
                return Err(Error::Matrix(
                    "lambda failed to bracket — matrix has no positive expected score".into(),
                ));
            }
        }
        let mut lambda_upper = lambda;

        while lambda_upper - lambda_lower > 1e-5 {
            lambda = (lambda_lower + lambda_upper) / 2.0;
            if sum_at(lambda) >= 1.0 {
                lambda_upper = lambda;
            } else {
                lambda_lower = lambda;
            }
        }
        Ok(Some(lambda))
    }
}

impl fmt::Display for SubstMatrix {
    /// Render in crossmatch layout, round-trippable through [`SubstMatrix::parse`].
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.has_freqs() {
            write!(f, "FREQS")?;
            for (i, &c) in self.alphabet.iter().enumerate() {
                if self.freqs[i] > 0.0 {
                    write!(f, " {} {}", c as char, self.freqs[i])?;
                }
            }
            writeln!(f)?;
        }
        if let (Some(go), Some(ge)) = (self.gap_open, self.gap_extend) {
            writeln!(f, "GAP {go} {ge}")?;
        }
        for &c in &self.alphabet {
            write!(f, "{:>4}", c as char)?;
        }
        writeln!(f)?;
        let size = self.alphabet.len();
        for i in 0..size {
            for j in 0..size {
                write!(f, "{:>4}", self.score_idx(i, j))?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// RepeatMasker `Matrices/crossmatch/14p35g.matrix`, verbatim.
    const M14P35G: &str = "\
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

    #[test]
    fn parses_crossmatch_matrix() {
        let m = SubstMatrix::parse(M14P35G).unwrap();
        assert_eq!(m.size(), 12);
        assert_eq!(m.alphabet(), b"ARGCYTKMSWNX");
        // Diagonal A/A and T/T.
        assert_eq!(m.score(b'A', b'A'), Some(8));
        assert_eq!(m.score(b'T', b'T'), Some(8));
        assert_eq!(m.score(b'C', b'C'), Some(12));
        assert_eq!(m.score(b'G', b'G'), Some(12));
    }

    #[test]
    fn score_is_subject_then_query_and_asymmetric() {
        let m = SubstMatrix::parse(M14P35G).unwrap();
        // Row G (index 2), column A (index 0) = -7.
        // Row A (index 0), column G (index 2) = -10.
        // If these were ever equal the test would not detect a swapped lookup.
        assert_eq!(m.score(b'G', b'A'), Some(-7));
        assert_eq!(m.score(b'A', b'G'), Some(-10));
    }

    #[test]
    fn lookup_is_case_insensitive_and_rejects_gaps() {
        let m = SubstMatrix::parse(M14P35G).unwrap();
        assert_eq!(m.score(b'a', b'a'), Some(8));
        assert_eq!(m.score(b'-', b'A'), None);
        assert_eq!(m.score(b'A', b' '), None);
    }

    #[test]
    fn freqs_land_on_the_right_symbols() {
        let m = SubstMatrix::parse(M14P35G).unwrap();
        let f = m.freqs();
        // Alphabet order is A R G C Y T ... so A=0, G=2, C=3, T=5.
        assert_eq!(f[0], 0.325);
        assert_eq!(f[2], 0.175);
        assert_eq!(f[3], 0.175);
        assert_eq!(f[5], 0.325);
        // R, Y, K, M, S, W, N, X carry no frequency.
        assert_eq!(f[1], 0.0);
    }

    #[test]
    fn lambda_is_derived_and_plausible() {
        let m = SubstMatrix::parse(M14P35G).unwrap();
        let lambda = m.lambda().expect("14p35g has FREQS");
        // Bracket generously; the exact value is pinned by the round-trip test
        // against Matrix.pm in tests/matrix_lambda.rs.
        assert!(lambda > 0.0 && lambda < 1.0, "lambda = {lambda}");
        // Verify it actually solves sum f_i f_j exp(lambda S_ij) = 1.
        let mut sum = 0.0;
        for (i, &fi) in m.freqs().iter().enumerate() {
            for (j, &fj) in m.freqs().iter().enumerate() {
                if fi > 0.0 && fj > 0.0 {
                    sum += fi * fj * (lambda * m.score_idx(i, j) as f64).exp();
                }
            }
        }
        assert!((sum - 1.0).abs() < 1e-3, "sum = {sum}");
    }

    #[test]
    fn no_freqs_means_no_lambda() {
        let text = "  A   C\n  1  -1\n -1   1\n";
        let m = SubstMatrix::parse(text).unwrap();
        assert!(m.lambda().is_none());
        assert!(!m.has_freqs());
    }

    #[test]
    fn parses_giri_gap_line_and_derives_gap_match_score() {
        let text = "GAP -25 -5\n  A   C\n  1  -1\n -1   1\n";
        let m = SubstMatrix::parse(text).unwrap();
        assert_eq!(m.gap_open(), Some(-25));
        assert_eq!(m.gap_extend(), Some(-5));
        // GIRI: gapMatchScore = -gapExt / 2 under integer division.
        assert_eq!(m.gap_match_score(), Some(2));
    }

    #[test]
    fn parses_labelled_rows_out_of_order() {
        // NCBI-style row labels; deliberately shuffled to prove placement is by
        // label rather than by line position.
        let text = "  A   C\nC -1   1\nA  1  -1\n";
        let m = SubstMatrix::parse(text).unwrap();
        assert_eq!(m.score(b'A', b'A'), Some(1));
        assert_eq!(m.score(b'C', b'C'), Some(1));
        assert_eq!(m.score(b'A', b'C'), Some(-1));
    }

    #[test]
    fn display_round_trips() {
        let m = SubstMatrix::parse(M14P35G).unwrap();
        let rendered = m.to_string();
        let m2 = SubstMatrix::parse(&rendered).unwrap();
        assert_eq!(m.alphabet(), m2.alphabet());
        assert_eq!(m.freqs(), m2.freqs());
        for i in 0..m.size() {
            for j in 0..m.size() {
                assert_eq!(m.score_idx(i, j), m2.score_idx(i, j));
            }
        }
    }

    #[test]
    fn rejects_wrong_row_width() {
        let text = "  A   C\n  1  -1\n -1\n";
        assert!(SubstMatrix::parse(text).is_err());
    }

    #[test]
    fn score_range_reports_extremes() {
        let m = SubstMatrix::parse(M14P35G).unwrap();
        assert_eq!(m.score_range(), (-30, 12));
    }
}
