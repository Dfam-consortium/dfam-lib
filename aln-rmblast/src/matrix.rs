//! Converting an [`aln_core::SubstMatrix`] into rmblast's `ScoreMatrix`.
//!
//! # The two matrices are transposes of one another
//!
//! RepeatMasker ships the same matrices twice, in two layouts:
//!
//! ```text
//! Matrices/crossmatch/14p35g.matrix      Matrices/ncbi/nt/14p35g.matrix
//! FREQS A 0.325 C 0.175 ...              # FREQS A 0.325 C 0.175 ...
//!   A   R   G   C  ...                       A   R   G   C  ...
//!   8   0 -10 -18  ...                   A   8   3  -7 -17  ...
//!   3   3  12 -17  ...                   R   0   3   2 -16  ...
//!  -7   2  12 -16  ...                   G -10  12  12 -16  ...
//! ```
//!
//! Read the first *column* of the crossmatch file and it is the first *row* of
//! the NCBI file.  They are transposes, because the two lineages index scores in
//! opposite orders:
//!
//! | | row | column |
//! |---|---|---|
//! | crossmatch / `SearchResult.pm::rescoreAlignment` / [`aln_core`] | subject (consensus) | query (genomic) |
//! | NCBI / rmblast `ScoreMatrix::score(q, s)` | query | subject |
//!
//! Both therefore end up with the *same* scoring semantics when each is fed its
//! own file.  But converting between the in-memory forms means transposing, and
//! getting that wrong is silent: the matrices are asymmetric, so a swapped
//! lookup changes scores without any error.
//! [`orientation_is_preserved`](#tests) pins it with a deliberately lopsided
//! matrix.
//!
//! # Two further quirks of rmblast's parser
//!
//! * Frequencies are only recognised on a **comment** line (`# FREQS ...`).  A
//!   bare `FREQS` line — which is what the crossmatch files carry — would fall
//!   through to the score-row parser and corrupt the matrix, so the rendering
//!   here always comments it.
//! * `X` maps to BLASTNA index 15, the gap slot, whose row and column are then
//!   overwritten with `INT_MIN/2`.  That is NCBI's behaviour and is preserved;
//!   it does not matter in practice because sequence encoding sends `X` to `N`.

use aln_core::SubstMatrix;
use aln_engine::EngineError;
use rmblast_lib::matrix::ScoreMatrix;

use crate::NAME;

/// Render `subst` in NCBI layout — i.e. transposed — and parse it back through
/// rmblast's own reader.
///
/// Going via the text format rather than poking `ScoreMatrix`'s array directly
/// means rmblast applies its own symbol mapping, gap-slot handling and lambda
/// estimation, exactly as it would for a file supplied on the command line.
pub fn to_rmblast(subst: &SubstMatrix) -> Result<ScoreMatrix, EngineError> {
    let text = render_ncbi(subst);
    ScoreMatrix::from_reader(subst.name(), std::io::Cursor::new(text.into_bytes())).map_err(|e| {
        EngineError::backend(NAME, format!("rmblast rejected the converted matrix: {e}"))
    })
}

/// Render an [`aln_core::SubstMatrix`] in NCBI nucleotide-matrix layout.
///
/// The output is the transpose of [`SubstMatrix`]'s own `Display`, carries row
/// labels, and comments the `FREQS` line.
pub fn render_ncbi(subst: &SubstMatrix) -> String {
    let alphabet = subst.alphabet();
    let mut out = String::new();

    if subst.has_freqs() {
        out.push_str("# FREQS");
        for (i, &c) in alphabet.iter().enumerate() {
            let f = subst.freqs()[i];
            if f > 0.0 {
                out.push_str(&format!(" {} {}", c as char, f));
            }
        }
        out.push('\n');
    }

    out.push_str("   ");
    for &c in alphabet {
        out.push_str(&format!("{:>4}", c as char));
    }
    out.push('\n');

    // NCBI row = query, column = subject.  aln-core is score_idx(subject, query).
    for (qi, &qc) in alphabet.iter().enumerate() {
        out.push_str(&format!("{:<3}", qc as char));
        for si in 0..alphabet.len() {
            out.push_str(&format!("{:>4}", subst.score_idx(si, qi)));
        }
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use rmblast_lib::encoding::IUPAC_TO_BLASTNA;

    const CROSSMATCH_14P35G: &str = "\
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

    fn blastna(c: u8) -> u8 {
        IUPAC_TO_BLASTNA[c as usize]
    }

    /// The load-bearing test: an asymmetric pair must survive the transpose.
    ///
    /// In 14p35g, `score(subject=G, query=A)` is -7 and `score(subject=A,
    /// query=G)` is -10.  After conversion rmblast must report the same values
    /// under its own `score(query, subject)` order.
    #[test]
    fn orientation_is_preserved() {
        let subst = SubstMatrix::parse(CROSSMATCH_14P35G).unwrap();
        assert_eq!(subst.score(b'G', b'A'), Some(-7), "aln-core: subject G, query A");
        assert_eq!(subst.score(b'A', b'G'), Some(-10), "aln-core: subject A, query G");

        let rm = to_rmblast(&subst).unwrap();
        // rmblast is score(query, subject).
        assert_eq!(rm.score(blastna(b'A'), blastna(b'G')), -7, "rmblast: query A, subject G");
        assert_eq!(rm.score(blastna(b'G'), blastna(b'A')), -10, "rmblast: query G, subject A");
    }

    #[test]
    fn every_unambiguous_cell_round_trips() {
        let subst = SubstMatrix::parse(CROSSMATCH_14P35G).unwrap();
        let rm = to_rmblast(&subst).unwrap();
        for &s in b"ACGT" {
            for &q in b"ACGT" {
                assert_eq!(
                    rm.score(blastna(q), blastna(s)),
                    subst.score(s, q).unwrap(),
                    "subject {} query {}",
                    s as char,
                    q as char
                );
            }
        }
    }

    #[test]
    fn iupac_codes_survive_too() {
        let subst = SubstMatrix::parse(CROSSMATCH_14P35G).unwrap();
        let rm = to_rmblast(&subst).unwrap();
        // R/Y/K/M/S/W/N all have distinct rows; X is deliberately excluded
        // because rmblast maps it onto the gap slot.
        for &s in b"ARGCYTKMSWN" {
            for &q in b"ARGCYTKMSWN" {
                assert_eq!(
                    rm.score(blastna(q), blastna(s)),
                    subst.score(s, q).unwrap(),
                    "subject {} query {}",
                    s as char,
                    q as char
                );
            }
        }
    }

    #[test]
    fn freqs_are_commented_so_rmblast_reads_them() {
        let subst = SubstMatrix::parse(CROSSMATCH_14P35G).unwrap();
        let text = render_ncbi(&subst);
        assert!(
            text.starts_with("# FREQS"),
            "a bare FREQS line would be parsed as a score row: {}",
            text.lines().next().unwrap()
        );
        let rm = to_rmblast(&subst).unwrap();
        assert!(rm.lambda > 0.0, "lambda needs frequencies; got {}", rm.lambda);
    }

    #[test]
    fn rendered_layout_matches_the_shipped_ncbi_file() {
        // The first data row of Matrices/ncbi/nt/14p35g.matrix, verbatim:
        //   A   8   3  -7 -17 -19 -21 -14  -4 -12  -6  -1 -30
        let subst = SubstMatrix::parse(CROSSMATCH_14P35G).unwrap();
        let text = render_ncbi(&subst);
        let row_a = text
            .lines()
            .find(|l| l.starts_with("A "))
            .expect("no A row rendered");
        let vals: Vec<i32> = row_a[3..]
            .split_whitespace()
            .map(|t| t.parse().unwrap())
            .collect();
        assert_eq!(vals, vec![8, 3, -7, -17, -19, -21, -14, -4, -12, -6, -1, -30]);
    }
}
